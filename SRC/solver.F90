module solver
#include <petsc/finclude/petscksp.h>
  use petscksp

! SOLVER: solver for elasto-dynamic equation
!         M.a = - K.d + F

  use problem_class
  use petsc_objects 
  use echo, only : ItInfo
  use stdio, only : IO_Abort
  use sources, only : SO_add
  use bc_gen , only : BC_apply, bc_apply_kind, &
                      bc_trans, bc_update_dfault, bc_nnodeActive_DYNFLT, &
                      bc_update_bcdv, bc_reset, bc_has_dynflt, bc_update_pre
  use fields_class, only: FIELD_SetVecFromField, FIELD_SetFieldFromVec, &
                      FIELDS_reset, FIELDS_Update_Pre

  implicit none
  private

  public :: solve, solve_GF_DYNFLT

contains

!=====================================================================

subroutine solve(pb, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp

  type(problem_type), intent(inout) :: pb  
  type(petsc_objects_type), intent(inout) :: petobj
  integer:: rank
  PetscErrorCode :: ierr
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)

  select case (pb%time%kind)
    case ('adaptive')
      call solve_adaptive(pb, petobj)
    case ('leapfrog')
      call solve_leapfrog(pb, petobj)
    case ('newmark')
      call solve_Newmark(pb, petobj)
    case ('HHT-alpha')
      call solve_HHT_alpha(pb, petobj)
    case default
      call solve_symplectic(pb, petobj)
  end select

end subroutine solve

!=====================================================================
! adaptive solver, a wrapper for dynamic and quasi-static solver
! depending on the 'isDynamic' flag
!
subroutine solve_adaptive(pb, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp
  use time_evol, only : TIME_broadcast
  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  logical::id1, id2, isw
  integer,parameter::IS_DYNFLT=6
  PetscErrorCode :: ierr
  double precision :: tmp(size(pb%fields%displ, 1), size(pb%fields%displ, 2))
  double precision :: tn, dt
  integer :: rank, i
  
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)
  
  id1 = pb%time%isdynamic
  tn  = pb%time%time - pb%time%dt

  if (pb%time%isDynamic) then
      ! dynamic, obtain perturbation

      pb%fields%displ = pb%fields%displ - pb%fields%dn

      if (rank==0) then
          call solve_dynamic(pb, petobj)
      end if
      
      ! add perturbation back to the background displacement
      pb%fields%displ = pb%fields%displ + pb%fields%dn
  else
      ! quasi-static
      call solve_quasi_static_petsc(pb, petobj)
      ! update background values for dynamic step
      pb%fields%dn = pb%fields%displ
  end if ! isdynamic
   
  if (rank==0) then
      if (.not. pb%time%fixdt) call update_adaptive_step(pb)
  end if

  ! broad cast constants in t
  call TIME_broadcast(pb%time)

  id2 = pb%time%isdynamic
  isw = id2 .and. (.not. id1)
  ! switch from static to dynamic
  if (isw) then
     pb%fields%veloc = 0d0
     pb%fields%veloc_pre = 0d0
  end if

end subroutine solve_adaptive

! update time step 
subroutine update_adaptive_step(pb)
  use bc_gen, only: BC_timestep
  type(problem_type), intent(inout) :: pb
  double precision :: vgmax
  vgmax = maxval(abs(pb%fields%veloc))

  ! rate-state time step following Lapusta et al., 2010
  call BC_timestep(pb%bc,pb%time, vgmax)
  
end subroutine update_adaptive_step

!
! a solver to compute the DYNFLT static Green's function
! from fault slip to fault shear w/o normal stress using petsc
!
! Algorithm:
!   1. do i = 1, number_of_active_nodes_on_dynflt 
!   2. For each active node, set slip of one node to 1 and others to zeros.
!   3. Solve global linear system to obtain solution, d, f
!   4. evaluate stresses on the fault and set values of one row of Green's function
!   5. end do !i
!

subroutine solve_GF_DYNFLT(pb, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp
  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  double precision, dimension(:,:), pointer :: d, f
  
  integer :: i, nnode

  integer :: IS_DIRNEU = 1, & 
             IS_KINFLT = 2, & 
             IS_DYNFLT = 6, &
             IS_DIRABS = 7
         
  integer :: IS_DIR    = 2, IS_NEU = 1
  PetscErrorCode :: ierr
  PetscScalar, pointer :: xx_d(:)
  PetscScalar, pointer :: xx_b(:)
  KSPConvergedReason :: reason
  logical :: print_time
  integer :: rank
  double precision :: start_time, end_time, elapse_time
  character(160):: outputString
  
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)
  print_time = rank==0 .and. mod(pb%time%it, ItInfo) == 0
  has_dynflt  = pb%has_dynflt
  
  if (rank==0 .and. (.not. has_dynflt)) then
      call IO_Abort('No DYNFLT defined! Stop!')
  end if

  if (rank==0) nnode = bc_nnodeActive_DYNFLT(pb%bc)
  call MPI_Bcast(nnode, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

  if (rank==0) then
      ! points to displacement
      d => pb%fields%displ
      f => pb%fields%accel
  end if

  do i = 1, nnode
    
    pb%time%pcg_iters = 0
    if (rank==0) then 

    ! set d to zero
    d  = 0d0 
    ! apply dirichlet boundary condition and kinematic boundary condition
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_DIR)
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRABS, IS_DIR)

    ! transform fields, apply half slip rate on side -1 and then transform back
    ! different from dynamic, update displacement as well.
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_KINFLT)
    end if
   
    if (rank==0) then
        ! change the this subroutine so that fault slip is updated
        ! on i-th active node on DYNFLT to 1 and others zero
!    call bc_update_dfault(pb%bc, d_pre, d, (v + v_pre)/2d0, dt)

    f=0d0
    ! apply newmann if there's any, used to solve dofs in the medium
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_NEU) 

    ! transform d to X*d, f to Xinv*f
    call BC_trans(pb%bc, d,  1)
    call BC_trans(pb%bc, f, -1)
    end if

    start_time = MPI_Wtime()
    call FIELD_SetVecFromField(petobj%d, d, ierr)
    call FIELD_SetVecFromField(petobj%b, f, ierr)
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time

    if (print_time) then 
    write(*, *) 'Elapse Time for setting d b from local vectors =', elapse_time, 's'
    end if

    ! scatter d to b
    start_time = MPI_Wtime()
    call VecScatterBegin(petobj%bd_fix_ctx,petobj%d,petobj%b,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(petobj%bd_fix_ctx,petobj%d,petobj%b,INSERT_VALUES,SCATTER_FORWARD,ierr)
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time

    ! solve the linear system for d
    start_time = MPI_Wtime()
    call KSPSolve(petobj%ksp, petobj%b, petobj%d, ierr)  
    CHKERRQ(ierr) 
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then 
    write(*, *) 'Elapse Time for KspSolve =', elapse_time, 's'
    end if
    
    ! return the number of iterations
    call KSPGetIterationNumber(petobj%ksp, pb%time%pcg_iters(i), ierr)

    ! return converged reason 
    call KSPGetConvergedReason(petobj%ksp, reason, ierr)
    pb%time%reasons(i) = reason

    if (reason<0 .and. rank==0) then
        write(*, *) "KSP solver diverges, with reason = ", reason
    end if

    ! copy the displacement from converged solution in petsc to d
    start_time = MPI_Wtime()
    call VecScatterBegin(petobj%d_ctx,petobj%d,petobj%d0,INSERT_VALUES,SCATTER_FORWARD,ierr);
    call VecScatterEnd(petobj%d_ctx,petobj%d,petobj%d0,INSERT_VALUES,SCATTER_FORWARD,ierr);

    if (rank==0) then
        call FIELD_SetFieldFromVec(d, petobj%d0, ierr)  
    end if 
    
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then
!    write(*, *) 'Elapse Time for copy petsc vec d to local d0 =', elapse_time, 's'
    end if

    if (rank==0) then

    ! transfer the d' to d
    call BC_trans(pb%bc, d,  -1)

    end if

    start_time = MPI_Wtime()

    ! compute stress on the fault and update one row of Green's function
    ! in BYNFLT_FAULT
    if (rank==0) then
    call compute_Fint(f, d, pb%fields%veloc, pb, .false.)

    ! bc_set_gf_iactive_dynflt(pb%bc, f, i)
!    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DYNFLT)
    end if

    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then
    write(*, *) 'Elapse Time for compute_Fint=', elapse_time, 's'
    end if

  end do !i
  
end subroutine solve_GF_DYNFLT

!=====================================================================
! a wrapper over all dynamic solver
subroutine solve_dynamic(pb, petobj)
    type(problem_type), intent(inout) :: pb
    type(petsc_objects_type), intent(inout) :: petobj
    select case (pb%time%kind_dyn)
      case ('leapfrog')
        call solve_leapfrog(pb, petobj)
      case ('newmark')
        call solve_Newmark(pb, petobj)
      case ('HHT-alpha')
        call solve_HHT_alpha(pb, petobj)
      case default
        call solve_symplectic(pb, petobj)
    end select
end subroutine solve_dynamic

! SOLVE: advance ONE time step
!        using single predictor-corrector Newmark (explicit)
!        in acceleration form
subroutine solve_Newmark(pb,petobj)

  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  
  double precision, dimension(:,:), pointer :: d,v,a,f
  double precision :: dt,beta,gamma

  d => pb%fields%displ
  v => pb%fields%veloc
  a => pb%fields%accel
  f => a ! NOTE: to save some memory

  dt = pb%time%dt
  beta = pb%time%beta
  gamma = pb%time%gamma
  
 !-- predictors
  d = d + dt*v + (0.5d0-beta)*dt*dt *a
  v = v + (1d0-gamma)*dt *a

 !-- internal forces
 ! WARNING: compute_Fint for plasticity assumes explicit displacement update
 !          which is not true when beta/=0
  call compute_Fint(f,d,v,pb)

 !-- Sources
  call SO_add(pb%src, pb%time%time, f)

 !-- Apply boundary conditions
  call BC_apply(pb%bc,pb%time,pb%fields,f)

 ! NOTE: if a source is an incident wave, it is not added during
 !       "call SO_add" but during "call BC_apply"
 !       Incident waves must be used with absorbing boundaries

 !-- Divide by mass matrix
  a = f*pb%rmass

 !--  Corrector 
  v = v + gamma*dt*a
  d = d + beta*dt*dt*a
  
end subroutine solve_Newmark

!=====================================================================
! Hilber-Hughes-Taylor alpha scheme
!
subroutine solve_HHT_alpha(pb,petobj)

  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  
  double precision, dimension(:,:), pointer :: d,v,a,f,d_alpha,v_alpha
  double precision :: t_alpha,tmp,dt,alpha,beta,gamma

  d => pb%fields%displ
  v => pb%fields%veloc
  a => pb%fields%accel
  f => a ! NOTE: to save some memory
  d_alpha => pb%fields%displ_alpha
  v_alpha => pb%fields%veloc_alpha

  dt = pb%time%dt
  alpha = pb%time%alpha
  beta = pb%time%beta
  gamma = pb%time%gamma
  
  d_alpha = d
  v_alpha = v
  d = d + dt*v + (0.5d0-beta)*dt*dt *a
  v = v + (1d0-gamma)*dt *a
  d_alpha = alpha*d + (1.d0-alpha)*d_alpha
  v_alpha = alpha*v + (1.d0-alpha)*v_alpha

  call compute_Fint(f,d_alpha,v_alpha,pb)
  t_alpha = pb%time%time +(alpha-1.d0)*dt
  call SO_add(pb%src, t_alpha, f)
  !call BC_apply(pb%bc,t_alpha,pb%fields,f)
  tmp = pb%time%time
  pb%time%time = t_alpha
  call BC_apply(pb%bc,pb%time,pb%fields,f)
  pb%time%time = tmp
  a = f*pb%rmass

  v = v + gamma*dt*a
  d = d + beta*dt*dt*a
  
end subroutine solve_HHT_alpha

!=====================================================================
!
! Second-order central difference (leap-frog)
!
! M*(v[n+1/2] - v[n-1/2])/dt = - K*d[n] - C*v[n] + B*t[n] + F(n)
! d[n+1] = d[n] + dt*v[n+1/2]
!
! also: a[n] = (v[n+1/2] - v[n-1/2])/dt

subroutine solve_leapfrog(pb,petobj)

  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  
  double precision, dimension(:,:), pointer :: d,v_mid,a,f

  d => pb%fields%displ
  v_mid => pb%fields%veloc
  a => pb%fields%accel
  f => a
  
  d = d + pb%time%dt * v_mid
  ! Maybe replace by K * d
  call compute_Fint(f,d,v_mid,pb) 
  call SO_add(pb%src, pb%time%time, f)
  call BC_apply(pb%bc,pb%time,pb%fields,f)

  a = pb%rmass * f
  v_mid = v_mid + pb%time%dt * a 
  
end subroutine solve_leapfrog

!=====================================================================
!
! Symplectic schemes
! WARNING: no boundary conditions implemented yet
! WARNING: use this only for pure elasticity
!
subroutine solve_symplectic(pb,petobj)

  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  
  double precision, dimension(:,:), pointer :: d,v,a,f
  double precision, dimension(:), pointer :: coa,cob
  double precision :: dt,t
  integer :: k

  d => pb%fields%displ
  v => pb%fields%veloc
  a => pb%fields%accel
  f => pb%fields%accel

  coa => pb%time%a
  cob => pb%time%b
  dt = pb%time%dt
  t = pb%time%time - dt
  
  do k = 1,pb%time%nstages
    d = d + dt*coa(k) * v
    call compute_Fint(f,d,v,pb) 
    t = t + dt*coa(k)
    call SO_add(pb%src, t, f)
    a = pb%rmass * f
    v = v + dt*cob(k) * a 
  enddo
  d = d + dt*coa(pb%time%nstages+1) * v
  
end subroutine solve_symplectic

!=====================================================================
! SOLVE: advance ONE time step
!        using quasi-static algorithm defined in Kaneko, et. al. 2011
!        in acceleration form
!
! Solve nonlinear RS BC and off-fault elasticity in two iterations.
!
! Currently, this solver only works for RS frictional faults and elasticity
!
! The elasticity is solved with (Jacobi-)preconditioned conjugate gradient
!
! The state and slip velocity of the fault is updated in an explicit
! but iterative way.
!
! Kinematic boundary condition must be modified
! 
! In Quasi-statics, only DIRNEU, KINFLT, DYNFLT BCs are applied.
!
! Symmetry assumption must be relaxed following Junpei Seki, 2017 
!
! bc_kind
! IS_EMPTY  = 0, &
! IS_DIRNEU = 1, &
! IS_KINFLT = 2, &
! IS_ABSORB = 3, &
! IS_PERIOD = 4, &
! IS_LISFLT = 5, &
! IS_DYNFLT = 6

subroutine solve_quasi_static_petsc(pb, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp
  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  double precision, dimension(:,:), pointer :: d, v, a, f
  double precision, dimension(pb%fields%npoin, pb%fields%ndof) :: d_pre, v_pre, fp
  
  ! vplate is built into the rate-state b.c.
  double precision :: dt, dtmax, tn 
  integer :: i, ndof, npoin,j, k
  logical :: has_dynflt = .false.

  integer :: IS_DIRNEU = 1, & 
             IS_KINFLT = 2, & 
             IS_DYNFLT = 6, &
             IS_DIRABS = 7
         
  integer :: IS_DIR    = 2, IS_NEU = 1
  PetscErrorCode :: ierr
  PetscScalar, pointer :: xx_d(:)
  PetscScalar, pointer :: xx_b(:)
  KSPConvergedReason :: reason
  logical :: isUpdateDMG
  logical :: print_time
  integer :: rank
  double precision :: start_time, end_time, elapse_time
  character(160):: outputString
  
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)
  print_time = rank==0 .and. mod(pb%time%it, ItInfo) == 0

  ndof  = pb%ndof
  npoin = pb%npoin

  isUpdateDMG = pb%isUpdateDMG
  has_dynflt  = pb%has_dynflt

  if (rank==0) then
      ! points to displacement
      ! save some memory by modifying in place
      d => pb%fields%displ
      v => pb%fields%veloc
      ! use acceleration and f to work with force
      a => pb%fields%accel
      f => pb%fields%accel

      dt = pb%time%dt
      tn = pb%time%time - dt! previous time step 
      
      ! use additional memory to keep the previous displacement
      ! use to update predictor for rsf faults
      d_pre = pb%fields%displ
      v_pre = pb%fields%veloc
      
      ! compute fp, which does not change over time step
      call compute_Fint_EP(fp, pb)
      
  end if
      ! reassemble the stiffness matrix and configure Ksp
      if (pb%time%isUpdateKsp) call UpdateKsp(pb, petobj)
  
  do j = 1, 2 ! for damage time step
    pb%time%pcg_iters = 0
    if (rank==0) then 
    ! update displacement to obtain a better initial guess for pcg solver
    ! update both d and pb%fields%disp
    
    d  = d_pre + dt * v_pre 

    ! apply dirichlet boundary condition and kinematic boundary condition
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_DIR)
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRABS, IS_DIR)

    ! transform fields, apply half slip rate on side -1 and then transform back
    ! different from dynamic, update displacement as well.
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_KINFLT)
    end if
   
    !mechanical balance without update damage variable
    ! start 2 passes
    do i = 1, 2
    ! update fault displacement for rsf 
    ! d(fault) = d_pre(fault) + dt * v(fault)
    ! v  = 0.5 * (vpre + v)
    !d = d_pre + (v + v_pre)/2d0 *dt

    if (rank==0) then
    call bc_update_dfault(pb%bc, d_pre, d, (v + v_pre)/2d0, dt)
    f=0d0
    
    ! apply newmann if there's any, used to solve dofs in the medium
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_NEU) 

    ! if there's damage or plastic strain, compute force from ep 
    f = f + fp

    ! transform d to X*d, f to Xinv*f
    call BC_trans(pb%bc, d,  1)
    call BC_trans(pb%bc, f, -1)
    end if

    start_time = MPI_Wtime()
    call FIELD_SetVecFromField(petobj%d, d, ierr)
    call FIELD_SetVecFromField(petobj%b, f, ierr)
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time

    if (print_time) then 
    write(*, *) 'Elapse Time for setting d b from local vectors =', elapse_time, 's'
    end if

    ! scatter d to b
    start_time = MPI_Wtime()
    call VecScatterBegin(petobj%bd_fix_ctx,petobj%d,petobj%b,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(petobj%bd_fix_ctx,petobj%d,petobj%b,INSERT_VALUES,SCATTER_FORWARD,ierr)
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time

    ! solve the linear system for d
    start_time = MPI_Wtime()
    call KSPSolve(petobj%ksp, petobj%b, petobj%d, ierr)  
    CHKERRQ(ierr) 
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then 
    write(*, *) 'Elapse Time for KspSolve =', elapse_time, 's'
    end if
    
    ! return the number of iterations
    call KSPGetIterationNumber(petobj%ksp, pb%time%pcg_iters(i), ierr)

    ! return converged reason 
    call KSPGetConvergedReason(petobj%ksp, reason, ierr)
    pb%time%reasons(i) = reason

    if (reason<0 .and. rank==0) then
        write(*, *) "KSP solver diverges, with reason = ", reason
    end if

    ! copy the displacement from converged solution in petsc to d
    start_time = MPI_Wtime()
    call VecScatterBegin(petobj%d_ctx,petobj%d,petobj%d0,INSERT_VALUES,SCATTER_FORWARD,ierr);
    call VecScatterEnd(petobj%d_ctx,petobj%d,petobj%d0,INSERT_VALUES,SCATTER_FORWARD,ierr);

    if (rank==0) then
        call FIELD_SetFieldFromVec(d, petobj%d0, ierr)  
    end if 
    
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then
!    write(*, *) 'Elapse Time for copy petsc vec d to local d0 =', elapse_time, 's'
    end if

    if (rank==0) then
    ! transfer the d' to d
    call BC_trans(pb%bc, d,  -1)

    ! update velocity for the off fault dofs
    v = (d - d_pre)/dt

    ! detect if there's dynflt fault boundary
    ! if not, break the loop and exit after the pcg solver
    end if

    if (.not. has_dynflt) exit
    ! apply boundary conditions including following: 
    !    1.  compute fault traction, compute state variable 
    !    2.  update slip rate, update fault velocity in fields%veloc 
    !
    ! velocity is modified inside this subroutine

    ! recompute the force use the updated total disp
    ! do not update damage state variable

    start_time = MPI_Wtime()
    if (rank==0) then
    call compute_Fint(f, d, pb%fields%veloc, pb, .false.)

    ! update rate-state fault
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DYNFLT)
    end if
    end_time   = MPI_Wtime()
    elapse_time = end_time - start_time
    if (print_time) then
    write(*, *) 'Elapse Time for compute_Fint=', elapse_time, 's'
    end if

  enddo !i 

  ! print the coordinates of the nodes that have high slip rate
!  if (rank==0) then
!      do k = 1, pb%grid%npoin
!         if (maxval(abs(v(k,:)))>1e-2) then
!             write(*, *) "Large velocity, v = ", v(k,:), " at coord = ", pb%grid%coord(:,k)
!         end if
!      end do
!  end if

  if (.not. isUpdateDMG) exit
!
  if (rank==0) then
  ! now compute the maximum time step to update the damage variable
  call compute_dtmax(d, pb, dtmax)
  
  if (dt<=dtmax .or. j==2) then
      ! update fault state variable
      call Update_DMG(d, pb)
      ! update dt
      ! if there's no rate state fault, update the time step
      if (.not. pb%time%fixdt) then
          if (.not. has_dynflt) then
              pb%time%dt = 0.9*dtmax
          end if
      end if
      exit
  else
      if (pb%time%fixdt) then
          call IO_Abort('Solver: fixed time step too large!!')
      end if

      ! change dt, t
      dt = dtmax*0.9
      pb%time%time = tn + dt
      pb%time%dt   = dt
      ! reset the fault state variable
      call bc_reset(pb%bc, dt)
      ! redo the mechanical balance
  end if
  end if ! rank==0
  enddo !j

  if (rank==0) then
      ! declare final slip values on the fault
      call bc_update_bcdv(pb%bc, d, pb%time%time)
      a   = 0d0
  end if
  
end subroutine solve_quasi_static_petsc

! compute the maximum allowed time step
! limited by updating the internal variable of 
! bulk material

subroutine compute_dtmax(d, pb, dtmax)
  use fields_class, only : FIELD_get_elem
  use mat_gen, only : MAT_dtmax

  type(problem_type), intent(inout) :: pb
  double precision, dimension(:,:) :: d
  double precision :: dtmax, dtmax_i
  integer :: e
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: dloc
  
  dtmax_i = 0d0
  dtmax   = huge(0d0)

  do e = 1,pb%grid%nelem
    dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))
    ! compute the maximum time step allowed 
    ! to update state variable in the bulk material
    call MAT_dtmax(dloc, pb%matwrk(e), pb%grid%ngll, pb%fields%ndof, dtmax_i)
    dtmax = min(dtmax, dtmax_i)
  end do

end subroutine

! update internal variable for damage models
subroutine Update_DMG(d, pb)
  use fields_class, only : FIELD_get_elem
  use mat_gen, only : MAT_Update_DMG

  type(problem_type), intent(inout) :: pb
  double precision, dimension(:,:) :: d
  integer :: e
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: dloc
  
  do e = 1,pb%grid%nelem
    dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))
    ! compute the maximum time step allowed 
    ! to update state variable in the bulk material
    call MAT_Update_DMG(dloc, pb%matwrk(e), pb%grid%ngll, & 
                        pb%fields%ndof, pb%time%dt, pb%time%isdynamic)
  end do

end subroutine

!=====================================================================
! f = - K * d 
! update: update the internal variables
subroutine compute_Fint(f, d, v, pb, update, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp
  use fields_class, only : FIELD_get_elem, FIELD_add_elem
  use mat_gen, only : MAT_Fint
  use constants, only : TINY_XABS
  
  double precision, dimension(:,:), intent(out) :: f
  double precision, dimension(:,:), intent(in) :: d,v
  type(problem_type), intent(inout) :: pb
  logical, optional :: update
  type(petsc_objects_type), optional :: petobj
  logical :: update_in, f_use_matrix

  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: dloc,vloc,floc
  double precision :: E_ep, E_el, sg(pb%fields%ndof+1), sgp(pb%fields%ndof+1)
  integer :: e

  f = 0d0
  pb%energy%E_el = 0d0
  pb%energy%sg   = 0d0
  pb%energy%sgp  = 0d0

  if (present(update)) then
      update_in = update
  else
      update_in = .true.
  end if


  do e = 1,pb%grid%nelem
    dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))
    vloc = FIELD_get_elem(v,pb%grid%ibool(:,:,e))

    ! skip the elements with zero v and d

    !if (all(abs(dloc)<TINY_XABS) .and. all(abs(vloc)<TINY_XABS)) cycle

    call MAT_Fint(floc,dloc,vloc,pb%matpro(e), pb%matwrk(e), & 
                   pb%grid%ngll,pb%fields%ndof,pb%time%dt, &
                   pb%grid, update_in, E_ep,E_el,sg,sgp, pb%time%isdynamic)
    call FIELD_add_elem(floc,f,pb%grid%ibool(:,:,e)) ! assembly

   ! total elastic energy change
    pb%energy%E_el = pb%energy%E_el +E_el
   ! cumulated plastic energy
    pb%energy%E_ep = pb%energy%E_ep +E_ep
   ! cumulated stress glut
    pb%energy%sg = pb%energy%sg + sg
    pb%energy%sgp = pb%energy%sgp + sgp

  enddo

!DEVEL: to parallelize this loop for multi-cores (OpenMP)
!DEVEL: reorder the elements to avoid conflict during assembly (graph coloring)
!DEVEL: loop on the colors, with sync at the end of each color
! do icol=1,size(colors)
!   do k = 1,colors(icol)%nelem  ! parallelize this loop
!     e = colors(icol)%elem(k)
!     ... compute Fint for element #e and assemble ...
!   enddo
! enddo

end subroutine compute_Fint

subroutine compute_Fint_EP(f, pb)
  use fields_class, only : FIELD_add_elem
  use mat_gen, only : MAT_Fint_EP

  double precision, dimension(:,:), intent(out) :: f
  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: floc
  integer :: e

  f = 0d0
  do e = 1,pb%grid%nelem
    call MAT_Fint_EP(floc, pb%matwrk(e), pb%grid%ngll,pb%fields%ndof,pb%time%isdynamic)
    call FIELD_add_elem(floc,f,pb%grid%ibool(:,:,e)) ! assembly
  enddo

end subroutine compute_Fint_EP

double precision function norml2(x)
  double precision, dimension(:,:), intent(in) :: x
  norml2 = sqrt( sum(x*x) )
end function norml2

end module solver
