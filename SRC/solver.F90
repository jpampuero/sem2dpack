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
  use bc_gen , only : BC_apply, bc_apply_kind, bc_select_kind, bc_set_kind, &
                      bc_select_fix, bc_set_fix_zero, bc_trans, bc_update_dfault, &
                      bc_update_bcdv, bc_reset, bc_has_dynflt, BC_D2S_ReInit
  use fields_class, only: FIELD_SetVecFromField, FIELD_SetFieldFromVec 

  implicit none
  private

  public :: solve

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
  integer :: rank
  
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)
  
  id1 = pb%time%isdynamic
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
      ! work with total displacement (relative to di) 
      if (pb%time%isUsePetsc) then
          call solve_quasi_static_petsc(pb, petobj)
      else
          call solve_quasi_static(pb)
      end if
      ! update background values for dynamic step
      pb%fields%dn = pb%fields%displ
  end if

  if (rank==0) then
      if (.not. pb%time%fixdt) call update_adaptive_step(pb)
  end if

  ! broad cast constants in t
  call TIME_broadcast(pb%time)

  id2 = pb%time%isdynamic
  isw = id2 .and. (.not. id1)
  if (isw) then
     ! keep only velocity on the boundaries as zero
  !   tmp=0d0
!     write(*, *)"HERE, switching!"
   !  call bc_select_kind(pb%bc, pb%fields%veloc,tmp, IS_DYNFLT)
   ! switch from static to dynamic
     pb%fields%veloc = 0d0
 !    write(*, *)"Done setting velocity!"
  end if

! switch from dynamic to static
  isw = id1 .and. (.not. id2)
  if (isw) then
      ! switch from dynamic to static
      pb%fields%displ = 0d0 
      pb%fields%accel = 0d0
      if (rank==0) then
          call BC_D2S_ReInit(pb%bc, pb%fields%veloc)
      end if
      ! reinitialize the fields
  end if
  
end subroutine solve_adaptive

! update time step 
subroutine update_adaptive_step(pb)
  use bc_gen, only: BC_timestep
  type(problem_type), intent(inout) :: pb

  ! rate-state time step following Lapusta et al., 2010
  call BC_timestep(pb%bc,pb%time)
end subroutine update_adaptive_step

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

subroutine solve_quasi_static(pb)
  type(problem_type), intent(inout) :: pb
  double precision, dimension(:,:), pointer :: d, v, a, f
  double precision, dimension(pb%fields%npoin, pb%fields%ndof) :: d_pre, d_fix, v_pre
  
  ! vplate is built into the rate-state b.c.
  double precision :: dt 
  integer :: i
  logical :: has_dynflt = .false.

  integer :: IS_DIRNEU = 1, & 
             IS_KINFLT = 2, & 
             IS_DYNFLT = 6, &
             IS_DIRABS = 7
         
  integer :: IS_DIR    = 2, IS_NEU = 1, flag = 0 

  ! points to displacement
  ! save some memory by modifying in place
  d => pb%fields%displ
  v => pb%fields%veloc
  ! use acceleration and f to work with force
  a => pb%fields%accel
  f => pb%fields%accel

  dt = pb%time%dt
  
  ! use additional memory to keep the previous displacement
  ! use to update predictor for rsf faults
  d_pre = pb%fields%displ
  
  ! update displacement to obtain a better initial guess for pcg solver
  ! update both d and pb%fields%disp
  d   = d + dt * v 

  ! apply dirichlet boundary condition and kinematic boundary condition
  ! set displacement to desired value and zero-out forcing f at DIR nodes
  call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_DIR)
  call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRABS, IS_DIR)
  v_pre = pb%fields%veloc

  ! transform fields, apply half slip rate on side -1 and then transform back
  ! different from dynamic, update displacement as well.
  call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_KINFLT)

  d_fix = 0.0d0


  ! check if there's dynflt boundary
  has_dynflt = bc_has_dynflt (pb%bc)
  
  ! start 2 passes
  pb%time%pcg_iters = 0

  do i = 1, 2

    ! update fault displacement for rsf 
    ! d(fault) = d_pre(fault) + dt * v(fault)
    ! v  = 0.5 * (vpre + v)
    !d = d_pre + (v + v_pre)/2d0 *dt
    call bc_update_dfault(pb%bc, d_pre, d, (v + v_pre)/2d0, dt)
    
    ! obtain d_fix, d'_fix rotated back to d_fix

    call bc_select_fix(pb%bc, d, d_fix)

    ! solve for forces created by d_fix, finding K_{21}d^f + K_{23}d^{dir} 
    call compute_Fint(f, d_fix, pb%fields%veloc, pb)
   
    ! apply newmann if there's any, used to solve dofs in the medium
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_NEU) 
    
    f = -f
   
    ! set the fixed dofs to 0
    call BC_trans(pb%bc, d,  1)
    call BC_set_fix_zero(pb%bc, d)
    call BC_trans(pb%bc, d, -1)
    
    ! (preconditioned) conjugate gradient method solver
    call pcg_solver(d, f, pb, flag)
    d = d + d_fix

    pb%time%pcg_iters(i) = flag 

    ! update velocity for the off fault dofs
    v = (d - d_pre)/dt

    ! detect if there's dynflt fault boundary
    ! if not, break the loop and exit after the pcg solver

    if (.not. has_dynflt) exit
    ! apply boundary conditions including following: 
    !    1.  compute fault traction, compute state variable 
    !    2.  update slip rate, update fault velocity in fields%veloc 
    !
    ! velocity is modified inside this subroutine

    ! recompute the force use the updated total disp
    call compute_Fint(f, d, pb%fields%veloc, pb)
    call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DYNFLT)

    ! reapply dichlet boundary condition to overwrite the fault tip node
    ! if there's a conflict between tip and Dirichlet bc.
    !
    ! Note this issue only occurs when simulating half domain 
    ! call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_DIR)
    ! call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRABS, IS_DIR)
  enddo
  
  ! declare final slip values on the fault
  call bc_update_bcdv(pb%bc, d, pb%time%time)
  a   = 0d0
  !a  = (v - v_pre)/dt ! a crude estimate of acceleration
  
end subroutine solve_quasi_static

! solve the mechanical balance 2 passes with dt
! compute dtmax to update the damage variable
! if dtmax<dt
!    dt = dtmax
!    resolve the mechanical balance, 2 passes
!
subroutine solve_quasi_static_petsc(pb, petobj)
#include <petsc/finclude/petscksp.h>
  use petscksp
  type(problem_type), intent(inout) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  double precision, dimension(:,:), pointer :: d, v, a, f
  double precision, dimension(pb%fields%npoin, pb%fields%ndof) :: d_pre, v_pre, fp
  
  ! vplate is built into the rate-state b.c.
  double precision :: dt, dtmax, tn 
  integer :: i, ndof, npoin,j
  logical :: has_dynflt = .false.

  integer :: IS_DIRNEU = 1, & 
             IS_KINFLT = 2, & 
             IS_DYNFLT = 6, &
             IS_DIRABS = 7
         
  integer :: IS_DIR    = 2, IS_NEU = 1
  PetscErrorCode :: ierr
  PetscScalar, pointer :: xx_d(:)
  PetscScalar, pointer :: xx_b(:)
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
  
  do j = 1, 2
    pb%time%pcg_iters = 0
    if (rank==0) then 
    ! update displacement to obtain a better initial guess for pcg solver
    ! update both d and pb%fields%disp
    
    d  = d_pre + dt * v_pre 

    ! apply dirichlet boundary condition and kinematic boundary condition
    ! set displacement to desired value and zero-out forcing f at DIR nodes
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

!    if (print_time) then 
!    write(*, *) 'Elapse Time for VecScatter fix dof =', elapse_time, 's'
!    end if

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
  !a  = (v - v_pre)/dt ! a crude estimate of acceleration
  
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

    call MAT_Fint(floc,dloc,vloc,pb%matwrk(e), & 
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

!=====================================================================
! NOT SUPPORTED
subroutine cg_solver(d, f, pb, tolerance)
  double precision, dimension(:,:), intent(inout) :: d
  double precision, dimension(:,:), intent(in) :: f
  double precision, intent(in) :: tolerance
  type(problem_type), intent(inout) :: pb

  ! internal variables
  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: r, p, K_p
  double precision, dimension(pb%fields%ndof, pb%fields%ndof) :: alpha_n, alpha_d, alpha, beta_n, beta_d, beta
  double precision :: norm_f,norm_r
  integer, parameter :: maxIterations=10000
  integer :: IS_DIRNEU = 1, & 
             IS_DYNFLT = 6
         
  integer :: IS_DIR    = 2, & 
             IS_NEU    = 1 
  integer :: it  
 
  call compute_Fint(K_p,p,pb%fields%veloc,pb)
  call BC_set_kind(pb%bc, K_p, 0.0d0, IS_DYNFLT, IS_NEU)
  call BC_set_kind(pb%bc, K_p, 0.0d0, IS_DIRNEU, IS_DIR)
  !call BC_set(pb%bc, K_p, 0.0d0, K_p) 
  
  ! initial residual
  r = f - K_p
  norm_f = norml2(f)

  p = r

  do it=1,maxIterations
    ! compute stiffness*p for later steps
    call compute_Fint(K_p,p,pb%fields%veloc,pb)
    call BC_set_kind(pb%bc, K_p, 0.0d0, IS_DYNFLT, IS_NEU)
    call BC_set_kind(pb%bc, K_p, 0.0d0, IS_DIRNEU, IS_DIR)
    !call BC_set(pb%bc, K_p, 0.0d0, K_p) 
    
    ! find step length
    alpha_n = sum(r*r)
    alpha_d = sum(p*K_p)

    alpha = alpha_n/alpha_d


    ! determine approximate solution
    d = d + alpha*p

    ! find the residual
    r = r - alpha*K_p

    ! test if within tolerance
    norm_r = norml2(r)
    if (norm_r/norm_f < tolerance) return
    if (norm_r/norm_f > 10.0d10) call IO_Abort('Conjugate Gradient does not converge')

    ! improve the step
    beta_n = sum(r*r)
    beta_d = alpha_n
    beta = alpha_n/beta_d

    ! search direction
    p = r + beta*p
  enddo
  
  call IO_Abort('Conjugate Gradient does not converge')

end subroutine cg_solver

double precision function norml2(x)
  double precision, dimension(:,:), intent(in) :: x
  norml2 = sqrt( sum(x*x) )
end function norml2


!=====================================================================
! Preconditioned conjugate gradient with Jacobi preconditioner
!
! Jacobi preconditioner is precomputed, inverse of diagonal of stiffness
! matrix. Global stiffness matrix is never assembled.
!
! d, f are untransformed as inputs
!

subroutine pcg_solver(d, f, pb, flag)

  double precision, dimension(:,:), intent(inout) :: d
  double precision, dimension(:,:), intent(inout) :: f
  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: r, p, K_p, z
  double precision :: alpha_n, alpha_d, alpha, beta_n, beta_d, beta
  double precision :: norm_f, norm_r, tolerance
  integer :: maxIterations
  ! one hardwired parameter to ensure stable division
  double precision, parameter :: eps_stable = 1.0d-15
  integer :: it, flag

  flag = -1

  tolerance     = pb%time%TolLin
  maxIterations = pb%time%MaxIterLin
  
  call bc_trans(pb%bc, f, -1) ! f -> f'
  call bc_trans(pb%bc, d,  1) ! d to d'

  ! From here on work with f', d', p', K' * p'
  call compute_kp_prime(K_p, d, pb)
  
  ! set the left hand side of fixed DOFs to zero
  call BC_set_fix_zero(pb%bc, K_p)
  call BC_set_fix_zero(pb%bc,   f)

  r = f - K_p

  norm_f = norml2(f) ! initial right hand side 

  ! obtain transformed preconditioner
!  invKDiag = pb%invKDiag
!  invKDiag = 1d0/invKDiag
!  call bc_trans(pb%bc, invKDiag, -1)
!  
!  ! might be zero at fault node 2, stablize before division
!  where (abs(invKDiag)<eps_stable) invKDiag = eps_stable
!  invKDiag = 1d0/invKDiag

  z = pb%invKDiagTrans * r ! z'
  p = z ! p'
  
  do it=1,maxIterations

    ! compute stiffness*p for later steps
    call compute_kp_prime(K_p, p, pb) 
    call BC_set_fix_zero(pb%bc, K_p)
    
    ! the following steps work with p', K' * p'
    ! find step length
    alpha_n = sum(z*r)
    alpha_d = sum(p*K_p)
    alpha = alpha_n/alpha_d

    ! determine approximate solution work with d'
    d = d + alpha*p

    ! find new residual
    r = r - alpha*K_p

    ! test if within tolerance
    norm_r = norml2(r)

    ! if |LHS-RHS|<tolerance*|RHS|

    if (norm_r/norm_f < tolerance) then 
        flag = it
        ! transform after convergence.
        ! convert d' to d, back transform
         call bc_trans(pb%bc, d, -1) 
        ! convert f' to f, forward transform
         call bc_trans(pb%bc, f,  1) 
        return
    end if
    if (norm_r/norm_f > 10.0d10) exit

    z = pb%invKDiag * r

    ! improve the step
    beta_n = sum(z*r)
    beta_d = alpha_n
    beta = beta_n/beta_d

    ! search direction
    p = z + beta*p
  enddo
  
  call IO_Abort('Preconditioned Conjugate Gradient does not converge')
  
end subroutine pcg_solver

! ==============================================================
! Compute K' * p' according to the transform of Junpei Seki, 2017 
!
subroutine compute_kp_prime(kp, pprime, pb)
  type(problem_type), intent(inout) :: pb
  double precision, dimension(:,:), intent(inout) :: pprime
  double precision, dimension(:,:), intent(inout) :: kp

  ! Transform P' to P
  ! P = X^{-1} * P'
  call bc_trans(pb%bc, pprime, -1) ! p' -> p 
  ! KP = K * P
  call compute_Fint(kp, pprime, pb%fields%veloc, pb) !K*p
  ! K'* P' = X^{-1} * KP = X^{-1} * K * X^{-1} P' = K' * P'
  call bc_trans(pb%bc, kp, -1)     ! K' * p'
  ! rotate P to P'
  call bc_trans(pb%bc, pprime, 1)  ! p -> p'

end subroutine compute_kp_prime

end module solver
