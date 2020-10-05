module solver

! SOLVER: solver for elasto-dynamic equation
!         M.a = - K.d + F

  use problem_class
  use stdio, only : IO_Abort
  use sources, only : SO_add
  use bc_gen , only : BC_apply, BC_set, bc_apply_kind, bc_select_kind, bc_set_kind, &
                      bc_select_fix, bc_set_fix_zero, bc_trans, bc_update_dfault, &
                      bc_has_dynflt, bc_update_bcdv

  implicit none
  private

  public :: solve

contains

!=====================================================================

subroutine solve(pb)

  type(problem_type), intent(inout) :: pb
  select case (pb%time%kind)
    case ('adaptive')
      call solve_adaptive(pb)
    case ('leapfrog')
      call solve_leapfrog(pb)
    case ('newmark')
      call solve_Newmark(pb)
    case ('HHT-alpha')
      call solve_HHT_alpha(pb)
    case default
      call solve_symplectic(pb)
  end select

end subroutine solve

!=====================================================================
! adaptive solver, a wrapper for dynamic and quasi-static solver
! depending on the 'isDynamic' flag
!
subroutine solve_adaptive(pb)
  type(problem_type), intent(inout) :: pb

  if (.not. pb%time%isDynamic) then
      ! quasi-static
      call solve_quasi_static(pb)
  else
      ! dynamic
      call solve_dynamic(pb)
  end if

  if (.not. pb%time%fixdt) call update_adaptive_step(pb)

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
subroutine solve_dynamic(pb)
    type(problem_type), intent(inout) :: pb
    select case (pb%time%kind_dyn)
      case ('leapfrog')
        call solve_leapfrog(pb)
      case ('newmark')
        call solve_Newmark(pb)
      case ('HHT-alpha')
        call solve_HHT_alpha(pb)
      case default
        call solve_symplectic(pb)
    end select
end subroutine solve_dynamic

! SOLVE: advance ONE time step
!        using single predictor-corrector Newmark (explicit)
!        in acceleration form
subroutine solve_Newmark(pb)

  type(problem_type), intent(inout) :: pb
  
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
subroutine solve_HHT_alpha(pb)

  type(problem_type), intent(inout) :: pb
  
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

subroutine solve_leapfrog(pb)

  type(problem_type), intent(inout) :: pb
  
  double precision, dimension(:,:), pointer :: d,v_mid,a,f

  d => pb%fields%displ
  v_mid => pb%fields%veloc
  a => pb%fields%accel
  f => a
                       
  d = d + pb%time%dt * v_mid
   
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
subroutine solve_symplectic(pb)

  type(problem_type), intent(inout) :: pb
  
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
             IS_DYNFLT = 6
         
  integer :: IS_DIR    = 2, IS_NEU = 1 

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
  v_pre = pb%fields%veloc
  
  ! update displacement to obtain a better initial guess for pcg solver
  ! update both d and pb%fields%disp
  d   = d + dt * v 

  ! apply dirichlet boundary condition and kinematic boundary condition
  ! set displacement to desired value and zero-out forcing f at DIR nodes
  call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_DIRNEU, IS_DIR)

  ! transform fields, apply half slip rate on side -1 and then transform back
  ! different from dynamic, update displacement as well.
  call bc_apply_kind(pb%bc, pb%time, pb%fields, f, IS_KINFLT)

  d_fix = 0.0d0


  ! check if there's dynflt boundary
  has_dynflt = bc_has_dynflt (pb%bc)
  
  ! start 2 passes

  do i = 1, 2

    write(*, *) "quasi-static solve, pass ", i
    ! update fault displacement for rsf 
    ! d(fault) = d_pre(fault) + dt * v(fault)
    ! v  = 0.5 * (vpre + v)
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
    call pcg_solver(d, f, pb)
    d = d + d_fix

    ! update velocity
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
  enddo
  
  ! declare final slip values on the fault
  call bc_update_bcdv(pb%bc, d, pb%time%time)
  a  = (v - v_pre)/dt ! a crude estimate of acceleration
  
end subroutine solve_quasi_static

!subroutine solve_quasi_static(pb)
!  type(problem_type), intent(inout) :: pb
!
!  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: d_pre, d_medium, d_fault, d, f
!  ! Note the d_fault and d_medium still have the same dimension of the d and d_pre.
!  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: v_pre, v_plate, v_f0
!  double precision, parameter :: tolerance = 10.0d-5
!  double precision :: plate_rate
!  integer :: i
! 
!  ! store initial 
!  ! the plate rate should be treated as dirichlet boundary condition
!  plate_rate = (2.0d-3)/(365*24*60*60) !DEVEL Trevor: ad hoc 
!  v_pre = pb%fields%veloc
!  d_pre = pb%fields%displ
! 
!  ! create field with plate velocity on fault, zeros in medium
!  v_f0 = 0d0
!  call BC_set(pb%bc, v_f0, plate_rate, v_plate)
!  ! Note that BC_set is only used for DYNFLT with RSF.
!  
!  ! add plate velocity to fault
!  ! veloc = v + vp on the fault plane 
!  ! veloc = v     off the fault plane
!
!  pb%fields%veloc = pb%fields%veloc + v_plate 
!
!  ! update the slip velocity
!  v_f0 = pb%fields%veloc
!
!  do i=1,2
!    ! correct velocity for improved estimate and
!    ! make prediction of all displacements 
!    pb%fields%veloc = 0.5*(v_f0 + pb%fields%veloc)
!    d = d_pre + pb%time%dt*pb%fields%veloc
!   
!    ! set displacements in medium to be zero
!   call BC_set(pb%bc, d, 0.0d0, d_medium)
!    d_fault = d - d_medium
!                  
!    ! solve for forces created by d_fault, finding K_{21}d^f 
!    ! (v_pre doesn't influence force, just energy, which isn't used)
!    call compute_Fint(f, d_fault, pb%fields%veloc, pb)
!    f = -f
!    call BC_set(pb%bc, f, 0.0d0, f) 
!
!    ! inputting the previous displacements as the initial guess into
!    ! (preconditioned) conjugate gradient method solver for the 
!    ! displacements in the medium
!    call BC_set(pb%bc, pb%fields%displ, 0.0d0, d_medium)
!    call pcg_solver(d_medium, f, pb, tolerance)
!    
!    ! combine displacements and calculate forces
!    d = d_fault + d_medium
!    call compute_Fint(f, d, pb%fields%veloc, pb)
!
!    ! apply boundary conditions 
!    call BC_apply(pb%bc, pb%time, pb%fields, f)
!  enddo
!
!  ! subtract plate velocity from fault for delta v 
!  pb%fields%veloc = pb%fields%veloc - v_plate
!
!  ! store final displacements
!  pb%fields%displ = d
!
!end subroutine solve_quasi_static

!=====================================================================
! f = - K * d 
subroutine compute_Fint(f,d,v,pb)

  use fields_class, only : FIELD_get_elem, FIELD_add_elem
  use mat_gen, only : MAT_Fint

  double precision, dimension(:,:), intent(out) :: f
  double precision, dimension(:,:), intent(in) :: d,v
  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: dloc,vloc,floc
  double precision :: E_ep, E_el, sg(3), sgp(3)
  integer :: e

  f = 0d0
  pb%energy%E_el = 0d0
  pb%energy%sg   = 0d0
  pb%energy%sgp  = 0d0

  do e = 1,pb%grid%nelem
    dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))
    vloc = FIELD_get_elem(v,pb%grid%ibool(:,:,e))
    call MAT_Fint(floc,dloc,vloc,pb%matpro(e),pb%matwrk(e), & 
                   pb%grid%ngll,pb%fields%ndof,pb%time%dt,pb%grid, &
                   E_ep,E_el,sg,sgp)
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

subroutine pcg_solver(d, f, pb)

  double precision, dimension(:,:), intent(inout) :: d
  double precision, dimension(:,:), intent(inout) :: f
  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: r, p, K_p, z
  double precision :: alpha_n, alpha_d, alpha, beta_n, beta_d, beta
  double precision :: norm_f, norm_r, tolerance
  integer :: maxIterations
  ! one hardwired parameter to ensure stable division
  double precision, parameter :: eps_stable = 1.0d-15
  integer :: it  

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
  pb%invKDiag = 1d0/pb%invKDiag
  call bc_trans(pb%bc, pb%invKDiag, -1)
  
  ! might be zero at fault node 2, stablize before division
  where (abs(pb%invKDiag)<eps_stable) pb%invKDiag = eps_stable
  pb%invKDiag = 1d0/pb%invKDiag

  z = pb%invKDiag * r ! z'
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
        print *, "PCG solver converges in ", it, " iterations." 
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
