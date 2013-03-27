module solver

! SOLVER: Newmark solver for elasto-dynamic equation
!         M.a = - K.d + F

  use problem_class
  use stdio, only : IO_Abort
  use sources, only : SO_add
  use bc_gen , only : BC_set, BC_zero

  implicit none
  private

  public :: solve

contains

!=====================================================================

subroutine solve(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%time%kind)
    case ('leapfrog')
      call solve_leapfrog(pb)
    case ('newmark')
      call solve_Newmark(pb)
    case ('HHT-alpha')
      call solve_HHT_alpha(pb)
    case ('quasi-static')
      call solve_quasi_static(pb)
    case default
      call solve_symplectic(pb)
  end select

  end subroutine solve

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
  call BC_set(pb%bc,pb%time%time,pb%fields,f)

 ! NOTE: if a source is an incident wave, it is not added during
 !       "call SO_add" but during "call BC_set"
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
  double precision :: t_alpha,dt,alpha,beta,gamma

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
  call BC_set(pb%bc,t_alpha,pb%fields,f)
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
  call BC_set(pb%bc,pb%time%time,pb%fields,f)

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
subroutine solve_quasi_static(pb)
  type(problem_type), intent(inout) :: pb
  !type(problem_type), pointer :: pb_pointer
  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: d, d_medium, d_fault, v, f
  double precision :: tolerance
  
  ! make prediction of all displacements 
  print*,'v',size(pb%fields%veloc,1),'x',size(pb%fields%veloc,2)
  v = pb%fields%veloc
  print*,'v',size(v,1),'x',size(v,2)
  d = pb%fields%displ + pb%time%dt*v
  
  ! set displacements in medium to be zero
  print*,'bczero'
  call BC_zero(pb%bc, pb%fields, d_medium)
  d_fault = d - d_medium
                
  ! solve for forces, finding K_{21}d^f
  print*,'computefint'
  call compute_Fint(f, d_fault, pb%fields%veloc, pb)
  f = -f

  ! use previous displacements as initial guess for displacements in medium 
  ! and compute resulting forces 
  print*,'bczero'
  call BC_zero(pb%bc, pb%fields, d_medium)
  
  ! input the initial guess and the function that computes K*d into the
  ! (preconditioned) conjugate gradient method solver
  tolerance = 10.0d-5
  print*,'cg_solver()'
  call cg_solver(d_medium, f, pb, tolerance)
  ! call pcg_solver(d_medium, compute_Fint, f, pb, tolerance)
  
  ! combine displacements and calculate forces
  d = d_fault + d_medium
  call compute_Fint(f, d, pb%fields%veloc, pb)
  !        Output Input       Input   Input

  ! apply boundary conditions
  call BC_set(pb%bc, pb%time%time, pb%fields, f)
  ! this will return the velocity 

  ! calculate final displacements 
  d = pb%fields%displ + 0.5d0*pb%time%dt*(v + pb%fields%veloc)

  

end subroutine solve_quasi_static

!=====================================================================
!
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
!
subroutine cg_solver(d, f, pb, tolerance)
  double precision, dimension(:,:), intent(inout) :: d
  double precision, dimension(:,:), intent(in) :: f
  double precision, intent(in) :: tolerance
  type(problem_type), intent(inout) :: pb

  ! internal variables
  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: r_new, r_old, p, K_p
  double precision, dimension(pb%fields%ndof, pb%fields%ndof) :: alpha_n, alpha_d, alpha, beta_d, beta
  double precision :: norm_f
  integer, parameter :: maxIterations=20000
  integer :: it  
 
  norm_f = norm2(f)
  r_new = f
  r_old = f
  p = r_new

  do it=1,maxIterations
    print*,it
    ! compute stiffness*p for later steps
    call compute_Fint(K_p,p,pb%fields%veloc,pb)
    ! find step length
    alpha_n = matmul(transpose(r_new),r_new)
    alpha_d = matmul(transpose(p),K_p)
    alpha = alpha_n/alpha_d
    print*,'alpha',alpha
    ! determine approximate solution
    d = d + matmul(K_p,alpha)
    ! test if within tolerance
    print*,'|r_new|',norm2(r_new)
    print*,'|f|',norm_f
    print*,'|r_new|/|f|',norm2(r_new)/norm_f
    if (norm2(r_new)/norm_f < tolerance) return
    ! find the residual
    r_new = r_new - matmul(K_p,alpha)
    ! improve the step
    beta_d = matmul(transpose(r_old),r_old)
    beta = alpha_n/beta_d
    print*,'beta',beta
    ! search direction
    p = r_new + matmul(p,beta)
    ! store old residual
    r_old = r_new  
  enddo
  
  call IO_Abort('Conjugate Gradient does not converge')

end subroutine cg_solver

! DEVEL: not finished yet!  - invKDiag needs to be integrated.
subroutine pcg_solver(d, f, pb, tolerance)
  double precision, dimension(:,:), intent(inout) :: d
  double precision, dimension(:,:), intent(in) :: f
  double precision, intent(in) :: tolerance
  type(problem_type), intent(inout) :: pb

  ! internal variables
  double precision, dimension(pb%fields%npoin,pb%fields%ndof) :: r_new, r_old, p, K_p, z_new, z_old
  double precision, dimension(pb%fields%ndof, pb%fields%ndof) :: alpha_n, alpha_d, alpha, beta_d, beta
  double precision :: norm_f
  integer, parameter :: maxIterations=4000
  integer :: it  

  norm_f = norm2(f)
  call compute_Fint(K_p,d,pb%fields%veloc,pb)
  r_new = f - K_p
  r_old = r_new 
  z_new = matmul(pb%invKDiag,r_new)
  z_old = z_new
  p = z_new

  do it=1,maxIterations
    print*,it
    ! compute stiffness*p for later steps
    call compute_Fint(K_p,p,pb%fields%veloc,pb)
    ! find step length
    alpha_n = matmul(transpose(r_new),z_new)
    alpha_d = matmul(transpose(p),K_p)
    alpha = alpha_n/alpha_d
    print*,'alpha',alpha
    ! determine approximate solution
    d = d + matmul(K_p,alpha)
    ! test if within tolerance
    print*,'|r_new|',norm2(r_new)
    print*,'|f|',norm_f
    print*,'|r_new|/|f|',norm2(r_new)/norm_f
    if (norm2(r_new)/norm_f < tolerance) return
    ! find the residual
    r_new = r_new - matmul(K_p,alpha)
    z_new = matmul(pb%invKDiag,r_new)
    ! improve the step
    beta_d = matmul(transpose(z_old),r_old)
    beta = alpha_n/beta_d
    print*,'beta',beta
    ! search direction
    p = r_new + matmul(p,beta)
    ! store old residual
    r_old = r_new  
    z_old = z_new
  enddo
  
  call IO_Abort('Conjugate Gradient does not converge')
  

end subroutine pcg_solver

end module solver
