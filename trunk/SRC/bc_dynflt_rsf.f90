module bc_dynflt_rsf

! BC_DYNFLT_RSF: rate and state dependent friction for dynamic faults
! DEVEL: module in progress ... add theta, write friction solver

  use distribution_cd
  use stdio, only: IO_abort

  implicit none
  private

  type rsf_input_type
    type(cd_type) :: dc, mus, a, b, Vstar, theta
  end type rsf_input_type

  type rsf_type
    private
    integer :: kind
    double precision, dimension(:), pointer :: dc=>null(), mus=>null(), a=>null(), b=>null(), &
                                               Vstar=>null(), theta=>null(), Tc=>null(), coeft=>null()
    double precision :: dt
    type(rsf_input_type) :: input
  end type rsf_type

  public :: rsf_type, rsf_read, rsf_init, rsf_mu, rsf_solver, rsf_qs_solver

contains

!---------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_RSF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Velocity and state dependent friction
! SYNTAX : &BC_DYNFLT_RSF kind, Dc | DcH, Mus | MusH , 
!                         a | aH, b | bH, Vstar | VstarH /
!          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
!          arguments with suffix H, if present, in the order listed above.
!
! ARG: kind     [int] [1] Type of rate-and-state friction law:
!                       1 = strong velocity-weakening at high speed
!                           as in Ampuero and Ben-Zion (2008)
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: a        [dble] [0.01d0] Direct effect coefficient
! ARG: b        [dble] [0.02d0] Evolution effect coefficient
! ARG: Vstar    [dble] [1d0] Characteristic or reference slip velocity
! ARG: theta    [dble] [1d0] State variable
!
! END INPUT BLOCK

! not implement yet:
!                       2 = logarithmic rate-and-state with aging state law
!                       3 = logarithmic rate-and-state with slip state law

! Read parameters from input file
  subroutine rsf_read(rsf,iin)

  use echo, only : echo_input,iout

  type(rsf_type), intent(out) :: rsf
  integer, intent(in) :: iin

  double precision :: Dc,MuS,a,b,Vstar,theta
  character(20) :: DcH,MuSH,aH,bH,VstarH,thetaH
  integer :: kind
  character(25) :: kind_txt

  NAMELIST / BC_DYNFLT_RSF / kind,Dc,MuS,a,b,Vstar,theta,DcH,MuSH,aH,bH,VstarH,thetaH

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  a = 0.01d0
  b = 0.02d0
  Vstar = 1d0
  theta = 0d0;
  DcH = ''
  MuSH = ''
  aH = ''
  bH = ''
  VstarH = ''
  thetaH = '';

  read(iin,BC_DYNFLT_RSF,END=300)
300 continue
  
  select case (kind)
    case(1); kind_txt = 'Strong velocity-weakening'
    case(2); kind_txt = 'Classical with aging law'
    case(3); kind_txt = 'Classical with slip law'
    case default; call IO_abort('BC_DYNFLT_RSF: invalid kind')
  end select
  rsf%kind = kind
  
  call DIST_CD_Read(rsf%input%Dc,Dc,DcH,iin,DcH)
  call DIST_CD_Read(rsf%input%MuS,MuS,MuSH,iin,MuSH)
  call DIST_CD_Read(rsf%input%a,a,aH,iin,aH)
  call DIST_CD_Read(rsf%input%b,b,bH,iin,bH)
  call DIST_CD_Read(rsf%input%Vstar,Vstar,VstarH,iin,VstarH)
  call DIST_CD_Read(rsf%input%theta,theta,thetaH,iin,thetaH)

  if (echo_input) write(iout,400) kind_txt,DcH,MuSH,aH,bH,VstarH,thetaH

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = Rate and State Dependent', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Direct effect coefficient . . . . . .(a) = ',A,&
            /5x,'  Evolution effect coefficient  . . . .(b) = ',A,&
            /5x,'  Velocity scale  . . . . . . . . .(Vstar) = ',A,&
            /5x,'  State variable  . . . . . . . . .(theta) = ',A)

  end subroutine rsf_read

!=====================================================================
! Initialize parameters
  subroutine rsf_init(rsf,coord,dt)

  type(rsf_type), intent(inout) :: rsf
  double precision, intent(in) :: coord(:,:),dt

  integer :: n
  
  call DIST_CD_Init(rsf%input%dc,coord,rsf%dc)
  call DIST_CD_Init(rsf%input%mus,coord,rsf%mus)
  call DIST_CD_Init(rsf%input%a,coord,rsf%a)
  call DIST_CD_Init(rsf%input%b,coord,rsf%b)
  call DIST_CD_Init(rsf%input%Vstar,coord,rsf%Vstar)
  call DIST_CD_Init(rsf%input%theta,coord,rsf%theta)
  
  n = size(coord,2)
  allocate(rsf%Tc(n))
  allocate(rsf%coeft(n))
                                              
  rsf%Tc = rsf%dc / rsf%Vstar
  rsf%coeft = exp(-dt/rsf%Tc)
  rsf%dt = dt
  end subroutine rsf_init

!=====================================================================
! Friction coefficient
  function rsf_mu(v,f) result(mu)

  double precision, dimension(:), intent(in) :: v
  type(rsf_type), intent(in) :: f
  double precision :: mu(size(v))

  select case(f%kind)
    case(1) 
      mu = f%mus +f%a*abs(v)/(abs(v)+f%Vstar) - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3) 
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      ! mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      mu = f%a*asinh(abs(v)/(2d0*f%Vstar)*exp((f%mus+f%b*log(f%Vstar*f%theta/f%Dc))/f%a))
  end select

  end function rsf_mu

!---------------------------------------------------------------------
! friction coefficient without the direct effect 
! (i.e. without the term that depends explicitly on slip velocity V)
  function rsf_mu_no_direct(f) result(mu)

  type(rsf_type), intent(in) :: f
  double precision :: mu(size(f%mus))

  select case(f%kind)
    case(1)
      mu = f%mus - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3) 
      !  Kaneko et al (2008) Eq. 12
      mu = f%mus + f%b*log(f%theta*f%Vstar/f%Dc) 
  end select

  end function rsf_mu_no_direct

!---------------------------------------------------------------------
!  Slip velocity at theta, tau
  function rsf_v(f,tau,sigma) result(v)

  type(rsf_type), intent(in) :: f
  double precision, dimension(:), intent(in) :: tau, sigma
  double precision, dimension(size(f%theta)) :: tmp, v

  !  Kaneko et al (2011) Eq. 14
  tmp = f%mus +f%b*log(f%Vstar*f%theta/f%Dc)
  tmp = 2d0*f%Vstar*exp(-tmp/f%a)
  v = sinh(tau/(-sigma*f%a))*tmp

  end function rsf_v

!=====================================================================
!       --------------- RSF_SOLVER subroutine -----------------
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!      3. update theta again, now using V=(Vnew+Vold)/2
!      4. solve again for Vnew
!
! Note: most often sigma is negative (compressive) 
!
  subroutine rsf_solver(v,tau_stick,sigma,f,Z)

  double precision, dimension(:), intent(inout) :: v
  double precision, dimension(:), intent(in) :: sigma,Z,tau_stick
  type(rsf_type), intent(inout) :: f

  double precision, dimension(size(v)) :: v_new,theta_new
  ! First pass: 
  theta_new = rsf_update_theta(f%theta,v,f)
  v_new = rsf_update_V(tau_stick, sigma, f, theta_new, Z)
  ! Second pass:
  theta_new = rsf_update_theta(f%theta,0.5d0*(v+v_new), f)
  v_new = rsf_update_V(tau_stick, sigma, f, theta_new, Z)
  
  ! store new velocity and state variable estimate in friction law 
  f%theta = theta_new 
  v = v_new

  end subroutine rsf_solver

!=====================================================================
!       --------------- RSF_QS_SOLVER subroutine -----------------
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!
! Note: most often sigma is negative (compressive) 
!
  subroutine rsf_qs_solver(v,tau,sigma,f)

  double precision, dimension(:), intent(inout) :: v
  double precision, dimension(:), intent(in) :: sigma,tau
  type(rsf_type), intent(inout) :: f
  
  f%theta = rsf_update_theta(f%theta,v,f)
  v = rsf_v(f, tau, sigma)

  end subroutine rsf_qs_solver

!---------------------------------------------------------------------
! Update state variable (theta) assuming slip velocity (v) is known

  function rsf_update_theta(theta,v,f) result(theta_new)

  double precision, dimension(:), intent(in) :: v,theta
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(v)) :: theta_new

  select case(f%kind)
    case(1) 
     ! exact integration assuming constant V over the timestep
     ! Tc = Dc/Vstar
     ! coeft = exp(-dt/Tc)
      theta_new = theta*f%coeft +f%Tc*abs(v)*(1d0-f%coeft)

    case(2) 
     ! Kaneko et al (2008) eq 19 - "Aging Law"
     ! theta_new = (theta-Dc/v)*exp(-v*dt/Dc) + Dc/v
      theta_new = f%Dc/abs(v)
      theta_new = (theta-theta_new)*exp(-f%dt/theta_new) + theta_new
    case(3) 
     ! Kaneko et al (2008) eq 20 - "Slip Law"
     ! theta_new = Dc/v *(theta*v/Dc)**exp(-v*dt/Dc)
      theta_new = f%Dc/abs(v)
      theta_new = theta_new *(theta/theta_new)**exp(-f%dt/theta_new)
  end select

  end function rsf_update_theta

!---------------------------------------------------------------------
! Update slip velocity assuming theta is known
! The constraints are
!   (abs(tau)-strength)*v = 0
!   abs(tau)-strength <= 0
!   sign(tau) = sign(v)
! where
!   strength = -sigma*( mu(theta) +a*v/(1+v) )
!   tau = tau_stick-Z*v
!
! Inherited from the SBIEM code BIMAT-PCSI
! WARNING: the SBIEM code assumed v>0
!          We should allow here for any sign of v
!          Exploit the fact that sign(tau)=sign(tau_stick) (because mu>0)

  function rsf_update_V(tau_stick,sigma,f,theta,Z) result(v)
   
  double precision, dimension(:), intent(in) :: tau_stick,sigma,theta,Z
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(tau_stick)) :: v
  double precision :: tmp(size(tau_stick)), tolerance, estimateLow, estimateHigh
  integer :: it

!  strength = -sigma*rsf_mu_no_direct(v,f) 
!  v = (tau_stick-strength)/Z

  select case(f%kind)
    case(1) 
      v = (tau_stick +sigma*rsf_mu_no_direct(f))/Z ! if v<0 will stop
      tmp = v -f%Vstar +sigma*f%a/Z
      v = 0.5d0*( tmp +sqrt(tmp*tmp +4d0*v*f%Vstar) )
      v = max(0d0,v)  ! arrest if v<0 
 
    case(2,3) 
     ! "Aging Law and Slip Law"
     ! Find each element's velocity:
     do it=1,size(tau_stick)
       !DEVEL: What are the accepted tolerances and bounds? User-input? 
       tolerance=0.001*f%a(it)*sigma(it) ! As used by Kaneko in MATLAB code
       estimateLow = -10.0 
       estimateHigh = 10.0
       v(it)=nr_solver(nr_fric_func_tau,estimateLow,estimateHigh,tolerance,f,it,theta(it),tau_stick(it),sigma(it),Z(it))
     enddo     
        
  end select

  end function rsf_update_V

!==================================================================
!         Newton-Raphson algorithm with bisection step         
! from  Numerical Recipes in Fortran 90 by Press et al.  (1996)

! N-R method, must provide a function (nr_fric_func) that gives
! the value of the function and the value of the derivative of the 
! function at a point and it finds a root to 0=nr_fric_func bounded 
! by [x1, x2] w/ error < x_acc

  function nr_solver(nr_fric_func_tau, xL, xR, x_acc, f, it, theta, tau_stick, sigma, Z) result(v)
  !WARNING: This is not a well tested portion !!!
 
  integer, parameter :: maxIteration=200
  integer :: is
  
  double precision, intent(in) :: xL, xR, x_acc
  double precision :: xLeft, xRight, x_est, v
  double precision :: dfunc_dx, dx, dx_old, func_x, f_high, f_low, x_high, x_low
  double precision :: temp
  
  ! Friction parameters:
  external nr_fric_func_tau
  type(rsf_type), intent(in) :: f
  double precision, intent(in) :: theta, tau_stick, sigma, Z
  integer, intent(in) :: it

  ! Find initial function estimates (dfunc_dx not needed yet)
  call nr_fric_func_tau(xL, f_low, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  call nr_fric_func_tau(xR, f_high, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  xLeft = xL
  xRight = xR  

  ! Ensure zero is bounded:
  do while (f_low*f_high>0)
    xLeft = xLeft*100
    xRight = xRight*100
    call nr_fric_func_tau(xLeft, f_low, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
    call nr_fric_func_tau(xRight, f_high, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  enddo
 
  ! Lucky guesses: 
  if (f_low==0)  then
    x_est=xLeft
    return
  else if (f_high==0) then
    x_est=xRight
    return
  ! Orient the search so that func_x(x_low) < 0
  else if (f_low<0) then 
    x_low=xLeft
    x_high=xRight
  else 
    x_high=xLeft
    x_low=xRight
  endif

  ! Initialize the guess and step sizes:
  x_est=.5d0*(xLeft+xRight)
  ! The stepsize before last
  dx_old=abs(xRight-xLeft)
  ! The last step:
  dx=dx_old
  
  call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)

  ! Loop over allowed iterations:
  do is=1,maxIteration
   
    ! Bisect if N-R out of range, or not decreasing fast enough:
    if( ((x_est-x_high)*dfunc_dx-func_x)*((x_est-x_low)*dfunc_dx-func_x)>0  .or. abs(2*func_x)>abs(dx_old*dfunc_dx) ) then
      dx_old=dx
      dx=0.5*(x_high-x_low)
      x_est=x_low+dx
      !  Check if change is negligible:
      if (x_low==x_est) return

    ! The Newton step is acceptable, move forward with algorithm:
    else
      dx_old=dx
      dx=func_x/dfunc_dx
      temp=x_est
      x_est=x_est-dx
      !  Check if change is negligible:
      if (temp==x_est) return
    endif
  
    ! Check convergence criterion:
    if (abs(dx)<abs(x_acc)) return
  
    ! Evaluate function with new estimate of x
    call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)

    ! Redefine the bounds for the next loop
    if (func_x<0) then
      x_low=x_est
    else 
      x_high=x_est
    endif
  enddo
  print*,'tau = ',x_est
  print*,'delta = ',dx,' > ',x_acc
  print*,'vel = ',v
  call IO_abort('NR_Solver has exceeded the maximum iterations (200)')
  
  end function nr_solver

!==================================================================
!         Newton-Raphson algorithm from Kaneko's MATLAB code
! uses NR algorithm from Kaneko's code used in 2011 paper, steps by
! delta to specified tolerance

  function nr_solver_Kaneko(nr_fric_func_tau, x_est, x_acc, f, theta, it, tau_stick, sigma, Z) result(v)
  !WARNING: This is not a well tested portion !!!
 
  integer, parameter :: maxIteration=10
  integer :: is
  
  double precision, intent(in) :: x_acc
  double precision, intent(inout) :: x_est
  double precision :: v
  double precision :: dfunc_dx, dx, func_x 
  double precision :: temp
  
  ! Friction parameters:
  external nr_fric_func_tau
  type(rsf_type), intent(in) :: f
  double precision, intent(in) :: tau_stick, sigma, Z, theta
  integer, intent(in) :: it

  ! Find initial function estimates 
  call nr_fric_func_tau(x_est, func_x, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  
  ! Loop over allowed iterations:
  do is=1,maxIteration
    dx=func_x/dfunc_dx
    temp=x_est
    x_est=x_est-dx
    !  Check if change is negligible:
    if (temp==x_est) return
    ! Check convergence criterion:
    if (abs(dx)<abs(x_acc)) then
      call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
      return
    endif
 
    if (abs(dx)>10**9) then
      print*,is
      print*,'v = ',v
      print*,'tau = ',x_est
      print*,'dx = ',dx
      call IO_abort('NR_Solver fails to converge')
    endif
    ! Evaluate function with new estimate of x
    call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  enddo
  
  print*,'v = ',v
  print*,'tau = ',x_est
  print*,'dx = ',dx
  call IO_abort('NR_Solver has exceeded the maximum iterations')
  
  end function nr_solver_Kaneko
!---------------------------------------------------------------------
!-----------Friction function for Newton Raphson method---------------
! This function returns the value of the friction function and its 
! derivative evaluated at a particular shear stress.  

subroutine nr_fric_func_tau(tau, func_tau, dfunc_dtau, v, f, theta, it, tau_stick, sigma, Z) 
  !WARNING: This is not a well tested portion!
  double precision, intent(out) :: func_tau, dfunc_dtau, v
  double precision, intent(in) ::  tau_stick,sigma,Z,tau,theta
  type(rsf_type), intent(in) :: f
  double precision :: dv_dtau, tmp
  integer, intent(in) :: it

      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      !  mu = mus +a*log(v/Vstar) + b*log(theta*Vstar/Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Lapusta et al. (2000))
      !  solved in terms of v
      tmp = f%mus(it) +f%b(it)*log( f%Vstar(it)*theta/f%Dc(it) )
      tmp = 2d0*f%Vstar(it)*exp(-tmp/f%a(it))
      v = sinh(tau/(-sigma*f%a(it)))*tmp
      func_tau = tau_stick - Z*v - tau
    
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used)
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      !  solved in terms of v, derivative wrt tau
      dv_dtau = cosh(tau/(-sigma*f%a(it)))*tmp/(-sigma*f%a(it))
      dfunc_dtau = -Z*dv_dtau - 1d0
  
end subroutine nr_fric_func_tau

subroutine nr_fric_func_v(v, func_v, dfunc_dv, f, it, tau_stick, sigma, Z) 
  !WARNING: This is not a well tested portion!
  double precision, intent(out) :: func_v, dfunc_dv
  double precision, intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(in) :: f
  double precision :: v, mu, dmu_dv
  integer, intent(in) :: it

      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      ! mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      mu = (f%mus(it)+f%b(it)*log(f%Vstar(it)*f%theta(it)/f%Dc(it)))/f%a(it)
      mu = f%a(it)*asinh(v/(2*f%Vstar(it))*exp(mu))
      func_v = tau_stick - Z*v - sigma*mu
      
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      !  dmu_dv = f%a*/v 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      dmu_dv = exp((f%mus(it)+f%b(it)*log(f%Vstar(it)*f%theta(it)/f%Dc(it)))/f%a(it))/(2*f%Vstar(it))
      dmu_dv = dmu_dv*f%a(it)/(sqrt(1+(dmu_dv*v)**2))
      dfunc_dv = -Z - sigma*dmu_dv
  
end subroutine nr_fric_func_v

end module bc_dynflt_rsf
