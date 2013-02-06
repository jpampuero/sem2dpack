module bc_dynflt_rsf

! BC_DYNFLT_RSF: rate and state dependent friction for dynamic faults
! DEVEL: module in progress ... add theta, write friction solver

  use distribution_cd
  use stdio, only: IO_abort

  implicit none
  private

  type rsf_input_type
    type(cd_type) :: dc, mus, a, b, Vstar
  end type rsf_input_type

  type rsf_type
    private
    integer :: kind
    double precision, dimension(:), pointer :: dc, mus, a, b, Vstar, theta, &
                                               Tc, coeft
    double precision :: dt
    type(rsf_input_type) :: input
  end type rsf_type

  public :: rsf_type, rsf_read, rsf_init, rsf_mu, rsf_solver

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

  double precision :: Dc,MuS,a,b,Vstar
  character(20) :: DcH,MuSH,aH,bH,VstarH
  integer :: kind
  character(20) :: kind_txt

  NAMELIST / BC_DYNFLT_RSF / kind,Dc,MuS,a,b,Vstar,DcH,MuSH,aH,bH,VstarH

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  a = 0.01d0
  b = 0.02d0
  Vstar = 1d0;
  DcH = ''
  MuSH = ''
  aH = ''
  bH = ''
  VstarH = '';

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

  if (echo_input) write(iout,400) kind_txt,DcH,MuSH,aH,bH,VstarH

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = rate and state dependent', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Direct effect coefficient . . . . . .(a) = ',A,&
            /5x,'  Evolution effect coefficient  . . . .(b) = ',A,&
            /5x,'  Velocity scale  . . . . . . . . .(Vstar) = ',A)

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

  n = size(coord,2)
  allocate(rsf%theta(n))
  allocate(rsf%Tc(n))
  allocate(rsf%coeft(n))
 !WARNING: theta initialization should be more general for Dieterich-Ruina rsf
!          Also needs option for input by user
  rsf%theta = 0d0 
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
      mu = f%mus +f%a*v/(v+f%Vstar) - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3) 
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      ! mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      mu = (f%a)*asinh((v/(2*f%Vstar))*exp((f%mus+b*log(f%Vstar*f%theta/f%Dc))/f%a))
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

!=====================================================================
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!      3. update theta again, now using V=(Vnew+Vold)/2
!      4. solve again for Vnew
  subroutine rsf_solver(v,tau_stick,sigma,f,Z)

  double precision, dimension(:), intent(inout) :: v
  double precision, dimension(:), intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(inout) :: f

  double precision, dimension(size(v)) :: v_new,theta_new

  theta_new = rsf_update_theta(f%theta,v,f)
  v_new = rsf_update_V(tau_stick, sigma, f, Z)
  theta_new = rsf_update_theta(f%theta,0.5d0*(v+v_new),f)
  v_new = rsf_update_V(tau_stick, sigma, f, Z)

  v = v_new
  f%theta = theta_new

  end subroutine rsf_solver

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

  function rsf_update_V(tau_stick,sigma,f,Z) result(v)
   
  double precision, dimension(:), intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(tau_stick)) :: v

  double precision :: tmp(size(tau_stick))

!  strength = -sigma*rsf_mu_no_direct(v,f) 
!  v = (tau_stick-strength)/Z

  select case(f%kind)
    case(1) 
      v = (tau_stick +sigma*rsf_mu_no_direct(f))/Z ! if v<0 will stop
      tmp = v -f%Vstar +sigma*f%a/Z
      v = 0.5d0*( tmp +sqrt(tmp*tmp +4d0*v*f%Vstar) )
      v = max(0d0,v)  ! arrest if v<0 
 
    case(2) 
      call IO_abort('rsf_solve: case 2 not implemented yet')
     !DEVEL: should call Newton-Raphson solver

    case(3) 
      call IO_abort('rsf_solve: case 3 not implemented yet')
     !DEVEL: should call Newton-Raphson solver

  end select

  end function rsf_update_V


!==================================================================
!         Newton-Raphson algorithm with bisection step         
! from  Numerical Recipes in Fortran 90 by Press et al.  (1996)
! finds a roots to 0=func_x bounded by [x1, x2] w/ error < x_acc

! generalized N-R method, must provide function (func_x) that gives
! the value of the function and the value of the derivative of the 
! function at a point 

  function nr_solver(func_x, x1, x2, x_acc)
 
  integer, parameter :: maxIteration=200
  real :: x_est, x_left, x_right, x_acc
  external :: func_x
  integer :: it
  real :: dfunc_dx, dx, dx_old, func_x, f_high, f_low, x_high, x_low
  real :: temp

  ! Find initial function estimates (dfunc_dx not needed yet)
  call func_x(x1, f_low, dfunc_dx)
  call func_x(x2, f_high, dfunc_dx)
 
  ! Safety checks:
  if ((f_low>0 .and. f_high>0) .or. (f_low<0 .and. f_high<0)) then
    call IO_abort('nr_solver - root must be bracketed! ')
 
  ! Lucky guesses: 
  if (f_low==0)  then
    x_est=x1
    return
  else if (f_high==0) then
    x_est=x2
    return
   
  ! Orient the search so that func_x(x_low) < 0
  else if (f_low<0) then 
    x_low=x1
    x_high=x2
  else 
    x_high=x1
    x_low=x2
  endif

  ! Initialize the guess and step sizes:
  x_est=.5*(x1+x2)
  ! The stepsize before last
  dx_old=abs(x2-x1)
  ! The last step:
  dx=dx_old
  
  call func_x(x_est,func_x,dfunc_dx)
  ! Loop over allowed iterations:
  do it=1,maxIteration
 
    ! Biset if N-R out of range, or not decreasing fast enough:
    if( (((x_est-x_high)*dfunc_dx-func_x)*((x_est-x_low)*dfunc_dx-func_x))>0  .or. abs(2*func_x)>abs(dx_old*dfunc_dx) ) then
      dxold=dx
      dx=0.5*(x_high-x_low)
      x_est=x_low+dx
      !  Check if change is negligible:
      if (x_low==x_est) return

    ! The Newton step is acceptable, move forward with algorithm:
    else
      dxold=dx
      dx=func_x/dfunc_dx
      temp=x_est
      x_est=x_est-dx
      !  Check if change is negligible:
      if (temp==x_est) return
    endif
  
    ! Check convergence criterion:
    if (abs(dx)<x_acc) return
  
    ! Evaluate function with new estimate of x
    call func_x(x_est,func_x,dfunc_dx)
    
    ! Redefine the bounds for the next loop
    if (f<0) then
      x_low=x_est
    else 
      x_high=x_est
    endif
  enddo
  call IO_abort('NR_Solver has exceeded the maximum iterations (200)')
  return
  end function

end module bc_dynflt_rsf
