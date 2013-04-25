module time_evol

  implicit none

  type timescheme_type
    !private !devel: main and solver need it public
    character(12) :: kind
    double precision :: dt,courant,time,total,alpha,beta,gamma,Omega_max
    double precision, dimension(:), pointer :: a,b 
    integer :: nt,nstages
  end type timescheme_type

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : TIME
! PURPOSE: Defines time integration scheme
! SYNTAX : &TIME kind, {Dt or Courant}, {NbSteps or TotalTime} /
!          Possibly followed by a TIME_XXXX block.
!
! ARG: kind      [char*10] ['leapfrog'] Type of scheme:
!                'newmark'       Explicit Newmark
!                'HHT-alpha'     Explicit HHT-alpha
!                'leapfrog'      Central difference
!                'symp_PV'       Position Verlet
!                'symp_PFR'      Position Forest-Ruth (4th order)
!                'symp_PEFRL'    Extended PFR (4th order)
!		 'quasi-static'  Quasi-static (Kaneko, et al, 2011)
! ARG: Dt        [dble] [none] Timestep (in seconds)
! ARG: Courant   [dble] [0.5d0] the maximum value of the Courant-Friedrichs-Lewy 
!                stability number (CFL), defined as
!                  CFL = Dt*wave_velocity/dx 
!                where dx is the distance between GLL nodes. Tipically CFL<= 0.5
! ARG: NbSteps   [int] [none] Total number of timesteps
! ARG: TotalTime [int] [none] Total duration (in seconds)
!
! NOTE   : The leap-frog scheme is recommended for dynamic faults. It is equivalent 
!          to the default Newmark scheme (beta=0, gamma=1/2). However it is 
!          faster and requires less memory.
!
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : TIME_NEWMARK
! GROUP  : TIME SCHEMES
! PURPOSE: Explicit Newmark time integration scheme 
! SYNTAX : &TIME_NEWMARK gamma, beta /
!
! ARG: beta     [dble] [0d0] First Newmark parameter.
!               If beta=0 the scheme is fully explicit (the update of
!               displacement depends only on the last value of acceleration),
!               otherwise it is a single-predictor-corrector scheme
! ARG: gamma    [dble] [0.5d0] Second Newmark parameter.
!               Second order requires gamma=1/2.
!
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : TIME_HHTA
! GROUP  : TIME SCHEMES
! PURPOSE: Explicit HHT-alpha time integration scheme, second order 
! SYNTAX : &TIME_HHTA alpha, rho /
!
! ARG: alpha    [dble] [0.5d0] Parameter in the HHT-alpha method. Values in [0,1].
!               Defined here as 1 + HHT's original definition of alpha.
!               When alpha=1 it reduces to second order explicit Newmark
!               (beta=0, gamma=0.5).
! ARG: rho      [dble] [0.5d0] Minimum damping factor for high frequencies.
!               Values in [0.5,1]. Rho=1 is non-dissipative.
!
! NOTE: We consider only second order schemes, for which alpha+gamma=3/2
!       If  alpha<1, Newmark's beta is related to the HHT parameters by
!         beta = 1 -alpha -rho^2*(rho-1)/[(1-alpha)*(1+rho)^3]
!       If alpha=1, we set rho=1 (beta=0, gamma=0.5)
! 
! NOTE: Dissipative schemes (rho<1) require slightly smaller Courant number
!       (0.56 for rho=0.5, compared to 0.6 for rho=1)
!
! NOTE:	This is an explicit version of the HHT-alpha scheme of
!         H.M. Hilber, T.J.R. Hughes and R.L. Taylor (1977) "Improved numerical 
!         dissipation for time integration algorithms in structural dynamics" 
!         Earthquake Engineering and Structural Dynamics, 5, 283-292
!       implemented with a slightly different definition of alpha (1+original).
!      	Its properties can be derived from the EG-alpha scheme of
!         G.M. Hulbert and J. Chung (1996) "Explicit time integration 
!         algorithms for structural dynamics with optimal numerical dissipation"
!         Comp. Methods Appl. Mech. Engrg. 137, 175-188
!       by setting alpha_m=0 and alpha=1-alpha_f.
!
! END INPUT BLOCK

!
!       . second order iff alpha+gamma=3/2 (from eq.21)
!
!       . let's call r=RhoP the spectral radius of the principal root
!         and RhoS that of the spurious root, 
!         in HHT-alpha: RhoS = (1-r)/(2*r) (from eq.24)
!
!       . to optimize high-frequency dissipation we require RhoS <= r
!         this can only be achieved when r>=0.5
!         (that's the limitation of HHT-alpha)
!
!       . beta= 1 -alpha -r^2*(r-1)/[(1-alpha)*(1+r)^3] (from eq.25)
!
!       . conservative schemes (r=1) are obtained when beta + alpha = 1
!
!       . the stability limit OmegaS=sqrt(-4*(1+r)^5/[6*r^2+9*r^4-15*r-34*r^3+r^5+1])
!         (from eq.28) ranges from approx. 1.87 at r=0.5 to 2.00 at r=1
!         The more dissipative schemes are also slightly more unstable
!
!       . the bifurcation limit OmegaB=sqrt[(1+r)^3/(2*r)] (from eq.27)
!         ranges from approx. 1.84 at r=0.5 to 2.00 at r=1
!

  subroutine TIME_read(t,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type(timescheme_type), intent(out) :: t

  double precision :: alpha,beta,gamma,dt,courant,TotalTime,rho &
                     ,xi,lambda,chi,theta
  integer :: NbSteps,n
  character(12) :: kind
  
  NAMELIST / TIME / kind,NbSteps,dt,courant,TotalTime
  NAMELIST / TIME_NEWMARK / beta,gamma
  NAMELIST / TIME_HHTA / alpha,rho
    
!-------------------------------------------------------------------------------

  kind = 'leapfrog'
  NbSteps  = 0 
  dt       = 0.d0
  courant  = 0.5d0
  TotalTime = 0.d0

  rewind(iin)
  read(iin,TIME,END = 100) 

  if (NbSteps < 0) call IO_abort('TIME: NbSteps must be positive')
  if (dt < 0.d0) call IO_abort('TIME: Dt must be positive')
  if (courant < 0.d0 .or. Courant > 0.6) &
    call IO_abort('TIME: Courant out of range [0,0.6]')
  if (TotalTime < 0.d0) call IO_abort('TIME: TotalTime must be positive')

 ! The user can set the total duration (secs) or the number of NbSteps.
  if (NbSteps*TotalTime /= 0.d0)  &
    call IO_abort('TIME: bad combination of settings, NbSteps or TotalTime')

  if (dt > 0.d0) then
    if (TotalTime > 0.d0) NbSteps = ceiling(TotalTime/dt)
    TotalTime = dt*NbSteps
  endif

  if (echo_input) then
    write(iout,200) 
    write(iout,201) kind
    if (NbSteps > 0) then
      write(iout,202) NbSteps
    else
      write(iout,203)
    endif
    if (dt > 0.d0) then
      write(iout,204) dt
    else
      write(iout,205)
      write(iout,206) courant
    endif
    if (TotalTime > 0) then
      write(iout,208) TotalTime
    else
      write(iout,209)
    endif
  endif

  t%nt      = NbSteps
  t%dt      = dt
  t%courant = courant
  t%total   = TotalTime
  t%kind    = kind

!-------------------------------------------------------------------------------

! old default was: alpha=1/2, beta=1/2, gamma=1, rho=1
! new default is equivalent to leapfrog: alpha=1, beta=0, gamma=1/2, rho=1.
  alpha    = 1d0
  beta     = 0d0
  gamma    = 0.5d0

  select case (kind)
   
   case ('quasi-static')
    t%Omega_max = huge(1d0)

   case ('leapfrog')
    t%Omega_max = 2d0

   case ('newmark')
    read(iin,TIME_NEWMARK, END = 101)
  
    if (beta < 0.d0 .or. beta > 1.d0) call IO_abort('TIME_NEWMARK: beta is out of range [0,1]')
    if (gamma < 0.d0 .or. gamma > 1.d0) call IO_abort('TIME_NEWMARK: gamma is out of range [0,1]')
  
101 continue
    if (echo_input) write(iout,300) beta,gamma
    t%Omega_max = sqrt(2d0/gamma)


   case ('HHT-alpha')

    alpha = 0.5d0
    rho = 0.5d0

    read(iin,TIME_HHTA, END = 102)

    if (alpha < 0.d0 .or. alpha > 1.d0) call IO_abort('TIME_HHTA: alpha is out of range [0,1]')
    if (rho<0.5d0 .or. rho>1.d0) call IO_abort('TIME_HHTA: rho is out of range [0.5,1]')

    gamma = 1.5d0-alpha
    if (alpha/=1.d0) then
      beta  = 1.d0 -alpha -rho**2*(rho-1.d0)/((1.d0-alpha)*(1.d0+rho)**3)
    else
      beta  = 0d0   ! back to explicit newmark
    endif

102 continue
    if (echo_input) write(iout,310) alpha,rho
    t%Omega_max = sqrt( -4d0*(1d0+rho)**5 &
                   /( rho**5 +9d0*rho**4 -34d0*rho**3 +6d0*rho**2 -15d0*rho +1d0 ) )

   case ('symp_PV')
    n = 1
    allocate(t%a(n+1), t%b(n))
    t%nstages = n
    t%a(1) = 0.5d0
    t%a(2) = t%a(1)
    t%b(1) = 1d0
    t%Omega_max = 2d0

   case ('symp_PFR')
    n = 3
    allocate(t%a(n+1), t%b(n))
    t%nstages = n
    theta = 1d0/(2d0-2d0**(1d0/3d0))
    t%a(1) = theta/2d0
    t%a(2) = (1d0-theta)/2d0
    t%a(3) = t%a(2)
    t%a(4) = t%a(1)
    t%b(1) = theta
    t%b(2) = 1d0-2d0*theta
    t%b(3) = t%b(1)
    t%Omega_max = 1.5734d0

   case ('symp_PEFRL')
    n = 4
    allocate(t%a(n+1), t%b(n))
    t%nstages = n
    xi = 0.1786178958448091d0 ;
    lambda = -0.2123418310626054d0 ;
    chi = -0.06626458266981849d0;
    t%a(1) = xi
    t%a(2) = chi
    t%a(3) = 1d0-2d0*(chi+xi)
    t%a(4) = t%a(2)
    t%a(5) = t%a(1)
    t%b(1) = 0.5d0-lambda
    t%b(2) = lambda
    t%b(3) = t%b(2)
    t%b(4) = t%b(1)
    t%Omega_max = 2.97633d0

   case default
    call IO_abort('TIME: unknown kind')

  end select

  t%alpha = alpha
  t%beta  = beta
  t%gamma = gamma

!-------------------------------------------------------------------------------

  return

  100 call IO_abort('TIME parameters not found')

  200 format(//' T i m e   i n t e g r a t i o n'/1x,31('='),/)
  201 format(5x,'Scheme. . . . . . . . . . . . . .(kind) = ',A)
  202 format(5x,'Number of time steps. . . . . (NbSteps) = ',I0)
  203 format(5x,'Number of time steps. . . . . (NbSteps) = will be set later')
  204 format(5x,'Time step increment . . . . . . . .(Dt) = ',EN12.3)
  205 format(5x,'Time step increment . . . . . . . .(Dt) = will be set later')
  206 format(5x,'Courant number. . . . . . . . (Courant) = ',F0.2)
  208 format(5x,'Total simulation duration . (TotalTime) = ',EN12.3)
  209 format(5x,'Total simulation duration . (TotalTime) = will be set later')

  300   format(///' N e w m a r k   p a r a m e t e r s '/1x,35('=')//5x, &
      'First integration parameter . . . . (beta) = ',F0.3,/5x, &
      'Second time integration parameter .(gamma) = ',F0.3)

  310   format(///' H H T - a l p h a   p a r a m e t e r s '/1x,39('=')//5x, &
      'Force collocation parameter . . . .(alpha) = ',F0.3,/5x, &
      'High-frequency damping factor . . . .(rho) = ',F0.3)

  end subroutine TIME_read



!=======================================================================
!
!  Set time integration parameters
!
  subroutine TIME_init(t,grid_cfl)

  use echo, only : echo_check,iout
  use stdio, only : IO_warning

  type(timescheme_type), intent(inout) :: t
  double precision, intent(in) :: grid_cfl

  double precision :: critical_CFL

 ! Check the Courant number or set the timestep:
  if (t%dt > 0.d0) then
    t%courant =  grid_cfl*t%dt 
  else
    t%dt      =  t%courant/grid_cfl
   ! Set the total duration or the number of timesteps
    if (t%total > 0.d0) t%nt = ceiling(t%total/t%dt)
    t%total = t%nt*t%dt
  endif

  if (echo_check) then

    write(iout,*) 
    write(iout,103) ' T i m e   s o l v e r' 
    write(iout,103) ' ====================='
    write(iout,*) 
    write(iout,101) '    Time step (secs)      = ',t%dt
    write(iout,104) '    Number of time steps  = ',t%nt
    write(iout,101) '    Total duration (secs) = ',t%total
    write(iout,101) '    Courant number        = ',t%courant
    write(iout,*) 
    write(iout,102) '    STABILITY:  CFL number               = ',grid_cfl*t%dt
    write(iout,*) 

  endif

! In 1D SEM, the maximum angular frequency of an element is 
!       wmax ~ 7/3 * c/dx *sqrt(D)
! where c  = wave velocity 
!       dx = minimal GLL spacing ~ 4 h/(ngll^2-1)
!       h  = element size
!       D  = dimension (here D=2)
! For stability we need:
!       dt*wmax < Omega_max
! where Omega_max depends on the time scheme (see above)
! The CFL number is defined as (see init.f90)
!       CFL = c*dt/dx
! so the stability condition is
!       CFL < Omega_max * 3/7/sqrt(D) 
! Example: 2D leapfrog CFL <~ 0.6
  critical_CFL = t%Omega_max * 3d0/7d0/sqrt(2d0) 

  if (t%courant > critical_CFL) then
    write(iout,*)
    write(iout,103) '*******************************************************'
    write(iout,102) '** WARNING: Courant number too high = ',t%courant,'   **'
    write(iout,102) '**          Numerical instability is expected !      **' 
    write(iout,102) '** Try a value smaller than ',critical_CFL,'             **'
    write(iout,102) '** or a timestep smaller than ', &
                    t%dt*critical_CFL/t%courant ,'           **'
    write(iout,103) '*******************************************************'
    write(iout,*)
    call IO_warning()
  endif

  return

  101 format(A,EN12.3)
  102 format(A,EN12.3,A)
  103 format(A)
  104 format(A,I0)

  end subroutine TIME_init

!=======================================================================
  logical function TIME_needsAlphaField(t) 
  type(timescheme_type), intent(in) :: t
  TIME_needsAlphaField = t%kind=='HHT-alpha'
  end function TIME_needsAlphaField
  
!=======================================================================
  double precision function TIME_getTimeStep(t)
  type(timescheme_type), intent(in) :: t
  TIME_getTimeStep = t%dt
  end function TIME_getTimeStep

!=======================================================================
  integer function TIME_getNbTimeSteps(t)
  type(timescheme_type), intent(in) :: t
  TIME_getNbTimeSteps = t%nt
  end function TIME_getNbTimeSteps

!=======================================================================
  double precision function TIME_getTime(t)
  type(timescheme_type), intent(in) :: t
  TIME_getTime = t%time
  end function TIME_getTime

!=======================================================================
! Coefficients of corrector phase (see solver.f90)
!  vnew = vpredictor + coefA2V *anew   
!  dnew = dpredictor + coefA2D *anew

  function TIME_getCoefA2D(t) result(c)

  use stdio, only : IO_abort

  type(timescheme_type), intent(in) :: t
  double precision :: c

  select case (t%kind)
    case ('newmark','HHT-alpha'); c = t%beta * t%dt**2
    case ('leapfrog','quasi-static'); c = 0d0
    case default
      c=0d0
      call IO_abort('TIME_getCoefA2D: unknown time scheme')
  end select

  end function TIME_getCoefA2D

!-----------------------------------------------------------------------
  function TIME_getCoefA2V(t) result(c)

  use stdio, only : IO_abort

  type(timescheme_type), intent(in) :: t
  double precision :: c

  select case (t%kind)
    case ('newmark','HHT-alpha'); c = t%gamma * t%dt
    case ('leapfrog','quasi-static'); c = t%dt
    case default
      c=0d0
      call IO_abort('TIME_getCoefA2D: unknown time scheme')
  end select

  end function TIME_getCoefA2V

!-----------------------------------------------------------------------
! Get the coefficient in 
!   v_rhs = v_rhs_pre + coefficient*a
! where v_rhs is the velocity on the r.h.s. of the equation M*a = ...
! and v_rhs_pre its predicted value
! 
  function TIME_getCoefA2Vrhs(t) result(c)

  use stdio, only : IO_abort

  type(timescheme_type), intent(in) :: t
  double precision :: c

  select case (t%kind)
    case ('newmark','HHT-alpha'); c = t%alpha * TIME_getCoefA2V(t)
    case ('leapfrog','quasi-static'); c = 0.5d0 * TIME_getCoefA2V(t)
 ! NOTE: for the leapfrog time scheme, v_rhs = v_(n+1)
 !	 but the velocity field is stored at n+1/2, 
 !       v_rhs = v_(n+1/2) + 1/2*dt*a_(n+1)
 !       v_rhs_pre = v_(n+1/2)
    case default
      c=0d0
      call IO_abort('TIME_getCoefA2D: unknown time scheme')
  end select

  end function TIME_getCoefA2Vrhs

end module time_evol
