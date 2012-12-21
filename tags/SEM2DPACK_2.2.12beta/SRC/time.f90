! SEM2DPACK version 2.2.12beta -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                             with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics Group
! ETH Hönggerberg HPP O 13.1
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
! 
! This software is freely available for scientific research purposes. 
! If you use this software in writing scientific papers include proper 
! attributions to its author, Jean-Paul Ampuero.
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
! 
module time_evol

  implicit none

  type timescheme_type
    character(10) :: kind
    double precision :: dt,courant,time,total,alpha,beta,gamma,rho,Omega_max
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
! SYNTAX : &TIME kind, NbSteps, Dt, Courant, TotalTime /
!
! ARG: kind     [char*10] ['leapfrog'] Type of scheme:
!			'newmark'	Newmark-alpha
!			'leapfrog'	Central difference
!			'symp_PV'		Position Verlet
!			'symp_PFR'	Position Forest-Ruth (4th order)
!			'symp_PEFRL'	Extended PFR (4th order)
! ARG: NbSteps  [int] [none] Number of timesteps to be performed
! ARG: Dt       [dble] [none] Amplitude of the timestep
! ARG: Courant  [dble] [0.5d0] Courant stability number: the maximum ratio
!               Dt*wave_velocity/dx where dx is the inter-GLL node distance
!               Tipically <= 0.5
! ARG: TotalTime[int] [none] Total duration (in seconds) of simulation
!
! NOTE   :      Not all combinations of parameters need to be set at once.
!               You can set the total duration (secs) or the number of steps.
!               You can set the timestep or the Courant number (or use default).
! 
! NOTE:		The leap-frog scheme is equivalent to the Newmark scheme
!		with alpha=1, beta=0, gamma=1/2.
!		However it is faster and requires less memory.
!		Dynamic faults require this scheme.
!
! END INPUT BLOCK

!		with alpha=1/2, beta=1/2, gamma=1.

! BEGIN INPUT BLOCK
!
! NAME   : TIME_NEWMARK
! PURPOSE: Parameters of the explicit Newmark or HHT-alpha time scheme 
! SYNTAX : &TIME_NEWMARK alpha|gamma, beta|rho /
!
! ARG: beta     [dble] [0.5d0] The algorithm is fully explicit if beta=0
!               otherwise it is a single-predictor-corrector scheme
! ARG: gamma    [dble] [1.d0] 
! ARG: alpha    [dble] [0.5d0] parameter in the Hilber-Hughes-Taylor method
!               Actually, here alpha = 1 + their original definition of alpha
! ARG: rho      [dble] [1.d0] high frequencies are damped by a factor>=rho. 
!               The default is non-dissipative. Dissipation is limited however 
!               to rho>=0.5 . For max dissipation you should work close to
!               the stability limit (Courant around 0.56 for rho=0.5).
!
! NOTE: For second order schemes only two parameters need to be set:
!       (alpha OR gamma) AND (beta OR rho)
!
! NOTE: Dissipative schemes (0.5<=rho<1) are slightly more unstable,
!       i.e. they require slightly smaller Courant number
!       (0.56 for rho=0.5, compared to 0.6 for rho=1)
!
! END INPUT BLOCK

! NOTE: Some properties of the explicit HHT-alpha scheme can be derived
!       from the more general EG-alpha scheme of
!       	G.M. Hulbert and J. Chung (1996) "Explicit time integration 
!		algorithms for structural dynamics with optimal numerical dissipation"
!		Comp. Methods Appl. Mech. Engrg. 137, 175-188
!
!       . our HHT-alpha = EG-alpha with alpha_m=0, alpha=1-alpha_f
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
  character(10) :: kind
  
  NAMELIST / TIME / kind,NbSteps,dt,courant,TotalTime
  NAMELIST / TIME_NEWMARK / alpha,beta,gamma,rho
    
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
      write(iout,'(A)') '     The number of steps will be set later'
    endif
    if (dt > 0.d0) then
      write(iout,204) dt
    else
      write(iout,'(A)') '     The timestep will be set later'
      write(iout,206) courant
    endif
    if (TotalTime > 0) then
      write(iout,208) TotalTime
    else
      write(iout,'(A)') '     The total duration will be set later'
    endif
  endif

  t%nt      = NbSteps
  t%dt      = dt
  t%courant = courant
  t%total   = TotalTime

!-------------------------------------------------------------------------------

  alpha    = 0.5d0
  beta     = 0.5d0
  gamma    = 1.d0
  rho      = 1.d0

  select case (kind)

   case ('leapfrog')
    t%Omega_max = 2d0

   case ('newmark')
    read(iin,TIME_NEWMARK, END = 101)
  
    if (alpha < 0.d0 .or. alpha > 1.d0) call IO_abort('TIME_NEWMARK: alpha is out of range [0,1]')
    if (beta < 0.d0 .or. beta > 1.d0) call IO_abort('TIME_NEWMARK: beta is out of range [0,1]')
    if (gamma < 0.d0 .or. gamma > 1.d0) call IO_abort('TIME_NEWMARK: gamma is out of range [0,1]')
    if (rho<0.5d0 .or. rho>1.d0) call IO_abort('TIME_NEWMARK: rho is out of range [0.5,1]')
  
    if (alpha/=0.5d0) gamma = 1.5d0-alpha
    if (gamma/=1.d0)  alpha = 1.5d0-gamma
    if (rho/=1.d0)    beta  = 1.d0 -alpha -rho**2*(rho-1.d0)/((1.d0-alpha)*(1.d0+rho)**3)

101 continue
    if (echo_input) write(iout,300) alpha,beta,gamma
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

  end select

  t%alpha = alpha
  t%beta  = beta
  t%gamma = gamma
  t%rho   = rho
  
  t%kind = kind


!-------------------------------------------------------------------------------

  return

  100 call IO_abort('TIME parameters not found')

  200 format(//' T i me   i n t e g r a t i o n'/1x,30('='),/)
  201 format(5x,'Scheme. . . . . . . . . . . . . .(kind) = ',A)
  202 format(5x,'Number of time steps. . . . . (NbSteps) = ',I0)
  204 format(5x,'Time step increment . . . . . . . .(Dt) = ',EN12.3)
  206 format(5x,'Courant number. . . . . . . . (Courant) = ',F0.2)
  208 format(5x,'Total simulation duration . (TotalTime) = ',EN12.3)

  300   format(///' N e w m a r k   p a r a m e t e r s '/1x,35('=')//5x, &
      'First integration parameter . . . .(alpha) = ',F0.3,/5x, &
      'Second integration parameter. . . . (beta) = ',F0.3,/5x, &
      'Third time integration parameter . (gamma) = ',F0.3)

  end subroutine TIME_read



!=======================================================================
!
!  Set time sequence parameters
!
  subroutine TIME_init(t,grid_cfl)

  use echo, only : echo_check,iout

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

  endif


  return

  101 format(A,EN12.3)
  102 format(A,EN12.3,A)
  103 format(A)
  104 format(A,I0)

  end subroutine TIME_init

end module time_evol
