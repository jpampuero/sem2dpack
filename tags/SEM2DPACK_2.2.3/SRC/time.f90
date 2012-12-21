! SEM2DPACK version 2.2.3 -- A Spectral Element Method tool for 2D wave propagation
!                            and earthquake source dynamics
! 
! Copyright (C) 2003 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics
! ETH Hönggerberg (HPP)
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 1 633 2197 (office)
! +41 1 633 1065 (fax)
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
    double precision :: dt,courant,time,total,alpha,beta,gamma
    double precision :: coef(4)
    integer :: nt,niter
  end type timescheme_type

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : TIME
! PURPOSE: Defines time integration scheme
! SYNTAX : &TIME NbSteps, Dt, Courant, TotalTime /
!
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
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : TIME_NEWMARK
! PURPOSE: Parameters of the explicit Newmark-alpha time scheme
! SYNTAX : &TIME_NEWMARK alpha|gamma, beta|rho /
!
! ARG: alpha    [dble] [0.5d0]
! ARG: beta     [dble] [0.5d0] 
! ARG: gamma    [dble] [1.d0] 
! ARG: rho      [dble] [1.d0] high frequencies are damped by a factor>=rho. 
!               The default is non-dissipative. Dissipation is limited however 
!               to rho>=0.5 . For max dissipation you should work close to
!               the stability limit (Courant around 0.6)
!
! NOTE: For second order schemes only two parameters need to be set:
!       (alpha OR gamma) AND (beta OR rho)
!
! END INPUT BLOCK

! ARG: niter    [int] [1] Number of iterations in predictor-multicorrector scheme
! NOTE: as of this version, niter>1 is not implemented

! NOTE: Some properties of the explicit Newmark-alpha as derived
!       from the analysis in HulChu96 of the more general EG-alpha scheme:
!
!       . eN-alpha = EG-alpha with alpha_m=0, alpha=1-alpha_f
!
!       . second order iff alpha+gamma=3/2 (from eq.21)
!
!       . let's call r=RhoP the spectral radius of the principal root
!         and RhoS that of the spurious root, 
!         in eN-alpha (alpha_m=0) we have RhoS = (1-r)/(2*r) (from eq.24)
!
!       . to optimize hf dissipation we require RhoS <= r
!         this can only be achieved when r>=0.5
!         (that's the limitation of eN-alpha)
!
!       . beta= 1 -alpha -r^2*(r-1)/[(1-alpha)*(1+r)^3] (from eq.25)
!
!       . this implies alpha+beta>1 (if alpha<1)
!
!       . the stability limit OmegaS=sqrt(-4*(1+r)^5/[6*r^2+9*r^4-15*r-34*r^3+r^5+1])
!         (from eq.28) ranges from approx. 1.87 at r=0.5 to 2.00 at r=1
!         The more dissipative schemes are also slightly more unstable
!
!       . the bifurcation limit OmegaB=sqrt[(1+r)^3/(2*r)] (from eq.27)
!         ranges from approx. 1.84 at r=0.5 to 2.00 at r=1

  subroutine TIME_read(t,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type(timescheme_type), intent(out) :: t

  double precision :: alpha,beta,gamma,dt,courant,TotalTime,rho
  integer :: NbSteps,niter
  
  NAMELIST / TIME / NbSteps,dt,courant,TotalTime
  NAMELIST / TIME_NEWMARK / alpha,beta,gamma,rho,niter
    
!-------------------------------------------------------------------------------

  NbSteps    = 0 
  dt       = 0.d0
  courant  = 0.5d0
  TotalTime    = 0.d0

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
    if (TotalTime > 0.d0) NbSteps = int(TotalTime/dt)
    TotalTime = dt*NbSteps
  endif

  if (echo_input) then
    write(iout,200) 
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
  niter    = 1

  read(iin,TIME_NEWMARK, END = 101)

  if (alpha < 0.d0 .or. alpha > 1.d0) call IO_abort('TIME_NEWMARK: alpha is out of range [0,1]')
  if (beta < 0.d0 .or. beta > 1.d0) call IO_abort('TIME_NEWMARK: beta is out of range [0,1]')
  if (gamma < 0.d0 .or. gamma > 1.d0) call IO_abort('TIME_NEWMARK: gamma is out of range [0,1]')
  if (rho<0.5d0 .or. rho>1.d0) call IO_abort('TIME_NEWMARK: rho is out of range [0.5,1]')
  if (niter < 1 ) call IO_abort('TIME_NEWMARK: niter must be positive')

  if (alpha/=0.5d0) gamma = 1.5d0-alpha
  if (gamma/=1.d0) alpha = 1.5d0-gamma
  if (rho/=1.d0) beta= 1.d0 -alpha -rho**2*(rho-1.d0)/((1.d0-alpha)*(1.d0+rho)**3)

101 continue
  if (echo_input) write(iout,300) alpha,beta,gamma,niter

  t%alpha = alpha
  t%beta  = beta
  t%gamma = gamma
  t%niter = niter

!-------------------------------------------------------------------------------

  return

  100 call IO_abort('TIME parameters not found')

  200 format(//' I t e r a t i o n   i n f o s '/1x,29('='),/)
  202 format(5x,'Number of time steps. . . . . (NbSteps) = ',I0)
  204 format(5x,'Time step increment . . . . . . . .(Dt) = ',EN12.3)
  206 format(5x,'Courant number. . . . . . . . (Courant) = ',F0.2)
  208 format(5x,'Total simulation duration . (TotalTime) = ',EN12.3)

  300   format(///' N e w m a r k   p a r a m e t e r s '/1x,35('=')//5x, &
      'First integration parameter . . . .(alpha) = ',F0.3,/5x, &
      'Second integration parameter. . . . (beta) = ',F0.3,/5x, &
      'Third time integration parameter . (gamma) = ',F0.3,/5x, &
      'Number of corrector iterations . . (niter) = ',I0)

  end subroutine TIME_read



!=======================================================================
!
!  Set time sequence parameters
!
  subroutine TIME_init(t)

  type(timescheme_type), intent(inout) :: t
  double precision, parameter :: one=1.d0,two=2.d0

!---- coefficients pour le schema de Newmark

  t%coef(1) = t%dt*t%dt*(one - two*t%beta)/two
  t%coef(2) = (one - t%gamma)*t%dt
  t%coef(3) = t%gamma*t%dt
  t%coef(4) = t%beta*t%dt*t%dt

  end subroutine TIME_init

end module time_evol
