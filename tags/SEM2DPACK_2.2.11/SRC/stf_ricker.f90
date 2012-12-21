! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module stf_ricker

  use constants, only: PI

  implicit none
  private

  type ricker_type
    private
    double precision :: f0,t0,ampli
  end type ricker_type

  public :: RICKER_type ,&
            RICKER_read ,&
            RICKER ,&
            RICKER_deriv ,&
            RICKER_int

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_RICKER
! GROUP  : SRC_TIMEFUNCTION
! PURPOSE: The Ricker wavelet is the second derivative of a gaussian.
! SYNTAX : &STF_RICKER ampli, f0, onset /
!
! ARG: ampli    [real] [1.] Signed amplitude of the central peak
! ARG: f0       [real >0] [0] Fundamental frequency (Hz).
!                 distribution: it has a peak at f0 and an exponential
!                 decay at high frequency. The cut-off high frequency is usually
!                 taken as fmax = 2.5 x f0. 
! ARG: onset    [real >1/f0] [0] Delay time (secs) with respect to the peak value. 
!
! NOTE   : The spectrum has a peak at f0 and decays exponentially at high 
!          frequencies. Beyond 2.5*f0 there is little energy, this is a 
!          recommended value for fmax.
! NOTE   : onset>1/f0 is needed to avoid a strong jump at t=0, which can cause
!          numerical oscillations. Ignore if using incident waves.
!
! END INPUT BLOCK

  subroutine RICKER_read(src,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(ricker_type), intent(out) :: src
  integer, intent(in) :: iin

  real :: f0,onset,ampli

  NAMELIST / STF_RICKER / f0,onset,ampli

  f0    = 0.
  onset = 0.
  ampli = 1.

  read(iin,STF_RICKER,END=100)

  if (f0 <= 0.) call IO_abort('RICKER_read: f0 must be positive')
  if (onset <= 1/f0) then
    !onset = 1./f0
    write(iout,*) '*** WARNING: RICKER_read: Onset time too small, ***'
!    write(iout,'(A,EN12.3)') '             will be reset to t0 = ',onset
  endif

  src%f0    = dble(f0)
  src%t0    = dble(onset)
  src%ampli = dble(ampli)

  if (echo_input) write(iout,200) f0,onset,ampli

  return
  
  100 call IO_abort('RICKER_read: STF_RICKER input block not found')
  200 format(5x, & 
     'Source time function . .(TimeFunction) = Ricker',/5x, &
     'Fundamental frequency (Hz) . . . .(f0) =',EN12.3,/5x, &
     'Time delay (s) . . . . . . . . (onset) =',EN12.3,/5x, &
     'Multiplying factor . . . . . . (ampli) =',EN12.3)
  
  end subroutine RICKER_read

!=====================================================================
  double precision function ricker(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker = - src%ampli * (1.d0-2.d0*arg)*exp(-arg)

  return
  end function ricker

!=====================================================================

  double precision function ricker_deriv(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker_deriv = (PI*src%f0)**2 *src%ampli *(t-src%t0)* (6.d0-4.d0*arg)*exp(-arg)

  return
  end function ricker_deriv
!
!=====================================================================

  double precision function ricker_int(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker_int = - src%ampli * (t-src%t0)*exp(-arg)

  return
  end function ricker_int

end module stf_ricker
