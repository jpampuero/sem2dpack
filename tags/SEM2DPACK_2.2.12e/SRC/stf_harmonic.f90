! SEM2DPACK version 2.2.12e -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! ETH H�nggerberg HPP O 13.1
! CH-8093 Z�rich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
! 
! This software is freely available for academic research purposes. 
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
module stf_harmonic

  implicit none
  private

  type STF_HARMONIC_type
    private
    double precision :: ampli, f0
  end type STF_HARMONIC_type

  public :: STF_HARMONIC_type, STF_HARMONIC_read, STF_HARMONIC_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_HARMONIC
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Harmonic source time function f(t) = ampli*sin(2*pi*t*f0)
! SYNTAX : &STF_HARMONIC ampli, f0 /
!
! ARG: ampli    [dble] [0d0] Amplitude
! ARG: f0       [dble] [0d0] Frequency
!
! END INPUT BLOCK

subroutine STF_HARMONIC_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_HARMONIC_type), intent(out) :: stf
  integer, intent(in) :: iin

  real :: ampli,f0
  NAMELIST / STF_HARMONIC / ampli,f0

  ampli = 0d0
  f0 = 0d0

  read(iin,STF_HARMONIC,END=500)

  if (f0 <= 0d0) call IO_abort('STF_HARMONIC_read: f0 must be positive')
  if (ampli == 0d0) call IO_abort('STF_HARMONIC_read: ampli must be non zero')

  stf%ampli = ampli
  stf%f0 = f0

  if (echo_input) write(iout,200) ampli,f0

  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Harmonic',/5x, &
     'Amplitude. . . . . . . . . . . . . . . =',EN12.3,/5x, &
     'Frequency. . . . . . . . . . . . . . . =',EN12.3)
  
  500 call IO_abort('STF_HARMONIC_read: input block STF_HARMONIC not found')

end subroutine STF_HARMONIC_read


!=====================================================================
function STF_HARMONIC_fun(stf,t) result(fun)

  use constants, only : PI

  type(STF_HARMONIC_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun

  fun = stf%ampli* sin(2d0*PI*t*stf%f0)

end function STF_HARMONIC_fun

end module stf_harmonic