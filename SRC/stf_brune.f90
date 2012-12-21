! SEM2DPACK version 2.3.6 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://web.gps.caltech.edu/~ampuero/
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
module stf_brune

! Brune (1970)'s assumes the seismic moment rate is 
!   Mo*wc^2*t*exp(-wc*t)
! This implies the seismic moment is 
!   Mo*( 1 - (1+wc*t)*exp(-wc*t) )

  implicit none
  private

  type STF_BRUNE_type
    private
    double precision :: ampli, fc
  end type STF_BRUNE_type

  public :: STF_BRUNE_type, STF_BRUNE_read, STF_BRUNE_fun

contains

!=====================================================================
!
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_BRUNE
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Brune (1970)'s model with omega-squared spectral fall-off:
!            stf(t) = ampli*( 1 - (1+2*pi*fc*t)*exp(-2*pi*fc*t) )
! SYNTAX : &STF_BRUNE ampli, fc /
!
! ARG: ampli    [dble] [1d0] Amplitude (usually the seismic moment)
! ARG: fc       [dble] [1d0] Corner frequency (Hz)
!
! END INPUT BLOCK

subroutine STF_BRUNE_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_BRUNE_type), intent(out) :: stf
  integer, intent(in) :: iin

  double precision :: ampli,fc
  NAMELIST / STF_BRUNE / ampli,fc

  ampli = 1d0
  fc = 1d0

  read(iin,STF_BRUNE,END=500)

  if (fc < 0d0) call IO_abort('STF_BRUNE_read: fc must be positive')

  stf%ampli = ampli
  stf%fc = fc

  if (echo_input) write(iout,200) ampli,fc
  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Brune',/5x, &
     'Amplitude. . . . . . . . . . . . . . . =',EN12.3,/5x, &
     'Corner frequency (Hz). . . . . . . . . =',EN12.3,/5x )
  
  500 call IO_abort('STF_BRUNE_read: input block STF_BRUNE not found')

end subroutine STF_BRUNE_read


!=====================================================================
function STF_BRUNE_fun(stf,t) result(fun)

  use constants, only : pi

  type(STF_BRUNE_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun, arg

  arg = 2*pi*stf%fc *max(t,0d0)
  fun = stf%ampli *( 1d0 - (1d0+arg)*exp(-arg) )

end function STF_BRUNE_fun

end module stf_brune
