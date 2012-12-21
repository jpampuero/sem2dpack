! SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module stf_gaussian

  implicit none
  private

  type STF_GAUSSIAN_type
    private
    double precision :: f0,ampli,onset
  end type STF_GAUSSIAN_type

  public :: STF_GAUSSIAN_type, STF_GAUSSIAN_read, STF_GAUSSIAN_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_GAUSSIAN
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Gaussian source time function  
!          stf(t) = ampli*exp[-(pi*f0*(t-onset))^2]
! SYNTAX : &STF_GAUSSIAN ampli, f0, onset /
!
! ARG: ampli    [dble] [1d0] Amplitude
! ARG: onset    [dble] [0d0] Delay time (s)
! ARG: f0       [dble] [1d0] Characteristic frequency bandwidth (Hz)  
!
! END INPUT BLOCK

subroutine STF_GAUSSIAN_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_GAUSSIAN_type), intent(out) :: stf
  integer, intent(in) :: iin

  real :: f0,onset,ampli
  
  NAMELIST / STF_GAUSSIAN / f0,onset,ampli

  f0=1d0
  onset = 0d0
  ampli = 1d0

  read(iin,STF_GAUSSIAN,END=500)

  if (f0 < 0d0) call IO_abort('STF_GAUSSIAN_read: f0 must be positive')

  stf%f0=f0
  stf%onset = onset
  stf%ampli = ampli

  if (echo_input) write(iout,200) f0,onset,ampli

  return

  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Gaussian',/5x, &
     'Characteristic frequency (Hz). . .(f0) =',EN12.3,/5x, &
     'Time delay (s) . . . . . . . . (onset) =',EN12.3,/5x, &
     'Amplitude. . . . . . . . . . . (ampli) =',EN12.3)
  
  500 call IO_abort('STF_GAUSSIAN_read: input block STF_GAUSSIAN not found')

end subroutine STF_GAUSSIAN_read


!=====================================================================
double precision function STF_GAUSSIAN_fun(stf,t)

  use constants, only: PI

  type(STF_GAUSSIAN_type), intent(in) :: stf
  double precision, intent(in) :: t

  double precision :: arg

  arg=PI*stf%f0*(t-stf%onset)
  arg=arg*arg
  STF_GAUSSIAN_fun = stf%ampli*exp(-arg)

end function STF_GAUSSIAN_fun

end module stf_gaussian
