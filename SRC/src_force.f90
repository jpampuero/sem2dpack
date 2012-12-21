! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module src_force

! SRC_FORCE: collocated force source

  implicit none
  private

! dir(ndof) = direction vector, only used in P-SV
! iglob     = global node index 
  type so_force_type
    private
    double precision :: dir(2)
    integer :: iglob
  end type so_force_type

  public :: so_force_type,FORCE_read,FORCE_init,FORCE_add

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : SRC_FORCE
! GROUP  : SOURCE MECHANISM
! PURPOSE: Point force source
! SYNTAX : &SRC_FORCE angle /
!
! ARG: angle    [dble] [0d0]	For P-SV, the angle of the applied force, 
!                  in degrees, counterclockwise from Z-UP, e.g.: 
!                  (90 points left, 180 points down)
!                  For SH, angle is ignored.
!
! END INPUT BLOCK

  subroutine FORCE_read(so,iin)

  use stdio, only : IO_abort
  use echo, only : echo_input,iout
  use constants, only: PI

  type(so_force_type), intent(out) :: so
  integer, intent(in) :: iin

  double precision :: angle
  NAMELIST / SRC_FORCE / angle
  
  angle  = 0.d0
  read(iin,SRC_FORCE,END=100)
  if (echo_input) write(iout,200) angle
  angle  = angle*PI/180.d0
  so%dir = (/ -sin(angle), cos(angle) /) ! counterclockwise / UP 

  return

  100 call IO_abort('FORCE_read: SRC_FORCE input block not found')
  200 format( &
     5x, 'Source Type. . . . . . . . . . . . . . = Collocated Force', &
     /5x,'If P-SV: counterclockwise angle / up . = ',F0.2)

  end subroutine FORCE_read

!=====================================================================
!
  subroutine FORCE_init(so,iglob)

  type(so_force_type), intent(inout) :: so
  integer, intent(in) :: iglob

  so%iglob = iglob

  end subroutine

!=====================================================================
!
  subroutine FORCE_add(so,ampli,MxA)

  type(so_force_type), intent(inout) :: so
  double precision, intent(in) :: ampli
  double precision, intent(inout) :: MxA(:,:)

  if ( size(MxA,2)==1 ) then
    MxA(so%iglob,:) = MxA(so%iglob,:) + ampli
  else
    MxA(so%iglob,:) = MxA(so%iglob,:) + so%dir *ampli
  endif

  end subroutine FORCE_add


end module src_force
