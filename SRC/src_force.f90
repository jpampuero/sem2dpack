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
module src_force

! SRC_FORCE: collocated force source

  implicit none
  private

  type double_pointer_type
    double precision, pointer :: P
  end type double_pointer_type

! dir(dof)   = direction vector
! field(dof) = points to the corresponding node of the field
!              where forces will be added
  type so_force_type
    private
    double precision :: dir(2)
    type(double_pointer_type), pointer :: field(:)
  end type so_force_type

  public :: so_force_type,FORCE_read,FORCE_init,FORCE_add

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : SRC_FORCE [source mechanism]
! PURPOSE: Point force source
! SYNTAX : &SRC_FORCE angle /
!
! ARG: angle    [dble] [0d0]   The angle of the force applied, in degrees,
!                              counterclockwise with respect to Z-UP: 90 points
!                              left, 180 points down
!
! END INPUT BLOCK

  subroutine FORCE_read(so,iin)

  use stdio, only : IO_abort
  use echo, only : echo_input,iout

  type(so_force_type), intent(out) :: so
  integer, intent(in) :: iin

  double precision, parameter :: pi = 3.141592653589793d0
  double precision :: angle
  NAMELIST / SRC_FORCE / angle
  
  angle  = 0.d0
  read(iin,SRC_FORCE,END=100)
  if (echo_input) write(iout,200) angle
  angle  = angle*pi/180.d0
  so%dir = (/ -sin(angle), cos(angle) /) ! counterclockwise / UP 

  return

  100 call IO_abort('FORCE_read: SRC_FORCE input block not found')
  200 format( &
     5x, 'Source Type. . . . . . . . . . . . . . = Collocated Force', &
     /5x,'Angle from vertical up . . . . . (deg) =',F0.2)

  end subroutine FORCE_read

!=====================================================================
!
  subroutine FORCE_init(so,field,iglob)

  type(so_force_type), intent(inout) :: so
  double precision, pointer  :: field(:,:)
  integer, intent(in) :: iglob

  integer :: d,dof

  dof = size(field,2)
  allocate(so%field(dof))
  do d=1,dof
    so%field(d)%P => field(iglob,d)
  enddo
 
  end subroutine

!=====================================================================
!
  subroutine FORCE_add(so,ampli)

  type(so_force_type), intent(inout) :: so
  double precision, intent(in) :: ampli

  integer :: d

  do d=1,size(so%field)
    so%field(d)%P = so%field(d)%P + so%dir(d) *ampli
  enddo

  end subroutine FORCE_add


end module src_force
