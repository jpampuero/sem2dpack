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
module distribution_gaussian

  use stdio, only: IO_abort
  implicit none
  private
 
  type gaussian_dist_type
    private
    double precision :: x_0=0d0, z_0=0d0,lx=1d0,lz=1d0,level_0=0d0,ampli=1d0
  end type gaussian_dist_type

  public :: gaussian_dist_type,read_gaussian_dist,generate_gaussian_dist ,&
            destroy_gaussian_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_GAUSSIAN
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Bell shaped (Gaussian) 2D distribution
! SYNTAX : &DIST_GAUSSIAN centered_at, length, offset, ampli /
!
! ARG: centered_at      [dble(2)] [none] Coordinates of the center point.
! ARG: length           [dble(2)] [none] Characteristic lengths on each axis.
! ARG: offset           [dble] [none]    Background level.    
! ARG: ampli            [dble] [none]    Amplitude from background.
!
! END INPUT BLOCK


  subroutine read_gaussian_dist (d, file)

  type(gaussian_dist_type) :: d
  integer , intent(in) :: file

  double precision :: centered_at(2),length(2),offset,ampli

  NAMELIST / DIST_GAUSSIAN / centered_at,length,offset,ampli

  centered_at = 0d0
  length = 1d0
  offset = 0d0
  ampli = 1d0

  read(file,DIST_GAUSSIAN,END=100)

  d%x_0 = centered_at(1)
  d%z_0 = centered_at(2)
  d%lx  = length(1)
  d%lz  = length(2)
  d%level_0 = offset
  d%ampli = ampli

  return

  100 call IO_abort('read_gaussian_dist: DIST_GAUSSIAN parameters missing')

  end subroutine read_gaussian_dist
!
!***********************************************************************
!

  subroutine generate_gaussian_dist(field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(gaussian_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  field = d%level_0 + d%ampli*exp(- ((coord(1,:)-d%x_0)/d%lx)**2 &
                                  - ((coord(2,:)-d%z_0)/d%lz)**2 )
 
  end subroutine generate_gaussian_dist

!
!***********************************************************************
!
!! gaussian_dist_type destructor
subroutine destroy_gaussian_dist(d)
  type(gaussian_dist_type), pointer :: d
  deallocate(d)
end subroutine destroy_gaussian_dist

end module distribution_gaussian
