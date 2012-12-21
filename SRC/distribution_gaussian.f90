! SEM2DPACK version 2.3.2 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! Phone: (626) 395-3429
! Fax  : (626) 564-0715
! 
! http://www.seismolab.caltech.edu
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
module distribution_gaussian

  use stdio, only: IO_abort
  implicit none
  private
 
  type gaussian_dist_type
    private
    double precision :: x_0=0d0,z_0=0d0,lx=1d0,lz=1d0,level_0=0d0,ampli=1d0
    integer :: order=1
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
! PURPOSE: Bell shaped Gaussian 2D distribution 
! SYNTAX : &DIST_GAUSSIAN centered_at, length, offset, ampli, order /
!
! ARG: centered_at      [dble(2)] [0,0] Coordinates of the center point.
! ARG: length           [dble(2)] [1]   Characteristic lengths on each axis.
! ARG: offset           [dble] [0]      Background level.    
! ARG: ampli            [dble] [1]      Amplitude from background.
! ARG: order            [dble] [1]      Exponent
!
! END INPUT BLOCK


  subroutine read_gaussian_dist (d, file)

  type(gaussian_dist_type) :: d
  integer , intent(in) :: file

  double precision :: centered_at(2),length(2),offset,ampli
  integer :: order

  NAMELIST / DIST_GAUSSIAN / centered_at,length,offset,ampli,order

  centered_at = 0d0
  length = 1d0
  offset = 0d0
  ampli = 1d0
  order = 1

  read(file,DIST_GAUSSIAN,END=100)

  d%x_0 = centered_at(1)
  d%z_0 = centered_at(2)
  d%lx  = length(1)
  d%lz  = length(2)
  d%level_0 = offset
  d%ampli = ampli
  d%order = order

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

  field = d%level_0 + d%ampli*exp(- ((coord(1,:)-d%x_0)/d%lx)**(2*d%order) &
                                  - ((coord(2,:)-d%z_0)/d%lz)**(2*d%order) )
 
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
