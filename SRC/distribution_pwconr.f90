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
module distribution_pwconr

  implicit none
  private
 
  type pwconr_dist_type
    private
    integer :: NumZon ! number of zones
    double precision :: RefPnt(2) ! reference point
    double precision, pointer :: RadZon(:)  &! ext.radius of zone(i)
                                ,ValZon(:)   ! value inside zone(i)
  end type pwconr_dist_type

  public :: pwconr_dist_type,read_pwconr_dist,generate_pwconr_dist ,&
            destroy_pwconr_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_PWCONR
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Piecewise constant radial (2D) distribution.
! SYNTAX : &DIST_PWCONR num, ref /
!             r(1)  ... ...  r(num-1)
!          v(1) v(2) ... v(num-1) v(num)
!
! ARG: num      [int] [none] Number of radial zones (including outermost)
! ARG: ref      [dble(2)] [(0d0,0d0)] Reference point: center of radial zones
! ARG: r        [dble(num-1)] [none] External radius of zones: 
!                first zone R < r(1), second r(1) =< R < r(2), ...
!                last r(num-1) =< R 
! ARG: v        [dble(num)] [none] Values inside each zone
!
! END INPUT BLOCK

  subroutine read_pwconr_dist (data, file)

  use stdio, only: IO_abort
  type(pwconr_dist_type), intent(out) :: data
  integer, intent(in) :: file
  integer :: num
  double precision :: ref(2)
  NAMELIST / DIST_PWCONR / num, ref

  ref = 0.d0
  read (file,DIST_PWCONR)
  if (num<2) call IO_abort('read_pwconr_dist: needs more than 2 zones (num)')
  data%NumZon = num
  data%RefPnt = ref

  allocate( data%RadZon(num-1) )
  read (file,*) data%RadZon
  allocate( data%ValZon(num) )
  read (file,*) data%ValZon
 
  end subroutine read_pwconr_dist
!
!***********************************************************************
!

  subroutine generate_pwconr_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(pwconr_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  double precision :: rad
  integer :: i,izone

  do i =1,size(field)
    rad = sqrt( (coord(1,i)-par%RefPnt(1))**2 +(coord(2,i)-par%RefPnt(2))**2 ) 
   !NOTE: when the loop goes to its end, izone=par%NumZon
    do izone=1,par%NumZon-1; if ( rad <= par%RadZon(izone) ) exit; end do
    field(i) = par%ValZon(izone)
  end do
 
  end subroutine generate_pwconr_dist

!
!***********************************************************************
!
!! pwconr_dist_type destructor
subroutine destroy_pwconr_dist(d)
  type(pwconr_dist_type), pointer :: d
  deallocate(d%RadZon,d%ValZon)
  deallocate(d)
end subroutine destroy_pwconr_dist

end module distribution_pwconr
