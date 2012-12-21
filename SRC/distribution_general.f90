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
module distribution_general

!---  List here distribution modules
  use distribution_order0
  use distribution_gaussian
  use distribution_spline
  use distribution_linear
  use distribution_gradient
  use distribution_pwconr

  use stdio, only: IO_abort

  implicit none
  private

  type distribution_type
    private
    integer :: kind = 0
!--- List here distribution types
    type (order0_dist_type)  , pointer :: order0 => null()
    type (gaussian_dist_type), pointer :: gaussian => null()
    type (spline_dist_type)  , pointer :: spline => null()
    type (linear_dist_type)  , pointer :: linear => null()
    type (gradient_dist_type), pointer :: gradient => null()
    type (pwconr_dist_type), pointer :: pwconr => null()
  end type distribution_type

!--- List here the tags, must be different for each distribution type
  integer, parameter :: tag_order0   = 1 &
                       ,tag_gaussian = 2 &
                       ,tag_spline   = 3 &
                       ,tag_linear   = 4 &
                       ,tag_gradient = 5 &
                       ,tag_pwconr   = 6

  character(10) :: dist_name(6) = (/ 'ORDER0    ','GAUSSIAN  ','SPLINE    ','LINEAR    ', &
                                     'GRADIENT  ','PWCONR    '/)

  public :: distribution_type,DIST_read,DIST_generate,DIST_destructor

contains

!=======================================================================
!! Reads the distribution parameters from input File
subroutine DIST_read(d,name,iin)

  type(distribution_type), intent(out) :: d
  integer , intent(in) :: iin
  character(*), intent(in) :: name

  integer :: i

  do i=1,size(dist_name)
    if (dist_name(i) == name) d%kind = i
  enddo
  if (d%kind == 0) call IO_abort('DIST_read: unknown distribution name')

  select case (d%kind)
    case(tag_order0)
      allocate(d%order0)
      call read_order0_dist (d%order0,iin)
    case(tag_gaussian)
      allocate(d%gaussian)
      call read_gaussian_dist (d%gaussian,iin)
    case(tag_spline)
      allocate(d%spline)
      call read_spline_dist(d%spline,iin)
    case(tag_linear)
      allocate(d%linear)
      call read_linear_dist(d%linear,iin)
    case(tag_gradient)
      allocate(d%gradient)
      call read_gradient_dist(d%gradient,iin)
    case(tag_pwconr)
      allocate(d%pwconr)
      call read_pwconr_dist(d%pwconr,iin)
    case default
      call IO_abort( 'DIST_read: illegal kind')
  end select

end subroutine DIST_read

!=======================================================================
!! Generates a field using distribution parameters
subroutine DIST_generate (field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(distribution_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field

!  print *," generating dist type =",d%kind
  select case (d%kind)
    case(tag_order0)
      call generate_order0_dist(field,coord,d%order0)
    case(tag_gaussian)
      call generate_gaussian_dist(field,coord,d%gaussian)
    case(tag_spline)
      call generate_spline_dist(field,coord,d%spline)
    case(tag_linear)
      call generate_linear_dist(field,coord,d%linear)
    case(tag_gradient)
      call generate_gradient_dist(field,coord,d%gradient)
    case(tag_pwconr)
      call generate_pwconr_dist(field,coord,d%pwconr)
    case default
      call IO_abort( 'DIST_generate: illegal kind')
  end select

end subroutine DIST_generate

!=======================================================================
!! distribution destructor
subroutine DIST_destructor(d)

  type(distribution_type), intent(inout) :: d

!  print *,"destroying dist type = ", d%kind
  select case (d%kind)
    case(tag_order0)
      call destroy_order0_dist(d%order0)
    case(tag_gaussian)
      call destroy_gaussian_dist(d%gaussian)
    case(tag_spline)
      call destroy_spline_dist(d%spline)
    case(tag_linear)
      call destroy_linear_dist(d%linear)
    case(tag_gradient)
      call destroy_gradient_dist(d%gradient)
    case(tag_pwconr)
      call destroy_pwconr_dist(d%pwconr)
    case default
      call IO_abort( 'DIST_destroy: illegal kind')
  end select

end subroutine DIST_destructor

end module distribution_general
