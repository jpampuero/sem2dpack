! SEM2DPACK version 2.3.3 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module distribution_order0

  implicit none
  private
 
  type order0_dist_type
    private
    integer :: x_nzones, z_nzones
    double precision, pointer, dimension(:) :: x_bound, z_bound
    double precision, pointer, dimension(:,:) :: val
  end type order0_dist_type

  public :: order0_dist_type,read_order0_dist,generate_order0_dist ,&
            destroy_order0_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_ORDER0
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Blockwise constant 2D distribution.
! SYNTAX : &DIST_ORDER0 xn, zn /
!          x(1) ...  x(xn-1)
!          z(1) ...  z(zn-1)
!          v(1,1)  ... v(xn,1)          
!            ...   ...   ...
!          v(1,zn) ... v(xn,zn)          
!
! ARG: xn       [int] [none] Number of zones along X
! ARG: zn       [int] [none] Number of zones along Z
! ARG: x        [dble(xn-1)] [none] Boundaries of X-zones: first zone X < x(1), 
!                second zone x(1) < X < x(2), ... , last zone x(xn-1) < X
! ARG: z        [dble(zn-1)] [none] Boundaries of Z-zones
! ARG: v        [dble(xn,zn)] [none] Values inside each zone
!
! END INPUT BLOCK

  subroutine read_order0_dist (data, file)

  type(order0_dist_type) :: data
  integer , intent(in) :: file

  integer :: xn,zn,i

  NAMELIST / DIST_ORDER0 / xn,zn

  read (file,DIST_ORDER0)

  data%x_nzones = xn
  data%z_nzones = zn

  allocate( data%x_bound(max(xn-1,1)) )
  allocate( data%z_bound(max(zn-1,1)) )
  if ( xn>1 ) read (file,*) data%x_bound
  if ( zn>1 ) read (file,*) data%z_bound

  allocate( data%val(xn,zn) )
  do i = 1, zn
    read (file,*) data%val(:,i) 
  end do
 
  end subroutine read_order0_dist
!
!***********************************************************************
!

  subroutine generate_order0_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(order0_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  integer :: nod

  do nod =1,size(field)
    field(nod) = par%val( zone(coord(1,nod), par%x_nzones, par%x_bound) &
                         ,zone(coord(2,nod), par%z_nzones, par%z_bound) )
  end do
 
  end subroutine generate_order0_dist

!
!***********************************************************************
!

  function zone (coord,nzones,bound)
 
  integer :: zone, nzones 
  double precision :: coord
  double precision :: bound(:)

  integer :: k

  zone = nzones
  if (nzones==1) return

  do k = 1, nzones - 1
    if ( coord < bound(k) ) exit
  end do
  zone = k

  end function zone

!
!***********************************************************************
!
!! order0_dist_type destructor
subroutine destroy_order0_dist(d)
  type(order0_dist_type), pointer :: d
  deallocate(d%x_bound,d%z_bound,d%val)
  deallocate(d)
end subroutine destroy_order0_dist

end module distribution_order0
