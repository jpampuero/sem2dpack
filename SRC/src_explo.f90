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
module src_explo

! SRC_EXPLO: explosion source

  implicit none
  private

  type double_pointer_type
    double precision, pointer :: P => null()
  end type double_pointer_type

  type so_explo_type
    private
    double precision, pointer :: coef(:,:) => null()
    type(double_pointer_type), pointer :: field(:,:) => null()
  end type so_explo_type

  public :: so_explo_type,EXPLO_read,EXPLO_init,EXPLO_add

contains

!=====================================================================
!
  subroutine EXPLO_read()
  
  use echo, only : echo_input,iout

  if (echo_input) write(iout,200)
  return
  200   format(5x,'Source Type. . . . . . . . . . . . . . = Explosion')

  end subroutine EXPLO_read

!=====================================================================
!-- Define the working arrays 
!   Spatial distribution: Dirac (schema en croix)
!     so%coef   = coefficients
!     so%field  = points to the corresponding nodes of the field
!                 where forces will be added
!
!   WARNING: the source must not be on an element edge
!
  subroutine EXPLO_init(so,igll,jgll,element,grid,field)

  use spec_grid, only : sem_grid_type,SE_inquire

  type(so_explo_type), intent(inout) :: so
  type(sem_grid_type), intent(in) :: grid
  double precision, pointer :: field(:,:)
  integer, intent(in) :: igll,jgll,element
  
  double precision :: a(grid%ngll,grid%ngll,grid%ndime)
  double precision :: jac_inv(2,2)
  integer :: d,i,j,ip,src_npoin,ngll

  src_npoin = 2*grid%ngll-1

  allocate(so%coef(src_npoin,grid%ndime))
  allocate(so%field(src_npoin,grid%ndime))

  !-- get the inverse jacobian
  call SE_inquire(grid,element=element,igll=igll,jgll=jgll &
                 ,inv_jacobian=jac_inv)

  a(:,jgll,1) = jac_inv(1,1) * grid%hprime(:,igll)
  a(:,jgll,2) = jac_inv(1,2) * grid%hprime(:,igll)
  a(igll,:,1) = jac_inv(2,1) * grid%hprime(:,jgll)
  a(igll,:,2) = jac_inv(2,2) * grid%hprime(:,jgll)
  a(igll,jgll,:) = ( jac_inv(1,:) + jac_inv(2,:) ) * grid%hprime(igll,jgll)

  ip = 0
  do j=1,grid%ngll
  do i=1,grid%ngll
    if (i /= igll .AND. j /= jgll) cycle
    ip = ip +1
    so%field(ip,1)%P => field(grid%ibool(i,j,element),1)
    so%field(ip,2)%P => field(grid%ibool(i,j,element),2)
    so%coef(ip,:) = a(i,j,:)
  enddo
  enddo

  end subroutine EXPLO_init

!=====================================================================

  subroutine EXPLO_add(so,ampli)

  type(so_explo_type), intent(inout) :: so
  double precision, intent(in) :: ampli

  integer :: i,j
    
  do j=1,size(so%field,2)
  do i=1,size(so%field,1)
    so%field(i,j)%P = so%field(i,j)%P + so%coef(i,j)*ampli
  enddo
  enddo

  end subroutine EXPLO_add


end module src_explo
