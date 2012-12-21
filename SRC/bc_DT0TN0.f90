! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bc_DT0TN0
! Boundary condition =
!   Null tangential displacement
!   Null normal stress
! Implemented for vertical flat boundary
  
  use bnd_grid, only : bnd_grid_type

  implicit none
  private

  type bc_DT0TN0_type
    private
    type(bnd_grid_type), pointer :: topo
  end type

  public :: BC_DT0TN0_type,BC_DT0TN0_read,BC_DT0TN0_init,BC_DT0TN0_set

contains

!=======================================================================
!
subroutine bc_DT0TN0_read(bc,iin)
  
  use echo , only: echo_input,iout

  type(bc_DT0TN0_type), pointer :: bc
  integer, intent(in) :: iin

  allocate(bc)
  if (echo_input) write(iout,200)

  return
  200 format(5x,'Type   = Null vertical displacement and horizontal stress')

end subroutine bc_DT0TN0_read


!=======================================================================
!
subroutine bc_DT0TN0_init(bc,tag,grid)

  use echo, only : echo_init,iout
  use spec_grid, only : sem_grid_type,BC_inquire

  type(bc_DT0TN0_type) , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: tag
  
  !-- bc%topo => grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  if (echo_init) write(iout,*) 'DT0TN0 boundary nodes = ',bc%topo%npoin

end subroutine bc_DT0TN0_init


!=======================================================================
!
subroutine bc_DT0TN0_set(bc,field)

  type(bc_DT0TN0_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)

  field(bc%topo%node,2) = 0.d0

end subroutine bc_DT0TN0_set

end module bc_DT0TN0
