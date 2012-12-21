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
module bc_periodic

  use spec_grid, only : bc_topo_type

  implicit none
  private

  type bc_periodic_type
    private
    type(bc_topo_type), pointer :: master,slave
  end type bc_periodic_type

  interface BC_PERIO_set
    module procedure BC_PERIO_set_field &
                    ,BC_PERIO_set_scal
  end interface

  public :: bc_periodic_type,BC_PERIO_read &
           ,BC_PERIO_set,BC_PERIO_init

contains

!=======================================================================

  subroutine BC_PERIO_read(bc,iin)

  use echo, only : echo_input,iout

  type (bc_periodic_type) :: bc
  integer, intent(in) :: iin

  if (echo_input) write(iout,200)
  
  return
  200 format(5x,'Type   = Periodic')

  end subroutine BC_PERIO_read


!=======================================================================
!
subroutine BC_PERIO_init(bc,tags,grid,rmass)

  use spec_grid, only : sem_grid_type,BC_inquire
  use stdio, only: IO_abort

  type(bc_periodic_type), intent(out) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(inout) :: rmass(:)
  integer, intent(in) :: tags(2)

  double precision :: shift1,shift2
  double precision, parameter :: tiny=1d-3

  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%master )
  call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%slave  )
  
  if (bc%master%nelem/=bc%slave%nelem) &
   call IO_abort('bc_perio_init: number of boundary elements do not match')

  if (bc%master%npoin/=bc%slave%npoin) &
   call IO_abort('bc_perio_init: number of nodes on boundaries do not match')

  shift1 = grid%coord(1,bc%master%bulk_node(1)) - grid%coord(1,bc%slave%bulk_node(1))
  shift2 = grid%coord(2,bc%master%bulk_node(1)) - grid%coord(2,bc%slave%bulk_node(1))
  if ( any(abs(grid%coord(1,bc%master%bulk_node)-grid%coord(1,bc%slave%bulk_node)-shift1)>tiny) &
   .OR.any(abs(grid%coord(2,bc%master%bulk_node)-grid%coord(2,bc%slave%bulk_node)-shift2)>tiny) )&
   call IO_abort('bc_perio_init: coordinates on boundaries do not match properly')

  call BC_PERIO_set_scal(bc,rmass,assemble='yes')

end subroutine BC_PERIO_init


!=======================================================================
!
  subroutine BC_PERIO_set_field(bc,field,assemble)

    type(bc_periodic_type), intent(in) :: bc
    character(*), intent(in) :: assemble
    double precision, intent(inout) :: field(:,:)

!------ Assemble in master nodes    
    if (assemble == 'yes') &
        field(bc%master%bulk_node,:) = field(bc%master%bulk_node,:) &
                                     + field(bc%slave%bulk_node,:)

!------ Redistribute to slave nodes
    field(bc%slave%bulk_node,:) = field(bc%master%bulk_node,:)

  end subroutine BC_PERIO_set_field


!=======================================================================
!
  subroutine BC_PERIO_set_scal(bc,scal,assemble)

    type(bc_periodic_type), intent(in) :: bc
    character(*), intent(in) :: assemble
    double precision, intent(inout) :: scal(:)
    integer :: i,iglob_master

!------ Assemble in master nodes    
    if (assemble == 'yes') scal(bc%master%bulk_node) = &
      scal(bc%master%bulk_node) + scal(bc%slave%bulk_node)

!------ Redistribute to slave nodes
    scal(bc%slave%bulk_node) = scal(bc%master%bulk_node)

    end subroutine BC_PERIO_set_scal


end module bc_periodic
