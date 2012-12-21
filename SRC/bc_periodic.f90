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
module bc_periodic

  use bnd_grid, only : bnd_grid_type

  implicit none
  private

  type bc_periodic_type
    private
    type(bnd_grid_type), pointer :: master,slave
  end type bc_periodic_type

  interface BC_PERIO_set
    module procedure BC_PERIO_set_field &
                    ,BC_PERIO_set_scal
  end interface

  public :: bc_periodic_type,BC_PERIO_read &
           ,BC_PERIO_set,BC_PERIO_init,BC_PERIO_intersects

contains

!=======================================================================

  subroutine BC_PERIO_read() !bc,iin)

  use echo, only : echo_input,iout

!  type (bc_periodic_type) :: bc
!  integer, intent(in) :: iin

  if (echo_input) write(iout,200)
  
  return
  200 format(5x,'Type   = Periodic')

  end subroutine BC_PERIO_read


!=======================================================================
!
subroutine BC_PERIO_init(bc,tags,grid,rmass)

  use spec_grid, only : sem_grid_type,BC_inquire
  use stdio, only: IO_abort
  use constants, only : TINY_XABS

  type(bc_periodic_type), intent(out) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(inout) :: rmass(:,:)
  integer, intent(in) :: tags(2)

  double precision :: shift1,shift2

  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%master )
  call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%slave  )
  
  if (bc%master%nelem/=bc%slave%nelem) &
   call IO_abort('bc_perio_init: number of boundary elements do not match')

  if (bc%master%npoin/=bc%slave%npoin) &
   call IO_abort('bc_perio_init: number of nodes on boundaries do not match')

  shift1 = grid%coord(1,bc%master%node(1)) - grid%coord(1,bc%slave%node(1))
  shift2 = grid%coord(2,bc%master%node(1)) - grid%coord(2,bc%slave%node(1))
  if ( any(abs(grid%coord(1,bc%master%node)-grid%coord(1,bc%slave%node)-shift1)>TINY_XABS) &
   .OR.any(abs(grid%coord(2,bc%master%node)-grid%coord(2,bc%slave%node)-shift2)>TINY_XABS) )&
   call IO_abort('bc_perio_init: coordinates on boundaries do not match properly')

  call BC_PERIO_set_field(bc,rmass)

end subroutine BC_PERIO_init


!=======================================================================
!
  subroutine BC_PERIO_set_field(bc,field)

    type(bc_periodic_type), intent(in) :: bc
    double precision, intent(inout) :: field(:,:)

!------ Assemble in master nodes    
    field(bc%master%node,:) = field(bc%master%node,:) + field(bc%slave%node,:)

!------ Redistribute to slave nodes
    field(bc%slave%node,:) = field(bc%master%node,:)

  end subroutine BC_PERIO_set_field


!=======================================================================
!
  subroutine BC_PERIO_set_scal(bc,scal)

    type(bc_periodic_type), intent(in) :: bc
    double precision, intent(inout) :: scal(:)

!------ Assemble in master nodes    
    scal(bc%master%node) = scal(bc%master%node) + scal(bc%slave%node)

!------ Redistribute to slave nodes
    scal(bc%slave%node) = scal(bc%master%node)

    end subroutine BC_PERIO_set_scal


!=======================================================================
!
  function BC_PERIO_intersects(bnd,perio) result(inter)

    type(bc_periodic_type), pointer :: perio
    type(bnd_grid_type), intent(in) :: bnd
    logical :: inter

    if (associated(perio)) then
      inter = any( perio%master%node==bnd%node(1) .or. perio%slave%node==bnd%node(1) ) &
        .and. any( perio%master%node==bnd%node(bnd%npoin) .or. perio%slave%node==bnd%node(bnd%npoin) )
    else
      inter = .false.
    endif

  end function BC_PERIO_intersects

end module bc_periodic
