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
module mesh_gen

! MESH_GEN is a driver for the finite element mesh generation.
! The methods implemented are:
!   - mesh generation by CARTesian topology
!   - interface to database in EMC2's ".ftq" format

  use mesh_cartesian
  use mesh_emc2

  implicit none
  private

  type mesh_type
    private
    integer :: kind
    type (mesh_cart_type), pointer :: cart
    type (mesh_emc2_type), pointer :: emc2
  end type mesh_type
  
  integer, parameter :: tag_cart = 1 &
                       ,tag_emc2 = 2

  public :: mesh_type, MESH_read, MESH_build

contains

!===========================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_DEF
! PURPOSE: Selects a method to import/generate a mesh.
! SYNTAX : &MESH_DEF method /
!
! ARG: method   [name] [none] Choices are: 'CARTESIAN' or 'EMC2'
!               
! END INPUT BLOCK

subroutine MESH_read(mesh,iin)

  use stdio, only : IO_abort

  type(mesh_type), intent(out) :: mesh
  integer, intent(in) :: iin

  character(10) :: method

  NAMELIST / MESH_DEF / method

  method = ' '
  rewind(iin)
  read(iin,MESH_DEF,END= 100)
  if (method == ' ') call IO_abort('mesh_read: you must set the "method" ')

  select case (method) 
    case('CARTESIAN') 
      mesh%kind = tag_cart
      allocate(mesh%cart)
      call CART_read(mesh%cart,iin)
    case('EMC2')
      mesh%kind = tag_emc2
      allocate(mesh%emc2)
      call EMC2_read(mesh%emc2,iin)
    case default
      call IO_abort('mesh_read: unknown "method" ')
  end select

  return

100 call IO_abort('mesh_read: MESH_DEF input block not found')

end subroutine MESH_read

!===========================================================================
! Build the mesh. On exit the following must be well defined :
!   - coorg   = coordinates of the control nodes
!   - knods   = control nodes indices for each element
!   - kmato   = domain index for each element
!   - bc_topo = boundary conditions topology
! Note that for some mesh methods some of these data were 
! already defined in a previous step (see MESH_read)

subroutine MESH_build(grid,mesh)

  use stdio, only: IO_new_unit
  use spec_grid, only : sem_grid_type

  type(sem_grid_type), intent(inout) :: grid
  type(mesh_type), intent(in) :: mesh

  integer ounit,e,i

  select case (mesh%kind)
    case(tag_cart)
      call CART_build(mesh%cart,grid)
    case(tag_emc2)
      call EMC2_build(mesh%emc2,grid)
  end select

  ounit = IO_new_unit()
  open(ounit,file='ElmtNodes_sem2d.tab')
  do e=1,grid%nelem
    write(ounit,*) grid%knods(:,e)
  enddo
  close(ounit)

  open(ounit,file='MeshNodesCoord_sem2d.tab')
  do i=1,grid%npgeo
    write(ounit,*) real(grid%coorg(:,i))
  enddo
  close(ounit)

end subroutine MESH_build

end module mesh_gen
