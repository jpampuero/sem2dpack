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
module mesh_gen

! MESH_GEN is a driver for the finite element mesh generation.
! The methods implemented are:
!   - mesh generation by CARTesian topology
!   - interface to database in EMC2's ".ftq" format

  use mesh_cartesian
  use mesh_layers
  use mesh_emc2
  use fem_grid, only : fem_grid_type

  implicit none
  private

  type mesh_type
    private
    integer :: kind
    type (mesh_cart_type), pointer :: cart
    type (mesh_layers_type), pointer :: layers
    type (fem_grid_type), pointer :: emc2
  end type mesh_type
  
  integer, parameter :: tag_cart = 1 &
                       ,tag_layers = 2 &
                       ,tag_emc2 = 3

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
! ARG: method   [name] [none] 'CARTESIAN', 'LAYERED' or 'EMC2'
!               The &MESH_DEF input block must be followed by a
!               &MESH_method input block
!               
! END INPUT BLOCK

subroutine MESH_read(mesh,iin)

  use stdio, only : IO_abort
  use echo, only : iout

  type(mesh_type), intent(out) :: mesh
  integer, intent(in) :: iin

  character(10) :: method

  NAMELIST / MESH_DEF / method

  method = ' '
  rewind(iin)
  read(iin,MESH_DEF,END= 100)
  if (method == ' ') call IO_abort('mesh_read: you must set the "method" ')
  write(iout,200) method

  select case (method) 
    case('CARTESIAN') 
      mesh%kind = tag_cart
      allocate(mesh%cart)
      call CART_read(mesh%cart,iin)
    case('LAYERED') 
      mesh%kind = tag_layers
      allocate(mesh%layers)
      call MESH_LAYERS_read(mesh%layers,iin)
    case('EMC2')
      mesh%kind = tag_emc2
      allocate(mesh%emc2)
      call EMC2_read(mesh%emc2,iin)
    case default
      call IO_abort('mesh_read: unknown "method" ')
  end select

  return

100 call IO_abort('mesh_read: MESH_DEF input block not found')

200 format(//,' M e s h   G e n e r a t i o n', &
      /1x,29('='),//5x, &
      'Method. . . . . . . . . . . . . . . . .(method) = ',a)

end subroutine MESH_read

!===========================================================================
! Build the mesh. On exit the following must be well defined :
!   - coord = coordinates of the control nodes
!   - knods = control nodes indices for each element
!   - tag   = domain index for each element
!   - bnds  = boundary conditions topology
! Note that for some mesh methods some of these data were 
! already defined in a previous step (see MESH_read)

subroutine MESH_build(grid,mesh)

  use stdio, only: IO_new_unit
  use fem_grid, only : fem_grid_type
  use echo, only : iout,echo_init, fmt1,fmtok

  type(fem_grid_type), intent(inout) :: grid
  type(mesh_type), intent(inout) :: mesh

  integer ounit,e,i
  character(10) :: fmt

  if (echo_init) then
    write(iout,*)
    write(iout,fmt1,advance='no') 'Defining the FEM mesh'
  endif

  select case (mesh%kind)
    case(tag_cart)
      call CART_build(mesh%cart,grid)
    case(tag_layers)
      call MESH_LAYERS_build(mesh%layers,grid)
    case(tag_emc2)
      call EMC2_build(mesh%emc2,grid)
  end select

  if (echo_init) then 
    write(iout,fmtok)
    write(iout,fmt1,advance='no') 'Saving node coordinates in file MeshNodesCoord_sem2d.tab'
  endif

  open(ounit,file='MeshNodesCoord_sem2d.tab')
  do i=1,grid%npoin
    write(ounit,*) real(grid%coord(:,i))
  enddo
  close(ounit)

  if (echo_init) then 
    write(iout,fmtok)
    write(iout,fmt1,advance='no') 'Saving element connectivity in file ElmtNodes_sem2d.tab'
  endif

  ounit = IO_new_unit()
  open(ounit,file='ElmtNodes_sem2d.tab')
  write(fmt,'(a,i1,a)') "(", grid%ngnod, "(x,i6))" 
  do e=1,grid%nelem
    write(ounit,fmt) grid%knods(:,e)
  enddo
  close(ounit)

  if (echo_init) write(iout,fmtok)

end subroutine MESH_build

end module mesh_gen
