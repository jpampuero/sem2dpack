module mesh_gen

! MESH_GEN is a driver for the finite element mesh generation.
! The methods implemented are:
!   - mesh generation by CARTesian topology
!   - interface to database in EMC2's ".ftq" format

  use mesh_cartesian
  use mesh_layers
  use mesh_emc2
  use mesh_mesh2d
  use fem_grid, only : fem_grid_type

  implicit none
  private

  type mesh_type
    private
    integer :: kind
    type (mesh_cart_type), pointer :: cart
    type (mesh_layers_type), pointer :: layers
    type (fem_grid_type), pointer :: emc2
    type (fem_grid_type), pointer :: mesh2d
  end type mesh_type
  
  integer, parameter :: tag_cart = 1 &
                       ,tag_layers = 2 &
                       ,tag_emc2 = 3 &
                       ,tag_mesh2d = 4

  public :: mesh_type, MESH_read, MESH_build, MESH_h

contains

!===========================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_DEF
! PURPOSE: Selects a method to import/generate a mesh.
! SYNTAX : &MESH_DEF method /
!          followed by a &MESH_method input block
!
! ARG: method   [name] [none] Meshing method name:
!               'CARTESIAN', 'LAYERED', 'EMC2', 'MESH2D'
!               
! END INPUT BLOCK

subroutine MESH_read(mesh,iin)

  use stdio, only : IO_abort
  use echo, only : iout, echo_input

  type(mesh_type), intent(out) :: mesh
  integer, intent(in) :: iin

  character(10) :: method

  NAMELIST / MESH_DEF / method

  method = ' '
  rewind(iin)
  read(iin,MESH_DEF,END= 100)
  if (method == ' ') call IO_abort('mesh_read: you must set the "method" ')
  if (echo_input) write(iout,200) method

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
    case('MESH2D')
      mesh%kind = tag_mesh2d
      allocate(mesh%mesh2d)
      call MESH2D_read(mesh%mesh2d,iin)
    case default
      call IO_abort('mesh_read: unknown method')
  end select

  return

100 call IO_abort('mesh_read: MESH_DEF input block not found')

200 format(//,' M e s h   G e n e r a t i o n', &
      /1x,29('='),//5x, &
      'Method  . . . . . . . . . . . . . . . .(method) = ',a)

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
    case(tag_mesh2d)
      call MESH2D_build(mesh%mesh2d,grid)
  end select

  if (echo_init) then 
    write(iout,fmtok)
    write(iout,fmt1,advance='no') 'Saving node coordinates in file MeshNodesCoord_sem2d.tab'
  endif

  ounit = IO_new_unit()
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
!  write(fmt,'(a,i1,a)') "(", grid%ngnod, "(1x,i6))" 
  write(fmt,'(a,i1,a)') "(", grid%ngnod, "(1x,i8))" 
  do e=1,grid%nelem
    write(ounit,fmt) grid%knods(:,e)
  enddo
  close(ounit)

  if (echo_init) write(iout,fmtok)


end subroutine MESH_build

subroutine MESH_h(mesh,h)
  use stdio, only : IO_abort

  type(mesh_type), intent(in) :: mesh
  double precision, intent(out) :: h

  select case (mesh%kind)
    case(tag_cart)
      call CART_h(mesh%cart,h)
    case(tag_layers)
      call IO_abort('MESH_h: Layers not implemented')
    case(tag_emc2)
      call IO_abort('MESH_h: EMC2 not implemented')
    case(tag_mesh2d)
      call IO_abort('MESH_h: Mesh2D not implemented')
  end select

end subroutine

end module mesh_gen
