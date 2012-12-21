! SEM2DPACK version 2.3.4 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module mesh_layers

! MESH_LAYERS: structured mesh generation for layered media

  use distribution_cd

!-- Renumbering of elements by Reverse Cuthill-McKee algorithm --
!   to improve data locality (optimize cache usage)
  use constants, only : OPT_RENUMBER, NDIME

  implicit none
  private

  type layer_type
    integer :: nez, tag
    double precision, pointer :: ztop(:) => null()
    type (cd_type) :: ztopin
  end type layer_type

  type mesh_layers_type
    private
    double precision :: xmin,xmax,zmin
    integer :: nlayer,ngnod,nz,nx,ezflt
    type (layer_type), pointer :: layer(:)=>null()
  end type mesh_layers_type

  public :: mesh_layers_type, MESH_LAYERS_read, MESH_LAYERS_build

contains

!=====================================================================
!
subroutine MESH_LAYERS_read(mesh,iin)

  use stdio, only : IO_abort, IO_new_unit
  use echo, only : echo_input, iout

  type(mesh_layers_type), intent(out) :: mesh
  integer, intent(in) :: iin

  integer :: i,iin2,ngnod
  character(50) :: filename
  character(2) :: eltype

  ! Read &MESH_LAYERED input block
  rewind(iin)
  call read_layered(iin,mesh,filename)

  ! Read nlayer &MESH_LAYER blocks or data from separate file
  ! Listed in input file from top to bottom
  ! but stored here from bottom to top (index i)
  ! Default tags are reversed: 1=top, nlayer=bottom 
  allocate(mesh%layer(mesh%nlayer))
  mesh%ngnod = 4

  if (filename=='') then
    do i=mesh%nlayer,1,-1
      call read_layer(mesh%layer(i),ngnod,iin, mesh%nlayer+1-i) 
      mesh%ngnod = max(mesh%ngnod, ngnod)
    enddo

  else
    iin2 = IO_new_unit()
    open(iin2,file=filename,status='old')
    do i=mesh%nlayer,1,-1
      call read_layer(mesh%layer(i),ngnod,iin2) 
      mesh%ngnod = max(mesh%ngnod, ngnod)
    enddo
    close(iin2)
  endif

  mesh%nz=sum(mesh%layer%nez) 

  if (echo_input) then
    write(eltype,'("Q",i1)') mesh%ngnod
    write(iout,200) mesh%nz,eltype
  endif

  return

200 format(/5x, &
      'Number of elements along Z. . . . . . . . . . . = ',i5,/5x, &
      'Element type. . . . . . . . . . . . . . . . . . = ',a)

end subroutine MESH_LAYERS_read

!--------------------------------------------------
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_LAYERED
! GROUP  : MESH_DEF
! PURPOSE: Structured mesh for layered medium 
!          with surface and interface topography. 
! SYNTAX : &MESH_LAYERED xlim,zmin,nx,file,nlayer,ezflt /
!
! ARG: xlim     [dble(2)] [none] X limits of the box (min and max)
! ARG: zmin     [dble] [none] bottom Z limit of the box 
! ARG: nx       [int] [none]  Number of elements along X direction
! ARG: file     [string] [''] Only for flat layers,
!                name of ASCII file containing layer parameters, 
!                one line per layer, listed from top to bottom, 
!                3 columns per line:
!                (1) vertical position of top boundary,
!                (2) number of elements along Z direction
!                (3) material tag
! ARG: nlayer   [int] [none]  Number of layers
!                If a file name is not given the layer parameters
!                must be given immediately after the &MESH_LAYERED block
!                by nlayer &MESH_LAYER input blocks,
!                one for each layer, listed from top to bottom.
! ARG: ezflt    [int][0] introduce a fault between the ezflt-th and the
!               (ezflt+1)-th element rows. If ezflt=0 (default), no fault
!               is introduced.
!
! NOTE: the following tags are automatically assigned to the boundaries: 
!               1       Bottom 
!               2       Right        
!               3       Top  
!               4       Left
!               5       Fault, bottom side
!               6       Fault, top side
!
! END INPUT BLOCK

subroutine read_layered(iin,mesh,file)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort, IO_file_length

  type(mesh_layers_type), intent(out) :: mesh
  integer, intent(in) :: iin
  character(50), intent(out) :: file

  double precision :: init_double,xlim(2),zmin
  integer :: nx,nlayer,ezflt

  NAMELIST / MESH_LAYERED /  xlim,zmin,nx,nlayer,file,ezflt

  init_double = huge(init_double)

  xlim = (/ 0.d0,init_double /)
  zmin = init_double
  nx = 0
  nlayer = 0
  file= ''
  ezflt = 0
  
  read(iin,MESH_LAYERED,END=100)

  if (xlim(2) == init_double) call IO_abort('MESH_LAYERS_read: you must set xlim')
  if (zmin == init_double) call IO_abort('MESH_LAYERS_read: you must set zmin')
  if (nx <= 0) call IO_abort('MESH_LAYERS_read: nx must be positive')
  if (nlayer <= 0 .and. file=='') &
    call IO_abort('MESH_LAYERS_read: nlayer or file must be set')
  if (ezflt<0) call IO_abort('MESH_LAYERS_read: ezflt must be positive')

  if (file/='') nlayer = IO_file_length(file)

  if (echo_input) then
    write(iout,200) xlim, zmin, nx, nlayer
    if (file/='') write(iout,210) trim(file)
    if (ezflt>0) write(iout,220) ezflt
  endif

  mesh%xmin = xlim(1)
  mesh%xmax = xlim(2)
  mesh%zmin = zmin
  mesh%nx   = nx
  mesh%nlayer = nlayer
  mesh%ezflt = ezflt

  return

100 call IO_abort('MESH_LAYERS_read: input block not found')

200 format(5x, &
      'Minimum X . . . . . . . . . . . . . . (xlim(1)) = ',EN12.3,/5x, &
      'Maximum X . . . . . . . . . . . . . . (xlim(2)) = ',EN12.3,/5x, &
      'Minimum Z . . . . . . . . . . . . . . . .(zmin) = ',EN12.3,/5x, &
      'Number of elements along X. . . . . . . . .(nx) = ',i5,/5x, &
      'Number of layers  . . . . . . . . . . .(nlayer) = ',i5)

210 format(5x, &
      'Layers read from file . . . . . . . . . .(file) = ',a)

220 format(5x, &
      'Fault on top of this element row  . . . (ezflt) = ',i0)

end subroutine read_layered

!--------------------------------------------------
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_LAYER
! GROUP  : MESH_DEF
! PURPOSE: Define mesh parameters for one layer
! SYNTAX : &MESH_LAYER nz, ztop|ztopH, tag /
!          followed by a DIST_XXXX block if ztopH is set
!
! ARG: nz       [int] [none]  Number of elements in layer along Z direction
! ARG: ztop     [dble] [none] Only for layers with flat top surface: 
!                vertical position of top boundary
! ARG: ztopH    [string] ['none'] Only for layers with irregular top boundary: 
!                name of distribution, 'LINEAR', 'SPLINE' or any other
!                1D distribution available through a DIST_XXXX block.
! ARG: tag      [int] [none]  Material tag
!                 If not given, a tag is automatically assigned to the layer, 
!                 sequentially numbered from top to bottom (top layer tag =1)
!
! END INPUT BLOCK

subroutine read_layer(layer,ngnod,iin,tagdef)
  
  use stdio, only : IO_abort
  use echo, only : echo_input, iout

  type(layer_type), intent(inout) :: layer
  integer, intent(out) :: ngnod
  integer, intent(in) :: iin
  integer, intent(in), optional :: tagdef

  double precision :: init_double,ztop
  integer :: nz,tag
  character(20) :: ztopH

  NAMELIST / MESH_LAYER / nz,ztop,ztopH,tag

  ngnod = 4
  init_double = huge(init_double)

  if (present(tagdef)) then ! read MESH_LAYER input blocks from Par.inp
    nz = 0
    ztop = init_double
    ztopH = ''
    tag  = tagdef 
    read(iin,MESH_LAYER)

    if (ztopH=='' .and. ztop==init_double) &
      call IO_abort('MESH_LAYERS_read: ztop or ztopH must be set')
    if (ztopH /= 'LINEAR' .and. ztopH/='') ngnod = 9

  else
    ztopH=''
    read(iin,*) ztop,nz,tag

  endif

  if (nz<=0) call IO_abort('MESH_LAYERS_read: nz null, negative or missing')
  layer%nez = nz

  call DIST_CD_Read(layer%ztopin,ztop,ztopH,iin,ztopH)

  if (tag<1) call IO_abort('MESH_LAYERS_read: tag null, negative or missing')
  layer%tag = tag

  if (echo_input) write(iout,300) tag,nz,ztopH

  return
  
300 format(/5x, &
      'Layer number. . . . . . . . . . . . . . . (tag) = ',i5,/5x, &
      'Number of elements along Z. . . . . . . . .(nz) = ',i5,/5x, &
      'Z of top surface. . . . . . . . (ztop or ztopH) = ',a)

end subroutine read_layer

!=====================================================================
! CART_BUILD:
!

subroutine MESH_LAYERS_build(mesh,grid)

  use fem_grid, only : fem_grid_type
  use memory_info
  use stdio, only : IO_abort
  use rcmlib, only : genrcm, perm_inverse
  use mesh_structured
  use utils, only : sub2ind

  type(mesh_layers_type), intent(inout) :: mesh
  type(fem_grid_type), intent(inout) :: grid

  double precision, pointer :: zbot(:),ztop(:),coord(:,:)
  integer :: nxp,nzp,i,j,k,ilast,ifirst,nj,jfirst,jlast,tag,idoubling
  integer :: jposlast,jposflt

  if (mesh%ngnod==4) then
    idoubling = 1
  else
    idoubling = 2
  endif
  nxp = idoubling*mesh%nx+1
  nzp = idoubling*mesh%nz+1
  if (mesh%ezflt>0) nzp=nzp+1
  grid%npoin = nxp*nzp
  grid%ngnod = mesh%ngnod
  grid%nelem = mesh%nx*mesh%nz
  grid%flat  = .false.

  ! Allocations
  allocate(grid%coord(NDIME,grid%npoin))
  allocate(grid%knods(grid%ngnod,grid%nelem))
  allocate(grid%tag(grid%nelem))
  call storearray('coorg',size(grid%coord),idouble)
  call storearray('knods',size(grid%knods),iinteg)
  call storearray('tag',size(grid%tag),iinteg)
  
 !-- Coordinates of control nodes

 ! generate surfaces
  allocate(coord(2,nxp))
  coord(1,:) = mesh%xmin +(mesh%xmax-mesh%xmin)/dble(nxp-1) *(/ (i, i=0,nxp-1) /)
  coord(2,:) = 0d0 ! not used
  do k=1,mesh%nlayer
    call DIST_CD_Init(mesh%layer(k)%ztopin,coord, mesh%layer(k)%ztop) 
    if (any(mesh%layer(k)%ztop<mesh%zmin)) &
      call IO_abort('MESH_LAYERS_build: zmin in MESH_LAYERED should be deeper than all layer interfaces')
  enddo

 ! bottom line of nodes
  ifirst = 1
  ilast = nxp
  grid%coord(1,ifirst:ilast) = coord(1,:)
  grid%coord(2,ifirst:ilast) = mesh%zmin
  ztop => grid%coord(2,ifirst:ilast)
  jposlast = 1

 ! for each layer: lines of nodes from bottom (not included) to top (included)
  do k=1,mesh%nlayer
    zbot => ztop 
    ztop => mesh%layer(k)%ztop
    nj = idoubling*mesh%layer(k)%nez
    jposflt = idoubling*mesh%ezflt+1 -jposlast
    do j=1,nj
      ifirst = ilast + 1
      ilast  = ilast + nxp
      grid%coord(1, ifirst:ilast ) = coord(1,:)
      grid%coord(2, ifirst:ilast ) = zbot+ dble(j)/dble(nj)*(ztop-zbot)
      if (mesh%ezflt>0 .and. j==jposflt) then
        ifirst = ilast + 1
        ilast  = ilast + nxp
        grid%coord(1, ifirst:ilast ) = coord(1,:)
        grid%coord(2, ifirst:ilast ) = zbot+ dble(j)/dble(nj)*(ztop-zbot)
      endif
    enddo
    jposlast = jposlast + nj
  enddo

  deallocate(coord)
  
 ! Domain tags, usually for material sets
  grid%tag = 0
  jlast=0
  do k=1,mesh%nlayer
    jfirst=jlast+1
    jlast=jfirst+mesh%layer(k)%nez -1
    tag = mesh%layer(k)%tag
    do j=jfirst,jlast
    do i=1,mesh%nx
      grid%tag( sub2ind(i,j,mesh%nx) ) = tag
    enddo
    enddo
  enddo
  if (any(grid%tag == 0)) call IO_abort('MESH_LAYERS_build: Domain tags not entirely set')

 ! Control nodes of each element
 ! Elements are sequentially numbered horizontally from bottom-left to top-right
  call MESH_STRUCTURED_connectivity(grid%knods,mesh%nx,mesh%nz,mesh%ngnod,mesh%ezflt)
  
 ! Boundary conditions
  if (mesh%ezflt>0) then
    allocate(grid%bnds(6))
  else
    allocate(grid%bnds(4))
  endif
  call MESH_STRUCTURED_boundaries(grid%bnds,mesh%nx,mesh%nz,mesh%ezflt)

 ! Renumber elements
  if (OPT_RENUMBER) call MESH_STRUCTURED_renumber(grid,mesh%nx,mesh%nz)
  
end subroutine  MESH_LAYERS_build

end module mesh_layers
