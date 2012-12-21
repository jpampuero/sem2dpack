! SEM2DPACK version 2.3.5 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! http://web.gps.caltech.edu/~ampuero/
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
module mesh_mesh2d

! MESH_MESH2D: an interface for finite element mesh database in .mesh2d format

  use constants, only : NDIME
  use fem_grid, only : fem_grid_type

  implicit none
  private

  public :: MESH2D_read, MESH2D_build

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_MESH2D
! GROUP  : MESH_DEF
! PURPOSE: Imports a mesh in mesh2d format 
!	   as defined by the PRE/mesh2d mesh generator tools for Matlab
! SYNTAX : &MESH_MESH2D file /
!
! ARG: file     [name] [none] Name of the MESH2D file, including suffix.
!                 The format of this file is:
!
!  "NEL NPEL NNOD NBC"
!  1 line with 4 integers: 
!  nb of elements, nodes per element, total nb of nodes, nb of boundaries
!  "NID X Y"
!  NNOD lines, one per node, with 1 integer and 2 reals: 
!  node id, x, y
!  "EID NODES TAG"
!  NEL lines, one per element, with NPEL+2 integers: 
!  element id, NPEL node ids, tag. 
!  "BCTAG NBEL"                                         |
!  2 integers: boundary tag, nb of boundary elements    |
!  "BID EID EDGE"                                       | repeat for each of
!  NBEL lines, one per boundary element, 3 integers:    | the NBC boundaries
!  boundary element id, bulk element id, edge id        |
!
! END INPUT BLOCK

subroutine MESH2D_read(mesh,iin)

  use memory_info
  use stdio, only : IO_new_unit,IO_abort
  use echo, only : echo_input, iout

  type(fem_grid_type), intent(out) :: mesh
  integer, intent(in) :: iin

  integer :: imesh2d,i,e,nb,idum
  character(50) :: file

  NAMELIST / MESH_MESH2D / file

  file = ' '
  read(iin,MESH_MESH2D,END=100)
  if (file == ' ') call IO_abort('MESH2D_read: you must provide a mesh2d file name')

  imesh2d = IO_new_unit()
  open(imesh2d,file=file,status='old')

  read(imesh2d,*)
  read(imesh2d,*) mesh%nelem, mesh%ngnod, mesh%npoin, nb
  if (echo_input) write(iout,200) file, mesh%npoin, mesh%nelem, mesh%ngnod, nb

  read(imesh2d,*)
  allocate( mesh%coord(NDIME,mesh%npoin) )
  call storearray('coorg',size(mesh%coord),idouble)
  do i=1,mesh%npoin
    read(imesh2d,*) idum,mesh%coord(:,i)
  enddo
 
  read(imesh2d,*)
  allocate (mesh%knods(mesh%ngnod,mesh%nelem) )
  allocate (mesh%tag(mesh%nelem))
  call storearray('knods',size(mesh%knods),iinteg)
  call storearray('tag',size(mesh%tag),iinteg)
  do e=1,mesh%nelem
    read(imesh2d,*) idum, mesh%knods(:,e),mesh%tag(e)
  enddo

  allocate(mesh%bnds(nb))
  do i=1,nb
    read(imesh2d,*)
    read(imesh2d,*) mesh%bnds(i)%tag, mesh%bnds(i)%nelem
    read(imesh2d,*)
    allocate( mesh%bnds(i)%elem( mesh%bnds(i)%nelem ))
    allocate( mesh%bnds(i)%edge( mesh%bnds(i)%nelem ))
    do e=1,mesh%bnds(i)%nelem
      read(imesh2d,*) idum, mesh%bnds(i)%elem(e), mesh%bnds(i)%edge(e)
    enddo
  enddo

  close(imesh2d)

  return

  100 call IO_abort('MESH2D_read: MESH_MESH2D input block not found')

200 format(5x, &
      'File name . . . . . . . . . . . . . . . .(file) = ',a,/5x, &
      'Number of nodes . . . . . . . . . . . . . . . . = ',i5,/5x, &
      'Number of elements. . . . . . . . . . . . . . . = ',i5,/5x, &
      'Number of nodes per element . . . . . . . . . . = ',i5,/5x, &
      'Number of boundaries  . . . . . . . . . . . . . = ',i5 )

end subroutine MESH2D_read


!=====================================================================
! MESH2D_BUILD: Define the model boundaries

subroutine MESH2D_build(mesh,grid)

  use fem_grid, only : fem_grid_type

  type(fem_grid_type), intent(in) :: mesh
  type(fem_grid_type), intent(inout) :: grid

  grid%npoin =  mesh%npoin
  grid%ngnod =  mesh%ngnod
  grid%nelem =  mesh%nelem
  grid%coord => mesh%coord
  grid%knods => mesh%knods
  grid%tag   => mesh%tag
  grid%bnds  => mesh%bnds

end subroutine MESH2D_build

end module mesh_mesh2d
