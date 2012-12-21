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
module mesh_cartesian

! MESH_CARTESIAN: generation of mesh with cartesian topology

  implicit none
  private

  type domain_type
    integer :: tag,ex(2),ez(2)
  end type domain_type

  type mesh_cart_type
    private
    logical :: FaultX
    double precision :: xmin,xmax,zmin,zmax
    integer :: ndom,nz,nx
    type (domain_type), pointer :: domains(:) 
  end type mesh_cart_type

  public :: mesh_cart_type,CART_read,CART_build

  integer, parameter :: ndime = 2, ngnod = 4

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_CART [mesh]
! PURPOSE: Rectangular box with structured mesh.
! SYNTAX : &CART xlim,zlim,nelem /
!
! ARG: xlim     [dble(2)] [none] X limits of the box (min and max)
! ARG: zlim     [dble(2)] [none] Z limits of the box (min and max)
! ARG: nelem    [int(2)] [none]  Number of elements along each direction
! ARG: FaultX   [log] [F] Cut the box in the middle by a horizontal fault
!                         If enabled, nelem(2) must be even
!
! NOTE: the tags for the boundaries are
!               1       Bottom 
!               2       Right        
!               3       Top  
!               4       Left
!               5       FaultBottom
!               6       FaultTop
!
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : MESH_CART_DOMAIN
! PURPOSE: Define a subdomain within a structured meshed box.
! SYNTAX : &MESH_CART_DOMAIN tag,ex,ez /
!
! ARG: tag      [int] [none] Tag number assigned to this domain. 
! ARG: ex       [int(2)] [none]  First and last element along the X direction.
! ARG: ez       [int(2)] [none]  First and last element along the Z direction.
!
! NOTE   : If you ignore this input block a single domain (tag=1) will span 
!          the whole box 
!
! END INPUT BLOCK


subroutine CART_read(mesh,iin)

  use stdio, only : IO_abort

  type(mesh_cart_type), intent(out) :: mesh
  integer, intent(in) :: iin

  double precision :: init_double,xlim(2),zlim(2)
  integer :: nelem(2),nx,nz,tag,ex(2),ez(2),n_domains,i
  logical :: FaultX

  NAMELIST / MESH_CART /  xlim,zlim,nelem,FaultX
  NAMELIST / MESH_CART_DOMAIN / tag,ex,ez

  init_double = huge(init_double)
  xlim = (/ 0.d0,init_double /)
  zlim = xlim
  nelem = 0
  FaultX = .false.
  
  rewind(iin)
  read(iin,MESH_CART,END=100)

  if (xlim(2) == init_double) call IO_abort('CART_read: you must set xlim')
  if (zlim(2) == init_double) call IO_abort('CART_read: you must set zlim')
  nx = nelem(1)
  nz = nelem(2)
  if (nx <= 0) call IO_abort('CART_read: nx must be positive')
  if (nz <= 0) call IO_abort('CART_read: nz must be positive')

  mesh%xmin = xlim(1)
  mesh%xmax = xlim(2)
  mesh%zmin = zlim(1)
  mesh%zmax = zlim(2)
  mesh%nx   = nx
  mesh%nz   = nz
  mesh%FaultX = FaultX

 ! Count the domains
  rewind(iin)
  n_domains = 0
  do 
    read(iin,MESH_CART_DOMAIN,END=110)
    n_domains = n_domains + 1
  enddo 
110 continue

  if (n_domains == 0) then
   ! A single domain spans the whole box
    mesh%ndom = 1
    allocate(mesh%domains(1))
    mesh%domains(1)%tag = 1
    mesh%domains(1)%ex  = (/ 1,nx /)
    mesh%domains(1)%ez  = (/ 1,nz /)

  else

    mesh%ndom = n_domains
    allocate(mesh%domains(n_domains))

    rewind(iin)
    do i=1,n_domains
      
      tag  = 0
      ex = 0
      ez = 0
      
      read(iin,MESH_CART_DOMAIN)

      if (tag<1) call IO_abort('CART_read: tag negative or missing')
      if (ex(1)<1 .or. ex(1)>nx) call IO_abort('CART_read: ex(1) out of bounds or missing')
      if (ex(2)<1 .or. ex(2)>nx) call IO_abort('CART_read: ex(2) out of bounds or missing')
      if (ez(1)<1 .or. ez(1)>nz) call IO_abort('CART_read: ez(1) out of bounds or missing')
      if (ez(2)<1 .or. ez(2)>nz) call IO_abort('CART_read: ez(2) out of bounds or missing')

      mesh%domains(i)%tag = tag
      mesh%domains(i)%ex  = ex
      mesh%domains(i)%ez  = ez

    enddo
  endif

  return
100 call IO_abort('CART_read: input block not found')

end subroutine CART_read

!=====================================================================
! CART_BUILD:
!

subroutine CART_build(mesh,grid)

  use spec_grid, only : sem_grid_type,edge_D,edge_R,edge_U,edge_L, & 
                        side_D,side_R,side_U,side_L
  use stdio, only : IO_abort

  type(mesh_cart_type), intent(in) :: mesh
  type(sem_grid_type), intent(inout) :: grid

  double precision, allocatable :: x(:),z(:)
  integer :: nxp,nzp,i,j,k,ilast,ifirst,idom

  integer, parameter :: fault_D = 5
  integer, parameter :: fault_U = 6

  nxp = mesh%nx+1 
  if (mesh%FaultX) then 
    nzp = mesh%nz+2
  else
    nzp = mesh%nz+1
  endif
  grid%npgeo = nxp*nzp
  grid%nelem = mesh%nx*mesh%nz
  grid%ndime = NDIME
  grid%ngnod = NGNOD

  ! Allocations
  allocate(grid%coorg(grid%ndime,grid%npgeo))
  allocate(grid%knods(grid%ngnod,grid%nelem))
  allocate(grid%tag(grid%nelem))
  
  ! Coordinates of control nodes
  allocate(x(nxp))
  allocate(z(nzp))
  x = mesh%xmin +(mesh%xmax-mesh%xmin)/dble(mesh%nx) *(/ (i, i=0,nxp-1) /)
  if (mesh%FaultX) then 
    z = mesh%zmin +(mesh%zmax-mesh%zmin)/dble(mesh%nz) &
               *(/ (j, j=0,nzp/2-1), (j, j=nzp/2-1,nzp-2) /)
  else
    z = mesh%zmin +(mesh%zmax-mesh%zmin)/dble(mesh%nz) *(/ (j, j=0,nzp-1) /) 
  endif
  ilast  = 0
  do j=1,nzp
    ifirst = ilast + 1
    ilast  = ilast + nxp
    grid%coorg(1, ifirst:ilast ) = x
    grid%coorg(2, ifirst:ilast ) = z(j)
  enddo
  deallocate(x,z)
  
 ! Domain tags, usually for material sets
  grid%tag = 0
  do idom=1,mesh%ndom
    do i=mesh%domains(idom)%ex(1),mesh%domains(idom)%ex(2)
    do j=mesh%domains(idom)%ez(1),mesh%domains(idom)%ez(2)
      grid%tag( CART_ij_to_index(i,j,mesh%nx) ) = mesh%domains(idom)%tag
    enddo
    enddo
  enddo
  if (any(grid%tag == 0)) call IO_abort('CART_build: Domain tags not entirely set')

 ! Control nodes of each element
  k = 0
  do j=1,mesh%nz
  do i=1,mesh%nx
    k = k + 1
    grid%knods(1,k) = CART_ij_to_index(i,j,nxp)
    grid%knods(2,k) = CART_ij_to_index(i+1,j,nxp)
    grid%knods(3,k) = CART_ij_to_index(i+1,j+1,nxp)
    grid%knods(4,k) = CART_ij_to_index(i,j+1,nxp)
  enddo
  enddo 
  if (mesh%FaultX) grid%knods(:,mesh%nx*mesh%nz/2+1:) = grid%knods(:,mesh%nx*mesh%nz/2+1:)+nxp
  
 ! Boundary conditions
  if (mesh%FaultX) then
    allocate(grid%bounds(6))
  else
    allocate(grid%bounds(4))
  endif

 ! Down
  grid%bounds(side_D)%tag = side_D
  grid%bounds(side_D)%nelem = mesh%nx
  allocate(grid%bounds(side_D)%bulk_element(mesh%nx))
  allocate(grid%bounds(side_D)%element_edge(mesh%nx))
  do i =1,mesh%nx
    grid%bounds(side_D)%bulk_element(i) = CART_ij_to_index(i,1,mesh%nx)
  enddo
  grid%bounds(side_D)%element_edge = edge_D

 ! Right
  grid%bounds(side_R)%tag = side_R
  grid%bounds(side_R)%nelem = mesh%nz
  allocate(grid%bounds(side_R)%bulk_element(mesh%nz))
  allocate(grid%bounds(side_R)%element_edge(mesh%nz))
  do j =1,mesh%nz
    grid%bounds(side_R)%bulk_element(j) = CART_ij_to_index(mesh%nx,j,mesh%nx)
  enddo
  grid%bounds(side_R)%element_edge = edge_R

 ! Up
  grid%bounds(side_U)%tag = side_U
  grid%bounds(side_U)%nelem = mesh%nx
  allocate(grid%bounds(side_U)%bulk_element(mesh%nx))
  allocate(grid%bounds(side_U)%element_edge(mesh%nx))
  do i =1,mesh%nx
    grid%bounds(side_U)%bulk_element(i) = CART_ij_to_index(i,mesh%nz,mesh%nx)
  enddo
  grid%bounds(side_U)%element_edge = edge_U

 ! Left
  grid%bounds(side_L)%tag = side_L
  grid%bounds(side_L)%nelem = mesh%nz
  allocate(grid%bounds(side_L)%bulk_element(mesh%nz))
  allocate(grid%bounds(side_L)%element_edge(mesh%nz))
  do j =1,mesh%nz
    grid%bounds(side_L)%bulk_element(j) = CART_ij_to_index(1,j,mesh%nx)
  enddo
  grid%bounds(side_L)%element_edge = edge_L

  if (mesh%FaultX) then

   ! Fault Up: 
    grid%bounds(fault_U)%tag = fault_U
    grid%bounds(fault_U)%nelem = mesh%nx
    allocate(grid%bounds(fault_U)%bulk_element(mesh%nx))
    allocate(grid%bounds(fault_U)%element_edge(mesh%nx))
    do i =1,mesh%nx
      grid%bounds(fault_U)%bulk_element(i) = CART_ij_to_index(i,mesh%nz/2+1,mesh%nx)
    enddo
    grid%bounds(fault_U)%element_edge = edge_D
    
   ! Fault Down
    grid%bounds(fault_D)%tag = fault_D
    grid%bounds(fault_D)%nelem = mesh%nx
    allocate(grid%bounds(fault_D)%bulk_element(mesh%nx))
    allocate(grid%bounds(fault_D)%element_edge(mesh%nx))
    do i =1,mesh%nx
      grid%bounds(fault_D)%bulk_element(i) = CART_ij_to_index(i,mesh%nz/2,mesh%nx)
    enddo
    grid%bounds(fault_D)%element_edge = edge_U

  endif

end subroutine CART_build

!=====================================================================
! In the following conversion function the numbering conventions are 
!   i=1:n 
!   j=1:
!   iglob = 1:

integer function CART_ij_to_index(i,j,n)

  integer, intent(in) :: i,j,n

  CART_ij_to_index = (j-1)*n + i

end function CART_ij_to_index

end module mesh_cartesian
