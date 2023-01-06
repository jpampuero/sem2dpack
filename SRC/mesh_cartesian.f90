module mesh_cartesian

! MESH_CARTESIAN: generation of a rectangular mesh

!-- Renumbering of elements by Reverse Cuthill-McKee algorithm --
!   to improve data locality (optimize cache usage)
  use constants, only : OPT_RENUMBER, NDIME

  implicit none
  private

  type domain_type
    integer :: tag,ex(2),ez(2)
  end type domain_type

  type mesh_cart_type
    private
    double precision :: xmin,xmax,zmin,zmax,splitD
    integer :: ndom,nz,nx,ezflt,fztag,fznz
    logical :: split
    type (domain_type), pointer :: domains(:) 
  end type mesh_cart_type

  public :: mesh_cart_type,CART_read,CART_build,CART_h

  integer, parameter :: NGNOD = 4

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_CART
! GROUP  : MESH_DEF
! PURPOSE: Rectangular box with structured mesh.
! SYNTAX : &MESH_CART xlim, zlim, nelem, ezflt,fztag, FaultX /
!
! ARG: xlim     [dble(2)] [none] X limits of the box (min and max)
! ARG: zlim     [dble(2)] [none] Z limits of the box (min and max)
! ARG: nelem    [int(2)] [none]  Number of elements along each direction
! ARG: ezflt    [int][0] introduce a horizontal fault between the ezflt-th
!                and the (ezflt+1)-th element rows. Rows are numbered from
!                bottom to top, starting at ezflt=1.
!                If ezflt=0, (default) no fault is introduced inside the box
!                (for symmetric problems a fault can still be set at an external boundary)
!                If ezflt=-1, a fault is introduced at/near the middle of the box
!                (ezflt is reset to int[nelem(2)/2])
! ARG: fztag    [int][0] fault zone tag for elements close to the fault
!                Useful to set a damping layer near the fault.
!                If ezflt=0, a fault is assumed at the bottom boundary
! ARG: split    [log][F] splits the bottom boundary into two segments.
! ARG: splitD   [dble(2)][none] X distance that splits bottom boundary
! ARG: fznz     [int][1] vertical size (number of elements) of near-fault layer
! ARG: FaultX   [log][F] Same as ezflt=-1. Obsolete (will be deprecated) 
!
! NOTE: the following tags are automatically assigned to the boundaries: 
!               1       Bottom or Bottom Right
!               2       Right        
!               3       Top  
!               4       Left
!               5       Fault, bottom side or Bottom Left
!               6       Fault, top side
!
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : MESH_CART_DOMAIN
! PURPOSE: Define a subdomain within a structured meshed box.
! SYNTAX : &MESH_CART_DOMAIN tag,ex,ez /
!
! ARG: tag      [int] [none] Tag number assigned to this domain. 
! ARG: ex       [int(2)] [none]	Horizontal index of the first and last elements.
!               The leftmost element column has horizontal index 1.
! ARG: ez       [int(2)] [none]	Vertical index of the first and last elements.
!               The bottom element row has vertical index 1.
!
! NOTE   : If you ignore this input block a single domain (tag=1) will span 
!          the whole box 
!
! END INPUT BLOCK


subroutine CART_read(mesh,iin)

  use stdio, only : IO_abort
  use echo, only : echo_input, iout

  type(mesh_cart_type), intent(out) :: mesh
  integer, intent(in) :: iin

  double precision :: init_double,xlim(2),zlim(2), splitD
  integer :: nelem(2),ezflt,nx,nz,tag,ex(2),ez(2),n_domains,i,fztag,fznz
  logical :: FaultX,split

  NAMELIST / MESH_CART /  xlim,zlim,nelem,FaultX,ezflt,split,fztag,fznz,splitD
  NAMELIST / MESH_CART_DOMAIN / tag,ex,ez

  init_double = huge(init_double)
  xlim = (/ 0.d0,init_double /)
  zlim = xlim
  nelem = 0
  FaultX = .false.
  ezflt = 0
  fztag = 0
  fznz = 1
  splitD = init_double
  split = .false.
  
  rewind(iin)
  read(iin,MESH_CART,END=100)

  if (xlim(2) == init_double) call IO_abort('CART_read: you must set xlim')
  if (zlim(2) == init_double) call IO_abort('CART_read: you must set zlim')
  nx = nelem(1)
  nz = nelem(2)
  if (nx <= 0) call IO_abort('CART_read: nelem(1) must be positive')
  if (nz <= 0) call IO_abort('CART_read: nelem(2) must be positive')
  if (FaultX) ezflt=-1
  if (ezflt >= nz) call IO_abort('CART_read: ezflt must be < nelem(2)')
  if (ezflt== -1) ezflt = nz/2
  if (ezflt <-1) call IO_abort('CART_read: ezflt must be >= -1')
  if (fztag<0) call IO_abort('MESH_LAYERS_read: fztag must be positive')
  if (fznz<1) call IO_abort('MESH_LAYERS_read: fznz must be strictly positive')
  if (split .and. ( splitD < xlim(1) .or. splitD > xlim(2) ) ) &
    call IO_abort('CART_read: splitD must be in xlim range')

  if (echo_input) then
    write(iout,200) xlim, zlim, nelem
    if (ezflt>0) write(iout,210) ezflt
    if (fztag>0) write(iout,230) fztag,fznz
    if (split) write(iout,250) splitD
  endif

  mesh%xmin = xlim(1)
  mesh%xmax = xlim(2)
  mesh%zmin = zlim(1)
  mesh%zmax = zlim(2)
  mesh%nx   = nx
  mesh%nz   = nz
  mesh%ezflt = ezflt
  mesh%fztag = fztag
  mesh%fznz = fznz
  mesh%split = split 
  mesh%splitD = splitD
 
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

      if (tag<1) call IO_abort('CART_read: tag null, negative or missing')
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

200 format(5x, &
      'Minimum X . . . . . . . . . . . . . . (xlim(1)) = ',1pe10.3,/5x, &
      'Maximum X . . . . . . . . . . . . . . (xlim(2)) = ',1pe10.3,/5x, &
      'Minimum Z . . . . . . . . . . . . . . (zlim(1)) = ',1pe10.3,/5x, &
      'Maximum Z . . . . . . . . . . . . . . (zlim(2)) = ',1pe10.3,/5x, &
      'Number of elements along X. . . . . .(nelem(1)) = ',I0,/5x, &
      'Number of elements along Z. . . . . .(nelem(2)) = ',I0)

210 format(5x, &
      'Fault on top of this element row  . . . (ezflt) = ',I0)

230 format(5x, &
      'Tag for elements in fault zone  . . . . (fztag) = ',I0,/5x, &
      'Vertical nb of elements in fault zone . .(fznz) = ',I0)
250 format(5x, &
      'Splitting the fault at . . . . . . . . (splitD) = ',1pe10.3,/5x)

end subroutine CART_read

!=====================================================================
! CART_BUILD:
!

subroutine CART_build(mesh,grid)

  use fem_grid, only : fem_grid_type
  use mesh_structured
  use memory_info
  use stdio, only : IO_abort
  use utils, only : sub2ind

  type(mesh_cart_type), intent(in) :: mesh
  type(fem_grid_type), intent(inout) :: grid

  double precision, allocatable :: x(:),z(:)
  integer :: nxp,nzp,i,j,j1,j2,ilast,ifirst,idom,splitN

  nxp = mesh%nx+1 
  if (mesh%ezflt>0) then 
    nzp = mesh%nz+2
  else
    nzp = mesh%nz+1
  endif
  grid%npoin = nxp*nzp
  grid%ngnod = NGNOD
  grid%nelem = mesh%nx*mesh%nz
  grid%flat  = .true.

  ! Allocations
  allocate(grid%coord(NDIME,grid%npoin))
  allocate(grid%knods(grid%ngnod,grid%nelem))
  allocate(grid%tag(grid%nelem))
  call storearray('coorg',size(grid%coord),idouble)
  call storearray('knods',size(grid%knods),iinteg)
  call storearray('tag',size(grid%tag),iinteg)
  
  ! Coordinates of control nodes
  allocate(x(nxp))
  allocate(z(nzp))
  x = mesh%xmin +(mesh%xmax-mesh%xmin)/dble(mesh%nx) *(/ (i, i=0,nxp-1) /)
  if (mesh%ezflt>0) then 
    z = mesh%zmin +(mesh%zmax-mesh%zmin)/dble(mesh%nz) &
               *(/ (j, j=0,mesh%ezflt), (j, j=mesh%ezflt,mesh%nz) /)
  else
    z = mesh%zmin +(mesh%zmax-mesh%zmin)/dble(mesh%nz) *(/ (j, j=0,mesh%nz) /) 
  endif
  ilast  = 0
  do j=1,nzp
    ifirst = ilast + 1
    ilast  = ilast + nxp
    grid%coord(1, ifirst:ilast ) = x
    grid%coord(2, ifirst:ilast ) = z(j)
  enddo
  deallocate(x,z)
  
 ! Domain tags, usually for material sets
  grid%tag = 0
  do idom=1,mesh%ndom
    do i=mesh%domains(idom)%ex(1),mesh%domains(idom)%ex(2)
      do j=mesh%domains(idom)%ez(1),mesh%domains(idom)%ez(2)
        grid%tag( sub2ind(i,j,mesh%nx) ) = mesh%domains(idom)%tag
      enddo
    enddo
  enddo
 ! tag the elements near the fault
 ! If ezflt=0 fault is at bottom
  if (mesh%fztag>0) then
    j1 = max(mesh%ezflt+1-mesh%fznz,1)
    j2 = min(mesh%ezflt+mesh%fznz,mesh%nz)
    do j=j1,j2
      do i=1,mesh%nx
        grid%tag( sub2ind(i,j,mesh%nx) ) = mesh%fztag
      enddo
    enddo
  endif
  if (any(grid%tag == 0)) call IO_abort('CART_build: Domain tags not entirely set')

 ! Control nodes of each element
 ! Elements are sequentially numbered horizontally from bottom-left to top-right
  call MESH_STRUCTURED_connectivity(grid%knods,mesh%nx,mesh%nz,grid%ngnod,mesh%ezflt)
  
 ! Boundary conditions
  splitN = 0
  if (mesh%ezflt>0) then
    allocate(grid%bnds(6))
  else
    if (mesh%split) splitN = floor( (mesh%splitD-mesh%xmin)/(mesh%xmax-mesh%xmin)*mesh%nx )
    if (splitN>0) then
      allocate(grid%bnds(5))  
    else 
      allocate(grid%bnds(4))
    endif
  endif
  call MESH_STRUCTURED_boundaries(grid%bnds,mesh%nx,mesh%nz,mesh%ezflt,splitN)

 ! Renumber elements
  if (OPT_RENUMBER) call MESH_STRUCTURED_renumber(grid,mesh%nx,mesh%nz)

end subroutine CART_build

subroutine CART_h(mesh,h)
  double precision, intent(out) :: h 
  type(mesh_cart_type), intent(in) :: mesh

  ! Devel Trevor: assumes fault on x boundary
  h = (mesh%xmax-mesh%xmin)/mesh%nx

end subroutine


end module mesh_cartesian
