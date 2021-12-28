module mesh_layers

! MESH_LAYERS: structured mesh generation for layered media

  use distribution_cd

!-- Renumbering of elements by Reverse Cuthill-McKee algorithm --
!   to improve data locality (optimize cache usage)
  use constants, only : OPT_RENUMBER, NDIME

  implicit none
  private

  type qc_spline_type
    integer :: N=0, kind=0
    double precision, pointer :: x(:)=>null(), y(:)=>null(), dy(:)=>null()
  end type qc_spline_type

  type layer_type
    integer :: nez=0, tag=0
    double precision, pointer :: top(:,:) => null()
    type (cd_type) :: ztopin
    type (qc_spline_type) :: qc_spline
  end type layer_type

  type mesh_layers_type
    private
    double precision :: xmin,xmax,zmin
    integer :: nlayer,ngnod,nz,nx,ezflt,fztag,fznz
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

  integer :: i,iin2,ngnod,nex_qc_spline
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

  nex_qc_spline = maxval(mesh%layer(:)%qc_spline%N)-1
  if (nex_qc_spline>0) then
    mesh%nx = nex_qc_spline
    do i=1,mesh%nlayer
      if (mesh%layer(i)%qc_spline%N>0 .and. mesh%layer(i)%qc_spline%N<nex_qc_spline+1) &
        call IO_abort('MESH_LAYERS_read: all qc_splines must have same data size')
    enddo
    write(iout,*) 
    write(iout,*) '    Reset nx to ',nex_qc_spline , ' (imposed by qc_spline)'
  endif

  mesh%nz=sum(mesh%layer%nez) 
  if (mesh%ezflt== -1) mesh%ezflt = mesh%nz/2
  if (mesh%ezflt >= mesh%nz) call IO_abort('MESH_LAYERS_read: ezflt must be < nz')

  if (echo_input) then
    write(eltype,'("Q",i1)') mesh%ngnod
    write(iout,200) mesh%nz,eltype
    if (mesh%ezflt>0) write(iout,220) mesh%ezflt
    if (mesh%fztag>0) write(iout,230) mesh%fztag,mesh%fznz
  endif

  return

200 format(/5x, &
      'Number of elements along Z. . . . . . . . . . . = ',i5,/5x, &
      'Element type. . . . . . . . . . . . . . . . . . = ',a)

220 format(5x, &
      'Fault on top of this element row  . . . (ezflt) = ',i0)

230 format(5x, &
      'Tag for elements in fault zone  . . . . (fztag) = ',I0,/5x, &
      'Vertical nb of elements in fault zone . .(fznz) = ',I0)

end subroutine MESH_LAYERS_read

!--------------------------------------------------
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_LAYERED
! GROUP  : MESH_DEF
! PURPOSE: Structured mesh for layered medium 
!          with surface and interface topography. 
! SYNTAX : &MESH_LAYERED xlim,zmin,nx,file,nlayer,ezflt,fztag /
!
! ARG: xlim     [dble(2)] [none] X limits of the box (min and max)
! ARG: zmin     [dble] [none] bottom Z limit of the box 
! ARG: nx       [int] [1]  Number of elements along the X direction.
!                Not needed if ztopH='QSPLINE' in a &MESH_LAYER block.
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
!                (ezflt+1)-th element rows, numbered from bottom to top. 
!                If ezflt=0 (default), no fault is introduced.
!                If ezflt=-1, a horizontal fault is introduced at/near the 
!                middle of the box: ezflt is reset to int[nelem(2)/2]
! ARG: fztag    [int][0] tag for elements near the fault
!                Useful to set a damping layer near the fault.
! ARG: fznz     [int][1] vertical size of near-fault layer
!                (half thickness in number of elements) 
!
! NOTE: the following tags are automatically assigned to the boundaries: 
!               1       Bottom 
!               2       Right        
!               3       Top  
!               4       Left
!               5       Fault, lower side
!               6       Fault, upper side
!
! END INPUT BLOCK

subroutine read_layered(iin,mesh,file)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort, IO_file_length

  type(mesh_layers_type), intent(out) :: mesh
  integer, intent(in) :: iin
  character(50), intent(out) :: file

  double precision :: init_double,xlim(2),zmin
  integer :: nx,nlayer,ezflt,fztag,fznz

  NAMELIST / MESH_LAYERED /  xlim,zmin,nx,nlayer,file,ezflt,fztag,fznz

  init_double = huge(init_double)

  xlim = init_double
  zmin = init_double
  nx = 1
  nlayer = 0
  file= ''
  ezflt = 0
  fztag = 0
  fznz = 1
  
  read(iin,MESH_LAYERED,END=100)

  if (xlim(1) == init_double) call IO_abort('MESH_LAYERS_read: you must set xlim')
  if (xlim(2) == init_double) call IO_abort('MESH_LAYERS_read: you must set xlim')
  if (zmin == init_double) call IO_abort('MESH_LAYERS_read: you must set zmin')
  if (nx <= 0) call IO_abort('MESH_LAYERS_read: nx must be positive')
  if (nlayer <= 0 .and. file=='') &
    call IO_abort('MESH_LAYERS_read: nlayer or file must be set')
  if (ezflt<-1) call IO_abort('MESH_LAYERS_read: ezflt must be >= -1')
  if (ezflt==0) fztag=0
  if (fztag<0) call IO_abort('MESH_LAYERS_read: fztag must be positive')
  if (fznz<1) call IO_abort('MESH_LAYERS_read: fznz must be strictly positive')

  if (file/='') nlayer = IO_file_length(file)

  if (echo_input) then
    write(iout,200) xlim, zmin, nx, nlayer
    if (file/='') write(iout,210) trim(file)
  endif

  mesh%xmin = xlim(1)
  mesh%xmax = xlim(2)
  mesh%zmin = zmin
  mesh%nx   = nx
  mesh%nlayer = nlayer
  mesh%ezflt = ezflt
  mesh%fztag = fztag
  mesh%fznz = fznz

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
! ARG: ztopH    [string] ['none'] Name of the type of spatial distribution to 
!                generate an irregular (non flat) top boundary. In general it is
!                one of the 1D distribution available through a DIST_XXXX block: 
!                  ztopH = 'LINEAR', or 
!                  ztopH = 'SPLINE', etc. 
!                There are two methods to generate a curve with a smooth normal, 
!                typically to guarantee smooth boundary conditions on curved faults.
!                The first method is based on quadratic splines and sometimes 
!                produces degenerated elements:
!                  ztopH='QSPLINE', followed by a &QC_SPLINE block
!                The second method is based on cubic splines and is more robust:
!                  ztopH='CSPLINE', followed by a &QC_SPLINE block
! ARG: tag      [int] [none]  Material tag
!                 If not given, a tag is automatically assigned to the layer, 
!                 sequentially numbered from top to bottom (top layer tag =1)
!
! NOTE: If ztopH='LINEAR' the mesh uses linearly deformed (Q4) elements,
!       otherwise it uses quadratically deformed (Q9) elements
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
    if (ztopH == 'LINEAR' .or. ztopH=='') then
      ngnod = 4
    elseif (ztopH == 'QSPLINE') then
      ngnod = 9
      layer%qc_spline%kind = 1
    elseif (ztopH == 'CSPLINE') then
      ngnod = 8
      layer%qc_spline%kind = 2
    else
      ngnod = 9
    endif

  else
    ztopH=''
    read(iin,*) ztop,nz,tag

  endif

  if (nz<=0) call IO_abort('MESH_LAYERS_read: nz null, negative or missing')
  layer%nez = nz

  if (ztopH=='QSPLINE' .or. ztopH=='CSPLINE') then
    call QC_SPLINE_read(layer%qc_spline,iin)
  else
    call DIST_CD_Read(layer%ztopin,ztop,ztopH,iin,ztopH)
  endif

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
!

subroutine MESH_LAYERS_build(mesh,grid)

  use fem_grid, only : fem_grid_type
  use memory_info
  use stdio, only : IO_abort
  use mesh_structured
  use utils, only : sub2ind

  type(mesh_layers_type), intent(inout) :: mesh
  type(fem_grid_type), intent(inout) :: grid

  double precision, pointer :: bot(:,:),top(:,:),coord(:,:),ztop(:)
  integer :: nxp,nzp,i,j,j1,j2,k,ilast,ifirst,nj,jfirst,jlast,tag,idoublingx,idoublingz
  integer :: jposlast,jposflt

  if (mesh%ngnod==4) then
    idoublingx = 1
    idoublingz = 1
  elseif (mesh%ngnod==9) then
    idoublingx = 2
    idoublingz = 2
  elseif (mesh%ngnod==8) then
    idoublingx = 3
    idoublingz = 1
  endif
  nxp = idoublingx*mesh%nx+1
  nzp = idoublingz*mesh%nz+1
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
  allocate(coord(NDIME,nxp))
  coord(1,:) = mesh%xmin +(mesh%xmax-mesh%xmin)/dble(nxp-1) *(/ (i, i=0,nxp-1) /)
  coord(2,:) = 0d0 ! not used
  do k=1,mesh%nlayer
    allocate(mesh%layer(k)%top(2,nxp))
    if (mesh%layer(k)%qc_spline%kind==1) then
      call QSPLINE_make( mesh%layer(k)%qc_spline, mesh%layer(k)%top(1,:), mesh%layer(k)%top(2,:))
    elseif (mesh%layer(k)%qc_spline%kind==2) then
      call CSPLINE_make( mesh%layer(k)%qc_spline, mesh%layer(k)%top(1,:), mesh%layer(k)%top(2,:))
    else
      mesh%layer(k)%top(1,:) = coord(1,:)
      ztop => mesh%layer(k)%top(2,:)
      call DIST_CD_Init(mesh%layer(k)%ztopin,coord, ztop)
    endif
    if (any(mesh%layer(k)%top(2,:)<mesh%zmin)) &
      call IO_abort('MESH_LAYERS_build: zmin in MESH_LAYERED should be deeper than all layer interfaces')
  enddo

 ! bottom line of nodes
  ifirst = 1
  ilast = nxp
  grid%coord(1,ifirst:ilast) = coord(1,:)
  grid%coord(2,ifirst:ilast) = mesh%zmin
  top => grid%coord(:,ifirst:ilast)
  jposlast = 1

 ! for each layer: lines of nodes from bottom (not included) to top (included)
  do k=1,mesh%nlayer
    bot => top 
    top => mesh%layer(k)%top
    nj = idoublingz*mesh%layer(k)%nez
    jposflt = idoublingz*mesh%ezflt+1 -jposlast
    do j=1,nj
      ifirst = ilast + 1
      ilast  = ilast + nxp
      grid%coord(:, ifirst:ilast ) = bot+ dble(j)/dble(nj)*(top-bot)
      if (j==jposflt) then
        ifirst = ilast + 1
        ilast  = ilast + nxp
        grid%coord(:, ifirst:ilast ) = bot+ dble(j)/dble(nj)*(top-bot)
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
 ! tag the elements near the fault
  if (mesh%fztag>0) then
    j1 = max(mesh%ezflt+1-mesh%fznz,1)
    j2 = min(mesh%ezflt+mesh%fznz,mesh%nz)
    do j=j1,j2
      do i=1,mesh%nx
        grid%tag( sub2ind(i,j,mesh%nx) ) = mesh%fztag
      enddo
    enddo
  endif
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

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : QC_SPLINE
! GROUP  : MESH_LAYER
! PURPOSE: Define the boundary of a layer using quadratic or cubic splines and
!          enforcing smooth (continuous) normal between elements, for instance
!          to guarantee smooth boundary conditions on curved faults.
!          
! SYNTAX : &QC_SPLINE file /
!
! ARG: file     [string] [''] Name of ASCII file containing information of
!                all the element vertex nodes lying on the boundary curve.
!                One line per node, ordered by increasing x, 3 columns per line:
!                  (1) x position
!                  (2) z position
!                  (3) derivative dz/dx of the curve at the node
!                All QC_SPLINE curves in a mesh must have the same number of nodes.
!                The parameter nx in &MESH_LAYERED is automatically reset
!                (nx = number of nodes in QC_SPLINE - 1)
!
! END INPUT BLOCK

subroutine QC_SPLINE_read(q,iin)

  use stdio, only : IO_file_length, IO_new_unit

  type (qc_spline_type), intent(inout) :: q
  integer, intent(in) :: iin

  character(50) :: file
  integer :: iunit, i

  NAMELIST / QC_SPLINE / file

  read(iin,QC_SPLINE)

  q%N = IO_file_length(file)
  allocate( q%x(q%N),q%y(q%N), q%dy(q%N) )
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  do i= 1,q%N
    read (iunit,*) q%x(i),q%y(i),q%dy(i)
  end do
  close(iunit)

end subroutine QC_SPLINE_read

!---------------------------------------------------------------------
subroutine QSPLINE_make(q,x,y)

  use stdio, only : IO_abort

  type (qc_spline_type), intent(in) :: q
  double precision, intent(out) :: x(:),y(:)

  double precision :: x1,x2,y1,y2,d1,d2
  integer :: i,k0,nspline

  nspline = 2*q%N-1  ! total number of nodes in the spline curve

  if (size(x)/=nspline .or. size(y)/=nspline) &
    call IO_abort('mesh_layers:QSPLINE_make: arguments have inconsistent size')

  do i=1,q%N-1

    x1=q%x(i)
    x2=q%x(i+1)
    y1=q%y(i)
    y2=q%y(i+1)
    d1=q%dy(i)
    d2=q%dy(i+1)

    k0 = 2*(i-1)

    x(k0+1) = x1
    y(k0+1) = y1

    if (d1==d2) then
      x(k0+2) = 0.5d0*(x1+x2)
      y(k0+2) = 0.5d0*(y1+y2)
    else
!      x(k0+2) = ( y2-y1 +(3d0*d1-d2)*x1*0.5d0 +(d1-3d0*d2)*x2*0.5d0 )/(d1-d2)*0.5d0
!      y(k0+2) = ( d2*d1*(x2-x1) +(3d0*d2-d1)*y1*0.5d0 + (d2-3d0*d1)*y2*0.5d0 )/(d2-d1)*0.5d0
      x(k0+2) = 0.5d0*(x1+x2) + 0.5d0*( y2-y1 +0.5d0*(d1+d2)*(x1 -x2) )/(d1-d2)
      y(k0+2) = 0.5d0*(y1+y2) + 0.5d0*( d2*d1*(x2-x1) +0.5d0*(d2+d1)*(y1-y2) )/(d2-d1)
    endif

    x(k0+3) = x1
    y(k0+3) = y1

  enddo

end subroutine QSPLINE_make

!---------------------------------------------------------------------
subroutine CSPLINE_make(c,x,y)

  use stdio, only : IO_abort

  type (qc_spline_type), intent(in) :: c
  double precision, intent(out) :: x(:),y(:)

  double precision :: x1,x2,y1,y2,yx1,yx2
  integer :: i,k0,nspline

  nspline = 3*(c%N-1)+1  ! total number of nodes in the spline curve

  if (size(x)/=nspline .or. size(y)/=nspline) &
    call IO_abort('mesh_layers:CSPLINE_make: arguments have inconsistent size')

 ! for each element compute y for the two interior nodes
  do i=1,c%N-1  

   ! constraints
    x1=c%x(i)
    y1=c%y(i)
    yx1=c%dy(i)
    x2=c%x(i+1)
    y2=c%y(i+1)
    yx2=c%dy(i+1)

    k0 = 3*(i-1)

    x(k0+1) = x1
    y(k0+1) = y1

    x(k0+2) = (2d0*x1+x2)/3d0
    y(k0+2) = ( (20d0*y1+7d0*y2)-(x2-x1)*(2d0*yx2-4d0*yx1) )/27d0

    x(k0+3) = (x1+2d0*x2)/3d0
    y(k0+3) = ( (7d0*y1+20d0*y2)-(x2-x1)*(4d0*yx2-2d0*yx1) )/27d0

    x(k0+4) = x2
    y(k0+4) = y2

  enddo

end subroutine CSPLINE_make

end module mesh_layers
