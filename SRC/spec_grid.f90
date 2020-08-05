module spec_grid

  use stdio, only : IO_abort
  use constants
  use fem_grid
  use bnd_grid

  implicit none
  private


!-----------------------------------------------------------------------
!-- Spectral element grid type
!
!------ Topology
!     nelem        = total number of elements
!     ngll         = number of Gauss-Lobato-Legendre nodes per unit segment
!     npoin        = total number of GLL nodes in the mesh
!     ibool        = local to global numbering for the GLL mesh
!                    (igll,jgll,element) -> bulk node index
!     coord        = coordinates of the GLL nodes
!
!------ Working data
!     shape        = control nodes to GLL nodes shape functions
!     dshape       = derivatives of shape functions
!     hprime       = lagrange derivatives on the unit square
!     hTprime      = transpose of hprime
!     wgll         = 1D GLL integration weights
!     wgll2        = 2D GLL integration weights
!     xgll         = GLL points on the unit segment
!
!------ Other
!     tag          = element -> domain tags
!     fmax         = Highest frequency to be resolved by the grid
!     W            = The seismogenic width. Infinity means 2D problem 
!                    and finite W means 2.5D problem (use for elastic material)

!
!
! NOTE: After initialization, the FEM mesh database is only used for plotting
!       purposes, so if no plots are needed you can save them on disk and 
!       clean it from memory.
!
! NOTE: COORD is used in the solver only for computation of the incident
!       wavefield, the boundary coordinates could be stored instead in the 
!       bc structure (boundary_condition).


  type sem_grid_type
    integer :: ngll=0,nelem=0,npoin=0
    type(fem_grid_type) :: fem
    double precision :: fmax=0d0, W = huge(1d0)
    double precision, pointer :: coord(:,:)      =>null(), &
                                 hprime(:,:)     =>null(), &
                                 hTprime(:,:)    =>null(), &
                                 wgll(:)         =>null(), &
                                 xgll(:)         =>null(), &
                                 wgll2(:,:)      =>null(), &
                                 shape(:,:,:)    =>null(), &
                                 dshape(:,:,:,:) =>null()
    integer, pointer :: ibool(:,:,:)  =>null(), &
                        tag(:) =>null()
    type (bnd_grid_type), pointer :: bounds(:) =>null()
  end type sem_grid_type

!-----------------------------------------------------------------------

  interface SE_Jacobian
    module procedure SE_Jacobian_e, SE_Jacobian_eij
  end interface SE_Jacobian

  interface SE_InverseJacobian
    module procedure SE_InverseJacobian_e, SE_InverseJacobian_eij
  end interface SE_InverseJacobian

  interface SE_VolumeWeights
    module procedure SE_VolumeWeights_e, SE_VolumeWeights_eij
  end interface SE_VolumeWeights

  interface SE_node_belongs_to
    module procedure SE_node_belongs_to_1, SE_node_belongs_to_2
  end interface SE_node_belongs_to

  ! (i,j) GLL for element corner nodes
  !icorner(1) = 1    ; jcorner(1) = 1
  !icorner(2) = ngll ; jcorner(2) = 1
  !icorner(3) = ngll ; jcorner(3) = ngll
  !icorner(4) = 1    ; jcorner(4) = ngll


  public :: sem_grid_type ,&
         SE_init ,&
         SE_init_interpol ,&
         SE_find_nearest_node ,&
         SE_node_belongs_to, &
         SE_find_point ,&
         SE_get_edge_nodes ,&
         SE_inquire , &
         SE_isFlat, &
         SE_firstElementTagged, &
         SE_VolumeWeights ,&
         SE_InverseJacobian ,&
         SE_Jacobian ,&
         SE_elem_coord, &
         BC_inquire ,&
         BC_tag_exists, &
         BC_get_normal_and_weights &
          , edge_D,edge_R,edge_U,edge_L
 
contains

!=======================================================================
!
  logical function SE_isFlat(se)
  type (sem_grid_type), intent(in) :: se
  SE_isFlat = se%fem%flat 
  end function SE_isFlat

!=======================================================================
!
  subroutine SE_init(se)

  use stdio, only : IO_new_unit
  type (sem_grid_type), intent(inout) :: se
  integer :: ounit

  call SE_init_gll(se) ! GLL quadrature points, weights 
                            ! and lagrange coefficients,
                            ! shape functions and their derivatives
  call SE_init_numbering(se) ! generate the global node numbering
  call SE_init_coord(se) ! get the coordinates of the GLL nodes,
                              ! and export the grid to a file
  call SE_BcTopoInit(se) ! initialize boundaries generic data

  ! export grid parameters
  ounit = IO_new_unit()
  open(ounit,file='grid_sem2d.hdr',status='replace')
  write(ounit,*) 'NELEM  NPGEO  NGNOD  NPOIN  NGLL'
  write(ounit,*) se%nelem, FE_GetNbNodes(se%fem), FE_GetNodesPerElement(se%fem), &
                 se%npoin, se%ngll
  close(ounit)

  end subroutine SE_init

!=======================================================================
!
! get coordinates and weights of the Gauss-Lobatto-Legendre points
! and derivatives of Lagrange polynomials  H_ij = h'_i(xgll(j))
  subroutine SE_init_gll(se)

  use gll, only : get_GLL_info, print_GLL

  type (sem_grid_type), intent(inout) :: se
  double precision :: xi,eta 
  integer :: ngll,i,j, ngnod

  ngll = se%ngll

  allocate(se%xgll(ngll))
  allocate(se%wgll(ngll))
  allocate(se%hprime(ngll,ngll))
  allocate(se%hTprime(ngll,ngll))
  allocate(se%wgll2(ngll,ngll))

  call get_GLL_info(ngll,se%xgll,se%wgll,se%hprime)
  call print_GLL(ngll,se%xgll,se%wgll,se%hprime)
  se%hTprime = transpose(se%hprime)
  
 ! wgll2(i,j) = wgll(i) * wgll(j)
  do j = 1,ngll
    se%wgll2(:,j) = se%wgll * se%wgll(j)
  enddo

!-----------------------------------------------------------------------
!-- set up the shape functions and their local derivatives
  ngnod = FE_GetNodesPerElement(se%fem)

  allocate(se%shape(ngnod,ngll,ngll))
  allocate(se%dshape(ngnod,NDIME,ngll,ngll))

  do j = 1,ngll
    eta = se%xgll(j)
    do i = 1,ngll
      xi = se%xgll(i)
      se%shape(:,i,j) = FE_getshape(xi,eta,ngnod)
      se%dshape(:,:,i,j) = FE_getdershape(xi,eta,ngnod)
    enddo
  enddo

  end subroutine SE_init_gll



!=======================================================================
!
!-- generate the global numbering
!
  subroutine SE_init_numbering(se)

  use memory_info
  use echo, only : echo_init,iout,fmt1,fmtok
  use stdio, only : IO_new_unit

  type(sem_grid_type), intent(inout) :: se

  integer, pointer :: ibool(:,:,:)
  integer, dimension(se%ngll,4) :: iedg,jedg,iedgR,jedgR
  integer, dimension(4) :: ivtx,jvtx
  integer :: i,j,k,e,n,ii,jj,ee,nn,ngll,npoin,iol,ounit
  integer, pointer :: ees(:),nns(:)

!-----------------------------------------------------------------------

  if (echo_init) then
    write(iout,*) 
    write(iout,'(a)') ' S p e c t r a l   e l e m e n t s   g r i d'
    write(iout,'(a)') ' ==========================================='
    write(iout,*) 
  endif

  ngll  = se%ngll
  se%nelem = FE_GetNbElements(se%fem)
  se%tag => se%fem%tag

! global numbering table
  allocate(se%ibool(ngll,ngll,se%nelem))
  call storearray('ibool',size(se%ibool),iinteg)
  ibool => se%ibool
  ibool = 0

!---- start numbering
  if (echo_init) write(iout,fmt1,advance='no') 'Numbering GLL points'
  npoin = 0

 ! GLL index (i,j) of edge nodes, counterclockwise
  do n = 1,4
    call SE_inquire(se,edge=n,itab=iedg(:,n),jtab=jedg(:,n))
  enddo
 ! reverse order for matching edges
  iedgR = iedg(ngll:1:-1,:)
  jedgR = jedg(ngll:1:-1,:)

 ! GLL index (i,j) for vertex nodes
  ivtx(1) = 1;    jvtx(1) = 1
  ivtx(2) = ngll; jvtx(2) = 1
  ivtx(3) = ngll; jvtx(3) = ngll
  ivtx(4) = 1;    jvtx(4) = ngll

  do e = 1,se%nelem


   !-- interior nodes are unique
    do j=2,ngll-1
    do i=2,ngll-1
      npoin = npoin + 1
      ibool(i,j,e) = npoin
    enddo
    enddo
  
   !-- edge nodes
    do n = 1,4
      if ( ibool(iedg(2,n),jedg(2,n),e)>0 ) cycle ! skip if already processed
      call FE_GetEdgeConn(se%fem,e,n,ee,nn)
      do k = 2,ngll-1
        npoin = npoin+1
        ibool(iedg(k,n),jedg(k,n),e) = npoin
        if (ee>0) ibool(iedgR(k,nn),jedgR(k,nn),ee) = npoin
      enddo
    enddo


   !-- vertex nodes
    do n = 1,4
      i = ivtx(n)
      j = jvtx(n)
      if ( ibool(i,j,e) > 0 ) cycle
      npoin = npoin +1
     ! scan Vertex_Conn for shared elements
      call FE_GetVertexConn(se%fem,e,n,ees,nns)
      do k= 1,size(ees)
        ii = ivtx(nns(k))
        jj = jvtx(nns(k))
        ibool(ii,jj,ees(k)) = npoin
      enddo
    enddo 

  enddo

  call FE_UnsetConnectivity(se%fem)

  se%npoin = npoin

  if (echo_init) then
    write(iout,fmtok)
    write(iout,100) 'Total number of elements. . . . . . . . = ',se%nelem
    write(iout,100) 'Total number of GLL points. . . . . . . = ',npoin
    write(iout,*)
    write(iout,fmt1,advance='no') 'Saving element/node table in binary file ibool_sem2d.dat'
  endif

  inquire( IOLENGTH=iol ) se%ibool(:,:,1)
  ounit = IO_new_unit()
  open(ounit,file='ibool_sem2d.dat',status='replace',access='direct',recl=iol)
  do e = 1,se%nelem
    write(ounit,rec=e) se%ibool(:,:,e)
  enddo
  close(ounit)
  if (echo_init) write(iout,fmtok)

  return

100 format(5X,A,I0)
  
  end subroutine SE_init_numbering

  
!=======================================================================
!
!  set the global nodal coordinates
!
  subroutine SE_init_coord(se)

  use echo, only : echo_check,echo_init,iout,fmt1,fmtok
  use memory_info
  use stdio, only : IO_new_unit

  type (sem_grid_type), intent(inout) :: se

  integer :: i,j,e,ounit,iol
  double precision, pointer :: coorg(:,:)

  allocate(se%coord(NDIME,se%npoin))
  call storearray('coord',size(se%coord),idouble)

!---- Coordinates of the global points 
  if (echo_init) write(iout,fmt1,advance='no') 'Defining nodes coordinates'
  do e = 1,se%nelem
    coorg => FE_GetElementCoord(se%fem,e)
    do j = 1,se%ngll
    do i = 1,se%ngll
      se%coord(:,se%ibool(i,j,e)) = matmul( coorg, se%shape(:,i,j) )
    enddo
    enddo
    deallocate(coorg)
  enddo
  if (echo_init) write(iout,fmtok)

  if (echo_check) then

!----  Save the grid in a text file
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Saving the grid coordinates (coord) in a text file'
    ounit = IO_new_unit()
    open(ounit,file='coord_sem2d.tab',status='unknown')
    write(ounit,*) se%npoin
    do i = 1,se%npoin
      write(ounit,*) i, se%coord(:,i)
    enddo
    close(ounit)
    write(iout,fmtok)
  endif
    
!----  Always save the grid in a binary file
  if (echo_init) write(iout,fmt1,advance='no') 'Saving the grid coordinates (coord) in a binary file'
  inquire( IOLENGTH=iol ) real(se%coord(:,1))
  ounit = IO_new_unit()
  open(ounit,file='coord_sem2d.dat',status='replace',access='direct',recl=iol)
  do i = 1,se%npoin
    write(ounit,rec=i) real(se%coord(:,i))
  enddo
  close(ounit)
  if (echo_init) write(iout,fmtok)

  end subroutine SE_init_coord



!=======================================================================
!
! lagrange interpolation functions at (xi,eta)
! interp(ngll*ngll) is then used as
!   interpolated_value = matmul(interp,nodal_values)
!   interpolated_vector(1:ndof) = matmul(interp,nodal_vectors(:,1:ndof))
!
  subroutine SE_init_interpol(xi,eta,interp,grid)

  use gll, only : hgll

  double precision, intent(in) :: xi,eta
  type (sem_grid_type), intent(in) :: grid
  double precision, intent(out) :: interp(grid%ngll*grid%ngll)

  double precision :: fi,fj
  integer :: i,j,k

  k=0
  do j=1,grid%ngll
    fj = hgll(j-1,eta,grid%xgll,grid%ngll)
    do i=1,grid%ngll
      k=k+1
      fi = hgll(i-1,xi,grid%xgll,grid%ngll)
      interp(k) = fi*fj
    enddo
  enddo

  end subroutine SE_init_interpol


!=====================================================================
!
  subroutine SE_find_nearest_node(coord_in,grid,iglob,coord,distance)

  type(sem_grid_type), intent(in)  :: grid
  double precision   , intent(in)  :: coord_in(:)
  double precision, optional, intent(out) :: coord(:),distance
  integer, intent(out) :: iglob

  double precision :: d2min,d2
  integer :: ip

  iglob = 0
  d2min = huge(d2min)
  do ip = 1,grid%npoin
    d2 = (coord_in(1)-grid%coord(1,ip))**2 + (coord_in(2)-grid%coord(2,ip))**2
    if (d2 <= d2min) then
      d2min  = d2
      iglob = ip
    endif
  enddo

  if (present(coord)) coord = grid%coord(:,iglob)
  if (present(distance)) distance = sqrt(d2min)

  end subroutine SE_find_nearest_node

!=====================================================================
! Find (e,i,j), element and local indices associated to a global node iglob
! Version 1: first element found
! Version 2: all elements
!
  subroutine SE_node_belongs_to_1(iglob,e,i,j,grid)

  integer, optional, intent(in) :: iglob
  integer, intent(out) :: i,j,e
  type(sem_grid_type), intent(in)  :: grid

  if (iglob>grid%npoin) call IO_abort('SE_node_belongs_to: iglob out of range')

  do e= 1,grid%nelem
  do j= 1,grid%ngll
  do i= 1,grid%ngll
    if (grid%ibool(i,j,e)==iglob) return
  enddo
  enddo
  enddo

  end subroutine SE_node_belongs_to_1

!---------------------------------------------------------------------

  subroutine SE_node_belongs_to_2(iglob,etab,itab,jtab,grid)

  integer, optional, intent(in) :: iglob
  integer, pointer, dimension(:) :: itab,jtab,etab
  type(sem_grid_type), intent(in)  :: grid

  integer :: i,j,e,nel
  ! WARNING: dimension should be larger than the maximum expected
  ! number of neighbour elements 
  integer, dimension(10) :: etmp,jtmp,itmp

  if (iglob>grid%npoin) call IO_abort('SE_node_belongs_to: iglob out of range')

  nel = 0
  do e= 1,grid%nelem
  do j= 1,grid%ngll
  do i= 1,grid%ngll
    if (grid%ibool(i,j,e)==iglob) then
      nel = nel +1
      etmp(nel) = e
      itmp(nel) = i
      jtmp(nel) = j
    endif
  enddo
  enddo
  enddo

  if (associated(etab)) deallocate(etab)
  if (associated(itab)) deallocate(itab)
  if (associated(jtab)) deallocate(jtab)
  allocate(etab(nel),itab(nel),jtab(nel))
  etab = etmp(1:nel)
  itab = itmp(1:nel)
  jtab = jtmp(1:nel)

  end subroutine SE_node_belongs_to_2

!=====================================================================
!
  subroutine SE_find_point(coord,grid,e,xi,eta,new_coord)

  double precision, intent(in)  :: coord(NDIME)
  type(sem_grid_type), intent(in) :: grid
  integer, intent(out) :: e
  double precision, intent(out) :: xi,eta,new_coord(NDIME)

  integer :: iglob,k,istat
  integer, pointer, dimension(:) :: itab,jtab,etab

  nullify(itab,jtab,etab)
  call SE_find_nearest_node(coord,grid,iglob)
  call SE_node_belongs_to(iglob,etab,itab,jtab,grid)
  do k = 1,size(etab)
    xi=grid%xgll(itab(k))
    eta=grid%xgll(jtab(k))
    e = etab(k)
    call FE_find_point(coord,grid%fem,e,xi,eta,new_coord,istat)
    if (istat==0) exit
  enddo
  
  deallocate(itab,jtab,etab)
  if (istat>0) call IO_abort('SE_find_point: point not found')

  end subroutine SE_find_point

!=====================================================================
! get the list of bulk indices of the EDGE nodes of an ELEMENT, 
! assume counterclockwise numbering for the output list NODES

  subroutine SE_get_edge_nodes(grid,element,edge,nodes)

  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: element,edge
  integer, intent(out) :: nodes(grid%ngll)

  select case (edge)
    case(edge_D); nodes = grid%ibool(:,1,element)
    case(edge_R); nodes = grid%ibool(grid%ngll,:,element)
    case(edge_U); nodes = grid%ibool(grid%ngll:1:-1,grid%ngll,element)
    case(edge_L); nodes = grid%ibool(1,grid%ngll:1:-1,element)
  end select

  end subroutine SE_get_edge_nodes


!=======================================================================
!
! Jacobian matrix  =  | dx/dxi   dx/deta |
!                     | dz/dxi   dz/deta |
!
  function SE_Jacobian_e(sem,e) result(jac)

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e
  double precision :: jac(2,2,sem%ngll,sem%ngll)

  double precision, pointer :: coorg(:,:)
  integer :: i,j

  coorg => FE_GetElementCoord(sem%fem,e)
  do j = 1,sem%ngll
  do i = 1,sem%ngll
    jac(:,:,i,j) = matmul(coorg, sem%dshape(:,:,i,j) )
  enddo
  enddo
  deallocate(coorg)

  end function SE_Jacobian_e

!-----------------------------------------------------------------------

  function SE_Jacobian_eij(sem,e,i,j) result(jac)

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e,i,j
  double precision :: jac(2,2)
  double precision, pointer :: coorg(:,:)

  coorg => FE_GetElementCoord(sem%fem,e)
  jac = matmul(coorg, sem%dshape(:,:,i,j) )
  deallocate(coorg)

  end function SE_Jacobian_eij


!=======================================================================
!
!-- Compute weights(i,j) = wgll(i)*wgll(j)*dvol(i,j) for one element
!
  function SE_VolumeWeights_e(sem,e) result(dvol)

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e
  double precision :: dvol(sem%ngll,sem%ngll)

  double precision :: jac(2,2,sem%ngll,sem%ngll)

  jac = SE_Jacobian_e(sem,e)
  dvol = jac(1,1,:,:)*jac(2,2,:,:) - jac(1,2,:,:)*jac(2,1,:,:)
  dvol = dvol*sem%wgll2

  end function SE_VolumeWeights_e

!-----------------------------------------------------------------------

  function SE_VolumeWeights_eij(sem,e,i,j) result(dvol)

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e,i,j
  double precision :: dvol

  double precision :: jac(2,2)

  jac = SE_Jacobian_eij(sem,e,i,j)
  dvol = ( jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1) )*sem%wgll2(i,j)

  end function SE_VolumeWeights_eij


!=======================================================================
!
! Inverse Jacobian matrix  =  | dxi/dx    dxi/dz |
!                             | deta/dx  deta/dz |
!
  function SE_InverseJacobian_e(sem,e) result(jaci)

  use utils, only : invert2

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e
  double precision :: jaci(2,2,sem%ngll,sem%ngll)

  double precision :: jac(2,2,sem%ngll,sem%ngll)
  integer :: i,j

  jac = SE_Jacobian_e(sem,e)
  do j = 1,sem%ngll
  do i = 1,sem%ngll
    jaci(:,:,i,j) = invert2(jac(:,:,i,j))
  enddo
  enddo

  end function SE_InverseJacobian_e

!----------------------------------------------------------------------

  function SE_InverseJacobian_eij(sem,e,i,j) result(jaci)

  use utils, only : invert2

  type(sem_grid_type), intent(in) :: sem
  integer, intent(in) :: e,i,j
  double precision :: jaci(2,2)

  double precision :: jac(2,2)

  jac = SE_Jacobian_eij(sem,e,i,j)
  jaci = invert2(jac)

  end function SE_InverseJacobian_eij


!=======================================================================
!
subroutine SE_BcTopoInit(grid)
  type(sem_grid_type), intent(inout) :: grid
  integer :: i
  grid%bounds => grid%fem%bnds
  do i=1,size(grid%bounds)
    call BC_set_bulk_node(grid%bounds(i),grid)
  enddo
end subroutine SE_BcTopoInit

!=======================================================================
!
subroutine BC_set_bulk_node(bc,grid)

  use generic_list, only : Link_Ptr_Type,Link_Type,List_Type &
     ,LI_Init_List,LI_Add_To_Head,LI_Get_Head &
     ,LI_Get_Next,LI_Associated,LI_Remove_Head

  use utils, only: drank

  type(bnd_grid_type), intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  
  type BC_Node_Data_Type
    integer :: bulk_node,bc_node
  end type BC_Node_Data_Type

  type BC_Node_Type
    TYPE(Link_Type) :: Link
    type(BC_Node_Data_Type), pointer :: Data
  end type BC_Node_Type

  type BC_Node_Ptr_Type
    type(BC_Node_Type), pointer :: P
  end type BC_Node_Ptr_Type

  integer :: bulk_node,bulk_node_vect(grid%ngll),i,n,kloc,bc_npoin,bc_inode
  logical :: node_at_corner,new_node
  type(List_Type)      :: BC_Nodes_List,BC_Corners_List
  TYPE(Link_Ptr_Type)  :: Link,Sublink
  TYPE(BC_Node_Ptr_Type)  :: BC_node,BC_corner

! For reordering:
  double precision, allocatable :: BcCoord(:,:)
  integer, allocatable :: Isort(:),Iback(:)
  double precision :: Lx,Lz

  bc%ngnod = grid%ngll
  allocate(bc%ibool(bc%ngnod,bc%nelem))

!-----------------------------------------------------------------------
!  Build the database as a linked list

  bc_npoin = 0
  call LI_Init_List(BC_Nodes_List)
  call LI_Init_List(BC_Corners_List)

  do n=1,bc%nelem

!   get edge nodes, counterclockwise by default
    call SE_get_edge_nodes(grid,bc%elem(n),bc%edge(n), bulk_node_vect)

    do kloc = 1,bc%ngnod

      bulk_node = bulk_node_vect(kloc)
      new_node = .true.

!  When the node is an element corner (vertex) 
!  check if it is already in the BC list.
!  Note that the search is done in a BC_corners sublist.

      node_at_corner = kloc==1 .or. kloc==bc%ngnod
      if (node_at_corner) then
        SubLink = LI_Get_Head(BC_Corners_List)
        do while(LI_Associated(SubLink))
          BC_corner = transfer(SubLink,BC_corner)
          if (BC_corner%P%Data%bulk_node == bulk_node) then
            new_node = .false.
            bc_inode = BC_corner%P%Data%bc_node
            exit
          endif
          Sublink = LI_Get_Next(Sublink)
        enddo
      endif

!  If it is a new node, add it to the list of BC nodes ...
      if (new_node) then 

        bc_npoin = bc_npoin +1
        bc_inode = bc_npoin
        allocate(BC_node%P)
        allocate(BC_node%P%data)
        BC_node%P%data%bulk_node = bulk_node
        BC_node%P%data%bc_node = bc_inode
        Link = transfer(BC_node,Link)
        call LI_Add_To_Head(Link,BC_Nodes_List)
         
!  ... and possibly to the sublist of BC corners.
        if (node_at_corner) then
          allocate(BC_corner%P)
          BC_corner%P%data => BC_node%P%data
          SubLink = transfer(BC_corner,SubLink)
          call LI_Add_To_Head(SubLink,BC_Corners_List)
        endif

      endif

     ! Set the [ (gll,bc_element) -> bc_node ] table
      bc%ibool(kloc,n) = bc_inode

    enddo
  enddo

 ! clean up the list of corners
  do 
    SubLink = LI_Remove_Head(BC_Corners_List)
    if (.not.LI_Associated(SubLink)) exit
    BC_corner = transfer(SubLink,BC_corner)
    deallocate(BC_corner%P)
  enddo

!-----------------------------------------------------------------------
!  Translate BC database from linked list to array storage

  bc%npoin = bc_npoin ! = LI_Get_Len(BC_Nodes_List)
  allocate(bc%node(bc%npoin)) 
  do i=1,bc%npoin
    Link = LI_Remove_Head(BC_Nodes_List)
    BC_node = transfer(Link,BC_node)
    bc%node(BC_node%P%data%bc_node) = BC_node%P%data%bulk_node
    deallocate(BC_node%P%data)
    deallocate(BC_node%P)
  enddo

!-----------------------------------------------------------------------
!  Sort BC nodes by increasing coordinate
!  Use the coordinate with the largest range
  allocate(BcCoord(2,bc%npoin))
  allocate(Isort(bc%npoin))
  BcCoord = grid%coord(:,bc%node)
  Lx = maxval(BcCoord(1,:)) - minval(BcCoord(1,:))
  Lz = maxval(BcCoord(2,:)) - minval(BcCoord(2,:))
  if (Lx>Lz) then
    call drank(BcCoord(1,:),Isort)
  else
    call drank(BcCoord(2,:),Isort)
  endif
!    write(51,'(I5,1X,I3)') ( bc%node(n),Isort(n),n=1,bc%npoin) 
  bc%node = bc%node(Isort)
!    write(51,'(I5)') bc%node
  allocate(Iback(bc%npoin))
  Iback(Isort) = (/ (i, i=1,bc%npoin) /) 
  do n=1,bc%nelem
    bc%ibool(:,n) = Iback(bc%ibool(:,n))
  enddo

  deallocate(BcCoord,Isort,Iback)

end subroutine BC_set_bulk_node


!=======================================================================
!
! example:
!
!  integer, pointer :: ptr(:)
!  ptr => SE_firstElementTagged(grid)
!  ...
!  deallocate(ptr)
!
  function SE_firstElementTagged(grid) result(elist)

  use utils, only : unique

  type (sem_grid_type), intent(in) :: grid
  integer, pointer :: elist(:)

  integer, pointer :: tags(:)
  integer :: ntags,k,e

  tags => unique(grid%tag)
  ntags = size(tags)

  allocate(elist(ntags))
  do k=1,ntags
    do e=1,grid%nelem
      if (grid%tag(e)==tags(k)) exit
    enddo
    elist(k) = e
  enddo

  deallocate(tags)
  
  end function SE_firstElementTagged

!=======================================================================
!
  subroutine SE_inquire(grid,element,edge &
                       ,itab,jtab,dim_t &
                       ,size_min,size_max)

  type (sem_grid_type), intent(in) :: grid
  integer, optional, intent(in) :: element,edge
  double precision, optional, intent(out) :: size_min,size_max
  integer, dimension(grid%ngll), optional, intent(out) :: itab,jtab
  integer, optional, intent(out) :: dim_t

  integer :: k
  double precision :: dmin,dmax

  !-- give the GLL local numbering ITAB,JTAB for a given EDGE
  !   assuming counterclockwise orientation
  if ( present(edge) .and. present(itab) .and. present(jtab) ) then
    select case(edge)
      case(edge_D); itab = (/ (k, k=1,grid%ngll) /)   ; jtab = 1
      case(edge_R); itab = grid%ngll                  ; jtab = (/ (k, k=1,grid%ngll) /)
      case(edge_U); itab = (/ (k, k=grid%ngll,1,-1) /); jtab = grid%ngll
      case(edge_L); itab = 1                          ; jtab = (/ (k, k=grid%ngll,1,-1) /)
    end select
  endif

  !-- give local dimension (xi,eta) of tangent direction to an EDGE
  if ( present(dim_t) ) then
    select case(edge)
      case(edge_D,edge_U); dim_t = 1
      case(edge_R,edge_L); dim_t = 2
    end select
  endif

  if (present(size_min) .or. present(size_max)) then
    call FE_GetElementSizes(dmax,dmin,element,grid%fem)
    if (present(size_min)) size_min = dmin
    if (present(size_max)) size_max = dmax
  endif
    
  end subroutine SE_inquire

!=======================================================================
!
  function SE_elem_coord(grid,e) result(ecoord)

  type (sem_grid_type), intent(in) :: grid
  integer, intent(in) :: e

  double precision :: ecoord(2,grid%ngll,grid%ngll)
  integer :: i,j

  do j=1,grid%ngll
  do i=1,grid%ngll
    ecoord(:,i,j) = grid%coord(:,grid%ibool(i,j,e))
  enddo
  enddo

  end function SE_elem_coord

!=======================================================================
!
  subroutine BC_inquire(bounds,tag,bc_topo_ptr)

  type (bnd_grid_type), pointer :: bounds(:),bc_topo_ptr
  integer, intent(in) :: tag
  
  integer :: i
 
 ! if no boundary exists with the requested tag a null pointer is returned
  nullify(bc_topo_ptr) 
  do i=1,size(bounds)
    if ( bounds(i)%tag == tag ) then
      bc_topo_ptr => bounds(i)
      return
    endif
  enddo

  end subroutine BC_inquire

!=======================================================================
!
  logical function BC_tag_exists(bounds,tag) result(exists)

  type (bnd_grid_type), intent(in) :: bounds(:)
  integer, intent(in) :: tag

  integer :: i
  
  exists = .false.
  do i=1,size(bounds)
    if ( bounds(i)%tag == tag ) then
      exists = .true.
      exit
    endif
  enddo

  end function BC_tag_exists

!=======================================================================
! Normal to a boundary, pointing out of the element
! Normals are assembled --> "average" normal between elements
  subroutine BC_get_normal_and_weights(bc,grid,NORM,W,periodic)

  type (bnd_grid_type), intent(inout) :: bc
  type (sem_grid_type), intent(in) :: grid
  double precision, intent(out) :: NORM(bc%npoin,2)
  double precision, intent(out) :: W(bc%npoin)
  logical, intent(in) :: periodic

  double precision :: DGlobDLoc(2,2),DxzDt(2),Tang(2),Jac1D,SignTang
  integer :: iGLLtab(bc%ngnod),jGLLtab(bc%ngnod)
  integer :: BcElmt,kGLL,BcNode,LocDimTanToEdge

  NORM = 0d0
  W = 0d0
  do BcElmt = 1,bc%nelem
    call SE_inquire(grid, edge=bc%edge(BcElmt) &
                   ,itab=iGLLtab, jtab=jGLLtab, dim_t=LocDimTanToEdge)
   ! LocDimTanToEdge = local dimension (1=xi, 2=eta) tangent to current edge                    
   ! SignTang enforces the counterclockwise convention for tangent vector
   ! (iGLLtab,jGLLtab) are the GLL indices of the nodes of boundary element BcElmt
   ! in the same order (counterclockwise) as they appear in bc%ibool(:,BcElmt)
    if (bc%edge(BcElmt)==edge_U .OR. bc%edge(BcElmt)==edge_L) then
      SignTang = -1d0
    else
      SignTang = 1d0
    endif
    do kGLL = 1,bc%ngnod
      DGlobDLoc = SE_Jacobian(grid,bc%elem(BcElmt),iGLLtab(kGLL),jGLLtab(kGLL))
      DxzDt = DGlobDLoc(:,LocDimTanToEdge)
      Jac1D = sqrt( DxzDt(1)**2 + DxzDt(2)**2 )
      Tang = SignTang*DxzDt/Jac1D
      BcNode = bc%ibool(kGLL,BcElmt)
                            ! assembled, outwards normal
      NORM(BcNode,:) = NORM(BcNode,:) + (/ Tang(2),-Tang(1) /)
      W(BcNode) = W(BcNode) + grid%wgll(kGLL)*Jac1D
    enddo
  enddo
  if (periodic) then
    NORM(1,:) = NORM(1,:) + NORM(bc%npoin,:)
    NORM(bc%npoin,:) = NORM(1,:)
    W(1) = W(1) + W(bc%npoin)
    W(bc%npoin) = W(1)
  endif
  do BcElmt = 1,bc%nelem ! fix interelement normals:
    BcNode = bc%ibool(1,BcElmt)
    NORM(BcNode,:) = NORM(BcNode,:) / sqrt( NORM(BcNode,1)**2 + NORM(BcNode,2)**2 )
    BcNode = bc%ibool(bc%ngnod,BcElmt)
    NORM(BcNode,:) = NORM(BcNode,:) / sqrt( NORM(BcNode,1)**2 + NORM(BcNode,2)**2 )
  enddo
  
  end subroutine BC_get_normal_and_weights


end module spec_grid
