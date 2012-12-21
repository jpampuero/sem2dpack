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
module spec_grid

!=======================================================================
!
!                       This module deals with quadrangle elements
!                       defined by 4 or 9 control nodes.
!                       The control nodes are defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         t         .
!                               .                   .
!                               8         9  s      6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                           Local coordinate system : s,t
!
!=======================================================================
 
  use stdio, only : IO_abort

  implicit none
  private

!-----------------------------------------------------------------------
!-- General topology for boundary condition :

!     nelem        = total number of bc elements              [input]
!     npoin        = total number of bc nodes (no redundant)
!     bulk_element = #bc_element --> #bulk_element            [input]
!     element_edge = #bc_element --> #edge in bulk element    [input]
!     bulk_node    = #bc_node    --> #bulk_node
!     ibool        = (#gll_node,#bc_element) --> #bc_node

  type bc_topo_type
    integer :: nelem=0,npoin=0,ngll=0,tag=0
    integer, pointer :: ibool(:,:)=>null()
    integer, dimension(:), pointer :: bulk_element=>null() &
                                     ,element_edge=>null() &
                                     ,bulk_node=>null()
  end type bc_topo_type

!-----------------------------------------------------------------------
!-- Q Spectral element grid type

!---- Finite element mesh (Q4 or Q9)
!     ndime        = physical dimension of the domain
!     npgeo        = total number of control nodes
!     coorg        = coordinates of the control nodes
!     nelem        = total number of elements
!     ngnod        = number of control nodes per element (4 or 9)
!     tag          = element -> domain tags
!     knods        = local to global numbering for the control nodes
!                    (ignod,element) -> control node index
!     flat         = flag for cartesian grid (allows simplifications)

!---- GLL spectral element mesh

!------ Topology
!     ngll         = number of Gauss-Lobato-Legendre nodes per unit segment
!     npoin        = total number of GLL nodes in the mesh
!     ibool        = local to global numbering for the GLL mesh
!                    (igll,jgll,element) -> bulk node index
!     coord        = coordinates of the GLL nodes

!------ Working data
!     shape        = control nodes to GLL nodes shape functions
!     xjaci        = jacobian inverted, at each node
!     weights      = volume integration weight, at each node
!     hprime       = lagrange derivatives on the unit Q
!     hTprime      = transpose of hprime
!     wgll         = GLL integration weights
!     xgll         = GLL points on the unit segement

!------ Checks
!     fmax         = Highest frequency to be resolved by the grid

! NOTE: XJACI is very large but usually it can be cleaned at the end of 
!       the initialization phase. It is needed for instance if energy is to be
!       computed.
!       After init-phase, the FEM mesh database is only used for plotting
!       purposes, so if no plots are needed you can save them on disk and 
!       clean it from memory.
!       COORD is used in the solver only for computation of the incident
!       wavefield.


  type sem_grid_type
    integer :: ngll,ndime,npgeo,ngnod,nelem,npoin
    double precision :: fmax
    double precision, pointer :: xjaci(:,:,:,:,:),weights(:,:,:) &
                                ,coorg(:,:),coord(:,:) &
                                ,hprime(:,:),hTprime(:,:) &
                                ,wgll(:),xgll(:),shape(:,:,:)
    integer, pointer :: ibool(:,:,:),knods(:,:),tag(:)
    type (bc_topo_type), pointer :: bounds(:)
    logical :: flat
  end type sem_grid_type

!-----------------------------------------------------------------------

  type interpol_type
    integer :: n
    double precision, pointer :: shape(:,:,:),flagrange(:,:),xi(:)
  end type interpol_type


  ! (i,j) GLL for element corner nodes
  !icorner(1) = 1    ; jcorner(1) = 1
  !icorner(2) = ngll ; jcorner(2) = 1
  !icorner(3) = ngll ; jcorner(3) = ngll
  !icorner(4) = 1    ; jcorner(4) = ngll

  !-- element edges tags Up, Down, Left, Right
  integer, parameter :: edge_D = 1
  integer, parameter :: edge_R = 2
  integer, parameter :: edge_U = 3
  integer, parameter :: edge_L = 4

  !-- domain side tags Up, Down, Left, Right
  integer, parameter :: side_D = 1
  integer, parameter :: side_R = 2
  integer, parameter :: side_U = 3
  integer, parameter :: side_L = 4

  !-- control nodes defining the element edges
  integer, dimension(4), parameter :: EdgeKnod1 = (/  1, 2, 3, 4 /) &
                                     ,EdgeKnod2   = (/  2, 3, 4, 1 /)
  
  !-- physical dimension of the domain
  integer, parameter :: ndime = 2

  public :: bc_topo_type ,&
         sem_grid_type ,&
         interpol_type ,&
         SE_init_numbering ,&
         SE_init_gll ,&
         SE_init_interpol ,&
         SE_init_coord ,&
         SE_find_nearest_node ,&
         SE_get_edge_nodes ,&
         SE_inquire ,&
         Q49_read ,&
         Q49_init ,&
         Q49_init_interpol ,&
         SE_BcTopoInit,&
         BC_inquire ,&
         BC_get_normal_and_weights ,&
         edge_D,edge_R,edge_U,edge_L, &
         side_D,side_R,side_U,side_L
 
contains


!=======================================================================
!
!-- generate the global numbering
!
  subroutine SE_init_numbering(se)

  use memory_info
  use echo, only : echo_init,iout,fmt1,fmtok
  use generic_list, only : Link_Ptr_Type,Link_Type,List_Type &
     ,LI_Init_List,LI_Add_To_Head,LI_Get_Head,LI_Get_Len &
     ,LI_Get_Next,LI_Associated,LI_Remove_Head,LI_Remove_Next,LI_Nullify
  use stdio, only : IO_new_unit

  type Edge_Data_Type
    integer :: bulk_element,element_edge,knods(2)
    integer, pointer :: ibool(:)
  end type Edge_Data_Type

  type Edge_Type
    type(Link_Type) :: Link
    type(Edge_Data_Type) :: Data
  end type Edge_Type

  type Edge_Ptr_Type
    type(Edge_Type), pointer :: P
  end type Edge_Ptr_Type

  type(sem_grid_type), intent(inout) :: se

  integer :: nelem,i,j,k,e,ee,npedge,npcorn,nedge,knods(2),ngll,iol,ounit

  type(List_Type)     :: Edges_List
  type(Link_Ptr_Type) :: Link,Link2,Link_pre
  type(Edge_Ptr_Type) :: Edge,Edge2

!-----------------------------------------------------------------------

  if (echo_init) write(iout,'(// " G L L   g r i d" / " ===============" /)' )

  nelem = se%nelem
  ngll  = se%ngll

! initialisation du tableau de numerotation globale
  allocate(se%ibool(ngll,ngll,nelem))
  call storearray('ibool',size(se%ibool),iinteg)
  se%ibool(:,:,:) = 0

  se%npoin  = 0
  npedge = 0
  npcorn = 0

!---- store the edges data in a linked list:
  if (echo_init) write(iout,fmt1,advance='no') 'Listing edges'
  call LI_Init_List(Edges_List)
  do e = nelem,1,-1
  do nedge   = 1,4
    allocate(Edge%P)
    Edge%P%Data%bulk_element = e
    Edge%P%Data%element_edge = nedge
    Edge%P%Data%knods(1)     = se%knods(EdgeKnod1(nedge),e)
    Edge%P%Data%knods(2)     = se%knods(EdgeKnod2(nedge),e)
    select case(nedge)
      case(edge_D) ; Edge%P%Data%ibool  => se%ibool(:,1,e)
      case(edge_R) ; Edge%P%Data%ibool  => se%ibool(ngll,:,e)
      case(edge_U) ; Edge%P%Data%ibool  => se%ibool(ngll:1:-1,ngll,e)
      case(edge_L) ; Edge%P%Data%ibool  => se%ibool(1,ngll:1:-1,e)
    end select
    Link = transfer(Edge,Link)
    call LI_Add_To_Head(Link,Edges_List)
  enddo
  enddo
  if (echo_init) write(iout,fmtok)

!---- start numbering
  if (echo_init) write(iout,fmt1,advance='no') 'Numbering interior points'
  Link = LI_Remove_Head(Edges_List)
  Edge = transfer(Link,Edge)
  ee = Edge%P%Data%bulk_element  ! current edge belongs to this element

  do e = 1,nelem

   !-- points inside an element are unique
    do j=2,ngll-1
    do i=2,ngll-1
      se%npoin = se%npoin + 1
      se%ibool(i,j,e) = se%npoin
    enddo
    enddo

   !-- points on element edges
   ! process current edge only if it belongs to current element
   do while (ee==e)

      ! numbering of the new edge
      npedge = npedge + ngll - 2 ! number of edge-no-corner nodes (just for info)
      do k = 1,ngll
        if (k==1 .or. k==ngll) then  ! skip corner if already processed
          if (Edge%P%Data%ibool(k) /=0 ) cycle
          npcorn = npcorn + 1  ! number of corner nodes (just for info)
        endif
        se%npoin  = se%npoin + 1
        Edge%P%Data%ibool(k)  = se%npoin  ! this ibool points to the global ibool
      enddo

      ! search the grid for a matching edge (matching control nodes),
      ! update its ibool and remove it from the list of edges 
      ! search also for matching vertices and fill their ibool

      knods = Edge%P%Data%knods
      call LI_Nullify(Link_pre)
      Link2 = LI_Get_Head(Edges_List)
      do while (LI_Associated(Link2))
        Edge2 = transfer(Link2,Edge2)
        if ( Edge2%P%Data%knods(1) == knods(2) &
        .AND.Edge2%P%Data%knods(2) == knods(1)  ) then ! matching edge
        ! the ordering of the control nodes is inverse for matching edges
          Edge2%P%Data%ibool = Edge%P%Data%ibool(ngll:1:-1)
          Link2 = LI_Remove_Next(Link_pre,Edges_List)
          deallocate(Edge2%P)
          if (LI_Associated(Link_pre)) then
            Link2 = LI_Get_Next(Link_pre)
          else
            Link2 = LI_Get_Head(Edges_List)
          endif
        else 
          ! matching first vertex
          if ( knods(1)==Edge2%P%Data%knods(1) ) Edge2%P%Data%ibool(1) = Edge%P%Data%ibool(1)
          Link_pre = Link2
          Link2 = LI_Get_Next(Link2)
        endif
      enddo

     ! pick next edge from list (and remove it from the list)
      Link = LI_Remove_Head(Edges_List)
      if (LI_Associated(Link)) then
        Edge = transfer(Link,Edge)
        ee = Edge%P%Data%bulk_element  ! current edge belongs to this element
      else
        ee = 0
      endif

    enddo

  enddo

  if (echo_init) then
    write(iout,fmtok)
    write(iout,*)
    write(iout,100) 'Nodes of the global mesh  . . . . = ',se%npoin
    write(iout,100) 'Interior points . . . . . . . . . = ',se%npoin-npedge-npcorn
    write(iout,100) 'Edge points (without corners) . . = ',npedge
    write(iout,100) 'Corner points . . . . . . . . . . = ',npcorn
    write(iout,*)
    write(iout,fmt1,advance='no') 'Saving node index table (ibool) in a binary file'
  endif

  inquire( IOLENGTH=iol ) se%ibool
  ounit = IO_new_unit()
  open(ounit,file='ibool_sem2d.dat',status='replace',access='direct',recl=iol)
  write(ounit,rec=1) se%ibool
  close(ounit)
  if (echo_init) write(iout,fmtok)

  return

100 format(5X,A,I0)
  
  end subroutine SE_init_numbering

  
!=======================================================================
!
  subroutine SE_init_gll(se)

  use gll, only : zwgljd,hdgll
  use memory_info

  type (sem_grid_type), intent(inout) :: se
  integer :: ngll,ip,i

  ngll = se%ngll

!----    set up coordinates of the Gauss-Lobatto-Legendre points
  allocate(se%xgll(ngll)); call storearray('xgll',ngll,idouble)
  allocate(se%wgll(ngll)); call storearray('wgll',ngll,idouble)
  call zwgljd(se%xgll,se%wgll,ngll,0.d0,0.d0)
!------ if nb of points is odd, the middle abscissa is exactly zero
  if(mod(ngll,2) /= 0) se%xgll((ngll-1)/2+1) = 0.d0

!-- compute hprime coefficients (derivatives of Lagrange polynomials)

  allocate(se%hprime(ngll,ngll))
  allocate(se%hTprime(ngll,ngll))
  call storearray('hprime',ngll*ngll,idouble)
  call storearray('hTprime',ngll*ngll,idouble)

  do ip=1,ngll
  do i=1,ngll
    se%hprime(ip,i)  = hdgll(ip-1,i-1,se%xgll,ngll)
    se%hTprime(i,ip) = se%hprime(ip,i)
  enddo
  enddo

  
  end subroutine SE_init_gll
  


!=======================================================================
!
!--- for interpolation
!
  subroutine SE_init_interpol(interp,grid,ipts)

  use memory_info
  use gll, only : hgll

  type (interpol_type), pointer :: interp
  type (sem_grid_type), intent(in) :: grid
  integer, intent(in) :: ipts

  integer :: i,k,ngll

  if (ipts == 0) return

  allocate(interp)
  interp%n = ipts
 
  !---- regular grid
  allocate( interp%xi(ipts) )
  interp%xi = 2.d0* (/ (i-1, i=1,ipts) /) /dble(ipts-1) - 1.d0

  !---- Lagrange interpolators
  
  allocate(interp%flagrange(grid%ngll,ipts))
  call storearray('interp%flagrange',size(interp%flagrange),idouble)

  do i=1,grid%ngll
  do k=1,ipts
    interp%flagrange(i,k) = hgll(i-1,interp%xi(k),grid%xgll,grid%ngll)
  enddo
  enddo

  call Q49_init_interpol(grid,interp)

  end subroutine SE_init_interpol


!=======================================================================
!
!  set the global nodal coordinates
!
  subroutine SE_init_coord(se)

  use echo, only : echo_check,echo_init,iout,fmt1,fmtok
  use memory_info
  use stdio, only : IO_new_unit

  type (sem_grid_type), intent(inout) :: se

  integer i,j,e,ounit,iol
  double precision :: dum_coorg(ndime,se%ngnod)

  allocate(se%coord(ndime,se%npoin))
  call storearray('coord',size(se%coord),idouble)

!---- Coordinates of the global points 
  if (echo_init) write(iout,fmt1,advance='no') 'Defining nodes coordinates'
  do e = 1,se%nelem
    dum_coorg = se%coorg(:,se%knods(:,e))
    do j = 1,se%ngll
    do i = 1,se%ngll
      se%coord(:,se%ibool(i,j,e)) = matmul( dum_coorg, se%shape(:,i,j) )
    enddo
    enddo
  enddo
  if (echo_init) write(iout,fmtok)

  if (echo_check) then

!----  Save the grid in a text file
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Saving the grid coordinates (coord) in a text file...'
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
  if (echo_init) write(iout,fmt1,advance='no') 'Saving the grid coordinates (coord) in a binary file...'
  inquire( IOLENGTH=iol ) se%coord
  iol=iol/2
  ounit = IO_new_unit()
  open(ounit,file='coord_sem2d.dat',status='replace',access='direct',recl=iol)
  write(ounit,rec=1) real(se%coord)
  close(ounit)
  if (echo_init) write(iout,fmtok)

  end subroutine SE_init_coord


!=====================================================================
!
  subroutine SE_find_nearest_node(coord_in,grid,iglob,igll,jgll &
                                 ,element,coord,inner,bound,distance)

  type(sem_grid_type), intent(in)  :: grid
  double precision   , intent(in)  :: coord_in(:)
  double precision, optional, intent(out) :: coord(:),distance
  integer, optional, intent(in)  :: bound
  integer, optional, intent(out) :: iglob,igll,jgll,element
  logical, optional  , intent(in)  :: inner

  double precision :: dmin,xs,zs,xp,zp,dist
  integer :: ip,i,j,e,ilo,jlo,ihi,jhi,iglob_loc,igll_loc,jgll_loc,element_loc
  logical :: search_inner

  if (present(bound)) &
    call IO_abort('SE_find_nearest_node: locate to boundary not implemented')

  ! coordonnees demandees
  xs = coord_in(1)
  zs = coord_in(2)

  ! eventuellement, on ne fait la recherche que sur l'interieur de l'element
  search_inner = .false.
  if(present(inner)) search_inner = inner

  if(search_inner) then
    ilo = 2
    jlo = 2
    ihi = grid%ngll-1
    jhi = grid%ngll-1
  else
    ilo = 1
    jlo = 1
    ihi = grid%ngll
    jhi = grid%ngll
  endif

  ! recherche du point de grille le plus proche
  dmin = huge(dmin)
  do e=1,grid%nelem
  do i=ilo,ihi
  do j=jlo,jhi

    ip = grid%ibool(i,j,e)
    xp = grid%coord(1,ip)
    zp = grid%coord(2,ip)
    dist = sqrt((xp-xs)**2 + (zp-zs)**2)

    if (dist < dmin) then
      dmin        = dist
      iglob_loc   = ip
      igll_loc    = i
      jgll_loc    = j
      element_loc = e
    endif

  enddo
  enddo
  enddo

  if (present(coord)) coord(:) = grid%coord(:,iglob_loc)
  if (present(igll)) igll = igll_loc
  if (present(jgll)) jgll = jgll_loc
  if (present(element)) element = element_loc
  if (present(iglob)) iglob = iglob_loc
  if (present(distance)) distance = dmin

  end subroutine SE_find_nearest_node


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
!--- Read Q-grid data from ASCII file
!
subroutine Q49_read(sem,iin)

  use echo, only : echo_init,iout
  use stdio, only : IO_new_unit

  type(sem_grid_type), intent(out) :: sem
  integer, intent(in) :: iin

  integer :: iin2
  character(50) :: file

  NAMELIST / QGRID / file
  
  file = 'HERE'
  rewind(iin)
  read(iin,QGRID,END=100)
  if (file == 'HERE') then
    iin2 = iin  
  else
    iin2 = IO_new_unit()
    open(iin2,file=file,status='old',action='read')
  endif

  call Q49_read_coorg(sem,iin2)

  call Q49_read_elements(sem,iin2)

  if (iin2 /= iin) close(iin2)
    
  if (echo_init) write(iout,200) sem%ndime,sem%npgeo,sem%flat &
                              ,sem%ngnod,sem%nelem

  return

  100 call IO_abort('Q49_read: QGRID data missing')

  200 format(//1x,'G r i d   P a r a m e t e r s   C o n t r o l   c a r d', &
  /1x,34('='),//5x,&
  'Number of space dimensions . . . . . . . . . (ndime) =',I0/5x, &
  'Number of spectral elements control nodes. . (npgeo) =',I0/5x, &
  'Flat topography or not. . . . . . . . . .(topoplane) =',L1/5x, &
  'Number of control nodes per element . . . . .(ngnod) =',I0/5x, &
  'Number of spectral elements . . . . . . . . .(nelem) =',I0)

end subroutine Q49_read


!=======================================================================
!
!  Read macroblocs nodal coordinates
!
subroutine Q49_read_coorg(grid,iin)

  use memory_info

  type(sem_grid_type), intent(inout) :: grid
  integer, intent(in) :: iin
  
  integer :: i,ip,npgeo,dim
  logical :: topoplane

  NAMELIST / COORG / dim,npgeo,topoplane

  dim         = 2
  npgeo       = -1
  topoplane   = .false. 

  rewind(iin)
  read(iin,COORG,END=100)

  if (ndime /= 2) call IO_abort('Input: only dim = 2 implemented')
  if (npgeo <= 0) call IO_abort('Parameter npgeo is null or missing')

  grid%npgeo = npgeo
  grid%flat  = topoplane

  allocate(grid%coorg(dim,npgeo)) 
  call storearray('grid%coorg',size(grid%coorg),idouble)

  do i = 1,npgeo
    read(iin,*) ip,grid%coorg(:,ip)
  enddo

  return

  100 call IO_abort('COORG parameters not found')
  
end subroutine Q49_read_coorg

!=======================================================================
!
!  Read elements topology and material set for
!  spectral elements bloc
! 
  subroutine Q49_read_elements(grid,iin)

  use memory_info 

  type(sem_grid_type), intent(inout) :: grid
  integer, intent(in) :: iin

  integer :: n,kmatoread,i,ngnod,nelem

  NAMELIST / ELEMENTS / ngnod,nelem

  ngnod = -1
  nelem = -1

  rewind(iin)
  read(iin,ELEMENTS,END=100)

  if (ngnod < 0) call IO_abort('Q49_read_elements: ngnod parameter missing')
  if (ngnod /= 4 .or. ngnod /= 9) &
    call IO_abort('Q49_read_elements: only Q4 and Q9 elements supported')
  grid%ngnod = ngnod

  if (nelem <= 0) call IO_abort('Q49_read_elements: nelem parameter null or missing')
  grid%nelem = nelem

  allocate(grid%knods(ngnod,nelem)); call storearray('knods',ngnod*nelem,iinteg)
  allocate(grid%tag(nelem)); call storearray('kmato',nelem,iinteg)

  do i = 1,nelem
    read(iin,*) n,grid%tag(n),grid%knods(:,n)
  enddo

  if (any(grid%tag == 0)) &
    call IO_abort('Q49_read_elements: material domains not completely set')

  return

100 call IO_abort('Q49_read_elements: ELEMENTS data not found')

  end subroutine Q49_read_elements



!=======================================================================
!
!-- set up the shape functions and their local derivatives
!-- and compute the jacobian matrix at the integration points
!
subroutine Q49_init(sem)

  use memory_info
  use echo, only : echo_init,iout,fmt1,fmtok

  type(sem_grid_type), intent(inout) :: sem

  double precision :: dershape(sem%ndime,sem%ngnod,sem%ngll,sem%ngll) &
                     ,coorg_e(sem%ndime,sem%ngnod),wgll2(sem%ngll,sem%ngll) &
                     ,dvolu(sem%ngll,sem%ngll)

  integer :: ngll,ngnod,nelem,ndime,l1,l2,e,i,j

!-----------------------------------------------------------------------

  ngnod = sem%ngnod
  nelem = sem%nelem
  ndime = sem%ndime
  ngll  = sem%ngll

!-----------------------------------------------------------------------
!-- set up the shape functions and their local derivatives

  if (echo_init) write(iout,fmt1,advance='no') 'Defining shape functions'
  allocate(sem%shape(ngnod,ngll,ngll))
  call storearray('shape',size(sem%shape),idouble)

!----    4-noded rectangular element
  if(ngnod  ==  4) then
    do l2 = 1,ngll
    do l1 = 1,ngll
      call Q4_getshape(sem%xgll(l1),sem%xgll(l2),sem%shape(:,l1,l2),dershape(:,:,l1,l2))
    enddo
    enddo

!--    9-noded rectangular element
  else if(ngnod  ==  9) then
    do l2 = 1,ngll
    do l1 = 1,ngll
      call Q9_getshape(sem%xgll(l1),sem%xgll(l2),sem%shape(:,l1,l2),dershape(:,:,l1,l2))
    enddo
    enddo

  else
    call IO_abort('Wrong number of control nodes')
  endif

  if (echo_init) write(iout,fmtok)

  
!-----------------------------------------------------------------------
!----    compute the jacobian matrix at the integration points
!----    compute the weighting for volume integrals

  if (echo_init) write(iout,fmt1,advance='no') 'Defining jacobian and weights'

  allocate(sem%xjaci(ndime,ndime,ngll,ngll,nelem))
  allocate(sem%weights(ngll,ngll,nelem))
  call storearray('xjaci',size(sem%xjaci),idouble)
  call storearray('weights',size(sem%weights),idouble)
  
 ! wgll2(i,j) = wgll(i) * wgll(j)
  do j = 1,ngll
    wgll2(:,j) = sem%wgll * sem%wgll(j)
  enddo

  do e = 1,nelem

    coorg_e = sem%coorg(:,sem%knods(:,e))
    do j = 1,ngll
    do i = 1,ngll
      call Q49_getjac(dvolu(i,j), sem%xjaci(:,:,i,j,e), coorg_e, dershape(:,:,i,j) )
    enddo
    enddo

    sem%weights(:,:,e) = wgll2 * dvolu(:,:)

  enddo

  if (echo_init) write(iout,fmtok)

end subroutine Q49_init


!=======================================================================
!
!---- set up the shape functions for the interpolated grid
!
subroutine Q49_init_interpol(sem,interp)

  use memory_info

  type(sem_grid_type), intent(in) :: sem
  type(interpol_type), intent(inout) :: interp
  
  integer :: l2,l1

  allocate(interp%shape(sem%ngnod,interp%n,interp%n))
  call storearray('interp%shape',size(interp%shape),idouble)

  do l2 = 1,interp%n
  do l1 = 1,interp%n
    if ( sem%ngnod == 4 ) then
      call Q4_getshape(interp%xi(l1),interp%xi(l2),shape = interp%shape(:,l1,l2))
    else
      call Q9_getshape(interp%xi(l1),interp%xi(l2),shape = interp%shape(:,l1,l2))
    endif  
  enddo
  enddo

end subroutine Q49_init_interpol


!=======================================================================
!
!! GET SHAPE FUNCTIONS AND/OR THEIR DERIVATIVES 
!! s,t = local coordinates

subroutine Q4_getshape(s,t,shape,dershape)

  double precision, intent(in) :: s,t
  double precision, intent(out), optional :: shape(4)
  double precision, intent(out), optional :: dershape(2,4)

  double precision :: sp,sm,tp,tm 
  double precision, parameter :: one=1d0,quart=0.25d0
  
  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one

  if (present(shape)) then

    shape(1) = quart * sm * tm
    shape(2) = - quart * sp * tm
    shape(3) = quart * sp * tp
    shape(4) = - quart * sm * tp

  endif

  if (present(dershape)) then

    dershape(1,1) = quart * tm
    dershape(1,2) = - quart * tm
    dershape(1,3) =  quart * tp
    dershape(1,4) = - quart * tp

    dershape(2,1) = quart * sm
    dershape(2,2) = - quart * sp
    dershape(2,3) =  quart * sp
    dershape(2,4) = - quart * sm

  endif

end subroutine Q4_getshape

!-----------------------------------------------------------------------

subroutine Q9_getshape(s,t,shape,dershape)

  double precision, intent(in) :: s,t
  double precision, intent(out), optional :: shape(9)
  double precision, intent(out), optional :: dershape(2,9)

  double precision :: sp,sm,tp,tm,s2,t2,ss,tt,st
  double precision, parameter :: two=2d0,one=1d0,half=0.5d0,quart=0.25d0

  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one
  s2 = s * two
  t2 = t * two
  ss = s * s
  tt = t * t
  st = s * t

  if (present(shape)) then

    shape(1) = quart * sm * st * tm
    shape(2) = quart * sp * st * tm
    shape(3) = quart * sp * st * tp
    shape(4) = quart * sm * st * tp

    shape(5) = half * tm * t * (one - ss)
    shape(6) = half * sp * s * (one - tt)
    shape(7) = half * tp * t * (one - ss)
    shape(8) = half * sm * s * (one - tt)

    shape(9) = (one - ss) * (one - tt)

  endif

  if (present(dershape)) then

    dershape(1,1) = quart * tm * t * (s2 - one)
    dershape(1,2) = quart * tm * t * (s2 + one)
    dershape(1,3) = quart * tp * t * (s2 + one)
    dershape(1,4) = quart * tp * t * (s2 - one)

    dershape(2,1) = quart * sm * s * (t2 - one)
    dershape(2,2) = quart * sp * s * (t2 - one)
    dershape(2,3) = quart * sp * s * (t2 + one)
    dershape(2,4) = quart * sm * s * (t2 + one)

    dershape(1,5) = -one  * st * tm
    dershape(1,6) =  half * (one - tt) * (s2 + one)
    dershape(1,7) = -one  * st * tp
    dershape(1,8) =  half * (one - tt) * (s2 - one)

    dershape(2,5) =  half * (one - ss) * (t2 - one)
    dershape(2,6) = -one  * st * sp
    dershape(2,7) =  half * (one - ss) * (t2 + one)
    dershape(2,8) = -one  * st * sm

    dershape(1,9) = -one * s2 * (one - tt)
    dershape(2,9) = -one * t2 * (one - ss)

  endif

end subroutine Q9_getshape


!
!=======================================================================
! 
!! GET JACOBIAN (INVERTED) AND ELEMENTAL VOLUME
!! i,j = local index 
!! coorg   = ordered coordinates of control nods (2, 4 or 9)
!! dershape= derivatives of shape function (2, 4 or 9)
subroutine Q49_getjac(dvolu,xjaci,coorg,dershape)

  double precision, intent(in) :: coorg(:,:),dershape(:,:)
  double precision, intent(out) :: dvolu,xjaci(2,2)

  double precision :: xjac2(2,2)

  xjac2 = matmul(coorg, transpose(dershape) )
  call invert2(xjac2,xjaci,dvolu)
  if (dvolu <= 0d0) call IO_abort('Q49_getjac: Jacobian undefined')

end subroutine Q49_getjac

!----------------------------------------------------------------------
!-- inversion and determinant of a 2 x 2 matrix
subroutine invert2(A,B,det)

  double precision, intent(in) :: A(2,2)
  double precision, intent(out) :: B(2,2)
  double precision, optional, intent(out) :: det

  double precision :: determinant

  determinant = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if (determinant == 0d0) call IO_abort( 'invert2: undefined inverse matrix')

  B(1,1) =   A(2,2)
  B(2,1) = - A(2,1)
  B(1,2) = - A(1,2)
  B(2,2) =   A(1,1)

  B = B/determinant

  if (present(det)) det = determinant 

end subroutine invert2
!
!=======================================================================
! 
!! GET GLOBAL COORDINATES FROM LOCAL
!! coorg(ndime, 4 or 9)
!! shape_loc(4 or 9) = shape functions ALREADY evaluated at X_loc
subroutine Q49_loc2glob(X_glob,coorg,shape_loc)

  double precision, intent(in) :: coorg(:,:),shape_loc(:)
  double precision, intent(out) :: X_glob(:)

  X_glob = matmul(coorg,shape_loc)

end subroutine Q49_loc2glob

!=======================================================================
!
subroutine SE_BcTopoInit(grid)
  type(sem_grid_type), intent(inout) :: grid
  integer :: i
  if (.not.associated(grid%bounds)) return
  do i=1,size(grid%bounds)
!    write(51,'("BC ",I2)') i 
    call BC_set_bulk_node(grid%bounds(i),grid)
  enddo
end subroutine SE_BcTopoInit

!=======================================================================
!
subroutine BC_set_bulk_node(bc,grid)

  use generic_list, only : Link_Ptr_Type,Link_Type,List_Type &
     ,LI_Init_List,LI_Add_To_Head,LI_Get_Head &
     ,LI_Get_Next,LI_Associated,LI_Get_Len,LI_Remove_Head

  use utils, only: drank

  type(bc_topo_type), intent(inout) :: bc
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

  bc%ngll = grid%ngll
  allocate(bc%ibool(bc%ngll,bc%nelem))

!-----------------------------------------------------------------------
!  Build the database as a linked list

  bc_npoin = 0
  call LI_Init_List(BC_Nodes_List)
  call LI_Init_List(BC_Corners_List)

  do n=1,bc%nelem

!   get edge nodes, counterclockwise by default
    call SE_get_edge_nodes(grid,bc%bulk_element(n),bc%element_edge(n) &
                          ,bulk_node_vect)

    do kloc = 1,bc%ngll

      bulk_node = bulk_node_vect(kloc)
      new_node = .true.

!  When the node is an element corner (vertex) 
!  check if it is already in the BC list.
!  Note that the search is done in a BC_corners sublist.

      node_at_corner = kloc==1 .or. kloc==bc%ngll
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
        allocate(BC_node%P); allocate(BC_node%P%data)
        BC_node%P%data%bulk_node = bulk_node
        BC_node%P%data%bc_node = bc_inode
        Link = transfer(BC_node,Link)
        call LI_Add_To_Head(Link,BC_Nodes_List)
         
!  ... and eventually to the sublist of BC corners.
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

!-----------------------------------------------------------------------
!  Translate BC database from linked list to array storage

  bc%npoin = bc_npoin ! = LI_Get_Len(BC_Nodes_List)
  allocate(bc%bulk_node(bc%npoin)) 
  do i=1,bc%npoin
    Link = LI_Remove_Head(BC_Nodes_List)
    BC_node = transfer(Link,BC_node)
    bc%bulk_node(BC_node%P%Data%bc_node) = BC_node%P%Data%bulk_node
    deallocate(BC_node%P)
  enddo

!-----------------------------------------------------------------------
!  Sort BC nodes by increasing coordinate
!  Use the coordinate with the largest range
  allocate(BcCoord(2,bc%npoin))
  allocate(Isort(bc%npoin))
  BcCoord = grid%coord(:,bc%bulk_node)
  Lx = maxval(BcCoord(1,:)) - minval(BcCoord(1,:))
  Lz = maxval(BcCoord(2,:)) - minval(BcCoord(2,:))
  if (Lx>Lz) then
    call drank(BcCoord(1,:),Isort)
  else
    call drank(BcCoord(2,:),Isort)
  endif
!    write(51,'(I5,1X,I3)') ( bc%bulk_node(n),Isort(n),n=1,bc%npoin) 
  bc%bulk_node = bc%bulk_node(Isort)
!    write(51,'(I5)') bc%bulk_node
  allocate(Iback(bc%npoin))
  Iback(Isort) = (/ (i, i=1,bc%npoin) /) 
  do n=1,bc%nelem
    bc%ibool(:,n) = Iback(bc%ibool(:,n))
  enddo

  deallocate(BcCoord,Isort,Iback)

end subroutine BC_set_bulk_node

!=======================================================================
!
  subroutine SE_get_element_info(dmax,dmin,e,grid)

    type (sem_grid_type), intent(in) :: grid
    double precision, intent(out) :: dmax,dmin 
    integer, intent(in) :: e

    double precision :: x0,z0,x1,z1,x2,z2,x3,z3,rdist(4)
    integer :: ngll

    ngll = grid%ngll

    x0 = grid%coord(1,grid%ibool(1,1,e))
    z0 = grid%coord(2,grid%ibool(1,1,e))
    x1 = grid%coord(1,grid%ibool(ngll,1,e))
    z1 = grid%coord(2,grid%ibool(ngll,1,e))
    x2 = grid%coord(1,grid%ibool(ngll,ngll,e))
    z2 = grid%coord(2,grid%ibool(ngll,ngll,e))
    x3 = grid%coord(1,grid%ibool(1,ngll,e))
    z3 = grid%coord(2,grid%ibool(1,ngll,e))

    rdist(1) = (x1-x0)**2 + (z1-z0)**2 
    rdist(2) = (x2-x1)**2 + (z2-z1)**2 
    rdist(3) = (x3-x2)**2 + (z3-z2)**2 
    rdist(4) = (x3-x0)**2 + (z3-z0)**2

    dmax = sqrt( maxval(rdist) )
    dmin = sqrt( minval(rdist) )

  end subroutine SE_get_element_info


!=======================================================================
!
  subroutine SE_inquire(grid,element,igll,jgll,edge &
                       ,inv_jacobian,jacobian,itab,jtab,dim_t &
                       ,size_min,size_max)

  type (sem_grid_type), intent(in) :: grid
  integer, optional, intent(in) :: element,igll,jgll,edge
  double precision, dimension(grid%ndime,grid%ndime), optional, intent(out) :: &
    inv_jacobian,jacobian
  double precision, optional, intent(out) :: size_min,size_max
  integer, dimension(grid%ngll), optional, intent(out) :: itab,jtab
  integer, optional, intent(out) :: dim_t

  integer :: k
  double precision :: dmin,dmax

  !-- JACOBIAN and its inverse INV_JACOBIAN at a GLL node defined
  !   DGlobDLoc                DLocDGlob
  !   by its ELEMENT and local numbering IGLL,JGLL
  if (present(inv_jacobian)) inv_jacobian = grid%xjaci(:,:,igll,jgll,element)
  if (present(jacobian)) call invert2(grid%xjaci(:,:,igll,jgll,element) &
                                         ,jacobian)

  !-- give the GLL local numbering ITAB,JTAB for a given EDGE
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
    call SE_get_element_info(dmax,dmin,element,grid)
    if (present(size_min)) size_min = dmin
    if (present(size_max)) size_max = dmax
  endif
    
  end subroutine SE_inquire

!=======================================================================
!
  subroutine BC_inquire(bounds,tag,bc_topo_ptr)

  type (bc_topo_type), pointer :: bounds(:),bc_topo_ptr
  integer, intent(in) :: tag
  
  integer :: i
 
  do i=1,size(bounds)
    if ( bounds(i)%tag == tag ) then
      bc_topo_ptr => bounds(i)
      return
    endif
  enddo

  end subroutine BC_inquire

!=======================================================================
! Normal to a boundary, pointing out of the element
! Normals are assembled --> "average" normal between elements
  subroutine BC_get_normal_and_weights(bc,grid,NORM,W,periodic)

  type (bc_topo_type), intent(inout) :: bc
  type (sem_grid_type), intent(in) :: grid
  double precision, intent(out) :: NORM(bc%npoin,2)
  double precision, intent(out) :: W(bc%npoin)
  logical, intent(in) :: periodic

  double precision :: DGlobDLoc(2,2),DxzDt(2),Tang(2),Jac1D,SignTang
  integer :: iGLLtab(bc%ngll),jGLLtab(bc%ngll)
  integer :: BcElmt,kGLL,BcNode,LocDimTanToEdge

  NORM = 0d0
  W = 0d0
  do BcElmt = 1,bc%nelem
    call SE_inquire(grid, edge=bc%element_edge(BcElmt) &
                   ,itab=iGLLtab, jtab=jGLLtab, dim_t=LocDimTanToEdge)
   ! SignTang enforces the counterclockwise convention for tangent vector
    if (bc%element_edge(BcElmt)==edge_U .OR. bc%element_edge(BcElmt)==edge_L) then
      SignTang = -1d0
    else
      SignTang = 1d0
    endif
    do kGLL = 1,bc%ngll
      call SE_inquire(grid,element=bc%bulk_element(BcElmt) &
                     ,igll=iGLLtab(kGLL),jgll=jGLLtab(kGLL) &
                     ,jacobian=DGlobDLoc)
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
    BcNode = bc%ibool(bc%ngll,BcElmt)
    NORM(BcNode,:) = NORM(BcNode,:) / sqrt( NORM(BcNode,1)**2 + NORM(BcNode,2)**2 )
  enddo
  
  end subroutine BC_get_normal_and_weights


end module spec_grid
