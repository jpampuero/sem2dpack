module mesh_emc2

! MESH_EMC2: an interface for finite element mesh database 
! in EMC2's .ftq format

  use constants, only : NDIME
  use fem_grid, only : fem_grid_type

  implicit none
  private

  public :: EMC2_read,EMC2_build

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_EMC2
! GROUP  : MESH_DEF
! PURPOSE: Imports a mesh from INRIA's EMC2 mesh generator in FTQ format
! SYNTAX : &MESH_EMC2 file /
!
! ARG: file     [name] [none] Name of the FTQ file, including suffix
!
! END INPUT BLOCK

subroutine EMC2_read(mesh,iin)

  use memory_info
  use stdio, only : IO_new_unit,IO_abort
  use echo, only : echo_input, iout

  type(fem_grid_type), intent(out) :: mesh
  integer, intent(in) :: iin

  integer :: iftq,triangles,quadrangles,four,i,e
  character(50) :: file

  NAMELIST / MESH_EMC2 / file

  file = ' '
  read(iin,MESH_EMC2,END=100)
  if (file == ' ') call IO_abort('EMC2_read: you must provide a ftq file name')

  iftq = IO_new_unit()
  open(iftq,file=file,status='old')

  read(iftq,*) mesh%npoin, mesh%nelem, triangles,quadrangles
  if (triangles /= 0 .OR. quadrangles /= mesh%nelem) &
    call IO_abort('EMC2_read: the mesh contains triangles')
  mesh%ngnod = 4

  if (echo_input) write(iout,200) file, mesh%npoin, mesh%nelem

  allocate (mesh%knods(mesh%ngnod,mesh%nelem) )
  allocate (mesh%tag(mesh%nelem))
  call storearray('knods',size(mesh%knods),iinteg)
  call storearray('tag',size(mesh%tag),iinteg)
  do e=1,mesh%nelem
    read(iftq,*) four, mesh%knods(:,e),mesh%tag(e)
  enddo

  allocate( mesh%coord(NDIME,mesh%npoin) , mesh%ref(mesh%npoin) )
  call storearray('coorg',size(mesh%coord),idouble)
  do i=1,mesh%npoin
    read(iftq,*) mesh%coord(:,i),mesh%ref(i)
  enddo
 
  close(iftq)

  return

  100 call IO_abort('EMC2_read: MESH_EMC2 input block not found')

200 format(5x, &
      'File name . . . . . . . . . . . . . . . .(file) = ',a,/5x, &
      'Number of nodes . . . . . . . . . . . . . . . . = ',i5,/5x, &
      'Number of elements. . . . . . . . . . . . . . . = ',i5)

end subroutine EMC2_read



!=====================================================================
! EMC2_BUILD: Define the model boundaries

subroutine EMC2_build(mesh,grid)

  use fem_grid, only : fem_grid_type
  use generic_list, only : Link_Ptr_Type,Link_Type,List_Type &
     ,LI_Init_List,LI_Add_To_Head,LI_Get_Head,LI_Get_Len &
     ,LI_Get_Next,LI_Associated,LI_Remove_Head
  use stdio, only : IO_abort

 !-- For boundary tags list:

  type Tag_Type
    type(Link_Type) :: Link
    integer :: Tag
  end type Tag_Type

  type Tag_Ptr_Type
    type(Tag_Type), pointer :: P
  end type Tag_Ptr_Type

 !-- For boundary elements lists:

  type BCE_Data_Type
    integer :: elem,edge
  end type BCE_Data_Type

  type BCE_Type
    type(Link_Type) :: Link
    type(BCE_Data_Type) :: Data
  end type BCE_Type

  type BCE_Ptr_Type
    type(BCE_Type), pointer :: P
  end type BCE_Ptr_Type

 !-- 
  type(fem_grid_type), intent(in) :: mesh
  type(fem_grid_type), intent(inout) :: grid

  type(List_Type), allocatable :: BC_List(:)
  type(List_Type) :: Tags_List
  type(Link_Ptr_Type) :: Link
  type(Tag_Ptr_Type) :: Tag_Ptr,New_Tag
  type(BCE_Ptr_Type) :: bce
  integer :: ref(mesh%ngnod),tab123(3),iplus(4),i,tag &
            ,nbounds,e,nbc,bc_ref1,bc_ref2,i_null,bc_extr(3) &
            ,ibc,bc_nelem
  logical :: not_in_list

  grid%npoin = mesh%npoin
  grid%ngnod = mesh%ngnod
  grid%nelem = mesh%nelem
  grid%coord => mesh%coord
  grid%knods => mesh%knods
  grid%tag => mesh%tag

  tab123 = (/ 1,2,3 /)
  iplus  = (/ 2,3,4,1 /)

 !-- Get the list of used boundary tags as a linked list, 
 !   do not consider tag=0 , reserved for bulk nodes

  call LI_Init_List(Tags_List)
  do i=1,mesh%npoin

    tag = mesh%ref(i)
    if (tag==0) cycle  ! skip if bulk node (default tag=0)

    not_in_list = .true.
    Link = LI_Get_Head(Tags_List)
    do while(LI_Associated(Link))
      Tag_Ptr = transfer(Link,Tag_Ptr)
      if (Tag_Ptr%P%tag == tag) then
        not_in_list = .false.
        exit
      endif
      Link = LI_Get_Next(Link)
    enddo

    if (not_in_list) then
      allocate(New_Tag%P)
      New_Tag%P%tag = tag
      Link = transfer(New_Tag,Link)
      call LI_Add_To_Head(Link,Tags_List)
    endif

  enddo

  !-- transfer to an array.

  nbounds = LI_Get_Len(Tags_List)
  allocate( grid%bnds(nbounds) )
  do i = 1,nbounds
    Link = LI_Remove_head(Tags_List)
    Tag_Ptr = transfer(Link,Tag_Ptr)
    grid%bnds(i)%tag = Tag_Ptr%P%tag
    deallocate(Tag_Ptr%P)
  enddo

  !-- Get the list of boundary elements for each boundary
  !   as a linked list

  allocate(BC_List(nbounds))
  do i = 1,nbounds
    call LI_Init_List(BC_List(i))
  enddo

  do e=1,grid%nelem

    ref = mesh%ref( mesh%knods(:,e) )
    nbc = max( 0 , COUNT(ref /= 0) -1 )
    if (nbc == 0) cycle  ! skip if no bc in this element

    select case(nbc)

      case(1) ! a single boundary
        do i=1,mesh%ngnod 
          bc_ref1 = ref(i)
          bc_ref2 = ref(iplus(i))
          if (bc_ref1 == 0 .or. bc_ref2 == 0) cycle
! At the tip of a crack, we expect the refs to be different, at least on one
! side of the crack. Convention in EMC2 will be: ref=-1 for the tip nodes
          bc_ref1 = max(bc_ref1,bc_ref2)

          call add_bce(element=e,edge=i,tag=bc_ref1)

        enddo

      case(2) ! two boundaries at corner
        i_null = minloc(ref,1)
        bc_extr = i_null+ tab123
        where(bc_extr > mesh%ngnod) bc_extr = bc_extr-mesh%ngnod

        call add_bce(element=e,edge=bc_extr(1),tag=ref(bc_extr(1)) )
        call add_bce(element=e,edge=bc_extr(2),tag=ref(bc_extr(3)) )
   
      case default
        call IO_abort('EMC2_build: unexpected configuration')

    end select

  enddo

  ! transfer BCs to array storage

  do ibc =1,nbounds
    bc_nelem = LI_Get_Len(BC_List(ibc))
    grid%bnds(ibc)%nelem = bc_nelem
    if ( bc_nelem == 0 ) call IO_abort('EMC2_build: a BC is empty')
    allocate( grid%bnds(ibc)%elem(bc_nelem) )
    allocate( grid%bnds(ibc)%edge(bc_nelem) )
    do i=1,bc_nelem
      Link = LI_Remove_Head(BC_List(ibc))
      bce = transfer(Link,bce)
      grid%bnds(ibc)%elem(i) = bce%P%data%elem
      grid%bnds(ibc)%edge(i) = bce%P%data%edge
      deallocate(bce%P)
    enddo
  enddo

contains

subroutine add_bce(element,edge,tag)

  integer, intent(in) :: element,edge,tag

  type (BCE_Ptr_Type) :: new_bce
  type(Link_Ptr_Type) :: Link
  integer :: ibc

  ! get the corresponding BC_List
  do ibc =1,nbounds
    if ( grid%bnds(ibc)%tag == tag ) exit
    if ( ibc==nbounds ) call IO_abort('tag not referenced')
  enddo

  ! add to the corresponding BC_List
  allocate(new_bce%P)
  new_bce%P%data%elem = element
  new_bce%P%data%edge = edge
  Link = transfer(new_bce,Link)
  call LI_Add_To_Head(Link,BC_List(ibc))

end subroutine add_bce

end subroutine EMC2_build

end module mesh_emc2
