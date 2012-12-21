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
module mesh_emc2

! MESH_EMC2: an interface for finite element mesh database 
! in EMC2's .ftq format

  implicit none
  private

  type mesh_emc2_type
    private
    double precision, pointer :: coorg(:,:)
    integer, pointer :: knods(:,:),tag(:),ref(:)
    integer :: npgeo,nelem,ndime,ngnod
  end type mesh_emc2_type

  integer, parameter :: ndime = 2

  public :: mesh_emc2_type,EMC2_read,EMC2_build

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MESH_EMC2 [mesh]
! PURPOSE: Imports a mesh from INRIA's EMC2 mesh generator in FTQ format
! SYNTAX : &MESH_EMC2 file /
!
! ARG: file     [name] [none] Name of the FTQ file, including suffix
!
! END INPUT BLOCK

subroutine EMC2_read(mesh,iin)

  use stdio, only : IO_new_unit,IO_abort

  type(mesh_emc2_type), intent(out) :: mesh
  integer, intent(in) :: iin

  integer :: iftq,triangles,quadrangles,four,i,e
  character(50) :: file

  NAMELIST / MESH_EMC2 / file

  file = ' '
  read(iin,MESH_EMC2,END=100)
  if (file == ' ') call IO_abort('EMC2_read: you must provide a ftq file name')
  
  iftq = IO_new_unit()
  open(iftq,file=file,status='old')

  read(iftq,*) mesh%npgeo, mesh%nelem, triangles,quadrangles
  if (triangles /= 0 .OR. quadrangles /= mesh%nelem) &
    call IO_abort('EMC2_read: the mesh contains triangles')
  mesh%ngnod = 4

  mesh%ndime = ndime

  allocate (mesh%knods(mesh%ngnod,mesh%nelem) )
  allocate (mesh%tag(mesh%nelem))
  do e=1,mesh%nelem
    read(iftq,*) four, mesh%knods(:,e),mesh%tag(e)
  enddo

  allocate( mesh%coorg(mesh%ndime,mesh%npgeo) , mesh%ref(mesh%npgeo) )
  do i=1,mesh%npgeo
    read(iftq,*) mesh%coorg(:,i),mesh%ref(i)
  enddo
 
  close(iftq)

  return

  100 call IO_abort('EMC2_read: MESH_EMC2 inout block not found')

end subroutine EMC2_read



!=====================================================================
! EMC2_BUILD: Define the model boundaries

subroutine EMC2_build(mesh,grid)

  use spec_grid, only : sem_grid_type
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
    integer :: bulk_element,element_edge
  end type BCE_Data_Type

  type BCE_Type
    type(Link_Type) :: Link
    type(BCE_Data_Type) :: Data
  end type BCE_Type

  type BCE_Ptr_Type
    type(BCE_Type), pointer :: P
  end type BCE_Ptr_Type

 !-- 
  type(mesh_emc2_type), intent(in) :: mesh
  type(sem_grid_type), intent(inout) :: grid

  type(List_Type), pointer :: BC_List(:)
  type(List_Type) :: Tags_List
  type(Link_Ptr_Type) :: Link
  type(Tag_Ptr_Type) :: Tag_Ptr,New_Tag
  type(BCE_Ptr_Type) :: bce
  integer :: ref(mesh%ngnod),tab123(3),iplus(4),i,tag &
            ,nbounds,e,nbc,bc_ref1,bc_ref2,i_null,bc_extr(3) &
            ,ibc,bc_nelem
  logical :: not_in_list

  grid%npgeo =  mesh%npgeo
  grid%nelem =  mesh%nelem
  grid%ndime =  mesh%ndime
  grid%ngnod =  mesh%ngnod
  grid%coorg => mesh%coorg
  grid%knods => mesh%knods
  grid%tag   => mesh%tag

  tab123 = (/ 1,2,3 /)
  iplus  = (/ 2,3,4,1 /)

 !-- Get the list of used boundary tags as a linked list, 
 !   do not consider tag=0 , reserved for bulk nodes

  call LI_Init_List(Tags_List)
  do i=1,mesh%npgeo

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
  allocate( grid%bounds(nbounds) )
  do i = 1,nbounds
    Link = LI_Remove_head(Tags_List)
    Tag_Ptr = transfer(Link,Tag_Ptr)
    grid%bounds(i)%tag = Tag_Ptr%P%tag
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
    grid%bounds(ibc)%nelem = bc_nelem
    if ( bc_nelem == 0 ) call IO_abort('EMC2_build: a BC is empty')
    allocate( grid%bounds(ibc)%bulk_element(bc_nelem) )
    allocate( grid%bounds(ibc)%element_edge(bc_nelem) )
    do i=1,bc_nelem
      Link = LI_Remove_Head(BC_List(ibc))
      bce = transfer(Link,bce)
      grid%bounds(ibc)%bulk_element(i) = bce%P%data%bulk_element
      grid%bounds(ibc)%element_edge(i) = bce%P%data%element_edge
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
    if ( grid%bounds(ibc)%tag == tag ) exit
    if ( ibc==nbounds ) call IO_abort('tag not referenced')
  enddo

  ! add to the corresponding BC_List
  allocate(new_bce%P)
  new_bce%P%data%bulk_element = element
  new_bce%P%data%element_edge = edge
  Link = transfer(new_bce,Link)
  call LI_Add_To_Head(Link,BC_List(ibc))

end subroutine add_bce

end subroutine EMC2_build

end module mesh_emc2
