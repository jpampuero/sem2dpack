! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                             with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics Group
! ETH Hönggerberg HPP O 13.1
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
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
module bc_gen

!-- Import specific boundary conditions modules
  use bc_abso
  use bc_periodic
  use bc_DT0TN0 
  use bc_DTTTN0 
  use bc_lsf
  use bc_swff

  implicit none
  private
  
!-- Object containing all boundary conditions
  type bc_type
    private
    integer :: tag(2)
    character(6) :: kind
!--- List here all bc types
    type(bc_DT0TN0_type)  , pointer :: DT0TN0
    type(bc_DTTTN0_type)  , pointer :: DTTTN0
    type(bc_abso_type)    , pointer :: abso
    type(bc_periodic_type), pointer :: perio
    type(bc_lsf_type)     , pointer :: lsf 
    type(bc_swff_type)     , pointer :: swff
  end type bc_type

  public :: bc_type,bc_read,bc_set,bc_init,bc_write

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DEF
! PURPOSE: Define a boundary condition
! SYNTAX : &BC_DEF tag, tags, kind /
!            followed eventually by &BC_XXXX blocks 
!
! ARG:  tag     [int] [none] A number assigned to the boundary. If you are
!               using SEM2D built-in structured mesher the conventions are:
!                       1       bottom
!                       2       right
!                       3       up
!                       4       left
!               If you are importing a mesh, you must use the tags assigned
!               to the boundaries during the mesh construction.
! ARG:  tags    [int(2)] [none] Two tags are needed for interfaces (split-node)
!               and for periodic boundaries.
! ARG:  kind    [char*6] [none] Type of boundary condition. The following are
!               implemented:
!               'DT0TN0', ' DTTTN0', 'ABSORB', 'PERIOD', 'LISFLT', 'SWFFLT'
!
! NOTE   : you must DEFINE FIRST ALL PERIODIC BOUNDARIES
!
! NOTE   : Some of the boundary conditions need additional data. See their
!          respective input blocks if any.
!
! END INPUT BLOCK

subroutine bc_read(bc,iunit)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  integer, intent(in) :: iunit
  
  integer :: i,nbc,tag,tags(2)
  character(6) :: kind

  NAMELIST / BC_DEF / tag,tags,kind

  !-- count the BCs
  rewind(iunit)
  nbc = 0
  do 
    read(iunit,BC_DEF,END=10) 
    nbc = nbc + 1
  enddo
  10 continue

  !-- leave if there is no BC
  if (nbc == 0) return

  !-- read the BCs definition
  allocate( bc(nbc) )
  rewind(iunit)
  do i = 1,nbc
    tag  = 0
    tags = 0
    kind = ' '
    read(iunit,BC_DEF)
    bc(i)%tag = 0
    if (tag /= 0) then
      bc(i)%tag(1) = tag
    elseif ( any(tags /= 0) ) then
      bc(i)%tag    = tags
    else
      call IO_abort('bc_read: tag null or not set')
    endif
    if (kind == ' ') call IO_abort('bc_read: kind not set')
    bc(i)%kind = kind
  enddo
 
  if (echo_input) write(iout,'(//1x,A,/1x,37("="),/)') &
    "B o u n d a r y   C o n d i t i o n s"

  !-- read the specific BC properties
  rewind(iunit)
  do i = 1,nbc
    if (echo_input) then
      if (bc(i)%tag(2) == 0) then
        write(iout,200) bc(i)%tag(1)
      else
        write(iout,201) bc(i)%tag
      endif
      write(iout,202) bc(i)%kind
    endif
    select case(bc(i)%kind)
      case('DT0TN0')
        allocate(bc(i)%DT0TN0)
        call BC_DT0TN0_read(bc(i)%DT0TN0,iunit)
      case('DTTTN0')
        allocate(bc(i)%DTTTN0)
        call BC_DTTTN0_read(bc(i)%DTTTN0,iunit)
      case('ABSORB')
        allocate(bc(i)%abso)
        call BC_ABSO_read(bc(i)%abso,iunit)
      case('PERIOD')
        allocate(bc(i)%perio)
        call BC_PERIO_read(bc(i)%perio,iunit)
      case('LISFLT')
        allocate(bc(i)%lsf)
        call BC_LSF_read(bc(i)%lsf,iunit)
      case('SWFFLT')
        allocate(bc(i)%swff)
        call BC_SWFF_read(bc(i)%swff,iunit)
      case default  
        call IO_abort('bc_read: unknown kind')
    end select
  enddo

  return

  200 format(/5x,'Boundary tag. . . . . . . . . . . . (tag) = ',I0)
  201 format(/5x,'Boundary tags . . . . . . . . . . .(tags) = ',I0,' and ',I0)
  202 format( 5x,'Boundary condition. . . . . . . . .(kind) = ',A)

end subroutine bc_read


!-----------------------------------------------------------------------
subroutine bc_init(bc,grid,elast,M,tim,fields,src)

  use spec_grid, only : sem_grid_type, BC_inquire
  use bnd_grid, only : bnd_grid_type
  use elastic, only : elast_type
  use time_evol, only: timescheme_type
  use fields_class, only: fields_type
  use sources, only: source_type
  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  type(sem_grid_type), intent(inout) :: grid
  type(elast_type), intent(in) :: elast
  type(timescheme_type), intent(in) :: tim
  double precision, intent(inout) :: M(:,:)
  type (fields_type), intent(inout) :: fields
  type (source_type), pointer :: src(:)

  type(bnd_grid_type), pointer :: bctp
  type(bc_periodic_type), pointer :: perio
  integer :: i,j
  
  if (.not. associated(bc)) return

 ! check if requested bc tags exist in the mesh
  do i = 1,size(bc)
  do j=1,2
  if (bc(i)%tag(j) /= 0) then
    call BC_inquire( grid%bounds, tag = bc(i)%tag(j), bc_topo_ptr = bctp)
    if (.not. associated(bctp)) call IO_abort('bc_init: a tag defined in a BC_DEF block is not present in the mesh')
  endif
  enddo
  enddo
  nullify(bctp)

 ! first initialize periodic boundaries
  nullify(perio)
  do i = 1,size(bc)
    if(bc(i)%kind=='PERIOD') then
      call BC_PERIO_init(bc(i)%perio,bc(i)%tag,grid,M)
      perio => bc(i)%perio
    endif
  enddo
  do i = 1,size(bc)
    select case(bc(i)%kind)
      case('DT0TN0')
        call BC_DT0TN0_init(bc(i)%DT0TN0,bc(i)%tag(1),grid)
      case('DTTTN0')
        call BC_DTTTN0_init(bc(i)%DTTTN0,bc(i)%tag(1),grid)
      case('ABSORB')
        call BC_ABSO_init(bc(i)%abso,bc(i)%tag(1),grid,elast,M,tim,src,perio)
      case('LISFLT')
        call BC_LSF_init(bc(i)%lsf,bc(i)%tag,grid,M(:,1),perio)
      case('SWFFLT')
        call BC_SWFF_init(bc(i)%swff,bc(i)%tag,grid,M(:,1),tim,fields%veloc,perio)
    end select
  enddo

 ! write initial data
  call BC_write(bc,fields,0)

end subroutine bc_init


!=======================================================================
!! Applies the boundary condition
subroutine bc_set(bc,time,fields,field,which,scal)

  use sources, only: source_type
  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  type (fields_type), intent(inout) :: fields
  double precision, intent(in) :: time
  character(*), intent(in), optional :: which
  double precision, dimension(:,:), intent(inout) :: field
  double precision, dimension(:), intent(inout), optional :: scal

  integer :: i

  if (.not. associated(bc)) return

  if ( present(which) ) then
    do i = 1,size(bc)
      if ( bc(i)%kind == which ) call bc_set_single(bc(i))
    enddo
 
  else
   ! apply first periodic, then absorbing, then the rest
   ! Sep 29 2006: to avoid conflict between ABSORB and DT0TN0
    do i = 1,size(bc)
      if ( bc(i)%kind == 'PERIOD') call bc_set_single(bc(i))
    enddo
    do i = 1,size(bc)
      if ( bc(i)%kind == 'ABSORB') call bc_set_single(bc(i))
    enddo
    do i = 1,size(bc)
      if ( bc(i)%kind /= 'PERIOD' .and. bc(i)%kind /='ABSORB') call bc_set_single(bc(i))
    enddo
    
  endif

contains

  subroutine bc_set_single(bc)

    type(bc_type), intent(in) :: bc

    select case(bc%kind)
      case('DT0TN0')
        call bc_DT0TN0_set(bc%DT0TN0,field)
      case('DTTTN0')
        call bc_DTTTN0_set(bc%DTTTN0,fields%veloc,time)
      case('ABSORB')
        call BC_ABSO_set(bc%abso,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case('PERIOD')
        call bc_perio_set(bc%perio,field)
        if (present(scal)) call bc_perio_set(bc%perio,scal)
      case('LISFLT')
        call BC_LSF_set(bc%lsf,fields%accel,fields%displ_alpha)
      case('SWFFLT')
        call BC_SWFF_set(bc%swff,fields%accel,fields%veloc,fields%displ)
    end select

  end subroutine bc_set_single
  
end subroutine bc_set


!=======================================================================
! Writes data for faults, and eventually other BCs
subroutine BC_write(bc,fs,itime)

  use fields_class, only: fields_type
  use echo, only : iout,echo_init, fmt1,fmtok

  type (fields_type), optional, intent(inout) :: fs
  type(bc_type), pointer :: bc(:)
  integer, intent(in) :: itime

  integer :: i

  if (.not. associated(bc)) return

  do i = 1,size(bc)
    if (bc(i)%kind == 'SWFFLT') then
      if (echo_init .and. itime==0) write(iout,fmt1,advance='no') 'Exporting initial boundary data'
      call BC_SWFF_write(bc(i)%swff,fs%displ,fs%veloc,itime)
      if (echo_init .and. itime==0) write(iout,fmtok)
    endif
  enddo

end subroutine BC_write

end module bc_gen
