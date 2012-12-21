! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bc_gen

!-- Import boundary conditions modules
  use bc_abso
  use bc_periodic
  use bc_dirneu
  use bc_DTTTN0 
  use bc_lsf
  use bc_swff

  implicit none
  private
  
!-- Object containing all boundary conditions
  type bc_type
    private
    integer :: tag(2)
    integer :: kind
!--- List here all bc types
    type(bc_dirneu_type)  , pointer :: dirneu
    type(bc_DTTTN0_type)  , pointer :: DTTTN0
    type(bc_abso_type)    , pointer :: abso
    type(bc_periodic_type), pointer :: perio
    type(bc_lsf_type)     , pointer :: lsf 
    type(bc_swff_type)    , pointer :: swff
  end type bc_type

  integer, parameter :: IS_DIRNEU = 1, &
                        IS_DTTTN0 = 2, &
                        IS_ABSORB = 3, &
                        IS_PERIOD = 4, &
                        IS_LISFLT = 5, &
                        IS_SWFFLT = 6

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
!               'DIRNEU', 'ABSORB', 'PERIOD', 'LISFLT', 'SWFFLT'
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

  if (echo_input) write(iout,'(//1x,A,/1x,37("="),/)') &
    "B o u n d a r y   C o n d i t i o n s"

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
      if (echo_input) write(iout,200) bc(i)%tag(1)
    elseif ( any(tags /= 0) ) then
      bc(i)%tag    = tags
      if (echo_input) write(iout,201) bc(i)%tag
    else
      call IO_abort('bc_read: tag null or not set')
    endif
    if (kind == ' ') call IO_abort('bc_read: kind not set')
    if (echo_input) write(iout,202) kind
 
   !-- read the BC properties
    select case(kind)
      case('DIRNEU')
        bc(i)%kind = IS_DIRNEU
        allocate(bc(i)%dirneu)
        call BC_DIRNEU_read(bc(i)%dirneu,iunit)
      case('DTTTN0')
        bc(i)%kind = IS_DTTTN0
        allocate(bc(i)%DTTTN0)
        call BC_DTTTN0_read(bc(i)%DTTTN0,iunit)
      case('ABSORB')
        bc(i)%kind = IS_ABSORB
        allocate(bc(i)%abso)
        call BC_ABSO_read(bc(i)%abso,iunit)
      case('PERIOD')
        bc(i)%kind = IS_PERIOD
        allocate(bc(i)%perio)
        call BC_PERIO_read(bc(i)%perio,iunit)
      case('LISFLT')
        bc(i)%kind = IS_LISFLT
        allocate(bc(i)%lsf)
        call BC_LSF_read(bc(i)%lsf,iunit)
      case('SWFFLT')
        bc(i)%kind = IS_SWFFLT
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
subroutine bc_init(bc,grid,mat,M,tim,fields,src)

  use spec_grid, only : sem_grid_type, BC_inquire
  use bnd_grid, only : bnd_grid_type
  use prop_mat, only : matpro_elem_type
  use time_evol, only: timescheme_type
  use fields_class, only: fields_type
  use sources, only: source_type
  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  type(sem_grid_type), intent(inout) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
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
    if(bc(i)%kind==IS_PERIOD) then
      call BC_PERIO_init(bc(i)%perio,bc(i)%tag,grid,M)
      perio => bc(i)%perio
    endif
  enddo
  do i = 1,size(bc)
    select case(bc(i)%kind)
      case(IS_DIRNEU)
        call BC_DIRNEU_init(bc(i)%dirneu,bc(i)%tag(1),grid,perio)
      case(IS_DTTTN0)
        call BC_DTTTN0_init(bc(i)%DTTTN0,bc(i)%tag(1),grid)
      case(IS_ABSORB)
        call BC_ABSO_init(bc(i)%abso,bc(i)%tag(1),grid,mat,M,tim,src,perio)
      case(IS_LISFLT)
        call BC_LSF_init(bc(i)%lsf,bc(i)%tag,grid,M(:,1),perio)
      case(IS_SWFFLT)
        call BC_SWFF_init(bc(i)%swff,bc(i)%tag,grid,M(:,1),tim,fields%veloc,perio)
    end select
  enddo

 ! write initial data
  call BC_write(bc,fields,0)

end subroutine bc_init


!=======================================================================
!! Applies the boundary condition
subroutine bc_set(bc,time,fields,field)

  use sources, only: source_type
  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, intent(in) :: time
  type (fields_type), intent(inout) :: fields
  double precision, dimension(:,:), intent(inout) :: field

  integer :: i

  if (.not. associated(bc)) return

 ! apply first periodic, then absorbing, then the rest
 ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_PERIOD) call bc_set_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_ABSORB) call bc_set_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB) call bc_set_single(bc(i))
  enddo
    
contains

  subroutine bc_set_single(bc)

    type(bc_type), intent(in) :: bc

    select case(bc%kind)
      case(IS_DIRNEU)
        call bc_DIRNEU_set(bc%dirneu,field,time)
      case(IS_DTTTN0)
        call bc_DTTTN0_set(bc%DTTTN0,fields%veloc,time)
      case(IS_ABSORB)
        call BC_ABSO_set(bc%abso,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case(IS_PERIOD)
        call bc_perio_set(bc%perio,field)
      case(IS_LISFLT)
        call BC_LSF_set(bc%lsf,fields%accel,fields%displ_alpha)
      case(IS_SWFFLT)
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
    if (bc(i)%kind == IS_SWFFLT) then
      if (echo_init .and. itime==0) write(iout,fmt1,advance='no') 'Exporting initial boundary data'
      call BC_SWFF_write(bc(i)%swff,fs%displ,fs%veloc,itime)
      if (echo_init .and. itime==0) write(iout,fmtok)
    endif
  enddo

end subroutine BC_write

end module bc_gen
