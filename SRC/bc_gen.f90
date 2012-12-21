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

  public :: bc_type,bc_read,bc_set,bc_init_periodic,bc_init_not_periodic,bc_write

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
  
  integer :: i,nbc,ibc,tag,tags(2)
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
 
  if (echo_input) write(iout,'(//1x,A,/1x,64("="),/)') &
    "B o u n d a r y   C o n d i t i o n s   C o n t r o l   c a r d"

  !-- read the specific BC properties
  rewind(iunit)
  do i = 1,nbc
    if (echo_input) then
      if (bc(i)%tag(2) == 0) then
        write(iout,200) bc(i)%tag(1)
      else
        write(iout,201) bc(i)%tag
      endif
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

  200 format(/4x,'Boundary Tag  = ',I0)
  201 format(/4x,'Boundary Tags = ',I0,' and ',I0)

end subroutine bc_read


!=======================================================================
!
subroutine bc_init_periodic(bc,grid,rmass)

  use spec_grid, only : sem_grid_type

  type(bc_type), pointer :: bc(:)
  type(sem_grid_type), intent(inout) :: grid
  double precision, intent(inout) :: rmass(:)

  integer :: i

  if (.not. associated(bc)) return
  do i = 1,size(bc)
    if(bc(i)%kind=='PERIOD') call BC_PERIO_init(bc(i)%perio,bc(i)%tag,grid,rmass)
  enddo

end subroutine bc_init_periodic

!-----------------------------------------------------------------------
subroutine bc_init_not_periodic(bc,grid,elast,rmass,tim)

  use spec_grid, only : sem_grid_type
  use elastic, only : elast_type
  use time_evol, only: timescheme_type

  type(bc_type), pointer :: bc(:)
  type(sem_grid_type), intent(inout) :: grid
  type(elast_type), intent(in) :: elast
  type(timescheme_type), intent(in) :: tim
  double precision, intent(inout) :: rmass(:)

  integer :: i
  
  if (.not. associated(bc)) return

  do i = 1,size(bc)
    select case(bc(i)%kind)
      case('DT0TN0')
        call BC_DT0TN0_init(bc(i)%DT0TN0,bc(i)%tag(1),grid)
      case('DTTTN0')
        call BC_DTTTN0_init(bc(i)%DTTTN0,bc(i)%tag(1),grid)
      case('ABSORB')
        call BC_ABSO_init(bc(i)%abso,bc(i)%tag(1),grid,elast)
      case('LISFLT')
        call BC_LSF_init(bc(i)%lsf,bc(i)%tag,grid,rmass)
      case('SWFFLT')
        call BC_SWFF_init(bc(i)%swff,bc(i)%tag,grid,rmass,tim)
    end select
  enddo

end subroutine bc_init_not_periodic


!=======================================================================
!! Applies the boundary condition
subroutine bc_set(bc,which,field,field2,scal,src,grid,fields,time,assemble)

  use sources, only: source_type
  use fields_class, only: fields_type
  use spec_grid, only : sem_grid_type

  type(bc_type), pointer :: bc(:)
  character(*), intent(in) :: which
  double precision, dimension(:,:), intent(inout), optional :: field,field2
  double precision, dimension(:), intent(inout), optional :: scal
  type (sem_grid_type), optional, intent(in) :: grid
  type (fields_type), optional, intent(inout) :: fields
  type(source_type), optional, pointer :: src
  double precision, intent(in), optional :: time
  character(*), intent(in), optional :: assemble

  integer :: i

  if (.not. associated(bc)) return

  if (which == 'ALL') then
    do i = 1,size(bc) 
      call bc_set_single(bc(i))
    enddo
 
  else
    do i = 1,size(bc)
      if ( bc(i)%kind == which ) call bc_set_single(bc(i))
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
        call BC_ABSO_set(bc%abso,grid,fields,src,time)
      case('PERIOD')
        if (present(field)) call bc_perio_set(bc%perio,field,assemble)
        if (present(scal)) call bc_perio_set(bc%perio,scal,assemble)
      case('LISFLT')
        call BC_LSF_set(bc%lsf,field,field2)
      case('SWFFLT')
        call BC_SWFF_set(bc%swff,fields%accel,fields%veloc,fields%displ)
    end select

  end subroutine bc_set_single
  
end subroutine bc_set


!=======================================================================
! Writes data for faults, and eventually other BCs
subroutine BC_write(bc,fs,itime)

  use fields_class, only: fields_type

  type (fields_type), optional, intent(inout) :: fs
  type(bc_type), pointer :: bc(:)
  integer, intent(in) :: itime

  integer :: i

  if (.not. associated(bc)) return

  do i = 1,size(bc)
    if (bc(i)%kind == 'SWFFLT') call BC_SWFF_write(bc(i)%swff,fs%displ,fs%veloc,itime)
  enddo

end subroutine BC_write

end module bc_gen
