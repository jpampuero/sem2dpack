! SEM2DPACK version 2.2.12beta -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module mat_gen

! MAT_GEN: handles material properties

  use prop_mat

  ! List here all material modules:
  use mat_mass
  use mat_elastic
  use mat_kelvin_voigt
  use mat_damage

  implicit none
  private

  public :: MAT_read, MAT_init_prop, MAT_init_work

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MATERIAL
! PURPOSE: Define elastic material properties of a tagged domain
! SYNTAX : &MATERIAL tag, mode /
!          Followed by material data, with format depending on the mode (see
!          below)
!
! ARG: tag      [int] [none] Number identifying a mesh domain 
! ARG: mode     [char*5] ['ISOTR'] Type of material and/or spatial distribution.
!               The following modes are implemented and this is their data
!               format:
!
!               'ISOTR' homogeneous isotropic elastic
!               One line, dble(3):
!               density, P-wave-velocity, S-wave-velocity
!
!               'ANISO' homogeneous anisotropic
!               One line, dble(5):
!               density, c11, c13, c33, c44
!
!               'XXXXX' isotropic with any 2D distribution
!               Three $DIST_XXXXX blocks: 
!               density, P-velocity, S-velocity
!
! NOTE   : two MATERIAL blocks can share the same domain tag,
!          for instance to assign elastic and plastic material properties
!          to the same domain
!
! END INPUT BLOCK

! Read properties of a two-dimensional
! isotropic or anisotropic linear elastic element
!
subroutine MAT_read(mat,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type (matpro_input_type), pointer :: mat(:)

  integer :: i,numat,ntags,tag
  character(10)  :: mode

 !-----------------------------------------------------------------------

  NAMELIST / MATERIAL / tag,mode

 ! numat = number of MATERIAL input blocks 
 ! ntags = number of materials (domains) 
  numat = 0
  ntags = 0
  rewind(iin)
  do 
    read(iin,MATERIAL,END=50) 
    ntags = max(tag,ntags)
    numat = numat+1
  enddo
  50 if (numat==0)  call IO_abort('MATERIAL parameters not found')
  if (echo_input) write(iout,100) ntags,numat

  allocate(mat(ntags))

 ! Read properties for each input block
  rewind(iin) 
  do i = 1,numat
    tag = 0
    mode = 'ISOTR'
    read(iin,MATERIAL)
    if (echo_input) write(iout,200) tag,mode
    if (tag<=0) call IO_abort('MAT_read: tag must be positive')
   !NOTE: in the subroutines called below make sure "mat" is "inout"
   !      so multiple properties can be appended to the same material
    select case (mode)
      case('KELVIN_VOIGT'); call MAT_KV_read(mat(tag),iin)
      case('ISOTR'); call ISOTR_read(mat(tag),iin)
      case('ANISO'); call ANISO_read(mat(tag),iin)
      case('DAMAGE'); call MAT_DMG_read(mat(tag),iin)
      case default; call HETE_read(mat(tag),mode,iin)
    end select
  enddo

  return

!---- formats
  100   format(//,' M a t e r i a l   s e t s :   2 D  e l a s t i c i t y', &
         /1x,54('='),//5x, &
         'Number of material sets . . . . . . . . . = ',i5/5x, &
         'Number of material inputs . . . . . . . . = ',i5)
  200   format(/5x, &
         'Material number . . . . . . . . . . (tag) = ',i5/5x, &
         'Type  . . . . . . . . . . . . . . .(mode) = ',a)

end subroutine MAT_read

!=======================================================================
! 
subroutine MAT_init_prop(mat_elem,mat_input,grid)

  use spec_grid, only : sem_grid_type, SE_elem_coord
  use echo, only : echo_init, iout, fmt1, fmtok
  use stdio, only : IO_abort

  type(matpro_elem_type), pointer :: mat_elem(:)
  type(matpro_input_type), intent(inout), target :: mat_input(:)
  type(sem_grid_type), intent(in) :: grid

  integer :: e,tag

  if (echo_init) then
    write(iout,*) 
    write(iout,'(a)') ' M a t e r i a l   p r o p e r t i e s'
    write(iout,'(a)') ' ====================================='
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Translating input velocity model'
  endif

  allocate(mat_elem(grid%nelem))

  do e=1,grid%nelem
    tag = grid%tag(e)
    if ( tag > size(mat_input) .or. tag<1 ) &
      call IO_abort('ELAST_init: element tag does not correspond to a material number')
    mat_elem(e)%input => mat_input(tag)
    call MAT_init_elem_prop(mat_elem(e), SE_elem_coord(grid,e))
  enddo
  if (echo_init) write(iout,fmtok)


end subroutine MAT_init_prop

!-----------------------------------------------------------------------
subroutine MAT_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_MASS_init_elem_prop(elem,ecoord)

  if (MAT_isElastic(elem)) then
    call ELAST_init_elem_prop(elem,ecoord)

  elseif (MAT_isKelvinVoigt(elem)) then
    call MAT_KV_init_elem_prop(elem,ecoord)

  elseif (MAT_isDamage(elem)) then
    call MAT_DMG_init_elem_prop(elem,ecoord)

  endif

end subroutine MAT_init_elem_prop

!=======================================================================
subroutine MAT_init_work(matwrk,matpro,matinp,grid,ndof,dt)

  use spec_grid, only : sem_grid_type

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elem_type), pointer :: matwrk(:)
  type(matpro_elem_type), intent(inout) :: matpro(grid%nelem)
  type(matpro_input_type), intent(inout) :: matinp(:)
  integer, intent(in) :: ndof
  double precision, intent(in) :: dt

  logical :: homog

  allocate(matwrk(grid%nelem))

  homog = size(matinp)==1

  call ELAST_init_work(matwrk,matpro,grid,ndof,homog)
  call MAT_KV_init_work(matwrk,matpro,grid%ngll,dt)
  call MAT_DMG_init_work(matwrk,matpro,grid,ndof)

end subroutine MAT_init_work

end module mat_gen
