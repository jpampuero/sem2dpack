! SEM2DPACK version 2.2.12d -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module mat_gen

! MAT_GEN: handles material properties

!! To create a new material (for instance a material called "user"):
!!  1. Create a mat_user.f90 module, following existing examples
!!     It should contain the following:
!!     . saved parameter isUser
!!     . function MAT_isUser : tells if an element is of material User
!!     . subroutine MAT_USER_read : reads material properties from main input file
!!     . subroutine MAT_USER_init_elem_prop : initializes material properties for one element
!!     . subroutine MAT_USER_init_work : initializes working data for all elements
!!     . subroutine MAT_USER_action : actions to be performed during the solver phase
!!       for instance:
!!       subroutine MAT_USER_f : computes contribution to internal forces from one element
!!  2. Modify the prop_mat.f90 module:
!!     . set MAX_NKIND larger or equal than the number of existing material types
!!     . add new working data to matwrk_elem_type
!!  3. Modify the mat_gen.f90 module: follow comments that start with "!!" 
!!  4. Modify init.f90: if needed, pass additional arguments to MAT_init_work
!!  5. Modify solver.f90: call MAT_USER_action where appropriate
!!  6. Modify Makefile.depend: add new dependencies
!!  7. Re-compile (make)


  use prop_mat

!! List here all material modules:
  use mat_mass
  use mat_elastic
  use mat_kelvin_voigt
  use mat_damage
!!  use mat_user  

  implicit none
  private

  public :: MAT_read, MAT_init_prop, MAT_init_work

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MATERIAL
! PURPOSE: Define the material type of a tagged domain
! SYNTAX : &MATERIAL tag,isElastic,isKelvinVoigt,isDamage /
!          Followed by MAT_XXXX input locks
!
! ARG: tag      [int] [none] Number identifying a mesh domain 
! ARG: isElastic        [log] [F] Elastic material (see MAT_ELAST)
! ARG: isKelvinVoigt    [log] [F] Kelvin-Voigt material (see MAT_KV)
! ARG: isDamage [log] [F] Damage material (see MAT_DMG)
!
! NOTE   : Multiple material type can be assigned to a domain.
!          The MAT_XXXX blocks must then be given in a specific order.
!          The following multiple combinations are allowed:  
!           . elastic then Kelvin-Voigt
!
! END INPUT BLOCK

subroutine MAT_read(mat,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type (matpro_input_type), pointer :: mat(:)

  integer :: i,numat,ntags,tag
!! add a flag isUser here
  logical :: isElastic,isKelvinVoigt,isDamage !! ,iUser

  NAMELIST / MATERIAL / tag,isElastic,isKelvinVoigt,isDamage !! ,isUser

 ! numat = number of MATERIAL input blocks
 ! ntags = number of materials (domains)
  numat = 0
  ntags = 0
  rewind(iin)
  do
    tag = 0
    read(iin,MATERIAL,END=50)
    ntags = max(tag,ntags)
    numat = numat+1
    if (tag<=0) call IO_abort('MAT_read: tag must be positive')
  enddo
  if (numat /= ntags) call IO_abort('MAT_read: inconsistent or missing tags')
  50 if (numat==0)  call IO_abort('MATERIAL parameters not found')
  if (echo_input) write(iout,100) numat

  allocate(mat(numat))

 ! Read properties for each input block
  rewind(iin)
  do i = 1,numat

    tag = 0
    isElastic = .false.
    isKelvinVoigt = .false.
    isDamage = .false.
    !! isUser = .false.

    read(iin,MATERIAL)

    if (echo_input) write(iout,200) tag,isElastic,isKelvinVoigt,isDamage !!,isUser

   !NOTE: in the subroutines called below "mat" must be "inout"
   !      so multiple properties can be appended to the same material
   !WARNING: there is no rewind in these subroutines, the order of inputs matters !
    if (isElastic) call MAT_ELAST_read(mat(tag),iin)
    if (isKelvinVoigt) call MAT_KV_read(mat(tag),iin)
    if (isDamage) call MAT_DMG_read(mat(tag),iin)
    !! if (isUser) call MAT_USER_read(mat(tag),iin)

  enddo

  return

  100   format(//,' M a t e r i a l   P r o p e r t i e s',/1x,37('='),//5x, &
         'Number of materials . . . . . . . . . . . = ',I0)

  200   format(/5x, &
         'Material index. . . . . . . . . . . (tag) = ',I0/5x, &
         'Is elastic. . . . . . . . . . . . . . . . = ',L1/5x, &
         'Is Kelving-Voigt. . . . . . . . . . . . . = ',L1/5x, &
         'Is damage . . . . . . . . . . . . . . . . = ',L1)
      !! 'Is user . . . . . . . . . . . . . . . . . = ',L1)

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
    write(iout,fmt1,advance='no') 'Translating input model'
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

  if (MAT_isElastic(elem)) call MAT_ELAST_init_elem_prop(elem,ecoord)
  if (MAT_isKelvinVoigt(elem)) call MAT_KV_init_elem_prop(elem,ecoord)
  if (MAT_isDamage(elem)) call MAT_DMG_init_elem_prop(elem,ecoord)
  !!if (MAT_isUser(elem)) call MAT_USER_init_elem_prop(elem,ecoord)

end subroutine MAT_init_elem_prop

!=======================================================================
subroutine MAT_init_work(matwrk,matpro,grid,ndof,dt)

  use spec_grid, only : sem_grid_type

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elem_type), pointer :: matwrk(:)
  type(matpro_elem_type), intent(inout) :: matpro(grid%nelem)
  integer, intent(in) :: ndof
  double precision, intent(in) :: dt

  allocate(matwrk(grid%nelem))

  call MAT_ELAST_init_work(matwrk,matpro,grid,ndof)
  call MAT_KV_init_work(matwrk,matpro,grid%ngll,dt)
  call MAT_DMG_init_work(matwrk,matpro,grid,ndof)
  !!call MAT_USER_init_work(matwrk,matpro,grid, ... )
  !!if required, pass additional arguments to MAT_init_work 

end subroutine MAT_init_work

end module mat_gen
