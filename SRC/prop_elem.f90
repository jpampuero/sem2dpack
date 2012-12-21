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
module prop_elem
! Material property of a spectral element
! Two possible databases:
!   1. homogeneous
!   2. heterogeneous

  implicit none
  private

! Structure "prop_elem_type" contains one (material) property 
! for one spectral element
! If constant property then storage = one float [homo], 
! else storage = ngll*ngll 
!
  type prop_elem_type
    private
    double precision :: homo = 0d0
    double precision, pointer :: hete(:,:) => null()
  end type prop_elem_type
  
  interface PROP_get
    module procedure PROP_get_1, PROP_get_2, PROP_get_ij_1, PROP_get_ij_2
  end interface PROP_get

  interface PROP_set
    module procedure PROP_set_cd1, PROP_set_cd2, PROP_set_val1, PROP_set_val2
  end interface PROP_set

  public :: prop_elem_type, PROP_get, PROP_set, PROP_isNull

contains

!=====================================================================

  function PROP_set_cd1(cd,coord) result(prop)

  use distribution_cd, only : cd_type, DIST_CD_Init, DIST_CD_isDist

  type(prop_elem_type) :: prop
  type(cd_type), intent(inout) :: cd
  double precision, intent(in) :: coord(:,:,:)

  double precision, pointer :: tmp(:,:)

  ! devel: could make DIST_CD_Init output argument of prop_elem_type
  ! call DIST_CD_Init(cd,coord,prop, keep_cd=.true.)

  call DIST_CD_Init(cd,coord,tmp, keep_cd=.true.)
  if (DIST_CD_isDist(cd)) then
    prop = PROP_set_val2(tmp)
  else
    prop = PROP_set_val1(tmp(1,1))
  endif
  deallocate(tmp)

  end function PROP_set_cd1

!------------------------------------------------

  function PROP_set_cd2(cd,coord) result(prop)

  use distribution_cd, only : cd_type

  type(cd_type), intent(inout) :: cd(:)
  double precision, intent(in) :: coord(:,:,:)

  type(prop_elem_type) :: prop(size(cd))

  integer :: k

  do k=1,size(cd)
    prop(k) = PROP_set_cd1(cd(k),coord)
  enddo

  end function PROP_set_cd2

!------------------------------------------------
  function PROP_set_val1(val) result(prop)

  type(prop_elem_type) :: prop
  double precision, intent(in) :: val

  prop%homo = val
  if (associated(prop%hete)) deallocate(prop%hete)

  end function PROP_set_val1

!------------------------------------------------
  function PROP_set_val2(val) result(prop)

  type(prop_elem_type) :: prop
  double precision, intent(in) :: val(:,:)

  prop%homo = 0d0
  if (associated(prop%hete)) deallocate(prop%hete)
  allocate(prop%hete(size(val,1),size(val,2)))
  prop%hete = val

  end function PROP_set_val2

!=====================================================================
  logical function PROP_isNull(prop)

  type(prop_elem_type), intent(in) :: prop

  PROP_isNull = ( prop%homo==0d0 .and. associated(prop%hete) )
  
  end function PROP_isNull

!=====================================================================

  function PROP_get_1(prop,n) result(eval)

  type(prop_elem_type), intent(in) :: prop
  integer, intent(in) :: n
  double precision :: eval(n,n)

  if (associated(prop%hete)) then
    eval = prop%hete
  else
    eval = prop%homo
  endif

  end function PROP_get_1

!------------------------------------------------

  function PROP_get_2(prop,n) result(eval)

  type(prop_elem_type) :: prop(:)
  integer, intent(in) :: n
  double precision :: eval(n,n,size(prop))

  integer :: k

  do k=1,size(prop)
    eval(:,:,k) = PROP_get_1(prop(k),n)
  enddo

  end function PROP_get_2

!------------------------------------------------
  function PROP_get_ij_1(prop,i,j) result(eval)

  type(prop_elem_type), intent(in) :: prop
  integer, intent(in) :: i,j
  double precision :: eval

  if (associated(prop%hete)) then
    eval = prop%hete(i,j)
  else
    eval = prop%homo
  endif

  end function PROP_get_ij_1

!------------------------------------------------

  function PROP_get_ij_2(prop,i,j) result(eval)

  type(prop_elem_type) :: prop(:)
  integer, intent(in) :: i,j
  double precision :: eval(size(prop))

  integer :: k

  do k=1,size(prop)
    eval(k) = PROP_get_ij_1(prop(k),i,j)
  enddo

  end function PROP_get_ij_2

end module prop_elem
