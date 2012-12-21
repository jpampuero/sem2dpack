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
module prop_mat

  use distribution_cd, only : cd_type
  use prop_elem
  use stdio, only : IO_abort

  implicit none
  private

  ! MAX_NKIND must be > number of material types implemented
  integer, parameter :: MAX_NKIND = 10 
  integer, save :: MAT_nkind=0


  type prop_input_link_type
    character(10) :: name=''
    type(cd_type) :: data
    type(prop_input_link_type), pointer :: next => null()
  end type prop_input_link_type

  ! kind = flags the type of material (many can be ON simultaneously)
  ! list = linked list of material properties
  type matpro_input_type
    logical, dimension(MAX_NKIND) :: kind = .false.
    type(prop_input_link_type), pointer :: list => null()
  end type matpro_input_type


  type prop_elem_link_type
    character(10) :: name=''
    type(prop_elem_type) :: data
    type(prop_elem_link_type), pointer :: next=>null()
  end type prop_elem_link_type

  type matpro_elem_type
    integer :: ngll
    type(matpro_input_type), pointer :: input => null()
    type(prop_elem_link_type), pointer :: list => null()
  end type matpro_elem_type


 ! working arrays for each element
 !! Add here all the stuff that is needed often by the solver
 !! Usually these are material properties (stored or pointer here 
 !! for faster access than lookup in the linked list of material properties),
 !! products of material properties and quadrature weights,
 !! state variables, etc
  type matwrk_elem_type

   !-- elastic
    double precision, pointer :: a(:,:,:) => null()
    double precision, pointer :: H(:,:)=>null(), Ht(:,:)=>null()

   !-- kelvin-voigt
    double precision, pointer :: eta(:,:) => null()

   !-- plasticity
    double precision, pointer, dimension(:,:,:) :: ep
    double precision, pointer, dimension(:,:) :: &
      mu,lambda, weights, dxi_dx,dxi_dy,deta_dx,deta_dy
    double precision, pointer, dimension(:) :: e0,s0

   !-- damage
    double precision, pointer, dimension(:,:) :: alpha,alpha_dot
    double precision :: xi_0,gamma_r,beta,Cd,Cv
    
  end type matwrk_elem_type


  interface MAT_isKind
    module procedure MAT_isKind_input, MAT_isKind_elem
  end interface MAT_isKind

  interface MAT_isProp
    module procedure MAT_isProp_input, MAT_isProp_elem
  end interface MAT_isProp

  interface MAT_setProp
    module procedure MAT_setProp_input, MAT_setProp_elem_from_values &
                   , MAT_setProp_elem_from_value, MAT_setProp_elem_from_input
  end interface MAT_setProp

  interface MAT_getProp
    module procedure MAT_getProp_e, MAT_getProp_ij, MAT_getProp_11
  end interface MAT_getProp

  public :: matpro_input_type, matpro_elem_type, matwrk_elem_type  &
          , MAT_isKind, MAT_setKind &
          , MAT_isProp, MAT_getProp, MAT_setProp

contains

!=======================================================================
  logical function MAT_newKind()

  if (MAT_nkind == MAX_NKIND) call IO_abort('MAT_newKind: exceeded MAX_NKIND')
  MAT_nkind = MAT_nkind + 1
  MAT_newKind = MAT_nkind

  end function MAT_newKind

!=======================================================================
! Inquires if "input" kind is "isThis"
  logical function MAT_isKind_input(input,isThis) 

  type(matpro_input_type), intent(in) :: input
  integer, intent(inout) :: isThis

  if (isThis<1 .or. isThis>MAT_nkind) &
    call IO_abort('MAT_isKind: argument out of bounds')
  MAT_isKind_input = input%kind(isThis)

  end function MAT_isKind_input

!-----------------------------------------------------------------------
  logical function MAT_isKind_elem(elem,isThis) 

  type(matpro_elem_type), intent(in) :: elem
  integer, intent(inout) :: isThis

  if (isThis<1 .or. isThis>MAT_nkind) &
    call IO_abort('MAT_isKind: argument out of bounds')
  MAT_isKind_elem = MAT_isKind_input(elem%input,isThis)

  end function MAT_isKind_elem

!=======================================================================
! Sets "input" kind to "isThis"
! If "isThis" is yet undefined: defines a new kind tag
  subroutine MAT_setKind(input,isThis)

  type(matpro_input_type), intent(inout) :: input
  integer, intent(inout) :: isThis

  if (isThis<1) isThis = MAT_newKind()
  if (isThis>MAT_nkind) call IO_abort('MAT_setKind: argument out of bounds')
  input%kind(isThis) = .true.

  end subroutine MAT_setKind


!=======================================================================
!
! creates a new (empty) input property 
! and returns a handle (pointer) to it
  function MAT_newProp_input(input) result(L)

  type(matpro_input_type), intent(inout), target :: input

  type(prop_input_link_type), pointer :: L

  if (.not.associated(input%list)) then
  ! create list head
    allocate(input%list)
    L => input%list
  else
  ! add a link to the tail of the list
    L => input%list
    do while(associated(L%next))
      L => L%next
    enddo
    allocate(L%next) 
    L => L%next
  endif
  
  end function MAT_newProp_input

!=======================================================================
!
  logical function MAT_isProp_input(name,input)

  character(*), intent(in) :: name
  type(matpro_input_type), intent(in) :: input

  MAT_isProp_input = associated( MAT_getProp_input(name,input) )

  end function MAT_isProp_input

!=======================================================================
!
  function MAT_getProp_input(name,input) result(L)

  character(*), intent(in) :: name
  type(matpro_input_type), intent(in) :: input

  type(prop_input_link_type), pointer :: L

  L => input%list
  do while(associated(L))
    if (L%name == name) return
    L => L%next
  enddo

  end function MAT_getProp_input

!=======================================================================
!
  subroutine MAT_setProp_input(input,property_name,value,distribution_name,iin,txt)

  use echo, only : echo_input, iout
  use distribution_cd, only : DIST_CD_Read

  type(matpro_input_type), intent(inout) :: input
  character(*), intent(in) :: property_name
  double precision, optional :: value
  character(*), optional :: distribution_name
  integer, optional :: iin
  character(*), intent(out), optional :: txt

  type(prop_input_link_type), pointer :: L

  L => MAT_getProp_input(property_name,input)

  if (.not.associated(L)) then
    L => MAT_newProp_input(input)
  elseif (echo_input) then
    write(iout,*) 'WARNING: MAT_setProp: overwrite'
  endif

  L%name = property_name
  call DIST_CD_Read(L%data,value,distribution_name,iin,txt)

  end subroutine MAT_setProp_input

!=======================================================================
!
! creates a new (empty) element property 
! and returns a handle (pointer) to it
  function MAT_newProp_elem(elem) result(L)

  type(matpro_elem_type), intent(inout), target :: elem

  type(prop_elem_link_type), pointer :: L

  if (.not.associated(elem%list)) then
  ! create list head
    allocate(elem%list)
    L => elem%list
  else
  ! add a link to the tail of the list
    L => elem%list
    do while(associated(L%next))
      L => L%next
    enddo
    allocate(L%next) 
    L => L%next
  endif
  
  end function MAT_newProp_elem

!=======================================================================
!
  logical function MAT_isProp_elem(name,elem)

  character(*), intent(in) :: name
  type(matpro_elem_type), intent(in) :: elem

  MAT_isProp_elem = associated( MAT_getProp_elem(name,elem) )

  end function MAT_isProp_elem

!=======================================================================
!
  function MAT_getProp_elem(name,elem) result(L)

  character(*), intent(in) :: name
  type(matpro_elem_type), intent(in) :: elem

  type(prop_elem_link_type), pointer :: L

  L => elem%list
  do while(associated(L))
    if (L%name == name) return
    L => L%next
  enddo

  end function MAT_getProp_elem

!=======================================================================
!
  subroutine MAT_setProp_elem_from_input(elem,property_name,ecoord)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: ecoord(:,:,:)
   
  type(prop_elem_link_type), pointer :: L
  type(prop_input_link_type), pointer :: L_input

  L => MAT_getProp_elem(property_name,elem)

  if (.not.associated(L)) then
    L => MAT_newProp_elem(elem)
    L%name = property_name
  elseif (echo_init) then
    write(iout,*) 'WARNING: MAT_setProp: overwrite'
  endif

  if (.not. associated(elem%input)) call IO_abort('MAT_setProp_elem: undefined element input')
  L_input => MAT_getProp_input(property_name,elem%input)
  if (.not. associated(L_input)) call IO_abort('MAT_setProp_elem: undefined input property')
  L%data = PROP_set(L_input%data,ecoord)

  end subroutine MAT_setProp_elem_from_input

!-----------------------------------------------------------------------
!
  subroutine MAT_setProp_elem_from_values(elem,property_name,val)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: val(:,:)
   
  type(prop_elem_link_type), pointer :: L

  L => MAT_getProp_elem(property_name,elem)

  if (.not.associated(L)) then
    L => MAT_newProp_elem(elem)
    L%name = property_name
  elseif (echo_init) then
    write(iout,*) 'WARNING: MAT_setProp: overwrite'
  endif

  L%data = PROP_set(val)

  end subroutine MAT_setProp_elem_from_values

!-----------------------------------------------------------------------
!
  subroutine MAT_setProp_elem_from_value(elem,property_name,val)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: val
   
  type(prop_elem_link_type), pointer :: L

  L => MAT_getProp_elem(property_name,elem)

  if (.not.associated(L)) then
    L => MAT_newProp_elem(elem)
    L%name = property_name
  elseif (echo_init) then
    write(iout,*) 'WARNING: MAT_setProp: overwrite'
  endif

  L%data = PROP_set(val)

  end subroutine MAT_setProp_elem_from_value

!=======================================================================
!
  subroutine MAT_getProp_e(val,elem,property_name)

  type(matpro_elem_type), intent(in) :: elem
  character(*), intent(in) :: property_name

  double precision, intent(out) :: val(:,:)

  type(prop_elem_link_type), pointer :: L

  L => MAT_getProp_elem(property_name,elem)
  if (associated(L)) then
    val = PROP_get(L%data,size(val,1))
  else
    call IO_abort('MAT_getProp: undefined property')
  endif
  
  end subroutine MAT_getProp_e

!-----------------------------------------------------------------------
!
  subroutine MAT_getProp_ij(val,elem,property_name,i,j)

  type(matpro_elem_type), intent(in) :: elem
  character(*), intent(in) :: property_name
  integer, intent(in) :: i,j

  double precision, intent(out) :: val

  type(prop_elem_link_type), pointer :: L

  L => MAT_getProp_elem(property_name,elem)
  if (associated(L)) then
    val = PROP_get(L%data,i,j)
  else
    call IO_abort('MAT_getProp: undefined property')
  endif
  
  end subroutine MAT_getProp_ij

!-----------------------------------------------------------------------
!
  subroutine MAT_getProp_11(val,elem,property_name)

  type(matpro_elem_type), intent(in) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(out) :: val

  call MAT_getProp(val,elem,property_name,1,1)
  
  end subroutine MAT_getProp_11

end module prop_mat
