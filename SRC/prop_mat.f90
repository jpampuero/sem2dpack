module prop_mat

  use distribution_cd, only : cd_type
  use prop_elem
  use stdio, only : IO_abort

  implicit none
  private

 !! MAX_NKIND must be > number of material types implemented
  integer, parameter :: MAX_NKIND = 10 
  integer, save :: MAT_nkind=0
  integer, parameter :: PROP_NAME_LEN = 8  ! max length of property names
                                            ! currently should be >= 6

  type prop_input_link_type
    character(PROP_NAME_LEN) :: name=''
    type(cd_type) :: data
    type(prop_input_link_type), pointer :: next => null()
  end type prop_input_link_type

  ! kind = flags a material attribute (material type or other)
  !        many flags can be ON simultaneously
  ! list = linked list of material properties
  type matpro_input_type
    private
    logical, dimension(MAX_NKIND) :: kind = .false.
    type(prop_input_link_type), pointer :: list => null()
  end type matpro_input_type

  type prop_elem_link_type
    character(PROP_NAME_LEN) :: name=''                !devel: PROP_NAME_LEN*1 bytes
    type(prop_elem_type) :: data                       !devel: (1 +1+6)*8 bytes
    type(prop_elem_link_type), pointer :: next=>null() !devel: 1*8 bytes
  end type prop_elem_link_type

  type matpro_elem_type
    private
    type(matpro_input_type), pointer :: input => null()
    type(prop_elem_link_type), pointer :: list => null()
  end type matpro_elem_type

  interface MAT_isKind
    module procedure MAT_isKind_input, MAT_isKind_elem
  end interface MAT_isKind

  interface MAT_setProp
    module procedure MAT_setProp_input, MAT_setProp_elem_from_values &
                   , MAT_setProp_elem_from_value, MAT_setProp_elem_from_input &
                   , MAT_setProp_elem_ptr_input
  end interface MAT_setProp

  interface MAT_getProp
    module procedure MAT_getProp_e, MAT_getProp_ij, MAT_getProp_11
  end interface MAT_getProp

  public :: matpro_input_type, matpro_elem_type &
          , MAT_isKind, MAT_setKind &
          , MAT_getProp, MAT_setProp &
          , MAT_isProp_input, MAT_isProp_elem

contains

!=======================================================================
  integer function MAT_newKind()

  if (MAT_nkind == MAX_NKIND) call IO_abort('MAT_newKind: exceeded MAX_NKIND')
  MAT_nkind = MAT_nkind + 1
  MAT_newKind = MAT_nkind

  end function MAT_newKind

!=======================================================================
! Inquires if "input" kind is "isThis"

  logical function MAT_isKind_input(input,isThis) 

  type(matpro_input_type), intent(in) :: input
  integer, intent(in) :: isThis

  if (isThis > 0 .and. isThis <= MAT_nkind) then
    MAT_isKind_input = input%kind(isThis)
  else
    MAT_isKind_input = .false.
  endif

  end function MAT_isKind_input

!-----------------------------------------------------------------------
  logical function MAT_isKind_elem(elem,isThis) 

  type(matpro_elem_type), intent(in) :: elem
  integer, intent(in) :: isThis

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
  logical function MAT_isProp_input(name,m)

  character(*), intent(in) :: name
  type(matpro_elem_type), intent(in) :: m

  MAT_isProp_input = associated( MAT_getProp_input(name,m%input) )

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

  subroutine  MAT_getProp_input_bis(name,input,L)

  character(*), intent(in) :: name
  type(matpro_input_type), intent(in) :: input
  type(prop_input_link_type), pointer :: L

  L => input%list
  do while(associated(L))
    if (L%name == name) return
    L => L%next
  enddo

  end subroutine MAT_getProp_input_bis
  
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

  if (len(property_name)>PROP_NAME_LEN) &
    call IO_abort('MAT_setProp_input: property_name is too long')

  L => MAT_getProp_input(property_name,input)

  if (.not.associated(L)) then
    L => MAT_newProp_input(input)
  elseif (echo_input) then
    write(iout,*) 'WARNING: MAT_setProp_input: overwrite'
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

  subroutine MAT_newProp_elem_bis(elem,L)

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
  
  end subroutine MAT_newProp_elem_bis

!=======================================================================
!
  logical function MAT_isProp_elem(name,elem)

  character(*), intent(in) :: name
  type(matpro_elem_type), intent(in) :: elem

  MAT_isProp_elem = associated( MAT_getProp_elem(name,elem) )

  end function MAT_isProp_elem

!=======================================================================
! Gets a property queried by its name
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

  subroutine MAT_getProp_elem_bis(name,elem,L)

  character(*), intent(in) :: name
  type(matpro_elem_type), intent(in) :: elem
  type(prop_elem_link_type), pointer :: L

  L => elem%list
  do while(associated(L))
    if (L%name == name) return
    L => L%next
  enddo

  end subroutine MAT_getProp_elem_bis

!=======================================================================
!
  subroutine MAT_setProp_elem_ptr_input(elem,input)

  type(matpro_elem_type), intent(inout) :: elem
  type(matpro_input_type), intent(in), target :: input

  elem%input => input

  end subroutine MAT_setProp_elem_ptr_input

!=======================================================================
!
  subroutine MAT_setProp_elem_from_input(elem,property_name,ecoord,memcount)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: ecoord(:,:,:)
  integer, intent(inout) :: memcount
   
  type(prop_elem_link_type), pointer :: L
  type(prop_input_link_type), pointer :: L_input

  if (len(property_name)>PROP_NAME_LEN) &
    call IO_abort('MAT_setProp_elem_from_input: property_name is too long')

  nullify(L)
  nullify(L_input)

!  L => MAT_getProp_elem(property_name,elem)
  call MAT_getProp_elem_bis(property_name,elem,L)

  if (.not.associated(L)) then
!    L => MAT_newProp_elem(elem)
    call MAT_newProp_elem_bis(elem,L)
    L%name = property_name
    memcount = memcount + size( transfer(L,(/0d0/)) ) 
  else 
    if (echo_init) write(iout,*) 'WARNING: MAT_setProp_elem_from_input: overwrite'
    memcount = memcount - PROP_size(L%data)
  endif

  if (.not. associated(elem%input)) call IO_abort('MAT_setProp_elem: undefined element input')
!  L_input => MAT_getProp_input(property_name,elem%input)
  call MAT_getProp_input_bis(property_name,elem%input,L_input)
  if (.not. associated(L_input)) call IO_abort('MAT_setProp_elem: undefined input property')
  L%data = PROP_set(L_input%data,ecoord)
  memcount = memcount + PROP_size(L%data)

  end subroutine MAT_setProp_elem_from_input

!-----------------------------------------------------------------------
!
  subroutine MAT_setProp_elem_from_values(elem,property_name,val,memcount)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: val(:,:)
  integer, intent(inout) :: memcount
   
  type(prop_elem_link_type), pointer :: L

  if (len(property_name)>PROP_NAME_LEN) &
    call IO_abort('MAT_setProp_elem_from_values: property_name is too long')

  L => MAT_getProp_elem(property_name,elem)

  if (.not.associated(L)) then
    L => MAT_newProp_elem(elem)
    L%name = property_name
    memcount = memcount + size( transfer(L,(/0d0/)) ) 
  else
    if (echo_init) write(iout,*) 'WARNING: MAT_setProp_elem_from_values: overwrite'
    memcount = memcount - PROP_size(L%data)
  endif

  L%data = PROP_set(val)
  memcount = memcount + PROP_size(L%data)

  end subroutine MAT_setProp_elem_from_values

!-----------------------------------------------------------------------
!
  subroutine MAT_setProp_elem_from_value(elem,property_name,val,memcount)

  use echo, only : echo_init, iout

  type(matpro_elem_type), intent(inout) :: elem
  character(*), intent(in) :: property_name
  double precision, intent(in) :: val
  integer, intent(inout) :: memcount
   
  type(prop_elem_link_type), pointer :: L

  if (len(property_name)>PROP_NAME_LEN) &
    call IO_abort('MAT_setProp_elem_from_value: property_name is too long')

  L => MAT_getProp_elem(property_name,elem)

  if (.not.associated(L)) then
    L => MAT_newProp_elem(elem)
    L%name = property_name
    memcount = memcount + size( transfer(L,(/0d0/)) ) 
  else
    if (echo_init) write(iout,*) 'WARNING: MAT_setProp_elem_from_value: overwrite'
    memcount = memcount - PROP_size(L%data)
  endif

  L%data = PROP_set(val)
  memcount = memcount + PROP_size(L%data)

  end subroutine MAT_setProp_elem_from_value

!=======================================================================
  subroutine MAT_getProp_e(val,elem,property_name)

! To DO: add optional ecoord, when prop has not been pre-evaluated
! we will generate it directly from input structure

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
