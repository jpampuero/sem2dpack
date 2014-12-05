module prop_elem
! Material property of a spectral element
! Two possible databases:
!   1. homogeneous
!   2. heterogeneous

  implicit none
  private

  type hete_type
    double precision, pointer :: v(:,:) => null()
  end type hete_type

! "prop_elem_type" contains one material property for one spectral element
! If the property is constant [homo] then storage size = 1, 
! else [hete] storage size = ngll*ngll 
!
! NOTE: To minimize memory usage, "hete" is a "hete_type" pointer,
!       which requires only 8 bytes (4 actual + 4 alignment) when hete%v is not allocated.
!       If instead we make "hete" a pointer to a rank 2 double precision 
!       array it takes 6*8 bytes.
  type prop_elem_type
    private
    double precision :: homo = 0d0
    type(hete_type), pointer :: hete => null()
  end type prop_elem_type
  
  interface PROP_get
    module procedure PROP_get_1, PROP_get_2, PROP_get_ij_1, PROP_get_ij_2
  end interface PROP_get

  interface PROP_set
    module procedure PROP_set_cd1, PROP_set_cd2, PROP_set_val1, PROP_set_val2
  end interface PROP_set

  public :: prop_elem_type, PROP_get, PROP_set, PROP_isNull, PROP_size

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

  nullify(tmp)
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
  if (associated(prop%hete)) then
    deallocate(prop%hete%v)
    deallocate(prop%hete)
  endif

  end function PROP_set_val1

!------------------------------------------------
  function PROP_set_val2(val) result(prop)

  type(prop_elem_type) :: prop
  double precision, intent(in) :: val(:,:)

  prop%homo = 0d0
  if (associated(prop%hete)) then
    deallocate(prop%hete%v)
  else
    allocate(prop%hete)
  endif
  allocate(prop%hete%v(size(val,1),size(val,2)))
  prop%hete%v = val

  end function PROP_set_val2

!=====================================================================
  logical function PROP_isNull(prop)

  type(prop_elem_type), intent(in) :: prop

  PROP_isNull = ( prop%homo==0d0 .and. .not.associated(prop%hete) )
  
  end function PROP_isNull

!=====================================================================

  function PROP_get_1(prop,n) result(eval)

  type(prop_elem_type), intent(in) :: prop
  integer, intent(in) :: n
  double precision :: eval(n,n)

  if (associated(prop%hete)) then
    eval = prop%hete%v
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
    eval = prop%hete%v(i,j)
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

!=====================================================================
! memory size (in double precision words) of allocated internal storage
  integer function PROP_size(prop)

  type(prop_elem_type), intent(in) :: prop

! not counted in: size( transfer(prop,(/0d0/)) ) 

  if (associated(prop%hete)) then
    PROP_size = size( transfer(prop%hete,(/0d0/))) + size(prop%hete%v)
  else
    PROP_size = 0
  endif

  end function PROP_size

end module prop_elem
