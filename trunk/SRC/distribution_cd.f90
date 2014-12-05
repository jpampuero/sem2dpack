module distribution_cd
! Constant or Distribution structure

  use distribution_general

  implicit none
  private

  type cd_type
    private
    double precision :: c=0d0
    type(distribution_type), pointer :: d => null()
  end type cd_type

  interface DIST_CD_Init
    module procedure DIST_CD_Init_1, DIST_CD_Init_2
  end interface DIST_CD_Init

  public :: cd_type, DIST_CD_Read, DIST_CD_Init, &
            DIST_CD_isDist, DIST_CD_isNull

contains

!=====================================================================

  subroutine DIST_CD_Read(CD,C,D,iin,txt)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort
  use utils, only : setopt

  type(cd_type), intent(inout) :: CD
  double precision, intent(in), optional :: C
  character(*), optional :: D
  integer, intent(in), optional :: iin
  character(*), optional :: txt

  character(10) :: lD

  if (.not.present(C) .and. .not.present(D)) &
    call IO_abort('DIST_CD_Read: missing C or D')
  if (present(D) .and. .not.present(iin)) &
    call IO_abort('DIST_CD_Read: missing D or IIN')

  if (.not.DIST_CD_isNull(CD)) then
    if (echo_input) write(iout,*) 'WARNING: DIST_CD_Read: overwrite'
    call DIST_CD_destroy(CD)
  endif
    
  call setopt(lD,'',D)

  if (lD=='') then
    CD%c = C
    if (present(txt)) write(txt,'(EN13.3)') C
  else
    allocate(CD%d)
    call DIST_read(CD%d,lD,iin)
    if (present(txt)) txt=lD
  endif

  end subroutine DIST_CD_Read

!=====================================================================

  subroutine DIST_CD_Init_1(CD,Coord,F, keep_cd)

  type(cd_type), intent(inout) :: CD
  double precision, intent(in) :: Coord(:,:) ! [2,n]
  double precision, pointer :: F(:) ! [n]
  logical, optional :: keep_cd

  logical :: destroy

  if (.not.associated(F)) allocate( F(size(Coord,2)) )

  if (present(keep_cd)) then
    destroy = .not. keep_cd
  else
    destroy = .true.
  endif
  
  if (DIST_CD_isDist(CD)) then
    call DIST_generate(F,Coord,CD%d)
    if (destroy) call DIST_destructor(CD%d)
  else
    F = CD%c
  endif

  end subroutine DIST_CD_Init_1

!---------------------------------------------------------------------

  subroutine DIST_CD_Init_2(CD,Coord,F, keep_cd)

  type(cd_type), intent(inout) :: CD
  double precision, intent(in) :: Coord(:,:,:) ! [2,n,n]
  double precision, pointer :: F(:,:) ! [n,n]
  logical, optional :: keep_cd

  logical :: destroy

  if (.not.associated(F)) allocate( F(size(Coord,2),size(Coord,3)) )

  if (present(keep_cd)) then
    destroy = .not. keep_cd
  else
    destroy = .true.
  endif

  if (DIST_CD_isDist(CD)) then
    call DIST_generate(F,Coord,CD%d)
    if (destroy) call DIST_destructor(CD%d)
  else
    F = CD%c
  endif

  end subroutine DIST_CD_Init_2

!=====================================================================
  logical function DIST_CD_isDist(cd)
  type(cd_type), intent(in) :: cd
  DIST_CD_isDist = associated(cd%d)
  end function DIST_CD_isDist

!=====================================================================

  logical function DIST_CD_isNull(cd)
  type(cd_type), intent(in) :: cd
  DIST_CD_isNull = (.not.associated(cd%d) .and. cd%c==0d0) 
  end function DIST_CD_isNull

!=====================================================================
  subroutine DIST_CD_destroy(cd)
  type(cd_type), intent(inout) :: cd
  cd%c = 0d0
  if (associated(cd%d)) then
    call DIST_destructor(cd%d)
    deallocate(cd%d)
  endif
  end subroutine DIST_CD_destroy

end module distribution_cd
