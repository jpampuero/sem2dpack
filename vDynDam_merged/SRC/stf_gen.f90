module stf_gen
  
!! To add a new source time function:
!!   1. Create an stf_XXX.f90 module following the template stf_user.f90
!!   2. Modify the current module following the instructions in lines starting by "!!"
!!   3. Modify the file Makefile.depend
!!   4. Re-compile

  use stf_ricker
  use stf_tabulated
  use butterworth_filter
  use stf_harmonic
  use stf_brune
  use stf_gaussian
  use stf_user
  !! use stf_XXX
  !! Add here your new stf_XXX module
  
  implicit none
  private

  type stf_type
    private
    integer :: kind = 0
    type (ricker_type), pointer :: ricker => null()
    type (butter_type), pointer :: butter => null()
    type (stf_tab_type), pointer :: tab => null()
    type (stf_harmonic_type), pointer :: harmonic => null()
    type (stf_brune_type), pointer :: brune => null()
    type (stf_gaussian_type), pointer :: gaussian => null()
    type (stf_user_type), pointer :: user => null()
    !! Add here your new stf_XXX_type pointer
  end type stf_type

  integer, parameter :: IS_RICKER   = 1 &
                       ,IS_BUTTER   = 2 &
                       ,IS_TAB      = 3 &   
                       ,IS_HARMONIC = 4 &
                       ,IS_BRUNE    = 5 &
                       ,IS_GAUSSIAN = 6 &
                       ,IS_USER     = 7
  !! add here a unique tag number for your new source time function

  public :: stf_type, STF_read, STF_get

contains

!=====================================================================
!
  subroutine STF_read(stfname,stf,iin)

  use stdio, only: IO_abort

  character(15), intent(in) :: stfname
  type(stf_type), intent(inout) :: stf
  integer, intent(in) :: iin

    select case (stfname)
      case('RICKER')
        stf%kind = IS_RICKER
        allocate(stf%ricker)
        call RICKER_read(stf%ricker,iin)

      case('BUTTERWORTH')
        stf%kind = IS_BUTTER
        allocate(stf%butter)
        call BUTTER_read(stf%butter,iin)

      case('TAB')
        stf%kind = IS_TAB
        allocate(stf%tab)
        call STF_TAB_read(stf%tab,iin)

      case('HARMONIC')
        stf%kind = IS_HARMONIC
        allocate(stf%harmonic)
        call STF_HARMONIC_read(stf%harmonic,iin)

      case('BRUNE')
        stf%kind = IS_BRUNE
        allocate(stf%brune)
        call STF_BRUNE_read(stf%brune,iin)

      case('GAUSSIAN')
        stf%kind = IS_GAUSSIAN
        allocate(stf%gaussian)
        call STF_GAUSSIAN_read(stf%gaussian,iin)

      case('USER')
        stf%kind = IS_USER
        allocate(stf%user)
        call STF_USER_read(stf%user,iin)

      !! add here an input block for your source time function
      !! following the models above

      case default
        call IO_abort('STF_read: unknown source time function ')
    end select

  end subroutine STF_read

!=====================================================================
! compute amplitude of source time function
!
  double precision function STF_get(stf,t)

  use stdio, only: IO_abort

  type(stf_type), intent(in) :: stf
  double precision , intent(in) :: t

  select case (stf%kind)
    case(IS_RICKER); STF_get = ricker(t,stf%ricker)
    case(IS_BUTTER); STF_get = butter(stf%butter,t)
    case(IS_TAB);    STF_get = STF_TAB_fun(stf%tab,t)
    case(IS_HARMONIC);   STF_get = STF_HARMONIC_fun(stf%harmonic,t)
    case(IS_BRUNE);  STF_get = STF_BRUNE_fun(stf%brune,t)
    case(IS_GAUSSIAN); STF_get = STF_GAUSSIAN_fun(stf%gaussian,t)
    case(IS_USER);   STF_get = STF_USER_fun(stf%user,t)
    !! add here a call to your new source time function
    case default
      STF_get = 0d0
      call IO_abort('STF_get: unknown source time function ')
  end select

  end function STF_get


end module stf_gen
