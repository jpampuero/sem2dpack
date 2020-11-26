module echo

! ECHO: control of message flow

  implicit none

  integer, parameter ::  iout = 6 !  standard output (screen)

  integer, save :: ItInfo
  logical, save :: echo_input,echo_init,echo_check,echo_run
  character(50), save :: title = ' '
  character(*), parameter :: fmt1='(5X,A," ...")', fmtok='("... [OK]")'

contains

  subroutine ECHO_set(verbose)

  character(4), intent(in) :: verbose

  echo_input     = verbose(1:1)=='1'
  echo_init      = verbose(2:2)=='1'
  echo_check     = verbose(3:3)=='1'
  echo_run       = verbose(4:4)=='1'

  end subroutine ECHO_set

!=====================================================================
! Prints a code banner.
! Gets date and time using f90 portable routines

  subroutine ECHO_banner(string1,iout)

  character*(*), intent(in) :: string1
  integer, intent(in) :: iout
  character(8)  :: datein
  character(10) :: timein
  character(16) :: dateprint
  character(8)  :: timeprint

  datein = ''
  timein = ''

  call date_and_time(datein,timein)

  dateprint = datein(7:8)//' - '//datein(5:6)//' - '//datein(1:4)
  timeprint = timein(1:2)//':'//timein(3:4)//':'//timein(5:6)

   write(iout,100) string1
   write(iout,101) title
   write(iout,102) dateprint,timeprint

  return

  100   format(//1x,79('-')/1x,79('-')/1x,a)
  101   format(1x,79('-')/1x,79('-')/1x,a50)
  102   format(1x,79('-')/,1x,79('-')/' D a t e : ',a16, &
         30x,' T i m e  : ',a8/1x,79('-'),/1x,79('-'))

  end subroutine ECHO_banner

end module echo
