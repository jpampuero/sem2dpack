! SEM2DPACK version 2.3.3 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://www.seismolab.caltech.edu
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
