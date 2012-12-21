! SEM2DPACK version 2.2.3 -- A Spectral Element Method tool for 2D wave propagation
!                            and earthquake source dynamics
! 
! Copyright (C) 2003 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics
! ETH Hönggerberg (HPP)
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 1 633 2197 (office)
! +41 1 633 1065 (fax)
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
module echo

! ECHO: control of message flow

  implicit none

  integer, save :: iout,ItInfo,ItSnapshots,ItSnapshot1
  logical, save :: echo_input,echo_init,echo_check,echo_run
  character(50), save :: title = ' '
  character(*), parameter :: fmt1='(5X,A," ...")', fmtok='("... [OK]")'

contains

!===============================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : ECHO
! PURPOSE: Parameters controlling runtime output
! SYNTAX : &ECHO Verbose,ItInfo,ItSnapshots,ItSnapshot1,OutputFile /
!
! ARG: Title            [word] [none] Title of the simulation
! ARG: Verbose          [char(4)] ['1101'] Verbose flags for input,initialization,
!                        check and solver phases. Example: '0001' is verbose only
!                        during solver.
! ARG: ItInfo           [int] [100] Frequency (in number of timesteps) at which
!                        solver echoes some basic information.
! ARG: ItSnapshots      [int] [100] Frequency (in number of timesteps) at which
!                        snapshots are dumped (usually PostScript) 
! ARG: ItSnaphot1       [int] [0]   Time step at which first snapshot is dumped
! ARG: OutputFile       [name] ['NONE'] Name of a file for echo instead of screen
!
! END INPUT BLOCK

  subroutine ECHO_read(iin)

  use stdio, only : IO_new_unit, IO_abort

  character(4) :: verbose
  integer, intent(in) :: iin
  character(50) :: OutputFile

  NAMELIST / ECHO / title,verbose,ItInfo,ItSnapshots,ItSnapshot1,OutputFile

  title          = ''
  verbose        = '1101'
  ItInfo       = 100
  ItSnapshots          = 100
  ItSnapshot1     = 0
  OutputFile        = 'NONE'

  rewind(iin)
  read(iin,ECHO,END=99)
99 continue

  echo_input     = verbose(1:1)=='1'
  echo_init      = verbose(2:2)=='1'
  echo_check     = verbose(3:3)=='1'
  echo_run       = verbose(4:4)=='1'

  if (OutputFile == 'NONE') then
    iout = 6 !  standard output (screen)
  else
    call ECHO_banner('Program  S P E C F E M : start',6)
    iout = IO_new_unit()
    open (iout,file=trim(OutputFile)//'_sem2d.txt',status='replace')
  endif

  call ECHO_banner('Program  S P E C F E M : start',iout)

  if (echo_input) write(iout,200) echo_input,echo_init,echo_check,echo_run &
                                 ,ItSnapshots,ItSnapshot1,ItInfo

  return

  200 format(//1x,'E c h o   C o n t r o l   c a r d',/1x,34('='),//5x,&
  'Echo during input phase. . . . . . . . . . . (Verbose(1)) = ',L1/ 5x, &
  'Echo during init phase . . . . . . . . . . . (Verbose(2)) = ',L1/ 5x, &
  'Echo during check phase. . . . . . . . . . . (Verbose(3)) = ',L1/ 5x, &
  'Echo during run phase. . . . . . . . . . . . (Verbose(4)) = ',L1/ 5x, &
  'Display frequency  . . . . . . . . . . . . .(ItSnapshots) = ',I0/ 5x, &
  'First display. . . . . . . . . . . . . . . .(ItSnapshot1) = ',I0/ 5x, &
  'Basic info output frequency. . . . . . . . . . . (ItInfo) = ',I0/ 5x)

  end subroutine ECHO_read

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
