! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! ETH H�nggerberg HPP O 13.1
! CH-8093 Z�rich
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
module stf_tabulated
! Tabulated source time function:
! read values from a file then spline-interpolate

  implicit none
  private

  type STF_TAB_type
    private
    double precision, pointer, dimension(:) :: t,v,v2
  end type STF_TAB_type

  public :: STF_TAB_type, STF_TAB_read, STF_TAB_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_TAB
! GROUP  : SRC_TIMEFUNCTION
! PURPOSE: Source time function spline-interpolated from values in a file 
! SYNTAX : &STF_TAB file /
!
! ARG: file     [string] ['stf.tab'] ASCII file containing the source time function,
!		 two columns per line:
!  		 (1) time
!		 (2) value
!
! NOTE   : time can be irregularly sampled
! NOTE   : assumes value=0 before min(time) and after max(time)
!
! END INPUT BLOCK

subroutine STF_TAB_read(stf,iin)

  use stdio, only: IO_abort, IO_new_unit, IO_file_length
  use echo , only : echo_input,iout
  use utils, only : spline

  type(STF_TAB_type), intent(out) :: stf
  integer, intent(in) :: iin

  integer :: N,i,iunit
  character(50) :: file
  NAMELIST / STF_TAB / file

  file = 'stf.tab'
  read(iin,STF_TAB,END=100)
100 continue

  N = IO_file_length(file)
  allocate( stf%t(N), stf%v(N) )
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  do i= 1,N
    read (iunit,*) stf%t(i), stf%v(i)
  end do
  close(iunit)

  allocate( stf%v2(N) )
  ! 0 derivative at extreme points of spline
  call spline(stf%t,stf%v,N,0d0,0d0,stf%v2)   

  if (echo_input) write(iout,200) trim(file)
  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = tabulated',/5x, &
     'File name . . . . . . . . . . . (file) =',A)
  
end subroutine STF_TAB_read


!=====================================================================
! could be made faster if assumes regular time sampling on input
function STF_TAB_fun(stf,t) result(fun)

  use utils, only : splint

  type(STF_TAB_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun

  integer :: N 

  N = size(stf%t)
  if (t>=stf%t(1) .and. t<=stf%t(N)) then
    call splint(stf%t,stf%v,stf%v2,N,t,fun)
  else
    fun = 0d0
  endif

end function STF_TAB_fun

end module stf_tabulated