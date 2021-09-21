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
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Source time function spline-interpolated from values in a file 
! SYNTAX : &STF_TAB file /
!
! ARG: file     [string] ['stf.tab'] ASCII file containing the source time function,
!               two columns: time and value. Time can be irregularly sampled and
!               must increase monotonically.
!
! NOTE   : assumes value(t<min(time))=value(min(time)) 
!          and value(t>max(time))=value(max(time))
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
  !NOTE: assumes 0 derivative at extreme points of spline
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

  double precision :: tbis
  integer :: N 

  N = size(stf%t)
  !NOTE: assumes value(t<t1) = value(t1) and value(t>tN) = value(tN) 
  tbis = max(t,stf%t(1))
  tbis = min(tbis,stf%t(N))
  call splint(stf%t,stf%v,stf%v2,N,tbis,fun)

end function STF_TAB_fun

end module stf_tabulated
