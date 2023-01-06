module src_force

! SRC_FORCE: collocated force source

  implicit none
  private

! dir(ndof) = direction vector, only used in P-SV
! iglob     = global node index 
  type so_force_type
    private
    double precision :: dir(2)
    integer :: iglob
  end type so_force_type

  public :: so_force_type,FORCE_read,FORCE_init,FORCE_add

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : SRC_FORCE
! GROUP  : SOURCE MECHANISM
! PURPOSE: Point force source
! SYNTAX : &SRC_FORCE angle /
!
! ARG: angle    [dble] [0d0]	For P-SV, the angle of the applied force, 
!                  in degrees, counterclockwise from Z-UP, e.g.: 
!                  90 points left, 180 points down
!                  For SH, angle is ignored and the SRC_FORCE block is not required.
!
! END INPUT BLOCK

  subroutine FORCE_read(so,iin)

  use stdio, only : IO_abort
  use echo, only : echo_input,iout
  use constants, only: PI

  type(so_force_type), intent(out) :: so
  integer, intent(in) :: iin

  double precision :: angle
  NAMELIST / SRC_FORCE / angle
  
  angle  = 0.d0
  read(iin,SRC_FORCE,END=100)
100 continue
  if (echo_input) write(iout,200) angle
  angle  = angle*PI/180.d0
  so%dir = (/ -sin(angle), cos(angle) /) ! counterclockwise / UP 

  return

!  100 call IO_abort('FORCE_read: SRC_FORCE input block not found')
  200 format( &
     5x, 'Source Type. . . . . . . . . . . . . . = Collocated Force', &
     /5x,'If P-SV: counterclockwise angle / up . = ',F0.2)

  end subroutine FORCE_read

!=====================================================================
!
  subroutine FORCE_init(so,iglob)

  type(so_force_type), intent(inout) :: so
  integer, intent(in) :: iglob

  so%iglob = iglob

  end subroutine

!=====================================================================
!
  subroutine FORCE_add(so,ampli,MxA)

  type(so_force_type), intent(inout) :: so
  double precision, intent(in) :: ampli
  double precision, intent(inout) :: MxA(:,:)

  if ( size(MxA,2)==1 ) then
    MxA(so%iglob,1) = MxA(so%iglob,1) + ampli
  else
    MxA(so%iglob,1) = MxA(so%iglob,1) + so%dir(1) *ampli
    MxA(so%iglob,2) = MxA(so%iglob,2) + so%dir(2) *ampli
  endif

  end subroutine FORCE_add


end module src_force
