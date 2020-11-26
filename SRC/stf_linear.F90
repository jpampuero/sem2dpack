module stf_linear

  implicit none
  private

  type STF_LINEAR_type
    private
    double precision :: intercept, slope
  end type STF_LINEAR_type

  public :: STF_LINEAR_type, STF_LINEAR_read, STF_LINEAR_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_LINEAR
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Linear source time function  
!          stf(t) = intercept + slope * t 
! SYNTAX : &STF_LINEAR intercept, slope /
!
! ARG: intercept    [dble] [0d0] 
! ARG: slope        [dble] [0d0]
!
! END INPUT BLOCK

subroutine STF_LINEAR_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_LINEAR_type), intent(out) :: stf
  integer, intent(in) :: iin

  double precision :: intercept, slope
  
  NAMELIST / STF_LINEAR / intercept,slope

  intercept = 0d0
  slope     = 0d0

  read(iin, STF_LINEAR, END=500)

  stf%intercept=intercept
  stf%slope = slope

  if (echo_input) write(iout,200) intercept, slope

  return

  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Linear',/5x, &
     'Intercept. . . . . . . . . (intercept) =',EN12.3,/5x, &
     'slope . . . . . . . . . . . .  (slope) =',EN12.3)
  
  500 call IO_abort('STF_LINEAR_read: input block STF_LINEAR not found')

end subroutine STF_LINEAR_read


!=====================================================================
double precision function STF_LINEAR_fun(stf,t)

  type(STF_LINEAR_type), intent(in) :: stf
  double precision, intent(in) :: t

  STF_LINEAR_fun = stf%intercept + stf%slope * t

end function STF_LINEAR_fun

end module stf_linear
