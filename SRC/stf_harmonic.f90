module stf_harmonic

  implicit none
  private

  type STF_HARMONIC_type
    private
    double precision :: ampli, f0
  end type STF_HARMONIC_type

  public :: STF_HARMONIC_type, STF_HARMONIC_read, STF_HARMONIC_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_HARMONIC
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Harmonic source time function f(t) = ampli*sin(2*pi*t*f0)
! SYNTAX : &STF_HARMONIC ampli, f0 /
!
! ARG: ampli    [dble] [0d0] Amplitude
! ARG: f0       [dble] [0d0] Frequency
!
! END INPUT BLOCK

subroutine STF_HARMONIC_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_HARMONIC_type), intent(out) :: stf
  integer, intent(in) :: iin

  real :: ampli,f0
  NAMELIST / STF_HARMONIC / ampli,f0

  ampli = 0d0
  f0 = 0d0

  read(iin,STF_HARMONIC,END=500)

  if (f0 <= 0d0) call IO_abort('STF_HARMONIC_read: f0 must be positive')
  if (ampli == 0d0) call IO_abort('STF_HARMONIC_read: ampli must be non zero')

  stf%ampli = ampli
  stf%f0 = f0

  if (echo_input) write(iout,200) ampli,f0

  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Harmonic',/5x, &
     'Amplitude. . . . . . . . . . . . . . . =',EN12.3,/5x, &
     'Frequency. . . . . . . . . . . . . . . =',EN12.3)
  
  500 call IO_abort('STF_HARMONIC_read: input block STF_HARMONIC not found')

end subroutine STF_HARMONIC_read


!=====================================================================
function STF_HARMONIC_fun(stf,t) result(fun)

  use constants, only : PI

  type(STF_HARMONIC_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun

  fun = stf%ampli* sin(2d0*PI*t*stf%f0)

end function STF_HARMONIC_fun

end module stf_harmonic
