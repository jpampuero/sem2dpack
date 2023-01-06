module stf_brune

! Brune (1970)'s assumes the seismic moment rate is 
!   Mo*wc^2*t*exp(-wc*t)
! This implies the seismic moment is 
!   Mo*( 1 - (1+wc*t)*exp(-wc*t) )

  implicit none
  private

  type STF_BRUNE_type
    private
    double precision :: ampli, fc
  end type STF_BRUNE_type

  public :: STF_BRUNE_type, STF_BRUNE_read, STF_BRUNE_fun

contains

!=====================================================================
!
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_BRUNE
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Brune (1970)'s model with omega-squared spectral fall-off:
!            stf(t) = ampli*( 1 - (1+2*pi*fc*t)*exp(-2*pi*fc*t) )
! SYNTAX : &STF_BRUNE ampli, fc /
!
! ARG: ampli    [dble] [1d0] Amplitude (usually the seismic moment)
! ARG: fc       [dble] [1d0] Corner frequency (Hz)
!
! END INPUT BLOCK

subroutine STF_BRUNE_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_BRUNE_type), intent(out) :: stf
  integer, intent(in) :: iin

  double precision :: ampli,fc
  NAMELIST / STF_BRUNE / ampli,fc

  ampli = 1d0
  fc = 1d0

  read(iin,STF_BRUNE,END=500)

  if (fc < 0d0) call IO_abort('STF_BRUNE_read: fc must be positive')

  stf%ampli = ampli
  stf%fc = fc

  if (echo_input) write(iout,200) ampli,fc
  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Brune',/5x, &
     'Amplitude. . . . . . . . . . . . . . . =',EN12.3,/5x, &
     'Corner frequency (Hz). . . . . . . . . =',EN12.3,/5x )
  
  500 call IO_abort('STF_BRUNE_read: input block STF_BRUNE not found')

end subroutine STF_BRUNE_read


!=====================================================================
function STF_BRUNE_fun(stf,t) result(fun)

  use constants, only : pi

  type(STF_BRUNE_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun, arg

  arg = 2*pi*stf%fc *max(t,0d0)
  fun = stf%ampli *( 1d0 - (1d0+arg)*exp(-arg) )

end function STF_BRUNE_fun

end module stf_brune
