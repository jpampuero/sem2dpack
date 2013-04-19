module stf_ricker

  use constants, only: PI

  implicit none
  private

  type ricker_type
    private
    double precision :: f0,t0,ampli
  end type ricker_type

  public :: RICKER_type ,&
            RICKER_read ,&
            RICKER ,&
            RICKER_deriv ,&
            RICKER_int

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_RICKER
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: The Ricker wavelet is the second derivative of a gaussian.
! SYNTAX : &STF_RICKER ampli, f0, onset /
!
! ARG: ampli    [real] [1.] Signed amplitude of the central peak
! ARG: f0       [real >0] [0] Fundamental frequency (Hz).
!                 distribution: it has a peak at f0 and an exponential
!                 decay at high frequency. The cut-off high frequency is usually
!                 taken as fmax = 2.5 x f0. 
! ARG: onset    [real >1/f0] [0] Delay time (secs) with respect to the peak value. 
!
! NOTE   : The spectrum has a peak at f0 and decays exponentially at high 
!          frequencies. Beyond 2.5*f0 there is little energy, this is a 
!          recommended value for fmax.
! NOTE   : onset>1/f0 is needed to avoid a strong jump at t=0, which can cause
!          numerical oscillations. Ignore if using incident waves.
!
! END INPUT BLOCK

  subroutine RICKER_read(src,iin)

  use stdio, only : IO_abort, IO_warning
  use echo , only : echo_input,iout

  type(ricker_type), intent(out) :: src
  integer, intent(in) :: iin

  real :: f0,onset,ampli

  NAMELIST / STF_RICKER / f0,onset,ampli

  f0    = 0.
  onset = 0.
  ampli = 1.

  read(iin,STF_RICKER,END=100)

  if (f0 <= 0.) call IO_abort('RICKER_read: f0 must be positive')
  if (onset < 1/f0) then
    !onset = 1./f0
    write(iout,*) '*** WARNING: RICKER_read: Onset time too small, ***'
!    write(iout,'(A,EN12.3)') '             will be reset to t0 = ',onset
    call IO_warning()
  endif

  src%f0    = dble(f0)
  src%t0    = dble(onset)
  src%ampli = dble(ampli)

  if (echo_input) write(iout,200) f0,onset,ampli

  return
  
  100 call IO_abort('RICKER_read: STF_RICKER input block not found')
  200 format(5x, & 
     'Source time function . . . . . . . . . = Ricker',/5x, &
     'Fundamental frequency (Hz) . . . .(f0) =',EN12.3,/5x, &
     'Time delay (s) . . . . . . . . (onset) =',EN12.3,/5x, &
     'Multiplying factor . . . . . . (ampli) =',EN12.3)
  
  end subroutine RICKER_read

!=====================================================================
  double precision function ricker(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker = - src%ampli * (1.d0-2.d0*arg)*exp(-arg)

  return
  end function ricker

!=====================================================================

  double precision function ricker_deriv(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker_deriv = (PI*src%f0)**2 *src%ampli *(t-src%t0)* (6.d0-4.d0*arg)*exp(-arg)

  return
  end function ricker_deriv
!
!=====================================================================

  double precision function ricker_int(t,src)

  type(ricker_type), intent(in) :: src
  double precision, intent(in) :: t

  double precision :: arg

  arg = PI*src%f0*(t-src%t0)
  arg = arg*arg
  ricker_int = - src%ampli * (t-src%t0)*exp(-arg)

  return
  end function ricker_int

end module stf_ricker
