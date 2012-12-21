module stf_gaussian

  implicit none
  private

  type STF_GAUSSIAN_type
    private
    double precision :: f0,ampli,onset
  end type STF_GAUSSIAN_type

  public :: STF_GAUSSIAN_type, STF_GAUSSIAN_read, STF_GAUSSIAN_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_GAUSSIAN
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Gaussian source time function  
!          stf(t) = ampli*exp[-(pi*f0*(t-onset))^2]
! SYNTAX : &STF_GAUSSIAN ampli, f0, onset /
!
! ARG: ampli    [dble] [1d0] Amplitude
! ARG: onset    [dble] [0d0] Delay time (s)
! ARG: f0       [dble] [1d0] Characteristic frequency bandwidth (Hz)  
!
! END INPUT BLOCK

subroutine STF_GAUSSIAN_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_GAUSSIAN_type), intent(out) :: stf
  integer, intent(in) :: iin

  real :: f0,onset,ampli
  
  NAMELIST / STF_GAUSSIAN / f0,onset,ampli

  f0=1d0
  onset = 0d0
  ampli = 1d0

  read(iin,STF_GAUSSIAN,END=500)

  if (f0 < 0d0) call IO_abort('STF_GAUSSIAN_read: f0 must be positive')

  stf%f0=f0
  stf%onset = onset
  stf%ampli = ampli

  if (echo_input) write(iout,200) f0,onset,ampli

  return

  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Gaussian',/5x, &
     'Characteristic frequency (Hz). . .(f0) =',EN12.3,/5x, &
     'Time delay (s) . . . . . . . . (onset) =',EN12.3,/5x, &
     'Amplitude. . . . . . . . . . . (ampli) =',EN12.3)
  
  500 call IO_abort('STF_GAUSSIAN_read: input block STF_GAUSSIAN not found')

end subroutine STF_GAUSSIAN_read


!=====================================================================
double precision function STF_GAUSSIAN_fun(stf,t)

  use constants, only: PI

  type(STF_GAUSSIAN_type), intent(in) :: stf
  double precision, intent(in) :: t

  double precision :: arg

  arg=PI*stf%f0*(t-stf%onset)
  arg=arg*arg
  STF_GAUSSIAN_fun = stf%ampli*exp(-arg)

end function STF_GAUSSIAN_fun

end module stf_gaussian
