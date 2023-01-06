module stf_gaussxsin

  implicit none
  private

  type STF_GAUSSXSIN_type
    private
    double precision :: f0,ampli,onset,std
  end type STF_GAUSSXSIN_type

  public :: STF_GAUSSXSIN_type, STF_GAUSSXSIN_read, STF_GAUSSXSIN_fun

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_GAUSSXSIN
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: Gaussian times sinus source time function  
!          stf(t) = ampli*exp[-1/2*((t-onset)/std)^2]*sin(f0*t)
! SYNTAX : &STF_GAUSSXSIN ampli, f0, onset, std /
!
! ARG: ampli    [dble] [1d0] Amplitude
! ARG: onset    [dble] [0d0] Delay time (s)
! ARG: f0       [dble] [1d0] Characteristic frequency bandwidth (Hz)
! ARG: std      [dble] [1d0] Standard deviation 
!
! END INPUT BLOCK

subroutine STF_GAUSSXSIN_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_GAUSSXSIN_type), intent(out) :: stf
  integer, intent(in) :: iin

  real :: f0,onset,ampli,std
  
  NAMELIST / STF_GAUSSXSIN / f0,onset,ampli,std

  f0=1d0
  onset = 0d0
  ampli = 1d0
  std= 1d0

  read(iin,STF_GAUSSXSIN,END=500)

  if (f0 < 0d0) call IO_abort('STF_GAUSSXSIN_read: f0 must be positive')

  stf%f0=f0
  stf%onset = onset
  stf%ampli = ampli
  stf%std = std

  if (echo_input) write(iout,200) f0,onset,ampli,std

  return

  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = Gaussian x sinus',/5x, &
     'Characteristic frequency (Hz). . .(f0) =',EN12.3,/5x, &
     'Time delay (s) . . . . . . . . (onset) =',EN12.3,/5x, &
     'Amplitude. . . . . . . . . . . (ampli) =',EN12.3,/5x, &
     'Standard deviation . . . . . . . (std) =',EN12.3)
  
  500 call IO_abort('STF_GAUSSXSIN_read: input block STF_GAUSSXSIN not found')

end subroutine STF_GAUSSXSIN_read


!=====================================================================
double precision function STF_GAUSSXSIN_fun(stf,t)

  use constants, only: PI

  type(STF_GAUSSXSIN_type), intent(in) :: stf
  double precision, intent(in) :: t

  double precision :: arg

  arg=(t-stf%onset)/stf%std
  arg=real(1)/real(2)*arg*arg
  STF_GAUSSXSIN_fun = stf%ampli*exp(-arg)*sin(stf%f0*(t-stf%onset)) !-stf%onset

end function STF_GAUSSXSIN_fun

end module stf_gaussxsin
