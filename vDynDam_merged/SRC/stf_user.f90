module stf_user

!! To add a new source time function:
!!   1. Modify the current template module following the instructions 
!!      in comment lines starting by "!!"
!!      or create a new stf_XXX.f90 module based on this template
!!   2. Modify the stf_gen.f90 module
!!   3. Modify the file Makefile.depend
!!   4. Re-compile

  implicit none
  private

  !! Define a data structure containing all the parameters
  !! needed to evaluate the source time function
  type STF_USER_type
    private
    double precision :: ampli, onset, par1, par2
    integer :: ipar1, ipar2
  end type STF_USER_type

  public :: STF_USER_type, STF_USER_read, STF_USER_fun

contains

!=====================================================================
!
!! Modify the documentation of the input block:
!
! BEGIN INPUT BLOCK
!
! NAME   : STF_USER
! GROUP  : SOURCE TIME FUNCTIONS
! PURPOSE: A template for user-supplied source time function.
!          File stf_user.f90 must be modified by the user to fit
!          special needs.
! SYNTAX : &STF_USER ampli, onset, par1, par2, ipar1, ipar2 /
!
! ARG: ampli    [dble] [1.] Amplitude
! ARG: onset    [dble] [0]  Delay time (secs)
! ARG: par1     [dble] [0]  Example parameter
! ARG: par1     [dble] [0]  Example parameter
! ARG: par1     [int] [0]  Example parameter
! ARG: par1     [int] [0]  Example parameter
!
!
! END INPUT BLOCK

subroutine STF_USER_read(stf,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  type(STF_USER_type), intent(out) :: stf
  integer, intent(in) :: iin

  !! Define type and place on a namelist the parameters you need to read from "Par.inp"
  real :: onset,ampli,par1,par2
  integer :: ipar1,ipar2
  NAMELIST / STF_USER / onset,ampli,par1,par2,ipar1,ipar2

  !! Set default values for the input parameters
  onset = 0d0
  ampli = 1d0
  par1 = 0d0
  par2 = 0d0
  ipar1 = 0
  ipar2 = 0

  read(iin,STF_USER,END=500)

  !! Control the input values
  if (ipar1 < 0) call IO_abort('STF_USER_read: ipar1 must be positive') ! an example

  !! Assign input values to data structure
  stf%onset = onset
  stf%ampli = ampli
  stf%par1  = par1
  stf%par2  = par2
  stf%ipar1 = ipar1
  stf%ipar2 = ipar2

  !! List here all input parameters
  !! and edit the message that will be written on screen during the input stage
  if (echo_input) write(iout,200) onset,ampli,par1,par2,ipar1,ipar2
  return
  200 format(5x, & 
     'Function Type. . . . . . . . . . . . . = User supplied',/5x, &
     'Time delay (s) . . . . . . . . . . . . =',EN12.3,/5x, &
     'Amplitude. . . . . . . . . . . . . . . =',EN12.3,/5x, &
     'Example parameter 1. . . . . . . . . . =',EN12.3,/5x, &
     'Example parameter 2. . . . . . . . . . =',EN12.3,/5x, &
     'Example parameter 3. . . . . . . . . . =',EN12.3,/5x, &
     'Example parameter 4. . . . . . . . . . =',EN12.3)
  
  500 call IO_abort('STF_USER_read: input block STF_USER not found')

end subroutine STF_USER_read


!=====================================================================
!! This function is called once every timestep, for each source
function STF_USER_fun(stf,t) result(fun)

  type(STF_USER_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun

!  integer :: seed = 12345       !! example of remanent variable
                                !! might be used for initialization
                                !! or for random number generator
  double precision :: arg

  !! evaluate your source time function 
  !! using the parameters in the structure stf
  arg = t-stf%onset
  fun = stf%ampli* sin(arg) + stf%par1 *arg**2 

  return
end function STF_USER_fun

end module stf_user
