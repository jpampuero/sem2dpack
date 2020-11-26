module stf_empty

  ! Empty source time function that gives 0
  !
  implicit none
  private

  type STF_empty_type
	! empty type, nothing
  end type STF_empty_type

  public :: STF_empty_type, STF_empty_read, STF_empty_fun

contains

subroutine STF_empty_read(stf,iin)

  use echo , only : echo_input,iout

  type(STF_empty_type), intent(out) :: stf
  integer, intent(in) :: iin
  
  ! do not read anything

  ! echo empty function defined 
  if (echo_input) write(iout, 200)
  return
  200 format(5x, &
     'Function Type. . . . . . . . . . . . . =  empty')

end subroutine STF_empty_read

!=====================================================================
!! This function is called once every timestep, for each source
function STF_empty_fun(stf, t) result(fun)

  ! return zero as empty source time function

  type(STF_empty_type), intent(in) :: stf
  double precision, intent(in) :: t
  double precision :: fun

  fun = 0.0d0 

  return
end function STF_empty_fun

end module stf_empty
