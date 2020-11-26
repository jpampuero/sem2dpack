module mxmlib
! Product of contiguously packed square matrices
! following Deville, Fischer and Mund (2002), p 388-389

  implicit none
  private

  public :: mxm

contains

  function mxm(A,B,n) result(C)

    integer, intent(in) :: n
    double precision, dimension(n,n), intent(in) :: A,B
    double precision, dimension(n,n) :: C

    select case (n)
      case (2); C = mxm2(A,B)
      case (3); C = mxm3(A,B)
      case (4); C = mxm4(A,B)
      case (5); C = mxm5(A,B)
      case (6); C = mxm6(A,B)
      case (7); C = mxm7(A,B)
      case (8); C = mxm8(A,B)
      case (9); C = mxm9(A,B)
      case (10); C = mxm10(A,B)
      case default; C = matmul(A,B)
    end select

  end function mxm

!-----------------------------------------------------------------------
  function mxm2(a,b) result(c)
    integer, parameter :: n=2
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j)
    enddo
    enddo
  end function mxm2

!-----------------------------------------------------------------------
  function mxm3(a,b) result(c)
    integer, parameter :: n=3
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j)
    enddo
    enddo
  end function mxm3

!-----------------------------------------------------------------------
  function mxm4(a,b) result(c)
    integer, parameter :: n=4
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j)
    enddo
    enddo
  end function mxm4

!-----------------------------------------------------------------------
  function mxm5(a,b) result(c)
    integer, parameter :: n=5
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j)
    enddo
    enddo
  end function mxm5

!-----------------------------------------------------------------------
  function mxm6(a,b) result(c)
    integer, parameter :: n=6
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j) &
             + a(i,6)*b(6,j)
    enddo
    enddo
  end function mxm6

!-----------------------------------------------------------------------
  function mxm7(a,b) result(c)
    integer, parameter :: n=7
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j) &
             + a(i,6)*b(6,j) &
             + a(i,7)*b(7,j)
    enddo
    enddo
  end function mxm7

!-----------------------------------------------------------------------
  function mxm8(a,b) result(c)
    integer, parameter :: n=8
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j) &
             + a(i,6)*b(6,j) &
             + a(i,7)*b(7,j) &
             + a(i,8)*b(8,j)
    enddo
    enddo
  end function mxm8

!-----------------------------------------------------------------------
  function mxm9(a,b) result(c)
    integer, parameter :: n=9
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j) &
             + a(i,6)*b(6,j) &
             + a(i,7)*b(7,j) &
             + a(i,8)*b(8,j) &
             + a(i,9)*b(9,j)
    enddo
    enddo
  end function mxm9

!-----------------------------------------------------------------------
  function mxm10(a,b) result(c)
    integer, parameter :: n=10
    double precision, dimension(n,n) :: a,b,c
    integer :: i,j
    do j=1,n
    do i=1,n
      c(i,j) = a(i,1)*b(1,j) &
             + a(i,2)*b(2,j) &
             + a(i,3)*b(3,j) &
             + a(i,4)*b(4,j) &
             + a(i,5)*b(5,j) &
             + a(i,6)*b(6,j) &
             + a(i,7)*b(7,j) &
             + a(i,8)*b(8,j) &
             + a(i,9)*b(9,j) &
             + a(i,10)*b(10,j)
    enddo
    enddo
  end function mxm10

end module mxmlib
