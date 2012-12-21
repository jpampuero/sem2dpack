! SEM2DPACK version 2.3.2 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-3429
! Fax  : (626) 564-0715
! 
! http://www.seismolab.caltech.edu
! 
! 
! This software is freely available for academic research purposes. 
! If you use this software in writing scientific papers include proper 
! attributions to its author, Jean-Paul Ampuero.
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
! 
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
