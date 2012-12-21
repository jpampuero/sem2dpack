! SEM2DPACK version 2.2.3 -- A Spectral Element Method tool for 2D wave propagation
!                            and earthquake source dynamics
! 
! Copyright (C) 2003 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics
! ETH Hönggerberg (HPP)
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 1 633 2197 (office)
! +41 1 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
! 
! This software is freely available for scientific research purposes. 
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
module distribution_linear
  
  use stdio, only: IO_abort, IO_new_unit
  implicit none
  private
 
  type linear_dist_type
    private
    double precision :: smoothing_length
    double precision, pointer :: X(:),val(:)
  end type linear_dist_type

  public :: linear_dist_type,read_linear_dist,generate_linear_dist ,&
            destroy_linear_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_LINEAR [distributions]
! PURPOSE: Piecewise linear 1D distribution along X.
! SYNTAX : &DIST_LINEAR file,length /
!
! ARG: file     [name] [none] Name of the ASCII file containing
!               the data to be interpolated in this format:
!                       n               (number of points to be read)
!                       x(1) v(1)       (coordinate and value)
!                       ...  ...
!                       x(n) v(n)
!               X must be sorted in increasing order
! ARG: length   [dble] [0]    Smoothing length (sliding average window)
!
! END INPUT BLOCK

  subroutine read_linear_dist (d, iin)

  type(linear_dist_type) :: d
  integer , intent(in) :: iin 

  double precision :: length
  integer :: iunit, i,N
  character(50) :: file

  NAMELIST / DIST_LINEAR / file,length

  length = 0d0
  file = ''
  read(iin,DIST_LINEAR, END = 100)

  d%smoothing_length = length

  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  read (iunit,*) N
  allocate( d%X(N),d%val(N) )
  do i= 1,N
    ! WARNING: X must be ordered
    read (iunit,*) d%x(i),d%val(i)
  end do
  close(iunit)

  return

  100 call IO_abort('read_linear_dist: DIST_LINEAR parameters missing')
!-- checkings:
  print *,"  Read ",file
  print *,"   X min max   ",minval(d%x),maxval(d%x)
  print *,"   VAL min max ",minval(d%val),maxval(d%val)

  end subroutine read_linear_dist

!
!***********************************************************************
!
  subroutine generate_linear_dist(field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(linear_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  integer :: unit,i 

    ! WARNING: interpolating along X
    if (d%smoothing_length > 0d0) then
      call smooth_interpol(d%X,d%val,coord(1,:),field,d%smoothing_length)
    else
      call interpol(d%X,d%val,coord(1,:),field)
    endif

    unit = IO_new_unit()
    open(unit,file='DistLinear_sem2d.tab',status='replace')
    do i= 1,size(field)
      write(unit,*) coord(1,i),field(i)
    enddo
    close(unit)

 
  end subroutine generate_linear_dist

!
!***********************************************************************
!
  subroutine interpol(xi,yi,xo,yo)

  use utils, only: hunt

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(out) :: yo(:)

  double precision :: eta
  integer :: i,ipos

  ipos = 1
  do i=1,size(xo)
    call hunt(xi,xo(i),ipos)
    eta = ( xo(i)-xi(ipos) ) / ( xi(ipos+1)-xi(ipos) )
    yo(i) = yi(ipos) + eta*( yi(ipos+1)-yi(ipos) ) 
  enddo

  end subroutine interpol

!
!***********************************************************************
! This is too complicated, better: smooth the underlying data at startup
! then just interpolate
!
  subroutine smooth_interpol(xi,yi,xo,yo,length)

  use utils, only: hunt

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(in) :: length
  double precision, intent(out) :: yo(:)

  double precision, allocatable :: sub_xi(:),sub_yi(:)
  double precision :: eta,sub_x0,sub_xN
  integer :: i,ipos0,iposN,sub_N

  ipos0 = 1
  iposN = 1

  do i=1,size(xo)

    sub_x0 = xo(i)-length/2d0
    call hunt(xi,sub_x0,ipos0)
    sub_xN = xo(i)+length/2d0
    call hunt(xi,sub_xN,iposN)

    sub_N = iposN - ipos0 +2
    allocate(sub_xi(sub_N),sub_yi(sub_N))

    sub_xi(1) = sub_x0
    eta = ( sub_xi(1)-xi(ipos0) ) / ( xi(ipos0+1)-xi(ipos0) )
    sub_yi(1) = yi(ipos0) + eta*( yi(ipos0+1)-yi(ipos0) ) 

    sub_xi(2:sub_N-1) = xi(ipos0+1:iposN) 
    sub_yi(2:sub_N-1) = yi(ipos0+1:iposN) 

    sub_xi(sub_N) = sub_xN
    eta = ( sub_xi(sub_N)-xi(iposN) ) / ( xi(iposN+1)-xi(iposN) )
    sub_yi(sub_N) = yi(iposN) + eta*( yi(iposN+1)-yi(iposN) ) 

    yo(i) = 0.5d0*SUM( (sub_yi(1:sub_N-1) + sub_yi(2:sub_N)) &
                    *(sub_xi(2:sub_N) - sub_xi(1:sub_N-1)) ) &
            / length

    deallocate(sub_xi,sub_yi)
  enddo

  end subroutine smooth_interpol

!
!***********************************************************************
!
!! linear_dist_type destructor
subroutine destroy_linear_dist(d)
  type(linear_dist_type), pointer :: d
  deallocate(d%x,d%val)
  deallocate(d)
end subroutine destroy_linear_dist

end module distribution_linear
