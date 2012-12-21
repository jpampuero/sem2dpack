! SEM2DPACK version 2.3.5 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://web.gps.caltech.edu/~ampuero/
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
module distribution_linear
  
  implicit none
  private
 
  type linear_dist_type
    private
    integer :: dim = 1
    double precision :: smoothing_length = 0d0
    double precision, pointer :: X(:)=>null(),val(:)=>null()
  end type linear_dist_type

  public :: linear_dist_type,read_linear_dist,generate_linear_dist ,&
            destroy_linear_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_LINEAR
! GROUP  : DISTRIBUTIONS_1D
! PURPOSE: Piecewise linear 1D distribution along X or Z.
! SYNTAX : &DIST_LINEAR n,dim,length /
!            followed immediately by the interpolation data, 
!            one line per point, two columns: position (X or Z), value
!            Position must be sorted in increasing order.
!          or
!          &DIST_LINEAR file,dim,length / 
!            and the interpolation data is read from a two-column file
!
! ARG: n        [int] [0] Number of points to be interpolated
! ARG: dim      [int] [1] Interpolate along X (dim=1) or along Z (dim=2)
! ARG: file     [name] [none] Name of the ASCII file containing the data
! ARG: length   [dble] [0] Smoothing length for sliding average window
!               No smoothing if length=0
!
! END INPUT BLOCK

  subroutine read_linear_dist (d, iin)

  use stdio, only: IO_abort, IO_new_unit, IO_file_length

  type(linear_dist_type) :: d
  integer , intent(in) :: iin 

  double precision :: length
  integer :: iunit, i,N,dim
  character(50) :: file

  NAMELIST / DIST_LINEAR / file,length,N,dim

  length = 0d0
  file = ''
  N = 0
  dim = 1
  read(iin,DIST_LINEAR, END = 100)

  d%smoothing_length = length
  if (dim<1 .or. dim>2) call IO_abort('read_linear_dist: dim must be 1 or 2')
  d%dim = dim

 ! WARNING: X must be ordered, not checked
  if (file=='') then
    allocate( d%X(N),d%val(N) )
    do i= 1,N
      read (iin,*) d%x(i),d%val(i)
    end do
    
  else
    N = IO_file_length(file)
    allocate( d%X(N),d%val(N) )
    iunit = IO_new_unit()
    open(iunit,file=file,status='old')
    do i= 1,N
      read (iunit,*) d%x(i),d%val(i)
    end do
    close(iunit)
  endif

  return

  100 call IO_abort('read_linear_dist: DIST_LINEAR parameters missing')

  end subroutine read_linear_dist

!
!***********************************************************************
!
  subroutine generate_linear_dist(field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(linear_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

    if (d%smoothing_length > 0d0) then
      call smooth_interpol(d%X,d%val,coord(d%dim,:),field,d%smoothing_length)
    else
      call interpol(d%X,d%val,coord(d%dim,:),field)
    endif

  end subroutine generate_linear_dist

!
!***********************************************************************
!
! Linear interpolation
! If the requested point is out of bounds extrapolate the two start/end values
  subroutine interpol(xi,yi,xo,yo)

  use utils, only: hunt

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(out) :: yo(:)

  double precision :: eta
  integer :: i,ipos,ni,no

  ni=size(xi)
  no=size(xo)
  ipos = 1
  do i=1,no
    call hunt(xi,xo(i),ipos)
    if (ipos==0) ipos=1
    if (ipos==ni) ipos=ni-1
    eta = ( xo(i)-xi(ipos) ) / ( xi(ipos+1)-xi(ipos) )
    yo(i) = yi(ipos) + eta*( yi(ipos+1)-yi(ipos) ) 
  enddo

  end subroutine interpol

!***********************************************************************
! This is too complicated
! Better solution: smooth the underlying data at startup then interpolate
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
