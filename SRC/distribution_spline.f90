! SEM2DPACK version 2.3.3 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module distribution_spline

  use stdio, only: IO_abort, IO_new_unit, IO_file_length
  implicit none
  private
 
  type spline_dist_type
    private
    double precision, pointer :: X(:),val(:)
  end type spline_dist_type

  public :: spline_dist_type,read_spline_dist,generate_spline_dist ,&
            destroy_spline_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_SPLINE
! GROUP  : DISTRIBUTIONS_1D
! PURPOSE: Spline interpolated 1D distribution along X.
! SYNTAX : &DIST_SPLINE file /
!
! ARG: file     [name] [none] Name of the ASCII file containing
!               the data to be interpolated, two columns per line: 
!               (1) X position, sorted in increasing order, and 
!               (2) data value at X
!
! END INPUT BLOCK

  subroutine read_spline_dist (data, iin)

  type(spline_dist_type) :: data
  integer , intent(in) :: iin

  character(50) :: file
  integer :: iunit, i,N

  NAMELIST / DIST_SPLINE / file

  read(iin,DIST_SPLINE)

  N = IO_file_length(file)
  allocate( data%X(N),data%val(N) )
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  do i= 1,N
    read (iunit,*) data%x(i),data%val(i)
  end do
  close(iunit)

  end subroutine read_spline_dist

!
!***********************************************************************
!
  subroutine generate_spline_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(spline_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  integer :: unit,i 

  if (size(coord,1) == 2) then
    ! WARNING: interpolating along X
    call interpol(par%X,par%val,coord(1,:),field)
    unit = IO_new_unit()
    open(unit,file='DistSpline_sem2d.tab',status='replace')
    do i= 1,size(field)
      write(unit,*) coord(1,i),field(i)
    enddo
    close(unit)
  else
    call IO_abort('generate_spline_distribution: not implemented for 3D')
  endif
 
  end subroutine generate_spline_dist

!
!***********************************************************************
!
  subroutine interpol(xi,yi,xo,yo)

  use utils, only : spline, splint

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(out) :: yo(:)

  double precision :: y2(size(xi))
  integer :: Ni,No,i

  Ni = size(xi)
  No = size(xo)
  ! WARNING: 0 derivative at extreme points of spline
  call spline(xi,yi,Ni,0d0,0d0,y2) 
  ! or natural spline:
  !call spline(xi,yi,size(xi),1d30,1d30,y2) 
  do i=1,No
    call splint(xi,yi,y2,Ni,xo(i),yo(i))
  enddo

  end subroutine interpol

!
!***********************************************************************
!
!! spline_dist_type destructor
subroutine destroy_spline_dist(d)
  type(spline_dist_type), pointer :: d
  deallocate(d%x,d%val)
  deallocate(d)
end subroutine destroy_spline_dist


end module distribution_spline
