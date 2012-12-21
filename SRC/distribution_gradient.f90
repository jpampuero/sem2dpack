! SEM2DPACK version 2.3.7 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module distribution_gradient
  
  implicit none
  private
 
  type gradient_dist_type
    private
    double precision :: grad,angle,valref
    double precision, pointer :: X(:),Y(:)
  end type gradient_dist_type

  public :: gradient_dist_type,read_gradient_dist,generate_gradient_dist ,&
            destroy_gradient_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_GRADIENT
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Constant gradient 2D distribution.
! SYNTAX : &DIST_GRADIENT file,valref ,grad,angle/
!
! ARG: file             [name] [none]    Name of the file containing the coordinates
!                        of the points defining the reference line.
!                        It is an ASCII file with 2 columns per line:
!                        (1) X position (in m) and
!                        (2) Z position (in m)
! ARG: valref           [dble] [none]    Value along the reference line
! ARG: grad             [dble >0] [none] Positive gradient (valref_units/meter)
! ARG: angle            [dble] [none]    Angle (degrees) between the vertical down 
!                        and the grad+ direction. Anticlockwise convention (grad+
!                        points down if 0, right if 90)
!
! NOTE   : Make sure the angle and ref-line are compatible. The code will
!          abort if the ref-line is too short: some points of the domain
!          cannot be projected to ref-line in the angle direction.
!
! END INPUT BLOCK

  subroutine read_gradient_dist (d, iin)

  use utils, only: dsort
  use stdio, only: IO_abort, IO_new_unit, IO_file_length

  type(gradient_dist_type) :: d
  integer , intent(in) :: iin 

  double precision, allocatable :: xread(:),yread(:)
  double precision :: grad,angle,cosa,sina,valref
  integer :: iunit, i,N
  character(50) :: file

  NAMELIST / DIST_GRADIENT / file,grad,angle,valref

  angle = 0
  read(iin,DIST_GRADIENT, END = 100)
  d%grad  = abs(grad)
  d%angle = angle*3.141592653589793d0/180.d0
  d%valref = valref

! Read the reference line
  N = IO_file_length(file)
  allocate( d%x(N),d%y(N) )
  allocate( xread(N),yread(N) )
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  do i= 1,N
    read (iunit,*) xread(i),yread(i)
  end do
  close(iunit)

! rotate frame
  cosa = cos(d%angle)
  sina = sin(d%angle)
  d%x = cosa*xread + sina*yread
  d%y =-sina*xread + cosa*yread
  deallocate(yread,xread)

! sort X in ascending order, carry Y
  call dsort(d%x,d%y) 

  return

  100 call IO_abort('read_gradient_dist: DIST_GRADIENT parameters missing')

  end subroutine read_gradient_dist

!
!***********************************************************************
!
  subroutine generate_gradient_dist(field, coord, d)

  use stdio, only: IO_abort

  double precision, intent(in), dimension(:,:) :: coord
  type(gradient_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  double precision, dimension(size(coord,2)) :: xrot,yrot,yref
  double precision :: cosa,sina

! rotate frame
  cosa = cos(d%angle)
  sina = sin(d%angle)
  xrot = cosa*coord(1,:) + sina*coord(2,:)
  yrot =-sina*coord(1,:) + cosa*coord(2,:)
  
  if (minval(xrot)<d%x(1) .or. maxval(xrot)>d%x(size(d%x)) ) &
    call IO_abort('generate_gradient_dist: reference line is too short for this angle')

  call interpol( d%x,d%y, xrot,yref ) ! get a reference point in the reference line
  field = d%valref + d%grad*( yref-yrot ) ! apply the gradient

  end subroutine generate_gradient_dist

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

!***********************************************************************
! gradient_dist_type destructor
subroutine destroy_gradient_dist(d)
  type(gradient_dist_type), pointer :: d
  deallocate(d%x,d%y)
  deallocate(d)
end subroutine destroy_gradient_dist

end module distribution_gradient
