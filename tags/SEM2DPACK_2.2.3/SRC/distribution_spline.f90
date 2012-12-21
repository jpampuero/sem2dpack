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
module distribution_spline

  use stdio, only: IO_abort, IO_new_unit
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
! NAME   : DIST_SPLINE [distributions]
! PURPOSE: Spline interpolated 1D distribution along X.
! SYNTAX : &DIST_SPLINE file /
!
! ARG: file     [name] [none] Name of the ASCII file containing
!               the data to be interpolated in this format:
!                       n               (number of points to be read)
!                       x(1) v(1)       (coordinate and value)
!                       ...  ...
!                       x(n) v(n)
!               X must be sorted in increasing order
!
! END INPUT BLOCK

  subroutine read_spline_dist (data, iin)

  type(spline_dist_type) :: data
  integer , intent(in) :: iin

  character(50) :: file
  integer :: iunit, i,N

  NAMELIST / DIST_SPLINE / file

  read(iin,DIST_SPLINE)

  allocate( data%X(N),data%val(N) )

  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  read (iunit,*) N
  do i= 1,N
    ! WARNING: X must be ordered
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

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(out) :: yo(:)

  double precision :: y2(size(xi))
  integer :: Ni,No,i

  Ni = size(xi)
  No = size(xo)
  ! WARNING: 0 derivative at extreme points of spline
  call spline(xi,yi,size(xi),0d0,0d0,y2) 
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

!
!***********************************************************************
!
!  The following routines are adapted from Numerical Recipes
!  for double precision, Fortran 90

  SUBROUTINE spline(x,y,n,yp1,ypn,y2)

  INTEGER, intent(in) :: n
  double precision, intent(in) :: yp1,ypn,x(n),y(n)
  double precision, intent(out) :: y2(n)

  INTEGER i,k
  double precision :: p,qn,sig,un
  double precision :: u(n)

  if (yp1 > .99d30) then
    y2(1)=0d0
    u(1)=0d0
  else
    y2(1)=-0.5d0
    u(1)=(3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2d0
    y2(i)=(sig-1d0)/p
    u(i)=(6d0*((y(i+1)-y(i)) &
         /(x(i+1)-x(i))-(y(i)-y(i-1)) &
         /(x(i)-x(i-1))) &
         /(x(i+1)-x(i-1)) &
         -sig*u(i-1))/p
  enddo
  if (ypn > .99e30) then
    qn=0d0
    un=0d0
  else
    qn=0.5d0
    un=(3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1d0)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  
  END SUBROUTINE spline

!---------------------------------------------------------------

  SUBROUTINE splint(xa,ya,y2a,n,x,y)

  INTEGER, intent(in) :: n
  double precision, intent(in) :: x,xa(n),y2a(n),ya(n)
  double precision, intent(out) :: y

  integer :: k,khi,klo
  double precision :: a,b,h

  klo=1
  khi=n
  do while (khi-klo > 1)
    k=(khi+klo)/2
    if(xa(k).gt.x)then
      khi=k
    else
      klo=k
    endif
  enddo
  h=xa(khi)-xa(klo)
  if (h == 0d0) &
    call IO_abort('distribution_spline.splint: input X must be distinct')
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
  
  END SUBROUTINE splint

end module distribution_spline
