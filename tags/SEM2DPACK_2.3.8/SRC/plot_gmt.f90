! SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module plot_gmt

  implicit none
  private

  public :: PLOT_GMT_init

contains
  
subroutine PLOT_GMT_init(ibool)

  use stdio, only : IO_new_unit
  use echo, only : echo_init,iout,fmt1,fmtok

  integer, intent(in) :: ibool(:,:,:)

  integer :: ngll,nelem,i,j,e,fid

  ngll = size(ibool,1)
  nelem = size(ibool,3)

  fid = IO_new_unit()
  open(unit=fid,file='grid_sem2d.gmt',status='unknown')
  if (echo_init) write(iout,fmt1,advance='no') "Write GMT triangulation file grid_sem2d.gmt"

  do e =1,nelem
  do j =1,ngll-1
  do i =1,ngll-1
    write(fid,*) ibool(i,j,e),ibool(i+1,j,e),ibool(i,j+1,e)
    write(fid,*) ibool(i+1,j,e),ibool(i+1,j+1,e),ibool(i,j+1,e)
  enddo
  enddo
  enddo
  
  close(fid)

  if (echo_init) write(iout,fmtok)

end subroutine PLOT_GMT_init

end module plot_gmt