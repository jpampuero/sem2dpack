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
module plot_visual3

  implicit none

contains

! routine sauvegarde fichier Visual3
  subroutine plotvisual3(field,grid,mat,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit
  use prop_mat, only : matpro_elem_type, MAT_getProp
  use echo, only : echo_run,iout,fmt1,fmtok

  type (sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  double precision    , intent(in) :: field(:,:)
  integer             , intent(in) :: it

  double precision rmaxnorm,zoffset,xdiff,zdiff
  double precision :: celem(grid%ngll,grid%ngll)
  integer :: icell,i,j,e,fid,ip
  character(50) :: name

  print *,'Entering Visual3 file generation...'

  fid = IO_new_unit()
  write(name,"('Visual3_',i5.5,'_sem2d.inp')") it
  if (echo_run) write(iout,fmt1,advance='no') "Dump Visual3 file "//trim(name)
  open(unit=fid,file=trim(name),status='unknown')

  ! dummy thickness = size of first element
  xdiff = grid%coord(1,grid%ibool(grid%ngll,1,1)) - grid%coord(1,grid%ibool(1,1,1))
  zdiff = grid%coord(2,grid%ibool(grid%ngll,1,1)) - grid%coord(2,grid%ibool(1,1,1))
  zoffset = sqrt(xdiff**2 + zdiff**2)

  ! number of nodes, cells, data per cell
  !write(fid,180) 2*grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),zoffset
  write(fid,180) size(field,2)*grid%npoin &
                 ,grid%nelem*(grid%ngll-1)*(grid%ngll-1),zoffset

  ! index and coordinates of nodes
  do ip=1,grid%npoin
    write(fid,200) ip,grid%coord(:,ip)
  enddo

  ! index and topology of the cells
  icell = 0
  do e=1,grid%nelem
    call MAT_getProp(celem,mat(e),'cp')

    do i=1,grid%ngll-1
    do j=1,grid%ngll-1

      icell = icell + 1

      write(fid,210) icell,celem(i,j), &
        grid%ibool(i,j,e),grid%ibool(i+1,j,e), &
        grid%ibool(i+1,j+1,e),grid%ibool(i,j+1,e)

    enddo
    enddo
  enddo

  ! normalized node data
  rmaxnorm = maxval(abs(field))
  do ip=1,grid%npoin
    write(fid,205) ip,field(ip,:)/rmaxnorm
  enddo

  close(fid)

  if (echo_run) write(iout,fmtok)

180 format(I0,1x,I0,1x,e12.5)
200 format(I0,1x,e12.5,1x,e12.5)
205 format(I0,1x,e12.5,1x,e12.5)
210 format(I0,1x,e12.5,1x,I0,1x,I0,1x,I0,1x,I0)

  end subroutine plotvisual3

end module plot_visual3
