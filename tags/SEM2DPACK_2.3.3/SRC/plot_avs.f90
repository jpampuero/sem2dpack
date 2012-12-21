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
module plot_avs

  implicit none

contains

! routine sauvegarde fichier AVS
  subroutine plotavs(field,grid,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit
  use echo, only : echo_run,iout,fmt1,fmtok

  type (sem_grid_type), intent(in) :: grid
  
  integer, intent(in) :: it
  double precision, intent(in) :: field(:,:)

  integer :: icell,i,j,e,fid,ip
  double precision :: rmaxnorm
  character(50) :: name

  fid = IO_new_unit()
  write(name,"('AVS_',I5.5,'_sem2d.inp')") it
  if (echo_run) write(iout,fmt1,advance='no') "Dump AVS file "//trim(name)
  open(unit=fid,file=trim(name),status='unknown')

  ! number of nodes, cells, data per cell
  write(fid,180) grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),1,0,0

  ! index and coordinates of nodes (3D coord = 0)
  do ip=1,grid%npoin
    write(fid,200) ip,grid%coord(:,ip)
  enddo

  ! index and topology of the cells
  icell = 0
  do e=1,grid%nelem
  do i=1,grid%ngll-1
  do j=1,grid%ngll-1

    icell = icell + 1
    write(fid,210) icell,grid%tag(e),grid%ibool(i,j+1,e), &
      grid%ibool(i,j,e),grid%ibool(i+1,j,e),grid%ibool(i+1,j+1,e)

  enddo
  enddo
  enddo

  ! dummy structure data vector and labels 
  write(fid,*) ' 1 1'
  write(fid,*) ' Label1, Label2'

  ! normalized node data (scalar)
  rmaxnorm = maxval(abs(field))
  if ( size(field,2)==1 ) then 
    do ip=1,grid%npoin
      write(fid,205) ip,sqrt(field(ip,1)**2 +field(ip,2)**2)/rmaxnorm
    enddo
  else
    do ip=1,grid%npoin
      write(fid,205) ip,field(ip,1)/rmaxnorm
    enddo
  endif

  close(fid)

  if (echo_run) write(iout,fmtok)

180 format(5(I0,1x))
200 format(I0,1x,e12.5,' 0. ',e12.5)
205 format(I0,1x,e12.5)
210 format(I0,1x,I0,' quad ',4(I0,1x))

  end subroutine plotavs

end module plot_avs
