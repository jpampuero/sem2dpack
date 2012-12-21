! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                             with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics Group
! ETH Hönggerberg HPP O 13.1
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
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
module plot_avs

  implicit none

contains

! routine sauvegarde fichier AVS
  subroutine plotavs(field,grid,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit

  type (sem_grid_type), intent(in) :: grid
  
  integer, intent(in) :: it
  double precision, intent(in) :: field(:,:)

  integer :: icell,i,j,e,iout,ip
  double precision :: rmaxnorm
  character(50) :: name

  print *,'Entering AVS file generation...'

  iout = IO_new_unit()
  write(name,"('AVS_',I5.5,'_sem2d.inp')") it
  open(unit=iout,file=trim(name),status='unknown')

  ! number of nodes, cells, data per cell
  write(iout,180) grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),1,0,0

  ! index and coordinates of nodes (3D coord = 0)
  do ip=1,grid%npoin
    write(iout,200) ip,grid%coord(:,ip)
  enddo

  ! index and topology of the cells
  icell = 0
  do e=1,grid%nelem
  do i=1,grid%ngll-1
  do j=1,grid%ngll-1

    icell = icell + 1
    write(iout,210) icell,grid%tag(e),grid%ibool(i,j+1,e), &
      grid%ibool(i,j,e),grid%ibool(i+1,j,e),grid%ibool(i+1,j+1,e)

  enddo
  enddo
  enddo

  ! dummy structure data vector and labels 
  write(iout,*) ' 1 1'
  write(iout,*) ' Label1, Label2'

  ! normalized node data (scalar)
  rmaxnorm = maxval(abs(field))
  if ( size(field,2)==1 ) then 
    do ip=1,grid%npoin
      write(iout,205) ip,sqrt(field(ip,1)**2 +field(ip,2)**2)/rmaxnorm
    enddo
  else
    do ip=1,grid%npoin
      write(iout,205) ip,field(ip,1)/rmaxnorm
    enddo
  endif

  print *,'Max norm in output AVS = ',rmaxnorm

  close(iout)

  print *,'End of AVS file generation...'

180 format(5(I0,1x))
200 format(I0,1x,e12.5,' 0. ',e12.5)
205 format(I0,1x,e12.5)
210 format(I0,1x,I0,' quad ',4(I0,1x))

  end subroutine plotavs

end module plot_avs
