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
module plot_avs

  implicit none

contains

! routine sauvegarde fichier AVS
  subroutine plotavs(displ,grid,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit
  use fields_class, only : FIELDS_get_max

  type (sem_grid_type), intent(in) :: grid
  
  integer, intent(in) :: it
  double precision, intent(in) :: displ(:,:)

  integer :: icell,i,j,ispec,iavsfile,ip
  double precision :: rmaxnorm
  character(50) :: name

  print *,'Entering AVS file generation...'

  !---- ouverture du fichier AVS
  iavsfile = IO_new_unit()
  write(name,"('AVS_',I5.5,'_sem2d.inp')") it
  open(unit=iavsfile,file=trim(name),status='unknown')

  ! nb de noeuds, de cellules, de donnees par cellule
  write(iavsfile,180) grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),1,0,0

  ! numero et coordonnees des points du maillage (3D fictif avec coordonnee nulle)
  do ip=1,grid%npoin
    write(iavsfile,200) ip,grid%coord(:,ip)
  enddo

  ! numero et topologie des cellules
  icell = 0
  do ispec=1,grid%nelem
  do i=1,grid%ngll-1
  do j=1,grid%ngll-1

    icell = icell + 1
    write(iavsfile,210) icell,grid%tag(ispec),grid%ibool(i,j+1,ispec), &
      grid%ibool(i,j,ispec),grid%ibool(i+1,j,ispec),grid%ibool(i+1,j+1,ispec)

  enddo
  enddo
  enddo

  ! structure data vector et labels bidons
  write(iavsfile,*) ' 1 1'
  write(iavsfile,*) ' Label1, Label2'

  ! donnees aux noeuds (norme du vecteur deplacement, normalisee a 1)
  rmaxnorm = FIELDS_get_max(displ)
  do ip=1,grid%npoin
    write(iavsfile,205) ip,sqrt(displ(ip,1)**2 +displ(ip,2)**2)/rmaxnorm
  enddo

  print *,'Max norme dans output AVS = ',rmaxnorm

  close(iavsfile)

  print *,'End of AVS file generation...'

180 format(5(I0,1x))
200 format(I0,1x,e12.5,' 0. ',e12.5)
205 format(I0,1x,e12.5)
210 format(I0,1x,I0,' quad ',4(I0,1x))

  end subroutine plotavs

end module plot_avs
