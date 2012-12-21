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
module plot_visual3

  implicit none

contains

! routine sauvegarde fichier Visual3
  subroutine plotvisual3(displ,grid,elast,it)

  use spec_grid, only : sem_grid_type
  use elastic, only : elast_type,ELAST_inquire
  use fields_class, only : FIELDS_get_max

  type (sem_grid_type), intent(in) :: grid
  type (elast_type)   , intent(in) :: elast
  double precision    , intent(in) :: displ(:,:)
  integer             , intent(in) :: it

  double precision rmaxnorm,zoffset,xdiff,zdiff
  double precision :: cploc
  integer :: icell,i,j,e,ivis3file,ip
  character(50) :: name

  print *,'Entering Visual3 file generation...'

  !---- ouverture du fichier Visual3
  ivis3file = 34
  write(name,"('Visual3_',i5.5,'_sem2d.inp')") it
  open(unit=ivis3file,file=trim(name),status='unknown')

  ! prendre comme epaisseur fictive la taille du premier element spectral
  xdiff = grid%coord(1,grid%ibool(grid%ngll,1,1)) - grid%coord(1,grid%ibool(1,1,1))
  zdiff = grid%coord(2,grid%ibool(grid%ngll,1,1)) - grid%coord(2,grid%ibool(1,1,1))
  zoffset = sqrt(xdiff**2 + zdiff**2)

  ! nb de noeuds, de cellules, de donnees par cellule
  write(ivis3file,180) 2*grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),zoffset

  ! numero et coordonnees des points du maillage
  do ip=1,grid%npoin
    write(ivis3file,200) ip,grid%coord(:,ip)
  enddo

  ! numero et topologie des cellules
  icell = 0
  do e=1,grid%nelem
  do i=1,grid%ngll-1
  do j=1,grid%ngll-1

    icell = icell + 1
    call ELAST_inquire(elast,i,j,e,cp=cploc)

    write(ivis3file,210) icell,cploc, &
        grid%ibool(i,j,e),grid%ibool(i+1,j,e), &
        grid%ibool(i+1,j+1,e),grid%ibool(i,j+1,e)

  enddo
  enddo
  enddo

  ! donnees aux noeuds (composantes du vecteur deplacement, normalisees a 1)
  rmaxnorm = FIELDS_get_max(displ)
  do ip=1,grid%npoin
    write(ivis3file,205) ip,displ(ip,:)/rmaxnorm
  enddo

  print *,'Max norme dans output Visual3 = ',rmaxnorm

  close(ivis3file)

  print *,'End of Visual3 file generation...'

180 format(I0,1x,I0,1x,e12.5)
200 format(I0,1x,e12.5,1x,e12.5)
205 format(I0,1x,e12.5,1x,e12.5)
210 format(I0,1x,e12.5,1x,I0,1x,I0,1x,I0,1x,I0)

  end subroutine plotvisual3

end module plot_visual3
