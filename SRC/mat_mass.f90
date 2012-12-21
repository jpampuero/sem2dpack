! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module mat_mass

  use prop_mat

  implicit none
  private

  public :: MAT_MASS_init, MAT_MASS_init_elem_prop

contains

!=====================================================================
  subroutine MAT_MASS_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_setProp(elem,'rho',ecoord)
  
  end subroutine MAT_MASS_init_elem_prop

!=====================================================================
!
! Compute and assemble the mass matrix
!
  subroutine MAT_MASS_init(mass,mat,grid,ndof)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_VolumeWeights
  use prop_mat
  use echo, only : iout,echo_init, fmt1,fmtok
  use fields_class, only : FIELD_add_elem

  double precision, pointer :: mass(:,:)
  type(sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  integer, intent(in) :: ndof
  
  double precision :: massloc(grid%ngll,grid%ngll)
  integer :: e

  if (echo_init) write(iout,fmt1,advance='no') 'Building the mass matrix'

  allocate(mass(grid%npoin,ndof))
  call storearray('mass',size(mass),idouble)
  mass = 0.d0

  do e = 1,grid%nelem
    call MAT_getProp(massloc,mat(e),'rho')
    massloc = massloc * SE_VolumeWeights(grid,e)
    call FIELD_add_elem( massloc, mass(:,1), grid%ibool(:,:,e) )
  enddo

  if (ndof==2) mass(:,2) = mass(:,1)

  if (echo_init) write(iout,fmtok)

  end subroutine MAT_MASS_init


end module mat_mass
