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
module problem_class

! PROBLEM_CLASS:  
!
  use mesh_gen    , only : mesh_type
  use spec_grid   , only : sem_grid_type
  use bc_gen      , only : bc_type
  use fields_class, only : fields_type
  use elastic     , only : elast_type
  use kelvin_voigt, only : kelvin_voigt_type
  use time_evol   , only : timescheme_type
  use sources     , only : source_type
  use receivers   , only : rec_type

  implicit none
  private

  type problem_type
    type(mesh_type)            :: mesh
    type(sem_grid_type)        :: grid
    type(bc_type)    , pointer :: bc(:) => null()
    type(fields_type)          :: fields
    double precision, pointer  :: rmass(:,:) => null()
    type(timescheme_type)      :: time
    type(elast_type)           :: elast
    type(kelvin_voigt_type), pointer :: kelvin_voigt => null()
    type(source_type), pointer :: src(:) => null()
    type(rec_type)   , pointer :: rec => null()
  end type problem_type

  public :: problem_type

end module problem_class
