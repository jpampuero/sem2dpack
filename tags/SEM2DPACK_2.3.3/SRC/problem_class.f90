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
module problem_class

! PROBLEM_CLASS:  
!
  use mesh_gen    , only : mesh_type
  use spec_grid   , only : sem_grid_type
  use prop_mat    , only : matpro_input_type, matpro_elem_type
  use mat_gen     , only : matwrk_elem_type
  use bc_gen      , only : bc_type
  use fields_class, only : fields_type
  use time_evol   , only : timescheme_type
  use sources     , only : source_type
  use receivers   , only : rec_type

  implicit none
  private

  type problem_type
   ! input data for the macro-mesh
    type(mesh_type) :: mesh
   ! data for the spectral element mesh
    type(sem_grid_type) :: grid
   ! input data for each material
    type(matpro_input_type), pointer :: matinp(:) 
   ! material properties for each element (slow access)
    type(matpro_elem_type), pointer :: matpro(:)
   ! working data for each element (rapid access)
    type(matwrk_elem_type), pointer :: matwrk(:)
   ! data for boundary conditions
    type(bc_type), pointer :: bc(:) => null()
   ! physical variables, global pointwise storage
    type(fields_type) :: fields
   ! inverse mass matrix (diagonal)
    double precision, pointer :: rmass(:,:) => null()
   ! time integration coefficients
    type(timescheme_type) :: time
   ! data for sources
    type(source_type), pointer :: src(:) => null()
   ! data for receivers
    type(rec_type), pointer :: rec => null()
  end type problem_type

  public :: problem_type

end module problem_class
