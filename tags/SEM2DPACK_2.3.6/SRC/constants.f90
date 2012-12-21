! SEM2DPACK version 2.3.6 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module constants
  
!-- Optimization settings ---------------------------------------------------------------

! Optimize the evaluation of internal elastic forces (K*d) for this value of NGLL 
  integer, parameter :: OPT_NGLL=5
! (NGLL = number of Gauss-Lobatto-Legendre nodes per element edge =polynomial degree+1)
! Method: declare NGLL statically to allow for compiler optimizations

! Renumber elements ? 
  logical, parameter :: OPT_RENUMBER = .true.  
! Renumbering improves data locality in an attempt to optimize cache usage
! Method: reverse Cuthill-McKee algorithm.
! This is only applied on structured grids 
! (most mesh generators for unstructured grids, like EMC2, have a renumbering tool)


      
!-- Expert settings: ---------------------------------------------------------------

! Compute energies or not ?
  logical, parameter :: COMPUTE_ENERGIES = .false.
! The evaluation of energies is expensive, so this flag is usually turned off (false)
! Energies are exported to file 'energy_sem2d.tab', 4 columns:
! time, cumulative plastic energy, kinetic energy, total change of elastic energy

! Compute stress glut or not ?
  logical, parameter :: COMPUTE_STRESS_GLUT = .false.
! If this flag is true the cumulative stress glut is computed 
! and exported to file 'stress_glut_sem2d.tab', 7 columns:
! time, sp_11, sp_22, sp_12, sd_11, sd_22, sd_12
! where sp = plasticity-related stress glut
! and sd = damage-related stress glut

! Tolerance (in meters) for comparing locations and spatial coordinates
  double precision, parameter :: TINY_XABS = 1d-6



!-- Other constants.  ---------------------------------------------------------------
! Do NOT modify anything below this line.
  integer, parameter :: NDIME = 2
  double precision, parameter :: PI = 3.141592653589793d0

end module constants
