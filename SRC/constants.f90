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
  double precision, parameter :: TINY_XABS = 1d-3



!-- Other constants.  ---------------------------------------------------------------
! Do NOT modify anything below this line.
  integer, parameter :: NDIME = 2
  double precision, parameter :: PI = 3.141592653589793d0

end module constants
