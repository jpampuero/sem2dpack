module problem_class
#include <petsc/finclude/petscksp.h>
    use petscksp

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
  use energy      , only : energy_type

  implicit none
  private

  type problem_type ! a huge structure that stores all the data needed to solve a problem
   ! input data for the macro-mesh
    type(mesh_type) :: mesh
   ! data for the spectral element mesh
    type(sem_grid_type) :: grid
   ! input data for each material
    type(matpro_input_type), pointer :: matinp(:) => null()
   ! material properties for each element (slow access)
    type(matpro_elem_type), pointer :: matpro(:) => null()
   ! working data for each element (rapid access)
    type(matwrk_elem_type), pointer :: matwrk(:) => null()
   ! data for boundary conditions
    type(bc_type), pointer :: bc(:) => null()
   ! energies over the whole domain
    type(energy_type) :: energy
   ! physical variables, global pointwise storage
    type(fields_type) :: fields
   ! inverse mass matrix (diagonal)
    double precision, pointer :: rmass(:,:) => null()

   ! do not need after using petsc 
   ! inverse stiffness matrix (diagonal)
    double precision, pointer :: invKdiag(:,:) => null()
    double precision, pointer :: invKdiagTrans(:,:) => null()

   ! petsc ksp linear solver
    KSP  :: ksp

   ! global stiffness matrix, [npoin*ndof, npoin*ndof]
   ! initialize subroutine in mat_gen 
    Mat :: K

   ! global transformation matrix, [npoin*ndof, npoin*ndof]
   ! initialize subroutine in bc_gen
    Mat :: X, Xinv

   ! global displacement, velocity, acceleration, [npoin*ndof, 1]
   ! initialize subroutine in fields
    Vec :: d, v, a, b 

   ! time integration coefficients
    type(timescheme_type) :: time
   ! data for sources
    type(source_type), pointer :: src(:) => null()
   ! data for receivers
    type(rec_type), pointer :: rec => null()
  end type problem_type

  public :: problem_type, destroyPetscStruct

  contains

  subroutine destroyPetscStruct(pb, ierr)
#include <petsc/finclude/petscksp.h>
     use petscksp
     implicit none
     type(problem_type)::pb
     PetscErrorCode :: ierr
     call MatDestroy(pb%K, ierr)
     call VecDestroy(pb%d, ierr)
     call VecDestroy(pb%v, ierr)
     call VecDestroy(pb%a, ierr)
     call VecDestroy(pb%b, ierr)
     call KSPDestroy(pb%ksp, ierr)
   end subroutine destroyPetscStruct

end module problem_class
