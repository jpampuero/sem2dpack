module petsc_objects
#include <petsc/finclude/petscksp.h>
    use petscksp

! A class that stores the petsc objects  
!
  implicit none
  private

  type petsc_objects_type ! a data structure that stores parallel matrix and vectors 
   ! petsc ksp linear solver
    KSP  :: ksp

   ! global stiffness matrix, [npoin*ndof, npoin*ndof]
   ! initialize subroutine in mat_gen 
    Mat :: K, MatA

   ! global transformation matrix, [npoin*ndof, npoin*ndof]
   ! initialize subroutine in bc_gen
    Mat :: X, Xinv

   ! global displacement, velocity, acceleration, [npoin*ndof, 1]
   ! initialize subroutine in fields
    Vec :: d, b 

  end type petsc_objects_type

  public :: petsc_objects_type, destroyPetscStruct

  contains

  subroutine destroyPetscStruct(po, ierr)
#include <petsc/finclude/petscksp.h>
     use petscksp
     implicit none
     type(petsc_objects_type)::po
     PetscErrorCode :: ierr
     call MatDestroy(po%K, ierr)
     call MatDestroy(po%MatA, ierr)
     call VecDestroy(po%d, ierr)
     call VecDestroy(po%b, ierr)
     call KSPDestroy(po%ksp, ierr)
   end subroutine destroyPetscStruct

end module petsc_objects
