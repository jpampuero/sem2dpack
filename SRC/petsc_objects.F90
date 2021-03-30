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
    Vec :: d, b, d0, b0 
    VecScatter :: d_ctx, b_ctx, bd_fix_ctx 
    IS :: isIndexDofFix

  end type petsc_objects_type
  double precision :: FillRatio = 1.01

  public :: petsc_objects_type, init_petsc_objects, destroyPetscStruct, UpdateKsp

  contains

  subroutine init_petsc_objects(petobj, pb, ierr)
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscsys.h>
  use petscksp
  use petscsys
  use fields_class
  use problem_class , only : problem_type
  use echo, only : iout,info=>echo_init, fmt1,fmtok
  use mat_gen, only : MAT_init_KG, MAT_AssembleK
  use bc_gen, only : BC_build_transform_mat

  type (problem_type), intent(in) :: pb
  type(petsc_objects_type), intent(inout) :: petobj
  integer :: ndof, ndim, ngll, npoin, rank
  PetscViewer :: viewer
  character(160)::ksptype
  character(160)::pctype
  PC            ::ksp_pc
  PetscErrorCode :: ierr
  real :: cputime0, cputime1
  integer:: start, fini

  ndof  = pb%ndof
  ndim  = 2
  ngll  = pb%ngll
  npoin = pb%npoin

  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)

 ! --------------------------------------------------------------------------
 ! initialize Petsc Vec objects
  call VecCreate(PETSC_COMM_WORLD, petobj%d, ierr)
  call VecSetSizes(petobj%d, PETSC_DECIDE, ndof * npoin,ierr)
  call VecSetBlockSize(petobj%d, ndof, ierr)
  call VecSetFromOptions(petobj%d, ierr)
  
  ! set initial values to 0
  call VecSet(petobj%d, 0d0, ierr)
  call VecDuplicate(petobj%d, petobj%b, ierr) ! right hand side
  
  ! Create scatter context and local sequential vectors
  call VecScatterCreateToZero(petobj%d, petobj%d_ctx, petobj%d0, ierr)
  call VecScatterCreateToZero(petobj%b, petobj%b_ctx, petobj%b0, ierr)

  ! Create IS
  call ISCreateGeneral(PETSC_COMM_WORLD, size(pb%IndexDofFix), &
             pb%IndexDofFix - 1, PETSC_COPY_VALUES, petobj%isIndexDofFix, ierr)

  ! Create scatter context from b to d in indexDofFix
  call VecScatterCreate(petobj%d, petobj%isIndexDofFix, petobj%b, &
             petobj%isIndexDofFix, petobj%bd_fix_ctx, ierr)

  if (rank==0) then
      if (info) then
        write(iout,*)
        write(iout,'(a)') '***********************************************'
        write(iout,'(a)') '*   Start Assemble global stiffness matrix    *' 
        write(iout,'(a)') '***********************************************'
        write(iout,*)
      endif
  end if

 ! --------------------------------------------------------------------------
 ! initialize assemble stiffness matrix
 ! this has to be done once for linear problem
 ! implemented in mat_gen
  call MAT_init_KG(petobj%K, ndof, npoin, ngll, ierr) 
  CHKERRA(ierr)
  call CPU_TIME(cputime0)
  call MAT_AssembleK(petobj%K, pb%matwrk, ndof, ngll, ndim, pb%grid%ibool, ierr)
  CHKERRA(ierr)

  ! print the ownership of matrix K
  call MatGetOwnershipRangeColumn(petobj%K, start, fini, ierr)

  call CPU_TIME(cputime1)
  cputime0 = cputime1-cputime0
  if (rank==0) then
  write(iout,'(/A,1(/2X,A,EN12.3),/)')   &
        '---  CPU TIME ESTIMATES (in seconds) :', &
        'CPU time for initial assembling stiffness matrix . .', cputime0
  
  if (info) then
    write(iout,*)
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*    Finish Assemble global stiffness matrix  *' 
    write(iout,'(a)') '***********************************************'
    write(iout,*)
  endif

  end if

 call CPU_TIME(cputime0)
 call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'sem2d_StiffnessMat',FILE_MODE_WRITE, viewer,ierr);CHKERRA(ierr)
 call MatView(petobj%K, viewer, ierr);CHKERRA(ierr)
 call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)
 call CPU_TIME(cputime1)
 cputime0 = cputime1-cputime0

 if (rank==0) then
 write(iout,'(/A,1(/2X,A,EN12.3),/)')   &
        '---  CPU TIME ESTIMATES (in seconds) :', &
        'CPU time for writing stiffness matrix . .', cputime0

  if (info) then
    write(iout,*)
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*   Done writing global stiffness matrix  *' 
    write(iout,'(a)') '***********************************************'
    write(iout,*)
  endif

  end if

 ! --------------------------------------------------------------------------
  ! reset initial displacement and velocity vector
  call FIELD_SetVecFromField(petobj%d, pb%fields%displ, ierr)
  CHKERRQ(ierr)

 ! --------------------------------------------------------------------------
 ! initialize the transformation matrix X and Xinv
 ! implemented in bc_gen
 call BC_build_transform_mat(pb%bc, petobj%X, petobj%Xinv, ndof, npoin, ierr)
 CHKERRQ(ierr)
 
 
 ! --------------------------------------------------------------------------
 ! initialize the Petsc KSP solver

 ! compute the A matrix: A = Xinv * K * Xinv
 ! the linear system is Xinv * K * Xinv * (X*d) = Xinv * f
 call CPU_TIME(cputime0)
 call MatMatMatMult(petobj%Xinv, petobj%K, petobj%Xinv, MAT_INITIAL_MATRIX, FillRatio, petobj%MatA, ierr)
 call CPU_TIME(cputime1)
 cputime0 = cputime1-cputime0
 
 call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'sem2d_MatA',FILE_MODE_WRITE, viewer,ierr);CHKERRA(ierr)
 call MatView(petobj%MatA, viewer, ierr);CHKERRA(ierr)
 call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)

! write(iout,'(/A,1(/2X,A,EN12.3),/)')   &
!        '---  CPU TIME ESTIMATES (in seconds) :', &
!        'CPU time for construct A matrix for KSP (linear system) . .', cputime0

 ! zero out the rows and columns of MatA for dirichlet boundary conditions
 ! right hand side vector is modified in the solver using fortran subroutines

 call MatZeroRows(petobj%MatA, size(pb%indexDofFix), pb%indexDofFix - 1, 1d0, petobj%d, petobj%b, ierr)

 !call MatZeroRowsColumns(pb%MatA, size(pb%indexDofFix), pb%indexDofFix - 1, 1d0, pb%d, pb%b, ierr)

 call KSPCreate(PETSC_COMM_WORLD, petobj%ksp, ierr)
 call KSPSetOperators(petobj%ksp,petobj%MatA, petobj%MatA, ierr)

 !   /*
 !    Set linear solver defaults for this problem (optional).
 !    - all these parameters could be specified at runtime via
 !      KSPSetFromOptions(). 
 ! */

 ! set the non-zero initial guess
 call KSPSetInitialGuessNonzero(petobj%ksp, PETSC_TRUE, ierr);

 ! set tolerance
  call KSPSetTolerances(petobj%ksp,pb%time%TolLin,1.d-50,&
                      PETSC_DEFAULT_REAL, pb%time%MaxIterLin,ierr)
  CHKERRQ(ierr)
 
 ! enabling setting solver option from run time flags
  call KSPSetFromOptions(petobj%ksp, ierr)
  call KSPGetType(petobj%ksp,ksptype, ierr);
  call KSPGetPC(petobj%ksp,ksp_pc, ierr);
  call PCGetType(ksp_pc, pctype, ierr)

  if (rank==0) then
      write(iout, *) "KSPTYPE:------------"
      write(iout, *) ksptype
      write(iout, *) "PCTYPE:------------"
      write(iout, *) pctype
  end if

  end subroutine init_petsc_objects

subroutine UpdateKsp(pb,petobj)
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscsys.h>
  use petscksp
  use petscsys
  use problem_class, only : problem_type
  use mat_gen, only : MAT_AssembleK
  type(problem_type)::pb
  type(petsc_objects_type)::petobj
  integer::ndof, ngll, ndim
  PetscErrorCode :: ierr

  ndof  = pb%ndof
  ndim  = 2
  ngll  = pb%ngll

  call MatZeroEntries(petobj%K, ierr)
  call MAT_AssembleK(petobj%K, pb%matwrk, ndof, ngll, ndim, pb%grid%ibool, ierr)

  call MatMatMatMult(petobj%Xinv, petobj%K, petobj%Xinv, MAT_REUSE_MATRIX, &
                    PETSC_DEFAULT_REAL, petobj%MatA, ierr)

  call MatZeroRows(petobj%MatA, size(pb%indexDofFix), pb%indexDofFix - 1, 1d0, petobj%d, petobj%b, ierr)
 ! call KSPCreate(PETSC_COMM_WORLD, pb%ksp, ierr)
  call KSPSetOperators(petobj%ksp, petobj%MatA, petobj%MatA, ierr)

end subroutine UpdateKsp

  subroutine destroyPetscStruct(po, ierr)
#include <petsc/finclude/petscksp.h>
     use petscksp
     implicit none
     type(petsc_objects_type)::po
     PetscErrorCode :: ierr
     call MatDestroy(po%K, ierr)
     call MatDestroy(po%MatA, ierr)
     call MatDestroy(po%X, ierr)
     call MatDestroy(po%Xinv, ierr)
     call KSPDestroy(po%ksp, ierr)
     call VecDestroy(po%d, ierr)
     call VecDestroy(po%b, ierr)
     call VecDestroy(po%d0, ierr)
     call VecDestroy(po%b0, ierr)
     call VecScatterDestroy(po%b_ctx, ierr)
     call VecScatterDestroy(po%d_ctx, ierr)
     call VecScatterDestroy(po%bd_fix_ctx, ierr)
     call ISDestroy(po%isIndexDofFix, ierr)
   end subroutine destroyPetscStruct

end module petsc_objects
