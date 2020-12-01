module bc_gen

!! To add a new boundary condition, say "bc_user" :
!!
!!   1. create a module bc_user.f90, following existing examples
!!     It should contain the following:
!!      . subroutine BC_USER_read : reads boundary condition parameters from main input file
!!      . subroutine BC_USER_init : initializes properties and data
!!      . subroutine BC_USER_apply : applies the boundary condition during the solver phase
!!   2. Modify bc_gen.f90 (this file) following the instructions
!!      and templates in the comments that start by "!!"
!!   3. Modify the file Makefile.depend
!!   4. Re-compile

!-- Import boundary conditions modules
  use bc_abso
  use bc_periodic
  use bc_dirneu
  use bc_KINFLT 
  use bc_lsf
  use bc_dynflt
  use bc_dirabs

!!  use bc_user
!! add your new module

  implicit none
  private
  
!-- Object containing all boundary conditions
  type bc_type
    private
    integer :: tag(2)
    integer :: kind
!--- List here all bc types
    type(bc_dirneu_type)  , pointer :: dirneu => null()
    type(bc_kinflt_type)  , pointer :: kinflt => null()
    type(bc_abso_type)    , pointer :: abso => null()
    type(bc_periodic_type), pointer :: perio => null()
    type(bc_lsf_type)     , pointer :: lsf  => null()
    type(bc_dynflt_type)  , pointer :: dynflt => null()
    type(bc_dirabs_type)  , pointer :: dirabs => null()
    !! type(bc_user_type) , pointer :: user
  end type bc_type

  integer, parameter :: IS_EMPTY  = 0, &
                        IS_DIRNEU = 1, &
                        IS_KINFLT = 2, &
                        IS_ABSORB = 3, &
                        IS_PERIOD = 4, &
                        IS_LISFLT = 5, &
                        IS_DYNFLT = 6, &
                        IS_DIRABS = 7
                        !! IS_USER = 8

  public :: bc_type,bc_read, bc_apply, bc_apply_kind, bc_init,bc_write, bc_timestep, & 
            bc_set_kind, bc_select_kind, bc_trans, bc_update_dfault, bc_set_fix_zero, &
            bc_select_fix, bc_has_dynflt, bc_update_bcdv, BC_build_transform_mat, &
            bc_GetIndexDofFix, bc_nDofFix

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DEF
! PURPOSE: Define a boundary condition
! SYNTAX : &BC_DEF tag, tags, kind /
!          possibly followed by &BC_kind blocks 
!
! ARG:  tag     [int] [none] A number assigned to the boundary. If you are
!               using SEM2D built-in structured mesher the conventions are:
!                       1       bottom
!                       2       right
!                       3       up
!                       4       left
!               If you are importing a mesh, you must use the tags assigned
!               to the boundaries during the mesh construction.
! ARG:  tags    [int(2)] [none] Two tags are needed for split-node interfaces (faults)
!               and for periodic boundaries.
! ARG:  kind    [char*6] [none] Type of boundary condition. The following are
!               implemented:
!               'DIRNEU', 'ABSORB', 'PERIOD', 'LISFLT', 'DYNFLT', 'KINFLT'
!               'DIRABS'
!
! NOTE   : Most of the boundary conditions need additional data, given
!          in a BC_kind input block of the BOUNDARY_CONDITIONS group
!          immediately following the BC_DEF block.
!
! END INPUT BLOCK

subroutine bc_read(bc,iunit)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  integer, intent(in) :: iunit
  
  integer :: i,j,nbc,tag,tags(2)
  character(6) :: kind

  NAMELIST / BC_DEF / tag,tags,kind

  !-- count the BCs
  rewind(iunit)
  nbc = 0
  do 
    read(iunit,BC_DEF,END=10) 
    nbc = nbc + 1
  enddo
  10 continue

  !-- leave if there is no BC
  if (nbc == 0) return

  if (echo_input) write(iout,'(//1x,A,/1x,37("="),/)') &
    "B o u n d a r y   C o n d i t i o n s"

  !-- read the BCs definition
  allocate( bc(nbc) )
  do i = 1,nbc

   ! read the i-th BC_DEF
   ! with rewind & loop to account for possible EOFs in BC_XXX_read
    rewind(iunit)
    do j = 1,i
      tag  = 0
      tags = 0
      kind = ' '
      read(iunit,BC_DEF)
    enddo

    bc(i)%tag = 0
    if (tag > 0) then
      bc(i)%tag(1) = tag
    elseif ( tags(1)>0 ) then
      bc(i)%tag    = tags
    else
      call IO_abort('bc_read: tag(s) are null or not set')
    endif
    if (echo_input) then
      if ( bc(i)%tag(2)==0 ) then
        write(iout,200) bc(i)%tag(1)
      else
        write(iout,201) bc(i)%tag
      endif
    endif

    if (kind == ' ') call IO_abort('bc_read: kind not set')
    if (echo_input) write(iout,202) kind
 
   !-- read the BC properties
    select case(kind)
      case('DIRNEU')
        bc(i)%kind = IS_DIRNEU
        allocate(bc(i)%dirneu)
        call BC_DIRNEU_read(bc(i)%dirneu,iunit)
      case('KINFLT')
        bc(i)%kind = IS_KINFLT
        allocate(bc(i)%kinflt)
        call BC_KINFLT_read(bc(i)%kinflt,iunit)
      case('ABSORB')
        bc(i)%kind = IS_ABSORB
        allocate(bc(i)%abso)
        call BC_ABSO_read(bc(i)%abso,iunit)
      case('PERIOD')
        bc(i)%kind = IS_PERIOD
        allocate(bc(i)%perio)
        call BC_PERIO_read() !bc(i)%perio,iunit)
      case('LISFLT')
        bc(i)%kind = IS_LISFLT
        allocate(bc(i)%lsf)
        call BC_LSF_read(bc(i)%lsf,iunit)
      case('DYNFLT')
        bc(i)%kind = IS_DYNFLT
        allocate(bc(i)%dynflt)
        call BC_DYNFLT_read(bc(i)%dynflt,iunit)
      case('DIRABS')
        bc(i)%kind = IS_DIRABS
        allocate(bc(i)%dirabs)
        call BC_DIRABS_read(bc(i)%dirabs,iunit)
!!      case('USERBC')
!!        bc(i)%kind = IS_USER
!!        allocate(bc(i)%user)
!!        call BC_USER_read(bc(i)%user,iunit)
      case default  
        call IO_abort('bc_read: unknown kind')
    end select
  enddo

  return

  200 format(/5x,'Boundary tag. . . . . . . . . . . . (tag) = ',I0)
  201 format(/5x,'Boundary tags . . . . . . . . . . .(tags) = ',I0,' and ',I0)
  202 format( 5x,'Boundary condition. . . . . . . . .(kind) = ',A)

end subroutine bc_read

!-----------------------------------------------------------------------
subroutine bc_init(bc,grid,mat,M,tim,src,d,v)

  use spec_grid, only : sem_grid_type, BC_tag_exists
  use bnd_grid, only : bnd_grid_type
  use prop_mat, only : matpro_elem_type
  use time_evol, only: timescheme_type
  use sources, only: source_type
  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  type(sem_grid_type), intent(inout) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  type(timescheme_type), intent(in) :: tim
  double precision, intent(inout) :: M(:,:)
  type (source_type), pointer :: src(:)
  double precision, dimension(:,:), intent(inout) :: d,v

  type(bc_periodic_type), pointer :: perio
  integer :: i,j
  
  if (.not. associated(bc)) return

 ! check if requested bc tags exist in the mesh
  do i = 1,size(bc)
    do j=1,2
      if (bc(i)%tag(j) /= 0) then
        if (.not. BC_tag_exists(grid%bounds,bc(i)%tag(j)) )  bc(i)%kind = IS_EMPTY
      endif
    enddo
  enddo

 ! first initialize periodic boundaries
  nullify(perio)
  do i = 1,size(bc)
    if(bc(i)%kind==IS_PERIOD) then
      call BC_PERIO_init(bc(i)%perio,bc(i)%tag,grid,M)
      perio => bc(i)%perio
    endif
  enddo
  do i = 1,size(bc)
    select case(bc(i)%kind)
      case(IS_EMPTY)
              ! DEVEL write warning
      case(IS_DIRNEU)
        call BC_DIRNEU_init(bc(i)%dirneu,bc(i)%tag(1),grid,perio)
      case(IS_KINFLT)
        write (* , *) 'Initialize KINFLT'
        call BC_KINFLT_init(bc(i)%kinflt,bc(i)%tag,grid,M,tim,perio)
        write (* , *) 'FINISH Initialize KINFLT'
      case(IS_ABSORB)
        call BC_ABSO_init(bc(i)%abso,bc(i)%tag(1),grid,mat,M,tim,src,perio)
      case(IS_DIRABS)
        call BC_DIRABS_init(bc(i)%dirabs,bc(i)%tag(1),grid,mat,M,tim,src,perio)
      case(IS_LISFLT)
        call BC_LSF_init(bc(i)%lsf,bc(i)%tag,grid,M,perio)
      case(IS_DYNFLT)
          ! Note rate and state fault does not start with zero velocity 
        call BC_DYNFLT_init(bc(i)%dynflt,bc(i)%tag,grid,M,tim,perio, v, mat)
!!      case(IS_USER)
!!        call BC_USER_init(bc(i)%user,bc(i)%tag, ... )
    end select
  enddo

 ! write initial data
  call BC_write(bc,tim,d,v)

end subroutine bc_init


!=======================================================================
!! Applies the boundary condition
subroutine bc_apply(bc,time,fields,force)

  use sources, only: source_type
  use fields_class, only: fields_type
  use time_evol, only : timescheme_type

  type(bc_type), pointer :: bc(:)
  type (fields_type), intent(inout) :: fields
  type(timescheme_type), intent(in) :: time
  double precision, dimension(:,:), intent(inout) :: force

  integer :: i

  if (.not. associated(bc)) return

 ! apply first periodic, then absorbing, then the rest
 ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_PERIOD) call bc_apply_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_ABSORB .or. bc(i)%kind == IS_DIRABS) call bc_apply_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB .and. &
        bc(i)%kind /= IS_DIRABS) call bc_apply_single(bc(i))
  enddo
    
contains

  subroutine bc_apply_single(bc)

    type(bc_type), intent(inout) :: bc

    select case(bc%kind)
      case(IS_DIRNEU)
        call bc_DIRNEU_apply(bc%dirneu, fields%displ, force, time)
      case(IS_KINFLT)
          ! fields%accel means Mxa
        call bc_KINFLT_apply(bc%kinflt,fields%accel,fields%veloc,fields%displ,time)
      case(IS_ABSORB)
        call BC_ABSO_apply(bc%abso,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case(IS_DIRABS)
        call BC_DIRABS_apply(bc%dirabs,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case(IS_PERIOD)
        call bc_perio_apply(bc%perio,force)
      case(IS_LISFLT)
          ! fields%accel means Mxa
        call BC_LSF_apply(bc%lsf,fields%accel,fields%displ_alpha)
      case(IS_DYNFLT)
        call BC_DYNFLT_apply(bc%dynflt,fields%accel,fields%veloc,fields%displ,time)
!!      case(IS_USER)
!!        call BC_USER_apply(bc%user, ... )
    end select

  end subroutine bc_apply_single
  
end subroutine bc_apply


!=======================================================================
!! Applies the boundary condition for a single boundary kind
!
! bc_kind =
! IS_EMPTY  = 0, &
! IS_DIRNEU = 1, &
! IS_KINFLT = 2, &
! IS_ABSORB = 3, &
! IS_PERIOD = 4, &
! IS_LISFLT = 5, &
! IS_DYNFLT = 6
!
subroutine bc_apply_kind(bc, time, fields, force, bc_kind, dirneu_kindin)

  use sources, only: source_type
  use fields_class, only: fields_type
  use time_evol, only : timescheme_type

  type(bc_type), pointer :: bc(:)
  type (fields_type), intent(inout) :: fields
  type(timescheme_type), intent(in) :: time
  double precision, dimension(:,:), intent(inout) :: force
  integer, intent(in) :: bc_kind 
  integer, optional, intent(in)::dirneu_kindin
  integer :: dirneu_kind
  integer :: IS_DIR    = 2, & 
             IS_NEU    = 1 

  integer :: i

  if (.not. associated(bc)) return

  dirneu_kind = IS_NEU 
  if (present(dirneu_kindin)) dirneu_kind = dirneu_kindin 

  ! apply a single boundary kind
  do i = 1,size(bc)
    if ( bc(i)%kind == bc_kind) call bc_apply_single_kind(bc(i))
  enddo

contains

  subroutine bc_apply_single_kind(bc)

    type(bc_type), intent(inout) :: bc
    ! Note at this point fields%accel has been modified by the pointer force
    ! So fields%accel really means Mxa.

    select case(bc%kind)
      case(IS_DIRNEU)
        call bc_DIRNEU_apply_kind(bc%dirneu, fields%displ, force, time, dirneu_kind)
      case(IS_DIRABS)
        call bc_DIRABS_apply_kind(bc%dirabs, fields%displ, force, time, dirneu_kind)
      case(IS_KINFLT)
          ! fields%accel means Mxa
        call bc_KINFLT_apply(bc%kinflt,fields%accel,fields%veloc,fields%displ,time)
      case(IS_ABSORB)
        call BC_ABSO_apply(bc%abso,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case(IS_PERIOD)
        call bc_perio_apply(bc%perio,force)
      case(IS_LISFLT)
          ! fields%accel means Mxa
        call BC_LSF_apply(bc%lsf,fields%accel,fields%displ_alpha)
      case(IS_DYNFLT)
        call BC_DYNFLT_apply(bc%dynflt,fields%accel,fields%veloc,fields%displ,time)
!!      case(IS_USER)
!!        call BC_USER_apply(bc%user, ... )
    end select

  end subroutine bc_apply_single_kind
  
end subroutine bc_apply_kind

!=======================================================================
function BC_nfaultnode(bc) result(n)
  type(bc_type), pointer :: bc(:)
  integer :: i, n

  if (.not. associated(bc)) return
  
  n = 0
  do i = 1,size(bc)
    select case(bc(i)%kind)
    case (IS_KINFLT)
        n = n + bc_kinflt_nnode(bc(i)%kinflt)
    case (IS_DYNFLT)
        n = n + bc_dynflt_nnode(bc(i)%dynflt)
    end select
  enddo

end function


!=======================================================================
! obtain all the fault nodes
subroutine BC_faultnode(bc, node1, node2, nnode)
  type(bc_type), pointer :: bc(:)
  integer, intent(in):: nnode
  integer, dimension(nnode), intent(inout):: node1, node2
  integer, dimension(:), allocatable:: node1_tmp, node2_tmp
  integer :: i, n, istart, iend

  if (.not. associated(bc)) return
  node1 = 0
  node2 = 0
  istart = 1
  
  do i = 1,size(bc)
    select case(bc(i)%kind)
    case (IS_KINFLT)
        n = bc_kinflt_nnode(bc(i)%kinflt)
        iend = istart + n - 1
        allocate(node1_tmp(n))
        allocate(node2_tmp(n))
        call bc_kinflt_node(bc(i)%kinflt, node1_tmp, node2_tmp)
        node1(istart:iend) = node1_tmp
        node2(istart:iend) = node2_tmp
        istart = iend + 1
        deallocate(node1_tmp)
        deallocate(node2_tmp)
    case (IS_DYNFLT)
        n = bc_dynflt_nnode(bc(i)%dynflt)
        iend = istart + n - 1
        allocate(node1_tmp(n))
        allocate(node2_tmp(n))
        call bc_dynflt_node(bc(i)%dynflt, node1_tmp, node2_tmp)
        node1(istart:iend) = node1_tmp
        node2(istart:iend) = node2_tmp
        istart = iend + 1
        deallocate(node1_tmp)
        deallocate(node2_tmp)
    end select

  enddo

end subroutine

!=======================================================================
! build the global transformation matrix
!
! u' = X * u and 
! u  = Xinv * u'
! X for each split node pair is:
!
!  0.5 * [-1 1]
!        [1  1]
!
! u'(1) = 0.5*(u(2)-u(1))
! u'(2) = 0.5*(u(2)+u(1))
!
!
! Backward transform
!
!  u(1) = u'(2) - u'(1)
!  u(2) = u'(2) + u'(1)
!
!  X'   = [-1, 1; 1, 1]
!
!
! if index1 = index2 (symmetric assumption)
! only node1 are stored
! u'(1) = -u(1)
! u(1)  = -u'(1)
!

subroutine BC_build_transform_mat(bc, X, Xinv, ndof, npoin, ierr)
#include <petsc/finclude/petscksp.h>
  use petscksp
  type(bc_type), pointer :: bc(:)
  Mat:: Xinv, X
  integer:: ndof, npoin, i, j, nnode, n1, n2
  integer,dimension(:), allocatable:: node1, node2
  integer, dimension(:), allocatable:: rows, cols
  double precision, dimension(:,:), allocatable:: vals, valsinv
  PetscErrorCode  :: ierr

  ! obtain all the fault nodes
  nnode = BC_nfaultnode(bc)
  allocate(node1(nnode), node2(nnode))
  call BC_faultnode(bc, node1, node2, nnode)

  call MatCreate(PETSC_COMM_WORLD, X, ierr) 
  CHKERRQ(ierr)
! set matrix option at runtime

  ! -mat_type mpiaij
  call MatSetFromOptions(X, ierr);CHKERRQ(ierr)
  call MatSetSizes(X, PETSC_DECIDE, PETSC_DECIDE, ndof*npoin, ndof*npoin, ierr)
  CHKERRQ(ierr)

! Multiple ways of allocating memory
  call MatMPIAIJSetPreallocation(X, 2, &
       PETSC_NULL_INTEGER, 2, PETSC_NULL_INTEGER,ierr);
  call MatSeqAIJSetPreallocation(X, 2, PETSC_NULL_INTEGER,ierr);
  call MatSetUp(X, ierr); CHKERRQ(ierr)

  call MatDuplicate(X, MAT_DO_NOT_COPY_VALUES, Xinv, ierr)
  call MatSetUp(Xinv, ierr)
  
  ! set diagonals, index is 0-based
  do i = 0, npoin*ndof - 1
      call MatSetValue(X, i, i, 1d0, INSERT_VALUES, ierr)
      call MatSetValue(Xinv, i, i, 1d0, INSERT_VALUES, ierr)
  end do

!  call MatDuplicate(X, MAT_COPY_VALUES, Xinv, ierr)

  ! create X and Xinv for each fault split node pair and each ndof
  do i = 1, nnode
      n1 = node1(i)
      n2 = node2(i)

      if (n1==n2) then
          ! symmetric node 2 is not stored
          allocate(rows(1), cols(1), vals(1,1), valsinv(1,1))
      else
          allocate(rows(2), cols(2), vals(2,2), valsinv(2,2))
      end if
      
     do j = 1, ndof
          if (n1==n2) then
             rows(1) = (n1 - 1) * ndof + j - 1 ! global index for node n1 and dof j
             cols    = rows
             vals(1,1) = -1d0
             valsinv(1,1) = -1d0
          else
             rows(1) = (n1 - 1) * ndof + j - 1 ! global index for node n1 and dof j
             rows(2) = (n2 - 1) * ndof + j - 1 ! global index for node n2 and dof j
             cols    = rows
             vals(1,1) = -0.5d0
             vals(1,2) =  0.5d0
             vals(2,:) =  0.5d0
             valsinv   = 2d0*vals
          end if
          call MatSetValues(X, size(rows), rows, size(cols), cols, vals, INSERT_VALUES, ierr)
          call MatSetValues(Xinv, size(rows), rows, size(cols), cols, valsinv, INSERT_VALUES, ierr)
     end do !j
     deallocate(rows, cols, vals, valsinv)
  end do !i

  deallocate(node1, node2)

  ! assemble matrix X and Xinv
  call MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyend(X, MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyBegin(Xinv, MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyend(Xinv, MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
end subroutine

!=======================================================================
! Writes data for faults, and possibly other BCs
subroutine BC_write(bc,time,d,v)
  use time_evol, only : timescheme_type
  use echo, only : iout,echo_init, fmt1,fmtok
  type(timescheme_type), intent(in) :: time
  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(in) :: d,v

  integer :: i

  if (.not. associated(bc)) return

  do i = 1,size(bc)
    if (bc(i)%kind == IS_DYNFLT .or. bc(i)%kind == IS_KINFLT) then
      if (echo_init .and. time%it==0) write(iout,fmt1,advance='no') 'Exporting initial boundary data'
      if (bc(i)%kind == IS_DYNFLT) then
        call BC_DYNFLT_write(bc(i)%dynflt, time, d, v)
      else
        call BC_KINFLT_write(bc(i)%kinflt, time, d, v)
      endif
      if (echo_init .and. time%it==0) write(iout,fmtok)
    endif
  enddo

  ! flush the output every ntflush time steps, default 100
  if (mod(time%nt, time%ntflush)==0) call flush()

end subroutine BC_write

! =================================================
! select the fields only along the one kind of bc 
subroutine bc_select_kind(bc, field_in, field_out, bc_kind, dirneu_kind_in)

  use stdio, only: IO_abort
  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(in) :: field_in
  double precision, dimension(:,:), intent(inout) :: field_out
  integer :: i
  integer, intent(in):: bc_kind
  integer, intent(in), optional:: dirneu_kind_in
  integer :: dirneu_kind
  integer :: IS_DIR    = 2, & 
             IS_NEU    = 1 

  dirneu_kind = IS_NEU 
  if (present(dirneu_kind_in)) dirneu_kind = dirneu_kind_in 
  
  if (.not. associated(bc)) return
  
  if (bc_kind /= IS_DYNFLT .or. bc_kind /= IS_DIRNEU) then
      call IO_abort('bc_select_kind: bc_kind must be IS_DYNFLT or IS_DIRNEU')
      return
  end if

  do i = 1,size(bc)
      if (bc(i)%kind /= bc_kind) continue
      select case (bc(i)%kind)
      case (IS_DYNFLT)
        call BC_DYNFLT_select(bc(i)%dynflt, field_in, field_out)
      case (IS_DIRNEU) 
        call BC_DIRNEU_select_kind(bc(i)%dirneu, field_in, field_out, dirneu_kind)
      case (IS_DIRABS) 
        call BC_DIRABS_select_kind(bc(i)%dirabs, field_in, field_out, dirneu_kind)
      end select
  enddo
end subroutine bc_select_kind

! ====================================================================
! select field on the nodes where degree of freedom is fixed
!
! first transform the field, d to d'
! select the following bc:
!       dirichlet boundaries,
!       kinematic fault boundary (side -1, node1)
!       dynflt boundary (side -1, node1) 
!
! undo the transform after the selection
!
! used only in quasi-statics
!

subroutine bc_select_fix(bc, field_in, field_out)

  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(inout) :: field_in
  double precision, dimension(:,:), intent(inout) :: field_out
  integer :: i, side
  integer :: IS_DIR    = 2 

  side = -1
!  side = 2

  ! transform field_in 
  call bc_trans(bc, field_in, 1)

  if (.not. associated(bc)) return

  do i = 1,size(bc)
      select case (bc(i)%kind)
      case (IS_DYNFLT)
          ! select dynflt side -1
        call BC_DYNFLT_select(bc(i)%dynflt, field_in, field_out, side)
      case (IS_DIRNEU) 
          ! select dirichlet 
        call BC_DIRNEU_select_kind(bc(i)%dirneu, field_in, field_out, IS_DIR)
      case (IS_DIRABS) 
          ! select dirichlet 
        call BC_DIRABS_select_kind(bc(i)%dirabs, field_in, field_out, IS_DIR)
      case (IS_KINFLT) 
        call BC_KINFLT_select(bc(i)%kinflt, field_in, field_out, side)
      end select
  enddo
  ! undo the transform
  call bc_trans(bc, field_in, -1)
  call bc_trans(bc, field_out, -1)

end subroutine bc_select_fix

! Return the index of fixed dofs
! Note the index is 1 based (fortran)
! dirichlet dofs in dirneu and dirabs
! -1 side dofs for dynflt and kinflt 
! (in the transformed problem)

subroutine bc_GetIndexDofFix(bc, indexDofFix, ndof)
  type(bc_type), pointer :: bc(:)
  integer, dimension(:), intent(inout) :: indexDofFix
  integer :: i, n, istart, ndof
  istart =1 

  do i = 1,size(bc)
      select case (bc(i)%kind)
      case (IS_DYNFLT)
        call BC_DYNFLT_AppendDofFix(bc(i)%dynflt, indexDofFix, istart, ndof) 
      case (IS_DIRNEU) 
        call BC_DIRNEU_AppendDofFix(bc(i)%dirneu, indexDofFix, istart, ndof) 
      case (IS_DIRABS) 
        call BC_DIRABS_AppendDofFix(bc(i)%dirabs, indexDofFix, istart, ndof) 
      case (IS_KINFLT) 
        call BC_KINFLT_AppendDofFix(bc(i)%kinflt, indexDofFix, istart, ndof)
      end select
  enddo

end subroutine

function bc_nDofFix(bc, ndof) result(n)
  type(bc_type), pointer :: bc(:)
  integer :: i, n, ndof

  n = 0
  do i = 1,size(bc)
      select case (bc(i)%kind)
      case (IS_DYNFLT)
        n = n + BC_DYNFLT_nDofFix(bc(i)%dynflt, ndof) 
      case (IS_DIRNEU) 
        n = n + BC_DIRNEU_nDofFix(bc(i)%dirneu, ndof) 
      case (IS_DIRABS) 
        n = n + BC_DIRABS_nDofFix(bc(i)%dirabs, ndof) 
      case (IS_KINFLT) 
        n = n + BC_KINFLT_nDofFix(bc(i)%kinflt, ndof) 
      end select
  enddo
end function 
!======================================================================
! set the field on the fixed degree of freedom to zero
!
! used only in quasi-static for static condensation
!
!       dirichlet boundaries,
!       kinematic fault boundary (side -1)
!       dynflt boundary (side -1) 
!

subroutine bc_set_fix_zero(bc, field)

  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(inout) :: field
  integer :: i, side
  integer :: IS_DIR  = 2 

  side = -1
  !side = 2 ! TEST
  ! transform
  !call bc_trans(bc, field, 1)

  if (.not. associated(bc)) return

  do i = 1,size(bc)
      select case (bc(i)%kind)
      case (IS_DYNFLT)
          ! select dynflt side -1
        call BC_DYNFLT_set(bc(i)%dynflt, field, 0d0, side)
      case (IS_DIRNEU) 
          ! select dirichlet 
        call BC_DIRNEU_set_kind(bc(i)%dirneu, field, 0d0, IS_DIR)
      case (IS_DIRABS) 
          ! select dirichlet 
        call BC_DIRABS_set_kind(bc(i)%dirabs, field, 0d0, IS_DIR)
      case (IS_KINFLT) 
        call BC_KINFLT_set(bc(i)%kinflt, field, 0d0, side)
      end select
  enddo
  !call bc_trans(bc, field, -1)

end subroutine bc_set_fix_zero


!=======================================================================
! Sets the field along the boundary to a specific value
! Only used in quasi-static solver
! use to fix values along certain boundary
!
! NOT USED
!subroutine bc_set(bc,field,input)
!
!  use fields_class, only: fields_type
!
!  type(bc_type), pointer :: bc(:)
!  double precision, intent(in) :: input
!  double precision, dimension(:,:), intent(inout) :: field
!
!  integer :: i
!
!  if (.not. associated(bc)) return
! ! apply first periodic, then absorbing, then the rest
! ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
!  do i = 1,size(bc)
!    if ( bc(i)%kind == IS_PERIOD) call bc_set_single(bc(i))
!  enddo
!  do i = 1,size(bc)
!    if ( bc(i)%kind == is_absorb) call bc_set_single(bc(i))
!  enddo
!  do i = 1,size(bc)
!    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB) call bc_set_single(bc(i))
!  enddo
!    
!contains
!
!  subroutine bc_set_single(bc)
!
!    type(bc_type), intent(inout) :: bc
!    ! DEVEL these other set functions will have to be created 
!    ! in their respective files.
!    select case(bc%kind)
!      case(IS_DIRNEU)
!        !call bc_DIRNEU_set(bc%dirneu,field_in,input,field_out)
!      case(IS_KINFLT)
!        !call bc_KINFLT_set(bc%kinflt,field_in,input,field_out)
!      case(IS_ABSORB)
!        !call BC_ABSO_set(bc%abso,field_in,input,field_out)
!      case(IS_PERIOD)
!        !call bc_perio_set(bc%perio,field_in,input,field_out)
!      case(IS_LISFLT)
!        !call BC_LSF_set(bc%lsf,field_in,input,field_out)
!      case(IS_DYNFLT)
!        call BC_DYNFLT_set(bc%dynflt,field,input)
!!!      case(IS_USER)
!!!        call BC_USER_set(bc%user,field_in,input,field_out)
!    end select
!
!  end subroutine bc_set_single
!  
!end subroutine bc_set

! set fields on the boundary for a specific kind
! Note only dirichlet and dynflt is supported
subroutine bc_set_kind(bc, field, input, bc_kind, dirneu_kind, side_in)

  use stdio, only: IO_abort
  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, intent(in) :: input
  double precision, dimension(:,:), intent(inout) :: field
  integer, intent(in):: bc_kind, dirneu_kind
  integer, intent(in), optional :: side_in 
  integer :: side

  integer :: i

  side  = 2 ! default, select all available sides
  if (present(side_in)) side = side_in

  if (.not. associated(bc)) return

  do i = 1, size(bc)
      if (bc(i)%kind /= bc_kind) continue
      select case (bc(i)%kind)
          case (IS_DYNFLT)
             call BC_DYNFLT_set(bc(i)%dynflt, field, input, side) 
          case (IS_KINFLT)
!		    Need to implement the kinematic faults
!             call BC_KINFLT_set(bc(i)%kinflt, field, input, side) 
          case (IS_DIRNEU)
             call BC_DIRNEU_set_kind(bc(i)%dirneu, field, input, dirneu_kind) 
          case default
             call IO_abort('bc_set_kind: unsupported bc kind') 
     end select
  enddo
    
end subroutine bc_set_kind

! ==============================================================
! Transform a field on the fault boundary
! direction = 1
! 
! fout(node1) = (fin(node2) - fin(node1))/2    Half jump
! fout(node2) = (fin(node2) + fin(node1))/2    Average
!
! If symmetry is assumed, then node2 is not stored
! fin(node2) = - fin(node1) 
!
! After transform:
! fout(node1) = - fin(node1)                 Half jump
!
! direction = -1 perform the inverse transformation
!
! This is used to impose slip on the fault
! for both dynflt/kinematic faults
!
! See master thesis of Junpei Seki
!

subroutine bc_trans(bc, field, direction)

  use stdio, only: IO_abort

  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(inout) :: field
  integer, intent(in) :: direction
  integer :: i

  if (direction/=1 .and. direction/= -1) then
     call IO_abort('bc_trans: direction must be 1 or -1')
  end if

  if (.not. associated(bc)) return
  
  do i = 1,size(bc)
    select case(bc(i)%kind)
    case (IS_KINFLT)
        call bc_kinflt_trans(bc(i)%kinflt, field, direction)
    case (IS_DYNFLT)
        call bc_dynflt_trans(bc(i)%dynflt, field, direction)
    end select
  enddo

end subroutine bc_trans

! ============================================
! update displacement only on the faults
! d(fault_node) = dpre(fault_node) + dt
! only used by dynflt in quasi-static mode

subroutine bc_update_dfault(bc, dpre, d, v, dt)
  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(in) :: dpre,v
  double precision, dimension(:,:), intent(inout) :: d
  double precision, intent(in) :: dt 
  integer :: i
  
  if (.not. associated(bc)) return
  
  do i = 1,size(bc)
    if (bc(i)%kind==IS_DYNFLT) then
      call bc_dynflt_update_disp(bc(i)%dynflt, dpre, d, v, dt)
    end if
  enddo

end subroutine bc_update_dfault

subroutine bc_update_bcdv(bc, d, time)
  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(in) :: d
  double precision :: time
  integer :: i
  
  if (.not. associated(bc)) return
  
  do i = 1,size(bc)
    if (bc(i)%kind==IS_DYNFLT) then
      call bc_dynflt_update_bcdv(bc(i)%dynflt, d, time)
    end if
  enddo

end subroutine bc_update_bcdv

subroutine bc_set_fault(bc, field, field_set, side_in)
  type(bc_type), pointer :: bc(:)
  double precision, dimension(:,:), intent(inout) :: field
  double precision, dimension(:,:), intent(in) :: field_set
  integer, intent(in), optional :: side_in 
  integer :: side
  integer :: i

  side  = 2 ! default, select all available sides
  if (present(side_in)) side = side_in
  
  if (.not. associated(bc)) return
  
  do i = 1,size(bc)
    if (bc(i)%kind==IS_DYNFLT) then
      call bc_dynflt_set_array(bc(i)%dynflt, field, field_set, side)
    end if
  enddo

end subroutine bc_set_fault

function bc_has_dynflt(bc) result(has_dynflt)
  type(bc_type), pointer, intent(in) :: bc(:)
  logical :: has_dynflt 
  integer :: i

  has_dynflt = .false.
  do i = 1, size(bc)
      if (bc(i)%kind==IS_DYNFLT) then
          has_dynflt=.true.
          return
      end if
  end do
end function bc_has_dynflt

subroutine bc_timestep(bc, time)

  use time_evol, only : timescheme_type
  type(bc_type), pointer, intent(in) :: bc(:)
  type(timescheme_type), intent(inout) :: time
  integer :: i

  if (.not. associated(bc)) return

  do i = 1, size(bc)
    select case(bc(i)%kind)
    case (IS_DYNFLT)
        call BC_DYNFLT_timestep(time, bc(i)%dynflt)
    end select
  enddo
    
end subroutine bc_timestep

end module bc_gen
