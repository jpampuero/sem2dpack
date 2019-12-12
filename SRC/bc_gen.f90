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
    !! type(bc_user_type) , pointer :: user
  end type bc_type

  integer, parameter :: IS_EMPTY  = 0, &
                        IS_DIRNEU = 1, &
                        IS_KINFLT = 2, &
                        IS_ABSORB = 3, &
                        IS_PERIOD = 4, &
                        IS_LISFLT = 5, &
                        IS_DYNFLT = 6
                        !! IS_USER = 7

  public :: bc_type,bc_read,bc_apply,bc_init,bc_write,bc_set,bc_timestep

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
  double precision, dimension(:,:), intent(in) :: d,v

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
        call BC_KINFLT_init(bc(i)%kinflt,bc(i)%tag(1),grid,M,tim,perio)
      case(IS_ABSORB)
        call BC_ABSO_init(bc(i)%abso,bc(i)%tag(1),grid,mat,M,tim,src,perio)
      case(IS_LISFLT)
        call BC_LSF_init(bc(i)%lsf,bc(i)%tag,grid,M,perio)
      case(IS_DYNFLT)
        call BC_DYNFLT_init(bc(i)%dynflt,bc(i)%tag,grid,M,tim,perio)
!!      case(IS_USER)
!!        call BC_USER_init(bc(i)%user,bc(i)%tag, ... )
    end select
  enddo

 ! write initial data
  call BC_write(bc,0,d,v)

end subroutine bc_init


!=======================================================================
!! Applies the boundary condition
subroutine bc_apply(bc,time,fields,field)

  use sources, only: source_type
  use fields_class, only: fields_type
  use time_evol, only : timescheme_type

  type(bc_type), pointer :: bc(:)
  type (fields_type), intent(inout) :: fields
  type(timescheme_type), intent(in) :: time
  double precision, dimension(:,:), intent(inout) :: field

  integer :: i

  if (.not. associated(bc)) return

 ! apply first periodic, then absorbing, then the rest
 ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_PERIOD) call bc_apply_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_ABSORB) call bc_apply_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB) call bc_apply_single(bc(i))
  enddo
    
contains

  subroutine bc_apply_single(bc)

    type(bc_type), intent(inout) :: bc

    select case(bc%kind)
      case(IS_DIRNEU)
        call bc_DIRNEU_apply(bc%dirneu,field,time)
      case(IS_KINFLT)
        call bc_KINFLT_apply(bc%kinflt,fields%accel,fields%veloc,time)
      case(IS_ABSORB)
        call BC_ABSO_apply(bc%abso,fields%displ_alpha,fields%veloc_alpha,fields%accel,time)
      case(IS_PERIOD)
        call bc_perio_apply(bc%perio,field)
      case(IS_LISFLT)
        call BC_LSF_apply(bc%lsf,fields%accel,fields%displ_alpha)
      case(IS_DYNFLT)
        call BC_DYNFLT_apply(bc%dynflt,fields%accel,fields%veloc,fields%displ,time)
!!      case(IS_USER)
!!        call BC_USER_apply(bc%user, ... )
    end select

  end subroutine bc_apply_single
  
end subroutine bc_apply


!=======================================================================
! Writes data for faults, and possibly other BCs
subroutine BC_write(bc,itime,d,v)

  use echo, only : iout,echo_init, fmt1,fmtok

  type(bc_type), pointer :: bc(:)
  integer, intent(in) :: itime
  double precision, dimension(:,:), intent(in) :: d,v

  integer :: i

  if (.not. associated(bc)) return

  do i = 1,size(bc)
    if (bc(i)%kind == IS_DYNFLT .or. bc(i)%kind == IS_KINFLT) then
      if (echo_init .and. itime==0) write(iout,fmt1,advance='no') 'Exporting initial boundary data'
      if (bc(i)%kind == IS_DYNFLT) then
        call BC_DYNFLT_write(bc(i)%dynflt,itime,d,v)
      else
        call BC_KINFLT_write(bc(i)%kinflt,itime)
      endif
      if (echo_init .and. itime==0) write(iout,fmtok)
    endif
  enddo

end subroutine BC_write

!=======================================================================
!! Sets the field along the boundary to a specific value
subroutine bc_set(bc,field_in,input,field_out)

  use sources, only: source_type
  use fields_class, only: fields_type

  type(bc_type), pointer :: bc(:)
  double precision, intent(in) :: input
  double precision, dimension(:,:), intent(in) :: field_in
  double precision, dimension(:,:), intent(out) :: field_out

  integer :: i

  if (.not. associated(bc)) return
 ! apply first periodic, then absorbing, then the rest
 ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_PERIOD) call bc_set_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind == is_absorb) call bc_set_single(bc(i))
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB) call bc_set_single(bc(i))
  enddo
    
contains

  subroutine bc_set_single(bc)

    type(bc_type), intent(inout) :: bc
    ! DEVEL these other set functions will have to be created 
    ! in their respective files.
    select case(bc%kind)
      case(IS_DIRNEU)
        !call bc_DIRNEU_set(bc%dirneu,field_in,input,field_out)
      case(IS_KINFLT)
        !call bc_KINFLT_set(bc%kinflt,field_in,input,field_out)
      case(IS_ABSORB)
        !call BC_ABSO_set(bc%abso,field_in,input,field_out)
      case(IS_PERIOD)
        !call bc_perio_set(bc%perio,field_in,input,field_out)
      case(IS_LISFLT)
        !call BC_LSF_set(bc%lsf,field_in,input,field_out)
      case(IS_DYNFLT)
        call BC_DYNFLT_set(bc%dynflt,field_in,input,field_out)
!!      case(IS_USER)
!!        call BC_USER_set(bc%user,field_in,input,field_out)
    end select

  end subroutine bc_set_single
  
end subroutine bc_set

!=======================================================================
!! Sets the field along the boundary to a specific value
subroutine bc_timestep(bc,time,hcell)

  use sources, only: source_type
  use fields_class, only: fields_type
  use time_evol, only : timescheme_type

  double precision, intent(in) :: hcell
  type(bc_type), pointer, intent(in) :: bc(:)
  type(timescheme_type), intent(inout) :: time

  integer :: i

  if (.not. associated(bc)) return
 ! apply first periodic, then absorbing, then the rest
 ! Sep 29 2006: to avoid conflict between ABSORB and DIRNEU
  do i = 1,size(bc)
    if ( bc(i)%kind == IS_PERIOD) call bc_timestep_single(bc(i),time,hcell)
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind == is_absorb) call bc_timestep_single(bc(i),time,hcell)
  enddo
  do i = 1,size(bc)
    if ( bc(i)%kind /= IS_PERIOD .and. bc(i)%kind /= IS_ABSORB) call bc_timestep_single(bc(i),time,hcell)
  enddo
    
contains

  subroutine bc_timestep_single(bc,time,hcell)

    type(bc_type), intent(in) :: bc
    double precision, intent(in) :: hcell
    type(timescheme_type), intent(inout) :: time
    ! DEVEL these other set functions will have to be created 
    ! in their respective files.
    select case(bc%kind)
      case(IS_DIRNEU,IS_ABSORB,IS_PERIOD)
        !do nothing
      case(IS_KINFLT)
        !not implemented yet
        !call bc_KINFLT_timestep(bc%kinflt,time)
      case(IS_LISFLT)
        !not implemented yet
        !call BC_LSF_timestep(bc%lisflt,time)
      case(IS_DYNFLT)
        call BC_DYNFLT_timestep(time,bc%dynflt,hcell)
!!      case(IS_USER)
!!        call BC_USER_timestep(bc%user,time)
    end select

  end subroutine bc_timestep_single
  
end subroutine bc_timestep

end module bc_gen
