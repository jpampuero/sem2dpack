module bc_dirabs

! Dirichlet (displacement) and/or absorbing conditions
! on vertical or horizontal boundaries
! depending on time stepping scheme.
! if isdynamic, then absorbing bc is applied
! if not isdynamic, then dirichlet bc is applied.
  
  use bnd_grid, only : bnd_grid_type
  use stf_gen
  use time_evol, only : timescheme_type
  use bc_dirneu
  use bc_abso
  use spec_grid
  use memory_info
  use sources

  implicit none
  private

  type bc_dirabs_type
    private
    type(bc_abso_type), pointer  :: bc_abso
    type(bc_dirneu_type), pointer  :: bc_dir
  end type

  public :: BC_DIRABS_type, BC_DIRABS_read, BC_DIRABS_init, & 
      BC_DIRABS_apply, BC_DIRABS_apply_kind, BC_DIRABS_set_kind, BC_DIRABS_select_kind

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DIRABS
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Dirichlet during static time stepping 
!          and absorbing boundary condition during dynamic
!          boundary conditions on vertical or horizontal boundaries
!
! SYNTAX : 
!          &BC_DIRNEU h, v, hsrc, vsrc /
!          possibly followed by one or two STF_XXXX blocks
!          &BC_ABSORB stacey, let_wave / 
!
! See BC_DIRNEU and BC_ABSORB for explanation of the input parameters
!
! END INPUT BLOCK

subroutine bc_DIRABS_read(bc,iin)
  
  use echo , only: echo_input,iout
  use stdio, only: IO_abort

  type(bc_DIRABS_type), pointer :: bc
  integer, intent(in) :: iin

  allocate(bc%bc_dir)
  allocate(bc%bc_abso)
  
  ! read dirneu boundary condition
  call bc_DIRNEU_read(bc%bc_dir,iin) 

  ! read absorbing boundary condition
  call bc_ABSO_read(bc%bc_abso,iin) 
  
  if (echo_input) write(iout, 200)

  200 format(5x,'Done reading DIRNEU and ABSORB mixed boundary.', &
            /5x)

end subroutine bc_DIRABS_read

!=======================================================================
subroutine BC_DIRABS_init(bc,tag,grid,mat,M,tim,src,perio)

  use spec_grid, only : sem_grid_type
  use prop_mat, only : matpro_elem_type
  use time_evol, only: timescheme_type 
  use bc_periodic, only : bc_periodic_type
  use stdio, only: IO_abort

  type(bc_DIRABS_type) , intent(inout) :: bc
  type(sem_grid_type), intent(in)    :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  integer            , intent(in)    :: tag
  type(timescheme_type), intent(in) :: tim
  double precision, intent(inout) :: M(:,:)
  type(source_type), pointer :: src(:)
  type(bc_periodic_type), pointer :: perio

  call bc_DIRNEU_init(bc%bc_dir,tag,grid,perio)
  call bc_ABSO_init(bc%bc_abso,tag,grid,mat,M,tim,src,perio)
end subroutine bc_DIRABS_init

!=======================================================================
! Apply the Dirichlet or Neumann boundary condition at quasi-static time
! Apply absorbing boundary condition at dynamic time
!
subroutine bc_DIRABS_apply(bc, D,V,MxA,time)
  type(timescheme_type), intent(in) :: time
  double precision, intent(inout) :: MxA(:,:)
  double precision, intent(inout) :: D(:,:)
  double precision, intent(in) :: V(:,:)

  type(bc_DIRABS_type), intent(in) :: bc
  if (time%isDynamic) then
      write(*,*) "BC_DIRABS: Apply absorbing bc"
      call bc_ABSO_apply(bc%bc_abso, D, V, MxA, time)
  else
      write(*,*) "BC_DIRABS: Apply dirichlet bc"
      call bc_DIRNEU_apply(bc%bc_dir, D, MxA, time)
  end if
 
end subroutine bc_DIRABS_apply

! only apply dirichlet or neumann boundary condition one at a time
! bc_kind = IS_NEUMANN=1 or IS_DIRICHLET=2
!

! apply/set by kinds are only used for quasi-static dirneu bcs
subroutine bc_DIRABS_apply_kind(bc, disp, force, time, bc_kind)

  type(bc_DIRABS_type), intent(in) :: bc
  type(timescheme_type), intent(in) :: time
  integer, intent(in) :: bc_kind
  double precision, intent(inout) :: disp(:,:)
  double precision, intent(inout) :: force(:,:)
      
  write(*,*) "BC_DIRABS: Apply dirichlet bc, kind ", bc_kind
  call  bc_DIRNEU_apply_kind(bc%bc_dir, disp, force, time, bc_kind)
 
end subroutine bc_DIRABS_apply_kind

subroutine bc_DIRABS_set_kind(bc, field, input, bc_kind)

  type(bc_DIRABS_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)
  double precision, intent(in) :: input
  integer, intent(in) :: bc_kind
  
  call bc_DIRNEU_set_kind(bc%bc_dir, field, input, bc_kind)
  
end subroutine bc_DIRABS_set_kind

subroutine bc_DIRABS_select_kind(bc, field_in, field_out, bc_kind)
  type(bc_DIRABS_type), intent(in) :: bc
  double precision, intent(in) :: field_in(:,:)
  double precision, dimension(:,:), intent(inout) :: field_out
  integer, intent(in) :: bc_kind
  
  call bc_DIRNEU_select_kind(bc%bc_dir, field_in, field_out, bc_kind)
end subroutine bc_DIRABS_select_kind

end module bc_DIRABS
