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
  use bc_absorb

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
! PURPOSE: Dirichlet (null displacement or time-dependent displacement) 
!          and/or Neumann (null or time-dependent traction) 
!          boundary conditions on vertical or horizontal boundaries
! SYNTAX : &BC_DIRABS h, v, hsrc, vsrc /
!          possibly followed by one or two STF_XXXX blocks
!
! ARG: h        [char]['N'] Boundary condition on the horizontal component
! ARG: v        [char]['N'] Boundary condition on the vertical component :
!                       'N' : Neumann 
!                       'D' : Dirichlet
! ARG: hstf     [name]['none'] Name of the source time function for a
!                time-dependent horizontal traction/displacement: 
!                'RICKER', 'TAB', 'USER', etc  (see STF_XXXX input blocks)
! ARG: vstf     [name]['none'] Same for the vertical component
!
! END INPUT BLOCK

subroutine bc_DIRABS_read(bc,iin)
  
  use echo , only: echo_input,iout
  use stdio, only: IO_abort

  type(bc_DIRABS_type), pointer :: bc
  integer, intent(in) :: iin

  character(1) :: h,v
  character(15) :: hstf,vstf
  character(10) :: htype,vtype

  NAMELIST / BC_DIRABS /  h,v,hstf,vstf

  h = 'N'
  v = 'N'
  hstf = 'EMPTY'
  vstf = 'EMPTY'

  read(iin,BC_DIRABS,END=100)

  allocate(bc)

  if (h=='N') then
    bc%kind(1) = IS_NEUMANN
    htype = 'Neumann'
  elseif (h=='D') then
    bc%kind(1) = IS_DIRICHLET
    htype = 'Dirichlet'
!    hstf ='none' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRABS_read: h must be N or D')
  endif
  if (echo_input) write(iout,200) htype,hstf

  ! Only allocate bc%hstf when hstf is not none
  allocate(bc%hstf)
  call STF_read(hstf,bc%hstf,iin)

  if (v=='N') then
    bc%kind(2) = IS_NEUMANN
    vtype = 'Neumann'
  elseif (v=='D') then
    bc%kind(2) = IS_DIRICHLET
    vtype = 'Dirichlet'
!    vstf ='none' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRABS_read: h must be N or D')
  endif
  if (echo_input) write(iout,300) vtype,vstf

  ! Only allocate bc%vstf when vstf is not none
  allocate(bc%vstf)
  call STF_read(vstf,bc%vstf,iin)

  return
  100 call IO_abort('bc_DIRABS_read: no BC_DIRABS block found')
  200 format(5x,'Horizontal component. . . . . . . . . (h) = ',A, &
            /5x,'         source time function . . .(hstf) = ',A)
  300 format(5x,'Vertical component. . . . . . . . . . (v) = ',A, &
            /5x,'         source time function . . .(vstf) = ',A)

end subroutine bc_DIRABS_read


!=======================================================================
!
subroutine bc_DIRABS_init(bc,tag,grid,perio)

  use spec_grid, only : sem_grid_type,BC_inquire, BC_get_normal_and_weights
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects
  use constants, only : TINY_XABS
  use stdio, only: IO_abort

  type(bc_DIRABS_type) , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: tag
  type(bc_periodic_type), pointer :: perio
  
  double precision, allocatable :: n(:,:)

  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  allocate( bc%B(bc%topo%npoin), n(bc%topo%npoin,2) )
  call BC_get_normal_and_weights(bc%topo, grid, n, bc%B, &
                                 BC_PERIO_intersects(bc%topo,perio) ) 
  if (.not. (associated(bc%hstf).or.associated(bc%vstf) )) deallocate(bc%B)

 ! check that the boundary is flat, vertical or horizontal:
  if ( .not. (all(abs(n(:,1))<TINY_XABS).or.all(abs(n(:,2))<TINY_XABS) )) &
    call IO_abort('BC_DIRABS_init: boundary is not vertical or horizontal')
  deallocate(n)

end subroutine bc_DIRABS_init

!=======================================================================
! Apply the Dirichlet or Neumann boundary condition 
!
! Dirichlet boundary condition: 
! 
! Direct overwrite boundary data and nullify the force
!
! set the force = 0, which then ==> acceleration = 0. 
! set     disp  = STF(time)
!
! For Neumann boundary condition:
! set force += B * T 
! 
! where T= STF(time)
!
! Modified by Chao, 08/12/2020
!

subroutine bc_DIRABS_apply(bc, disp, force,time)

  type(bc_DIRABS_type), intent(in) :: bc
  type(timescheme_type), intent(in) :: time
  double precision, intent(inout) :: disp(:,:)
  double precision, intent(inout) :: force(:,:)

  select case(bc%kind(1))
      case(IS_DIRICHLET)
          ! Dirichle
          force(bc%topo%node,1) = 0.d0
          if (associated(bc%hstf)) then
            disp(bc%topo%node,1) = STF_get(bc%hstf,time%time)
          endif
      case(IS_NEUMANN)
          ! Neumann
          force(bc%topo%node,1) = force(bc%topo%node,1) + STF_get(bc%hstf,time%time)*bc%B
  end select

  if (size(force,2)==1) return

  select case(bc%kind(2))
      case(IS_DIRICHLET)
          ! Dirichle
          force(bc%topo%node,2) = 0.d0
          if (associated(bc%vstf)) then
            disp(bc%topo%node,2) = STF_get(bc%vstf,time%time)
          endif
      case(IS_NEUMANN)
          ! Neumann
          force(bc%topo%node,2) = force(bc%topo%node,2) + STF_get(bc%vstf,time%time)*bc%B
  end select
 
end subroutine bc_DIRABS_apply

! only apply dirichlet or neumann boundary condition one at a time
! bc_kind = IS_NEUMANN=1 or IS_DIRICHLET=2
!
subroutine bc_DIRABS_apply_kind(bc, disp, force,time, bc_kind)

  type(bc_DIRABS_type), intent(in) :: bc
  type(timescheme_type), intent(in) :: time
  integer, intent(in) :: bc_kind
  double precision, intent(inout) :: disp(:,:)
  double precision, intent(inout) :: force(:,:)
  if (bc%kind(1) == bc_kind) then
      select case(bc%kind(1))
          case(IS_DIRICHLET)
              ! Dirichle
              force(bc%topo%node,1) = 0.d0
              if (associated(bc%hstf)) then
                disp(bc%topo%node,1) = STF_get(bc%hstf,time%time)
              endif
          case(IS_NEUMANN)
              ! Neumann
              force(bc%topo%node,1) = force(bc%topo%node,1) + STF_get(bc%hstf,time%time)*bc%B
      end select
  end if

  if (size(force,2)==1) return

  if (bc%kind(2) == bc_kind) then
      select case(bc%kind(2))
          case(IS_DIRICHLET)
              ! Dirichle
              force(bc%topo%node,2) = 0.d0
              if (associated(bc%vstf)) then
                disp(bc%topo%node,2) = STF_get(bc%vstf,time%time)
              endif
          case(IS_NEUMANN)
              ! Neumann
              force(bc%topo%node,2) = force(bc%topo%node,2) + STF_get(bc%vstf,time%time)*bc%B
      end select
  end if
 
end subroutine bc_DIRABS_apply_kind

! only apply dirichlet or neumann boundary condition one at a time
! bc_kind = IS_NEUMANN=1 or IS_DIRICHLET=2
!

subroutine bc_DIRABS_set_kind(bc, field, input, bc_kind)

  type(bc_DIRABS_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)
  double precision, intent(in) :: input

  integer, intent(in) :: bc_kind
  
  if (bc%kind(1) == bc_kind) then
      field(bc%topo%node, 1) = input
  end if

  if (size(field,2)==1) return

  if (bc%kind(2) == bc_kind) then
	  field(bc%topo%node, 2) = input
  end if

end subroutine bc_DIRABS_set_kind

subroutine bc_DIRABS_select_kind(bc, field_in, field_out, bc_kind)

  type(bc_DIRABS_type), intent(in) :: bc
  double precision, intent(in) :: field_in(:,:)
  double precision, dimension(:,:), intent(inout) :: field_out

  integer, intent(in) :: bc_kind
  
  if (bc%kind(1) == bc_kind) then
	  field_out(bc%topo%node, 1) = field_in(bc%topo%node, 1)
  end if

  if (size(field_in,2)==1) return

  if (bc%kind(2) == bc_kind) then
	  field_out(bc%topo%node, 2) = field_in(bc%topo%node, 2)
  end if

end subroutine bc_DIRABS_select_kind

end module bc_DIRABS
