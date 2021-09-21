module bc_dirneu

! Dirichlet (null displacement) and/or Neumann (null traction) conditions
! on vertical or horizontal boundaries
! possibly time-dependent
  
  use bnd_grid, only : bnd_grid_type
  use stf_gen
  use time_evol, only : timescheme_type

  implicit none
  private

  type bc_dirneu_type
    private
    integer :: kind(2)
    type(bnd_grid_type), pointer :: topo =>null()
    type (stf_type), pointer  :: hstf =>null(), vstf =>null() 
    double precision, pointer :: B(:) =>null() 
  end type

  integer,parameter :: IS_NEUMANN=1, IS_DIRICHLET=2

  public :: BC_DIRNEU_type, BC_DIRNEU_read, BC_DIRNEU_init, BC_DIRNEU_apply

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DIRNEU
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Dirichlet (null displacement) 
!          and/or Neumann (null or time-dependent traction) 
!          boundary conditions on vertical or horizontal boundaries
! SYNTAX : &BC_DIRNEU h, v, hsrc, vsrc /
!          possibly followed by one or two STF_XXXX blocks
!
! ARG: h        [char]['N'] Boundary condition on the horizontal component
! ARG: v        [char]['N'] Boundary condition on the vertical component :
!                       'N' : Neumann 
!                       'D' : Dirichlet
! ARG: hsrc     [name]['none'] Name of the source time function for a
!                time-dependent horizontal traction: 
!                'RICKER', 'TAB', 'USER', etc  (see STF_XXXX input blocks)
! ARG: vsrc     [name]['none'] Same for the vertical component
!
! END INPUT BLOCK

subroutine bc_DIRNEU_read(bc,iin)
  
  use echo , only: echo_input,iout
  use stdio, only: IO_abort

  type(bc_DIRNEU_type), pointer :: bc
  integer, intent(in) :: iin

  character(1) :: h,v
  character(15) :: hstf,vstf
  character(10) :: htype,vtype

  NAMELIST / BC_DIRNEU /  h,v,hstf,vstf

  h = 'N'
  v = 'N'
  hstf = 'none'
  vstf = 'none'

  read(iin,BC_DIRNEU,END=100)

  allocate(bc)

  if (h=='N') then
    bc%kind(1) = IS_NEUMANN
    htype = 'Neumann'
  elseif (h=='D') then
    bc%kind(1) = IS_DIRICHLET
    htype = 'Dirichlet'
    hstf ='none' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRNEU_read: h must be N or D')
  endif
  if (echo_input) write(iout,200) htype,hstf
  if (hstf/='none') then
    allocate(bc%hstf)
    call STF_read(hstf,bc%hstf,iin)
  endif

  if (v=='N') then
    bc%kind(2) = IS_NEUMANN
    vtype = 'Neumann'
  elseif (v=='D') then
    bc%kind(2) = IS_DIRICHLET
    vtype = 'Dirichlet'
    vstf ='none' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRNEU_read: h must be N or D')
  endif
  if (echo_input) write(iout,300) vtype,vstf
  if (vstf/='none') then
    allocate(bc%vstf)
    call STF_read(vstf,bc%vstf,iin)
  endif

  return
  100 call IO_abort('bc_DIRNEU_read: no BC_DIRNEU block found')
  200 format(5x,'Horizontal component. . . . . . . . . (h) = ',A, &
            /5x,'         source time function . . .(hstf) = ',A)
  300 format(5x,'Vertical component. . . . . . . . . . (v) = ',A, &
            /5x,'         source time function . . .(vstf) = ',A)

end subroutine bc_DIRNEU_read


!=======================================================================
!
subroutine bc_DIRNEU_init(bc,tag,grid,perio)

  use spec_grid, only : sem_grid_type,BC_inquire, BC_get_normal_and_weights
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects
  use constants, only : TINY_XABS
  use stdio, only: IO_abort

  type(bc_DIRNEU_type) , intent(inout) :: bc
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
    call IO_abort('BC_DIRNEU_init: boundary is not vertical or horizontal')
  deallocate(n)

end subroutine bc_DIRNEU_init


!=======================================================================
!
subroutine bc_DIRNEU_apply(bc,field,time)

  type(bc_DIRNEU_type), intent(in) :: bc
  type(timescheme_type), intent(in) :: time
  double precision, intent(inout) :: field(:,:)


  if (bc%kind(1)==IS_DIRICHLET) then
    field(bc%topo%node,1) = 0.d0
  elseif (associated(bc%hstf)) then
    field(bc%topo%node,1) = field(bc%topo%node,1) + STF_get(bc%hstf,time%time)*bc%B
  endif

  if (size(field,2)==1) return
 
  if (bc%kind(2)==IS_DIRICHLET) then
    field(bc%topo%node,2) = 0.d0
  elseif (associated(bc%vstf)) then
    field(bc%topo%node,2) = field(bc%topo%node,2) + STF_get(bc%vstf,time%time)*bc%B
  endif

end subroutine bc_DIRNEU_apply

end module bc_DIRNEU
