! SEM2DPACK version 2.3.5 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bc_dirneu

! Dirichlet (null displacement) and/or Neumann (null traction) conditions
! on vertical or horizontal boundaries
! possibly time-dependent
  
  use bnd_grid, only : bnd_grid_type
  use stf_gen

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

  public :: BC_DIRNEU_type, BC_DIRNEU_read, BC_DIRNEU_init, BC_DIRNEU_set

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
! ARG: hsrc     [name]['null'] Name of the source time function for a
!                time-dependent horizontal traction: 
!                'RICKER', 'TAB', 'USER', etc  (see STF_XXXX input blocks)
! ARG: vsrc     [name]['null'] Same for the vertical component
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
  hstf = 'null'
  vstf = 'null'

  read(iin,BC_DIRNEU,END=100)

  allocate(bc)

  if (h=='N') then
    bc%kind(1) = IS_NEUMANN
    htype = 'Neumann'
  elseif (h=='D') then
    bc%kind(1) = IS_DIRICHLET
    htype = 'Dirichlet'
    hstf ='null' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRNEU_read: h must be N or D')
  endif
  if (echo_input) write(iout,200) htype,hstf
  if (hstf/='null') then
    allocate(bc%hstf)
    call STF_read(hstf,bc%hstf,iin)
  endif

  if (v=='N') then
    bc%kind(2) = IS_NEUMANN
    vtype = 'Neumann'
  elseif (v=='D') then
    bc%kind(2) = IS_DIRICHLET
    vtype = 'Dirichlet'
    vstf ='null' ! time-dependent only for Neumann
  else
    call IO_abort('bc_DIRNEU_read: h must be N or D')
  endif
  if (echo_input) write(iout,300) vtype,vstf
  if (vstf/='null') then
    allocate(bc%vstf)
    call STF_read(vstf,bc%vstf,iin)
  endif

  return
  100 call IO_abort('bc_DIRNEU_read: no BC_DIRNEU block found')
  200 format(5x,'Horizontal component. . . . . . . . . (h) = ',A, &
            /5x,'           source time function . .(hstf) = ',A)
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
subroutine bc_DIRNEU_set(bc,field,t)

  type(bc_DIRNEU_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)
  double precision, intent(in) :: t

  if (bc%kind(1)==IS_DIRICHLET) then
    field(bc%topo%node,1) = 0.d0
  elseif (associated(bc%hstf)) then
    field(bc%topo%node,1) = field(bc%topo%node,1) + STF_get(bc%hstf,t)*bc%B
  endif

  if (size(field,1)==1) return
 
  if (bc%kind(2)==IS_DIRICHLET) then
    field(bc%topo%node,2) = 0.d0
  elseif (associated(bc%vstf)) then
    field(bc%topo%node,2) = field(bc%topo%node,2) + STF_get(bc%vstf,t)*bc%B
  endif

end subroutine bc_DIRNEU_set

end module bc_DIRNEU
