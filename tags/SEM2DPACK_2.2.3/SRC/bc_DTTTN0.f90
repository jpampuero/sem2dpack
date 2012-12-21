! SEM2DPACK version 2.2.3 -- A Spectral Element Method tool for 2D wave propagation
!                            and earthquake source dynamics
! 
! Copyright (C) 2003 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics
! ETH Hönggerberg (HPP)
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 1 633 2197 (office)
! +41 1 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
! 
! This software is freely available for scientific research purposes. 
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
module bc_DTTTN0

! Applies a horizontal displacement rate history on a boundary
! combined with free vertical traction.

! The original purpose is to apply a kinematic slip model on an inplane
! fault at a flat horizontal boundary of the model (top or bottom).

! The kinematic model is a ramp-like velocity time function defined by 
!   rupture time
!   final slip
!   rise time
! given at control nodes.
! This time function is smoothed by a cosine tapper of
! characteristic time scale "t_smooth"

! WARNING: the boundary is on bottom of model
 
  use spec_grid, only : bc_topo_type
  use distribution_general

  implicit none
  private

  type input_type
    type(distribution_type) :: rupture_time,slip,rise_time
  end type

  type bc_DTTTN0_type
    private
    double precision, dimension(:), pointer :: rupture_time &
                                              ,slip &
                                              ,rise_time &
                                              ,slip_rate
    type(input_type) :: input
    type(bc_topo_type), pointer :: topo
    double precision :: t_smooth
  end type 

  public :: BC_DTTTN0_type,BC_DTTTN0_read,BC_DTTTN0_init,BC_DTTTN0_set

contains

!=======================================================================
!
subroutine bc_DTTTN0_read(bc,iunit)

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit

  type(bc_DTTTN0_type), pointer :: bc
  integer, intent(in) :: iunit

  double precision :: t_smooth
  integer :: i,iin
  character(50) :: file
  character(10) :: distribution

  NAMELIST / BC_DTTTN0 / t_smooth,file
  NAMELIST / SRC_RUPTIME / distribution
  NAMELIST / SRC_SLIP / distribution
  NAMELIST / SRC_RISTIME / distribution

  t_smooth    = 0.d0
  file        = ' '

  read(iunit,BC_DTTTN0,END=100)

  allocate(bc)
  bc%t_smooth   = t_smooth
  
  if (file /= 'HERE') then
    iin = IO_new_unit()
    open(iin,file=file)
  else
    iin = iunit
  endif

  distribution = ' '
  read(iin,SRC_RUPTIME)
  call DIST_read(bc%input%rupture_time,distribution,iin)

  rewind(iin)
  distribution = ' '
  read(iin,SRC_SLIP)
  call DIST_read(bc%input%slip,distribution,iin)

  rewind(iin)
  distribution = ' '
  read(iin,SRC_RISTIME)
  call DIST_read(bc%input%rise_time,distribution,iin)

  if (file /= 'HERE') close(iin)
 
  100 return

end subroutine bc_DTTTN0_read

!=======================================================================
!
subroutine bc_DTTTN0_init(bc,tag,grid)

  use echo, only : echo_init,iout
  use spec_grid, only : sem_grid_type,BC_inquire

  type(bc_DTTTN0_type), intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: tag
 
  double precision, allocatable :: bc_coord(:,:)

  !-- bc%topo => grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  if (echo_init) write(iout,*) "  DTTTN0 boundary nodes = ",bc%topo%npoin

  allocate( bc_coord(grid%ndime,bc%topo%npoin) )
  bc_coord = grid%coord(:,bc%topo%bulk_node)

  allocate(bc%rupture_time(bc%topo%npoin))
  call DIST_generate(bc%rupture_time,bc_coord,bc%input%rupture_time)
  call DIST_destructor(bc%input%rupture_time)

  allocate(bc%slip(bc%topo%npoin))
  call DIST_generate(bc%slip,bc_coord,bc%input%slip)
  call DIST_destructor(bc%input%slip)

  allocate(bc%rise_time(bc%topo%npoin))
  call DIST_generate(bc%rise_time,bc_coord,bc%input%rise_time)
  call DIST_destructor(bc%input%rise_time)

  where (bc%rise_time /= 0d0)
    bc%slip_rate = bc%slip/bc%rise_time
  elsewhere
    bc%slip_rate = 0d0
  end where
  
end subroutine bc_DTTTN0_init


!=======================================================================
!
subroutine bc_DTTTN0_set(bc,vfield,time)

  type(bc_DTTTN0_type), intent(in) :: bc
  double precision, intent(out) :: vfield(:,:)
  double precision, intent(in) :: time

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, pointer :: v
  double precision :: t1,t2,w
  integer :: i

  w = 0.25d0*pi/bc%t_smooth

  do i=1,bc%topo%npoin
    t1 = bc%rupture_time(i) - bc%t_smooth
    t2 = bc%rupture_time(i)+bc%rise_time(i) + bc%t_smooth
   
    if (bc%slip(i) == 0d0 .or. time <= t1 .or. time >= t2 ) then
      v = 0d0
    else
      if ( time < bc%rupture_time(i)+bc%t_smooth ) then
        v = bc%slip_rate(i)*( 1d0-cos((time-t1)*w) )
      else if (time > bc%rupture_time(i)+bc%rise_time(i)-bc%t_smooth) then
        v = bc%slip_rate(i)*( 1d0-cos((time-t2)*w) )
      else
        v = bc%slip_rate(i)
      endif
    endif

    ! WARNING: bc on bottom
    vfield(bc%topo%bulk_node(i),1) = v

  enddo

end subroutine bc_DTTTN0_set

end module bc_DTTTN0
