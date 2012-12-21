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
module sources

  use src_explo
  use src_force
  use src_wave

  use ricker_functions
  use butterworth_filter

  implicit none
  private

  integer, parameter :: ndime = 2

  type src_timef_type
    private
    integer :: kind = 0
    type (ricker_type), pointer :: ricker => null()
    type (butter_type), pointer :: butter => null()
  end type src_timef_type

  type src_mechanism_type
    private
    integer :: kind = 0
    type (so_explo_type), pointer :: explo => null()
    type (so_force_type), pointer :: force => null()
    type (src_wave_type), pointer :: wave  => null()
  end type src_mechanism_type

  type source_type
    private
    double precision :: coord(ndime) = 0.d0
    integer :: relocate_to = 0
    type (src_timef_type) :: time
    type (src_mechanism_type) :: mech
  end type source_type

  integer, parameter :: tag_explo  = 1 &
                       ,tag_force  = 2 &
                       ,tag_wave   = 3

  integer, parameter :: tag_ricker = 1 &
                       ,tag_butter = 2

  public :: source_type,SO_check,SO_read,SO_init,SO_add &
           ,SO_WAVE_veloc,SO_WAVE_accel,SO_inquire

contains


!=====================================================================
!
  subroutine SO_inquire(src,coord,is_wave)

  type(source_type), pointer :: src
  double precision, intent(out), optional :: coord(ndime)
  logical, intent(out), optional :: is_wave

  if ( present(coord) )  coord = src%coord
  if ( present(is_wave) ) is_wave = src%mech%kind == tag_wave
  
  end subroutine SO_inquire


!=====================================================================
!
! Read and store source functions

! BEGIN INPUT BLOCK
!
! NAME   : SRC_DEF
! PURPOSE: Define the sources.
! SYNTAX : &SRC_DEF   TimeFunction,mechanism,coord,relocate_to /
!          Immediately followed by SRC_TIMEFUNCTION then SRC_MECHANISM blocks
!
! ARG: TimeFunction    [name] [none]   The name of the source time function:
!                        'RICKER' or 'BUTTERWORTH'
! ARG: mechanism        [name] [none]   The name of the source mechanism:
!                        'EXPLOSION','WAVE' or 'FORCE'
! ARG: coord            [dble] [huge]   Location of the source (m). 
! ARG: relocate_to      [int] [0]       The coordinates can be relocated (projected)
!                        to a boundary if you give here one of the "tags" read in
!                        the DEF_BC input blocks. Usually 1 for the free surface.

! END INPUT BLOCK

  subroutine SO_read(so,iin)

  use stdio, only : IO_abort
  use echo , only : echo_input,iout

  integer, intent(in) :: iin
  type(source_type), pointer :: so

  double precision :: coord(ndime),init_double
  integer :: relocate_to
  character(15) :: TimeFunction,mechanism

  NAMELIST / SRC_DEF / TimeFunction,mechanism,coord,relocate_to

  init_double = huge(init_double)

!-- read general parameters

  TimeFunction = ' '
  mechanism     = ' '
  coord         = init_double
  relocate_to   = 0

  rewind(iin)
  read(iin,SRC_DEF,END=200)

  if (TimeFunction == ' ')  &
    call IO_abort('SO_read: You must set the "TimeFunction" ')
  if (mechanism == ' ') &
    call IO_abort('SO_read: You must set the "mechanism" ')
  if (any(coord == init_double)) &
    call IO_abort('SO_read: You must set "coord" ')
  if (relocate_to /= 0 .and. mechanism == 'EXPLOSION') &
    call IO_abort('SO_read: cannot relocate explosion on a boundary')

  if(echo_input) then
    write(iout,100) coord
    if (relocate_to /= 0) write(iout,110) relocate_to
  endif

!-- store data 
  allocate(so)
  so%relocate_to = relocate_to
  so%coord       = coord 
  nullify(so%time%butter,so%time%ricker)
  
!-- read time function parameters
  select case (TimeFunction)
    case('BUTTERWORTH')
      so%time%kind = tag_butter
      allocate(so%time%butter)
      call BUTTER_read(so%time%butter,iin)
    case('RICKER')
      so%time%kind = tag_ricker
      allocate(so%time%ricker)
      call RICKER_read(so%time%ricker,iin)
    case default
      call IO_abort('SO_read: unknown TimeFunction ')
  end select

!-- read mechanism
  select case (mechanism)
    case('EXPLOSION')
      so%mech%kind = tag_explo
      allocate(so%mech%explo)
      call EXPLO_read()
    case('FORCE')
      so%mech%kind = tag_force
      allocate(so%mech%force)
      call FORCE_read(so%mech%force,iin)
    case('WAVE')
      so%mech%kind = tag_wave
      allocate(so%mech%wave)
      call WAVE_read(so%mech%wave,iin)
    case default
      call IO_abort('SO_read: unknown Mechanism ')
  end select

 return

!---- formats and escapes

  100   format(//,' S o u r c e  F u n c t i o n s',/1x,28('='),//5x, &
     'X-position (meters). . . . . . . . . . =',EN12.3,/5x, &
     'Y-position (meters). . . . . . . . . . =',EN12.3)
  110   format(5x, &
     'Relocated on boundary tagged by. . . . =',I0)
  200   nullify(so) ; return

         
end subroutine SO_read

!=====================================================================
!
  subroutine SO_init(so,grid,KD,fields,elast)

  use spec_grid, only : sem_grid_type, SE_find_nearest_node
  use fields_class, only : fields_type
  use elastic, only: elast_type

  type(source_type), intent(inout) :: so
  type(sem_grid_type), intent(in)  :: grid
  type(fields_type), intent(inout) :: fields
  type(elast_type), intent(in) :: elast
  double precision, pointer :: KD(:,:)
  
  double precision :: coord(ndime)
  integer :: iglob,i,j,e

  if (so%relocate_to /= 0) then
    !-- find the nearest node on a BC
    call SE_find_nearest_node(so%coord,grid,iglob,i,j,e,coord &
                             ,bound = so%relocate_to)
  else
    call SE_find_nearest_node(so%coord,grid,iglob,i,j,e,coord &
                           ,inner = so%mech%kind == tag_explo)
    !-- if explosion source the node must be inside the e
  endif
   
  select case (so%mech%kind)
    case(tag_explo)
      call EXPLO_init(so%mech%explo,i,j,e,grid,KD)
    case(tag_force)
      call FORCE_init(so%mech%force,KD,iglob)
    case(tag_wave)
      call WAVE_init(so%mech%wave,so%time%ricker,fields,grid,elast &
                   ,i,j,e,so%coord)
  end select

  so%coord = coord
  
  end subroutine SO_init

!=====================================================================
!
  subroutine SO_check(so,dt,nt) 

  use stdio, only : IO_new_unit,IO_abort
  
  type(source_type), intent(inout) :: so
  double precision, intent(in) :: dt
  integer, intent(in) :: nt
 
  integer :: it,fid
  double precision :: ampli

  fid = IO_new_unit()  
  open(unit=fid,file='SourcesTime_sem2d.tab',status='unknown')
  do it=1,nt
    write(fid,*) real(it*dt),real( SO_get_F(so,it*dt) )
  enddo
  close(fid)

  end subroutine SO_check


!=====================================================================
!
  subroutine SO_add(so,t)

  type(source_type), intent(inout) :: so
  double precision, intent(in) :: t
  
  double precision :: ampli

  ampli = SO_get_F(so,t)
  call SO_add_F(so,ampli)

  end subroutine SO_add


!=====================================================================
!
  subroutine SO_add_F(so,ampli)

  type (source_type), intent(inout) :: so
  double precision, intent(in) :: ampli

  select case (so%mech%kind)
    case(tag_explo)
      call EXPLO_add(so%mech%explo,ampli)
    case(tag_force)
      call FORCE_add(so%mech%force,ampli)
  end select

  end subroutine SO_add_F


!=====================================================================
! compute amplitude of source time function
!
  double precision function SO_get_F(so,t)

  type(source_type), intent(in) :: so
  double precision , intent(in) :: t

  select case (so%time%kind)
    case(tag_ricker)
      SO_get_F = ricker(t,so%time%ricker)
    case(tag_butter)
      SO_get_F = butter(so%time%butter,t)
  end select

  end function SO_get_F


!=====================================================================
! For incident plane wave

function SO_wave_veloc(time,coord,src)

  type (source_type), intent(in) :: src
  double precision, intent(in) :: coord(ndime),time
  double precision :: SO_wave_veloc(ndime)

  SO_wave_veloc = WAVE_veloc(time,coord,src%mech%wave,src%time%ricker)

end function SO_wave_veloc

!-----------------------------------------------------------------------------

function SO_wave_accel(time,coord,src)

  type (source_type), intent(in) :: src
  double precision, intent(in) :: coord(ndime),time
  double precision :: SO_wave_accel(ndime)

  SO_wave_accel = WAVE_accel(time,coord,src%mech%wave,src%time%ricker)
  
end function SO_wave_accel
  

endmodule sources
