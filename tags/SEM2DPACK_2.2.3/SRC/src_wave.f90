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
module src_wave

! Incident plane wave
! The velocity time function is a ricker

  use ricker_functions

  implicit none 
  private

  integer, parameter :: ndime = 2

  type src_wave_type
    private
    double precision,dimension(ndime) :: pol_vect,k,ref_coord 
    double precision :: c
    character :: phase
  end type src_wave_type

  public :: src_wave_type,WAVE_read,WAVE_init,WAVE_veloc,WAVE_accel

contains

!----------------------------------------------------------------------
!
! BEGIN INPUT BLOCK
!
! NAME   : SRC_WAVE [source mechanism]
! PURPOSE: Incident plane wave
! SYNTAX : &SRC_WAVE angle /
!
! ARG: angle    [dble] [90d0]   Incidence angle, in degrees, counterclockwise with 
!                               respect to X+ : 45 comes from bottom left, 
!                               135 comes from bottom right
! ARG: phase    [char] ['S']    Phase: 'S' or 'P'
!
! END INPUT BLOCK

!
subroutine WAVE_read(src,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_abort

  type (src_wave_type), intent(out) :: src
  integer, intent(in) :: iin

  double precision, parameter :: pi = 3.141592653589793d0
  double precision :: angle,angle_rad,cosa,sina
  character :: phase  
  
  NAMELIST / SRC_WAVE / angle,phase

  angle = 90d0
  phase = 'S'

  read(iin,SRC_WAVE,END=200)

  if (phase == ' ') call IO_abort('WAVE_read: Phase is a mandatory parameter')

  angle_rad = pi*angle/180.d0
  cosa = cos(angle_rad)
  sina = sin(angle_rad)
  src%k = (/ cosa ,sina /)

  select case (phase)
    case('P'); src%pol_vect = (/ cosa ,sina /)
    case('S'); src%pol_vect = (/ -sina,cosa /)
    case default
      call IO_abort('WAVE_read: Phase must be S or P') 
  end select

  src%phase = phase
  
  if (echo_input) write(iout,100) angle,phase

  return

 100 format(//,5x, &
     'Source Type. . . . . . . . . . . . . . = Plane Wave',/5x, &
     'Incidence angle (deg/horizon)  . . . . =',F0.1,/5x, &
     'Polarization . . . . . . . . . . . . . =',1x,A,/5x)
 200 call IO_abort('WAVE_read: SRC_WAVE input block not found')

end subroutine WAVE_read


!-------------------------------------------------------------------
! Initializes the incident wave
!
subroutine WAVE_init(wave,timef,fields,grid,elast,isrc,jsrc,esrc,ref_coord)

  use echo, only: echo_init,iout
  use ricker_functions, only : ricker_type,ricker_int,ricker,ricker_deriv
  use fields_class, only : fields_type
  use spec_grid, only : sem_grid_type
  use elastic, only : elast_type,ELAST_inquire

  type (src_wave_type), intent(inout) :: wave
  type (ricker_type)  , intent(in)    :: timef
  type (fields_type)  , intent(out)   :: fields
  type (sem_grid_type), intent(in)    :: grid
  type (elast_type)   , intent(in)    :: elast
  double precision    , intent(in)    :: ref_coord(grid%ndime) ! coord of reference point
  integer             , intent(in)    :: isrc,jsrc,esrc ! nearest node to the reference point
  double precision :: cp,cs,phase
  integer :: ip

  call ELAST_inquire(elast,isrc,jsrc,esrc,cp=cp,cs=cs)
  
  select case (wave%phase)
    case('P'); wave%c = cp
    case('S'); wave%c = cs
  end select
  
  wave%k         = wave%k/wave%c
  wave%ref_coord = ref_coord

  if (echo_init) then
    write(iout,'(2X,A,EN12.3,A)') 'Incident wave phase velocity = ',wave%c,' m/s' 
    write(iout,'(2X,A,2EN12.3)') 'Reference point = ',wave%ref_coord
  endif

  do ip=1,grid%npoin
    phase = - DOT_PRODUCT( grid%coord(:,ip) - wave%ref_coord , wave%k )
    fields%displ(ip,:) = wave%pol_vect * ricker_int(phase,timef)
    fields%veloc(ip,:) = wave%pol_vect * ricker(phase,timef)
    fields%accel(ip,:) = wave%pol_vect * ricker_deriv(phase,timef)
  enddo
  
end subroutine WAVE_init


!----------------------------------------------------------------------------

function wave_veloc(time,coord,wave,timef)

  type (src_wave_type), intent(in) :: wave
  type (ricker_type), intent(in) :: timef
  double precision, intent(in) :: coord(ndime),time
  double precision :: wave_veloc(ndime)

  double precision :: phase

  phase = time - DOT_PRODUCT( coord(:) - wave%ref_coord, wave%k ) 
  wave_veloc = wave%pol_vect *ricker(phase,timef)  

end function wave_veloc

!-----------------------------------------------------------------------------

function wave_accel(time,coord,wave,timef)

  type (src_wave_type), intent(in) :: wave
  type (ricker_type), intent(in) :: timef
  double precision, intent(in) :: coord(ndime),time
  double precision :: wave_accel(ndime)

  double precision :: phase

  phase = time - DOT_PRODUCT( coord(:) - wave%ref_coord , wave%k ) 
  wave_accel = wave%pol_vect *ricker_deriv(phase,timef)  

end function wave_accel


end module src_wave
