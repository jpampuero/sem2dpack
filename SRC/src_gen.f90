! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                             with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics Group
! ETH Hönggerberg HPP O 13.1
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
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
module sources

  use src_moment
  use src_force
  use src_wave
  use stf_gen

  use constants, only: NDIME

  implicit none
  private

  type src_mechanism_type
    private
    integer :: kind = 0
    type (so_moment_type), pointer :: moment => null()
    type (so_force_type), pointer :: force => null()
    type (src_wave_type), pointer :: wave  => null()
  end type src_mechanism_type

  type source_type
    private
    double precision :: coord(NDIME) 
    double precision :: tdelay
    type (stf_type), pointer  :: stf
    type (src_mechanism_type) :: mech
  end type source_type

  integer, parameter :: tag_moment = 1 &
                       ,tag_force  = 2 &
                       ,tag_wave   = 3

  public :: source_type,SO_check,SO_read,SO_init,SO_add &
           ,SO_WAVE_get_VT,SO_inquire

contains


!=====================================================================
!
  subroutine SO_inquire(src,coord,is_wave)

  type(source_type) :: src
  double precision, intent(out), optional :: coord(NDIME)
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
! SYNTAX : &SRC_DEF stf, mechanism, coord, file /
!          followed by one SOURCE MECHANISM block (SRC_XXXX) 
!          and one SOURCE TIME FUNCTION block (STF_XXXX)
!
! ARG: stf        [name] [none] Name of the source time function:
!                  'RICKER', 'TAB' or 'USER'
! ARG: mechanism  [name] [none] Name of the source mechanism:
!                  'FORCE', 'EXPLOSION', 'DOUBLE_COUPLE', 'MOMENT' or 'WAVE'
! ARG: coord      [dble] [huge] Location of the source (m). 
! ARG: file       [string] ['none'] Station coordinates and delay times can
!                  be read from an ASCII file, with 3 columns per line:
!                  (1) X position (in m),
!                  (2) Z position (in m) and
!                  (3) time delay (in seconds)
!
! END INPUT BLOCK

  subroutine SO_read(so,iin,ndof)

  use stdio, only : IO_abort,IO_new_unit, IO_file_length
  use echo , only : echo_input,iout

  integer, intent(in) :: iin,ndof
  type(source_type), pointer :: so(:)

  double precision :: coord(NDIME),init_double
  integer :: i,n, iin2
  character(15) :: stf,mechanism
  character(50) :: file

  NAMELIST / SRC_DEF / stf,mechanism,coord,file

  init_double = huge(init_double)

!-- read general parameters

  stf = ' '
  mechanism    = ' '
  coord        = init_double
  file         = 'none'

  rewind(iin)
  read(iin,SRC_DEF,END=200)

  if (stf == ' ')  &
    call IO_abort('SO_read: stf not set')
  if (mechanism == ' ') &
    call IO_abort('SO_read: mechanism not set')
  if (any(coord == init_double) .and. file == 'none') &
    call IO_abort('SO_read: coord not set')

  if (echo_input) then
    if (file == 'none') then
      write(iout,100) coord
    else
      write(iout,105) file
    endif
  endif

!-- store data 
  if (file == 'none') then
    n=1
    allocate(so(1))
    so(1)%coord = coord
    so(1)%tdelay = 0d0

  else
    n = IO_file_length(file)
    allocate(so(n))
    iin2 = IO_new_unit()
    open(iin2,file=file)
    do i=1,n
      read(iin2 ,*) so(i)%coord(:), so(i)%tdelay
    enddo
    close(iin2)
  endif

  do i=1,n
  
   !-- read source time function parameters
    if (i==1) then
      allocate(so(i)%stf)
      call STF_read(stf,so(i)%stf,iin)
    else
     ! NOTE: all sources share the same source time function
      so(i)%stf => so(1)%stf
    endif

   !-- read mechanism
    if (i==1) rewind(iin)
    select case (mechanism)
    case('MOMENT','EXPLOSION','DOUBLE_COUPLE')
        so(i)%mech%kind = tag_moment
        allocate(so(i)%mech%moment)
        if (i==1) then
          call SRC_MOMENT_read(so(i)%mech%moment,iin,mechanism,ndof)
        else
          so(i)%mech%moment = so(1)%mech%moment
        endif
      case('FORCE')
        so(i)%mech%kind = tag_force
        allocate(so(i)%mech%force)
        if (i==1) then
          call FORCE_read(so(i)%mech%force,iin)
        else
          so(i)%mech%force = so(1)%mech%force
        endif
      case('WAVE')
        so(i)%mech%kind = tag_wave
        allocate(so(i)%mech%wave)
        if (i==1) then
          call WAVE_read(so(i)%mech%wave,iin)
        else
          so(i)%mech%wave = so(1)%mech%wave
        endif
      case default
        call IO_abort('SO_read: unknown mechanism ')
    end select

  enddo

 return

!---- formats and escapes

  100   format(//,' S o u r c e s',/1x,13('='),//5x, &
     'X-position (meters). . . . .(coord(1)) = ',EN12.3,/5x, &
     'Y-position (meters). . . . .(coord(2)) = ',EN12.3)
  105   format(//,' S o u r c e s',/1x,13('='),//5x, &
     'Source location file . . . . . .(file) = ',A)
  200   nullify(so) ; return

         
end subroutine SO_read

!=====================================================================
!
  subroutine SO_init(so,grid,mat)

  use spec_grid, only : sem_grid_type, SE_find_nearest_node
  use echo, only: iout,echo_init
  use prop_mat, only : matpro_elem_type

  type(source_type), intent(inout) :: so(:)
  type(sem_grid_type), intent(in)  :: grid
  type(matpro_elem_type), intent(in)  :: mat(:)
  
  double precision :: coord(NDIME),dist,distmax
  integer :: iglob,k

  distmax = 0.d0
  if (echo_init) write(iout,200)

  do k=1,size(so)

    call SE_find_nearest_node(so(k)%coord,grid,iglob,coord,dist)
    
    select case (so(k)%mech%kind)
      case(tag_moment)
        call SRC_MOMENT_init(so(k)%mech%moment,iglob,grid)
      case(tag_force)
        call FORCE_init(so(k)%mech%force,iglob)
      case(tag_wave)
        call WAVE_init(so(k)%mech%wave,grid,mat,coord)
    end select
  
    if (echo_init) then
      write(iout,150) k,so(k)%coord,coord,dist
      distmax = max(dist,distmax)
    endif
    
    so(k)%coord = coord
  
  enddo

  if (echo_init) write(iout,160) distmax

 150   format(2x,i7,5(1x,EN12.3))
 160   format(/2x,'Maximum distance between requested and real =',EN12.3)
 200  format(//1x,'S o u r c e s'/1x,13('=')// &
  ' Sources have been relocated to the nearest GLL node'// &
  '   Source  x-requested  z-requested   x-obtained   z-obtained     distance'/)

  end subroutine SO_init

!=====================================================================
!
! WARNING: output only for first source
  subroutine SO_check(so,dt,nt) 

  use stdio, only : IO_new_unit,IO_abort
  
  type(source_type), intent(inout) :: so(:)
  double precision, intent(in) :: dt
  integer, intent(in) :: nt
 
  integer :: it,fid

  fid = IO_new_unit()  
  open(unit=fid,file='SourcesTime_sem2d.tab',status='unknown')
  do it=1,nt
    write(fid,*) real(it*dt),real( STF_get(so(1)%stf,it*dt-so(1)%tdelay) )
  enddo
  close(fid)

  end subroutine SO_check


!=====================================================================
!
! adds point sources (not incident wave)
  subroutine SO_add(so,t, MxA)

  type(source_type), pointer :: so(:)
  double precision, intent(in) :: t
  double precision, intent(inout) :: MxA(:,:)
  
  double precision :: ampli
  integer :: k

  if (.not. associated(so)) return

  do k=1,size(so)

   ! get source amplitude
    ampli = STF_get(so(k)%stf,t-so(k)%tdelay)
  
   ! add source term
    select case (so(k)%mech%kind)
      case(tag_moment)
        call SRC_MOMENT_add(so(k)%mech%moment,ampli,MxA)
      case(tag_force)
        call FORCE_add(so(k)%mech%force,ampli,MxA)
    end select

  enddo

  end subroutine SO_add


!=====================================================================
! incident plane wave
subroutine SO_WAVE_get_VT( Vin, Tin, time,coord,n,src)

  double precision, intent(out) :: Vin(:,:), Tin(:,:)
  double precision, intent(in) :: coord(:,:),n(:,:),time
  type (source_type), intent(in) :: src

  double precision :: ampli,phase
  integer :: k

  call SRC_WAVE_get_VT(Vin,Tin, n,src%mech%wave)
  do k=1,size(coord,2)
    phase = SRC_WAVE_get_phase(src%mech%wave,time,coord(:,k)) -src%tdelay
    ampli = STF_get(src%stf, phase )
    Vin(k,:) = ampli*Vin(k,:)
    Tin(k,:) = ampli*Tin(k,:)
  enddo

end subroutine SO_WAVE_get_VT

endmodule sources
