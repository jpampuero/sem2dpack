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
module receivers

  use stdio, only: IO_abort
  use constants, only : NDIME

  implicit none
  private

  type rec_type
    private
    double precision, pointer :: coord(:,:),field(:,:)
    real, pointer :: sis(:,:,:)
    integer, pointer :: iglob(:)
    double precision :: tsamp
    integer :: nx,nt,isamp,irepr
    logical :: AtNode
    double precision, pointer :: interp(:,:)
    integer, pointer :: einterp(:)
    character :: SeisField
  end type rec_type

  public :: rec_type,REC_read,REC_init,REC_store,REC_write,REC_inquire

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : REC_LINE
! PURPOSE: Defines a line of receivers
! SYNTAX : &REC_LINE number,isamp,field,first,last,file,AtNode,irepr /
!
! ARG: number   [int] [0] Number of stations in the line
! ARG: isamp    [int] [1] Sampling stride (in number of timesteps). Note that
!                for stability reasons the timestep can be very small.
! ARG: field    [char] ['V'] The field in the seismogram:
!                               'D'     displacement
!                               'V'     velocity
!                               'A'     acceleration
! ARG: first    [dble(2)] Receivers can be located along a line,
!                this is the position (x,z) of the first receiver
! ARG: last     [dble(2)] Position (x,z) of the last receiver,
!                other receivers will be located with regular spacing
!                between First and Last.
! ARG: file     [name] ['none'] Station positions can instead be read 
!                from an ASCII file, with 2 columns per line:
!                (1) X position (in m) and
!                (2) Z position (in m)
! ARG: AtNode	[log] [T] Relocate the stations at the nearest GLL node
! ARG: irepr    [char] ['D'] Abscissa for the seismic multitrace plot:
!                               'X' Horizontal position
!                               'Z' Depth
!                               'D' Distance to the first station
!
! NOTE   : to locate receivers at the free surface set their vertical position 
!          above the free surface and AtNode=T
!
! END INPUT BLOCK

  subroutine REC_read(rec,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit, IO_file_length

  type(rec_type), pointer :: rec
  integer, intent(in) :: iin

  double precision :: first(NDIME),last(NDIME),init_double
  integer :: i,irec,iin2,isamp,number
  logical :: AtNode
  character(50) :: file
  character :: field,irepr
  
  NAMELIST / REC_LINE / number,isamp,field,first,last,AtNode,file,irepr

  init_double = huge(init_double) ! set to an unlikely value

  number     = 0
  isamp      = 1
  field      = 'V'
  first      = init_double
  last       = init_double
  AtNode     = .true.
  file       = 'none'
  irepr      = 'D'

  rewind(iin)
  read(iin,REC_LINE,END = 200) 
  allocate(rec)

  if (number <= 0) call IO_abort('REC_read: "number" is null or missing')

  rec%isamp  = isamp
  rec%SeisField = field
  rec%AtNode = AtNode

  if ( field /='D' .and. field /='V ' .and. field /='A') &
    call IO_abort('REC_read: parameter field has wrong value [D,V,A]')

  select case(irepr)
    case('X'); rec%irepr  = 1
    case('Z'); rec%irepr  = 2
    case('D'); rec%irepr  = 3
    case default; call IO_abort('REC_read: unknown "irepr" [X,Z,D]')
  end select
  
  if (echo_input) write(iout,100) number,isamp,field,irepr

  if (file == 'none') then

    if ( any(first == init_double) ) &
      call IO_abort('REC_read: you must set "first" station coordinates')
    if ( any(last == init_double) ) &
      call IO_abort('REC_read: you must set "last" station coordinates')
    rec%nx   = number
    allocate(rec%coord(NDIME,rec%nx))
    if (number>1) then
      do i = 1,rec%nx
        rec%coord(:,i) = first + (i-1)/dble(rec%nx-1) *(last-first)
      enddo
    else
      rec%coord(:,1) = first    
    endif

  else
    rec%nx = IO_file_length(file)
    allocate(rec%coord(NDIME,rec%nx))
    iin2 = IO_new_unit()
    open(iin2,file=file,status='old')
    do i=1,rec%nx
      read(iin2 ,*) rec%coord(:,i) 
    enddo
    close(iin2)
  endif

  return

  100 format(//1x,'R e c e i v e r s', &
  /1x,17('='),//5x, &
  'Number of receivers . . . . . . . . . . . . (number) = ',I0/5x, &
  'Subsampling for seismograms recording . . . .(isamp) = ',I0/5x, &
  'Field recorded. . . . . . . . . . . . . . . .(field) = ',A/5x, &
  'Axis of the seismogram plot . . . . . . . . .(irepr) = ',A)

  200 return

  end subroutine REC_read

!=====================================================================
!
!!  Initialize
  subroutine REC_init(rec,grid,time,fields)

    use fields_class, only : fields_type
    use time_evol, only : timescheme_type
    use spec_grid, only : sem_grid_type
    use memory_info
    use echo, only: iout,info=>echo_init
    use stdio, only: IO_new_unit

    type(rec_type), intent(inout) :: rec
    type(sem_grid_type), intent(in) :: grid
    type(fields_type), intent(in) :: fields
    type(timescheme_type), intent(in) :: time

    integer :: ounit,irec

    call REC_posit(rec,grid)

    rec%tsamp = time%dt*rec%isamp

    rec%nt = time%nt/rec%isamp +1
    if (rec%nt == 0) call IO_abort('REC_init: zero samples')
      
    allocate(rec%sis(rec%nt,rec%nx,fields%ndof))
    rec%sis = 0.
    call storearray('rec%sis',size(rec%sis),idouble)

    select case(rec%SeisField)
      case('D'); rec%field => fields%displ
      case('V'); rec%field => fields%veloc
      case('A'); rec%field => fields%accel
      case default
        call IO_abort('REC_init: Wrong field code for seismograms display')
    end select

  if (info) then
    write(iout,*)
    write(iout,100) 'Sampling rate (Hz)        = ',1d0/rec%tsamp
    write(iout,100) 'Sampling timestep (secs)  = ',rec%tsamp
    write(iout,200) 'Total number of samples   = ',rec%nt
    write(iout,200) 'Number of receivers       = ',rec%nx
    write(iout,*)
  endif

! Save headers
  ounit = IO_new_unit()
  open(ounit,file='SeisHeader_sem2d.hdr',status='replace')
  write(ounit,*) 'DT NSAMP NSTA'
  write(ounit,*) real(rec%tsamp),rec%nt,rec%nx
  write(ounit,*) 'XSTA ZSTA'
  do irec=1,rec%nx
    write(ounit,*) rec%coord(:,irec)
  enddo
  close(ounit)

100 format(2X,A,EN12.3)
200 format(2X,A,I0)

  end subroutine REC_init

!=====================================================================
!
!!  Real receivers position
  subroutine REC_posit(rec,grid)

  use echo, only: iout,info=>echo_init
  use spec_grid, only : sem_grid_type,SE_find_nearest_node,SE_init_interpol,SE_find_point
  use memory_info

  type(rec_type), intent(inout) :: rec
  type(sem_grid_type), intent(in) :: grid

  integer :: iglob_tmp(rec%nx)
  double precision :: distmax,xs,zs,xp,zp,dist, xi,eta, new_coord(NDIME)
  integer :: n,ip,ipoint,irec


  distmax = 0.d0

  if (rec%AtNode) then

    if (info) write(iout,200)
    iglob_tmp = 0
    irec = 0
    do n=1,rec%nx
      call SE_find_nearest_node(rec%coord(:,n),grid,ipoint,distance = dist)
     ! do not keep it if already in list
      if ( irec > 1 .and. any(iglob_tmp(:irec) == ipoint)) cycle
      irec = irec+1
      iglob_tmp(irec) = ipoint
      if (info) then
        write(iout,150) irec,rec%coord(:,n),grid%coord(:,ipoint),dist
        distmax = max(dist,distmax)
      endif
    enddo
    if (irec < rec%nx) then 
      rec%nx = irec
      deallocate(rec%coord)
      allocate(rec%coord(NDIME,rec%nx))
    endif
    call storearray('rec%coord',size(rec%coord),idouble)
    allocate(rec%iglob(rec%nx))
    call storearray('rec%iglob',size(rec%iglob),isngl)
    rec%iglob = iglob_tmp(:rec%nx)
    rec%coord(:,:) = grid%coord(:,rec%iglob)

  else

    if (info) write(iout,300)
    allocate(rec%interp(grid%ngll*grid%ngll,rec%nx))
    allocate(rec%einterp(rec%nx))
    do n=1,rec%nx
      call SE_find_point(rec%coord(:,n),grid,rec%einterp(n),xi,eta,new_coord)
      call SE_init_interpol(xi,eta,rec%interp(:,n),grid)
      dist = sqrt( (rec%coord(1,n)-new_coord(1))**2 + (rec%coord(2,n)-new_coord(2))**2 )
      if (info) then
        write(iout,150) n,rec%coord(:,n),new_coord,dist
        distmax = max(dist,distmax)
      endif
      rec%coord(:,n)=new_coord
    enddo

  endif

  if (info) write(iout,160) distmax

 150   format(2x,i7,5(1x,EN12.3))
 160   format(/2x,'Maximum distance between asked and real =',EN12.3)
 200  format(//1x,'R e c e i v e r s'/1x,17('=')// &
  ' Receivers have been relocated to the nearest GLL node'// &
  ' Receiver  x-requested  z-requested   x-obtained   z-obtained   distance'/)

 300  format(//1x,'R e c e i v e r s'/1x,17('=')// &
  ' Receiver            x            z   element           xi          eta'/)
 310   format(2x,i7,1x,EN12.3,1x,EN12.3,1x,i9,1x,EN12.3,EN12.3)

  end subroutine REC_posit


!=====================================================================
!
!!  Store the seismograms
  subroutine REC_store(rec,it,grid)

  use spec_grid, only : sem_grid_type

  type(rec_type), intent(inout) :: rec
  integer, intent(in) :: it
  type(sem_grid_type), intent(in) :: grid

  integer :: itsis,n,i,k,j,iglob
  double precision, allocatable :: vloc(:,:)

  if ( mod(it,rec%isamp) /= 0 ) return

  itsis = it/rec%isamp +1
  if (itsis > rec%nt) call IO_abort('receivers.REC_store: storage is full')

  if (rec%AtNode) then
    rec%sis(itsis,:,:) = rec%field(rec%iglob,:)
  else
    allocate(vloc(grid%ngll*grid%ngll,size(rec%field,2)))
    do n=1,rec%nx
      k=1
      do j=1,grid%ngll
      do i=1,grid%ngll
        iglob = grid%ibool(i,j,rec%einterp(n))
        vloc(k,:) = rec%field(iglob,:)
        k=k+1
      enddo
      enddo
      rec%sis(itsis,n,:) = matmul(rec%interp(:,n),vloc)
    enddo
    deallocate(vloc)
  endif
    
  end subroutine REC_store



!=====================================================================
!
!!  Export the seismograms
  subroutine REC_write(rec,iout)

  use stdio, only: IO_new_unit
  double precision, parameter :: factorxsu = 3.5d0

  type(rec_type), intent(inout) :: rec
  integer, intent(in) :: iout

  double precision :: xval(rec%nx),xref,zref
  integer :: iol,ounit
  character(30) :: ylabel
  character, parameter :: posvars(3) = (/ 'X','Z','D' /)
  character :: posvar

!---- binary data -------------------------------------------------------
  write(iout,*) 'Storing sismos data (SEP format) ...'
  INQUIRE( IOLENGTH=iol ) rec%sis(:,:,1)
  ounit = IO_new_unit()

  if ( size(rec%sis,3)==1) then

    open(ounit,file='Uy_sem2d.dat',status='replace',access='direct',recl=iol)
    write(ounit,rec=1) rec%sis(:,:,1)
    close(ounit)

  else

    open(ounit,file='Ux_sem2d.dat',status='replace',access='direct',recl=iol)
    write(ounit,rec=1) rec%sis(:,:,1)
    close(ounit)
  
    open(ounit,file='Uz_sem2d.dat',status='replace',access='direct',recl=iol)
    write(ounit,rec=1) rec%sis(:,:,2)
    close(ounit)

  endif

!=== Seismic Unix scripts: =================================================

  posvar = posvars(rec%irepr)

  select case(rec%irepr)
    case(1,2) ! recepteurs suivant coordonnee X ou Z
      xval = rec%coord(rec%irepr,:)
    case(3) ! recepteurs en distance
      xref = rec%coord(1,1)
      zref = rec%coord(2,1)
      xval = sqrt((rec%coord(1,:) - xref)**2 + &
              (rec%coord(2,:) - zref)**2)
  end select

  select case(rec%SeisField)
    case('D'); ylabel='Displacement (m)'
    case('V'); ylabel='Velocity (m/s)'
    case('A'); ylabel='Acceleration (m/s^2)'
  end select

! station "coordinates" for plots
  open(ounit,file='x2_sem2d.tab',status='unknown')
  write(ounit,'(1000(f0.2:","))') xval
               !this number > max nb stations ever
  close(ounit)

!-- Xwindow -------------------------------------------------------------------
  open(ounit,file='Xline_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set x2 = `cat x2_sem2d.tab`'
  write(ounit,100) 'xwigb',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Ux component" < Ux_sem2d.dat'
  write(ounit,100) 'xwigb',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Uz component" < Uz_sem2d.dat'
  close(ounit)

!-- PostScript ---------------------------------------------------------------
  open(ounit,file='PSline_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set x2 = `cat x2_sem2d.tab`'
  write(ounit,100) 'pswigp',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Ux component" < Ux_sem2d.dat >! UxPoly_sem2d.ps'
  write(ounit,100) 'pswigp',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Uz component" < Uz_sem2d.dat >! UzPoly_sem2d.ps'
  close(ounit)

!-- one trace for Xwindow --------------------------------------------

  open(ounit,file='Xtrace_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set nt = ',rec%nt
  write(ounit,*) 'set dt = ',real(rec%tsamp)
  write(ounit,*) '@ trace=10000'
  write(ounit,*) 'while ($trace > -1)'
  write(ounit,*) 'echo Trace number ?'
  write(ounit,*) 'set rep=$<'
  write(ounit,*) '@ trace = $rep'
  write(ounit,*) 'echo Trace asked : $trace'
  write(ounit,*) '# traces commencent a zero dans format SEP'
  write(ounit,*) '@ septrace = $trace - 1'
  write(ounit,*) 'subset < Ux_sem2d.dat n1=$nt', &
   ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
   ' label1="Time (s)" label2="',trim(ylabel),'" n=$nt style=normal d1=$dt', &
   ' title="Ux component trace "$trace &'
  write(ounit,*) 'subset < Uz_sem2d.dat n1=$nt', &
   ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
   ' label1="Time (s)" label2="',trim(ylabel),'" n=$nt style=normal d1=$dt', &
   ' title="Uz component trace "$trace &'
  write(ounit,*) 'end'
  close(ounit)

!-- one trace for PostScript --------------------------------------------

  open(ounit,file='PStrace_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set nt = ',rec%nt
  write(ounit,*) 'set dt = ',real(rec%tsamp)
  write(ounit,*) '@ trace=10000'
  write(ounit,*) 'while ($trace > -1)'
  write(ounit,*) 'echo Trace number ?'
  write(ounit,*) 'set rep=$<'
  write(ounit,*) '@ trace = $rep'
  write(ounit,*) 'echo Trace asked : $trace'
  write(ounit,*) '# traces commencent a zero dans format SEP'
  write(ounit,*) '@ septrace = $trace - 1'
  write(ounit,*) 'rm -f UxTrace{$trace}_sem2d.ps UzTrace{$trace}_sem2d.ps'
  write(ounit,*) 'subset < Ux_sem2d.dat n1=$nt', &
   '" if2s=$septrace n2s=1 | psgraph label1="Time (s)" label2="',trim(ylabel),&
   ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
   ' title="Ux component trace "$trace > UxTrace{$trace}_sem2d.ps'
  write(ounit,*) 'subset < Uz_sem2d.dat n1=$nt', &
   '" if2s=$septrace n2s=1 | psgraph label1="Time (s)" label2="',trim(ylabel),&
   ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
   ' title="Uz component trace "$trace > UzTrace{$trace}_sem2d.ps'
  write(ounit,*) 'end'
  close(ounit)

 100  format(A,' xcur=',F0.2,' n1=',I0,' d1=',F0.8,' label1="Time (s)" label2="',A, &
             ' (m)" n2=',I0,' x2=$x2',A )
 110  format('#!/bin/csh -f')
  
  end subroutine REC_write


!=====================================================================
!
  subroutine REC_inquire(rec,coord,isamp)

  type(rec_type), intent(in) :: rec
  double precision, pointer, optional :: coord(:,:)
  integer, optional :: isamp

  if (present(coord))  coord => rec%coord
  if (present(isamp))  isamp = rec%isamp

  end subroutine REC_inquire

end module receivers
