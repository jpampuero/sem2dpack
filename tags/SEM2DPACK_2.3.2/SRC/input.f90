! SEM2DPACK version 2.3.2 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! Phone: (626) 395-3429
! Fax  : (626) 564-0715
! 
! http://www.seismolab.caltech.edu
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
module input

  implicit none
  private

  public :: read_main

contains

!=========================================================================
!
  subroutine read_main(pb,iexec,input_file)

  use problem_class, only : problem_type
  use stdio, only : IO_new_unit,IO_abort

  use time_evol, only : TIME_read
  use mat_gen, only : MAT_read
  use bc_gen, only : BC_read
  use mesh_gen, only : MESH_read
  use sources, only : SO_read
  use receivers, only : REC_read
  use plot_gen, only : read_plot_gen

  type(problem_type), intent(out) :: pb
  character(*), intent(in) :: input_file
  integer, intent(out) :: iexec

  integer :: iin

!-----------------------------------------------------------------------

  iin  = IO_new_unit()
  open (iin,file=input_file)
  
  call read_gen(iexec,pb%fields%ndof,pb%grid%ngll,pb%grid%fmax,iin)

 !---- mesh generation parameters     
  call MESH_read(pb%mesh,iin)

 !---- timescheme settings
  call TIME_read(pb%time,iin)

 !---- material properties
  call MAT_read(pb%matinp,iin)

 !---- boundary conditions properties     
  call BC_read(pb%bc,iin) 
 
 !---- Collocated forces or pressure sources
  call SO_read(pb%src,iin, pb%fields%ndof)

 !---- receivers 
  call REC_read(pb%rec,iin)

 !---- plots settings
  call read_plot_gen(iin, pb%fields%ndof)

 !---- close input file
  close(iin)

 end subroutine read_main

!===============================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : GENERAL
! PURPOSE: General parameters
! SYNTAX : &GENERAL iexec, ngll, fmax, title, verbose, itInfo /
!
! ARG: iexec    [int] [0] Run level:
!                       0 = just check
!                       1 = solve
! ARG: ngll     [int] [9] Number of GLL nodes per edge on each spectral element
!                ( polynomial order +1 ). Usually 5 to 9.
! ARG: fmax     [dble] [0.d0] Maximum frequency to be well resolved. Mandatory.
!                This is a target frequency, the code will check if it is
!                compatible with the mesh and issue a warning if not. To
!                improve the resolution for a given fmax you must increase ngll 
!                (but you will have to use shorter timesteps) or refine/redesign 
!                the mesh.
! ARG: ndof     [int] [2] Number of degrees of freedom per node
!                       1 = SH waves, anti-plane
!                       2 = P-SV waves, in-plane
! ARG: title    [word] [none] Title of the simulation
! ARG: verbose  [char(4)] ['1101'] Print progress information during each phase:
!                       verbose(1) = input phase
!                       verbose(2) = initialization phase
!                       verbose(3) = check phase
!                       verbose(4) = solver phase
!                Example: '0001' is verbose only during solver.
! ARG: itInfo   [int] [100] Frequency (in number of timesteps) for printing
!                progress information during the solver phase.
!
! END INPUT BLOCK

  subroutine read_gen(iexec,ndof,ngll,fmax,iin)

  use echo
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  integer  :: iexec,ndof,ngll
  double precision :: fmax
  character(10) :: iexecname
  character(4) :: verbose

  NAMELIST / GENERAL / iexec,ngll,ndof,fmax,title,verbose,itInfo

  iexec = 0
  ndof = 2
  ngll = 9
  fmax = 0.d0
  title   = ''
  verbose = '1101'
  itInfo  = 100

  rewind(iin)
  read(iin,GENERAL,END=100)

  if (ndof>2 .or. ndof<1) call IO_abort('GENERAL input block: ndof must be 1 or 2 (SH or P-SV)')
  if (ngll <= 0) call IO_abort('GENERAL input block: ngll must be positive')
  if (fmax <= 0.d0) call IO_abort('GENERAL input block: fmax null or missing')
  if (itInfo<=0) call IO_abort('GENERAL input block: itInfo must be positive')

  if (iexec==0) then
    iexecname = 'check'
  else
    iexecname = 'solve'
  endif
  
  call ECHO_set(verbose)

  if (echo_input) then
    write(iout,*)
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*            I n p u t   p h a s e            *'
    write(iout,'(a)') '***********************************************'
    write(iout,*)
    write(iout,200) iexecname,ngll,ndof,fmax, &
      echo_input,echo_init,echo_check,echo_run,itInfo
  endif

  return

100   call IO_abort('Input: GENERAL parameters not found')
200   format(//1x,'G e n e r a l   P a r a m e t e r s', &
  /1x,35('='),//5x,&
  'Execution mode . . . . . . . . . . . . . . . (iexec) = ',a/ 5x, &
  'Number of nodes per edge . . . . . . . . . . .(ngll) = ',I0/5x, &
  'Number of d.o.f per node . . . . . . . . . . .(ndof) = ',I0/5x, &
  'Highest frequency to be resolved . . . . . . .(fmax) = ',EN12.3/5x, &
  'Print progress information during ',/5x, &
  '            input phase  . . . . . . . .(verbose(1)) = ',L1/ 5x, &
  '            initialization phase . . . .(verbose(2)) = ',L1/ 5x, &
  '            checking phase . . . . . . .(verbose(3)) = ',L1/ 5x, &
  '            solver phase . . . . . . . .(verbose(4)) = ',L1/ 5x, &
  'Frequency for solver progress information  .(itInfo) = ',I0/ 5x)

  end subroutine read_gen

end module input
