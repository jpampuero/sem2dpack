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
  use echo, only : ECHO_read,echo_input,iout
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
  
  call ECHO_read(iin)

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
! SYNTAX : &GENERAL iexec,ngll,fmax /
!
! ARG: iexec    [int] [0] Run level:
!                       0 = just check
!                       1 = solve
! ARG: ngll     [int] [9] Number of GLL nodes per edge on each spectral element
!               ( polynomial order +1 ). Usually 5 to 9.
! ARG: fmax     [dble] [0.d0] Maximum frequency to be well resolved. Mandatory.
!               This is a target frequency, the code will check if it is
!               compatible with the mesh and eventually issue a warning. To
!               improve the resolution for a given fmax you must increase ngll 
!               (but you will have to use shorter timesteps) or refine/redesign 
!               the mesh.
! ARG: ndof     [int] [2] Number of degrees of freedom per node
!                       1 = SH waves, anti-plane
!                       2 = P-SV waves, in-plane
!
! END INPUT BLOCK

  subroutine read_gen(iexec,ndof,ngll,fmax,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  integer  :: iexec,ndof,ngll
  double precision :: fmax
  character(10) :: iexecname

  NAMELIST / GENERAL / iexec,ngll,ndof,fmax

  iexec = 0
  ndof = 2
  ngll = 9
  fmax = 0.d0

  rewind(iin)
  read(iin,GENERAL,END=100)

  if (ndof>2 .or. ndof<1) call IO_abort('GENERAL: ndof must be 1 or 2 (SH or P-SV)')
  if (ngll <= 0) call IO_abort('GENERAL: ngll must be positive')
  if (fmax <= 0.d0) call IO_abort('GENERAL: fmax null or missing')

  if (iexec==0) then
    iexecname = 'check'
  else
    iexecname = 'solve'
  endif
  
  if (echo_input) write(iout,200) iexecname,ngll,ndof,fmax

  return

100   call IO_abort('Input: GENERAL parameters not found')
200   format(//1x,'G e n e r a l   P a r a m e t e r s', &
  /1x,35('='),//5x,&
  'Execution mode . . . . . . . . . . . . . . . (iexec) = ',a/ 5x, &
  'Number of nodes per edge . . . . . . . . . . .(ngll) = ',I0/5x, &
  'Number of d.o.f per node . . . . . . . . . . .(ndof) = ',I0/5x, &
  'Highest frequency to be resolved . . . . . . .(fmax) = ',EN12.3)

  end subroutine read_gen

end module input
