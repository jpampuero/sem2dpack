! SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
  program main

  use problem_class
  use echo, only : iout,ECHO_banner,ItInfo,title
  use input, only : read_main
  use init, only : init_main
  use solver, only : solve
  use plot_gen, only : PLOT_FIELD
  use receivers, only : REC_store,REC_write
  use bc_gen, only : BC_write
  use energy, only : energy_compute, energy_write, stress_glut_write
  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  implicit none

  type(problem_type) :: pb
  real :: cputime0, cputime1, cputime2,cputime3
  integer :: it,iexec,nt
  integer, parameter :: NT_CHECK=10

  call CPU_TIME(cputime0)

  call ECHO_banner('Program  S E M 2 D P A C K : start',iout)

!*************  i n p u t   p h a s e  ***************

  call read_main(pb,iexec,input_file='Par.inp')

!*************  i n i t i a l i z a t i o n   p h a s e  ***************

  call init_main(pb)

  pb%time%time = 0.d0

  !-- store seismograms at time = 0
  if(associated(pb%rec) .and.iexec>0) call REC_store(pb%rec,0,pb%grid)

  call PLOT_FIELD(pb,0,title,iout)

!*********************  s o l v e r   p h a s e  ***********************

    write(iout,*)
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*           S o l v e r   p h a s e           *'
    write(iout,'(a)') '***********************************************'
    write(iout,*)

  if (iexec==0) then
    nt = NT_CHECK
  else
    write(iout,200) 0,0d0,maxval(abs(pb%fields%veloc)),maxval(abs(pb%fields%displ))
    nt = pb%time%nt
  endif

  call CPU_TIME( cputime1 )
  cputime0 = cputime1-cputime0

  Time_Loop: do it=1,nt
    pb%time%time = it*pb%time%dt
    
    call solve(pb)
    
  !-- CPU time info -----------------------------------------------------

    if (it == NT_CHECK) then
      call CPU_TIME(cputime2)
      cputime2=(cputime2-cputime1)/dble(NT_CHECK)
      write(iout,'(/A,5(/2X,A,EN12.3),/)')   &
        '---  CPU TIME ESTIMATES (in seconds) :', &
        'CPU time for initialization . .', cputime0 ,&
        'CPU time per timestep . . . . .', cputime2,&
        'Total solver CPU time . . . . .', cputime2*pb%time%nt ,&
        '                 (mins) . . . .', cputime2*pb%time%nt/60. ,&
        '                 (hours). . . .', cputime2*pb%time%nt/3600.

    endif

 !--- Intermediate OUTPUTS -------------------------------------------

    if (mod(it,ItInfo) == 0)  &
      write(iout,200) it,pb%time%time,maxval(abs(pb%fields%veloc)),maxval(abs(pb%fields%displ))

    if (iexec==0) cycle

    !-- snapshot outputs
    call PLOT_FIELD(pb,it,title,iout)

    !-- store seismograms
    if(associated(pb%rec)) call REC_store(pb%rec,it,pb%grid)

    !-- write data for faults, and possibly other BCs
    call BC_write(pb%bc,it,pb%fields%displ,pb%fields%veloc)
  
    !-- export energies
    if (COMPUTE_ENERGIES) then
      call energy_compute(pb%energy,pb%matpro,pb%matwrk,pb%grid,pb%fields)
      call energy_write(pb%energy,pb%time%time)
    endif

    if (COMPUTE_STRESS_GLUT) call stress_glut_write(pb%energy,pb%time%time)

  !------------------------------------------------------------------------

  enddo Time_Loop



!*****************  g l o b a l   o u t p u t   p h a s e  **************

  if(iexec == 0) then
   !-- if data check mode then stop ---------------------------------
    write(iout,*) 
    write(iout,*) '**********************************'
    write(iout,*) '* Aborting, data check mode only *'
    write(iout,*) '**********************************'
    write(iout,*) 

  else

   !-- CPU TIME INFO
    call CPU_TIME(cputime3)      
    cputime3 = cputime3-cputime1
    write(iout,'(//A,5(/2X,A,EN12.3),/)')   &
        '---  CPU TIME INFORMATION (in seconds) :', &
        'CPU time for initialization . .', cputime0 ,&
        'CPU time per timestep . . . . .', cputime3/pb%time%nt ,&
        'Total solver CPU time . . . . .', cputime3,&
        '                 (mins) . . . .', cputime3/60. ,&
        '                 (hours). . . .', cputime3/3600.
  
   !  save seismograms
    if(associated(pb%rec)) call REC_write(pb%rec,iout)
  endif

!********************* e x i t   p h a s e *********************

  call ECHO_banner('Program  S E M 2 D P A C K :  end',iout)

  stop
200 format("Timestep #",I8,"  t = ",EN11.3,"  vmax = ",EN11.3,"  dmax = ",EN11.3)

  end program main