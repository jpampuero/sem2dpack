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
  program main

  use problem_class
  use echo, only : iout,ECHO_banner,ItInfo,ItSnapshots,ItSnapshot1,title
  use input, only : read_main
  use init, only : init_main
  use solver, only : solve
  use plot_gen, only : PLOT_FIELD
  use fields_class, only : FIELDS_get_max
  use receivers, only : REC_store,REC_write
  use bc_gen, only : BC_write

  implicit none

  type(problem_type) :: pb
  real :: cputime0, cputime1, cputime2,cputime3
  integer :: it,iexec
  character(10) :: tagSnapshots

  call CPU_TIME(cputime0)

  call read_main(pb,iexec,input_file='Par.inp')
  call init_main(pb)

  call CPU_TIME( cputime1 )
  cputime0 = cputime1-cputime0

!*********************  s o l v e r   p h a s e ***********************

  pb%time%time = 0.d0

 !-- write initial data for faults, and eventually other BCs
  write(iout,'(5X,A)',advance='no') 'Exporting initial boundary data ...'
  call BC_write(pb%bc,pb%fields,0)
  write(iout,'(A)') '... [OK]'

  if (iexec == 0) call PLOT_FIELD(pb,'test',0,title)

  write(iout,400)
  call CPU_TIME( cputime1 )

  Time_Loop: do it=1,pb%time%nt

    pb%time%time = it*pb%time%dt
    
    call solve(pb)

    
  !-- CPU time info -----------------------------------------------------

    if (it == 1) then
      call CPU_TIME(cputime2)
      cputime2=cputime2-cputime1
      write(iout,'(/A,5(/2X,A,EN12.3),/)')   &
        '---  CPU TIME ESTIMATES (in seconds) :', &
        'CPU time for initialization . .', cputime0 ,&
        'CPU time per timestep . . . . .', cputime2,&
        'Total solver CPU time . . . . .', cputime2*pb%time%nt ,&
        '                 (mins) . . . .', cputime2*pb%time%nt/60. ,&
        '                 (hours). . . .', cputime2*pb%time%nt/3600.

  !-- if data check mode then stop ---------------------------------
      if(iexec == 0) then
        write(iout,*) 
        write(iout,*) '**********************************'
        write(iout,*) '* Aborting, data check mode only *'
        write(iout,*) '**********************************'
        write(iout,*) 
        call ECHO_banner('Program  S P E C F E M :  end',iout)
        stop
      endif

    endif

 !--- Intermediate OUTPUTS -------------------------------------------

    if (mod(it-1,ItInfo) == 0)  write(iout,'(A,I0,A,EN10.1,A,EN12.3)') & 
      'Timestep ',it,'  t = ',pb%time%time,'  dmax = ',FIELDS_get_max(pb%fields%displ)

    if (mod(it-ItSnapshot1,ItSnapshots) == 0 .and. it >= ItSnapshot1) then
      write(iout,'(/"Snapshot at timestep = ",I0)') it
      write(tagSnapshots,'(i3.3)') (it-ItSnapshot1)/ItSnapshots
      call PLOT_FIELD(pb,tagSnapshots,it,title)
      write(iout,*)
    endif

    !-- store seismograms
    if(associated(pb%rec)) call REC_store(pb%rec,it)

    !-- write data for faults, and eventually other BCs
    call BC_write(pb%bc,pb%fields,it)
  
  !------------------------------------------------------------------------

  enddo Time_Loop



!*****************  g l o b a l   o u t p u t   p h a s e  **************


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
  
 !  sauvegarder sismogrammes en fin de simulation
  if(associated(pb%rec)) call REC_write(pb%rec,iout)

!********************* e x i t   p h a s e *********************

  call ECHO_banner('Program  S P E C F E M :  end',iout)
  stop

  !  formats
 400  format(/1x,41('=')/,' =  T i m e   e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end program main
