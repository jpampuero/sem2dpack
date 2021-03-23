program main
#include <petsc/finclude/petscksp.h>
  use petscksp
  
  use problem_class
  use petsc_objects
  use echo, only : iout,ECHO_banner,ItInfo,title
  use input, only : read_main
  use init, only : init_main
!  use solver, only : solve
  use plot_gen, only : PLOT_FIELD
  use receivers, only : REC_store,REC_write
  use bc_gen, only : BC_write
  use energy, only : energy_compute, energy_write, stress_glut_write
  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  implicit none

  type(problem_type) :: pb
  type(petsc_objects_type) :: petobj
  real :: cputime0, cputime1, cputime2,cputime3
  integer :: it, iexec, rank
  integer, parameter :: NT_CHECK=10
  PetscErrorCode :: ierr
  CHARACTER(len=32) :: inputfile

!------------- Start Petsc session ---------------
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CHKERRQ(ierr)

  ! get rank
  call MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr)

  ! allow different input file names
  if (iargc()>0) then
      call getarg(1, inputfile)
  else
      inputfile='Par.inp'
  endif

  call CPU_TIME(cputime0)

  if (rank==0) call ECHO_banner('Program  S E M 2 D P A C K : start', iout)
!*************  i n p u t   p h a s e  **************

  call read_main(pb,iexec,input_file=inputfile)

!*************  i n i t i a l i z a t i o n   p h a s e  ***************
  call init_main(pb, petobj, ierr)

  CHKERRQ(ierr)
!----------------------------------------------------

  pb%time%time = 0.d0

  !-- store seismograms at time = 0
  if (rank==0) then
     if(associated(pb%rec) .and.iexec>0) call REC_store(pb%rec,0,pb%grid)
  end if

  if (rank==0) call PLOT_FIELD(pb,0,title,iout)

!*********************  s o l v e r   p h a s e  ***********************

  if (rank==0) then
      write(iout,*)
      write(iout,'(a)') '***********************************************'
      write(iout,'(a)') '*           S o l v e r   p h a s e           *'
      write(iout,'(a)') '***********************************************'
      write(iout,*)
      write(iout,200) 0,0d0,maxval(abs(pb%fields%veloc)),maxval(abs(pb%fields%displ))
  end if

  call CPU_TIME( cputime1 )
  cputime0 = cputime1-cputime0

  it = 1
  pb%time%time = pb%time%time + pb%time%dt

! comment out to test the input 
  do while (pb%time%time < pb%time%total .and. pb%time%EQNum <= pb%time%EQNumMax) 

    ! save the time step number
    pb%time%it = it
 !   
 !   call solve(pb, petobj)

 ! !-- CPU time info -----------------------------------------------------

 !   if (it == NT_CHECK) then
 !     call CPU_TIME(cputime2)
 !     cputime2=(cputime2-cputime1)/dble(NT_CHECK)
 !     write(iout,'(/A,5(/2X,A,EN12.3),/)')   &
 !       '---  CPU TIME ESTIMATES (in seconds) :', &
 !       'CPU time for initialization . .', cputime0 ,&
 !       'CPU time per timestep . . . . .', cputime2,&
 !       'Total solver CPU time . . . . .', cputime2*pb%time%nt ,&
 !       '                 (mins) . . . .', cputime2*pb%time%nt/60. ,&
 !       '                 (hours). . . .', cputime2*pb%time%nt/3600.
 !     if (iexec==0) exit
 !   endif

 !!--- Intermediate OUTPUTS -------------------------------------------

 !   if (mod(it,ItInfo) == 0) then
 !     write(iout,200) it,pb%time%time,maxval(abs(pb%fields%veloc)),maxval(abs(pb%fields%displ))
 !     if (pb%time%kind == 'adaptive') write(iout,*) 'dt = ', pb%time%dt
 !     if (.not. pb%time%isdynamic) then 
 !         write(iout,*) & 
 !         'PCG iteration numbers (first/second passes) = ', pb%time%pcg_iters
 !         else
 !         write(iout,*) & 
 !         'Max NR iteration numbers (first/second passes) = ', pb%time%nr_iters
 !     end if
 !   endif

 !   if (iexec>0) then

 !     !-- snapshot outputs
 !     call PLOT_FIELD(pb,it,title,iout)

 !     !-- store seismograms
!!      if(associated(pb%rec) .and. pb%time%kind .ne. 'quasi-static') call REC_store(pb%rec,it,pb%grid)
 !     if(associated(pb%rec)) call REC_store(pb%rec,it,pb%grid)

 !     !-- write data for faults, and possibly other BCs
 !     call BC_write(pb%bc,pb%time,pb%fields%displ,pb%fields%veloc)
 !   
 !     !-- export energies
 !     if (COMPUTE_ENERGIES) then
 !       call energy_compute(pb%energy,pb%matpro,pb%matwrk,pb%grid,pb%fields)
 !       call energy_write(pb%energy,pb%time%time)
 !     endif

 !     if (COMPUTE_STRESS_GLUT) call stress_glut_write(pb%energy,pb%time%time)

 !   endif
 ! !------------------------------------------------------------------------
 !   if (pb%time%kind=='adaptive') then
 !     write(pb%time%ou_time, *) pb%time%it, pb%time%dt, pb%time%time
 !   end if

    it = it + 1
    pb%time%time = pb%time%time + pb%time%dt

  end do

!*****************  g l o b a l   o u t p u t   p h a s e  **************
  pb%time%nt = it

  if (iexec == 0) then
   !-- if data check mode then stop ---------------------------------
    write(iout,*) 
    write(iout,*) '**********************************'
    write(iout,*) '* Aborting, data check mode only *'
    write(iout,*) '**********************************'
    write(iout,*) 

  else
   if (rank==0) then
   !-- CPU TIME INFO
    call CPU_TIME(cputime3)      
    cputime3 = cputime3-cputime1
    write(iout,'(//A,5(/2X,A,EN12.3),/)')   &
        '---  CPU TIME INFORMATION (in seconds) :', &
        'CPU time for initialization . .', cputime0 ,&
        'CPU time per timestep . . . . .', cputime3/dble(pb%time%nt) ,&
        'Total solver CPU time . . . . .', cputime3,&
        '                 (mins) . . . .', cputime3/60. ,&
        '                 (hours). . . .', cputime3/3600.
  
   !  save seismograms
     if(associated(pb%rec)) call REC_write(pb%rec,iout)
   end if

  endif

!********************* e x i t   p h a s e *********************

  if (rank==0) then
    call ECHO_banner('Program  S E M 2 D P A C K :  end',iout)
  end if

  !---------- End Petsc session and delete Petsc objects --------
  call destroyPetscStruct(petobj, ierr)
  CHKERRQ(ierr)
  call PetscFinalize(ierr)
  !------------------- Petsc session End-------------------------

  !stop
200 format("Timestep #",I8,"  t = ",EN12.3,"  vmax = ",EN12.3,"  dmax = ",EN12.3)
  end program main
