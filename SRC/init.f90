module init

! INIT:  i n i t i a l i z a t i o n   p h a s e

  implicit none
  private

  public :: init_main

contains

!=====================================================================
! INIT_MAIN:
!

subroutine init_main(pb,InitFile)

  use problem_class, only : problem_type
  use mesh_gen, only : MESH_build
  use spec_grid, only : SE_init
  use bc_gen, only : BC_init
  use mat_mass, only : MAT_MASS_init
  use mat_gen, only : MAT_init_prop, MAT_init_work, MAT_diag_stiffness_init
  use mat_elastic, only : MAT_IsElastic
  use time_evol, only : TIME_init, TIME_needsAlphaField, TIME_getTimeStep, TIME_getNbTimeSteps
  use plot_gen, only : PLOT_init
  use receivers, only : REC_init,REC_inquire
  use sources, only : SO_init,SO_check
  use fields_class, only : FIELDS_init,FIELDS_read
  use memory_info
  use echo, only : iout,info=>echo_init, fmt1,fmtok
  use energy, only : energy_init, stress_glut_init
  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  type(problem_type), intent(inout) :: pb
  character(50), intent(in), optional :: InitFile
  
  double precision :: grid_cfl
  integer :: recsamp,i,j
  logical :: init_cond

  init_cond = present(InitFile)

  if (info) then
    write(iout,*)
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*   I n i t i a l i z a t i o n   p h a s e   *'
    write(iout,'(a)') '***********************************************'
    write(iout,*)
  endif

  call MESH_build(pb%grid%fem,pb%mesh) ! finite elements mesh
  call SE_init(pb%grid) ! spectral elements mesh
  
 ! initialize material properties database
  call MAT_init_prop(pb%matpro, pb%matinp, pb%grid)

  call PLOT_init(pb%grid)
  call CHECK_grid(pb%grid,pb%matpro,pb%fields%ndof,grid_cfl)
  call TIME_init(pb%time,grid_cfl) ! define time evolution coefficients
 
 ! define work arrays and data
  call MAT_init_work(pb%matwrk,pb%matpro,pb%grid,pb%fields%ndof,TIME_getTimeStep(pb%time))

 ! initialise fields
  if (info) write(iout,fmt1,advance='no') 'Initializing kinematic fields'
  call FIELDS_init(pb%fields,pb%grid%npoin, TIME_needsAlphaField(pb%time))
  if(init_cond) call FIELDS_read(pb%fields,InitFile) ! read from a file
  if (info) then
    write(iout,fmtok)
    write(iout,'(7X,A,EN12.3)') 'Max displ = ',maxval(abs(pb%fields%displ))
    write(iout,'(7X,A,EN12.3)') 'Max veloc = ',maxval(abs(pb%fields%veloc))
    write(iout,*)
  endif

 ! build the mass matrix
  call MAT_MASS_init(pb%rmass,pb%matpro,pb%grid,pb%fields%ndof)

 ! initialize physical parameters of all boundary conditions
  if (info) write(iout,fmt1,advance='no') 'Defining boundary conditions'
  call BC_init(pb%bc,pb%grid,pb%matpro,pb%rmass,pb%time,pb%src,pb%fields%displ,pb%fields%veloc)
  if (info) write(iout,fmtok)

 ! define position of receivers and allocate database
  if (associated(pb%rec)) then
    if (info) write(iout,fmt1,advance='no') 'Initializing receivers'
    call REC_init(pb%rec,pb%grid,pb%time,pb%fields)
    if (info) write(iout,fmtok)
  endif

 ! initialize sources
  if (associated(pb%src)) then
    if (info) write(iout,fmt1,advance='no') 'Initializing sources'
    call SO_init(pb%src,pb%grid,pb%matpro) 
    if (associated(pb%rec)) then
      call REC_inquire(pb%rec, isamp=recsamp )
    else
      recsamp=1
    endif
    call SO_check(pb%src, TIME_getTimeStep(pb%time)*recsamp, TIME_getNbTimeSteps(pb%time)/recsamp)
    if (info) write(iout,fmtok)
  endif

 !  list memory usage
  call MEMO_echo()

 ! invert the mass matrix
 !   pb%rmass = 1.d0 / pb%rmass
 ! NOTE: the line above crashes with segmentation fault if pb%rmass exceeds the stack size
 !       Under bash the stack size can be extended by the command: ulimit -s unlimited
 !       We rather expand the loop:
  do j=1,pb%fields%ndof
  do i=1,pb%grid%npoin
    pb%rmass(i,j) = 1.d0 / pb%rmass(i,j)
  enddo
  enddo

 ! macroscopic outputs
  if (COMPUTE_ENERGIES) call energy_init(pb%energy)
  if (COMPUTE_STRESS_GLUT) call stress_glut_init(pb%energy)

 ! required for quasi-static solution
 ! DEVEL only implemented in elastic material              
  if (pb%time%kind =='quasi-static') then
    write(iout,'(a)') '***********************************************'
    write(iout,'(a)') '*    Finding diagonal of stiffness matrix     *'
    write(iout,'(a)') '***********************************************'
    call MAT_diag_stiffness_init(pb%invKDiag,pb%grid,pb%fields,pb%matwrk)
  endif

end subroutine init_main

!=======================================================================
!
! Checks:
!
! 1. the resolution of the spectral element GRID
!    for a target frequency FMAX
!    The ELASTic properties (VP and VS) are used.
!
! 2. the Courant number

! mode = 1 : SH
!      = 2 : P-SV
  subroutine CHECK_grid(grid,mat,mode,max_c_dx)

  use echo, only : echo_check,iout,fmt1,fmtok
  use stdio, only : IO_new_unit, IO_warning
  use spec_grid, only : sem_grid_type,SE_inquire,SE_elem_coord
  use prop_mat, only : matpro_elem_type, MAT_getProp
  use plot_postscript, only : PLOT_PS

  type(sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  integer, intent(in) :: mode
  double precision, intent(out) :: max_c_dx

  double precision, allocatable :: check1(:),check2(:)
  double precision :: ecoord(2,grid%ngll,grid%ngll) &
                     ,celem(grid%ngll,grid%ngll),csmin,min_cs &
                     ,x0,z0,x1,z1,x2,z2 &
                     ,rdistmax,rdist1,rdist2 &
                     ,rsizemin,rsizemax &
                     ,rlamdaSmin,ratiomax,rlambmin
  integer :: ngll,e,i,j,c1unit,c2unit

  ngll = grid%ngll

  if (echo_check) then

    allocate(check1(grid%nelem),check2(grid%nelem))

    write(iout,*) 
    write(iout,103) ' M e s h   p r o p e r t i e s'
    write(iout,103) ' ============================='
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Checking mesh'

  endif

  rsizemin   = huge(rsizemin)
  rsizemax   = 0d0
  max_c_dx   = 0d0
  rlamdaSmin = huge(rlamdaSmin)
  min_cs = huge(min_cs)

  do e=1,grid%nelem

    call SE_inquire(grid,element=e,size_max=rdistmax)
    call MAT_getProp(celem,mat(e),'cs')
    csmin = minval(celem)
    min_cs = min(csmin,min_cs)
    rlambmin = csmin/rdistmax
    rlamdaSmin = min(rlamdaSmin,rlambmin)

    ratiomax = 0.d0

    if (mode==2) call MAT_getProp(celem,mat(e),'cp')

    ecoord = SE_elem_coord(grid,e)

    do j=1,ngll-1
    do i=1,ngll-1
    
      x0 = ecoord(1,i,j)
      z0 = ecoord(2,i,j)
      x1 = ecoord(1,i+1,j)
      z1 = ecoord(2,i+1,j)
      x2 = ecoord(1,i,j+1)
      z2 = ecoord(2,i,j+1)

      rdist1 = sqrt((x1-x0)**2 + (z1-z0)**2)
      rdist2 = sqrt((x2-x0)**2 + (z2-z0)**2)

      rsizemin = min(rsizemin,rdist1)
      rsizemin = min(rsizemin,rdist2)
      rsizemax = max(rsizemax,rdist1)
      rsizemax = max(rsizemax,rdist2)

      ratiomax = max(ratiomax, celem(i,j)/min(rdist1,rdist2) )

    enddo
    enddo

    max_c_dx = max(max_c_dx,ratiomax)

    if (echo_check) then
      check1(e) = rlambmin ! resolution: nodes per min wavelength
      check2(e) = ratiomax ! stability: Courant number (CFL) assuming dt=1
    endif

  enddo

  if (echo_check) then
   
    write(iout,fmtok)
    write(iout,101) '    Max mesh size = ',rsizemax
    write(iout,101) '    Min mesh size = ',rsizemin
    write(iout,101) '    Ratio max/min = ',rsizemax/rsizemin
    write(iout,*) 
    write(iout,101) '    RESOLUTION: nodes per min wavelength = ',(ngll-1)*rlamdaSmin/grid%fmax
    write(iout,102) '                for maximum frequency   = ',grid%fmax, ' Hz'
    write(iout,102) '                    minimum wavelength  = ',min_cs/grid%fmax, ' m'
    write(iout,*) 

    check1 = check1*(ngll-1)/grid%fmax

    c1unit = IO_new_unit()
    open(c1unit,file="Resolution_sem2d.tab")
    do e=1,grid%nelem
      write(c1unit,'(e10.3)') check1(e)
    enddo
    close(c1unit)

    c2unit = IO_new_unit()
    open(c2unit,file="Stability_sem2d.tab")
    do e=1,grid%nelem
      write(c2unit,'(e10.3)') check2(e)
    enddo
    close(c2unit)

    call PLOT_PS(file="Resolution_sem2d.ps",efield=check1,grid=grid,mat=mat &
             ,stitle='Resolution: nodes per minimum wavelength')

    call PLOT_PS(file="Stability_sem2d.ps",efield=check2,grid=grid,mat=mat &
             ,stitle='Stability: CFL number')

    deallocate(check1,check2)
  endif

  if (ngll*rlamdaSmin/grid%fmax < 5.d0) then
    write(iout,*)
    write(iout,103) '*******************************************************'
    write(iout,102) '** WARNING: Low resolution at fmax = ',grid%fmax,' Hz **'
    write(iout,102) '**          Degraded accuracy is expected !          **' 
    write(iout,103) '** Try a finer mesh OR higher polynomial degree OR   **'
    write(iout,103) '**      use this mesh with lower frequencies         **'
    write(iout,103) '*******************************************************'
    write(iout,*)
    call IO_warning()
  endif

  return
 
  101 format(A,EN12.3)
  102 format(A,EN12.3,A)
  103 format(A)

  end subroutine CHECK_grid



end module init
