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
  use spec_grid, only : SE_init_gll,SE_init_numbering,Q49_init &
                       ,SE_init_coord,SE_init_interpol,SE_BcTopoInit
  use bc_gen, only : BC_init_periodic,BC_init_not_periodic
  use elastic, only : ELAST_init_mass,ELAST_init
  use time_evol, only : TIME_init
  use plot_postscript, only : POST_PS_init
  use receivers, only : REC_init,REC_inquire
  use sources, only : SO_init,SO_check
  use fields_class, only : FIELDS_init,FIELDS_read,FIELDS_get_max
  use memory_info
  use echo, only : iout,info=>echo_init, fmt1,fmtok
  use stdio, only : IO_new_unit

  type(problem_type), intent(inout) :: pb
  character(50), intent(in), optional :: InitFile
  
!  integer :: ngll,nspec
  integer :: recsamp,ounit
  logical :: init_cond

  init_cond = present(InitFile)

  if (info) write(iout,fmt1,advance='no') 'Defining the FEM mesh'
  call MESH_build(pb%grid,pb%mesh) ! build the Finite Elements mesh
  if (info) write(iout,fmtok)

  call SE_init_gll(pb%grid)  ! get GLL quadrature points, weights 
                              ! and lagrange coefficients
  call SE_init_numbering(pb%grid) ! generate the global node numbering
  call Q49_init(pb%grid) ! shape functions and jacobian

  call SE_init_coord(pb%grid) ! get the coordinates of the GLL nodes,
                                     ! and export the grid to a file

 ! export grid parameters
  ounit = IO_new_unit()
  open(ounit,file='grid_sem2d.hdr',status='replace')
  write(ounit,*) 'NELEM  NPGEO  NGNOD  NPOIN  NGLL'
  write(ounit,*) pb%grid%nelem, pb%grid%npgeo, pb%grid%ngnod, pb%grid%npoin, pb%grid%ngll
  close(ounit)
  
  if (info) write(iout,fmt1,advance='no') 'Defining boundaries'
  call SE_BcTopoInit(pb%grid) ! initialize boundaries generic data
  if (info) write(iout,fmtok)

 ! define arrays for the computation of elastic internal forces
  call ELAST_init(pb%elast,pb%grid)

  call POST_PS_init(pb%grid) ! For display
 ! verifier le maillage, la stabilite et le nb de points par lambda
  call CHECK_grid(pb%grid,pb%elast,pb%time)
  call TIME_init(pb%time) ! define time evolution coefficients

 ! build the mass matrix for spectral elements
  if (info) write(iout,fmt1,advance='no') 'Defining and inverting the mass matrix'
  call ELAST_init_mass(pb%rmass,pb%grid,pb%elast)
  call BC_init_periodic(pb%bc,pb%grid,pb%rmass)
  pb%rmass = 1.d0 / pb%rmass ! invert the mass matrix
  if (info) write(iout,fmtok)

!   deallocate(pb%grid%weights) 
 ! initialize physical parameters of all boundary conditions
 ! WARNING: this will eventually apply periodicity to the mass matrix
 !          => it is safer to declare first all periodic boundaries
  if (info) write(iout,fmt1,advance='no') 'Defining physical properties of boundaries'
  call BC_init_not_periodic(pb%bc,pb%grid,pb%elast,pb%rmass,pb%time) 
  if (info) write(iout,fmtok)

 ! initialise fields
  if (info) write(iout,fmt1,advance='no') 'Initializing kinematic fields'
  call FIELDS_init(pb%fields,pb%grid%npoin)
  ! eventually read initial conditions from a file
  if(init_cond) call FIELDS_read(pb%fields,InitFile)
  if (info) write(iout,fmtok)
 ! afficher le max du deplacement initial
  if (info) then
    write(iout,'(7X,A,EN12.3)') 'Max displ = ',FIELDS_get_max(pb%fields%displ)
    write(iout,'(7X,A,EN12.3)') 'Max veloc = ',FIELDS_get_max(pb%fields%veloc)
  endif

 ! define position of receivers and allocate database
  if (associated(pb%rec)) then
    if (info) write(iout,fmt1,advance='no') 'Initializing receivers'
    call REC_init(pb%rec,pb%grid,pb%time,pb%fields)
    if (info) write(iout,fmtok)
  endif

 ! initialize sources
  if (associated(pb%src)) then
    if (info) write(iout,fmt1,advance='no') 'Initializing sources'
    call SO_init(pb%src,pb%grid,pb%fields%accel,pb%fields,pb%elast) 
    call REC_inquire(pb%rec, isamp=recsamp )
    call SO_check(pb%src, pb%time%dt*recsamp, pb%time%nt/recsamp)
    if (info) write(iout,fmtok)
  endif

!  deallocate(pb%grid%xjaci) ! needed by explosion sources 

! declare working space, only if ELAST_KD1 is used in the solver
!  ngll = pb%grid%ngll
!  nspec = pb%grid%nelem
!  call storearray('dUx_dxi',ngll*ngll*nspec,idouble)
!  call storearray('dUz_dxi',ngll*ngll*nspec,idouble)
!  call storearray('dUx_deta',ngll*ngll*nspec,idouble)
!  call storearray('dUz_deta',ngll*ngll*nspec,idouble)
!  call storearray('Uxloc',ngll*ngll*nspec,idouble)
!  call storearray('Uzloc',ngll*ngll*nspec,idouble)
!  call storearray('t1',ngll*ngll*nspec,idouble)
  
 !  list a short directory after input phase
  call MEMO_echo()

end subroutine init_main

!=======================================================================
!
! Checks:
!
! 1. the resolution of the spectral element GRID
!    for a target frequency FMAX
!    The ELASTic properties (VP and VS) are used.
!
! 2. if timestep was set by the user: the Courant number is checked 
!    if a target Courant number was set by the user: the timestep is fixed

  subroutine CHECK_grid(grid,elast,time)

  use echo, only : echo_check,iout
  use stdio, only : IO_new_unit
  use spec_grid, only : sem_grid_type,SE_inquire
  use elastic, only : elast_type,ELAST_csmin,ELAST_inquire
  use time_evol, only : timescheme_type
  use plot_postscript, only : PLOT_PS

  type(sem_grid_type), intent(in) :: grid
  type(elast_type)   , intent(in) :: elast
  type(timescheme_type), intent(inout) :: time

  double precision, allocatable :: check1(:),check2(:)
  double precision :: cploc,csmin &
                     ,x0,z0,x1,z1,x2,z2,x3,z3 &
                     ,rdistmin,rdistmax,rdist1,rdist2 &
                     ,rsizemin,rsizemax,max_cp_dx &
                     ,rlamdaSmin,ratiomax,rlambmin,rlambmax &
                     ,dhuge &
                     ,cpcor(4),cscor(4),rhocor(4)
  integer :: e,i,j,c1unit,c2unit,cpunit,csunit,rhounit

  if (echo_check) then
    allocate(check1(grid%nelem),check2(grid%nelem))
    write(iout,'(5X,A)',advance='no') 'Exporting model ...'

    cpunit = IO_new_unit()
    open(cpunit,file='Cp_sem2d.tab')
    csunit = IO_new_unit()
    open(csunit,file='Cs_sem2d.tab')
    rhounit = IO_new_unit()
    open(rhounit,file='Rho_sem2d.tab')

    do e=1,grid%nelem
      call ELAST_inquire(elast,1,1,e,cp=cpcor(1),cs=cscor(1),rho=rhocor(1))
      call ELAST_inquire(elast,grid%ngll,1,e,cp=cpcor(2),cs=cscor(2),rho=rhocor(2))
      call ELAST_inquire(elast,grid%ngll,grid%ngll,e,cp=cpcor(3),cs=cscor(3),rho=rhocor(3))
      call ELAST_inquire(elast,1,grid%ngll,e,cp=cpcor(4),cs=cscor(4),rho=rhocor(4))
      write(cpunit,*) cpcor
      write(csunit,*) cscor
      write(rhounit,*) rhocor
    enddo

    close(cpunit)
    close(csunit)
    close(rhounit)

    write(iout,'(A)') '... [OK]'

  endif

  if (echo_check) write(iout,'(5X,A)',advance='no') 'Checking mesh ...'
  dhuge = huge(dhuge)

  rsizemin     =  dhuge
  rsizemax     = -dhuge
  max_cp_dx    = -dhuge
  rlamdaSmin   =  dhuge

  do e=1,grid%nelem

    call SE_inquire(grid,element=e,size_max=rdistmax,size_min=rdistmin)
    call ELAST_csmin(elast,e,csmin)
    rlambmin = csmin/rdistmax
    rlamdaSmin = min(rlamdaSmin,rlambmin)

    ratiomax = 0.d0

    do j=1,grid%ngll-1
    do i=1,grid%ngll-1
    
      x0 = grid%coord(1,grid%ibool(i,j,e))
      z0 = grid%coord(2,grid%ibool(i,j,e))
      x1 = grid%coord(1,grid%ibool(i+1,j,e))
      z1 = grid%coord(2,grid%ibool(i+1,j,e))
      x2 = grid%coord(1,grid%ibool(i,j+1,e))
      z2 = grid%coord(2,grid%ibool(i,j+1,e))

      rdist1 = sqrt((x1-x0)**2 + (z1-z0)**2)
      rdist2 = sqrt((x2-x0)**2 + (z2-z0)**2)

      rsizemin = min(rsizemin,rdist1)
      rsizemin = min(rsizemin,rdist2)
      rsizemax = max(rsizemax,rdist1)
      rsizemax = max(rsizemax,rdist2)

      call ELAST_inquire(elast,i,j,e,cp=cploc)
      ratiomax = max(ratiomax, cploc / min(rdist1,rdist2))

    enddo
    enddo

    max_cp_dx = max(max_cp_dx,ratiomax)

    if (echo_check) then
      check1(e) = rlambmin ! resolution: nodes per min wavelength
      check2(e) = ratiomax ! stability: Courant number (CFL) assuming dt=1
    endif

  enddo
  if (echo_check) write(iout,'(A)') '... [OK]'

 ! Check the Courant number or set the timestep:
  if (time%dt > 0.d0) then
    time%courant =  max_cp_dx*time%dt 
  else
    time%dt      =  time%courant/max_cp_dx
   ! Set the total duration or the number of timesteps
    if (time%total > 0.d0) time%nt = int(time%total/time%dt)
    time%total = time%nt*time%dt
  endif

  if (echo_check) then
   
    write(iout,*) 
    write(iout,103) '==================================================='
    write(iout,103) '=          Time solver settings                   =' 
    write(iout,103) '==================================================='
    write(iout,*) 
    write(iout,101) '    Time step (secs)      = ',time%dt
    write(iout,104) '    Number of time steps  = ',time%nt
    write(iout,101) '    Total duration (secs) = ',time%total
    write(iout,101) '    Courant number        = ',time%courant
    write(iout,*) 

    write(iout,*) 
    write(iout,103) '==================================================='
    write(iout,103) '=          Mesh properties checkboard             ='
    write(iout,103) '==================================================='
    write(iout,*) 
    write(iout,101) '    Max mesh size = ',rsizemax
    write(iout,101) '    Min mesh size = ',rsizemin
    write(iout,101) '    Ratio max/min = ',rsizemax/rsizemin
    write(iout,*) 
    write(iout,102) '    STABILITY:  CFL number               = ',max_cp_dx*time%dt
    write(iout,*) 
    write(iout,102) '    RESOLUTION: nodes per min wavelength = ',(grid%ngll-1)*rlamdaSmin/grid%fmax
    write(iout,*) 

    check1 = check1*(grid%ngll-1)/grid%fmax
    check2 = check2*time%dt

    c1unit = IO_new_unit()
    open(c1unit,file="Resolution_sem2d.tab")
    c2unit = IO_new_unit()
    open(c2unit,file="Stability_sem2d.tab")
    do e=1,grid%nelem
      write(c1unit,*) check1(e)
      write(c2unit,*) check2(e)
    enddo
    close(c1unit)
    close(c2unit)
    call PLOT_PS(efield=check1,file="Resolution_sem2d.ps" &
             ,grid=grid,elast=elast &
             ,stitle='Resolution: nodes per minimum wavelength')
    call PLOT_PS(efield=check2,file="Stability_sem2d.ps" &
             ,grid=grid,elast=elast &
             ,stitle='Stability: CFL number')
    deallocate(check1,check2)
  endif

  if (grid%ngll*rlamdaSmin/grid%fmax < 5.d0) then
    write(iout,103) '***************************************************************'
    write(iout,102) '** WARNING: poor grid resolution at fmax = ',grid%fmax,' Hz **'
    write(iout,103) '** Try a finer mesh OR higher polynomial degree OR   **'
    write(iout,103) '**      use this mesh with lower frequencies         **'
    write(iout,103) '***************************************************************'
  endif

  if (time%courant > 0.55d0) then
    write(iout,103) '***************************************************************'
    write(iout,102) '** WARNING: Courant number is too high = ',time%courant,' ******'
    write(iout,102) '** Try a timestep smaller than ',time%dt*0.55d0/time%courant &
                   ,'       **'
    write(iout,103) '***************************************************************'

  endif
    
  return
 
  101 format(A,EN12.3)
  102 format(A,EN12.3,A)
  103 format(A)
  104 format(A,I0)

  end subroutine CHECK_grid



end module init
