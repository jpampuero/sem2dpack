! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
  use spec_grid, only : SE_init
  use fem_grid, only : FE_GetNodesPerElement,FE_GetNbNodes
  use bc_gen, only : BC_init
  use elastic, only : ELAST_init_mass,ELAST_init
  use kelvin_voigt, only : KV_init
  use time_evol, only : TIME_init
  use plot_postscript, only : POST_PS_init
  use receivers, only : REC_init,REC_inquire
  use sources, only : SO_init,SO_check
  use fields_class, only : FIELDS_init,FIELDS_read
  use memory_info
  use echo, only : iout,info=>echo_init, fmt1,fmtok

  type(problem_type), intent(inout) :: pb
  character(50), intent(in), optional :: InitFile
  
  double precision :: grid_cfl
  integer :: recsamp
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
  
 ! define arrays for the computation of elastic internal forces
  call ELAST_init(pb%elast,pb%grid,pb%fields%ndof)

  call POST_PS_init(pb%grid) ! For display
  call CHECK_grid(pb%grid,pb%elast,pb%fields%ndof,grid_cfl)
  call TIME_init(pb%time,grid_cfl) ! define time evolution coefficients

 ! Kelvin-Voigt coefficient may be normalized by dt
  call KV_init(pb%kelvin_voigt,pb%grid%coord,pb%time%dt,pb%fields%ndof)

 ! initialise fields
  if (info) write(iout,fmt1,advance='no') 'Initializing kinematic fields'
  call FIELDS_init(pb%fields,pb%grid%npoin, pb%time%kind=='newmark')
  if(init_cond) call FIELDS_read(pb%fields,InitFile) ! read from a file
  if (info) then
    write(iout,fmtok)
    write(iout,'(7X,A,EN12.3)') 'Max displ = ',maxval(abs(pb%fields%displ))
    write(iout,'(7X,A,EN12.3)') 'Max veloc = ',maxval(abs(pb%fields%veloc))
  endif

 ! build the mass matrix
  if (info) write(iout,fmt1,advance='no') 'Building the mass matrix'
  call ELAST_init_mass(pb%rmass,pb%grid,pb%elast,pb%fields%ndof)
  if (info) write(iout,fmtok)

 ! initialize physical parameters of all boundary conditions
 ! For faults: may reset the initial velocity
  if (info) write(iout,fmt1,advance='no') 'Defining boundary conditions'
  call BC_init(pb%bc,pb%grid,pb%elast,pb%rmass,pb%time,pb%fields,pb%src)
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
    call SO_init(pb%src,pb%grid,pb%elast) 
    if (associated(pb%rec)) then
      call REC_inquire(pb%rec, isamp=recsamp )
    else
      recsamp=1
    endif
    call SO_check(pb%src, pb%time%dt*recsamp, pb%time%nt/recsamp)
    if (info) write(iout,fmtok)
  endif

 !  list memory usage
  call MEMO_echo()

 ! invert the mass matrix
 ! NOTE: crashes here with segmentation fault if pb%rmass exceeds stack size
 !       FIX: (bash) ulimit -s unlimited
  pb%rmass = 1.d0 / pb%rmass

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
  subroutine CHECK_grid(grid,elast,mode,max_c_dx)

  use echo, only : echo_check,iout,fmt1,fmtok
  use stdio, only : IO_new_unit
  use spec_grid, only : sem_grid_type,SE_inquire
  use elastic, only : elast_type,ELAST_inquire
  use plot_postscript, only : PLOT_PS

  type(sem_grid_type), intent(in) :: grid
  type(elast_type)   , intent(in) :: elast
  integer, intent(in) :: mode
  double precision, intent(out) :: max_c_dx

  double precision, allocatable :: check1(:),check2(:)
  double precision :: celem(grid%ngll,grid%ngll),csmin &
                     ,x0,z0,x1,z1,x2,z2 &
                     ,rdistmax,rdist1,rdist2 &
                     ,rsizemin,rsizemax &
                     ,rlamdaSmin,ratiomax,rlambmin &
                     ,dhuge &
                     ,cpcor(4),cscor(4),rhocor(4)
  integer :: e,i,j,c1unit,c2unit,cpunit,csunit,rhounit

  if (echo_check) then

    allocate(check1(grid%nelem),check2(grid%nelem))
    write(iout,fmt1,advance='no') 'Exporting model'

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
      write(cpunit,'(4(e11.4,1x))') cpcor
      write(csunit,'(4(e11.4,1x))') cscor
      write(rhounit,'(4(e11.4,1x))') rhocor
    enddo

    close(cpunit)
    close(csunit)
    close(rhounit)

    write(iout,fmtok)

    write(iout,*) 
    write(iout,103) ' M e s h   p r o p e r t i e s'
    write(iout,103) ' ============================='
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Checking mesh'

  endif

  dhuge = huge(dhuge)

  rsizemin   = dhuge
  rsizemax   = 0d0
  max_c_dx   = 0d0
  rlamdaSmin = dhuge

  do e=1,grid%nelem

    call SE_inquire(grid,element=e,size_max=rdistmax)
    call ELAST_inquire(elast,e,cs=celem)
    csmin = minval(celem)
    rlambmin = csmin/rdistmax
    rlamdaSmin = min(rlamdaSmin,rlambmin)

    ratiomax = 0.d0

    if (mode==2) call ELAST_inquire(elast,e,cp=celem)

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
    write(iout,102) '    RESOLUTION: nodes per min wavelength = ',(grid%ngll-1)*rlamdaSmin/grid%fmax
    write(iout,*) 

    check1 = check1*(grid%ngll-1)/grid%fmax

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

    call PLOT_PS(file="Resolution_sem2d.ps",efield=check1 &
             ,grid=grid,elast=elast &
             ,stitle='Resolution: nodes per minimum wavelength')

    call PLOT_PS(file="Stability_sem2d.ps",efield=check2 &
             ,grid=grid,elast=elast &
             ,stitle='Stability: CFL number')

    deallocate(check1,check2)
  endif

  if (grid%ngll*rlamdaSmin/grid%fmax < 5.d0) then
    write(iout,*)
    write(iout,103) '*******************************************************'
    write(iout,102) '** WARNING: Low resolution at fmax = ',grid%fmax,' Hz **'
    write(iout,102) '**          Degraded accuracy is expected !          **' 
    write(iout,103) '** Try a finer mesh OR higher polynomial degree OR   **'
    write(iout,103) '**      use this mesh with lower frequencies         **'
    write(iout,103) '*******************************************************'
    write(iout,*)
  endif

  return
 
  101 format(A,EN12.3)
  102 format(A,EN12.3,A)
  103 format(A)
  104 format(A,I0)

  end subroutine CHECK_grid



end module init
