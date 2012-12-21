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
module plot_postscript

  use spec_grid 

  implicit none
  private

  double precision, parameter :: centim = 28.5d0
  logical, parameter :: legendes=.true.  !--- ecrire legendes ou non

  integer, parameter :: maxcolors = 20
  real, dimension(3,maxcolors), save :: RGB

  logical, save :: interpol,vectors,mesh,symbols, &
                   cpmodel,boundaries,usletter,color,numbers

  type(interpol_type), pointer , save :: interp

  integer, save :: isubsamp,DisplayPts
                   
  double precision, save :: cutvect,scalex,scalez,ArrowHeadAngle,ArrowHBRatio, &
      sizex,sizez,PageX0,PageZ0,rapp_page, &
      sizemax,xmin,zmin,xmax,zmax,ScaleField

  real, save :: usoffset

  integer, parameter :: ndime = 2
  
  public :: PLOT_PS,POST_PS_read,POST_PS_init

contains

!=======================================================================
!
!
! BEGIN INPUT BLOCK
!
! NAME   : PLOTS_POSTCRIPT [plots]
! PURPOSE: Preferences for PostScript snapshots
! SYNTAX : &PLOTS_POSTSCRIPT vectors, mesh, CpModel, color,
!               isubsamp, boundaries, symbols, numbers,
!               ScaleField, CutVect, SizeMax, 
!               ArrowHeadAngle, ArrowHBRatio,
!               Interpol, DisplayPts,
!               ScaleX, ScaleZ, USLetter, PageX0, PageZ0 /
!
! ARG: vectors          [log] [F] Plots a vectorial field (with arrows)
! ARG: mesh             [log] [F] Plots the mesh on background
! ARG: CpModel          [log] [F] Plots the P-velocity model as background
! ARG: color            [log] [F] Color output
! ARG: isubsamp         [int] [0] Subsampling of the GLL nodes for the
!                                 output of velocity model. The default uses 
!                                 vertex and middle point values only. 
!                                 isubsamp=3 samples every 3 GLL point, etc
! ARG: boundaries       [log] [T] Colors every tagged boundary
! ARG: symbols          [log] [T] Plots symbols for sources and receivers
! ARG: numbers          [log] [F] Plots the element numbers
! ARG: ScaleField       [dble] [0d0] Scale factor for the field amplitude. The
!                       default scales each snapshot by its maximum amplitude.
!                       This is not convenient for comparing snapshots or 
!                       running movies
! ARG: CutVect          [dble] [1.d0] Minimum arrow length as a percentage of 
!                       scale factor
! ARG: SizeMax          [dble] [1.d0] Maximum arrow length (in cm)
! ARG: ArrowHeadAngle   [dble] [20d0] Angle between arrow body and head (degrees)
! ARG: ArrowHBRatio     [dble] [0.4d0] Arrow head to body ratio 
! ARG: Interpol         [log] [T] Interpolate field. Subsampling allows smaller 
!                       snapshot files.
! ARG: DisplayPts       [log] [3] Number of points/element/direction for 
!                       interpolation. The default subsamples using only vertex 
!                       and middle point values
! ARG: ScaleX           [dble] [1d0] Scale factor in the X direction
! ARG: ScaleZ           [dble] [1d0] Scale factor in the Z direction
! ARG: USLetter         [log] [T] US letter format or French A4
! ARG: PageX0           [dble] [2.2d0] X-origin in the page (cm)
! ARG: PageZ0           [dble] [2.2d0] Z-origin in the page (cm)
!               
! END INPUT BLOCK

  subroutine POST_PS_read(iin)

  use echo, only : iout,echo_input

  integer, intent(in) :: iin

  NAMELIST / PLOTS_POSTSCRIPT / vectors,numbers &
                         ,cpmodel,isubsamp,color &
                         ,boundaries &
                         ,symbols &
                         ,mesh,cutvect,sizemax,ArrowHeadAngle &
                         ,interpol,DisplayPts &
                         ,scalex,scalez,ScaleField &
                         ,ArrowHBRatio,usletter,PageX0,PageZ0

  vectors     = .false. 
  numbers       = .false. 
  cpmodel    = .false.
  isubsamp      = 2
  color         = .false.
  boundaries    = .true.
  symbols      = .true.
  mesh          = .false.
  cutvect       = 1.d0
  sizemax       = 1.d0
  ArrowHeadAngle         = 20.d0
  interpol      = .true.
  DisplayPts      = 3
  scalex        = 1.d0
  scalez        = 1.d0
  ScaleField   = 0.d0 
  ArrowHBRatio       = 0.4d0
  usletter      = .true.
  PageX0        = 2.2d0
  PageZ0        = 2.2d0
     
  rewind(iin)
  read(iin,PLOTS_POSTSCRIPT,END=100)
100 continue
 
  cutvect = cutvect / 100.d0

  if (echo_input) write(iout,200) mesh,numbers &
                         ,cpmodel,isubsamp,color &
                         ,boundaries &
                         ,symbols &
                         ,vectors,100.d0*cutvect,sizemax,ArrowHeadAngle &
                         ,scalex,scalez,ScaleField &
                         ,ArrowHBRatio,PageX0,PageZ0,usletter &
                         ,interpol,DisplayPts
  return

  200 format(//1x,'P o s t S c r i p t   C o n t r o l   c a r d',/1x,34('='),//5x, &
  'Plot mesh . . . . . . . . . . . . . . . . . . (mesh) = ',L1/ 5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',L1/ 5x, &
  'Plot model  . . . . . . . . . . . . . . (cpmodel) = ',L1/ 5x, &
  'Subsampling for velocity model display  . (isubsamp) = ',I0/5x, &
  'Color display . . . . . . . . . . . . . . . .(color) = ',L1/ 5x, &
  'Plot boundaries . . . . . . . . . . . . (boundaries) = ',L1/ 5x, &
  'Plot symbols  . . . . . . . . . . . . . . (symbols) = ',L1/ 5x, &
  'Plot vector fields  . . . . . . . . . . . .(vectors) = ',L1/ 5x, &
  'Percentage of cut for vector plots. . . . .(cutvect) = ',F0.2/5x, &
  'Max size of arrows. . . . . . . . . . . . .(sizemax) = ',F0.2/5x, &
  'Angle of vector arrows. . . . . . . . . . . .(ArrowHeadAngle) = ',F0.2/5x, &
  'X-Scaling . . . . . . . . . . . . . . . . . (scalex) = ',F0.2/5x, &
  'Z-Scaling . . . . . . . . . . . . . . . . . (scalez) = ',F0.2/5x, &
  'A-Scaling . . . . . . . . . . . . . . .(ScaleField) = ',F0.2/5x, &
  'Head to body ratio for arrows . . . . . . .(ArrowHBRatio) = ',F0.2/5x, &
  'X origin. . . . . . . . . . . . . . . . . . (PageX0) = ',F0.2/5x, &
  'Z origin. . . . . . . . . . . . . . . . . . (PageZ0) = ',F0.2/5x, &
  'US letter format or French A4 . . . . . . (usletter) = ',L1/5x, &
  'Interpolate vector field  . . . . . . . . (interpol) = ',L1/5x, &
  'Points per edge for interpolation . . . . (DisplayPts) = ',I0)


  end subroutine POST_PS_read

!=======================================================================
! 
  subroutine POST_PS_init(grid)

  use spec_grid, only : sem_grid_type,SE_init_interpol

  type(sem_grid_type), intent(in) :: grid

  !-- taille de la fenetre de display Postscript en pourcentage de la feuille
  double precision, parameter :: rpercentx = 70.0d0, rpercentz = 77.0d0

  !-- definition de la palette de couleur

  RGB(:,1)  = (/ 1.00, 0.00, 0.00 /) ! red
  RGB(:,2)  = (/ 0.00, 0.00, 1.00 /) ! blue
  RGB(:,3)  = (/ 0.93, 0.51, 0.93 /) ! violet
  RGB(:,4)  = (/ 0.73, 0.33, 0.83 /) ! medium orchid
  RGB(:,5)  = (/ 0.60, 0.20, 0.80 /) ! dark orchid
  RGB(:,6)  = (/ 0.54, 0.17, 0.89 /) ! blue violet
  RGB(:,7)  = (/ 0.42, 0.35, 0.80 /) ! slate blue
  RGB(:,8)  = (/ 1.00, 0.08, 0.58 /) ! deep pink
  RGB(:,9)  = (/ 0.12, 0.56, 1.00 /) ! dodger blue
  RGB(:,10) = (/ 0.00, 0.81, 0.82 /) ! dark turquoise
  RGB(:,11) = (/ 0.25, 0.88, 0.82 /) ! turquoise
  RGB(:,12) = (/ 0.20, 0.80, 0.20 /) ! lime green
  RGB(:,13) = (/ 0.00, 1.00, 0.50 /) ! spring green
  RGB(:,14) = (/ 0.50, 1.00, 0.00 /) ! chartreuse
  RGB(:,15) = (/ 0.68, 1.00, 0.18 /) ! green yellow
  RGB(:,16) = (/ 1.00, 1.00, 0.00 /) ! yellow
  RGB(:,17) = (/ 1.00, 0.98, 0.80 /) ! lemon chiffon
  RGB(:,18) = (/ 1.00, 0.84, 0.00 /) ! gold
  RGB(:,19) = (/ 1.00, 0.89, 0.71 /) ! mocassin
  RGB(:,20) = (/ 1.00, 0.85, 0.73 /) ! peach puff

  ! papier A4 ou US letter
  if(usletter) then
    usoffset = 1.75
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.
    sizex = 29.7d0
    sizez = 21.d0
  endif

  !-- grid extrema
  xmax=maxval(grid%coord(1,:))
  xmin=minval(grid%coord(1,:))
  zmax=maxval(grid%coord(2,:))
  zmin=minval(grid%coord(2,:))

  ! ratio taille page/taille domaine physique
  rapp_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin))/100.d0

 ! weigths for interpolation to more regular grid
  if (interpol) call SE_init_interpol(interp,grid,DisplayPts) 

 
  end subroutine POST_PS_init


!=======================================================================
! routine affichage postscript
!
  subroutine PLOT_PS(vfield,efield,grid,elast,stitle,file &
                    ,it_in,time_in,src,rec)

  use stdio, only : IO_new_unit
  use echo, only : echo_run,iout
  use elastic, only : elast_type
  use fields_class, only : FIELDS_get_max
  use sources, only : source_type,SO_inquire
  use receivers, only : rec_type,REC_inquire

  type(sem_grid_type), intent(in) :: grid
  type(elast_type)   , intent(in) :: elast
  character(*)       , intent(in) :: stitle
  double precision   , intent(in), optional :: vfield(grid%npoin,2)
  double precision   , intent(in), optional :: efield(grid%nelem)
  character(*)       , intent(in), optional :: file 
  integer            , intent(in), optional :: it_in
  double precision   , intent(in), optional :: time_in
  type(source_type)  , optional, pointer :: src
  type(rec_type)     , optional, pointer :: rec

 ! hauteur des numeros de domaine en CM
  double precision, parameter :: height = 0.25d0 

  double precision :: maxfield,time
  integer :: psunit,it
  character(50) :: name

  if (present(time_in) .and. present(it_in)) then
    time = time_in
    it = it_in
  else
    time = 0d0
    it = 0
  endif
  if (present(vfield)) maxfield = FIELDS_get_max(vfield)
  if (present(efield)) maxfield = maxval(abs(efield))

! call system('cp .ps_header vect....ps')
! open(psunit,file=vect...ps,position='append')

  if (present(vfield)) write(name,'( "Snapshot_",I5.5,"_sem2d.ps" )') it
  if (present(file)) name = file

  psunit = IO_new_unit()
  open(unit=psunit,file=name,status='unknown')
  if (echo_run) write(iout,'("Dump PostScript ",A," ...")',advance='no') trim(name)

  call plot_header()
  
  write(psunit,*) '0 setgray'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.5 CM SCSF'

  if (legendes) call plot_legend()

! close(psunit)
! call system('cat .ps_model >> vect...ps' )
! open(psunit,file=vect...ps,position='append')

  write(psunit,*) '%'
  write(psunit,*) scalex,' ',scalez,' scale'
  write(psunit,*) '%'

  if (cpmodel)     call plot_model()
  if (vectors)      call plot_mesh()
  call plot_boundaries()

  if (present(vfield) .and. maxfield > 0.d0) call plot_vect()
  if (present(efield)) call plot_efield()

! close(psunit)
! call system('cat .ps_trailer >> vect...ps' )

  if(cpmodel) then ! Sources and receivers in color ?
    write(psunit,*) 'Colreceiv'
  else
    write(psunit,*) '0 setgray'
  endif

  if (present(src) .and. associated(src)) call plot_sources()
  if (present(rec) .and. associated(rec)) call plot_receivers()

  write(psunit,*) '%'
  write(psunit,*) 'grestore'
  write(psunit,*) 'showpage'

  close(psunit)
  if (echo_run) write(iout,'(A)') '... [OK]'

  return
 
  contains 


!-----------------------------------------------------------------------
! Plot headers

  subroutine plot_header()

  write(psunit,10) stitle
  write(psunit,*) '/CM {28.5 mul} def'
  write(psunit,*) '/LR {rlineto} def'
  write(psunit,*) '/LT {lineto} def'
  write(psunit,*) '/L {lineto} def'
  write(psunit,*) '/MR {rmoveto} def'
  write(psunit,*) '/MV {moveto} def'
  write(psunit,*) '/M {moveto} def'
  write(psunit,*) '/MK {mark} def'
  write(psunit,*) '/ST {stroke} def'
  write(psunit,*) '/CP {closepath} def'
  write(psunit,*) '/RG {setrgbcolor} def'
  write(psunit,*) '/GF {gsave fill grestore} def'
  write(psunit,*) '/GG {0 setgray ST} def'
  write(psunit,*) '/GC {Colmesh ST} def'
  write(psunit,*) '/RF {setrgbcolor fill} def'
  write(psunit,*) '/SF {setgray fill} def'
  write(psunit,*) '/GS {gsave} def'
  write(psunit,*) '/GR {grestore} def'
  write(psunit,*) '/SLW {setlinewidth} def'
  write(psunit,*) '/SCSF {scalefont setfont} def'
  write(psunit,*) '% differents symboles utiles'
  write(psunit,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(psunit,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(psunit,*) 'CP fill} def'
  write(psunit,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(psunit,*) 'CP fill} def'
  write(psunit,*) '/Cross {GS 0.05 CM SLW'
  write(psunit,*) 'GS 3 3 MR -6. -6. LR ST GR'
  write(psunit,*) 'GS 3 -3 MR -6. 6. LR ST GR'
  write(psunit,*) '0.01 CM SLW} def'
  write(psunit,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(psunit,*) '/Losange {GS 0.05 CM SLW 0 4.2 MR'
  write(psunit,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(psunit,*) 'GR 0.01 CM SLW} def'
  write(psunit,*) '%'
  write(psunit,*) '% niveaux de gris pour le modele de vitesse'
  write(psunit,*) '/BK {setgray fill} def'
  write(psunit,*) '% version noir et blanc'
  write(psunit,*) '%/BK {pop 1 setgray fill} def'
  write(psunit,*) '%'
  write(psunit,*) '% magenta pour les vecteurs deplacement'
  write(psunit,*) '/Colvects {0.01 CM SLW 1. 0. 1. RG} def'
  write(psunit,*) '% version noir et blanc'
  write(psunit,*) '%/Colvects {0.01 CM SLW 0. setgray} def'
  write(psunit,*) '%'
  write(psunit,*) '% chartre pour le maillage des macroblocs'
  write(psunit,*) '/Colmesh {0.02 CM SLW 0.5 1. 0. RG} def'
  write(psunit,*) '% version noir et blanc'
  write(psunit,*) '%/Colmesh {0.02 CM SLW 0. setgray} def'
  write(psunit,*) '%'
  write(psunit,*) '% cyan pour les sources et recepteurs'
  write(psunit,*) '/Colreceiv {0. 1. 1. RG} def'
  write(psunit,*) '% version noir et blanc'
  write(psunit,*) '%/Colreceiv {0. setgray} def'
  write(psunit,*) '%'
  write(psunit,*) '% macro dessin fleche'
  write(psunit,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(psunit,*) '% macro dessin contour elements'
  write(psunit,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(psunit,*) '%'
  write(psunit,*) '.01 CM SLW'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.35 CM SCSF'
  write(psunit,*) '%'
  write(psunit,*) '/vshift ',-height/2,' CM def'
  write(psunit,*) '/Rshow { currentpoint stroke MV'
  write(psunit,*) 'dup stringwidth pop neg vshift MR show } def'
  write(psunit,*) '/Cshow { currentpoint stroke MV'
  write(psunit,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(psunit,*) '/fN {/Helvetica-Bold findfont ',height,' CM SCSF} def'
  write(psunit,*) '%'
  write(psunit,*) 'gsave newpath 90 rotate'
  write(psunit,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(psunit,*) '%'

  return

10   format('%!PS-Adobe-2.0',/,'%%',/,'%%Title: ',A,/, &
        '%%Creator: SEM2DPACK Version 2.2',/, &
        '%%Author: Jean-Paul Ampuero',/, &
        '%%BoundingBox: 0 0 612 792',/,'%%')

  end subroutine plot_header


!-----------------------------------------------------------------------
! Plot legends

  subroutine plot_legend()

  write(psunit,*) '24. CM 1.2 CM MV'
  write(psunit,610) usoffset,it
  write(psunit,*) '%'

  write(psunit,*) '24. CM 1.95 CM MV'
  write(psunit,600) usoffset,time
  write(psunit,*) '%'
  write(psunit,*) '24. CM 2.7 CM MV'
  write(psunit,640) usoffset,maxfield
  write(psunit,*) '%'
  write(psunit,*) '24. CM 3.45 CM MV'
  write(psunit,620) usoffset,cutvect*100.d0

  write(psunit,*) '%'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.6 CM SCSF'
  if (color) write(psunit,*) '.4 .9 .9 setrgbcolor'
  write(psunit,*) '11 CM 1.1 CM MV'
  write(psunit,*) '(X axis) show'
  write(psunit,*) '%'
  write(psunit,*) '1.4 CM 9.5 CM MV'
  write(psunit,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(psunit,*) '(Y axis) show'
  write(psunit,*) 'grestore'
  write(psunit,*) '%'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.7 CM SCSF'
  if (color) write(psunit,*) '.8 0 .8 setrgbcolor'
  write(psunit,*) '25.35 CM 18.9 CM MV'
  write(psunit,*) usoffset,' CM 2 div neg 0 MR'
  write(psunit,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(psunit,*) '(',stitle,') show'
  write(psunit,*) 'grestore'
  write(psunit,*) '26.45 CM 18.9 CM MV'
  write(psunit,*) usoffset,' CM 2 div neg 0 MR'
  write(psunit,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(psunit,*) '(Elastic Wave 2D - Spectral Elements Method) show'
  write(psunit,*) 'grestore'

  return

 600  format(F0.3,' neg CM 0 MR (Time =',EN12.3,' s) show')
 610  format(F0.3,' neg CM 0 MR (Time step = ',I0,') show')
 620  format(F0.3,' neg CM 0 MR (Cut =',F0.2,' \%) show')
 640  format(F0.3,' neg CM 0 MR (Max =',EN12.3,') show')

  end subroutine plot_legend

!-----------------------------------------------------------------------
! Draw the velocity model in background

  subroutine plot_model

  use elastic, only : elast_type,ELAST_inquire,ELAST_cpminmax

  double precision :: vpmax,vpmin,cploc,x1
  integer :: i,j,ispec
  logical :: heterogeneous

  call ELAST_cpminmax(elast,vpmin,vpmax)

  heterogeneous = (vpmax-vpmin)/vpmin > 0.02d0
  if (.not. heterogeneous) x1 = 0.5d0

  do ispec=1,grid%nelem
  do i=1,grid%ngll-isubsamp,isubsamp
  do j=1,grid%ngll-isubsamp,isubsamp

    write(psunit,500) point_scaled( grid%coord(:,grid%ibool(i,j,ispec)) )
    write(psunit,499) point_scaled( grid%coord(:,grid%ibool(i+isubsamp,j,ispec)) )
    write(psunit,499) point_scaled( grid%coord(:,grid%ibool(i+isubsamp,j+isubsamp,ispec)) )
    write(psunit,499) point_scaled( grid%coord(:,grid%ibool(i,j+isubsamp,ispec)) )

    if (heterogeneous) then
      call ELAST_inquire(elast,i,j,ispec,cp=cploc)
      x1 = (cploc-vpmin)/(vpmax-vpmin)
      x1 = min( x1*0.7 + 0.2 , 1.d0 ) ! rescaler pour eviter gris trop sombre
      x1 = 1.d0 - x1 ! inverser echelle : blanc = vpmin, gris = vpmax
    endif

    write(psunit,604) x1

  enddo
  enddo
  enddo

 499  format(F0.2,1x,F0.2,' L')
 500  format(F0.2,1x,F0.2,' M')
 604  format('CP ',F0.2,' BK')

  end subroutine plot_model

!-----------------------------------------------------------------------
!-- Draw spectral element mesh

  subroutine plot_mesh

  double precision :: coorg(ndime,grid%ngnod) &
                     ,coord(ndime,grid%ngll,grid%ngll) &
                     ,point(ndime)
  integer :: e,i,j,is,ir,imat,icol

  write(psunit,*) '%'
  write(psunit,*) '% spectral element mesh'
  write(psunit,*) '%'

  do e=1,grid%nelem

    write(psunit,*) '% elem ',e
    write(psunit,*) 'MK'

    if (grid%ngnod == 4) then ! tracer des droites si elements Q4
   
     ! get the coordinates of the control nodes
      coorg = grid%coorg(:,grid%knods(:,e)) 

      do i=1,grid%ngnod
        write(psunit,601) point_scaled( coorg(:,i) )
      enddo
      write(psunit,601) point_scaled( coorg(:,1) )
    
    else ! tracer des courbes si elements Q9

     ! get the coordinates of the GLL points
      do j=1,grid%ngll
      do i=1,grid%ngll
        coord(:,i,j) = grid%coord(:,grid%ibool(i,j,e))
      enddo
      enddo

      is=1
      do ir=1,grid%ngll
        write(psunit,601) point_scaled( coord(:,ir,is) )
      enddo

      ir=grid%ngll
      do is=2,grid%ngll
        write(psunit,601) point_scaled( coord(:,ir,is) )
      enddo
        
      is=grid%ngll
      do ir=grid%ngll-1,1,-1
        write(psunit,601) point_scaled( coord(:,ir,is) )
      enddo
        
      ir=1
      do is=grid%ngll-1,2,-1
        write(psunit,601) point_scaled( coord(:,ir,is) )
      enddo
      
    endif

    write(psunit,*) 'CO'

    if (color) then

     ! Use a different color for each material set
      imat = grid%tag(e)
      icol = mod(imat - 1,maxcolors) + 1
      write(psunit,600) RGB(:,icol)

    endif

    if(cpmodel) then
      write(psunit,*) 'GC'
    else
      write(psunit,*) 'GG'
    endif

    ! write the element number, the group number and the
    ! material number inside the element
    if (numbers) then

      point = 0.25d0 * SUM( grid%coorg(:,grid%knods(1:4,e)) , dim =2 )
      point = point_scaled(point)
      if (color) write(psunit,*) '1 setgray'
      write(psunit,500) point
      write(psunit,502) e !-- write element number

    endif

  enddo

 500  format(F0.2,1x,F0.2,' M')
 502  format('fN (',I0,') Cshow')
 600  format(F0.2,1x,F0.2,1x,F0.2,' RG GF')
 601  format(F0.2,1x,F0.2)

  end subroutine plot_mesh

!-----------------------------------------------------------------------
!--  Draw the boundaries
   
  subroutine plot_boundaries

  double precision :: point1(ndime),point2(ndime)
  integer :: i,bce,ideb,iend

  write(psunit,*) '%'
  write(psunit,*) '% boundaries of the mesh'
  write(psunit,*) '%'

  write(psunit,*) '0.05 CM SLW'

  do i = 1,size(grid%bounds(:))
  ! set color cyclically from the palette
    write(psunit,'( 3(F0.2,1x),"RG" )') RGB(:, mod(grid%bounds(i)%tag-1,maxcolors)+1 )
  ! draw a straight segment for each boundary element
    do bce = 1,grid%bounds(i)%nelem
      ideb   = grid%bounds(i)%bulk_node( grid%bounds(i)%ibool(1,bce) )
      point1 = point_scaled( grid%coord(:,ideb) )
      iend   = grid%bounds(i)%bulk_node( grid%bounds(i)%ibool(grid%ngll,bce) )
      point2 = point_scaled( grid%coord(:,iend) )
      write(psunit,'( 2(2(F0.2,1x),A) )') point1,'M ',point2,'L ST'
    enddo
  enddo

  write(psunit,*) '0 setgray'
  write(psunit,*) '0.01 CM SLW'

  end subroutine plot_boundaries


!---------------------------------------------------------------------
! Print the displacement vector field
!
  subroutine plot_vect()

  double precision :: coorg(ndime,grid%ngnod),point(ndime) &
                     ,Uinterp(ndime),factor
  integer :: e,i,j,k,l,ipoin

  write(psunit,*) '%'
  write(psunit,*) '% vector field'
  write(psunit,*) '%'

 ! scale the vector field
  if (ScaleField > 0.d0) then
    factor = 1.d0/ScaleField
  else
    factor = 1.d0/maxfield
  endif

 ! fleches en couleur si modele de vitesse en background
  if(cpmodel) then
    write(psunit,*) 'Colvects'
  else
    write(psunit,*) '0 setgray'
  endif

  if (interpol) then  !-- interpolate the vector field

    do e=1,grid%nelem

      coorg = grid%coorg(:,grid%knods(:,e))

      do j=1,DisplayPts
      do i=1,DisplayPts

        point = matmul( coorg, interp%shape(:,i,j) )

        Uinterp = 0.d0
        do l=1,grid%ngll
        do k=1,grid%ngll
          Uinterp = Uinterp &
           + vfield(grid%ibool(k,l,e),:)*interp%flagrange(k,i)*interp%flagrange(l,j)
        enddo
        enddo

        call plot_single_vect(factor*Uinterp,point)

      enddo
      enddo
    enddo

  else !-- Plot vectors at mesh nodes

    do ipoin=1,grid%npoin
      call plot_single_vect(factor*vfield(ipoin,:), grid%coord(:,ipoin) )
    enddo

  endif

  write(psunit,*) '0 setgray'


  end subroutine plot_vect


!----------------------------------------------------------------------
! Plots scalar values defined element wise
  subroutine plot_efield

  double precision :: min_efield,max_efield,scale
  real :: x1
  integer :: e

  write(psunit,*) '% Plot an element wise scalar field'

  min_efield = minval(efield)
  max_efield = maxval(efield)
  if (max_efield-min_efield<abs(min_efield)*1d-10) then !WARNING: ugly test
    write(psunit,*) '% WARNING: constant field, skip plot'
    return    
  else
    scale=1d0/(max_efield-min_efield)
  endif

  do e=1,grid%nelem
    write(psunit,500) point_scaled( grid%coorg(:,grid%knods(1,e)) )
    write(psunit,499) point_scaled( grid%coorg(:,grid%knods(2,e)) )
    write(psunit,499) point_scaled( grid%coorg(:,grid%knods(3,e)) )
    write(psunit,499) point_scaled( grid%coorg(:,grid%knods(4,e)) )

    x1 = (efield(e)-min_efield)*scale
    !x1 = min( x1*0.7 + 0.2 , 1. ) ! rescale to avoid too dark a gray
    x1 = 1. - x1 ! invert the scale: white = min, gray = max

    write(psunit,604) x1

  enddo

 499  format(F0.2,1x,F0.2,' L')
 500  format(F0.2,1x,F0.2,' M')
 604  format('CP ',F0.2,' BK')
    
  end subroutine plot_efield

!----------------------------------------------------------------------
! Places a vector VECT at location ORIGIN
  subroutine plot_single_vect(VECT,ORIGIN)

    double precision, dimension(ndime), intent(in) :: VECT,ORIGIN

    double precision :: d,d1,d2,dummy,theta,thetaup,thetadown &
                       ,x2,z2,xa,za,xb,zb,convert
    
   ! pi conversion
    convert = 3.141592653589793d0/180.d0

    x2 = VECT(1)*sizemax
    z2 = VECT(2)*sizemax

    d = sqrt(x2*x2 + z2*z2)

! ignorer si vecteur trop petit
    if (d < cutvect*sizemax) return

    d1 = d * ArrowHBRatio
    d2 = d1 * cos(ArrowHeadAngle*convert)

    dummy = x2/d
    dummy = min(dummy,0.9999d0)
    dummy = max(dummy,-0.9999d0)
    theta = acos(dummy)

    if(z2 < 0.d0) theta = 360.d0*convert - theta
    thetaup = theta - ArrowHeadAngle*convert
    thetadown = theta + ArrowHeadAngle*convert

! tracer le vecteur proprement dit
    x2 = x2 * centim
    z2 = z2 * centim
    xa = -d2*cos(thetaup)
    za = -d2*sin(thetaup)
    xa = xa * centim
    za = za * centim
    xb = -d2*cos(thetadown)
    zb = -d2*sin(thetadown)
    xb = xb * centim
    zb = zb * centim
    write(psunit,700) xb,zb,xa,za,x2,z2,point_scaled(ORIGIN)

  700  format(8(F0.1,1x),'F') 
    
  end subroutine plot_single_vect

!------------------------------------------------------------------
!
  function point_scaled(coord)

    double precision, intent(in) :: coord(ndime)
    double precision :: point_scaled(ndime)

    point_scaled(1) = (coord(1)-xmin)*rapp_page + PageX0
    point_scaled(2) = (coord(2)-zmin)*rapp_page + PageZ0
    point_scaled    = point_scaled * centim

  end function point_scaled

!------------------------------------------------------------------
!
  subroutine plot_sources

  double precision :: coord(ndime)

  call SO_inquire(src,coord=coord)
  write(psunit,510) point_scaled(coord)
  if (symbols) then
    write(psunit,*) 'Cross'
  else
    write(psunit,*) '(S1) show'
  endif

  510  format(F0.1,1x,F0.1,' M')
  end subroutine plot_sources

!------------------------------------------------------------------
!
  subroutine plot_receivers
  
  double precision, pointer :: posrec(:,:)
  integer :: nrec,i
  !-- limite pour afficher des points a la place des recepteurs
  integer, parameter :: ndots = 10

  call REC_inquire(rec,coord=posrec)
  nrec = size(posrec,2)
  write(psunit,*) '% debut ligne recepteurs'
  do i=1,nrec
    write(psunit,510) point_scaled(posrec(:,i))
    if (symbols) then
      if(nrec > ndots.and.i /= 1.and.i /= nrec) then
        write(psunit,*) 'VDot'
      else
        write(psunit,*) 'Losange'
      endif
    else
      write(psunit,*) '(R',i,') show'
    endif
  enddo
  write(psunit,*) '% fin ligne recepteurs'

  510  format(F0.1,1x,F0.1,' M')
  end subroutine plot_receivers

end subroutine PLOT_PS

end module plot_postscript
