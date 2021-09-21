module plot_postscript

  use spec_grid 
  use constants, only : NDIME

  implicit none
  private

  double precision, parameter :: centim = 28.5d0
  integer, parameter :: maxcolors = 20
  real, dimension(3,maxcolors), save :: RGB
  logical, save :: legend=.true., interpol=.false., vectors=.false. &
                  ,mesh=.false. , symbols=.true., boundaries=.true. &
                  ,numbers=.false., color=.true.
  integer, save :: set_background=0

  double precision, allocatable, save :: se_interp(:,:,:), fe_interp(:,:,:)

  integer, save :: isubsamp=3, DisplayPts=3
  double precision, save :: ScaleField = 0.d0 
                   
  logical, parameter :: usletter=.true. ! Page format: US letter or A4

  double precision, save :: sizex,sizez,rapp_page,xmin,zmin,xmax,zmax
  character(3), save :: version='2.3'

  real, save :: usoffset

  public :: PLOT_PS, PLOT_PS_read

contains

!=======================================================================
!
!
! BEGIN INPUT BLOCK
!
! NAME   : SNAP_PS
! GROUP  : SNAPSHOT_OUTPUTS
! PURPOSE: Preferences for PostScript snapshots
! SYNTAX : &SNAP_PS vectors, mesh, background, color,
!               isubsamp, boundaries, symbols, numbers, legend,
!               ScaleField, Interpol, DisplayPts /
!
! ARG: vectors          [log] [F] Plots a vectorial field with arrows
! ARG: mesh             [log] [F] Plots the mesh on background
! ARG: background       [char] [''] Filled background, only for vector plots:
!                                   ''   none 
!                                   'P'  P-velocity model
!                                   'S'  S-velocity model
!                                   'T'  domains 
! ARG: isubsamp         [int] [3] Subsampling of the GLL nodes for the
!                                 output of velocity model. 
!                                 The default samples every 3 GLL points.
! ARG: boundaries       [log] [T] Colors every tagged boundary
! ARG: symbols          [log] [T] Plots symbols for sources and receivers
! ARG: numbers          [log] [F] Plots the element numbers
! ARG: legend           [log] [T] Writes legends
! ARG: color            [log] [T] Color output
! ARG: ScaleField       [dble] [0d0] Fixed amplitude scale (saturation),
!                       convenient for comparing snapshots and making movies. 
!                       The default scales each snapshot by its maximum amplitude
! ARG: Interpol         [log] [F] Interpolate field on a regular subgrid 
!                       inside each element
! ARG: DisplayPts       [log] [3] Size of interpolation subgrid inside each 
!                       element is DisplayPts*DisplayPts. The default plots at 
!                       vertices, mid-edges and element center.
!               
! END INPUT BLOCK

  subroutine PLOT_PS_read(iin)

  use echo, only : iout,echo_input

  integer, intent(in) :: iin
  character :: background
  character(10) :: bg_name

  NAMELIST / SNAP_PS / vectors,numbers &
                         ,background,isubsamp,color &
                         ,boundaries,symbols,legend,mesh &
                         ,interpol,DisplayPts,ScaleField

  background = ''
     
  rewind(iin)
  read(iin,SNAP_PS,END=100)
100 continue
 
  if (vectors) then
    select case(background)
      case('P'); set_background=1; bg_name='P model'
      case('S'); set_background=2; bg_name='S model'
      case('T'); set_background=3; bg_name='domains'
      case default;   set_background=0; bg_name='none'
    end select
  else
    set_background=0
    bg_name = 'none'
  endif

  if (echo_input) write(iout,200) mesh,numbers &
                         ,bg_name,isubsamp,color &
                         ,boundaries,symbols,legend &
                         ,vectors,ScaleField,interpol,DisplayPts
  return

  200 format(//1x,'P o s t S c r i p t   O u t p u t s',/1x,35('='),//5x, &
  'Plot mesh . . . . . . . . . . . . . . . . . . (mesh) = ',L1/ 5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',L1/ 5x, &
  'Background fill . . . . . . . . . . . . (background) = ',A/ 5x, &
  'Subsampling for velocity model display  . (isubsamp) = ',I0/5x, &
  'Color display . . . . . . . . . . . . . . . .(color) = ',L1/ 5x, &
  'Plot boundaries . . . . . . . . . . . . (boundaries) = ',L1/ 5x, &
  'Plot symbols  . . . . . . . . . . . . . . .(symbols) = ',L1/ 5x, &
  'Write legends . . . . . . . . . . . . . . . (legend) = ',L1/ 5x, &
  'Plot vector fields  . . . . . . . . . . . .(vectors) = ',L1/ 5x, &
  'Amplitude-Scaling . . . . . . . . . . . (ScaleField) = ',F0.2/5x, &
  'Interpolate vector field  . . . . . . . . (interpol) = ',L1/5x, &
  'Points per edge for interpolation . . . (DisplayPts) = ',I0)


  end subroutine PLOT_PS_read

!=======================================================================
! NOTE: indexed color scales can also be defined in Level 2 PostScript by
!       [/Indexed /DeviceRGB 255 <... ... ...> ] setcolorspace 

  subroutine PLOT_PS_init(grid)

  use spec_grid, only : sem_grid_type,SE_init_interpol
  use fem_grid, only : FE_getshape,Fe_GetNodesPerElement

  type(sem_grid_type), intent(in) :: grid

  !-- figure size as % of page
  double precision, parameter :: rpercentx = 70.0d0, rpercentz = 77.0d0

  double precision :: xi , eta
  integer :: i,j,ngnod

  !-- color palette

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

  ! paper format A4 or US letter
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

  ! page size / physical domain size
  rapp_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin))/100.d0

 ! weigths for interpolation (regular subgrids on each element)
  if (interpol) then
    ngnod = FE_GetNodesPerElement(grid%fem)
    allocate(se_interp(grid%ngll*grid%ngll,DisplayPts,DisplayPts))
    allocate(fe_interp(ngnod,DisplayPts,DisplayPts))
    do j=1,DisplayPts
      eta = 2.d0* dble(j-1)/dble(DisplayPts-1) - 1.d0
      do i=1,DisplayPts
        xi = 2.d0* dble(i-1)/dble(DisplayPts-1) - 1.d0
        call SE_init_interpol(xi,eta,se_interp(:,i,j),grid) 
        fe_interp(:,i,j) = FE_getshape(xi,eta,ngnod)
      enddo
    enddo
  endif

 
  end subroutine PLOT_PS_init


!=======================================================================
! PostScript plot manager
!
! Cases:
!   vfield is present and is vectorial: 
!       if requested ('vectors') plot as arrows and
!       fill the elements with efield, velocity model, domain tags or none
!       else plot the amplitude as color cells
!   vfield is present and is scalar: color cells
!   efield is present: if vfield is present, use as background,
!       else plot as color cells
!
  subroutine PLOT_PS(file,vfield,efield,grid,mat,stitle &
                    ,it_in,time_in,src,rec,comp)

  use stdio, only : IO_new_unit
  use echo, only : echo_run,iout,fmt1,fmtok
  use prop_mat, only : matpro_elem_type
  use sources, only : source_type,SO_inquire
  use receivers, only : rec_type,REC_inquire
  use fem_grid, only : FE_GetNodesPerElement, FE_GetElementCoord

  character(*)       , intent(in) :: file 
  type(sem_grid_type), intent(in) :: grid
  double precision   , intent(in), optional :: vfield(:,:)
  double precision   , intent(in), optional :: efield(grid%nelem)
  type(matpro_elem_type), intent(in) :: mat(:)
  character(*)       , intent(in) :: stitle
  integer            , intent(in), optional :: it_in
  double precision   , intent(in), optional :: time_in
  type(source_type)  , optional, pointer :: src(:)
  type(rec_type)     , optional, pointer :: rec
  integer            , optional :: comp

 ! height of the domain tags/numbers, in centimeters
  double precision, parameter :: height = 0.25d0 

  double precision :: maxfield,time
  integer :: psunit,it,opt_comp
  logical :: nodal_field,elem_field,vector_field
  logical, save :: initialized = .false.

  if (.not.initialized) then
    call PLOT_PS_init(grid)
    initialized = .true.
  endif

  if (present(time_in) .and. present(it_in)) then
    time = time_in
    it = it_in
  else
    time = 0d0
    it = 0
  endif

  if (present(comp)) then
    opt_comp = comp
  else
    opt_comp = 0
  endif

  vector_field = .false.
  nodal_field = present(vfield)
  if (nodal_field) then
    maxfield = maxval(abs(vfield))
    vector_field = (size(vfield,2)==2)
  endif

  elem_field = present(efield)
  if (elem_field) maxfield = maxval(abs(efield))

  psunit = IO_new_unit()
  open(unit=psunit,file=trim(file),status='unknown')
  if (echo_run) write(iout,fmt1,advance='no') "Dump PostScript "//trim(file)

  call plot_header()

  if (legend) call plot_legend()

  ! Scale factor in the X and Z direction
  write(psunit,*) '%'
  write(psunit,*) '1.0 1.0 scale'
  write(psunit,*) '%'

  if (nodal_field) then

    if (vector_field) then
     !-- plot a nodal vector field with arrows
      if (vectors) then
        if (elem_field) then
          call plot_efield()
        else 
          call plot_model(set_background)
        endif
        if (mesh) call plot_mesh(color,numbers, set_background==3 )
        call plot_boundaries()
        call plot_vect()
     !-- plot the amplitude of a nodal vector field
      else
        if (opt_comp<1 .or. opt_comp>2) then
          call plot_scal( sqrt(vfield(:,1)*vfield(:,1)+vfield(:,2)*vfield(:,2)) )
        else
          call plot_scal( vfield(:,opt_comp) )
        endif
        if (mesh) call plot_mesh(color,numbers, .false.)
        call plot_boundaries()
      endif

   !-- plot a nodal scalar field as colored GLL cells
    else
      call plot_scal(vfield(:,1))
      if (mesh) call plot_mesh(color,numbers, .false.)
      call plot_boundaries()
    endif

 !-- plot an element scalar field as colored elements
  else if (elem_field) then
    call plot_efield()
    if (mesh) call plot_mesh(color,numbers, .false. )
    call plot_boundaries()
  endif


  if (symbols) then 
    if(color) then ! Sources and receivers in color ?
      write(psunit,*) 'Colreceiv'
    else
      write(psunit,*) '0 setgray'
    endif
    if (present(src) .and. associated(src)) call plot_sources()
    if (present(rec) .and. associated(rec)) call plot_receivers()
  endif

  write(psunit,*) '%'
  write(psunit,*) 'grestore'
  write(psunit,*) 'showpage'

  close(psunit)
  if (echo_run) write(iout,fmtok)

  return
 
  contains 


!-----------------------------------------------------------------------
! Plot headers

  subroutine plot_header()

  write(psunit,10) stitle,version
  write(psunit,*) '/CM {28.5 mul} def'  ! convert from cm to points (1/72 inch)
  write(psunit,*) '/L {lineto} def'
  write(psunit,*) '/LR {rlineto} def'
  write(psunit,*) '/M {moveto} def'
  write(psunit,*) '/MR {rmoveto} def'
  write(psunit,*) '/MK {mark} def' ! start array construction
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
  write(psunit,*) '%---- symbols '
  write(psunit,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(psunit,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR CP fill} def'
  write(psunit,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR CP fill} def'
  write(psunit,*) '/Cross {GS 0.05 CM SLW GS 3 3 MR -6. -6. LR ST GR'
  write(psunit,*) 'GS 3 -3 MR -6. 6. LR ST GR 0.01 CM SLW} def'
  write(psunit,*) '/SmallLine {M 0.07 CM 0 rlineto} def'
  write(psunit,*) '/Losange {GS 0.05 CM SLW 0 4.2 MR -3 -4.2 LR 3 -4.2 LR'
  write(psunit,*) '3 4.2 LR CP ST GR 0.01 CM SLW} def'
  write(psunit,*) '%---- color settings'
  write(psunit,*) '% vector fields in magenta'
  write(psunit,*) '/Colvects {0.01 CM SLW 1. 0. 1. RG} def'
  write(psunit,*) '% element mesh in chartre'
  write(psunit,*) '/Colmesh {0.02 CM SLW 0.5 1. 0. RG} def'
  write(psunit,*) '% source and receivers in cyan'
  write(psunit,*) '/Colreceiv {0. 1. 1. RG} def'
  write(psunit,*) '%---- macros'
  write(psunit,*) '% arrow'
  write(psunit,*) '/F {M LR gsave LR ST grestore LR ST} def'
  write(psunit,*) '% element contour'
  write(psunit,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(psunit,*) '% filled quad cell'
  write(psunit,*) '/FQ {M L L L CP SF} def'
  write(psunit,*) '/FQC {M L L L CP RF} def'
  write(psunit,*) '%'
  write(psunit,*) '.01 CM SLW'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.35 CM SCSF'
  write(psunit,*) '%'
  write(psunit,*) '/vshift ',-height/2,' CM def'
  write(psunit,*) '/Rshow { currentpoint stroke M'
  write(psunit,*) 'dup stringwidth pop neg vshift MR show } def'
  write(psunit,*) '/Cshow { currentpoint stroke M'
  write(psunit,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(psunit,*) '/fN {/Helvetica-Bold findfont ',height,' CM SCSF} def'
  write(psunit,*) '%'
  write(psunit,*) 'gsave newpath 90 rotate'
  write(psunit,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(psunit,*) '%'

  return

10   format('%!PS-Adobe-2.0',/,'%%',/,'%%Title: ',A,/, &
        '%%Creator: SEM2DPACK Version ',A,/, &
        '%%Author: Jean-Paul Ampuero',/, &
        '%%BoundingBox: 0 0 612 792',/,'%%')

  end subroutine plot_header


!-----------------------------------------------------------------------
! Plot legends

  subroutine plot_legend()

  write(psunit,*) '0 setgray'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.5 CM SCSF'
  write(psunit,*) '24. CM 1.2 CM M'
  write(psunit,610) usoffset,it
  write(psunit,*) '%'
  write(psunit,*) '24. CM 1.95 CM M'
  write(psunit,600) usoffset,time
  write(psunit,*) '%'
  write(psunit,*) '24. CM 2.7 CM M'
  write(psunit,640) usoffset,maxfield
  write(psunit,*) '%'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.6 CM SCSF'
  if (color) write(psunit,*) '.4 .9 .9 RG'
  write(psunit,*) '11 CM 1.1 CM M'
  write(psunit,*) '(X) show'
  write(psunit,*) '%'
  write(psunit,*) '1.4 CM 9.5 CM M'
  write(psunit,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(psunit,*) '(Z) show'
  write(psunit,*) 'grestore'
  write(psunit,*) '%'
  write(psunit,*) '/Times-Roman findfont'
  write(psunit,*) '.7 CM SCSF'
  if (color) write(psunit,*) '.8 0 .8 RG'
  write(psunit,*) '25.35 CM 18.9 CM M'
  write(psunit,*) usoffset,' CM 2 div neg 0 MR'
  write(psunit,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(psunit,*) '(',stitle,') show'
  write(psunit,*) 'grestore'
  write(psunit,*) '26.45 CM 18.9 CM M'
  write(psunit,*) usoffset,' CM 2 div neg 0 MR'
  write(psunit,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(psunit,*) '(SEM2DPACK '//version//' - Spectral Element Method) show'
  write(psunit,*) 'grestore'

  return

 600  format(F0.3,' neg CM 0 MR (Time =',EN12.3,' s) show')
 610  format(F0.3,' neg CM 0 MR (Time step = ',I0,') show')
 640  format(F0.3,' neg CM 0 MR (Max =',EN12.3,') show')

  end subroutine plot_legend

!-----------------------------------------------------------------------
! Draw the velocity model in background
! INPUT PS:  1=P, 2=S, default=no_plot

  subroutine plot_model(PS)

  use prop_mat, only : MAT_getProp
  integer, intent(in) :: PS   

  double precision :: cmax,cmin,celem(grid%ngll,grid%ngll),c
  integer :: i,ip,j,jp,e

  if (PS<1 .or. PS>2) return

  cmin =  huge(cmin)
  cmax = -huge(cmax)
  do e=1,grid%nelem
    if (PS==1) then
      call MAT_getProp(celem,mat(e),'cp')
    else
      call MAT_getProp(celem,mat(e),'cs')
    endif
    cmin = min(cmin, minval(celem))
    cmax = max(cmax, maxval(celem))
  enddo

 ! quit if quasi-homogeneous velocity model, variation < 1%
  if ((cmax-cmin)/(cmin+cmax) < 0.02d0) then
    write(psunit,*) '%'
    write(psunit,*) '% no background : delta_v/v = ', 2d0*(cmax-cmin)/(cmax+cmin)
    write(psunit,*) '%'
  endif

  write(psunit,*) '%'
  write(psunit,*) '% background element fill'
  write(psunit,*) '%'

  do e=1,grid%nelem

    if (PS==1) then
      call MAT_getProp(celem,mat(e),'cp')
    else
      call MAT_getProp(celem,mat(e),'cs')
    endif

    do i=1,grid%ngll-1,isubsamp
      ip = min(grid%ngll,i+isubsamp)
      
      do j=1,grid%ngll-1,isubsamp
        jp = min(grid%ngll,j+isubsamp)
 
        write(psunit,500) point_scaled( grid%coord(:,grid%ibool(i,j,e)) )
        write(psunit,499) point_scaled( grid%coord(:,grid%ibool(ip,j,e)) )
        write(psunit,499) point_scaled( grid%coord(:,grid%ibool(ip,jp,e)) )
        write(psunit,499) point_scaled( grid%coord(:,grid%ibool(i,jp,e)) )
    
        c = (celem(i,j)-cmin)/(cmax-cmin)
        c = min( c*0.7 + 0.2 , 1.d0 ) ! rescale to avoid dark gray
        c = 1.d0 - c ! inverse scale: white = cmin, gray = cmax
    
        write(psunit,604) c
 
      enddo
    enddo

  enddo

 499  format(F0.2,1x,F0.2,' L')
 500  format(F0.2,1x,F0.2,' M')
 604  format('CP ',F0.2,' SF')

  end subroutine plot_model

!-----------------------------------------------------------------------
!-- Draw spectral element mesh

  subroutine plot_mesh(colors,elem_numbers,fill_domains)

  logical, intent(in) :: colors,elem_numbers,fill_domains
  double precision :: coord(NDIME,grid%ngll,grid%ngll) &
                     ,point(NDIME)
  double precision, pointer :: coorg(:,:)
  integer :: e,i,j,is,ir,imat,icol,ngnod

  write(psunit,*) '%'
  write(psunit,*) '% spectral element mesh'
  write(psunit,*) '%'

  ngnod = FE_GetNodesPerElement(grid%fem)

  do e=1,grid%nelem

    write(psunit,*) '% elem ',e
    write(psunit,*) 'MK'

    if (ngnod == 4) then ! linear element shape
   
     ! get the coordinates of the control nodes
      coorg => FE_GetElementCoord(grid%fem,e)
      do i=1,ngnod
        write(psunit,601) point_scaled( coorg(:,i) )
      enddo
      write(psunit,601) point_scaled( coorg(:,1) )
      deallocate(coorg)
    
    else ! curved element shape

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

    if (fill_domains) then
     ! Use a different color for each domain
      imat = grid%tag(e)
      icol = mod(imat - 1,maxcolors) + 1
      write(psunit,600) RGB(:,icol)
    endif

    if (colors) then
      write(psunit,*) 'GC'
    else
      write(psunit,*) 'GG'
    endif

    ! write the element number
    if (elem_numbers) then

      coorg => FE_GetElementCoord(grid%fem,e)
      point = 0.25d0 * SUM( coorg(:,1:4) , dim =2 )
      deallocate(coorg)
      point = point_scaled(point)
      if (colors) write(psunit,*) '1 setgray'
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

  double precision :: point1(NDIME),point2(NDIME)
  integer :: i,bce,ideb,iend

  write(psunit,*) '%'
  write(psunit,*) '% mesh boundaries'
  write(psunit,*) '%'

  write(psunit,*) '0.05 CM SLW'

  do i = 1,size(grid%bounds(:))
    write(psunit,*) '% boundary tag ',grid%bounds(i)%tag
  ! set color cyclically from the palette
    write(psunit,'( 3(F0.2,1x),"RG" )') RGB(:, mod(grid%bounds(i)%tag-1,maxcolors)+1 )
  ! draw a straight segment for each boundary element
    do bce = 1,grid%bounds(i)%nelem
      ideb   = grid%bounds(i)%node( grid%bounds(i)%ibool(1,bce) )
      point1 = point_scaled( grid%coord(:,ideb) )
      iend   = grid%bounds(i)%node( grid%bounds(i)%ibool(grid%ngll,bce) )
      point2 = point_scaled( grid%coord(:,iend) )
      write(psunit,'( 2(2(F0.2,1x),A) )') point1,'M ',point2,'L ST'
    enddo
  enddo

  write(psunit,*) '0 setgray'
  write(psunit,*) '0.01 CM SLW'

  end subroutine plot_boundaries


!---------------------------------------------------------------------
! Plot a vector field defined at each node
!
  subroutine plot_vect()

  double precision :: point(NDIME),Uinterp(NDIME),factor
  double precision, pointer :: coorg(:,:)
  double precision, allocatable :: vloc(:,:)
  integer :: e,i,j,k,ipoin

  if (maxfield>0d0) then
    write(psunit,*) '%'
    write(psunit,*) '% vector field'
    write(psunit,*) '%'
  else
    write(psunit,*) '%'
    write(psunit,*) '% empty vector field'
    write(psunit,*) '%'
    return
  endif

 ! scale the vector field
  if (ScaleField > 0.d0) then
    factor = 1.d0/ScaleField
  else
    factor = 1.d0/maxfield
  endif

  if(color) then
    write(psunit,*) 'Colvects'
  else
    write(psunit,*) '0 setgray'
  endif

  if (interpol) then  !-- interpolate the vector field

    allocate(vloc(grid%ngll*grid%ngll,size(vfield,2)))
    do e=1,grid%nelem

      coorg => FE_GetElementCoord(grid%fem,e)

      k=1
      do i=1,grid%ngll
      do j=1,grid%ngll
        vloc(k,:) = vfield(grid%ibool(i,j,e),:)
        k=k+1
      enddo
      enddo

      do j=1,DisplayPts
      do i=1,DisplayPts
        point = matmul( coorg, fe_interp(:,i,j) )
        Uinterp = matmul(se_interp(:,i,j),vloc)
        call plot_single_vect(factor*Uinterp,point)
      enddo
      enddo

      deallocate(coorg)

    enddo
    deallocate(vloc)

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

  double precision :: min_efield,max_efield,scale,x1
  double precision, pointer :: coorg(:,:)
  integer :: e

  min_efield = minval(efield)
  max_efield = maxval(efield)
  scale = max_efield - min_efield
  if (scale < abs(min_efield+max_efield)*1d-3) then
    write(psunit,*) '%'
    write(psunit,*) '% empty element-wise scalar field'
    write(psunit,*) '%'
    return    
  else
    write(psunit,*) '%'
    write(psunit,*) '% element-wise scalar field'
    write(psunit,*) '%'
    scale=1d0/scale
  endif

  do e=1,grid%nelem

    x1 = (efield(e)-min_efield)*scale
    !x1 = min( x1*0.7 + 0.2 , 1. ) ! rescale to avoid too dark a gray
    x1 = 1.d0 - x1 ! invert the scale: white = min, gray = max

    coorg => FE_GetElementCoord(grid%fem,e)
    write(psunit,'(9(F0.2,1x),"FQ")') x1 &
                     ,point_scaled( coorg(:,1) ) &
                     ,point_scaled( coorg(:,2) ) &
                     ,point_scaled( coorg(:,3) ) &
                     ,point_scaled( coorg(:,4) )
    deallocate(coorg)

  enddo
    
  end subroutine plot_efield

!----------------------------------------------------------------------
! Plots nodal scalar values 
! Actually, the average of the 4 neighbouring GLL nodes
  subroutine plot_scal(s)

  double precision, intent(in) :: s(:)

  double precision, allocatable :: spl(:,:), cpl(:,:,:)
  double precision :: stmp(grid%ngll*grid%ngll)
  double precision, pointer :: coorg(:,:)
  double precision :: val, scale
  integer :: e,i,j,k,iglob,npl

  if (maxfield>0d0) then
    write(psunit,*) '%'
    write(psunit,*) '% scalar field'
    write(psunit,*) '%'
  else
    write(psunit,*) '%'
    write(psunit,*) '% empty scalar field'
    write(psunit,*) '%'
    return
  endif

  if (ScaleField>0d0) then
    scale = 1d0/ScaleField
  else
    scale = 1d0/maxval(abs(s))
  endif

  if (interpol) then
    npl = DisplayPts
  else
    npl = grid%ngll
  endif
  allocate(spl(npl,npl), cpl(npl,npl,2) )

  do e=1,grid%nelem

    if (interpol) then
      coorg => FE_GetElementCoord(grid%fem,e)
      k=0
      do j=1,grid%ngll
      do i=1,grid%ngll
        k=k+1
        stmp(k) = s(grid%ibool(i,j,e))
      enddo
      enddo
      do j=1,npl
      do i=1,npl
        cpl(i,j,:) = matmul( coorg, fe_interp(:,i,j) )
        spl(i,j) = sum( se_interp(:,i,j)*stmp)
      enddo
      enddo
      deallocate(coorg)

    else
      do j=1,grid%ngll
      do i=1,grid%ngll
        iglob = grid%ibool(i,j,e)
        cpl(i,j,:) = grid%coord(:,iglob)
        spl(i,j) = s(iglob)
      enddo
      enddo
    endif
    
    do j=1,npl-1
    do i=1,npl-1
      val = 0.25d0 *(spl(i,j)+spl(i+1,j)+spl(i+1,j+1)+spl(i,j+1)) *scale
      val = min(val,1d0)
      val = max(val,-1d0)
     ! at this point, val in [-1:1]

      if (color) then ! min = -1 = blue, max = 1 = red
        if (abs(val)<0.01d0) cycle   ! skip if white
        write(psunit,'(11(F0.2,1x),"FQC")') &
          1d0+min(val,0d0),1d0-abs(val),1d0-max(val,0d0) &
          ,point_scaled( cpl(i,j,:) ) &
          ,point_scaled( cpl(i+1,j,:) ) &
          ,point_scaled( cpl(i+1,j+1,:) ) &
          ,point_scaled( cpl(i,j+1,:) )

      else            ! min = 1 = white, max = 0 = gray
        if (val< -0.99d0) cycle   ! skip if white
        val = 0.5d0*(1d0-val) 
        write(psunit,'(9(F0.2,1x),"FQ")') val &
          ,point_scaled( cpl(i,j,:) ) &
          ,point_scaled( cpl(i+1,j,:) ) &
          ,point_scaled( cpl(i+1,j+1,:) ) &
          ,point_scaled( cpl(i,j+1,:) )

      endif

    enddo
    enddo
  enddo

  deallocate(spl, cpl )

  end subroutine plot_scal

!----------------------------------------------------------------------
! Places a vector VECT at location ORIGIN
  subroutine plot_single_vect(VECT,ORIGIN)

    double precision, dimension(NDIME), intent(in) :: VECT,ORIGIN

    double precision, parameter :: cutvect = 0.01d0 ! min arrow length / scale factor
    double precision, parameter :: sizemax = 1d0 ! Maximum arrow length (in cm)
    double precision, parameter :: ArrowHeadAngle = 20.d0 ! Angle between arrow body and head (degrees)
    double precision, parameter :: ArrowHBRatio = 0.4d0 ! Arrow head to body ratio

    double precision :: d,d1,d2,dummy,theta,thetaup,thetadown &
                       ,x2,z2,xa,za,xb,zb,convert
    
   ! pi conversion
    convert = 3.141592653589793d0/180.d0

    x2 = VECT(1)*sizemax
    z2 = VECT(2)*sizemax

    d = sqrt(x2*x2 + z2*z2)

! ignore small vectors
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

! plot vector
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

    double precision, intent(in) :: coord(NDIME)
    double precision :: point_scaled(NDIME)
   
   ! origin in the page (cm)
    double precision, parameter :: PageX0 = 2.2d0, PageZ0 = 2.2d0

    point_scaled(1) = (coord(1)-xmin)*rapp_page + PageX0
    point_scaled(2) = (coord(2)-zmin)*rapp_page + PageZ0
    point_scaled    = point_scaled * centim

  end function point_scaled

!------------------------------------------------------------------
!
  subroutine plot_sources

  double precision :: coord(NDIME)
  integer :: i

  write(psunit,*) '% sources'
  do i=1,size(src)
    call SO_inquire(src(i),coord=coord)
    write(psunit,510) point_scaled(coord)
    write(psunit,*) 'Cross'
  enddo

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
  write(psunit,*) '% start receiver line'
  do i=1,nrec
    write(psunit,510) point_scaled(posrec(:,i))
    if(nrec > ndots.and.i /= 1.and.i /= nrec) then
      write(psunit,*) 'VDot'
    else
      write(psunit,*) 'Losange'
    endif
  enddo
  write(psunit,*) '% end receiver line'

  510  format(F0.1,1x,F0.1,' M')
  end subroutine plot_receivers

end subroutine PLOT_PS

end module plot_postscript
