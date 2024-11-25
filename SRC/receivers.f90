module receivers

  use stdio, only: IO_abort
  use constants, only : NDIME

  implicit none
  private

  type rec_type
    private
    double precision, pointer :: coord(:,:),field(:,:)
    real, pointer :: sis(:,:,:)
    integer, pointer :: iglob(:)
    double precision :: tsamp
    integer :: nx,nt,isamp,irepr
    logical :: AtNode
    double precision, pointer :: interp(:,:)
    integer, pointer :: einterp(:)
    character :: SeisField

    ! Elif (2020)
    ! extra stations for stress-strain data
    integer :: Enx
    logical :: AtNodeE
    double precision, pointer :: Ecoord(:,:)
    double precision, pointer :: fieldsig(:,:)
    real, pointer :: sigma(:,:,:)
    integer, pointer :: Eiglob(:)
    double precision, pointer :: ELinterp(:,:)
    integer, pointer :: ELeinterp(:)    

  end type rec_type

  public :: rec_type,REC_read,REC_init,REC_store,REC_write,REC_inquire

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : REC_LINE
! PURPOSE: Defines a line of receivers
! SYNTAX : If single receiver line: 
!            &REC_LINE number,first,last,AtNode,isamp,field,irepr /
!          If receiver locations from file:
!            &REC_LINE file,AtNode,isamp,field,irepr /
!
! ARG: number   [int] [0] Number of stations in the line
! ARG: first    [dble(2)] Receivers can be located along a line,
!                this is the position (x,z) of the first receiver
! ARG: last     [dble(2)] Position (x,z) of the last receiver,
!                other receivers will be located with regular spacing
!                between First and Last.
! ARG: file     [name] ['none'] Station positions can instead be read 
!                from an ASCII file, with 2 columns: X and Z (in meters)
! ARG: AtNode   [log] [T] Relocate the stations at the nearest GLL node
! ARG: isamp    [int] [1] Sampling stride (in number of timesteps). Note that
!                for stability reasons the timestep can be very small.
! ARG: field    [char] ['V'] The field in the seismogram:
!                               'D'     displacement
!                               'V'     velocity
!                               'A'     acceleration
! ARG: irepr    [char] ['D'] Abscissa for the seismic multitrace plot:
!                               'X' Horizontal position
!                               'Z' Depth
!                               'D' Distance to the first station
!
! ARG: extra    [log] [F] Write out stress-strain curves (not for elasticity)
!
!
! NOTE   : to locate receivers at the free surface set their vertical position 
!          above the free surface and AtNode=T
!
! END INPUT BLOCK

  subroutine REC_read(rec,iin)

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit, IO_file_length

  type(rec_type), pointer :: rec
  integer, intent(in) :: iin

  double precision :: first(NDIME),last(NDIME),init_double
  integer :: i,iin2,isamp,number
  logical :: AtNode, extra
  character(50) :: file
  character :: field,irepr
  
  NAMELIST / REC_LINE / number,isamp,field,first,last,AtNode,file,irepr,extra

  init_double = huge(init_double) ! set to an unlikely value

  number     = 0
  isamp      = 1
  field      = 'V'
  first      = init_double
  last       = init_double
  AtNode     = .true.
  file       = 'none'
  irepr      = 'D'
  extra      = .false.

  rewind(iin)
  read(iin,REC_LINE,END = 200) 
  allocate(rec)

  if (number < 0) call IO_abort('REC_read: "number" must be positive')

  rec%isamp  = isamp
  rec%SeisField = field
  rec%AtNode = AtNode

  if ( field /='D' .and. field /='V ' .and. field /='A') &
    call IO_abort('REC_read: parameter field has wrong value [D,V,A]')

  select case(irepr)
    case('X'); rec%irepr  = 1
    case('Z'); rec%irepr  = 2
    case('D'); rec%irepr  = 3
    case default; call IO_abort('REC_read: unknown "irepr" [X,Z,D]')
  end select

  if (file == 'none') then

    if ( any(first == init_double) ) &
      call IO_abort('REC_read: you must set "first" station coordinates')
    if ( any(last == init_double) ) &
      call IO_abort('REC_read: you must set "last" station coordinates')
    rec%nx   = number
    if (echo_input) write(iout,100) rec%nx,first,last
    allocate(rec%coord(NDIME,rec%nx))
    if (number>1) then
      do i = 1,rec%nx
        rec%coord(:,i) = first + (i-1)/dble(rec%nx-1) *(last-first)
      enddo
    else
      rec%coord(:,1) = first    
    endif

  else
    rec%nx = IO_file_length(file)
    if (echo_input) write(iout,110) rec%nx,file
    allocate(rec%coord(NDIME,rec%nx))
    iin2 = IO_new_unit()
    open(iin2,file=file,status='old')
    do i=1,rec%nx
      read(iin2 ,*) rec%coord(:,i) 
    enddo
    close(iin2)
  endif

  if (echo_input) write(iout,120) AtNode,isamp,field,irepr


  ! Extra receivers for stress-strain
  rec%Enx = 0
  if (extra) &
  call REC_extra_read(rec,iin)


  return

  100 format(//1x,'R e c e i v e r s', &
  /1x,17('='),//5x, &
  'Number of receivers . . . . . . . . . . . . (number) = ',I0/5x, &
  'First receiver X location . . . . . . . . (first(1)) = ',EN12.3/5x, &
  'First receiver Z location . . . . . . . . (first(2)) = ',EN12.3/5x, &
  'Last receiver X location. . . . . . . . . .(last(1)) = ',EN12.3/5x, &
  'Last receiver Z location. . . . . . . . . .(last(2)) = ',EN12.3)

  110 format(//1x,'R e c e i v e r s', &
  /1x,17('='),//5x, &
  'Number of receivers . . . . . . . . . . . . (number) = ',I0/5x, &
  'Receiver locations file name. . . . . . . . . (file) = ',A)

  120 format(/5x, &
  'Relocate to the nearest GLL node. . . . . . (AtNode) = ',L1/5x, &
  'Subsampling for seismograms recording . . . .(isamp) = ',I0/5x, &
  'Field recorded. . . . . . . . . . . . . . . .(field) = ',A/5x, &
  'Axis of the seismogram plot . . . . . . . . .(irepr) = ',A)

  200 return

  end subroutine REC_read
!=====================================================================
!

subroutine REC_extra_read(rec,iin)
!
! BEGIN INPUT BLOCK
!
! NAME   : REC_LINEX
! PURPOSE: Defines a line of receivers
! SYNTAX : If single receiver line: 
!            &REC_LINEX numberE,firstE,lastE,AtNodeE /
!          If receiver locations from file:
!            &REC_LINEX fileE,AtNodeE /
!
! ARG: numberE   [int] [0] Number of stations in the line
! ARG: firstE    [dble(2)] Receivers can be located along a line,
!                this is the position (x,z) of the first receiver
! ARG: lastE     [dble(2)] Position (x,z) of the last receiver,
!                other receivers will be located with regular spacing
!                between First and Last.
! ARG: fileE     [name] ['none'] Station positions can instead be read 
!                from an ASCII file, with 2 columns: X and Z (in meters)
! ARG: AtNodeE   [log] [T] Relocate the stations at the nearest GLL node
!
!
! END INPUT BLOCK

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit, IO_file_length

  type(rec_type), pointer :: rec
  integer, intent(in) :: iin

  integer :: numberE,i,iin2
  double precision :: firstE(NDIME),lastE(NDIME),init_double
  character(50) :: fileE
  logical :: AtNodeE

  NAMELIST / REC_LINEX / numberE,firstE,lastE,AtNodeE,fileE

  init_double = huge(init_double) ! set to an unlikely value
  numberE     = 0
  firstE      = init_double
  lastE       = init_double
  AtNodeE     = .true.
  fileE       = 'none'

  read(iin,REC_LINEX,END = 200) 
 
  if (numberE < 0) call IO_abort('REC_extra_read: "numberE" must be positive')
  rec%AtNodeE = AtNodeE

   if (fileE == 'none') then
    if ( any(firstE == init_double) ) &
      call IO_abort('REC_extra_read: you must set "firstE" station coordinates')
    if ( any(lastE == init_double) ) &
      call IO_abort('REC_extra_read: you must set "lastE" station coordinates')
    rec%Enx   = numberE
    if (echo_input) write(iout,130) rec%Enx,firstE,lastE
    allocate(rec%Ecoord(NDIME,rec%Enx))
    if (numberE>1) then
      do i = 1,rec%Enx
        rec%Ecoord(:,i) = firstE + (i-1)/dble(rec%Enx-1) *(lastE-firstE)
      enddo
    else
      rec%Ecoord(:,1) = firstE    
    endif
  else
    rec%Enx = IO_file_length(fileE)
    if (echo_input) write(iout,140) rec%Enx,fileE
    allocate(rec%Ecoord(NDIME,rec%Enx))
    iin2 = IO_new_unit()
    open(iin2,file=fileE,status='old')
    do i=1,rec%Enx
      read(iin2 ,*) rec%Ecoord(:,i) 
    enddo
    close(iin2)
  endif

  if (echo_input) write(iout,150) AtNodeE


  130 format(//1x,'E x t r a  R e c e i v e r s', &
  /1x,28('='),//5x, &
  'Number of extra receivers . . . . . . . . .(numberE) = ',I0/5x, &
  'First receiver X location . . . . . . . . (first(1)) = ',EN12.3/5x, &
  'First receiver Z location . . . . . . . . (first(2)) = ',EN12.3/5x, &
  'Last receiver X location. . . . . . . . . .(last(1)) = ',EN12.3/5x, &
  'Last receiver Z location. . . . . . . . . .(last(2)) = ',EN12.3)

  140 format(//1x,'E x t r a  R e c e i v e r s', &
  /1x,28('='),//5x, &
  'Number of extra receivers . . . . . . . . .(numberE) = ',I0/5x, &
  'Extra receiver locations file name . . . . . (fileE) = ',A)

  150 format(/5x, &
  'Relocate to the nearest GLL node. . . . . . (AtNode) = ',L1/5x)

  200 return

end subroutine REC_extra_read
!=====================================================================
!=====================================================================
!
!!  Initialize
  subroutine REC_init(rec,grid,time,fields)

    use fields_class, only : fields_type
    use time_evol, only : timescheme_type, TIME_getTimeStep, TIME_getNbTimeSteps
    use spec_grid, only : sem_grid_type
    use memory_info
    use echo, only: iout,info=>echo_init
    use stdio, only: IO_new_unit

    type(rec_type), intent(inout) :: rec
    type(sem_grid_type), intent(in) :: grid
    type(fields_type), intent(in) :: fields
    type(timescheme_type), intent(in) :: time

    integer :: ounit,irec

    call REC_posit(rec,grid)

    rec%tsamp = TIME_getTimeStep(time)*rec%isamp

    rec%nt = TIME_getNbTimeSteps(time)/rec%isamp +1
    if (rec%nt == 0) call IO_abort('REC_init: zero samples')
      
    allocate(rec%sis(rec%nt,rec%nx,fields%ndof))
    rec%sis = 0.
    call storearray('rec%sis',size(rec%sis),idouble)

    ! extra
    if (rec%Enx > 0) then
      allocate(rec%sigma(rec%nt,rec%Enx,3))
      rec%sigma = 0d0
      call storearray('rec%sigma',size(rec%sigma),idouble)        
    endif

    ! Extra receivers for stress-strain curves
    if (rec%Enx > 0) &
    rec%fieldsig => fields%sigeps


    select case(rec%SeisField)
      case('D'); rec%field => fields%displ
      case('V'); rec%field => fields%veloc
      case('A'); rec%field => fields%accel
      case default
        call IO_abort('REC_init: Wrong field code for seismograms display')
    end select

  if (info) then
    write(iout,*)
    write(iout,100) 'Sampling rate (Hz)        = ',1d0/rec%tsamp
    write(iout,100) 'Sampling timestep (secs)  = ',rec%tsamp
    write(iout,200) 'Total number of samples   = ',rec%nt
    write(iout,200) 'Number of receivers       = ',rec%nx
    write(iout,*)
  endif

! Save headers
  ounit = IO_new_unit()
  open(ounit,file='SeisHeader_sem2d.hdr',status='replace')
  write(ounit,*) 'DT NSAMP NSTA'
  write(ounit,*) real(rec%tsamp),rec%nt,rec%nx
  write(ounit,*) 'XSTA ZSTA'
  do irec=1,rec%nx
    write(ounit,*) rec%coord(:,irec)
  enddo

  ! Additional headers
  if (rec%Enx > 0) then
    write(ounit,*) rec%Enx
    write(ounit,*) 'XSTA-extra ZSTA-extra'
    do irec=1,rec%Enx
      write(ounit,*) rec%Ecoord(:,irec)
    enddo
  endif

  close(ounit)

100 format(2X,A,EN12.3)
200 format(2X,A,I0)

  end subroutine REC_init

!=====================================================================
!
!!  Real receivers position
  subroutine REC_posit(rec,grid)

  use echo, only: iout,info=>echo_init
  use spec_grid, only : sem_grid_type,SE_find_nearest_node,SE_init_interpol,SE_find_point
  use memory_info

  type(rec_type), intent(inout) :: rec
  type(sem_grid_type), intent(in) :: grid

  integer :: iglob_tmp(rec%nx)
  double precision :: distmax,dist, xi,eta, new_coord(NDIME)
  integer :: n,ipoint,irec


  distmax = 0.d0

  if (rec%AtNode) then

    if (info) write(iout,200)
    iglob_tmp = 0
    irec = 0
    do n=1,rec%nx
      call SE_find_nearest_node(rec%coord(:,n),grid,ipoint,distance = dist)
     ! do not keep it if already in list
      if ( irec > 1 .and. any(iglob_tmp(:irec) == ipoint)) cycle
      irec = irec+1
      iglob_tmp(irec) = ipoint
      if (info) then
        write(iout,150) irec,rec%coord(:,n),grid%coord(:,ipoint),dist
        distmax = max(dist,distmax)
      endif
    enddo
    if (irec < rec%nx) then 
      rec%nx = irec
      deallocate(rec%coord)
      allocate(rec%coord(NDIME,rec%nx))
    endif
    call storearray('rec%coord',size(rec%coord),idouble)
    allocate(rec%iglob(rec%nx))
    call storearray('rec%iglob',size(rec%iglob),isngl)
    rec%iglob = iglob_tmp(:rec%nx)
    rec%coord(:,:) = grid%coord(:,rec%iglob)

  else

    if (info) write(iout,300)
    allocate(rec%interp(grid%ngll*grid%ngll,rec%nx))
    allocate(rec%einterp(rec%nx))
    do n=1,rec%nx
      call SE_find_point(rec%coord(:,n),grid,rec%einterp(n),xi,eta,new_coord)
      call SE_init_interpol(xi,eta,rec%interp(:,n),grid)
      dist = sqrt( (rec%coord(1,n)-new_coord(1))**2 + (rec%coord(2,n)-new_coord(2))**2 )
      if (info) then
        write(iout,150) n,rec%coord(:,n),new_coord,dist
        distmax = max(dist,distmax)
      endif
      rec%coord(:,n)=new_coord
    enddo

  endif

  if (info) write(iout,160) distmax

  ! For extra receivers
  if (rec%Enx > 0) &
  call REC_extra_init(rec,grid)


 150   format(2x,i7,5(1x,EN12.3))
 160   format(/2x,'Maximum distance between asked and real =',EN12.3)
 200  format(//1x,'R e c e i v e r s'/1x,17('=')// &
  ' Receivers have been relocated to the nearest GLL node'// &
  ' Receiver  x-requested  z-requested   x-obtained   z-obtained   distance'/)

 300  format(//1x,'R e c e i v e r s'/1x,17('=')// &
  ' Receiver  x-requested  z-requested   x-obtained   z-obtained   distance'/)

  end subroutine REC_posit
!=====================================================================
!

subroutine REC_extra_init(rec,grid)

  use echo, only: iout,info=>echo_init
  use spec_grid, only : sem_grid_type,SE_find_nearest_node,SE_init_interpol,SE_find_point
  use memory_info

  type(rec_type), intent(inout) :: rec
  type(sem_grid_type), intent(in) :: grid

  integer :: iglob_tmp(rec%Enx)
  double precision :: distmax,dist, xi,eta, new_coord(NDIME)
  integer :: n,ipoint,irec

  distmax = 0.d0

  ! Save at closest nodes
  if (rec%AtNodeE) then
    if (info) write(iout,200)
    iglob_tmp = 0
    irec = 0
    do n=1,rec%Enx
      call SE_find_nearest_node(rec%Ecoord(:,n),grid,ipoint,distance = dist)
     ! do not keep it if already in list
      if ( irec > 1 .and. any(iglob_tmp(:irec) == ipoint)) cycle
      irec = irec+1
      iglob_tmp(irec) = ipoint
      if (info) then
        write(iout,150) irec,rec%Ecoord(:,n),grid%coord(:,ipoint),dist
        distmax = max(dist,distmax)
      endif
    enddo
    if (irec < rec%Enx) then 
      rec%Enx = irec
      deallocate(rec%Ecoord)
      allocate(rec%Ecoord(NDIME,rec%Enx))
    endif
    call storearray('rec%Ecoord',size(rec%Ecoord),idouble)
    allocate(rec%Eiglob(rec%Enx))
    call storearray('rec%Eiglob',size(rec%Eiglob),isngl)
    rec%Eiglob = iglob_tmp(:rec%Enx)
    rec%Ecoord(:,:) = grid%coord(:,rec%Eiglob)
  else
  ! Interpolate at specified coordinates
    if (info) write(iout,300)
    allocate(rec%ELinterp(grid%ngll*grid%ngll,rec%Enx))
    allocate(rec%ELeinterp(rec%Enx))
    do n=1,rec%Enx
      call SE_find_point(rec%Ecoord(:,n),grid,rec%ELeinterp(n),xi,eta,new_coord)
      call SE_init_interpol(xi,eta,rec%ELinterp(:,n),grid)
      dist = sqrt( (rec%Ecoord(1,n)-new_coord(1))**2 + (rec%Ecoord(2,n)-new_coord(2))**2 )
      if (info) then
        write(iout,150) n,rec%Ecoord(:,n),new_coord,dist
        distmax = max(dist,distmax)
      endif
      rec%Ecoord(:,n)=new_coord
    enddo
  endif

  if (info) write(iout,160) distmax

 150   format(2x,i7,5(1x,EN12.3))
 160   format(/2x,'Maximum distance between asked and real =',EN12.3)

 200  format(//1x,'E x t r a   R e c e i v e r s'/1x,17('=')// &
  ' Receivers have been relocated to the nearest GLL node'// &
  ' Receiver  x-requested  z-requested   x-obtained   z-obtained   distance'/)

 300  format(//1x,'R e c e i v e r s'/1x,17('=')// &
  ' Receiver  x-requested  z-requested   x-obtained   z-obtained   distance'/)

end subroutine REC_extra_init
!=====================================================================

!=====================================================================
!
!!  Store the seismograms
  subroutine REC_store(rec,it,grid)

  use spec_grid, only : sem_grid_type

  type(rec_type), intent(inout) :: rec
  integer, intent(in) :: it
  type(sem_grid_type), intent(in) :: grid

  integer :: itsis,n,i,k,j,e,iglob
  double precision, allocatable :: vloc(:,:), sigloc(:,:)

  if ( mod(it,rec%isamp) /= 0 ) return

  itsis = it/rec%isamp +1
  if (itsis > rec%nt) call IO_abort('receivers.REC_store: storage is full')

  if (rec%AtNode) then
    rec%sis(itsis,:,:) = rec%field(rec%iglob,:)
  else
    allocate(vloc(grid%ngll*grid%ngll,size(rec%field,2)))
    do n=1,rec%nx
      e = rec%einterp(n)
      k=1
      do j=1,grid%ngll
      do i=1,grid%ngll
        iglob = grid%ibool(i,j,e)
        vloc(k,:) = rec%field(iglob,:)
        k=k+1
      enddo
      enddo
      rec%sis(itsis,n,:) = matmul(rec%interp(:,n),vloc)
    enddo
    deallocate(vloc)
  endif


! Extra receivers
  if (rec%Enx > 0) then
    if (rec%AtNodeE) then
      rec%sigma(itsis,:,:) = rec%fieldsig(rec%Eiglob,:)
    else
      allocate(sigloc(grid%ngll*grid%ngll,size(rec%fieldsig,2)))
      do n=1,rec%Enx
        e = rec%ELeinterp(n)
        k=1
        do j=1,grid%ngll
        do i=1,grid%ngll
          iglob = grid%ibool(i,j,e)
          sigloc(k,:) = rec%fieldsig(iglob,:)
          k=k+1
        enddo
        enddo
        rec%sigma(itsis,n,:) = matmul(rec%ELinterp(:,n),sigloc)
      enddo
      deallocate(sigloc)
    endif
  endif


    
  end subroutine REC_store



!=====================================================================
!
!!  Export the seismograms
  subroutine REC_write(rec,iout)

  use stdio, only: IO_new_unit
  double precision, parameter :: factorxsu = 3.5d0

  type(rec_type), intent(inout) :: rec
  integer, intent(in) :: iout

  double precision :: xval(rec%nx),xref,zref
  integer :: iol,ounit,k
  character(30) :: ylabel
  character, parameter :: posvars(3) = (/ 'X','Z','D' /)
  character :: posvar

!---- binary data -------------------------------------------------------
  write(iout,*) 'Storing seismograms (SEP format) ...'
  INQUIRE( IOLENGTH=iol ) rec%sis(:,1,1)
  ounit = IO_new_unit()

  if ( size(rec%sis,3)==1) then

    open(ounit,file='Uy_sem2d.dat',status='replace',access='direct',recl=iol)
    do k=1,size(rec%sis,2)
      write(ounit,rec=k) rec%sis(:,k,1)
    enddo
    close(ounit)

  else

    open(ounit,file='Ux_sem2d.dat',status='replace',access='direct',recl=iol)
    do k=1,size(rec%sis,2)
      write(ounit,rec=k) rec%sis(:,k,1)
    enddo
    close(ounit)
  
    open(ounit,file='Uz_sem2d.dat',status='replace',access='direct',recl=iol)
    do k=1,size(rec%sis,2)
      write(ounit,rec=k) rec%sis(:,k,2)
    enddo
    close(ounit)

  endif


! ----Exporting stress-strain (only) in binary data ----------------------------
  if (rec%Enx > 0) then
    
    open(ounit,file='Strain_sem2d.dat',status='replace',access='direct',recl=iol)
    ! Write out strain
    do k=1,size(rec%sigma,2)
      write(ounit,rec=k) rec%sigma(:,k,1)
    enddo
    close(ounit)

    open(ounit,file='Stress_sem2d.dat',status='replace',access='direct',recl=iol)
    ! Write out stress
    do k=1,size(rec%sigma,2)
      write(ounit,rec=k) rec%sigma(:,k,2)
    enddo
    close(ounit)

    open(ounit,file='Nonlinear_surfaces_sem2d.dat',status='replace',access='direct',recl=iol)
    ! Write out stress
    do k=1,size(rec%sigma,2)
      write(ounit,rec=k) rec%sigma(:,k,3)
    enddo
    close(ounit)

  endif


!=== Seismic Unix scripts: =================================================

  posvar = posvars(rec%irepr)

  select case(rec%irepr)
    case(1,2) ! recepteurs suivant coordonnee X ou Z
      xval = rec%coord(rec%irepr,:)
    case(3) ! recepteurs en distance
      xref = rec%coord(1,1)
      zref = rec%coord(2,1)
      xval = sqrt((rec%coord(1,:) - xref)**2 + &
              (rec%coord(2,:) - zref)**2)
  end select

  select case(rec%SeisField)
    case('D'); ylabel='Displacement (m)'
    case('V'); ylabel='Velocity (m/s)'
    case('A'); ylabel='Acceleration (m/s^2)'
  end select

! station "coordinates" for plots
  open(ounit,file='x2_sem2d.tab',status='unknown')
  write(ounit,'(1000(f0.2:","))') xval
               !this number > max nb stations ever
  close(ounit)

!-- Xwindow -------------------------------------------------------------------
  open(ounit,file='Xline_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set x2 = `cat x2_sem2d.tab`'
  write(ounit,100) 'xwigb',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Ux component" < Ux_sem2d.dat'
  write(ounit,100) 'xwigb',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Uz component" < Uz_sem2d.dat'
  close(ounit)

!-- PostScript ---------------------------------------------------------------
  open(ounit,file='PSline_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set x2 = `cat x2_sem2d.tab`'
  write(ounit,100) 'pswigp',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Ux component" < Ux_sem2d.dat >! UxPoly_sem2d.ps'
  write(ounit,100) 'pswigp',factorxsu,rec%nt,rec%tsamp,posvar,rec%nx,&
    ' title="Uz component" < Uz_sem2d.dat >! UzPoly_sem2d.ps'
  close(ounit)

!-- one trace for Xwindow --------------------------------------------

  open(ounit,file='Xtrace_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set nt = ',rec%nt
  write(ounit,*) 'set dt = ',real(rec%tsamp)
  write(ounit,*) '@ trace=10000'
  write(ounit,*) 'while ($trace > -1)'
  write(ounit,*) 'echo Trace number ?'
  write(ounit,*) 'set rep=$<'
  write(ounit,*) '@ trace = $rep'
  write(ounit,*) 'echo Trace asked : $trace'
  write(ounit,*) '# traces commencent a zero dans format SEP'
  write(ounit,*) '@ septrace = $trace - 1'
  write(ounit,*) 'subset < Ux_sem2d.dat n1=$nt', &
   ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
   ' label1="Time (s)" label2="',trim(ylabel),'" n=$nt style=normal d1=$dt', &
   ' title="Ux component trace "$trace &'
  write(ounit,*) 'subset < Uz_sem2d.dat n1=$nt', &
   ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
   ' label1="Time (s)" label2="',trim(ylabel),'" n=$nt style=normal d1=$dt', &
   ' title="Uz component trace "$trace &'
  write(ounit,*) 'end'
  close(ounit)

!-- one trace for PostScript --------------------------------------------

  open(ounit,file='PStrace_sem2d.csh',status='unknown')
  write(ounit,110)
  write(ounit,*) 'set nt = ',rec%nt
  write(ounit,*) 'set dt = ',real(rec%tsamp)
  write(ounit,*) '@ trace=10000'
  write(ounit,*) 'while ($trace > -1)'
  write(ounit,*) 'echo Trace number ?'
  write(ounit,*) 'set rep=$<'
  write(ounit,*) '@ trace = $rep'
  write(ounit,*) 'echo Trace asked : $trace'
  write(ounit,*) '# traces commencent a zero dans format SEP'
  write(ounit,*) '@ septrace = $trace - 1'
  write(ounit,*) 'rm -f UxTrace{$trace}_sem2d.ps UzTrace{$trace}_sem2d.ps'
  write(ounit,*) 'subset < Ux_sem2d.dat n1=$nt', &
   '" if2s=$septrace n2s=1 | psgraph label1="Time (s)" label2="',trim(ylabel),&
   ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
   ' title="Ux component trace "$trace > UxTrace{$trace}_sem2d.ps'
  write(ounit,*) 'subset < Uz_sem2d.dat n1=$nt', &
   '" if2s=$septrace n2s=1 | psgraph label1="Time (s)" label2="',trim(ylabel),&
   ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
   ' title="Uz component trace "$trace > UzTrace{$trace}_sem2d.ps'
  write(ounit,*) 'end'
  close(ounit)

 100  format(A,' xcur=',F0.2,' n1=',I0,' d1=',F0.8,' label1="Time (s)" label2="',A, &
             ' (m)" n2=',I0,' x2=$x2',A )
 110  format('#!/bin/csh -f')
  
  end subroutine REC_write


!=====================================================================
!
  subroutine REC_inquire(rec,coord,isamp)

  type(rec_type), intent(in) :: rec
  double precision, pointer, optional :: coord(:,:)
  integer, optional :: isamp

  if (present(coord))  coord => rec%coord
  if (present(isamp))  isamp = rec%isamp

  end subroutine REC_inquire

end module receivers
