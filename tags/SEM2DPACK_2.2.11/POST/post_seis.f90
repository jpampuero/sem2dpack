program post_seis

  implicit none

  real :: dt
  integer :: answer,it,nt,nsta,trace,irecl,i
  real, allocatable :: Udata(:,:),x(:),z(:)
  character(30) :: file_name
  character :: comp
  logical :: data_loaded

  data_loaded = .false.

  !-- Read control parameters
  open(10,file='SeisHeader_sem2d.hdr')
  read(10,*) 
  read(10,*) dt,nt,nsta
  read(10,*) 
  allocate(x(nsta),z(nsta))
  do i=1,nsta
    read(10,*) x(i),z(i)
  enddo
  close(10)

  allocate( Udata(nt,nsta) )
  INQUIRE( IOLENGTH=irecl ) Udata

  do 

    write(*,*) 
    write(*,*) 'MENU'
    write(*,*) '1. Exit'
    write(*,*) '2. Read seismograms'
    write(*,*) '3. Export one station to ASCII'
    write(*,*) '4. Export all stations to ASCII'
    write(*,*) '5. Plot all stations'
    write(*,*)
   
    read(*,*) answer
    if (answer>3 .and. .not.data_loaded) then
      write(*,*) 'No data loaded yet! First select 2 or 3.'
      cycle
    endif

    select case(answer)

      case(1)
        write(*,*) 'Bye'
        stop

      case(2)
        write(*,"(A)",advance='no') '  Component (x,y,z) = '
        read(*,*) comp
        open(unit=11,file='U'//comp//'_sem2d.dat',status='old',access='direct',recl=irecl) 
        read(11,rec=1) Udata
        close(11)
        data_loaded = .true.

      case(3)
        write(*,"(A)",advance='no') '  Trace number = '
        read(*,*) trace
        write(*,"(A)",advance='no') '  Name of output data file (ASCII) = '
        read(*,*) file_name
        open(unit=11,file=file_name,status='replace')
        write(11,"(A,I0)") '# Trace at seismogram ', trace
        do it = 1,nt
          write(11,*) (it-1)*dt, Udata(it,trace)
        enddo
        close(11)

      case(4)
        write(*,"(A)",advance='no') '  Name of output data file (ASCII) = '
        read(*,*) file_name
        open(unit=11,file=file_name,status='replace')
        do it = 1,nt
          write(11,'(1000(e14.7,1x))') (it-1)*dt,Udata(it,:)
               !this number > max nb stations ever
        enddo
        close(11)

      case(5) 
        call plot_all_traces(Udata,dt,x,z)

    end select
    
  enddo

contains

! PLOT ALL TRACES, as a function of X (horizontal axis)
  subroutine plot_all_traces(udata,dt,x,z)

  real, intent(in) :: udata(:,:),dt,x(:),z(:)

  integer :: ix,nt,offset_mode
  character(30) :: file_name
  real :: dx,umax,scal
  real :: offset(size(udata,2)),factor 
  character :: ylabel

  nt   = size(udata,1)
  nsta = size(udata,2)

  write(*,"(A)",advance='no') '  Prefix of output file (prefix.ps)  = '
  read(*,*) file_name
  write(*,"(A)",advance='no') '  Offset by (1) X, (2) Z or (3) distance = '
  read(*,*) offset_mode
  write(*,"(A)",advance='no') '  Amplitude scaling (fraction of mean offset)  = '
  read(*,*) factor

  select case(offset_mode)
    case(1) 
      offset = x
      ylabel = 'X'
    case(2)
      offset = z
      ylabel = 'Z'
    case default
      offset = sqrt((x-x(1))**2+(z-z(1))**2)
      ylabel = 'D'
  end select

  dx   = real( abs(offset(nsta)-offset(1)) )/nsta 
  umax = maxval(abs(Udata))
  scal = factor*dx/umax

  open(unit=11,file='tmp.ascii',status='replace')
  do ix = 1,nsta
    do it = 1,nt
      write(11,*) (it-1)*dt,offset(ix) + scal*udata(it,ix)
    enddo
    write(11,*)
    write(11,*)
  enddo
  close(11)

 ! GNUPLOT script
  open(13,file='plot_traces.gnu',status='replace')
  write(13,*) 'set term postscript ; set output "',trim(file_name),'.ps" '
  write(13,*) 'set xlabel "Time (secs)" ; set ylabel "',ylabel,' (m)" '
  write(13,*) 'plot "tmp.ascii" notitle w l'
!  write(13,*) 'plot [0:',dt*nt,'][',minval(offset)-factor*dx,':',&
!              maxval(offset)+factor*dx,'] "tmp.ascii" notitle w l'
  close(13)

  call system('gnuplot plot_traces.gnu ; rm -f plot_traces.gnu tmp.ascii')

  end subroutine plot_all_traces

end program post_seis
