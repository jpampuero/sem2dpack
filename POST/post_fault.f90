program post_process_fault

  implicit none

  integer, parameter :: ndime = 2

  real :: dt,dummy
  integer :: answer,tag,i,it,ix,npoin,ndat,nt,trace,trace_init,trace_end,trace_stride
  real, allocatable :: Udata(:,:,:),coord(:,:),ref_val(:,:)
  character(30) :: dat_file,hdr_file,file_name
  logical :: data_loaded

  allocate(coord(1,ndime))
  allocate(Udata(1,1,1) )
  data_loaded = .false.

  do 

    write(*,*) 
    write(*,*) 'FAULT MENU'
    write(*,*) '1. Exit'
    write(*,*) '2. Read fault data'
    write(*,*) '3. Extract data at single X'
    write(*,*) '4. Extract data at multiple X'
    write(*,*) '5. Extract data at single T'
    write(*,*) '6. Extract data at multiple T'
    write(*,*) '7. Convert to GNUPLOT format'
    write(*,*) '8. Add initial stress'
!    write(*,*) '9. Normalize at each X'
    write(*,*) 

    read(*,*) answer
    if (answer>2 .and. .not.data_loaded) then
      write(*,*) 'No data loaded yet! First select 2.'
      cycle
    endif

    select case(answer)

      case(1) 
        write(*,*) 'Bye'
        stop

      case(2)

        write(*,"(A)",advance='no') 'Fault tag (1-99)  = '
        read(*,*) tag

       ! Read control parameters
        write(hdr_file,'("Flt",I2.2,"_sem2d.hdr")') tag
        open(10,file=hdr_file,status='old')
        read(10,*) 
        read(10,*) npoin,ndat,nt,dt
        read(10,*) 
        read(10,*) 
        deallocate(coord)
        allocate(coord(npoin,ndime))
        do i=1,npoin
          read(10,*) coord(i,:)
        enddo
        close(10)

       ! Read fault results
        deallocate(Udata)
        allocate( Udata(npoin,ndat,nt) )
        write(dat_file,'("Flt",I2.2,"_sem2d.dat")') tag
        open(unit=11,file=dat_file,status='old',form='unformatted')
        do it=1,nt
        do i=1,ndat
          read(11) Udata(:,i,it)
        enddo
        enddo
        close(11)

        data_loaded = .true.

      case(3,4)

        write(*,*) 'About fault coords:'
        write(*,*) '  Nodes . . . . ',npoin
        write(*,*) '  Xmin Xmax . . ',minval(coord(:,1)),'  ',maxval(coord(:,1))
        write(*,*) '  Zmin Zmax . . ',minval(coord(:,2)),'  ',maxval(coord(:,2))

        if (answer == 3) then
          write(*,"(A)",advance='no') '  Node number = '
          read(*,*) trace_init
          trace_end = trace_init
          trace_stride = 1
        else
          write(*,"(A)",advance='no') '  Begin at Node number = '
          read(*,*) trace_init
          write(*,"(A)",advance='no') '  End at Node number   = '
          read(*,*) trace_end
          write(*,"(A)",advance='no') '  Node stride          = '
          read(*,*) trace_stride
        endif

        write(*,"(A)",advance='no') '  Name of output data file (ASCII) = '
        read(*,*) file_name
        open(unit=11,file=file_name,status='replace')
        write(11,"(2A)") '# Source data file = ',dat_file
        write(11,"(A)") '# Fields: 1=time 2=slip 3=slip_rate 4=shear_stress 5=normal_stress 6=strength'
        do trace = trace_init,trace_end,trace_stride
          write(11,"(A,I0,A,F0.3)") '# Fault data at node number ', trace,' , X = ',coord(trace,1)
          do it = 2,nt
            write(11,'(10(e14.7,1x))') (it-1)*dt,Udata(trace,:,it)
                       !this number > max nb fields ever
          enddo
          !write(11,"(A)") '&'
          write (11,'(X)')
        enddo
        close(11)

      case(5,6)

        write(*,*) 'About time axis:'
        write(*,*) '  Samples . . ',nt
        write(*,*) '  DT .  . . . ',dt
        write(*,*) '  Tmax  . . . ',nt*dt

        if (answer == 5) then
          write(*,"(A)",advance='no') '  Time step = '
          read(*,*) trace_init
          trace_end = trace_init
          trace_stride = 1
        else
          write(*,"(A)",advance='no') '  Begin at time step = '
          read(*,*) trace_init
          write(*,"(A)",advance='no') '  End at time step   = '
          read(*,*) trace_end
          write(*,"(A)",advance='no') '  Time step stride   = '
          read(*,*) trace_stride
        endif

        write(*,"(A)",advance='no') '  Name of output data file (ASCII) = '
        read(*,*) file_name
        open(unit=11,file=file_name,status='replace')
        write(11,"(2A)") '# Source data file = ',dat_file
        write(11,"(A)") '# Fields: 1=x 2=z 3=slip 4=slip_rate 5=shear_stress 6=normal_stress 7=strength'
        do trace = trace_init,trace_end,trace_stride
          write(11,"(A,I0,A,F0.3)") '# Fault data at time step ', trace,' , T = ',trace*dt
          do ix = 1,npoin
            write(11,'(10(e14.7,1x))') coord(ix,:), Udata(ix,:,trace)
                       !this number > max nb fields ever
          enddo
          !write(11,"(A)") '&'
          write (11,'(X)')
        enddo
        close(11)

      case(7)

        write(*,"(A)",advance='no') 'Name of output data file (ASCII) = '
        read(*,*) file_name
        open(unit=11,file=file_name,status='replace')
        write(11,"(2A)") '# Source data file = ',dat_file
        write(11,"(A)") '# Fields: 1=x 2=z 3=time 4=slip 5=slip_rate 6=shear_stress 7=normal_stress 8=strength'
        do it = 2,nt
          do ix = 1,npoin
            write(11,'(10(e14.7,1x))')  coord(ix,:), (it-1)*dt,Udata(ix,:,it)
                       !this number > max nb fields ever
          enddo
          write (11,'(X)')
        enddo
        close(11)

      case(8)
        do it=2,nt
          Udata(:,3:4,it) = Udata(:,3:4,it)+Udata(:,3:4,1)
        enddo

      case(9)

        write(*,"(A)",advance='no') 'Name of reference data file (ASCII) ='
        read(*,*) file_name
        open(unit=11,file=file_name,status='old')
        allocate(ref_val(npoin,ndat))
        do ix = 1,npoin
          read(11,*) dummy,dummy,ref_val(ix,:)
        enddo
        close(11)
        do it=1,nt
          Udata(:,:,it) = Udata(:,:,it)/ref_val
        enddo
        deallocate(ref_val)

    end select
    
  enddo

end program post_process_fault
