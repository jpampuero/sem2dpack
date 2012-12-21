module plot_gen

  use plot_visual3
  use plot_avs
  use plot_postscript
  use plot_gmt

  implicit none
  private

  integer, parameter :: nb_fields = 7
  logical, save :: selected_fields(nb_fields)
  character(nb_fields), parameter :: field_names='DVAESdc'

  integer, parameter :: nb_comps = 4
  logical, save :: selected_comps(nb_comps)
  character(nb_comps), parameter :: comp_names='xyza'

  integer, save :: ITD,IT1
  logical, save :: plot_snap,visual3,avs,ps,bin,gmt

  public :: read_plot_gen,PLOT_FIELD, PLOT_init

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : SNAP_DEF
! GROUP  : SNAPSHOT_OUTPUTS
! PURPOSE: Set preferences for exporting snapshots
! SYNTAX : &SNAP_DEF it1, itd, fields, components, bin, visual3, avs, ps, gmt /
!          Followed by a &SNAP_PS block if ps=T.
!
! ARG: it1      [int] [0]   Time step of first snapshot output
! ARG: itd      [int] [100] Number of timesteps between snapshots
! ARG: fields   [char*] ['V'] fields to export in snapshots (the prefix of the 
!                output file names is given in parenthesis):
!                 'D'     displacements (dx,dy,dz,da)
!                 'V'     velocity (vx,vy,vz,va)
!                 'A'     acceleration (ax,ay,az,aa)
!                 'E'     strain (e11,e22,e12,e23,e13)
!                 'S'     stress (s11,s22,s12,s33,e13,e23)
!                 'd'     divergence rate (dvx/dx + dvz/dz)
!                 'c'     curl rate (dvx/dz - dvz/dx)
! ARG: components [char*] ['ya'] components for PostScript outputs:
!                 in P-SV: 'x','z' and/or 'a' (amplitude). 'y' is ignored 
!                 in SH:   'y' only. Other values are ignored. 
! ARG: ps       [log] [T] PostScript (see &SNAP_PS input block)
! ARG: gmt      [log] [F] output triangulation file grid_sem2d.gmt
!                 to be used in "pscontour -T" of the General Mapping Tool (GMT)
! ARG: avs      [log] [F] AVS (only for D,V and A fields)
! ARG: visual3  [log] [F] Visual3 (only for D,V and A fields)
! ARG: bin      [log] [T] binary
!               
! NOTE   : E and S fields are exported only as binary.
!
! END INPUT BLOCK

  subroutine read_plot_gen(iin,ndof)
  
  use echo, only : iout,echo_input

  integer, intent(in) :: iin,ndof

  character(nb_fields) :: fields
  character(nb_comps) :: components
  integer :: i

  NAMELIST / SNAP_DEF / visual3,avs,gmt,ps,bin,fields,components,ITD,IT1

  ITD        = 100
  IT1        = 0
  visual3    = .false.
  avs        = .false.
  gmt        = .false.
  ps         = .true.
  bin        = .true.
  fields     = 'V'
  components = 'ya'
  
  rewind(iin)
  read(iin,SNAP_DEF,END=100) 
100 continue

  do i=1,nb_fields
    selected_fields(i) = scan(fields,field_names(i:i)) > 0
  enddo
  
  do i=1,nb_comps
    selected_comps(i) = scan(components,comp_names(i:i)) > 0
  enddo
  
  if (ndof==1) then
    selected_comps(1)=.false.
    selected_comps(2)=.true.
    selected_comps(3)=.false.
    selected_comps(4)=.false.
    selected_fields(6:7) = .false. ! no div/curl for SH
  else
    selected_comps(2)=.false.
  endif

  plot_snap = (visual3 .or. avs .or. ps .or. bin) &
            .and. any(selected_comps) .and. any(selected_fields)

  if (echo_input) write(iout,200) IT1,ITD, &
                                  ps,gmt,avs,visual3,bin, &
                                  selected_fields(:),selected_comps(:)

  if (ps) call PLOT_PS_read(iin)

  return

  200 format(//1x,'S n a p s h o t   O u t p u t s',/1x,31('='),//5x,&
  'Timestep of first snapshot output  . . . . . . (it1) = ',I0/ 5x, &
  'Number of timesteps between snapshots. . . . . (itd) = ',I0/ 5x, &
  'Save results in PS file or not . . . . . . . . .(ps) = ',L1/5x, &
  'Save grid triangulation for GMT. . . . . . . . (gmt) = ',L1/5x, &
  'Save results in AVS file or not. . . . . . . . (avs) = ',L1/5x, &
  'Save results in Visual3 file or not. . . . (visual3) = ',L1/5x, &
  'Save results in binary file or not . . . . . . (bin) = ',L1/5x, &
  'Selected fields :',/5x, &
  '  Displacement . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Velocity . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Acceleration . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Strain . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Stress . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Divergence . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Curl . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  'Selected components for PostScript snapshots :',/5x, &
  '  X  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Y  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Z  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Amplitude  . . . . . . . . . . . . . . . . . . . . = ',L1)

  end subroutine read_plot_gen


!=======================================================================

  subroutine PLOT_FIELD(pb,it,stitle,iout)

  use problem_class
  use stdio, only : IO_rw_field, IO_new_unit
  use fields_class, only : FIELD_get_elem
  use mat_gen, only : MAT_strain, MAT_stress_dv, MAT_write, MAT_divcurl
  use time_evol, only : TIME_getTime

  type(problem_type), intent(in) :: pb
  integer, intent(in) :: it
  character(*), intent(in) :: stitle
  integer, intent(in), optional :: iout

  double precision, pointer :: f(:,:)
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: d,v
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof+1) :: stress,strain
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,2) :: divcurl
  integer :: iunit(20)
  character(10) :: tag
  character(30) :: fnames(20)
  character(50) :: ps_file_name
  character :: fchar
  integer :: ndof,ngll,i,e,nssf,k,k0_stress,k_div,k_curl,iol
  integer :: comp_indx(nb_comps)

  if (.not.plot_snap) return
  if (mod(it-IT1,ITD) > 0 .or. it<IT1) return

  write(tag,'(i3.3)') (it-IT1)/ITD

  comp_indx = (/ 1,1,2,0 /)

  if (present(iout)) write(iout,300) it

  do i = 1,3

    if (.not.selected_fields(i)) cycle

   ! Select the displayed field
    select case(i)
      case(1)
        f => pb%fields%displ
        fchar = 'd'
      case(2)
        f => pb%fields%veloc
        fchar = 'v'
      case(3)
        f => pb%fields%accel
        fchar = 'a'
    end select
  
    if (pb%fields%ndof==1) then
      write(fnames(1),'(A,"y_",A)') fchar,trim(tag)
    else
      write(fnames(1),'(A,"x_",A)') fchar,trim(tag) 
      write(fnames(2),'(A,"z_",A)') fchar,trim(tag)
    endif

    if (ps) then
      do k=1,nb_comps
      if (selected_comps(k)) then
        write(ps_file_name,'( A,A,"_",A,"_sem2d.ps" )') fchar,comp_names(k:k),trim(tag)
        call PLOT_PS(file=trim(ps_file_name), vfield=f,comp=comp_indx(k) &
                  ,it_in=it,time_in=TIME_getTime(pb%time) &
                  ,grid=pb%grid,mat=pb%matpro,stitle=stitle &
                  ,src=pb%src,rec=pb%rec)
      endif
      enddo
    endif
    if (avs) call plotavs(f,pb%grid,it)
    if (visual3) call plotvisual3(f,pb%grid,pb%matpro,it)
    if (bin) call IO_rw_field(f,fnames(1:2),'w')

  enddo

 !-- other snapshots in binary output files
 !   strain, stress, div, curl

  nssf = 0  ! number of snapshot fields
  k0_stress = 0
  k_div = 0
  k_curl = 0

  if (pb%fields%ndof==1) then

    if (selected_fields(4)) then
      fnames(1) = "e13"
      fnames(2) = "e23"
      nssf=2
    endif
    if (selected_fields(5)) then
      k0_stress = nssf
      fnames(k0_stress+1) = "s13"
      fnames(k0_stress+2) = "s23"
      nssf=nssf+2
    endif

  else

    if (selected_fields(4)) then
      fnames(1) = "e11"
      fnames(2) = "e22"
      fnames(3) = "e12"
      nssf=3
    endif
    if (selected_fields(5)) then
      k0_stress = nssf
      fnames(k0_stress+1) = "s11"
      fnames(k0_stress+2) = "s22"
      fnames(k0_stress+3) = "s12"
      nssf=nssf+3
    endif
    if (selected_fields(6)) then
      k_div = nssf+1
      fnames(k_div) = "div"
      nssf=nssf+1
    endif
    if (selected_fields(7)) then
      k_curl = nssf+1
      fnames(k_curl) = "curl"
      nssf=nssf+1
    endif

  endif

  if (nssf>0 .and. bin) then

    ndof = pb%fields%ndof
    ngll = pb%grid%ngll

    inquire( IOLENGTH=iol ) real(d(:,:,1))
    do k=1,nssf
      iunit(k) = IO_new_unit()
      open(unit=iunit(k),file=trim(fnames(k))//'_'//trim(tag)//'_sem2d.dat' &
          ,status='replace',access='direct',recl=iol)
    enddo

    do e=1,pb%grid%nelem

      d = FIELD_get_elem(pb%fields%displ,pb%grid%ibool(:,:,e))
      v = FIELD_get_elem(pb%fields%veloc,pb%grid%ibool(:,:,e))

      if (selected_fields(4)) then
        strain = MAT_strain(d,pb%matwrk(e),pb%grid,e,ngll,ndof)
        do k=1,ndof+1
          write(iunit(k),rec=e) real(strain(:,:,k))
        enddo
      endif
      
      if (selected_fields(5)) then
        call MAT_stress_dv(stress,d,v,pb%matwrk(e),pb%matpro(e),pb%grid,e,ngll,ndof)
        do k=1,ndof+1
          write(iunit(k0_stress+k),rec=e) real(stress(:,:,k))
        enddo
      endif

      if (selected_fields(6) .or. selected_fields(7)) then
        divcurl = MAT_divcurl(v,pb%matwrk(e),pb%grid,e,ngll,ndof)
        !divcurl = MAT_divcurl(d,pb%matwrk(e),pb%grid,e,ngll,ndof)
        if (selected_fields(6)) write(iunit(k_div),rec=e) real(divcurl(:,:,1))
        if (selected_fields(7)) write(iunit(k_curl),rec=e) real(divcurl(:,:,2))
      endif

    enddo

    do k=1,nssf
      close(iunit(k))
    enddo

  endif
  
  call MAT_write(pb%matwrk,pb%matpro,tag)

  if (present(iout)) write(iout,*)

  return 
300 format(/"Snapshot at timestep = ",I0)

  end subroutine PLOT_FIELD
  
!=======================================================================
  subroutine PLOT_init(grid)

  use spec_grid, only : sem_grid_type

  type (sem_grid_type), intent(inout) :: grid

  if (gmt) call PLOT_GMT_init(grid%ibool)
  
  end subroutine PLOT_init

end module plot_gen
