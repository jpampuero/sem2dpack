module bc_kinflt

! WARNING: This module is still in development

  use bnd_grid, only : bnd_grid_type
  use distribution_cd
  use time_evol, only : timescheme_type
  use stf_gen

  implicit none
  private

  type input_type
    type(cd_type) :: trup,slipr,tris
  end type

  type bc_kinflt_type
    private
    double precision, dimension(:), pointer :: trup,slipr,tris,B,Z,invM1
    double precision, dimension(:,:), pointer:: T,V
    double precision :: CoefA2V
    type(input_type) :: input
    type(bnd_grid_type), pointer :: topo
    type(stf_type) :: stf
   ! for outputs:
    double precision :: ot1,odt
    integer :: comp,oit,oitd,ounit,oix1,oixn,oixd
  end type 

  public :: BC_kinflt_type, BC_kinflt_read, BC_kinflt_init, BC_kinflt_apply, BC_kinflt_write

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_KINFLT
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Prescribe a kinematic source (the spatio-temporal distribution of slip rate)
!          on a finite fault
! STATUS : The current implementation has the following features and restrictions:
!            . the fault is a flat horizontal boundary at the bottom of the model
!            . prescribed (x,t)-dependent horizontal or vertical velocity
!            . traction is free in the other component 
!            . the source time function is the same everywhere (but variably shifted and scaled)
!            . the rupture time, final slip and rise time can be spatially variable
! SYNTAX : &BC_KINFLT comp, trup|trupH, slipr|sliprH, tris|trisH, stf, ot1, otd, oxi /
!          followed, in order, by:
!          1. &DIST_XXX blocks (from the DISTRIBUTIONS group) for arguments
!             with suffix H, if present, in the order listed above.
!          2. one SOURCE TIME FUNCTION block (&STF_XXXX)
! 
! ARG: comp     [int] [1] component to apply velocity, 1=horizontal, 2=vertical
! ARG: trup     [dble] [0d0] Rupture time (s)
! ARG: slipr    [dble] [0d0] Slip rate (m/s)
! ARG: tris     [dble] [0d0] Rise time (s)
! ARG: stf      [name] [none] Name of the source time function for slip rate:
!                'RICKER', 'TAB', 'HARMONIC', 'BRUNE', 'USER', etc.
!                See the SOURCE TIME FUNCTION blocks.
! ARG: otd      [dble] [0d0] Time lag between outputs (in seconds)
!                Internally adjusted to the nearest multiple of the timestep
!                Its value can be found in the output file FltXX_sem2d.hdr
!                The default internally resets otd = timestep
! ARG: ot1      [dble] [0d0] Time of first output (in seconds)
!                Internally adjusted to the nearest multiple of the timestep
!                Its value can be found in the output file FltXX_sem2d.hdr
! ARG: oxi      [int(3)] [(1,huge,1)] First, last node and stride for output
!                The default resets oxi(2) = last fault node
! 
!
! END INPUT BLOCK

subroutine bc_kinflt_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_kinflt_type), intent(inout) :: bc
  integer, intent(in) :: iin

  character(20) :: trupH,sliprH,trisH,dt_txt,oxi2_txt
  character(15) :: stf
  double precision :: trup,slipr,tris,ot1, otd
  integer :: oxi(3),comp

  NAMELIST / BC_KINFLT / comp, trup, trupH, slipr, sliprH, tris, trisH, stf &
                        ,ot1, otd, oxi

  comp = 1
  trup = 0d0
  trupH = ''
  tris = 0d0
  trisH = ''
  slipr = 0d0
  sliprH = ''
  stf = ''

  ot1 = 0.d0
  otd = 0.d0
  oxi(1) = 1
  oxi(2) = huge(iin)
  oxi(3) = 1

  read(iin,BC_KINFLT,END=100)

  bc%comp = comp
  bc%ot1 = ot1
  bc%odt = otd
  bc%oix1 = oxi(1)
  bc%oixn = oxi(2) 
  bc%oixd = oxi(3) 

  call DIST_CD_Read(bc%input%trup,trup,trupH,iin,trupH)
  call DIST_CD_Read(bc%input%tris,tris,trisH,iin,trisH)
  call DIST_CD_Read(bc%input%slipr,slipr,sliprH,iin,sliprH)
  call STF_read(stf,bc%stf,iin)
  
  if (echo_input) then
    if (otd==0d0) then
      write(dt_txt,'(A)') 'dt'
    else
      write(dt_txt,'(EN13.3)') otd
    endif
    if (oxi(2)==huge(iin)) then
      write(oxi2_txt,'(A)') 'Last fault node'
    else
      write(oxi2_txt,'(I0)') oxi(2)
    endif
    write(iout,400) comp, trupH, trisH, sliprH,stf &
                   ,ot1,dt_txt,oxi(1),oxi2_txt,oxi(3)
  endif
 
  return

  100 call IO_abort('BC_KINFLT_read: BC_KINFLT input block not found')
  400 format(5x,'Kinematic slip source', &
            /5x,'  Component . . . . . . . . . . . . (comp) = ',I0,&
            /5x,'  Rupture time. . . . . . . . . . . (trup) = ',A,&
            /5x,'  Rise time . . . . . . . . . . . . (trup) = ',A,&
            /5x,'  Slip rate . . . . . . . . . . . .(slipr) = ',A,&
            /5x,'  Source time function. . . . . . . .(stf) = ',A, &
            /5x,'  Output first time . . . . . . . . .(ot1) = ',EN13.3,&
            /5x,'       time step  . . . . . . . . . .(otd) = ',A,&
            /5x,'       first node . . . . . . . . (oxi(1)) = ',I0,&
            /5x,'       last node  . . . . . . . . (oxi(2)) = ',A,&
            /5x,'       node stride  . . . . . . . (oxi(3)) = ',I0)

end subroutine bc_kinflt_read

!=======================================================================
!
subroutine bc_kinflt_init(bc,tag,grid,M,time,perio)

  use echo, only : echo_init,iout
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use constants, only : NDIME
  use stdio, only: IO_new_unit
  use time_evol, only: timescheme_type, TIME_getTimeStep, TIME_getNbTimeSteps, &
                       TIME_getCoefA2V
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_kinflt_type), intent(inout) :: bc
  integer, intent(in) :: tag
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: M(:,:)
  type(timescheme_type), intent(in) :: time
  type(bc_periodic_type), pointer :: perio
 
  double precision, dimension(:,:), allocatable :: bc_coord, normal
  double precision :: dt
  integer :: ndof,NDAT,NSAMP,onx,i,npoin, hunit
  character(25) :: oname,hname

  ndof = size(M,2)

 ! bc%topo => grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  npoin = bc%topo%npoin
  if (echo_init) write(iout,*) "  kinflt boundary nodes = ",npoin

  allocate( bc_coord(NDIME,npoin) )
  bc_coord = grid%coord(:,bc%topo%node)
  call DIST_CD_Init(bc%input%trup,bc_coord,bc%trup)
  call DIST_CD_Init(bc%input%tris,bc_coord,bc%tris)
  call DIST_CD_Init(bc%input%slipr,bc_coord,bc%slipr)

  allocate( normal(npoin,NDIME) )
  allocate( bc%B(npoin) )
  call BC_get_normal_and_weights(bc%topo,grid,normal,bc%B, BC_PERIO_intersects(bc%topo,perio) )

  bc%CoefA2V = TIME_getCoefA2V(time)

  allocate(bc%invM1(npoin))
  bc%invM1 = 1d0/M(bc%topo%node,1)

  allocate(bc%Z(npoin))
  bc%Z = 0.5d0/(bc%CoefA2V * bc%B * bc%invM1 )

  allocate(bc%T(npoin,ndof))
  allocate(bc%V(npoin,ndof))
  bc%T = 0d0
  bc%V = 0d0

!-- Open output files -----
! Set output parameters

  bc%oix1 = max(bc%oix1, 1)
  bc%oixn = min(bc%oixn, npoin)

! NOTE: file names limited to tags(1)<100
  write(oname,'("Flt",I2.2,"_sem2d.dat")') tag
  bc%ounit = IO_new_unit()
  !DEVEL: use direct access, and recl
  !inquire( IOLENGTH=iol ) real( bc%V(bc%oix1:bc%oixn:bc%oixd,1) )
  !open(bc%ounit,file=oname,status='replace',access='direct',recl=iol)
  open(bc%ounit,file=oname,status='replace',form='unformatted')

 ! adjust output timestep to the nearest multiple of dt:
  dt = TIME_getTimeStep(time)
  bc%oitd = max(1,nint(bc%odt/dt))
  bc%odt = dt * dble(bc%oitd)
  bc%odt = dt * bc%oitd
  bc%oit = nint(bc%ot1/dt)

 ! HEADER FILE FORMAT:
 !      Name            FltXX_sem2d.hdr where XX=tags(1) of the BC_DYNFLT
 !                      input block, the tag of the first side of the fault.
 !      Line 1          Name of size/other variables
 !      Line 2          Value of size/other variables
 !      Line 3          Name of data fields, separated by ":"
 !      Line 4          Name of coordinate axis
 !      Line 5:end      Two column table of nodes coordinates
 !
 ! OUTPUT FILE FORMAT:
 !      At each time (lag DELT, total NSAMP), NDAT data lines with
 !      NPTS columns, one per fault node, are written.
 !
 ! NOTE: when the binary file is written one line at a time
 !      a record tag is also written at the beginning and end of line,
 !      => number of columns in output file = nb of data columns +2
 !
  write(hname,'("Flt",I2.2,"_sem2d")') tag
  NDAT = ndof*2
  NSAMP = (TIME_getNbTimeSteps(time) -bc%oit)/bc%oitd +1
  hunit = IO_new_unit()
  open(hunit,file=trim(hname)//'.hdr',status='replace')
  write(hunit,*) 'NPTS NDAT NSAMP DELT'
  onx = (bc%oixn-bc%oix1)/bc%oixd +1
  write(hunit,*) onx,NDAT,NSAMP,bc%odt
  if (ndof==1) then
    write(hunit,'(A)') ' Vx:Tx'
  else
    write(hunit,'(A)') ' Vx:Vz:Tx:Tz'
  endif
  write(hunit,*) 'XPTS ZPTS'
  do i=bc%oix1,bc%oixn,bc%oixd
    write(hunit,*) bc_coord(:,i)
  enddo
  close(hunit)

  deallocate(bc_coord)

end subroutine bc_kinflt_init


!=======================================================================
!
! Subroutine bc_kinflt_apply 
!
! PURPOSE:	compute the fault boundary term B*T 
! 		with fault tractions T satisfying kinematic source boundary conditions,
! 		and add it to the internal forces
!
! INPUT: 	velocities (V_pre) and internal forces (MxA_pre) 
! 		predicted or partially updated by the time solver. 
!
! OUTPUT:	Updated fault velocities (bc%V) and tractions (bc%T)
!   			T = T_stick - Z*V
! 		where
!   			T_stick = Z * (dV_pre + coef_a2v*dA_pre)
!		Updated internal forces:
!			MxA = MxA_pre + B*T
!
! SEE ALSO:	bc_dynflt::bc_dynflt_apply
!   
subroutine bc_kinflt_apply(bc,MxA,v,time)

  type(bc_kinflt_type), intent(inout) :: bc
  double precision, intent(in) :: v(:,:)
  double precision, intent(inout) :: MxA(:,:)
  type(timescheme_type), intent(in) :: time

  integer :: i,ic

 ! WARNING: boundary is horizontal

 ! prescribed velocity (kinematic source)
  ic = bc%comp
  do i=1,bc%topo%npoin
    bc%V(i,ic) = bc%slipr(i)*STF_get(bc%stf,(time%time-bc%trup(i))/bc%tris(i))
  enddo
 ! T_stick, only x component
  bc%T(:,ic) = bc%Z * ( -2d0*v(bc%topo%node,ic) -2d0*bc%CoefA2V*bc%invM1*MxA(bc%topo%node,ic) -bc%V(:,ic))
  MxA(bc%topo%node,ic) = MxA(bc%topo%node,ic) + bc%B*bc%T(:,ic)
 ! Tz component is assumed null

end subroutine bc_kinflt_apply


!=====================================================================
!
! OUTPUT FORMAT: 
! For P-SV:	at each time, 4 data lines, one column per fault node:
! 		Vx, Vz, Tx, Tz
! For SH:	at each time, 2 data lines, one column per fault node:
! 		Vx, Tx
!
  subroutine BC_KINFLT_write(bc,itime)

  type(bc_kinflt_type), intent(inout) :: bc
  integer, intent(in) :: itime

  if ( itime < bc%oit ) return

  !DEVEL: use direct access and rec

  if (size(bc%V,2)==1) then
    write(bc%ounit) real( bc%V(bc%oix1:bc%oixn:bc%oixd,1) )
    write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,1) )
  else
    write(bc%ounit) real( bc%V(bc%oix1:bc%oixn:bc%oixd,1) )
    write(bc%ounit) real( bc%V(bc%oix1:bc%oixn:bc%oixd,2) )
    write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,1) )
    write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,2) )
  endif

  bc%oit = bc%oit + bc%oitd

  end subroutine BC_KINFLT_write

end module bc_kinflt
