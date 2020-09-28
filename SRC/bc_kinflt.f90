module bc_kinflt
  
  ! The new kinematic fault module
  !
  ! Precribe fault dislocation, general fault geometry
  !

  use bnd_grid, only: bnd_grid_type
  use distribution_cd
  use time_evol, only : timescheme_type
  use stf_gen

  implicit none
  private

  type input_type
    type(cd_type) :: trup1, slipr1, tris1, trup2, slipr2, tris2
  end type

  type bc_kinflt_type
    private
    integer :: npoin
    double precision, dimension(:), pointer :: trup1, slipr1, tris1 
    double precision, dimension(:), pointer :: trup2, slipr2, tris2 
    double precision, dimension(:, :), pointer:: n1=>null(), B=>null(), &
      invM1=>null(),invM2=>null(),Z=>null(), T=>null(), V=>null(),&
      D=>null(), coord=>null()
    
    double precision :: CoefA2V, CoefA2D
    type(input_type) :: input

    integer, dimension(:), pointer :: node1=>null(), node2=>null()
    type(bnd_grid_type), pointer :: bc1 => null(), bc2 => null()
    type(stf_type), dimension(2) :: stf ! source time functions, 2 components

   ! for outputs:
    double precision :: ot1, odt,ot,odtD, odtS
    integer :: oit,oitd,ounit,oix1,oixn,oixd,ou_pot,ou_time
    integer :: otdD, otdS
    logical :: osides
  end type

  public :: BC_kinflt_type, BC_kinflt_read, BC_kinflt_init, BC_kinflt_apply, BC_kinflt_write, &
            BC_kinflt_select, BC_kinflt_trans, BC_kinflt_set

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_KINFLT
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Prescribe a kinematic source (the spatio-temporal distribution of slip rate)
!          on a finite fault. 
! 'slip rate' is prescribed in fault tangential (1) and fault normal (2) directions 
!
! STATUS : The current implementation:
!            . the source time function is the same everywhere (but variably shifted and scaled)
!                    and different for horizontal and vertical components.
!            . if stf is not specified for a certain component, it is assumed to be zero  
!            . the rupture time, slip rate and rise time can be spatially variable
!
!   slip velocity = (slip_rate) * STF((t - rupture_time)/rise_time)
!
! SYNTAX : &BC_KINFLT stf1, trup1|trup1H, slipr1|slipr1H, tris1|tris1H, 
!                     stf2, trup2|trup2H, slipr2|slipr2H, tris2|tris2H, 
!                     ot1, otd, oxi, osides/
!
!          followed, in order, by:
!          1. &DIST_XXX blocks (from the DISTRIBUTIONS group) for arguments
!             with suffix H, if present, in the order listed above.
!          2. Up to two SOURCE TIME FUNCTION block (&STF_XXXX), if present
!             . 0 STF block, both stf1 and stf2 are empty type 
!             . 1 STF block, stf1 defined, stf2 are type 
! 
! Source parameters for component 1: tangential
! ARG: stf1      [name] [none] Name of the source time function for slip rate:
!                'RICKER', 'TAB', 'HARMONIC', 'BRUNE', 'USER', etc.
!                See the SOURCE TIME FUNCTION blocks.
! ARG: trup1     [dble] [0d0] Rupture time (s)
! ARG: slipr1    [dble] [0d0] Slip rate (m/s)
! ARG: tris1     [dble] [0d0] Rise time (s)
!
! Source parameters for component 2: normal
! ARG: stf2      [name] [none] Name of the source time function for slip rate:
!                'RICKER', 'TAB', 'HARMONIC', 'BRUNE', 'USER', etc.
!                See the SOURCE TIME FUNCTION blocks.
! ARG: trup2     [dble] [0d0] Rupture time (s)
! ARG: slipr2    [dble] [0d0] Slip rate (m/s)
! ARG: tris2     [dble] [0d0] Rise time (s)
!
! ARG: otd       [dble] [0d0] Time lag between outputs (in seconds)
!                Internally adjusted to the nearest multiple of the timestep
!                Its value can be found in the output file FltXX_sem2d.hdr
!                The default internally resets otd = timestep
! ARG: ot1      [dble] [0d0] Time of first output (in seconds)
!                Internally adjusted to the nearest multiple of the timestep
!                Its value can be found in the output file FltXX_sem2d.hdr
!
! ARG: otdD     [dble] time lag between outputs dynamic 
! ARG: otdS     [dble] time lag between outputs static 

! ARG: oxi      [int(3)] [(1,huge,1)] First, last node and stride for output
!                The default resets oxi(2) = last fault node
! 
! ARG: osides   [log] [F] Export displacement and velocities on each side
!                of the fault
!
! END INPUT BLOCK


  subroutine BC_KINFLT_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_kinflt_type), intent(out) :: bc
  integer, intent(in) :: iin

  double precision :: trup1, slipr1, tris1, &
                      trup2, slipr2, tris2, ot1, otd, & 
                      otdD, otdS
  character(20) :: trup1H, slipr1H, tris1H, & 
                   trup2H, slipr2H, tris2H, dt_txt, oxi2_txt
  
  character(15) :: stf1, stf2 

  integer :: i,oxi(3) 
  logical :: osides

  NAMELIST / BC_KINFLT /  &
      stf1, trup1, trup1H, slipr1, slipr1H, tris1, tris1H, &
      stf2, trup2, trup2H, slipr2, slipr1H, tris2, tris2H, &
      ot1, otd, oxi, osides, otdD, otdS 

! default parameters:
! component 1
  trup1 = 0d0
  trup1H = ''
  tris1 = 1.0d0
  tris1H = ''
  slipr1 = 0d0
  slipr1H = ''
  stf1 = 'EMPTY'

! component 2 
  trup2 = 0d0
  trup2H = ''
  tris2 = 1.0d0
  tris2H = ''
  slipr2 = 0d0
  slipr2H = ''
  stf2 = 'EMPTY'

  ot1 = 0.d0
  otd = 0.d0

  otdD = 0.d0
  otdS = 0.d0

  oxi(1) = 1
  oxi(2) = huge(i)
  oxi(3) = 1
  osides = .false.

! read parameters from the namelist
  read(iin, BC_KINFLT, END=100)

  bc%ot1 = ot1
  bc%odt = otd
  bc%ot  = bc%ot1

  bc%otdD = otdD
  bc%otdS = otdS

  bc%oix1 = oxi(1)
  bc%oixn = oxi(2) 
  bc%oixd = oxi(3) 
  bc%osides = osides

  ! read distribution if there's any or set fields into a constant
  call DIST_CD_Read(bc%input%trup1,trup1,trup1H,iin,trup1H)
  call DIST_CD_Read(bc%input%slipr1,slipr1,slipr1H,iin,slipr1H)
  call DIST_CD_Read(bc%input%tris1,tris1,tris1H,iin,tris1H)
  
  call DIST_CD_Read(bc%input%trup2,trup2,trup2H,iin,trup2H)
  call DIST_CD_Read(bc%input%slipr2,slipr2,slipr2H,iin,slipr2H)
  call DIST_CD_Read(bc%input%tris2,tris2,tris2H,iin,tris2H)
  
  call STF_read(stf1, bc%stf(1), iin)
  call STF_read(stf2, bc%stf(2), iin)

  if (echo_input) then
     if (otdD+otdS>0) then
        write(dt_txt, '(A)') 'adaptive'
     else
        if (otd==0d0) then
          write(dt_txt,'(A)') 'dt'
        else
          write(dt_txt,'(EN13.3)') otd
         end if
    endif

    if (oxi(2)==huge(i)) then
      write(oxi2_txt,'(A)') 'Last fault node'
    else
      write(oxi2_txt,'(I0)') oxi(2)
    endif
    write(iout,200) stf1, trup1H, tris1H, slipr1H, stf2, trup2H, tris2H, slipr2H,& 
                    ot1, dt_txt, otdD, otdS, & 
					oxi(1), oxi2_txt, oxi(3), osides
  endif

  return
  100 call IO_abort('BC_KINFLT_read: BC_KINFLT input block not found')
  200 format(5x,'Kinematic fault boundary:', &
            /5x,'  Component 1 (tangential)',&
            /5x,'  Source time function. . . . . . . .(stf1)  = ',A, &
            /5x,'  Rupture time . . . . . . . . . . . (trup1) = ',A,&
            /5x,'  Rise time . . . . . . . . . . . . .(tris1) = ',A,&
            /5x,'  Slip rate . . . . . . . . . . . .(slipr1)  = ',A,&
            /5x,'  Component 2 (normal)',&
            /5x,'  Source time function. . . . . . . .(stf2)  = ',A, &
            /5x,'  Rupture time . . . . . . . . . . . (trup2) = ',A,&
            /5x,'  Rise time . . . . . . . . . . . . (tris2)  = ',A,&
            /5x,'  Slip rate . . . . . . . . . . . .(slipr2)  = ',A,&
            /5x,'  Output first time . . . . . . . . .(ot1)   = ',EN13.3,&
            /5x,'       time step  . . . . . . . . . .(otd)   = ',A,&
            /5x,'       time step dynamic. . . . . . (otdD)   = ',EN16.8,&
            /5x,'       time step static . . . . . . (otdS)   = ',EN16.8,&
            /5x,'       first node . . . . . . . . (oxi(1))   = ',I0,&
            /5x,'       last node  . . . . . . . . (oxi(2))   = ',A,&
            /5x,'       node stride  . . . . . . . (oxi(3))   = ',I0,&
            /5x,'       data from each fault side. (osides)   = ',L1)

  end subroutine BC_KINFLT_read

!=====================================================================
!
  subroutine BC_KINFLT_init(bc,tags,grid,M,time,perio)
  
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort,IO_new_unit
  use time_evol, only: timescheme_type, TIME_getTimeStep, TIME_getNbTimeSteps, &
                       TIME_getCoefA2V, TIME_getCoefA2D
  use constants, only : TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_kinflt_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: M(:,:)
  integer, intent(in) :: tags(2)
  type(timescheme_type), intent(in) :: time
  type(bc_periodic_type), pointer :: perio

  double precision, pointer :: tmp_n1(:,:), tmp_B(:)
  double precision :: dt
  integer :: i,j,k,hunit,npoin,NSAMP,NDAT,ndof,onx,ounit
  character(25) :: oname,hname
  logical :: two_sides, adapt_time
  
  ndof = size(M,2)

! if fault is on a domain boundary and symmetry is required, 
! bc2 and invM2 are not defined
  two_sides = tags(2)>0

!================== READ/PROCESS IN BOUNDARY NODES =========================

! bc1 --> grid%bounds(i) boundary corresponding to TAG
  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%bc1 )
  if ( two_sides ) then
    call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%bc2 )

    if (bc%bc1%nelem/=bc%bc2%nelem) &
     call IO_abort('bc_kinflt_init: number of boundary elements do not match')

    if (bc%bc1%npoin/=bc%bc2%npoin) &
     call IO_abort('bc_kinflt_init: number of nodes on boundaries do not match')
  endif

! nodes that are common to both sides (non-split nodes) are sticky nodes
! they must be deleted from the fault boundary
  npoin = bc%bc1%npoin
  if ( two_sides ) then
    do k=1,bc%bc1%npoin
      if (bc%bc1%node(k)==bc%bc2%node(k))  npoin = npoin-1
    enddo
    if (npoin<bc%bc1%npoin) then
      allocate( bc%node1(npoin) )
      allocate( bc%node2(npoin) )
      j = 0
      do k=1,bc%bc1%npoin
        if (bc%bc1%node(k)/=bc%bc2%node(k)) then
          j=j+1
          bc%node1(j) = bc%bc1%node(k)
          bc%node2(j) = bc%bc2%node(k)
        endif
      enddo
    else
      bc%node1 => bc%bc1%node
      bc%node2 => bc%bc2%node
    endif
  else
    bc%node1 => bc%bc1%node
  endif
  bc%npoin = npoin

  allocate(bc%coord(2,npoin))
  bc%coord = grid%coord(:,bc%node1)
  if ( two_sides ) then
    if ( any(abs(bc%coord(1,:)-grid%coord(1,bc%node2))>TINY_XABS) &
      .OR.any(abs(bc%coord(2,:)-grid%coord(2,bc%node2))>TINY_XABS) )&
      call IO_abort('bc_kinflt_init: coordinates on boundaries do not match properly')
  endif

  call DIST_CD_Init(bc%input%trup1,bc%coord,bc%trup1)
  call DIST_CD_Init(bc%input%tris1,bc%coord,bc%tris1)
  call DIST_CD_Init(bc%input%slipr1,bc%coord,bc%slipr1)
  call DIST_CD_Init(bc%input%trup2,bc%coord,bc%trup2)
  call DIST_CD_Init(bc%input%tris2,bc%coord,bc%tris2)
  call DIST_CD_Init(bc%input%slipr2,bc%coord,bc%slipr2)

! NOTE: the mesh being conformal, the weights B=GLL_weights*jac1D are equal on both
!       sides of the fault. 
  allocate( bc%n1(npoin,2) )
  allocate( bc%B(npoin,ndof) )
  if (npoin<bc%bc1%npoin) then
    allocate( tmp_n1(bc%bc1%npoin,2) )
    allocate( tmp_B(bc%bc1%npoin) ) ! assembled[ GLL_weights * jac1D ]
    call BC_get_normal_and_weights(bc%bc1,grid,tmp_n1,tmp_B, BC_PERIO_intersects(bc%bc1,perio) )
    j = 0
    do k=1,bc%bc1%npoin
      if (bc%bc1%node(k)/=bc%bc2%node(k)) then
        j=j+1
        bc%B(j,1) = tmp_B(k)
        bc%n1(j,:) = tmp_n1(k,:)
      endif
    enddo
    deallocate(tmp_n1,tmp_B)
  else
    call BC_get_normal_and_weights(bc%bc1,grid,bc%n1,bc%B(:,1), BC_PERIO_intersects(bc%bc1,perio) )
  endif
  bc%B(:,ndof) = bc%B(:,1)

! Coefficients of corrector phase (see solver.f90)
!  vnew = vpredictor + coefA2V *anew   
!  dnew = dpredictor + coefA2D *anew
  bc%CoefA2V = TIME_getCoefA2V(time)
  bc%CoefA2D = TIME_getCoefA2D(time)

! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1

! Inverse of mass matrix
  allocate(bc%invM1(npoin,ndof))
  bc%invM1 = 1d0/M(bc%node1,:)
  if ( two_sides ) then
    allocate(bc%invM2(npoin,ndof))
    bc%invM2 = 1d0/M(bc%node2,:)
  endif

! Fault impedance, Z in :  Trac=T_Stick-Z*dV
!   Z = 1/( B1/M1 + B2/M2 )
! T_Stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity) 
! NOTE: same Bi on both sides, see note above
  allocate(bc%Z(npoin,ndof))
  if ( two_sides ) then
    bc%Z = 1.d0/(bc%CoefA2V * bc%B *( bc%invM1 + bc%invM2 ))
  else
    bc%Z = 0.5d0/(bc%CoefA2V * bc%B * bc%invM1 )
  endif

  allocate(bc%T(npoin,2))
  allocate(bc%D(npoin,ndof))
  allocate(bc%V(npoin,ndof))

  bc%T = 0d0
  bc%D = 0d0
  bc%V = 0d0 

!-- Open output files -----
! Set output parameters

  bc%oix1 = max(bc%oix1, 1)
  bc%oixn = min(bc%oixn, npoin)
  adapt_time =  time%kind=='adaptive'

! NOTE: file names limited to tags(1)<100
  write(oname,'("Flt",I2.2,"_sem2d.dat")') tags(1)
  bc%ounit = IO_new_unit()
  !DEVEL: use direct access, and recl
  !inquire( IOLENGTH=iol ) real( bc%D(bc%oix1:bc%oixn:bc%oixd,1) )
  !open(bc%ounit,file=oname,status='replace',access='direct',recl=iol)
  open(bc%ounit,file=oname,status='replace',form='unformatted')

  ! File to store the output time information if adaptive time stepping
  if (adapt_time) then
      write(oname,'("Flt",I2.2,"_time_sem2d.tab")') tags(1)
      bc%ou_time = IO_new_unit()
      open(bc%ou_time, file=oname, status='replace')
  end if

 ! adjust output timestep to the nearest multiple of dt:
  dt = TIME_getTimeStep(time)

  if (.not. adapt_time) then
      bc%oitd = max(1,nint(bc%odt/dt))
      bc%odt = dt * dble(bc%oitd)
      bc%odt = dt * bc%oitd
      bc%oit = nint(bc%ot1/dt)
      bc%ot  = dble(bc%oit) * dt
  end if

 ! HEADER FILE FORMAT:
 !      Name            FltXX_sem2d.hdr where XX=tags(1) of the BC_KINFLT
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

 !
 ! WARNING: 
 ! Header file is different if adaptive time stepping is used
 !
 ! NSAMP and DELT is not known before the end of simulation
 !
 ! And should be estimated from the FltXX_time_sem2d.tab file
 !
 
  write(hname,'("Flt",I2.2,"_sem2d")') tags(1)
  NDAT = 4
  if (bc%osides) NDAT = NDAT + 4*ndof
  NSAMP = (TIME_getNbTimeSteps(time) -bc%oit)/bc%oitd +1
  hunit = IO_new_unit()

  onx = (bc%oixn-bc%oix1)/bc%oixd +1
  open(hunit,file=trim(hname)//'.hdr',status='replace')

  if (adapt_time) then
      write(hunit,*) 'NPTS NDAT'
      write(hunit,*) onx, NDAT
  else
      NSAMP = (TIME_getNbTimeSteps(time) -bc%oit)/bc%oitd +1
      write(hunit,*) 'NPTS NDAT NSAMP DELT'
      write(hunit,*) onx,NDAT,NSAMP,bc%odt
  end if

  if (bc%osides) then 
    if (ndof==1) then
      write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:D1t:D2t:V1t:V2t'
    else
      write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:D1t:D1n:D2t:D2n:V1t:V1n:V2t:V2n'
    endif
  else
    write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress'
  endif
  write(hunit,*) 'XPTS ZPTS'
  do i=bc%oix1,bc%oixn,bc%oixd
    write(hunit,*) bc%coord(:,i)
  enddo
  close(hunit)

  end subroutine BC_KINFLT_init

!=====================================================================
!
! Apply kinematic boundary condition, different in dynamic and static
!
! Computes the boundary term B*T, 
! with fault tractions T satisfying kinematic source boundary conditions, 
! and adds it to the internal forces.
!
! On entry to this subroutine, displacements (D), velocities (V)
! and internal forces (MxA) have values predicted or partially updated
! by the time solver. 
!
! The subroutine computes the fault tractions T that satisfy
! the contact and friction conditions and the discrete relation
!   T = T_stick - Z*dV
! where Z is the discrete fault impedance
!   Z = 1/( B1/M1 + B2/M2 ) / coef_a2v
! (B and M are diagonal boundary and mass marices 
! and coef_a2v depends on the time integration scheme),
! dV is the slip rate and T_stick is the traction that would prevail 
! if the fault was suddenly stress free:
!   T_stick = Z * dV_free
!           = Z * (dV_predictor + coef_a2v*dA_free)
! where: 
!   dV_free = slip velocity if the fault was suddenly stress free
!   dA_free = slip acceleration if the fault was suddenly stress free
!             (obtained from MxA_predictor)
!
! On exit, the new values should be
!   MxA = MxA_predictor + B*T
!   D = D_predictor + coef_a2d*MxA/M
!   V = V_predictor + coef_a2v*MxA/M
! where the coefficients coef_xxx depend on the time scheme.
! However, this subroutine updates only MxA. The update of D and V is done 
! elsewhere (see solver.f90).
!
! Convention:   T = sigma*n1
!               V = v2-v1
! => T and V have same sign
! 1: minus side
! 2: plus  side
!

subroutine BC_KINFLT_apply(bc,MxA,V,D,time)
      ! MxA: force
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use constants, only : PI
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_kinflt_type), intent(inout) :: bc
  double precision, intent(inout) :: V(:,:)
  double precision, intent(inout) :: D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  if (.not. time%isDynamic) then
      call BC_KINFLT_apply_quasi_static(bc,MxA,V,D,time)
  else
      call BC_KINFLT_apply_dynamic(bc,MxA,V,D,time)
  end if

end subroutine BC_KINFLT_apply

! ==============================================================================
!
!  Subroutine bc_kinflt_apply_quasi_static
!
! Apply kinematic condition in quasi-statics 
! following Seki, 2017; Mekosh and Raefsky, 1981.
!
! Transform V into half slip and average velocity
! overwrite the half slip by kinematic prescribed values and back transform
!
! Update displacement accordingly 
!

subroutine BC_KINFLT_apply_quasi_static(bc,MxA,V,D,time)
      ! MxA: force or -K*d
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_kinflt_type), intent(inout) :: bc
  double precision, intent(inout) :: V(:,:), D(:,:), MxA(:,:)
  integer :: ndof, i

  ndof = size(MxA,2)

  ! Evaluate slip rate
  do i=1,bc%npoin
      bc%V(i,1) = bc%slipr1(i)*STF_get(bc%stf(1),(time%time-bc%trup1(i))/bc%tris1(i))
  enddo
  
  if (ndof==2) then
      do i=1,bc%npoin
          bc%V(i,2) = bc%slipr2(i)*STF_get(bc%stf(2),(time%time-bc%trup2(i))/bc%tris2(i))
      enddo
  end if

  ! update slip
  bc%D = bc%D + time%dt * bc%V

  ! rotate D, V back into x, y coordinate 
  if (ndof==2) then
      bc%V = rotate(bc, bc%V, -1)
      bc%D = rotate(bc, bc%D, -1)
  end if

  ! transform global velocity, and displacement
  call BC_KINFLT_trans(bc, V, 1)
  call BC_KINFLT_trans(bc, D, 1)

  ! only update node1, stores half slip
  ! node2 now stores the average
  V(bc%node1, :) = bc%V/2.0d0
  D(bc%node1, :) = bc%D/2.0d0

  ! transform back
  call BC_KINFLT_trans(bc, V, -1)
  call BC_KINFLT_trans(bc, D, -1)
  
  ! rotate V,D back to fault tangential/normal coordinate
  if (ndof==2) then
      bc%V = rotate(bc, bc%V, 1)
      bc%D = rotate(bc, bc%D, 1)
  end if

  ! apply symmetry to normal stress when needed
  if (.not.associated(bc%bc2) .or. ndof==1) bc%T(:,2)=0d0 

  ! fault tangential traction is not updated, unknown at this point
  ! only known after the quasi-static solve

end subroutine BC_KINFLT_apply_quasi_static

!=======================================================================
!
! Subroutine bc_kinflt_apply_dynamic 
!
! PURPOSE:  compute the fault boundary term B*T 
!       with fault tractions T satisfying kinematic source boundary conditions,
!       and add it to the internal forces
!
! INPUT:    velocities (V_pre) and internal forces (MxA_pre) 
!       predicted or partially updated by the time solver. 
!
! OUTPUT:   Updated fault velocities (bc%V) and tractions (bc%T)
!               T = T_stick - Z*V
!       where
!               T_stick = Z * (dV_pre + coef_a2v*dA_pre)
!       Updated internal forces:
!           MxA = MxA_pre + B*T
!
! SEE ALSO: bc_dynflt::bc_dynflt_apply
!   
! Kinematic faults are treated as traction boundary condition
!
! The traction is computed using prescribed slip, this traction 
! then automatically leads to the relative velocity between 
! two split nodes 
!
! In this sense, the dynamic kinematic boundary condition is weakly
! enforeced.

subroutine BC_KINFLT_apply_dynamic(bc,MxA,V,D,time)
      ! MxA: force
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use constants, only : PI
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_kinflt_type), intent(inout) :: bc
  double precision, intent(inout) :: V(:,:),D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%npoin,size(V,2)) :: T, dV, dA
  integer :: ndof, i

  ndof = size(MxA,2)
  
  ! get predicted values
  dV = get_jump(bc,V) ! dV_predictor
  dA = get_weighted_jump(bc,MxA) ! dA_free

  ! T_stick
  T(:, 1:ndof) = bc%Z * ( dV + bc%CoefA2V*dA )

  ! Evaluate slip rate
  do i=1,bc%npoin
      bc%V(i,1) = bc%slipr1(i)*STF_get(bc%stf(1),(time%time-bc%trup1(i))/bc%tris1(i))
  enddo
  
  if (ndof==2) then
      do i=1,bc%npoin
          bc%V(i,2) = bc%slipr2(i)*STF_get(bc%stf(2),(time%time-bc%trup2(i))/bc%tris2(i))
      enddo
  end if

  ! WARNING
  ! update slip, may not be consistent with all explicit time schemes
  bc%D = bc%D + time%dt * bc%V
  
  ! rotate V back into x, y coordinate 
  if (ndof==2) bc%V = rotate(bc, bc%V, -1)
  
  ! update traction in x, y coordinate
  T(:, 1:ndof) = T(:, 1:ndof) - bc%Z * bc%V(:, 1:ndof)
  
  ! apply symmetry to normal stress when needed
  if (.not.associated(bc%bc2) .or. ndof==1) T(:, 2)=0d0 

  ! Save tractions
  bc%T = T
  
  ! Rotate tractions back to (x,z) frame 
  if (ndof==2) T = rotate(bc,T,-1)

  ! Add boundary term B*T to M*a
  MxA(bc%node1,:) = MxA(bc%node1,:) + bc%B * T(:,1:ndof)
  if (associated(bc%bc2)) MxA(bc%node2,:) = MxA(bc%node2,:) - bc%B*T(:,1:ndof)

  end subroutine BC_KINFLT_apply_dynamic

!---------------------------------------------------------------------
  function get_jump (bc, v) result(dv)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(in) :: v(:,:)
  double precision :: dv(bc%npoin,size(v,2))

  if (associated(bc%bc2)) then
    dv = v(bc%node2,:)-v(bc%node1,:)
  else
    dv = -2d0*v(bc%node1,:)
  endif

  end function get_jump

!---------------------------------------------------------------------
  function get_weighted_jump (bc,f) result(da)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(in) :: f(:,:)
  double precision :: da(bc%npoin,size(f,2))

  if (associated(bc%bc2)) then
    da = bc%invM2*f(bc%node2,:)-bc%invM1*f(bc%node1,:)
  else
    da = -2d0*bc%invM1*f(bc%node1,:)
  endif

  end function get_weighted_jump

!---------------------------------------------------------------------
  function rotate(bc,v,fb) result(vr)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(in) :: v(bc%npoin,2)
  integer, intent(in) :: fb
  double precision :: vr(bc%npoin,2)

 ! forward rotation
  if (fb==1) then
    vr(:,1) = bc%n1(:,2)*v(:,1) - bc%n1(:,1)*v(:,2)
    vr(:,2) = bc%n1(:,1)*v(:,1) + bc%n1(:,2)*v(:,2)
 ! backward rotation
  else
    vr(:,1) = bc%n1(:,2)*v(:,1) + bc%n1(:,1)*v(:,2)
    vr(:,2) =-bc%n1(:,1)*v(:,1) + bc%n1(:,2)*v(:,2)
  endif

  end function rotate


!=====================================================================
! OUTPUT FORMAT: at each time, 4 data lines, one column per fault node
! 1: Slip
! 2: Slip rate
! 3: Shear stress 
! 4: Normal stress 
!
  subroutine BC_KINFLT_write(bc,time,d,v)

  use time_evol, only : timescheme_type
  type(bc_kinflt_type), intent(inout) :: bc
  type(timescheme_type), intent(in) :: time 
  double precision, dimension(:,:), intent(in) :: d,v


  ! reset ot for each step that a switch in time scheme takes place
  if (time%switch) bc%ot = time%time
  if ( time%time < bc%ot ) return

  !DEVEL: use direct access, and rec
  !D,V,T are fault tangential and fault normal directions
  write(bc%ounit) real( bc%D(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%V(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,2) )

  if (bc%osides) then
    call export_side(bc,get_side(bc,d,1))
    call export_side(bc,get_side(bc,d,2))
    call export_side(bc,get_side(bc,v,1))
    call export_side(bc,get_side(bc,v,2))
  endif

  ! update the next output time step index
  if (.not. time%kind=='adaptive') then
      bc%ot = bc%ot + bc%odt
  else
      ! adaptive time stepping
      if (time%isDynamic) then
          bc%ot = bc%ot + bc%odtD
      else
          bc%ot = bc%ot + bc%odtS
      end if
      ! write time information if adaptive time
      write(bc%ou_time, *)  time%it, time%dt, time%time, &
                        time%EQNum, time%isDynamic, time%switch
  end if

  end subroutine BC_KINFLT_write

  !----------
  function get_side(bc,d,side) result(delta)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, dimension(:,:) :: d
  integer, intent(in) :: side
  double precision :: delta(bc%npoin,size(d,2))

  if (side==1) then
    delta = d(bc%node1,:)
  elseif (associated(bc%bc2)) then
    delta = d(bc%node2,:)
  else  
    delta = -d(bc%node1,:)
  endif

  end function get_side

  !----------
  subroutine export_side(bc,delta)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, dimension(:,:) :: delta

  if (size(delta,2)==1) then 
    write(bc%ounit) real( delta(bc%oix1:bc%oixn:bc%oixd,1) )
  else
    delta = rotate(bc,delta,1)
    write(bc%ounit) real( delta(bc%oix1:bc%oixn:bc%oixd,1) )
    write(bc%ounit) real( delta(bc%oix1:bc%oixn:bc%oixd,2) )
  endif

  end subroutine export_side

!=====================================================================
! set values along boundary nodes with scalar input

subroutine BC_KINFLT_set(bc, field, input, side_in)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)
  double precision, intent(in) :: input
  integer, intent(in), optional :: side_in 
  integer :: side 
  
  side  = 2 ! default, select all available sides
  if (present(side_in)) side = side_in

  ! side = -1, node1
  ! side =  1, node2
  ! side =  2, node1 + node2

  select case(side)
  case (-1)
      field(bc%node1,:) = input
      return
  case (1)
      if (associated(bc%bc2)) then
          field(bc%node2, :) = input
      end if
      return
  case (2)
      field(bc%node1,:) = input
      if (associated(bc%bc2)) then
          field(bc%node2, :) = input
      end if
  end select

end subroutine BC_KINFLT_set

! ======================================================
! set values along boundary nodes with array input
subroutine BC_KINFLT_set_array(bc, field, input, side_in)

  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(inout) :: field(:,:)
  double precision, intent(in) :: input(:,:)
  integer, intent(in), optional :: side_in 
  integer :: side 
  
  side  = 2 ! default, select all available sides
  if (present(side_in)) side = side_in

  ! side = -1, node1
  ! side =  1, node2
  ! side =  2, node1 + node2

  select case(side)
  case (-1)
      field(bc%node1,:) = input(bc%node1,:)
      return
  case (1)
      if (associated(bc%bc2)) then
          field(bc%node2, :) = input(bc%node2,:)
      end if
      return
  case (2)
      field(bc%node1,:) = input(bc%node1,:)
      if (associated(bc%bc2)) then
          field(bc%node2, :) = input(bc%node2,:)
      end if
  end select

end subroutine BC_KINFLT_set_array
 
! select fields only on dynamic fault boundary and zero other components
subroutine BC_KINFLT_select(bc, field_in, field_out, side_in)
   
  use stdio, only: IO_abort
  
  type(bc_kinflt_type), intent(in) :: bc
  double precision, intent(in) :: field_in(:,:)
  double precision, dimension(:,:), intent(inout) :: field_out
  integer, intent(in), optional :: side_in 
  integer :: side 

  side  = 2 ! default, select all available sides
  if (present(side_in)) side = side_in

  ! side = -1, node1
  ! side =  1, node2
  ! side =  2, node1 + node2

  select case(side)
  case(-1)
      field_out(bc%node1,:) = field_in(bc%node1, :)
      return
  case(1)
      if (associated(bc%bc2)) then
      field_out(bc%node2,:) = field_in(bc%node2, :)
      end if
      return
  case(2)
      field_out(bc%node1,:) = field_in(bc%node1, :)
      if (associated(bc%bc2)) then
          field_out(bc%node2, :) = field_in(bc%node2, :)
      end if
  end select
  
  end subroutine BC_KINFLT_select

! ==============================================================
! Transform a field on the fault boundary
! direction = 1
! 
! fout(node1) = (fin(node2) - fin(node1))/2    Half jump
! fout(node2) = (fin(node2) + fin(node1))/2    Average
!
! After the transformation, node1 stores half jump 
! and node2 stores the average between two split nodes
!
! If symmetry is assumed, then node2 is not stored
! fin(node2) = - fin(node1) 
!
! After transform:
! fout(node1) = - fin(node1)                 Half jump
!
! direction = -1 perform the inverse transformation
!
! This is used to impose slip on the fault
!
! See master thesis of Junpei Seki
!

subroutine BC_KINFLT_trans(bc, field, direction)
  ! note that field is directly modified
  type(bc_kinflt_type), intent(in) :: bc
  double precision, dimension(:,:), intent(inout) :: field
  double precision, dimension(:,:), pointer :: tmp1, tmp2
  integer, intent(in) :: direction
  integer :: i, Nbcnode, ndof

  ndof    = size(field, 2)
  Nbcnode = size(bc%node1, 1)

  allocate(tmp1(Nbcnode, ndof)) 
  allocate(tmp2(Nbcnode, ndof)) 

  tmp1    = field(bc%node1, :)  
  tmp2    = 0d0

  if (associated(bc%bc2)) then
        tmp2 = field(bc%node2, :)
  end if

  select case(direction)
  case (1)
      ! forward transform
      if (associated(bc%bc2)) then
          ! two sided fault
          field(bc%node1, :) = 0.5d0 * (tmp2 - tmp1)
          field(bc%node2, :) = 0.5d0 * (tmp2 + tmp1)
      else
          ! one sided fault
          field(bc%node1,:) = - field(bc%node1, :)
      end if
  case (-1)
      ! backward transform
      if (associated(bc%bc2)) then
          ! two sided fault
          field(bc%node1, :) = tmp2 - tmp1 
          field(bc%node2, :) = tmp2 + tmp1 
      else
          ! one sided fault
          field(bc%node1,:) = - field(bc%node1, :)
      end if
  end select

  deallocate(tmp1, tmp2)
end subroutine BC_KINFLT_trans


end module bc_kinflt
  
