module bc_dynflt
! slip weakening friction fault

  use bnd_grid, only: bnd_grid_type
  use distribution_cd
  use bc_dynflt_normal
  use bc_dynflt_swf
  use bc_dynflt_rsf
  use bc_dynflt_twf

  implicit none
  private

  type bc_dynflt_input_type
    type(cd_type) :: T,N,Sxx,Sxy,Sxz,Syz,Szz,cohesion,V
  end type bc_dynflt_input_type

  type bc_dynflt_type
    private
    integer :: npoin
    integer, dimension(:), pointer :: node1=>null(), node2=>null()
    double precision :: CoefA2V,CoefA2D
    double precision, dimension(:,:), pointer:: n1=>null(),B=>null(), &
      invM1=>null(),invM2=>null(),Z=>null(),T0=>null(),T=>null(),V=>null(),&
      D=>null(),coord=>null()
    double precision, dimension(:), pointer:: MU=>null(),cohesion=>null()
    type(swf_type), pointer :: swf => null()
    type(rsf_type), pointer :: rsf => null()
    type(twf_type), pointer :: twf => null()
    logical :: allow_opening
    type(normal_type) :: normal
    type(bnd_grid_type), pointer :: bc1 => null(), bc2 => null()
    type(bc_dynflt_input_type) :: input
   ! for outputs:
    double precision :: ot1,odt
    integer :: oit,oitd,ounit,oix1,oixn,oixd, ou_pot
    logical :: osides
  end type bc_dynflt_type

  public :: BC_DYNFLT_type, BC_DYNFLT_read, BC_DYNFLT_init, BC_DYNFLT_apply, BC_DYNFLT_write, BC_DYNFLT_set, BC_DYNFLT_timestep

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT
! GROUP  : BOUNDARY_CONDITION, DYNAMIC_FAULT
! PURPOSE: Dynamic fault with friction
! SYNTAX : &BC_DYNFLT friction, cohesion|cohesionH, opening, Tn|TnH, Tt|TtH,
!                     Sxx|SxxH, Sxy|SxyH, Sxz|SxzH, Syz|SyzH, Szz|SzzH
!                     ot1, otd, oxi, osides /
!          followed, in order, by:
!          1. &DIST_XXX blocks (from the DISTRIBUTIONS group) for arguments
!             with suffix H, if present, in the order listed above.
!          2. &BC_DYNFLT_SWF, &BC_DYNFLT_TWF or &BC_DYNFLT_RSF block(s) 
!             (if absent, default values are used)
!          3. &BC_DYNFLT_NOR block (if absent, default values are used)
!
! ARG: friction [name(2)] ['SWF',''] Friction law type:
!                  SWF = slip weakening friction
!                  TWF = time weakening friction
!                  RSF = rate and state dependent friction
!                Some friction types can be combined. E.g. to set the 
!                friction coefficient to the minimum of SWF and TWF, set 
!                  friction='SWF','TWF'
! ARG: cohesion [dble] [0d0] part of the strength that is not proportional to 
!                normal stress. It must be positive or zero.
! ARG: opening  [log] [T] Allow fault opening instead of tensile normal stress
! ARG: Tn       [dble] [0d0] Initial normal traction (positive = tensile)
! ARG: Tt       [dble] [0d0] Initial tangent traction 
!                (positive antiplane: y>0; positive inplane: right-lateral slip)
! ARG: Sxx      [dble] [0d0] Initial stress sigma_xx
! ARG: Sxy      [dble] [0d0] Initial stress sigma_xy
! ARG: Sxz      [dble] [0d0] Initial stress sigma_xz
! ARG: Syz      [dble] [0d0] Initial stress sigma_yz
! ARG: Szz      [dble] [0d0] Initial stress sigma_zz
! ARG: otd      [dble] [0.d0] Time lag between outputs (in seconds).
!                Internally adjusted to the nearest multiple of the timestep.
!                Its value can be found in the output file FltXX_sem2d.hdr.
!                The default internally resets otd = timestep
! ARG: ot1      [dble] [0.d0] Time of first output (in seconds).
!                Internally adjusted to the nearest multiple of the timestep.
!                Its value can be found in the output file FltXX_sem2d.hdr
! ARG: oxi      [int(3)] [(1,huge,1)] First, last node and stride for output.
!                The default resets oxi(2) = last fault node
! ARG: osides   [log] [F] Export displacement and velocities on each side
!                of the fault
! ARG: V        [dble] [1d-12] Initial velocity (needed for RSF)
!
! NOTE: The initial stress can be set as a stress tensor (Sxx,etc), as
!       initial tractions on the fault plane (Tn and Tt) or as the sum of both.
!
! NOTE: We recommend to use dynamic faults with the leapfrog time scheme
!       and a layer of Kelvin-Voigt damping material near the fault.
!
! END INPUT BLOCK


  subroutine BC_DYNFLT_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_dynflt_type), intent(out) :: bc
  integer, intent(in) :: iin
  double precision :: cohesion,Tt,Tn,ot1,otd,Sxx,Sxy,Sxz,Syz,Szz,V
  character(20) :: TtH,TnH, SxxH,SxyH,SxzH,SyzH,SzzH &
                  ,dt_txt,oxi2_txt, cohesionH, VH
  character(3) :: friction(2)
  integer :: i,oxi(3)
  logical :: opening,osides

  NAMELIST / BC_DYNFLT /  Tt,Tn,Sxx,Sxy,Sxz,Syz,Szz &
                         ,TtH,TnH,SxxH,SxyH,SxzH,SyzH,SzzH &
                         ,ot1,otd,oxi,osides, friction, opening &
                         ,cohesion, cohesionH, V, VH

  Tt = 0d0
  Tn = 0d0
  Sxx = 0d0
  Sxy = 0d0
  Sxz = 0d0
  Syz = 0d0
  Szz = 0d0

  TtH = ''
  TnH = ''
  SxxH = ''
  SxyH = ''
  SxzH = ''
  SyzH = ''
  SzzH = ''

  V = 1d-12
  VH = ''

  ot1 = 0.d0
  otd = 0.d0
  oxi(1) = 1
  oxi(2) = huge(i)
  oxi(3) = 1
  osides = .false.

  friction(1) = 'SWF'
  friction(2) = ''

  cohesion = 0d0
  cohesionH = ''
  
  opening = .true.

  read(iin,BC_DYNFLT,END=100)

  bc%ot1 = ot1
  bc%odt = otd
  bc%oix1 = oxi(1)
  bc%oixn = oxi(2) 
  bc%oixd = oxi(3) 
  bc%osides = osides

  call DIST_CD_Read(bc%input%cohesion,cohesion,cohesionH,iin,cohesionH)
  call DIST_CD_Read(bc%input%N,Tn,TnH,iin,TnH)
  call DIST_CD_Read(bc%input%T,Tt,TtH,iin,TtH)
  call DIST_CD_Read(bc%input%Sxx,Sxx,SxxH,iin,SxxH)
  call DIST_CD_Read(bc%input%Sxy,Sxy,SxyH,iin,SxyH)
  call DIST_CD_Read(bc%input%Sxz,Sxz,SxzH,iin,SxzH)
  call DIST_CD_Read(bc%input%Syz,Syz,SyzH,iin,SyzH)
  call DIST_CD_Read(bc%input%Szz,Szz,SzzH,iin,SzzH)
  call DIST_CD_Read(bc%input%V,V,VH,iin,VH)

  bc%allow_opening = opening

  if (echo_input) then
    if (otd==0d0) then
      write(dt_txt,'(A)') 'dt'
    else
      write(dt_txt,'(EN13.3)') otd
    endif
    if (oxi(2)==huge(i)) then
      write(oxi2_txt,'(A)') 'Last fault node'
    else
      write(oxi2_txt,'(I0)') oxi(2)
    endif
    write(iout,200) TnH,TtH,SxxH,SxyH,SxzH,SyzH,SzzH,VH, & 
                    cohesionH,opening,ot1,dt_txt,oxi(1),oxi2_txt,oxi(3),osides
  endif

  do i=1,2
    select case (friction(i))
    case ('SWF')
      allocate(bc%swf)
      call swf_read(bc%swf,iin)
    case ('RSF')
      allocate(bc%rsf)
      call rsf_read(bc%rsf,iin)
    case ('TWF')
      allocate(bc%twf)
      call twf_read(bc%twf,iin)
    case ('')
    case default
      call IO_abort('BC_DYNFLT: invalid friction')
    end select
  enddo

  call normal_read(bc%normal,iin)

  return
  100 call IO_abort('BC_DYNFLT_read: BC_DYNFLT input block not found')
  200 format(5x,'Initial traction normal. . . . . . . .(Tn) = ',A,&
            /5x,'                 tangent . . . . . . .(Tt) = ',A,&
            /5x,'Initial stress xx . . . . . . . . . .(Sxx) = ',A,&
            /5x,'               xy . . . . . . . . . .(Sxy) = ',A,&
            /5x,'               xz . . . . . . . . . .(Sxz) = ',A,&
            /5x,'               yz . . . . . . . . . .(Syz) = ',A,&
            /5x,'               zz . . . . . . . . . .(Szz) = ',A,&
            /5x,'Initial velocity on boundary . . . . . (V) = ',A,&
            /5x,'Cohesion  . . . . . . . . . . . (cohesion) = ',A,&
            /5x,'Allow opening . . . . . . . . . .(opening) = ',L1,&
            /5x,'Output first time . . . . . . . . . .(ot1) = ',EN13.3,&
            /5x,'       time step  . . . . . . . . . .(otd) = ',A,&
            /5x,'       first node . . . . . . . . (oxi(1)) = ',I0,&
            /5x,'       last node  . . . . . . . . (oxi(2)) = ',A,&
            /5x,'       node stride  . . . . . . . (oxi(3)) = ',I0,&
            /5x,'       data from each fault side. (osides) = ',L1)

  end subroutine BC_DYNFLT_read

!=====================================================================
!
  subroutine BC_DYNFLT_init(bc,tags,grid,M,time,perio)
  
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort,IO_new_unit
  use time_evol, only: timescheme_type, TIME_getTimeStep, TIME_getNbTimeSteps, &
                       TIME_getCoefA2V, TIME_getCoefA2D
  use constants, only : TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_dynflt_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: M(:,:)
  integer, intent(in) :: tags(2)
  type(timescheme_type), intent(in) :: time
  type(bc_periodic_type), pointer :: perio

  double precision, dimension(:), pointer :: Tt0,Tn0,Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz,nx,nz,V
  double precision, pointer :: tmp_n1(:,:), tmp_B(:)
  double precision :: dt
  integer :: i,j,k,hunit,npoin,NSAMP,NDAT,ndof,onx,ounit
  character(25) :: oname,hname
  logical :: two_sides
  
  nullify(Tt0,Tn0,Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz,nx,nz,V)

  ndof = size(M,2)

! if fault is on a domain boundary and symmetry is required, 
! bc2 and invM2 are not defined
  two_sides = tags(2)>0

! bc1 --> grid%bounds(i) boundary corresponding to TAG
  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%bc1 )
  if ( two_sides ) then
    call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%bc2 )

    if (bc%bc1%nelem/=bc%bc2%nelem) &
     call IO_abort('bc_dynflt_init: number of boundary elements do not match')

    if (bc%bc1%npoin/=bc%bc2%npoin) &
     call IO_abort('bc_dynflt_init: number of nodes on boundaries do not match')
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
      call IO_abort('bc_dynflt_init: coordinates on boundaries do not match properly')
  endif

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

  dt = TIME_getTimeStep(time)

  if (associated(bc%swf)) then
    call swf_init(bc%swf,bc%coord,dt)
  elseif (associated(bc%rsf)) then
    call rsf_init(bc%rsf,bc%coord,dt)
  endif

  allocate(bc%T(npoin,2))
  allocate(bc%D(npoin,ndof))
  allocate(bc%V(npoin,ndof))
  bc%T = 0d0
  bc%D = 0d0
  !bc%V = 0d0  ! DEVEL: RSF can't have zero slip velocity !
  call DIST_CD_Init(bc%input%V,bc%coord,V)
  bc%V(:,1)=V
  deallocate(V)

! Initial stress
  allocate(bc%T0(npoin,2))
  call DIST_CD_Init(bc%input%T,bc%coord,Tt0)
  call DIST_CD_Init(bc%input%N,bc%coord,Tn0)
  call DIST_CD_Init(bc%input%Sxx,bc%coord,Sxx)
  call DIST_CD_Init(bc%input%Sxy,bc%coord,Sxy)
  call DIST_CD_Init(bc%input%Sxz,bc%coord,Sxz)
  call DIST_CD_Init(bc%input%Syz,bc%coord,Syz)
  call DIST_CD_Init(bc%input%Szz,bc%coord,Szz)
  nx => bc%n1(:,1)
  nz => bc%n1(:,2)
  allocate(Tx(npoin))
  allocate(Ty(npoin))
  allocate(Tz(npoin))
  Tx = Sxx*nx + Sxz*nz
  Ty = Sxy*nx + Syz*nz
  Tz = Sxz*nx + Szz*nz
  if (ndof==1) then ! SH
    bc%T0(:,1) = Tt0 + Ty
  else ! P-SV
    bc%T0(:,1) = Tt0 + Tx*nz - Tz*nx
  endif 
  bc%T0(:,2) = Tn0 + Tx*nx + Tz*nz
  deallocate(Tt0,Tn0,Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz)

! Initial friction
  allocate(bc%MU(npoin))
  if (associated(bc%swf)) then
    bc%MU = swf_mu(bc%swf)
    if (associated(bc%twf)) bc%MU = min( bc%MU,twf_mu(bc%twf,bc%coord,0d0,bc%D(:,1)) )
  elseif (associated(bc%rsf)) then
    bc%MU = rsf_mu(bc%V(:,1),bc%rsf)
    if (associated(bc%twf)) bc%MU = min( bc%MU,twf_mu(bc%twf,bc%coord,0d0,bc%D(:,1)) )
  elseif (associated(bc%twf)) then
    bc%MU = twf_mu(bc%twf,bc%coord,0d0,bc%D(:,1))
  endif
 
! cohesion
  call DIST_CD_Init(bc%input%cohesion,bc%coord,bc%cohesion)
  if (any(bc%cohesion < 0d0)) call IO_abort('bc_dynflt_init: cohesion must be positive')

  ! Normal stress response
  call normal_init(bc%normal,dt,bc%T0(:,2))

!-- Open output files -----
! Set output parameters

  bc%oix1 = max(bc%oix1, 1)
  bc%oixn = min(bc%oixn, npoin)

! NOTE: file names limited to tags(1)<100
  write(oname,'("Flt",I2.2,"_sem2d.dat")') tags(1)
  bc%ounit = IO_new_unit()
  !DEVEL: use direct access, and recl
  !inquire( IOLENGTH=iol ) real( bc%D(bc%oix1:bc%oixn:bc%oixd,1) )
  !open(bc%ounit,file=oname,status='replace',access='direct',recl=iol)
  open(bc%ounit,file=oname,status='replace',form='unformatted')

  write(oname,'("Flt",I2.2,"_potency_sem2d.tab")') tags(1)
  bc%ou_pot = IO_new_unit()
  open(bc%ou_pot,file=oname,status='replace')

  write(oname,'("Flt",I2.2,"_init_sem2d.tab")') tags(1)
  ounit = IO_new_unit()
  open(ounit,file=oname,status='replace')
  do i=bc%oix1,bc%oixn,bc%oixd
    write(ounit,*) bc%T0(i,1),bc%T0(i,2),bc%MU(i)
  enddo
  close(ounit)

 ! adjust output timestep to the nearest multiple of dt:
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
 !      Stress fields are relative to initial values and
 !      the first NDAT lines contain the initial conditions.
 !
 ! NOTE: when the binary file is written one line at a time
 !      a record tag is also written at the beginning and end of line,
 !      => number of columns in output file = nb of data columns +2
 !
  write(hname,'("Flt",I2.2,"_sem2d")') tags(1)
  NDAT = 5
  if (bc%osides) NDAT = NDAT + 4*ndof
  NSAMP = (TIME_getNbTimeSteps(time) -bc%oit)/bc%oitd +1
  hunit = IO_new_unit()

  onx = (bc%oixn-bc%oix1)/bc%oixd +1
  open(hunit,file=trim(hname)//'.hdr',status='replace')
  write(hunit,*) 'NPTS NDAT NSAMP DELT'
  write(hunit,*) onx,NDAT,NSAMP,bc%odt
  if (bc%osides) then 
    if (ndof==1) then
      write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D2t:V1t:V2t'
    else
      write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D1n:D2t:D2n:V1t:V1n:V2t:V2n'
    endif
  else
    write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction'
  endif
  write(hunit,*) 'XPTS ZPTS'
  do i=bc%oix1,bc%oixn,bc%oixd
    write(hunit,*) bc%coord(:,i)
  enddo
  close(hunit)

! Sample SU script subset|b2a
 ! NOTE: how to convert to a human-readable ASCII table
 !      the Slip at time 10*DELT (assuming first output at t=0), 
 !      using Seismic Unix's subset and b2a :
 !              subset <Flt05_sem2d.hdr n1=(NPTS+2) n2=NDAT n3=NSAMP \
 !              if1s=1 n1s=NPTS ix2s=0 ix3s=10 \
 !              | b2a n1=NPTS > slip_t10.tab
 !
  open(hunit,file=trim(hname)//'.csh',status='replace')
  write(hunit,*) '#!/bin/csh -f'
  write(hunit,*) '#USAGE: '//trim(hname)//'.csh node field'
  write(hunit,*) '@ ifs = $2 - 1'
  write(hunit,'(A,I0,A,I0,A,I0,A)') 'subset <'//trim(hname)//'.dat n1=',onx+2, &
    ' n2=',NDAT,' n3=',NSAMP,' ix1s=$1 ix2s=$ifs | b2a n1=1 > '//trim(hname)//'.$1.$ifs.tab'
  close(hunit)

  end subroutine BC_DYNFLT_init

!=====================================================================
! Computes the boundary term B*T, 
! with fault tractions T satisfying contact and friction conditions,
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
!
! NOTE: possible periodicity does not need to be enforced at this level
!       because it is assumed that MxA and D are already periodic
!
! WARNING: The opening conditions are incomplete
!          there can be a problem if the fault opens and closes again
!          There should be no update of the friction state variable if there is opening
!
! WARNING: Need compatibility conditions (additional constraints)
!          when the fault reaches the free surface
!

!! DEVEL: This subroutine needs significant modifications and refactoring
!! for rate-and-state friction: a two-loop solver for second order accuracy

  subroutine BC_DYNFLT_apply(bc,MxA,V,D,time)

  use stdio, only: IO_abort
  use constants, only : PI
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_dynflt_type), intent(inout) :: bc
  double precision, intent(in) :: V(:,:),D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%npoin) :: strength
  double precision, dimension(bc%npoin,2) :: T
  double precision, dimension(bc%npoin,size(V,2)) :: dD,dV,dA
  integer :: ndof

  ndof = size(MxA,2)

! get predicted values
  dD = get_jump(bc,D) ! dD_predictor
  dV = get_jump(bc,V) ! dV_predictor
  dA = get_weighted_jump(bc,MxA) ! dA_free

! T_stick
  T(:,1:ndof) = bc%Z * ( dV + bc%CoefA2V*dA )

! rotate to fault frame (tangent,normal) if P-SV
  if (ndof==2) then
    dD = rotate(bc,dD,1)
    dV = rotate(bc,dV,1)
    dA = rotate(bc,dA,1)
    T  = rotate(bc,T,1)      
  endif

! apply symmetry to normal stress when needed
  if (.not.associated(bc%bc2) .or. ndof==1) T(:,2)=0d0 

! add initial stress
  T = T + bc%T0

! Solve for normal stress (negative is compressive)
  ! Opening implies free stress
  if (bc%allow_opening) T(:,2) = min(T(:,2),0.d0) 
  ! Update normal stress variables
  call normal_update(bc%normal,T(:,2),dV(:,1)) 
 !DEVEL: maybe need here a second loop to obtain second order

 ! Update friction and shear stress
 !WARNING: during opening the friction state variable should not evolve

 !-- velocity and state dependent friction 
  if (associated(bc%rsf)) then
    if (time%kind=='quasi-static') then
      call rsf_qs_solver(bc%V(:,1),T(:,1), normal_getSigma(bc%normal), bc%rsf)
    else
      call rsf_solver(bc%V(:,1), T(:,1), normal_getSigma(bc%normal), bc%rsf, bc%Z(:,1))
    endif
    bc%MU = rsf_mu(bc%V(:,1), bc%rsf)
                                        
   !DEVEL combined with time-weakening
   !DEVEL WARNING: slip rate is updated later, but theta is not

   ! superimposed time-weakening
    if (associated(bc%twf)) bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,time%time,bc%D(:,1)) )

    strength = - bc%MU * normal_getSigma(bc%normal)
                                         
    T(:,1) = sign( strength, T(:,1))

  else
   !-- slip weakening
    if (associated(bc%swf)) then
     ! For time schemes that update displacements explicitly 
     ! (i.e. displacement at time n+1 is independent of velocity and acceleration at time n+1)
     ! update the slip-dependent friction coefficient using the updated slip.
     ! Otherwise, use the slip from the previous time step (one-timestep delay)
      if (bc%CoefA2D==0d0) then
        call swf_update_state(dD(:,1),dV(:,1),bc%swf)
      else
        call swf_set_state(bc%D(:,1), bc%swf)
      endif
      bc%MU = swf_mu(bc%swf)
  
     ! superimposed time-weakening
      if (associated(bc%twf)) bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,time%time,bc%D(:,1)) )

   !-- pure time-weakening
    elseif (associated(bc%twf)) then
      bc%MU = twf_mu(bc%twf,bc%coord,time%time,bc%D(:,1))
    endif

   ! Update strength
    strength = bc%cohesion - bc%MU * normal_getSigma(bc%normal)
                                         
   ! Solve for shear stress
    T(:,1) = sign( min(abs(T(:,1)),strength), T(:,1))
                                  
  endif

! Subtract initial stress
  T = T - bc%T0

! Save tractions
  bc%T = T

! Rotate tractions back to (x,z) frame 
  if (ndof==2) T = rotate(bc,T,-1)

! Add boundary term B*T to M*a
  MxA(bc%node1,:) = MxA(bc%node1,:) + bc%B*T(:,1:ndof)
  if (associated(bc%bc2)) MxA(bc%node2,:) = MxA(bc%node2,:) - bc%B*T(:,1:ndof)

! Update slip and slip rate, in fault frame
  dA = dA - bc%T(:,1:ndof)/(bc%Z*bc%CoefA2V)
  bc%D = dD + bc%CoefA2D*dA
  bc%V = dV + bc%CoefA2V*dA

  end subroutine BC_DYNFLT_apply

!---------------------------------------------------------------------
  function get_jump (bc,v) result(dv)

  type(bc_dynflt_type), intent(in) :: bc
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

  type(bc_dynflt_type), intent(in) :: bc
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

  type(bc_dynflt_type), intent(in) :: bc
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
! OUTPUT FORMAT: at each time, 5 data lines, one column per fault node
! 1: Slip
! 2: Slip rate
! 3: Shear stress 
! 4: Normal stress 
! 5: Friction coefficient
!
  subroutine BC_DYNFLT_write(bc,itime,d,v)

  type(bc_dynflt_type), intent(inout) :: bc
  integer, intent(in) :: itime
  double precision, dimension(:,:), intent(in) :: d,v

  write(bc%ou_pot,'(6D24.16)') BC_DYNFLT_potency(bc,d), BC_DYNFLT_potency(bc,v)

  if ( itime < bc%oit ) return

  !DEVEL: use direct access, and rec
  write(bc%ounit) real( bc%D(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%V(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,1) )
  write(bc%ounit) real( bc%T(bc%oix1:bc%oixn:bc%oixd,2) )
  write(bc%ounit) real( bc%MU(bc%oix1:bc%oixn:bc%oixd) )

  if (bc%osides) then
    call export_side(bc,get_side(bc,d,1))
    call export_side(bc,get_side(bc,d,2))
    call export_side(bc,get_side(bc,v,1))
    call export_side(bc,get_side(bc,v,2))
  endif

  bc%oit = bc%oit + bc%oitd

  end subroutine BC_DYNFLT_write

  !----------
  function get_side(bc,d,side) result(delta)

  type(bc_dynflt_type), intent(in) :: bc
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

  type(bc_dynflt_type), intent(in) :: bc
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
  ! Sets the displacement on the boundary to a specific value

  subroutine BC_DYNFLT_set(bc,field_in,input,field_out)
   
  use stdio, only: IO_abort
  
  type(bc_dynflt_type), intent(in) :: bc
  double precision, intent(in) :: field_in(:,:)
  double precision, intent(in) :: input
  double precision, dimension(:,:), intent(out) :: field_out
  
  field_out = field_in
  field_out(bc%node1,:) = input

  end subroutine BC_DYNFLT_set

!=====================================================================
  function BC_DYNFLT_potency(bc,d) result(p)

  type(bc_dynflt_type), intent(in) :: bc
  double precision, dimension(:,:), intent(in) :: d
  double precision, dimension(1+size(d,2)) :: p

  double precision :: delta(bc%npoin,size(d,2))

  delta = get_jump(bc,d)

  if (size(d,2)==2) then
   ! p11, p22, p12
    p(1) = sum( bc%n1(:,1)*delta(:,1)*bc%B(:,1) )
    p(2) = sum( bc%n1(:,2)*delta(:,2)*bc%B(:,1) )
    p(3) = 0.5d0*sum( (bc%n1(:,1)*delta(:,2)+bc%n1(:,2)*delta(:,1))*bc%B(:,1) )
  
  else
   ! p13, p23
    p(1) = 0.5d0*sum( bc%n1(:,1)*delta(:,1)*bc%B(:,1) )
    p(2) = 0.5d0*sum( bc%n1(:,2)*delta(:,1)*bc%B(:,1) )

  endif

  end function BC_DYNFLT_potency

!=====================================================================
  ! Adapts timestep for quasi-static solver

  subroutine BC_DYNFLT_timestep(time,bc,hcell)
 
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(inout) :: time
  type(bc_dynflt_type), intent(in) :: bc
  double precision, intent(in) :: hcell

  if (associated(bc%rsf)) then 
    call rsf_timestep(time,bc%rsf,bc%V(:,1),normal_getSigma(bc%normal),hcell)
  endif  

  end subroutine

end module bc_dynflt
  
