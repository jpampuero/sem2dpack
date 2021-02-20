module bc_dynflt
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
    double precision:: hnodeFix
    integer, dimension(:), pointer :: node1=>null(), node2=>null()
    double precision :: CoefA2V,CoefA2D
    double precision, dimension(:,:), pointer:: n1=>null(),B=>null(), &
      invM1=>null(),invM2=>null(),Z=>null(),T0=>null(),T=>null(),V=>null(),&
      D=>null(),coord=>null(), Tp=>null(), Dp=>null()
    double precision, dimension(:), pointer:: MU=>null(),cohesion=>null()
    double precision, dimension(:), pointer:: hnode=>null(), MuStar=>null() 
    type(swf_type), pointer :: swf => null()
    type(rsf_type), pointer :: rsf => null()
    type(twf_type), pointer :: twf => null()
    logical :: allow_opening
    type(normal_type) :: normal
    type(bnd_grid_type), pointer :: bc1 => null(), bc2 => null()
    type(bc_dynflt_input_type) :: input

   ! for outputs:
    double precision :: ot1, odt, odtD, odtS, ot
    integer :: oit,oitd,ounit,oix1,oixn,oixd,ou_pot,ou_time
    logical :: osides
  end type bc_dynflt_type

  public :: BC_DYNFLT_type, BC_DYNFLT_read, BC_DYNFLT_init,  & 
            BC_DYNFLT_apply, BC_DYNFLT_write, BC_DYNFLT_set, &
            BC_DYNFLT_select, BC_DYNFLT_timestep, BC_DYNFLT_trans, &
            BC_DYNFLT_set_array, BC_DYNFLT_update_disp,BC_DYNFLT_update_BCDV, &
            BC_DYNFLT_nnode, BC_DYNFLT_node, BC_DYNFLT_AppendDofFix, BC_DYNFLT_nDofFix,&
            BC_DYNFLT_reset

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT
! GROUP  : BOUNDARY_CONDITION, DYNAMIC_FAULT
! PURPOSE: Dynamic fault with friction
! SYNTAX : &BC_DYNFLT friction, cohesion|cohesionH, opening, Tn|TnH, Tt|TtH,
!                     Sxx|SxxH, Sxy|SxyH, Sxz|SxzH, Syz|SyzH, Szz|SzzH
!                     ot1, otd, oxi, osides , otdS, otdD, V/
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
! ARG: hnodeFix [dble] [0d0] user defined fixed hnode (adaptive time stepping)
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
! ARG: otdD     [dble] [0d0] Time lag between outputs for dynamic 
! ARG: otdS     [dble] [0d0] Time lag between outputs for static
! ARG: oxi      [int(3)] [(1,huge,1)] First, last node and stride for output.
!                The default resets oxi(2) = last fault node
!
! ARG: osides   [log] [F] Export displacement and velocities on each side
!                of the fault
! ARG: V        [dble] [1d-12] Initial slip velocity (needed for RSF)
!
!
! NOTE:
!     if otdD=0, otdS=0, then uniform time step is used
!     output steps are set using ot1 and otd, internally adjusted.
!
!     if adaptive time step is used, ot1, otdD, otdS are used to set 
!     output steps and are not internally adjusted as time step is
!     changing every step. (otdS may still be internally adjustable
!     if uniform step is used for dynamic case) 
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
  double precision :: otdD, otdS
  character(20) :: TtH,TnH, SxxH,SxyH,SxzH,SyzH,SzzH &
                  ,dt_txt,oxi2_txt, cohesionH, VH
  character(3) :: friction(2)
  integer :: i,oxi(3) 
  logical :: opening, osides
  double precision:: hnodeFix

  NAMELIST / BC_DYNFLT /  Tt,Tn,Sxx,Sxy,Sxz,Syz,Szz &
                         ,TtH,TnH,SxxH,SxyH,SxzH,SyzH,SzzH &
                         ,ot1,otd,oxi,osides, friction, opening &
                         ,cohesion, cohesionH, V, VH, &
                         otdD, otdS,hnodeFix

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

  V    = 1d-12
  VH   = ''

  ot1 = 0.d0
  otd = 0.d0
  otdD = 0.d0
  otdS = 0.d0

  oxi(1) = 1
  oxi(2) = huge(i)
  oxi(3) = 1
  osides = .false.

  friction(1) = 'SWF'
  friction(2) = ''

  cohesion = 0d0
  cohesionH = ''
  hnodeFix = 0d0
  
  opening = .true.

  read(iin,BC_DYNFLT,END=100)

  bc%ot1 = ot1
  bc%ot  = ot1
  bc%odt = otd
  bc%odtD = otdD
  bc%odtS = otdS

  bc%oix1 = oxi(1)
  bc%oixn = oxi(2) 
  bc%oixd = oxi(3) 
  bc%osides = osides
  bc%hnodeFix = hnodeFix

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
      if (otdD+otdS>0d0) then
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
    write(iout,200) TnH,TtH,SxxH,SxyH,SxzH,SyzH,SzzH,VH, & 
                    cohesionH,opening,ot1,dt_txt,&
                    otdD, otdS, oxi(1),oxi2_txt,oxi(3),osides 
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
            /5x,'Initial slip velocity . . . . . . . . .(V) = ',A,&
            /5x,'Cohesion  . . . . . . . . . . . (cohesion) = ',A,&
            /5x,'Allow opening . . . . . . . . . .(opening) = ',L1,&
            /5x,'Output first time . . . . . . . . . .(ot1) = ',EN16.8,&
            /5x,'       time step  . . . . . . . . . .(otd) = ',A,&
            /5x,'       time step in dynamic. . . . .(otdD) = ', EN16.8,&
            /5x,'       time step in static . . . . .(otdS) = ', EN16.8,&
            /5x,'       first node . . . . . . . . (oxi(1)) = ',I0,&
            /5x,'       last node  . . . . . . . . (oxi(2)) = ',A,&
            /5x,'       node stride  . . . . . . . (oxi(3)) = ',I0,&
            /5x,'       data from each fault side. (osides) = ',L1)

  end subroutine BC_DYNFLT_read

!=====================================================================
!
  subroutine BC_DYNFLT_init(bc,tags,grid,M,time,perio,vg,mat)

  use prop_mat, only : matpro_elem_type
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort,IO_new_unit
  use time_evol, only: timescheme_type, TIME_getTimeStep, TIME_getNbTimeSteps, &
                       TIME_getCoefA2V, TIME_getCoefA2D
  use constants, only : TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_dynflt_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  double precision, intent(in) :: M(:,:)
  double precision, intent(inout) :: vg(:,:)
  integer, intent(in) :: tags(2)
  type(timescheme_type), intent(in) :: time
  type(bc_periodic_type), pointer :: perio

  double precision, dimension(:), pointer :: Tt0,Tn0,Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz,nx,nz,V
  double precision, pointer :: tmp_n1(:,:), tmp_B(:), tmp_mustar(:)
  double precision :: dt
  integer :: i,j,k,hunit,npoin,NSAMP,NDAT,ndof,onx,ounit
  character(25) :: oname,hname
  logical :: two_sides, adapt_time
  double precision :: hnode_left, hnode_right
  double precision, allocatable:: theta(:) ! state variable
  double precision, allocatable:: rsf_a(:) ! rsf a parameter
  double precision, allocatable:: rsf_b(:) ! rsf b parameter
  
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

! compute temperary mustar
  allocate(tmp_mustar(bc%bc1%npoin))
  call BC_DYNFLT_Compute_MuStar(bc, grid, mat, tmp_mustar, ndof)

! nodes that are common to both sides (non-split nodes) are sticky nodes
! they must be deleted from the fault boundary
! remove the first and last nodes
  npoin = bc%bc1%npoin
  if ( two_sides ) then
    do k=1,bc%bc1%npoin
      if (bc%bc1%node(k)==bc%bc2%node(k))  npoin = npoin-1
    enddo
    if (npoin<bc%bc1%npoin) then
      allocate( bc%node1(npoin) )
      allocate( bc%node2(npoin) )
      allocate( bc%MuStar(npoin))
      j = 0
      do k=1,bc%bc1%npoin
        if (bc%bc1%node(k)/=bc%bc2%node(k)) then
          j=j+1
          bc%node1(j) = bc%bc1%node(k)
          bc%node2(j) = bc%bc2%node(k)
          bc%MuStar(j) = tmp_mustar(k)
        endif
      enddo
    deallocate(tmp_mustar)
    else
      bc%node1 => bc%bc1%node
      bc%node2 => bc%bc2%node
      bc%MuStar => tmp_mustar
    endif
  else
    bc%node1 => bc%bc1%node
    bc%MuStar => tmp_mustar
  endif
  bc%npoin = npoin

  allocate(bc%coord(2,npoin))
  bc%coord = grid%coord(:,bc%node1)
  if ( two_sides ) then
    if ( any(abs(bc%coord(1,:)-grid%coord(1,bc%node2))>TINY_XABS) &
      .OR.any(abs(bc%coord(2,:)-grid%coord(2,bc%node2))>TINY_XABS) )&
      call IO_abort('bc_dynflt_init: coordinates on boundaries do not match properly')
  endif

  ! compute node spacing for each node 
  ! allocate space for hnode

  allocate(bc%hnode(npoin))
  bc%hnode     = 0.0d0 

  if (bc%hnodeFix>TINY_XABS) then
      bc%hnode=bc%hnodeFix
  else 
      do i = 1, bc%npoin
         if (i==1) then
             hnode_left  = huge(0d0)
             hnode_right = sqrt(sum((bc%coord(:, i) - bc%coord(:,i+1))**2.0d0))
         elseif (i==bc%npoin) then
             hnode_left  = sqrt(sum((bc%coord(:, i-1) - bc%coord(:,i))**2.0d0))
             hnode_right = huge(0d0) 
         else
             hnode_left  = sqrt(sum((bc%coord(:, i-1) - bc%coord(:,i))**2.0d0))
             hnode_right = sqrt(sum((bc%coord(:, i) - bc%coord(:,i+1))**2.0d0))
         end if
         bc%hnode(i)   = min(hnode_left, hnode_right) 
      end do
  end if

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

  allocate(bc%T(npoin,2))
  allocate(bc%D(npoin,ndof))
  allocate(bc%Dp(npoin,ndof))
  allocate(bc%V(npoin,ndof))

  bc%T = 0d0
  bc%D = 0d0
  bc%Dp = 0d0
  bc%V = 0d0  ! initialize 
  call DIST_CD_Init(bc%input%V, bc%coord, V)

  bc%V(:,1) = V
  deallocate(V)

! Initial stress
  allocate(bc%T0(npoin,2))
  allocate(bc%Tp(npoin,2))
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
  
  bc%Tp = bc%T0

  dt =  time%dt 

  ! Normal stress response
  call normal_init(bc%normal, dt, bc%T0(:,2))

! Initial friction
  allocate(bc%MU(npoin))

  if (associated(bc%swf)) then
    call swf_init(bc%swf,bc%coord,dt)
    bc%MU = swf_mu(bc%swf)
    if (associated(bc%twf)) then 
        bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,0d0) )
    end if
  elseif (associated(bc%rsf)) then
    ! bc%MU = rsf_mu(bc%V(:,1),bc%rsf)
    ! Let mu be computed directly from initial stress
    bc%Mu = abs(bc%T0(:, 1)/bc%T0(:, 2))
    if (associated(bc%twf)) then
        bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,0d0) )
    end if
    call rsf_init(bc%rsf,bc%coord, dt, bc%Mu, abs(bc%V(:, 1)))

  elseif (associated(bc%twf)) then
    bc%MU = twf_mu(bc%twf,bc%coord,0d0)
  endif
 
! cohesion
  call DIST_CD_Init(bc%input%cohesion,bc%coord,bc%cohesion)
  if (any(bc%cohesion < 0d0)) call IO_abort('bc_dynflt_init: cohesion must be positive')


  ! update global velocity at the fault nodes
  ! rotate bc%V to x,z frame 
  bc%V  = rotate(bc, bc%V, -1)
  
  ! transform global velocity
  call BC_DYNFLT_trans(bc, vg, 1) 

  ! update fault slip velocity
  vg(bc%node1, :) = bc%V/2d0
  
  ! rotate back 
  bc%V  = rotate(bc, bc%V, 1)
  
  ! if rsf fault, subtract plate rate 
  if (associated(bc%rsf)) then
      vg(bc%node1, 1) = vg(bc%node1, 1) - rsf_vplate(bc%rsf)/2d0 
  end if
  
  call BC_DYNFLT_trans(bc, vg, -1) 

!-- Open output files -----
! Set output parameters

  bc%oix1 = max(bc%oix1, 1)
  bc%oixn = min(bc%oixn, npoin)

  adapt_time =  (time%kind=='adaptive') & 
                .and. (.not. time%fixdt) 

! NOTE: file names limited to tags(1)<100

  ! Binary file that stores data for outputing time steps
  write(oname,'("Flt",I2.2,"_sem2d.dat")') tags(1)
  bc%ounit = IO_new_unit()
  open(bc%ounit,file=oname,status='replace',form='unformatted')
  
  ! File to store the potency
  write(oname,'("Flt",I2.2,"_potency_sem2d.tab")') tags(1)
  bc%ou_pot = IO_new_unit()
  open(bc%ou_pot,file=oname,status='replace')
  
  ! File to store the output time information if adaptive time stepping
  ! write into binary files
  if (adapt_time) then
      write(oname,'("Flt",I2.2,"_time_sem2d.tab")') tags(1)
      bc%ou_time = IO_new_unit()
      open(bc%ou_time, file=oname, status='replace')
  end if

  ! write the initial stress and friction
  write(oname,'("Flt",I2.2,"_init_sem2d.tab")') tags(1)
  ounit = IO_new_unit()

  ! output initial parameters including state variable
  open(ounit,file=oname,status='replace')
  allocate(theta(bc%npoin))
  allocate(rsf_a(bc%npoin))
  allocate(rsf_b(bc%npoin))
  
  if (associated(bc%rsf)) then
      ! get state
      theta = rsf_get_theta(bc%rsf)
      rsf_a = rsf_get_a(bc%rsf)
      rsf_b = rsf_get_b(bc%rsf)
      do i=bc%oix1,bc%oixn,bc%oixd
        write(ounit,*) bc%T0(i,1),bc%T0(i,2),bc%MU(i), theta(i), bc%V(i, 1), &
                       rsf_a(i), rsf_b(i)                     
      enddo
  else
      do i=bc%oix1,bc%oixn,bc%oixd
        write(ounit,*) bc%T0(i,1),bc%T0(i,2),bc%MU(i)
      enddo
  end if
  close(ounit)

  deallocate(theta, rsf_a, rsf_b)

 ! output times 
 ! Only adjust times when using uniform time step
 ! Actual time should be read in _time_sem2d.tab
  if (.not. adapt_time) then
      bc%oitd = max(1,nint(bc%odt/dt))
      bc%odt = dt * dble(bc%oitd)
      bc%odt = dt * bc%oitd
      bc%oit = nint(bc%ot1/dt)
      bc%ot  = dt * dble(bc%oit) 
  end if

 ! HEADER FILE FORMAT:
 !      Name            FltXX_sem2d.hdr where XX=tags(1) of the BC_DYNFLT
 !                      input block, the tag of the first side of the fault.
 !      Line 1          Name of size/other variables
 !      Line 2          Value of size/other variables
 !      Line 3          Name of data fields, separated by ":"
 !      Line 4          Name of coordinate axis
 !      Line 5:end      Two column table of nodes coordinates
 !
 ! TIME FILE FORMAT:
 !      Name            FltXX_time_sem2d.tab where XX=tags(1) of the BC_DYNFLT
 !                      input block, the tag of the first side of the fault.
 !                      fields: 'it', 'dt', 'time', 'EqNum', 'isDynamic' 
 !      Line 1:end      7 columns: it, dt, time, EQNum, isDynamic, switch, isEQ
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

 !
 ! WARNING: 
 ! Header file is different if adaptive time stepping is used
 !
 ! NSAMP and DELT is not known before the end of simulation
 !
 ! And should be estimated from the FltXX_time_sem2d.tab file
 !
  
  NDAT = 5
  if (bc%osides) NDAT = NDAT + 4*ndof
  ! output state variable if rsf is associated
  if (associated(bc%rsf)) NDAT = NDAT + 1

  onx = (bc%oixn-bc%oix1)/bc%oixd +1

  write(hname,'("Flt",I2.2,"_sem2d")') tags(1)
  hunit = IO_new_unit()
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
    if (associated(bc%rsf)) then
        if (ndof==1) then
          write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D2t:V1t:V2t:theta'
        else
          write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D1n:D2t:D2n:V1t:V1n:V2t:V2n:theta'
        endif
    else
        if (ndof==1) then
          write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D2t:V1t:V2t'
        else
          write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:D1t:D1n:D2t:D2n:V1t:V1n:V2t:V2n'
        endif
    end if ! associated bc%rsf

  else
    ! if do not output sides
    if (associated(bc%rsf)) then
        write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction:theta'
    else
        write(hunit,'(A)') ' Slip:Slip_Rate:Shear_Stress:Normal_Stress:Friction'
    end if

  endif ! bc%osides

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
! 1: minus side
! 2: plus  side
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
!

subroutine BC_DYNFLT_apply(bc,MxA,V,D,time)
      ! MxA: force
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use constants, only : PI
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_dynflt_type), intent(inout) :: bc
  double precision, intent(inout) :: V(:,:)
  double precision, intent(in) :: D(:,:)
  double precision, intent(inout) :: MxA(:,:)
  if (.not. time%isDynamic) then
      call BC_DYNFLT_apply_quasi_static(bc,MxA,V,D,time)
  else
      call BC_DYNFLT_apply_dynamic(bc,MxA,V,D,time)
  end if

end subroutine BC_DYNFLT_apply

! reset the state parameter of rsf
subroutine BC_DYNFLT_reset(bc, dt)
  type(bc_dynflt_type), intent(inout) :: bc
  double precision :: dt
  call rsf_reset(bc%rsf, dt)
end subroutine BC_DYNFLT_reset

! ==============================================================================
! quasi-static kaneko 2011.
! Only works with classic rate-state friction
! 1. obtain fault traction
! 2. update state variable
! 3. update fault slip velocity
! 4. update fault velocity

subroutine BC_DYNFLT_apply_quasi_static(bc,MxA,V,D,time)
      ! MxA: force or -K*d
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use time_evol, only : timescheme_type

  type(timescheme_type), intent(in) :: time
  type(bc_dynflt_type), intent(inout) :: bc
  double precision, intent(inout) :: V(:,:)
  double precision, intent(in) :: D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%npoin, 2) :: T
  double precision, dimension(bc%npoin, size(V, 2)) :: dV, dF 
  integer :: ndof, i

  ndof = size(MxA,2)

  ! compute jump in force
  dF = get_jump(bc, MxA)     ! jump in force -Kd 
  dF = -dF                   ! correct sign, dF = [Kd]

  dV     = get_jump(bc, V)   ! slip velocity  

! compute fault traction
  do i =1,ndof
      ! Conforming mesh is assumed here
      T(:, i) = -dF(:, i)/(2d0*bc%B(:,1))
  enddo

! rotate to fault frame (tangent,normal) if P-SV
  if (ndof==2) then
    dF  = rotate(bc,dF,1)
    T   = rotate(bc,T,1)      
    dV  = rotate(bc,dV,1)      
  endif
  
  ! add plate rate to the fault slip
  if (associated(bc%rsf)) then
      dV(:, 1) = dV(:, 1) + rsf_vplate(bc%rsf)
  else
      call IO_abort('BC_DYNFLT_apply_quasi_static: only rsf is implemented!!')
  end if

! apply symmetry to normal stress when needed
  if (.not.associated(bc%bc2) .or. ndof==1) T(:,2)=0d0 
! add initial stress
  
  T(:,:) = T(:,:) + bc%T0(:,:)

! Solve for normal stress (negative is compressive)
  ! Opening implies free stress
  if (bc%allow_opening) T(:,2) = min(T(:,2),0.d0) 

  ! Update normal stress variables
  ! Note this regularization is negligible in quasi-static
  call normal_update(bc%normal, T(:,2), dV(:,1)) 

 !-- velocity and state dependent friction 
 ! compute new slip velocity

  if (associated(bc%rsf)) then
    call rsf_qs_solver(dV(:,1), T(:,1), normal_getSigma(bc%normal), bc%rsf)
  else
    call IO_abort('BC_DYNFLT_apply_quasi_static: only rsf is implemented!!') 
  endif

! save total stress as Tn
  bc%Tp = T

! Subtract initial stress
  T = T - bc%T0

! Subtract plate rate
  if (associated(bc%rsf)) then
    dV(:, 1) = dV(:, 1) - rsf_vplate(bc%rsf)
  else
    call IO_abort('BC_DYNFLT_apply_quasi_static: only rsf is implemented!!')
  end if
  
! Save tractions in fault local frame
! This traction is relative to T0
  bc%T = T

! Rotate tractions back to (x,z) frame 
  if (ndof==2) T  = rotate(bc, T , -1)
  if (ndof==2) dV = rotate(bc, dV, -1)

! Update fault velocity following Junpei Seki, 2017

! forward transform global velocity
  call BC_DYNFLT_trans(bc, V, 1)

! update half slip velocity by to the current prediction 
  V(bc%node1, :) = 0.5d0 * dV
  
! transform back, global velocity 
  call BC_DYNFLT_trans(bc, V, -1)

  bc%V = dV
  if (associated(bc%rsf)) bc%V(:,1) = bc%V(:,1) + rsf_vplate(bc%rsf)
  !bc%D = bc%D + bc%V * time%dt

end subroutine BC_DYNFLT_apply_quasi_static

subroutine BC_DYNFLT_update_BCDV(bc, D, time)
  type(bc_dynflt_type), intent(inout) :: bc
  double precision, intent(in) :: D(:,:)
  double precision :: time 

  bc%D  = get_jump(bc,D);
  bc%D  = rotate(bc,bc%D, 1)
  bc%Dp = bc%D
  if (associated(bc%rsf)) bc%D(:,1)=bc%D(:,1) + rsf_vplate(bc%rsf)*time 

end subroutine BC_DYNFLT_update_BCDV

subroutine BC_DYNFLT_apply_dynamic(bc,MxA,V,D,time)
      ! MxA: force
      ! V  : velocity
      ! D  : displacement

  use stdio, only: IO_abort
  use constants, only : PI
  use time_evol, only : timescheme_type, TIME_getCoefA2V, TIME_getCoefA2D

  type(timescheme_type), intent(in) :: time
  type(bc_dynflt_type), intent(inout) :: bc
  double precision, intent(in) :: V(:,:),D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%npoin) :: strength
  double precision, dimension(bc%npoin,2) :: T
  double precision, dimension(bc%npoin,size(V,2)) :: dD,dV,dA,dF
  integer :: ndof, i

  ndof = size(MxA,2)
  
  bc%CoefA2V = TIME_getCoefA2V(time)
  bc%CoefA2D = TIME_getCoefA2D(time)

! get predicted values
  dD = get_jump(bc,D) ! dD_predictor
  dV = get_jump(bc,V) ! dV_predictor
  dF = get_jump(bc,MxA) ! dF_predictor

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
! T is the stress from perturbation, must add background stress
! in dynamic stepping, background stress is relative to previous quasi-static
! solution Tn
  T = T + bc%Tp

! apply backslip
!  if (associated(bc%rsf)) dV(:,1) = dV(:,1) + rsf_vplate(bc%rsf)

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
    call rsf_solver(bc%V(:,1), T(:,1), normal_getSigma(bc%normal), bc%rsf, bc%Z(:,1), time)
    bc%MU = rsf_mu(bc%V(:,1), bc%rsf)

   !DEVEL combined with time-weakening
   !DEVEL WARNING: slip rate is updated later, but theta is not

   ! superimposed time-weakening
    if (associated(bc%twf)) bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,time%time) )

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
      if (associated(bc%twf)) bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,time%time) )

   !-- pure time-weakening
    elseif (associated(bc%twf)) then
      bc%MU = twf_mu(bc%twf,bc%coord,time%time)
    endif

   ! Update strength
    strength = bc%cohesion - bc%MU * normal_getSigma(bc%normal)

   ! Solve for shear stress
    T(:,1) = sign( min(abs(T(:,1)),strength), T(:,1))
                                  
  endif

! Subtract initial stress
  T = T - bc%Tp

! Save tractions
  bc%T = T

! Rotate tractions back to (x,z) frame 
  if (ndof==2) T = rotate(bc,T,-1)

! Add boundary term B*T to M*a
  MxA(bc%node1,:) = MxA(bc%node1,:) + bc%B*T(:,1:ndof)
  if (associated(bc%bc2)) MxA(bc%node2,:) = MxA(bc%node2,:) - bc%B*T(:,1:ndof)

  dA(:, 1:ndof) = dA(:, 1:ndof) - bc%T(:,1:ndof)/(bc%Z*bc%CoefA2V)
  
  ! note bc%V and bc%D store the total slip velocity and slip
  ! dD is perturbation
  bc%D = dD + bc%CoefA2D*dA
  bc%V = dV + bc%CoefA2V*dA

  if (associated(bc%rsf)) then
      bc%V(:,1)=bc%V(:,1) + rsf_vplate(bc%rsf)
      bc%D(:,1)=bc%D(:,1) + bc%Dp(:,1) + rsf_vplate(bc%rsf)*time%time
  end if

  ! eventually still save relative stress to T0
  bc%T = bc%T + bc%Tp - bc%T0 

  end subroutine BC_DYNFLT_apply_dynamic


function BC_DYNFLT_nDofFix(bc, ndof) result(n)
    type(bc_dynflt_type), intent(in) :: bc
    integer:: n, ndof
    n = size(bc%node1) * ndof
end function

! append index of dof to indexFixDof
subroutine BC_DYNFLT_AppendDofFix(bc, indexFixDof, istart, ndof)
  type(bc_dynflt_type), intent(in) :: bc
  integer, dimension(:), intent(inout) :: indexFixDof
  integer :: istart, ndof, i, iend, nnode_i

  nnode_i = size(bc%node1)
  do i = 1, ndof
      iend = istart + nnode_i - 1
      indexFixDof(istart: iend) = (bc%node1 - 1)*ndof + i
      istart = iend + 1
  end do

end subroutine BC_DYNFLT_AppendDofFix

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
  subroutine BC_DYNFLT_write(bc,time,d,v)
  use time_evol, only : timescheme_type
  type(bc_dynflt_type), intent(inout) :: bc
  type(timescheme_type), intent(in) :: time
  double precision, dimension(:,:), intent(in) :: d,v
  Logical :: adapt_time

  write(bc%ou_pot,'(6D24.16)') BC_DYNFLT_potency(bc,d), BC_DYNFLT_potency(bc,v)

  ! reset ot for each step that a switch in time scheme takes place
  if (time%switch) bc%ot = time%time

  if ( time%time < bc%ot ) return

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

  call export_theta(bc) ! export state variable

  ! update the next output time step index
  adapt_time = (time%kind=='adaptive') &
               .and. (.not.time%fixdt)
  if (.not. adapt_time) then
      bc%ot = bc%ot + bc%odt
  else
      ! adaptive time stepping 
      if (time%isdynamic) then
          bc%ot = bc%ot + bc%odtD 
      else
          bc%ot = bc%ot + bc%odtS 
      end if
      ! write time information if adaptive time into bindary
      write(bc%ou_time, *)  time%it, time%dt, time%time, & 
                        time%EQNum, time%isDynamic, time%switch, time%isEQ  
  end if
  end subroutine BC_DYNFLT_write

  subroutine export_theta(bc)
  ! export the state variable for rsf faults
  type(bc_dynflt_type), intent(in) :: bc
  double precision, dimension(bc%npoin) :: theta
  if (associated(bc%rsf)) then
      theta = rsf_get_theta(bc%rsf)
      write(bc%ounit) real( theta(bc%oix1:bc%oixn:bc%oixd) )
  end if

  end subroutine export_theta

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
! set values along boundary nodes with scalar input

subroutine BC_DYNFLT_set(bc, field, input, side_in)

  type(bc_dynflt_type), intent(in) :: bc
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

end subroutine BC_DYNFLT_set

! ======================================================
! set values along boundary nodes with array input
subroutine BC_DYNFLT_set_array(bc, field, input, side_in)

  type(bc_dynflt_type), intent(in) :: bc
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

end subroutine BC_DYNFLT_set_array
 
! select fields only on dynamic fault boundary and zero other components
subroutine BC_DYNFLT_select(bc, field_in, field_out, side_in)
   
  use stdio, only: IO_abort
  
  type(bc_dynflt_type), intent(in) :: bc
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
  
  end subroutine BC_DYNFLT_select

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
subroutine BC_DYNFLT_trans(bc, field, direction)
  ! note that field is directly modified
  type(bc_dynflt_type), intent(in) :: bc
  double precision, dimension(:,:), intent(inout) :: field
  double precision, dimension(size(bc%node1, 1),size(field, 2)) :: tmp1, tmp2
  integer, intent(in) :: direction
  integer :: i, Nbcnode, ndof

  ndof    = size(field, 2)
  Nbcnode = size(bc%node1, 1)

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

end subroutine BC_DYNFLT_trans

function BC_DYNFLT_nnode(bc) result(n)
  type(bc_dynflt_type), intent(in) :: bc
  integer :: n
  n = size(bc%node1, 1)
end function BC_DYNFLT_nnode

subroutine BC_DYNFLT_node(bc, node1, node2)
  type(bc_dynflt_type), intent(in) :: bc
  integer, dimension(size(bc%node1, 1)), intent(inout) :: node1, node2
  node1 = bc%node1 
  if (associated(bc%bc2)) then
      node2 = bc%node2
  else
      node2 = bc%node1
  end if
end subroutine BC_DYNFLT_node

! update displacement on the fault nodes 
! only used by quasi-static
subroutine BC_DYNFLT_update_disp(bc, dpre, d, v, dt)
  type(bc_dynflt_type), intent(in) :: bc  
  double precision, dimension(:,:), intent(in) :: dpre, v
  double precision, dimension(:,:) :: d
  double precision, intent(in) :: dt 
  integer :: i
  
  d(bc%node1, :) = dpre(bc%node1, :) + dt * v(bc%node1, :)
  if (associated(bc%bc2)) then
	  ! two sided fault
	  d(bc%node2, :) = dpre(bc%node2, :) + dt * v(bc%node2, :) 
  end if

end subroutine BC_DYNFLT_update_disp

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

  subroutine BC_DYNFLT_timestep(time, bc)
 
    use time_evol, only : timescheme_type
    type(timescheme_type), intent(inout) :: time
    type(bc_dynflt_type), intent(in) :: bc
  
    if (associated(bc%rsf)) then 
      call rsf_timestep(time, bc%rsf, bc%V(:,1), &
  				     -normal_getSigma(bc%normal), bc%hnode, bc%MuStar)
    endif 

  end subroutine
  
! compute the mu_star in Lapusta 2000 for adaptive time stepping 
! mu_star = mu for mode 3, ndof=1
! mu_star = mu/(1-nu) for mode 2, ndof=2
! mu is shear modulus, not confused with the frictional coefficient mu.
!
! For each fault node, compute its mu_star, maximum between two split nodes
!
! if a node is at the boundary of two elements, take either one 
! A more rigorious way is to take the larger one among the two elements.
!
  subroutine BC_DYNFLT_Compute_MuStar(bc, grid, mat, MuStar, ndof)
	use prop_mat, only  : matpro_elem_type, MAT_getProp
	use time_evol, only : timescheme_type 
	use constants, only : TINY_XABS
	use stdio, only: IO_abort
    use spec_grid, only : sem_grid_type, SE_inquire

	type(bc_dynflt_type) , intent(inout) :: bc
	type(sem_grid_type), intent(in) :: grid
	type(matpro_elem_type), intent(in) :: mat(:)
	double precision, intent(inout) :: MuStar(:)
    integer, intent(in) :: ndof

    double precision :: MuStar2
    double precision :: rho, cp, cs, nu
    integer :: i,j,k,e,bck,ngll,bc_nelem,&
               bc_npoin,ebulk,itab(grid%ngll),jtab(grid%ngll)
    
    MuStar    = 0.0d0
	ngll      = grid%ngll
	bc_nelem  = bc%bc1%nelem
	bc_npoin  = bc%bc1%npoin

    ! bc1 boundary grid 1
    do e=1, bc_nelem
	  ebulk = bc%bc1%elem(e)
	  call SE_inquire(grid, edge=bc%bc1%edge(e) &
                     ,itab=itab, jtab=jtab)
	  
	  do k = 1,ngll

	    i = itab(k)
	    j = jtab(k)

	    call MAT_getProp(rho,mat(ebulk),'rho',i,j)
	    call MAT_getProp(cp,mat(ebulk),'cp',i,j)
	    call MAT_getProp(cs,mat(ebulk),'cs',i,j)
	    
	    ! Poisson ratio 
	    nu = 0.5d0*((cp/cs)**2d0-2d0)/((cp/cs)**2d0-1d0)		  
	    !mu star1
	    bck = bc%bc1%ibool(k,e)
	    MuStar(bck) = cs**2*rho
        if (ndof==2) MuStar(bck)=MuStar(bck)/(1d0-nu)
	  enddo
    end do
     
     !bc 2
	if (associated(bc%bc2)) then
      do e=1, bc_nelem
		ebulk = bc%bc2%elem(e)
		call SE_inquire(grid, edge=bc%bc2%edge(e) &
                       ,itab=itab, jtab=jtab)
		
		do k = 1,ngll

		  i = itab(k)
		  j = jtab(k)

		  call MAT_getProp(rho,mat(ebulk),'rho',i,j)
		  call MAT_getProp(cp,mat(ebulk),'cp',i,j)
		  call MAT_getProp(cs,mat(ebulk),'cs',i,j)
		  
		  ! Poisson ratio 
		  nu = 0.5d0*((cp/cs)**2d0-2d0)/((cp/cs)**2d0-1d0)		  
	      !mu star1
		  bck = bc%bc1%ibool(k,e)
		  MuStar2 = cs**2*rho
          if (ndof==2) MuStar2=MuStar2/(1d0-nu)
		  MuStar(bck) = max(MuStar(bck), MuStar2)
		enddo
      end do
	 end if
  end subroutine BC_DYNFLT_Compute_MuStar

end module bc_dynflt
  
