module bc_abso

! Absorbing boundary conditions P1 and P3 of
! R. Stacey, BSSA Vol.78, No.6, pp. 2089-2097, Dec 1988
!
! P1. Clayton and Engquist, first order accurate:
!   V_t = - Cs * U_t,n
!   V_n = - Cp * U_n,n
!
! P3. Stacey, second order accurate:
!   V_t = - Cs * U_t,n - (Cp-Cs) * U_n,t
!   V_n = - Cp * U_n,n - (Cp-Cs) * U_t,t
!
! Implemented as a boundary term (virtual work of tractions),
! involving only tangential derivatives :
!
! P1. T_t = - rho*Cs * V_t
!     T_n = - rho*Cp * V_n
!
! P3. T_t = - rho*Cs * V_t + rho*Cs*(2*Cs-Cp) * U_n,t
!     T_n = - rho*Cp * V_n - rho*Cs*(2*Cs-Cp) * U_t,t
!
! NOTE: this form of P3.2 comes from e.g. Casadei et al. (2002)
!     T_n = - (lambda+2*mu)/Cp * V_n + (lambda*Cs+2*mu*(Cs-Cp))/Cp * U_t,t
!   with lambda = mu/Cs^2 *(Cp^2-2*Cs^2) 
!   It corresponds to the first two terms in eq.11 of Kausel (1992)
!

  use spec_grid
  use memory_info
  use bnd_grid, only : bnd_grid_type
  use sources
  use time_evol, only : timescheme_type

  implicit none
  private
  
  type bc_abso_type
    private
    double precision, dimension(:,:,:), pointer :: K => null()
    double precision, dimension(:,:), pointer :: C => null(),Ht => null(),n => null(),coord => null()
    double precision, dimension(:), pointer :: B => null()
    type(bnd_grid_type), pointer :: topo => null()
    type(source_type), pointer :: wav=>null() ! incident wave
    logical :: stacey,periodic,let_wave,is_flat
  end type

  public :: BC_ABSO_type, BC_ABSO_read, BC_ABSO_init, BC_ABSO_apply
    
contains

!=====================================================================
!   
! BEGIN INPUT BLOCK
!
! NAME   : BC_ABSORB
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Absorbing boundary
! SYNTAX : &BC_ABSORB stacey, let_wave /
!
! ARG: stacey   [log] [F] Apply Stacey absorbing conditions for P-SV.
!                Higher order than Clayton-Engquist (the default).
! ARG: let_wave [log] [T] Allow incident waves across this boundary 
!                if mechanism='WAVE' in &SRC_DEF
!
! NOTE   : Stacey conditions and incident waves are only implemented for 
!          vertical and horizontal boundaries
!
! END INPUT BLOCK

  subroutine BC_ABSO_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_abso_type), intent(out) :: bc
  integer, intent(in) :: iin

  logical :: stacey,let_wave
  character(20) :: abs_type

  NAMELIST / BC_ABSORB / stacey,let_wave

  stacey  = .false.
  let_wave = .true.

  read(iin,BC_ABSORB,END=100)
100 continue
  
  bc%stacey = stacey
  bc%let_wave = let_wave

  if (echo_input) then 
    if (stacey) then
      abs_type = 'Stacey'
    else
      abs_type = 'Clayton-Engquist'
    endif
    write(iout,200) abs_type,let_wave
  endif

  return

  200 format(5x, 'Type of absorbing boundary. . . .(stacey) = ',A &
            /5x, 'Allow incident wave . . . . . .(let_wave) = ',L1 )

  end subroutine BC_ABSO_read



!=====================================================================
!
!--- Definition of the working coefficients
!
  subroutine BC_ABSO_init(bc,tag,grid,mat,M,tim,src,perio)

  use prop_mat, only : matpro_elem_type, MAT_getProp
  use time_evol, only: timescheme_type, TIME_getCoefA2Vrhs
  use constants, only : NDIME, TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects
  use stdio, only: IO_abort

  type(bc_abso_type) , intent(inout) :: bc
  type(sem_grid_type), intent(in)    :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  integer            , intent(in)    :: tag
  type(timescheme_type), intent(in) :: tim
  double precision, intent(inout) :: M(:,:)
  type(source_type), pointer :: src(:)
  type(bc_periodic_type), pointer :: perio

  double precision :: CoefIntegr
  double precision, dimension(NDIME,NDIME) :: xjac
  double precision :: rho,c(2)
  integer :: i,j,k,e,bck,LocDimTan,GeoDimTan,GeoDimNor &
            ,ngll,ndof,bc_nelem,bc_npoin,ebulk,itab(grid%ngll),jtab(grid%ngll)
  logical :: is_wave

  !-- bc%topo => grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  bc%periodic = BC_PERIO_intersects(bc%topo,perio)

  ndof      = size(M,2)
  ngll      = grid%ngll
  bc_nelem  = bc%topo%nelem
  bc_npoin  = bc%topo%npoin

  allocate(bc%B(bc_npoin), bc%n(bc_npoin,2))
  call BC_get_normal_and_weights(bc%topo,grid,bc%n,bc%B, bc%periodic)

 ! check that the boundary is flat, vertical or horizontal:
 ! if horizontal:
  if (all( abs(bc%n(:,1))<TINY_XABS) ) then
    bc%is_flat = .true.
    GeoDimTan = 1 ! Coordinate component (x,z) tangent to domain side
    GeoDimNor = 2 ! Coordinate component (x,z) normal to domain side
 ! if vertical:
  elseif (all( abs(bc%n(:,2))<TINY_XABS) ) then
    bc%is_flat = .true.
    GeoDimTan = 2
    GeoDimNor = 1
 ! if not flat:
  else
    bc%is_flat = .false.
    GeoDimTan = 2
    GeoDimNor = 1
!    call IO_abort('BC_ABSO_init: boundary is not vertical or horizontal')
  endif
  
  allocate(bc%C(bc_npoin,ndof))
  call storearray('bc%C',bc_npoin*ndof,idouble)
  bc%C = 0d0

! NOTE: Stacey only implemented for flat boundaries, for inplane mode 
  bc%stacey = bc%stacey .and. (ndof==2) .and. bc%is_flat
  if (bc%stacey) then
    allocate(bc%K(ngll,ndof,bc_nelem))
    call storearray('bc%K',ngll*ndof*bc_nelem,idouble)
  endif

  do e =1,bc_nelem

    ebulk = bc%topo%elem(e)

    call SE_inquire(grid, edge=bc%topo%edge(e) &
                   ,itab=itab, jtab=jtab, dim_t=LocDimTan)
   ! LocDimTan = local dimension (xi,eta) tangent to the edge

    do k = 1,ngll

      i = itab(k)
      j = jtab(k)

      call MAT_getProp(rho,mat(ebulk),'rho',i,j)
      call MAT_getProp(c(GeoDimNor),mat(ebulk),'cp',i,j)
      call MAT_getProp(c(GeoDimTan),mat(ebulk),'cs',i,j)

      xjac  = SE_Jacobian(grid,ebulk,i,j)        ! xjac  = DGlobDLoc
      CoefIntegr = grid%wgll(k)*sqrt( xjac(1,LocDimTan)**2 + xjac(2,LocDimTan)**2 ) 

      bck = bc%topo%ibool(k,e)
      if (ndof==1) then
        bc%C(bck,1)  = bc%C(bck,1) + rho*c(GeoDimTan)*CoefIntegr
      else
        bc%C(bck,:)  = bc%C(bck,:) + rho*c*CoefIntegr
      endif

      if (bc%stacey) then
! P3. T_t = - rho*Cs * V_t + mu/Cs*(2*Cs-Cp) * U_n,t
!     T_n = - rho*Cp * V_n - mu/Cs*(2*Cs-Cp) * U_t,t
        xjac = SE_InverseJacobian(grid,ebulk,i,j) ! xjac = DLocDGlob
        bc%K(k,1,e)  = CoefIntegr * xjac(LocDimTan,GeoDimTan) &
          * rho*c(GeoDimTan)*(2d0*c(GeoDimTan)-c(GeoDimNor))
      endif

    enddo
  enddo

  if (bc%stacey) then
    bc%K(:,2,:)  = bc%K(:,1,:)
    bc%K(:,GeoDimTan,:) = - bc%K(:,GeoDimTan,:)
    bc%Ht => grid%hTprime
  endif

! NOTE: periodic only possible for flat boundaries
  if (bc%periodic .and. bc%is_flat) then
    bc%C(1,:) = bc%C(1,:) + bc%C(bc_npoin,:)
    bc%C(bc_npoin,:) = bc%C(1,:)
  endif

 ! Modify the mass matrix for implicit treatment of C*v :
 ! NOTE: Only for flat boundaries. 
 !       Non-flat boundaries are usually circular boundaries far away and combined with mesh coarsening,
!        so the elements are usually large and the explicit scheme is stable.
 !
 !   M*a_(n+1) = -K*d_pre_(n+alpha) -C*v_(n+alpha)
 ! with  v_(n+alpha) = v_pre_(n+alpha) + coefA2Vrhs*a_(n+1)
 ! and coefA2Vrhs depends on the time scheme
 !
 ! => (M+coefA2Vrhs*C)*a_(n+1) = -K*d_pre_(n+alpha) -C*v_pre_(n+alpha)
 !
  if (bc%is_flat) M(bc%topo%node,:) = M(bc%topo%node,:)  + TIME_getCoefA2Vrhs(tim)*bc%C

 ! search for a wave source, pick the first found
 ! NOTE: incident waves only implemented for flat boundaries
  if (.not. bc%is_flat) bc%let_wave=.false.
  if (bc%let_wave .and. associated(src)) then
    do i=1,size(src)
      call SO_inquire(src(i),is_wave=is_wave)
      if (is_wave ) then
        bc%wav => src(i)
        exit
      endif
    enddo
  endif
 ! if a wave source was found:
  if (associated(bc%wav)) then
    allocate(bc%coord(2,bc_npoin)) 
    bc%coord = grid%coord(:,bc%topo%node)
  else
    deallocate(bc%B)
    if (bc%is_flat) deallocate(bc%n)
  endif

  end subroutine BC_ABSO_init



!=====================================================================
!
! NOTE: for Stacey conditions we assume vertical or horizontal boundaries
!       otherwise C and K are 2x2 matrices
!
!  D => fields%displ_alpha
!  V => fields%veloc_alpha
!  MxA => fields%accel
!
! If there is an incident wave the total field = incident + diffracted
!     T = T_i + T_d
!     V = V_i + V_d
! and the absorbing conditions apply only to the diffracted field
!     T = T_i - C*V_d
!       = T_i - C*( V-V_i )

  subroutine BC_ABSO_apply(bc,D,V,MxA,time)

  type(bc_abso_type) , intent(in)    :: bc
  type(timescheme_type), intent(in) :: time
  double precision, intent(inout) :: MxA(:,:)
  double precision, intent(in) :: D(:,:),V(:,:)

  double precision, dimension(:,:), allocatable :: KxD, Vin, Tin, Vn
  integer, pointer :: nodes(:)
  integer :: e,i,k(bc%topo%ngnod),ndof

  nodes => bc%topo%node
  ndof=size(MxA,2)

  if (associated(bc%wav)) then
    allocate( Vin(bc%topo%npoin,ndof), Tin(bc%topo%npoin,ndof) )
    call SO_WAVE_get_VT( Vin, Tin, time%time,bc%coord,bc%n,bc%wav)
    do i=1,ndof
      MxA(nodes,i) = MxA(nodes,i) - bc%C(:,i)*( V(nodes,i) - Vin(:,i) ) +bc%B*Tin(:,i)
    enddo
    deallocate(Vin,Tin)
  else
    if (bc%is_flat .or. ndof==1) then
      MxA(nodes,:) = MxA(nodes,:) - bc%C*V(nodes,:)
    else
      allocate( Vn(bc%topo%npoin,ndof) )
      Vn(:,2) = ( V(nodes,1)*bc%n(:,1) + V(nodes,2)*bc%n(:,2) ) ! tmp dot product V.n
      Vn(:,1) = Vn(:,2)*bc%n(:,1)
      Vn(:,2) = Vn(:,2)*bc%n(:,2)
      MxA(nodes,1) = MxA(nodes,1) - bc%C(:,1)*Vn(:,1) - bc%C(:,2)*( V(nodes,1)-Vn(:,1) )
      MxA(nodes,2) = MxA(nodes,2) - bc%C(:,1)*Vn(:,2) - bc%C(:,2)*( V(nodes,2)-Vn(:,2) )
      deallocate(Vn)
    endif
  endif

  if (bc%stacey) then
    allocate( KxD(bc%topo%npoin,2) )
    KxD = 0d0
    do e=1,bc%topo%nelem
      k = bc%topo%ibool(:,e)
      KxD(k,:) = KxD(k,:) + bc%K(:,:,e) * matmul( bc%Ht, D(nodes(k),:) )
    enddo
    if (bc%periodic) then
      KxD(1,:) = KxD(1,:)+KxD(bc%topo%npoin,:)
      KxD(bc%topo%npoin,:) = KxD(1,:)
    endif 
    MxA(nodes,:) = MxA(nodes,:) - KxD 
    deallocate(KxD)
  endif

  end subroutine BC_ABSO_apply

end module bc_abso
