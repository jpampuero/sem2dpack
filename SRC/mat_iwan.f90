module mat_iwan

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private


  type matwrk_iwan_type
    private
    integer :: Nspr
    double precision :: mu, lambda, gref, Kmod, ni, K0
    double precision, pointer, dimension(:,:,:) :: eps => null() 
    double precision, pointer, dimension(:,:,:) :: sig => null() 
    double precision, pointer, dimension(:,:,:) :: Rgll => null() 
    double precision, pointer, dimension(:,:,:) :: CNinvgll => null()     
    double precision, pointer, dimension(:) :: R => null()
    double precision, pointer, dimension(:) :: CNinv => null()
    double precision, pointer, dimension(:,:,:) :: S => null() 
    double precision, pointer, dimension(:,:,:,:) :: Sa => null() 
    double precision, pointer, dimension(:,:,:) :: F => null()
    integer, pointer, dimension(:,:) :: aktif => null()  

    double precision, pointer, dimension(:,:) :: Ws => null()
    double precision, pointer, dimension(:,:) :: pS => null()
    double precision, pointer, dimension(:,:) :: pS0 => null()
    double precision, pointer, dimension(:,:) :: pSb => null()
    double precision, pointer, dimension(:,:) :: T => null()
    logical, pointer, dimension(:,:) :: drysoil => null()


    double precision, pointer, dimension(:,:) :: Peff0, Gm0
    double precision, pointer, dimension(:,:) :: Gact     !!!
  
    double precision, pointer, dimension(:) :: ews => null()
    double precision, pointer, dimension(:) :: ewp => null()
    double precision, pointer, dimension(:,:,:,:) :: zipt => null()



  end type matwrk_iwan_type


  integer, save :: isIwan = 0
  integer, save :: isIai  = 0
  integer, save :: isVEP  = 0
  integer, save :: isOverburden = 0


  ! for memory report
  integer, save:: MAT_IWAN_mempro = 0
  integer, save:: MAT_IWAN_memwrk = 0


  double precision, save:: MAT_IWAN_WaterTable = -9d13


  public :: matwrk_iwan_type &
          , MAT_isIwan, MAT_isOverburden, MAT_IWAN_read &
          , MAT_IWAN_init_elem_prop &
          , MAT_IWAN_init_elem_work, MAT_IWAN_stress, MAT_IWAN_initial_stress &
          , MAT_IWAN_mempro, MAT_IWAN_memwrk



contains



!=======================================================================
  logical function MAT_isIwan(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isIwan = MAT_isKind(m,isIwan)
  end function MAT_isIwan
!=======================================================================
  logical function MAT_isOverburden(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isOverburden = MAT_isKind(m,isOverburden)
  end function MAT_isOverburden
!=======================================================================
  logical function MAT_isIai(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isIai = MAT_isKind(m,isIai)
  end function MAT_isIai
!=======================================================================
  logical function MAT_isVEP(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isVEP = MAT_isKind(m,isVEP)
  end function MAT_isVEP
!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MAT_IWAN
! GROUP  : MATERIALS
! PURPOSE: Set material properties for soil nonlinearity 
!          after the MPII model of Iwan(1967)
!
! SYNTAX 1: Viscoelastoplastic/Elastoplastic model with gref given 
!           and without pore-pressure development, in other words,
!           Pressure-independent model (Material type 1)
!           &MAT_IWAN cp,cs,rho,Nspr,gref,VEPMOD /
!
! SYNTAX 2: Viscoelastoplastic/Elastoplastic model with no gref given 
!           and with/out pore pressure development 
!           Pressure-dependent model, total/effective stress analysis (Material types 2-3)
!           &MAT_IWAN cp,cs,rho,Nspr,phi_f,WT,K0,cohesion,IAIMOD,VEPMOD /
!
!
! ARG: cp       [dble][0d0]  P wave velocity
! ARG: cs       [dble][0d0]  S wave velocity
! ARG: rho      [dble][0d0]  density
! ARG: Nspr     [int] [0]    number of Iwan springs 
! ARG: gref     [dble][-1d0] reference strain
! ARG: phi_f    [dble][-1d0] angle of failure line
! ARG: WT       [dble][0d0]  water table depth
! ARG: K0       [dble][1d0]  coefficient of Earth at rest
! ARG: cohesion [dble][0d0]  cohesion
! ARG: IAIMOD   [log] [  ]   excess-pore pressure development 
! ARG: VEPMOD   [log] [  ]   visco-elastoplasticity (referred to Liu& Archuleta, 2006)
!
!
! NOTE: WT must be defined only once in the first material block of types 2/3, 
! otherwise it is overwritten.
!
!
! END INPUT BLOCK


  subroutine MAT_IWAN_read(input,iin)

  use echo, only : echo_input, iout
  use constants, only: PI

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: cp,cs,rho,Nspr,gref,phi_f,WT,K0,cohesion,dum
  logical  :: IAIMOD, VEPMOD

  NAMELIST / MAT_IWAN / cp,cs,rho,Nspr,gref,phi_f,WT,IAIMOD,K0,cohesion,VEPMOD

  call MAT_setKind(input,isIwan)

  cp       = 0d0
  cs       = 0d0
  rho      = 0d0
  Nspr     = 0
  gref     = -1d0
  phi_f    = -1d0
  dum      = huge(dum) 
  WT       = dum
  K0       = 1d0
  IAIMOD   = .False.
  VEPMOD   = .False.
  cohesion = 0d0

  

  read(iin, MAT_IWAN, END=100)


  ! Setting water table level globally
  if (WT .NE. dum) &
  MAT_IWAN_WaterTable = WT


  write(iout,200) cp,cs,rho,Nspr,IAIMOD,VEPMOD
  ! Material type 2-3
  if (gref < 0d0) then
    write(iout,300) MAT_IWAN_WaterTable,K0,cohesion,phi_f
    call MAT_setKind(input,isOverburden)     
  ! Material type 1 
  else
    write(iout,400) gref 
  endif

  ! convert to radians
  phi_f = phi_f* PI/ 180d0


  if (Nspr == 0) &
  call IO_abort('MAT_IWAN.MAT_IWAN_read: Spring number zero')

  if (gref < 0d0 .AND. dsin(phi_f) < 0d0)&
  call IO_abort('MAT_IWAN.MAT_IWAN_read: Specify positive gref or failure-line angle')

  if (gref > 0d0 .AND. dsin(phi_f) > 0d0)&
  call IO_abort('MAT_IWAN.MAT_IWAN_read: Either gref or failure-line angle must be specified')

  if (MAT_IWAN_WaterTable == dum .AND. dsin(phi_f) > 0d0)&
  call IO_abort('MAT_IWAN.MAT_IWAN_read: Specify water table depth')

  if (IAIMOD .AND. gref > 0d0)&
  call IO_abort('MAT_IWAN.MAT_IWAN_read: Effective stress analysis allowed only for pressure-dependent models')

  if (IAIMOD .AND. dsin(phi_f) < 0d0)&
  call IO_abort('MAT_IWAN.MAT_IWAN_read: failure-line angle required for IAI model')



  ! Reading Iai model parameters
  if (IAIMOD) &
  call MAT_IAI_READ(input,iin)
  
  if (VEPMOD) &
  call MAT_VEP_READ(input,iin)



  ! Setting properties
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'Nspr',Nspr)


  ! Material type 2-3
  if (gref < 0d0) then
    call MAT_setProp(input,'m1', dsin(phi_f))
    call MAT_setProp(input,'cos_phif', dcos(phi_f))
    call MAT_setprop(input,'K0', K0)
    call MAT_setprop(input,'cohesion', cohesion) 
  ! Material type 1
  else
    call MAT_setProp(input,'gref',gref)
  endif    




  return

  100 call IO_abort('MAT_IWAN_read: MAT_IWAN input block not found')
  
  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Number of Iwan mechanisms . . . . .(Nspr) =',EN12.3,/5x, &
    'Front saturation model. . . . . .(IAIMOD) =',L3/5x, &
    'Viscoelastoplasticity . . . . . .(VEPMOD) =',L3/5x)

  300   format(5x, &
    'Water table level. . . . . . . . . . (WT) =',EN12.3,/5x, &  
    'Coefficient of Earth at rest . . . . (K0) =',EN12.3,/5x, &   
    'Cohesion . . . . . . . . . . . (cohesion) =',EN12.3,/5x, &      
    'Angle of failure line. . . . . . (phi_f) =',EN12.3,/5x)

  400   format(5x, &
    'Reference strain. . . . . . . . . .(gref) =',EN12.3,/5x)

  end subroutine MAT_IWAN_read

!=======================================================================

subroutine MAT_IAI_READ(input,iin)

! BEGIN INPUT BLOCK
!
! NAME   : MAT_IAI
! GROUP  : MATERIALS
! PURPOSE: Set material properties for 
!          liquefaction front model (Iai et al., 1990)
!   
! SYNTAX : &MAT_IAI  phi_p,p1,p2,s1,w1 /
!
! ARG: phi_p    [dble][-1d0] 
! ARG: p1       [dble][-1d0] 
! ARG: p2       [dble][-1d0] 
! ARG: S1       [dble][-1d0] 
! ARG: w1       [dble][-1d0] 
!
!
! END INPUT BLOCK

  use echo, only : echo_input, iout
  use constants, only: PI

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: p1,p2,S1,w1,phi_p

  NAMELIST / MAT_IAI / p1,p2,S1,w1,phi_p

  call MAT_setKind(input,isIai)
  p1    = -1d0 
  p2    = -1d0
  S1    = -1d0 
  w1    = -1d0 
  phi_p = -1d0

  read(iin, MAT_IAI, END=110)
  write(iout,210) phi_p,p1,p2,S1,w1

  ! convert to radians
  phi_p = phi_p* PI/ 180d0

  if (dsin(phi_p)<0d0 .OR. p1<0d0 .OR. p2<0d0 .OR. S1<0d0 .OR. w1<0d0 ) &
  call IO_abort('MAT_IWAN.MAT_IAI_read: Missing Iai model parameters')

  call MAT_setProp(input,'m2', dsin(phi_p))
  call MAT_setProp(input,'p1', p1)
  call MAT_setProp(input,'p2', p2)
  call MAT_setProp(input,'S1', S1)
  call MAT_setProp(input,'w1', w1)
  
  return 


  110 call IO_abort('MAT_IWAN.MAT_IAI_READ:MAT_IAI input block not found')

  210   format(5x, &
    'Angle of phase-transf. line . . . (phi_p) =',EN12.3,/5x, &
    'p1 . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x, &
    'p2 . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x, &
    'S1 . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x, &
    'w1 . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x)

end subroutine MAT_IAI_READ
!=======================================================================

subroutine MAT_VEP_READ(input,iin)

! BEGIN INPUT BLOCK
!
! NAME   : MAT_VEP
! GROUP  : MATERIALS
! PURPOSE: Set viscoelasticity properties for 
!          (Liu and Archuleta, 2006) model to be used in
!          visco-elastoplasticity
!   
! SYNTAX : &MAT_VEP  Qp,Qs,fr /
!
! ARG: Qp       [dble][-1d0] P wave quality factor
! ARG: Qs       [dble][-1d0] S wave quality factor
! ARG: fr       [dble][-1d0] Reference frequency
!
!
! END INPUT BLOCK

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: Qp,Qs,fr

  NAMELIST / MAT_VEP / Qp,Qs,fr

  call MAT_setKind(input,isVEP)

  Qp = -1d0
  Qs = -1d0
  fr = -1d0

  read(iin, MAT_VEP, END=120)
  write(iout,220) Qp,Qs,fr

  if (Qp<0d0 .OR. Qs<0d0 .OR. fr<0d0 ) &
  call IO_abort('MAT_IWAN.MAT_VEP_read: Missing viscoelasticity parameters')

  call MAT_setProp(input,'Qp',Qp)
  call MAT_setProp(input,'Qs',Qs)
  call MAT_setProp(input,'fr',fr)

  return 

  120 call IO_abort('MAT_IWAN.MAT_VEP_READ:MAT_VEP input block not found')

  220   format(5x, &
    'Qp . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x, &
    'Qs . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x, &
    'fr . . . . . . . . . . . . . . . . . . . .=',EN12.3,/5x)

end subroutine MAT_VEP_READ
!=======================================================================

  subroutine MAT_IWAN_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:) 
  
  integer :: Nspr
  double precision :: cp,cs,rho,dNspr,gref,mu,ni,Kmod
  double precision :: Qp,Qs,fr,visMs,visMp


  call MAT_setProp(elem,'cp',ecoord,MAT_IWAN_mempro)!
  call MAT_setProp(elem,'cs',ecoord,MAT_IWAN_mempro)!
  call MAT_setProp(elem,'Nspr',ecoord,MAT_IWAN_mempro)! 

  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_getProp(dNspr,elem,'Nspr') !dNspr=dble(Nspr)
  Nspr=int(dNspr)                     !Attention 



  ! Material types 2-3
  if (MAT_isOverburden(elem)) then
    call MAT_setProp(elem,'m1',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'K0',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'cohesion',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'cos_phif',ecoord,MAT_IWAN_mempro)!  
  ! Material type 1
  else
    call MAT_setProp(elem,'gref',ecoord,MAT_IWAN_mempro)!
  endif


  ! ELASTOPLASTICITY
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_IWAN_mempro)
  call MAT_setProp(elem,'mu',rho*cs*cs,MAT_IWAN_mempro)
  


  ! VISCO-ELASTOPLASTICITY
  if (MAT_isVEP(elem)) then
    call MAT_setProp(elem,'Qp',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'Qs',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'fr',ecoord,MAT_IWAN_mempro)!  
  endif



  ! PORE-PRESSURE EFFECT
  if (MAT_isIai(elem)) then
    call MAT_setProp(elem,'m2',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'p1',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'p2',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'S1',ecoord,MAT_IWAN_mempro)!  
    call MAT_setProp(elem,'w1',ecoord,MAT_IWAN_mempro)!  
  endif


  ! Elastic properties
  call MAT_getProp(mu,elem,'mu')

  ni = (cp*cp)/(cs*cs)  
  ni = (ni-2d0)/(2d0*(ni-1d0))
  call MAT_setprop(elem,'ni',ni,MAT_IWAN_mempro)

  Kmod = 2d0*mu*(1d0+ni)          ! E (Young modulus) 
  Kmod = Kmod/(3d0*(1d0-2d0*ni))  ! K (Bulk modulus)

  call MAT_setprop(elem,'Kmod',Kmod,MAT_IWAN_mempro)

  end subroutine MAT_IWAN_init_elem_prop

!=======================================================================

  subroutine MAT_IWAN_init_elem_work(m,p,ngll,ndof,sigmid,siginit,e,grid)

  use spec_grid, only : sem_grid_type, SE_elem_coord


  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  integer, intent(in) :: ndof 
  double precision, optional, intent(in) :: sigmid
  double precision, optional, intent(in) :: siginit(ngll,ngll)
  integer, intent(in),optional :: e 
  type(sem_grid_type), intent(in),optional :: grid

  double precision :: dNspr

 
  call MAT_getProp(dNspr,p,'Nspr')
  m%Nspr=int(dNspr)
  call MAT_getProp(m%lambda,p,'lambda')
  call MAT_getProp(m%mu,p,'mu')
  call MAT_getProp(m%ni,p,'ni')
  call MAT_getProp(m%Kmod,p,'Kmod')

  !Initiate strain and stress matrices
  allocate(m%eps(ngll,ngll,6))
  allocate(m%sig(ngll,ngll,6))
  m%eps = 0d0
  m%sig = 0d0

  ! Iwan parameters
  allocate(m%S(ngll,ngll,6))
  allocate(m%F(ngll,ngll,m%Nspr))
  allocate(m%Sa(ngll,ngll,m%Nspr,6))
  allocate(m%aktif(ngll,ngll))
  allocate(m%Gm0  (ngll,ngll))
  allocate(m%Rgll(ngll,ngll,m%Nspr))
  allocate(m%CNinvgll(ngll,ngll,m%Nspr-1))
  allocate(m%Peff0(ngll,ngll))
  allocate(m%Gact(ngll,ngll))       
  m%aktif = -1
  m%S     = 1d-3
  m%Sa    = 0d0
  m%F     = 0d0  
  

  ! don't change this order *
  if (MAT_isVEP(p)) &
  call MAT_IWAN_init_elem_vepmod(m,p,ngll,ndof)


  ! don't change this order *
  select case (MAT_isOverburden(p))
    case (.True.) 
      call MAT_IWAN_init_elem_overburden(m,p,ngll,siginit,sigmid,grid,e)
    case (.False.)
      call MAT_IWAN_init_elem_no_overburden(m,p,ngll)
  end select


  ! Memory envanter 
  MAT_IWAN_memwrk = MAT_IWAN_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%eps) &
                   + size(m%sig) &
                   + size(m%S) &
                   + size(m%F) &
                   + size(m%Sa) &
                   + size(m%aktif) &
                   + size(m%Rgll)     &
                   + size(m%CNinvgll) &
                   + size(m%Peff0)    &
                   + size(m%Gm0)  &
                   + size(m%Gact)                      


  end subroutine MAT_IWAN_init_elem_work
!============================================================================
!

subroutine MAT_IWAN_init_elem_vepmod(m,p,ngll,ndof)

  use mat_visla, only: MAT_VISLA_module

  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  integer, intent(in) :: ndof 

  double precision :: Qp,Qs,fr,visMs,visMp


  allocate(m%ews(8))
  allocate(m%ewp(8))   
  allocate(m%zipt(ngll,ngll,ndof+1,8))
  m%zipt = 0d0

  call MAT_getProp(Qp,p,'Qp')
  call MAT_getProp(Qs,p,'Qs')
  call MAT_getProp(fr,p,'fr')

  call MAT_VISLA_module(Qs,fr,m%mu,visMs,m%ews)
  call MAT_VISLA_module(Qp,fr,(m%lambda+2d0*m%mu),visMp,m%ewp)

  m%mu   = visMs
  m%Kmod = 2d0*visMs*(1d0+m%ni)            ! E (Young modulus) 
  m%Kmod = m%Kmod/(3d0*(1d0-2d0*m%ni))     ! K (Bulk modulus)

  m%lambda = visMp- 2d0*visMs

  MAT_IWAN_memwrk = MAT_IWAN_memwrk &
                 + size( transfer(m, (/ 0d0 /) )) &
                 + size(m%ews)     &
                 + size(m%ewp)     &
                 + size(m%zipt)     


end subroutine MAT_IWAN_init_elem_vepmod
!============================================================================

subroutine MAT_IWAN_init_elem_overburden(m,p,ngll,siginit,sigmid,grid,e)

  use spec_grid, only : sem_grid_type, SE_elem_coord

  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  double precision, optional, intent(in) :: sigmid
  double precision, optional, intent(in) :: siginit(ngll,ngll)
  type(sem_grid_type), intent(in),optional :: grid
  integer, intent(in),optional :: e 

  double precision :: gref(ngll,ngll), m1, K0, cohesion, cos_phif
  double precision :: ecoord(2,ngll,ngll), WT
  integer :: i,j


  WT = MAT_IWAN_WaterTable

  call MAT_getProp(K0,p,'K0')
  call MAT_getProp(m1,p,'m1')
  call MAT_getProp(cohesion,p,'cohesion')
  call MAT_getProp(cos_phif,p,'cos_phif')
  call MAT_getProp(m%K0,p,'K0')

  ! Initial confining stress
  m%sig(:,:,1) = siginit* K0    ! xx
  m%sig(:,:,6) = siginit* K0    ! yy
  m%sig(:,:,2) = siginit        ! zz

 
  m%Peff0 = siginit* (1d0+ 2d0*K0)/3d0
  m%Gm0   = m%mu * (abs(m%Peff0/(sigmid* (1d0+ 2d0*K0)/3d0)))**0.5

  !
  m%Gact  = m%Gm0

  m%S(:,:,1) = max(1d-3, m%sig(:,:,1)-m%Peff0)
  m%S(:,:,6) = max(1d-3, m%sig(:,:,6)-m%Peff0)
  m%S(:,:,2) = max(1d-3, m%sig(:,:,2)-m%Peff0)


  do i=1,ngll
    do j=1,ngll
      gref(i,j) = (m1* m%Peff0(i,j)+cohesion*cos_phif)/ m%Gm0(i,j)
      call MAT_IWAN_backbone_elem(gref(i,j),m%Nspr,m%Gm0(i,j),m%Rgll(i,j,:),m%CNinvgll(i,j,:))
    enddo
  enddo

  ! Liquefaction front model (Model type 3)
  if (MAT_isIai(p)) then
    call MAT_IWAN_init_shear_work(m,p,ngll,e)

    allocate(m%drysoil(ngll,ngll))
    m%drysoil = .False.
    ecoord    = SE_elem_coord(grid,e)

    do i=1,ngll
    do j=1,ngll
        ! Soil above water-table level
        if (abs(WT) > abs(ecoord(2,i,j))) &
        m%drysoil(i,j) = .True. 
    enddo
    enddo

    MAT_IWAN_memwrk = MAT_IWAN_memwrk &
                     + size( transfer(m, (/ 0d0 /) )) &
                     + size(m%drysoil)     
  endif


end subroutine MAT_IWAN_init_elem_overburden

!============================================================================
!

subroutine MAT_IWAN_init_elem_no_overburden(m,p,ngll)

  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll

  double precision :: gref(ngll,ngll)
  integer :: i,j

  m%Gm0   = m%mu  
  m%Gact  = m%mu

  do i=1,ngll
    do j=1,ngll
      call MAT_getProp(gref(i,j),p,'gref') 
      call MAT_IWAN_backbone_elem(gref(i,j),m%Nspr,m%Gm0(i,j),m%Rgll(i,j,:),m%CNinvgll(i,j,:))
    enddo
  enddo

end subroutine MAT_IWAN_init_elem_no_overburden
!============================================================================
!

subroutine MAT_IWAN_init_shear_work(m,p,ngll,e)

  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  integer, intent(in),optional :: e   

  integer :: i,j
  double precision :: m1,S1,w1,p1,p2, Wn,r, m4,m3,m2,a,b,c,w0

  call MAT_getProp(m1,p,'m1')
  call MAT_getProp(m2,p,'m2')
  call MAT_getProp(S1,p,'S1')
  call MAT_getProp(w1,p,'w1')
  call MAT_getProp(p1,p,'p1')
  call MAT_getProp(p2,p,'p2')

  allocate(m%Ws  (ngll,ngll))
  allocate(m%pS  (ngll,ngll))
  allocate(m%pS0 (ngll,ngll))
  allocate(m%pSb (ngll,ngll))
  allocate(m%T(ngll,ngll))

  m%pS  = 1d0 
  m%pS0 = 1d0 


  do i=1,ngll
    do j=1,ngll

      m%Gact(i,j) =  m%Gm0(i,j)

      Wn = (m%Peff0(i,j)* m1)**2/ (2d0* m%Gm0(i,j))
      r  = (m%sig(i,j,2)- m%sig(i,j,1))/ 2d0    
      r  = r/ m%Peff0(i,j) 

      m4 = 0d0
      m3 = 0.67d0* m2


      if (r .le. m3  .OR.  .not. MAT_isOverburden(p)) then
        m%Ws(i,j)  = 0d0
        m%pS0(i,j) = 1d0
        m%pSb(i,j) = 0.4d0
      else
        m4 = 1d0- (m2-m3)/ m1
        a  = m4*m4*m1*m1- m2*m2- 2d0*m3*m3+ 2d0*m2*m3
        b  = 2d0*r*m3- 2d0*m1*m1*m4
        c  = m1*m1- r*r
        
        m%pS0(i,j) = (-b- sqrt(b*b- 4d0*a*c))/ (2d0*a) 

        if (m%pS0(i,j) .le. S1) &
        m%pS0(i,j) = S1 

        if (m%pS0(i,j) .ge. 0.4) then
          m%pSb(i,j) = 0.4d0
        else
          m%pSb(i,j) = m%pS0(i,j)
        endif

        if (m%pS0(i,j) .ge. 0.4) then
          w0 = w1* ((1d0-m%pS0(i,j))/ 0.6d0)**(1d0/p1)

        else if (m%pS0(i,j) .le. S1) then
          w0 = 1d0
        else
          w0 = w1* ((m%pS0(i,j)-S1)/ (0.4d0-S1))**(-1d0/p2)
        endif
        m%Ws(i,j) = w0* Wn

      endif
    enddo
  enddo


  ! Memory report
  MAT_IWAN_memwrk = MAT_IWAN_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%Ws) &
                   + size(m%pS) &
                   + size(m%pS0) &
                   + size(m%pSb) &
                   + size(m%T) 


  end subroutine MAT_IWAN_init_shear_work
!============================================================================
! Constitutive law

  subroutine MAT_IWAN_stress(ndof,m,ngll,dt,de,sigsys,p)

!  subroutine MAT_IWAN_stress(ndof,m,ngll,dt,de,sigeps,sature,sigsys,p)

  use mat_visla, only: MAT_VISLA_strain, MAT_VISLA_strain2

  integer, intent(in) :: ndof, ngll
  double precision, intent(in) :: dt
  double precision, intent(in) :: de(ngll,ngll,ndof+1) 
  type (matwrk_iwan_type), intent(inout) :: m
!  double precision, intent(out) :: sigeps(ngll,ngll,2)
!  double precision, intent(out) :: sature(ngll,ngll,3)
  double precision, intent(out) :: sigsys(ngll,ngll,ndof+1)
  type(matpro_elem_type), intent(in) :: p

  double precision :: lambda
  integer :: i,j, index
  double precision :: deps(ngll,ngll,6), dsig(ngll,ngll,6),dplastic(ngll,ngll,6)
  double precision :: mu,Kmod,Gact,sig1,sig2,coeff, lambda2mu, ni_viscous


! Strain increment
  deps = 0d0
  index = -ndof*ndof+ 4
  do i=1,ndof+1
    deps(:,:,i+index) = de(:,:,i)
  enddo


! Viscoelasticity applies to strain increment        
  if (MAT_isVEP(p)) then
    if (ndof==2) then  !P-SV waves
      do i=1,ngll
        do j=1,ngll
          mu         = m%mu           !visMs  
          lambda     = m%lambda
          lambda2mu  = lambda+ 2* mu  !visMp

          call MAT_VISLA_strain2(m%zipt(i,j,1,:),dt,m%ewp,m%ews,de(i,j,1),de(i,j,2),lambda2mu,mu,sig1) !xx
          call MAT_VISLA_strain2(m%zipt(i,j,2,:),dt,m%ewp,m%ews,de(i,j,2),de(i,j,1),lambda2mu,mu,sig2) !zz
          coeff = 1d0 / (lambda2mu *lambda2mu - lambda* lambda)
          deps(i,j,1) = coeff * (lambda2mu * sig1 - lambda * sig2)            !xx
          deps(i,j,2) = coeff * (lambda2mu * sig2 - lambda * sig1)            !zz
          deps(i,j,3) = MAT_VISLA_strain(m%zipt(i,j,3,:),dt,m%ews,de(i,j,3))  !xz
        enddo
      enddo
    else               !SH waves
      do i=1,ngll
        do j=1,ngll
          deps(i,j,4) = MAT_VISLA_strain(m%zipt(i,j,1,:),dt,m%ews,de(i,j,1))  !xy
          deps(i,j,5) = MAT_VISLA_strain(m%zipt(i,j,2,:),dt,m%ews,de(i,j,2))  !yz
        enddo
      enddo
    endif    
  endif


  m%eps = m%eps+ deps 
!  if (ndof==1)   sigeps(:,:,1) = m%eps(:,:,5)   !yz component to write out    ndof=1  
!  if (ndof==2)   sigeps(:,:,1) = m%eps(:,:,3)   !xz component to write out    ndof=2


  ! Stress increment !
  dsig = 0d0
  do i=1,ngll
    do j=1,ngll
      mu     = m%Gm0(i,j)   
      Gact   = m%Gact(i,j)

      Kmod   = 2d0* m%Gm0(i,j)* (1d0+m%ni)      ! E (Young modulus) 
      lambda = Kmod* m%ni/ (1d0+m%ni)/ (1d0- 2d0*m%ni)
      Kmod   = Kmod/(3d0*(1d0- 2d0*m%ni))       ! K (Bulk modulus)

      call MAT_IWAN_stress_inc(deps(i,j,:),dsig(i,j,:),m%aktif(i,j),m%S(i,j,:)&
                        ,m%Nspr,m%Sa(i,j,:,:),m%F(i,j,:), m%Rgll(i,j,:),m%CNinvgll(i,j,:) &
                        ,mu,lambda,Kmod,dplastic(i,j,:),Gact)        

    enddo
  enddo
  m%sig = m%sig+ dsig
!  if (ndof==1)   sigeps(:,:,2) = m%sig(:,:,5)   !yz component to write out    ndof=1  
!  if (ndof==2)   sigeps(:,:,2) = m%sig(:,:,3)   !xz component to write out    ndof=2

  index = -ndof*ndof+ 4
  do i=1,ndof+1
    sigsys(:,:,i) = m%sig(:,:,i+index)
  enddo      


  if (MAT_isOverburden(p)  .and. ndof==2 ) then
      sigsys(:,:,1) = m%sig(:,:,1)- m%Peff0* m%K0* 3d0/(1d0+ 2d0* m%K0)   !xx
      sigsys(:,:,2) = m%sig(:,:,2)- m%Peff0* 3d0/(1d0+ 2d0* m%K0)         !zz  
  endif


  if ( MAT_isIai(p) ) then
    call MAT_IWAN_shear_work(m,p,ngll,dplastic)

!    ! Write out Iai-model-related changes
!    sature(:,:,1) = m%Gact(:,:)                 ! Current G modulus
!    sature(:,:,2) = m%T/ m%Peff0                ! Current normalized-deviatoric-stress
!    sature(:,:,3) = m%pS                        ! Current normalized-effective-stress
  endif


  end subroutine MAT_IWAN_stress

!============================================================================
!

  subroutine MAT_IWAN_backbone_elem(gref,Nspr,mu,R,CNinv)

  integer, intent(in) :: Nspr
  double precision, intent(in) :: gref,mu
  !double precision, intent(out) :: R(Nspr), CNinv(Nspr-1)
  double precision, intent(inout) :: R(:), CNinv(:)

  double precision :: x0,xu,dx, gama(Nspr),G(Nspr),summy
  integer :: i

  x0 = -6d0
  xu = log10(0.1d0)
  dx = (xu-x0)/(Nspr-1)

  gama  = 0d0
  G     = 0d0

  R     = 0d0
  CNinv = 0d0

  do i=1,Nspr
    gama(i) = 10d0**(x0+dx*(i-1))
    G(i)    = 1d0/ (1d0+ abs(gama(i)/gref))
    R(i)    = G(i)* mu* gama(i)
  enddo
! 

  summy = 0d0
  do i=1,Nspr-1
    CNinv(i) = (gama(i+1)/2d0- gama(i)/2d0)/(R(i+1)-R(i))  &
                - 0.5d0/mu - summy
    summy = summy+ CNinv(i)
  enddo   

  end subroutine MAT_IWAN_backbone_elem
 
!============================================================================

  subroutine MAT_IWAN_stress_inc(deps,dsig,aktif,S,Nspr,Sa,F,R,CNinv,mu,lambda,Kmod,dplastic,Gact)

  use lu

  double precision, intent(in)  :: deps(:)
  double precision, intent(out) :: dsig(:)
  integer, intent(inout) :: aktif
  double precision, intent(inout)  :: S(:)
  integer, intent(in) :: Nspr
  double precision, intent(inout)  :: Sa(:,:) 
  double precision, intent(inout)  :: F(:)
  double precision, intent(in)  :: R(:)
  double precision, intent(inout)  :: CNinv(:)  
  double precision, intent(in)  :: mu,lambda,Kmod
  double precision, intent(inout)  :: dplastic(:)
  double precision, intent(in)  :: Gact

  double precision ::  dS(6), dF(Nspr), Ed(6,6), Esd(6,6), de(6)
  double precision :: depsm, dsigm
  integer :: surface, errorflag, start,k,j
 
  integer :: D, INDX(6)


!   ! Elasticity test !
!   ! By commenting everything below the first three lines
!   ! in order to test the routine numerically
!   dsig  = 0d0
!   dsigm = (deps(1)+deps(2)+deps(6))* Kmod
!   dsig = MAT_IWAN_elastic(mu,Gact,lambda,deps)


!
  ! Nonlinearity
  de    = 0d0
  dsig  = 0d0 
  dsigm = (deps(1)+deps(2)+deps(6))* Kmod

  ! First time step
  if (aktif == -1) then
    dsig = MAT_IWAN_elastic(mu,Gact,lambda,deps)
    dS    = dsig
    dS(1) = dsig(1)- dsigm
    dS(2) = dsig(2)- dsigm
    dS(6) = dsig(6)- dsigm   
  
    aktif = 0
    S     = S+ dS
    Sa    = 0d0 
    dplastic = 0d0

    return
  endif


  depsm = (deps(1)+deps(2)+deps(6))/3d0
  ! Incremental deviatoric strain
  de(1) = deps(1)- depsm
  de(2) = deps(2)- depsm
  de(3) = deps(3)
  de(4) = deps(4)
  de(5) = deps(5)
  de(6) = deps(6)- depsm


  if (aktif == 0) then
    F(1)  = MAT_IWAN_surface(S,Sa(1,:))
    dF(1) = MAT_IWAN_dsurface(S,Sa(1,:),de)
  endif

  surface = 0
  if (aktif .GT. 0) then
    do j =1,aktif     

      ! New centers
      do k=1,6
        Sa(j,k) = S(k)- R(j)/sqrt(F(j))* (S(k)-Sa(j,k))
      enddo

      F(j)  = MAT_IWAN_surface(S,Sa(j,:))
      dF(j) = MAT_IWAN_dsurface(S,Sa(j,:),de)

      if ( (dF(j) .GE. 0d0)  .AND. (F(j) .GE. R(j)**2) ) &
      surface = surface+ 1  
    enddo
  endif

  if ( (dF(1) .GE. 0d0)  .AND. (F(1) .LT. R(1)**2) ) then
    dsig = MAT_IWAN_elastic(mu,Gact,lambda,deps)
    dS    = dsig
    dS(1) = dsig(1)- dsigm
    dS(2) = dsig(2)- dsigm
    dS(6) = dsig(6)- dsigm   
  
    S = S+ dS
    dplastic = 0d0
    return
  endif

  ! Ed computation
  Ed    = 0d0
  start = 1

  Ed(1,1) = 0.5d0/Gact
  Ed(2,2) = 0.5d0/Gact
  Ed(3,3) = 0.5d0/Gact
  Ed(4,4) = 0.5d0/Gact
  Ed(5,5) = 0.5d0/Gact
  Ed(6,6) = 0.5d0/Gact
 
  if (surface .gt. 0) &
  call MAT_IWAN_Ematris(start,Nspr,Ed,CNinv,S,Sa,F,surface)  


  do j = surface+1,Nspr-1
      F(j)  = MAT_IWAN_surface(S,Sa(j,:))
      dF(j) = MAT_IWAN_dsurface(S,Sa(j,:),de)    

      if ( (dF(j) .GE. 0d0)  .AND. (F(j) .GE. R(j)**2) ) then
        surface = surface+ 1  
        start   = surface
        call MAT_IWAN_Ematris(start,Nspr,Ed,CNinv,S,Sa,F,surface)
      else
        EXIT
      endif
   enddo 

!   Esd = 0d0
!   call MAT_IWAN_FINDInv(Ed,Esd,6,errorflag)
!   dS  = MATMUL(Esd,de)               
  
  ! Alternative for inversion
  call LUDCMP(Ed, 6, INDX, D, errorflag)
  if (errorflag .ne. 0) stop 'not invertible matrix'
  call LUBKSB(Ed, 6, INDX, de) ! solve Ed·x = de (de is used as input/ouput)
  dS = de

  
  S = S+ dS
  dsig = dS
  dsig(1) = dsig(1)+ dsigm
  dsig(2) = dsig(2)+ dsigm
  dsig(6) = dsig(6)+ dsigm


  ! Plastic strain-increment
  dplastic    = 0d0
  dplastic(1) = deps(1)- dsig(1)/ (lambda+ 2d0*mu)           !xx
  dplastic(2) = deps(2)- dsig(2)/ (lambda+ 2d0*mu)           !zz
  dplastic(3) = 2d0*deps(3)- dsig(3)/ Gact                   !xz 
  dplastic(4) = 2d0*deps(4)- dsig(4)/ Gact                   !xy  
  dplastic(5) = 2d0*deps(5)- dsig(5)/ Gact                   !yz
  dplastic(6) = deps(6)- dsig(6)/ (lambda+ 2d0*mu)           !yy

  !
  aktif = max (1, surface)

  end subroutine MAT_IWAN_stress_inc
!=======================================================================
! Incremental elastic computation of 
! stress matrix for given strain matrix

  function MAT_IWAN_elastic(mu,Gact,lambda,deps) result(dsig)

  double precision, intent(in) :: mu, lambda, Gact
  double precision, intent(in) :: deps(:)
  double precision :: dsig(6)

  double precision :: depsvol

  depsvol = deps(1)+deps(2)+deps(6)

  dsig(1) = 2d0*mu*deps(1)+ lambda*depsvol
  dsig(2) = 2d0*mu*deps(2)+ lambda*depsvol
  dsig(3) = 2d0*Gact*deps(3)
  dsig(4) = 2d0*Gact*deps(4)
  dsig(5) = 2d0*Gact*deps(5)
  dsig(6) = 2d0*mu*deps(6)+ lambda*depsvol

  end function MAT_IWAN_elastic
!=======================================================================
! Iwan surface(s) computation 

  function MAT_IWAN_surface(S,Sa)  result(F)
  
  double precision, intent(in) :: S(:)
  double precision, intent(in) :: Sa(:)
  double precision :: F
 
  F = 0d0
  F = 0.5d0* ((S(1)-Sa(1))**2+ (S(2)-Sa(2))**2+ &
                  2d0* (S(3)-Sa(3))**2 + &
                  2d0* (S(4)-Sa(4))**2 + &
                  2d0* (S(5)-Sa(5))**2 &
                  + (S(6)-Sa(6))**2)


  end function MAT_IWAN_surface  
!=======================================================================
! Iwan surface(s) movement computation 

  function MAT_IWAN_dsurface(S,Sa,de)  result(dF)
  
  double precision, intent(in) :: S(:)
  double precision, intent(in) :: Sa(:)
  double precision, intent(in) :: de(:)
  double precision :: dF
 
  dF = 0d0
  dF = 0.5d0* ((S(1)-Sa(1))*de(1)+ (S(2)-Sa(2))*de(2)+ &
                  2d0* (S(3)-Sa(3))*de(3) + &
                  2d0* (S(4)-Sa(4))*de(4) + &
                  2d0* (S(5)-Sa(5))*de(5) &
                  + (S(6)-Sa(6))*de(6))

  end function MAT_IWAN_dsurface  
!=======================================================================
! Iwan plasticity matrix
!

  subroutine MAT_IWAN_Ematris (start,Nspr,Ed,CNinv,S1,Sa1,F1,aktif)

  integer, intent(IN) :: start
  integer, intent(IN) :: Nspr
  double precision, intent(inout)  :: Ed(:,:)
  double precision, INTENT(IN)     :: CNinv(:)
  double precision, INTENT(IN)     :: S1(:)
  double precision, INTENT(IN)     :: Sa1(:,:)
  double precision, INTENT(INOUT)  :: F1 (:)
  integer, intent(IN) :: aktif

  integer :: j,m,k
  double precision :: ss(6)

  ss(1) = 1d0
  ss(2) = 1d0
  ss(3) = 2d0
  ss(4) = 2d0
  ss(5) = 2d0
  ss(6) = 1d0

  j = start
  do while (j .lt. aktif+1)
    do m = 1,6
      do k = 1,6    
        Ed(m,k) = Ed(m,k)+ CNinv(j)* ss(k)* (S1(m)-Sa1(j,m))&
                  *(S1(k)-Sa1(j,k))/ (2.0* F1(j))
      enddo
    enddo
    j = j+1
  enddo

end subroutine MAT_IWAN_Ematris
!=======================================================================

subroutine MAT_IWAN_initial_stress(mat_elem,grid,ntags,sigmid,siginit)

  ! NOTE: While meshing, z coordinate is very important for this part
  ! The origin z=0 point should not be in the middle because everything
  ! below is calculated based on absolute values of depth and WT depth.
  ! it could be modified later !

  use spec_grid, only : sem_grid_type, SE_elem_coord

  type(matpro_elem_type), intent(in) :: mat_elem(:)
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: ntags
  double precision, intent(inout)   :: sigmid(ntags)
  double precision, intent(inout)   :: siginit(grid%ngll,grid%ngll,grid%nelem)

  integer :: e, i,j, tag, count, n
  double precision, dimension(:), allocatable :: zmin,sigmax,zmax
  double precision :: ecoord(2,grid%ngll,grid%ngll)
  double precision :: rho, minim, maxim, WT, ekle


  ! Water-table level
  WT = MAT_IWAN_WaterTable

  allocate(zmin(ntags),zmax(ntags),sigmax(ntags))
  ! Reference depth for each domain
  ! finding min and max depths of each layer
  do i=1,ntags
    count = 0
    do e=1,grid%nelem
      tag = grid%tag(e) 
      if (tag == i) then
        ecoord  = SE_elem_coord(grid,e)
        if (count == 0) then 
          zmin(i) = minval(abs(ecoord(2,:,:)))
          zmax(i) = maxval(abs(ecoord(2,:,:)))
        else
          zmin(i) = min(zmin(i), minval(abs(ecoord(2,:,:))))
          zmax(i) = max(zmax(i), maxval(abs(ecoord(2,:,:))))
        endif
        count = count+1
      endif
    enddo
  enddo


  ! ! Relative initial stresses
  ! do e=1,grid%nelem
  !   ecoord  = SE_elem_coord(grid,e)
  !   call MAT_getProp(rho,mat_elem(e),'rho')

  !   do i=1,grid%ngll
  !     do j=1,grid%ngll

  !       ! Water table level - bulk density correction
  !       if (abs(WT) <  abs(ecoord(2,i,j))) then
  !         siginit(i,j,e) = 9.8*(rho-1d3)* (abs(ecoord(2,i,j))-zmin(grid%tag(e)))
  !       else
  !         siginit(i,j,e) = 9.8*(rho)* (abs(ecoord(2,i,j))-zmin(grid%tag(e)))
  !       endif

  !       ! TO CHANGE LATER - ATTENTION !
  !       ! how to understand the surface level !
  !       ! Assuming surface coordinate (z) = 0d0
  !       if (abs(ecoord(2,i,j)) == 0d0) &
  !       siginit(i,j,e) = siginit(i,j-1,e)

  !     enddo
  !   enddo  
  ! enddo




! Quick and dirt modification for Nepal simulations
  ! Relative initial stresses
  do e=1,grid%nelem
    ecoord  = SE_elem_coord(grid,e)
    call MAT_getProp(rho,mat_elem(e),'rho')

    do i=1,grid%ngll
      do j=1,grid%ngll

        siginit(i,j,e) = 9.8*(rho)* (abs(ecoord(2,i,j))-zmin(grid%tag(e)))

        ! fixing minimum depth to 1 m.
        if (ecoord(2,i,j) > -0.5d0  ) &
        siginit(i,j,e) = 9.8* rho* 0.5d0

      enddo
    enddo  
  enddo
!


  ! IF WATER TABLE IS INSIDE THE DOMAIN
  do n=1,ntags
    if ( abs(WT) > zmin(n)  .AND.  abs(WT)< zmax(n) ) then
      
      ekle = 1d3* 9.8d0* (abs(WT)-zmin(n))

      do e=1,grid%nelem
      tag = grid%tag(e) 

      if (tag == n) then
        ecoord  = SE_elem_coord(grid,e)
        do i=1,grid%ngll
        do j=1,grid%ngll
          if ( abs(ecoord(2,i,j)) >  abs(WT) ) then
            siginit(i,j,e) = siginit(i,j,e)+ ekle           
          endif
        enddo
        enddo
      endif
      enddo
    endif
  enddo


  ! Max initial stress for each domain   !!! CORRECTED
  ! to be changed for irregular layer interfaces!!!
  sigmax = 0d0
  do i=1,ntags-1

    do e=1,grid%nelem
      if (grid%tag(e) == i) then
      sigmax(i+1) = max(sigmax(i+1), maxval(siginit(:,:,e)))
      endif
    enddo
    
    sigmax(i+1) = sigmax(i+1)+ sigmax(i)
  enddo

  do e=1,grid%nelem
    do i=1,grid%ngll
      do j=1,grid%ngll
        siginit(i,j,e) = siginit(i,j,e)+ sigmax(grid%tag(e))
      enddo
    enddo
  enddo

  
  ! Midlayer stress for each domain
  do i=1,ntags
    count = 0
    do e=1,grid%nelem

      if (grid%tag(e) == i) then
        if (count == 0) then
          minim = minval(siginit(:,:,e))
          maxim = maxval(siginit(:,:,e))
        else
          minim = min(minim, minval(siginit(:,:,e)))
          maxim = max(maxim, maxval(siginit(:,:,e)))
        endif
        count = count+ 1
      endif
    enddo

    if (i == 1) then 
      sigmid(i) = 0.5d0* (maxim)
    else
      sigmid(i) = 0.5d0* (minim+maxim)
    endif
  enddo

  deallocate(zmin,sigmax)
  return

end  subroutine MAT_IWAN_initial_stress
!=======================================================================
  
  subroutine MAT_IWAN_shear_work(m,p,ngll,dplastic)


  type(matwrk_iwan_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  double precision, intent(in) :: dplastic(ngll,ngll,6)

  integer :: i,j,k
  double precision :: m1,m2,S1,w1,p1,p2
  double precision :: sigmap(3),dWs, coef(6),Told,r,ssl,CoR,w,Wn
  double precision :: rr2,mm3,rr3,S2
  double precision :: gref,sigmax,delta,cohesion,cos_phif


  call MAT_getProp(m1,p,'m1')
  call MAT_getProp(m2,p,'m2')
  call MAT_getProp(S1,p,'S1')
  call MAT_getProp(w1,p,'w1')
  call MAT_getProp(p1,p,'p1')
  call MAT_getProp(p2,p,'p2')
  call MAT_getProp(cohesion,p,'cohesion')
  call MAT_getProp(cos_phif,p,'cos_phif')

  coef    = 0d0
  coef(1) = 1d0
  coef(2) = 1d0
  coef(6) = 1d0
 
  
  do i=1,ngll
    do j=1,ngll

      call MAT_IWAN_principal_stress(m%sig(i,j,:), sigmap)
      Told     = m%T(i,j)
      m%T(i,j) = (sigmap(1)- sigmap(3))/2d0

      ! Plastic shear work increment
      dWs  = 0d0
      do k=1,6
        dWs = dWs+ (m%sig(i,j,k)-coef(k)*(m%sig(i,j,1)+m%sig(i,j,2)+m%sig(i,j,6))/3d0)*&
                    dplastic(i,j,k)
    enddo

    ! Previous deviatoric ratio
    r  = Told/ m%Peff0(i,j)   

    ssl = 0d0
    ssl = 0.4d0+ (m%pSb(i,j)- 0.4d0)*m%pS0(i,j)/ m%pSb(i,j)

    ! Correction of shear work
    if (m%pS(i,j) .ge. ssl) then
        if (r/m%pS0(i,j) .le. 0.67d0*m2) then
          CoR = 1d0
        else
          CoR = (m1- r/m%pS(i,j))/ (m1- 0.67d0*m2)
        endif
    else
        if (r .le. ssl*0.67d0*m2) then
          CoR = 1d0
        else
          CoR = (ssl*m1- r)/ (ssl*(m1- 0.67d0*m2) )
        endif
    endif

    if (dWs  .gt. 0d0)&
    dWs = CoR* dWs

    !
    dWs = max(0d0, dWs)

    m%Ws(i,j) = m%Ws(i,j)+ dWs
      
    ! Normalized shear work
    Wn = (m%Peff0(i,j)* m1)**2/ (2d0* m%Gm0(i,j))
    w  = m%Ws(i,j)/ Wn

    ! Variable S0
    if (w   .le.  0d0)    then
      m%pS0(i,j) =   1d0
    elseif  (w   .le.  w1)   then
        m%pS0(i,j) =   1d0- 0.6d0* ((w/w1)**p1)
    elseif  (w   .gt.  w1)   then
        m%pS0(i,j) =   (0.4d0- S1)* ((w1/w)**p2)+ S1
    
    endif

    ! Actual deviatoric ratio
    r   =   m%T(i,j)/ m%Peff0(i,j)

    rr2     =   m2* m%pS0(i,j)
    mm3     =   m2* 0.67d0
    rr3     =   mm3* m%pS0(i,j)
    S2      =   m%pS0(i,j)- (rr2-rr3)/ m1

    ! VARIABLE S
    if (r .le. rr3) then
        m%pS(i,j)  = m%pS0(i,j)
    else
        m%pS(i,j)  = S2+ sqrt((m%pS0(i,j)- S2)**2+ ((r- rr3)/ m1)**2)
    endif

    ! Soil above water table level (dry soil)
    if (m%drysoil(i,j)) then
      m%pS0(i,j) = 1d0
      m%pS(i,j)  = 1d0
    endif
    

    gref = (m1* m%Peff0(i,j)+cohesion*cos_phif)/ m%Gm0(i,j)   

    ! [Iai et al. 1990 - Eqns 86-91]
    delta =  0d0 
    if (m%pS0(i,j) .gt.  m%pSb(i,j)) then
      sigmax      = m%Peff0(i,j)* m1* m%pS(i,j)
    else
      delta  = (m1-m2)* (m%pSb(i,j)-m%pS0(i,j))* (0.4d0/m%pSb(i,j))* m%Peff0(i,j)
      sigmax = m%Peff0(i,j)* m1* m%pS(i,j)+ delta
      gref   = gref/ (m%pS0(i,j)/m%pSb(i,j))
    endif

    ! Updated shear modulus
    m%Gact(i,j) = sigmax/ gref 

    call MAT_IWAN_backbone_elem(gref,m%Nspr,m%Gact(i,j),m%Rgll(i,j,:),m%CNinvgll(i,j,:))
    enddo
  enddo

  end subroutine MAT_IWAN_shear_work

!============================================================================
!

  subroutine MAT_IWAN_principal_stress(sigma,sigmap)

    double precision,    intent(IN)      :: sigma(:)
    double precision,    intent(INOUT)   :: sigmap(:)

    double precision :: y(3),TG,a,b,c,d,delta,x1,x2,x3,gamma
    integer  :: l,i,j

    if (sigma(3) == 0d0 .AND. sigma(4)== 0d0 .AND. sigma(5) == 0d0)   then
      y(1) = sigma(1)   !xx
      y(2) = sigma(6)   !yy
      y(3) = sigma(2)   !zz

      do i = 1,2     
        l = i
        do j = i,3
          if (y(l) .gt. y(j) )   l = j
        enddo

        if(l.NE.i) then
          TG   = y(i)
          y(i) = y(l)
          y(l) = TG
        end if
      end do 

      do i = 1,3
        sigmap(i) = y(4-i)
      end do
      return
    endif

    ! Coefficients of the eqn of 3rd power
    a = -1d0
    b = sigma(1)+ sigma(2)+ sigma(6)
    c = sigma(3)**2 + sigma(4)**2+ sigma(5)**2 &
        - sigma(1)*sigma(6)&
        - sigma(2)*sigma(6)&
        - sigma(2)*sigma(1)
    d = sigma(1)*sigma(2)*sigma(6) &
        + 2d0*sigma(3)*sigma(4)*sigma(5) &
        - sigma(1)*sigma(3)**2 &      
        - sigma(6)*sigma(5)**2 &      
        - sigma(2)*sigma(4)**2        

    delta = b**2- 3d0*a*c
    if (delta .ne. 0d0) then
        gamma = (9d0*a*b*c- 2d0*b**3- 27d0*(a**2)*d)/(2d0* sqrt(abs(delta**3)))
    endif

    if (delta > 0d0) then
      if (abs(gamma) .le. 1d0) then
        x1 = (2d0* sqrt(delta)*(cos(acos(gamma)/3d0))-b)/(3d0*a)
        x2 = (2d0* sqrt(delta)* cos(acos(gamma)/3d0 - 2d0* 4d0* atan(1d0)/3d0)-b)/(3d0*a)
        x3 = (2d0* sqrt(delta)* cos(acos(gamma)/3d0 + 2d0* 4d0* atan(1d0)/3d0)-b)/(3d0*a)
      else
        return
      endif
    else if (delta .ge. 0d0) then
      return
    endif

    y(1) = x1
    y(2) = x2
    y(3) = x3

    do i = 1,2     
      l = i
      do j = i,3
        if (y(l) .gt. y(j) )   l = j
      enddo

      if(l.NE.i) then
        TG   = y(i)
        y(i) = y(l)
        y(l) = TG
      end if
    end do 

    do i = 1,3
      sigmap(i) = y(4-i)
    end do

    return

  end subroutine MAT_IWAN_principal_stress
!

end module mat_iwan





