module mat_dmg3
!
! Damage rheology for mode 3 (antiplane shear) 
! following the formulation by Vladimir Lyakhovsky, as in Lyakhovsky et al (1997)
! drop the interaction between i1 and i2 and focus only on antiplane strain components
!
! . Stress-strain relation:
!
!     sigma_ij = 2*mu*e_ij
!
!     where sigma_ij = total stress
!           e_ij = total elastic strain
!
! . Evolution of elastic moduli:
!
!     lambda = lambda_0 constant
!     mu     = mu_0 + mu_r * alpha
!
!     where alpha = damage state variable, in [0,1]
!           mu_r  = damage modulus, mu_0 + mu_r is minimal 
!                   shear modulus when alpha = 1 
!
!  By choosing |mu_r|<mu_0, convexity is guaranteed
!
!
! . Damage evolution:
!
!   i2>i2_cr, damage increase:
!
!        dalpha/dt = Cd*i2 (positive)
!
!        where i2 is e_{31}^2 + e_{32}^2 
!              i2_cr = 0.5*[(-f0*sigma_0+c0)/mu_0]^2
!              sigma_0 is the initial mean stress [negative]
!              f0,c0 initial static friction, cohesion 
!
!        . Damage-related plasticity:
!
!        if dalpha/dt > 0, plastic strain rate = dep_ij/dt = tij * Cv *dalpha/dt
!
!          where tij = deviatoric stress
!                 Cv = R/mu_0 with R=O(1)
!          [Cv is on the order of 1/mu_0] 
!
!   i2<=i2_cr, damage decrease, healing:
!
!       .Damage law similar to Lyakhovsky et al., 1997 (logarithmic healing)
!
!           dalpha/dt = C1 * exp(alpha/C2) * e_cmp^2 
!
!           alpha = alpha_n - C2*ln(10)*log10(1-C1/C2*exp(alpha_n/C2)*e_cmp^2*dt)
!
!           where C1<0, C2>0 are constants
!                 e_cmp = compaction strain ~-sigma_0/K0
!                 value of C1 and C2 are estimated according to rate state experiment
!
!                C2*ln(10)~A, C1~-B*C2*exp(-1/C2)/e_cmp^2
!                given, A~0.01, B~1-2 s^-1, e_cmp~10^-2 
!                C2=0.05, C1~10^-6; or C2 = 0.03, C1~10^-12.
!                C1 is very sensitive to C2
!
!
!       .empirical healing law, exponential healing
!           alpha  = alpha_0 - (alpha_0 - alpha_0*Rp) * (1 - exp(-t0/Th))
!           where t0 is the time to alpha_0, Th is the healing time
!                 alpha0 is the initial damage level entering for current healing period
!                 Rp is the ratio of permenant damage
!           Note both alpha_0, t0 can be evolves over time can be spatially dependent.
!

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
 !-- damage
  type matwrk_dmg3_type
    private
!    double precision, pointer, dimension(:,:) :: mu => null(), lambda => null()
    double precision :: mu
    double precision, pointer, dimension(:,:,:) :: ep      => null()
    double precision, pointer, dimension(:,:,:) :: e0      => null()
    double precision, pointer, dimension(:, :)  :: alpha   => null()
    double precision, pointer, dimension(:, :)  :: e_cmp   => null()
    double precision, pointer, dimension(:, :)  :: i2_cr   => null()
    double precision, pointer, dimension(:, :)  :: alpha_0 => null()
    double precision, pointer, dimension(:, :)  :: t_0     => null()
    double precision :: mu_0 = 0d0, mu_r = 0d0 
    double precision :: C1 = 0d0, C2 = 0d0, Cd = 0d0, Cv = 0d0, Th = 0d0
  end type matwrk_dmg3_type

  integer, save :: isDamage3 = 0

  ! for memory report
  integer, save :: MAT_DMG3_mempro = 0
  integer, save :: MAT_DMG3_memwrk = 0

  public :: matwrk_dmg3_type &
          , MAT_isDamage3, MAT_anyDamage3, MAT_DMG3_read, MAT_DMG3_init_elem_prop &
          , MAT_DMG3_init_elem_work, MAT_DMG3_stress, MAT_DMG3_export &
          , MAT_DMG3_mempro, MAT_DMG3_memwrk
  
contains

!=======================================================================
  logical function MAT_isDamage3(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isDamage3 = MAT_isKind(m,isDamage3)
  end function MAT_isDamage3

!=======================================================================
  logical function MAT_anyDamage3()
  MAT_anyDamage3 = (isDamage3>0)
  end function MAT_anyDamage3

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_DAMAGE3
! GROUP  : MATERIALS
! PURPOSE: Set material properties for simplified damage rheology 
!          of mode 3, antiplane shear problem.
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: f0       [dble][0d0] internal friction coefficient
! ARG: c0       [dble][0d0] internal cohesion
! ARG: Sm0      [dble][0d0] initial value of mean stress [negative] 
! ARG: alpha    [dble][0d0] initial value of damage variable
! ARG: Cd       [dble][0d0] damage evolution coefficient
! ARG: C1       [dble][0d0] damage evolution coefficient [healing]
! ARG: C2       [dble][0d0] damage evolution coefficient [healing]
! ARG: Rp       [dble][0d0] ratio of permanent damage [healing] 
! ARG: Th       [dble][0d0] time scale for exponential recovery [healing]
! ARG: R        [dble][0d0] damage-related plasticity coefficient Cv
!                 normalized by the inverse of the intact shear modulus
! ARG: e0       [dble(2)][0d0] initial total strain (31, 32)
! ARG: ep       [dble(2)][0d0] initial plastic strain (31, 32)
!
! END INPUT BLOCK
!
! SYNTAX : &MAT_DAMAGE cp, cs, rho, f0, c0, Sm0, Sm0H, e1H, e2H, 
!                      alpha, Cd, C1, C2, Rp, Th, R, e0, ep /
!

  subroutine MAT_DMG_read(input,iin)

  use distribution_cd, only : cd_type
  use echo, only : echo_input, iout
  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho,cp, cs, f0, c0, alpha, Cd, C1, C2
  double precision :: Rp, Th, R, ep(2)
  double precision :: e0(2) 
  cd_type :: cd_e1, cd_e2, cd_sm0, cd_alpha,

  NAMELIST / MAT_DAMAGE / cp,cs,rho,phi,alpha,C,beta,R,e0,ep
  
  cp    = 0d0
  cs    = 0d0
  rho   = 0d0
  Cd    = 0d0
  R     = 0d0
  beta  = 0d0
  e0    = 0d0
  ep    = 0d0
  phi   = 0d0
  alpha = 0d0

  read(iin, MAT_DAMAGE, END=100)

  write(iout,200) cp,cs,rho,phi,alpha,Cd,beta,R,e0,ep

  call MAT_setKind(input,isDamage)
  
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'Cd',Cd)
  call MAT_setProp(input,'R',R)
  call MAT_setProp(input,'beta',beta)
  call MAT_setProp(input,'phi',phi)
  call MAT_setProp(input,'alpha',alpha)

  call MAT_setProp(input,'e31_0',e0(1))
  call MAT_setProp(input,'e32_0',e0(2))

  call MAT_setProp(input,'e31_p',ep(1))
  call MAT_setProp(input,'e32_p',ep(2))

  return

  100 call IO_abort('MAT_DMG_read: MAT_DAMAGE input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Yield angle . . . . . . . . . . . . (phi) =',EN12.3,/5x, &
    'Initial damage variable . . . . . (alpha) =',EN12.3,/5x, &
    'Damage evolution coefficient. . . . .(Cd) =',EN12.3,/5x, &
    'Damage evolution exponent . . . . .(beta) =',EN12.3,/5x, &
    'Damage-related plasticity coefficient (R) =',EN12.3,/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x, &
    '                     12 . . . . . (e0(3)) =',EN12.3,/5x, &
    'Initial plastic strain 11 . . . . (ep(1)) =',EN12.3,/5x, &
    '                       22 . . . . (ep(2)) =',EN12.3,/5x, &
    '                       12 . . . . (ep(3)) =',EN12.3)

  end subroutine MAT_DMG_read

!=======================================================================
  subroutine MAT_DMG_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision :: rho,cp,cs

  call MAT_setProp(elem,'cp',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'cs',ecoord,MAT_DMG_mempro)

  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_DMG_mempro)
  call MAT_setProp(elem,'mu',rho*cs*cs,MAT_DMG_mempro)

  call MAT_setProp(elem,'Cd',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'R',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'beta',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'phi',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'alpha',ecoord,MAT_DMG_mempro)

  call MAT_setProp(elem,'e11_0',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'e22_0',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'e12_0',ecoord,MAT_DMG_mempro)

  call MAT_setProp(elem,'e11_p',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'e22_p',ecoord,MAT_DMG_mempro)
  call MAT_setProp(elem,'e12_p',ecoord,MAT_DMG_mempro)

  end subroutine MAT_DMG_init_elem_prop

!-----------------------------------------------------------------------
  subroutine MAT_DMG_init_elem_work(m,p,n)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights

  type(matwrk_dmg_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: n

  double precision :: lambda,mu,phi,R,rg,i1,i2,xi

!  allocate(m%lambda(n,n)) ! NOTE: assume uniform intact elastic moduli
!  allocate(m%mu(n,n))     ! NOTE: assume uniform intact elastic moduli
  allocate(m%alpha(n,n))
  allocate(m%e0(3))
  allocate(m%ep(n,n,3))
  allocate(m%s0(3))

 ! elastic coefficients for intact material (zero damage)
  call MAT_getProp(lambda,p,'lambda')
  call MAT_getProp(mu,p,'mu')
  m%lambda = lambda
  m%mu = mu

 ! damage variable and its time derivative
  call MAT_getProp(m%alpha,p,'alpha')

 ! xi_0 = damage yield threshold on strain
 !        derived from the internal friction angle phi
  call MAT_getProp(phi,p,'phi')
  m%xi_0 = xi_zero_2d(phi,lambda,mu)

 ! gamma_r = damage influence coefficient
 !           such that convexity is lost at xi=xi_0 and alpha=1
  m%gamma_r = gamma_r_2d(m%xi_0,lambda,mu)
!  m%gamma_r = 0d0 !test: no damage

 ! Cd = damage evolution coefficient 
  call MAT_getProp(m%Cd,p,'Cd')

 ! beta = damage evolution exponent 
  call MAT_getProp(m%beta,p,'beta')

 ! Cv = R/mu = coefficient for damage-related plasticity 
  call MAT_getProp(R,p,'R')
  m%Cv = R / mu

 ! initial strain
  call MAT_getProp(m%e0(1),p,'e11_0')
  call MAT_getProp(m%e0(2),p,'e22_0')
  call MAT_getProp(m%e0(3),p,'e12_0')

 ! initial plastic strain
  call MAT_getProp(m%ep(:,:,1),p,'e11_p')
  call MAT_getProp(m%ep(:,:,2),p,'e22_p')
  call MAT_getProp(m%ep(:,:,3),p,'e12_p')

 ! initial stress
  mu = mu + m%xi_0 * m%gamma_r * m%alpha(1,1)
  rg = m%gamma_r * (m%alpha(1,1)**(1d0+m%beta)) / (1d0+m%beta)
  call compute_stress(m%s0,m%e0-m%ep(1,1,:),lambda,mu,rg, i1,i2,xi)

  MAT_DMG_memwrk = MAT_DMG_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%alpha) &
                   + size(m%ep) &
                   + size(m%e0) &
                   + size(m%s0)
!                   + size(m%lambda) &
!                   + size(m%mu) &

  end subroutine MAT_DMG_init_elem_work

! compute the critical I2 beyond which damage starts
!-----------------------------------------------------------------------
  double precision function compute_i2_cr(f0, mu0, Sm0) result(i2_cr)
      ! allow heterogeneous Sm0

  double precision, intent(in) :: f0, mu_0, Sm0(:,:)
  double precision, dimension(size(Sm0,1), size(Sm0,2)) :: i2_cr

  i2_cr = (f0*Sm0/sqrt(2.0d0)/mu_0)**2d0 
  
  end function compute_i2_cr

!- compute the compaction strain for healing
!-----------------------------------------------------------------------
  double precision function compute_ecmp(Sm0, lambda0, mu0) result(e_cmp)
  double precision, intent(in) :: Sm0(:, :), lambda0, mu0, K0
  double precision, dimension(size(Sm0,1), size(Sm0,2)) :: e_cmp

  ! compute bulk modulus
  K0 = lambda0 + 2d0/3d0*mu0

  ! compute the compaction strain
  e_cmp = abs(Sm0/K0)
  end function xi_zero_2d


!=======================================================================
! Constitutive law:
! compute relative stresses from given strains
! and (if requested) update damage variable and plastic strain
!
! subroutine MAT_DMG_stress(s,e,m,ngll,update,dt)
  subroutine MAT_DMG_stress(s,etot,m,ngll,update,dt,E_ep,E_el,sg,sgp)

  use utils, only: positive_part
  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  integer, intent(in) :: ngll
  double precision, intent(in) :: etot(ngll,ngll,3)
  double precision, intent(out) :: s(ngll,ngll,3)
  type (matwrk_dmg_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt
  double precision, intent(out) :: E_ep(ngll,ngll)
  double precision, intent(out) :: E_el(ngll,ngll)
  double precision, intent(out) :: sg(ngll,ngll,3)
  double precision, intent(out) :: sgp(ngll,ngll,3)

  double precision :: e(ngll,ngll,3)
  double precision, dimension(3) :: eij,sij

  double precision :: i1_0, i2_0
  double precision, dimension(ngll,ngll) :: i1, i2, rl,rm,rg, xi, sm,dalpha
  double precision, dimension(ngll,ngll,3) :: dep
  integer :: i,j

 !-- total strain
  e(:,:,1) = etot(:,:,1) + m%e0(1)
  e(:,:,2) = etot(:,:,2) + m%e0(2)
  e(:,:,3) = etot(:,:,3) + m%e0(3)

 !-- elastic strain
  e = e - m%ep

 !-- damaged elastic moduli
  rl = m%lambda
  rm = m%mu + m%xi_0 * m%gamma_r * m%alpha
!  rl = m%lambda(1,1)
!  rm = m%mu(1,1) + m%xi_0 * m%gamma_r * m%alpha
  rg = m%gamma_r * (m%alpha**(1d0+m%beta)) / (1d0+m%beta) ! Hamiel et al (2004) eq 3

 !-- compute stresses and two strain invariants
  do j=1,ngll
  do i=1,ngll
    eij = e(i,j,:)
    call compute_stress(sij,eij,rl(i,j),rm(i,j),rg(i,j),i1(i,j),i2(i,j),xi(i,j))
    s(i,j,:) = sij
  enddo
  enddo

  if (update) then

    if (.not.present(dt)) &
      call IO_abort('mat_damage:MAT_DMG_stress: update requested but argument dt is absent')

   !-- damage evolution
    if (m%beta==0d0) then
      dalpha = dt*m%Cd*i2*positive_part( xi - m%xi_0 )
    else
     ! from Hamiel et al (2004) eq 7
      dalpha = dt*m%Cd*i2*positive_part( xi * m%alpha**m%beta  - m%xi_0 ) 
    endif
    m%alpha = m%alpha + dalpha
  
   !-- plasticity update
   ! Damage-related viscosity, if alpha_dot > 0
   ! Plastic strain rate = deij/dt = tij * Cv *dalpha/dt
   ! where tij = deviatoric stress
    sm = 0.5d0*(s(:,:,1)+s(:,:,2)) ! mean stress
    dalpha = m%Cv * positive_part(dalpha)
    dep(:,:,1) = (s(:,:,1) - sm)*dalpha
    dep(:,:,2) = (s(:,:,2) - sm)*dalpha
    dep(:,:,3) = s(:,:,3)*dalpha
    m%ep = m%ep + dep

    if (COMPUTE_ENERGIES) then
     ! increment of plastic energy dissipation
      E_ep = s(:,:,1)*dep(:,:,1) + s(:,:,2)*dep(:,:,2) + 2d0*s(:,:,3)*dep(:,:,3)
     ! total elastic energy change
     ! WARNING: assumes zero initial plastic strain
      i1_0 = m%e0(1) + m%e0(2)
      i2_0 = m%e0(1)*m%e0(1) +m%e0(2)*m%e0(2) +2d0*m%e0(3)*m%e0(3)
      E_el = 0.5d0*(rl*i1*i1 - m%lambda*i1_0*i1_0) + ( rm*i2 - m%mu*i2_0 ) - rg*i1*sqrt(i2)
    else
      E_ep = 0d0
      E_el = 0d0
    endif

    if (COMPUTE_STRESS_GLUT) then
      rm = 2d0*m%mu
     ! damage components
      sg(:,:,1) = s(:,:,1) - (rl+rm)*e(:,:,1) - rl*e(:,:,2)
      sg(:,:,2) = s(:,:,2) - rl*e(:,:,1) - (rl+rm)*e(:,:,2)
      sg(:,:,3) = s(:,:,3) - rm*e(:,:,3)
     ! plastic components
      sgp(:,:,1) = - rm*m%ep(:,:,1)
      sgp(:,:,2) = - rm*m%ep(:,:,2)
      sgp(:,:,3) = - rm*m%ep(:,:,3)
    else
      sg = 0d0
      sgp = 0d0
    endif

  endif

 !-- relative stresses
  s(:,:,1) = s(:,:,1) - m%s0(1)
  s(:,:,2) = s(:,:,2) - m%s0(2)
  s(:,:,3) = s(:,:,3) - m%s0(3)

  end subroutine MAT_DMG_stress

!-------------------------------------------------------------------
! Stress-strain relation:
!   i1 and i2 = strain invariants
!   sigma_ij = ( lambda*i1 - gamma*sqrt(i2) )*delta_ij 
!             + ( 2*mu - gamma*i1/sqrt(i2) )*e_ij
!
  subroutine compute_stress(s,e,rl,rm,rg,i1,i2,xi)

  double precision, intent(in) :: e(3),rl,rm,rg
  double precision, intent(out) :: s(3),i1,i2,xi

  double precision :: si2,two_mue,p,q,d

 !-- invariants of the elastic strain tensor (2D plane strain)
  i2 = e(1)*e(1) +e(2)*e(2)

  two_mue = 2d0*rm-rg*xi
  s = two_mue*e

  end subroutine compute_stress

!=======================================================================
! export output data

  function MAT_DMG_export(m) result(dat)

  type(matwrk_dmg_type), intent(in) :: m
  real :: dat(size(m%alpha,1),size(m%alpha,2),3)  ! 4 = 1(alpha)+size(ep,3)

  dat(:,:,1) = real(m%alpha)
  dat(:,:,2:3) = real(m%ep)
  
  end function MAT_DMG_export

end module mat_dmg3
