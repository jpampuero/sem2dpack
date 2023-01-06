module mat_damage
!
! Damage rheology following the formulation by Vladimir Lyakhovsky
! as in Lyakhovsky et al (1997)
! See also Hamiel et al (2004) for the parameter beta
!
! . Stress-strain relation:
!
!     sigma_ij = ( lambda*i1 - gamma*sqrt(i2) )*delta_ij + ( 2*mu - gamma*i1/sqrt(i2) )*e_ij
!
!     where sigma_ij = total stress
!           e_ij = total elastic strain
!           i1 and i2 = strain invariants
!
!
! . Evolution of elastic moduli:
!
!     lambda = lambda_0 constant
!     mu = mu_0 + gamma_r * xi_0 * alpha
!     gamma = gamma_r * alpha^(1+beta) / (1+beta)
!
!     where alpha = damage state variable, in [0,1]
!           xi_0 = strain yield threshold (<0)
!
!
! . Damage evolution:
!
!     dalpha/dt = Cd*i2*positive_part[ xi * alpha^beta  - xi_0 ]
!
!     where xi = i1/sqrt(i2)
!
!
! . Damage-related plasticity:
!
!     if dalpha/dt > 0, plastic strain rate = dep_ij/dt = tij * Cv *dalpha/dt
!
!     where tij = deviatoric stress
!           Cv = R/mu with R=O(1)

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
 !-- damage
  type matwrk_dmg_type
    private
!    double precision, pointer, dimension(:,:) :: mu => null(), lambda => null()
    double precision :: mu, lambda
    double precision, pointer, dimension(:,:,:) :: ep => null()
    double precision, pointer, dimension(:) :: e0=>null(), s0=>null()
    double precision, pointer, dimension(:,:) :: alpha=>null()
    double precision :: xi_0=0d0,gamma_r=0d0,beta=0d0,Cd=0d0,Cv=0d0
  end type matwrk_dmg_type

  integer, save :: isDamage = 0

  ! for memory report
  integer, save :: MAT_DMG_mempro = 0
  integer, save :: MAT_DMG_memwrk = 0

  public :: matwrk_dmg_type &
          , MAT_isDamage, MAT_anyDamage, MAT_DMG_read, MAT_DMG_init_elem_prop &
          , MAT_DMG_init_elem_work, MAT_DMG_stress, MAT_DMG_export &
          , MAT_DMG_mempro, MAT_DMG_memwrk
  
contains

!=======================================================================
  logical function MAT_isDamage(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isDamage = MAT_isKind(m,isDamage)
  end function MAT_isDamage

!=======================================================================
  logical function MAT_anyDamage()
  MAT_anyDamage = (isDamage>0)
  end function MAT_anyDamage

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_DAMAGE
! GROUP  : MATERIALS
! PURPOSE: Set material properties for the damage rheology of 
!          Lyakhovsky, Ben-Zion and Agnon (J. Geophys. Res. 1997) 
!          and Hamiel et al (Geophys. J. Int. 2004)
! SYNTAX : &MAT_DAMAGE cp,cs,rho,phi,alpha,Cd,R,e0,ep /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: phi      [dble][0d0] internal friction angle
! ARG: alpha    [dble][0d0] initial value of damage variable
! ARG: Cd       [dble][0d0] damage evolution coefficient
! ARG: R        [dble][0d0] damage-related plasticity coefficient Cv
!                 normalized by the inverse of the intact shear modulus
! ARG: e0       [dble(3)][0d0] initial total strain (11, 22 and 12)
! ARG: ep       [dble(3)][0d0] initial plastic strain (11, 22 and 12)
!
! END INPUT BLOCK

! DEVEL: still testing beta>0
! SYNTAX : &MAT_DAMAGE cp,cs,rho,phi,alpha,Cd,beta,R,e0,ep /
! ARG: beta     [dble][0d0] damage evolution exponent

  subroutine MAT_DMG_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho,cp,cs,phi,alpha,Cd,beta,R,e0(3),ep(3)

  NAMELIST / MAT_DAMAGE / cp,cs,rho,phi,alpha,Cd,beta,R,e0,ep
  
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

  call MAT_setProp(input,'e11_0',e0(1))
  call MAT_setProp(input,'e22_0',e0(2))
  call MAT_setProp(input,'e12_0',e0(3))

  call MAT_setProp(input,'e11_p',ep(1))
  call MAT_setProp(input,'e22_p',ep(2))
  call MAT_setProp(input,'e12_p',ep(3))

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

!-----------------------------------------------------------------------
  elemental double precision function xi_zero_3d(phi,lambda,mu)

  use constants, only : PI

  double precision, intent(in) :: phi,lambda,mu
  double precision :: q

  q = sin(phi* PI/180d0)
  q = q / (1d0- q/3d0)
  xi_zero_3d = -sqrt(3d0) / sqrt(2d0*q*q*(lambda/mu +2d0/3d0)**2 +1d0)

  end function xi_zero_3d
!-----------------------------------------------------------------------
  elemental double precision function xi_zero_2d(phi,lambda,mu)

  use constants, only : PI

  double precision, intent(in) :: phi,lambda,mu
  double precision :: q

  q = sin(phi* PI/180d0)
  xi_zero_2d = -sqrt(2d0) / sqrt(q*q*(lambda/mu +1d0)**2 +1d0)

  end function xi_zero_2d

!-----------------------------------------------------------------------
  elemental double precision function gamma_r_2d(xi0,lambda,mu)

  double precision, intent(in) :: xi0,lambda,mu
  double precision :: q,p

  q = 2d0*(mu+lambda)/(2d0-xi0**2)
  p = 0.5d0*xi0*(q + lambda)
  gamma_r_2d = p + sqrt( p*p + 2d0*mu*q )

  end function gamma_r_2d

!-----------------------------------------------------------------------
  elemental double precision function gamma_r_3d(xi0,lambda,mu)

  double precision, intent(in) :: xi0,lambda,mu
  double precision :: q,p

  q = (2d0*mu+3d0*lambda)/(3d0-xi0**2)
  p = 0.5d0*xi0*(q + lambda)
  gamma_r_3d = p + sqrt( p*p + 2d0*mu*q )

  end function gamma_r_3d

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
  i1 = e(1) + e(2)
  i2 = e(1)*e(1) +e(2)*e(2) +2d0*e(3)*e(3)
  si2 = sqrt(i2)
  if (si2<1d-10) then
    xi = 0d0
  else
    xi = i1/si2
  endif

  two_mue = 2d0*rm-rg*xi

 !-- absolute stress tensor
  s(1) = rl*i1 - rg*si2 + two_mue*e(1)
  s(2) = rl*i1 - rg*si2 + two_mue*e(2)
  s(3) = two_mue*e(3)

!return !DEBUG: include this line to turn OFF convexity checks, at your own risk

 !-- check if the damage is above the critical value (loss of convexity)
 ! 2d plane strain:
  p = -(4d0*rm+2d0*rl-3d0*rg*xi) 
  q = two_mue*two_mue + two_mue*(2d0*rl-rg*xi) + rg*(rl*xi-rg)*(2d0-xi*xi)
 ! 3d:
  !p = -(4d0*rm+3d0*rl-3d0*rg*xi)
  !q = two_mue*two_mue + two_mue*(3d0*rl-rg*xi) + rg*(rl*xi-rg)*(3d0-xi*xi)
  d = p*p/4d0 - q
  if (d <= 0d0) call IO_abort('mat_damage:elastic: discriminant < 0')
  if ( p/2d0+sqrt(d) >= 0d0 ) call IO_abort('MAT_DMG: damage exceeded critical value (1st type)')
  if ( two_mue <= 0d0 ) call IO_abort('MAT_DMG: damage exceeded critical value (2nd type)')

  end subroutine compute_stress

!=======================================================================
! export output data

  function MAT_DMG_export(m) result(dat)

  type(matwrk_dmg_type), intent(in) :: m
  real :: dat(size(m%alpha,1),size(m%alpha,2),4)  ! 4 = 1(alpha)+size(ep,3)

  dat(:,:,1) = real(m%alpha)
  dat(:,:,2:4) = real(m%ep)
  
  end function MAT_DMG_export

end module mat_damage
