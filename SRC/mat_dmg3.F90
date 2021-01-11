module mat_dmg3
!
! Damage rheology for only mode 3 (antiplane shear) 
! following the formulation by Vladimir Lyakhovsky, as in Lyakhovsky et al (1997)
! but modified by Chao Liang. Here, I drop the interaction between i1 and i2 and 
! focus only on antiplane strain components
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
!     mu     = mu_0 + mu_r * mu_0 * alpha
!
!     where alpha = damage state variable, in [0,1]
!           mu_r  = (-1, 0], (1 + mu_r)*mu_0 is minimal 
!                   shear modulus when alpha = 1 
!
!  By choosing |mu_r| < 1, convexity is guaranteed
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
!              f0, c0 initial static friction, cohesion 
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
    double precision, pointer, dimension(:,:,:) :: ep      => null()
    double precision, pointer, dimension(:)     :: e0      => null()
    double precision, pointer, dimension(:)     :: s0      => null()
    double precision, pointer, dimension(:, :)  :: alpha   => null()
    double precision, pointer, dimension(:, :)  :: e_cmp   => null()
    double precision, pointer, dimension(:, :)  :: i2_cr   => null()
    double precision, pointer, dimension(:, :)  :: alpha0 => null()
    double precision, pointer, dimension(:, :)  :: t0     => null()
    double precision :: mu0 = 0d0, mu_r = 0d0, Cd = 0d0, Cv=0d0 
    double precision :: C1 = 0d0, C2 = 0d0, Th = 0d0, Rp = 0d0
    character(5) :: healLaw = 'EXP'
  end type matwrk_dmg3_type

  integer, save :: isDmg3 = 0

  ! for memory report
  integer, save :: MAT_DMG3_mempro = 0
  integer, save :: MAT_DMG3_memwrk = 0

  public :: matwrk_dmg3_type &
          , MAT_isDmg3, MAT_anyDmg3, MAT_DMG3_read, MAT_DMG3_init_elem_prop &
          , MAT_DMG3_init_elem_work, MAT_DMG3_stress, MAT_DMG3_export &
          , MAT_DMG3_mempro, MAT_DMG3_memwrk
  
contains

!=======================================================================
  logical function MAT_isDmg3(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isDmg3 = MAT_isKind(m,isDmg3)
  end function MAT_isDmg3

!=======================================================================
  logical function MAT_anyDmg3()
  MAT_anyDmg3 = (isDmg3>0)
  end function MAT_anyDmg3

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_DMG3
! GROUP  : MATERIALS
! PURPOSE: Set material properties for simplified damage rheology 
!          of mode 3, antiplane shear problem.
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: mu_r     [dble][0d0] damage ratio for mu (-1, 0]
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
! SYNTAX : &MAT_DMG3 alpha|alphaH, Sm0|Sm0H, e0, ep 
!                    cp, cs, rho, mu_r, f0, c0, Cd, C1, C2, Rp, Th, R /
!
! Note that: alphaH, Sm0H can be spatial distributions
!

  subroutine MAT_DMG3_read(input,iin)

  use echo, only : echo_input, iout
  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho, cp, cs, f0, c0, alpha, Cd, C1, C2
  double precision :: mu_r, Rp, Th, R, ep(2), e0(2), sm0
  character(20):: alphaH, sm0H, healLaw

  NAMELIST / MAT_DAMAGE / alpha, sm0, e0, ep, alphaH, sm0H, mu_r, & 
                          cp, cs, rho, f0, c0, Cd, C1, C2, R, Rp, Th
  
  ! undamaged elastic property
  cp    = 0d0
  cs    = 0d0
  rho   = 0d0

  ! initial total/plastic strain
  e0    = 0d0
  ep    = 0d0

  ! static friction strength
  f0    = 0d0
  c0    = 0d0
  sm0   = 0d0
  sm0H  = ''

  ! damage parameters
  alpha = 0d0
  alphaH = ''
  mu_r  = 0d0
  Cd    = 0d0
  R     = 0d0 ! viscous relaxation

  ! parameters for healing
  C1    = 0d0
  C2    = 0d0
  Rp    = 0d0
  Th    = 0d0

  ! read parameters
  read(iin, MAT_DAMAGE, END=100)

  if (Rp>0 .and. Th>0) then
      healLaw = 'Exponential'
  else
      healLaw = 'Logarithmic'
  end if

  call MAT_setKind(input,isDmg3)

  ! undamaged properties
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  
  call MAT_setProp(input,'e31_0',e0(1))
  call MAT_setProp(input,'e32_0',e0(2))
  call MAT_setProp(input,'e31_p',ep(1))
  call MAT_setProp(input,'e32_p',ep(2))

  ! yielding criteria
  call MAT_setProp(input,'f0',f0)
  call MAT_setProp(input,'c0',c0)
  call MAT_setProp(input,'sm0',sm0,sm0H,iin,sm0H)

  ! damage parameters
  call MAT_setProp(input,'alpha',alpha,alphaH,iin,alphaH)
  call MAT_setProp(input,'Cd', Cd)
  call MAT_setProp(input,'mu_r',mu_r)
  call MAT_setProp(input,'R', R)

  ! healing parameters
  ! Logarithmic healing parameters
  call MAT_setProp(input,'C1', C1)
  call MAT_setProp(input,'C2', C2)

  ! exponential healing parameters
  call MAT_setProp(input,'Rp', Rp)
  call MAT_setProp(input,'Th', Th)

  write(iout,200) cp, cs, rho, f0, c0, &
                  sm0H, alphaH, Cd, mu_r, R, &
                  C1, C2, Rp, Th, healLaw, e0, ep

  return

  100 call IO_abort('MAT_DMG3_read: MAT_DMG3 input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Yield friction . . . . . . . . . . . (f0) =',EN12.3,/5x, &
    'Yield cohesion . . . . . . . . . . . (c0) =',EN12.3,/5x, &
    'Mean confining stress . . . . . . . (sm0) =', A,/5x, &
    'Initial damage variable . . . . . (alpha) =', A,/5x, &
    'Damage evolution coefficient. . . . .(Cd) =',EN12.3,/5x, &
    'Damage modulus ratio . . . . . . . (mu_r) =',EN12.3,/5x, &
    'Damage-related plasticity coefficient (R) =',EN12.3,/5x, &
    'Damage healing coefficient . . . . . (C1) =',EN12.3,/5x, &
    'Damage healing coefficient . . . . . (C2) =',EN12.3,/5x, &
    'Damage healing coefficient . . . . . (Rp) =',EN12.3,/5x, &
    'Damage healing coefficient . . . . . (Th) =',EN12.3,/5x, &
    'Healing law. . . . . . . . . . . . . . .  =', A ,/5x, &
    'Initial total strain 31 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     32 . . . . . (e0(2)) =',EN12.3,/5x, &
    'Initial plastic strain 31 . . . . (ep(1)) =',EN12.3,/5x, &
    '                       32 . . . . (ep(2)) =',EN12.3)

  end subroutine MAT_DMG3_read

! Initialize elemental properties
!=======================================================================
  subroutine MAT_DMG3_init_elem_prop(elem,ecoord)
  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)
  double precision :: rho, cp, cs

  call MAT_setProp(elem,'cp',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'cs',ecoord,MAT_DMG3_mempro)
  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_DMG3_mempro)
  call MAT_setProp(elem,'mu',rho*cs*cs,MAT_DMG3_mempro)

  call MAT_setProp(elem,'e31_0',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'e32_0',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'e31_p',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'e32_p',ecoord,MAT_DMG3_mempro)

  call MAT_setProp(elem,'f0',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'c0',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'sm0',ecoord,MAT_DMG3_mempro)

  call MAT_setProp(elem,'alpha',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'Cd',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'mu_r',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'R',ecoord,MAT_DMG3_mempro)

  call MAT_setProp(elem,'C1',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'C2',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'Rp',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'Th',ecoord,MAT_DMG3_mempro)

  end subroutine MAT_DMG3_init_elem_prop

!-----------------------------------------------------------------------
! construct working data type
!======================================================================
  subroutine MAT_DMG3_init_elem_work(m,p,ngll)

    type(matwrk_dmg3_type), intent(inout) :: m
    type(matpro_elem_type), intent(in) :: p
    integer, intent(in) :: ngll
    double precision :: R, f0, c0, sm0(ngll, ngll), mu0, lambda0, mu
   
    allocate(m%alpha(ngll,ngll))
    allocate(m%e0(2))
    allocate(m%s0(2))
    allocate(m%ep(ngll,ngll,2))
    allocate(m%i2_cr(ngll,ngll))
    
   
    ! elastic coefficients for intact material (zero damage)
    call MAT_getProp(mu0,p,'mu')
    call MAT_getProp(lambda0,p,'lambda')
    m%mu0 = mu0
   
    ! constants for damage, healing, yielding
    call MAT_getProp(m%Cd, p, 'Cd')
    call MAT_getProp(m%C1, p, 'C1')
    call MAT_getProp(m%C2, p, 'C2')
    call MAT_getProp(m%mu_r,p,'mu_r')
   
    ! Cv = R/mu = coefficient for damage-related plasticity 
     call MAT_getProp(R,p,'R')
     m%Cv = R / m%mu0
   
    ! damage variable 
    call MAT_getProp(m%alpha, p, 'alpha')

    ! exp healing
    call MAT_getProp(m%Rp,p,'Rp')
    call MAT_getProp(m%Th,p,'Th')
    
    if (m%Rp>0 .and. m%Th>0) then
        m%healLaw = 'EXP'
       allocate(m%alpha0(ngll, ngll))
       allocate(m%t0(ngll, ngll))
       m%alpha0 = m%alpha
       m%t0     = 0d0 
    else
        m%healLaw = 'LOG'
        allocate(m%e_cmp(ngll,ngll)) 
    end if
   
    ! i2_cr = damage yield threshold for i2
    call Mat_getProp(sm0, p, 'sm0')
    call Mat_getProp(f0, p, 'f0')
    call Mat_getProp(c0, p, 'c0')
   
    m%i2_cr = compute_i2_cr(f0, mu0, c0, sm0) 
    if (m%healLaw=='LOG') m%e_cmp = compute_ecmp(sm0, lambda0, mu0) 
   
    ! initial strain
    call MAT_getProp(m%e0(1), p, 'e31_0')
    call MAT_getProp(m%e0(2), p, 'e32_0')
   
    ! initial plastic strain
    call MAT_getProp(m%ep(:,:,1), p, 'e31_p')
    call MAT_getProp(m%ep(:,:,2), p, 'e32_p')
   
    ! initial antiplane shear stress, uniform
    mu = m%mu0 + m%mu0 * m%mu_r * m%alpha(1,1)
    m%s0 = 2d0 * mu * (m%e0 - m%ep(1,1,:)) 
    
    ! track memory usage
    MAT_DMG3_memwrk = MAT_DMG3_memwrk &
                      + size( transfer(m, (/ 0d0 /) )) &
                      + size(m%alpha) &
                      + size(m%ep) &
                      + size(m%e0) &
                      + size(m%s0) &
                      + size(m%i2_cr)
    
    if (associated(m%alpha0)) MAT_DMG3_memwrk = MAT_DMG3_memwrk + size(m%alpha0)
    if (associated(m%t0)) MAT_DMG3_memwrk = MAT_DMG3_memwrk + size(m%t0)
    if (associated(m%e_cmp)) MAT_DMG3_memwrk = MAT_DMG3_memwrk + size(m%e_cmp)

end subroutine MAT_DMG3_init_elem_work

! compute the critical I2 beyond which damage starts
! allow sm0 be spatially variable
!-----------------------------------------------------------------------
  function compute_i2_cr(f0, mu0, c0, sm0) result(i2_cr)
      double precision, intent(in) :: f0, mu0, c0, sm0(:,:)
      double precision, dimension(size(sm0,1), size(sm0,2)) :: i2_cr

      i2_cr = 0.5d0* ((-f0*sm0+c0)/mu0)**2d0 
      
  end function compute_i2_cr

!- compute the compaction strain for healing
!-----------------------------------------------------------------------
  function compute_ecmp(sm0, lambda0, mu0) result(e_cmp)
      double precision, intent(in) :: sm0(:, :), lambda0, mu0
      double precision, dimension(size(sm0,1), size(sm0,2)) :: e_cmp
      double precision :: K0

      ! compute bulk modulus
      K0 = lambda0 + 2d0/3d0*mu0

      ! compute the compaction strain
      e_cmp = abs(Sm0/K0)
  end function compute_ecmp

!=======================================================================
! Constitutive law:
! compute relative stresses from given strains
! and (if requested) update damage variable and plastic strain
!
! subroutine MAT_DMG3_stress(s,e,m,ngll,update,dt)
  subroutine MAT_DMG3_stress(s,etot,m,ngll,update,dt,E_ep,E_el,sg,sgp)

  use utils, only: positive_part
  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  integer, intent(in) :: ngll
  double precision, intent(in) :: etot(ngll,ngll,2)
  double precision, intent(out) :: s(ngll,ngll,2)
  type (matwrk_dmg3_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt
  double precision, intent(out) :: E_ep(ngll,ngll)
  double precision, intent(out) :: E_el(ngll,ngll)
  double precision, intent(out) :: sg(ngll,ngll,2)
  double precision, intent(out) :: sgp(ngll,ngll,2)

  double precision :: e(ngll,ngll,2)
  double precision, dimension(2) :: eij,sij

  double precision :: i2_0
  double precision, dimension(ngll,ngll) :: i1, i2, rl,rm,rg, xi, sm,dalpha
  double precision, dimension(ngll,ngll,2) :: dep
  integer :: i,j

 !!-- total strain
 ! e(:,:,1) = etot(:,:,1) + m%e0(1)
 ! e(:,:,2) = etot(:,:,2) + m%e0(2)

 !!-- elastic strain
 ! e = e - m%ep

 !!-- damaged elastic moduli
 ! rl = m%lambda
 ! rm = m%mu + m%xi_0 * m%gamma_r * m%alpha
!!  rl = m%lambda(1,1)
!!  rm = m%mu(1,1) + m%xi_0 * m%gamma_r * m%alpha
 ! rg = m%gamma_r * (m%alpha**(1d0+m%beta)) / (1d0+m%beta) ! Hamiel et al (2004) eq 3

 !!-- compute stresses and two strain invariants
 ! do j=1,ngll
 ! do i=1,ngll
 !   eij = e(i,j,:)
 !   call compute_stress(sij,eij,rl(i,j),rm(i,j),rg(i,j),i1(i,j),i2(i,j),xi(i,j))
 !   s(i,j,:) = sij
 ! enddo
 ! enddo

!!------- Update damage variables and plastic strain --------------- 
 ! if (update) then

 !   if (.not.present(dt)) &
 !     call IO_abort('mat_dmg3:MAT_DMG3_stress: update requested but argument dt is absent')

 !  !-- damage evolution
!!    if (m%beta==0d0) then
 !     dalpha = dt*m%Cd*i2*positive_part( xi - m%xi_0 )
!!    else
 !    ! from Hamiel et al (2004) eq 7
!!      dalpha = dt*m%Cd*i2*positive_part( xi * m%alpha**m%beta  - m%xi_0 ) 
!!    endif
 !   m%alpha = m%alpha + dalpha
 ! 
 !  !-- plasticity update
 !  ! Damage-related viscosity, if alpha_dot > 0
 !  ! Plastic strain rate = deij/dt = tij * Cv *dalpha/dt
 !  ! where tij = deviatoric stress

 !   dalpha = m%Cv * positive_part(dalpha)
 !   dep(:,:,1) = s(:,:,1)*dalpha
 !   dep(:,:,2) = s(:,:,2)*dalpha
 !   m%ep = m%ep + dep

 !   if (COMPUTE_ENERGIES) then
 !    ! increment of plastic energy dissipation
 !     E_ep = s(:,:,1)*dep(:,:,1) + s(:,:,2)*dep(:,:,2) + 2d0*s(:,:,3)*dep(:,:,3)
 !    ! total elastic energy change
 !    ! WARNING: assumes zero initial plastic strain
 !     i2_0 = m%e0(1)*m%e0(1) +m%e0(2)*m%e0(2) 
!!      E_el = 0.5d0*(rl*i1*i1 - m%lambda*i1_0*i1_0) + ( rm*i2 - m%mu*i2_0 ) - rg*i1*sqrt(i2)
 !   else
 !     E_ep = 0d0
 !     E_el = 0d0
 !   endif

 !   if (COMPUTE_STRESS_GLUT) then
 !     rm = 2d0*m%mu
 !    ! damage components
 !     sg(:,:,1) = s(:,:,1) - (rl+rm)*e(:,:,1) - rl*e(:,:,2)
 !     sg(:,:,2) = s(:,:,2) - rl*e(:,:,1) - (rl+rm)*e(:,:,2)
 !     sg(:,:,3) = s(:,:,3) - rm*e(:,:,3)
 !    ! plastic components
 !     sgp(:,:,1) = - rm*m%ep(:,:,1)
 !     sgp(:,:,2) = - rm*m%ep(:,:,2)
 !     sgp(:,:,3) = - rm*m%ep(:,:,3)
 !   else
 !     sg = 0d0
 !     sgp = 0d0
 !   endif

 ! endif

 !!-- relative stresses
 ! s(:,:,1) = s(:,:,1) - m%s0(1)
 ! s(:,:,2) = s(:,:,2) - m%s0(2)

  end subroutine MAT_DMG3_stress

!-------------------------------------------------------------------
!
  function compute_stress(e, mu) result(s)
      double precision, intent(in) :: e(2), mu
      double precision :: s(2)
      s = 2d0*mu*e
  end function compute_stress

!=======================================================================
! export output data

  function MAT_DMG3_export(m) result(dat)

      type(matwrk_dmg3_type), intent(in) :: m
      real :: dat(size(m%alpha,1),size(m%alpha,2),3)  ! 4 = 1(alpha)+size(ep,3)

      dat(:,:,1) = real(m%alpha)
      dat(:,:,2:3) = real(m%ep)
  
  end function MAT_DMG3_export

end module mat_dmg3
