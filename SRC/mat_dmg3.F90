module mat_dmg3
!
! Damage rheology for only mode 3 (antiplane shear) 
! following the formulation by Vladimir Lyakhovsky, as in Lyakhovsky et al (1997)
! but modified by Chao Liang, drop the interaction between i1 and i2 and 
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
!     mu     = mu_0 - mu_r * mu_0 * alpha
!
!     where alpha = damage state variable, in [0,1]
!           mu_r  = [0, 1), (1 - mu_r)*mu_0 is minimal 
!                   shear modulus when alpha = 1 
!
!  By choosing |mu_r| < 1, convexity is guaranteed
!
!
! . Damage evolution:
!
!   i2>i2_cr, damage increase:
!
!        dalpha/dt = Cd*(i2-i2_cr) (positive)
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
!           dalpha/dt = C1 * exp(alpha/C2) * (i2-i2_cr) 
!
!           alpha = alpha_n - C2*ln(10)*log10(1+C1/C2*exp(alpha_n/C2)*(i2-i2_cr)*dt)
!
!           where C1>0, C2>0 are constants
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
    double precision, pointer, dimension(:, :)  :: i2_cr   => null()
    double precision, pointer, dimension(:, :)  :: alpha0 => null()
    double precision, pointer, dimension(:, :)  :: Cd => null()
    double precision :: mu0 = 0d0, mu_r = 0d0, Cv=0d0 
    double precision :: C1 = 0d0, C2 = 0d0, Th = 0d0, Rp = 0d0
    double precision :: R  = 0d0, da_max=0d0, dtmax=1d10
    character(5) :: healLaw = 'EXP'
  end type matwrk_dmg3_type

  integer, save :: isDmg3 = 0
  character(10), save :: healLaw = 'EXP'
  logical :: update_in=.true.
  logical :: fix_cd=.true.
  double precision :: edot0=1.d-4, Cd0=1.d1, Cdm=0.d0
  double precision :: Cdmin = 0.1d0, Cdmax = 1.0d4 
  integer :: EqStart = 1
  ! implement strain rate dependent Cd
  !log10(Cd/Cd0) = 1 + Cdm * log(edot/edot0)

  ! for memory report
  integer, save :: MAT_DMG3_mempro = 0
  integer, save :: MAT_DMG3_memwrk = 0

  public :: matwrk_dmg3_type &
          , MAT_isDmg3, MAT_anyDmg3, MAT_DMG3_read, MAT_DMG3_init_elem_prop &
          , MAT_DMG3_init_elem_work, MAT_DMG3_stress, MAT_DMG3_export &
          , MAT_DMG3_mempro, MAT_DMG3_memwrk, MAT_DMG3_dtmax,MAT_DMG3_Cijkl &
          , MAT_DMG3_stress_ep,MAT_isUpdateDmg3, MAT_DMG3_Set_Cd
  
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
  
  logical function MAT_isUpdateDmg3()
  MAT_isUpdateDmg3 = update_in
  end function MAT_isUpdateDmg3

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
! ARG: sm0      [dble][0d0] initial confining mean stress [positive] 
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
! ARG: da_max   [dble] 1d-2, max damage per step 
! ARG: healLaw  [character(10)] 'EXP','LOG','NONE'
! ARG: dtmax    [dble] maximum time step limit for controlling time step
! ARG: update   [logic] update damage or not
! 
! Strain rate dependent Cd
! ARG: fix_cd  [logic] use fixe value cd
! ARG: edot0   [dble] reference strain rate
! ARG: Cdm     [dble] exponent
! ARG: Cdmin   [dble] 0.1
! ARG: Cdmax   [dble] 1d4
!
! log10(Cd/Cd0) = 1 + Cdm * log10(edot/edot0)
!
! END INPUT BLOCK
!
! SYNTAX : &MAT_DMG3 alpha|alphaH, Sm0|Sm0H, e0, ep 
!                    cp, cs, rho, mu_r, f0, c0, Cd, 
!                    C1, C2, Rp, Th, R , da_max, update
!                    edot0, Cd0, Cdm /
!
! Note that: alphaH, Sm0H can be spatial distributions
!

  subroutine MAT_DMG3_read(input,iin)

  use echo, only : echo_input, iout
  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho, cp, cs, f0, c0, alpha, Cd, C1, C2, dtmax
  double precision :: mu_r, Rp, Th, R, ep(2), e0(2), sm0,da_max
  character(20):: alphaH, sm0H
  logical :: update

  NAMELIST / MAT_DMG3 / alpha, sm0, e0, ep, alphaH, sm0H, mu_r, & 
                          cp, cs, rho, f0, c0, Cd, C1, C2, R, Rp, Th, &
                          da_max, healLaw, dtmax, update, edot0, Cdm, &
                          fix_cd, Cdmin, Cdmax, EqStart
  
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
  da_max= 1d-2
  dtmax = 1d10
  update = .true.

  ! read parameters
  read(iin, MAT_DMG3, END=100)
  Cd0   = Cd

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
  call MAT_setProp(input,'da_max', da_max)
  call MAT_setProp(input,'dtmax', dtmax)

  ! healing parameters
  ! Logarithmic healing parameters
  call MAT_setProp(input,'C1', C1)
  call MAT_setProp(input,'C2', C2)

  ! exponential healing parameters
  call MAT_setProp(input,'Rp', Rp)
  call MAT_setProp(input,'Th', Th)
  update_in = update

  write(iout,200) cp, cs, rho, f0, c0, &
                  sm0H, alphaH, Cd, mu_r, R, da_max, &
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
    'Maximum damage per step allowed .(da_max) =',EN12.3,/5x, &
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
  call MAT_setProp(elem,'da_max',ecoord,MAT_DMG3_mempro)
  call MAT_setProp(elem,'dtmax',ecoord,MAT_DMG3_mempro)

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
    allocate(m%alpha0(ngll, ngll))
    allocate(m%Cd(ngll, ngll))
   
    ! elastic coefficients for intact material (zero damage)
    call MAT_getProp(mu0,p,'mu')
    call MAT_getProp(lambda0,p,'lambda')
    m%mu0 = mu0
   
    ! constants for damage, healing, yielding
    call MAT_getProp(m%Cd, p, 'Cd')
    call MAT_getProp(m%C1, p, 'C1')
    call MAT_getProp(m%C2, p, 'C2')
    call MAT_getProp(m%mu_r,p,'mu_r')
    call MAT_getProp(m%da_max,p,'da_max')
   
    ! Cv = R/mu = coefficient for damage-related plasticity 
     call MAT_getProp(R,p,'R')
     m%Cv = R / m%mu0
   
    ! damage variable 
    call MAT_getProp(m%alpha, p, 'alpha')

    ! exp healing
    call MAT_getProp(m%Rp,p,'Rp')
    call MAT_getProp(m%Th,p,'Th')

    m%alpha0 = m%alpha*m%Rp
    
    m%healLaw = healLaw
   
    ! i2_cr = damage yield threshold for i2
    call Mat_getProp(sm0, p, 'sm0')
    call Mat_getProp(f0, p, 'f0')
    call Mat_getProp(c0, p, 'c0')
   
    m%i2_cr = compute_i2_cr(f0, mu0, c0, sm0) 
   
    ! initial strain
    call MAT_getProp(m%e0(1), p, 'e31_0')
    call MAT_getProp(m%e0(2), p, 'e32_0')
   
    ! initial plastic strain
    call MAT_getProp(m%ep(:,:,1), p, 'e31_p')
    call MAT_getProp(m%ep(:,:,2), p, 'e32_p')
   
    ! initial antiplane shear stress, uniform
    mu = m%mu0 - m%mu0 * m%mu_r * m%alpha(1,1)
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

end subroutine MAT_DMG3_init_elem_work

! compute the critical I2 beyond which damage starts
! allow sm0 be spatially variable
!-----------------------------------------------------------------------
  function compute_i2_cr(f0, mu0, c0, sm0) result(i2_cr)
      double precision, intent(in) :: f0, mu0, c0, sm0(:,:)
      double precision, dimension(size(sm0,1), size(sm0,2)) :: i2_cr

      i2_cr = 0.5d0* ((f0*sm0+c0)/mu0)**2d0 
      
  end function compute_i2_cr

!
! Compute the Cijkl with damage (scale shear modulus)
! ===========================================================
!
subroutine MAT_DMG3_Cijkl(Cijkl, isCijklZero, m, ndof, ngll)
   implicit none
   integer :: i, j, k , l, ndim, ndof, ngll
   double precision :: Cijkl(ndof, 2, ndof, 2, ngll, ngll), v
   logical :: isCijklZero(ndof, 2, ndof, 2)
   type (matwrk_dmg3_type), intent(in) :: m
   double precision :: mu(ngll, ngll)

   ndim = 2
   Cijkl = 0d0
   isCijklZero = .true.

   if (.not. ndof==1) call IO_abort('ndof must be 1 for mode 3 damage model') 

   mu = m%mu0 - m%mu0*m%mu_r*m%alpha

   Cijkl(1,1,1,1, :, :) = mu
   Cijkl(1,2,1,2, :, :) = mu
   isCijklZero(1,1,1,1) = .false.
   isCijklZero(1,2,1,2) = .false.

end subroutine MAT_DMG3_Cijkl

!=======================================================================
! Constitutive law:
! compute relative stresses from given strains
! and (if requested) update damage variable and plastic strain
!
! compute stress and update or not the damage variable
!

subroutine MAT_DMG3_stress(s, etot, m, ngll, update, dt)
  use echo, only : iout
  use utils, only: positive_part
  integer, intent(in) :: ngll
  double precision, intent(in)  :: etot(ngll,ngll,2)
  double precision, intent(out) :: s(ngll,ngll,2)
  type (matwrk_dmg3_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt
  double precision :: e(ngll,ngll,2), mu(ngll, ngll)
  double precision, dimension(ngll,ngll) :: i2
  double precision, dimension(2) :: dep
  double precision :: ap, da, t0, an,cd
  double precision :: edot_eff(ngll,ngll)
  integer :: i,j
  ! TESTING ONLY!!
  double precision :: e_cr, s_cr

  da     = 0d0

  e_cr = sqrt(maxval(m%i2_cr))
  s_cr = 2*m%mu0*e_cr 

 !!-- total strain
  e(:,:,1) = etot(:,:,1) + m%e0(1)
  e(:,:,2) = etot(:,:,2) + m%e0(2)
 !!-- elastic strain
  e  = e - m%ep

  i2 = e(:,:,1)**2d0 + e(:,:,2)**2d0 
  mu = m%mu0 - m%mu0*m%mu_r*m%alpha 

  ! total stress
  do i = 1, 2
    s(:, :, i)  = 2d0*mu*e(:,:,i)
  end do

!!------- Update damage variables and plastic strain --------------- 
  if (update .and. update_in) then

    if (.not.present(dt)) &
      call IO_abort('MAT_DMG3_stress: update requested but argument dt is absent')

    ! determine the maximum allowed time step 
    do i =1, ngll
    do j =1, ngll
      da = 0d0
      cd = m%Cd(i,j)
      ! make cd strain rate dependent

      !-- damage evolution
      if (i2(i,j)>m%i2_cr(i,j)) then
        ! damage increase
        da           = dt*cd*(i2(i,j) - m%i2_cr(i,j))
        m%alpha(i,j) = m%alpha(i,j) + da

        ! increase the permenent damage
        if (m%healLaw == 'EXP') then
            m%alpha0(i,j) = m%alpha0(i,j) + da*m%Rp
        end if

        if (m%alpha(i, j) > 1d0) then 
!            write(*, *) "Cd = ", cd, "da = ", da  
!            call IO_abort('Damage exceeds 1, simulation Stop!!')
            m%alpha(i,j) = 1d0
        end if
        !-- plasticity update
        ! Damage-related viscosity, if alpha_dot > 0
        ! Plastic strain rate = deij/dt = tij * Cv *dalpha/dt
        ! where tij = deviatoric stress
        da           = m%Cv * da
        dep(1)       = s(i,j,1)*da
        dep(2)       = s(i,j,2)*da
        m%ep(i,j,:)  = m%ep(i,j,:) + dep
      else
        ! damage healing, two laws: EXP or LOG
        select case (m%healLaw)
          case ('EXP')
            ! exponential healing
            ap = m%alpha0(i,j) ! permenent damage
            an = m%alpha(i,j) 
            m%alpha(i,j) = ap + (an-ap) * exp(-dt/m%Th) 
          case ('LOG')
            an = m%alpha(i,j)
            da = -m%C2*log(1 - m%C1/m%C2*exp(an/m%C2)*(i2(i,j)-m%i2_cr(i,j))*dt)
            m%alpha(i,j) = an + da 
            ! damage must be positive
            m%alpha(i,j) = max(m%alpha(i,j), 0d0)
          case ('NONE')
          ! do nothing, no healing
        end select
      end if ! damage evolve
    end do !j
    end do !i
  end if !update

 !-- relative stresses
  s(:,:,1) = s(:,:,1) - m%s0(1)
  s(:,:,2) = s(:,:,2) - m%s0(2)

  end subroutine MAT_DMG3_stress

! compute frictitious stress due to plasticity
! relative stress:
! sij = 2*mu*e_ij + s_ij_ep
!

subroutine MAT_DMG3_stress_ep(s, m, ngll)
  use echo, only : iout
  integer, intent(in) :: ngll
  double precision, intent(out) :: s(ngll,ngll,2)
  type (matwrk_dmg3_type), intent(inout) :: m
  double precision :: e(ngll,ngll,2), mu(ngll, ngll)
  integer :: i,j

  mu = m%mu0 - m%mu0*m%mu_r*m%alpha 
 !!-- elastic strain
  do i = 1, 2
    s(:,:,i)  = 2d0*mu*(m%e0(i) - m%ep(:,:,i)) - m%s0(i)
  end do

end subroutine MAT_DMG3_stress_ep

subroutine MAT_DMG3_Set_Cd(edot, m, ngll, ndof, isdynamic)
  use echo, only : iout
  integer, intent(in) :: ngll, ndof
  double precision, intent(in) :: edot(ngll, ngll, ndof + 1)
  type (matwrk_dmg3_type), intent(inout) :: m
  double precision:: edot_eff(ngll, ngll)
  logical, optional::isdynamic
  integer :: i, j

  ! hard coded to only damage in the dynamic steping
  if (.not. present(isdynamic)) then 
      isdynamic= .true.
  end if

  if (.not. isdynamic) then
      m%Cd = 0.d0 
      return 
  end if

  if (fix_cd) then
     m%Cd = Cd0 
     return
  end if

 !!-- elastic strain
  edot_eff = sqrt(0.5d0*(edot(:,:,1)**2d0 + edot(:,:,2)**2d0))

  do i = 1, ngll
      do j=1, ngll
          m%Cd(i,j) = compute_cd(edot_eff(i,j))
      end do
  end do

end subroutine MAT_DMG3_Set_Cd

function compute_cd(edot) result(cd_out)
    double precision::edot, cd_out
    ! Cd0, Cdm, edot0 are saved constants for the module
    if (abs(edot)<1.0d-12) edot= 1.0d-12
    cd_out = Cd0 * 10.0**(1d0 + Cdm * log10(edot/edot0))
    cd_out = max(cd_out, Cdmin)
    cd_out = min(cd_out, Cdmax)
end function compute_cd

subroutine MAT_DMG3_dtmax(s, etot, m, ngll, dtmax)
  integer, intent(in) :: ngll
  double precision, intent(in)  :: etot(ngll,ngll,2)
  double precision, intent(out) :: s(ngll,ngll,2)
  type (matwrk_dmg3_type), intent(inout) :: m
  double precision :: dtmax
  double precision :: e(ngll,ngll,2), mu(ngll, ngll)
  double precision, dimension(ngll,ngll) :: i2 
  double precision :: dtmax_ij, di2, an, tcr, ap, da
  integer :: i,j

  dtmax=huge(0d0)

  if (.not. update_in) return

 !!-- total strain
  e(:,:,1) = etot(:,:,1) + m%e0(1)
  e(:,:,2) = etot(:,:,2) + m%e0(2)

 !!-- elastic strain
  e  = e - m%ep

  i2 = e(:,:,1)**2d0 + e(:,:,2)**2d0 

  mu = m%mu0 - m%mu0*m%mu_r*m%alpha 

  ! total stress
  do i = 1, 2
    s(:,:,i)  = 2d0*mu*e(:,:,i)
  end do
  
  dtmax_ij = 0d0

  do i =1, ngll
    do j =1, ngll
   !-- damage evolution
   di2 = i2(i,j) - m%i2_cr(i,j)
    if (di2>0) then
      ! damage increase
      dtmax_ij    = min(m%da_max/(m%Cd(i,j)*di2), m%da_max/(m%Cd(i,j)*di2)/m%R)
    else
      ! maximum time step allowed from healing 
      select case(m%healLaw)
      case ('EXP')
         ap = m%alpha0(i, j)
         an = m%alpha(i, j)
         da = an-ap
         if (da<=m%da_max) then 
             dtmax_ij = huge(0d0)
         else 
             dtmax_ij = abs(log(da-m%da_max/da)*m%Th) 
         end if
      case ('LOG')
         an  = m%alpha(i, j) 
         tcr = m%C2/m%C1/(-di2)*exp(-an/m%C2)
         dtmax_ij = (exp(m%da_max/m%C2)-1)*tcr 
         if (an<=m%da_max) dtmax_ij=huge(0d0)
      case ('NONE')
      end select
    end if ! damage evolve
      dtmax       = min(dtmax, dtmax_ij)
    end do !j
  end do !i

  dtmax = min(dtmax, m%dtmax)

  end subroutine MAT_DMG3_dtmax

!=======================================================================
! export output data

  function MAT_DMG3_export(m) result(dat)

      type(matwrk_dmg3_type), intent(in) :: m
      real :: dat(size(m%alpha,1),size(m%alpha,2),3)  ! 4 = 1(alpha)+size(ep,3)

      dat(:,:,1) = real(m%alpha)
      dat(:,:,2:3) = real(m%ep)
  
  end function MAT_DMG3_export

end module mat_dmg3
