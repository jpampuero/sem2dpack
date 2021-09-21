module mat_plastic
! plane strain Coulomb plasticity following Andrews (2005)

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
 !-- elasto-visco-plasticity
  type matwrk_plast_type
    private
!    double precision, pointer, dimension(:,:) :: mu => null(), lambda => null()
    double precision :: mu, lambda
    double precision, pointer, dimension(:,:,:) :: ep => null()
    double precision, pointer, dimension(:) :: e0=>null(), s0=>null()
    double precision :: yield_mu=0d0, yield_co=0d0, vp_factor=1d0
  end type matwrk_plast_type

  integer, save :: isPlastic = 0

  ! for memory report
  integer, save :: MAT_PLAST_mempro = 0
  integer, save :: MAT_PLAST_memwrk = 0

  public :: matwrk_plast_type &
          , MAT_isPlastic, MAT_anyPlastic, MAT_PLAST_read, MAT_PLAST_init_elem_prop &
          , MAT_PLAST_init_elem_work, MAT_PLAST_stress, MAT_PLAST_export &
          , MAT_PLAST_mempro, MAT_PLAST_memwrk
  
contains

!=======================================================================
  logical function MAT_isPlastic(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isPlastic = MAT_isKind(m,isPlastic)
  end function MAT_isPlastic

!=======================================================================
  logical function MAT_anyPlastic()
  MAT_anyPlastic = (isPlastic>0)
  end function MAT_anyPlastic

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_PLASTIC
! GROUP  : MATERIALS
! PURPOSE: Set material properties for elasto-plastic material
!          with Mohr-Coulomb yield criterion
!          and non-dilatant (null volumetric plastic strain)
! SYNTAX : &MAT_PLASTIC cp,cs,rho,phi,coh,Tv,e0 /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: phi      [dble][0d0] internal friction angle
! ARG: coh      [dble][0d0] cohesion
! ARG: Tv       [dble][0d0] visco-plastic relaxation time
! ARG: e0       [dble(3)][0d0] initial total strain (11, 22 and 12)
!
! END INPUT BLOCK

  subroutine MAT_PLAST_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho,cp,cs,phi,coh,Tv,e0(3)

  NAMELIST / MAT_PLASTIC / cp,cs,rho,phi,coh,Tv,e0
  
  cp  = 0d0
  cs  = 0d0
  rho = 0d0
  e0  = 0d0
  phi = 0d0
  coh = 0d0
  Tv  = 0d0

  read(iin, MAT_PLASTIC, END=100)

  write(iout,200) cp,cs,rho,phi,coh,Tv,e0

  call MAT_setKind(input,isPlastic)
  
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'phi',phi)
  call MAT_setProp(input,'coh',coh)
  call MAT_setProp(input,'Tv',Tv)

  call MAT_setProp(input,'e11_0',e0(1))
  call MAT_setProp(input,'e22_0',e0(2))
  call MAT_setProp(input,'e12_0',e0(3))


  return

  100 call IO_abort('MAT_PLAST_read: MAT_PLASTIC input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Internal friction angle . . . . . . (phi) =',EN12.3,/5x, &
    'Cohesion. . . . . . . . . . . . . . (coh) =',EN12.3,/5x, &
    'Visco-plastic timescale . . . . . . .(Tv) =',EN12.3,/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x, &
    '                     12 . . . . . (e0(3)) =',EN12.3)

  end subroutine MAT_PLAST_read

!=======================================================================
  subroutine MAT_PLAST_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision :: rho,cp,cs

  call MAT_setProp(elem,'cp',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'cs',ecoord,MAT_PLAST_mempro)

  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_PLAST_mempro)
  call MAT_setProp(elem,'mu',rho*cs*cs,MAT_PLAST_mempro)

  call MAT_setProp(elem,'phi',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'coh',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'Tv',ecoord,MAT_PLAST_mempro)

  call MAT_setProp(elem,'e11_0',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e22_0',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e12_0',ecoord,MAT_PLAST_mempro)

  end subroutine MAT_PLAST_init_elem_prop

!=======================================================================
  subroutine MAT_PLAST_init_elem_work(m,p,ngll,dt)

  use constants, only : PI

  type(matwrk_plast_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  double precision, intent(in) :: dt

  double precision :: phi,Tv,lambda,two_mu

 ! elastic coefficients
!  allocate(m%lambda(n,n)) ! NOTE: uniform elastic moduli 
!  allocate(m%mu(n,n))     ! NOTE: uniform elastic moduli
  call MAT_getProp(m%lambda,p,'lambda')
  call MAT_getProp(m%mu,p,'mu')

 ! Coulomb yield function Y = max_shear_stress - cohesion*cos(phi) + mean_stress*sin(phi)
 ! where phi = internal friction angle
 ! The yield condition Y=0 is equivalent to  
 ! max_over_all_planes[ shear_stress + tan(phi)*normal_stress ] = cohesion
  call MAT_getProp(phi,p,'phi')
  phi = PI/180d0 *phi
  call MAT_getProp(m%yield_co,p,'coh')
  m%yield_co = m%yield_co * cos(phi)  ! store cohesion*cos(phi)
  m%yield_mu = sin(phi)

 ! visco-plastic coefficient
  call MAT_getProp(Tv,p,'Tv')
  if (Tv>0d0) then 
    m%vp_factor = 1d0-exp(-dt/Tv)
  else
    m%vp_factor = 1d0
  endif

 ! initial strain
  allocate(m%e0(3))
  call MAT_getProp(m%e0(1),p,'e11_0')
  call MAT_getProp(m%e0(2),p,'e22_0')
  call MAT_getProp(m%e0(3),p,'e12_0')

 ! initial plastic strain = 0
  allocate(m%ep(ngll,ngll,3))
  m%ep = 0d0

 ! initial stress
  allocate(m%s0(3))
  lambda = m%lambda
  two_mu = 2d0*m%mu
  m%s0(1) = (lambda+two_mu)*m%e0(1) + lambda*m%e0(2)
  m%s0(2) = lambda*m%e0(1) + (lambda+two_mu)*m%e0(2)
  m%s0(3) = two_mu*m%e0(3)

  MAT_PLAST_memwrk = MAT_PLAST_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%ep) &
                   + size(m%e0) &
                   + size(m%s0)
!                   + size(m%lambda) &
!                   + size(m%mu) &

  end subroutine MAT_PLAST_init_elem_work


!=======================================================================
! Constitutive law:
! compute relative stresses from given strains
! and (if requested) update plastic strain
  subroutine MAT_PLAST_stress(s,etot,m,ngll,update,E_ep,E_el,sg)

  use constants, only : COMPUTE_ENERGIES, COMPUTE_STRESS_GLUT

  integer, intent(in) :: ngll
  double precision, intent(in) :: etot(ngll,ngll,3)
  double precision, intent(out) :: s(ngll,ngll,3)
  type (matwrk_plast_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, intent(out) :: E_ep(ngll,ngll)
  double precision, intent(out) :: E_el(ngll,ngll)
  double precision, intent(out) :: sg(ngll,ngll,3)

  double precision :: lambda, two_mu, i1_0, i2_0
  double precision, dimension(ngll,ngll,3) :: e,sd,sd_tr,dep
  double precision, dimension(ngll,ngll) :: sm,tau,Y,i1,i2 ,factor

 ! relative elastic strain
  e = etot - m%ep

  lambda = m%lambda
  two_mu = 2d0*m%mu

 ! if requested to update internal variables
  if (update) then

   ! trial absolute elastic strain (assumes the whole strain increment is elastic)
    e(:,:,1) = e(:,:,1) + m%e0(1)
    e(:,:,2) = e(:,:,2) + m%e0(2)
    e(:,:,3) = e(:,:,3) + m%e0(3)
  
   ! trial stress (assumes elastic increment)
    s(:,:,1) = (lambda+two_mu)*e(:,:,1) + lambda*e(:,:,2)
    s(:,:,2) = lambda*e(:,:,1) + (lambda+two_mu)*e(:,:,2)
    s(:,:,3) = two_mu*e(:,:,3)
  
   ! maximum shear stress
    tau = sqrt( 0.25d0*(s(:,:,1)-s(:,:,2))**2 +s(:,:,3)**2 )
  
   ! mean stress
    sm = 0.5d0*( s(:,:,1) + s(:,:,2) )
  
   ! Coulomb yield stress
    Y = m%yield_co - m%yield_mu *sm 
  
   ! visco-plastic update of deviatoric stresses
   ! see Andrews (2005), corrected by Duan and Day (2008)
   ! This is equivalent to 
   !   stress_rate = c:total_strain_rate - (stress-yield)/Tv
   ! (Simo and Hughes, "Computational Inelasticity", 1998, eq 1.7.15)
   ! Trial (elastic) deviatoric stresses
    sd_tr(:,:,1)  = s(:,:,1) - sm
    sd_tr(:,:,2)  = s(:,:,2) - sm
    sd_tr(:,:,3)  = s(:,:,3)
   ! Here vp_factor = 1-exp(-dt/Tv)
    factor = 1d0 - max( 1d0-Y/tau, 0d0) * m%vp_factor
   ! if Tv=0 (no viscosity) then factor=min(Y/tau,1)
    sd(:,:,1) = factor * sd_tr(:,:,1)
    sd(:,:,2) = factor * sd_tr(:,:,2)
    sd(:,:,3) = factor * sd_tr(:,:,3)
  
   ! update plastic strain
    dep = (sd_tr - sd)/two_mu
    m%ep = m%ep + dep
  
   ! recompose stress (mean + deviatoric)
    s(:,:,1) = sd(:,:,1) + sm
    s(:,:,2) = sd(:,:,2) + sm
    s(:,:,3) = sd(:,:,3)
  
    if (COMPUTE_ENERGIES) then
     ! increment of plastic energy dissipation
      E_ep = s(:,:,1)*dep(:,:,1) + s(:,:,2)*dep(:,:,2) + 2d0*s(:,:,3)*dep(:,:,3)
     ! total elastic energy change
      i1_0 = m%e0(1) + m%e0(2)
      i2_0 = m%e0(1)*m%e0(1) +m%e0(2)*m%e0(2) +2d0*m%e0(3)*m%e0(3)
      e = e - dep ! update absolute elastic strain
      i1 = e(:,:,1) + e(:,:,2)
      i2 = e(:,:,1)*e(:,:,1) +e(:,:,2)*e(:,:,2) +2d0*e(:,:,3)*e(:,:,3)
      E_el = 0.5d0* ( lambda*( i1*i1-i1_0*i1_0 ) + two_mu*(i2-i2_0) )
    else
      E_ep = 0d0
      E_el = 0d0
    endif
  
    if (COMPUTE_STRESS_GLUT) then
      ! no lambda terms because ep is deviatoric
      sg(:,:,1) = - two_mu*m%ep(:,:,1)  
      sg(:,:,2) = - two_mu*m%ep(:,:,2)
      sg(:,:,3) = - two_mu*m%ep(:,:,3)
    else
      sg = 0d0
    endif
  
   ! relative stress
    s(:,:,1) = s(:,:,1) - m%s0(1)
    s(:,:,2) = s(:,:,2) - m%s0(2)
    s(:,:,3) = s(:,:,3) - m%s0(3)

  else ! if no update, just compute relative stress
    s(:,:,1) = (lambda+two_mu)*e(:,:,1) + lambda*e(:,:,2)
    s(:,:,2) = lambda*e(:,:,1) + (lambda+two_mu)*e(:,:,2)
    s(:,:,3) = two_mu*e(:,:,3)
  endif

  end subroutine MAT_PLAST_stress

!=======================================================================
! export output data

  function MAT_PLAST_export(m) result(dat)

  type(matwrk_plast_type), intent(in) :: m
  real :: dat(size(m%ep,1),size(m%ep,2),size(m%ep,3))

  dat = real(m%ep)
  
  end function MAT_PLAST_export


end module mat_plastic
