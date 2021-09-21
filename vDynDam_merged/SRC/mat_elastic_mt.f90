!=======================================================================
!***********************   ELASTICITY   ****************************
!This module accounts for Linear Elasticity in an homogeneous isotropic medium
!				    Written by Marion THOMAS
!					Last Modified May 2017
!=======================================================================
module mat_elastic_mt

!Modules to be used
  use prop_mat
  use stdio, only : IO_abort
  use utils_dyn_damage

  implicit none
  private

!-----------------------------------------------------------------------
! Define a new derived type corresponding to the material: dynamic damage
!-----------------------------------------------------------------------
  type matwrk_elasM_type
    private
    double precision :: mu, lambda
    double precision, pointer, dimension(:,:,:) :: e => null()
    double precision, pointer, dimension(:,:,:) :: s => null()
    double precision, pointer, dimension(:) :: e0=>null()
    double precision, pointer, dimension(:) :: s0=>null()
    double precision, pointer, dimension(:,:) :: alpha=>null(),erate=>null()
    double precision, pointer, dimension(:,:) :: invI=>null(),invII=>null()
  end type matwrk_elasM_type

  ! Defined in the first call to MAT_setKind
  integer, save :: isElasticM = 0
  integer, save :: isAntiplane = 0 ! needed for optimizations

  ! for memory report
  integer, save :: MAT_ELASM_mempro = 0
  integer, save :: MAT_ELASM_memwrk = 0

  public :: matwrk_elasM_type &
          , MAT_isElasticM, MAT_isAntiplane, MAT_anyElasticM, MAT_ELASM_read &
          , MAT_ELASM_init_elem_prop, MAT_ELASM_init_elem_work &
          , MAT_ELASM_stress, MAT_ELASM_export &
          , MAT_ELASM_memwrk, MAT_ELASM_mempro

contains

!=======================================================================
! Tell if an element is of material "ElasticM" and Anti-Plane Strain
!=======================================================================
  logical function MAT_isElasticM(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isElasticM = MAT_isKind(m,isElasticM)
  end function MAT_isElasticM
!-----------------------------------------------------------------------
  logical function MAT_isAntiplane(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isAntiplane = MAT_isKind(m,isAntiplane)
  end function MAT_isAntiplane
!-----------------------------------------------------------------------
  logical function MAT_anyElasticM()
  MAT_anyElasticM = (isElasticM>0)
  end function MAT_anyElasticM

!=======================================================================
! Read material properties form main input file
!=======================================================================
subroutine MAT_ELASM_read(input,iin,ndof)
!-----------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : MAT_ELASTICM
! GROUP  : MATERIALS
! PURPOSE: Set material properties for a linear elastic medium
! SYNTAX : For isotropic material:
!           &MAT_ELASTIC rho, cp, cs, alpha, e0 /
!  
! ARG: cp       [dble][0d0] P wave velocity (m/s)
! ARG: cs       [dble][0d0] S wave velocity (m/s)
! ARG: rho      [dble][0d0] density (kg/m^3)
! ARG: alpha    [dble][0d0] state variable
! ARG: e0       [dble(3)][0d0] initial total strain if plane strain (11, 22 and 12)
! ARG: e0       [dble(2)][0d0] initial total strain if antiplane (13 and 23)

! END INPUT BLOCK
!-----------------------------------------------------------------------

  !modules
  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  !input/ouput variables
  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin, ndof

  !local variables
  double precision :: Kmod,poiss,rho,alpha,cp,cs,mu,lambda,e0(3)
  double precision :: invI,invII

  NAMELIST / MAT_ELASTICM / rho,cp,cs,e0,alpha

  !-----------------------------------------------------------------------  
  !initial value
  !-----------------------------------------------------------------------  
  cp    = 0d0
  cs    = 0d0
  rho   = 0d0
  e0    = 0d0
  alpha = 0d0
  invI=0d0
  invII=0d0
  !-----------------------------------------------------------------------  
  !Get the values from the input file (par.inp)
  !-----------------------------------------------------------------------  
  read(iin, MAT_ELASTICM, END=100)
  call MAT_setKind(input,isElasticM)

  !Isotropic and Homegeneous material: get cs, cp and rho form the input file (par.inp)
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'alpha',alpha)
  call MAT_setProp(input,'invI',invI)
  call MAT_setProp(input,'invII',invII)
  !-----------------------------------------------------------------------   
  !Define others variables: mu, Lambda, poisson ratio, Kmod
  !-----------------------------------------------------------------------  
  if (cp>0d0 .and. cs>0d0 .and. rho>0d0) then
    mu   = rho*cs*cs
    Kmod  = rho*cp*cp
    lambda  = rho*(cp*cp - 2d0*cs*cs)
    poiss = 0.5d0*lambda/(Kmod-mu)  !(cp*cp-2d0*cs*cs)/(cp*cp-cs*cs)
    if (poiss < 0.d0 .or. poiss > 0.5d0) call IO_abort('Poisson''s ratio out of range !')
    if (echo_input) write(iout,200) cp,cs,rho,alpha,poiss,lambda,mu, &
                   lambda+2d0*mu/3d0, 2d0*mu*(1d0+poiss)  ! Kvol,young  
  !-----------------------------------------------------------------------  
  !assign initial value
  !-----------------------------------------------------------------------  
   call MAT_setProp(input,'lambda',lambda)
    call MAT_setProp(input,'mu',mu)
    if (ndof==1) then !SH antiplane
        call MAT_setKind(input,isAntiplane) 
        call MAT_setProp(input,'e13_0',e0(1))
        call MAT_setProp(input,'e23_0',e0(2))
        if (echo_input) write(iout,250) e0
    else !P-SV Plane strain/ inplane
        call MAT_setProp(input,'e11_0',e0(1))
        call MAT_setProp(input,'e22_0',e0(2))
        call MAT_setProp(input,'e12_0',e0(3))
        if (echo_input) write(iout,300) e0
    endif  
  else
    call IO_abort('MAT_ELASTM_read: incomplete input')
  endif

  return
  100 call IO_abort('MAT_ELASM_read: MAT_ELASTICM input block not found')
  !-----------------------------------------------------------------------  
  ! simulation update info
  !-----------------------------------------------------------------------  
  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'State variable. . . . . . . . . . (alpha) =',EN12.3,/5x, &
    'Poisson''s ratio . . . . . . . . . . . . .=',EN12.3,/5x, &
    'First Lame parameter Lambda . . . . . . . =',EN12.3,/5x, &
    'Second Lame parameter Mu. . . . . . . . . =',EN12.3,/5x, &
    'Bulk modulus K. . . . . . . . . . . . . . =',EN12.3,/5x, &
    'Young''s modulus E . . . . . . . . . . . .=',EN12.3)

  250   format(5x, &
    '                                            antiplane',/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x) 

  300   format(5x, &
    '                                            plane strain',/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x, &
    '                     12 . . . . . (e0(3)) =',EN12.3)

end subroutine MAT_ELASM_read

!=======================================================================
! Initialise material properties for one element
!=======================================================================
  subroutine MAT_ELASM_init_elem_prop(elem,ecoord)

  !input/ouput variables
  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  !-----------------------------------------------------------------------  
  ! assign value for each element
  !-----------------------------------------------------------------------  
  call MAT_setProp(elem,'cp',ecoord,MAT_ELASM_mempro)
  call MAT_setProp(elem,'cs',ecoord,MAT_ELASM_mempro)   
  call MAT_setProp(elem,'alpha',ecoord,MAT_ELASM_mempro)
  call MAT_setProp(elem,'lambda',ecoord,MAT_ELASM_mempro)
  call MAT_setProp(elem,'mu',ecoord,MAT_ELASM_mempro)
  call MAT_setProp(elem,'invI',ecoord,MAT_ELASM_mempro)   
  call MAT_setProp(elem,'invII',ecoord,MAT_ELASM_mempro)   

  !SH antiplane
  if (MAT_isKind(elem,isAntiplane)) then 
    call MAT_setProp(elem,'e13_0',ecoord,MAT_ELASM_mempro)
    call MAT_setProp(elem,'e23_0',ecoord,MAT_ELASM_mempro)
  !P-SV Plane strain/ inplane
  else 
    call MAT_setProp(elem,'e11_0',ecoord,MAT_ELASM_mempro)
    call MAT_setProp(elem,'e22_0',ecoord,MAT_ELASM_mempro)
    call MAT_setProp(elem,'e12_0',ecoord,MAT_ELASM_mempro)
  endif  

  end subroutine MAT_ELASM_init_elem_prop

!=======================================================================
! Assign material properties for one element
!=======================================================================
 subroutine MAT_ELASM_init_elem_work(m,p,n,ndof)

  !input/ouput variables
  type(matwrk_elasM_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: n, ndof

  !local variables
  double precision :: lambda,two_mu
  
 !-----------------------------------------------------------------------  
 ! elastic coefficients for homogeneous istropic material
 !-----------------------------------------------------------------------  
 call MAT_getProp(m%lambda,p,'lambda')
 call MAT_getProp(m%mu,p,'mu')
 !-----------------------------------------------------------------------  
 ! State variable
 !-----------------------------------------------------------------------  
  allocate(m%alpha(n,n))
  call MAT_getProp(m%alpha,p,'alpha')
  allocate(m%invI(n,n))
  call MAT_getProp(m%invI,p,'invI')
  allocate(m%invII(n,n))
  call MAT_getProp(m%invII,p,'invII')
 !-----------------------------------------------------------------------  
 ! initial strain
 !-----------------------------------------------------------------------
  allocate(m%e(n,n,ndof+1))
  allocate(m%erate(n,n))
  allocate(m%e0(ndof+1))
  if (ndof==1) then !SH antiplane
    call MAT_getProp(m%e(:,:,1),p,'e13_0')
    call MAT_getProp(m%e(:,:,2),p,'e23_0')
    call MAT_getProp(m%e0(1),p,'e13_0')
    call MAT_getProp(m%e0(2),p,'e23_0')
  else !P-SV plane strain
    call MAT_getProp(m%e(:,:,1),p,'e11_0')
    call MAT_getProp(m%e(:,:,2),p,'e22_0')
    call MAT_getProp(m%e(:,:,3),p,'e12_0')
    call MAT_getProp(m%e0(1),p,'e11_0')
    call MAT_getProp(m%e0(2),p,'e22_0')
    call MAT_getProp(m%e0(3),p,'e12_0')
 endif
 m%erate(:,:) = 0.0d0
 !-----------------------------------------------------------------------  
 ! initial stress
 !-----------------------------------------------------------------------  
  allocate(m%s0(3))
  allocate(m%s(n,n,ndof+1))
  lambda = m%lambda
  two_mu = 2d0*m%mu
  if (ndof==1) then !SH antiplane
      m%s0(1) = two_mu*m%e0(1)! sigma_13
      m%s0(2) = two_mu*m%e0(2)! sigma_23
      m%s(:,:,1) = m%s0(1)
      m%s(:,:,2) = m%s0(2)
   else !PSV
      m%s0(1) = (lambda+two_mu)*m%e0(1) + lambda*m%e0(2)! sigma_11
      m%s0(2) = lambda*m%e0(1) + (lambda+two_mu)*m%e0(2)! sigma_22
      m%s0(3) = two_mu*m%e0(3)! sigma_12
      m%s(:,:,1) = m%s0(1)
      m%s(:,:,2) = m%s0(2)
      m%s(:,:,3) = m%s0(3)
  endif
  !-----------------------------------------------------------------------
  ! assign memory
  !-----------------------------------------------------------------------
  MAT_ELASM_memwrk = MAT_ELASM_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%alpha) &
 				   + size(m%invI) &
				   + size(m%invII) &
                   + size(m%s) &
                   + size(m%e) &
                   + size(m%erate) &                   
                   + size(m%e0) &
                   + size(m%s0)
          
end subroutine MAT_ELASM_init_elem_work
  
!=======================================================================  
! MAR_USER_actions: Action to perform during the solver phase
! Constitutive law: compute relative stresses from given strains
!=======================================================================
  subroutine MAT_ELASM_stress(s,etot,m,ngll,ndof,update,dt)

  !input/ouput variables
  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: etot(ngll,ngll,ndof+1)
  double precision, intent(out) :: s(ngll,ngll,ndof+1)
  type (matwrk_elasM_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt

  !local variables
  double precision :: e(ngll,ngll,ndof+1)
  double precision :: erate(ngll,ngll,ndof+1)
  double precision, dimension(ngll,ngll) :: dalpha
  double precision, dimension(ngll,ngll) :: invItemp,invIItemp
  double precision, dimension(ndof+1) :: eij(ndof+1),sij(ndof+1)
  double precision :: alphaij
  double precision :: stress(3,3)
  double precision :: sigma,tau
  integer :: i,j
  double precision :: test1,epsi,gamma,nu,strain(3,3)
  
 !-----------------------------------------------------------------------  
 ! elastic coefficients
 !-----------------------------------------------------------------------  
  double precision :: mu,lambda
  lambda = m%lambda
  mu = m%mu
 !-----------------------------------------------------------------------  
 ! Update strain
 !-----------------------------------------------------------------------  
  e(:,:,1) = etot(:,:,1) + m%e0(1)
  e(:,:,2) = etot(:,:,2) + m%e0(2)
  if (ndof==2) then !PSV 
    e(:,:,3) = etot(:,:,3) + m%e0(3)
  endif  
 !-----------------------------------------------------------------------
 !	Get Saved Variables from previous time-step
 !-----------------------------------------------------------------------
  invItemp(:,:)= m%invI
  invIItemp(:,:)= m%invII
 !-----------------------------------------------------------------------  
 ! compute stresses and two strain invariants
 !-----------------------------------------------------------------------  
  do j=1,ngll
  	do i=1,ngll
    eij = e(i,j,:)

    !antiplane
    if (ndof==1) then !SH
      sij(1) = 2d0*mu*eij(1) ! sigma_13
      sij(2) = 2d0*mu*eij(2) ! sigma_23

    !plane strain
    else !PSV 
      sij(1) = (lambda+2d0*mu)*eij(1) + lambda*eij(2) ! sigma_11
      sij(2) = lambda*eij(1) + (lambda+2d0*mu)*eij(2) ! sigma_22
      sij(3) = 2d0*mu*eij(3)  ! sigma_12
    endif  

    s(i,j,:) = sij
	call vec2mat(stress,sij,ndof)
	call stress_invariant(stress,sigma,tau)
	invIItemp(i,j) = tau
	invItemp(i,j) = sigma				
  	enddo
  enddo
  !-----------------------------------------------------------------------
  !	Compute the strain rate
  !-----------------------------------------------------------------------
 ! if requested to update internal variables
  if (update) then

    !update strain
	m%erate(:,:)=(m%e(:,:,1)-e(:,:,1))/dt

    !-- compute the state variables
    do j=1,ngll
      do i=1,ngll
        eij = e(i,j,:)
        call state_variable(eij,alphaij,ndof)
        dalpha(i,j) = alphaij
      enddo
    enddo
  !-----------------------------------------------------------------------
  !	Update variables to for next times step and output requests
  !-----------------------------------------------------------------------
   m%alpha =  m%erate(:,:)
   m%invI=invItemp
   m%invII=invIItemp
  !update strain
   m%e = e

  endif
  !-----------------------------------------------------------------------
  ! relative stresses
  !-----------------------------------------------------------------------
  s(:,:,1) = s(:,:,1) - m%s0(1)
  s(:,:,2) = s(:,:,2) - m%s0(2)
  if (ndof==2) then !
      s(:,:,3) = s(:,:,3) - m%s0(3)
  endif
  m%s = s

  end subroutine MAT_ELASM_stress
  
!=======================================================================
! state variable evolution
!   alpha = e11+e22+e12
!=======================================================================
  subroutine state_variable(e,alpha,ndof)

  !input/ouput variables
  integer, intent(in) :: ndof
  double precision, intent(in) :: e(ndof+1)
  double precision, intent(out) :: alpha
  
  !-----------------------------------------------------------------------
  ! compute alpha
  !-----------------------------------------------------------------------
  if (ndof==1) then !SH
    alpha = e(1)+e(2)    
  else !PSV
    alpha = e(1)+e(2)+e(3)
  endif  

  end subroutine state_variable

!=======================================================================
! Extra data to be written on output binary file
! export output data
!=======================================================================
  function MAT_ELASM_export(m,ndof) result(dat)

  !input/ouput variables
  integer, intent(in) :: ndof
  type(matwrk_elasM_type), intent(in) :: m
  real :: dat(size(m%e,1),size(m%e,2),3)

  !-----------------------------------------------------------------------
  !assign values
  !-----------------------------------------------------------------------
  dat(:,:,1) = real(m%alpha)
  dat(:,:,2) = real(m%invI)
  dat(:,:,3) = real(m%invII)

  end function MAT_ELASM_export

end module mat_elastic_mt
