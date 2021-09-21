!=======================================================================
!***********************   DYNAMIC DAMAGE   ****************************
!				Written by Marion THOMAS & Harsha S. BHAT
!					Based on Bhat et al. (2012)
!					Last Modified May 2017
!=======================================================================
module mat_dyn_damage

 !Modules to be used
  use prop_mat
  use stdio, only : IO_abort
  use constants, only : PI
  use utils_dyn_damage
  
  implicit none
  private
  
!-----------------------------------------------------------------------
! Define a new derived type corresponding to the material: dynamic damage
!-----------------------------------------------------------------------
  type matwrk_dyndmg_type
    private
    double precision :: lambda
    double precision :: fs,nu,mu,rho,cr,Nv,beta,omega,KicSS,a0,alpha,D0,vm,fitKicD
    double precision, pointer, dimension(:,:,:) :: e => null()
    double precision, pointer, dimension(:,:,:) :: erate => null()
    double precision, pointer, dimension(:) :: e0=>null(), s0=>null()
    double precision, pointer, dimension(:,:) :: D=>null()
    double precision, pointer, dimension(:,:) :: css=>null()
    double precision, pointer, dimension(:,:) :: cps=>null()
    double precision, pointer, dimension(:,:) :: v=>null()
    double precision, pointer, dimension(:,:) :: KI=>null()
    double precision, pointer, dimension(:,:) :: R=>null()
    double precision, pointer, dimension(:,:) :: Rcum=>null()
!    double precision, pointer, dimension(:,:) :: A=>null(),B=>null(),C=>null()
    double precision, pointer, dimension(:,:) :: invI=>null(),invII=>null()
  end type matwrk_dyndmg_type
  
  !constants
  double precision, private :: cD0,calpha,cvm
  double precision, private :: cfs,cnu,cmu,crho,ccr,cfitKicD
  double precision, private :: cNv,comega,cbeta,cKicSS,ca0
  integer :: nn

  ! Variables to save for the RKF45 method
  double precision :: str_rate(3,3), strain0(3,3)
  double precision :: KIsav,start_time,end_time,Rsav

  ! Defined in the first call to MAT_setKind
  integer, save :: isDynDamage = 0
  integer, save :: isAntiplane = 0 ! needed for optimizations

  ! for memory report
  integer, save :: MAT_DyDMG_mempro = 0
  integer, save :: MAT_DyDMG_memwrk = 0

  public :: matwrk_dyndmg_type &
          , MAT_isDynDamage, MAT_isAntiplane, MAT_anyDynDamage, MAT_DyDMG_read &
          , MAT_DyDMG_init_elem_prop, MAT_DyDMG_init_elem_work &
          , MAT_DyDMG_stress, MAT_DyDMG_export &
          , MAT_DyDMG_memwrk, MAT_DyDMG_mempro

contains

!=======================================================================
! Tell if element is of material type "DynDamage" and Anti-Plane Strain
!=======================================================================
  logical function MAT_isDynDamage(m)
	  type(matpro_elem_type), intent(in) :: m
	  MAT_isDynDamage = MAT_isKind(m,isDynDamage)
  end function MAT_isDynDamage
!-----------------------------------------------------------------------
  logical function MAT_isAntiplane(m)
	  type(matpro_elem_type), intent(in) :: m
	  MAT_isAntiplane = MAT_isKind(m,isAntiplane)
  end function MAT_isAntiplane
!-----------------------------------------------------------------------
  logical function MAT_anyDynDamage()
 	 MAT_anyDynDamage = (isDynDamage>0)
  end function MAT_anyDynDamage
  
!=======================================================================
! Read material properties from main input file
!=======================================================================
subroutine MAT_DyDMG_read(input,iin,ndof)
!-----------------------------------------------------------------------
	! fs        static friction coefficient
	! rho       density (kg/m^3)
	! vm        Branching speed (m/s)
	! cs        S wave velocity (m/s)
	! cp        P wave velocity (m/s)
	! Nv        volume density of cracks
	! beta      Ashby and Sammis(1990) factor (usually 0.1)
	! omega     Cracks factor (usually 2.0)
	! KicSS     quasi-static fracture toughness
	! a0        initial radius of cracks
	! fitKicD   to compute the dynamic initation toughness
	! e0        initial total strain if plane strain (11, 22 and 12)
	! e0        initial total strain if antiplane (13 and 23)
!-----------------------------------------------------------------------

	  !modules
	  use echo, only : echo_input, iout
	  use stdio, only : IO_abort

	  !input/ouput variables
	  type (matpro_input_type), intent(inout) :: input
	  integer, intent(in) :: iin, ndof

	  !local variables
	  double precision :: fs,nu,mu,rho,lambda,cr,cp,cs,Nv,beta,omega,KicSS,a0,fitKicD,e0(3)
	  double precision :: KI,D,D0,alpha,phi,sigma,tau,R,vm,css,cps,v
	  double precision :: invI,invII,Rcum
!	  double precision :: A,B,C

	  NAMELIST / MAT_DYNDAMAGE / fs,nu,mu,rho,vm,cp,cs,Nv,beta,omega,KicSS,a0,e0,fitKicD

 	  !-----------------------------------------------------------------------  
	  !initial value
 	  !-----------------------------------------------------------------------  
	  fs      = 0d0
	  rho     = 0d0
	  cp      = 0d0
	  cs      = 0d0
	  vm      = 0d0
	  Nv      = 0d0
	  beta    = 0d0
	  omega   = 0d0
	  KicSS   = 0d0
	  a0      = 0d0
	  fitKicD = 0d0
	  e0      = 0d0
 	  !-----------------------------------------------------------------------  
	  !Get the values from the input file (par.inp)
 	  !-----------------------------------------------------------------------  
	  read(iin, MAT_DYNDAMAGE, END=100)
	  call MAT_setKind(input,isDynDamage)
  
  	  call MAT_setProp(input,'fs',fs)
	  call MAT_setProp(input,'rho',rho)
	  call MAT_setProp(input,'cp',cp)
	  call MAT_setProp(input,'cs',cs)
	  call MAT_setProp(input,'vm',vm)
	  call MAT_setProp(input,'Nv',Nv)
	  call MAT_setProp(input,'beta',beta)
	  call MAT_setProp(input,'omega',Omega)
	  call MAT_setProp(input,'KicSS',KicSS)
	  call MAT_setProp(input,'a0',a0)
	  call MAT_setProp(input,'fitKicD',fitKicD)
  
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
 	  !-----------------------------------------------------------------------  
	  !Define others variables: shear modulus, lame parameter, poisson ratio
 	  !-----------------------------------------------------------------------  
	  if (cp>0d0 .and. cs>0d0 .and. rho>0d0) then
		mu = rho*cs*cs !shear Modulus
		lambda = rho*(cp*cp - 2d0*cs*cs) !lame parameters
		nu = 0.5d0-cs*cs/(2*cp*cp-2*cs*cs) !Poisson ratio
		cr = ((0.87 + 1.12*nu)/(1+nu))*cs
		if (nu < 0.d0 .or. nu > 0.5d0) call IO_abort('Poisson''s ratio out of range !')
		if (echo_input) write(iout,200) fs,Nv,KicSS,fitKicD,cr,cp,cs,vm,rho,nu,mu,lambda, &
					   lambda+2d0*mu/3d0, 2d0*mu*(1d0+nu)  ! Kvol,young  
		call MAT_setProp(input,'lambda',lambda)
		call MAT_setProp(input,'mu',mu)
	    call MAT_setProp(input,'nu',nu)
	    call MAT_setProp(input,'cr',cr)
	  else
		call IO_abort('MAT_DyDMG_read: incomplete input')
	  endif
 	  !-----------------------------------------------------------------------  
	  !Define damage variables
 	  !-----------------------------------------------------------------------  
	  phi=0.5d0*atan(1.0d0/fs) 
	  alpha=cos(phi)
	  D0=((4d0/3d0)*PI*Nv*(alpha*a0)**3)
	  D=D0
	  css=cs
	  cps=cp
	  v=0d0
	  KI=0d0 
	  R=1d0
	  Rcum=1d0
	  invI=0d0
	  invII=0d0
!	  A=0d0
!	  B=0d0
!	  C=0d0
 	  !-----------------------------------------------------------------------  
	  !assign initial value
 	  !-----------------------------------------------------------------------  
	  call MAT_setProp(input,'phi',phi)
	  call MAT_setProp(input,'alpha',alpha)
	  call MAT_setProp(input,'D0',D0)
	  call MAT_setProp(input,'D',D)
	  call MAT_setProp(input,'css',css)
	  call MAT_setProp(input,'cps',cps)
	  call MAT_setProp(input,'v',v)
	  call MAT_setProp(input,'KI',KI)
	  call MAT_setProp(input,'R',R)
	  call MAT_setProp(input,'Rcum',Rcum)
	  call MAT_setProp(input,'invI',invI)
	  call MAT_setProp(input,'invII',invII)
!	  call MAT_setProp(input,'A',A)
!	  call MAT_setProp(input,'B',B)
!	  call MAT_setProp(input,'C',C)
	  if (echo_input) write(iout,400) D0,D,phi*180.0d0/PI,alpha

  return
  100 call IO_abort('MAT_DyDMG_read: MAT_DYNDAMAGE input block not found')
  !-----------------------------------------------------------------------  
  ! simulation update info
  !-----------------------------------------------------------------------  
  200   format(5x, &
    'static friction coefficient . . . . .(fs) =',EN12.3,/5x, &
    'Volume density of cracks. . . . . . .(Nv) =',EN12.3,/5x, &
    'Quasi-static fracture toughness . (KicSS) =',EN12.3,/5x, &
    'For KicD . . . . . . . . . . . . (fitKicD) =',EN12.3,/5x, &
    'Rayleigh-wave velocity . . . . . . . (cr) =',EN12.3,/5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Branching speed . . . . . . . . . . .(vm) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Poisson''s ratio . . . . . . . . . . (nu) =',EN12.3,/5x, &
    'Shear modulus . . . . . . . . . . . .(mu) =',EN12.3,/5x, &
    'First Lame parameter Lambda . . .(lambda) =',EN12.3,/5x, &
    'Bulk modulus . . . . . . . . . . . . .(K) =',EN12.3,/5x, &
    'Young''s modulus . . . . . . . . . . .(E) =',EN12.3)

  250   format(5x, &
    '                                            antiplane',/5x, &
    'Initial total strain 13 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     23 . . . . . (e0(2)) =',EN12.3,/5x) 

  300   format(5x, &
    '                                            plane strain',/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x, &
    '                     12 . . . . . (e0(3)) =',EN12.3)
    
  400   format(5x, &
    '                                     damage variable',/5x, &
    'Initial Damage variable . . . . . . .(D0) =',EN12.3,/5x, &
    'Damage variable. . . . . . . . . . . .(D) =',EN12.3,/5x, &
    'angle to sigma1. . . . . . . . . . .(phi) =',EN12.3,/5x, &
    'projection of a0 to sigma1 . . . .(alpha) =',EN12.3)
  
end subroutine MAT_DyDMG_read

!=======================================================================
! Initialise material properties for one element
!=======================================================================
  subroutine MAT_DyDMG_init_elem_prop(elem,ecoord)

      !input/ouput variables
	  type(matpro_elem_type), intent(inout) :: elem
	  double precision, intent(in) :: ecoord(:,:,:)

 	  !-----------------------------------------------------------------------  
	  ! assign value for each element
 	  !-----------------------------------------------------------------------  
	  call MAT_setProp(elem,'fs',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'nu',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'mu',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'rho',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'cr',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'cp',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'Nv',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'beta',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'omega',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'KicSS',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'a0',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'phi',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'alpha',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'D0',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'D',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'vm',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'fitKicD',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'lambda',ecoord,MAT_DyDMG_mempro)
	  call MAT_setProp(elem,'cs',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'css',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'cps',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'v',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'KI',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'R',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'Rcum',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'invI',ecoord,MAT_DyDMG_mempro)   
	  call MAT_setProp(elem,'invII',ecoord,MAT_DyDMG_mempro)   
!	  call MAT_setProp(elem,'A',ecoord,MAT_DyDMG_mempro)   
!	  call MAT_setProp(elem,'B',ecoord,MAT_DyDMG_mempro)   
!	  call MAT_setProp(elem,'C',ecoord,MAT_DyDMG_mempro)   
 
	  !SH antiplane
	  if (MAT_isKind(elem,isAntiplane)) then 
		call MAT_setProp(elem,'e13_0',ecoord,MAT_DyDMG_mempro)
		call MAT_setProp(elem,'e23_0',ecoord,MAT_DyDMG_mempro)
	  !P-SV Plane strain/ inplane
	  else 
		call MAT_setProp(elem,'e11_0',ecoord,MAT_DyDMG_mempro)
		call MAT_setProp(elem,'e22_0',ecoord,MAT_DyDMG_mempro)
		call MAT_setProp(elem,'e12_0',ecoord,MAT_DyDMG_mempro)
	  endif  

  end subroutine MAT_DyDMG_init_elem_prop

!=======================================================================
! Assign material properties for one element
!=======================================================================
 subroutine MAT_DyDMG_init_elem_work(m,p,n,ndof)

    !input/ouput variables
	type(matwrk_dyndmg_type), intent(inout) :: m
	type(matpro_elem_type), intent(in) :: p
	integer, intent(in) :: n, ndof

	!local variables
	double precision :: mu,nu
	
	!-----------------------------------------------------------------------  
	! Constant variables for damage material
	!-----------------------------------------------------------------------
		call MAT_getProp(m%fs,p,'fs')
		call MAT_getProp(m%nu,p,'nu')
		call MAT_getProp(m%mu,p,'mu')
		call MAT_getProp(m%rho,p,'rho')
		call MAT_getProp(m%cr,p,'cr')
		call MAT_getProp(m%Nv,p,'Nv')
		call MAT_getProp(m%beta,p,'beta')
		call MAT_getProp(m%omega,p,'omega')
		call MAT_getProp(m%a0,p,'a0')
		call MAT_getProp(m%KicSS,p,'KicSS')
		call MAT_getProp(m%alpha,p,'alpha')
		call MAT_getProp(m%D0,p,'D0')
		call MAT_getProp(m%vm,p,'vm')
		call MAT_getProp(m%fitKicD,p,'fitKicD')
		call MAT_getProp(m%lambda,p,'lambda')
	!-----------------------------------------------------------------------
	! Damage variable
	!-----------------------------------------------------------------------
		allocate(m%D(n,n))
		call MAT_getProp(m%D,p,'D')
	!-----------------------------------------------------------------------
	! Stress Intensity factor
	!-----------------------------------------------------------------------
		allocate(m%KI(n,n))
		call MAT_getProp(m%KI,p,'KI')
!		allocate(m%A(n,n))
!		call MAT_getProp(m%A,p,'A')
!        allocate(m%B(n,n))
!		call MAT_getProp(m%B,p,'B')
!        allocate(m%C(n,n))
!		call MAT_getProp(m%C,p,'C')
	!-----------------------------------------------------------------------  
	! Regime
	!-----------------------------------------------------------------------
		allocate(m%R(n,n))
		call MAT_getProp(m%R,p,'R')
		allocate(m%Rcum(n,n))
		call MAT_getProp(m%Rcum,p,'Rcum')
        allocate(m%invI(n,n))
		call MAT_getProp(m%invI,p,'invI')
        allocate(m%invII(n,n))
		call MAT_getProp(m%invII,p,'invII')
	!-----------------------------------------------------------------------  
	! Apparent S and P waves
	!-----------------------------------------------------------------------
		allocate(m%css(n,n))
		call MAT_getProp(m%css,p,'css')
		allocate(m%cps(n,n))
		call MAT_getProp(m%cps,p,'cps')
	!-----------------------------------------------------------------------  
	! Instantaneous wing-crack speed v=dl/dt
	!-----------------------------------------------------------------------
		allocate(m%v(n,n))
		call MAT_getProp(m%v,p,'v')
	!-----------------------------------------------------------------------
	! initial strain
	!-----------------------------------------------------------------------
		allocate(m%e(n,n,ndof+1))
		allocate(m%erate(n,n,ndof+1))
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
		m%erate(:,:,:) = 0.0d0
	!----------------------------------------------------------------------- 
	! initial stress (assume no Damage)
	!-----------------------------------------------------------------------
		nu = m%nu
		mu = m%mu
		allocate(m%s0(3))
		if (ndof==1) then !SH antiplane
		  m%s0(1) = 2.0*mu*m%e0(1)! sigma_13
		  m%s0(2) = 2.0*mu*m%e0(2)! sigma_23
		else !PSV
		  m%s0(1) = 2.0*mu*m%e0(1) + ((2.0*mu*nu)/(1.0-2.0*nu))*(m%e0(2)+m%e0(1))! sigma_11
		  m%s0(2) = 2.0*mu*m%e0(2) + ((2.0*mu*nu)/(1.0-2.0*nu))*(m%e0(2)+m%e0(1))! sigma_22
		  m%s0(3) = 2.0*mu*m%e0(3)! sigma_12
		endif
	!-----------------------------------------------------------------------
	! assign memory
	!-----------------------------------------------------------------------
		MAT_DyDMG_memwrk = MAT_DyDMG_memwrk &
					   + size( transfer(m, (/ 0d0 /) )) &
					   + size(m%D) &
					   + size(m%css) &
					   + size(m%cps) &
					   + size(m%v) &
					   + size(m%KI) &
					   + size(m%R) &
					   + size(m%Rcum) &
					   + size(m%invI) &
					   + size(m%invII) &
!					   + size(m%A) &
!					   + size(m%B) &
!					   + size(m%C) &
					   + size(m%e) &
					   + size(m%erate) &
					   + size(m%e0) &
					   + size(m%s0)
	!-----------------------------------------------------------------------
	!	To set up the constant variables use in damage update
	!-----------------------------------------------------------------------  
		call ini_dyndmg_cst(m)   
end subroutine MAT_DyDMG_init_elem_work 

!=======================================================================
! Main Subroutine for Constitutive Update
! Compute Stresses from given Strains and state variables
!=======================================================================
  subroutine MAT_DyDMG_stress(s,etot,m,ngll,ndof,update,dt,tt)

  !input/ouput variables
  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: etot(ngll,ngll,ndof+1)
  double precision, intent(out) :: s(ngll,ngll,ndof+1)
  type (matwrk_dyndmg_type), intent(inout) :: m
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt
  double precision, optional, intent(in) :: tt

  !local variables
  double precision, dimension(ngll,ngll,ndof+1) :: e,erate,eprevious
  double precision, dimension(ngll,ngll) :: D,KItemp,Rtemp,Rctemp
  double precision, dimension(ngll,ngll) :: csstemp,cpstemp,vtemp
  double precision, dimension(ngll,ngll) :: invItemp,invIItemp
!  double precision, dimension(ngll,ngll) :: Atemp,Btemp,Ctemp
  double precision, dimension(ndof+1) :: eij,sij,erij
  double precision :: itime,etime,iti,eti
  double precision :: alphaij
  double precision :: strain(3,3),stress(3,3), strain_rate(3,3)
  double precision :: sigma,tau,A,B,C,O,R,Dsav
  integer :: i,j
  double precision ::KIdot,KicD,dldt
  double precision :: mu,lambda
  double precision :: css,cps
  double precision :: dumcss,dumcps,dumepsi,dumgamma,dumnu,test1
  integer ::ii,jj
  
  !for the Runge-Kutta-Fehlberg method
  integer ( kind = 4 ), parameter :: neqn = 1
  double precision :: y(neqn)
  double precision :: yp(neqn)
  double precision :: dDdt(neqn)
  double precision :: abserr,relerr
  integer ( kind = 4 ) flag
  
  !-----------------------------------------------------------------------
  !	elastic constants
  !-----------------------------------------------------------------------
	  lambda = m%lambda
	  mu = m%mu
  !-----------------------------------------------------------------------
  !	define error tolerances for ODE solver
  !-----------------------------------------------------------------------
	  abserr = 1.0d-9
	  relerr = 1.0d-9  
  !-----------------------------------------------------------------------  
  !	get Damage related constants
  !-----------------------------------------------------------------------
	  call ini_dyndmg_cst(m)
  !----------------------------------------------------------------------- 
  !	Update strain
  !-----------------------------------------------------------------------
	  e(:,:,1) = etot(:,:,1) + m%e0(1)
	  e(:,:,2) = etot(:,:,2) + m%e0(2)
		  if (ndof==2) then !PSV
			e(:,:,3) = etot(:,:,3) + m%e0(3)
		  endif  
  !-----------------------------------------------------------------------
  !	Compute the strain rate
  !-----------------------------------------------------------------------
	  if (update) then
		  m%erate(:,:,1)=(-m%e(:,:,1)+e(:,:,1))/dt
		  m%erate(:,:,2)=(-m%e(:,:,2)+e(:,:,2))/dt
	      if (ndof==2) then !PSV
			m%erate(:,:,3)=(-m%e(:,:,3)+e(:,:,3))/dt
  		  endif  
	      if (.not.present(dt)) &
		  call IO_abort('mat_damage:MAT_DMG_stress: update requested but argument dt is absent')
  !-----------------------------------------------------------------------
  ! Start and End time for RKF45 ODE integration	
  !-----------------------------------------------------------------------
	  	iti = tt
	  	eti = tt+dt
	  endif
  !-----------------------------------------------------------------------
  !	Get Saved Variables from previous time-step
  !-----------------------------------------------------------------------
	  D(:,:)=m%D
	  KItemp(:,:)=m%KI
	  erate =m%erate
	  csstemp(:,:)= m%css
	  cpstemp(:,:)= m%cps
	  vtemp(:,:)= m%v
	  Rtemp(:,:)= m%R
	  Rctemp(:,:)= m%Rcum
	  invItemp(:,:)= m%invI
	  invIItemp(:,:)= m%invII
!	  Atemp(:,:)= m%A
!	  Btemp(:,:)= m%B
!	  Ctemp(:,:)= m%C
	  eprevious=m%e
	  m%e = e
  !-----------------------------------------------------------------------
  ! CORE PHILOSOPHY OF Finite Element and Spectral Element Methods
  !  
  !	At each integration point calculate stress at t+dt given 
  ! strain at t+dt and stress, Damage at t i.e.
  ! Compute sigma(t+dt) based on epsilon(t+dt), sigma(t) and D(t)
  !
  !-----------------------------------------------------------------------  
	  do j=1,ngll
		  do i=1,ngll 
			  !-----------------------------------------------------------------------	
			  !	previous Damage
			  !-----------------------------------------------------------------------	
				y(1)=D(i,j)
				Dsav=D(i,j)
				flag = 1
			  !-----------------------------------------------------------------------	
			  !	strain and strain rate
			  !-----------------------------------------------------------------------	
				eij = eprevious(i,j,:)
				erij = erate(i,j,:)
				call vec2mat(strain,eij,ndof)
				call vec2mat(strain_rate,erij,ndof)
			  !-----------------------------------------------------------------------
			  !	Update Variables !F90 bizzareness
			  !-----------------------------------------------------------------------	
				itime=iti
			    etime=eti  
				call damage_var(strain, strain_rate, KItemp(i,j), itime,etime, Rtemp(i,j)) 
			  !-----------------------------------------------------------------------		
			  !	Damage evolution (Runge-Kutta-Fehlberg method)
			  !	Solve dD/dt = derivs() from t = tt to t = tt+dt with D(tt) = D(i,j)	
			  !-----------------------------------------------------------------------	   
				call r8_rkf45(derivs, neqn, y(1), yp(1), itime,etime, relerr, abserr, flag)
			  !-----------------------------------------------------------------------
			  !	Updated Damage, D
			  !-----------------------------------------------------------------------	  
				D(i,j)=y(1)   
			  !-----------------------------------------------------------------------
			  !	Sanity Check	
			  !-----------------------------------------------------------------------      
				if (D(i,j).ge.0.99d0) then
				  D(i,j)=0.99d0
				else if (D(i,j).lt.cD0) then
				  D(i,j)=cD0
				end if  
			  !-----------------------------------------------------------------------	
			  !	Compute the KI parameters needed
			  !-----------------------------------------------------------------------
				call KI_parameters(D(i,j),A,B,C,O)
!				Atemp(i,j) = A
!				Btemp(i,j) = B
!				Ctemp(i,j) = C
			  !-----------------------------------------------------------------------	
			  !	Do a final check of the regime and Compute stress
			  !-----------------------------------------------------------------------
				eij =m%e(i,j,:)
				call vec2mat(strain,eij,ndof)
				call stress_damage(R,strain,A,B,C,O,stress,D(i,j),dumcss,dumcps,Rtemp(i,j))
				if (D(i,j).gt.m%D(i,j)) then
				   if (dumcss.le.m%css(i,j)) then
					   csstemp(i,j) = dumcss
				   end if
				   if (dumcps.le.m%cps(i,j)) then
					   cpstemp(i,j) = dumcps
				   end if
				end if
				Rtemp(i,j) = R
				if (Rctemp(i,j).lt.R) then
					Rctemp(i,j) = R
				end if
				call mat2vec(stress,sij,ndof)
				s(i,j,:) = sij
				call stress_invariant(stress,sigma,tau)
				invIItemp(i,j) = tau
				invItemp(i,j) = sigma				
			  !-----------------------------------------------------------------------	
			  !	Compute the stress intensity factor
			  !-----------------------------------------------------------------------
				call compute_KI(A,B,C,O,sigma,tau,D(i,j),KItemp(i,j),R)	  
			  !-----------------------------------------------------------------------	
			  !	Compute v=dl/dt crack speed
			  !-----------------------------------------------------------------------
				call derivs(end_time,D(i,j),dDdt(1))
				vtemp(i,j)=dDdt(1)/((3.0*D(i,j)**(2.0/3.0)*cD0**(1.0/3.0))/(calpha*ca0))
			  enddo
	  enddo  
  !-----------------------------------------------------------------------
  !	Update variables to for next times step and output requests
  !-----------------------------------------------------------------------
	  if (update) then
		m%KI=KItemp
		m%D=D
		m%css=csstemp
		m%cps=cpstemp
		m%v=vtemp
		m%R=Rtemp
		m%Rcum=Rctemp
		m%invI=invItemp
		m%invII=invIItemp
!		m%A=Atemp
!		m%B=Btemp
!		m%C=Ctemp
	  endif
  !-----------------------------------------------------------------------
  !	relative stresses
  !-----------------------------------------------------------------------
	  s(:,:,1) = s(:,:,1) - m%s0(1)
	  s(:,:,2) = s(:,:,2) - m%s0(2)
		  if (ndof==2) then !
			  s(:,:,3) = s(:,:,3) - m%s0(3)
		  endif
  end subroutine MAT_DyDMG_stress
  
!=======================================================================
! Extra data to be written to output binary file
!=======================================================================
  function MAT_DyDMG_export(m,ndof) result(dat)
  ! dat has to be the size of ngll,ngll,#of output variables and you need to change the matlab
  ! file "sem2d_snapshot_read.m" consequently in "elseif ~isempty(strfind(fname_pre,'elm_'))"

  !input/ouput variables
  integer, intent(in) :: ndof
  type(matwrk_dyndmg_type), intent(in) :: m
  real :: dat(size(m%D,1),size(m%D,2),8)
  
  !-----------------------------------------------------------------------
  !assign values
  !-----------------------------------------------------------------------
  dat(:,:,1) = real(m%D)
  dat(:,:,2) = real(m%css)
  dat(:,:,3) = real(m%cps)
  dat(:,:,4) = real(m%v)
  dat(:,:,5) = real(m%R)
  dat(:,:,6) = real(m%invI)
  dat(:,:,7) = real(m%invII)
  dat(:,:,8) = real(m%Rcum)
!  dat(:,:,8) = real(m%KI)
!  dat(:,:,9) = real(m%A)
!  dat(:,:,10) = real(m%B)
!  dat(:,:,11) = real(m%C)

  !real :: dat(size(m%D,1),size(m%D,2),ndof+ndof+3)
  !dat(:,:,1) = real(m%D)
  !dat(:,:,2:ndof+2) = real(m%e)
  !dat(:,:,ndof+3:ndof+5) = real(m%s)

  end function MAT_DyDMG_export

!=======================================================================
! Assign Constants used in the Damage model
!=======================================================================
  subroutine ini_dyndmg_cst(m)
!-----------------------------------------------------------------------
  !	fs		: static friction coefficient
  !	nu		: Poisson ratio
  !	mu 		: shear Modulus (Pa)
  !	rho		: density (kg/m^3)
  !	cr 		: Rayleigh wave velocity (m/s)
  !	Nv 		: volume density of cracks (#/m^3)
  !	beta	: Ashby and Sammis (1990) factor (usually 0.1)
  !	omega	: Cracks factor (usually 2.0)
  !	KicSS	: Quasi-static fracture toughness
  !	a0 		: initial radius of penny shape cracks
  !	alpha	: projection of the crack radius in a vertical plane parallel to 
  !       		the direction of sigma1
  !	D0   	: Initial Damage parameter
  !	vm		: Branching speed
  !	fitKicD	: Curve fit parameter to compute the dynamic initation toughness
!-----------------------------------------------------------------------

  !input/ouput variables
  type(matwrk_dyndmg_type), intent(in) :: m

  !-----------------------------------------------------------------------
  !assign constant
  !-----------------------------------------------------------------------
  nn=3
  cfs=m%fs
  cnu=m%nu 
  cmu=m%mu 
  crho=m%rho 
  ccr=m%cr 
  cNv=m%Nv 
  cbeta=m%beta 
  comega=m%omega 
  ca0=m%a0 
  cKicSS=m%KicSS 
  calpha=m%alpha 
  cD0=m%D0 
  cvm=m%vm 
  cfitKicD=m%fitKicD

  end subroutine ini_dyndmg_cst  
  
!=======================================================================
! Variables that need to be saved for the Damage model (Fortran QUIRCK)
!=======================================================================
 subroutine damage_var(strainINI,strain_R,KIcomp,xit,xet,Rcomp)
!-----------------------------------------------------------------------
  !	strainINI	: initial strain for the RKF loop
  !	strain_R  	: strain_rate
  !	KIcomp		: stress intensity factor
  !	xit         : time at the beginning of the ODE integration
  !	xet       	: output time for the ODE integration
!-----------------------------------------------------------------------
  
  !input/ouput variables
  double precision :: strainINI(3,3),strain_R(3,3),KIcomp,xit,xet,Rcomp

  !save variables
  strain0=strainINI
  str_rate = strain_R
  KIsav=KIcomp
  Rsav=Rcomp
  start_time=xit
  end_time=xet

 end subroutine damage_var
 
!=======================================================================
! Compute the Damage Dependent Constants used to calculate K_I
!=======================================================================
  subroutine KI_parameters(D,A,B,C,O)
!-----------------------------------------------------------------------
  !INPUT VARIABLES
  !----------------
  !	D	: Damage parameter
  !OUTPUT VARIABLES
  !----------------
  !	A,B,C,O	: are functions used to compute the stress intensity factor KI
  ! Note "O" is "E" in Bhat et al. (2012)
!-----------------------------------------------------------------------
  !modules
  use constants, only : PI
  implicit none
  
  !input/ouput variables
  double precision, intent(inout) :: D
  double precision, intent(out) :: A,B,C,O

  !local variables
  double precision :: c1,c2,c3
  double precision :: alpha, D0,beta,omega,fs
  
  !constants
  alpha=calpha
  D0=cD0
  beta=cbeta
  omega=comega
  fs=cfs
  
  !-----------------------------------------------------------------------
  !	Sanity Check	
  !-----------------------------------------------------------------------      
	  if (D.ge.0.99d0) then
		D = 0.99d0
	  else if (D.lt.D0) then
		D = D0
	  end if
  !-----------------------------------------------------------------------
  !	Damage Dependent Constants
  !-----------------------------------------------------------------------      
	  c1 = sqrt(1.0d0-alpha**2)/((PI*alpha**(3.0d0/2.0d0))* &
	  			((D/D0)**(1.0d0/3.0d0)-1.0d0+beta/alpha)**(3.0d0/2.0d0))
	  c2 = (sqrt(1.0d0-alpha**2)/alpha**2)*((D0**(2.0d0/3.0d0))/(1.0d0-D**(2.0d0/3.0d0)))
	  c3 = ((2.0d0*sqrt(alpha))/PI)*sqrt((D/D0)**(1.0d0/3.0d0)-1.0d0)
  !-----------------------------------------------------------------------
  !	Damage Dependent Constants
  !-----------------------------------------------------------------------      
	  B = c1+c2*c3
	  A = fs*c1+(1.0d0+fs*c2)*c3
	  C = A+omega*sqrt(alpha*(D/D0)**(1./3.))
	  O = B*C/sqrt((C**2-A**2)) 

  end subroutine KI_parameters
    
!=======================================================================
! check the regime and compute the stress intensity factor KI.
!=======================================================================
  subroutine compute_KI(A,B,C,O,sigma,tau,D,KI,R)
!-----------------------------------------------------------------------
  !INPUT VARIABLES
  !----------------
  !	A,B,C,O	: are functions used to compute the stress intensity factor KI
  ! Note "O" is "E" in Bhat et al. (2012)
  !	sigma: first invariant of the stress tensor
  !	tau	: second invariant of the deviatoric stress tensor
  !	D	: Damage parameter
  !	R	: Deformation Regime 
  !OUTPUT VARIABLES
  !----------------
  !	KI	: stress intensity factor
!-----------------------------------------------------------------------
  !modules used
  use constants, only : PI
  implicit none
  
  !input/ouput variables
  double precision, intent(in) :: A,B,C,O, sigma,tau,D,R
  double precision, intent(out) :: KI
 
  !local variables
  double precision :: KI1,a0
  
  !constants
  a0=ca0

  !-----------------------------------------------------------------------  
  !	Scaling Factor for KI
  !-----------------------------------------------------------------------  
	  KI1 = sqrt(PI*a0)
  !-----------------------------------------------------------------------  
  !	Check the Regime and compute the stress intensity factor KI
  !-----------------------------------------------------------------------  
	  if ((R-1.0).eq.0d0) then
		KI = 0.0d0
	  else if ((R-2.0).eq.0d0) then
		KI = KI1*(A*sigma+B*tau)   
	  else if ((R-3.0).eq.0d0) then
		KI = KI1*sqrt((C**2)*(sigma**2)+(O**2)*(tau**2))
	  else 
		call IO_abort('No regime!! HELP!')
	  end if
   
  end subroutine compute_KI  
  
!=======================================================================
! strain-stress relationship based on Bhat et al. (2012)
!=======================================================================
  subroutine stress_damage(R,strain,A,B,C,O,stress,D,css,cps,Rs)
!-----------------------------------------------------------------------
  !INPUT VARIABLES
  !----------------
  !	strain	: strain tensor
  !	A,B,C,O	: are functions used to compute the stress intensity factor KI
  ! Note "O" is "E" in Bhat et al. (2012)
  !	D		: Damage
  ! Rs		: Previous Regime
  !OUTPUT VARIABLES
  !----------------
  !	R		: Deformation Regime
  !	stress	: stress tensor
  ! css   	: Apparent S wave velocity
  ! cps   	: Apparent P wave velocity
!-----------------------------------------------------------------------
  use constants, only : PI
  use stdio, only : IO_abort
  use utils_dyn_damage
  implicit none
  
  !input/ouput variables
  integer, parameter :: n=3
  double precision, intent(in) :: strain(n,n),A,B,C,O,D,Rs
  double precision, intent(out) ::  R, stress(n,n),css,cps

  !local variables
  double precision :: epsi, gamma, delta(n,n)!, 
  double precision :: Dvar,gammacst,cst1,cst2,cst3,A1,C1,B1,O1
  double precision :: fac,test1,test2
  integer :: i,j
  double precision :: nu,mu,alpha,D0,rho
  double precision :: s(3),e(3),lambda
  double precision :: mus
  double precision :: lambdas
  
  !Constants
  nu=cnu
  mu=cmu
  alpha=calpha
  D0=cD0
  rho=crho
  lambda=(2.0d0*mu*nu)/(1.0d0-2.0d0*nu)
  
  !-----------------------------------------------------------------------  
  !	compute strain invariants and kronecker delta function
  !-----------------------------------------------------------------------  
	  call strain_invariant(strain,epsi,gamma)
	  call kronecker(delta,n)
  !-----------------------------------------------------------------------  
  !	factor useful for calculations Regime II and III
  !-----------------------------------------------------------------------    
	  Dvar=sqrt((PI*D0*(1-nu))/(alpha**3))
  !-----------------------------------------------------------------------  
  ! Check Regime
  !-----------------------------------------------------------------------
	  test1=-3.0*gamma*B*(1.0-2.0*nu)/(2.0*A*(1.0+nu)) 
	  test2=3.0*(1.0-2.0*nu)/(2.0*mu*(1.0+nu))+C**2*Dvar**2/(2*mu)
	  test2=test2/(1.0/mu+B**2*C**2*Dvar**2/(2*mu*(C**2-A**2)))
	  test2=test2*A*B/(C**2-A**2)*gamma
	  !-----------------------------------------------------------------------
	  !	Regime 1 
	  !if (epsi.le.test1) then
	  !-----------------------------------------------------------------------
	  !if (epsi.le.test1 .AND. D==D0) then
	  if (epsi.le.test1) then
	    R=1.0d0
	  !-----------------------------------------------------------------------
	  !	Regime 2 
	  !else if  (epsi.gt.test1 .AND. D.gt.D0 .AND. epsi.lt.test2) then
	  !-----------------------------------------------------------------------
	  else if  (epsi.gt.test1 .AND. epsi.le.test2) then
		R=2.0d0
	  !-----------------------------------------------------------------------  
	  !	Regime 3 
	  !  elseif  (epsi.gt.test1 .AND. D.gt.D0 .AND. epsi.gt.test2) then
	  !-----------------------------------------------------------------------
	  else if  (epsi.gt.test1 .AND. epsi.gt.test2) then
		R=3.0d0  	
	  else 
		!call IO_abort('No regime!! HELP!')
		R=Rs
		!print*, 'we are in no regime mode', Rs, Rsav
	  end if
  !------------------------------------------------------------------------
  !	Compute the stress based on the Regime
  !------------------------------------------------------------------------  
	  !------------------------------------------------------------------------  
	  !	Regime 1 : stress
	  !------------------------------------------------------------------------
	  	if ((R-1.0).eq.0d0) then
	   		!stress
	   		stress=2.0*mu*(strain+(nu/(1.0-2.0*nu))*epsi*delta)
	   		!Apparent shear and lame modulus
	   		mus     = mu
	   		lambdas = lambda
	  !------------------------------------------------------------------------
	  !	Regime 2 : stress
	  !------------------------------------------------------------------------
	  	else if ((R-2.0).eq.0d0) then
			A1 = A*Dvar
			B1 = B*Dvar
			gammacst = (3.0*(1.0-2.0*nu))/(2.0*(1.0+nu)) + &
					   (3.0*(1.0-2.0*nu)*B1**2)/(4.0*(1.0+nu)) + (A1**2)/2.0
			cst1 = (3.0*(1.0-2.0*nu))/(1.0+nu) + A1**2 - (A1*B1*epsi)/gamma
			cst2 = (3.0*nu)/(1.0+nu) + (B1**2.0)/2.0 - (A1**2)/3.0 + (A1*B1*epsi)/(3.0*gamma)
			cst3 = -(A1*B1)/2.0
	   		!stress
			stress = (mu/gammacst)*(cst1*strain + cst2*epsi*delta + cst3*gamma*delta)
	   		!Apparent shear and lame modulus
	   		mus     = (1.0/2.0)*(mu/gammacst)*(3.0*(1.0-2.0*nu)/(1+nu) + A1**2)
	   		lambdas = (mu/gammacst)*((3.0*nu)/(1.0+nu) + (B1**2.0)/2.0 - (A1**2)/3.0)
	  !------------------------------------------------------------------------
	  !Regime 3 : stress
	  !------------------------------------------------------------------------
	  	else if ((R-3.0).eq.0d0) then
			C1 = C*Dvar
			O1 = O*Dvar
			cst1 = ((3.0*(1-2*nu))/(1+nu)+C1**2)**(-1)
			cst2 = (2.0+O1**2)**(-1)
	   		!stress
			stress = mu*(4.0*cst2*strain + (2*cst1-4.0*cst2/3.0)*epsi*delta);
	   		!Apparent shear and lame modulus
	   		mus     = 2.0*mu*cst2
	   		lambdas = mu*(2.0*cst1-4.0*cst2/3.0)
	  !------------------------------------------------------------------------  
	  	else 
			call IO_abort('Regime does not exist, cannot compute the stress')
	  	end if
  !------------------------------------------------------------------------
  !	Compute the apparent S and P waves velocity
  !------------------------------------------------------------------------  
	  	css = sqrt(mus/rho)
	  	cps = sqrt((lambdas+2*mus)/rho)
	  	
end subroutine stress_damage

!=======================================================================
! Newton Raphson solver for crack rupture velocity
!=======================================================================
  subroutine  rupture_velocity(x,KI)
!-----------------------------------------------------------------------
  !INPUT VARIABLES
  !----------------
  !	x		: rupture speed
  !	KicSS	: quasi-static fracture toughness
  !	KI		: stress intensity factor
  !	vm		: branching speed
  !	cr		: Rayleigh waves speed
  !	fvalue	: function for which we want to find the root
  !	dfvalue	: derivative of function f
  !	x0		: initial guess for the root
  !	cst		: csontant value you might need in function f
  !OUTPUT VARIABLES
  !----------------
  !	x 	: root (crack velocity)
!-----------------------------------------------------------------------
implicit none

  !input/ouput variables
  double precision, intent(in) :: KI
  double precision, intent(out) :: x 

  !local variables
  double precision :: fvalue,dfvalue, x0, a00, b00, c00
  integer :: iter ! current iteration number
  real, parameter :: tolerance = 1.0e-3 ! tolerance for near-zero
  integer, parameter :: itermx = 100000 ! maximum number of iterations
  double precision :: KicSS,cr,vm
  
  !Get constant value
  KicSS=cKicSS
  cr=ccr
  vm=cvm

  !Definition of local variables
  a00 = 1.0d0/(vm**5)
  b00 = KI/(KicSS*cr)
  c00 = (1.0d0-KI/KicSS)
  iter = 0

  !------------------------------------------------------------------------
  !	initial guess
  !------------------------------------------------------------------------
  x0=vm*0.5d0
  fvalue=a00*x0**5+b00*x0+c00 !Crack evolution
  dfvalue=5.0d0*a00*x0**4+b00 !derivative
  x=x0
  !------------------------------------------------------------------------
  !	loop until root found or maximum iterations reached
  !	Newton-Raphson method
  !------------------------------------------------------------------------
	  do
		x = x - fvalue / dfvalue 		! update x by newton-raphson formula
		fvalue = a00*x**5+b00*x+c00 	! update value of function
		dfvalue = 5.0d0*a00*x**4+b00
		iter = iter + 1 				! update iteration number
		if (abs(fvalue).lt.tolerance .or. iter.ge.itermx) then
		  exit
		end if
	  end do
  !------------------------------------------------------------------------
  !	Check that v is in the range [0, vm]
  !------------------------------------------------------------------------  
	  if (abs(fvalue).gt.tolerance .or. iter.ge.itermx) then
		x = 0.0d0
	  elseif (x.gt.vm) then
		x = vm
	  elseif (x.lt.0.0d0) then
		x = 0.0d0
	  end if
 
end subroutine rupture_velocity

!=======================================================================
! Solve for damage evolution
! dD/dt = F(t,D,sigma,epsilon,KI,KICD,dl/dt)
! derivs outputs dD/dt
!=======================================================================
  subroutine derivs(x,y,dydx)
!-----------------------------------------------------------------------
  !	INPUT/OUTPUT VAIRABLES
  !----------------------
  !	x	: independent variable (time)
  !	y	: the solution vector at x of the dependent varibale (damage)
  !	dydx	: The current value of the derivative the dependent variable
!-----------------------------------------------------------------------
	implicit none

	!input/ouput variables
	integer ( kind = 4 ), parameter :: neqns = 1
	double precision, intent(in) :: x,y(neqns)
	double precision, intent(out) :: dydx(neqns)

	!local variables
	integer :: i,j
	double precision :: strain(nn,nn), stress(nn,nn),eige(nn), eigv(nn,nn)
	double precision :: A,B,C,O,R,KI,sigma,tau,KIdot,KicD,dldt,D
	double precision :: KicSS,fitKicD,D0,alpha,a0
	double precision :: css,cps
   
	!constants
	KicSS=cKicSS
	fitKicD=cfitKicD
	D0=cD0
	alpha=calpha
	a0=ca0
    
	!-----------------------------------------------------------------------
	!	Initial Condition to solve the ODE
	!-----------------------------------------------------------------------
		D=y(1)
	!-----------------------------------------------------------------------
	!	Compute the KI parameters
	!-----------------------------------------------------------------------
		call KI_parameters(D,A,B,C,O)
	!-----------------------------------------------------------------------
	!	Update strain 
	!	Linear interpolation assuming constant strain rate between t and t+dt
	!	strain(t+c*dt) = strain(t) + c*dt*strainrate(t)
	!	c is in the range [0,1] and is chosen by the adaptive RKF45 solver
	!-----------------------------------------------------------------------
		do i=1,nn
		  do j=1,nn
			strain(i,j)=strain0(i,j)+(x-start_time)*str_rate(i,j)!
		  enddo
		enddo
	!-----------------------------------------------------------------------
	!	Update stress
	!-----------------------------------------------------------------------
		call stress_damage(R,strain,A,B,C,O,stress,D,css,cps,Rsav)
		call stress_invariant(stress,sigma,tau)
	!-----------------------------------------------------------------------
	!	Compute the stress intensity factor
	!-----------------------------------------------------------------------
		call compute_KI(A,B,C,O,sigma,tau,D,KI,R)
	!-----------------------------------------------------------------------
	!	Compute the stress intensity factor rate
	!	for Dynamic Initiation Toughness K_{IC}^{D}
	! 	Approximated as (K(t+dt)-K(t))/dt instead of
	!	(K(t+c*dt)-K(t))/(c*dt) because the latter is quite noisy
	!-----------------------------------------------------------------------
		if ((x-start_time).gt.0.0d0) then
!		  KIdot=(KI-KIsav)/(end_time-start_time) 
		  KIdot=(KI-KIsav)/(x-start_time)
		!-----------------------------------------------------------------------
		!	If K(t) < K(t-dt) then dK/dt = 0 
		!-----------------------------------------------------------------------
		elseif (KI.lt.KIsav) then
		  KIdot=0.0d0
		end if
	!-----------------------------------------------------------------------
	!	Extra cautious step to ensure dK/dt >= 0 
	!-----------------------------------------------------------------------
		if (KIdot.lt.0.0d0) then
		  KIdot=0.0d0
		endif
	!-----------------------------------------------------------------------
	!	Compute the Dynamic Initiation Toughness K_{IC}^{D}
	!-----------------------------------------------------------------------
		KicD = (1.0d0 + KIdot/fitKicD)*KicSS; 
	!-----------------------------------------------------------------------
	!	Compute dl/dt (or v_r), the micro-crack rupture velocity
	!-----------------------------------------------------------------------
		if (KI.lt.KicD) then
		  dldt=0.0d0
		else 
		  call rupture_velocity(dldt, KI)
		end if
	!-----------------------------------------------------------------------
	!	Damage Evolution Law (dD/dt)
	!-----------------------------------------------------------------------
		dydx(1)=((3.0*y(1)**(2.0/3.0)*D0**(1.0/3.0))/(alpha*a0))*dldt
	!-----------------------------------------------------------------------
	!	More Sanity Checks
	!-----------------------------------------------------------------------
		if (y(1).ge.0.99d0) then
			dydx(1)=0.0d0
		else if (y(1).lt.D0) then
			dydx(1)=0.0d0
		end if
	!-----------------------------------------------------------------------
	return

	end subroutine derivs

end module mat_dyn_damage
