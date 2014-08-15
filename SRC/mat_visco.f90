module mat_visco
! Added by Yihe (2012)
! viscoelastic medium following Mozco (2004)
  
  use attenuation
  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private

  type matwrk_visco_type
    private
    integer :: Nbody
    double precision :: mu, lambda
    double precision, pointer, dimension(:,:) :: theta => null() ! anelastic memory variables
    double precision, pointer, dimension(:,:,:,:) :: el => null() ! Anelastic function in loop
    double precision, pointer, dimension(:,:,:) :: etot_old => null() ! Store strain at last time step
!    double precision, pointer, dimension(:) :: e0 => null(), s0 => null()
    double precision, pointer, dimension(:) :: wbody => null() ! central frequencies of viscoelastic mechanisms
  end type matwrk_visco_type

  integer, save :: isVisco = 0

  ! for memory report
  integer, save:: MAT_VISCO_mempro = 0
  integer, save:: MAT_VISCO_memwrk = 0

  public :: matwrk_visco_type &
          , MAT_isVisco, MAT_anyVisco, MAT_VISCO_read, MAT_VISCO_init_elem_prop &
          , MAT_VISCO_init_elem_work, MAT_VISCO_stress & 
          , MAT_VISCO_mempro, MAT_VISCO_memwrk

contains

!=======================================================================
  logical function MAT_isVisco(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isVisco = MAT_isKind(m,isVisco)
  end function MAT_isVisco

!=======================================================================
  logical function MAT_anyVisco()
  MAT_anyVisco = (isVisco>0)
  end function MAT_anyVisco

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_VISCO
! GROUP  : MATERIALS
! PURPOSE: Set material properties for in viscoelastic medium and allow 
!          attenuation in a certain frequency band
! SYNTAX : &MAT_PLASTIC cp,cs,rho,phi,coh,Tv,e0 /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: QP       [dble][0d0] attenuation quality factor of P wave
! ARG: QS       [dble][0d0] attenuation quality factor of S wave
! ARG: Nbody    [int][0] number of viscoelastic mechanisms
! ARG: f0       [dble][0d0] central frequency
! ARG: fmin     [dble][0d0] minimum frequency
! ARG: fmax     [dble][0d0] maximum frequency
!
! END INPUT BLOCK

  subroutine MAT_VISCO_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: cp,cs,rho,QP,QS,Nbody,f0,fmin,fmax

  NAMELIST / MAT_VISCO / cp,cs,rho,QP,QS,Nbody,f0,fmin,fmax


  cp = 0d0
  cs = 0d0
  rho = 0d0
  QP = 0d0
  QS = 0d0
  Nbody = 0
  f0 = 0d0
  fmin = 0d0
  fmax = 0d0

  read(iin, MAT_VISCO, END=100)
  write(iout,200) cp,cs,rho,QP,QS,Nbody,f0,fmin,fmax

  call MAT_setKind(input,isVisco)

  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'QP',QP)
  call MAT_setProp(input,'QS',QS)
  call MAT_setProp(input,'Nbody',Nbody)
  call MAT_setProp(input,'f0',f0)
  call MAT_setProp(input,'fmin',fmin)
  call MAT_setProp(input,'fmax',fmax)

  return

  100 call IO_abort('MAT_VISCO_read: MAT_VISCO input block not found')
  
  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'P-wave Q . . . . . . .  . . . . . . .(QP) =',EN12.3,/5x, &
    'S-wave Q. . . . . . . . . . . . . . .(QS) =',EN12.3,/5x, &
    'Number of visco mechanisms . . . .(Nbody) =',EN12.3,/5x, &
    'Central frequency . . . . . . . . . .(f0) =',EN12.3,/5x, &
    'Minimum frequency . . . . . . . . .(fmin) =',EN12.3,/5x, &
    'Maximum frequency . . . . . . . . .(fmax) =',EN12.3)  

  end subroutine MAT_VISCO_read

!=======================================================================
  subroutine MAT_VISCO_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:) !don't need eccord in the following lines
  
  integer :: i

  character(len=2) :: ctemp !The number of mechanisms < 100

  double precision :: cp,cs,rho,QP,QS,f0,fmin,fmax,mu_inf,lambda_inf, dNbody

  integer :: Nbody

  double precision, pointer, dimension(:,:) :: theta => null()
 
  double precision, pointer, dimension(:) :: wbody => null()
     
  call MAT_setProp(elem,'cp',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'cs',ecoord,MAT_VISCO_mempro)!
    
  call MAT_setProp(elem,'QP',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'QS',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'Nbody',ecoord,MAT_VISCO_mempro)!
    
  call MAT_setProp(elem,'f0',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'fmin',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'fmax',ecoord,MAT_VISCO_mempro)!
    
  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
    
  call MAT_getProp(QP,elem,'QP')
  call MAT_getProp(QS,elem,'QS')
  
  !dNbody=dble(Nbody)
  call MAT_getProp(dNbody,elem,'Nbody')
    
  call MAT_getProp(f0,elem,'f0')
  call MAT_getProp(fmin,elem,'fmin')
  
  call MAT_getProp(fmax,elem,'fmax')
  
!!!! Pay attention 
  Nbody=int(dNbody)

  allocate(theta(Nbody,3))
  allocate(wbody(Nbody))
  print *,'CHECK3'
  call get_attenuation(theta,wbody,mu_inf,lambda_inf,cp,cs,rho,QP,QS,Nbody,f0,fmin,fmax)
!!!!
  print *,'CHECK4'
!!!! debug
 ! print *, 'Nbody is', Nbody
 ! print *, 'wbody is', wbody
 ! print *, 'mu_inf is', mu_inf
 ! print *, 'lambda_inf is', lambda_inf
 ! print *, 'theta is', theta
!!!!
  do i=1,Nbody
     write(ctemp,'(i2)') i
     !call MAT_setProp(elem,'theta1'//trim(adjustl(ctemp)),theta(i,1),MAT_VISCO_mempro)
     call MAT_setProp(elem,'theta1'//trim(ctemp),theta(i,1),MAT_VISCO_mempro)
     call MAT_setProp(elem,'theta2'//trim(ctemp),theta(i,2),MAT_VISCO_mempro)
     call MAT_setProp(elem,'theta3'//trim(ctemp),theta(i,3),MAT_VISCO_mempro)
     call MAT_setProp(elem,'wbody'//trim(ctemp),wbody(i),MAT_VISCO_mempro)
  end do
      
  call MAT_setProp(elem,'mu',mu_inf,MAT_VISCO_mempro)
  call MAT_setProp(elem,'lambda',lambda_inf,MAT_VISCO_mempro)

  end subroutine MAT_VISCO_init_elem_prop

!=======================================================================
  subroutine MAT_VISCO_init_elem_work(m,p,ngll,dt)

  type(matwrk_visco_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  double precision, intent(in) :: dt
  integer :: i
  double precision :: dNbody
  character(len=2) :: ctemp !The number of mechanisms < 100

  call MAT_getProp(dNbody,p,'Nbody')
  m%Nbody=int(dNbody)

  call MAT_getProp(m%lambda,p,'lambda')
  call MAT_getProp(m%mu,p,'mu')

  allocate(m%theta(m%Nbody,3))
  allocate(m%wbody(m%Nbody))

  do i=1,m%Nbody
     write(ctemp,'(i2)') i
     call MAT_getProp(m%theta(i,1),p,'theta1'//trim(ctemp))
     call MAT_getProp(m%theta(i,2),p,'theta2'//trim(ctemp))
     call MAT_getProp(m%theta(i,3),p,'theta3'//trim(ctemp))
     call MAT_getProp(m%wbody(i),p,'wbody'//trim(ctemp))
  end do

  !Initiate anelastic function el
  allocate(m%el(ngll,ngll,m%Nbody,3))
  m%el = 0d0

  !Initiate old strain at time 0
  allocate(m%etot_old(ngll,ngll,3))
  m%etot_old = 0d0

  MAT_VISCO_memwrk = MAT_VISCO_memwrk+size(transfer(m,(/0d0/)))+size(m%el)

  end subroutine MAT_VISCO_init_elem_work

!============================================================================

! Constitutive law

  subroutine MAT_VISCO_stress(s,etot,m,ngll,dt)

  integer, intent(in) :: ngll
  double precision, intent(in) :: dt
  double precision, intent(in) :: etot(ngll,ngll,3)
  double precision, intent(out):: s(ngll,ngll,3)
  type (matwrk_visco_type), intent(inout) :: m
  
  double precision :: lambda, two_mu, RK_factor
  integer :: i
!  double precision, pointer, dimension(:,:,:) :: el_old
  double precision, dimension(ngll,ngll,3):: s_an(ngll,ngll,3)

  lambda = m%lambda
  two_mu = 2d0*m%mu

!  allocate(el_old(ngll,ngll,m%Nbody,3))
!  el_old=m%el

! Update anelastic function el(m) using el(m-1) and etot(m-1)
  do i=1,m%Nbody
     RK_factor=m%wbody(i)*dt-(m%wbody(i)**2)*(dt**2)/2d0+(m%wbody(i)**3)*(dt**3)/6d0  &
               -(m%wbody(i)**4)*(dt**4)/24d0
     m%el(:,:,i,1)=m%el(:,:,i,1)+RK_factor*(m%etot_old(:,:,1)-m%el(:,:,i,1))
     m%el(:,:,i,2)=m%el(:,:,i,2)+RK_factor*(m%etot_old(:,:,2)-m%el(:,:,i,2))
     m%el(:,:,i,3)=m%el(:,:,i,3)+RK_factor*(m%etot_old(:,:,3)-m%el(:,:,i,3))   
  end do

! Update strain to time step m
  m%etot_old = etot

! Calculate anelastic part of stress at time step m
  s_an=0d0
  
  do i=1,m%Nbody
     s_an(:,:,1)=s_an(:,:,1)+m%theta(i,1)*m%el(:,:,i,1)+m%theta(i,2)*m%el(:,:,i,2)
     s_an(:,:,2)=s_an(:,:,2)+m%theta(i,2)*m%el(:,:,i,1)+m%theta(i,1)*m%el(:,:,i,2)
     s_an(:,:,3)=s_an(:,:,3)+m%theta(i,3)*m%el(:,:,i,3)
  end do

! Calculate the total stress at time step m

  s(:,:,1) = (lambda+two_mu)*etot(:,:,1) + lambda*etot(:,:,2)-s_an(:,:,1)
  s(:,:,2) = lambda*etot(:,:,1) + (lambda+two_mu)*etot(:,:,2)-s_an(:,:,2)
  s(:,:,3) = two_mu*etot(:,:,3)-s_an(:,:,3)

  end subroutine MAT_VISCO_stress

!===============================================================================

!  function MAT_VISCO_export(m) result(dat)
  
!  type(matwrk_visco_type), intent(in) :: m
!  real :: dat(size(m%el,1),size(m%el,2),size(m%el,3),size(m%el,4))

!  dat=real(m%el)

!  end function MAT_VISCO_export


!===============================================================================



!=====================================================================================
end module mat_visco














































































