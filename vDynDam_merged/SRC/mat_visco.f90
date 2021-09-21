module mat_visco
! Added by Yihe Huang (2012)
! viscoelastic medium following Moczo (2004)
! Modified by J.-P. Ampuero (2015)
  
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
          , MAT_isVisco, MAT_VISCO_read, MAT_VISCO_init_elem_prop &
          , MAT_VISCO_init_elem_work, MAT_VISCO_stress & 
          , MAT_VISCO_mempro, MAT_VISCO_memwrk

contains

!=======================================================================
  logical function MAT_isVisco(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isVisco = MAT_isKind(m,isVisco)
  end function MAT_isVisco

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_VISCO
! GROUP  : MATERIALS
! PURPOSE: Set material properties for viscoelastic medium with
!          approximately constant quality factor in a prescribed frequency band
! SYNTAX : &MAT_VISCO cp,cs,rho,QP,QS,Nbody,fmin,fmax /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: QP       [dble][0d0] attenuation quality factor of P waves
! ARG: QS       [dble][0d0] attenuation quality factor of S waves
! ARG: Nbody    [int][0] number of viscoelastic mechanisms
! ARG: fmin     [dble][0d0] minimum frequency for constant Q
! ARG: fmax     [dble][0d0] maximum frequency for constant Q
!
! NOTE   : For Nbody=3, constant Q with less than 5% error can be achieved
!          over a maximum bandwidth fmax/fmin ~ 100
!
! END INPUT BLOCK

  subroutine MAT_VISCO_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: cp,cs,rho,QP,QS,Nbody,fmin,fmax

  NAMELIST / MAT_VISCO / cp,cs,rho,QP,QS,Nbody,fmin,fmax

  call MAT_setKind(input,isVisco)

  cp = 0d0
  cs = 0d0
  rho = 0d0
  QP = 0d0
  QS = 0d0
  Nbody = 0
  fmin = 0d0
  fmax = 0d0

  read(iin, MAT_VISCO, END=100)
  write(iout,200) cp,cs,rho,QP,QS,Nbody,fmin,fmax

  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'QP',QP)
  call MAT_setProp(input,'QS',QS)
  call MAT_setProp(input,'Nbody',Nbody)
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
    'Minimum frequency . . . . . . . . .(fmin) =',EN12.3,/5x, &
    'Maximum frequency . . . . . . . . .(fmax) =',EN12.3)  

  end subroutine MAT_VISCO_read

!=======================================================================
  subroutine MAT_VISCO_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:) !don't need eccord in the following lines
  
  integer :: i,Nbody
  character(len=2) :: ctemp !The number of mechanisms < 100
  double precision :: cp,cs,rho,QP,QS,fmin,fmax,mu_inf,lambda_inf, dNbody
  double precision, dimension(:,:), allocatable :: theta
  double precision, dimension(:), allocatable :: wbody 
     
  call MAT_setProp(elem,'cp',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'cs',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'QP',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'QS',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'Nbody',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'fmin',ecoord,MAT_VISCO_mempro)!
  call MAT_setProp(elem,'fmax',ecoord,MAT_VISCO_mempro)!
    
  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_getProp(QP,elem,'QP')
  call MAT_getProp(QS,elem,'QS')
  call MAT_getProp(dNbody,elem,'Nbody') !dNbody=dble(Nbody)
  Nbody=int(dNbody) !!!! Pay attention 
  call MAT_getProp(fmin,elem,'fmin')
  call MAT_getProp(fmax,elem,'fmax')

  allocate(theta(Nbody,3))
  allocate(wbody(Nbody))
  call get_attenuation(theta,wbody,mu_inf,lambda_inf,cp,cs,rho,QP,QS,Nbody,fmin,fmax)
  do i=1,Nbody
    write(ctemp,'(i2)') i
    !call MAT_setProp(elem,'theta1'//trim(adjustl(ctemp)),theta(i,1),MAT_VISCO_mempro)
    call MAT_setProp(elem,'theta1'//trim(ctemp),theta(i,1),MAT_VISCO_mempro)
    call MAT_setProp(elem,'theta2'//trim(ctemp),theta(i,2),MAT_VISCO_mempro)
    call MAT_setProp(elem,'theta3'//trim(ctemp),theta(i,3),MAT_VISCO_mempro)
    call MAT_setProp(elem,'wbody'//trim(ctemp),wbody(i),MAT_VISCO_mempro)
  end do
  deallocate(theta,wbody)
      
  call MAT_setProp(elem,'mu',mu_inf,MAT_VISCO_mempro)
  call MAT_setProp(elem,'lambda',lambda_inf,MAT_VISCO_mempro)

  end subroutine MAT_VISCO_init_elem_prop

!=======================================================================
  subroutine MAT_VISCO_init_elem_work(m,p,ngll)

  type(matwrk_visco_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
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
  
  double precision, dimension(ngll,ngll,3):: s_an(ngll,ngll,3)
  double precision :: lambda, two_mu, RK_factor
  integer :: i

  lambda = m%lambda
  two_mu = 2d0*m%mu

! Update anelastic function el(m) using el(m-1) and etot(m-1)
  do i=1,m%Nbody
    RK_factor=m%wbody(i)*dt-(m%wbody(i)*dt)**2/2d0+(m%wbody(i)*dt)**3/6d0  &
             -(m%wbody(i)*dt)**4/24d0
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

!=======================================================================
  subroutine get_attenuation(theta,wbody,mu_inf,lambda_inf,cp,cs,rho,QP,QS,Nbody,fmin,fmax)

  use constants, only : PI

  integer, intent(in) :: Nbody
  double precision, intent(in) :: cp, cs, rho, QP, QS, fmin, fmax
  double precision, intent(out) :: theta(Nbody,3), wbody(Nbody)
  double precision, intent(out) :: mu_inf, lambda_inf

  integer :: i,j,Nfrequency
  double precision :: w0, wmin, wmax, &
                      RP1, RP2, RP, RS1, RS2, RS, mu, lambda ! to calculate unrelaxed moduli
  double precision, dimension(2*Nbody-1,Nbody) :: AP, AS  ! (Nfrequency*Nbody) Matrices for P and S
  double precision, dimension(2*Nbody-1) :: QPInv,QSInv,w ! (Nfrequency) 1/Q and frequencies used for inversion
  double precision, dimension(Nbody) :: Y_alpha, Y_beta   ! Anelastic coefficients for P and S
  double precision :: WM(Nbody), VM(Nbody,Nbody)

!Calculate anelastic functions Y_alpha and Y_beta

  Nfrequency = 2*Nbody-1

  w0=2.0d0*PI*(fmin*fmax)**0.5d0
  wmin=2.0d0*PI*fmin
  wmax=2.0d0*PI*fmax

  if(Nbody>1) then
    do i=1,Nfrequency
      w(i)=exp(log(wmin)+(i-1)*(log(wmax)-log(wmin))/(Nfrequency-1))
    end do
  else
    w(:)=w0
  end if
  
  do j=1,Nbody
    wbody(j)=w(2*j-1)
  end do

  QPInv=1.0d0/QP
  QSInv=1.0d0/QS

  do i=1,Nfrequency
    do j=1,Nbody
      AP(i,j)=(wbody(j)*w(i)+wbody(j)**2/QP)/(wbody(j)**2+w(i)**2)
      AS(i,j)=(wbody(j)*w(i)+wbody(j)**2/QS)/(wbody(j)**2+w(i)**2)
    end do
  end do
  
  WM=0d0
  VM=0d0

  call svdcmp(AP,Nfrequency,Nbody,WM,VM)
  call svbksb(AP,WM,VM,Nfrequency,Nbody,QPInv,Y_alpha)

  call svdcmp(AS,Nfrequency,Nbody,WM,VM)
  call svbksb(AS,WM,VM,Nfrequency,Nbody,QSInv,Y_beta)


!Calculate unrelaxed modulus mu_inf and lambda_inf

  RP1=1d0
  RP2=0d0
  RS1=1d0
  RS2=0d0

  do j=1,Nbody
    RP1=RP1-Y_alpha(j)/(1d0+(w0/wbody(j))**2)
    RP2=RP2+Y_alpha(j)*(w0/wbody(j))/(1d0+(w0/wbody(j))**2)
    RS1=RS1-Y_beta(j)/(1d0+(w0/wbody(j))**2)
    RS2=RS2+Y_beta(j)*(w0/wbody(j))/(1d0+(w0/wbody(j))**2)
  end do

  RP=sqrt(RP1**2+RP2**2)
  RS=sqrt(RS1**2+RS2**2)

  mu=rho*cs*cs
  lambda=rho*(cp*cp-2d0*cs*cs)

  mu_inf=mu*(RS+RS1)/(2*RS**2)
  lambda_inf=(lambda+2d0*mu)*(RP+RP1)/(2*RP**2)-2d0*mu_inf


!Calculate Y+, Y- and 2d0*mu*Y_beta

  do j=1,Nbody
    theta(j,1)=(lambda_inf+2d0*mu_inf)*Y_alpha(j)
    theta(j,2)=(lambda_inf+2d0*mu_inf)*Y_alpha(j)-2d0*mu_inf*Y_beta(j)
    theta(j,3)=2d0*mu_inf*Y_beta(j)
  end do

  end subroutine get_attenuation
!=====================================================================================

  SUBROUTINE svbksb(u,w,v,m,n,b,x) 
  INTEGER nmax !Added Yihe
  PARAMETER (nmax=100) 
  INTEGER m,n !Added Yihe
  DOUBLE PRECISION u(m,n),w(n),v(n,n),b(m),x(n),tnmax(nmax),s !changed from b(n) to b(m)
  INTEGER i,j,jj 
  do j=1,n 
    s=0.0 
    if (w(j)/=0.0) then 
      do i=1,m 
        s=s+u(i,j)*b(i) 
      end do 
      s=s/w(j) 
    end if 
    tnmax(j)=s 
  end do 
  do j=1,n 
    s=0.0 
    do jj=1,n 
      s=s+v(j,jj)*tnmax(jj) 
    end do 
    x(j)=s 
  end do 
  END SUBROUTINE svbksb 

!====================================================================

  SUBROUTINE svdcmp(a,m,n,w,v) 
  INTEGER nmax
  PARAMETER(nmax=100)
  INTEGER m,n !Added Yihe 
  DOUBLE PRECISION a(m,n),w(n),v(n,n),rv1(nmax) 
  INTEGER i,j,l,k,nm,its 
  DOUBLE PRECISION g,s,scale1,anorm,f,h,c,x,y,z  !need to add pythag if used
   !if(m<n) pause 'you must augment a with extra zero rows.' 
   ! m > or = n
  g=0.0 
  scale1=0.0 
  anorm=0.0 
  do i=1,n   !25
    l=i+1 
    rv1(i)=scale1*g 
    g=0.0 
    s=0.0 
    scale1=0.0 
    if(i.le.m) then !changed by Yihe
      do k=i,m     !11
        scale1=scale1+abs(a(k,i)) 
      end do       !11
      if(scale1/=0.0) then 
        do k=i,m     !12
          a(k,i)=a(k,i)/scale1 
          s=s+a(k,i)*a(k,i) 
        end do       !12
        f=a(i,i) 
        g=-sign(sqrt(s),f) 
        h=f*g-s 
        a(i,i)=f-g 
        if (i/=n) then    !in book version of Recipe
          do j=l,n     !15
            s=0.0 
            do k=i,m !13
              s=s+a(k,i)*a(k,j) 
            end do !13
            f=s/h 
            do k=i,m !14
              a(k,j)=a(k,j)+f*a(k,i) 
            end do !14
          end do !15
        endif      !in book version of Recipe
        do k=i,m !16
          a(k,i)=scale1*a(k,i) 
        end do !16
      endif 
    endif 
    w(i)=scale1*g 
    g=0.0 
    s=0.0 
    scale1=0.0 
    if ((i.le.m).and.(i/=n)) then !changed by Yihe
      do k=l,n !17
        scale1=scale1+abs(a(i,k)) 
      end do !17
      if (scale1/=0.0) then 
        do k=l,n !18
          a(i,k)=a(i,k)/scale1 
          s=s+a(i,k)*a(i,k) 
        end do !18
        f=a(i,l) 
        g=-sign(sqrt(s),f) 
        h=f*g-s 
        a(i,l)=f-g 
        do k=l,n !19
          rv1(k)=a(i,k)/h 
        end do !19
        if (i/=m) then      !in book version of Recipe
          do j=l,m !23
            s=0.0 
            do k=l,n !21
              s=s+a(j,k)*a(i,k) 
            end do !21
            do k=l,n !22
              a(j,k)=a(j,k)+s*rv1(k) 
            end do !22
          end do !23
        endif      !in book version of Recipe
        do k=l,n !24
          a(i,k)=scale1*a(i,k) 
        end do !24
      endif 
    endif 
    anorm=max(anorm,(abs(w(i))+abs(rv1(i)))) 
  end do !25
  do i=n,1,-1 !32
    if (i<n) then 
      if (g/=0.0) then 
        do j=l,n !26
          v(j,i)=(a(i,j)/a(i,l))/g 
        end do !26
        do j=l,n !29
          s=0.0 
          do k=l,n !27
             s=s+a(i,k)*v(k,j) 
          end do !27
          do k=l,n !28
            v(k,j)=v(k,j)+s*v(k,i) 
          end do !28
        end do !29
      endif 
      do j=l,n !31
        v(i,j)=0.0 
        v(j,i)=0.0 
      end do !31
    endif 
    v(i,i)=1.0 
    g=rv1(i) 
    l=i 
  end do !32
  do i=n,1,-1 !39 assume m>n, or otherwise i=min(m,n),1,-1 
    l=i+1 
    g=w(i) 
    if (i<n) then      !in book version of Recipe
      do j=l,n !33
        a(i,j)=0.0 
      end do !33
    endif         !in book version of Recipe
    if (g/=0.0) then 
      g=1.0/g 
      if (i/=n) then       !in book version of Recipe
        do j=l,n !36
          s=0.0 
          do k=l,m !34
            s=s+a(k,i)*a(k,j) 
          end do !34
          f=(s/a(i,i))*g 
          do k=i,m !35
            a(k,j)=a(k,j)+f*a(k,i) 
          end do !35
        end do !36
      endif        !in book version of Recipe 
      do j=i,m !37
        a(j,i)=a(j,i)*g 
      end do !37
    else 
      do j=i,m !38
        a(j,i)=0.0 
      end do !38
    endif 
    a(i,i)=a(i,i)+1.0 
  end do !39
  do k=n,1,-1 !49
    do its=1,30 !48
      do l=k,1,-1 !41
        nm=l-1 
        if((abs(rv1(l))+anorm)==anorm) exit !in recipe: goto 2
        if((abs(w(nm))+anorm)==anorm) exit !in recipe: goto 1
      end do !41
      if (abs(rv1(l))+anorm/=anorm) then !!!!! not in recipe
!!!!!!!!!!!!!! Below in recipe: 1
        c=0.0 
        s=1.0 
        do i=l,k !43
          f=s*rv1(i)
! not in book version of recipe: rv1(i)=c*rv1(i) !added by Yihe 
          if((abs(f)+anorm)/=anorm) then     !in book version of Recipe 
            g=w(i) 
            h=sqrt(f*f+g*g)     !in book version of Recipe 
            w(i)=h 
            h=1.0/h 
            c=(g*h) 
            s=-(f*h) 
            do j=1,m !42
              y=a(j,nm) 
              z=a(j,i) 
              a(j,nm)=(y*c)+(z*s) 
              a(j,i)=-(y*s)+(z*c) 
            end do !42
          endif      !in book version of Recipe
        end do !43
      end if  !!!!! not in recipe
!!!!!!!!!!!!!! Below in recipe: 2
      z=w(k) 
      if(l==k) then 
        if(z<0.0) then 
          w(k)=-z 
          do j=1,n !44
            v(j,k)=-v(j,k) 
          end do !44
        endif 
      exit ! in recipe: goto 3 (exit do 48)
      endif 
      if(its==30) stop 'no convergence in 30 iterations' 
      x=w(l) 
      nm=k-1 
      y=w(nm) 
      g=rv1(nm) 
      h=rv1(k) 
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y) 
      g=sqrt(f*f+1.0) 
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x 
      c=1.0 
      s=1.0 
      do j=l,nm !47
        i=j+1 
        g=rv1(i) 
        y=w(i) 
        h=s*g 
        g=c*g 
        z=sqrt(f*f+h*h) 
        rv1(j)=z 
        c=f/z 
        s=h/z 
        f=(x*c)+(g*s) 
        g=-(x*s)+(g*c) 
        h=y*s 
        y=y*c 
        do nm=1,n !45 
          x=v(nm,j) 
          z=v(nm,i) 
          v(nm,j)=(x*c)+(z*s) 
          v(nm,i)=-(x*s)+(z*c) 
        end do !45
        z=sqrt(f*f+h*h) 
        w(j)=z 
        if(z/=0.0) then 
          z=1.0/z 
          c=f*z 
          s=h*z 
        endif 
        f=(c*g)+(s*y) 
        x=-(s*g)+(c*y) 
        do nm=1,m !46 
          y=a(nm,j) 
          z=a(nm,i) 
          a(nm,j)=(y*c)+(z*s) 
          a(nm,i)=-(y*s)+(z*c) 
        end do !46
      end do !47
      rv1(l)=0.0 
      rv1(k)=f 
      w(k)=x 
    end do !48
! in recipe: 3 continue
  end do !49
  END SUBROUTINE svdcmp 

end module mat_visco
