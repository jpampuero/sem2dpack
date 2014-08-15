module attenuation
! Added by Yihe (2012)

  implicit none
  contains


  subroutine get_attenuation(theta,wbody,mu_inf,lambda_inf,cp,cs,rho,QP,QS,Nbody,f0,fmin,fmax)

  use constants, only : PI

  implicit none

  integer :: i,j,Nbody,Nfrequency

  double precision :: mu_inf, lambda_inf, cp, cs, rho, QP, QS, f0, fmin, fmax, w0, wmin, wmax

  double precision :: theta(Nbody,3), wbody(Nbody)

  double precision, pointer, dimension(:,:) :: AP => null(), AS => null() ! (Nfrequency*Nbody) Matrices for P and S

  double precision, pointer, dimension(:) :: QPInv => null(), QSInv => null() ! (Nfrequency) 1/Q

  double precision, pointer, dimension(:) :: w => null() ! (Nfrequency) frequencies used for inversion

  double precision :: Y_alpha(Nbody), Y_beta(Nbody) ! Anelastic function for P and S

  double precision, pointer, dimension(:) :: WM => null()

  double precision, pointer, dimension(:,:) :: VM => null()

  double precision :: RP1, RP2, RP, RS1, RS2, RS, mu, lambda ! Used to calculate unrelaxted modulus

  intent(in) :: cp, cs, rho, QP, QS, Nbody, f0, fmin, fmax

  intent(out) :: theta, wbody, mu_inf, lambda_inf

!Calculate anelastic functions Y_alpha and Y_beta

  Nfrequency = 2*Nbody-1

  allocate(AP(Nfrequency,Nbody))
  allocate(AS(Nfrequency,Nbody))
  allocate(QPInv(Nfrequency))
  allocate(QSInv(Nfrequency))
  allocate(w(Nfrequency))

  w0=2.0d0*PI*f0
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

  do i=1,Nfrequency
     QPInv(i)=1.0d0/QP
     QSInv(i)=1.0d0/QS
     do j=1,Nbody
        AP(i,j)=(wbody(j)*w(i)+wbody(j)**2/QP)/(wbody(j)**2+w(i)**2)
        AS(i,j)=(wbody(j)*w(i)+wbody(j)**2/QS)/(wbody(j)**2+w(i)**2)
     end do
  end do
  
  allocate(WM(Nbody))
  allocate(VM(Nbody,Nbody))

  WM=0d0
  VM=0d0

  call svdcmp(AP,Nfrequency,Nbody,WM,VM)
  call svbksb(AP,WM,VM,Nfrequency,Nbody,QPInv,Y_alpha)

  call svdcmp(AS,Nfrequency,Nbody,WM,VM)
  call svbksb(AS,WM,VM,Nfrequency,Nbody,QSInv,Y_beta)

!========================================================================================

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

!=====================================================================================

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
      if(its==30) pause 'no convergence in 30 iterations' 
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

end module attenuation
