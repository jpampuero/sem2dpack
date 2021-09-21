module gll

  implicit none
  private

  interface print_GLL
    module procedure print_GLL_a, print_GLL_b
  end interface print_GLL

  public :: get_GLL_info, print_GLL, hgll,hdgll,zwgljd

contains

!=======================================================================
! 
! Get coordinates and weights of the Gauss-Lobatto-Legendre points
! and derivatives of Lagrange polynomials  H_ij = h'_i(xgll(j))

  subroutine get_GLL_info(n,x,w,H)

  integer, intent(in) :: n
  double precision, intent(out) :: w(n),x(n),H(n,n)

  integer :: ix,ip

  call zwgljd(x,w,n,0.d0,0.d0)
  !if n is odd make sure the middle point is exactly zero:
  if(mod(n,2) /= 0) x((n-1)/2+1) = 0.d0

  do ix=1,n
  do ip=1,n
    H(ip,ix)  = hdgll(ip-1,ix-1,x,n) ! in hdggl, indices start at 0
  enddo
  enddo

  end subroutine get_GLL_info
!
!=======================================================================
!
! Print GLL information to a text file

  subroutine print_GLL_a(n,x,w,H)

  use stdio, only: IO_new_unit

  integer, intent(in) :: n
  double precision, intent(in) :: w(n),x(n),H(n,n)

  integer :: ip,iout
  character(15) :: fmt

 ! fmt = (n(x,F0.16))
  write(fmt,'("(",I0,"(x,F0.16))")') n

  iout = IO_new_unit()
  open(unit=iout,file='gll_sem2d.tab')
  write(iout,fmt) x
  write(iout,fmt) w
  do ip=1,n
    write(iout,fmt) H(ip,:)
  enddo
  close(iout)

  end subroutine print_GLL_a

!-----------------------------------------------------------------------

  subroutine print_GLL_b(n)

  integer, intent(in) :: n

  double precision :: w(n),x(n),H(n,n)

  call get_GLL_info(n,x,w,H)
  call print_GLL_a(n,x,w,H)

  end subroutine print_GLL_b

      
!=======================================================================
!
!     E n d w 1 :
!     ---------
!
!=======================================================================
!
  double precision function endw1 (n,alpha,beta)
  

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0, &
        three=3.d0,four=4.d0

  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3

  integer i
!
!-----------------------------------------------------------------------
!
  f3 = zero
  apb   = alpha+beta
  if (n == 0) then
   endw1 = zero
   return
  endif
  f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw1 = f1
   return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw1 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw1  = f3

  return
  end function endw1

  double precision function endw2 (n,alpha,beta)
!
!=======================================================================
!
!     E n d w 2 :
!     ---------
!
!=======================================================================
!
  

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0, &
      three=3.d0,four=4.d0

  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3

  integer i

!
!-----------------------------------------------------------------------
!
  apb   = alpha+beta
  f3 = zero
  if (n == 0) then
   endw2 = zero
   return
  endif
  f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw2 = f1
   return
  endif
  fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw2 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw2  = f3
  return
  end function endw2

  double precision function gammaf (x)
!
!=======================================================================
!
!     G a m m a f :
!     -----------
!
!=======================================================================
!

  double precision x
  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: half=0.5d0,one=1.d0,two=2.d0

  gammaf = one

  if (x == -half) gammaf = -two*sqrt(pi)
  if (x ==  half) gammaf =  sqrt(pi)
  if (x ==  one ) gammaf =  one
  if (x ==  two ) gammaf =  one
  if (x ==  1.5d0) gammaf =  sqrt(pi)/2.d0
  if (x ==  2.5d0) gammaf =  1.5d0*sqrt(pi)/2.d0
  if (x ==  3.5d0) gammaf =  2.5d0*1.5d0*sqrt(pi)/2.d0
  if (x ==  3.d0 ) gammaf =  2.d0
  if (x ==  4.d0 ) gammaf = 6.d0
  if (x ==  5.d0 ) gammaf = 24.d0
  if (x ==  6.d0 ) gammaf = 120.d0

  return
  end function gammaf

  double precision FUNCTION HDGLL (I,j,ZGLL,NZ)
!-------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrangian interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j).
!
!-------------------------------------------------------------

  

  integer i,j,nz
  double precision zgll(0:nz-1)

  integer idegpoly
  double precision rlegendre1,rlegendre2,rlegendre3

    idegpoly = nz - 1
  if ((i == 0).and.(j == 0)) then
          hdgll = - dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
  else if ((i == idegpoly).and.(j == idegpoly)) then
          hdgll = dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
  else if (i == j) then
          hdgll = 0.d0
  else
         rlegendre1 = pnleg(zgll(j),idegpoly)
         rlegendre2 = pndleg(zgll(j),idegpoly)
         rlegendre3 = pnleg(zgll(i),idegpoly)
  hdgll = rlegendre1 / (rlegendre3*(zgll(j)-zgll(i))) &
    + (1.d0-zgll(j)*zgll(j))*rlegendre2/(dble(idegpoly)* &
    (dble(idegpoly)+1.d0)*rlegendre3* &
    (zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  return
  end FUNCTION hdgll

!=====================================================================

  double precision FUNCTION HGLL (I,Z,ZGLL,NZ)
!-------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant L through
!     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
!
!-------------------------------------------------------------

  

  integer i,nz
  double precision z
  double precision ZGLL(0:nz-1)

  integer n
  double precision EPS,DZ,ALFAN

  EPS = 1.d-5
  DZ = Z - ZGLL(I)
  IF (ABS(DZ)  <  EPS) THEN
   HGLL = 1.d0
   RETURN
  ENDIF
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+1.d0)
  HGLL = - (1.d0-Z*Z)*PNDLEG(Z,N)/ (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))
  RETURN
  end function hgll
!
!=====================================================================

  subroutine jacg (xjac,np,alpha,beta)
!
!=======================================================================
!
!     J a c g : Compute np Gauss points, which are the zeros of the
!               Jacobi polynomial with parameter alpha and beta.
!
!=======================================================================
!
!     Note :
!     ----
!          .Alpha and Beta determines the specific type of gauss points.
!                  .alpha = beta =  0.0  ->  Legendre points
!                  .alpha = beta = -0.5  ->  Chebyshev points
!
!=======================================================================
!
  

  integer np
  double precision alpha,beta
  double precision xjac(np)

  integer k,j,i,jmin,jm,n
  double precision xlast,dth,x,x1,x2,recsum,delx,xmin,swap
  double precision p,pd,pm1,pdm1,pm2,pdm2

  integer, parameter :: kstop = 10
  double precision, parameter :: zero = 0.d0, eps = 1.0d-12

!
!-----------------------------------------------------------------------
!
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  xlast = 0.d0
  n   = np-1
  dth = 4.d0*datan(1.d0)/(2.d0*dble(n)+2.d0)
  p = 0.d0
  pd = 0.d0
  jmin = 0
  do 40 j=1,np
   if (j == 1) then
      x = cos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
   else
      x1 = cos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
      x2 = xlast
      x  = (x1+x2)/2.d0
   endif
   do 30 k=1,kstop
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
      recsum = 0.d0
      jm = j-1
      do 29 i=1,jm
         recsum = recsum+1.d0/(x-xjac(np-i+1))
 29         continue
      delx = -p/(pd-recsum*p)
      x    = x+delx
      if (abs(delx) < eps) goto 31
 30      continue
 31      continue
   xjac(np-j+1) = x
   xlast        = x
 40   continue
  do 200 i=1,np
   xmin = 2.d0
   do 100 j=i,np
      if (xjac(j) < xmin) then
         xmin = xjac(j)
         jmin = j
      endif
 100     continue
   if (jmin /= i) then
      swap = xjac(i)
      xjac(i) = xjac(jmin)
      xjac(jmin) = swap
   endif
 200  continue
  return
  end subroutine jacg
!
!=====================================================================

  subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp, &
                     bet,x)
!
!=======================================================================
!
!     J a c o b f : Compute the Jacobi polynomial and its derivative
!     -----------   of degree n at x.
!
!=======================================================================
  

  double precision poly,pder,polym1,pderm1,polym2,pderm2,alp,bet,x
  integer n

  double precision apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave
  integer k

  apb  = alp+bet
  poly = 1.d0
  pder = 0.d0
  psave = 0.d0
  pdsave = 0.d0

  if (n  ==  0) return

  polyl = poly
  pderl = pder
  poly  = (alp-bet+(apb+2.d0)*x)/2.d0
  pder  = (apb+2.d0)/2.d0
  if (n  ==  1) return
  do 20 k=2,n
   dk = dble(k)
   a1 = 2.d0*dk*(dk+apb)*(2.d0*dk+apb-2.d0)
   a2 = (2.d0*dk+apb-1.d0)*(alp**2-bet**2)
   b3 = (2.d0*dk+apb-2.d0)
   a3 = b3*(b3+1.d0)*(b3+2.d0)
   a4 = 2.d0*(dk+alp-1.d0)*(dk+bet-1.d0)*(2.d0*dk+apb)
   polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
   pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
   psave  = polyl
   pdsave = pderl
   polyl  = poly
   poly   = polyn
   pderl  = pder
   pder   = pdern
 20   continue
  polym1 = polyl
  pderm1 = pderl
  polym2 = psave
  pderm2 = pdsave
!
  return
  end subroutine jacobf
!
!=====================================================================

  double precision FUNCTION PNDLEG (Z,N)
!-------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!-------------------------------------------------------------
  

  double precision z
  integer n

  double precision P1,P2,P1D,P2D,P3D,FK,P3
  integer k

  P1   = 1.d0
  P2   = Z
  P1D  = 0.d0
  P2D  = 1.d0
  P3D  = 1.d0
  DO 10 K = 1, N-1
   FK  = dble(K)
   P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
   P3D = ((2.d0*FK+1.d0)*P2 + (2.d0*FK+1.d0)*Z*P2D - FK*P1D) &
                          /(FK+1.d0)
   P1  = P2
   P2  = P3
   P1D = P2D
   P2D = P3D
 10   CONTINUE
  PNDLEG = P3D
  RETURN
  end function pndleg
!
!=====================================================================

  double precision FUNCTION PNLEG (Z,N)
!-------------------------------------------------------------
!
!     Compute the value of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!-------------------------------------------------------------
  

  double precision z
  integer n

  double precision P1,P2,P3,FK
  integer k

  P1   = 1.d0
  P2   = Z
  P3   = P2
  DO 10 K = 1, N-1
   FK  = dble(K)
   P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
   P1  = P2
   P2  = P3
 10   CONTINUE
  PNLEG = P3
  RETURN
  end function pnleg
!
!=====================================================================

  double precision function pnormj (n,alpha,beta)
!
!=======================================================================
!
!     P n o r m j
!     -----------
!
!=======================================================================
!
  

  double precision alpha,beta
  integer n

  double precision one,two,dn,const,prod,dindx,frac
  integer i

  one   = 1.d0
  two   = 2.d0
  dn    = dble(n)
  const = alpha+beta+one
  if (n <= 1) then
   prod   = gammaf(dn+alpha)*gammaf(dn+beta)
   prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
   pnormj = prod * two**const/(two*dn+const)
   return
  endif
  prod  = gammaf(alpha+one)*gammaf(beta+one)
  prod  = prod/(two*(one+const)*gammaf(const+one))
  prod  = prod*(one+alpha)*(two+alpha)
  prod  = prod*(one+beta)*(two+beta)
  do 100 i=3,n
   dindx = dble(i)
   frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
   prod  = prod*frac
 100  continue
  pnormj = prod * two**const/(two*dn+const)

  return
  end function pnormj
!
!=====================================================================

  subroutine zwgjd(z,w,np,alpha,beta)
!
!=======================================================================
!
!     Z w g j d : Generate np Gauss-Jacobi points and weights
!                 associated with Jacobi polynomial of degree n = np-1
!
!=======================================================================
!
!     Note : Coefficients alpha and beta must be greater than -1.
!     ----
!
!=======================================================================
!
  

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision z(np),w(np)
  double precision alpha,beta

  integer n,np1,np2,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef
!
!-----------------------------------------------------------------------
!
  pd = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n    = np-1
  apb  = alpha+beta
  p    = zero
  pdm1 = zero

  if (np <= 0) stop 'Minimum number of Gauss points is 1'

  if ((alpha <= -one).or.(beta <= -one)) &
    stop 'Alpha and Beta must be greater than -1'

  if (np == 1) then
   z(1) = (beta-alpha)/(apb+two)
   w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
   return
  endif

  call jacg (z,np,alpha,beta)

  np1   = n+1
  np2   = n+2
  dnp1  = dble(np1)
  dnp2  = dble(np2)
  fac1  = dnp1+alpha+beta+one
  fac2  = fac1+dnp1
  fac3  = fac2+one
  fnorm = pnormj(np1,alpha,beta)
  rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
  do i=1,np
    call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
    w(i) = -rcoef/(p*pdm1)
  enddo

  return
  end subroutine zwgjd
!
!=====================================================================

  subroutine zwgljd (z,w,np,alpha,beta)
!
!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!=======================================================================
!
!     Note : alpha and beta coefficients must be greater than -1.
!     ----
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================
!
  

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision alpha,beta
  double precision z(np), w(np)

  integer n,nm1,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision alpg,betg
!
!-----------------------------------------------------------------------
!
  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero

  if (np <= 1) stop 'Minimum number of Gauss-Lobatto points is 2'

  if ((alpha <= -one).or.(beta <= -one)) &
    stop 'Alpha and Beta must be greater than -1'

  if (nm1 > 0) then
   alpg  = alpha+one
   betg  = beta+one
   call zwgjd (z(2),w(2),nm1,alpg,betg)
  endif
  z(1)  = - one
  z(np) =  one
  do  110 i=2,np-1
   w(i) = w(i)/(one-z(i)**2)
  110 continue
  call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1 (n,alpha,beta)/(two*pd)
  call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2 (n,alpha,beta)/(two*pd)

  return
  end subroutine zwgljd

end module gll
