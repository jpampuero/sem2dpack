! This is a set of basic tools, adapted from different open sources
module utils

  implicit none
  private

  interface setopt
    module procedure double_setopt, integer_setopt, string_setopt
  end interface setopt

  public :: setopt, invert2, unique, drank, dsort, hunt &
          , sub2ind, sub2ind_b, spline, splint, heaviside, positive_part

contains

!=====================================================================
elemental function heaviside(x)

  double precision, intent(in) :: x
  double precision :: heaviside

  if (x>=0d0) then
    heaviside = 1d0
  else
    heaviside = 0d0
  endif

end function heaviside

!=====================================================================
elemental function positive_part(x)

  double precision, intent(in) :: x
  double precision :: positive_part

  if (x>=0d0) then
    positive_part = x
  else
    positive_part = 0d0
  endif

end function positive_part

!=====================================================================
subroutine double_setopt(value,deflt,opt)

  double precision, intent(out) :: value
  double precision, intent(in) :: deflt
  double precision, optional :: opt

  if( present(opt)) then
    value=opt
  else
    value=deflt
  endif

end subroutine double_setopt

!----------------------------------------------------------------------
subroutine integer_setopt(value,deflt,opt)

  integer, intent(out) :: value
  integer, intent(in) :: deflt
  integer, optional :: opt
  
  if( present(opt)) then
    value=opt
  else
    value=deflt
  endif

end subroutine integer_setopt

!----------------------------------------------------------------------
subroutine string_setopt(value,deflt,opt)

  character(*), intent(out) :: value
  character(*), intent(in) :: deflt
  character(*), optional :: opt
  
  if( present(opt)) then
    value=opt
  else
    value=deflt
  endif

end subroutine string_setopt


!=====================================================================
!-- inversion of a 2 x 2 matrix
function invert2(A) result(B)

  use stdio, only : IO_abort

  double precision, intent(in) :: A(2,2)
  double precision :: B(2,2)

  double precision :: determinant

  determinant = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if (determinant <= 0d0) call IO_abort('SE_InverseJacobian: undefined Jacobian')

  B(1,1) =   A(2,2)
  B(2,1) = - A(2,1)
  B(1,2) = - A(1,2)
  B(2,2) =   A(1,1)

  B = B/determinant

end function invert2

!=====================================================================
! subscript (i,j) to index (iglob)
! numbering conventions :
!   i=1:n 
!   j=1:p
!   iglob = 1:

integer function sub2ind(i,j,n)
  integer, intent(in) :: i,j,n
  sub2ind = (j-1)*n + i
end function sub2ind

!--------------------------------------------------------------------- 
! returns 0 if out of box
integer function sub2ind_b(i,j,n,p)
  integer, intent(in) :: i,j,n,p
  if (i<=n .and. i>=1 .and. j<=p .and. j>=1) then 
    sub2ind_b= (j-1)*n + i
  else
    sub2ind_b=0 
  endif
end function sub2ind_b

!-----------------------------------------------------------------------
function unique(k) result(ku)

integer, intent(in) :: k(:)
integer, pointer :: ku(:)

integer :: ks(size(k)),n,p,i

call drank(dble(k),ks)
ks = k(ks)

n = 1
p = ks(1)
do i=2,size(ks)
  if (ks(i)>p) then
    n=n+1
    p = ks(i)
  endif
enddo

allocate( ku(n) )
n = 1
ku(n) = ks(1)
do i=2,size(ks)
  if (ks(i)>ku(n)) then
    n = n+1
    ku(n) = ks(i)
  endif
enddo

end function unique


!-----------------------------------------------------------------------
!From Michel Olagnon's ORDERPACK 2.0
!http://www.fortran-2000.com/rank/index.html
!
!Subroutine D_mrgrnk (XDONT, IRNGT)
Subroutine drank(XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
!      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      double precision, Dimension (:), Intent (In) :: XDONT !jpampuero
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
!      Real (kind=kdp) :: XVALA, XVALB
      double precision :: XVALA, XVALB !jpampuero
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
!End Subroutine D_mrgrnk
End Subroutine drank

!-----------------------------------------------------------------------
! DSORT adapted from SLATEC
! Converted to f90 using f2f.perl
! Modified:     + to be used always with DY
!               + only increasing order
!               + no checks / error messages
! Jean-Paul Ampuero     dim oct 19 17:15:29 EDT 2003


! ECK DSORT
!    SUBROUTINE DSORT (DX, DY, N, KFLAG)
    SUBROUTINE DSORT (DX, DY)
!***BEGIN PROLOGUE  DSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
! an auxiliary array.  The array may be sorted in increasing
! or decreasing order.  A slightly modified QUICKSORT
! algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
! Wisniewski, J. A., (SNLA)
!***DESCRIPTION

! DSORT sorts array DX and optionally makes the same interchanges in
! array DY.  The array DX may be sorted in increasing order or
! decreasing order.  A slightly modified quicksort algorithm is used.

! Description of Parameters
! DX - array of values to be sorted   (usually abscissas)
! DY - array to be (optionally) carried along
! N  - number of values in array DX to be sorted
! KFLAG - control parameter
! =  2  means sort DX in increasing order and carry DY along.
! =  1  means sort DX in increasing order (ignoring DY)
! = -1  means sort DX in decreasing order (ignoring DY)
! = -2  means sort DX in decreasing order and carry DY along.

!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
! for sorting with minimal storage, Communications of
! the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
! 761101  DATE WRITTEN
! 761118  Modified to use the Singleton quicksort algorithm.  (JAW)
! 890531  Changed all specific intrinsics to generic.  (WRB)
! 890831  Modified array declarations.  (WRB)
! 891009  Removed unreferenced statement labels.  (WRB)
! 891024  Changed category.  (WRB)
! 891024  REVISION DATE from Version 3.2
! 891214  Prologue converted to Version 4.0 format.  (BAB)
! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
! 901012  Declared all variables; changed X,Y to DX,DY; changed
! code to parallel SSORT. (M. McClain)
! 920501  Reformatted the REFERENCES section.  (DWL, WRB)
! 920519  Clarified error messages.  (DWL)
! 920801  Declarations section rebuilt and code restructured to use
! IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  DSORT
! .. Array Arguments ..
    double precision :: DX(:), DY(:)
! .. Local Scalars ..
    double precision :: R, T, TT, TTY, TY
    INTEGER :: I, IJ, J, K, L, M, NN
! .. Local Arrays ..
    INTEGER :: IL(21), IU(21)
! .. External Subroutines ..
!    EXTERNAL XERMSG
! .. Intrinsic Functions ..
!    INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  DSORT
    NN = size(DX)

! Alter array DX to get decreasing order if needed
!
!    IF (KFLAG <= -1) THEN
!        DO 10 I=1,NN
!            DX(I) = -DX(I)
!        10 end do
!    ENDIF

! Sort DX and carry DY along

    M = 1
    I = 1
    J = NN
    R = 0.375D0

    110 IF (I == J) GO TO 150
    IF (R <= 0.5898437D0) THEN
        R = R+3.90625D-2
    ELSE
        R = R-0.21875D0
    ENDIF

    120 K = I

! Select a central element of the array and save it in location T

    IJ = I + INT((J-I)*R)
    T = DX(IJ)
    TY = DY(IJ)

! If first element of array is greater than T, interchange with T

    IF (DX(I) > T) THEN
        DX(IJ) = DX(I)
        DX(I) = T
        T = DX(IJ)
        DY(IJ) = DY(I)
        DY(I) = TY
        TY = DY(IJ)
    ENDIF
    L = J

! If last element of array is less than T, interchange with T

    IF (DX(J) < T) THEN
        DX(IJ) = DX(J)
        DX(J) = T
        T = DX(IJ)
        DY(IJ) = DY(J)
        DY(J) = TY
        TY = DY(IJ)
    
    ! If first element of array is greater than T, interchange with T
    
        IF (DX(I) > T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
        ENDIF
    ENDIF

! Find an element in the second half of the array which is smaller
! than T

    130 L = L-1
    IF (DX(L) > T) GO TO 130

! Find an element in the first half of the array which is greater
! than T

    140 K = K+1
    IF (DX(K) < T) GO TO 140

! Interchange these elements

    IF (K <= L) THEN
        TT = DX(L)
        DX(L) = DX(K)
        DX(K) = TT
        TTY = DY(L)
        DY(L) = DY(K)
        DY(K) = TTY
        GO TO 130
    ENDIF

! Save upper and lower subscripts of the array yet to be sorted

    IF (L-I > J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M+1
    ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M+1
    ENDIF
    GO TO 160

! Begin again on another portion of the unsorted array

    150 M = M-1
    IF (M == 0) GO TO 190
    I = IL(M)
    J = IU(M)

    160 IF (J-I >= 1) GO TO 120
    IF (I == 1) GO TO 110
    I = I-1

    170 I = I+1
    IF (I == J) GO TO 150
    T = DX(I+1)
    TY = DY(I+1)
    IF (DX(I) <= T) GO TO 170
    K = I

    180 DX(K+1) = DX(K)
    DY(K+1) = DY(K)
    K = K-1
    IF (T < DX(K)) GO TO 180
    DX(K+1) = T
    DY(K+1) = TY
    GO TO 170

! Clean up

!    190 IF (KFLAG <= -1) THEN
!        DO 200 I=1,NN
!            DX(I) = -DX(I)
!        200 end do
!    ENDIF
190    RETURN
    end SUBROUTINE DSORT

!-----------------------------------------------------------------------
!  HUNT is adapted from Numerical Recipes for double precision, F90
!
!  jlo = index in xx right before x
!      = 0 or length(xx) if x is out of xx bounds 
!
subroutine hunt(xx,x,jlo)
  integer, intent(inout) :: jlo
  double precision, intent(in) :: x,xx(:)
  integer :: inc,jhi,jm,n
  logical :: ascnd
  n = size(xx)
  ascnd=xx(n)>xx(1)
  if(jlo<=0.or.jlo>n)then
    jlo=0
    jhi=n+1
  else
   inc=1
   if(x>=xx(jlo).eqv.ascnd)then
     do
        jhi=jlo+inc
        if(jhi>n)then
          jhi=n+1
          exit
        else if(x>=xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
        else
          exit
        endif
     enddo
   else
    jhi=jlo
    do
        jlo=jhi-inc
        if(jlo<1)then
          jlo=0
          exit
        else if(x<xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
        else
          exit
        endif
    enddo
   endif
  endif
  do while(jhi-jlo/=1)
    jm=(jhi+jlo)/2
    if(x>xx(jm).eqv.ascnd)then
        jlo=jm
    else
        jhi=jm
    endif
  enddo
 
end subroutine hunt

!***********************************************************************
!  SPLINE INTERPOLATION
!  adapted from Numerical Recipes for double precision, Fortran 90

  SUBROUTINE spline(x,y,n,yp1,ypn,y2)

  INTEGER, intent(in) :: n
  double precision, intent(in) :: yp1,ypn,x(n),y(n)
  double precision, intent(out) :: y2(n)

  INTEGER i,k
  double precision :: p,qn,sig,un
  double precision :: u(n)

  if (yp1 > .99d30) then
    y2(1)=0d0
    u(1)=0d0
  else
    y2(1)=-0.5d0
    u(1)=(3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2d0
    y2(i)=(sig-1d0)/p
    u(i)=(6d0*((y(i+1)-y(i)) &
         /(x(i+1)-x(i))-(y(i)-y(i-1)) &
         /(x(i)-x(i-1))) &
         /(x(i+1)-x(i-1)) &
         -sig*u(i-1))/p
  enddo
  if (ypn > .99e30) then
    qn=0d0
    un=0d0
  else
    qn=0.5d0
    un=(3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1d0)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  
  END SUBROUTINE spline

!---------------------------------------------------------------
!  adapted from Numerical Recipes for double precision, Fortran 90
!  and to speed up table lookup

  SUBROUTINE splint(xa,ya,y2a,n,x,y,dydx)

  INTEGER, intent(in) :: n
  double precision, intent(in) :: x,xa(n),y2a(n),ya(n)
  double precision, intent(out) :: y
  double precision, intent(out), optional :: dydx

!  integer :: k
  integer :: khi
  integer, save :: klo=1
  double precision :: a,b,h

! original : find bracketing indices bissection
!  klo=1
!  khi=n
!  do while (khi-klo > 1)
!    k=(khi+klo)/2
!    if(xa(k).gt.x)then
!      khi=k
!    else
!      klo=k
!    endif
!  enddo

! modified: 
  call hunt(xa,x,klo)
  khi=klo+1

  h=xa(khi)-xa(klo)
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h*h)/6d0
  
  if (present(dydx)) then
    dydx = -ya(klo)/h + ya(khi)/h + (-(2*a*a-1)*y2a(klo)+(2*b*b-1)*y2a(khi))*h/6d0
  endif

  END SUBROUTINE splint

end module utils
