
! Low-pass Butterworth recursive filter, order 1 to 4
! Adapted from J-P Moreau, see http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tfilters_f90.txt
module butterworth_filter

  implicit none

! BUTTERWORTH MODULE HANDLE:

  type butter_handle_type
    integer :: NSections=0
    double precision, pointer :: C(:,:)=>null(), D(:,:)=>null()
  end type

!        C........: Table[1..5,1..NSections] of filter coefficients      
!                   calculated by BUTTERWORTH subroutine      
!        NSections: Number of required 2nd order sections (integer)      
!                   = n/2     for n even (n=order of filter)             
!                   = (n+1)/2 for n odd                                  
!                   calculated by BUTTERWORTH subroutine      
!        D........: Table[1..2,1..NSections] of coefficients defining    
!                   the filter memory, initialized by INIT subroutine    

  type butter_type
    type(butter_handle_type) :: H
  end type butter_type

contains


!=======================================================================
subroutine BUTTER_read(butter,iin)

  use stdio, only : IO_abort

  !type (butter_type), intent(out) :: butter
  type (butter_type) :: butter
  integer, intent(in) :: iin
  
  call IO_abort('BUTTER_read: not implemented')
  return

end subroutine BUTTER_read


!=======================================================================
double precision function butter(b,t)

  use stdio, only : IO_abort

  type (butter_type), intent(in) :: b
  double precision, intent(in) :: t

  call IO_abort('BUTTER: not implemented')
  butter = 0.d0

end function butter

!=======================================================================
!! A simple driver follows. 
!! For more intensive use you must decompose the calls in your program.
subroutine butter_single(signal,filtered,Fc,ndata,iorder,dt)

  integer, intent(in) :: ndata,iorder
  double precision, intent(in) :: signal(ndata),dt,Fc
  double precision, intent(out) :: filtered(ndata)

  type (butter_handle_type) :: H
  integer :: i

  !1. Calculate the filter coefficients
  !   for given cut-off frequency (Fc)
  !           , sampling period (dt) 
  !         and order (iorder)
  call Butterworth(Fc,dt,iorder,H)

  !2. Initialize filter memory
  !   Must be called for each data set before filtering
  call Init(H)

  !3. Recursively call Butterworth filter
  do i=1,ndata
    call Filter(filtered(i),signal(i),H)
  end do

end subroutine butter_single

!=======================================================================
!***********************************************************************
!*          Filtering a signal F(t) by Butterworth method              *
!*             (removing frequencies greater then Fc)                  *
!* ------------------------------------------------------------------- *
!* Calling mode:   Filter(Xs,Xd,H);                                    *
!* ------------------------------------------------------------------- *
!* INPUTS:                                                             *
!* -------                                                             *
!*      Xd.......: current value of input signal (double)              *
!* ------------------------------------------------------------------- *
!* OUTPUTS:                                                            *
!* -------                                                             *
!*      Xs.........: current value of filtered signal (double)         *
!* ------------------------------------------------------------------- *
!* INOUT:                                                              *
!* ------                                                              *
!*      H..........: Module handle                                     *
!* ------------------------------------------------------------------- *
!* Reference                                                           *
!* ---------                                                           *
!*  "Lawrence R.Rabiner et Bernard Gold                                *
!*   Theory and application of digital processing.                     *
!*   Prentice Hall Inc., EnglewoodclIFfs,NEW JERSEY,1975."             *
!*                                                                     *
!*             Module version by J-P Ampuero, ampuero@ipgp.jussieu.fr  *
!*                                   F90 version by J-P Moreau, Paris  *
!*               from Fortran version by J-P Dumont / Dang Trong Tuan  *
!***********************************************************************
Subroutine Filter(Xs,Xd,H) 

  type (butter_handle_type), intent(inout) :: H
  double precision, intent(in) ::  Xd
  double precision, intent(out) ::  Xs

  double precision :: x,xerr,y
  integer :: ii

    x=Xd
    do ii=1, H%NSections
      xerr=x+H%C(1,ii)*H%D(1,ii)+H%C(2,ii)*H%D(2,ii)
      y=H%C(5,ii)*(xerr +H%C(3,ii)*H%D(1,ii)+H%C(4,ii)*H%D(2,ii))
      H%D(2,ii)=H%D(1,ii)
      H%D(1,ii)=xerr
      x=y
    end do
    Xs=x
    
End subroutine

!**************************************************************************
!*                       INIT FILTER PROCEDURE                            *
!* ---------------------------------------------------------------------- *
!* The filter response is initialized to stationnary value for a constant *
!* input signal value.                                                    *
!*                                                                        *
!* Calling mode:   INIT(H,Xdc);                                           *
!* ---------------------------------------------------------------------- *
!* OPTIONAL INPUT:                                                        *
!* ---------------                                                        *
!*        Xdc......: constant input value (double)                        *
!* ---------------------------------------------------------------------- *
!* INOUT:                                                                 *
!* ------                                                                 *
!*        H........: Module handle                                        *
!**************************************************************************
Subroutine Init(H,Xdc)

  type (butter_handle_type), intent(inout) :: H
  double precision, optional, intent(in) :: Xdc

  double precision :: dc,Csum 
  integer :: j,ii

  if (present(Xdc)) then
    dc=Xdc
    DO j=1, H%NSections
      H%D(2,j)=dc/(1.d0-H%C(1,j)-H%C(2,j))
      H%D(1,j)=H%D(2,j)
      Csum=0.d0
      do ii=1, 4
        Csum=Csum + H%C(ii,j)
      end do
      dc=H%C(5,j)*(dc+H%D(2,j)*Csum)
    END DO

  else
    H%D = 0.d0
  endif
    
END subroutine

!**********************************************************************
!*          Calculates the Butterworth filter coefficients            *
!* ------------------------------------------------------------------ *
!*  Calling mode:   Butterworth(Fc,Ts,n,C,NSections,Tg);              *
!* ------------------------------------------------------------------ *
!*  INPUTS:                                                           *
!*  ------                                                            *
!*         Fc.......: Cut off frequency                               *
!*         dt.......: Sampling time of input signal                   *
!*         N........: Order of filter (1 to 4)                        *
!* ------------------------------------------------------------------ *
!*  OUTPUTS:                                                          *
!*  -------                                                           *
!*         H........: Module handle                                   *
!*         Tgd......: Group delay in seconds (optional)               *
!**********************************************************************
Subroutine Butterworth(Fc,dt,N,H,Tgd)

  double precision, intent(in) :: Fc,dt
  integer, intent(in) :: N
  type(butter_handle_type) :: H
  double precision, optional, intent(out) :: Tgd

  double precision :: arg,m,Omega,OmegaSq,Pi,temp,Zero,ONE,TWO,HALF &
           ,Rep,W0,W1,Tg
  integer :: i

  Zero = 0.d0
  ONE  = 1.d0
  TWO  = 2.d0
  HALF = 0.5d0

  Pi = 3.1415926535d0
  arg=Pi*dt*Fc
  If (abs(arg) > 2.d0*Pi) THEN
    m=INT(arg/2.d0/Pi)
    arg=arg-(m*2.d0*Pi)
  END IF
  Omega= tan(arg)
  OmegaSq=Omega*Omega
  select case(N)
    case(1)
      H%NSections = 1
      temp = zero
    case(2)
      H%NSections = 1
      temp = half
    case(3)
      H%NSections = 2
      temp = zero
    case(4)
      H%NSections = 2
      temp = half
    case default
      stop 'Butterworth: oder must be in [1,4]'
  end select

  if (associated(H%C)) deallocate(H%C) 
  allocate( H%C(5,H%NSections) )

  Tg=Zero
  If (N>1) THEN
    DO i=1, N/2
      Rep=Omega*COS(Pi*(i-temp)/N)
      Tg=Tg+dt*Rep/OmegaSq
      W0=TWO*Rep
      W1=ONE + W0+OmegaSq
      H%C(1,i)=-TWO*(OmegaSq-ONE)/W1
      H%C(2,i)=-(ONE-W0+OmegaSq)/W1
      H%C(3,i)=TWO
      H%C(4,i)=ONE
      H%C(5,i)=OmegaSq/W1
    end do
  end if
  If (temp.eq.Zero) THEN
    H%C(1,H%Nsections)=(ONE-Omega)/(ONE+Omega)
    H%C(2,H%NSections)= Zero
    H%C(3,H%NSections)= ONE
    H%C(4,H%NSections)= Zero
    H%C(5,H%NSections)= Omega/(ONE+Omega)
    Tg= Tg+dt/(TWO*Omega)
  END IF
  if (present(Tgd)) Tgd = Tg

  if (associated(H%D)) deallocate(H%D)
  allocate( H%D(2,H%NSections) )

END subroutine !of Butterworth()

end module butterworth_filter
