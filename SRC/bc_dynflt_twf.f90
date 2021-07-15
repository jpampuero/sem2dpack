module bc_dynflt_twf

! BC_DYNFLT_TWF: Time weakening friction for dynamic faults
! Usually for Andrews nucleation procedure:
! Minimum rupture velocity V imposed through time-dependent weakening
! L is the (minimum) width of the imposed rupture front
! Only for t < T.
!

  implicit none
  private

  type twf_type
    private
    integer :: kind
    double precision :: X,Z,mus,mud,mu0,L,V,T,Dc
  end type twf_type

  public :: twf_type, twf_read, twf_mu

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_TWF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Time weakening friction for dynamic faults 
!          with prescribed rupture speed.
! SYNTAX : &BC_DYNFLT_TWF kind, MuS, MuD, Mu0, X, Z, V, L, T /
!
! ARG: kind     [int] [1] Type of time-weakening history:
!               1 = expansion at constant speed V up to time T
!               2 = expansion at decreasing speed then contraction
!               3 = expansion as a pulse with length=L
!                   as in Andrews and Ben-Zion (JGR 1997, eqs 2 and 3)
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
! ARG: Mu0      [dble] [0.6d0] Friction coefficient at the hypocenter at time=0
! ARG: X,Z      [dble] [0d0] Position of hypocenter
! ARG: V        [dble] [1d3] Rupture propagation speed (initial speed if kind=2)
! ARG: L        [dble] [1d0] Size of weakening zone
! ARG: T        [dble] [huge] Total duration
!
! NOTE   : Time-weakening is usually applied as an artificial nucleation procedure.
!          The maximum size of the nucleation region is 2*V*T if kind=1, V*T/2 if kind=2
!
! END INPUT BLOCK

  subroutine twf_read(tw,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(twf_type), intent(out) :: tw
  integer, intent(in) :: iin

  double precision :: mus,mud,mu0,X,Z,L,V,T,Dc
  integer :: kind

  NAMELIST / BC_DYNFLT_TWF / kind,mus,mud,mu0,X,Z,L,V,T,Dc

  kind = 1
  mus = 0.6d0
  mud = 0.5d0
  mu0 = 0.6d0
  X = 0d0
  Z = 0d0
  V = 1d3
  L = 1d0
  T = huge(1d0)
  Dc= huge(1d0)

  read(iin,BC_DYNFLT_TWF,END=100)
100 continue

  if (kind<1) call IO_abort('BC_SWFF_init: kind must be > 0')
  if (kind>3) call IO_abort('BC_SWFF_init: kind must be < 4')
  if (mus<0d0) call IO_abort('BC_SWFF_init: MuS must be positive in BC_DYNFLT_TWF input block')
  if (mud<0d0) call IO_abort('BC_SWFF_init: MuD must be positive in BC_DYNFLT_TWF input block')
  if (mu0<0d0) call IO_abort('BC_SWFF_init: Mu0 must be positive in BC_DYNFLT_TWF input block')
  if (L<=0d0) call IO_abort('BC_SWFF_init: L must be positive in BC_DYNFLT_TWF input block')
  if (V<=0d0) call IO_abort('BC_SWFF_init: V must be positive in BC_DYNFLT_TWF input block')
  if (T<=0d0) call IO_abort('BC_SWFF_init: T must be positive in BC_DYNFLT_TWF input block')
  if (Dc<=0d0) call IO_abort('BC_SWFF_init: Dc must be positive in BC_DYNFLT_TWF input block')

  tw%kind = kind
  tw%mus = mus
  tw%mud = mud
  tw%mu0 = mu0
  tw%X = X
  tw%Z = Z
  tw%L = L
  tw%T = T
  tw%V = V
  tw%Dc= Dc

  if (echo_input) write(iout,200) kind,mus,mud,mu0,X,Z,V,L,T,Dc

  return
  200 format(5x,'Friction law  . . . . . . . . . . . . . .  = time weakening', &
            /5x,'  Type of weakening history . . . . (kind) = ',I0, &
            /5x,'  Static friction coefficient . . . .(MuS) = ',EN13.3, &
            /5x,'  Dynamic friction coefficient  . . .(MuD) = ',EN13.3, &
            /5x,'  Initial friction coefficient  . . .(Mu0) = ',EN13.3, &
            /5x,'  Hypocenter position X . . . . . . . .(X) = ',EN13.3, &
            /5x,'  Hypocenter position Z . . . . . . . .(Z) = ',EN13.3, &
            /5x,'  Propagation speed . . . . . . . . . .(V) = ',EN13.3, &
            /5x,'  Size of weakening zone  . . . . . . .(L) = ',EN13.3, &
            /5x,'  Duration  . . . . . . . . . . . . . .(T) = ',EN13.3, &
            /5x,'  Turn-off slip  . . . . . . . . . . .(Dc) = ',EN13.3 )

  end subroutine twf_read


!=====================================================================
  function twf_mu(tw,coord,time,d) result(mu)

  type(twf_type), intent(in) :: tw
  double precision, intent(in) :: coord(:,:),time
  double precision :: mu(size(coord,2))
  double precision, dimension(:), intent(in) :: d

  integer :: i
  double precision :: r(size(coord,2)),t
  double precision, parameter :: VERY_LARGE_VALUE = huge(1d0)

 ! compute the position of the front (where mu=mus)
 ! the time shift sets mu=mu0 at the hypocenter
  if (tw%kind==1) then
    t = time + (tw%mus-tw%mu0)*tw%L/((tw%mus-tw%mud)*tw%V)
!    t = min(t,tw%T)  ! version 1 (old): total weakening persists after nucleation
    if (t> tw%T) t=0d0 ! version 2 (new): reset the time-weakening coefficient to its static value after nucleation is over
    r = tw%V*t
  elseif (tw%kind==2) then
    t = time+ 0.5d0*tw%T*( 1d0-sqrt( 1d0-4d0*(tw%mus-tw%mu0)*tw%L/((tw%mus-tw%mud)*tw%T*tw%V) ) )
    t = min(t,tw%T)
    r = tw%V*t*(1d0-t/tw%T) 
  endif
  
  if (tw%kind==1 .or. tw%kind==2) then
    ! relative position of fault node with respect to the front
    r = sqrt( (coord(1,:)-tw%X)*(coord(1,:)-tw%X) + (coord(2,:)-tw%Z)*(coord(2,:)-tw%Z) ) -r

    ! friction coefficient
    ! version A (old): mu = mus beyond the nucleation front
    ! mu = tw%mus + (tw%mus-tw%mud)/tw%L *r
    ! mu = max( mu, tw%mud )
    ! mu = min( mu, tw%mus )
   
    ! friction coefficient
    ! version B (new): mu keeps growing linearly beyond the nucleation front up to distance L, 
    !                  then jumps to a huge value
    do i=1,size(r)
      if (r(i)< -tw%L) then
        mu(i) = tw%mud
      elseif ( r(i) <= tw%L ) then  
     !NOTE: replace by "if r(i)<=0d0" to strongly enforce the position of the front
        mu(i) = tw%mus + (tw%mus-tw%mud)/tw%L *r(i)
      else
        mu(i) = VERY_LARGE_VALUE
      endif
    enddo

  else  ! New added cohesive nucleation
    t = time
    r = sqrt((coord(1,:)-tw%X)*(coord(1,:)-tw%X) + (coord(2,:)-tw%Z)*(coord(2,:)-tw%Z))
    do i=1,size(r)
      if (r(i) <= tw%V*tw%T .and. d(i) <= tw%Dc) then
        if (r(i) < tw%V*t - tw%L ) then
          mu(i) = tw%mud
        elseif (r(i) >= tw%V*t - tw%L .and. r(i) <= tw%V*t) then
          mu(i) = tw%mus + (tw%mus-tw%mud)/tw%L * (r(i)-tw%V*t)
        else
          mu(i) = VERY_LARGE_VALUE
        endif
      elseif (r(i) <= tw%V*tw%T .and. d(i) > tw%Dc) then 
          mu(i) = VERY_LARGE_VALUE
      endif
    enddo
  endif

  
  end function twf_mu

end module bc_dynflt_twf
