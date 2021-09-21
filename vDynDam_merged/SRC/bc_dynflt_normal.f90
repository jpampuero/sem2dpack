module bc_dynflt_normal

! BC_DYNFLT_NOR: normal stress response for dynamic faults 

  implicit none
  private

  type normal_type
    private
    integer :: kind
    double precision, dimension(:), pointer:: sigma
    double precision :: T,L,V,coef
  end type normal_type

  public :: normal_type, normal_read, normal_init, normal_update, normal_getSigma

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_NOR
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Normal stress response for dynamic faults.
! SYNTAX : &BC_DYNFLT_NOR kind, V, L, T /
!
! ARG: kind     [int] [1] Type of normal stress response:
!                       0 = shear strength is independent of normal stress
!                           (the cohesive strength is set as the product of
!                           friction coefficient and initial normal stress)
!                       1 = Coulomb 
!                       2 = Prakash-Clifton with regularizing time scale
!                       3 = Prakash-Clifton with regularizing length scale
! ARG: T        [dble] [1d0] Regularization time scale if kind=2
! ARG: V        [dble] [1d0] Characteristic velocity if kind=3
! ARG: L        [dble] [1d0] Regularization length scale if kind=3
!
! END INPUT BLOCK

  subroutine normal_read(n,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(normal_type), intent(out) :: n
  integer, intent(in) :: iin

  double precision :: L,V,T
  integer :: kind

  NAMELIST / BC_DYNFLT_NOR / kind,L,V,T

  kind = 1
  L = 1d0
  V = 1d0
  T = 1d0

  read(iin,BC_DYNFLT_NOR,END=100)
100 continue

  if (kind>3 .or. kind<0) call IO_abort('BC_SWFF_init: invalid kind in BC_DYNFLT_NOR input block')
  if (L<=0d0) call IO_abort('BC_SWFF_init: L must be positive in BC_DYNFLT_NOR input block')
  if (V<=0d0) call IO_abort('BC_SWFF_init: V must be positive in BC_DYNFLT_NOR input block')
  if (T<=0d0) call IO_abort('BC_SWFF_init: T must be positive in BC_DYNFLT_NOR input block')

  n%kind = kind
  n%L = L
  n%T = T
  n%V = V

  if (echo_input) then
    select case (n%kind)
    case (0) 
      write(iout,200) 'Cohesion'
    case (1) 
      write(iout,200) 'Coulomb'
    case (2) ! modified Prakash Clifton law #1
      write(iout,210) 'Prakash-Clifton with time-scale',T
    case (3) ! modified Prakash Clifton law #2
      write(iout,220) 'Prakash-Clifton with length-scale',L,V
    end select
  endif

  return
  200 format(5x,'Normal stress law . . . . . . . . . (kind) = ',A)
  210 format(5x,'Normal stress law . . . . . . . . . (kind) = ',A,&
            /5x,'  Time scale  . . . . . . . . . . . . .(T) = ',EN13.3)
  220 format(5x,'Normal stress law . . . . . . . . . (kind) = ',A,&
            /5x,'  Length scale  . . . . . . . . . . . .(L) = ',EN13.3,&
            /5x,'  Velocity scale  . . . . . . . . . . .(V) = ',EN13.3)

  end subroutine normal_read

!=====================================================================

  subroutine normal_init(n,dt,sigma_0)

  type(normal_type), intent(inout) :: n
  double precision, intent(in) :: dt
  double precision, intent(in) :: sigma_0(:)

  allocate(n%sigma(size(sigma_0)))

  select case (n%kind)
    case (0,1) ! cohesion or Coulomb friction
      continue
    case (2) ! modified Prakash Clifton law #1
      n%coef = exp(-dt/n%T)
    case (3) ! modified Prakash Clifton law #2
      n%coef = dt/n%L
  end select

  n%sigma = sigma_0

  end subroutine normal_init

!=====================================================================
  subroutine normal_update(n,Tn,V)

  type(normal_type), intent(inout) :: n
  double precision, intent(in) :: Tn(:),V(:)

  select case (n%kind)
    case (0) ! cohesion
      continue
    case (1) ! Coulomb friction
      n%sigma = Tn
    case (2) ! modified Prakash Clifton law #1
      n%sigma = Tn + n%coef *(n%sigma - Tn)
                       ! coef = exp(-dt/sigma_T)
    case (3) ! modified Prakash Clifton law #2
      n%sigma = Tn + exp(-(abs(V)+n%V)*n%coef) *(n%sigma - Tn)
                                               ! coef = dt/L
  end select

  end subroutine normal_update

!=====================================================================
  function normal_getSigma(n) result(sigma)

  type(normal_type), intent(in) :: n
  double precision, pointer :: sigma(:)

  sigma => n%sigma

  end function normal_getSigma

end module bc_dynflt_normal
