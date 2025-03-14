module bc_dynflt_swf

! BC_DYNFLT_SWF: slip weakening friction for dynamic faults

  use distribution_cd

  implicit none
  private

  type swf_input_type
    type(cd_type) :: dc, mus, mud, alpha,p
  end type swf_input_type

  type swf_type
    private
    integer :: kind
    double precision :: dt
    logical :: healing
    double precision, dimension(:), pointer :: dc=>null(), mus=>null(), mud=>null(), theta=>null(), p=>null()
    double precision, dimension(:), pointer :: alpha=>null()
    type(swf_input_type) :: input
  end type swf_type

 !slip velocity threshold for healing
 !WARNING: not very robust
  double precision, parameter :: V_HEALING = 1d-14

  public :: swf_type, swf_read, swf_init, swf_mu, &
            swf_update_state, swf_set_state

contains

!---------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_SWF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Slip-weakening friction
! SYNTAX : &BC_DYNFLT_SWF Dc | DcH, MuS | MuSH , MuD | MuDH, healing /
!          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
!          arguments with suffix H, if present, in the order listed above.
!
! ARG: kind     [int] [1] Type of slip weakening function:
!                       1 = linear
!                       2 = exponential
!                       3 = power-law
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
! ARG: p        [dble] [3d0] Exponent of power-law slip-weakening
! ARG: alpha    [dble] [0.0d0] Roughness drag coefficient (normalized by Tn)
! ARG: healing  [log] [F] Instantaneous healing upon slip arrest
!               Healing is currently valid only with the leapfrog time scheme
!
! END INPUT BLOCK

! Read parameters from input file
  subroutine swf_read(swf,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(swf_type), intent(out) :: swf
  integer, intent(in) :: iin

  double precision :: Dc,MuS,MuD,alpha,p
  character(20) :: DcH,MuSH,MuDH,alphaH,pH
  integer :: kind
  character(20) :: kind_txt
  logical :: healing

  NAMELIST / BC_DYNFLT_SWF / kind,Dc,MuS,MuD,alpha,p,DcH,MuSH,MuDH,alphaH,pH,healing

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  MuD = 0.5d0
  alpha = 0.0d0
  p=3d0
  DcH = ''
  MuSH = ''
  MuDH = ''
  alphaH = ''
  pH=''
  healing = .false.

  read(iin,BC_DYNFLT_SWF,END=300)
300 continue
  
  select case (kind)
    case(1); kind_txt = 'Linear'
    case(2); kind_txt = 'Exponential'
    case(3); kind_txt = 'Power-law'
    case default; call IO_abort('BC_DYNFLT_SWF: invalid kind')
  end select
  swf%kind = kind
  swf%healing = healing
  
  call DIST_CD_Read(swf%input%Dc,Dc,DcH,iin,DcH)
  call DIST_CD_Read(swf%input%MuS,MuS,MuSH,iin,MuSh)
  call DIST_CD_Read(swf%input%MuD,MuD,MuDH,iin,MuDH)
  call DIST_CD_Read(swf%input%alpha,alpha,alphaH,iin,alphaH)
  call DIST_CD_Read(swf%input%p,p,pH,iin,pH)
  if (echo_input) write(iout,400) kind_txt,DcH,MuSH,MuDH,alphaH,pH,healing

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = slip weakening', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Dynamic friction coefficient  . . .(MuD) = ',A,&
            /5x,'  Roughness drag coefficient  . . .(alpha) = ',A,&
            /5x,'  Power-law exponent  . . . . . . . . .(p) = ',A,&
            /5x,'  Instantaneous healing . . . . .(healing) = ',L1)

  end subroutine swf_read

!=====================================================================
! Initialize parameters
  subroutine swf_init(swf,coord,dt)

  type(swf_type), intent(inout) :: swf
  double precision, intent(in) :: coord(:,:)
  double precision, intent(in) :: dt

  call DIST_CD_Init(swf%input%dc,coord,swf%dc)
  call DIST_CD_Init(swf%input%mus,coord,swf%mus)
  call DIST_CD_Init(swf%input%mud,coord,swf%mud)
  call DIST_CD_Init(swf%input%p,coord,swf%p)
  call DIST_CD_Init(swf%input%alpha,coord,swf%alpha)

  allocate( swf%theta(size(coord,2)) )
  swf%theta = 0d0 

  swf%dt = dt

  end subroutine swf_init

!=====================================================================
! Friction coefficient
  function swf_mu(f) result(mu)

  type(swf_type), intent(in) :: f
  double precision :: mu(size(f%theta))

 !-- linear slip weakening:
  if (f%kind==1) then
    mu = f%mus -(f%mus-f%mud)*min(f%theta/f%dc,1d0)
  elseif (f%kind==2) then 
 !-- exponential slip weakening:
    mu = f%mud -(f%mud-f%mus)*exp(-f%theta/f%dc)
 !--  power-law slip weakening:
  elseif (f%kind==3) then  
    mu = f%mud + (f%mus-f%mud) / (1d0+f%theta/f%dc)**f%p
  endif

  mu = mu + f%alpha * f%theta

  end function swf_mu

!=====================================================================
  subroutine swf_update_state(d,v,f)

  double precision, dimension(:), intent(in) :: d,v
  type(swf_type), intent(inout) :: f

  integer :: k

  if (f%healing) then
   !WARNING: valid only for leapfrog. For other schemes we should use dD-bc%D
    f%theta = f%theta +abs(v)*f%dt
    do k=1,size(v)
    if (abs(v(k))<V_HEALING) then
      f%theta(k) = 0d0
    endif
    enddo

  else
    f%theta = abs(d)
  endif

  end subroutine swf_update_state


!=====================================================================
  subroutine swf_set_state(d,f)

  double precision, dimension(:), intent(in) :: d
  type(swf_type), intent(inout) :: f

  f%theta = abs(d)

  end subroutine swf_set_state


end module bc_dynflt_swf
