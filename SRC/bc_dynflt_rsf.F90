module bc_dynflt_rsf

! BC_DYNFLT_RSF: rate and state dependent friction for dynamic faults
! DEVEL: module in progress ... add theta, write friction solver

  use distribution_cd
  use stdio, only: IO_abort

  implicit none
  private

  type rsf_input_type
    ! plate rate is added
    type(cd_type) :: dc, mus, a, b, Vstar, theta, vplate 
  end type rsf_input_type

  type rsf_type
    private
    integer :: kind
    double precision, dimension(:), pointer :: dc=>null(),& 
                     mus=>null(), a=>null(), b=>null(), &
                     Vstar=>null(), theta=>null(), Tc=>null(),&
                     coeft=>null(), vplate=>null(), theta_pre=>null()
    double precision :: dt
    double precision :: dtScale
    double precision :: vmaxPZ ! maximum slip rate to control the process zone size
    integer :: iter, NRMaxIter, StepLock, minStep
    double precision :: vmaxD2S, vmaxS2D, NRTol, vEQ, tEqPrev, minGap
    double precision :: vmaxD2S_VS
    type(rsf_input_type) :: input
  end type rsf_type

  public :: rsf_type, rsf_read, rsf_init, rsf_mu, & 
            rsf_solver, rsf_qs_solver, rsf_timestep, & 
            rsf_vplate, rsf_get_theta, rsf_get_a, &
            rsf_get_b, rsf_get_NRTol, rsf_reset,rsf_vmaxD2S,& 
            rsf_vmaxS2D,rsf_update_theta_pre 

contains

!---------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_RSF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Velocity and state dependent friction
! SYNTAX : &BC_DYNFLT_RSF kind, Dc | DcH, Mus | MusH , 
!                         a | aH, b | bH, Vstar | VstarH 
!                         vmaxS2D, vmaxD2S /
!          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
!          arguments with suffix H, if present, in the order listed above.
!
! ARG: kind     [int] [1] Type of rate-and-state friction law:
!                       1 = strong velocity-weakening at high speed
!                           as in Ampuero and Ben-Zion (2008)
!                       2 = logarithmic rate-and-state with aging state law
!                       3 = logarithmic rate-and-state with slip state law
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: a        [dble] [0.01d0] Direct effect coefficient
! ARG: b        [dble] [0.02d0] Evolution effect coefficient
! ARG: Vstar    [dble] [1d0] Characteristic or reference slip velocity
! ARG: theta    [dble] [1d0] State variable
! ARG: vplate   [dble] [1d0] plate rate used for 'back-slip' loading
! ARG: vmaxS2D  [dble] [5d-3] max slip velocity for switch from static to dynamic 
! ARG: vmaxD2S  [dble] [2d-3] max slip velocity for switch from dynamic to static
! ARG: vEQ      [dble] [10d-3] threshold slip velocity for earthquake event 
! ARG: minGap   [dble] [10 s] minimum time gap between two earthquake events
! ARG: NRTol    [dble] [1.0d-4] relative tolerance for NR solver 
! ARG: NRMaxIter [Int] [200] maximum interation number for NR solver 
! ARG: minStep [Int]  [50] minimum step before another switch 
!
! END INPUT BLOCK

! Read parameters from input file
  subroutine rsf_read(rsf,iin)

  use echo, only : echo_input,iout

  type(rsf_type), intent(out) :: rsf
  integer, intent(in) :: iin

  double precision :: Dc,MuS,a,b,Vstar,theta,vplate,vmaxS2D,vmaxD2S, vEQ
  double precision :: NRTol,minGap,dtScale,vmaxPZ, vmaxD2S_VS 
  integer :: NRMaxIter, minStep 
  character(20) :: DcH,MuSH,aH,bH,VstarH,thetaH,vplateH
  integer :: kind
  character(25) :: kind_txt

  NAMELIST / BC_DYNFLT_RSF / kind,Dc,MuS,a,b,Vstar,theta,&
           DcH,MuSH,aH,bH,VstarH,thetaH, vplate, vplateH,&
           vmaxS2D, vmaxD2S, NRMaxIter, NRTol, vEQ, minGap, &
           dtScale, vmaxPZ, minStep,vmaxD2S_VS

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  a = 0.01d0
  b = 0.02d0
  Vstar = 1d0
  theta = 0d0
  vplate = 0d0
  DcH = ''
  MuSH = ''
  aH = ''
  bH = ''
  VstarH = ''
  thetaH = ''
  vplateH = ''
  vmaxS2D = 5d-3 ! 5 mm/s
  vmaxD2S = 2d-3 ! 2 mm/s
  vmaxD2S_VS = 1d-2 ! 10 mm/s for velocity strengthening region 
  vEQ     = 10d-3 ! 10 mm/s
  minGap  = 10 ! 10 s
  NRMaxIter = 2000 
  minStep  = 50
  NRTol     = 1.0d-4
  dtScale = 1d0
  vmaxPZ  = 1d10 ! by default a very large value, no control on process zone size

  read(iin,BC_DYNFLT_RSF,END=300)
300 continue
  
  select case (kind)
    case(1); kind_txt = 'Strong velocity-weakening'
    case(2); kind_txt = 'Classical with aging law'
    case(3); kind_txt = 'Classical with slip law'
    case default; call IO_abort('BC_DYNFLT_RSF: invalid kind')
  end select
  rsf%kind = kind
  
  call DIST_CD_Read(rsf%input%Dc,Dc,DcH,iin,DcH)
  call DIST_CD_Read(rsf%input%MuS,MuS,MuSH,iin,MuSH)
  call DIST_CD_Read(rsf%input%a,a,aH,iin,aH)
  call DIST_CD_Read(rsf%input%b,b,bH,iin,bH)
  call DIST_CD_Read(rsf%input%Vstar,Vstar,VstarH,iin,VstarH)
  call DIST_CD_Read(rsf%input%theta,theta,thetaH,iin,thetaH)
  call DIST_CD_Read(rsf%input%vplate,vplate,vplateH,iin,vplateH)

  rsf%vmaxS2D = vmaxS2D
  rsf%vmaxD2S = vmaxD2S
  rsf%vmaxD2S_VS = vmaxD2S_VS
  rsf%vEQ     = vEQ
  rsf%minGap  = minGap
  rsf%NRTol   = NRTol
  rsf%NRMaxIter = NRMaxIter
  rsf%minStep = minStep
  rsf%StepLock = 0
  rsf%teqPrev = -huge(1d0) 
  rsf%dtScale = dtScale
  rsf%vmaxPZ  = vmaxPZ

  if (echo_input) write(iout,400) kind_txt,DcH, &
               MuSH,aH,bH,VstarH,thetaH,vplateH,& 
               vmaxS2D,vmaxD2S, vEQ, NRMaxIter, minStep, NRTol

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = Rate and State Dependent', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Direct effect coefficient . . . . . .(a) = ',A,&
            /5x,'  Evolution effect coefficient  . . . .(b) = ',A,&
            /5x,'  Velocity scale  . . . . . . . . .(Vstar) = ',A,&
            /5x,'  State variable  . . . . . . . . .(theta) = ',A,&
            /5x,'  Plate rate  . . . . . . . . . . (vplate) = ',A,&
            /5x,'  V threshold (Static to Dyanmic) (vmaxS2D) = ',EN12.3,&
            /5x,'  V threshold (Dyanmic to Static) (vmaxD2S) = ',EN12.3,&
            /5x,'  V threshold (earthquake event) . . .(vEQ) = ',EN12.3,&
            /5x,'  NR Max iteration . . . . . . .(NRMaxIter) = ',I4,&
            /5x,'  Minimum locked steps . . . . . .(minStep) = ',I4,&
            /5x,'  NR tolerance . . . . . . . . . . . (NRTol)= ',EN12.3)

  end subroutine rsf_read

!=====================================================================
! Initialize parameters
! TO DO: Change the initialization method for rsf
! Initialising the state variable by specifing initial slip velocity
! and initial stress (or Mu)
!

  subroutine rsf_init(rsf, coord, dt, mu, v)

  type(rsf_type), intent(inout) :: rsf
  double precision, intent(in) :: coord(:,:),dt
  double precision, intent(in) :: mu(:), v(:)
  integer :: n
  
  call DIST_CD_Init(rsf%input%dc,coord,rsf%dc)
  call DIST_CD_Init(rsf%input%mus,coord,rsf%mus)
  call DIST_CD_Init(rsf%input%a,coord,rsf%a)
  call DIST_CD_Init(rsf%input%b,coord,rsf%b)
  call DIST_CD_Init(rsf%input%Vstar,coord,rsf%Vstar)
  call DIST_CD_Init(rsf%input%theta,coord,rsf%theta)
  call DIST_CD_Init(rsf%input%vplate,coord,rsf%vplate)

  ! reinitialize state variable theta using initial mu and v
  call rsf_init_theta(rsf, mu, v)
  n = size(coord,2)

  allocate(rsf%Tc(n))
  allocate(rsf%coeft(n))
  allocate(rsf%theta_pre(n))
  rsf%theta_pre = rsf%theta
                                              
  rsf%Tc = rsf%dc / rsf%Vstar
  ! if adaptive stepping is activated, 
  ! dt and coeft needs to be updated
  ! however coeft and Tc seems only been used
  ! in the kind = 1 case.
  rsf%coeft = exp(-dt/rsf%Tc)
  rsf%dt = dt
  rsf%iter = 0
  end subroutine rsf_init

  subroutine rsf_init_theta(f, mu, v)
      type(rsf_type), intent(inout) :: f
      double precision, intent(in) :: mu(:), v(:)
      double precision, dimension(size(f%theta)) :: tmp
      double precision:: tmp_i
      integer :: i
      
      ! initialize
!      f%theta = exp((mu - f%mus - f%a * log(v/f%Vstar))/f%b) & 
!               * f%Dc / f%Vstar
      
      tmp = f%a/f%b * log(2*f%Vstar/v * sinh(mu/f%a))
      f%theta = f%Dc/f%Vstar * exp(tmp - f%mus/f%b)
      tmp_i   = 0d0

      do i = 1, size(f%theta)
          if (f%theta(i)>huge(0d0) .or. tmp(i)>huge(0d0)) then
              ! expand the sinh and cancel the factor a
              tmp_i = f%a(i)/f%b(i) * log(f%Vstar(i)/v(i))
              tmp_i = mu(i)/f%b(i) + tmp_i
              f%theta(i) = f%Dc(i)/f%Vstar(i) * exp(tmp_i- f%mus(i)/f%b(i))
          end if
      end do

  end subroutine rsf_init_theta

! reset the state variable and update to new dt
  subroutine rsf_reset(f, dt)
    type(rsf_type), intent(inout) :: f
    double precision :: dt
   
    f%theta = f%theta_pre
    f%dt    = dt
    f%iter  = 0
    if (f%kind==1) f%coeft = exp(-f%dt/f%Tc)
  
  end subroutine rsf_reset

  function rsf_get_theta(f) result(theta)
      type(rsf_type), intent(in) :: f
      double precision, dimension(size(f%theta)):: theta
      theta = f%theta
  end function !rsf_get_theta 
  
  function rsf_get_a(f) result(a)
      type(rsf_type), intent(in) :: f
      double precision, dimension(size(f%a)):: a
      a = f%a
  end function !rsf_get_a 
  
  function rsf_get_b(f) result(b)
      type(rsf_type), intent(in) :: f
      double precision, dimension(size(f%b)):: b
      b = f%b
  end function !rsf_get_b 
  
  function rsf_get_NRTol(f) result(NRTol)
      type(rsf_type), intent(in) :: f
      double precision:: NRTol
      NRTol = f%NRTol
  end function !rsf_get_NRTol 

!=====================================================================
! Friction coefficient
! v : slip velocity

  function rsf_mu(v,f) result(mu)

      double precision, dimension(:), intent(in) :: v
      type(rsf_type), intent(in) :: f
      double precision :: mu(size(v)), tmp(size(v))
      integer :: i

      select case(f%kind)
        case(1) 
          mu = f%mus +f%a*abs(v)/(abs(v)+f%Vstar) - f%b*f%theta/(f%theta+f%Dc) 
        case(2,3) 
          !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
          ! mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
          !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
          ! expand asinh(x) = log(x+sqrt(x.^2+1))
          ! mu = a*log(x+sqrt(x.^2+1)), 
          ! where x = exp((f0 + b*log(v0*theta/dc))/a + log(v/(2*v0)))

          tmp  = (f%mus+f%b*log(f%Vstar*f%theta/f%Dc))/f%a + log(abs(v)/(2d0*f%Vstar)) 

          ! if exp(tmp) is too large, it causes trouble! expand sinh
          mu = f%a*asinh(exp(tmp))
          do i = 1, size(v)
              if (mu(i)>huge(0d0) .or. exp(tmp(i))>1d50) then
                  ! asinh(exp(tmp)) ~ log(2*exp(tmp)) = tmp + log(2)
                  mu(i) = f%a(i)*(tmp(i) + log(2d0))
              end if
          end do
      end select

  end function rsf_mu
  
  function rsf_vmaxD2S(f) result(vmaxd2s)
      double precision :: vmaxd2s
      type(rsf_type), intent(in) :: f
      vmaxd2s = f%vmaxd2s
  end function
  
  function rsf_vmaxS2D(f) result(vmaxs2d)
      double precision :: vmaxs2d
      type(rsf_type), intent(in) :: f
      vmaxs2d = f%vmaxs2d
  end function

!---------------------------------------------------------------------
! friction coefficient without the direct effect 
! (i.e. without the term that depends explicitly on slip velocity V)
  function rsf_mu_no_direct(f) result(mu)

  type(rsf_type), intent(in) :: f
  double precision :: mu(size(f%mus))

  select case(f%kind)
    case(1)
      mu = f%mus - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3) 
      !  Kaneko et al (2008) Eq. 12
      mu = f%mus + f%b*log(f%theta*f%Vstar/f%Dc) 
  end select

  end function rsf_mu_no_direct

!---------------------------------------------------------------------
!  Slip velocity at theta, tau (only used in quasi-static simulation)
!

  function rsf_v(f,tau,sigma,theta) result(v)

  type(rsf_type), intent(in) :: f
  double precision, dimension(:), intent(in) :: tau, sigma, theta
  double precision, dimension(size(theta)) :: tmp, v
  double precision:: tmp1, tmp2
  integer::i

  !  Kaneko et al (2011) Eq. 14
!  tmp = f%mus +f%b*log(f%Vstar*theta/f%Dc)
!  tmp = 2d0*f%Vstar*exp(-tmp/f%a)
!  v = sinh(abs(tau)/(-sigma*f%a))*tmp
!  v = sign(v, tau)

  tmp = -(f%mus +f%b*log( f%Vstar*theta/f%Dc))/f%a
  v = f%Vstar * (exp(tau/(-sigma*f%a)+tmp) - exp(tau/(sigma*f%a) + tmp))
  v = sign(v, tau)

  do i = 1, size(theta)
      if (abs(v(i))>1d0) then
          write(*, *) "error in rsf_v, abs(v)>1d0, i, v(i) = ", i, v(i)
          write(*, *) "tau(i), theta(i), a(i), b(i), dc(i) = ", tau(i), theta(i), &
                      f%a(i), f%b(i), f%dc(i)
          tmp1 = f%mus(i) +f%b(i)*log(f%Vstar(i)*theta(i)/f%Dc(i)) 
          tmp2 = 2d0*f%Vstar(i)*exp(-tmp1/f%a(i))
          write(*, *) "tmp1=", tmp1 
          write(*, *) "tmp2=", tmp2
          write(*, *) "sinh(abs(tau)/(sigma*a))", sinh(abs(tau(i))/(-sigma(i)*f%a(i)))
          write(*, *) "v = sinh(abs(tau)/(sigma*a))*tmp2 = ", sinh(abs(tau(i))/(-sigma(i)*f%a(i)))*tmp2
      end if
  end do

  end function rsf_v

!=====================================================================
!       --------------- RSF_SOLVER subroutine -----------------
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!      3. update theta again, now using V=(Vnew+Vold)/2
!      4. solve again for Vnew
!
! Note: most often sigma is negative (compressive) 
!
! add control on slip velocity
!
!
  subroutine rsf_solver(v, tau_stick, sigma, f, Z, time, isSetmaxV)
  use time_evol, only: timescheme_type
  type(timescheme_type):: time
  double precision, dimension(:), intent(inout) :: v
  logical, dimension(size(v)), intent(out) :: isSetmaxV
  double precision, dimension(:), intent(in) :: sigma,Z,tau_stick
  type(rsf_type), intent(inout) :: f

  double precision, dimension(size(v)) :: v_new,theta_new
  integer :: it, iter,  stat_i

  isSetmaxV = .false.

  time%solver_converge_stat = 1

  ! First pass: 
  theta_new = rsf_update_theta(f%theta,v,f)
  
  call rsf_update_V(tau_stick, sigma, f, theta_new, Z, v_new, time%nr_iters(1), stat_i)
  if (stat_i<0) then
      write(*, *) "rsf_solver, pass 1 diverge!"
      time%solver_converge_stat = stat_i
  end if

  ! limit the slip velocity to v_new
  do it = 1, size(v_new)
      if (abs(v_new(it))>f%vmaxPZ) then
          v_new(it) = sign(f%vmaxPZ, v_new(it))
      end if
  end do

  theta_new = rsf_update_theta(f%theta,0.5d0*(v+v_new), f)
  call rsf_update_V(tau_stick, sigma, f, theta_new, Z, v_new, time%nr_iters(2), stat_i)
  
  ! limit the slip velocity to v_new
  do it = 1, size(v_new)
      if (abs(v_new(it))>f%vmaxPZ) then
          v_new(it) = sign(f%vmaxPZ, v_new(it))
          isSetmaxV(it) = .true.
      end if
  end do

  if (stat_i<0) then
      write(*, *) "rsf_solver, pass 2 diverge!"
      time%solver_converge_stat = stat_i
  end if
  
  ! store new velocity and state variable estimate in friction law 
  f%theta = theta_new 
  v = v_new

  end subroutine rsf_solver

!=====================================================================
!       --------------- RSF_QS_SOLVER subroutine -----------------
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!
! Note: most often sigma is negative (compressive) 
!
  subroutine rsf_qs_solver(v,tau,sigma,f, time)
  use time_evol, only: timescheme_type
  type(timescheme_type):: time
  double precision, dimension(:), intent(in) :: sigma,tau
  double precision, dimension(:), intent(inout) :: v
  double precision, dimension(size(v)) :: theta, vnew
  type(rsf_type), intent(inout) :: f
  integer::i
 
  time%solver_converge_stat = 1

  f%theta_pre = f%theta
  theta = rsf_update_theta(f%theta_pre,v,f)
  vnew  = rsf_v(f, tau, sigma, theta)

  do i = 1, size(v)
     if (abs(v(i))>1d0) then
         ! static solver produces unphysical slip solution
         write(*, *) "rsf_qs_solver, unphysical v(it)>1, it, v(it)", i, vnew(i)
         write(*, *) "theta_pre(i), v_in(i) = ", f%theta_pre(i), v(i)
         exit
     end if
  end do

  v  = vnew

  ! increment the counter until 2 (passes)
  f%iter = f%iter + 1

  if (f%iter==2) then
      ! save theta after 2 passes
      f%theta = theta
      ! reset the counter to 0
      f%iter  = 0

      ! second pass, check convergence
      ! solution is higher than 1d-1
      do i = 1, size(v)
         if (abs(v(i))>1d0) then
             ! static solver produces unphysical slip solution
             write(*, *) "rsf_qs_solver, unphysical static slip rate, it, v(it)", i, v(i)
             time%solver_converge_stat = -3
             ! if solver does not converge then do not update theta
             f%theta = f%theta_pre
             exit
         end if
      end do
  endif
  if (time%solver_converge_stat<-3) call IO_Abort("rsf_qs_solver diverge! stop!")

  end subroutine rsf_qs_solver

  subroutine rsf_update_theta_pre(f)
      type(rsf_type), intent(inout) :: f
      f%theta_pre = f%theta
  end subroutine 

!---------------------------------------------------------------------
! Update state variable (theta) assuming slip velocity (v) is known

  function rsf_update_theta(theta,v,f) result(theta_new)

  double precision, dimension(:), intent(in) :: v,theta
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(v)) :: theta_new
!  integer::i

  select case(f%kind)
    case(1) 
     ! exact integration assuming constant V over the timestep
     ! Tc = Dc/Vstar, precomputed
     ! coeft = exp(-dt/Tc), precomputed
      theta_new = theta*f%coeft +f%Tc*abs(v)*(1d0-f%coeft)
    case(2) 
     ! Kaneko et al (2008) eq 19 - "Aging Law"
     where (abs(v)*f%dt/f%Dc < 1d-8)
         !using results from taylor series expansion when abs(v)*f%dt/f%Dc too small
         theta_new  = theta*exp(-abs(v)*f%dt/f%Dc) + f%dt*(1-0.5d0*abs(v)*f%dt/f%Dc)
     elsewhere
         theta_new = theta*exp(-abs(v)*f%dt/f%Dc)
         theta_new = theta_new + exp(log(f%Dc/abs(v)) + log(1-exp(-abs(v)*f%dt/f%Dc))) 
     end where

    case(3) 
     ! Kaneko et al (2008) eq 20 - "Slip Law"
     ! theta_new = Dc/v *(theta*v/Dc)**exp(-v*dt/Dc)
      theta_new = f%Dc/abs(v)
      theta_new = theta_new *(theta/theta_new)**exp(-f%dt/theta_new)
  end select
  
!  do i = 1, size(theta)
!      if (theta_new(i)<1d-20 .or. theta(i)<1d-20) then
!          write(*, *) "rsf_update_theta, it", i
!          write(*, *) "theta_new<1d-20, theta, theta_new, v", theta(i), theta_new(i),v(i)
!          write(*, *) "theta_new<1d-20, log10(theta, theta_new, v)", &
!          log10(theta(i)), log10(theta_new(i)),log10(v(i))
!      end if
!  end do

  end function rsf_update_theta

!---------------------------------------------------------------------
! Update slip velocity assuming theta is known
! The constraints are
!   (abs(tau)-strength)*v = 0
!   abs(tau)-strength <= 0
!   sign(tau) = sign(v)
! where
!   strength = -sigma*( mu(theta) +a*v/(1+v) )
!   tau = tau_stick-Z*v
!
! Inherited from the SBIEM code BIMAT-PCSI
! WARNING: the SBIEM code assumed v>0
!          We should allow here for any sign of v
!          Exploit the fact that sign(tau)=sign(tau_stick) (because mu>0)

  subroutine rsf_update_V(tau_stick,sigma,f,theta, Z, v, iter, stat) 
   
  double precision, dimension(:), intent(in) :: tau_stick,sigma,theta,Z
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(tau_stick)) :: v
  double precision :: tmp(size(tau_stick)), tolerance, estimateLow, estimateHigh
  double precision :: mintau(size(tau_stick)), maxtau(size(tau_stick))
  integer :: it, iter, iter_i, stat_i, stat!,i

  stat   = 1 ! converge normally
  maxtau = tau_stick
  iter_i = 0

  ! compute the lower bound initial guess
  tmp    = 2d0*f%Vstar/(-f%a*sigma)*exp(-(f%mus+f%b*log(f%Vstar*f%theta/f%dc))/f%a)
  mintau = tau_stick/(1d0 + Z*tmp)
  mintau = 0d0

  select case(f%kind)
    case(1) 
      v = (tau_stick +sigma*rsf_mu_no_direct(f))/Z ! if v<0 will stop
      tmp = v -f%Vstar +sigma*f%a/Z
      v = 0.5d0*( tmp +sqrt(tmp*tmp +4d0*v*f%Vstar) )
      v = max(0d0,v)  ! arrest if v<0 
 
    case(2,3) 
     ! "Aging Law and Slip Law"
     ! Find each element's velocity:
     do it=1,size(tau_stick)
       !DEVEL: What are the accepted tolerances and bounds? User-input? 
       tolerance=f%NRTol*f%a(it)*sigma(it) ! As used by Kaneko in MATLAB code
       estimateLow  = min(mintau(it), maxtau(it))
       estimateHigh = max(mintau(it), maxtau(it))
       stat_i       = 1
        
!       if (abs(estimateLow-estimateHigh) < abs(tolerance)) then
!           ! initial guess extremely close to tau_stick, fault is almost locked
!           v(it) = (tau_stick(it) - (estimateLow + estimateHigh)/2d0)/Z(it) 
!           v(it) = sign(v(it), tau_stick(it))
!       else
           !the initial guess assume zero increment of slip velocity
           
           ! try nr_solver first
           call nr_solver(nr_fric_func_tau,estimateLow,estimateHigh,&
               tolerance,f,it,theta(it),tau_stick(it),sigma(it),Z(it), v(it), iter_i, stat_i)

           if (stat_i <0) then
               call bs_solver(nr_fric_func_tau,estimateLow,estimateHigh,&
                   tolerance,f,it,theta(it),tau_stick(it),sigma(it),Z(it), v(it), iter_i, stat_i)
           end if
           iter = max(iter, iter_i)

!       end if

       if (stat_i<0) then
           stat = stat_i
           write(*, *) "estimatehigh, estimatelow", estimateHigh, estimateLow
           select case (stat_i)
               case (-1)
               write(*, *) "nr_solver diverge with maximum iteration number reached, it = ", it
               case (-2)
               write(*, *) "nr_solver diverge with NAN results, it = ", it
               case (-3)
               write(*, *) "nr_solver diverge, bs_solver reach maximum iteration number"
           end select
       end if
     enddo     
  end select
  
  if (stat<0) call IO_Abort("nr_solver diverge!")

  end subroutine rsf_update_V

! use bisection to solve the rsf
! when a is too small, it is not very stable to use newton raphson iteration
subroutine bs_solver(nr_fric_func_tau, xL, xR, x_acc, f, it, theta, tau_stick, sigma, Z, v, is, stat)
  !WARNING: This is not a well tested portion !!!
  integer, intent(out) :: is
  double precision, intent(in) :: xL, xR, x_acc
  double precision :: xLeft, xRight, x_est
  double precision, intent(out) :: v
  double precision :: dfunc_dx, dx, dx_old, func_x, f_high, f_low, x_high, x_low
  double precision :: temp
  integer :: stat ! solution status
  
  ! Friction parameters:
  external nr_fric_func_tau
  type(rsf_type), intent(in) :: f
  double precision, intent(in) :: theta, tau_stick, sigma, Z
  integer, intent(in) :: it
  double precision:: df, f_pre

  stat  = 1

  ! Find initial function estimates (dfunc_dx not needed yet)
  call nr_fric_func_tau(xL, f_low, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  call nr_fric_func_tau(xR, f_high, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  xLeft = xL
  xRight = xR  

  ! Lucky guesses: 
  if (f_low==0)  then
    x_est=xLeft
    return
  else if (f_high==0) then
    x_est=xRight
    return
  ! Orient the search so that func_x(x_low) < 0
  else if (f_low<0) then 
    x_low=xLeft
    x_high=xRight
  else 
    x_high=xLeft
    x_low=xRight
  endif

  f_pre = 0d0
  df    = 0d0

  do is=1,f%NRMaxIter
      x_est = (x_low+x_high)/2d0

      ! evaluate the function
      call nr_fric_func_tau(x_est, func_x, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)

      df = func_x - f_pre

      if (func_x>=0) then
          x_high = x_est
      else
          x_low  = x_est
      end if
      
      if (abs(x_low-x_high)<abs(x_acc) .and. abs(df)<abs(x_acc)) then
          ! evaluate the objective function one last time
          x_est = (x_low+x_high)/2d0
          call nr_fric_func_tau(x_est, func_x, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
          exit
      end if
      f_pre = func_x

      if (is==f%NRMaxIter) then 
          write(*, *) "maximum bisect iteration reached, maxiter=", f%NRMaxIter
          write(*, *) "xlow, xhigh", x_low, x_high
          stat = -3
      end if

  end do
  
end subroutine

!==================================================================
!         Newton-Raphson algorithm with bisection step         
! from  Numerical Recipes in Fortran 90 by Press et al.  (1996)

! N-R method, must provide a function (nr_fric_func) that gives
! the value of the function and the value of the derivative of the 
! function at a point and it finds a root to 0=nr_fric_func bounded 
! by [x1, x2] w/ error < x_acc

  subroutine nr_solver(nr_fric_func_tau, xL, xR, x_acc, f, it, theta, tau_stick, sigma, Z, v, is, stat)
  !WARNING: This is not a well tested portion !!!
  integer, intent(out) :: is
  double precision, intent(in) :: xL, xR, x_acc
  double precision :: xLeft, xRight, x_est
  double precision, intent(out) :: v
  double precision :: dfunc_dx, dx, dx_old, func_x, f_high, f_low, x_high, x_low
  double precision :: temp, df, f_pre
  integer :: stat ! solution status

  !stat =  1 normal converge
  !stat = -1 diverge with max iteration number reached
  !stat = -2 diverge with NAN
  
  ! Friction parameters:
  external nr_fric_func_tau
  type(rsf_type), intent(in) :: f
  double precision, intent(in) :: theta, tau_stick, sigma, Z
  integer, intent(in) :: it

! for debugging purposes, check input values
!  if (theta<1d-18) then
!      write(*, *) "nr_solver input, it", it
!      write(*, *) "theta<1d-18, theta=", theta
!      write(*, *) "theta<1d-18, log10(theta)=", log10(theta)
!  end if

  ! Find initial function estimates (dfunc_dx not needed yet)
  call nr_fric_func_tau(xL, f_low, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  call nr_fric_func_tau(xR, f_high, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  xLeft = xL
  xRight = xR  

  ! Ensure zero is bounded:
  do while (f_low*f_high>0)
    xLeft = xLeft/2.0
    xRight = xRight*2
    call nr_fric_func_tau(xLeft, f_low, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
    call nr_fric_func_tau(xRight, f_high, dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  enddo
 
  ! Lucky guesses: 
  if (f_low==0)  then
    x_est=xLeft
    return
  else if (f_high==0) then
    x_est=xRight
    return
  ! Orient the search so that func_x(x_low) < 0
  else if (f_low<0) then 
    x_low=xLeft
    x_high=xRight
  else 
    x_high=xLeft
    x_low=xRight
  endif

  ! Initialize the guess and step sizes:
  x_est=.5d0*(xLeft+xRight)
  ! The stepsize before last
  dx_old=abs(xRight-xLeft)
  ! The last step:
  dx=dx_old
  
  call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
  f_pre = 0d0

  ! Loop over allowed iterations:
  do is=1,f%NRMaxIter

    df = func_x - f_pre
    f_pre = func_x
   
    ! Bisect if N-R out of range, or not decreasing fast enough:
    if( ((x_est-x_high)*dfunc_dx-func_x)*((x_est-x_low)*dfunc_dx-func_x)>0 .or. abs(2*func_x)>abs(dx_old*dfunc_dx) &
        .or. abs(dfunc_dx)>huge(0d0) .or. abs(func_x)>huge(0d0)) then
      dx_old=dx
      dx=0.5*(x_high-x_low)
      x_est=x_low+dx
      !  Check if change is negligible:
      if (x_low==x_est) return

    ! The Newton step is acceptable, move forward with algorithm:
    else
      dx_old=dx
      dx=func_x/dfunc_dx
      temp=x_est
      x_est=x_est-dx
      !  Check if change is negligible:
      if (temp==x_est) return
    endif
  
    ! Check convergence criterion:
    if (abs(dx)<abs(x_acc) .and. abs(df) < abs(x_acc))  then
        ! evaluate the function and update v, one last time
        call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)
        ! converge as normal
        stat = 1
        return
    end if
    
    if (isnan(x_est) .or. isnan(dx) .or. isnan(v) & 
        .or. isnan(sigma) .or. isnan(x_acc) .or. isnan(func_x) .or. isnan(dfunc_dx)) then
      print*,'tau = ',x_est
      print*,'vel = ',v
      print*,'delta = ',dx,' > ',x_acc
      print*,'func_x=', func_x
      print*,'dfunc_dx=', dfunc_dx
      print*,'theta', theta
      print*,'log10(theta)', log10(theta)
      print*,'it', it
      print*,'tau_stick', tau_stick
      print*,'sigma', sigma
      print*, 'xL', xL
      print*, 'xR', xR
      print*, 'Z', Z
      stat = -2 ! diverge with NAN
      return
    end if
    
    ! Evaluate function with new estimate of x
    call nr_fric_func_tau(x_est,func_x,dfunc_dx, v, f, theta, it, tau_stick, sigma, Z)

    ! Redefine the bounds for the next loop
    if (func_x<0) then
      x_low=x_est
    else 
      x_high=x_est
    endif
  enddo
  print*,'tau = ',x_est
  print*,'delta = ',dx,' > ',x_acc
  print*,'vel = ',v
  stat = -1
  ! maximum iteration number reached!
  ! call IO_abort('NR_Solver has exceeded the maximum iterations')
  
  end subroutine nr_solver

!---------------------------------------------------------------------
!-----------Friction function for Newton Raphson method---------------
! This function returns the value of the friction function and its 
! derivative evaluated at a particular shear stress.  

subroutine nr_fric_func_tau(tau, func_tau, dfunc_dtau, v, f, theta, it, tau_stick, sigma, Z) 
  !WARNING: This is not a well tested portion!
  double precision, intent(out) :: func_tau, dfunc_dtau, v
  double precision, intent(in) ::  tau_stick,sigma,Z,tau,theta
  type(rsf_type), intent(in) :: f
  double precision :: dv_dtau, tmp
  integer, intent(in) :: it

      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      !  mu = mus +a*log(v/Vstar) + b*log(theta*Vstar/Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Lapusta et al. (2000))
      !  solved in terms of v
      !tmp = (f%mus(it) +f%b(it)*log( f%Vstar(it)*theta/f%Dc(it)))
!      tmp = 2d0*f%Vstar(it)*exp(-tmp/f%a(it))
      !v = sinh(tau/(-sigma*f%a(it)))*tmp

      tmp = -(f%mus(it) +f%b(it)*log( f%Vstar(it)*theta/f%Dc(it)))/f%a(it)
      v = f%Vstar(it) * (exp(tau/(-sigma*f%a(it))+tmp) - exp(tau/(sigma*f%a(it)) + tmp))
      v = sign(v, tau_stick)
      func_tau = tau_stick - Z*v - tau
    
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used)
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      !  solved in terms of v, derivative wrt tau
      !dv_dtau = cosh(tau/(-sigma*f%a(it)))*tmp/(-sigma*f%a(it))
      dv_dtau = f%Vstar(it) * (exp(tau/(-sigma*f%a(it))+tmp) + exp(tau/(sigma*f%a(it)) + tmp)) /(-sigma*f%a(it))
      dfunc_dtau = -Z*dv_dtau - 1d0
end subroutine nr_fric_func_tau

subroutine nr_fric_func_v(v, func_v, dfunc_dv, f, it, tau_stick, sigma, Z) 
  !WARNING: This is not a well tested portion!
  double precision, intent(out) :: func_v, dfunc_dv
  double precision, intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(in) :: f
  double precision :: v, mu, dmu_dv
  integer, intent(in) :: it

      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      ! mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      mu = (f%mus(it)+f%b(it)*log(f%Vstar(it)*f%theta(it)/f%Dc(it)))/f%a(it)
      mu = f%a(it)*asinh(v/(2*f%Vstar(it))*exp(mu))
      func_v = tau_stick - Z*v - sigma*mu
      
      !  Kaneko et al. (2008) Eq. 12 (this is unphysical, so Eq. 15 is used):
      !  dmu_dv = f%a*/v 
      !  Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
      dmu_dv = exp((f%mus(it)+f%b(it)*log(f%Vstar(it)*f%theta(it)/f%Dc(it)))/f%a(it))/(2*f%Vstar(it))
      dmu_dv = dmu_dv*f%a(it)/(sqrt(1+(dmu_dv*v)**2))
      dfunc_dv = -Z - sigma*dmu_dv
  
end subroutine nr_fric_func_v

!---------------------------------------------------------------------
!-----------Adapts timestep for quasi-static algorithm----------------
! limiting time step from fault
!
! I follow the matlab code to pass in the average node distance
! dt in rsf type must also be updated
!
! this subroutine also determines if an earthquake is occurring

subroutine rsf_timestep(time,f,v,sigma,hnode,mu_star,vgmax)
  use time_evol, only : timescheme_type
  
  use constants, only: PI

  type(rsf_type), intent(inout) :: f
  double precision, intent(in) :: hnode(:) 
  double precision, intent(in) :: mu_star(:) 
  type(timescheme_type), intent(inout) :: time
  double precision, dimension(:), intent(in) :: v,sigma
  double precision ::vgmax, vmax_VW, vmax_VS

  double precision :: k, xi, chi, dti, tmp, max_timestep,vmax
  integer :: it

  ! determine if the timestepping scheme needs to be changed
  ! using the maximum slip velocity on the fault

  vmax  = maxval(abs(v)) 
  vmax_VW = maxval(abs(v), f%a<f%b) ! vmax for velocity weakeaning
  vmax_VS = maxval(abs(v), f%a>f%b) ! vmax for velocity strengtening

  dti   = 0.0d0

  if (vmax>f%vEQ) then
      if ((.not. time%isEQ) .and. &
          (time%time-f%tEqPrev)>f%minGap .and. (.not. time%EQStart)) then
          time%EQStart = .true.
          write(*,*) "Start another earthquake ... "
          write(*,*) "max sliprate:", vmax
          ! this is a new earthquake
          time%isEQ = .true.
      end if
  else
      if (time%isEQ .and. time%isDynamic .and. time%EQStart) then
          time%EQNum  = time%EQNum + 1
          time%EQStart= .false.
          write(*,*) "Finish EQNum: ", time%EQNum
          write(*,*) "max sliprate:", vmax
          f%tEqPrev   = time%time
      end if
      time%isEQ = .false.
  end if

  f%StepLock = max(f%StepLock - 1, 0)
  
  if ((time%isDynamic .and. vmax_VW<f%vmaxD2S .and. vmax_VS<f%vmaxD2S_VS .and. f%StepLock==0) .or. &
      ((.not. time%isDynamic) .and. (vmax_VW<f%vmaxS2D)) .or. & 
      ((.not. time%isDynamic) .and. (vmax_VW>f%vmaxS2D) .and. (f%StepLock>0))) then
      ! three cases of being static:
      ! 1: previous step dynamic + vmax < vmaxD2S + steplock==0 (normal switch)
      ! 2: previous step static + vmax < vmaxS2D (keep being static)
      ! 3: previous step static + vmax > vmaxS2D + steplock>0 (forced being static for minStep)

      if (time%isDynamic) then
          time%switch = .true.
          write(*,*) "max sliprate in velocity weakening:", vmax_VW
          write(*,*) "max sliprate in velocity strengthening:", vmax_VS
          write(*,*) "max global velocity:", vgmax
          write(*, *) "Switching from dynamic to static, EQNum = ", time%EQNum
          ! reset steplock to minStep
          f%StepLock = f%minStep
      else
          time%switch = .false.
      end if

      ! stay in static
      time%isDynamic = .false.

      ! calculate time step
      max_timestep = time%dtev_max
      
      xi = 0.5d0 ! xi critical

      do it=1,size(v)
        k = (PI/4d0)* mu_star(it)/hnode(it) ! cell (two nodes) stiffness 
        ! Determine xi, as in Lapusta et al, 2000
        tmp = k*f%Dc(it)/(f%a(it)*sigma(it))
        chi = (0.25d0)*(tmp - (f%b(it) - f%a(it))/f%a(it))**2 - tmp ! Eq. 15c
        if (chi >= 0.0d0) then
          xi = f%a(it)*sigma(it)/(k*f%Dc(it) - sigma(it)*(f%b(it) - f%a(it))) ! Eq. 15a
        else
          xi = 1 - sigma(it)*(f%b(it)-f%a(it))/(k*f%Dc(it)) ! Eq. 15b
        endif
        xi = min(xi, 0.5d0) ! Eq. 15a,b
        dti = xi*f%Dc(it)/abs(v(it)) * f%dtScale
        if (dti > 0.0d0 .and. dti < max_timestep) max_timestep = dti
      enddo
      
      ! ---------- Limit the time step increase by a factor ---------
      if (max_timestep>time%dt*time%dt_incf) max_timestep = time%dt*time%dt_incf

      time%dt = max_timestep 
      f%dt    = max_timestep
  else
      ! dynamic stepping
      ! if previous state was static switch to dynamic
      ! update number of earthquakes occurred so far 
      
      if (.not. time%isDynamic) then
          time%switch = .true.
          ! switch to Dynamic and minimal time step
          time%isDynamic = .true.
          write(*, *) "Switching from static to dynamic, EQNum = ", time%EQNum
          write(*,*) "max sliprate in velocity weakening:", vmax_VW
          write(*,*) "max sliprate in velocity strengthening:", vmax_VS
          ! reset steplock to minstep
          f%StepLock = f%minStep
      else
          time%switch    = .false.
      end if
      time%dt     = time%dt_min
      f%dt        = time%dt_min
  end if

  if (f%kind==1) f%coeft = exp(-f%dt/f%Tc)

end subroutine

  function rsf_vplate(f) result(vplate)
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(f%vplate)) :: vplate
      vplate = f%vplate
  end function 



end module bc_dynflt_rsf
