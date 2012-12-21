! SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://web.gps.caltech.edu/~ampuero/
! 
! This software is freely available for academic research purposes. 
! If you use this software in writing scientific papers include proper 
! attributions to its author, Jean-Paul Ampuero.
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
! 
module bc_dynflt_rsf

! BC_DYNFLT_RSF: rate and state dependent friction for dynamic faults
! DEVEL: module in progress ... add theta, write friction solver

  use distribution_cd
  use stdio, only: IO_abort

  implicit none
  private

  type rsf_input_type
    type(cd_type) :: dc, mus, a, b, Vstar
  end type rsf_input_type

  type rsf_type
    private
    integer :: kind
    double precision, dimension(:), pointer :: dc, mus, a, b, Vstar, theta, &
                                               Tc, coeft
    double precision :: dt
    type(rsf_input_type) :: input
  end type rsf_type

  public :: rsf_type, rsf_read, rsf_init, rsf_mu, rsf_solver

contains

!---------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_RSF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Velocity and state dependent friction
! SYNTAX : &BC_DYNFLT_RSF kind, Dc | DcH, Mus | MusH , 
!                         a | aH, b | bH, Vstar | VstarH /
!          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
!          arguments with suffix H, if present, in the order listed above.
!
! ARG: kind     [int] [1] Type of rate-and-state friction law:
!                       1 = strong velocity-weakening at high speed
!                           as in Ampuero and Ben-Zion (2008)
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: a        [dble] [0.01d0] Direct effect coefficient
! ARG: b        [dble] [0.02d0] Evolution effect coefficient
! ARG: Vstar    [dble] [1d0] Characteristic or reference slip velocity
!
! END INPUT BLOCK

! not implement yet:
!                       2 = logarithmic rate-and-state with aging state law
!                       3 = logarithmic rate-and-state with slip state law

! Read parameters from input file
  subroutine rsf_read(rsf,iin)

  use echo, only : echo_input,iout

  type(rsf_type), intent(out) :: rsf
  integer, intent(in) :: iin

  double precision :: Dc,MuS,a,b,Vstar
  character(20) :: DcH,MuSH,aH,bH,VstarH
  integer :: kind
  character(20) :: kind_txt

  NAMELIST / BC_DYNFLT_RSF / kind,Dc,MuS,a,b,Vstar,DcH,MuSH,aH,bH,VstarH

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  a = 0.01d0
  b = 0.02d0
  Vstar = 1d0;
  DcH = ''
  MuSH = ''
  aH = ''
  bH = ''
  VstarH = '';

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

  if (echo_input) write(iout,400) kind_txt,DcH,MuSH,aH,bH,VstarH

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = rate and state dependent', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Direct effect coefficient . . . . . .(a) = ',A,&
            /5x,'  Evolution effect coefficient  . . . .(b) = ',A,&
            /5x,'  Velocity scale  . . . . . . . . .(Vstar) = ',A)

  end subroutine rsf_read

!=====================================================================
! Initialize parameters
  subroutine rsf_init(rsf,coord,dt)

  type(rsf_type), intent(inout) :: rsf
  double precision, intent(in) :: coord(:,:),dt

  integer :: n

  call DIST_CD_Init(rsf%input%dc,coord,rsf%dc)
  call DIST_CD_Init(rsf%input%mus,coord,rsf%mus)
  call DIST_CD_Init(rsf%input%a,coord,rsf%a)
  call DIST_CD_Init(rsf%input%b,coord,rsf%b)
  call DIST_CD_Init(rsf%input%Vstar,coord,rsf%Vstar)

  n = size(coord,2)
  allocate(rsf%theta(n))
  allocate(rsf%Tc(n))
  allocate(rsf%coeft(n))
 !WARNING: theta initialization should be more general for Dieterich-Ruina rsf
!          Also needs option for input by user
  rsf%theta = 0d0 
  rsf%Tc = rsf%dc / rsf%Vstar
  rsf%coeft = exp(-dt/rsf%Tc)
  rsf%dt = dt

  end subroutine rsf_init

!=====================================================================
! Friction coefficient
  function rsf_mu(v,f) result(mu)

  double precision, dimension(:), intent(in) :: v
  type(rsf_type), intent(in) :: f
  double precision :: mu(size(v))

  select case(f%kind)
    case(1); mu = f%mus +f%a*v/(v+f%Vstar) - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3); mu = f%mus +f%a*log(v/f%Vstar) + f%b*log(f%theta*f%Vstar/f%Dc) 
  end select

  end function rsf_mu

!---------------------------------------------------------------------
! friction coefficient without the direct effect 
! (i.e. without the term that depends explicitly on slip velocity V)
  function rsf_mu_no_direct(f) result(mu)

  type(rsf_type), intent(in) :: f
  double precision :: mu(size(f%mus))

  select case(f%kind)
    case(1); mu = f%mus - f%b*f%theta/(f%theta+f%Dc) 
    case(2,3); mu = f%mus + f%b*log(f%theta*f%Vstar/f%Dc) 
  end select

  end function rsf_mu_no_direct

!=====================================================================
! Two passes:
!      1. update theta using Vold from the previous time step
!      2. solve for Vnew
!      3. update theta again, now using V=(Vnew+Vold)/2
!      4. solve again for Vnew
  subroutine rsf_solver(v,tau_stick,sigma,f,Z)

  double precision, dimension(:), intent(inout) :: v
  double precision, dimension(:), intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(inout) :: f

  double precision, dimension(size(v)) :: v_new,theta_new

  theta_new = rsf_update_theta(f%theta,v,f)
  v_new = rsf_update_V(tau_stick, sigma, f, Z)
  theta_new = rsf_update_theta(f%theta,0.5d0*(v+v_new),f)
  v_new = rsf_update_V(tau_stick, sigma, f, Z)

  v = v_new
  f%theta = theta_new

  end subroutine rsf_solver

!---------------------------------------------------------------------
! Update state variable (theta) assuming slip velocity (v) is known

  function rsf_update_theta(theta,v,f) result(theta_new)

  double precision, dimension(:), intent(in) :: v,theta
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(v)) :: theta_new

  select case(f%kind)
    case(1) 
     ! exact integration assuming constant V over the timestep
     ! Tc = Dc/Vstar
     ! coeft = exp(-dt/Tc)
      theta_new = theta*f%coeft +f%Tc*abs(v)*(1d0-f%coeft)

    case(2) 
     ! Kaneko et al (2008) eq 19
     ! theta_new = (theta-Dc/v)*exp(-v*dt/Dc) + Dc/v
      theta_new = f%Dc/abs(v)
      theta_new = (theta-theta_new)*exp(-f%dt/theta_new) + theta_new
    case(3) 
     ! Kaneko et al (2008) eq 20
     ! theta_new = Dc/v *(theta*v/Dc)**exp(-v*dt/Dc)
      theta_new = f%Dc/abs(v)
      theta_new = theta_new *(theta/theta_new)**exp(-f%dt/theta_new)
  end select

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

  function rsf_update_V(tau_stick,sigma,f,Z) result(v)
   
  double precision, dimension(:), intent(in) :: tau_stick,sigma,Z
  type(rsf_type), intent(in) :: f
  double precision, dimension(size(tau_stick)) :: v

  double precision :: tmp(size(tau_stick))

!  strength = -sigma*rsf_mu_no_direct(v,f) 
!  v = (tau_stick-strength)/Z
  v = (tau_stick +sigma*rsf_mu_no_direct(f))/Z ! if v<0 will stop

  select case(f%kind)
    case(1) 
      tmp = v -f%Vstar +sigma*f%a/Z
      v = 0.5d0*( tmp +sqrt(tmp*tmp +4d0*v*f%Vstar) )
      v = max(0d0,v)  ! arrest if v<0 
 
    case(2) 
      call IO_abort('rsf_solve: case 2 not implemented yet')
     !DEVEL: should call Newton-Raphson solver

    case(3) 
      call IO_abort('rsf_solve: case 3 not implemented yet')
     !DEVEL: should call Newton-Raphson solver

  end select

  end function rsf_update_V

end module bc_dynflt_rsf
