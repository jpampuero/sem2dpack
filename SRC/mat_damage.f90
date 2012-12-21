! SEM2DPACK version 2.2.12beta -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                             with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics Group
! ETH Hönggerberg HPP O 13.1
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 44 633 2197 (office)
! +41 44 633 1065 (fax)
! 
! http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
! 
! 
! This software is freely available for scientific research purposes. 
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
module mat_damage
! Damage (formulation by Vladimir Lyakhovsky)

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
  integer, save :: isDamage = 0

  public :: MAT_isDamage, MAT_DMG_read, MAT_DMG_init_elem_prop &
          , MAT_DMG_init_work, MAT_DMG_f
  
contains

!=======================================================================
  logical function MAT_isDamage(m)

  type(matpro_elem_type), intent(in) :: m

  if (isDamage>0) then
    MAT_isDamage = MAT_isKind(m,isDamage)
  else
    MAT_isDamage = .false.
  endif

  end function MAT_isDamage

!=======================================================================
  subroutine MAT_DMG_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho,cp,cs,beta,Cd,Cv,R,e0(4),ep(4),phi,alpha 

  NAMELIST / MAT_DAMAGE / rho,cp,cs,beta,Cd,R,e0,ep,phi,alpha
  
  rho = 0d0
  cp = 0d0
  cs = 0d0
  Cd = 0d0
  R = 0d0
  beta = 0d0
  e0 = 0d0
  ep = 0d0
  phi = 0d0
  alpha = 0d0

  rewind(iin)
  read(iin, MAT_DAMAGE, END=100)

  call MAT_setKind(input,isDamage)
  
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'Cd',Cd)
  call MAT_setProp(input,'R',R)
  call MAT_setProp(input,'beta',beta)
  call MAT_setProp(input,'phi',phi)
  call MAT_setProp(input,'alpha',alpha)

  call MAT_setProp(input,'e11_0',e0(1))
  call MAT_setProp(input,'e22_0',e0(2))
  call MAT_setProp(input,'e33_0',e0(3))
  call MAT_setProp(input,'e12_0',e0(4))

  call MAT_setProp(input,'e11_p',ep(1))
  call MAT_setProp(input,'e22_p',ep(2))
  call MAT_setProp(input,'e33_p',ep(3))
  call MAT_setProp(input,'e12_p',ep(4))

  return
  100 call IO_abort('MAT_DMG_read: MAT_DAMAGE input block not found')

  end subroutine MAT_DMG_read

!=======================================================================
  subroutine MAT_DMG_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision :: rho,cp,cs

  call MAT_setProp(elem,'cp',ecoord)
  call MAT_setProp(elem,'cs',ecoord)

  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs))
  call MAT_setProp(elem,'mu',rho*cs*cs)

  call MAT_setProp(elem,'Cd',ecoord)
  call MAT_setProp(elem,'R',ecoord)
  call MAT_setProp(elem,'beta',ecoord)
  call MAT_setProp(elem,'phi',ecoord)
  call MAT_setProp(elem,'alpha',ecoord)

  call MAT_setProp(elem,'e11_0',ecoord)
  call MAT_setProp(elem,'e22_0',ecoord)
  call MAT_setProp(elem,'e33_0',ecoord)
  call MAT_setProp(elem,'e12_0',ecoord)

  call MAT_setProp(elem,'e11_p',ecoord)
  call MAT_setProp(elem,'e22_p',ecoord)
  call MAT_setProp(elem,'e33_p',ecoord)
  call MAT_setProp(elem,'e12_p',ecoord)

  end subroutine MAT_DMG_init_elem_prop

!=======================================================================
  subroutine MAT_DMG_init_work(matwrk,mat_elem,grid,ndof)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights &
                      , SE_elem_coord, SE_isFlat
  use echo, only : echo_init,iout,fmt1,fmtok

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elem_type), intent(inout) :: matwrk(grid%nelem)
  type(matpro_elem_type), intent(in) :: mat_elem(grid%nelem)
  integer, intent(in) :: ndof

  double precision :: xjaci(2,2,grid%ngll,grid%ngll)
  double precision :: lambda,mu,tmp,rg
  integer :: e,nelem,ngll

  if (isDamage==0) return

  if (echo_init) &
    write(iout,fmt1,advance='no') 'Defining work arrays for damage'

  if (ndof/=2) call IO_abort('MAT_DMG_init: ndof must be 2 (P-SV)')

  ngll = grid%ngll
  nelem = 0
  do e=1,grid%nelem
  if (MAT_isDamage(mat_elem(e))) then

    nelem = nelem +1

    matwrk(e)%H => grid%hprime
    matwrk(e)%Ht => grid%hTprime

    allocate(matwrk(e)%dxi_dx(ngll,ngll))
    allocate(matwrk(e)%dxi_dy(ngll,ngll))
    allocate(matwrk(e)%deta_dx(ngll,ngll))
    allocate(matwrk(e)%deta_dy(ngll,ngll))
    xjaci = SE_InverseJacobian(grid,e)
    matwrk(e)%dxi_dx  = xjaci(1,1,:,:)
    matwrk(e)%dxi_dy  = xjaci(1,2,:,:)
    matwrk(e)%deta_dx = xjaci(2,1,:,:)
    matwrk(e)%deta_dy = xjaci(2,2,:,:)

    allocate(matwrk(e)%weights(ngll,ngll))
    matwrk(e)%weights = SE_VolumeWeights(grid,e)

    allocate(matwrk(e)%lambda(ngll,ngll))
    allocate(matwrk(e)%mu(ngll,ngll))
    call MAT_getProp(lambda,mat_elem(e),'lambda')
    call MAT_getProp(mu,mat_elem(e),'mu')
    matwrk(e)%lambda = lambda
    matwrk(e)%mu = mu

    allocate(matwrk(e)%alpha(ngll,ngll))
    call MAT_getProp(matwrk(e)%alpha,mat_elem(e),'alpha')

    allocate(matwrk(e)%alpha_dot(ngll,ngll))
    matwrk(e)%alpha_dot = 0d0

    call MAT_getProp(tmp,mat_elem(e),'phi')
    matwrk(e)%xi_0 = xi_zero(tmp,lambda,mu)

    matwrk(e)%gamma_r = gamma_r(matwrk(e)%xi_0,lambda,mu)
!    matwrk(e)%gamma_r = 0d0 !test: no damage

    call MAT_getProp(matwrk(e)%beta,mat_elem(e),'beta')

    call MAT_getProp(matwrk(e)%Cd,mat_elem(e),'Cd')

    call MAT_getProp(tmp,mat_elem(e),'R')
    matwrk(e)%Cv = tmp / mu

    allocate(matwrk(e)%e0(4))
    call MAT_getProp(matwrk(e)%e0(1),mat_elem(e),'e11_0')
    call MAT_getProp(matwrk(e)%e0(2),mat_elem(e),'e22_0')
    call MAT_getProp(matwrk(e)%e0(3),mat_elem(e),'e33_0')
    call MAT_getProp(matwrk(e)%e0(4),mat_elem(e),'e12_0')

    allocate(matwrk(e)%s0(4))
    mu = mu + matwrk(e)%xi_0 * matwrk(e)%gamma_r * matwrk(e)%alpha(1,1)
    rg = matwrk(e)%gamma_r * matwrk(e)%alpha(1,1)**(1d0+matwrk(e)%beta)
    call compute_stress(matwrk(e)%s0,matwrk(e)%e0,lambda,mu &
                       ,rg, tmp,tmp)

    allocate(matwrk(e)%ep(ngll,ngll,4))
    call MAT_getProp(matwrk(e)%ep(:,:,1),mat_elem(e),'e11_p')
    call MAT_getProp(matwrk(e)%ep(:,:,2),mat_elem(e),'e22_p')
    call MAT_getProp(matwrk(e)%ep(:,:,3),mat_elem(e),'e33_p')
    call MAT_getProp(matwrk(e)%ep(:,:,4),mat_elem(e),'e12_p')

  endif
  enddo

  call storearray('matwrk:damage',nelem*(ngll*ngll*13+4*2+5),idouble)
  
  end subroutine MAT_DMG_init_work

!-----------------------------------------------------------------------
  elemental double precision function xi_zero(phi,lambda,mu)

  use constants, only : PI

  double precision, intent(in) :: phi,lambda,mu
  double precision :: q

  q = sin(phi* PI/180d0)
  q = q / (1d0- q/3d0)
  xi_zero = -sqrt(3d0) / sqrt(2d0*q*q*(lambda/mu +2d0/3d0)**2 +1d0)

  end function xi_zero

!-----------------------------------------------------------------------
  elemental double precision function gamma_r(xi0,lambda,mu)

  double precision, intent(in) :: xi0,lambda,mu
  double precision :: q,p

  q = (2d0*mu+3d0*lambda)/(3d0-xi0**2)
  p = 0.5d0*xi0*(q + lambda)
  gamma_r = p + sqrt( p*p + 2d0*mu*q )

  end function gamma_r

!=======================================================================
  subroutine MAT_DMG_f(f,d,m,ngll,dt)

  use utils, only: heaviside

  integer, intent(in) :: ngll
  double precision, intent(out):: f(ngll,ngll,2)
  double precision, intent(in) :: d(ngll,ngll,2)
  type (matwrk_elem_type), intent(inout) :: m
  double precision, intent(in) :: dt

  double precision, dimension(ngll,ngll,4) :: e,s
  double precision, dimension(ngll,ngll) :: &
    dUx_dxi,dUy_dxi,dUx_deta,dUy_deta, &
    i2, rl,rm,rg, xi, tmp1,tmp2, tau,vsdt
  double precision :: a2,dtp
  integer :: i,j

 !-- local gradient
  dUx_dxi  = MATMUL( m%Ht, d(:,:,1) )
  dUy_dxi  = MATMUL( m%Ht, d(:,:,2) )
  dUx_deta = MATMUL( d(:,:,1), m%H )
  dUy_deta = MATMUL( d(:,:,2), m%H )

 !-- total strain 
  e(:,:,1) = m%e0(1) + dUx_dxi * m%dxi_dx + dUx_deta * m%deta_dx
  e(:,:,2) = m%e0(2) + dUy_dxi * m%dxi_dy + dUy_deta * m%deta_dy
  e(:,:,3) = m%e0(3)
  e(:,:,4) = m%e0(4) + 0.5d0*( dUx_dxi * m%dxi_dy + dUx_deta * m%deta_dy  &
                             + dUy_dxi * m%dxi_dx + dUy_deta * m%deta_dx  )

 !-- elastic strain
  e = e - m%ep

 !-- effective moduli
  rl = m%lambda
  rm = m%mu + m%xi_0 * m%gamma_r * m%alpha
  rg = m%gamma_r * m%alpha**(1d0+m%beta)

  do j=1,ngll
  do i=1,ngll
    call compute_stress(s(i,j,:),e(i,j,:),rl(i,j),rm(i,j),rg(i,j) &
                       ,i2(i,j),xi(i,j))
  enddo
  enddo


 !-- damage evolution, only if xi > xi_0
  m%alpha_dot = m%Cd*i2*( (1d0+m%beta)*m%alpha**m%beta *xi - m%xi_0 ) &
              * heaviside(xi-m%xi_0)
  m%alpha = m%alpha + dt * m%alpha_dot

 !-- plasticity update
 ! Damage-related viscosity, if alpha_dot > 0
 ! Plastic strain rate = deij/dt = tij * Cv *dalpha/dt
 ! where tij = deviatoric stress
  tau = (s(:,:,1)+s(:,:,2)+s(:,:,3))/3d0
  vsdt = m%Cv * m%alpha_dot *dt
  m%ep(:,:,1) = m%ep(:,:,1) + (s(:,:,1)-tau) * vsdt
  m%ep(:,:,2) = m%ep(:,:,2) + (s(:,:,2)-tau) * vsdt
  m%ep(:,:,3) = m%ep(:,:,3) + (s(:,:,3)-tau) * vsdt
  m%ep(:,:,4) = m%ep(:,:,4) + s(:,:,4) * vsdt

 !-- relative stresses
  do i=1,4
    s(:,:,i) = s(:,:,i) - m%s0(i)
  enddo

 !-- internal forces
  tmp1 = -m%weights * ( m%dxi_dx  * s(:,:,1) + m%dxi_dy  * s(:,:,4) )
  tmp2 = -m%weights * ( m%deta_dx * s(:,:,1) + m%deta_dy * s(:,:,4) )
  f(:,:,1) = MATMUL(m%H,tmp1) + MATMUL(tmp2,m%Ht)

  tmp1 = -m%weights * ( m%dxi_dx  * s(:,:,4) + m%dxi_dy  * s(:,:,2) )
  tmp2 = -m%weights * ( m%deta_dx * s(:,:,4) + m%deta_dy * s(:,:,2) )
  f(:,:,2) = MATMUL(m%H,tmp1) + MATMUL(tmp2,m%Ht)

 !--
  a2 = maxval(m%alpha_dot)*dt
  dtp = 1d-7 * dt/ maxval(vsdt)

  end subroutine MAT_DMG_f

!-------------------------------------------------------------------
  subroutine compute_stress(s,e,rl,rm,rg,i2,xi)

  double precision, intent(in) :: e(4),rl,rm,rg
  double precision, intent(out) :: s(4),i2,xi

  double precision :: i1,si2,two_mue,p,q,d

 !-- invariants of the elastic strain tensor
  i1 = e(1) + e(2) + e(3)
  i2 = e(1)*e(1) +e(2)*e(2) +e(3)*e(3) +2d0*e(4)*e(4)
  si2 = sqrt(i2)
  if (si2<1d-10) then
    xi =0d0
  else
    xi = i1/si2
  endif

 !-- check if the damage is above the critical value (loss of convexity)
  two_mue = 2d0*rm-rg*xi
  if ( two_mue <= 0d0 ) &
    call IO_abort('MAT_DMG: damage exceeded critical value (2nd type)')
  p = -(4d0*rm+3d0*rl-3d0*rg*xi)
  q = two_mue*two_mue + two_mue*(3d0*rl-rg*xi) + rg*(rl*xi-rg)*(3d0-xi*xi)
  d = p*p/4d0 - q
  if (d <= 0d0) call IO_abort('mat_damage:elastic: discriminant < 0')
  if ( p/2d0+sqrt(d) >= 0d0 ) &
    call IO_abort('MAT_DMG: damage exceeded critical value (1st type)')

 !-- absolute stress tensor
  s(1) = rl*i1 - rg*si2 + two_mue*e(1)
  s(2) = rl*i1 - rg*si2 + two_mue*e(2)
  s(3) = rl*i1 - rg*si2 + two_mue*e(3)
  s(4) = two_mue*e(4)

  end subroutine compute_stress

end module mat_damage
