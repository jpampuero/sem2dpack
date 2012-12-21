! SEM2DPACK version 2.3.3 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! http://www.seismolab.caltech.edu
! 
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
module mat_plastic
! plane strain Coulomb plasticity following Andrews (2005)

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
 !-- elasto-visco-plasticity
  type matwrk_plast_type
    private
!    double precision, pointer, dimension(:,:) :: mu => null(), lambda => null()
    double precision :: mu, lambda
    double precision, pointer, dimension(:,:,:) :: ep => null()
    double precision, pointer, dimension(:) :: e0=>null(), s0=>null()
    double precision :: yield_mu=0d0, yield_co=0d0, vp_factor=1d0
  end type matwrk_plast_type

  integer, save :: isPlastic = 0

  ! for memory report
  integer, save :: MAT_PLAST_mempro = 0
  integer, save :: MAT_PLAST_memwrk = 0

  public :: matwrk_plast_type &
          , MAT_isPlastic, MAT_anyPlastic, MAT_PLAST_read, MAT_PLAST_init_elem_prop &
          , MAT_PLAST_init_elem_work, MAT_PLAST_stress, MAT_PLAST_export &
          , MAT_PLAST_mempro, MAT_PLAST_memwrk
  
contains

!=======================================================================
  logical function MAT_isPlastic(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isPlastic = MAT_isKind(m,isPlastic)
  end function MAT_isPlastic

!=======================================================================
  logical function MAT_anyPlastic()
  MAT_anyPlastic = (isPlastic>0)
  end function MAT_anyPlastic

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_PLASTIC
! GROUP  : MATERIALS
! PURPOSE: Set material properties for elasto-plastic material
!          with Mohr-Coulomb yield criterion
!          and non-dilatant (null volumetric plastic strain)
! SYNTAX : &MAT_PLASTIC cp,cs,rho,phi,coh,Tv,e0,ep /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: phi      [dble][0d0] internal friction angle
! ARG: coh      [dble][0d0] cohesion
! ARG: Tv       [dble][0d0] Maxwellian visco-plastic timescale
! ARG: e0       [dble(4)][0d0] initial total strain (11, 22 and 12)
! ARG: ep       [dble(4)][0d0] initial plastic strain (11, 22 and 12)
!
! END INPUT BLOCK

  subroutine MAT_PLAST_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: rho,cp,cs,phi,coh,Tv,e0(3),ep(3)

  NAMELIST / MAT_PLASTIC / cp,cs,rho,phi,coh,Tv,e0,ep
  
  cp  = 0d0
  cs  = 0d0
  rho = 0d0
  e0  = 0d0
  ep  = 0d0
  phi = 0d0
  coh = 0d0
  Tv  = 0d0

  read(iin, MAT_PLASTIC, END=100)

  write(iout,200) cp,cs,rho,phi,coh,Tv,e0,ep

  call MAT_setKind(input,isPlastic)
  
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'phi',phi)
  call MAT_setProp(input,'coh',coh)
  call MAT_setProp(input,'Tv',Tv)

  call MAT_setProp(input,'e11_0',e0(1))
  call MAT_setProp(input,'e22_0',e0(2))
  call MAT_setProp(input,'e12_0',e0(3))

  call MAT_setProp(input,'e11_p',ep(1))
  call MAT_setProp(input,'e22_p',ep(2))
  call MAT_setProp(input,'e12_p',ep(3))

  return

  100 call IO_abort('MAT_PLAST_read: MAT_PLASTIC input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Internal friction angle . . . . . . (phi) =',EN12.3,/5x, &
    'Cohesion. . . . . . . . . . . . . . (coh) =',EN12.3,/5x, &
    'Visco-plastic timescale . . . . . . .(Tv) =',EN12.3,/5x, &
    'Initial total strain 11 . . . . . (e0(1)) =',EN12.3,/5x, &
    '                     22 . . . . . (e0(2)) =',EN12.3,/5x, &
    '                     12 . . . . . (e0(3)) =',EN12.3,/5x, &
    'Initial plastic strain 11 . . . . (ep(1)) =',EN12.3,/5x, &
    '                       22 . . . . (ep(2)) =',EN12.3,/5x, &
    '                       12 . . . . (ep(3)) =',EN12.3)

  end subroutine MAT_PLAST_read

!=======================================================================
  subroutine MAT_PLAST_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision :: rho,cp,cs

  call MAT_setProp(elem,'cp',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'cs',ecoord,MAT_PLAST_mempro)

  call MAT_getProp(rho,elem,'rho')
  call MAT_getProp(cp,elem,'cp')
  call MAT_getProp(cs,elem,'cs')
  call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_PLAST_mempro)
  call MAT_setProp(elem,'mu',rho*cs*cs,MAT_PLAST_mempro)

  call MAT_setProp(elem,'phi',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'coh',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'Tv',ecoord,MAT_PLAST_mempro)

  call MAT_setProp(elem,'e11_0',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e22_0',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e12_0',ecoord,MAT_PLAST_mempro)

  call MAT_setProp(elem,'e11_p',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e22_p',ecoord,MAT_PLAST_mempro)
  call MAT_setProp(elem,'e12_p',ecoord,MAT_PLAST_mempro)

  end subroutine MAT_PLAST_init_elem_prop

!=======================================================================
  subroutine MAT_PLAST_init_elem_work(m,p,ngll,dt)

  use constants, only : PI

  type(matwrk_plast_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll
  double precision, intent(in) :: dt

  double precision :: phi,Tv,lambda,two_mu,e_el(3)

 ! elastic coefficients
!  allocate(m%lambda(n,n)) ! NOTE: uniform elastic moduli 
!  allocate(m%mu(n,n))     ! NOTE: uniform elastic moduli
  call MAT_getProp(m%lambda,p,'lambda')
  call MAT_getProp(m%mu,p,'mu')

 ! Coulomb yield function Y = max_shear_stress - cohesion*cos(phi) + mean_stress*sin(phi)
 ! where phi = internal friction angle
 ! The yield condition Y=0 is equivalent to  
 ! max_over_all_planes[ shear_stress + tan(phi)*normal_stress ] = cohesion
  call MAT_getProp(phi,p,'phi')
  phi = PI/180d0 *phi
  call MAT_getProp(m%yield_co,p,'coh')
  m%yield_co = m%yield_co * cos(phi)  ! store cohesion*cos(phi)
  m%yield_mu = sin(phi)

 ! visco-plastic coefficient
  call MAT_getProp(Tv,p,'Tv')
  if (Tv>0d0) then 
    m%vp_factor = 1d0-exp(-dt/Tv)
  else
    m%vp_factor = 1d0
  endif

 ! initial strain
  allocate(m%e0(3))
  call MAT_getProp(m%e0(1),p,'e11_0')
  call MAT_getProp(m%e0(2),p,'e22_0')
  call MAT_getProp(m%e0(3),p,'e12_0')

 ! initial plastic strain
  allocate(m%ep(ngll,ngll,3))
  call MAT_getProp(m%ep(:,:,1),p,'e11_p')
  call MAT_getProp(m%ep(:,:,2),p,'e22_p')
  call MAT_getProp(m%ep(:,:,3),p,'e12_p')

 ! initial stress
  allocate(m%s0(3))
  e_el = m%e0 - m%ep(1,1,:)
  lambda = m%lambda
  two_mu = 2d0*m%mu
  m%s0(1) = (lambda+two_mu)*e_el(1) + lambda*e_el(2)
  m%s0(2) = lambda*e_el(1) + (lambda+two_mu)*e_el(2)
  m%s0(3) = two_mu*e_el(3)

  MAT_PLAST_memwrk = MAT_PLAST_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%ep) &
                   + size(m%e0) &
                   + size(m%s0)
!                   + size(m%lambda) &
!                   + size(m%mu) &

  end subroutine MAT_PLAST_init_elem_work


!=======================================================================
! Constitutive law:
! compute relative stresses from given strains
! and (if requested) update plastic strain
  subroutine MAT_PLAST_stress(s,e,m,ngll,update)

  integer, intent(in) :: ngll
  double precision, intent(inout) :: e(ngll,ngll,3)
  double precision, intent(out) :: s(ngll,ngll,3)
  type (matwrk_plast_type), intent(inout) :: m
  logical, intent(in) :: update

  double precision :: lambda, two_mu
  double precision, dimension(ngll,ngll,3) :: sd,sd_tr
  double precision, dimension(ngll,ngll) :: sm,tau,Y
  integer :: i,j

 ! elastic strain
  e = e - m%ep

 ! trial (elastic) stresses 
  lambda = m%lambda
  two_mu = 2d0*m%mu
  s(:,:,1) = (lambda+two_mu)*e(:,:,1) + lambda*e(:,:,2)
  s(:,:,2) = lambda*e(:,:,1) + (lambda+two_mu)*e(:,:,2)
  s(:,:,3) = two_mu*e(:,:,3)

 ! continue only if request to update internal variables
  if (.not. update) return

 ! total stresses
  s(:,:,1) = s(:,:,1) + m%s0(1)
  s(:,:,2) = s(:,:,2) + m%s0(2)
  s(:,:,3) = s(:,:,3) + m%s0(3)

 ! mean stress
  sm = 0.5d0*( s(:,:,1) + s(:,:,2) )

 ! trial (elastic) deviatoric stresses
  sd_tr(:,:,1)  = s(:,:,1) - sm
  sd_tr(:,:,2)  = s(:,:,2) - sm
  sd_tr(:,:,3)  = s(:,:,3)

 ! maximum shear stress
  tau = sqrt( (0.5d0*(sd_tr(:,:,1)-sd_tr(:,:,2)))**2 +sd_tr(:,:,3)**2 )

 ! Coulomb yield stress
  Y = m%yield_co - m%yield_mu *sm 

 ! test yield and update deviatoric stress
  do j=1,ngll
  do i=1,ngll
    if (tau(i,j)>Y(i,j)) then 
     ! visco-plastic update of deviatoric stresses
     ! where vp_factor = 1-exp(-dt/Tvisc)
     ! sd(i,j,:) = Y(i,j)/tau(i,j)*m%vp_factor* sd_tr(i,j,:) ! Andrews (2005)
      sd(i,j,:) = ( 1d0-(1d0-Y(i,j)/tau(i,j))*m%vp_factor )* sd_tr(i,j,:) ! correction by Duan (2008)
    else
      sd(i,j,:) = sd_tr(i,j,:)
    endif
  enddo
  enddo

 ! update plastic strain
  m%ep = m%ep + (sd_tr - sd)/two_mu

 ! recompose stresses (mean + deviatoric)
  s(:,:,1) = sd(:,:,1) + sm
  s(:,:,2) = sd(:,:,2) + sm
  s(:,:,3) = sd(:,:,3)

 ! relative stresses
  s(:,:,1) = s(:,:,1) - m%s0(1)
  s(:,:,2) = s(:,:,2) - m%s0(2)
  s(:,:,3) = s(:,:,3) - m%s0(3)

  end subroutine MAT_PLAST_stress

!=======================================================================
! export output data

  function MAT_PLAST_export(m) result(dat)

  type(matwrk_plast_type), intent(in) :: m
  real :: dat(size(m%ep,1),size(m%ep,2),size(m%ep,3))

  dat = real(m%ep)
  
  end function MAT_PLAST_export


end module mat_plastic
