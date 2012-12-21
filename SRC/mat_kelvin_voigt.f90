! SEM2DPACK version 2.3.4 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module mat_kelvin_voigt

  use prop_mat

  implicit none
  private

 !-- kelvin-voigt
  type matwrk_kv_type
    private
    double precision, pointer :: eta(:,:) => null()
  end type matwrk_kv_type
  
  integer, save :: isKelvinVoigt = 0
  logical, save :: NormalizeByDT=.true.

  ! for memory report
  integer, save :: MAT_KV_mempro = 0
  integer, save :: MAT_KV_memwrk = 0

  public :: matwrk_kv_type &
          , MAT_isKelvinVoigt, MAT_KV_read, MAT_KV_init_elem_prop &
          , MAT_KV_init_elem_work, MAT_KV_add_etav &
          , MAT_KV_mempro, MAT_KV_memwrk

contains

!=======================================================================
  logical function MAT_isKelvinVoigt(m)

  type(matpro_elem_type), intent(in) :: m

  MAT_isKelvinVoigt = MAT_isKind(m,isKelvinVoigt)

  end function MAT_isKelvinVoigt

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_KV
! GROUP  : MATERIALS
! PURPOSE: Set material properties for Kelvin-Voigt viscosity 
!          Adds a damping term C*v = K*eta*v, where eta is a viscous time
!          This produces attenuation with frequency-dependent quality factor
!          Q(f) = 1/(eta*2*pi*f)
!          This is often useful as an artificial damping layer in fault zones
!          to control high-frequency numerical artifacts, setting eta=0.1*dt 
!          and a layer thickness of 4 to 5 GLL nodes.
! SYNTAX : &MAT_KV eta|etaH, ETAxDT /
!	   followed by a DIST_XXX input block if etaH is present.
!
! ARG: eta	[dble][0d0] Viscosity coefficient
! ARG: ETAxDT	[log][T] If eta is given in units of dt (timestep)
!
! END INPUT BLOCK

! NOTE: courant_kv < courant_elastic / sqrt( 1+ 2*eta/dtc)
!       where dtc is the critical timestep for elastic medium with eta=0
!       (see maybe in Hughes's book)
 
  subroutine MAT_KV_read(input,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: eta
  character(20) :: etaH
  logical :: ETAxDT

  NAMELIST / MAT_KV / eta, etaH, ETAxDT

  eta = 0d0
  etaH = ' '
  ETAxDT = .true.

  read(iin, MAT_KV, END=100)

  call MAT_setKind(input,isKelvinVoigt)

  call MAT_setProp(input,'eta',eta,etaH,iin,etaH)
  NormalizeByDT = ETAxDT

  if (echo_input) write(iout,200) etaH,ETAxDT

  return
  
  100 call IO_abort('MAT_KV_read: MAT_KV input block not found')

  200   format(5x, &
    'Viscosity coefficient . . . . .(eta/etaH) = ',a,/5x, &
    'Normalized by dt. . . . . . . . .(ETAxDT) = ',L1)

  end subroutine MAT_KV_read

!=======================================================================
  subroutine MAT_KV_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_setProp(elem,'eta',ecoord,MAT_KV_mempro)

  end subroutine MAT_KV_init_elem_prop

!=======================================================================
  subroutine MAT_KV_init_elem_work(matwrk,matpro,ngll,dt)

  type(matwrk_kv_type), intent(inout) :: matwrk
  type(matpro_elem_type), intent(in) :: matpro
  double precision, intent(in) :: dt
  integer, intent(in) :: ngll

  if (.not. MAT_isKelvinVoigt(matpro)) return
  allocate( matwrk%eta(ngll,ngll) )
  call MAT_getProp(matwrk%eta, matpro,'eta')
  if (NormalizeByDT) matwrk%eta = dt * matwrk%eta

  MAT_KV_memwrk = MAT_KV_memwrk &
                + size( transfer(matwrk, (/ 0d0 /) )) &
                + size(matwrk%eta)

  end subroutine MAT_KV_init_elem_work

!=======================================================================

  subroutine MAT_KV_add_etav(d,v,m,ngll,ndof)

  integer, intent(in) :: ngll,ndof
  double precision, intent(inout) :: d(ngll,ngll,ndof)
  double precision, intent(in) :: v(ngll,ngll,ndof)
  type(matwrk_kv_type), intent(in) :: m

  integer :: i

  do i=1,ndof
    d(:,:,i) = d(:,:,i) + m%eta * v(:,:,i)
  enddo

  end subroutine MAT_KV_add_etav

end module mat_kelvin_voigt
