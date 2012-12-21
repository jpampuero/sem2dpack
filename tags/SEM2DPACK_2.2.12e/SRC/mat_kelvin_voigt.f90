! SEM2DPACK version 2.2.12e -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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

  integer, save :: isKelvinVoigt = 0
  logical, save :: NormalizeByDT=.true.

  public :: MAT_isKelvinVoigt, MAT_KV_read, MAT_KV_init_elem_prop &
          , MAT_KV_init_work, MAT_KV_add_etav

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
!	   possibly followed by a DIST_XXX input block.
!
! ARG: eta	[dble][0d0] Viscosity coefficient
! ARG: etaH	[char*][''] If eta is distributed non uniformly
!		  give here the name of the distribution (see DIST_XXX)
! ARG: ETAxDT	[log][T] If eta is given in units of dt (timestep)
!
! END INPUT BLOCK

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

  call MAT_setProp(elem,'eta',ecoord)

  end subroutine MAT_KV_init_elem_prop

!=======================================================================
  subroutine MAT_KV_init_work(matwrk,matpro,ngll,dt)

  use memory_info

  type(matwrk_elem_type), intent(inout) :: matwrk(:)
  type(matpro_elem_type), intent(in) :: matpro(:)
  double precision, intent(in) :: dt
  integer, intent(in) :: ngll

  integer :: nelem,e
 
  if (isKelvinVoigt==0) return

  nelem = 0
  do e=1,size(matwrk)
    if (.not. MAT_isKelvinVoigt(matpro(e))) cycle
    nelem = nelem +1
    allocate(matwrk(e)%eta(ngll,ngll))
    call MAT_getProp(matwrk(e)%eta, matpro(e),'eta')
    if (NormalizeByDT) matwrk(e)%eta = dt * matwrk(e)%eta
  enddo

  if (nelem>0) call storearray('eta',nelem*ngll*ngll,idouble)

  end subroutine MAT_KV_init_work

!=======================================================================

  subroutine MAT_KV_add_etav(fout,v,m)

  double precision, intent(inout) :: fout(:,:,:)
  double precision, intent(in) :: v(:,:,:)
  type(matwrk_elem_type), intent(in) :: m

  integer :: i

  do i=1,size(fout,3) 
    fout(:,:,i) = fout(:,:,i) + m%eta * v(:,:,i)
  enddo

  end subroutine MAT_KV_add_etav

end module mat_kelvin_voigt
