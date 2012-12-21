! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module kelvin_voigt

  use distribution_general, only: cd_type

  implicit none
  private

  type kv_input_type
    logical :: NormalizeByDT=.true.
    type(cd_type) :: eta
  end type kv_input_type

  type kelvin_voigt_type
    private
    double precision, pointer :: eta(:,:)=>null()
    type(kv_input_type) :: input
  end type kelvin_voigt_type

  public :: kelvin_voigt_type, KV_read, KV_init, KV_etamul

contains

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : KELVIN_VOIGT
! PURPOSE: Define Kelvin-Voigt viscosity properties (whole domain)
!          i.e. add damping term C*v = K*eta*v
!          (eta is a viscous time)
! SYNTAX : &KELVIN_VOIGT eta|etaH, ETAxDT /
!
! ARG: eta	[dble][0d0] Viscosity coefficient
! ARG: etaH	[char*][] If eta is distributed non uniformly
!		  give here the nameof the distribution,
!		  followed by a DIST_XXX input block.
! ARG: ETAxDT	[log][T] If eta is given in units of dt (timestep)
!
! NOTE: useful as artificial damping layer in fault zones to control
!	high frequency noise. Set eta=0.1*dt and a thickness of 4-5 GLL nodes.
!
! END INPUT BLOCK

  subroutine KV_read(kv,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort
  use distribution_general, only: DIST_Read_CD

  type (kelvin_voigt_type), pointer :: kv
  integer, intent(in) :: iin

  double precision :: eta
  character(20) :: etaH
  logical :: ETAxDT

  NAMELIST / KELVIN_VOIGT / eta, etaH, ETAxDT

  nullify(kv)

  eta = 0d0
  etaH = ' '
  ETAxDT = .true.

  rewind(iin)
  read(iin, KELVIN_VOIGT, END=100)

  allocate(kv)
  call DIST_Read_CD(eta,etaH,kv%input%eta,iin) 
  kv%input%NormalizeByDT = ETAxDT

  if (echo_input) write(iout,200) etaH,ETAxDT

  100 return
  200   format(//,' M a t e r i a l   s e t s :   K e l v i n - V o i g t', &
         /1x,54('='),//5x, &
         'Viscosity coefficient . . . . . . . . . . =',a,/5x, &
         'Normalized by dt. . . . . . . . . . . . . =',L1)

  end subroutine KV_read

!=======================================================================

  subroutine KV_init(kv,coord,dt,ndof)

  use memory_info
  use distribution_general, only: DIST_Init_CD

  type (kelvin_voigt_type), pointer :: kv
  double precision, intent(in) :: coord(:,:),dt
  integer, intent(in) :: ndof

  integer :: npoin,n
  double precision, pointer :: eta(:)
  
  if (.not. associated(kv)) return

  npoin = size(coord,2)
  allocate( kv%eta(npoin,ndof) )
  call storearray('eta',npoin*ndof,idouble)

  eta => kv%eta(:,1)
  call DIST_Init_CD(kv%input%eta,coord,eta)
  
  if (ndof==2) kv%eta(:,2) = eta
  if (kv%input%NormalizeByDT) kv%eta = dt*kv%eta
  
! TO DO:
! define local node numbers for the viscous layer
! so KV_etamul will work only on those
!  kv%ibool = unique( ibool(:,:,e_visc) )

  end subroutine KV_init

!=======================================================================

  function KV_etamul(kv,v) result(eta_v)

  type (kelvin_voigt_type) :: kv
  double precision, intent(in) :: v(:,:)
  double precision :: eta_v(size(v,1),size(v,2))

  eta_v = kv%eta*v

! --> add_eta_v
! i => kv%ibool
! d(i,:) = d(i,:) + kv%eta*v(i,:)

  end function KV_etamul

end module kelvin_voigt
