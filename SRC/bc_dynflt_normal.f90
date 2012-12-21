! SEM2DPACK version 2.3.2 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! Phone: (626) 395-3429
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

  if (kind>3 .or. kind<1) call IO_abort('BC_SWFF_init: invalid kind in BC_DYNFLT_NOR input block')
  if (L<=0d0) call IO_abort('BC_SWFF_init: L must be positive in BC_DYNFLT_NOR input block')
  if (V<=0d0) call IO_abort('BC_SWFF_init: V must be positive in BC_DYNFLT_NOR input block')
  if (T<=0d0) call IO_abort('BC_SWFF_init: T must be positive in BC_DYNFLT_NOR input block')

  n%kind = kind
  n%L = L
  n%T = T
  n%V = V

  if (echo_input) then
    select case (n%kind)
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
    case (1) ! Coulomb friction
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
