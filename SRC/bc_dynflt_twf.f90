! SEM2DPACK version 2.3.6 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bc_dynflt_twf

! BC_DYNFLT_TWF: Time weakening friction for dynamic faults
! Usually for Andrews nucleation procedure:
! Minimum rupture velocity V imposed through time-dependent weakening
! L is the (minimum) width of the imposed rupture front
! Only for t < T.
!

  implicit none
  private

  type twf_type
    private
    double precision :: X,Z,mus,mud,mu0,L,V,T
  end type twf_type

  public :: twf_type, twf_read, twf_mu

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_TWF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Time weakening friction for dynamic faults 
!          with prescribed rupture speed.
! SYNTAX : &BC_DYNFLT_TWF MuS, MuD, Mu0, X, Z, V, L, T /
!
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
! ARG: Mu0      [dble] [0.6d0] Friction coefficient at the hypocenter at time=0
! ARG: X,Z      [dble] [0d0] Position of hypocenter
! ARG: V        [dble] [1d3] Rupture propagation speed
! ARG: L        [dble] [1d0] Size of weakening zone 
! ARG: T        [dble] [huge] Total duration
!
! END INPUT BLOCK

  subroutine twf_read(tw,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(twf_type), intent(out) :: tw
  integer, intent(in) :: iin

  double precision :: mus,mud,mu0,X,Z,L,V,T

  NAMELIST / BC_DYNFLT_TWF / mus,mud,mu0,X,Z,L,V,T

  mus = 0.6d0
  mud = 0.5d0
  mu0 = 0.6d0
  X = 0d0
  Z = 0d0
  V = 1d3
  L = 1d0
  T = huge(1d0)

  read(iin,BC_DYNFLT_TWF,END=100)
100 continue

  if (mus<0d0) call IO_abort('BC_SWFF_init: MuS must be positive in BC_DYNFLT_TWF input block')
  if (mud<0d0) call IO_abort('BC_SWFF_init: MuD must be positive in BC_DYNFLT_TWF input block')
  if (mu0<0d0) call IO_abort('BC_SWFF_init: Mu0 must be positive in BC_DYNFLT_TWF input block')
  if (L<=0d0) call IO_abort('BC_SWFF_init: L must be positive in BC_DYNFLT_TWF input block')
  if (V<=0d0) call IO_abort('BC_SWFF_init: V must be positive in BC_DYNFLT_TWF input block')
  if (T<=0d0) call IO_abort('BC_SWFF_init: T must be positive in BC_DYNFLT_TWF input block')

  tw%mus = mus
  tw%mud = mud
  tw%mu0 = mu0
  tw%X = X
  tw%Z = Z
  tw%L = L
  tw%T = T
  tw%V = V

  if (echo_input) write(iout,200) mus,mud,mu0,X,Z,V,L,T

  return
  200 format(5x,'Friction law  . . . . . . . . . . . . . .  = time weakening', &
            /5x,'  Static friction coefficient . . . .(MuS) = ',EN13.3, &
            /5x,'  Dynamic friction coefficient  . . .(MuD) = ',EN13.3, &
            /5x,'  Initial friction coefficient  . . .(Mu0) = ',EN13.3, &
            /5x,'  Hypocenter position X . . . . . . . .(X) = ',EN13.3, &
            /5x,'  Hypocenter position Z . . . . . . . .(Z) = ',EN13.3, &
            /5x,'  Propagation speed . . . . . . . . . .(V) = ',EN13.3, &
            /5x,'  Size of weakening zone  . . . . . . .(L) = ',EN13.3, &
            /5x,'  Duration  . . . . . . . . . . . . . .(T) = ',EN13.3 )

  end subroutine twf_read


!=====================================================================
  function twf_mu(tw,coord,time) result(mu)

  type(twf_type), intent(in) :: tw
  double precision, intent(in) :: coord(:,:),time
  double precision :: mu(size(coord,2))

  double precision :: r(size(coord,2)),t

  t = min(time,tw%T)
  r = sqrt( (coord(1,:)-tw%X)*(coord(1,:)-tw%X) + (coord(2,:)-tw%Z)*(coord(2,:)-tw%Z) )
  mu = tw%mu0 - (tw%mus-tw%mud)*(tw%V*t-r)/tw%L
  mu = max( mu, tw%mud )
  
  end function twf_mu

end module bc_dynflt_twf
