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
module bc_dynflt_swf

! BC_DYNFLT_SWF: slip weakening friction for dynamic faults

  use distribution_cd

  implicit none
  private

  type swf_input_type
    type(cd_type) :: dc, mus, mud
  end type swf_input_type

  type swf_type
    private
    integer :: kind
    double precision, dimension(:), pointer :: dc, mus, mud
    type(swf_input_type) :: input
  end type swf_type

  public :: swf_type, swf_read, swf_init, swf_mu

contains

!---------------------------------------------------------------------
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_SWF
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Slip-weakening friction
! SYNTAX : &BC_DYNFLT_SWF Dc | DcH, MuS | MuSH , MuD | MuDH /
!          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
!          arguments with suffix H, if present, in the order listed above.
!
! ARG: kind     [int] [1] Type of slip weakening function:
!                       1 = linear
!                       2 = exponential
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
!
! END INPUT BLOCK

! Read parameters from input file
  subroutine swf_read(swf,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(swf_type), intent(out) :: swf
  integer, intent(in) :: iin

  double precision :: Dc,MuS,MuD
  character(20) :: DcH,MuSH,MuDH
  integer :: kind
  character(20) :: kind_txt

  NAMELIST / BC_DYNFLT_SWF / kind,Dc,MuS,MuD,DcH,MuSH,MuDH

  kind = 1
  Dc = 0.5d0
  MuS = 0.6d0
  MuD = 0.5d0
  DcH = ''
  MuSH = ''
  MuDH = ''

  read(iin,BC_DYNFLT_SWF,END=300)
300 continue
  
  select case (kind)
    case(1); kind_txt = 'Linear'
    case(2); kind_txt = 'Exponential'
    case default; call IO_abort('BC_DYNFLT_SWF: invalid kind')
  end select
  swf%kind = kind
  
  call DIST_CD_Read(swf%input%Dc,Dc,DcH,iin,DcH)
  call DIST_CD_Read(swf%input%MuS,MuS,MuSH,iin,MuSH)
  call DIST_CD_Read(swf%input%MuD,MuD,MuDH,iin,MuDH)

  if (echo_input) write(iout,400) kind_txt,DcH,MuSH,MuDH

  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = slip weakening', &
            /5x,'  Type of weakening . . . . . . . . (kind) = ',A,&
            /5x,'  Critical slip . . . . . . . . . . . (Dc) = ',A,&
            /5x,'  Static friction coefficient . . . .(MuS) = ',A,&
            /5x,'  Dynamic friction coefficient  . . .(MuD) = ',A)

  end subroutine swf_read

!=====================================================================
! Initialize parameters
  subroutine swf_init(swf,coord)

  type(swf_type), intent(inout) :: swf
  double precision, intent(in) :: coord(:,:)

  call DIST_CD_Init(swf%input%dc,coord,swf%dc)
  call DIST_CD_Init(swf%input%mus,coord,swf%mus)
  call DIST_CD_Init(swf%input%mud,coord,swf%mud)

  end subroutine swf_init

!=====================================================================
! Friction coefficient
  function swf_mu(s,f) result(mu)

  double precision, dimension(:), intent(in) :: s
  type(swf_type), intent(in) :: f
  double precision :: mu(size(s))

 !-- linear slip weakening:
  if (f%kind==1) then
    mu = f%mus -(f%mus-f%mud)/f%dc *abs(s)
    mu = max( mu, f%mud)
  else 
 !-- exponential slip weakening:
    mu = f%mus -(f%mus-f%mud)*exp(-abs(s)/f%dc)
  endif

  end function swf_mu

end module bc_dynflt_swf
