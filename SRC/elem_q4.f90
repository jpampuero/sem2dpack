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
module elem_Q4

!=======================================================================
!
!                       This module deals with quadrilateral elements
!                       defined by 4 control nodes as follows:
!
!                                     edge 3
!
!                               4 . . . . . . . . 3
!                               .                 .
!                               .                 .
!                               .                 .
!                   edge 4      .                 .     edge 2
!                               .                 .
!                               .                 .
!                               .                 .
!                               1 . . . . . . . . 2
!
!                                     edge 1
!
!                                                         t
!                       The local coordinate system is :  .
!                       with (s,t) in [-1,1]              . . s
!
!=======================================================================

  implicit none
  private

  double precision, parameter :: one=1d0,quart=0.25d0

  public :: Q4_getshape,Q4_getdershape

contains

function Q4_getshape(s,t) result(shape)

  double precision, intent(in) :: s,t
  double precision :: shape(4)

  double precision :: sp,sm,tp,tm 
  
  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one

  shape(1) = quart * sm * tm
  shape(2) = - quart * sp * tm
  shape(3) = quart * sp * tp
  shape(4) = - quart * sm * tp

end function Q4_getshape

!-----------------------------------------------------------------------
function Q4_getdershape(s,t) result(dershape)

  double precision, intent(in) :: s,t
  double precision :: dershape(4,2)

  double precision :: sp,sm,tp,tm 
  
  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one

  dershape(1,1) = quart * tm
  dershape(2,1) = - quart * tm
  dershape(3,1) =  quart * tp
  dershape(4,1) = - quart * tp

  dershape(1,2) = quart * sm
  dershape(2,2) = - quart * sp
  dershape(3,2) =  quart * sp
  dershape(4,2) = - quart * sm

end function Q4_getdershape

end module elem_Q4
