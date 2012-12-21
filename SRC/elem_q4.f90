! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module elem_Q4

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
