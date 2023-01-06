module elem_Q9

!=======================================================================
!
!                       This module deals with quadrilateral elements
!                       defined by 9 control nodes as follows:
!
!                                       edge 3
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .                   .
!                               .                   .
!                   edge 4      8         9         6     edge 2
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                                       edge 1
!
!                                                         t
!                       The local coordinate system is :  .
!                       with (s,t) in [-1,1]              . . s
!
!=======================================================================

  implicit none
  private

  double precision, parameter :: two=2d0,one=1d0,half=0.5d0,quart=0.25d0

  public :: Q9_getshape,Q9_getdershape

contains

function Q9_getshape(s,t) result(shape)

  double precision, intent(in) :: s,t
  double precision :: shape(9)

  double precision :: sp,sm,tp,tm,ss,tt,st

  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one
  ss = s * s
  tt = t * t
  st = s * t

  shape(1) = quart * sm * st * tm
  shape(2) = quart * sp * st * tm
  shape(3) = quart * sp * st * tp
  shape(4) = quart * sm * st * tp

  shape(5) = half * tm * t * (one - ss)
  shape(6) = half * sp * s * (one - tt)
  shape(7) = half * tp * t * (one - ss)
  shape(8) = half * sm * s * (one - tt)

  shape(9) = (one - ss) * (one - tt)

end function Q9_getshape

!-----------------------------------------------------------------------
function Q9_getdershape(s,t) result(dershape)

  double precision, intent(in) :: s,t
  double precision :: dershape(9,2)

  double precision :: sp,sm,tp,tm,s2,t2,ss,tt,st

  sp = s + one
  sm = s - one
  tp = t + one
  tm = t - one
  ss = s * s
  tt = t * t
  st = s * t
  s2 = s * two
  t2 = t * two

  dershape(1,1) = quart * tm * t * (s2 - one)
  dershape(2,1) = quart * tm * t * (s2 + one)
  dershape(3,1) = quart * tp * t * (s2 + one)
  dershape(4,1) = quart * tp * t * (s2 - one)

  dershape(1,2) = quart * sm * s * (t2 - one)
  dershape(2,2) = quart * sp * s * (t2 - one)
  dershape(3,2) = quart * sp * s * (t2 + one)
  dershape(4,2) = quart * sm * s * (t2 + one)

  dershape(5,1) = -one  * st * tm
  dershape(6,1) =  half * (one - tt) * (s2 + one)
  dershape(7,1) = -one  * st * tp
  dershape(8,1) =  half * (one - tt) * (s2 - one)

  dershape(5,2) =  half * (one - ss) * (t2 - one)
  dershape(6,2) = -one  * st * sp
  dershape(7,2) =  half * (one - ss) * (t2 + one)
  dershape(8,2) = -one  * st * sm

  dershape(9,1) = -one * s2 * (one - tt)
  dershape(9,2) = -one * t2 * (one - ss)

end function Q9_getdershape

end module elem_Q9
