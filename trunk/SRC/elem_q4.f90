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
