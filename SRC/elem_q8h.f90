module elem_Q8h

!=======================================================================
!
!                       This module deals with quadrilateral elements
!                       defined by 8 control nodes as follows:
!
!                                     edge 3
!
!                               4 . . 7 . . 8 . . 3
!                               .                 .
!                               .                 .
!                               .                 .
!                   edge 4      .                 .     edge 2
!                               .                 .
!                               .                 .
!                               .                 .
!                               1 . . 5 . . 6 . . 2
!
!                                     edge 1
!
!                                                         t
!                       The local coordinate system is :  .
!                       with (s,t) in [-1,1]              . . s
!
!=======================================================================
!
! Written by Yihe Huang (Caltech Seismolab) February 2010 

  implicit none
  private

  double precision, parameter :: one=1d0,two=2d0,three=3d0,eight=8d0,nine=9d0,eighteen=18d0,tseven=27d0,thtwo=32d0

  public :: Q8h_getshape,Q8h_getdershape

contains

function Q8h_getshape(s,t) result(shape)

  double precision, intent(in) :: s,t
  double precision :: shape(8)

  double precision :: sm,sp,tm,tp,ss,tt,ths

  sm = one - s  
  sp = one + s
  tm = one - t
  tp = one + t
  ss = s * s
  tt = t * t
  ths = s * three

  shape(1) = sm * tm * (nine * ss - one) / thtwo
  shape(2) = sp * tm * (nine * ss - one) / thtwo
  shape(3) = sp * tp * (nine * ss - one) / thtwo
  shape(4) = sm * tp * (nine * ss - one) / thtwo
  shape(5) = nine * (one - ss) * (one - ths) * tm / thtwo
  shape(6) = nine * (one - ss) * (one + ths) * tm / thtwo
  shape(7) = nine * (one - ss) * (one - ths) * tp / thtwo
  shape(8) = nine * (one - ss) * (one + ths) * tp / thtwo

end function Q8h_getshape

!-----------------------------------------------------------------------
function Q8h_getdershape(s,t) result(dershape)

  double precision, intent(in) :: s,t
  double precision :: dershape(8,2)

  double precision :: sp,sm,tp,tm,ss,tt,s2,ths

  sm = one - s  
  sp = one + s
  tm = one - t
  tp = one + t
  ss = s * s
  tt = t * t
  s2 = s * two
  ths = s * three

  dershape(1,1) = tm * (-tseven * ss + eighteen * s + one) / thtwo
  dershape(2,1) = tm * ( tseven * ss + eighteen * s - one) / thtwo
  dershape(3,1) = tp * ( tseven * ss + eighteen * s - one) / thtwo
  dershape(4,1) = tp * (-tseven * ss + eighteen * s + one) / thtwo

  dershape(1,2) = -sm * (nine * ss - one) / thtwo
  dershape(2,2) = -sp * (nine * ss - one) / thtwo
  dershape(3,2) =  sp * (nine * ss - one) / thtwo
  dershape(4,2) =  sm * (nine * ss - one) / thtwo

  dershape(5,1) = nine * tm * ( nine * ss - s2 - three) / thtwo
  dershape(6,1) = nine * tm * (-nine * ss - s2 + three) / thtwo
  dershape(7,1) = nine * tp * ( nine * ss - s2 - three) / thtwo
  dershape(8,1) = nine * tp * (-nine * ss - s2 + three) / thtwo

  dershape(5,2) = -nine * (one - ss) * (one - ths) / thtwo
  dershape(6,2) = -nine * (one - ss) * (one + ths) / thtwo
  dershape(7,2) =  nine * (one - ss) * (one - ths) / thtwo
  dershape(8,2) =  nine * (one - ss) * (one + ths) / thtwo


end function Q8h_getdershape

end module elem_Q8h
