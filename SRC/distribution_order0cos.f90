module distribution_order0cos

  implicit none
  private
 
  type order0cos_dist_type
    private
    integer :: x_nzones, dim
    double precision, pointer, dimension(:) :: x_bound
    double precision, pointer, dimension(:) :: L
    double precision, pointer, dimension(:) :: val
  end type order0cos_dist_type

  public :: order0cos_dist_type,read_order0cos_dist,generate_order0cos_dist ,&
            destroy_order0cos_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_ORDER0COS
! GROUP  : DISTRIBUTIONS_1D
! PURPOSE: Piecewise constant 1D distribution with cosine taper transition
! SYNTAX : &DIST_ORDER0COS xn, dim /
!          x(1) ...  x(xn-1)
!          L(1)  ... L(xn-1)          
!          v(1)  ... v(xn)          
!
! ARG: xn       [int] [none] Number of zones along X
! ARG: dim      [int] [1] x: dim=1, z: dim=2
! ARG: x        [dble(xn-1)] [none] Boundaries of X-zones: 
!                first zone  X < x(1), 
!                second zone x(1) < X < x(2), ... 
!                last zone   x(xn-1) < X
! ARG: L        [dble(xn-1)] [none] transition length 
! ARG: v        [dble(xn)] [none] Values inside each zone
!
!
!    Between Zone i and Zone i + 1:
!          x < xi - Li/2: v(x) = v(i)
!          x > xi + Li/2: v(x) = v(i+1)
!             otherwise : 
!         v(x) = (v(i) + v(i+1))/2 - (v(i)-v(i+1))/2 * sin(2*pi*(x-xi)/2L)
!
! END INPUT BLOCK

  subroutine read_order0cos_dist (data, file)

  type(order0cos_dist_type) :: data
  integer , intent(in) :: file
  integer :: xn, dim

  NAMELIST / DIST_ORDER0COS / xn, dim

  dim = 1

  read (file, DIST_ORDER0COS)

  data%x_nzones = xn
  data%dim      = dim

  allocate( data%x_bound(max(xn-1,1)) )
  allocate( data%L(max(xn-1,1)) )

  if ( xn>1 ) read (file,*) data%x_bound
  if ( xn>1 ) read (file,*) data%L

  allocate( data%val(xn))
  read (file,*) data%val 
 
  end subroutine read_order0cos_dist
!
!***********************************************************************
!

  subroutine generate_order0cos_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(order0cos_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  integer :: nod, zones(2)

  do nod =1,size(field)
      ! obtain the zone
      zones = zone(coord(par%dim, nod), & 
                  par%x_nzones, par%x_bound, par%L)
      ! evaluate the function
      field(nod) = feval(coord(par%dim, nod), &
                         par%x_bound(zones(1)),& 
                         par%val(zones(1)),    &
                         par%val(zones(2)),    &
                         par%L(zones(1)))
  end do
 
  end subroutine generate_order0cos_dist

!
!***********************************************************************
!

  function zone(coord, nzones, bound, L)
 
  integer :: nzones 
  integer :: zone(2) 
  double precision :: coord
  double precision :: bound(:)
  double precision :: L(:)

  integer :: k

  zone = nzones

  if (nzones==1) return

  do k = 1, nzones
    if ( coord < bound(k) ) exit
  end do

  if (coord < (bound(k)-L(k)/2d0)) then
      zone(1) = max(k-1, 1)
      zone(2) = k
  else
      zone(1) = k 
      zone(2) = min(k+1, nzones)
  end if

  end function zone

  function feval(x, xi, vm, vp, L) result(v)
      double precision:: x, xi, vm, vp, L
      double precision:: v
      double precision, parameter :: PI = 3.141592653589793d0

      if (x<(xi-L/2.0d0)) then
          v = vm
          return
      end if

      if (x>(xi+L/2.0d0)) then
          v = vp
          return
      end if

      ! x>=xi-L/2 and x<=xi+L/2 
      v = (vm+vp)/2.0d0 - (vm-vp)/2.0d0*sin(PI*(x-xi)/L) 

  end function feval
!
!***********************************************************************
!
!! order0cos_dist_type destructor
subroutine destroy_order0cos_dist(d)
  type(order0cos_dist_type), pointer :: d
  deallocate(d%x_bound, d%L, d%val)
  deallocate(d)
end subroutine destroy_order0cos_dist

end module distribution_order0cos
