module distribution_pwconr

  implicit none
  private
 
  type pwconr_dist_type
    private
    integer :: NumZon ! number of zones
    double precision :: RefPnt(2) ! reference point
    double precision, pointer :: RadZon(:)  &! ext.radius of zone(i)
                                ,ValZon(:)   ! value inside zone(i)
  end type pwconr_dist_type

  public :: pwconr_dist_type,read_pwconr_dist,generate_pwconr_dist ,&
            destroy_pwconr_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_PWCONR
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Piecewise constant radial (2D) distribution.
!          This distribution defines a set of annular zones, centered
!          at an arbitrary reference point, and assigns constant values 
!          within each zone.
! SYNTAX : &DIST_PWCONR num, ref /
!             r(1)  ... ...  r(num-1)
!          v(1) v(2) ... v(num-1) v(num)
!
! ARG: num      [int] [none] Number of annular zones (including inner and exterior)
! ARG: ref      [dble(2)] [(0d0,0d0)] Reference point: center of radial zones
! ARG: r        [dble(num-1)] [none] External radius of zones:
!                first zone R <= r(1), 
!                second r(1) < R <= r(2), ...
!                last r(num-1) < R 
! ARG: v        [dble(num)] [none] Value inside each zone
!
! END INPUT BLOCK

  subroutine read_pwconr_dist (data, file)

  use stdio, only: IO_abort
  type(pwconr_dist_type), intent(out) :: data
  integer, intent(in) :: file
  integer :: num
  double precision :: ref(2)
  NAMELIST / DIST_PWCONR / num, ref

  ref = 0.d0
  read (file,DIST_PWCONR)
  if (num<2) call IO_abort('read_pwconr_dist: needs more than 2 zones (num)')
  data%NumZon = num
  data%RefPnt = ref

  allocate( data%RadZon(num-1) )
  read (file,*) data%RadZon
  allocate( data%ValZon(num) )
  read (file,*) data%ValZon
 
  end subroutine read_pwconr_dist
!
!***********************************************************************
!

  subroutine generate_pwconr_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(pwconr_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  double precision :: rad
  integer :: i,izone

  do i =1,size(field)
    rad = sqrt( (coord(1,i)-par%RefPnt(1))**2 +(coord(2,i)-par%RefPnt(2))**2 ) 
   !NOTE: when the loop goes to its end, izone=par%NumZon
    do izone=1,par%NumZon-1; if ( rad <= par%RadZon(izone) ) exit; end do
    field(i) = par%ValZon(izone)
  end do
 
  end subroutine generate_pwconr_dist

!
!***********************************************************************
!
!! pwconr_dist_type destructor
subroutine destroy_pwconr_dist(d)
  type(pwconr_dist_type), pointer :: d
  deallocate(d%RadZon,d%ValZon)
  deallocate(d)
end subroutine destroy_pwconr_dist

end module distribution_pwconr
