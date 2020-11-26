module distribution_gaussian

  use stdio, only: IO_abort
  implicit none
  private
 
  type gaussian_dist_type
    private
    double precision :: x_0=0d0,z_0=0d0,lx=1d0,lz=1d0,level_0=0d0,ampli=1d0
    integer :: order=1
  end type gaussian_dist_type

  public :: gaussian_dist_type,read_gaussian_dist,generate_gaussian_dist ,&
            destroy_gaussian_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_GAUSSIAN
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Bell shaped Gaussian 2D distribution 
! SYNTAX : &DIST_GAUSSIAN centered_at, length, offset, ampli, order /
!
! ARG: centered_at      [dble(2)] [0,0] Coordinates of the center point.
! ARG: length           [dble(2)] [1]   Characteristic lengths on each axis.
! ARG: offset           [dble] [0]      Background level.    
! ARG: ampli            [dble] [1]      Amplitude from background.
! ARG: order            [int] [1]       Exponent
!
! END INPUT BLOCK


  subroutine read_gaussian_dist (d, file)

  type(gaussian_dist_type) :: d
  integer , intent(in) :: file

  double precision :: centered_at(2),length(2),offset,ampli
  integer :: order

  NAMELIST / DIST_GAUSSIAN / centered_at,length,offset,ampli,order

  centered_at = 0d0
  length = 1d0
  offset = 0d0
  ampli = 1d0
  order = 1

  read(file,DIST_GAUSSIAN,END=100)

  d%x_0 = centered_at(1)
  d%z_0 = centered_at(2)
  d%lx  = length(1)
  d%lz  = length(2)
  d%level_0 = offset
  d%ampli = ampli
  d%order = order

  return

  100 call IO_abort('read_gaussian_dist: DIST_GAUSSIAN parameters missing')

  end subroutine read_gaussian_dist
!
!***********************************************************************
!

  subroutine generate_gaussian_dist(field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(gaussian_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  !field = d%level_0 + d%ampli*exp(- ((coord(1,:)-d%x_0)/d%lx)**(2d0*d%order) &
  !                                - ((coord(2,:)-d%z_0)/d%lz)**(2d0*d%order) )
  field = d%level_0 + d%ampli*exp(- ( ((coord(1,:)-d%x_0)/d%lx)**2d0 &
                                     +((coord(2,:)-d%z_0)/d%lz)**2d0  )**d%order )
 
  end subroutine generate_gaussian_dist

!
!***********************************************************************
!
!! gaussian_dist_type destructor
subroutine destroy_gaussian_dist(d)
  type(gaussian_dist_type), pointer :: d
  deallocate(d)
end subroutine destroy_gaussian_dist

end module distribution_gaussian
