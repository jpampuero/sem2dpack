module distribution_elliptic1D
  
  implicit none
  private
 
  type elliptic1D_dist_type
    private
    integer :: dim = 1
    double precision :: x0, L
  end type elliptic1D_dist_type

  public :: elliptic1D_dist_type, read_elliptic1D_dist, generate_elliptic1D_dist ,&
            destroy_elliptic1D_dist

  contains

!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_ELLIPTIC1D
! GROUP  : DISTRIBUTIONS_1D
! PURPOSE: Elliptical distribution in 1D as:
!              |x-x0|<=L:    d = sqrt{1 - ((x-x0)/L)^2} 
!              |x-x0|> L:    d = 0 
!
! SYNTAX : &DIST_ELLIPTIC1D dim, x0, L /
!
! ARG: dim      [int] [1]   Distribution along x (dim=1) or along z (dim=2) 
! ARG: x0       [dbl] [0d0] center point of the patch 
! ARG: L        [int] [0d0] half length of the patch 
!
! NOTE   : the distribution must be 1D along either x or z direction. 
!
! NOTE   : mainly used to build the elliptical slip profile of a 1D fault under 
!          constant stress drop, to verify the implementation of kinflt module
!
! END INPUT BLOCK

  subroutine read_elliptic1D_dist(d, iin)
  use stdio, only: IO_abort
  use constants, only: tiny_xabs
  type(elliptic1D_dist_type) :: d
  integer , intent(in) :: iin

  double precision :: x0, L
  integer :: dim

  NAMELIST / DIST_ELLIPTIC1D / dim, x0, L

  L  = 0d0
  x0 = 0d0
  dim = 1
  
  read(iin, DIST_ELLIPTIC1D, END = 100)

  if (dim<1 .or. dim>2) call IO_abort('read_elliptic1D_dist: dim must be 1 or 2')
  if (L<0d0 .or. abs(L)<tiny_xabs) then
	 call IO_abort('read_elliptic1D_dist: L is negative or L is too small')
  end if

  d%dim = dim
  d%L   = L
  d%x0  = x0
  return

  100 call IO_abort('read_elliptic1D_dist: DIST_ELLIPTIC1D parameters missing')

  end subroutine read_elliptic1D_dist

!
!***********************************************************************
!
  subroutine generate_elliptic1D_dist(field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord  ! ndime*npoin
  type(elliptic1D_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  where (abs(coord(d%dim,:) - d%x0)<=d%L)
	field = sqrt(1d0 - ((coord(d%dim,:) - d%x0)/d%L)**2) 
  elsewhere
    field = 0d0
  end where
             
  end subroutine generate_elliptic1D_dist

!***********************************************************************
! elliptic1D_dist_type destructor
subroutine destroy_elliptic1D_dist(d)
  type(elliptic1D_dist_type), pointer :: d
  deallocate(d)
end subroutine destroy_elliptic1D_dist

end module distribution_elliptic1D
