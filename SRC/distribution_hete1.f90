module distribution_hete1
  
  implicit none
  private
 
  type hete1_dist_type
    private
    integer :: nx,nz
    double precision :: x0,z0,dx,dz
    double precision, pointer :: val(:,:)
  end type hete1_dist_type

  public :: hete1_dist_type, read_hete1_dist, generate_hete1_dist ,&
            destroy_hete1_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_HETE1
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Linear interpolation of values from a regular 2D grid.
! SYNTAX : &DIST_HETE1 file, col /
!
! ARG: file     [name] [none] Name of the file containing the definition
!               of the regular grid and values at grid points.
!               The format of this ASCII file is:
!                 Line 1 :  ncol nx nz x0 z0 dx dz
!                   ncol  = [int] number of data columns 
!                   nx,nz = [2*int] number of nodes along x and z
!                   x0,z0 = [2*dble] bottom-left corner 
!                   dx,dz = [2*dble] spacing along x and z
!                 Line 2 to nx*nz+1 : [ncol*dble] values at grid points
!                   listed from left to right (x0 to x0+nx*dx), 
!                   then from bottom to top (z0 to z0+nz*dx)
! ARG: col      [int] [1] Column of the file to be read
!
! NOTE   : The same file can contain values for (ncol) different properties,
!          (e.g. rho, vp, vs) but each DIST_HETE1 block will read only one.
!
! NOTE   : Even if the original model domain has an irregular shape, 
!          the regular grid where input values are defined must be rectangular
!          and large enough to contain the whole model domain. 
!          The regular grid possibly contains buffer areas with dummy values. 
!          These dummy values should be assigned carefully (not random nor zero)
!          because SEM2D might use them during nearest-neighbor interpolation.
!
! END INPUT BLOCK

  subroutine read_hete1_dist (d, iin)

  use stdio, only: IO_abort, IO_new_unit

  type(hete1_dist_type) :: d
  integer , intent(in) :: iin 

  double precision, allocatable :: vread(:)
  integer :: iunit, i,j,ncol, col
  character(50) :: file

  NAMELIST / DIST_HETE1 / file,col

  file =''
  col = 1
  read(iin,DIST_HETE1, END = 100)

! Read the file
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  read (iunit,*) ncol,d%nx,d%nz,d%x0,d%z0,d%dx,d%dz
  allocate( d%val(d%nx,d%nz) )
  allocate( vread(ncol) )
  do j= 1,d%nz
  do i= 1,d%nx
    read (iunit,*) vread
    d%val(i,j) = vread(col)
  end do
  end do
  close(iunit)

  deallocate(vread)

  return

  100 call IO_abort('read_hete1_dist: DIST_HETE1 parameters missing')

  end subroutine read_hete1_dist

!
!***********************************************************************
!
! Linear interpolation between nearest-neighbor grid points
!
! A smoothing kernel (a weighting function) could be implemented here
! involving more neighbours
!
  subroutine generate_hete1_dist(field, coord, d)

  use stdio, only: IO_abort
  use constants, only: TINY_XABS

  double precision, intent(in), dimension(:,:) :: coord  ! ndime*npoin
  type(hete1_dist_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field 

  double precision :: xi,eta
  integer :: i,j,k,ip,jp

  if (    minval(coord(1,:))<d%x0-TINY_XABS  &
     .or. maxval(coord(1,:))>d%x0+(d%nx-1)*d%dx+TINY_XABS  &
     .or. minval(coord(2,:))<d%z0-TINY_XABS  &
     .or. maxval(coord(2,:))>d%z0+(d%nz-1)*d%dz+TINY_XABS ) then
    call IO_abort('generate_hete1_dist: input grid is too small')
  endif

  do k=1,size(coord,2)

    xi  = (coord(1,k)-d%x0)/d%dx
    i = floor(xi)
    if (i>=d%nx-1) i = d%nx-2
    if (i<0) i = 0
    xi  = xi -dble(i)
    i = i + 1

    eta = (coord(2,k)-d%z0)/d%dz
    j = floor(eta)
    if (j>=d%nz-1) j = d%nz-2
    if (j<0) j = 0
    eta = eta -dble(j) 
    j = j + 1

    ip=i+1
    jp=j+1

    field(k) = (1d0-xi)*(1d0-eta)* d%val(  i,  j)   &
             +       xi*(1d0-eta)* d%val( ip,  j)   &
             +       xi*     eta * d%val( ip, jp)   &
             + (1d0-xi)*     eta * d%val(  i, jp)

  enddo

  end subroutine generate_hete1_dist


!***********************************************************************
! hete1_dist_type destructor
subroutine destroy_hete1_dist(d)
  type(hete1_dist_type), pointer :: d
  deallocate(d%val)
  deallocate(d)
end subroutine destroy_hete1_dist

end module distribution_hete1
