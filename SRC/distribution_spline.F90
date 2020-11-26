module distribution_spline

  use stdio, only: IO_abort, IO_new_unit, IO_file_length
  implicit none
  private
 
  type spline_dist_type
    private
    integer :: dim = 1
    double precision, pointer :: X(:)=>null(),val(:)=>null()
  end type spline_dist_type

  public :: spline_dist_type,read_spline_dist,generate_spline_dist ,&
            destroy_spline_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_SPLINE
! GROUP  : DISTRIBUTIONS_1D
! PURPOSE: Spline interpolated 1D distribution along X or Z.
! SYNTAX : &DIST_SPLINE file,dim /
!
! ARG: file     [name] [none] Name of the ASCII file containing
!               the interpolation data, one line per point, two columns: 
!               one line per point, two columns: 
!               position (X or Z), value
! ARG: dim      [int] [1] Interpolate along X (dim=1) or along Z (dim=2)
!
! END INPUT BLOCK

  subroutine read_spline_dist (data, iin)

  use utils, only: dsort

  type(spline_dist_type) :: data
  integer , intent(in) :: iin

  character(50) :: file
  integer :: iunit, i,N,dim

  NAMELIST / DIST_SPLINE / file,dim

  dim = 1

  read(iin,DIST_SPLINE)

  if (dim<1 .or. dim>2) call IO_abort('read_spline_dist: dim must be 1 or 2')
  data%dim = dim
  
  N = IO_file_length(file)
  allocate( data%X(N),data%val(N) )
  iunit = IO_new_unit()
  open(iunit,file=file,status='old')
  do i= 1,N
    read (iunit,*) data%x(i),data%val(i)
  end do
  close(iunit)

! sort X in ascending order, carry VAL
  call dsort(data%x,data%val) 

  end subroutine read_spline_dist

!
!***********************************************************************
!
  subroutine generate_spline_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(spline_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 

  integer :: unit,i 

  call interpol(par%X,par%val,coord(par%dim,:),field)

  unit = IO_new_unit()
  open(unit,file='DistSpline_sem2d.tab',status='replace')
  do i= 1,size(field)
    write(unit,*) coord(par%dim,i),field(i)
  enddo
  close(unit)
 
  end subroutine generate_spline_dist

!
!***********************************************************************
!
  subroutine interpol(xi,yi,xo,yo)

  use utils, only : spline, splint

  double precision, intent(in), dimension(:) :: xi,yi,xo
  double precision, intent(out) :: yo(:)

  double precision :: y2(size(xi))
  integer :: Ni,No,i

  Ni = size(xi)
  No = size(xo)
  ! WARNING: 0 derivative at extreme points of spline
  call spline(xi,yi,Ni,0d0,0d0,y2) 
  ! or natural spline:
  !call spline(xi,yi,size(xi),1d30,1d30,y2) 
  do i=1,No
    call splint(xi,yi,y2,Ni,xo(i),yo(i))
  enddo

  end subroutine interpol

!
!***********************************************************************
!
!! spline_dist_type destructor
subroutine destroy_spline_dist(d)
  type(spline_dist_type), pointer :: d
  deallocate(d%x,d%val)
  deallocate(d)
end subroutine destroy_spline_dist


end module distribution_spline
