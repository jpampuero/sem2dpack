! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! This software is freely available for scientific research purposes. 
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
module distribution_general

! TO ADD A NEW DISTRIBUTION :
!
! 1. Write a distribution module (e.g. "distribution_new"):
!    Take as a template any of the existing distribution_*.f90 modules.
!    A good starting point is distribution_hete1.f90
!
! 2. Modify the current module (distribution_general.f90):
!    Follow the template of the "new" distribution in the commented lines,
!    replacing "new" by your own distribution name.
!
! 3. Modify the file Makefile.depend :
!    Add your distribution_new.o to the list of dependencies of distribution_general.o
!    
! 4. Recompile (make)

  use distribution_order0
  use distribution_gaussian
  use distribution_spline
  use distribution_linear
  use distribution_gradient
  use distribution_pwconr
  use distribution_hete1
!  use distribution_new

  use stdio, only: IO_abort

  implicit none
  private

  type distribution_type
    private
    integer :: kind = 0
!--- List here distribution types
    type (order0_dist_type)  , pointer :: order0 => null()
    type (gaussian_dist_type), pointer :: gaussian => null()
    type (spline_dist_type)  , pointer :: spline => null()
    type (linear_dist_type)  , pointer :: linear => null()
    type (gradient_dist_type), pointer :: gradient => null()
    type (pwconr_dist_type)  , pointer :: pwconr => null()
    type (hete1_dist_type)   , pointer :: hete1 => null()
!    type (new_dist_type)   , pointer :: new => null()
  end type distribution_type

 ! CD = Constant or Distribution
  type cd_type
    double precision :: c=0d0
    type(distribution_type), pointer :: d => null()
  end type cd_type

!--- List here the tags, must be different for each distribution type
  integer, parameter :: tag_order0   = 1 &
                       ,tag_gaussian = 2 &
                       ,tag_spline   = 3 &
                       ,tag_linear   = 4 &
                       ,tag_gradient = 5 &
                       ,tag_pwconr   = 6 &
                       ,tag_hete1    = 7
!                       ,tag_new = next_number

!                  dist_name(number_of_tags_above)
  character(10) :: dist_name(7) = (/ 'ORDER0    ','GAUSSIAN  ','SPLINE    ','LINEAR    ', &
                                     'GRADIENT  ','PWCONR    ','HETE1     ' /)  
!                                    'NEW       '

  public :: distribution_type,DIST_read,DIST_generate,DIST_destructor &
            ,cd_type,DIST_read_CD,DIST_Init_CD

contains

!=======================================================================
!! Reads the distribution parameters from input File
subroutine DIST_read(d,name,iin)

  type(distribution_type), intent(out) :: d
  integer , intent(in) :: iin
  character(*), intent(in) :: name

  integer :: i

  do i=1,size(dist_name)
    if (dist_name(i) == name) d%kind = i
  enddo
  if (d%kind == 0) call IO_abort('DIST_read: unknown distribution name')

  select case (d%kind)
    case(tag_order0)
      allocate(d%order0)
      call read_order0_dist (d%order0,iin)
    case(tag_gaussian)
      allocate(d%gaussian)
      call read_gaussian_dist (d%gaussian,iin)
    case(tag_spline)
      allocate(d%spline)
      call read_spline_dist(d%spline,iin)
    case(tag_linear)
      allocate(d%linear)
      call read_linear_dist(d%linear,iin)
    case(tag_gradient)
      allocate(d%gradient)
      call read_gradient_dist(d%gradient,iin)
    case(tag_pwconr)
      allocate(d%pwconr)
      call read_pwconr_dist(d%pwconr,iin)
    case(tag_hete1)
      allocate(d%hete1)
      call read_hete1_dist(d%hete1,iin)
!    case(tag_new)
!      allocate(d%new)
!      call read_new_dist(d%new,iin)
    case default
      call IO_abort( 'DIST_read: illegal kind')
  end select

end subroutine DIST_read

!=======================================================================
!! Generates a field using distribution parameters
subroutine DIST_generate (field, coord, d)

  double precision, intent(in), dimension(:,:) :: coord
  type(distribution_type), intent(in) :: d

  double precision, intent(out), dimension(:) :: field

  select case (d%kind)
    case(tag_order0)
      call generate_order0_dist(field,coord,d%order0)
    case(tag_gaussian)
      call generate_gaussian_dist(field,coord,d%gaussian)
    case(tag_spline)
      call generate_spline_dist(field,coord,d%spline)
    case(tag_linear)
      call generate_linear_dist(field,coord,d%linear)
    case(tag_gradient)
      call generate_gradient_dist(field,coord,d%gradient)
    case(tag_pwconr)
      call generate_pwconr_dist(field,coord,d%pwconr)
    case(tag_hete1)
      call generate_hete1_dist(field,coord,d%hete1)
!    case(tag_new)
!      call generate_new_dist(field,coord,d%new)
    case default
      call IO_abort( 'DIST_generate: illegal kind')
  end select

end subroutine DIST_generate

!---- matrix input/output
!subroutine DIST_generate_2 (field, coord, d)
!
!  double precision, intent(in), dimension(:,:) :: coord
!  type(distribution_type), intent(in) :: d
!  double precision, intent(out), dimension(:) :: field
!end subroutine DIST_generate

!=======================================================================
!! distribution destructor
subroutine DIST_destructor(d)

  type(distribution_type), intent(inout) :: d

  select case (d%kind)
    case(tag_order0)
      call destroy_order0_dist(d%order0)
    case(tag_gaussian)
      call destroy_gaussian_dist(d%gaussian)
    case(tag_spline)
      call destroy_spline_dist(d%spline)
    case(tag_linear)
      call destroy_linear_dist(d%linear)
    case(tag_gradient)
      call destroy_gradient_dist(d%gradient)
    case(tag_pwconr)
      call destroy_pwconr_dist(d%pwconr)
    case(tag_hete1)
      call destroy_hete1_dist(d%hete1)
!    case(tag_new)
!      call destroy_new_dist(d%new)
    case default
      call IO_abort( 'DIST_destroy: illegal kind')
  end select

end subroutine DIST_destructor


!=====================================================================

  subroutine DIST_Read_CD(C,D,CD,iin)

  double precision, intent(in) :: C
  character(20), intent(inout) :: D
  type(cd_type), intent(out) :: CD
  integer, intent(in) :: iin

  if (D/='') then
    allocate(CD%d)
    call DIST_read(CD%d,D,iin)
  else
    CD%c = C
    write(D,'(EN13.3)') C
  endif

  end subroutine DIST_Read_CD

!---------------------------------------------------------------------

  subroutine DIST_Init_CD(CD,Coord,F)

  type(cd_type), intent(inout) :: CD
  double precision, intent(in) :: Coord(:,:)
  double precision, pointer :: F(:)

  if (.not.associated(F)) allocate( F(size(Coord,2)) )
  if (associated(CD%d)) then
    call DIST_generate(F,Coord,CD%d)
    call DIST_destructor(CD%d)
  else
    F = CD%c
  endif

  end subroutine DIST_Init_CD


end module distribution_general
