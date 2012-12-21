! SEM2DPACK version 2.2.12e -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
! This software is freely available for academic research purposes. 
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
module fields_class
! Fields (scalar or vector) with nodal storage

  implicit none
  private

  type fields_type
    integer :: ndof,npoin
    double precision, dimension(:,:), pointer ::  displ,veloc,accel, &
                                                  displ_alpha,veloc_alpha
  end type fields_type

  interface FIELD_get_elem
    module procedure FIELD_get_elem_1, FIELD_get_elem_2
  end interface FIELD_get_elem

  interface FIELD_add_elem
    module procedure FIELD_add_elem_1, FIELD_add_elem_2
  end interface FIELD_add_elem

  interface FIELD_strain_elem
    module procedure FIELD_strain_elem_1, FIELD_strain_elem_2
  end interface FIELD_strain_elem

  public :: fields_type, FIELDS_read, FIELDS_init &
          , FIELD_get_elem, FIELD_add_elem, FIELD_strain_elem

contains

!===================================================================================
!
subroutine FIELDS_read(fields,file)

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit

  type(fields_type), intent(inout) :: fields
  character(50), intent(in) :: file

  integer :: n,inump,iunit
  
  if (echo_input) write(iout,'(/,"Reading initial fields from external file",/)')

  iunit = IO_new_unit()
  open(iunit,file=trim(file),status='old')
  do n = 1,fields%npoin
    read(iunit,*) inump, fields%displ(inump,:), fields%veloc(inump,:),fields%accel(inump,:)
  enddo
  close(iunit)
  
end subroutine FIELDS_read


!===================================================================================
!
subroutine FIELDS_init(fields,npoin,alpha)

  type(fields_type), intent(out) :: fields
  integer, intent(in) :: npoin
  logical, intent(in) :: alpha

  fields%npoin = npoin
  call FIELD_init(fields%displ,'displ')
  call FIELD_init(fields%veloc,'veloc')
  call FIELD_init(fields%accel,'accel')
  if (alpha) then
    call FIELD_init(fields%displ_alpha,'displ_alpha')
    call FIELD_init(fields%veloc_alpha,'veloc_alpha')
  else
    fields%displ_alpha => fields%displ
    fields%veloc_alpha => fields%veloc
  endif

contains

   !-- single field:
    subroutine FIELD_init(field,name)
  
    use memory_info
  
    double precision, pointer :: field(:,:)
    character(*), intent(in) :: name
  
    allocate(field(npoin,fields%ndof))
    call storearray(name,size(field),idouble)
    field = 0.d0
  
    end subroutine FIELD_init

end subroutine FIELDS_init

!==============================================================================
! Assemble element contribution to global field
subroutine FIELD_add_elem_1(fin,Fout,ibool)
  
  double precision, intent(in) :: fin(:,:)
  double precision, intent(inout) :: Fout(:)
  integer, intent(in) :: ibool(:,:)

  integer :: i,j,k,ngll

  ngll = size(ibool,1)
  do j=1,ngll
  do i=1,ngll
    k = ibool(i,j)
    Fout(k) = Fout(k) + fin(i,j)
  enddo
  enddo

end subroutine FIELD_add_elem_1

!-----------------------------------------------------------------------------
subroutine FIELD_add_elem_2(fin,Fout,ibool)
  
  double precision, intent(in) :: fin(:,:,:)
  double precision, intent(inout) :: Fout(:,:)
  integer, intent(in) :: ibool(:,:)

  integer :: i,j,k,ngll

  ngll = size(ibool,1)
  do j=1,ngll
  do i=1,ngll
    k = ibool(i,j)
    Fout(k,:) = Fout(k,:) + fin(i,j,:)
  enddo
  enddo

end subroutine FIELD_add_elem_2

!=============================================================================
! Get element contribution from a global field, in local element storage 
function FIELD_get_elem_1(Fin,ibool) result(fout)

  double precision, intent(in) :: Fin(:)
  integer, intent(in) :: ibool(:,:)

  double precision :: fout(size(ibool,1),size(ibool,2))

  integer :: i,j,k,ngll

  ngll = size(ibool,1)
  do j=1,ngll
  do i=1,ngll
    k = ibool(i,j)
    fout(i,j) = Fin(k)
  enddo
  enddo

end function FIELD_get_elem_1

!-----------------------------------------------------------------------------
function FIELD_get_elem_2(Fin,ibool) result(fout)

  double precision, intent(in) :: Fin(:,:)
  integer, intent(in) :: ibool(:,:)

  double precision :: fout(size(ibool,1),size(ibool,2),size(Fin,2))

  integer :: i,j,k,ngll

  ngll = size(ibool,1)
  do j=1,ngll
  do i=1,ngll
    k = ibool(i,j)
    fout(i,j,:) = Fin(k,:)
  enddo
  enddo

end function FIELD_get_elem_2

!=============================================================================
function FIELD_strain_elem_1(Uloc,grid,e) result(eij)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian

  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: Uloc(grid%ngll,grid%ngll)
  integer             , intent(in) :: e

  double precision :: eij(grid%ngll,grid%ngll,2)
  double precision, dimension(grid%ngll,grid%ngll) :: dU_dxi, dU_deta
  double precision, dimension(2,2,grid%ngll,grid%ngll), target :: xjaci
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dz, deta_dx, deta_dz


!-- Local gradient
    dU_dxi  = MATMUL( grid%hTprime, Uloc )
    dU_deta = MATMUL( Uloc, grid%hprime )

!-- Jacobian matrix
    xjaci = SE_InverseJacobian(grid,e)
    dxi_dx  => xjaci(1,1,:,:)
    dxi_dz  => xjaci(1,2,:,:)
    deta_dx => xjaci(2,1,:,:)
    deta_dz => xjaci(2,2,:,:)

!-- Strain 
    eij(:,:,1) = 0.5d0*( dU_dxi*dxi_dx + dU_deta*deta_dx ) ! e13 = dUy_dx
    eij(:,:,2) = 0.5d0*( dU_dxi*dxi_dz + dU_deta*deta_dz ) ! e23 = dUy_dz

end function FIELD_strain_elem_1

!-----------------------------------------------------------------------------
function FIELD_strain_elem_2(Uloc,grid,e) result(eij)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian

  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: Uloc(grid%ngll,grid%ngll,2)
  integer             , intent(in) :: e

  double precision :: eij(grid%ngll,grid%ngll,3)
  double precision, dimension(grid%ngll,grid%ngll) :: &
    dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  double precision, dimension(2,2,grid%ngll,grid%ngll), target :: xjaci
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dz, deta_dx, deta_dz

!-- Local gradient
    dUx_dxi  = MATMUL( grid%hTprime, Uloc(:,:,1) )
    dUz_dxi  = MATMUL( grid%hTprime, Uloc(:,:,2) )
    dUx_deta = MATMUL( Uloc(:,:,1), grid%hprime )
    dUz_deta = MATMUL( Uloc(:,:,2), grid%hprime )

!-- Jacobian matrix
    xjaci = SE_InverseJacobian(grid,e)
    dxi_dx  => xjaci(1,1,:,:)
    dxi_dz  => xjaci(1,2,:,:)
    deta_dx => xjaci(2,1,:,:)
    deta_dz => xjaci(2,2,:,:)

!-- Strain 
    eij(:,:,1) = dUx_dxi*dxi_dx + dUx_deta*deta_dx  ! e11
    eij(:,:,2) = dUz_dxi*dxi_dz + dUz_deta*deta_dz  ! e22
    eij(:,:,3) = 0.5d0*( dUx_dxi*dxi_dz + dUx_deta*deta_dz  &
                       + dUz_dxi*dxi_dx + dUz_deta*deta_dx  ) ! e12


end function FIELD_strain_elem_2

end module fields_class
