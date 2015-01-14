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

  public :: fields_type, FIELDS_read, FIELDS_init &
          , FIELD_get_elem, FIELD_add_elem, FIELD_strain_elem, FIELD_divcurl_elem

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
!    if (name == 'displ') then 
!      field = 0.5d-3 ! WARNING: For testing benchmark against Kaneko 2008
!    else 
    field = 0.d0
!    endif
    
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
function FIELD_strain_elem(Uloc,ngll,ndof,grid,e) result(eij)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian
  use mxmlib

  integer, intent(in) :: ngll,ndof
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: Uloc(ngll,ngll,ndof)
  integer             , intent(in) :: e

  double precision :: eij(ngll,ngll,ndof+1)

  double precision, dimension(ngll,ngll) :: dUy_dxi, dUy_deta,dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  double precision, dimension(2,2,ngll,ngll), target :: xjaci
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dz, deta_dx, deta_dz

!-- Local gradient
  if (ndof==1) then
    dUy_dxi  = mxm( grid%hTprime, Uloc(:,:,1), ngll )
    dUy_deta = mxm( Uloc(:,:,1), grid%hprime, ngll )
  else
    dUx_dxi  = mxm( grid%hTprime, Uloc(:,:,1), ngll )
    dUz_dxi  = mxm( grid%hTprime, Uloc(:,:,2), ngll )
    dUx_deta = mxm( Uloc(:,:,1), grid%hprime, ngll )
    dUz_deta = mxm( Uloc(:,:,2), grid%hprime, ngll )
  endif

!-- Jacobian matrix
  xjaci = SE_InverseJacobian(grid,e)
  dxi_dx  => xjaci(1,1,:,:)
  dxi_dz  => xjaci(1,2,:,:)
  deta_dx => xjaci(2,1,:,:)
  deta_dz => xjaci(2,2,:,:)

!-- Strain 
  if (ndof==1) then
    eij(:,:,1) = 0.5d0*( dUy_dxi*dxi_dx + dUy_deta*deta_dx ) ! e13
    eij(:,:,2) = 0.5d0*( dUy_dxi*dxi_dz + dUy_deta*deta_dz ) ! e23
  else
    eij(:,:,1) = dUx_dxi*dxi_dx + dUx_deta*deta_dx  ! e11
    eij(:,:,2) = dUz_dxi*dxi_dz + dUz_deta*deta_dz  ! e22
    eij(:,:,3) = 0.5d0*( dUx_dxi*dxi_dz + dUx_deta*deta_dz  &
                       + dUz_dxi*dxi_dx + dUz_deta*deta_dx  ) ! e12
  endif

end function FIELD_strain_elem


!=============================================================================
! divergence and curl, only for P-SV (ndof=2)
function FIELD_divcurl_elem(Uloc,ngll,grid,e) result(eij)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian
  use mxmlib

  integer, intent(in) :: ngll
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: Uloc(ngll,ngll,2)
  integer             , intent(in) :: e

  double precision :: eij(ngll,ngll,2)

  double precision, dimension(ngll,ngll) :: dUx_dxi, dUx_deta, dUy_dxi, dUy_deta
  double precision, dimension(2,2,ngll,ngll), target :: xjaci
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dy, deta_dx, deta_dy

!-- Local gradient
  dUx_dxi  = mxm( grid%hTprime, Uloc(:,:,1), ngll )
  dUy_dxi  = mxm( grid%hTprime, Uloc(:,:,2), ngll )
  dUx_deta = mxm( Uloc(:,:,1), grid%hprime, ngll )
  dUy_deta = mxm( Uloc(:,:,2), grid%hprime, ngll )

!-- Jacobian matrix
  xjaci = SE_InverseJacobian(grid,e)
  dxi_dx  => xjaci(1,1,:,:)
  dxi_dy  => xjaci(1,2,:,:)
  deta_dx => xjaci(2,1,:,:)
  deta_dy => xjaci(2,2,:,:)

 ! div = dUx/dx + dUy/dy
  eij(:,:,1) = dUx_dxi * dxi_dx + dUx_deta * deta_dx &
             + dUy_dxi * dxi_dy + dUy_deta * deta_dy

 ! curl = dUx/dy - dUy/dx
  eij(:,:,2) = dUx_dxi * dxi_dy + dUx_deta * deta_dy &
             - dUy_dxi * dxi_dx - dUy_deta * deta_dx  

end function FIELD_divcurl_elem

end module fields_class
