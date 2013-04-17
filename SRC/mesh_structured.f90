module mesh_structured

  implicit none
  private

  public :: MESH_STRUCTURED_connectivity, MESH_STRUCTURED_boundaries &
           ,MESH_STRUCTURED_renumber

contains

! Set connectivity table: indices of nodes for each element
subroutine MESH_STRUCTURED_connectivity(knods,nx,nz,ngnod,ezflt)

  use stdio, only : IO_abort
  use utils, only : sub2ind

  integer, intent(out) :: knods(ngnod,nx*nz)
  integer, intent(in) :: nx,nz,ngnod
  integer, intent(in), optional :: ezflt

  integer :: ie,je,i,j,k,nxp

 ! Elements are sequentially numbered horizontally from bottom-left to top-right
  if (ngnod==4) then
    nxp = nx+1 
    k = 0
    do j=1,nz
    do i=1,nx
      k = k + 1
      knods(1,k) = sub2ind(i,j,nxp)
      knods(2,k) = sub2ind(i+1,j,nxp)
      knods(3,k) = sub2ind(i+1,j+1,nxp)
      knods(4,k) = sub2ind(i,j+1,nxp)
    enddo
    enddo 

  elseif (ngnod==8) then
    nxp = 3*nx+1 
    k = 0
    do je=1,nz
    do ie=1,nx
      k = k + 1
      i = 3*(ie-1)+1
      j = je
      knods(1,k) = sub2ind(i,j,nxp)
      knods(2,k) = sub2ind(i+3,j,nxp)
      knods(3,k) = sub2ind(i+3,j+1,nxp)
      knods(4,k) = sub2ind(i,j+1,nxp)
      knods(5,k) = sub2ind(i+1,j,nxp)
      knods(6,k) = sub2ind(i+2,j,nxp)
      knods(7,k) = sub2ind(i+1,j+1,nxp)
      knods(8,k) = sub2ind(i+2,j+1,nxp)
    enddo
    enddo 

  elseif (ngnod==9) then
    nxp = 2*nx+1 
    k = 0
    do je=1,nz
    do ie=1,nx
      k = k + 1
      i = 2*(ie-1)+1
      j = 2*(je-1)+1
      knods(1,k) = sub2ind(i,j,nxp)
      knods(2,k) = sub2ind(i+2,j,nxp)
      knods(3,k) = sub2ind(i+2,j+2,nxp)
      knods(4,k) = sub2ind(i,j+2,nxp)
      knods(5,k) = sub2ind(i+1,j,nxp)
      knods(6,k) = sub2ind(i+2,j+1,nxp)
      knods(7,k) = sub2ind(i+1,j+2,nxp)
      knods(8,k) = sub2ind(i,j+1,nxp)
      knods(9,k) = sub2ind(i+1,j+1,nxp)
    enddo
    enddo 
  else
    call IO_abort('MESH_STRUCTURED_connectivity: ngnod must be 4 or 9')
  endif

  if (present(ezflt)) then
    if (ezflt>0) knods(:,nx*ezflt+1:) = knods(:,nx*ezflt+1:)+nxp
  endif

end subroutine MESH_STRUCTURED_connectivity

!=====================================================================
subroutine MESH_STRUCTURED_boundaries(bnds,nx,nz,ezflt,splitN)

  use fem_grid, only : edge_D,edge_R,edge_U,edge_L, & 
                       side_D,side_R,side_U,side_L
  use bnd_grid, only : bnd_grid_type
  use utils, only : sub2ind

  type(bnd_grid_type) :: bnds(:)
  integer, intent(in) :: nx,nz
  integer, intent(in), optional :: ezflt
  integer, intent(in), optional :: splitN

  integer :: i,j
  integer, parameter :: fault_D = 5
  integer, parameter :: fault_U = 6
  logical :: split_bottom

  split_bottom = .false.
  if (present(splitN)) then
    if (splitN>0) split_bottom = .true.
  endif

  if (split_bottom) then
    ! Down Left
    bnds(fault_D)%tag = fault_D
    bnds(fault_D)%nelem = splitN
    allocate(bnds(fault_D)%elem(splitN))
    allocate(bnds(fault_D)%edge(splitN))
    do i = 1,splitN
      bnds(fault_D)%elem(i) = sub2ind(i,1,nx)
    enddo
    bnds(fault_D)%edge = edge_D
 
    ! Down Right
    bnds(side_D)%tag = side_D
    bnds(side_D)%nelem = nx - splitN
    allocate(bnds(side_D)%elem((nx-splitN)))
    allocate(bnds(side_D)%edge((nx-splitN)))
    do i = 1, nx-splitN
      bnds(side_D)%elem(i) = sub2ind(splitN+i,1,nx)
    enddo
    bnds(side_D)%edge = edge_D

  else
    ! Down
    bnds(side_D)%tag = side_D
    bnds(side_D)%nelem = nx
    allocate(bnds(side_D)%elem(nx))
    allocate(bnds(side_D)%edge(nx))
    do i =1,nx
      bnds(side_D)%elem(i) = sub2ind(i,1,nx)
    enddo
    bnds(side_D)%edge = edge_D
  endif

 ! Right
  bnds(side_R)%tag = side_R
  bnds(side_R)%nelem = nz
  allocate(bnds(side_R)%elem(nz))
  allocate(bnds(side_R)%edge(nz))
  do j =1,nz
    bnds(side_R)%elem(j) = sub2ind(nx,j,nx)
  enddo
  bnds(side_R)%edge = edge_R

 ! Up
  bnds(side_U)%tag = side_U
  bnds(side_U)%nelem = nx
  allocate(bnds(side_U)%elem(nx))
  allocate(bnds(side_U)%edge(nx))
  do i =1,nx
    bnds(side_U)%elem(i) = sub2ind(i,nz,nx)
  enddo
  bnds(side_U)%edge = edge_U

 ! Left
  bnds(side_L)%tag = side_L
  bnds(side_L)%nelem = nz
  allocate(bnds(side_L)%elem(nz))
  allocate(bnds(side_L)%edge(nz))
  do j =1,nz
    bnds(side_L)%elem(j) = sub2ind(1,j,nx)
  enddo
  bnds(side_L)%edge = edge_L

  if (present(ezflt)) then
  if (ezflt>0) then
   ! Fault Up: 
    bnds(fault_U)%tag = fault_U
    bnds(fault_U)%nelem = nx
    allocate(bnds(fault_U)%elem(nx))
    allocate(bnds(fault_U)%edge(nx))
    do i =1,nx
      bnds(fault_U)%elem(i) = sub2ind(i,ezflt+1,nx)
    enddo
    bnds(fault_U)%edge = edge_D

   ! Fault Down
    bnds(fault_D)%tag = fault_D
    bnds(fault_D)%nelem = nx
    allocate(bnds(fault_D)%elem(nx))
    allocate(bnds(fault_D)%edge(nx))
    do i =1,nx
      bnds(fault_D)%elem(i) = sub2ind(i,ezflt,nx)
    enddo
    bnds(fault_D)%edge = edge_U

  endif
  endif

end subroutine MESH_STRUCTURED_boundaries

!=====================================================================
!-- Renumbering of elements by Reverse Cuthill-McKee algorithm --
!   to improve data locality (optimize cache usage)
! perm: new index --> old index
! perm_inv: old index --> new index

subroutine MESH_STRUCTURED_renumber(grid,nx,nz)

  use fem_grid, only : fem_grid_type, FE_reorder
  use rcmlib, only : genrcm, perm_inverse
  use utils, only : sub2ind_b

  type(fem_grid_type), intent(inout) :: grid
  integer, intent(in) :: nx,nz
  
  integer, allocatable :: adj(:), adj_row(:)
  integer, allocatable :: perm(:), perm_inv(:)
  integer :: eadj(8),nadj,i,j,k,e,nelem

  nelem = nx*nz
  allocate(perm(nelem)) 
  allocate(perm_inv(nelem)) 

  ! Create the element adjacency structure 
  allocate(adj_row(nelem+1))
  nadj =  8*(nz-2)*(nx-2)+10*(nx+nz-4)+12
  allocate(adj(nadj))
  e =0
  nadj=0
  do j=1,nz
  do i=1,nx
    e=e+1
    adj_row(e)=nadj+1
   !
   ! For element e(i,j) :
   !   . . . . . . .
   !   . . 6 7 8 . .
   !   . . 4 e 5 . .
   !   . . 1 2 3 . .
   !   . . . . . . .
   !
    eadj(1) = sub2ind_b(i-1,j-1,nx,nz) ! returns 0 if out of box
    eadj(2) = sub2ind_b(i  ,j-1,nx,nz)
    eadj(3) = sub2ind_b(i+1,j-1,nx,nz)
    eadj(4) = sub2ind_b(i-1,j  ,nx,nz)
    eadj(5) = sub2ind_b(i+1,j  ,nx,nz)
    eadj(6) = sub2ind_b(i-1,j+1,nx,nz)
    eadj(7) = sub2ind_b(i  ,j+1,nx,nz)
    eadj(8) = sub2ind_b(i+1,j+1,nx,nz)
    do k=1,8
      if (eadj(k)>0) then
        nadj=nadj+1
        adj(nadj)=eadj(k)
      endif
    enddo
  enddo
  enddo
  adj_row(e+1)=nadj+1

  ! 2. Reverse Cuthill-McKee and permutation table
  perm=0
  call genrcm( nelem, nadj, adj_row, adj, perm)
  deallocate(adj_row,adj)
  perm_inv=0
  call perm_inverse ( nelem, perm, perm_inv )

  ! 3. Reorder grid data
  call FE_reorder(grid,perm,perm_inv)

  deallocate(perm,perm_inv)

end subroutine MESH_STRUCTURED_renumber

end module mesh_structured
