module mat_mass

  use prop_mat

  implicit none
  private

  ! for memory report
  integer, save :: MAT_MASS_mempro = 0

  public :: MAT_MASS_init, MAT_MASS_init_elem_prop, MAT_MASS_mempro

contains

!=====================================================================
  subroutine MAT_MASS_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_setProp(elem,'rho',ecoord,MAT_MASS_mempro)
  
  end subroutine MAT_MASS_init_elem_prop

!=====================================================================
!
! Compute and assemble the mass matrix
!
  subroutine MAT_MASS_init(mass,mat,grid,ndof)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_VolumeWeights
  use prop_mat
  use echo, only : iout,echo_init, fmt1,fmtok
  use fields_class, only : FIELD_add_elem

  double precision, pointer :: mass(:,:)
  type(sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  integer, intent(in) :: ndof
  
  double precision :: massloc(grid%ngll,grid%ngll)
  integer :: e

  if (echo_init) write(iout,fmt1,advance='no') 'Building the mass matrix'

  allocate(mass(grid%npoin,ndof))
  call storearray('mass',size(mass),idouble)
  mass = 0.d0

  do e = 1,grid%nelem
    call MAT_getProp(massloc,mat(e),'rho')
    massloc = massloc * SE_VolumeWeights(grid,e)
    call FIELD_add_elem( massloc, mass(:,1), grid%ibool(:,:,e) )
  enddo

  if (ndof==2) mass(:,2) = mass(:,1)

  if (echo_init) write(iout,fmtok)

  end subroutine MAT_MASS_init


end module mat_mass
