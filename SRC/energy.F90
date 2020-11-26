module energy

  implicit none
  private

  type energy_type ! stores macroscopic quantities derived from calculated fields
    double precision :: E_ep = 0d0 ! cumulative plastic energy
    double precision :: E_k = 0d0  ! kinetic energy
    double precision :: E_el = 0d0  ! elastic energy
    double precision :: sg(3) = 0d0  ! damage stress glut
    double precision :: sgp(3) = 0d0  ! plastic stress glut 
    integer :: iout_E = 0          ! output unit for energies
    integer :: iout_SG = 0          ! output unit for stress glut
  end type energy_type

  public :: energy_type, energy_init, energy_compute, energy_write, &
            stress_glut_init, stress_glut_write
  
contains

!=======================================================================
  subroutine energy_init(energy)

  use stdio, only : IO_new_unit

  type(energy_type), intent(inout) :: energy

  energy%iout_E = IO_new_unit()
  open(energy%iout_E,file='energy_sem2d.tab',status='replace')

  end subroutine energy_init

!=======================================================================
  subroutine stress_glut_init(energy)

  use stdio, only : IO_new_unit

  type(energy_type), intent(inout) :: energy

  energy%iout_SG = IO_new_unit()
  open(energy%iout_SG,file='stress_glut_sem2d.tab',status='replace')

  end subroutine stress_glut_init

!=======================================================================
  subroutine energy_compute(energy,matpro,matwrk,grid,fields)

  use prop_mat    , only : matpro_elem_type
  use mat_gen     , only : matwrk_elem_type
  use spec_grid   , only : sem_grid_type, SE_VolumeWeights
  use fields_class, only : fields_type,FIELD_get_elem
  use prop_mat, only : MAT_getProp

  type(energy_type), intent(inout) :: energy
  type(matpro_elem_type), intent(in) :: matpro(:)
  type(matwrk_elem_type), intent(in) :: matwrk(:)
  type(sem_grid_type), intent(in) :: grid
  type(fields_type), intent(in) :: fields

  double precision :: v(grid%ngll,grid%ngll,fields%ndof)
  double precision, dimension(grid%ngll,grid%ngll) :: rho, v2, w
  integer :: e

  energy%E_k = 0d0

  do e=1,grid%nelem

    v = FIELD_get_elem(fields%veloc,grid%ibool(:,:,e))
    if (fields%ndof==1) then
      v2 = v(:,:,1)*v(:,:,1)
    else
      v2 = v(:,:,1)*v(:,:,1) + v(:,:,2)*v(:,:,2)
    endif
    call MAT_getProp(rho,matpro(e),'rho')
    if (associated(matwrk(e)%derint)) then
      w = matwrk(e)%derint%weights
    else
      w = SE_VolumeWeights(grid,e)
    endif
    energy%E_k = energy%E_k + sum(w*rho*v2) 

  enddo

  energy%E_k = 0.5d0*energy%E_k
  
  end subroutine energy_compute

!=======================================================================
  subroutine energy_write(energy,time)

  type(energy_type), intent(in) :: energy
  double precision, intent(in) :: time

  write(energy%iout_E,'(4D24.16)') time, energy%E_ep, energy%E_k, energy%E_el
  
  end subroutine energy_write

!=======================================================================
  subroutine stress_glut_write(energy,time)

  type(energy_type), intent(in) :: energy
  double precision, intent(in) :: time

  write(energy%iout_SG,'(7D24.16)') time, energy%sgp, energy%sg
  
  end subroutine stress_glut_write

end module energy
