! SEM2DPACK version 2.3.4 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://www.seismolab.caltech.edu
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
module energy

  implicit none
  private

  type energy_type ! stores macroscopic quantities derived from calculated fields
    double precision :: E_ep = 0d0 ! cumulative plastic energy
    double precision :: E_k = 0d0  ! kinetic energy
    double precision :: E_el = 0d0  ! elastic energy
    integer :: iout_E = 0          ! output unit for energies
  end type energy_type

  public :: energy_type, energy_init, energy_compute, energy_write
  
contains

!=======================================================================
  subroutine energy_init(energy)

  use stdio, only : IO_new_unit

  type(energy_type), intent(out) :: energy

  energy%iout_E = IO_new_unit()
  open(energy%iout_E,file='energy_sem2d.tab',status='replace')

  end subroutine energy_init

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

end module energy
