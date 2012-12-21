! SEM2DPACK version 2.2.3 -- A Spectral Element Method tool for 2D wave propagation
!                            and earthquake source dynamics
! 
! Copyright (C) 2003 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! ETH Zurich (Swiss Federal Institute of Technology)
! Institute of Geophysics
! Seismology and Geodynamics
! ETH Hönggerberg (HPP)
! CH-8093 Zürich
! Switzerland
! 
! ampuero@erdw.ethz.ch
! +41 1 633 2197 (office)
! +41 1 633 1065 (fax)
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
module solver

! SOLVER: Newmark solver for elasto-dynamic equation
!         M.a = - K.d + F

  use problem_class

  implicit none
  private

  public :: solve

contains

!=====================================================================
! SOLVE: advance ONE time step
!        using single predictor-corrector Newmark-alpha (explicit)
!        in acceleration form

subroutine solve(pb)

  use elastic, only : ELAST_KD => ELAST_KD2
  use sources, only : SO_add
  use bc_gen , only : BC_set

  type(problem_type), intent(inout) :: pb
  
  double precision, pointer :: tmp(:,:)

  tmp => pb%fields%tmp
  tmp = pb%fields%displ

 !-- predictors
  pb%fields%displ = pb%fields%displ + pb%time%dt*pb%fields%veloc
  if (pb%time%coef(1)/=0d0) pb%fields%displ = pb%fields%displ + pb%time%coef(1)*pb%fields%accel

  pb%fields%veloc = pb%fields%veloc + pb%time%coef(2)*pb%fields%accel

 !-- compute - K*D(i,n+alpha) 
 !   store the result in pb%fields%accel
  tmp = pb%time%alpha*pb%fields%displ + (1.d0-pb%time%alpha)*tmp 
  !print *,'Entering elast_kd'
  call ELAST_KD( pb%elast, pb%grid, tmp, pb%fields%accel)
  !print *,'done'

 !-- Sources
  if (associated(pb%src)) call SO_add(pb%src,pb%time%time+(pb%time%alpha-1.d0)*pb%time%dt)

 !-- Some BCs
  !print *,'Entering BC_set'
  call BC_set(pb%bc,which='PERIOD',field=pb%fields%accel,assemble='yes')
  call BC_set(pb%bc,which='DT0TN0',field=pb%fields%accel)
  call BC_set(pb%bc,which='LISFLT',field=pb%fields%accel,field2=tmp)
  call BC_set(pb%bc,which='SWFFLT',fields=pb%fields)
  !print *,'done'

 !-- Divide by mass matrix
  pb%fields%accel(:,1) =  pb%fields%accel(:,1)*pb%rmass
  pb%fields%accel(:,2) =  pb%fields%accel(:,2)*pb%rmass

 !--  apply absorbing boundary conditions
 ! NOTE: I don't like the implementation of this ABSORB
 !       Will replace it by a classical Clayton-Engquist
 !       so I can use which='ALL' above
  !print *,'Entering BC_set for ABSORB'
  call BC_set(pb%bc,which='ABSORB',src=pb%src,grid=pb%grid,fields=pb%fields &
             ,time=pb%time%time)
  !print *,'done'

 !--  Corrector 
  pb%fields%veloc = pb%fields%veloc + pb%time%coef(3)*pb%fields%accel
  pb%fields%displ = pb%fields%displ + pb%time%coef(4)*pb%fields%accel
  
end subroutine

end module solver
