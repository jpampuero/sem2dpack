! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module solver

! SOLVER: Newmark solver for elasto-dynamic equation
!         M.a = - K.d + F

  use problem_class
  use mat_elastic, only : ELAST_KD,MAT_isElastic
  use mat_kelvin_voigt, only : MAT_KV_add_etav,MAT_isKelvinVoigt
  use mat_damage, only : MAT_DMG_f,MAT_isDamage
  use sources, only : SO_add
  use bc_gen , only : BC_set
  use fields_class, only : FIELD_get_elem, FIELD_add_elem


  implicit none
  private

  public :: solve

contains

!=====================================================================

subroutine solve(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%time%kind)
    case ('leapfrog')
      call solve_leapfrog(pb)
    case ('newmark')
      call solve_Newmark_alpha(pb)
    case default
      call solve_symplectic(pb)
  end select

end subroutine solve


!=====================================================================
! SOLVE: advance ONE time step
!        using single predictor-corrector Newmark-alpha (explicit)
!        in acceleration form
subroutine solve_Newmark_alpha(pb)

  type(problem_type), intent(inout) :: pb
  
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: &
    dloc, floc
  double precision, dimension(:,:), pointer :: d,v,a,f,d_alpha,v_alpha
  double precision :: t_alpha,dt,alpha,beta,gamma
  integer :: e

  d => pb%fields%displ
  v => pb%fields%veloc
  a => pb%fields%accel
  f => a
  d_alpha => pb%fields%displ_alpha
  v_alpha => pb%fields%veloc_alpha

  dt = pb%time%dt
  alpha = pb%time%alpha
  beta = pb%time%beta
  gamma = pb%time%gamma
  
 !-- predictors
  d_alpha = d
  v_alpha = v
  d = d + dt*v + (0.5d0-1d0*beta)*dt*dt *a
  v = v + (1d0-gamma)*dt *a
  d_alpha = alpha*d + (1.d0-alpha)*d_alpha
  v_alpha = alpha*v + (1.d0-alpha)*v_alpha

  f = 0d0
  do e = 1,pb%grid%nelem

    dloc = FIELD_get_elem(d_alpha, pb%grid%ibool(:,:,e))

    if (MAT_isKelvinVoigt(pb%matpro(e))) &
      call MAT_KV_add_etav( dloc, FIELD_get_elem(v_alpha,pb%grid%ibool(:,:,e)) &
                      , pb%matwrk(e) )

   !-- compute - K*D(n+alpha) 
    if (MAT_isElastic(pb%matpro(e))) call ELAST_KD( floc, dloc, pb%matwrk(e) )

   ! NOTE: f=>a, the result is actually stored in a
    call FIELD_add_elem(floc,f,pb%grid%ibool(:,:,e))

  enddo

 !-- Sources
  t_alpha = pb%time%time +(alpha-1.d0)*dt
  call SO_add(pb%src, t_alpha, f)

 !-- Apply boundary conditions
  call BC_set(pb%bc,t_alpha,pb%fields,f)

 ! NOTE: if a source is an incident wave, it is not added during
 !       "call SO_add" but during "call BC_set"
 !       Incident waves must be used with absorbing boundaries

 !-- Divide by mass matrix
  a = f*pb%rmass

 !--  Corrector 
  v = v + gamma*dt*a
  d = d + beta*dt*dt*a
  
end subroutine solve_Newmark_alpha


!=====================================================================
!
! Second-order central difference (leap-frog)
!
! M*(v[n+1/2] - v[n-1/2])/dt = - K*d[n] - C*v[n] + B*t[n] + F(n)
! d[n+1] = d[n] + dt*v[n+1/2]
!
! also: a[n] = (v[n+1/2] - v[n-1/2])/dt

subroutine solve_leapfrog(pb)

  type(problem_type), intent(inout) :: pb
  
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: &
    dloc, floc
  double precision, dimension(:,:), pointer :: d,v_mid,a,f
  integer :: e

  d => pb%fields%displ
  v_mid => pb%fields%veloc
  a => pb%fields%accel
  f => a

  d = d + pb%time%dt * v_mid
  f = 0d0
  do e = 1,pb%grid%nelem

    dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))

    if (MAT_isKelvinVoigt(pb%matpro(e))) &
      call MAT_KV_add_etav( dloc, FIELD_get_elem(v_mid,pb%grid%ibool(:,:,e)) &
                      , pb%matwrk(e) )

    if (MAT_isElastic(pb%matpro(e))) call ELAST_KD(floc,dloc,pb%matwrk(e))

    if (MAT_isDamage(pb%matpro(e))) &
      call MAT_DMG_f(floc,dloc,pb%matwrk(e),pb%grid%ngll,pb%time%dt)

    call FIELD_add_elem(floc,f,pb%grid%ibool(:,:,e))

  enddo

  call SO_add(pb%src, pb%time%time, f)
  call BC_set(pb%bc,pb%time%time,pb%fields,f)

  a = pb%rmass * f
  v_mid = v_mid + pb%time%dt * a 
  
end subroutine solve_leapfrog


!=====================================================================
!
! Symplectic schemes
! WARNING: no boundary conditions implemented yet
!
subroutine solve_symplectic(pb)

  type(problem_type), intent(inout) :: pb
  
  double precision, dimension(pb%grid%ngll,pb%grid%ngll,pb%fields%ndof) :: &
    dloc, floc
  double precision, dimension(:,:), pointer :: d,v,a,f
  double precision, dimension(:), pointer :: coa,cob
  double precision :: dt,t
  integer :: k,e

  d => pb%fields%displ
  v => pb%fields%veloc
  a => pb%fields%accel
  f => pb%fields%accel

  coa => pb%time%a
  cob => pb%time%b
  dt = pb%time%dt
  t = pb%time%time - dt
  
  do k = 1,pb%time%nstages
    d = d + dt*coa(k) * v
    f = 0d0
    do e = 1,pb%grid%nelem
      dloc = FIELD_get_elem(d,pb%grid%ibool(:,:,e))
      if (MAT_isElastic(pb%matpro(e))) call ELAST_KD( floc, dloc, pb%matwrk(e) )
      call FIELD_add_elem(floc,f,pb%grid%ibool(:,:,e))
    enddo
    t = t + dt*coa(k)
    call SO_add(pb%src, t, f)
    a = pb%rmass * f
    v = v + dt*cob(k) * a 
  enddo
  d = d + dt*coa(pb%time%nstages+1) * v
  
  
end subroutine solve_symplectic


end module solver
