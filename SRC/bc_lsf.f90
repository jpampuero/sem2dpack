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
module bc_lsf
! Linear-slip fault

  use bnd_grid, only: bnd_grid_type

  implicit none
  private

  type bc_lsf_type
! kind = 1 only tangent displacement discontinuity
!        2 only normal ...
!        3 both
    integer :: kind
!    double precision, pointer :: kt(:),kn(:)
    double precision :: kt,kn
    double precision, pointer :: normal(:,:),weights(:),MassRatio1(:),MassRatio2(:)
    type(bnd_grid_type), pointer :: bc1,bc2
  end type bc_lsf_type

  public :: BC_LSF_type, BC_LSF_read, BC_LSF_init, BC_LSF_set

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_LSF 
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Linear slip fault, a displacement discontinuity interface
!          where stress and disp.discont. are linearly related
! SYNTAX : &BC_LSF Ktang,Knorm /
!
! ARG: Ktang    [dble] [Inf] Tangential stiffness
! ARG: Ctang    [dble] [0d0] Tangential compliance
! ARG: Knorm    [dble] [Inf] Normal stiffness
! ARG: Cnorm    [dble] [0d0] Normal compliance
!
! NOTE: for each component you can set K _or_ C, but _not_both_
!       
! NOTE: if one of the C=0d0 or K=Inf (the default) then
!       no displacement discontinuity is allowed for that component
!       (transparent),
!       if K=0d0 the fault is a free stress boundary for that component
!       In summary the fault can behave as:
!              -1       transparent T&N (Tangent and Normal)
!               0       stress free T&N
!               1       linear-slip/free T, transparent N 
!               2       transparent T, linear-slip/free N
!               3       linear-slip/free T&N
!
! END INPUT BLOCK

  subroutine BC_LSF_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_lsf_type), intent(out) :: bc
  integer, intent(in) :: iin
  double precision :: Ktang,Knorm,InitDble,Ctang,Cnorm

  NAMELIST / BC_LSF /  Ktang,Knorm,Ctang,Cnorm

  Ctang = 0d0
  Cnorm = 0d0
  InitDble = huge(InitDble) ! an implausible initial value
  Ktang = InitDble          ! the default is : welded contact
  Knorm = InitDble

  read(iin,BC_LSF,END=100)
  
  if (Ctang/=0d0) Ktang= 1d0/Ctang ! Compliance input overwrites
  if (Cnorm/=0d0) Knorm= 1d0/Cnorm ! stiffness input

  bc%kind = -1  ! default: welded contact (continuous displ and sress)
  if (Ktang/=InitDble) then
    bc%kind = 1
    bc%kt   = Ktang
  endif
  if (Knorm/=InitDble) then
    bc%kind = bc%kind+2
    bc%kn   = Knorm
  endif
  if (Knorm==0d0.and.Ktang==0d0) bc%kind=0 ! free boundaries

  if (echo_input) write(iout,200) bc%kind,bc%kt,bc%kn

  return
  100 call IO_abort('BC_LSF_read: BC_LSF input block not found')
  200 format(5x,'Type   = Linear slip fault', &
            /5x,'  Kind                 =',I0,&
            /5x,'  Tangential stiffness =',EN12.3,&
            /5x,'  Normal stiffness     =',EN12.3)

  end subroutine BC_LSF_read


!=====================================================================
!
  subroutine BC_LSF_init(bc,tags,grid,M,perio)
  
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort
  use constants, only : TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_lsf_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: M(:)
  integer, intent(in) :: tags(2)
  type(bc_periodic_type), pointer :: perio

  double precision, pointer :: tmp(:)
  
  if (bc%kind==0) return ! free boundaries

! bc1 --> grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%bc1 )
  call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%bc2 )

  if (bc%bc1%nelem/=bc%bc2%nelem) &
   call IO_abort('bc_lsf_init: number of boundary elements do not match')

  if (bc%bc1%npoin/=bc%bc2%npoin) &
   call IO_abort('bc_lsf_init: number of nodes on boundaries do not match')

  if ( any(abs(grid%coord(1,bc%bc1%node)-grid%coord(1,bc%bc2%node))>TINY_XABS) &
   .OR.any(abs(grid%coord(2,bc%bc1%node)-grid%coord(2,bc%bc2%node))>TINY_XABS) )&
   call IO_abort('bc_lsf_init: coordinates on boundaries do not match properly')

  allocate( bc%normal(bc%bc1%npoin,2) )
  allocate( bc%weights(bc%bc1%npoin) ) ! assembled[ GLL_weights * jac1D ]
  call BC_get_normal_and_weights(bc%bc1,grid,bc%normal,bc%weights, BC_PERIO_intersects(bc%bc1,perio))
! NOTE: the mesh being conformal, the [GLL_weights*jac1D] are equal on both
!       sides of the fault. 

 !The following are needed to compute the traction 
 !as if the fault was stuck (no displ discontinuity)
 !Tstick = ( -K2*d2/M2 + K1*d1/M1 )/( B1/M1 + B2/M2 )
 !Here we store: bc%MassRatioi = 1/Mi /( B1/M1 + B2/M2 )
  allocate(bc%MassRatio1(bc%bc1%npoin))
  allocate(bc%MassRatio2(bc%bc2%npoin))
! NOTE: B1=B2, same GLL_weights*jac1D on both sides, see note above
! NOTE: assuming periodic boundaries have been already initialized 
!       (M and weights respect periodicity)
  tmp => bc%MassRatio2
  tmp = bc%weights*( 1d0/M(bc%bc1%node)+1d0/M(bc%bc2%node) )
  tmp = 1d0/tmp
  bc%MassRatio1 = 1d0/M(bc%bc1%node)*tmp
  bc%MassRatio2 = 1d0/M(bc%bc2%node)*tmp

  end subroutine BC_LSF_init


!=====================================================================
! Computes the boundary term B*T with T=k*D
! Traction = Fault_stiffness * displacement_discontinuity
!
! Convention:   T = sigma*n1
!               D = d2-d1
! => T and D have same sign
!
! NOTE: this is an EXPLICIT scheme, 
!       we use the PREDICTED displacement to compute the slip
!
! NOTE: eventual periodicity does not need to be enforced at this level
!       because it is assumed that MxA and D are already periodic
!
  subroutine BC_LSF_set(bc,MxA,D)

  type(bc_lsf_type), intent(in)    :: bc
  double precision, intent(in) :: D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%bc1%npoin) :: Tt,Tn
  double precision, dimension(bc%bc1%npoin,2) :: BxT,Tstick,dD
  integer, pointer :: i1(:),i2(:)

  if (bc%kind ==0)  return ! free boundaries: do nothing

  i1 => bc%bc1%node
  i2 => bc%bc2%node

! STEP 1: get discontinuity of input fields across the fault
 !Displacement discontinuity = Face2 - Face1
  dD = D(i2,:) - D(i1,:)
 !Traction as if the fault was sticking (no accel discontinuity)
 !Tstick = ( -K2*d2/M2 + K1*d1/M1 )/( B1/M1 + M2/B2 )
 !on entry -K*d is stored in MxA
 !       , 1/Mi /( B1/M1 + M2/B2 ) is stored in bc%MassRatioi
  Tstick(:,1)= bc%MassRatio2*MxA(i2,1) - bc%MassRatio1*MxA(i1,1)
  Tstick(:,2)= bc%MassRatio2*MxA(i2,2) - bc%MassRatio1*MxA(i1,2)

! STEP 2: Get tractions in fault frame (Tt,Tn)
  select case(bc%kind)

  case(-1) ! welded contact in both components
    Tt = bc%normal(:,2)*Tstick(:,1) - bc%normal(:,1)*Tstick(:,2)
    Tn = bc%normal(:,1)*Tstick(:,1) + bc%normal(:,2)*Tstick(:,2)

  case(1)
  ! Linear slip law for tangent component
  ! NOTE: this is an EXPLICIT scheme 
  ! If Newmark: uses the PREDICTED slip to compute the tractions
  ! If leapfrog: second order ok
    Tt = bc%normal(:,2)*dD(:,1) - bc%normal(:,1)*dD(:,2) ! = tangent slip
    Tt = bc%kt * Tt ! = tangent traction
  ! Normal accel is continuous
    Tn = bc%normal(:,1)*Tstick(:,1) + bc%normal(:,2)*Tstick(:,2)

  case(2)
  ! Linear slip law for normal component
    Tn = bc%normal(:,1)*dD(:,1) + bc%normal(:,2)*dD(:,2)
    Tn = bc%kn * Tn
  ! Tangent accel is continuous
    Tt = bc%normal(:,2)*Tstick(:,1) - bc%normal(:,1)*Tstick(:,2)

  case(3)
  ! Rotate slip to fault frame (tangent,normal)
    Tt = bc%normal(:,2)*dD(:,1) - bc%normal(:,1)*dD(:,2)
    Tn = bc%normal(:,1)*dD(:,1) + bc%normal(:,2)*dD(:,2)
  ! Apply linear slip law
    Tt = bc%kt * Tt
    Tn = bc%kn * Tn

  end select

! STEP 3: Add boundary term B*T to M*a
 !Rotate tractions back to original frame and apply boundary weights
  BxT(:,1) = bc%weights*(  bc%normal(:,2)*Tt + bc%normal(:,1)*Tn )
  BxT(:,2) = bc%weights*( -bc%normal(:,1)*Tt + bc%normal(:,2)*Tn )
 !add
  MxA(i1,:) = MxA(i1,:) + BxT
  MxA(i2,:) = MxA(i2,:) - BxT

  end subroutine BC_LSF_set

  end module bc_lsf
