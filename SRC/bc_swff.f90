! SEM2DPACK version 2.2.12beta -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bc_swff
! slip weakening friction fault

  use bnd_grid, only: bnd_grid_type
  use distribution_cd

  implicit none
  private

  type bc_swff_input_type
    type(cd_type) :: Dc,MuS,MuD,T,N,Sxx,Sxy,Sxz,Syz,Szz
  end type bc_swff_input_type

  type bc_swff_type
    private
    logical :: UpdateFrictionFirst
    double precision :: CoefA2V,CoefA2D,OutputT1,OutputDt
    integer :: OutputNextTime,OutputNdt,OutputUnit,OutputIFir &
              ,OutputILas,OutputILag,OutputNpoin
    double precision, pointer :: n1(:,:)
    double precision, dimension(:), pointer:: &
      B,invM1=>null(),invM2=>null() ,Z,Tt0,Tn0 &
     ,MU,Tt,Tn,Dc,MuStatic,MuDynamic
    type(bnd_grid_type), pointer :: bc1 => null(), bc2 => null()
    type(bc_swff_input_type) :: input
  end type bc_swff_type

  public :: BC_SWFF_type, BC_SWFF_read, BC_SWFF_init, BC_SWFF_set, BC_SWFF_write

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_SWFFLT
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Slip weakening friction fault
! SYNTAX : &BC_SWFFLT Dc | DcHet, MuS | MuSHet , MuD | MuDHet, 
!                     Tn | TnHet, Tt | TtHet,
!                     Sxx | SxxHet, Sxy | SxyHet, Sxz | SxzHet, 
!		      Syz | SyzHet, Szz | SzzHet
!                     FirstOutput, DtOutput, IxOut /
!             followed eventually by distribution input blocks &DIST_XXX
!             for Dc,MuS,MuD,Tn and/or Tt (the order is important)
!
! NOTE: for better results, use dynamic faults with the leapfrog time scheme
!       and with a layer of damping material (Kelvin-Voigt) near the fault.
!
! Friction law:
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
!
! Initial stress, can be a superposition of tractions and background stress:
! ARG: Tn       [dble] [-100d6] Normal traction (positive = tensile)
! ARG: Tt       [dble] [55d6] Tangential traction (positive antiplane: y>0)
! ARG: Sxx      [dble] [0d0] sigma_xx
! ARG: Sxy      [dble] [0d0] sigma_xy
! ARG: Sxz      [dble] [0d0] sigma_xz
! ARG: Syz      [dble] [0d0] sigma_yz
! ARG: Szz      [dble] [0d0] sigma_zz
!
! NOTE: arguments with the suffix "Het" are used to give
!       friction and initial stress parameters non uniform values.
!	For instance, DcHet='GAUSSIAN' followed by a DIST_GAUSSIAN block
!	sets a gaussian distribution of Dc.
!	Several heterogeneous distributions are available, 
!       See DIST_XXX for their syntax.
!
! For outputs in FltXX_sem2d.dat:
! ARG: DtOutput [dble] [0.d0] Time lag between outputs (in seconds)
!               Default resets DtOutput = global timestep
! ARG: FirstOutput [dble] [0.d0] Start output at this time
! ARG: IxOut    [int(3)] [(1,huge,1)] First node, last node and stride
!               Default resets Ixout(2)= last point
!
! NOTE: DtOutput is internally adjusted to the nearest multiple 
!       of the global timestep
!
! END INPUT BLOCK

  subroutine BC_SWFF_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_swff_type), intent(out) :: bc
  integer, intent(in) :: iin
  double precision :: Dc,MuS,MuD,Tt,Tn,FirstOutput,DtOutput &
                     ,Sxx,Sxy,Sxz,Syz,Szz
  character(20) :: DcHet,MuSHet,MuDHet,TtHet,TnHet &
                  ,SxxHet,SxyHet,SxzHet,SyzHet,SzzHet
  integer :: i,IxOut(3)

  NAMELIST / BC_SWFFLT /  Dc,MuS,MuD,Tt,Tn &
                         ,DcHet,MuSHet,MuDHet,TtHet,TnHet &
                         ,Sxx,Sxy,Sxz,Syz,Szz &
                         ,SxxHet,SxyHet,SxzHet,SyzHet,SzzHet &
                         ,FirstOutput,DtOutput,IxOut

  Dc = 0.5d0
  MuS = 0.6d0
  MuD = 0.5d0
  Tt = 55d6
  Tn = -100.d6
  Sxx = 0d0
  Sxy = 0d0
  Sxz = 0d0
  Syz = 0d0
  Szz = 0d0

  DcHet = ''
  MuSHet = ''
  MuDHet = ''
  TtHet = ''
  TnHet = ''
  SxxHet = ''
  SxyHet = ''
  SxzHet = ''
  SyzHet = ''
  SzzHet = ''

  FirstOutput = 0.d0
  DtOutput = 0.d0
  IxOut(1) = 1
  IxOut(2) = huge(i)
  IxOut(3) = 1

  read(iin,BC_SWFFLT,END=100)
  
  bc%OutputT1 = FirstOutput
  bc%OutputDt = DtOutput
  bc%OutputIFir = IxOut(1)
  bc%OutputILas = IxOut(2) 
  bc%OutputILag = IxOut(3) 

  call DIST_CD_Read(bc%input%Dc,Dc,DcHet,iin,DcHet)
  call DIST_CD_Read(bc%input%MuS,MuS,MuSHet,iin,MuSHet)
  call DIST_CD_Read(bc%input%MuD,MuD,MuDHet,iin,MuDHet)
  call DIST_CD_Read(bc%input%N,Tn,TnHet,iin,TnHet)
  call DIST_CD_Read(bc%input%T,Tt,TtHet,iin,TtHet)
  call DIST_CD_Read(bc%input%Sxx,Sxx,SxxHet,iin,SxxHet)
  call DIST_CD_Read(bc%input%Sxy,Sxy,SxyHet,iin,SxyHet)
  call DIST_CD_Read(bc%input%Sxz,Sxz,SxzHet,iin,SxzHet)
  call DIST_CD_Read(bc%input%Syz,Syz,SyzHet,iin,SyzHet)
  call DIST_CD_Read(bc%input%Szz,Szz,SzzHet,iin,SzzHet)

  if (echo_input) then
    write(iout,200) DcHet,MuSHet,MuDHet,TnHet,TtHet &
                   ,SxxHet,SxyHet,SxzHet,SyzHet,SzzHet &
                   ,FirstOutput,DtOutput,IxOut
    write(iout,*) 'NOTE: FirstOutput and DtOutput may be modified later'
    if (DtOutput==0d0) &
     write(iout,*) 'NOTE: DtOutput will be reset to the global timestep'
    write(iout,*) '      Final values can be found in FltXX_sem2d.hdr'
    if (IxOut(2)==huge(i)) &
     write(iout,*) 'NOTE: Last output node, IxOut(2), will be reset later'
  endif

  return
  100 call IO_abort('BC_SWFF_read: BC_SWFFLT input block not found')
  200 format(5x,'Type   = Slip weakening fault', &
            /5x,'  Critical slip         = ',A,&
            /5x,'  Static friction       = ',A,&
            /5x,'  Dynamic friction      = ',A,&
            /5x,'  Initial traction T_n  = ',A,&
            /5x,'  Initial traction T_t  = ',A,&
            /5x,'  Initial stress S_xx   = ',A,&
            /5x,'  Initial stress S_xy   = ',A,&
            /5x,'  Initial stress S_xz   = ',A,&
            /5x,'  Initial stress S_yz   = ',A,&
            /5x,'  Initial stress S_zz   = ',A,&
            /5x,'  First output time     = ',EN13.3,&
            /5x,'  Output time stride    = ',EN13.3,&
            /5x,'  First output node     = ',I0,&
            /5x,'  Last output node      = ',I0,&
            /5x,'  Output node stride    = ',I0)

  end subroutine BC_SWFF_read


!=====================================================================
!
  subroutine BC_SWFF_init(bc,tags,grid,M,time,v,perio)
  
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort,IO_new_unit
  use time_evol, only: timescheme_type
  use constants, only : TINY_XABS
  use bc_periodic, only : bc_periodic_type,BC_PERIO_intersects

  type(bc_swff_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: M(:)
  integer, intent(in) :: tags(2)
  type(timescheme_type), intent(in) :: time
  double precision, intent(inout) :: v(:,:)
  type(bc_periodic_type), pointer :: perio

  double precision, allocatable :: FltCoord(:,:)
  double precision, dimension(:), pointer :: Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz,nx,nz
  integer :: i,HdrUnit,npoin,NSAMP,NDAT
  character(15) :: OutputFileName,HdrFileName
  logical :: two_sides
  
! if fault is on a domain boundary and symmetry is required, 
! bc2 and invM2 are not defined
  two_sides = tags(2)>0

! bc1 --> grid%bounds(i) boundary corresponding to TAG
  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%bc1 )
  if ( two_sides ) then
    call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%bc2 )

    if (bc%bc1%nelem/=bc%bc2%nelem) &
     call IO_abort('bc_swff_init: number of boundary elements do not match')

    if (bc%bc1%npoin/=bc%bc2%npoin) &
     call IO_abort('bc_swff_init: number of nodes on boundaries do not match')
  endif
  npoin = bc%bc1%npoin

  allocate(FltCoord(2,npoin))
  FltCoord = grid%coord(:,bc%bc1%node)
  if ( two_sides ) then
    if ( any(abs(FltCoord(1,:)-grid%coord(1,bc%bc2%node))>TINY_XABS) &
      .OR.any(abs(FltCoord(2,:)-grid%coord(2,bc%bc2%node))>TINY_XABS) )&
      call IO_abort('bc_swff_init: coordinates on boundaries do not match properly')
  endif

! NOTE: the mesh being conformal, the weights B=GLL_weights*jac1D are equal on both
!       sides of the fault. 
  allocate( bc%n1(npoin,2) )
  allocate( bc%B(npoin) ) ! assembled[ GLL_weights * jac1D ]
  call BC_get_normal_and_weights(bc%bc1,grid,bc%n1,bc%B, BC_PERIO_intersects(bc%bc1,perio) )

! Update coefficients of Newmark scheme
! vnew = vpredictor +gamma*dt*anew   See solver.f90 
! dnew = dpredictor +beta*dt^2*anew   See solver.f90
  bc%CoefA2V = time%gamma * time%dt
  bc%CoefA2D = time%beta * time%dt**2
  bc%UpdateFrictionFirst = time%kind == 'leapfrog'
!  bc%UpdateFrictionFirst = .true. !WARNING: testing

! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1
  allocate(bc%invM1(npoin))
  bc%invM1 = 1d0/M(bc%bc1%node)
  if ( two_sides ) then
    allocate(bc%invM2(npoin))
    bc%invM2 = 1d0/M(bc%bc2%node)
  endif

! Fault impedance, Z in :  Trac=T_Stick-Z*dV
!   Z = 1/( B1/M1 + B2/M2 )
! T_Stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity) 
! NOTE: same Bi on both sides, see note above
  allocate(bc%Z(npoin))
  if ( two_sides ) then
    bc%Z = 1.d0/(bc%CoefA2V * bc%B *( bc%invM1 + bc%invM2 ))
  else
    bc%Z = 0.5d0/(bc%CoefA2V * bc%B * bc%invM1 )
  endif

! Friction parameters
  call DIST_CD_Init(bc%input%Dc,FltCoord,bc%Dc)
  call DIST_CD_Init(bc%input%MuS,FltCoord,bc%MuStatic)
  call DIST_CD_Init(bc%input%MuD,FltCoord,bc%MuDynamic)

! Initial stress
  call DIST_CD_Init(bc%input%T,FltCoord,bc%Tt0)
  call DIST_CD_Init(bc%input%N,FltCoord,bc%Tn0)
  call DIST_CD_Init(bc%input%Sxx,FltCoord,Sxx)
  call DIST_CD_Init(bc%input%Sxy,FltCoord,Sxy)
  call DIST_CD_Init(bc%input%Sxz,FltCoord,Sxz)
  call DIST_CD_Init(bc%input%Syz,FltCoord,Syz)
  call DIST_CD_Init(bc%input%Szz,FltCoord,Szz)
  nx => bc%n1(:,1)
  nz => bc%n1(:,2)
  allocate(Tx(npoin))
  allocate(Ty(npoin))
  allocate(Tz(npoin))
  Tx = Sxx*nx + Sxz*nz
  Ty = Sxy*nx + Syz*nz
  Tz = Sxz*nx + Szz*nz
  bc%Tn0 = bc%Tn0 + Tx*nx + Tz*nz
  if (size(v,2)==1) then ! SH
    bc%Tt0 = bc%Tt0 + Ty
  else ! P-SV
    bc%Tt0 = bc%Tt0 + Tx*nz - Tz*nx
  endif 
  deallocate(Sxx,Sxy,Sxz,Syz,Szz,Tx,Ty,Tz)

! Initial friction
  allocate(bc%MU(npoin))
  bc%MU = bc%MuStatic
 
! BETA :
! Initial velocities
! To guarantee second order when initial stress drop is abrupt, 
! the initial velocity must be set analytically:
!   impedance*slip_rate = max( initial_stress - initial_strength, 0)
! where impedance = 0.5*rho*cs
!       slip rate = v_2 - v_1
!                 = -2*v_1 if symmetry
!  if (size(v,2)==1) then ! SH
!    v(bc%bc1%node,1) = -sign( max( abs(bc%Tt0)+bc%MU*bc%Tn0, 0d0), bc%Tt0 ) &
!                         / (2670d0*3464d0)
!    if (two_sides) v(bc%bc2%node,1) = -v(bc%bc1%node,1)
!  else ! P-SV
!! WARNING: not implemented yet
!  endif


! Output arrays
  bc%OutputIFir = max(bc%OutputIFir, 1)
  bc%OutputILas = min(bc%OutputILas, npoin)
  bc%OutputNpoin = (bc%OutputILas-bc%OutputIFir)/bc%OutputILag +1
  allocate(bc%Tt(bc%OutputNpoin))
  allocate(bc%Tn(bc%OutputNpoin))
  bc%Tt=bc%Tt0(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  bc%Tn=bc%Tn0(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  if (size(bc%Tn0(bc%OutputIFir:bc%OutputILas:bc%OutputILag))/=bc%OutputNpoin) &
    stop 'couille avec npoin'

!-- Open output files -----
! Set output parameters

! NOTE: file names limited to tags(1)<100
  write(OutputFileName,'("Flt",I2.2,"_sem2d.dat")') tags(1)
  bc%OutputUnit = IO_new_unit()
  open(bc%OutputUnit,file=OutputFileName,status='replace',form='unformatted')

 ! adjust output timestep to the nearest multiple of dt:
  if (bc%OutputDt==0d0) then
    bc%OutputNdt=1
  else
    bc%OutputNdt = nint(bc%OutputDt/time%dt)
  endif
  bc%OutputDt = time%dt * bc%OutputNdt
  bc%OutputNextTime = nint(bc%OutputT1/time%dt)

 ! HEADER FILE FORMAT:
 !      Name            FltXX_sem2d.hdr where XX=tags(1) of the BC_SWFFLT
 !                      input block, the tag of the first side of the fault.
 !      Line 1          Name of size/other variables
 !      Line 2          Value of size/other variables
 !      Line 3          Name of data fields, separated by ":"
 !      Line 4          Name of coordinate axis
 !      Line 5:end      Two column table of nodes coordinates
 !
 ! OUTPUT FILE FORMAT:
 !      At each time (lag DELT, total NSAMP), NDAT data lines with
 !      NPTS columns, one per fault node, are written.
 !      Stress fields are relative to initial values and
 !      the first NDAT lines contain the initial conditions.
 !
 ! NOTE: when the binary file is written one line at a time
 !      a record tag is also written at the beginning and end of line,
 !      => number of columns in output file = nb of data columns +2
 !
  write(HdrFileName,'("Flt",I2.2,"_sem2d")') tags(1)
  NDAT = 5
  NSAMP = (time%nt -bc%OutputNextTime)/bc%OutputNdt +1
  HdrUnit = IO_new_unit()

  open(HdrUnit,file=trim(HdrFileName)//'.hdr',status='replace')
  write(HdrUnit,*) 'NPTS NDAT NSAMP DELT'
  write(HdrUnit,*) bc%OutputNpoin,NDAT,NSAMP,bc%OutputDt
                  
  write(HdrUnit,*) 'Slip:Slip Rate:Shear Stress:Normal Stress:Friction'
  write(HdrUnit,*) 'XPTS ZPTS'
  do i=bc%OutputIFir,bc%OutputILas,bc%OutputILag
    write(HdrUnit,*) FltCoord(:,i)
  enddo
  close(HdrUnit)
  deallocate(FltCoord)

! Sample SU script subset|b2a
 ! NOTE: how to convert to a human-readable ASCII table
 !      the Slip at time 10*DELT (assuming first output at t=0), 
 !      using Seismic Unix's subset and b2a :
 !              subset <Flt05_sem2d.hdr n1=(NPTS+2) n2=NDAT n3=NSAMP \
 !              if1s=1 n1s=NPTS ix2s=0 ix3s=10 \
 !              | b2a n1=NPTS > slip_t10.tab
 !
  open(HdrUnit,file=trim(HdrFileName)//'.csh',status='replace')
  write(HdrUnit,*) '#!/bin/csh -f'
  write(HdrUnit,*) '#USAGE: '//trim(HdrFileName)//'.csh node field'
  write(HdrUnit,*) '@ ifs = $2 - 1'
  write(HdrUnit,'(A,I0,A,I0,A,I0,A)') 'subset <'//trim(HdrFileName)//'.dat n1=',bc%OutputNpoin+2, &
    ' n2=',NDAT,' n3=',NSAMP,' ix1s=$1 ix2s=$ifs | b2a n1=1 > '//trim(HdrFileName)//'.$1.$ifs.tab'
  close(HdrUnit)

  end subroutine BC_SWFF_init


!=====================================================================
! Computes the boundary term B*T with T satisfying a friction law
!
! Convention:   T = sigma*n1
!               D = d2-d1
! => T and D have same sign
!
! NOTE: this is an EXPLICIT scheme FOR FRICTION, 
!       we use the displacement at the PREVIOUS timestep 
!       to compute the friction coefficient
!
! NOTE: possible periodicity does not need to be enforced at this level
!       because it is assumed that MxA and D are already periodic
!
  subroutine BC_SWFF_set(bc,MxA,V,D)

  type(bc_swff_type) :: bc
  double precision, intent(in) :: V(:,:),D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  if (size(V,2) == 1) then
    call BC_SWFF_set_SH(bc,MxA(:,1),V(:,1),D(:,1))
  else
    call BC_SWFF_set_PSV(bc,MxA,V,D)
  endif

  end subroutine BC_SWFF_set

!---------------------------------------------------------------------

  subroutine BC_SWFF_set_SH(bc,MxA,V,D)

  type(bc_swff_type) :: bc
  double precision, intent(in) :: V(:),D(:)
  double precision, intent(inout) :: MxA(:)

  double precision, dimension(bc%bc1%npoin) :: T,dD
  integer, pointer :: i1(:),i2(:)
  logical :: two_sides

  two_sides = associated(bc%bc2)
  i1 => bc%bc1%node

! A test traction is computed as if the fault was sticking at n+1 (no slip rate)
!           = Z * SlipRateFree 
!	    = Z * (SlipRate_predicted + coef*SlipAccelFree)
! where: 
!   Z = fault impedance matrix (diagonal)
!   SlipRateFree, SlipAccelFree = slip velocity and acceleration 
!     				  as if the fault was stress free
!   coef = depends on the time integration scheme
!
  if (two_sides) then
    i2 => bc%bc2%node
    T = bc%Z* ( V(i2)-V(i1) + bc%CoefA2V*( bc%invM2*MxA(i2)-bc%invM1*MxA(i1) ))
    dD = D(i2) - D(i1)
  else
    T  = -2d0*bc%Z* ( V(i1) + bc%CoefA2V * bc%invM1*MxA(i1) )
    dD = -2d0*D(i1)
  endif

! Update slip weakening (only for leapfrog scheme)
  if (bc%UpdateFrictionFirst) bc%MU = friction(dD,bc%MuStatic,bc%MuDynamic,bc%Dc)

! Coulomb friction, applied to the absolute stress
  T = T + bc%Tt0
  T = sign( min(abs(T), -bc%Tn0*bc%MU), T)
  T = T - bc%Tt0

! Add boundary term B*T to M*a
  MxA(i1) = MxA(i1) + bc%B*T
  if (two_sides) MxA(i2) = MxA(i2) - bc%B*T

! Store relative tractions for output
  bc%Tt = T(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  bc%Tn = 0d0

! Newmark scheme: update friction for next time step
! Introduces a delay of dt, but keeps the friction solver fully explicit
  if (.not. bc%UpdateFrictionFirst) then
    if (two_sides) then
      dD = dD + bc%CoefA2D *( bc%invM2*MxA(i2)-bc%invM1*MxA(i1) )
    else
      dD = dD -2d0*bc%CoefA2D * bc%invM1*MxA(i1)
    endif
    bc%MU = friction(dD, bc%MuStatic, bc%MuDynamic, bc%Dc)
  endif

  end subroutine BC_SWFF_set_SH


!---------------------------------------------------------------------
! WARNING: the opening conditions are incomplete
!          there can be a problem if the fault opens and closes again

  subroutine BC_SWFF_set_PSV(bc,MxA,V,D)

  type(bc_swff_type) :: bc
  double precision, intent(in) :: V(:,:),D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%bc1%npoin) :: Slip,Tt,Tn
  double precision, dimension(bc%bc1%npoin,2) :: BxT,T_Stick,dD,dA

  integer, pointer :: i1(:),i2(:)
  logical :: two_sides

  two_sides = associated(bc%bc2)

! Test traction
  i1 => bc%bc1%node
  if (two_sides) then
    i2 => bc%bc2%node
    dD = D(i2,:) - D(i1,:)
    dA(:,1) = bc%invM2*MxA(i2,1) - bc%invM1*MxA(i1,1)
    dA(:,2) = bc%invM2*MxA(i2,2) - bc%invM1*MxA(i1,2)
    T_Stick = V(i2,:)-V(i1,:) + bc%CoefA2V * dA
  else
    dD = - 2d0*D(i1,:) ! used only to compute tangential discontinuity (slip)
    dA(:,1) = - 2d0*bc%invM1*MxA(i1,1)
    dA(:,2) = - 2d0*bc%invM1*MxA(i1,2)
    T_Stick = -2d0*V(i1,:) + bc%CoefA2V * dA
  endif
  T_Stick(:,1) = bc%Z * T_Stick(:,1)
  T_Stick(:,2) = bc%Z * T_Stick(:,2)

 ! Update slip weakening:
!WARNING: there should be no update if there is opening, 
!         if (Norm>0) return
!         should be implemented here if required (unlikely)
!
  Slip = bc%n1(:,2)*dD(:,1) - bc%n1(:,1)*dD(:,2)
  bc%MU = friction(Slip, bc%MuStatic, bc%MuDynamic, bc%Dc)

! Rotate tractions to fault frame (Tt,Tn)
  Tt = bc%n1(:,2)*T_Stick(:,1) - bc%n1(:,1)*T_Stick(:,2)
  if (two_sides) then
    Tn = bc%n1(:,1)*T_Stick(:,1) + bc%n1(:,2)*T_Stick(:,2)
  else
    Tn = 0d0    ! if symmetry: costant normal stress
  endif

! Add initial stress
  Tt = Tt +bc%Tt0
  Tn = Tn +bc%Tn0

 !Opening-->free stress
  Tn = min(Tn,0.d0) ! negative normal stress is compressive
 !Coulomb friction
  Tt = sign( min(abs(Tt),-Tn*bc%MU), Tt)

! Subtract initial stress
  Tt = Tt -bc%Tt0
  Tn = Tn -bc%Tn0

! Rotate tractions back to (x,z) frame and apply boundary weights
  BxT(:,1) = bc%B*(  bc%n1(:,2)*Tt + bc%n1(:,1)*Tn )
  BxT(:,2) = bc%B*( -bc%n1(:,1)*Tt + bc%n1(:,2)*Tn )

! Add boundary term B*T to M*a
  MxA(i1,:) = MxA(i1,:) + BxT
  if (two_sides) MxA(i2,:) = MxA(i2,:) - BxT

! Store tractions for output
  bc%Tt = Tt(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  bc%Tn = Tn(bc%OutputIFir:bc%OutputILas:bc%OutputILag)

  end subroutine BC_SWFF_set_PSV



!=====================================================================
! OUTPUT FORMAT: at each time, 5 data lines, one column per fault node
! 1: Slip
! 2: Slip rate
! 3: Shear stress 
! 4: Normal stress 
! 5: Friction coefficient
!
  subroutine BC_SWFF_write(bc,displ,veloc,itime)

  type(bc_swff_type) :: bc
  double precision, intent(in) :: Veloc(:,:),Displ(:,:)
  integer, intent(in) :: itime

  double precision, pointer :: n1(:,:)
  integer, pointer :: i1(:),i2(:)

  if (itime<bc%OutputNextTime) return
  bc%OutputNextTime = bc%OutputNextTime + bc%OutputNdt

  if (associated(bc%bc2)) then

    i1 => bc%bc1%node(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
    i2 => bc%bc2%node(bc%OutputIFir:bc%OutputILas:bc%OutputILag)

    if ( size(displ,2)==1 ) then
      write(bc%OutputUnit) real(displ(i2,1)-displ(i1,1))
      write(bc%OutputUnit) real(veloc(i2,1)-veloc(i1,1))
    else
      n1 => bc%n1(bc%OutputIFir:bc%OutputILas:bc%OutputILag,:)
      write(bc%OutputUnit) real( n1(:,2)*(displ(i2,1)-displ(i1,1)) &
                               - n1(:,1)*(displ(i2,2)-displ(i1,2)) )
      write(bc%OutputUnit) real( n1(:,2)*(veloc(i2,1)-veloc(i1,1)) &
                               - n1(:,1)*(veloc(i2,2)-veloc(i1,2)) )
    endif

  else

    i1 => bc%bc1%node(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
    if ( size(displ,2)==1 ) then
      write(bc%OutputUnit) -real(2d0*displ(i1,1))
      write(bc%OutputUnit) -real(2d0*veloc(i1,1))
    else
      n1 => bc%n1(bc%OutputIFir:bc%OutputILas:bc%OutputILag,:)
      write(bc%OutputUnit) -2.0*real( n1(:,2)*displ(i1,1) &
                                - n1(:,1)*displ(i1,2) )
      write(bc%OutputUnit) -2.0*real( n1(:,2)*veloc(i1,1) &
                                - n1(:,1)*veloc(i1,2) )
    endif
  endif

  write(bc%OutputUnit) real( bc%Tt )
  write(bc%OutputUnit) real( bc%Tn )
  write(bc%OutputUnit) real( bc%MU(bc%OutputIFir:bc%OutputILas:bc%OutputILag) )

  end subroutine BC_SWFF_write



!=====================================================================
 !-- linear slip weakening:
  function friction(s,mus,mud,dc) result(mu)

  double precision, dimension(:), intent(in) :: s,mus,mud,dc
  double precision :: mu(size(s))

  mu = mus -(mus-mud)/dc *abs(s)
  mu = max( mu, mud)

  end function friction

!---------------------------------------------------------------------
 !-- exponential slip weakening:
  function friction_exponential(s,mus,mud,dc) result(mu)

  double precision, dimension(:), intent(in) :: s,mus,mud,dc
  double precision :: mu(size(s))

  mu = mus -(mus-mud)*exp(-abs(s)/dc)

  end function friction_exponential

  end module bc_swff
