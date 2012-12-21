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
module bc_swff
! slip weakening friction fault

  use spec_grid, only: bc_topo_type
  use distribution_general, only: distribution_type

  implicit none
  private

  type bc_swff_type
    private
    logical :: periodic
    double precision :: CoefA2V,CoefA2D,OutputT1,OutputDt &
                       ,Dc_const,MuS_const,MuD_const,InitStrT_const,InitStrN_const
    integer :: OutputNextTime,OutputNdt,OutputUnit,OutputIFir,OutputILas,OutputILag,OutputNpoin
    double precision, pointer :: n1(:,:)
    double precision, dimension(:), pointer:: B,RMass1,RMass2,Z,InitStrT,InitStrN &
                                             ,MU,StrT,StrN,Dc,MuStatic,MuDynamic
    type(bc_topo_type), pointer :: bc1,bc2
    type(distribution_type), pointer :: Dc_dist,MuS_dist,MuD_dist,InitStrT_dist,InitStrN_dist
  end type bc_swff_type

  public :: BC_SWFF_type, BC_SWFF_read, BC_SWFF_init, BC_SWFF_set, BC_SWFF_write

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_SWFFLT [boundary condition]
! PURPOSE: Slip weakening friction fault
! SYNTAX : &BC_SWFFLT Dc|DcHet, MuS|MuSHet ,MuD|MuDHet, 
!                     StrN|StrNHet, StrT|StrTHet,
!                     Periodic,
!                     FirstOutput, DtOutput, IxOut /
!             followed eventually by distribution input blocks &DIST_XXX
!             for Dc,MuS,MuD,StrN and/or StrT (the order is important)
!
! Friction law:
! ARG: Dc       [dble] [0.5d0] Critical slip 
! ARG: MuS      [dble] [0.6d0] Static friction coefficient
! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
! ARG: DcHet    [dble] [none] Heterogeneous distribution for critical slip 
! ARG: MuSHet   [dble] [none] Het. dist. for static friction coefficient
! ARG: MuDHet   [dble] [none] Het. dist. for dynamic friction coefficient
!
! Initial stress:
! ARG: StrT     [dble] [5.5d6] Constant shear stress
! ARG: StrN     [dble] [10.d6] Constant normal stress
! ARG: StrTHet  [name] [none] Het. dist. for shear stress
! ARG: StrNHet  [name] [none] Het. dist. for normal stress
!
! ARG: Periodic [log] [F]  Enable if the fault crosses a periodic boundary
!
! For outputs in FltXX_sem2d.dat:
! ARG: DtOutput [dble] [0.d0] Time lag between outputs (in seconds)
!               Default resets DtOutput = global timestep
! ARG: FirstOutput [dble] [0.d0] Start output at this time
! ARG: IxOut    [int(3)] [(1,huge,1)] First node, last node and stride
!               Default resets Ixout(2)= last point
!
! NOTE: several heterogeneous distributions are available, 
!       see DIST_XXX for their syntax
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
  double precision :: Dc,MuS,MuD,StrT,StrN,FirstOutput,DtOutput
  character(10) :: DcHet,MuSHet,MuDHet,StrTHet,StrNHet
  logical :: Periodic
  integer :: i,IxOut(3)

  NAMELIST / BC_SWFFLT /  Dc,MuS,MuD,StrT,StrN &
                         ,DcHet,MuSHet,MuDHet,StrTHet,StrNHet &
                         ,Periodic,FirstOutput,DtOutput,IxOut

  Dc = 0.5d0
  MuS = 0.6d0
  MuD = 0.5d0
  StrT = 5.5d6
  StrN = 10.d6
  DcHet = ''
  MuSHet = ''
  MuDHet = ''
  StrTHet = ''
  StrNHet = ''
  FirstOutput = 0.d0
  DtOutput = 0.d0
  Periodic = .false.
  IxOut(1) = 1
  IxOut(2) = huge(i)
  IxOut(3) = 1

  read(iin,BC_SWFFLT,END=100)
  
  bc%periodic = Periodic
  bc%OutputT1 = FirstOutput
  bc%OutputDt = DtOutput
  bc%OutputIFir = IxOut(1)
  bc%OutputILas = IxOut(2) 
  bc%OutputILag = IxOut(3) 

  call ReadConstOrDist(Dc,DcHet,bc%Dc_const,bc%Dc_dist,iin)
  call ReadConstOrDist(MuS,MuSHet,bc%MuS_const,bc%MuS_dist,iin)
  call ReadConstOrDist(MuD,MuDHet,bc%MuD_const,bc%MuD_dist,iin)
  call ReadConstOrDist(StrN,StrNHet,bc%InitStrN_const,bc%InitStrN_dist,iin)
  call ReadConstOrDist(StrT,StrTHet,bc%InitStrT_const,bc%InitStrT_dist,iin)

  if (echo_input) then
    write(iout,200) DcHet,MuSHet,MuDHet,StrNHet,StrTHet,Periodic &
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
            /5x,'  Initial normal stress = ',A,&
            /5x,'  Initial shear stress  = ',A,&
            /5x,'  Periodic fault        = ',L1,&
            /5x,'  First output time     = ',EN12.3,&
            /5x,'  Output time stride    = ',EN12.3,&
            /5x,'  First output node     = ',I0,&
            /5x,'  Last output node      = ',I0,&
            /5x,'  Output node stride    = ',I0)

  end subroutine BC_SWFF_read

!---------------------------------------------------------------------
  subroutine ReadConstOrDist(Cin,Din,Cout,Dout,iin)

  use distribution_general, only: DIST_read

  double precision, intent(in) :: Cin
  character(10), intent(inout) :: Din
  double precision, intent(out) :: Cout
  type(distribution_type), pointer :: Dout
  integer, intent(in) :: iin

  if (Din/='') then
    allocate(Dout)
    call DIST_read(Dout,Din,iin)
  else
    Cout = Cin
    write(Din,'(F10.3)') Cin
  endif

  end subroutine ReadConstOrDist

!=====================================================================
!
  subroutine BC_SWFF_init(bc,tags,grid,RMass,time)
  
  use spec_grid, only : sem_grid_type,BC_inquire,BC_get_normal_and_weights
  use stdio, only: IO_abort,IO_new_unit
  use time_evol, only: timescheme_type

  type(bc_swff_type)  , intent(inout) :: bc
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(in) :: RMass(:)
  integer, intent(in) :: tags(2)
  type(timescheme_type), intent(in) :: time

  double precision, parameter :: tiny=1d-3
  double precision, allocatable :: FltCoord(:,:)
  integer :: i,HdrUnit,npoin,NSAMP,NDAT
  character(15) :: OutputFileName,HdrFileName
  
! bc1 --> grid%bounds(i) boundary corresponding to TAG
  call BC_inquire( grid%bounds, tag = tags(1), bc_topo_ptr = bc%bc1 )
  call BC_inquire( grid%bounds, tag = tags(2), bc_topo_ptr = bc%bc2 )

  if (bc%bc1%nelem/=bc%bc2%nelem) &
   call IO_abort('bc_swff_init: number of boundary elements do not match')

  if (bc%bc1%npoin/=bc%bc2%npoin) &
   call IO_abort('bc_swff_init: number of nodes on boundaries do not match')
  npoin = bc%bc1%npoin

  allocate(FltCoord(2,npoin))
  FltCoord = grid%coord(:,bc%bc1%bulk_node)
  if ( any(abs(FltCoord(1,:)-grid%coord(1,bc%bc2%bulk_node))>tiny) &
   .OR.any(abs(FltCoord(2,:)-grid%coord(2,bc%bc2%bulk_node))>tiny) )&
   call IO_abort('bc_swff_init: coordinates on boundaries do not match properly')

! NOTE: the mesh being conformal, the weights B=GLL_weights*jac1D are equal on both
!       sides of the fault. 
  allocate( bc%n1(npoin,2) )
  allocate( bc%B(npoin) ) ! assembled[ GLL_weights * jac1D ]
  call BC_get_normal_and_weights(bc%bc1,grid,bc%n1,bc%B,bc%periodic)

! Update coefficients of Newmark scheme
! vnew = vpredictor +coef3*anew   See time.f90 module
! dnew = dpredictor +coef4*anew   See time.f90 module
  bc%CoefA2V = time%coef(3)
  bc%CoefA2D = time%coef(4)

! RMassi = 1/Mi
! Needed in DaccelFree = -K2*d2/M2 + K1*d1/M1
  allocate(bc%RMass1(npoin))
  allocate(bc%RMass2(npoin))
  bc%RMass1 = RMass(bc%bc1%bulk_node)
  bc%RMass2 = RMass(bc%bc2%bulk_node)

! Fault impedance, Z in :  Trac=TracStick-Z*Dveloc
!   Z = 1/( B1/M1 + B2/M2 )
! TracStick = Z*Vfree traction as if the fault was stuck (no displ discontinuity) 
! NOTE: same Bi on both sides, see note above
  allocate(bc%Z(npoin))
  bc%Z = 1.d0/(bc%CoefA2V*bc%B*( bc%RMass1 + bc%RMass2 ))

! Friction parameters
  call InitConstOrDist(bc%Dc_const,bc%Dc_dist,FltCoord,bc%Dc)
  call InitConstOrDist(bc%MuS_const,bc%MuS_dist,FltCoord,bc%MuStatic)
  call InitConstOrDist(bc%MuD_const,bc%MuD_dist,FltCoord,bc%MuDynamic)

! Initial stress
  call InitConstOrDist(bc%InitStrT_const,bc%InitStrT_dist,FltCoord,bc%InitStrT)
  call InitConstOrDist(bc%InitStrN_const,bc%InitStrN_dist,FltCoord,bc%InitStrN)

! Initial friction
  allocate(bc%MU(npoin))
  bc%MU = bc%MuStatic
 
! Output arrays
  bc%OutputIFir = max(bc%OutputIFir, 1)
  bc%OutputILas = min(bc%OutputILas, npoin)
  bc%OutputNpoin = (bc%OutputILas-bc%OutputIFir)/bc%OutputILag +1
  allocate(bc%StrT(bc%OutputNpoin))
  allocate(bc%StrN(bc%OutputNpoin))
  bc%StrT=bc%InitStrT(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  bc%StrN=bc%InitStrN(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  if (size(bc%InitStrN(bc%OutputIFir:bc%OutputILas:bc%OutputILag))/=bc%OutputNpoin) &
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

!---------------------------------------------------------------------
  subroutine InitConstOrDist(Cin,Din,Coord,Fout)

  use distribution_general, only: DIST_generate, DIST_destructor

  double precision, intent(in) :: Cin,Coord(:,:)
  type(distribution_type), pointer :: Din
  double precision, pointer :: Fout(:)

  allocate( Fout(size(Coord,2)) )
  if (associated(Din)) then
    call DIST_generate(Fout,Coord,Din)
    call DIST_destructor(Din)
  else
    Fout = Cin
  endif

  end subroutine InitConstOrDist

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
! NOTE: eventual periodicity does not need to be enforced at this level
!       because it is assumed that MxAccel and Displ are already periodic
!
  subroutine BC_SWFF_set(bc,MxAccel,Veloc,Displ)

  type(bc_swff_type) :: bc
  double precision, intent(in) :: Veloc(:,:),Displ(:,:)
  double precision, intent(inout) :: MxAccel(:,:)

  double precision, dimension(bc%bc1%npoin) :: Tang,Norm
  double precision, dimension(bc%bc1%npoin,2), target :: tmp1,tmp2
  double precision, dimension(:,:), pointer :: BxTrac,TracStick,Ddispl,Dveloc,Daccel,DaccelFree

  integer, pointer :: i1(:),i2(:)

! NOTE: using pointers here to save memory. Not critical, just good practice.
!       But beware, if you forget it you will scratch data
  Dveloc => tmp1
  DaccelFree => tmp2
  TracStick => tmp1
  BxTrac => tmp2
  Ddispl => tmp1
  Daccel => tmp2

  i1 => bc%bc1%bulk_node
  i2 => bc%bc2%bulk_node

! velocity discontinuity = Face2 - Face1
  Dveloc = Veloc(i2,:) - Veloc(i1,:)
! acceleration discontinuity as if the fault was stress free
  DaccelFree(:,1)= bc%RMass2*MxAccel(i2,1) - bc%RMass1*MxAccel(i1,1)
  DaccelFree(:,2)= bc%RMass2*MxAccel(i2,2) - bc%RMass1*MxAccel(i1,2)
! Traction as if the fault was stuck (no velocity discontinuity)
! TracStick = Z * DvelocFree = Z * (Dveloc+coef*DaccelFree)
  TracStick(:,1) = bc%Z * (Dveloc(:,1) + bc%CoefA2V*DaccelFree(:,1) )
  TracStick(:,2) = bc%Z * (Dveloc(:,2) + bc%CoefA2V*DaccelFree(:,2) )

! Rotate tractions to fault frame (Tang,Norm)
! Add initial stress
  Tang =  bc%n1(:,2)*TracStick(:,1) + bc%n1(:,1)*TracStick(:,2) +bc%InitStrT
  Norm = -bc%n1(:,1)*TracStick(:,1) + bc%n1(:,2)*TracStick(:,2) +bc%InitStrN


!-- BEGIN SOLVE FAULT CONSTITUTIVE LAW ----------------------------
 !Opening-->free stress
  Norm = min(Norm,0.d0) ! negative normal stress is compressive
 !Coulomb friction
 ! NOTE: Friction coefficient MU is updated at the end of this subroutine
 !       The friction solver is explicit: MU is evaluated using 
 !       the slip at the previous timestep
  Tang = sign( min(abs(Tang),-Norm*bc%MU), Tang)
!-- END SOLVE FAULT CONSTITUTIVE LAW ------------------------------

! Subtract initial stress
  Tang = Tang -bc%InitStrT
  Norm = Norm -bc%InitStrN

! Store tractions for output
  bc%StrT = Tang(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  bc%StrN = Norm(bc%OutputIFir:bc%OutputILas:bc%OutputILag)

! Rotate tractions back to (x,z) frame and apply boundary weights
  BxTrac(:,1) = bc%B*( bc%n1(:,2)*Tang - bc%n1(:,1)*Norm )
  BxTrac(:,2) = bc%B*( bc%n1(:,1)*Tang + bc%n1(:,2)*Norm )

! Add boundary term B*T to M*a
  MxAccel(i1,:) = MxAccel(i1,:) + BxTrac
  MxAccel(i2,:) = MxAccel(i2,:) - BxTrac


!-- BEGIN UPDATE FRICTION LAW ------------------------------
! NOTE: update will be used in NEXT timestep
 !Update displacement discontinuity:
  Ddispl = Displ(i2,:) - Displ(i1,:)
  Daccel(:,1)= bc%RMass2*MxAccel(i2,1) - bc%RMass1*MxAccel(i1,1)
  Daccel(:,2)= bc%RMass2*MxAccel(i2,2) - bc%RMass1*MxAccel(i1,2)
  Ddispl = Ddispl + bc%CoefA2D*Daccel
  !WARNING: there should be no update if there is opening, if (Norm>0) return
  !         should be implemented here if the physical situation appears (unlikely)
  Tang = bc%n1(:,2)*Ddispl(:,1) + bc%n1(:,1)*Ddispl(:,2) ! = tangent slip
 !Update slip weakening:
 !NOTE: to implement other slip weakening friction laws just modify next line
  bc%MU = max( bc%MuStatic -(bc%MuStatic-bc%MuDynamic)*abs(Tang)/bc%Dc , bc%MuDynamic) ! linear
!  bc%MU = bc%MuStatic -(bc%MuStatic-bc%MuDynamic)*exp(-Tang/bc%Dc) ! exponential

!-- END UPDATE FRICTION LAW --------------------------------

  end subroutine BC_SWFF_set


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

  i1 => bc%bc1%bulk_node(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  i2 => bc%bc2%bulk_node(bc%OutputIFir:bc%OutputILas:bc%OutputILag)
  n1 => bc%n1(bc%OutputIFir:bc%OutputILas:bc%OutputILag,:)

  write(bc%OutputUnit) real( n1(:,2)*(displ(i2,1)-displ(i1,1)) &
                           + n1(:,1)*(displ(i2,2)-displ(i1,2)) )
  write(bc%OutputUnit) real( n1(:,2)*(veloc(i2,1)-veloc(i1,1)) &
                           + n1(:,1)*(veloc(i2,2)-veloc(i1,2)) )
  write(bc%OutputUnit) real( bc%StrT )
  write(bc%OutputUnit) real( bc%StrN )
  write(bc%OutputUnit) real( bc%MU(bc%OutputIFir:bc%OutputILas:bc%OutputILag) )

  end subroutine BC_SWFF_write

  end module bc_swff
