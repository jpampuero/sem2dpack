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
module bc_abso

! Absorbing boundaries using first order derivatives
! Reference: R. Stacey, BSSA Vol.78, No.6, pp. 2089-2097, Dec 1988

! P1. Clayton and Engquist, first order:
!   V_t = - Cs * U_t,n
!   V_n = - Cp * U_n,n

! P3. Stacey, second order:
!   V_t = - Cs * U_t,n + (Cp-Cs) * U_n,t
!   V_n = - Cp * U_n,n + (Cp-Cs) * U_t,t

! Note: this is implemented as a constraint on acceleration
! after deriving once wrt time
! A different implementation would be as a boundary term (virtual work
! of tractions)
! P1. T_t = - mu/Cs * V_t
!     T_n = - mu/Cp * V_n
! P3. see Kausel, BSSA 1992, eq.7

  use spec_grid
  use memory_info

  implicit none
  private
  
  type bc_abso_type
    private
    double precision, dimension(:,:), pointer :: &
      Cx_DetaDn,Cx_DxiDn,Cx_DetaDt,Cx_DxiDt &
     ,Cz_DetaDn,Cz_DxiDn,Cz_DetaDt,Cz_DxiDt
    double precision, dimension(:), pointer :: cdamp
    type(bc_topo_type), pointer :: topo
    integer :: side
    logical :: stacey,Periodic
  end type

  integer, parameter :: ndime = 2

  public :: BC_ABSO_type, BC_ABSO_read, BC_ABSO_init, BC_ABSO_set
    
contains

!=====================================================================
!   
! BEGIN INPUT BLOCK
!
! NAME   : BC_ABSORB [boundary condition]
! PURPOSE: Absorbing boundary
! SYNTAX : &BC_ABSORB side,stacey,periodic /
!
! ARG: side     [char] [none] Which side of the model corresponds to this
!               boundary:       'U'     Up,top
!                               'D'     Down,bottom
!                               'L'     Left
!                               'R'     Right
! ARG: stacey   [log] [T] Apply or not Stacey absorbing conditions. If F, we
!                will simply use Clayton and Engquist.
! ARG: Periodic [log] [F]  Enable if the fault crosses a periodic boundary
!
! NOTE   : Only implemented for vertical and horizontal boundaries.
!
! END INPUT BLOCK

  subroutine BC_ABSO_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_abso_type), intent(out) :: bc
  integer, intent(in) :: iin

  integer :: i
  logical :: stacey,periodic
  character :: side
  character(20) :: abs_type

  NAMELIST / BC_ABSORB / side,stacey,periodic

  side = ''
  stacey  = .true.
  periodic = .false.

  read(iin,BC_ABSORB,END=100)
  
  select case (side)
  case('U'); bc%side =  side_U
  case('D'); bc%side =  side_D
  case('L'); bc%side =  side_L
  case('R'); bc%side =  side_R
  case default; call IO_abort('BC_ABSO_read: side must be U,D,R or L')
  end select

  bc%stacey = stacey
  bc%periodic = periodic

  if (echo_input) then 
    if (stacey) then
      abs_type = 'Stacey'
    else
      abs_type = 'Clayton-Engquist'
    endif
    write(iout,200) abs_type,periodic
  endif

  return

  100 call IO_abort('BC_ABSO_read: no BC_ABSORB input block found') 
  200 format(6x,'Type = Absorbing boundary ',A &
            /6x,'Periodic =',L1 )

  end subroutine BC_ABSO_read



!=====================================================================
!
!--- Definition of the working coefficients
!
  subroutine BC_ABSO_init(bc,tag,grid,elast)

  use elastic, only : elast_type, ELAST_inquire

  type(bc_abso_type) , intent(inout) :: bc
  type(sem_grid_type), intent(in)    :: grid
  type(elast_type)   , intent(in)    :: elast
  integer            , intent(in)    :: tag

  double precision, allocatable :: a15(:,:)
  double precision :: CoefIntegr,CoefStacey,sign
  double precision, target :: ijac(ndime,ndime),jac(ndime,ndime)
  double precision, pointer :: DxiDt,DxiDn,DetaDt,DetaDn
  double precision :: c(2)
  integer :: i,j,k,e,LocDimTanToEdge,GeoDimTanToSide,GeoDimNorToSide &
            ,ngll,bc_nelem,ebulk,itab(grid%ngll),jtab(grid%ngll)

  !-- bc%topo => grid%bounds(i) corresponding to TAG
  call BC_inquire( grid%bounds, tag = tag, bc_topo_ptr = bc%topo )

  ngll         = grid%ngll
  bc_nelem     = bc%topo%nelem

  allocate(bc%Cx_DetaDn(ngll,bc_nelem))
  allocate(bc%Cx_DxiDn(ngll,bc_nelem))
  allocate(bc%Cz_DetaDn(ngll,bc_nelem))
  allocate(bc%Cz_DxiDn(ngll,bc_nelem))

  call storearray('bc%Cx_DetaDn',ngll*bc_nelem,idouble)
  call storearray('bc%Cx_DxiDn',ngll*bc_nelem,idouble)
  call storearray('bc%Cz_DetaDn',ngll*bc_nelem,idouble)
  call storearray('bc%Cz_DxiDn',ngll*bc_nelem,idouble)

  if (bc%stacey) then

    allocate(bc%Cx_DxiDt(ngll,bc_nelem))
    allocate(bc%Cx_DetaDt(ngll,bc_nelem))
    bc%Cz_DxiDt  => bc%Cx_DxiDt
    bc%Cz_DetaDt => bc%Cx_DetaDt
    call storearray('bc%Cx_DxiDt',ngll*bc_nelem,idouble)
    call storearray('bc%Cx_DetaDt',ngll*bc_nelem,idouble)

  endif

  allocate(a15(ngll,bc_nelem)) 

  do e =1,bc_nelem

    ebulk = bc%topo%bulk_element(e)

    call SE_inquire(grid, edge=bc%topo%element_edge(e) &
                   ,itab=itab, jtab=jtab, dim_t=LocDimTanToEdge)
   ! LocDimTanToEdge = local dimension (xi,eta) tangent to the edge

    select case(bc%side)
      case(side_U,side_D)
        GeoDimTanToSide = 1 ! Geographical dimension (x,z) tangent to domain side
        GeoDimNorToSide = 2 ! Geographical dimension (x,z) normal to domain side
      case(side_L,side_R)
        GeoDimTanToSide = 2
        GeoDimNorToSide = 1
    end select
    select case(bc%side)
      case(side_U,side_R); sign = -1.d0
      case(side_D,side_L); sign =  1.d0
    end select

    do k = 1,ngll

      i = itab(k)
      j = jtab(k)

      call ELAST_inquire(elast,i,j,ebulk, cp=c(GeoDimNorToSide), cs=c(GeoDimTanToSide))

     ! ijac = DLocDGlob
     ! jac  = DGlobDLoc
      call SE_inquire(grid,element=ebulk,igll=i,jgll=j, inv_jacobian=ijac, jacobian=jac)

      DxiDt   => ijac(1,GeoDimTanToSide)
      DxiDn   => ijac(1,GeoDimNorToSide)
      DetaDt  => ijac(2,GeoDimTanToSide)
      DetaDn  => ijac(2,GeoDimNorToSide)

     ! the differential length has this simplified expression for
     ! vertical or horizontal edges:
      CoefIntegr       = grid%wgll(k)*abs(jac(GeoDimTanToSide,LocDimTanToEdge))

      a15(k,e)         = CoefIntegr

      bc%Cx_DxiDn(k,e)  = sign*c(1)*CoefIntegr*DxiDn
      bc%Cx_DetaDn(k,e) = sign*c(1)*CoefIntegr*DetaDn

      bc%Cz_DxiDn(k,e)  = sign*c(2)*CoefIntegr*DxiDn
      bc%Cz_DetaDn(k,e) = sign*c(2)*CoefIntegr*DetaDn
      
      if (bc%stacey) then
        CoefStacey = c(GeoDimNorToSide)-c(GeoDimTanToSide)  ! Cp - Cs
        bc%Cx_DxiDt(k,e)  = sign*CoefIntegr*DxiDt*CoefStacey
        bc%Cx_DetaDt(k,e) = sign*CoefIntegr*DetaDt*CoefStacey
      endif

    enddo
  enddo

!--- assemble damping diagonal matrix
  allocate(bc%cdamp(bc%topo%npoin))
  call storearray('bc%cdamp',bc%topo%npoin,idouble)
  bc%cdamp = 0.d0
  do e=1,bc_nelem
    bc%cdamp( bc%topo%ibool(:,e) ) = bc%cdamp( bc%topo%ibool(:,e) ) + a15(:,e)
  enddo
  if (bc%periodic) then
    bc%cdamp(1) = bc%cdamp(1) + bc%cdamp(bc%topo%npoin)
    bc%cdamp(bc%topo%npoin) = bc%cdamp(1)
  endif

  deallocate(a15)

  end subroutine BC_ABSO_init



!=====================================================================
!
  subroutine BC_ABSO_set(bc,grid,fields,src,time)

  use sources, only : source_type,SO_WAVE_veloc,SO_WAVE_accel,SO_inquire
  use fields_class, only : fields_type

  type(bc_abso_type) , intent(in)    :: bc
  type(sem_grid_type), intent(in)    :: grid
  type(fields_type)  , intent(inout) :: fields
  type(source_type)  , pointer       :: src
  double precision   , intent(in)    :: time

  double precision, dimension(bc%topo%ngll,bc%topo%nelem) ::  &
    dUx_dxi,dUz_dxi,dUx_deta,dUz_deta
  double precision :: Uloc(bc%topo%ngll,bc%topo%ngll,ndime) &
                     ,Floc(bc%topo%ngll,bc%topo%nelem,ndime) &
                     ,Fglob(bc%topo%npoin,ndime)
  integer :: i,j,k,iglob,d,itab(bc%topo%ngll),jtab(bc%topo%ngll)
  logical :: anywave

  if (associated(src)) then
    call SO_inquire(src,is_wave=anywave)
  else
    anywave = .false.
  endif

  do k=1,bc%topo%nelem
    
    do j=1,bc%topo%ngll
    do i=1,bc%topo%ngll
      iglob = grid%ibool(i,j,bc%topo%bulk_element(k))
      Uloc(i,j,:) = fields%veloc(iglob,:)
      if (anywave) Uloc(i,j,:) = Uloc(i,j,:) - SO_WAVE_veloc(time,grid%coord(:,iglob),src)
    enddo
    enddo

    call SE_inquire(grid,edge = bc%topo%element_edge(k),itab=itab,jtab=jtab)
    do i= 1,bc%topo%ngll
      dUx_dxi(i,k)  = dot_product( Uloc(:,jtab(i),1), grid%hprime(:,itab(i)) )
      dUz_dxi(i,k)  = dot_product( Uloc(:,jtab(i),2), grid%hprime(:,itab(i)) )
      dUx_deta(i,k) = dot_product( Uloc(itab(i),:,1), grid%hprime(:,jtab(i)) )
      dUz_deta(i,k) = dot_product( Uloc(itab(i),:,2), grid%hprime(:,jtab(i)) )
    enddo

  enddo

  if (bc%stacey) then
! P3. Stacey, second order:
!   V_t = - Cs * U_t,n + (Cp-Cs) * U_n,t
!   V_n = - Cp * U_n,n + (Cp-Cs) * U_t,t
    Floc(:,:,1) = bc%Cx_DetaDn*dUx_deta + bc%Cx_DxiDn *dUx_dxi &
                + bc%Cx_DxiDt *dUz_dxi  + bc%Cx_DetaDt*dUz_deta
    Floc(:,:,2) = bc%Cz_DetaDn*dUz_deta + bc%Cz_DxiDn *dUz_dxi &
                + bc%Cz_DxiDt *dUx_dxi  + bc%Cz_DetaDt*dUx_deta
  else
! P1. Clayton and Engquist, first order:
!   V_t = - Cs * U_t,n
!   V_n = - Cp * U_n,n
    Floc(:,:,1) = bc%Cx_DetaDn*dUx_deta + bc%Cx_DxiDn*dUx_dxi
    Floc(:,:,2) = bc%Cz_DetaDn*dUz_deta + bc%Cz_DxiDn*dUz_dxi
  endif

 ! assemblage des contributions
  Fglob = 0.d0 
  do k=1,bc%topo%nelem
    Fglob(bc%topo%ibool(:,k),:) = Fglob(bc%topo%ibool(:,k),:) + Floc(:,k,:)
  enddo

  fields%accel(bc%topo%bulk_node,1) = Fglob(:,1) / bc%cdamp 
  fields%accel(bc%topo%bulk_node,2) = Fglob(:,2) / bc%cdamp 

  if (anywave) then
    do i= 1,bc%topo%npoin
      iglob = bc%topo%bulk_node(i)
      fields%accel(iglob,:) = fields%accel(iglob,:) + SO_WAVE_accel(time,grid%coord(:,iglob),src)
    enddo
  endif

  end subroutine BC_ABSO_set

end module bc_abso
