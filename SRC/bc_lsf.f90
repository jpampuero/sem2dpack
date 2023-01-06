module bc_lsf
! Linear-slip fault

  use bnd_grid, only: bnd_grid_type

  implicit none
  private

  type bc_lsf_type
    double precision :: kt,kn
    double precision, dimension(:,:), pointer :: normal,B,W1,W2
    type(bnd_grid_type), pointer :: bc1,bc2
  end type bc_lsf_type

  double precision, parameter :: inf = huge(inf)

  public :: BC_LSF_type, BC_LSF_read, BC_LSF_init, BC_LSF_apply

contains

!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_LSF 
! GROUP  : BOUNDARY_CONDITION
! PURPOSE: Linear slip fault, a displacement discontinuity interface
!          where stress and slip are linearly related
! SYNTAX : &BC_LSF Ktang | Ctang, Knorm | Cnorm /
!
! ARG: Ktang    [dble] [Inf] Tangential stiffness
! ARG: Ctang    [dble] [0d0] Tangential compliance
! ARG: Knorm    [dble] [Inf] Normal stiffness
! ARG: Cnorm    [dble] [0d0] Normal compliance
!
! NOTE: For each component:
!       You can set K _or_ C, but _not_both_
!       If C=0d0 or K=Inf then no discontinuity is allowed (transparent)
!       If K=0d0 the fault is free stress boundary
!
! END INPUT BLOCK

  subroutine BC_LSF_read(bc,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(bc_lsf_type), intent(out) :: bc
  integer, intent(in) :: iin
  double precision :: Ktang,Knorm,Ctang,Cnorm
  character(13) :: kt_txt, kn_txt

  NAMELIST / BC_LSF /  Ktang,Knorm,Ctang,Cnorm

  Ctang = 0d0
  Cnorm = 0d0
  Ktang = inf          ! the default is : welded contact
  Knorm = inf

  read(iin,BC_LSF,END=100)
  
  if (Ctang/=0d0) Ktang= 1d0/Ctang ! Compliance input overwrites stiffness input
  if (Ktang/=inf) then
    bc%kt   = Ktang
    write(kt_txt,'(EN12.3)') Ktang
  else
    kt_txt = 'Inf'
  endif

  if (Cnorm/=0d0) Knorm= 1d0/Cnorm 
  if (Knorm/=inf) then
    bc%kn   = Knorm
    write(kn_txt,'(EN12.3)') Knorm
  else
    kn_txt = 'Inf'
  endif

  if (echo_input) write(iout,200) kt_txt, kn_txt

  return
  100 call IO_abort('BC_LSF_read: BC_LSF input block not found')
  200 format(5x,'Type   = Linear slip fault', &
            /5x,'  Tangential stiffness = ',A,&
            /5x,'  Normal stiffness     = ',A)

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
  double precision, intent(in) :: M(:,:)
  integer, intent(in) :: tags(2)
  type(bc_periodic_type), pointer :: perio

  integer :: ndof 

  ndof = size(M,2)

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
  allocate( bc%B(bc%bc1%npoin,ndof) ) ! assembled[ GLL_weights * jac1D ]
  call BC_get_normal_and_weights(bc%bc1,grid,bc%normal,bc%B(:,1), BC_PERIO_intersects(bc%bc1,perio))
! NOTE: the mesh being conformal, the [GLL_weights*jac1D] are equal on both
!       sides of the fault. 
  if (ndof==2) bc%B(:,2) = bc%B(:,1)

 !The following are needed to compute the traction 
 !as if the fault was stuck (no displ discontinuity)
 !Tstick = ( -K2*d2/M2 + K1*d1/M1 )/( B1/M1 + B2/M2 )
 !Here we store: bc%Wi = 1/Mi /( B1/M1 + B2/M2 )
  allocate(bc%W1(bc%bc1%npoin,ndof))
  allocate(bc%W2(bc%bc2%npoin,ndof))
! NOTE: B1=B2, same GLL_weights*jac1D on both sides, see note above
! NOTE: assuming periodic boundaries have been already initialized 
!       (M and B respect periodicity)
  bc%W1 = 1d0/( bc%B*( 1d0+M(bc%bc1%node,:)/M(bc%bc2%node,:) ))
  bc%W2 = 1d0/( bc%B*( M(bc%bc2%node,:)/M(bc%bc1%node,:)+1d0 ))

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
!       If Newmark: we use the PREDICTED displacement to compute the slip
!       If leapfrog: second order ok
!
! NOTE: possible periodicity does not need to be enforced at this level
!       because it is assumed that MxA and D are already periodic
!
  subroutine BC_LSF_apply(bc,MxA,D)

  type(bc_lsf_type), intent(in) :: bc
  double precision, intent(in) :: D(:,:)
  double precision, intent(inout) :: MxA(:,:)

  double precision, dimension(bc%bc1%npoin,size(D,2)) :: T,dD

 ! Traction as if the fault was welded (no slip acceleration)
 !   Tstick = ( -K2*d2/M2 + K1*d1/M1 )/( B1/M1 + B2/M2 )
 ! On entry MxA = -K*d 
 !          bc%Wi = 1/Mi /( B1/M1 + B2/M2 ) 
  T = bc%W2*MxA(bc%bc2%node,:) - bc%W1*MxA(bc%bc1%node,:)

 ! slip = displacement discontinuity = Face2 - Face1
  dD = D(bc%bc2%node,:) - D(bc%bc1%node,:)

  if (size(D,2)==2) then

   ! rotate to fault frame (Tt,Tn)
    T = rotate(bc,T,1)
    dD = rotate(bc,dD,1)
 
   ! apply the linear slip fault condition
    if (bc%kt<inf) T(:,1) = bc%kt * dD(:,1)
    if (bc%kn<inf) T(:,2) = bc%kn * dD(:,2)

   ! rotate tractions back to original frame
    T = rotate(bc,T,-1)

  else
    if (bc%kt<inf) T(:,1) = bc%kt * dD(:,1)

  endif
  
 ! Add boundary term B*T to M*a
  MxA(bc%bc1%node,:) = MxA(bc%bc1%node,:) + bc%B*T
  MxA(bc%bc2%node,:) = MxA(bc%bc2%node,:) - bc%B*T

  end subroutine BC_LSF_apply


!---------------------------------------------------------------------
  function rotate(bc,v,fb) result(vr)

  type(bc_lsf_type), intent(in) :: bc
  double precision, intent(in) :: v(bc%bc1%npoin,2)
  integer, intent(in) :: fb
  double precision :: vr(bc%bc1%npoin,2)

  if (fb==1) then
    vr(:,1) = bc%normal(:,2)*v(:,1) - bc%normal(:,1)*v(:,2)
    vr(:,2) = bc%normal(:,1)*v(:,1) + bc%normal(:,2)*v(:,2)
  else
    vr(:,1) = bc%normal(:,2)*v(:,1) + bc%normal(:,1)*v(:,2)
    vr(:,2) =-bc%normal(:,1)*v(:,1) + bc%normal(:,2)*v(:,2)
  endif

  end function rotate


  end module bc_lsf
