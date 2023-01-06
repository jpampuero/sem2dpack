module src_moment

! SRC_MOMENT: moment tensor sources
! + special cases: explosion, double-couple, tensile

  implicit none
  private

  type so_moment_type
    private
    double precision, pointer :: M(:,:)=>null()
    double precision, dimension(:,:,:), pointer :: coef_xi=>null(), coef_eta=>null()
    integer, dimension(:,:), pointer :: iglob_xi=>null(), iglob_eta=>null()
  end type so_moment_type

  public :: so_moment_type,SRC_MOMENT_read,SRC_MOMENT_init,SRC_MOMENT_add

contains

!=====================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : SRC_DOUBLE_COUPLE
! GROUP  : SOURCE MECHANISM
! PURPOSE: Define a double-couple source
! SYNTAX : &SRC_DOUBLE_COUPLE  dip /
!
! ARG: dip           [dble] [90] Dip angle, in degrees, clockwise 
!                    from the positive X direction
!                     
! NOTE   : Sign convention: if the source amplitude is positive the right block
!          moves up (positive Z direction) in PSV and forward (positive Y 
!          direction) in SH.
!
! NOTE   : The source time function gives the cumulative seismic moment Mo(t), 
!          NOT the seismic moment rate.
!
! NOTE   : The seismic moment Mo must be rescaled because a 2D point source is 
!          equivalent to a 3D line source. A proper scaling is obtained by 
!          dividing the original 3D moment by the characteristic size of the 
!          rupture area in the off-plane dimension. An approximate scaling for
!          a fault area with aspect ratio close to unity is
!            Mo_2D = (Mo_3D/dtau)^2/3 * dtau
!          where dtau is the stress drop (typically a few MPa).
!
! END INPUT BLOCK


! BEGIN INPUT BLOCK
!
! NAME   : SRC_MOMENT
! GROUP  : SOURCE MECHANISM
! PURPOSE: Define a moment tensor source
! SYNTAX : &SRC_MOMENT Mxx,Mxz,Mzx,Mzz /
!          &SRC_MOMENT Myx,Myz /
!
! ARG: Mxx,Mxz,Mzx,Mzz [dble] [0] Tensor components for PSV
! ARG: Myx,Myz         [dble] [0] Tensor components for SH
!                     
! END INPUT BLOCK

  subroutine SRC_MOMENT_read(so,iin,mode,ndof)
  
  use echo, only : echo_input,iout
  use stdio, only : IO_abort
  use constants, only : PI

  type(so_moment_type), intent(inout) :: so
  integer, intent(in) :: iin,ndof
  character(*), intent(in) :: mode

  double precision :: Mxx,Mxz,Mzx,Mzz , Myx,Myz, dip, n(2),r(2)

  NAMELIST / SRC_MOMENT / Mxx,Mxz,Mzx,Mzz , Myx,Myz
!  NAMELIST / SRC_EXPLOSION / 
  NAMELIST / SRC_DOUBLE_COUPLE / dip

  Mxx = 0d0
  Mxz = 0d0
  Mzx = 0d0
  Mzz = 0d0
  Myx = 0d0
  Myz = 0d0
  dip = 90d0

! define moment tensor components
  allocate(so%M(2,ndof))
  select case (mode)

    case('EXPLOSION') 
      if (echo_input) write(iout,200)
      if (ndof==2) then
        so%M(1,1) = 1d0 
        so%M(1,2) = 0d0 
        so%M(2,1) = 0d0 
        so%M(2,2) = 1d0 
      else
        call IO_abort('SRC_MOMENT_read: explosion only allowed in PSV (ndof=2)')
      endif

    case('DOUBLE_COUPLE') 
      read(iin,SRC_DOUBLE_COUPLE,END=100)
      if (echo_input) write(iout,210) dip
      dip = dip*PI/180d0
      ! n=fault_normal vector, points from left block to right block
      n(1) = sin(dip)
      n(2) = cos(dip)
      if (ndof==2) then !PSV
        ! r=slip_vector, if ampli>0 the block on the right moves upward
        r(1) = -cos(dip)
        r(2) =  sin(dip)
        so%M(1,1) = 2d0*r(1)*n(1)
        so%M(1,2) = r(1)*n(2) + r(2)*n(1)
        so%M(2,1) = so%M(1,2)
        so%M(2,2) = 2d0*r(2)*n(2)
      else !SH
        ! right block moves forward (y>0 direction)
        so%M(1,1) = n(1)  
        so%M(2,1) = n(2)
      endif

    case('MOMENT')
      read(iin,SRC_MOMENT,END=101)
      if (ndof==2) then !PSV
        if (echo_input) write(iout,220) Mxx,Mxz,Mzx,Mzz
        so%M(1,1) = Mxx
        so%M(1,2) = Mxz
        so%M(2,1) = Mzx
        so%M(2,2) = Mzz
      else !SH
        if (echo_input) write(iout,230) Myx,Myz
        so%M(1,1) = Myx
        so%M(2,1) = Myz
      endif

    case default
      call IO_abort('SRC_MOMENT_read: unknown type')
  end select


  return

 100 call IO_abort('SRC_MOMENT_read: SRC_DOUBLE_COUPLE input block not found')
 101 call IO_abort('SRC_MOMENT_read: SRC_MOMENT input block not found')
 200 format(5x, &
     'Source Type. . . . . . . . (mechanism) = explosion')
 210 format(5x, &
     'Source Type. . . . . . . . (mechanism) = double couple',/5x, &
     'Dip angle. . . . . . . . . . . . (dip) = ',F0.2)
 220 format(5x, &
     'Source Type. . . . . . . . (mechanism) = moment tensor',/5x, &
     'Mxx . . . .  . . . . . . . . . . . . . = ',EN12.3,/5x, &
     'Mxz . . . .  . . . . . . . . . . . . . = ',EN12.3,/5x, &
     'Mzx . . . .  . . . . . . . . . . . . . = ',EN12.3,/5x, &
     'Mzz . . . .  . . . . . . . . . . . . . = ',EN12.3)
 230 format(5x, &
     'Source Type. . . . . . . . (mechanism) = moment tensor',/5x, &
     'Myx . . . .  . . . . . . . . . . . . . = ',EN12.3,/5x, &
     'Myz . . . .  . . . . . . . . . . . . . = ',EN12.3)

  end subroutine SRC_MOMENT_read

!=====================================================================
!-- Define the working arrays 
!   Spatial distribution: Dirac (cross GLL stencil)
!     so%coef_  = coefficients
!     so%iglob_ = global node indices
!
  subroutine SRC_MOMENT_init(so,iglob,grid)

  use spec_grid, only : sem_grid_type,SE_InverseJacobian,SE_node_belongs_to

  type(so_moment_type), intent(inout) :: so
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: iglob
  
  double precision :: jac_inv(2,2), G(2,2)
  integer :: ndof,nel,k,e,i,j
  integer, pointer, dimension(:) :: itab,jtab,etab

  call SE_node_belongs_to(iglob,etab,itab,jtab,grid)
  nel = size(etab)
  ndof = size(so%M,2)

  allocate(so%iglob_xi(grid%ngll,nel))
  allocate(so%iglob_eta(grid%ngll,nel))
  allocate(so%coef_xi(grid%ngll,ndof,nel))
  allocate(so%coef_eta(grid%ngll,ndof,nel))

  do k=1,nel

    e = etab(k)
    i = itab(k)
    j = jtab(k)

    so%iglob_xi(:,k)  = grid%ibool(:,j,e)
    so%iglob_eta(:,k) = grid%ibool(i,:,e)
  
    jac_inv = SE_InverseJacobian(grid,e,i,j)
  
    if (ndof==2) then ! PSV
      G = matmul(so%M,transpose(jac_inv))
      so%coef_xi(:,1,k)  = G(1,1) * grid%hprime(:,i) ! Mxx*dxi_dx  + Mxz*dxi_dz
      so%coef_eta(:,1,k) = G(1,2) * grid%hprime(:,j) ! Mxx*deta_dx + Mxz*deta_dz
      so%coef_xi(:,2,k)  = G(2,1) * grid%hprime(:,i) ! Mzx*dxi_dx  + Mzz*dxi_dz
      so%coef_eta(:,2,k) = G(2,2) * grid%hprime(:,j) ! Mzx*deta_dx + Mzz*deta_dz
  
    else ! SH
      G(:,1) = matmul(jac_inv,so%M(:,1))
      so%coef_xi(:,1,k)  = G(1,1) * grid%hprime(:,i) ! Myx*dxi_dx  + Myz*dxi_dz
      so%coef_eta(:,1,k) = G(2,1) * grid%hprime(:,j) ! Myx*deta_dx + Myz*deta_dz
  
    endif
  enddo

  deallocate(etab,itab,jtab)

  end subroutine SRC_MOMENT_init

!=====================================================================

  subroutine SRC_MOMENT_add(so,ampli,MxA)

  type(so_moment_type), intent(inout) :: so
  double precision, intent(in) :: ampli
  double precision, intent(inout) :: MxA(:,:)

  integer :: nel,k

  nel = size(so%iglob_xi,2)
  do k=1,nel
    MxA(so%iglob_xi(:,k) ,:) = MxA(so%iglob_xi(:,k),:)  + ampli * so%coef_xi(:,:,k) 
    MxA(so%iglob_eta(:,k),:) = MxA(so%iglob_eta(:,k),:) + ampli * so%coef_eta(:,:,k) 
  enddo

  end subroutine SRC_MOMENT_add


end module src_moment
