module mat_elastic
! Linear elasticity
! Isotropic or transverse anisotropic
  
  use prop_mat

  implicit none
  private
  
 !-- purely elastic
  type matwrk_elast_type
    private
    double precision, pointer :: a(:,:,:) => null(), beta(:,:) => null()
  end type matwrk_elast_type

  ! defined in the first call to MAT_setKind
  integer, save :: isElastic = 0
  integer, save :: isElasticIsotropic = 0
  integer, save :: isElasticAnisotropic = 0
  integer, save :: isElasticHomogeneous = 0 ! needed for optimizations

  ! for memory report
  integer, save :: MAT_ELAST_mempro = 0
  integer, save :: MAT_ELAST_memwrk = 0

  public :: matwrk_elast_type &
          , MAT_isElastic, MAT_isElasticHomogeneous &
          , MAT_ELAST_read &
          , MAT_ELAST_init_elem_prop, MAT_ELAST_init_elem_work &
          , MAT_ELAST_f, MAT_ELAST_stress &
          , MAT_ELAST_memwrk, MAT_ELAST_mempro &
          , MAT_ELAST_add_25D_f


contains

!=======================================================================
  logical function MAT_isElastic(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isElastic = MAT_isKind(m,isElastic)
  end function MAT_isElastic

!-----------------------------------------------------------------------
  logical function MAT_isElasticHomogeneous(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isElasticHomogeneous = MAT_isKind(m,isElasticHomogeneous)
  end function MAT_isElasticHomogeneous

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_ELASTIC
! GROUP  : MATERIALS
! PURPOSE: Set material properties for a linear elastic medium
! SYNTAX : For isotropic material:
!           &MAT_ELASTIC rho|rhoH, cp|cpH, cs|csH /
!          For transverse anisotropy with vertical symmetry axis:
!           &MAT_ELASTIC rho|rhoH, c11|c11H, c13|c13H, c33|c33H, c55|c55H, c66|c66H /
!          Followed by one DIST_XXXX blocks for each argument present with suffix H,
!          in the same order as listed above.
!
! ARG: cp       [dble][0d0] P wave velocity (m/s)
! ARG: cs       [dble][0d0] S wave velocity (m/s)
! ARG: rho      [dble][0d0] density (kg/m^3)
! ARG: c11,c13,c33,c55,c66  [dble][0d0] anisotropic elastic moduli (Pa)
!
! END INPUT BLOCK

  subroutine MAT_ELAST_read(input,iin)
  
  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type(matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: Kmod,poiss,rho,cp,cs,mu,lambda
  double precision :: c11,c13,c33,c55,c66
  character(20) :: rhoH,cpH,csH,c11H,c13H,c33H,c55H,c66H

  NAMELIST / MAT_ELASTIC / rho,cp,cs, rhoH,cpH,csH, c11,c13,c33,c55,c66, c11H,c13H,c33H,c55H,c66H

  call MAT_setKind(input,isElastic)

  rho = 0d0
  cp = 0d0
  cs = 0d0
  rhoH = ''
  cpH = ''
  csH = ''
  c11 = 0d0
  c13 = 0d0
  c33 = 0d0
  c55 = 0d0
  c66 = 0d0
  c11H = ''
  c13H = ''
  c33H = ''
  c55H = ''
  c66H = ''

  read(iin,MAT_ELASTIC,END=100)

  if (rho<=0d0 .and. rhoH=='') call IO_abort('MAT_ELAST_read: undefined density (rho)')
  call MAT_setProp(input,'rho',rho,rhoH,iin,rhoH)

 ! Isotropic material: cp, cs (m/s)
  if ( (cp>0d0 .or. cpH/='') .and. (cs>0d0 .or. csH/='') ) then

    call MAT_setKind(input,isElasticIsotropic)

    call MAT_setProp(input,'cp',cp,cpH,iin,cpH)
    call MAT_setProp(input,'cs',cs,csH,iin,csH)

    if (echo_input) write(iout,150) cpH,csH,rhoH

   ! if uniform properties
    if (cp>0d0 .and. cs>0d0 .and. rho>0d0) then
      call MAT_setKind(input,isElasticHomogeneous)
      mu   = rho*cs*cs
      Kmod  = rho*cp*cp
      lambda  = rho*(cp*cp - 2d0*cs*cs)
      poiss = 0.5d0*lambda/(Kmod-mu)  !(cp*cp-2d0*cs*cs)/(cp*cp-cs*cs)
      if (poiss < 0.d0 .or. poiss > 0.5d0) call IO_abort('Poisson''s ratio out of range !')
      if (echo_input) write(iout,200) poiss,lambda,mu, &
                     lambda+2d0*mu/3d0, 2d0*mu*(1d0+poiss)  ! Kvol,young  
      call MAT_setProp(input,'lambda',lambda)
      call MAT_setProp(input,'mu',mu)
    endif

 ! Anisotropic material: c11, c13, c33, c55, c66 (Pa)
  elseif ( ((c11>0d0 .or. c11H/='') .and. (c13>0d0 .or. c13H/='') .and. &
            (c33>0d0 .or. c33H/='') .and. (c55>0d0 .or. c55H/='')) .or. &
           ((c55>0d0 .or. c55H/='') .and. (c66>0d0 .or. c66H/='')) ) then

    call MAT_setKind(input,isElasticAnisotropic)

    call MAT_setProp(input,'c11',c11,c11H,iin,c11H)
    call MAT_setProp(input,'c13',c13,c13H,iin,c13H)
    call MAT_setProp(input,'c33',c33,c33H,iin,c33H)
    call MAT_setProp(input,'c55',c55,c55H,iin,c55H)
    call MAT_setProp(input,'c66',c66,c66H,iin,c66H)

    if (echo_input) write(iout,250) c11H,c13H,c33H,c55H,c66H,rhoH

   ! if uniform properties
    if (c11>0d0 .and. c13>0d0 .and. c33>0d0 .and. c55>0d0 .and. c66>0d0 .and. rho>0d0) then
      call MAT_setKind(input,isElasticHomogeneous)
      if (echo_input) write(iout,300) sqrt(c33/rho),sqrt(c11/rho),sqrt(c55/rho),sqrt(c66/rho)
    endif
  
  else
    call IO_abort('MAT_ELAST_read: incomplete input')
  endif

  return

  100 call IO_abort('MAT_ELAST_read: MAT_ELASTIC input block not found')

  150   format(5x, &
    '                                            Isotropic',/5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',A,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',A,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',A)

  200   format(5x, &
    'Poisson''s ratio . . . . . . . . . . . . . =',EN12.3,/5x, &
    'First Lame parameter Lambda . . . . . . . =',EN12.3,/5x, &
    'Second Lame parameter Mu. . . . . . . . . =',EN12.3,/5x, &
    'Bulk modulus K. . . . . . . . . . . . . . =',EN12.3,/5x, &
    'Young''s modulus E . . . . . . . . . . . . =',EN12.3)

  250   format(5x, &
    '                                            Anisotropic',/5x, &
    'c11 coefficient (Pascal). . . . . . (c11) =',A,/5x, &
    'c13 coefficient (Pascal). . . . . . (c13) =',A,/5x, &
    'c33 coefficient (Pascal). . . . . . (c33) =',A,/5x, &
    'c55 coefficient (Pascal). . . . . . (c55) =',A,/5x, &
    'c66 coefficient (Pascal). . . . . . (c66) =',A,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',A)

  300   format(5x, &
    'Velocity of qP along vertical axis. . . . =',EN12.3,/5x, &
    'Velocity of qP along horizontal axis. . . =',EN12.3,/5x, &
    'Velocity of qSV and vertical qSH  . . . . =',EN12.3,/5x, &
    'Velocity of qSH along horizontal axis . . =',EN12.3)

  end subroutine MAT_ELAST_read


!=======================================================================
! Initialize material properties
!
  subroutine MAT_ELAST_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision, dimension(size(ecoord,2),size(ecoord,3)) :: cp,cs,rho,c11,c55
  
  if (MAT_isKind(elem,isElasticIsotropic)) then

    call MAT_setProp(elem,'cp',ecoord,MAT_ELAST_mempro)
    call MAT_setProp(elem,'cs',ecoord,MAT_ELAST_mempro)

    if (MAT_isProp_input('lambda',elem)) then
      call MAT_setProp(elem,'lambda',ecoord,MAT_ELAST_mempro)
      call MAT_setProp(elem,'mu',ecoord,MAT_ELAST_mempro)

    else
      call MAT_getProp(rho,elem,'rho')
      call MAT_getProp(cp,elem,'cp')
      call MAT_getProp(cs,elem,'cs')
      call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs),MAT_ELAST_mempro)
      call MAT_setProp(elem,'mu',rho*cs*cs,MAT_ELAST_mempro)
    endif

  else
    call MAT_setProp(elem,'c13',ecoord,MAT_ELAST_mempro)
    call MAT_setProp(elem,'c55',ecoord,MAT_ELAST_mempro)
    call MAT_setProp(elem,'c66',ecoord,MAT_ELAST_mempro)
    call MAT_setProp(elem,'c11',ecoord,MAT_ELAST_mempro)
    call MAT_setProp(elem,'c33',ecoord,MAT_ELAST_mempro)

   ! a rough idea of wave speeds for timestep setting, ABCs and plots
    call MAT_getProp(rho,elem,'rho')
    call MAT_getProp(c11,elem,'c11')
    call MAT_getProp(c55,elem,'c55')
    call MAT_setProp(elem,'cp',sqrt(c11/rho),MAT_ELAST_mempro)
    call MAT_setProp(elem,'cs',sqrt(c55/rho),MAT_ELAST_mempro)
  
  endif
      
  end subroutine MAT_ELAST_init_elem_prop

!=======================================================================
!
!  Define arrays a1 to a10 for the computation of elastic internal forces
!
! ndof = 1 : SH
!      = 2 : P-SV
  subroutine MAT_ELAST_init_elem_work(matwrk,matpro,grid,e,ndof,flat_grid)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights
  use stdio, only : IO_abort

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elast_type), intent(inout) :: matwrk
  type(matpro_elem_type), intent(in) :: matpro
  integer, intent(in) :: e,ndof
  logical, intent(in) :: flat_grid

  integer :: nelast

 ! if box mesh ignore the crossed terms
  if (flat_grid) then
    if (ndof==1) then
      nelast = 2
    else 
      nelast = 6
    endif
  else
    if (ndof==1) then
      nelast = 3
    else 
      nelast = 10
    endif
  endif

  allocate(matwrk%a(grid%ngll,grid%ngll,nelast))
  !NOTE: temporary array is created for arguments #4 and #5 
  !      but does not seem to hurt (no memory leak)
  call MAT_ELAST_init_a( matwrk%a, grid%ngll, nelast &
    , SE_InverseJacobian(grid,e), SE_VolumeWeights(grid,e), matpro )

  MAT_ELAST_memwrk = MAT_ELAST_memwrk &
                   + size( transfer(matwrk, (/ 0d0 /) )) &
                   + size(matwrk%a)

  if (grid%W < huge(1d0)) then
     allocate(matwrk%beta(grid%ngll,grid%ngll))
     call MAT_ELAST_init_25D(matwrk%beta,grid%ngll &
                   , SE_VolumeWeights(grid,e),grid%W,matpro,ndof)
  endif


  end subroutine MAT_ELAST_init_elem_work

!---------------------------------------------------------------------
  subroutine MAT_ELAST_init_a(a,ngll,nelast,xjaci,weights,mat)

  integer, intent(in) :: ngll,nelast
  double precision, intent(out) :: a(ngll,ngll,nelast)
  double precision, intent(in) :: xjaci(2,2,ngll,ngll)
  double precision, intent(in) :: weights(ngll,ngll)
  type(matpro_elem_type), intent(in) :: mat

  double precision, dimension(ngll,ngll) :: DxiDx,DxiDz,DetaDx,DetaDz
  double precision, dimension(ngll,ngll) :: Kx,Kz,la,mu,mux,muz
  integer :: k

  DxiDx  = xjaci(1,1,:,:)
  DxiDz  = xjaci(1,2,:,:)
  DetaDx = xjaci(2,1,:,:)
  DetaDz = xjaci(2,2,:,:)
  
  if (MAT_isKind(mat,isElasticIsotropic)) then
    call MAT_getProp(la, mat,'lambda')
    call MAT_getProp(mu, mat,'mu')
    mux = mu
    muz = mu
    Kx = la + 2d0*mu
    Kz = Kx
  else
    call MAT_getProp(la, mat,'c13')
    call MAT_getProp(mu, mat,'c55')
    muz = mu
    call MAT_getProp(mux, mat,'c66')
    call MAT_getProp(Kx, mat,'c11')
    call MAT_getProp(Kz, mat,'c33')
  endif

  select case(nelast)

  case(2) ! SH flat
    a(:,:,1) = mux * DxiDx*DxiDx
    a(:,:,2) = muz * DetaDz*DetaDz
    
  case(3) ! SH general
    a(:,:,1) = mux * DxiDx*DxiDx   + muz * DxiDz*DxiDz 
    a(:,:,2) = mux * DetaDx*DetaDx + muz * DetaDz*DetaDz 
    a(:,:,3) = mux * DxiDx*DetaDx  + muz * DxiDz*DetaDz 

  case(6) ! P-SV flat
    a(:,:,1) = Kx * DxiDx*DxiDx
    a(:,:,2) = la * DxiDx*DetaDz
    a(:,:,3) = Kz * DetaDz*DetaDz
    a(:,:,4) = mu * DetaDz*DetaDz
    a(:,:,5) = mu * DxiDx*DetaDz
    a(:,:,6) = mu * DxiDx*DxiDx

  case(10) ! P-SV general
    a(:,:,1) = Kx * DxiDx*DxiDx   + mu * DxiDz*DxiDz
    a(:,:,2) = la * DxiDx*DetaDz  + mu * DxiDz*DetaDx
    a(:,:,3) = Kz * DetaDz*DetaDz + mu * DetaDx*DetaDx
    a(:,:,4) = Kx * DetaDx*DetaDx + mu * DetaDz*DetaDz
    a(:,:,5) = la * DxiDz*DetaDx  + mu * DxiDx*DetaDz
    a(:,:,6) = Kz * DxiDz*DxiDz   + mu * DxiDx*DxiDx
    a(:,:,7)  = Kx * DxiDx*DetaDx  + mu * DxiDz*DetaDz
    a(:,:,8)  = (la+mu) * DxiDx*DxiDz
    a(:,:,9)  = (la+mu) * DetaDx*DetaDz
    a(:,:,10) = Kz * DxiDz*DetaDz  + mu * DxiDx*DetaDx

  end select

  do k=1,nelast
    a(:,:,k) = - weights*a(:,:,k)
  enddo

  end subroutine MAT_ELAST_init_a

!==============================================
  subroutine MAT_ELAST_init_25D(beta,ngll,dvol,W,mat,ndof)

  integer, intent(in) :: ngll,ndof
  double precision, intent(out) :: beta(ngll,ngll)
  double precision, dimension(ngll,ngll), intent(in) :: dvol
  double precision, intent(in) :: W
  type(matpro_elem_type), intent(in) :: mat

  double precision, dimension(ngll,ngll) :: mu, lambda, nu

  call MAT_getProp(mu,mat,'mu')
  call MAT_getProp(lambda,mat,'lambda')

  nu(:,:) = lambda(:,:) / (lambda(:,:) + mu(:,:)) / 2.0
  if(ndof == 1) then
      beta(:,:) = dvol(:,:) * mu(:,:) * (4.D0*DATAN(1.D0)/W)**2
  else
      beta(:,:) = dvol(:,:) * mu(:,:) * (4.D0*DATAN(1.D0)*(1-nu(:,:))/W)**2
  endif

end subroutine MAT_ELAST_init_25D

!=======================================================================
!
! Computes the elastic internal forces term = -K*displ - displ/(beta*W) 
! in a SEM grid using the coefficients in elast.
! On output the result is stored in the field KD (scratched)
!
! Number of multiplications per GLL node ( = total / ngll^2*nelem) :
!                       SH              P-SV
!       flat            4*ngll +2       8*ngll +8
!       general         4*ngll +4       8*ngll +16
!
  subroutine MAT_ELAST_f(f,d,m,H,Ht,ngll,ndof)

  use constants, only : OPT_NGLL 

  integer, intent(in) :: ngll,ndof
  double precision, intent(out):: f(ngll,ngll,ndof)
  double precision, intent(in) :: d(ngll,ngll,ndof)
  double precision, dimension(ngll,ngll), intent(in) :: H,Ht
  type (matwrk_elast_type), intent(in) :: m

  integer :: nelast

  nelast = size(m%a,3)

  if ( ndof==1 ) then
    if (OPT_NGLL==ngll) then
      f(:,:,1) = ELAST_KD2_SH(d(:,:,1),m%a,nelast,H,Ht)
    else
      f(:,:,1) = ELAST_KD1_SH(d(:,:,1),m%a,nelast,ngll,H,Ht)
    endif
  else
    if (OPT_NGLL==ngll) then
      f = ELAST_KD2_PSV(d,m%a,nelast,H,Ht)
    else
      f = ELAST_KD1_PSV(d,m%a,nelast,ngll,H,Ht)
    endif
  endif

  end subroutine MAT_ELAST_f

!=======================================================================
! Compute the forces acted on "crustal plane" approximated by Lapusta and Rice
!
subroutine MAT_ELAST_add_25D_f(f,d,elast,ngll,ndof)
  integer, intent(in) :: ngll,ndof
  double precision, intent(out):: f(ngll,ngll,ndof)
  double precision, intent(in) :: d(ngll,ngll,ndof)
  type (matwrk_elast_type), intent(in) :: elast

  integer :: k

    do k=1,ndof
       f(:,:,k) = f(:,:,k) - elast%beta(:,:) * d(:,:,k)
    enddo

end subroutine MAT_ELAST_add_25D_f


!----------------------------------------------------------------

  function ELAST_KD1_PSV(displ,a,nelast,ngll,H,Ht) result(f)

  use mxmlib 

  integer, intent(in) :: ngll,nelast
  double precision, dimension(ngll,ngll,nelast), intent(in) :: a
  double precision, dimension(ngll,ngll), intent(in) :: H,Ht
  double precision, dimension(ngll,ngll,2), intent(in) :: displ

  double precision, dimension(ngll,ngll,2) :: f
  double precision, dimension(ngll,ngll) :: tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta

!-- Local gradient
    dUx_dxi  = mxm( Ht, displ(:,:,1), ngll)
    dUz_dxi  = mxm( Ht, displ(:,:,2), ngll)
    dUx_deta = mxm( displ(:,:,1), H, ngll )
    dUz_deta = mxm( displ(:,:,2), H, ngll )

!-- Elementwise forces

   if (size(a,3)==6) then

    tmp = a(:,:,1)*dUx_dxi + a(:,:,2)*dUz_deta
    f(:,:,1) = mxm( H, tmp, ngll )

    tmp = a(:,:,4)*dUx_deta + a(:,:,5)*dUz_dxi
    f(:,:,1) = f(:,:,1) + mxm( tmp, Ht, ngll )

    tmp = a(:,:,5)*dUx_deta + a(:,:,6)*dUz_dxi
    f(:,:,2) = mxm( H, tmp, ngll )

    tmp = a(:,:,2)*dUx_dxi + a(:,:,3)*dUz_deta
    f(:,:,2) = f(:,:,2) + mxm( tmp, Ht, ngll )

   else

    tmp = a(:,:,1)*dUx_dxi + a(:,:,7)*dUx_deta &
        + a(:,:,8)*dUz_dxi + a(:,:,2)*dUz_deta
    f(:,:,1) = mxm( H, tmp, ngll )

    tmp = a(:,:,7)*dUx_dxi + a(:,:,4)*dUx_deta &
        + a(:,:,5)*dUz_dxi + a(:,:,9)*dUz_deta
    f(:,:,1) = f(:,:,1) + mxm( tmp, Ht, ngll )

    tmp = a(:,:,8)*dUx_dxi + a(:,:,5)*dUx_deta &
        + a(:,:,6)*dUz_dxi + a(:,:,10)*dUz_deta
    f(:,:,2) = mxm( H, tmp, ngll )

    tmp = a(:,:,2)*dUx_dxi + a(:,:,9)*dUx_deta &
        + a(:,:,10)*dUz_dxi + a(:,:,3)*dUz_deta
    f(:,:,2) = f(:,:,2) + mxm( tmp, Ht, ngll )

   endif

  end function ELAST_KD1_PSV

!----------------------------------------------------------------

  function ELAST_KD1_SH(displ,a,nelast,ngll,H,Ht) result(f)

  use mxmlib 

  integer, intent(in) :: ngll,nelast
  double precision, dimension(ngll,ngll,nelast), intent(in) :: a
  double precision, dimension(ngll,ngll), intent(in) :: H,Ht,displ

  double precision, dimension(ngll,ngll) :: f
  double precision, dimension(ngll,ngll) :: tmp, dU_dxi, dU_deta


!-- Local gradient
    dU_dxi  = mxm( Ht, displ, ngll)
    dU_deta = mxm( displ, H, ngll)

!-- Elementwise forces

    if (size(a,3)==2) then

      tmp = a(:,:,1)*dU_dxi
      f = mxm( H, tmp, ngll)

      tmp = a(:,:,2)*dU_deta
      f = f + mxm( tmp, Ht, ngll)
    
    else

      tmp = a(:,:,1)*dU_dxi + a(:,:,3)*dU_deta
      f = mxm( H, tmp, ngll)

      tmp = a(:,:,3)*dU_dxi + a(:,:,2)*dU_deta
      f = f + mxm( tmp, Ht, ngll)

    endif

  end function ELAST_KD1_SH


!=======================================================================
! Version 2: OPT_NGLL declared statically, allows for compiler optimizations
! 
  function ELAST_KD2_PSV(displ,a,nelast,H,Ht) result(f)

  use constants, only : OPT_NGLL 

  integer, intent(in) :: nelast
  double precision, intent(in) :: a(OPT_NGLL,OPT_NGLL,nelast) &
    , displ(OPT_NGLL,OPT_NGLL,2), H(OPT_NGLL,OPT_NGLL), Ht(OPT_NGLL,OPT_NGLL)

  double precision, dimension(OPT_NGLL,OPT_NGLL,2) :: f
  double precision, dimension(OPT_NGLL,OPT_NGLL) :: tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta

!-- Local gradient
    dUx_dxi  = My_MATMUL( Ht, displ(:,:,1))
    dUz_dxi  = My_MATMUL( Ht, displ(:,:,2))
    dUx_deta = My_MATMUL( displ(:,:,1), H )
    dUz_deta = My_MATMUL( displ(:,:,2), H )

!-- Elementwise forces
    if (nelast==6) then

      tmp = a(:,:,1)*dUx_dxi + a(:,:,2)*dUz_deta
      f(:,:,1) = My_MATMUL( H, tmp )
  
      tmp = a(:,:,4)*( dUx_deta + dUz_dxi )
      f(:,:,1) = f(:,:,1) + My_MATMUL( tmp, Ht )

      tmp = a(:,:,5)*dUx_deta + a(:,:,6)*dUz_dxi
      f(:,:,2) = My_MATMUL( H, tmp )
  
      tmp = a(:,:,2)*dUx_dxi + a(:,:,3)*dUz_deta
      f(:,:,2) = f(:,:,2) + My_MATMUL( tmp, Ht )

    elseif (nelast==10) then

      tmp = a(:,:,1)*dUx_dxi + a(:,:,7)*dUx_deta &
          + a(:,:,8)*dUz_dxi + a(:,:,2)*dUz_deta
      f(:,:,1) = My_MATMUL( H, tmp )
  
      tmp = a(:,:,7)*dUx_dxi + a(:,:,4)*dUx_deta &
          + a(:,:,5)*dUz_dxi + a(:,:,9)*dUz_deta
      f(:,:,1) = f(:,:,1) + My_MATMUL( tmp, Ht )
  
      tmp = a(:,:,8)*dUx_dxi + a(:,:,5)*dUx_deta &
          + a(:,:,6)*dUz_dxi + a(:,:,10)*dUz_deta
      f(:,:,2) = My_MATMUL( H, tmp )
  
      tmp = a(:,:,2)*dUx_dxi + a(:,:,9)*dUx_deta &
          + a(:,:,10)*dUz_dxi + a(:,:,3)*dUz_deta
      f(:,:,2) = f(:,:,2) + My_MATMUL( tmp, Ht )

    endif

  end function ELAST_KD2_PSV


!----------------------------------

  function ELAST_KD2_SH(displ,a,nelast,H,Ht) result(f)

  use constants, only : OPT_NGLL 

  integer, intent(in) :: nelast
  double precision, intent(in) :: a(OPT_NGLL,OPT_NGLL,nelast) &
    , displ(OPT_NGLL,OPT_NGLL), H(OPT_NGLL,OPT_NGLL), Ht(OPT_NGLL,OPT_NGLL)

  double precision, dimension(OPT_NGLL,OPT_NGLL) :: f
  double precision, dimension(OPT_NGLL,OPT_NGLL) :: tmp, dU_dxi, dU_deta

!-- Local gradient
  dU_dxi  = My_MATMUL( Ht, displ)
  dU_deta = My_MATMUL( displ, H)

!-- Elementwise forces

  if (nelast==2) then

    tmp = a(:,:,1)*dU_dxi
    f = My_MATMUL( H, tmp )

    tmp = a(:,:,2)*dU_deta
    f = f + My_MATMUL( tmp, Ht )
    
  elseif (nelast==3) then

    tmp = a(:,:,1)*dU_dxi + a(:,:,3)*dU_deta
    f = My_MATMUL( H, tmp )

    tmp = a(:,:,3)*dU_dxi + a(:,:,2)*dU_deta
    f = f + My_MATMUL( tmp, Ht )

  endif


  end function ELAST_KD2_SH


!----------------------------------
  function My_MATMUL(A,B) result(C)

    use constants, only : OPT_NGLL 

    double precision, dimension(OPT_NGLL,OPT_NGLL), intent(in) :: A,B
    double precision, dimension(OPT_NGLL,OPT_NGLL) :: C

    double precision :: Cij
    integer :: i,j,k

    do j=1,OPT_NGLL
    do i=1,OPT_NGLL
      Cij = 0d0
      do k=1,OPT_NGLL
        Cij = Cij + A(i,k)*B(k,j)
      enddo 
      C(i,j) = Cij
    enddo
    enddo

  end function My_MATMUL


!=======================================================================
!
  subroutine MAT_ELAST_stress(s,e,m,ngll,ndof)

  use spec_grid, only : sem_grid_type

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: e(ngll,ngll,ndof+1)
  double precision, intent(out) :: s(ngll,ngll,ndof+1)
  type (matpro_elem_type), intent(in) :: m

  if (MAT_isKind(m,isElasticIsotropic)) then
    call stress_isotropic
  else
    call stress_anisotropic
  endif

  contains

!-----------------------------------------------------------------------
    subroutine stress_isotropic
  
    double precision, dimension(ngll,ngll) :: mu,lambda

    if (ndof==1) then !SH
      call MAT_getProp(mu,m,'mu')
      s(:,:,1) = 2d0*mu*e(:,:,1) ! sigma_13
      s(:,:,2) = 2d0*mu*e(:,:,2) ! sigma_23
  
    else !PSV
      call MAT_getProp(lambda,m,'lambda')
      call MAT_getProp(mu,m,'mu')
      s(:,:,1) = (lambda+2d0*mu)*e(:,:,1) + lambda*e(:,:,2) ! sigma_11
      s(:,:,2) = lambda*e(:,:,1) + (lambda+2d0*mu)*e(:,:,2) ! sigma_22
      s(:,:,3) = 2d0*mu*e(:,:,3)  ! sigma_12
    endif

    end subroutine stress_isotropic

!-----------------------------------------------------------------------
    subroutine stress_anisotropic
  
    use stdio, only : IO_abort

    double precision, dimension(ngll,ngll) :: c11,c13,c33,c55

    if (ndof==2) then
      call MAT_getProp(c11,m,'c11')
      call MAT_getProp(c13,m,'c13')
      call MAT_getProp(c33,m,'c33')
      call MAT_getProp(c55,m,'c55')
      s(:,:,1) = c11*e(:,:,1) + c13*e(:,:,2)
      s(:,:,2) = c13*e(:,:,1) + c33*e(:,:,2)
      s(:,:,3) = 2d0*c55*e(:,:,3)
    else
      call IO_abort('MAT_ELAST_stress: anisotropy requires ndof=2')
    endif

    end subroutine stress_anisotropic

  end subroutine MAT_ELAST_stress

end module mat_elastic
