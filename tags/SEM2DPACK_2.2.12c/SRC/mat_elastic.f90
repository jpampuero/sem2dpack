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
module mat_elastic
! Linear elasticity
! Isotropic or transverse anisotropic

  use prop_mat

  implicit none
  private
  
  ! defined in the first call to MAT_setKind
  integer, save :: isElastic = 0
  integer, save :: isElasticIsotropic = 0
  integer, save :: isElasticAnisotropic = 0
  integer, save :: isElasticHomogeneous = 0 ! needed for optimizations

  public :: MAT_isElastic &
          , MAT_ELAST_read &
          , MAT_ELAST_init_elem_prop, MAT_ELAST_init_work &
          , ELAST_KD, ELAST_strain_stress

contains

!=======================================================================
  logical function MAT_isElastic(m)

  type(matpro_elem_type), intent(in) :: m

  MAT_isElastic = MAT_isKind(m,isElastic)

  end function MAT_isElastic

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_ELASTIC
! GROUP  : MATERIALS
! PURPOSE: Set material properties for a linear elastic medium
! SYNTAX : &MAT_ELASTIC rho|rhoH, cp|cpH, cs|csH /   if isotropic
!          possibly followed by DIST_XXXX blocks, in the order: rho,cp,cs
!          or &MAT_ELASTIC rho, c11,c13,c33,c44 /  if anisotropic
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: cpH,csH,rhoH     [name][''] name of non uniform distribution in the
!                        DISTRIBUTIONS_2D group, to set non uniform values
! ARG: c11,c13,c33,c44  [dble][0d0] anisotropic elastic moduli 
!
! END INPUT BLOCK

  subroutine MAT_ELAST_read(input,iin)
  
  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type(matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: Kmod,poiss,rho,cp,cs,mu,lambda
  double precision :: c11,c13,c33,c44
  character(20) :: rhoH,cpH,csH

  NAMELIST / MAT_ELASTIC / rho,cp,cs, rhoH,cpH,csH, c11,c13,c33,c44

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
  c44 = 0d0

  read(iin,MAT_ELASTIC,END=100)

  if (rho<=0d0 .and. rhoH=='') call IO_abort('MAT_ELAST_read: undefined density (rho)')
  call MAT_setProp(input,'rho',rho,rhoH,iin,rhoH)

! Isotropic material: cp, cs (m/s)
! uniform or general properties
!
! Old notations:
!   coef(1) = lambda
!   coef(2) = mu
!   coef(3) = Kmod = lambda + 2*mu
!   coef(4) = Kmod
!
  if ( (cp>0d0 .or. cpH/='') .and. (cs>0d0 .or. csH/='') ) then

    call MAT_setKind(input,isElasticIsotropic)
    call MAT_setProp(input,'cp',cp,cpH,iin,cpH)
    call MAT_setProp(input,'cs',cs,csH,iin,csH)

   ! if uniform properties
    if (cp>0d0 .and. cs>0d0 .and. rho>0d0) then
      call MAT_setKind(input,isElasticHomogeneous)
      mu   = rho*cs*cs
      Kmod  = rho*cp*cp
      lambda  = rho*(cp*cp - 2d0*cs*cs)
      poiss = 0.5d0*lambda/(Kmod-mu)  !(cp*cp-2d0*cs*cs)/(cp*cp-cs*cs)
      if (poiss < 0.d0 .or. poiss > 0.5d0) call IO_abort('Poisson''s ratio out of range !')
      if(echo_input) write(iout,200) cp,cs,rho,poiss,lambda,mu, &
                     lambda+2d0*mu/3d0, 2d0*mu*(1d0+poiss)  ! Kvol,young  
      call MAT_setProp(input,'lambda',lambda)
      call MAT_setProp(input,'mu',mu)
      call MAT_setProp(input,'Kmod',Kmod)
    endif

! anisotropic material: c11, c13, c33, c44 (Pa)
! uniform properties
!
! Old notations:
!   coef(1) = c13
!   coef(2) = c44
!   coef(3) = c11
!   coef(4) = c33
!
  elseif (rho>0d0 .and. c11>0d0 .and. c13>0 .and. c33>0 .and. c44>0) then

    call MAT_setKind(input,isElasticAnisotropic)
    call MAT_setKind(input,isElasticHomogeneous)

   ! a rough idea of wave speeds for timestep setting, ABCs and plots
    cp  = sqrt(c11/rho)
    cs  = sqrt(c44/rho)
  
    if(echo_input) write(iout,300) c11,c13,c33,c44,rho, &
      sqrt(c33/rho),sqrt(c11/rho),sqrt(c44/rho),sqrt(c44/rho)
  
    call MAT_setProp(input,'cp',cp)
    call MAT_setProp(input,'cs',cs)
    call MAT_setProp(input,'c13',c13)
    call MAT_setProp(input,'c44',c44)
    call MAT_setProp(input,'c11',c11)
    call MAT_setProp(input,'c33',c33)

  else
    call IO_abort('MAT_ELAST_read: incomplete input')
  endif

  return

  100 call IO_abort('MAT_ELAST_read: MAT_ELASTIC input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Poisson''s ratio . . . . . . . . . . . . . =',EN12.3,/5x, &
    'First Lame parameter Lambda . . . . . . . =',EN12.3,/5x, &
    'Second Lame parameter Mu. . . . . . . . . =',EN12.3,/5x, &
    'Bulk modulus K. . . . . . . . . . . . . . =',EN12.3,/5x, &
    'Young''s modulus E . . . . . . . . . . . . =',EN12.3)

  300   format(5x, &
    'c11 coefficient (Pascal). . . . . . (c11) =',EN12.3,/5x, &
    'c13 coefficient (Pascal). . . . . . (c13) =',EN12.3,/5x, &
    'c33 coefficient (Pascal). . . . . . (c33) =',EN12.3,/5x, &
    'c44 coefficient (Pascal). . . . . . (c44) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Velocity of qP along vertical axis. . . . =',EN12.3,/5x, &
    'Velocity of qP along horizontal axis. . . =',EN12.3,/5x, &
    'Velocity of qSV along vertical axis . . . =',EN12.3,/5x, &
    'Velocity of qSV along horizontal axis . . =',EN12.3)

  end subroutine MAT_ELAST_read


!=======================================================================
! Initialize material properties
!
  subroutine MAT_ELAST_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  double precision, dimension(size(ecoord,2),size(ecoord,3)) :: cp,cs,rho
  
  call MAT_setProp(elem,'cp',ecoord)
  call MAT_setProp(elem,'cs',ecoord)

  if (MAT_isKind(elem,isElasticIsotropic)) then

    if (MAT_isProp('lambda',elem%input)) then
      call MAT_setProp(elem,'lambda',ecoord)
      call MAT_setProp(elem,'mu',ecoord)
      call MAT_setProp(elem,'Kmod',ecoord)

    else
      call MAT_getProp(rho,elem,'rho')
      call MAT_getProp(cp,elem,'cp')
      call MAT_getProp(cs,elem,'cs')
      call MAT_setProp(elem,'lambda',rho*(cp*cp-2d0*cs*cs))
      call MAT_setProp(elem,'mu',rho*cs*cs)
      call MAT_setProp(elem,'Kmod',rho*cp*cp)
    endif

  else
    call MAT_setProp(elem,'c13',ecoord)
    call MAT_setProp(elem,'c44',ecoord)
    call MAT_setProp(elem,'c11',ecoord)
    call MAT_setProp(elem,'c33',ecoord)
  endif
      
  end subroutine MAT_ELAST_init_elem_prop

!=======================================================================
!
!  Define arrays a1 to a10 for the computation of elastic internal forces
!
! ndof = 1 : SH
!      = 2 : P-SV
  subroutine MAT_ELAST_init_work(matwrk,mat_elem,grid,ndof)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights &
                      , SE_elem_coord, SE_isFlat, SE_single_tag
  use stdio, only : IO_abort
  use echo, only : echo_init,iout,fmt1,fmtok

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elem_type), intent(inout) :: matwrk(grid%nelem)
  type(matpro_elem_type), intent(in) :: mat_elem(grid%nelem)
  integer, intent(in) :: ndof

  integer :: e,nelast,nelem_a,nelem

  if (isElastic==0) return

  if (echo_init) &
    write(iout,fmt1,advance='no') 'Defining work arrays for elasticity'

  if (ndof<1 .or. ndof>2) &
    call IO_abort('MAT_ELAST_init_work: ndof must be 1 or 2 (SH or P-SV)')

  nelem_a = grid%nelem

 ! if box mesh ignore the crossed terms
  if (SE_isFlat(grid)) then
    if (ndof==1) then
      nelast = 2
    else 
      nelast = 6
    endif
   ! if box mesh and homogeneous medium: store only one element
    if ( MAT_isKind(mat_elem(1),isElasticHomogeneous) &
       .and. SE_single_tag(grid) )  nelem_a = 1
  else
    if (ndof==1) then
      nelast = 3
    else 
      nelast = 10
    endif
  endif

  nelem = 0
  do e=1,nelem_a
    if (.not. MAT_isElastic(mat_elem(e))) cycle
    nelem = nelem +1
    allocate(matwrk(e)%a(grid%ngll,grid%ngll,nelast))
    call MAT_ELAST_init_elem_work( matwrk(e)%a, grid%ngll, nelast &
      , SE_InverseJacobian(grid,e), SE_VolumeWeights(grid,e) &
      , mat_elem(e) )
    matwrk(e)%H => grid%hprime
    matwrk(e)%Ht => grid%hTprime
  enddo

  call storearray('matwrk%a',nelem*grid%ngll*grid%ngll*nelast,idouble)

  if (nelem_a==1) then
    do e=2,grid%nelem
      matwrk(e)%a => matwrk(1)%a
      matwrk(e)%H => matwrk(1)%H
      matwrk(e)%Ht => matwrk(1)%Ht
    enddo
  endif

  if (echo_init) write(iout,fmtok)
   
  end subroutine MAT_ELAST_init_work


!---------------------------------------------------------------------
  subroutine MAT_ELAST_init_elem_work(a,ngll,nelast,xjaci,weights,mat)

  integer, intent(in) :: ngll,nelast
  double precision, intent(out) :: a(ngll,ngll,nelast)
  double precision, intent(in) :: xjaci(2,2,ngll,ngll)
  double precision, intent(in) :: weights(ngll,ngll)
  type(matpro_elem_type), intent(in) :: mat

  double precision, dimension(ngll,ngll) :: DxiDx,DxiDz,DetaDx,DetaDz
  double precision, dimension(ngll,ngll) :: Kx,Kz,la,mu
  integer :: k

  DxiDx  = xjaci(1,1,:,:)
  DxiDz  = xjaci(1,2,:,:)
  DetaDx = xjaci(2,1,:,:)
  DetaDz = xjaci(2,2,:,:)
  
  if (MAT_isKind(mat,isElasticIsotropic)) then
    call MAT_getProp(la, mat,'lambda')
    call MAT_getProp(mu, mat,'mu')
    call MAT_getProp(Kx, mat,'Kmod')
    Kz = Kx
  else
    call MAT_getProp(la, mat,'c13')
    call MAT_getProp(mu, mat,'c44')
    call MAT_getProp(Kx, mat,'c11')
    call MAT_getProp(Kz, mat,'c33')
  endif

  select case(nelast)

  case(2) ! SH flat
    a(:,:,1) = mu * DxiDx*DxiDx
    a(:,:,2) = mu * DetaDz*DetaDz
    
  case(3) ! SH general
    a(:,:,1) = mu *( DxiDx*DxiDx + DxiDz*DxiDz )
    a(:,:,2) = mu *( DetaDx*DetaDx + DetaDz*DetaDz )
    a(:,:,3) = mu *( DxiDx*DetaDx + DxiDz*DetaDz )

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

  end subroutine MAT_ELAST_init_elem_work


!=======================================================================
!
! Computes the elastic internal forces term = -K*displ 
! in a SEM grid using the coefficients in elast.
! On output the result is stored in the field KD (scratched)
!
! Number of multiplications per GLL node ( = total / ngll^2*nelem) :
!                       SH              P-SV
!       flat            4*ngll +2       8*ngll +8
!       general         4*ngll +4       8*ngll +16
!
  subroutine ELAST_KD(f,d,m)

  use constants, only : OPT_NGLL 

  double precision, intent(out):: f(:,:,:)
  double precision, intent(in) :: d(:,:,:)
  type (matwrk_elem_type), intent(in) :: m

  integer :: ngll,ndof,nelast

  ngll = size(d,1)
  ndof = size(d,3)
  nelast = size(m%a,3)

  if ( ndof==1 ) then
    if (OPT_NGLL==ngll) then
      f(:,:,1) = ELAST_KD2_SH(d(:,:,1),m%a,nelast,m%H,m%Ht)
    else
      f(:,:,1) = ELAST_KD1_SH(d(:,:,1),m%a,nelast,ngll,m%H,m%Ht)
    endif
  else
    if (OPT_NGLL==ngll) then
      f = ELAST_KD2_PSV(d,m%a,nelast,m%H,m%Ht)
    else
      f = ELAST_KD1_PSV(d,m%a,nelast,ngll,m%H,m%Ht)
    endif
  endif

  end subroutine ELAST_KD

!----------------------------------------------------------------

  function ELAST_KD1_PSV(displ,a,nelast,ngll,H,Ht) result(f)

  integer, intent(in) :: ngll,nelast
  double precision, dimension(ngll,ngll,nelast), intent(in) :: a
  double precision, dimension(ngll,ngll), intent(in) :: H,Ht
  double precision, dimension(ngll,ngll,2), intent(in) :: displ

  double precision, dimension(ngll,ngll,2) :: f
  double precision, dimension(ngll,ngll) :: tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta

!-- Local gradient
    dUx_dxi  = MATMUL( Ht, displ(:,:,1))
    dUz_dxi  = MATMUL( Ht, displ(:,:,2))
    dUx_deta = MATMUL( displ(:,:,1), H )
    dUz_deta = MATMUL( displ(:,:,2), H )

!-- Elementwise forces

   if (size(a,3)==6) then

    tmp = a(:,:,1)*dUx_dxi + a(:,:,2)*dUz_deta
    f(:,:,1) = MATMUL( H, tmp )

    tmp = a(:,:,4)*dUx_deta + a(:,:,5)*dUz_dxi
    f(:,:,1) = f(:,:,1) + MATMUL( tmp, Ht )

    tmp = a(:,:,5)*dUx_deta + a(:,:,6)*dUz_dxi
    f(:,:,2) = MATMUL( H, tmp )

    tmp = a(:,:,2)*dUx_dxi + a(:,:,3)*dUz_deta
    f(:,:,2) = f(:,:,2) + MATMUL( tmp, Ht )

   else

    tmp = a(:,:,1)*dUx_dxi + a(:,:,7)*dUx_deta &
        + a(:,:,8)*dUz_dxi + a(:,:,2)*dUz_deta
    f(:,:,1) = MATMUL( H, tmp )

    tmp = a(:,:,7)*dUx_dxi + a(:,:,4)*dUx_deta &
        + a(:,:,5)*dUz_dxi + a(:,:,9)*dUz_deta
    f(:,:,1) = f(:,:,1) + MATMUL( tmp, Ht )

    tmp = a(:,:,8)*dUx_dxi + a(:,:,5)*dUx_deta &
        + a(:,:,6)*dUz_dxi + a(:,:,10)*dUz_deta
    f(:,:,2) = MATMUL( H, tmp )

    tmp = a(:,:,2)*dUx_dxi + a(:,:,9)*dUx_deta &
        + a(:,:,10)*dUz_dxi + a(:,:,3)*dUz_deta
    f(:,:,2) = f(:,:,2) + MATMUL( tmp, Ht )

   endif

  end function ELAST_KD1_PSV

!----------------------------------------------------------------

  function ELAST_KD1_SH(displ,a,nelast,ngll,H,Ht) result(f)

  integer, intent(in) :: ngll,nelast
  double precision, dimension(ngll,ngll,nelast), intent(in) :: a
  double precision, dimension(ngll,ngll), intent(in) :: H,Ht,displ

  double precision, dimension(ngll,ngll) :: f

  double precision, dimension(ngll,ngll) :: tmp, dU_dxi, dU_deta


!-- Local gradient
    dU_dxi  = MATMUL( Ht, displ)
    dU_deta = MATMUL( displ, H )

!-- Elementwise forces

    if (size(a,3)==2) then

      tmp = a(:,:,1)*dU_dxi
      f = MATMUL( H, tmp )

      tmp = a(:,:,2)*dU_deta
      f = f + MATMUL( tmp, Ht )
    
    else

      tmp = a(:,:,1)*dU_dxi + a(:,:,3)*dU_deta
      f = MATMUL( H, tmp )

      tmp = a(:,:,3)*dU_dxi + a(:,:,2)*dU_deta
      f = f + MATMUL( tmp, Ht )

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
  subroutine ELAST_strain_stress(mat,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type

  type (matpro_elem_type)   , intent(in) :: mat(:)
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  if (.not.(store_strain .or. store_stress)) return
  if ( size(displ,2)==1) then
    call ELAST_strain_stress_SH(mat,grid,displ(:,1),store_strain,store_stress,dataout)
  else
    call ELAST_strain_stress_PSV(mat,grid,displ,store_strain,store_stress,dataout)
  endif

  end subroutine ELAST_strain_stress

!-----------------------------------------------------------------------
!
!*** devel: make these element-wise
!           code element loop in caller

  subroutine ELAST_strain_stress_SH(mat,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type
  use fields_class, only : FIELD_strain_elem, FIELD_get_elem

  type (matpro_elem_type), intent(in) :: mat(:)
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  double precision, dimension(grid%ngll,grid%ngll,2) :: eij
  double precision, dimension(grid%ngll,grid%ngll) :: mu,dloc
  integer :: is1,e

  if (store_strain) then
    is1 = 2
  else
    is1 = 0
  endif

  do e=1,grid%nelem

    dloc = FIELD_get_elem(displ, grid%ibool(:,:,e))
    eij = FIELD_strain_elem(dloc,grid,e)

    if (store_strain) then
      dataout(:,:,e,1) = eij(:,:,1)
      dataout(:,:,e,2) = eij(:,:,2)
    endif

!-- Stress, isotropic medium, antiplane strain
    if (store_stress) then
      call MAT_getProp(mu, mat(e),'mu')
      dataout(:,:,e,is1+1) = 2d0*mu*eij(:,:,1) ! sigma_13
      dataout(:,:,e,is1+2) = 2d0*mu*eij(:,:,2) ! sigma_23
    endif

  enddo

  end subroutine ELAST_strain_stress_SH

!-----------------------------------------------------------------------
!
  subroutine ELAST_strain_stress_PSV(mat,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type
  use fields_class, only : FIELD_strain_elem, FIELD_get_elem

  type (matpro_elem_type), intent(in) :: mat(:)
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  double precision, dimension(grid%ngll,grid%ngll,3) :: eij
  double precision, dimension(grid%ngll,grid%ngll,2) :: dloc
  double precision, dimension(grid%ngll,grid%ngll) :: mu, lambda, s33
  integer :: is1,e

  if (store_strain) then
    is1 = 3
  else
    is1 = 0
  endif

  do e=1,grid%nelem

    dloc = FIELD_get_elem(displ, grid%ibool(:,:,e))
    eij = FIELD_strain_elem(dloc,grid,e)

    if (store_strain) then
      dataout(:,:,e,1) = eij(:,:,1)
      dataout(:,:,e,2) = eij(:,:,2)
      dataout(:,:,e,3) = eij(:,:,3)
    endif

!-- Stress, isotropic medium, plane strain
    if (store_stress) then
      call MAT_getProp(lambda,mat(e),'lambda')
      call MAT_getProp(mu,mat(e),'mu')
      s33 = lambda*( eij(:,:,1) + eij(:,:,2) )
      dataout(:,:,e,is1+1) = s33 +2d0*mu*eij(:,:,1) ! sigma_11
      dataout(:,:,e,is1+2) = s33 +2d0*mu*eij(:,:,2) ! sigma_22
      dataout(:,:,e,is1+3) = s33                    ! sigma_33
      dataout(:,:,e,is1+4) = 2d0*mu*eij(:,:,3)      ! sigma_12
    endif

  enddo

  end subroutine ELAST_strain_stress_PSV


end module mat_elastic
