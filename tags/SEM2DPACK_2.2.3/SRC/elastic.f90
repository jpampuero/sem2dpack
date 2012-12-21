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
module elastic
! Elastic properties database
! For each element there is two possible databases:
!   1. homogeneous
!   2. heterogeneous

  use distribution_general, only: distribution_type

  implicit none
  private
  
  type elast_type
    private
    double precision, pointer :: a(:,:,:,:) => null() ! coefs for internal forces
    type(mat_elem_type), pointer :: elem(:) => null() ! element number --> material properties
    type(mat_input_type), pointer :: input(:) => null() ! material index --> input properties
  end type elast_type

  type mat_input_type
    integer :: kind ! 0= homogeneous element
                    ! 1= heterogeneous element, pointwise definition
    type(matpro_type), pointer :: homo => null()  ! elementwise
    type(hete_input_type), pointer :: hete => null() 
  end type mat_input_type

  type mat_elem_type
    type(matpro_type), pointer :: homo, hete(:,:) => null()
  end type mat_elem_type

  type matpro_type  ! material properties
    double precision :: cp=0d0,cs=0d0,rho=0d0,coef(4)=0d0
  end type matpro_type

  type hete_input_type  ! heterogeneous domain, info needed to construct the model
    type(distribution_type) :: cp,cs,rho
  end type hete_input_type

  interface ELAST_inquire
    module procedure ELAST_inquire_node, ELAST_inquire_element
  end interface ELAST_inquire

  logical, save :: ireadmodel

  public :: elast_type,ELAST_read,ELAST_init,ELAST_init_mass,ELAST_inquire &
           ,ELAST_cpminmax,ELAST_csmin,ELAST_KD1,ELAST_KD2, ELAST_strain_stress

contains


!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MATERIAL
! PURPOSE: Define material properties of a tagged domain
! SYNTAX : &MATERIAL tag, mode /
!          Followed by material data, with format depending on the mode (see
!          below)
!
! ARG: tag      [int] [none] The number assigned to a domain during mesh generation
! ARG: mode     [char*5] ['ISOTR'] Type of material and/or spatial distribution.
!               The following modes are implemented and this is their data
!               format:
!
!               'ISOTR' homogeneous isotropic elastic
!               One line, dble(3):
!               density, P-wave-velocity, S-wave-velocity
!
!               'ANISO' homogeneous anisotropic
!               One line, dble(5):
!               density, c11, c13, c33, c44
!
!               'GRADI' isotropic with constant gradient
!               Three $DIST_GRADIENT blocks: 
!               density, P-velocity, S-velocity
!
! END INPUT BLOCK

! Read properties of a two-dimensional
! isotropic or anisotropic linear elastic element
!
  subroutine ELAST_read(elast,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort
  use memory_info

  integer, intent(in) :: iin
  type (elast_type), intent(out) :: elast

  integer :: i,n,indic,numat,ntags,tag
  character(5)  :: mode

 !-----------------------------------------------------------------------

  NAMELIST / MATERIAL / tag,mode

 ! Count material sets 
  numat = 0
  ntags = 0
  rewind(iin)
  do 
    read(iin,MATERIAL,END=50) 
    ntags = max(tag,ntags)
    numat = numat+1
  enddo
  50 if (numat==0)  call IO_abort('MATERIAL parameters not found')
  if (echo_input) write(iout,100) numat

 ! Read properties for each set 
  allocate(elast%input(ntags))
  rewind(iin) 
  do i = 1,numat
    mode = 'ISOTR'
    read(iin,MATERIAL)
    if (echo_input) write(iout,200) tag
    select case(mode)
     
      case('ISOTR') ! Isotropic material: RHO, VP and VS given 
        elast%input(tag)%kind = 0
        allocate(elast%input(tag)%homo)
        call ISOTR_read(elast%input(tag)%homo,iin)

      case('ANISO') ! anisotropic material: RHO,c11, c13, c33 et c44 given (Pa)
        elast%input(tag)%kind = 0
        allocate(elast%input(tag)%homo)
        call ANISO_read(elast%input(tag)%homo,iin)
       
      case('GRADI') ! isotropic with vertical gradients
        elast%input(tag)%kind = 1
        allocate(elast%input(tag)%hete)
        call GRADI_read(elast%input(tag)%hete,iin)

      case default 
        call IO_abort('Improper value while reading material sets')

    end select

  enddo

  return

!---- formats
  100   format(//,' M a t e r i a l   s e t s :   2 D  e l a s t i c i t y', &
         /1x,54('='),//5x, &
         'Number of material sets . . . . . . . . . =',i5)
  200   format(//5x,'------------------------',/5x, &
         'Material set number . . . . . . . . . . . =',i5)

  end subroutine ELAST_read

!=======================================================================
! Homogeneous block, user input
! Isotropic :  lambda, mu, K (= lambda + 2*mu), zero
  subroutine ISOTR_read(homo,iin)
  
  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type(matpro_type), intent(out) :: homo
  integer, intent(in) :: iin

  double precision :: Kmod,Kvol,young,poiss,denst,cp,cs,mu,lambda,cs2,cp2

  read(iin ,*) denst,cp,cs
  cp2 = cp*cp
  cs2 = cs*cs
  mu   = denst*cs2
  Kmod  = denst*cp2 
  lambda  = denst*(cp2 - 2d0*cs2)

  poiss = 0.5d0*(cp2-2d0*cs2)/(cp2-cs2)
  if (poiss < 0.d0 .or. poiss > 0.5d0) call IO_abort('Poisson''s ratio out of range !')

  if(echo_input) then
    Kvol  = lambda + 2d0*mu/3.d0
    young = 2d0*mu*(1d0+poiss)
    write(iout,200) cp,cs,denst,poiss,lambda,mu,Kvol,young
  endif

  homo%coef(1) = lambda
  homo%coef(2) = mu
  homo%coef(3) = Kmod
  homo%coef(4) = Kmod
  homo%cp   = cp
  homo%cs   = cs
  homo%rho  = denst

  200   format(5x, &
         '-- Isotropic material --',/5x, &
         'P-wave velocity. . . . . . . . . . . (cp) =',1pe15.8,/5x, &
         'S-wave velocity. . . . . . . . . . . (cs) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',1pe15.8,/5x, &
         'Poisson''s ratio . . . . . . . . . (poiss) =',1pe15.8,/5x, &
         'First Lame parameter Lambda. . . . (alam) =',1pe15.8,/5x, &
         'Second Lame parameter Mu. . . . . . (amu) =',1pe15.8,/5x, &
         'Bulk modulus K . . . . . . . . . . (Kvol) =',1pe15.8,/5x, &
         'Young''s modulus E . . . . . . . . (young) =',1pe15.8)

  end subroutine ISOTR_read

!=======================================================================
! Homogeneous block, user input
! Transverse anisotropic :  c11, c13, c33, c44
!
  subroutine ANISO_read(homo,iin)

  use echo, only : echo_input, iout

  type(matpro_type), intent(out) :: homo
  integer, intent(in) :: iin

  double precision :: denst,c11,c13,c33,c44,cp,cs

  read(iin ,*) denst, c11, c13, c33, c44

 ! juste une idee des proprietes
  cp  = sqrt(c11/denst)
  cs  = sqrt(c44/denst)

  if(echo_input) write(iout,300) c11,c13,c33,c44,denst, &
  sqrt(c33/denst),sqrt(c11/denst),sqrt(c44/denst),sqrt(c44/denst)

  homo%coef(1) = c13
  homo%coef(2) = c44
  homo%coef(3) = c11
  homo%coef(4) = c33
  homo%cp     = cp
  homo%cs     = cs
  homo%rho    = denst

 300   format(5x, &
         '-- Transverse anisotropic material --',/5x, &
         'c11 coefficient (Pascal). . . . . . (c11) =',1pe15.8,/5x, &
         'c13 coefficient (Pascal). . . . . . (c13) =',1pe15.8,/5x, &
         'c33 coefficient (Pascal). . . . . . (c33) =',1pe15.8,/5x, &
         'c44 coefficient (Pascal). . . . . . (c44) =',1pe15.8,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',1pe15.8,/5x, &
         'Velocity of qP along vertical axis. . . . =',1pe15.8,/5x, &
         'Velocity of qP along horizontal axis. . . =',1pe15.8,/5x, &
         'Velocity of qSV along vertical axis . . . =',1pe15.8,/5x, &
         'Velocity of qSV along horizontal axis . . =',1pe15.8)

  end subroutine ANISO_read

!=======================================================================
! Heterogeneous block, vertical gradient user input
  subroutine GRADI_read(hete,iin)

  use distribution_general, only : DIST_read

  type(hete_input_type), intent(out) :: hete
  integer, intent(in) :: iin

  call DIST_read(hete%rho,'GRADIENT',iin)
  call DIST_read(hete%cp,'GRADIENT',iin)
  call DIST_read(hete%cs,'GRADIENT',iin)

  end subroutine GRADI_read


!=======================================================================
!
!  Define arrays a1 to a10 for the computation of elastic internal forces 
!
  subroutine ELAST_init(elast,grid)

  use memory_info
  use spec_grid
  use stdio, only : IO_abort
  use echo, only : echo_init,iout,fmt1,fmtok

  type(elast_type), intent(inout), target :: elast
  type(sem_grid_type), intent(in) :: grid

  double precision, dimension(grid%ngll,grid%ngll) :: &
    DxiDx,DxiDz,DetaDx,DetaDz, &
    DxiDx_DxiDx,DxiDx_DetaDz,DxiDz_DxiDz,DxiDz_DetaDx, &
    DetaDx_DetaDx,DetaDz_DetaDz,DxiDx_DetaDx,DxiDz_DetaDz
  double precision :: Kx1,Kz1,la1,mu1
  double precision, dimension(grid%ngll,grid%ngll) :: Kx2,Kz2,la2,mu2
  integer :: k,e,nelast,tag

  if (echo_init) write(iout,fmt1,advance='no') 'Translating input velocity model'
  allocate(elast%elem(grid%nelem))
  do e=1,grid%nelem
    tag = grid%tag(e)
    if ( elast%input(tag)%kind == 0 ) then
      elast%elem(e)%homo => elast%input(tag)%homo
    else
      allocate( elast%elem(e)%hete(grid%ngll,grid%ngll) )
      call HETE_init_elem( elast%elem(e)%hete, elast%input(tag)%hete, grid%ibool(:,:,e), grid%coord)
    endif
  enddo
  if (echo_init) write(iout,fmtok)

  if (echo_init) write(iout,fmt1,advance='no') 'Defining elasticity work arrays'
  nelast = 10
  if (grid%flat) nelast = 6
  allocate(elast%a(grid%ngll,grid%ngll,grid%nelem,nelast))
  call storearray('elast%a',size(elast%a),idouble)
  
  do e=1,grid%nelem
   
    DxiDx    = grid%xjaci(1,1,:,:,e)
    DxiDz    = grid%xjaci(1,2,:,:,e)
    DetaDx   = grid%xjaci(2,1,:,:,e)
    DetaDz   = grid%xjaci(2,2,:,:,e)
    
    DxiDx_DxiDx   = DxiDx*DxiDx
    DxiDx_DetaDz  = DxiDx*DetaDz
    DxiDz_DxiDz   = DxiDz*DxiDz
    DxiDz_DetaDx  = DxiDz*DetaDx
    DetaDx_DetaDx = DetaDx*DetaDx
    DetaDz_DetaDz = DetaDz*DetaDz

    if ( .not.associated(elast%elem(e)%hete) ) then

      la1 = elast%elem(e)%homo%coef(1) ! c13
      mu1 = elast%elem(e)%homo%coef(2) ! c44
      Kx1 = elast%elem(e)%homo%coef(3) ! c11
      Kz1 = elast%elem(e)%homo%coef(4) ! c33
  
      elast%a(:,:,e,1) = Kx1 * DxiDx_DxiDx   + mu1 * DxiDz_DxiDz
      elast%a(:,:,e,2) = la1 * DxiDx_DetaDz  + mu1 * DxiDz_DetaDx
      elast%a(:,:,e,3) = Kz1 * DetaDz_DetaDz + mu1 * DetaDx_DetaDx
      elast%a(:,:,e,4) = Kx1 * DetaDx_DetaDx + mu1 * DetaDz_DetaDz
      elast%a(:,:,e,5) = la1 * DxiDz_DetaDx  + mu1 * DxiDx_DetaDz
      elast%a(:,:,e,6) = Kz1 * DxiDz_DxiDz   + mu1 * DxiDx_DxiDx
  
      if (nelast==10) then
        DxiDx_DetaDx  = DxiDx * DetaDx
        DxiDz_DetaDz  = DxiDz * DetaDz
        elast%a(:,:,e,7)  = Kx1 * DxiDx_DetaDx  + mu1 * DxiDz_DetaDz
        elast%a(:,:,e,8)  = (la1+mu1) * DxiDx * DxiDz
        elast%a(:,:,e,9)  = (la1+mu1) * DetaDx * DetaDz
        elast%a(:,:,e,10) = Kz1 * DxiDz_DetaDz  + mu1 * DxiDx_DetaDx
      endif

    else  
      la2 = elast%elem(e)%hete(:,:)%coef(1) ! c13
      mu2 = elast%elem(e)%hete(:,:)%coef(2) ! c44
      Kx2 = elast%elem(e)%hete(:,:)%coef(3) ! c11
      Kz2 = elast%elem(e)%hete(:,:)%coef(4) ! c33
  
      elast%a(:,:,e,1) = Kx2 * DxiDx_DxiDx   + mu2 * DxiDz_DxiDz
      elast%a(:,:,e,2) = la2 * DxiDx_DetaDz  + mu2 * DxiDz_DetaDx
      elast%a(:,:,e,3) = Kz2 * DetaDz_DetaDz + mu2 * DetaDx_DetaDx
      elast%a(:,:,e,4) = Kx2 * DetaDx_DetaDx + mu2 * DetaDz_DetaDz
      elast%a(:,:,e,5) = la2 * DxiDz_DetaDx  + mu2 * DxiDx_DetaDz
      elast%a(:,:,e,6) = Kz2 * DxiDz_DxiDz   + mu2 * DxiDx_DxiDx
  
      if (nelast==10) then
        DxiDx_DetaDx  = DxiDx * DetaDx
        DxiDz_DetaDz  = DxiDz * DetaDz
        elast%a(:,:,e,7)  = Kx2 * DxiDx_DetaDx  + mu2 * DxiDz_DetaDz
        elast%a(:,:,e,8)  = (la2+mu2) * DxiDx * DxiDz
        elast%a(:,:,e,9)  = (la2+mu2) * DetaDx * DetaDz
        elast%a(:,:,e,10) = Kz2 * DxiDz_DetaDz  + mu2 * DxiDx_DetaDx
      endif
    endif  

  enddo

  do k=1,nelast
    elast%a(:,:,:,k) = - grid%weights*elast%a(:,:,:,k)
                ! weights(i,j,e) = wgll(i) * wgll(j) * dvolu(e,i,j)
  enddo
  if (echo_init) write(iout,fmtok)
   
  end subroutine ELAST_init

!=====================================================================
  subroutine HETE_init_elem(matpro,matin,ibool,gcoord)

  use distribution_general, only : DIST_generate

  type(matpro_type), intent(out) :: matpro(:,:)
  type(hete_input_type), intent(in) :: matin
  integer, intent(in) :: ibool(:,:)
  double precision, intent(in) :: gcoord(:,:)

  double precision :: estuff(size(ibool)),ecoord(2,size(ibool))
  integer :: sha(2)
  
  sha = shape(ibool)
  ecoord = gcoord(:, pack(ibool,mask=.true.) )
  call DIST_generate( estuff, ecoord, matin%cp ) 
  matpro%cp = reshape(estuff,sha)
  call DIST_generate( estuff, ecoord, matin%cs ) 
  matpro%cs = reshape(estuff,sha)
  call DIST_generate( estuff, ecoord, matin%rho ) 
  matpro%rho= reshape(estuff,sha)

  matpro%coef(1) = matpro%rho*(matpro%cp**2-2d0*matpro%cs**2)
  matpro%coef(2) = matpro%rho*matpro%cs**2
  matpro%coef(3) = matpro%rho*matpro%cp**2
  matpro%coef(4) = matpro%coef(3)
      
  end subroutine HETE_init_elem

!=====================================================================
!
! Allocate and compute the mass matrix by summing the contribution of each point
!
  subroutine ELAST_init_mass(mass,grid,elast)

  use memory_info
  use spec_grid, only : sem_grid_type

  double precision, pointer :: mass(:)
  type(sem_grid_type), intent(in) :: grid
  type(elast_type), intent(in) :: elast
  
  double precision :: massloc(grid%ngll,grid%ngll)
  integer :: e,i,j,iglob

  allocate(mass(grid%npoin))
  call storearray('mass',size(mass),idouble)
  mass = 0.d0

  do e = 1,grid%nelem
    if ( .not.associated(elast%elem(e)%hete) ) then
      massloc = grid%weights(:,:,e) * elast%elem(e)%homo%rho 
    else
      massloc = grid%weights(:,:,e) * elast%elem(e)%hete%rho 
    endif
    do j=1,grid%ngll
    do i=1,grid%ngll
      iglob = grid%ibool(i,j,e)
      mass(iglob) = mass(iglob) + massloc(i,j)
    enddo
    enddo
  enddo

  end subroutine ELAST_init_mass




!=======================================================================
!
  subroutine ELAST_inquire_node(elast,i,j,e,rho,coefs,cp,cs)

    type(elast_type), intent(in) :: elast
    integer, intent(in) :: e,i,j
    double precision, optional, intent(out) :: coefs(4),rho,cp,cs
    
    type (matpro_type), pointer :: mato

    if ( .not.associated(elast%elem(e)%hete) ) then
      mato => elast%elem(e)%homo
    else
      mato => elast%elem(e)%hete(i,j)
    endif  

    if (present(rho))   rho   = mato%rho
    if (present(coefs)) coefs = mato%coef
    if (present(cp))    cp    = mato%cp
    if (present(cs))    cs    = mato%cs

  end subroutine ELAST_inquire_node

!=======================================================================
!
  subroutine ELAST_inquire_element(elast,e,rho,coefs,cp,cs)

    type(elast_type), intent(in) :: elast
    integer, intent(in) :: e
    double precision, optional, intent(out) :: coefs(:,:,:),rho(:,:),cp(:,:),cs(:,:)
    
    type (matpro_type), pointer :: mato1,mato2(:,:)

    if ( .not.associated(elast%elem(e)%hete) ) then
      mato1 => elast%elem(e)%homo
      if (present(rho))   rho   = mato1%rho
      if (present(coefs)) then
        coefs(:,:,1) = mato1%coef(1)
        coefs(:,:,2) = mato1%coef(2)
        coefs(:,:,3) = mato1%coef(3)
        coefs(:,:,4) = mato1%coef(4)
      endif
      if (present(cp))    cp    = mato1%cp
      if (present(cs))    cs    = mato1%cs
    else
      mato2 => elast%elem(e)%hete(:,:)
      if (present(rho))   rho   = mato2%rho
      if (present(coefs)) then
        coefs(:,:,1) = mato2%coef(1)
        coefs(:,:,2) = mato2%coef(2)
        coefs(:,:,3) = mato2%coef(3)
        coefs(:,:,4) = mato2%coef(4)
      endif
      if (present(cp))    cp    = mato2%cp
      if (present(cs))    cs    = mato2%cs
    endif  


  end subroutine ELAST_inquire_element

!=======================================================================
  subroutine ELAST_cpminmax(elast,cpmin,cpmax)

    type(elast_type), intent(in) :: elast
    double precision, intent(out) :: cpmin,cpmax

    integer :: i,e

    cpmax = -huge(cpmax)
    cpmin = huge(cpmin)

    do i=1,size(elast%input)
    if (elast%input(i)%kind == 0) then
      cpmax = max(cpmax, elast%input(i)%homo%cp )
      cpmin = min(cpmin, elast%input(i)%homo%cp )
    endif  
    enddo

    do e=1,size(elast%elem)
    if ( associated(elast%elem(e)%hete) ) then
      cpmax = max(cpmax,maxval(elast%elem(e)%hete%cp))
      cpmin = min(cpmin,minval(elast%elem(e)%hete%cp))
    endif      
    enddo

  end subroutine ELAST_cpminmax


!=======================================================================
  subroutine ELAST_csmin(elast,e,csmin)

    type(elast_type), intent(in) :: elast
    integer, intent(in) :: e
    double precision, intent(out) :: csmin

    if ( associated(elast%elem(e)%hete) ) then
      csmin = minval( elast%elem(e)%hete%cs )
    else
      csmin = elast%elem(e)%homo%cs
    endif

  end subroutine ELAST_csmin


!=======================================================================
!
! Computes the elastic internal forces term = -K.displ 
! in a SEM grid using the coeffcieints in elast.
! On output the result is stored in the field KD (scratched)
!
! Version 1 : using large arrays
! 
  subroutine ELAST_KD1(elast,grid,displ,KD)

  use spec_grid, only : sem_grid_type

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: KD(:,:)

  double precision, dimension(grid%ngll,grid%ngll,grid%nelem) :: &
    Uxloc, Uzloc, tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  integer :: i,j,k,iglob

!-- Elementwise storage of DISPL field

  do k=1,grid%nelem
  do j=1,grid%ngll
  do i=1,grid%ngll
    iglob = grid%ibool(i,j,k)
    Uxloc(i,j,k) = displ(iglob,1)
    Uzloc(i,j,k) = displ(iglob,2)
  enddo
  enddo
  enddo

!-- Local gradient

  do k=1,grid%nelem
    dUx_dxi(:,:,k)  = MATMUL( grid%hTprime, Uxloc(:,:,k) )
    dUz_dxi(:,:,k)  = MATMUL( grid%hTprime, Uzloc(:,:,k) )
    dUx_deta(:,:,k) = MATMUL( Uxloc(:,:,k), grid%hprime )
    dUz_deta(:,:,k) = MATMUL( Uzloc(:,:,k), grid%hprime )
  enddo

!-- Elementwise forces

  if (grid%flat) then

    tmp = elast%a(:,:,:,1)*dUx_dxi + elast%a(:,:,:,2)*dUz_deta
    do k=1,grid%nelem
      Uxloc(:,:,k) = MATMUL( grid%hprime, tmp(:,:,k) )
    enddo

    tmp = elast%a(:,:,:,4)*dUx_deta + elast%a(:,:,:,5)*dUz_dxi
    do k=1,grid%nelem
      Uxloc(:,:,k) = Uxloc(:,:,k) + MATMUL( tmp(:,:,k), grid%hTprime )
    enddo

    tmp = elast%a(:,:,:,5)*dUx_deta + elast%a(:,:,:,6)*dUz_dxi
    do k=1,grid%nelem
      Uzloc(:,:,k) = MATMUL( grid%hprime, tmp(:,:,k) )
    enddo

    tmp = elast%a(:,:,:,2)*dUx_dxi + elast%a(:,:,:,3)*dUz_deta
    do k=1,grid%nelem
      Uzloc(:,:,k) = Uzloc(:,:,k) + MATMUL( tmp(:,:,k), grid%hTprime )
    enddo

  else

    tmp = elast%a(:,:,:,1)*dUx_dxi + elast%a(:,:,:,7)*dUx_deta &
        + elast%a(:,:,:,8)*dUz_dxi + elast%a(:,:,:,2)*dUz_deta
    do k=1,grid%nelem
      Uxloc(:,:,k) = MATMUL( grid%hprime, tmp(:,:,k) )
    enddo

    tmp = elast%a(:,:,:,7)*dUx_dxi + elast%a(:,:,:,4)*dUx_deta &
        + elast%a(:,:,:,5)*dUz_dxi + elast%a(:,:,:,9)*dUz_deta
    do k=1,grid%nelem
      Uxloc(:,:,k) = Uxloc(:,:,k) + MATMUL( tmp(:,:,k), grid%hTprime )
    enddo

    tmp = elast%a(:,:,:,8)*dUx_dxi + elast%a(:,:,:,5)*dUx_deta &
        + elast%a(:,:,:,6)*dUz_dxi + elast%a(:,:,:,10)*dUz_deta
    do k=1,grid%nelem
      Uzloc(:,:,k) = MATMUL( grid%hprime, tmp(:,:,k) )
    enddo

    tmp = elast%a(:,:,:,2)*dUx_dxi + elast%a(:,:,:,9)*dUx_deta &
        + elast%a(:,:,:,10)*dUz_dxi + elast%a(:,:,:,3)*dUz_deta
    do k=1,grid%nelem
      Uzloc(:,:,k) = Uzloc(:,:,k) + MATMUL( tmp(:,:,k), grid%hTprime )
    enddo

  endif

!-- Pointwise assemble of internal forces

  KD = 0.d0

  do k=1,grid%nelem
  do j=1,grid%ngll
  do i=1,grid%ngll
    iglob = grid%ibool(i,j,k)
    KD(iglob,1) = KD(iglob,1) + Uxloc(i,j,k)
    KD(iglob,2) = KD(iglob,2) + Uzloc(i,j,k)
  enddo
  enddo
  enddo

  end subroutine ELAST_KD1


!=======================================================================
! Version 2: element by element 
! Performs better than version 1
! 
  subroutine ELAST_KD2(elast,grid,displ,KD)

  use spec_grid, only : sem_grid_type

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: KD(:,:)

  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uxloc, Uzloc, tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  integer :: i,j,k,iglob

  KD = 0.d0
  do k=1,grid%nelem

!-- Elementwise storage of DISPL field
    do j=1,grid%ngll
    do i=1,grid%ngll
      iglob = grid%ibool(i,j,k)
      Uxloc(i,j) = displ(iglob,1)
      Uzloc(i,j) = displ(iglob,2)
    enddo
    enddo

!-- Local gradient
    dUx_dxi  = MATMUL( grid%hTprime, Uxloc )
    dUz_dxi  = MATMUL( grid%hTprime, Uzloc )
    dUx_deta = MATMUL( Uxloc, grid%hprime )
    dUz_deta = MATMUL( Uzloc, grid%hprime )

!-- Elementwise forces

   if (grid%flat) then

!-- For PML, define PHIba = PHI_a * dUb_da (convolution)
!   Assume "box" grid: xi//x and eta//z
!    call PML_UpdateStrain(PML(k),dUx_dxi,dUz_dxi,dUx_deta,dUz_deta,grid%ngll)

    tmp = elast%a(:,:,k,1)*dUx_dxi + elast%a(:,:,k,2)*dUz_deta
    Uxloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,k,4)*dUx_deta + elast%a(:,:,k,5)*dUz_dxi
    Uxloc = Uxloc + MATMUL( tmp, grid%hTprime )

    tmp = elast%a(:,:,k,5)*dUx_deta + elast%a(:,:,k,6)*dUz_dxi
    Uzloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,k,2)*dUx_dxi + elast%a(:,:,k,3)*dUz_deta
    Uzloc = Uzloc + MATMUL( tmp, grid%hTprime )

!-- For PML, define PHIfa = PHI_a * KDa (convolution)
!    call PML_UpdateStress(PML(k),Uxloc,Uzloc,grid%ngll)

   else

    tmp = elast%a(:,:,k,1)*dUx_dxi + elast%a(:,:,k,7)*dUx_deta &
        + elast%a(:,:,k,8)*dUz_dxi + elast%a(:,:,k,2)*dUz_deta
    Uxloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,k,7)*dUx_dxi + elast%a(:,:,k,4)*dUx_deta &
        + elast%a(:,:,k,5)*dUz_dxi + elast%a(:,:,k,9)*dUz_deta
    Uxloc = Uxloc + MATMUL( tmp, grid%hTprime )

    tmp = elast%a(:,:,k,8)*dUx_dxi + elast%a(:,:,k,5)*dUx_deta &
        + elast%a(:,:,k,6)*dUz_dxi + elast%a(:,:,k,10)*dUz_deta
    Uzloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,k,2)*dUx_dxi + elast%a(:,:,k,9)*dUx_deta &
        + elast%a(:,:,k,10)*dUz_dxi + elast%a(:,:,k,3)*dUz_deta
    Uzloc = Uzloc + MATMUL( tmp, grid%hTprime )

   endif

!-- Assemble internal forces
    do j=1,grid%ngll
    do i=1,grid%ngll
      iglob = grid%ibool(i,j,k)
      KD(iglob,1) = KD(iglob,1) + Uxloc(i,j)
      KD(iglob,2) = KD(iglob,2) + Uzloc(i,j)
    enddo
    enddo

  enddo


  end subroutine ELAST_KD2


!=======================================================================
!
  subroutine ELAST_strain_stress(elast,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uxloc, Uzloc, tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta, &
    dxi_dx, dxi_dz, deta_dx, deta_dz, e11, e22, e12, s33
  double precision, dimension(grid%ngll,grid%ngll,4), target :: coefs
  double precision, pointer :: two_mu(:,:), lambda(:,:)
  integer :: is1,iglob,i,j,e

  if (.not.(store_strain .or. store_stress)) return

  if (store_strain) then
    is1 = 3
  else
    is1 = 0
  endif

  do e=1,grid%nelem

!-- Elementwise storage of DISPL field
    do j=1,grid%ngll
    do i=1,grid%ngll
      iglob = grid%ibool(i,j,e)
      Uxloc(i,j) = displ(iglob,1)
      Uzloc(i,j) = displ(iglob,2)
    enddo
    enddo

!-- Local gradient
    dUx_dxi  = MATMUL( grid%hTprime, Uxloc )
    dUz_dxi  = MATMUL( grid%hTprime, Uzloc )
    dUx_deta = MATMUL( Uxloc, grid%hprime )
    dUz_deta = MATMUL( Uzloc, grid%hprime )

!-- Jacobian matrix
    dxi_dx    = grid%xjaci(1,1,:,:,e)
    dxi_dz    = grid%xjaci(1,2,:,:,e)
    deta_dx   = grid%xjaci(2,1,:,:,e)
    deta_dz   = grid%xjaci(2,2,:,:,e)

!-- Strain 
    e11 = dUx_dxi*dxi_dx + dUx_deta*deta_dx
    e22 = dUz_dxi*dxi_dz + dUz_deta*deta_dz
    e12 = 0.5d0*( dUx_dxi*dxi_dz + dUx_deta*deta_dz  &
                       + dUz_dxi*dxi_dx + dUz_deta*deta_dx  )

    if (store_strain) then
      dataout(:,:,e,1) = e11
      dataout(:,:,e,2) = e22
      dataout(:,:,e,3) = e12
    endif

    if (store_stress) then

!-- Elastic moduli
      call ELAST_inquire(elast,e,coefs=coefs)
      lambda => coefs(:,:,1)
      two_mu => coefs(:,:,2)
      two_mu = 2d0*two_mu

!-- Stress, isotropic medium, plane strain
      s33 = lambda*(e11+e22)
      dataout(:,:,e,is1+1) = s33 +two_mu*e11     ! sigma_11
      dataout(:,:,e,is1+2) = s33 +two_mu*e22     ! sigma_22
      dataout(:,:,e,is1+3) = s33                 ! sigma_33
      dataout(:,:,e,is1+4) = two_mu*e12          ! sigma_12

    endif

  enddo

  end subroutine ELAST_strain_stress


end module elastic
