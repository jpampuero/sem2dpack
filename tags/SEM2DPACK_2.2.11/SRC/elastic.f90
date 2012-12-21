! SEM2DPACK version 2.2.11 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
    integer :: ndof
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
    type(matpro_type), pointer :: homo => null(), hete(:,:) => null()
  end type mat_elem_type

  type matpro_type  ! material properties
    double precision :: cp=0d0,cs=0d0,rho=0d0,coef(4)=0d0
  end type matpro_type

  type hete_input_type  ! heterogeneous domain, info needed to construct the model
    type(distribution_type) :: cp,cs,rho
  end type hete_input_type

  interface ELAST_inquire
    module procedure ELAST_inquire_node, ELAST_inquire_element !, &
!      ELAST_inquire_domain
  end interface ELAST_inquire

  public :: elast_type,ELAST_read,ELAST_init,ELAST_init_mass,ELAST_inquire &
           ,ELAST_cpminmax,ELAST_csminmax &
           ,ELAST_KD, ELAST_strain_stress, ELAST_inquire_domain

contains


!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MATERIAL
! PURPOSE: Define elastic material properties of a tagged domain
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
!               'XXXXX' isotropic with any 2D distribution
!               Three $DIST_XXXXX blocks: 
!               density, P-velocity, S-velocity
!
!
! END INPUT BLOCK

! Read properties of a two-dimensional
! isotropic or anisotropic linear elastic element
!
  subroutine ELAST_read(elast,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type (elast_type), intent(inout) :: elast

  integer :: i,numat,ntags,tag
  character(10)  :: mode

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
    if (echo_input) write(iout,200) tag,mode
    select case(mode)
     
      case('ISOTR') ! Isotropic material: RHO, VP and VS given 
        elast%input(tag)%kind = 0
        allocate(elast%input(tag)%homo)
        call ISOTR_read(elast%input(tag)%homo,iin)

      case('ANISO') ! anisotropic material: RHO,c11, c13, c33 et c44 given (Pa)
        elast%input(tag)%kind = 0
        allocate(elast%input(tag)%homo)
        call ANISO_read(elast%input(tag)%homo,iin)
       
      case default ! isotropic with arbitrary distribution
        elast%input(tag)%kind = 1
        allocate(elast%input(tag)%hete)
        call HETE_read(elast%input(tag)%hete,mode,iin)

!      case default 
!        call IO_abort('Improper value while reading material sets')

    end select

  enddo

  return

!---- formats
  100   format(//,' M a t e r i a l   s e t s :   2 D  e l a s t i c i t y', &
         /1x,54('='),//5x, &
         'Number of material sets . . . . . . . . . = ',i5)
  200   format(/5x, &
         'Material number . . . . . . . . . . (tag) = ',i5/5x, &
         'Type  . . . . . . . . . . . . . . .(mode) = ',a)

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
         'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
         'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',EN12.3,/5x, &
         'Poisson''s ratio . . . . . . . . . (poiss) =',EN12.3,/5x, &
         'First Lame parameter Lambda . . . .(alam) =',EN12.3,/5x, &
         'Second Lame parameter Mu. . . . . . (amu) =',EN12.3,/5x, &
         'Bulk modulus K. . . . . . . . . . .(Kvol) =',EN12.3,/5x, &
         'Young''s modulus E. . . . . . . . (young) =',EN12.3)

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
         'c11 coefficient (Pascal). . . . . . (c11) =',EN12.3,/5x, &
         'c13 coefficient (Pascal). . . . . . (c13) =',EN12.3,/5x, &
         'c33 coefficient (Pascal). . . . . . (c33) =',EN12.3,/5x, &
         'c44 coefficient (Pascal). . . . . . (c44) =',EN12.3,/5x, &
         'Mass density. . . . . . . . . . . (denst) =',EN12.3,/5x, &
         'Velocity of qP along vertical axis. . . . =',EN12.3,/5x, &
         'Velocity of qP along horizontal axis. . . =',EN12.3,/5x, &
         'Velocity of qSV along vertical axis . . . =',EN12.3,/5x, &
         'Velocity of qSV along horizontal axis . . =',EN12.3)

  end subroutine ANISO_read

!=======================================================================
! Heterogeneous block
  subroutine HETE_read(hete,mode,iin)

  use distribution_general, only : DIST_read

  type(hete_input_type), intent(out) :: hete
  character(*), intent(in) :: mode
  integer, intent(in) :: iin

  call DIST_read(hete%rho,mode,iin)
  call DIST_read(hete%cp,mode,iin)
  call DIST_read(hete%cs,mode,iin)

  end subroutine HETE_read


!=======================================================================
!
!  Define arrays a1 to a10 for the computation of elastic internal forces 
!
! mode = 1 : SH
!      = 2 : P-SV
  subroutine ELAST_init(elast,grid,mode)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights
  use stdio, only : IO_abort
  use echo, only : echo_init,iout,fmt1,fmtok
  use constants, only : OPT_HOMOG

  type(elast_type), intent(inout), target :: elast
  type(sem_grid_type), intent(in) :: grid
  integer, intent(in) :: mode

  double precision, dimension(grid%ngll,grid%ngll,4), target :: coefs
  double precision, dimension(2,2,grid%ngll,grid%ngll), target :: xjaci
  double precision, dimension(grid%ngll,grid%ngll) :: weights
  double precision, dimension(:,:), pointer :: DxiDx,DxiDz,DetaDx,DetaDz
  double precision, dimension(:,:), pointer :: Kx,Kz,la,mu
  integer :: k,e,nelast,tag,nelem_a

  if (echo_init) then
    write(iout,*) 
    write(iout,'(a)') ' M a t e r i a l   p r o p e r t i e s'
    write(iout,'(a)') ' ====================================='
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Translating input velocity model'
  endif

  allocate(elast%elem(grid%nelem))
  do e=1,grid%nelem
    tag = grid%tag(e)
    if ( tag > size(elast%input) .or. tag<1 ) &
      call IO_abort('ELAST_init: tag in mesh does not correspond to a defined material')
    if ( elast%input(tag)%kind == 0 ) then
      elast%elem(e)%homo => elast%input(tag)%homo
    else
      allocate( elast%elem(e)%hete(grid%ngll,grid%ngll) )
      call HETE_init_elem( elast%elem(e)%hete, elast%input(tag)%hete, grid%ibool(:,:,e), grid%coord)
    endif
  enddo
  if (echo_init) write(iout,fmtok)

  if (echo_init) write(iout,fmt1,advance='no') 'Defining elasticity work arrays'
  if (mode==1) then
    nelast = 3
    if (grid%flat) nelast = 2
  elseif (mode==2) then
    nelast = 10
    if (grid%flat) nelast = 6
  else
    call IO_abort('ELAST_init: mode must be 1 or 2 (SH or P-SV)')
  endif
  elast%ndof = mode
  if (OPT_HOMOG) then
    nelem_a = 1
  else
    nelem_a = grid%nelem
  endif
  allocate(elast%a(grid%ngll,grid%ngll,nelast,nelem_a))
  call storearray('elast%a',size(elast%a),idouble)
  
  do e=1,nelem_a 
   
    xjaci = SE_InverseJacobian(grid,e)
    DxiDx  => xjaci(1,1,:,:)
    DxiDz  => xjaci(1,2,:,:)
    DetaDx => xjaci(2,1,:,:)
    DetaDz => xjaci(2,2,:,:)
    
    call ELAST_inquire(elast,e,coefs=coefs)
    la => coefs(:,:,1) ! c13
    mu => coefs(:,:,2) ! c44
    Kx => coefs(:,:,3) ! c11
    Kz => coefs(:,:,4) ! c33

    select case(nelast)

    case(2) ! SH flat
      elast%a(:,:,1,e) = mu * DxiDx*DxiDx
      elast%a(:,:,2,e) = mu * DetaDz*DetaDz
      
    case(3) ! SH general
      elast%a(:,:,1,e) = mu *( DxiDx*DxiDx + DxiDz*DxiDz )
      elast%a(:,:,2,e) = mu *( DetaDx*DetaDx + DetaDz*DetaDz )
      elast%a(:,:,3,e) = mu *( DxiDx*DetaDx + DxiDz*DetaDz )

    case(6) ! P-SV flat
      elast%a(:,:,1,e) = Kx * DxiDx*DxiDx
      elast%a(:,:,2,e) = la * DxiDx*DetaDz
      elast%a(:,:,3,e) = Kz * DetaDz*DetaDz
      elast%a(:,:,4,e) = mu * DetaDz*DetaDz
      elast%a(:,:,5,e) = mu * DxiDx*DetaDz
      elast%a(:,:,6,e) = mu * DxiDx*DxiDx
  
    case(10) ! P-SV general
      elast%a(:,:,1,e) = Kx * DxiDx*DxiDx   + mu * DxiDz*DxiDz
      elast%a(:,:,2,e) = la * DxiDx*DetaDz  + mu * DxiDz*DetaDx
      elast%a(:,:,3,e) = Kz * DetaDz*DetaDz + mu * DetaDx*DetaDx
      elast%a(:,:,4,e) = Kx * DetaDx*DetaDx + mu * DetaDz*DetaDz
      elast%a(:,:,5,e) = la * DxiDz*DetaDx  + mu * DxiDx*DetaDz
      elast%a(:,:,6,e) = Kz * DxiDz*DxiDz   + mu * DxiDx*DxiDx
      elast%a(:,:,7,e)  = Kx * DxiDx*DetaDx  + mu * DxiDz*DetaDz
      elast%a(:,:,8,e)  = (la+mu) * DxiDx*DxiDz
      elast%a(:,:,9,e)  = (la+mu) * DetaDx*DetaDz
      elast%a(:,:,10,e) = Kz * DxiDz*DetaDz  + mu * DxiDx*DetaDx

    end select

    weights = SE_VolumeWeights(grid,e)
    do k=1,nelast
      elast%a(:,:,k,e) = - weights*elast%a(:,:,k,e)
    enddo

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

  ! DIST_generate work with vector input/outputs
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
  subroutine ELAST_init_mass(mass,grid,elast,ndof)

  use memory_info
  use spec_grid, only : sem_grid_type, SE_VolumeWeights

  double precision, pointer :: mass(:,:)
  type(sem_grid_type), intent(in) :: grid
  type(elast_type), intent(in) :: elast
  integer, intent(in) :: ndof
  
  double precision :: massloc(grid%ngll,grid%ngll)
  integer :: e,i,j,iglob

  allocate(mass(grid%npoin,ndof))
  call storearray('mass',size(mass),idouble)
  mass = 0.d0

  do e = 1,grid%nelem
    call ELAST_inquire(elast,e,rho=massloc)
    massloc = massloc * SE_VolumeWeights(grid,e)
    do j=1,grid%ngll
    do i=1,grid%ngll
      iglob = grid%ibool(i,j,e)
      mass(iglob,1) = mass(iglob,1) + massloc(i,j)
    enddo
    enddo
  enddo

  if (ndof==2) mass(:,2) = mass(:,1)

  end subroutine ELAST_init_mass


!=======================================================================
!
  subroutine ELAST_inquire_node(elast,i,j,e,rho,coefs,cp,cs,mu,lambda)

    type(elast_type), intent(in) :: elast
    integer, intent(in) :: e,i,j
    double precision, optional, intent(out) :: coefs(4),rho,cp,cs,mu,lambda
    
    type (matpro_type), pointer :: mato

    if ( .not.associated(elast%elem(e)%hete) ) then
      mato => elast%elem(e)%homo
    else
      mato => elast%elem(e)%hete(i,j)
    endif  

    if (present(rho))   rho   = mato%rho
    if (present(coefs)) coefs = mato%coef
    if (present(lambda)) lambda = mato%coef(1)
    if (present(mu))    mu    = mato%coef(2)
    if (present(cp))    cp    = mato%cp
    if (present(cs))    cs    = mato%cs

  end subroutine ELAST_inquire_node

!=======================================================================
!
  subroutine ELAST_inquire_element(elast,e,rho,coefs,cp,cs,mu,lambda)

    type(elast_type), intent(in) :: elast
    integer, intent(in) :: e
    double precision, optional, dimension(:,:,:), intent(out) :: coefs
    double precision, optional, dimension(:,:), intent(out) :: &
      rho,cp,cs,mu,lambda
    
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
      if (present(lambda)) lambda = mato1%coef(1)
      if (present(mu))    mu    = mato1%coef(2)
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
      if (present(lambda)) lambda = mato2%coef(1)
      if (present(mu))    mu    = mato2%coef(2)
      if (present(cp))    cp    = mato2%cp
      if (present(cs))    cs    = mato2%cs
    endif  


  end subroutine ELAST_inquire_element


!=======================================================================
!
  subroutine ELAST_inquire_domain(elast,ndof)

    type(elast_type), intent(in) :: elast
    integer, intent(out) :: ndof

    ndof = elast%ndof

  end subroutine ELAST_inquire_domain

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
  subroutine ELAST_csminmax(elast,csmin,csmax)

    type(elast_type), intent(in) :: elast
    double precision, intent(out) :: csmin,csmax

    integer :: i,e

    csmax = -huge(csmax)
    csmin = huge(csmin)

    do i=1,size(elast%input)
    if (elast%input(i)%kind == 0) then
      csmax = max(csmax, elast%input(i)%homo%cs )
      csmin = min(csmin, elast%input(i)%homo%cs )
    endif  
    enddo

    do e=1,size(elast%elem)
    if ( associated(elast%elem(e)%hete) ) then
      csmax = max(csmax,maxval(elast%elem(e)%hete%cs))
      csmin = min(csmin,minval(elast%elem(e)%hete%cs))
    endif      
    enddo

  end subroutine ELAST_csminmax


!=======================================================================
!
! Computes the elastic internal forces term = -K.displ 
! in a SEM grid using the coefficients in elast.
! On output the result is stored in the field KD (scratched)
!
! Number of multiplications per GLL node ( = total / ngll^2*nelem) :
!                       SH              P-SV
!       flat            4*ngll +2       8*ngll +8
!       general         4*ngll +4       8*ngll +16
!
  subroutine ELAST_KD(elast,grid,displ,KD)

  use spec_grid, only : sem_grid_type
  use constants, only : OPT_NGLL 

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: KD(:,:)

  if ( elast%ndof==1 ) then
    if (OPT_NGLL==grid%ngll) then
      call ELAST_KD2_SH(displ(:,1),KD(:,1),elast%a,grid%Hprime,grid%Htprime,grid%ibool &
                       ,size(elast%a,3),grid%nelem,grid%npoin)
    else
      call ELAST_KD1_SH(elast,grid,displ(:,1),KD(:,1))
    endif
  else
    if (OPT_NGLL==grid%ngll) then
      call ELAST_KD2_PSV(displ,KD,elast%a,grid%Hprime,grid%Htprime,grid%ibool &
                       ,size(elast%a,3),grid%nelem,grid%npoin)
    else
      call ELAST_KD1_PSV(elast,grid,displ,KD)
    endif
  endif
 
  end subroutine ELAST_KD

!----------------------------------------------------------------

  subroutine ELAST_KD1_PSV(elast,grid,displ,KD)

  use spec_grid, only : sem_grid_type

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: KD(:,:)

  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uxloc, Uzloc, tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  integer :: i,j,e,k

  KD = 0.d0
  do e=1,grid%nelem

!-- Elementwise storage of DISPL field
    do j=1,grid%ngll
    do i=1,grid%ngll
      k = grid%ibool(i,j,e)
      Uxloc(i,j) = displ(k,1)
      Uzloc(i,j) = displ(k,2)
    enddo
    enddo

!-- Local gradient
    dUx_dxi  = MATMUL( grid%hTprime, Uxloc )
    dUz_dxi  = MATMUL( grid%hTprime, Uzloc )
    dUx_deta = MATMUL( Uxloc, grid%hprime )
    dUz_deta = MATMUL( Uzloc, grid%hprime )

!-- Elementwise forces

   if (grid%flat) then

    tmp = elast%a(:,:,1,e)*dUx_dxi + elast%a(:,:,2,e)*dUz_deta
    Uxloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,4,e)*dUx_deta + elast%a(:,:,5,e)*dUz_dxi
    Uxloc = Uxloc + MATMUL( tmp, grid%hTprime )

    tmp = elast%a(:,:,5,e)*dUx_deta + elast%a(:,:,6,e)*dUz_dxi
    Uzloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,2,e)*dUx_dxi + elast%a(:,:,3,e)*dUz_deta
    Uzloc = Uzloc + MATMUL( tmp, grid%hTprime )

   else

    tmp = elast%a(:,:,1,e)*dUx_dxi + elast%a(:,:,7,e)*dUx_deta &
        + elast%a(:,:,8,e)*dUz_dxi + elast%a(:,:,2,e)*dUz_deta
    Uxloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,7,e)*dUx_dxi + elast%a(:,:,4,e)*dUx_deta &
        + elast%a(:,:,5,e)*dUz_dxi + elast%a(:,:,9,e)*dUz_deta
    Uxloc = Uxloc + MATMUL( tmp, grid%hTprime )

    tmp = elast%a(:,:,8,e)*dUx_dxi + elast%a(:,:,5,e)*dUx_deta &
        + elast%a(:,:,6,e)*dUz_dxi + elast%a(:,:,10,e)*dUz_deta
    Uzloc = MATMUL( grid%hprime, tmp )

    tmp = elast%a(:,:,2,e)*dUx_dxi + elast%a(:,:,9,e)*dUx_deta &
        + elast%a(:,:,10,e)*dUz_dxi + elast%a(:,:,3,e)*dUz_deta
    Uzloc = Uzloc + MATMUL( tmp, grid%hTprime )

   endif

!-- Assemble internal forces
    do j=1,grid%ngll
    do i=1,grid%ngll
      k = grid%ibool(i,j,e)
      KD(k,1) = KD(k,1) + Uxloc(i,j)
      KD(k,2) = KD(k,2) + Uzloc(i,j)
    enddo
    enddo

  enddo


  end subroutine ELAST_KD1_PSV

!----------------------------------------------------------------

  subroutine ELAST_KD1_SH(elast,grid,displ,KD)

  use spec_grid, only : sem_grid_type
  use constants, only : OPT_HOMOG

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  double precision    , intent(in) :: displ(:)
  double precision    , intent(out):: KD(:)

  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uloc, hprime,htprime, tmp, dU_dxi, dU_deta
  integer :: i,j,e,ea,k

  hprime = grid%hprime
  htprime = grid%hTprime

  KD = 0.d0
  do e=1,grid%nelem

!-- Elementwise storage of DISPL field
    do j=1,grid%ngll
    do i=1,grid%ngll
      k = grid%ibool(i,j,e)
      Uloc(i,j) = displ(k)
    enddo
    enddo

!-- Local gradient
    dU_dxi  = MATMUL( hTprime, Uloc )
    dU_deta = MATMUL( Uloc, hprime )

!-- Elementwise forces
   if (OPT_HOMOG) then
     ea = 1
   else
     ea = e
   endif

   if (grid%flat) then

    tmp = elast%a(:,:,1,ea)*dU_dxi
    Uloc = MATMUL( hprime, tmp )

    tmp = elast%a(:,:,2,ea)*dU_deta
    Uloc = Uloc + MATMUL( tmp, hTprime )
    
   else

    tmp = elast%a(:,:,1,e)*dU_dxi + elast%a(:,:,3,e)*dU_deta
    Uloc = MATMUL( hprime, tmp )

    tmp = elast%a(:,:,3,e)*dU_dxi + elast%a(:,:,2,e)*dU_deta
    Uloc = Uloc + MATMUL( tmp, hTprime )

   endif

!-- Assemble internal forces
    do j=1,grid%ngll
    do i=1,grid%ngll
      k = grid%ibool(i,j,e)
      KD(k) = KD(k) + Uloc(i,j)
    enddo
    enddo

  enddo

  end subroutine ELAST_KD1_SH


!=======================================================================
! Version 2: OPT_NGLL declared statically, allows for compiler optimizations
! 
  subroutine ELAST_KD2_PSV(displ,KD,a,H,Ht,ibool &
                          ,nelast,nelem,npoin)

  use constants, only : OPT_NGLL 

  integer, intent(in) :: nelast,nelem,npoin, ibool(OPT_NGLL,OPT_NGLL,nelem)
  double precision, intent(in) :: a(OPT_NGLL,OPT_NGLL,nelast,nelem) &
                  , displ(npoin,2), H(OPT_NGLL,OPT_NGLL), Ht(OPT_NGLL,OPT_NGLL)
  double precision, intent(out):: KD(npoin,2)

  double precision, dimension(OPT_NGLL,OPT_NGLL) :: Uxloc, Uzloc &
                  , tmp, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
  integer :: i,j,e,k

  KD = 0.d0

  do e=1,nelem

!-- Elementwise storage of DISPL field
    do j=1,OPT_NGLL
    do i=1,OPT_NGLL
      k = ibool(i,j,e)
      Uxloc(i,j) = displ(k,1)
      Uzloc(i,j) = displ(k,2)
    enddo
    enddo

!-- Local gradient
    dUx_dxi  = My_MATMUL( Ht, Uxloc )
    dUz_dxi  = My_MATMUL( Ht, Uzloc )
    dUx_deta = My_MATMUL( Uxloc, H )
    dUz_deta = My_MATMUL( Uzloc, H )

!-- Elementwise forces
    if (nelast==6) then

      tmp = a(:,:,1,e)*dUx_dxi + a(:,:,2,e)*dUz_deta
      Uxloc = My_MATMUL( H, tmp )
  
      tmp = a(:,:,4,e)*( dUx_deta + dUz_dxi )
      Uxloc = Uxloc + My_MATMUL( tmp, Ht )

      tmp = a(:,:,5,e)*dUx_deta + a(:,:,6,e)*dUz_dxi
      Uzloc = My_MATMUL( H, tmp )
  
      tmp = a(:,:,2,e)*dUx_dxi + a(:,:,3,e)*dUz_deta
      Uzloc = Uzloc + My_MATMUL( tmp, Ht )

    elseif (nelast==10) then

      tmp = a(:,:,1,e)*dUx_dxi + a(:,:,7,e)*dUx_deta &
          + a(:,:,8,e)*dUz_dxi + a(:,:,2,e)*dUz_deta
      Uxloc = My_MATMUL( H, tmp )
  
      tmp = a(:,:,7,e)*dUx_dxi + a(:,:,4,e)*dUx_deta &
          + a(:,:,5,e)*dUz_dxi + a(:,:,9,e)*dUz_deta
      Uxloc = Uxloc + My_MATMUL( tmp, Ht )
  
      tmp = a(:,:,8,e)*dUx_dxi + a(:,:,5,e)*dUx_deta &
          + a(:,:,6,e)*dUz_dxi + a(:,:,10,e)*dUz_deta
      Uzloc = My_MATMUL( H, tmp )
  
      tmp = a(:,:,2,e)*dUx_dxi + a(:,:,9,e)*dUx_deta &
          + a(:,:,10,e)*dUz_dxi + a(:,:,3,e)*dUz_deta
      Uzloc = Uzloc + My_MATMUL( tmp, Ht )

    endif

!-- Assemble internal forces
    do j=1,OPT_NGLL
    do i=1,OPT_NGLL
      k = ibool(i,j,e)
      KD(k,1) = KD(k,1) + Uxloc(i,j)
      KD(k,2) = KD(k,2) + Uzloc(i,j)
    enddo
    enddo

  enddo

  end subroutine ELAST_KD2_PSV


!----------------------------------

  subroutine ELAST_KD2_SH(displ,KD,a,H,Ht,ibool &
                          ,nelast,nelem,npoin)

  use constants, only : OPT_NGLL 

  integer, intent(in) :: nelast,nelem,npoin, ibool(OPT_NGLL,OPT_NGLL,nelem)
  double precision, intent(in) :: a(OPT_NGLL,OPT_NGLL,nelast,nelem) &
                  , displ(npoin), H(OPT_NGLL,OPT_NGLL), Ht(OPT_NGLL,OPT_NGLL)
  double precision, intent(out):: KD(npoin)

  double precision, dimension(OPT_NGLL,OPT_NGLL) :: Uloc, tmp, dU_dxi, dU_deta
  integer :: i,j,e,k

  KD = 0.d0

  do e=1,nelem

!-- Elementwise storage of DISPL field
    do j=1,OPT_NGLL
    do i=1,OPT_NGLL
      k = ibool(i,j,e)
      Uloc(i,j) = displ(k)
    enddo
    enddo

!-- Local gradient
    dU_dxi  = My_MATMUL( Ht, Uloc )
    dU_deta = My_MATMUL( Uloc, H)

!-- Elementwise forces

   if (nelast==2) then

    tmp = a(:,:,1,e)*dU_dxi
    Uloc = MATMUL( H, tmp )

    tmp = a(:,:,2,e)*dU_deta
    Uloc = Uloc + MATMUL( tmp, Ht )
    
   elseif (nelast==3) then

    tmp = a(:,:,1,e)*dU_dxi + a(:,:,3,e)*dU_deta
    Uloc = MATMUL( H, tmp )

    tmp = a(:,:,3,e)*dU_dxi + a(:,:,2,e)*dU_deta
    Uloc = Uloc + MATMUL( tmp, Ht )

   endif

    
!-- Assemble internal forces
    do j=1,OPT_NGLL
    do i=1,OPT_NGLL
      k = ibool(i,j,e)
      KD(k) = KD(k) + Uloc(i,j)
    enddo
    enddo

  enddo

  end subroutine ELAST_KD2_SH


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
  subroutine ELAST_strain_stress(elast,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  if (.not.(store_strain .or. store_stress)) return
  if ( elast%ndof==1) then
    call ELAST_strain_stress_SH(elast,grid,displ(:,1),store_strain,store_stress,dataout)
  else
    call ELAST_strain_stress_PSV(elast,grid,displ,store_strain,store_stress,dataout)
  endif

  end subroutine ELAST_strain_stress

!-----------------------------------------------------------------------
!
  subroutine ELAST_strain_stress_SH(elast,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uloc, dU_dxi, dU_deta, e13, e23
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dz, deta_dx, deta_dz
  double precision, dimension(grid%ngll,grid%ngll,4), target :: coefs
  double precision, dimension(2,2,grid%ngll,grid%ngll), target :: xjaci
  double precision, pointer :: two_mu(:,:)
  integer :: is1,i,j,e

  if (store_strain) then
    is1 = 2
  else
    is1 = 0
  endif

  do e=1,grid%nelem

!-- Elementwise storage of DISPL field
    do j=1,grid%ngll
    do i=1,grid%ngll
      Uloc(i,j) = displ( grid%ibool(i,j,e) )
    enddo
    enddo

!-- Local gradient
    dU_dxi  = MATMUL( grid%hTprime, Uloc )
    dU_deta = MATMUL( Uloc, grid%hprime )

!-- Jacobian matrix
    xjaci = SE_InverseJacobian(grid,e)
    dxi_dx  => xjaci(1,1,:,:)
    dxi_dz  => xjaci(1,2,:,:)
    deta_dx => xjaci(2,1,:,:)
    deta_dz => xjaci(2,2,:,:)

!-- Strain 
    e13 = 0.5d0*( dU_dxi*dxi_dx + dU_deta*deta_dx ) ! = dUy_dx
    e23 = 0.5d0*( dU_dxi*dxi_dz + dU_deta*deta_dz ) ! = dUy_dz

    if (store_strain) then
      dataout(:,:,e,1) = e13
      dataout(:,:,e,2) = e23
    endif

    if (store_stress) then

!-- Shear modulus
      call ELAST_inquire(elast,e,coefs=coefs)
      two_mu => coefs(:,:,2)
      two_mu = 2d0*two_mu

!-- Stress, isotropic medium, antiplane strain
      dataout(:,:,e,is1+1) = two_mu*e13     ! sigma_13
      dataout(:,:,e,is1+2) = two_mu*e23     ! sigma_23

    endif

  enddo

  end subroutine ELAST_strain_stress_SH

!-----------------------------------------------------------------------
!
  subroutine ELAST_strain_stress_PSV(elast,grid,displ,store_strain,store_stress,dataout)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian

  type (elast_type)   , intent(in) :: elast
  type (sem_grid_type), intent(in) :: grid
  logical             , intent(in) :: store_strain,store_stress
  double precision    , intent(in) :: displ(:,:)
  double precision    , intent(out):: dataout(:,:,:,:)
  
  double precision, dimension(grid%ngll,grid%ngll) :: &
    Uxloc, Uzloc, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta, &
    e11, e22, e12, s33
  double precision, dimension(:,:), pointer :: dxi_dx, dxi_dz, deta_dx, deta_dz
  double precision, dimension(grid%ngll,grid%ngll,4), target :: coefs
  double precision, dimension(2,2,grid%ngll,grid%ngll), target :: xjaci
  double precision, pointer :: two_mu(:,:), lambda(:,:)
  integer :: is1,iglob,i,j,e

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
    xjaci = SE_InverseJacobian(grid,e)
    dxi_dx  => xjaci(1,1,:,:)
    dxi_dz  => xjaci(1,2,:,:)
    deta_dx => xjaci(2,1,:,:)
    deta_dz => xjaci(2,2,:,:)

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

  end subroutine ELAST_strain_stress_PSV


end module elastic
