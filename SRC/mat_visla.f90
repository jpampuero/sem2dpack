module mat_visla

  use prop_mat

  implicit none
  private

  ! Viscoelasticity by Liu and Archuleta (2006)
  type matwrk_visla_type
    private
    double precision :: Ms, Mp
    double precision, pointer, dimension(:) :: ws => null()
    double precision, pointer, dimension(:) :: wp => null()
    double precision, pointer, dimension(:,:,:) :: sig => null()
    double precision, pointer, dimension(:,:,:,:) :: zipt => null()
  end type matwrk_visla_type

  integer, save :: isViscoLA = 0

  ! for memory report
  integer, save :: MAT_visla_mempro = 0
  integer, save :: MAT_visla_memwrk = 0

  public :: matwrk_visla_type, MAT_isViscoLA, MAT_VISLA_read, &
                   MAT_VISLA_init_elem_prop, MAT_VISLA_init_elem_work, &
                   MAT_VISLA_stress, MAT_VISLA_module, MAT_VISLA_strain, &
                   MAT_VISLA_strain2, MAT_visla_mempro, MAT_visla_memwrk

contains
!=======================================================================

logical function MAT_isViscoLA(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isViscoLA = MAT_isKind(m,isViscoLA)
end function MAT_isViscoLA
!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MAT_VISLA
! GROUP  : MATERIALS
! PURPOSE: Set viscoelasticity properties 
!          by Liu and Archuleta (2006)
!
! SYNTAX : &MAT_VISLA cp,cs,rho,Qp,Qs,fr/
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: Qp       [dble][0d0] P quality factor                                   
! ARG: Qs       [dble][0d0] S quality factor                           
! ARG: fr       [dble][0d0] reference frequency
!
!
!
! END INPUT BLOCK

subroutine MAT_VISLA_read(input,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: cp,cs,rho,Qp, Qs, fr, mu, lambda

  NAMELIST / MAT_VISLA / cp,cs,rho,Qp, Qs, fr

  cp    = 0d0
  cs    = 0d0
  rho   = 0d0
  Qp    = 0d0
  Qs    = 0d0
  fr    = 0d0

  read(iin, MAT_VISLA, END=100)
  call MAT_setKind(input,isViscoLA)

  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)  

  mu = rho*cs*cs
  lambda  = rho*(cp*cp - 2d0*cs*cs)

  ! Elif (10/19)
  ! used for incident wave condition
  call MAT_setProp(input,'lambda',lambda)
  call MAT_setProp(input,'mu',mu)

  call MAT_setProp(input,'Qp',Qp)
  call MAT_setProp(input,'Qs',Qs)
  call MAT_setProp(input,'fr',fr)

  if (echo_input) write(iout,200)  cp,cs,rho,Qp,Qs,fr

  return


  100 call IO_abort('MAT_VISLA_read: MAT_VISLA input block not found')

  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'P wave quality factor . . . . . . . .(Qp) =',EN12.3,/5x, &
    'S wave quality factor . . . . . . . .(Qs) =',EN12.3,/5x, &
    'Reference frequency . . . . . . . . .(fr) =',EN12.3,/5x)

end subroutine MAT_VISLA_read
!=======================================================================

subroutine MAT_VISLA_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_setProp(elem,'cp',ecoord,MAT_VISLA_mempro)!
  call MAT_setProp(elem,'cs',ecoord,MAT_VISLA_mempro)!
  call MAT_setProp(elem,'Qp',ecoord,MAT_VISLA_mempro)
  call MAT_setProp(elem,'Qs',ecoord,MAT_VISLA_mempro)
  call MAT_setProp(elem,'fr',ecoord,MAT_VISLA_mempro)

  ! Elif (10/19)
  ! used for incident wave condition
  call MAT_setProp(elem,'mu',ecoord,MAT_VISLA_mempro)
  call MAT_setProp(elem,'lambda',ecoord,MAT_VISLA_mempro)


end subroutine MAT_VISLA_init_elem_prop
!=======================================================================

subroutine MAT_VISLA_init_elem_work(m,p,ngll,ndof)

  type(matwrk_visla_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ngll,ndof

  double precision :: rho,cp,cs
  double precision :: Qp,Qs,fr
  integer :: i

  if (.not. MAT_isViscoLA(p)) return

  allocate(m%sig(ngll,ngll,ndof+1))
  m%sig = 0d0

  allocate(m%zipt(ngll,ngll,ndof+1,8))
  m%zipt = 0d0

  allocate(m%ws(8))
  allocate(m%wp(8))

  call MAT_getProp(Qp,p,'Qp')
  call MAT_getProp(Qs,p,'Qs')
  call MAT_getProp(fr,p,'fr')

  call MAT_getProp(rho,p,'rho')
  call MAT_getProp(cp,p,'cp')
  call MAT_getProp(cs,p,'cs')

  call MAT_VISLA_module(Qs,fr,rho*cs*cs,m%Ms,m%ws)
  call MAT_VISLA_module(Qp,fr,rho*cp*cp,m%Mp,m%wp)

  MAT_VISLA_memwrk = MAT_VISLA_memwrk &
                + size( transfer(m, (/ 0d0 /) )) &
                + size(m%ws) &
                + size(m%wp) &
                + size(m%zipt) &
                + size(m%sig)


end subroutine MAT_VISLA_init_elem_work
!=======================================================================

subroutine MAT_VISLA_stress(m,p,ngll,ndof,s,e,de,dt)

  type (matwrk_visla_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  integer, intent(in) :: ndof, ngll
  double precision, intent(out) :: s(ngll,ngll,ndof+1)
  double precision, intent(in) :: e(ngll,ngll,ndof+1)  
  double precision, intent(in) :: de(ngll,ngll,ndof+1)
  double precision, intent(in) :: dt

  integer :: i,j
  double precision :: deps(3), depsvol,sig1,sig2,coeff

  if (ndof==2) then     ! P-SV waves
    do i=1,ngll
      do j=1,ngll
        ! modification Celine 07/11/2016 xx et zz et MAT_VISLA_strain2
        call MAT_VISLA_strain2(m%zipt(i,j,1,:),dt,m%wp,m%ws,de(i,j,1),de(i,j,2),m%Mp,m%Ms,sig1) !xx
        call MAT_VISLA_strain2(m%zipt(i,j,2,:),dt,m%wp,m%ws,de(i,j,2),de(i,j,1),m%Mp,m%Ms,sig2) !zz
        coeff = 1d0 / (m%Mp * m%Mp - (m%Mp - 2d0 * m%Ms) * (m%Mp - 2d0 * m%Ms))
        deps(1) = de(i,j,1) - coeff * (m%Mp * sig1 - (m%Mp - 2d0 * m%Ms) * sig2)
        deps(2) = de(i,j,2) - coeff * (m%Mp * sig2 - (m%Mp - 2d0 * m%Ms) * sig1)
        deps(3) = MAT_VISLA_strain(m%zipt(i,j,3,:),dt,m%ws,de(i,j,3)) !xz

        m%sig(i,j,1) = m%sig(i,j,1) + sig1 * dt
        m%sig(i,j,2) = m%sig(i,j,2) + sig2 * dt
        m%sig(i,j,3) = m%sig(i,j,3)+ m%Ms* dt* deps(3)*2d0
        ! fin modification Celine 07/11/2016 xx et zz et MAT_VISLA_strain2

        ! Write out stres-strain
!        sigeps(i,j,1) = e(i,j,3)      !xz
!        sigeps(i,j,2) = m%sig(i,j,3)  !xz
      enddo
    enddo
  else                  ! SH waves
    do i=1,ngll
      do j=1,ngll
        deps(1) = MAT_VISLA_strain(m%zipt(i,j,1,:),dt,m%ws,de(i,j,1))  !xy
        deps(2) = MAT_VISLA_strain(m%zipt(i,j,2,:),dt,m%ws,de(i,j,2))  !yz

        m%sig(i,j,1) = m%sig(i,j,1)+ m%Ms* dt* deps(1)*2d0
        m%sig(i,j,2) = m%sig(i,j,2)+ m%Ms* dt* deps(2)*2d0

        ! Write out stres-strain
!        sigeps(i,j,1) = e(i,j,2)     !yz
!        sigeps(i,j,2) = m%sig(i,j,2) !yz

      enddo
    enddo
  endif

  s = m%sig


end subroutine MAT_VISLA_stress
!=======================================================================

subroutine MAT_VISLA_module(Q,fr,Melast,M,w)

  double precision, intent(in) :: Q
  double precision, intent(in) :: fr
  double precision, intent(in) :: Melast
  double precision, intent(out) :: M
  double precision, intent(out) :: w(8)


  real :: alpha(8),beta(8),tau(8), chi, w2(8)
  integer :: i
  complex :: A, A1(8)


  alpha = (/1.66958e-2,3.81644e-2,9.84666e-3,-1.36803e-2,&
        -2.85125e-2,-5.37309e-2,-6.65035e-2,-1.33696e-1/)

  beta  = (/8.98758e-2,6.84635e-2,9.67052e-2,1.20172e-1,&
        1.30728e-1,1.38746e-1,1.40705e-1,2.14647e-1/)

  tau   = (/1.72333e-3,1.80701e-3,5.38887e-3,1.99322e-2,&
        8.49833e-2,4.09335e-1,2.05951,13.2629/)


  chi  = (3.071+ 1.433*(real(Q)**(-1.158))* &
        log(real(Q)/ 5.0))/(1.0+ 0.415*real(Q))

  do i=1,8
    w2(i)  = chi* (chi*alpha(i)+ beta(i))
  enddo

  do i=1,8
    A1(i) = 1.0+ cmplx(0.0,1.0)*2.0*(4.0*atan(1.0))*real(fr)*tau(i)
    A1(i) = w2(i)/ A1(i)
  enddo

  A = A1(1)
  do i=2,8
    A = A+ A1(i)
  enddo

  w = dble(w2)
  M = Melast/ dble(cabs(1.0-A))

end subroutine MAT_VISLA_module

!=======================================================================

function MAT_VISLA_strain(zipt,dt,w,de)  result(deps)

  double precision, intent(inout) :: zipt(:)
  double precision, intent(in) :: dt
  double precision, intent(in) :: w(8)
  double precision, intent(in) :: de
  double precision  :: deps

  integer :: i
  double precision :: zipto,tau(8)


  zipto = 0d0

  tau   = (/1.72333d-3,1.80701d-3,5.38887d-3,1.99322d-2,&
        8.49833d-2,4.09335d-1,2.05951d0,13.2629d0/)

  do i=1,8
    zipt(i)  = zipt(i)* exp(-dt/tau(i))+ w(i)* (1d0-exp(-dt/tau(i)))*de
    zipto = zipto+ zipt(i)
  enddo

  deps = de- zipto

end function MAT_VISLA_strain
!=======================================================================
! ajout Celine

subroutine MAT_VISLA_strain2(zipt,dt,wp,ws,de1,de2,Mp,Ms,sig)

  double precision, intent(inout) :: zipt(8)
  double precision, intent(in) :: dt
  double precision, intent(in) :: Mp, Ms
  double precision, intent(in) :: wp(8), ws(8)
  double precision, intent(in) :: de1,de2
  double precision, intent(out) :: sig

  integer :: i
  double precision :: tau(8),zipto

  zipto = 0d0
  tau   = (/1.72333d-3,1.80701d-3,5.38887d-3,1.99322d-2,&
        8.49833d-2,4.09335d-1,2.05951d0,13.2629d0/)

  do i=1,8
    zipt(i)  = zipt(i)* exp(-dt/tau(i)) &
              + (wp(i) * Mp - 2. * ws(i) * Ms)* (1d0-exp(-dt/tau(i)))*(de1 + de2) &
              + 2. * ws(i)* (1d0-exp(-dt/tau(i)))*de1 * Ms
    zipto = zipto + zipt(i)
  enddo

  sig = ( Mp - 2. * Ms)* (de1 + de2)  + 2. * de1 * Ms - zipto

end subroutine MAT_VISLA_strain2


end module mat_visla
