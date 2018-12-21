module mat_cpml

  use prop_mat
  use stdio, only : IO_abort

  implicit none
  private
  
  type matwrk_cpml_type
    private
    integer :: np
    double precision :: Ap, lambda2mu,mu
    double precision, pointer, dimension(:,:) :: DumpS1 => null()     
    double precision, pointer, dimension(:,:) :: DumpS2 => null() 
    double precision, pointer, dimension(:,:) :: DumpP1 => null()     
    double precision, pointer, dimension(:,:) :: DumpP2 => null()     
    double precision, pointer, dimension(:,:,:) :: sig => null() 
  end type matwrk_cpml_type  


  integer, save :: isCPML = 0
  integer, save :: isDownward = 0

  ! for memory report
  integer, save:: MAT_CPML_mempro = 0
  integer, save:: MAT_CPML_memwrk = 0


  public :: matwrk_cpml_type &
                   ,MAT_isCPML,MAT_CPML_mempro,MAT_CPML_memwrk, MAT_CPML_read &
                   ,MAT_CPML_init_elem_prop, MAT_CPML_init_elem_work, MAT_CPML_stress

contains

!=======================================================================
  
  logical function MAT_isCPML(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isCPML = MAT_isKind(m,isCPML)
  end function MAT_isCPML
!=======================================================================
  
  logical function MAT_isDownward(m)
  type(matpro_elem_type), intent(in) :: m
  MAT_isDownward = MAT_isKind(m,isDownward)
  end function MAT_isDownward
!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MAT_CPML
! GROUP  : MATERIALS
! PURPOSE: Set absorbing boundary layer properties 
!          (Classical-Perfectly Matched Layers)
!          by Festa and Vilotte (2005)
!
! SYNTAX : &MAT_CPML cp,cs,rho,np,Ap,dir /
!
! ARG: cp       [dble][0d0] P wave velocity
! ARG: cs       [dble][0d0] S wave velocity
! ARG: rho      [dble][0d0] density
! ARG: np       [dble][0d0] attenuation power                                   
! ARG: Ap       [dble][0d0] attenuation coefficient                             
! ARG: dir      [char]['D'] attenuation direction
!                      'D': Downward
!                      'U': Upward 
!
!
! NOTE: Attenuation is done only along the z-axis for 
! upward or downward (default) direction
! Other options could be added later  
!
! END INPUT BLOCK


  subroutine MAT_CPML_read(input,iin)

  use echo, only : echo_input, iout

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: cp,cs,rho,np,Ap
  character(1) :: dir

  NAMELIST / MAT_CPML / cp,cs,rho,np,Ap,dir

  call MAT_setKind(input,isCPML)
  cp      = 0d0
  cs      = 0d0
  rho     = 0d0
  np      = 0d0
  Ap      = 0d0
  dir     = 'D'

  read(iin, MAT_CPML, END=100)

  if (dir .ne. 'D' .AND.  dir .ne. 'U') &
  call IO_abort('MAT_CPML_read: attenuation direction dir must be set')

  if (cp .le. 0d0 .OR. cs .le. 0d0 .OR. rho .le. 0d0) &
  call IO_abort('MAT_CPML_read: velocity and/or density must be set')

  if (np .le. 0d0 .OR. Ap .le. 0d0 ) &
  call IO_abort('MAT_CPML_read: attenuation coefficients must be set')

  ! Setting properties
  call MAT_setProp(input,'cp',cp)
  call MAT_setProp(input,'cs',cs)
  call MAT_setProp(input,'rho',rho)
  call MAT_setProp(input,'np',np)
  call MAT_setProp(input,'Ap',Ap)

  if (dir == 'D')  call MAT_setKind(input,isDownward)

  write(iout,200) cp,cs,rho,np,Ap,dir

  return

  100 call IO_abort('MAT_CPML_read: MAT_CPML input block not found')
  200   format(5x, &
    'P-wave velocity . . . . . . . . . . .(cp) =',EN12.3,/5x, &
    'S-wave velocity . . . . . . . . . . .(cs) =',EN12.3,/5x, &
    'Mass density. . . . . . . . . . . . (rho) =',EN12.3,/5x, &
    'Absorbing layer attenuation . . . . .(np) =',EN12.3,/5x, &
    'Attenuation coefficient . . . . . . .(Ap) =',EN12.3,/5x, &
    'Attenuation direction . . . . . . . (dir) = ',a)


  end subroutine MAT_CPML_read
!=======================================================================

  subroutine MAT_CPML_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:) 
  
  call MAT_setProp(elem,'rho',ecoord,MAT_CPML_mempro)!
  call MAT_setProp(elem,'cp',ecoord,MAT_CPML_mempro)!
  call MAT_setProp(elem,'cs',ecoord,MAT_CPML_mempro)!
  call MAT_setProp(elem,'np',ecoord,MAT_CPML_mempro)! 
  call MAT_setProp(elem,'Ap',ecoord,MAT_CPML_mempro)! 


  end subroutine MAT_CPML_init_elem_prop
!=======================================================================

  subroutine MAT_CPML_init_elem_work(m,p,grid,ngll,ndof,e,dt,coef1,coef2)

  use spec_grid, only : sem_grid_type,SE_elem_coord,SE_VolumeWeights

  type(matwrk_cpml_type), intent(inout) :: m
  type(matpro_elem_type), intent(in) :: p
  type(sem_grid_type), intent(in) :: grid  
  integer, intent(in) :: ngll
  integer, intent(in) :: ndof  
  integer, intent(in) :: e  
  double precision, intent(in) :: dt
  double precision, dimension(ngll,ngll,ndof), intent(out) :: coef1,coef2

  integer :: i,j,np
  double precision :: w(ngll),wp(ngll),Id(ngll), xi(ngll),dz,ri
  double precision :: ecoord(2,grid%ngll,grid%ngll)
  double precision :: cp,cs,Ap,rho,dnp

  allocate(m%DumpS1(ngll,ngll))
  allocate(m%DumpS2(ngll,ngll))
  allocate(m%DumpP1(ngll,ngll))
  allocate(m%DumpP2(ngll,ngll))
  m%DumpS1 = 1d0
  m%DumpS2 = 1d0
  m%DumpP1 = 1d0
  m%DumpP2 = 1d0 

  allocate(m%sig(ngll,ngll,ndof+1))
  m%sig    = 0d0

  call MAT_getProp(rho,p,'rho')
  call MAT_getProp(cp,p,'cp')
  call MAT_getProp(cs,p,'cs')
  call MAT_getProp(Ap,p,'Ap')
  call MAT_getProp(dnp,p,'np')
  np=int(dnp)

  m%lambda2mu = rho*cp*cp
  m%mu = rho*cs*cs


  ecoord  = SE_elem_coord(grid,e)
  dz = maxval(abs(ecoord(2,:,:)))- minval(abs(ecoord(2,:,:)))
  xi = grid%xgll


  ! Downward attenuation 
  if (MAT_isDownward(p)) then
    do j=1,ngll
      ri    = 0.5d0* (1d0+ xi(1+ngll-j)) * dble(ngll-1)
      w(j)  = CPML_attenuation(ri,cs,ngll,dz,np,Ap)
      wp(j) = CPML_attenuation(ri,cp,ngll,dz,np,Ap)
    enddo
  else
  ! Upward attenuation    
    do j=1,ngll
      ri    = 0.5d0* (1d0+ xi(j)) * dble(ngll-1)
      w(j)  = CPML_attenuation(ri,cs,ngll,dz,np,Ap)
      wp(j) = CPML_attenuation(ri,cp,ngll,dz,np,Ap)
    enddo
  endif

  Id = 1d0  
  do i=1,ngll
    m%DumpS2(i,:) = Id+ 0.5d0* dt* w 
    m%DumpS2(i,:) = 1d0/ m%DumpS2(i,:)
    m%DumpS1(i,:) = (Id- 0.5d0*dt*w)* m%DumpS2(i,:)
  
    m%DumpP2(i,:) = Id+ 0.5d0* dt* wp 
    m%DumpP2(i,:) = 1d0/ m%DumpP2(i,:)
    m%DumpP1(i,:) = (Id- 0.5d0*dt*wp)* m%DumpP2(i,:)
  enddo
  
  coef1(:,:,1) = m%DumpS1
  coef2(:,:,1) = m%DumpS2
   
  if (ndof == 2) then
    coef1(:,:,2) = m%DumpP1
    coef2(:,:,2) = m%DumpP2
  endif

  ! Memory envanter
  MAT_CPML_memwrk = MAT_CPML_memwrk &
                   + size( transfer(m, (/ 0d0 /) )) &
                   + size(m%DumpS1) &
                   + size(m%DumpS2) &
                   + size(m%DumpP1) &
                   + size(m%DumpP2) &                   
                   + size(m%sig)                  

  end subroutine MAT_CPML_init_elem_work
!=======================================================================
! 

  function CPML_attenuation(dist,veloc,ngll,dz,np,Ap)  result(w)
  
  double precision, intent(in) :: dist,veloc,dz,Ap
  integer, intent(in) :: ngll, np
  double precision :: w
 
  w = 1d0/ dz
  w = w* veloc* Ap* (dist/ dble(ngll-1))**np

  end function CPML_attenuation  
!=======================================================================

  subroutine MAT_CPML_stress(m,ngll,ndof,s,de,dt)

  type (matwrk_cpml_type), intent(inout) :: m
  integer, intent(in) :: ndof, ngll
  double precision, intent(out) :: s(ngll,ngll,ndof+1)  
  double precision, intent(in) :: de(ngll,ngll,ndof+1) 
  double precision, intent(in) :: dt

  integer :: i,j

  ! P-SV waves
  if (ndof == 2) then
        do i=1,ngll
        do j=1,ngll
        m%sig(i,j,1) = m%DumpP1(i,j)*m%sig(i,j,1)+ m%DumpP2(i,j)*dt*de(i,j,1)*m%lambda2mu  
        m%sig(i,j,2) = m%DumpP1(i,j)*m%sig(i,j,2)+ m%DumpP2(i,j)*dt*de(i,j,2)*m%lambda2mu  
        m%sig(i,j,3) = m%DumpS1(i,j)*m%sig(i,j,3)+ m%DumpS2(i,j)*dt*de(i,j,3)*m%mu*2d0
        enddo
        enddo
  else
  ! SH waves
        do i=1,ngll
        do j=1,ngll
          m%sig(i,j,1) = m%DumpS1(i,j)*m%sig(i,j,1)+ m%DumpS2(i,j)*dt*de(i,j,1)*m%mu*2d0
          m%sig(i,j,2) = m%DumpS1(i,j)*m%sig(i,j,2)+ m%DumpS2(i,j)*dt*de(i,j,2)*m%mu*2d0
        enddo
        enddo
  endif 

  s = m%sig

  end subroutine MAT_CPML_stress


end module mat_cpml
