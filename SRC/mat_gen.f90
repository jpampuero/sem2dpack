module mat_gen

! MAT_GEN: handles material properties

!! To create a new material (for instance a material called "User"):
!!  1. Create a mat_user.f90 module, following existing examples
!!     It should contain the following:
!!     . saved parameter isUser 
!!     . type matwrk_user_type : user working data type containing
!!        all the stuff that is needed often by the solver.
!!        Usually these are material properties (copy or pointer) 
!!        for faster access than lookup in the linked list of material properties,
!!        products of material properties and quadrature weights, state variables, etc
!!     . function MAT_isUser : tells if an element is of material User
!!     . subroutine MAT_USER_read : reads material properties from main input file
!!     . subroutine MAT_USER_init_elem_prop : initializes material properties for one element
!!     . subroutine MAT_USER_init_elem_work : initializes working data for one element
!!     . subroutine MAT_USER_action : actions to be performed during the solver phase
!!        in calls to MAT_Fint (see module mat_gen.f90). The most typical is a
!!        subroutine MAT_USER_stress that computes stresses from strains for one element
!!        and updates internal variables if requested
!!     . subroutine MAT_USER_export: extracts data to be written on output binary files
!!  2. Modify the prop_mat.f90 module: set MAX_NKIND larger or equal than 
!!     the total number of existing material attributes (types and sub-types)
!!  3. Modify the mat_gen.f90 module: follow comments that start with "!!" 
!!  4. Modify init.f90: if needed, pass additional arguments to MAT_init_work
!!  5. Modify Makefile.depend: add new dependencies
!!  6. Re-compile (make)

  use prop_mat
  use memory_info

!! List here all material modules:
  use mat_mass
  use mat_elastic
  use mat_kelvin_voigt
  use mat_damage
  use mat_plastic
  use mat_visco
!!  use mat_user  

  implicit none
  private

 !-- derivation and integration stuff
  type derint_type
    double precision, pointer, dimension(:,:) :: H=>null(), Ht=>null(), weights => null(), &
      dxi_dx => null(), dxi_dy => null(), deta_dx => null(), deta_dy => null()
  end type derint_type


 ! working arrays for each element
  type matwrk_elem_type
    type(derint_type)      , pointer :: derint =>null()
    type(matwrk_elast_type), pointer :: elast=>null() 
    type(matwrk_kv_type)   , pointer :: kv=>null() 
    type(matwrk_plast_type), pointer :: plast=>null()
    type(matwrk_dmg_type)  , pointer :: dmg=>null()
    type(matwrk_visco_type), pointer :: visco=>null()
!!    type(matwrk_user_type)  , pointer :: user=>null()
  end type matwrk_elem_type

  interface MAT_strain
    module procedure MAT_strain_gen,MAT_strain_m
  end interface MAT_strain

  ! for memory report
  integer, save :: MAT_DERINT_memwrk = 0

  public :: MAT_read, MAT_init_prop, MAT_init_work, MAT_write, MAT_Fint &
           ,MAT_stress_dv
  public :: matwrk_elem_type, &
            matwrk_elast_type, matwrk_plast_type, matwrk_dmg_type, matwrk_kv_type, &
            derint_type, &
            MAT_set_derint, MAT_strain, MAT_divcurl, MAT_forces, &
            MAT_diag_stiffness_init

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : MATERIAL
! PURPOSE: Define the material type of a tagged domain
! SYNTAX : &MATERIAL tag, kind /
!          followed by one or two MAT_XXXX input blocks.
!
! ARG: tag      [int] [none]    Number identifying a mesh domain 
! ARG: kind     [name(2)] ['ELAST','']  Material types:
!               'ELAST', 'DMG','PLAST', 'KV' 
!
! NOTE   : Some combinations of material kinds can be assigned to the same domain.
!          Any material type can be combined with 'KV', for instance:
!            &MATERIAL tag=1, kind='ELAST','KV' /
!            followed by a &MAT_ELAST block and a &MAT_KV block
!          sets an elastic material with Kelvin-Voigt damping. 
!
! END INPUT BLOCK

subroutine MAT_read(mat,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  integer, intent(in) :: iin
  type (matpro_input_type), pointer :: mat(:)

  integer :: i,j,k,numat,ntags,tag
  character*6 :: kind(2) 
  character*20 :: name(2) 

  NAMELIST / MATERIAL / tag,kind

 ! numat = number of MATERIAL input blocks
 ! ntags = number of materials (domains)
  numat = 0
  ntags = 0
  rewind(iin)
  do
    tag = 0
    read(iin,MATERIAL,END=50)
    ntags = max(tag,ntags)
    numat = numat+1
    if (tag<=0) call IO_abort('MAT_read: tag must be positive')
  enddo
  if (numat /= ntags) call IO_abort('MAT_read: inconsistent or missing tags')
  50 if (numat==0)  call IO_abort('MAT_read: MATERIAL block not found')
  if (echo_input) write(iout,100) numat

  allocate(mat(numat))

 ! Read properties for each input block
  rewind(iin)
  do i = 1,numat

    tag = 0
    kind(1) = 'ELAST'
    kind(2) = ''

    read(iin,MATERIAL)

    if (echo_input) then
      write(iout,200) tag
      do k=1,2
        select case(kind(k))
          case ('ELAST'); name(k) = 'Elastic'
          case ('PLAST'); name(k) = 'Plastic'
          case ('DMG')  ; name(k) = 'Damage'
          case ('KV')   ; name(k) = 'Kelvin-Voigt'
          case ('VISCO'); name(k) = 'Visco'
          !! case ('USER'); name(k) = 'User'
          case default  ; name(k) = kind(k)
        end select
      enddo
      if (kind(2)=='') then
        write(iout,210) trim(name(1))
      else
        write(iout,220) trim(name(1)), trim(name(2))
      endif
    endif

   !NOTE: in the MAT_xxx_read subroutines "mat" must be "inout"
   !      so multiple properties can be appended to the same material
   !NOTE: there is no rewind in these subroutines
    do k=1,2
      select case(kind(k))
        case ('ELAST'); call MAT_ELAST_read(mat(tag),iin)
        case ('PLAST'); call MAT_PLAST_read(mat(tag),iin)
        case ('DMG')  ; call MAT_DMG_read(mat(tag),iin)
        case ('KV')   ; call MAT_KV_read(mat(tag),iin)
        case ('VISCO'); call MAT_VISCO_read(mat(tag),iin)
        !! case ('USER'); call MAT_USER_read(mat(tag),iin)
        case ('')     ; continue
        case default  ; call IO_abort('MAT_read: unkown kind')
      end select

     ! reposition the file right after the i-th MATERIAL block
      rewind(iin)
      do j=1,i
        tag = 0
        kind(1) = 'ELAST'
        kind(2) = ''
        read(iin,MATERIAL)
      enddo

    enddo

  enddo

  return

  100   format(//,' M a t e r i a l   P r o p e r t i e s',/1x,37('='),//5x, &
         'Number of materials . . . . . . . . . . . = ',I0)

  200   format(/5x,'Material index. . . . . . . . . . . (tag) = ',I0)
  210   format(5x,'Material type . . . . . . . . . . .(kind) = ',A)
  220   format(5x,'Material type . . . . . . . . . . .(kind) = ',A,' and ',A)

end subroutine MAT_read

!=======================================================================
! 
subroutine MAT_init_prop(mat_elem,mat_input,grid)

  use spec_grid, only : sem_grid_type, SE_elem_coord
  use echo, only : echo_init, iout, fmt1, fmtok
  use stdio, only : IO_abort, IO_new_unit
  use memory_info

  type(matpro_elem_type), pointer :: mat_elem(:)
  type(matpro_input_type), intent(in), target :: mat_input(:)
  type(sem_grid_type), intent(in) :: grid

  double precision :: celem(grid%ngll,grid%ngll)
  integer :: e,tag
  integer :: ngll,cpunit,csunit,rhounit,cpunit2,csunit2,rhounit2,iol

  if (echo_init) then
    write(iout,*) 
    write(iout,'(a)') ' M a t e r i a l   p r o p e r t i e s'
    write(iout,'(a)') ' ====================================='
    write(iout,*) 
    write(iout,fmt1,advance='no') 'Translating input model'
  endif

  allocate(mat_elem(grid%nelem))

  do tag=1,size(mat_input)
    if (all(grid%tag<tag)) write(iout,*) 'WARNING: material ',tag,' is not assigned to any element'
  enddo

  do e=1,grid%nelem
    tag = grid%tag(e)
    if ( tag > size(mat_input) .or. tag<1 ) &
      call IO_abort('ELAST_init: element tag does not correspond to a material number')
    call MAT_setProp(mat_elem(e),mat_input(tag))
    call MAT_init_elem_prop(mat_elem(e), SE_elem_coord(grid,e))
  enddo
  if (echo_init) write(iout,fmtok)

 ! report memory allocations
  call storearray('matpro', size(transfer(mat_elem(1),(/0/)))*grid%nelem,iinteg)
  call storearray('matpro:mass',MAT_MASS_mempro,idouble)
  call storearray('matpro:elast',MAT_ELAST_mempro,idouble)
  call storearray('matpro:kv',MAT_KV_mempro,idouble)
  call storearray('matpro:dmg',MAT_DMG_mempro,idouble)
  call storearray('matpro:plast',MAT_PLAST_mempro,idouble)
  call storearray('matpro:visco',MAT_VISCO_mempro,idouble)
  !!call storearray('matpro:user',MAT_USER_mempro,idouble)

  ! export velocity and density model
  if (echo_init) write(iout,fmt1,advance='no') 'Exporting model'

  inquire( IOLENGTH=iol ) real(celem)

  cpunit = IO_new_unit()
  open(cpunit,file='Cp_sem2d.tab')
  cpunit2 = IO_new_unit()
  open(cpunit2,file='Cp_sem2d.dat',status='replace',access='direct',recl=iol)

  csunit = IO_new_unit()
  open(csunit,file='Cs_sem2d.tab')
  csunit2 = IO_new_unit()
  open(csunit2,file='Cs_sem2d.dat',status='replace',access='direct',recl=iol)

  rhounit = IO_new_unit()
  open(rhounit,file='Rho_sem2d.tab')
  rhounit2 = IO_new_unit()
  open(rhounit2,file='Rho_sem2d.dat',status='replace',access='direct',recl=iol)

  ngll = grid%ngll

  do e=1,grid%nelem

    call MAT_getProp(celem,mat_elem(e),'cp')
    write(cpunit,100) celem(1,1),celem(ngll,1),celem(ngll,ngll),celem(1,ngll)
    write(cpunit2,rec=e) real(celem)

    call MAT_getProp(celem,mat_elem(e),'cs')
    write(csunit,100) celem(1,1),celem(ngll,1),celem(ngll,ngll),celem(1,ngll)
    write(csunit2,rec=e) real(celem)

    call MAT_getProp(celem,mat_elem(e),'rho')
    write(rhounit,100) celem(1,1),celem(ngll,1),celem(ngll,ngll),celem(1,ngll)
    write(rhounit2,rec=e) real(celem)

  enddo

  close(cpunit)
  close(csunit)
  close(rhounit)
  close(cpunit2)
  close(csunit2)
  close(rhounit2)

  if (echo_init) write(iout,fmtok)
  
  return

  100 format(4(e11.4,1x))

end subroutine MAT_init_prop

!-----------------------------------------------------------------------
subroutine MAT_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_MASS_init_elem_prop(elem,ecoord)

  if (MAT_isElastic(elem)) call MAT_ELAST_init_elem_prop(elem,ecoord)
  if (MAT_isPlastic(elem)) call MAT_PLAST_init_elem_prop(elem,ecoord)
  if (MAT_isKelvinVoigt(elem)) call MAT_KV_init_elem_prop(elem,ecoord)
  if (MAT_isDamage(elem)) call MAT_DMG_init_elem_prop(elem,ecoord)
  if (MAT_isVisco(elem)) call MAT_VISCO_init_elem_prop(elem,ecoord)
  !!if (MAT_isUser(elem)) call MAT_USER_init_elem_prop(elem,ecoord)

end subroutine MAT_init_elem_prop

!=======================================================================
subroutine MAT_init_work(matwrk,matpro,grid,ndof,dt)

  use spec_grid, only : sem_grid_type, SE_isFlat, SE_firstElementTagged, SE_elem_coord
  use echo, only : echo_init,iout,fmt1,fmtok
  use stdio, only : IO_abort
  use memory_info

  type(sem_grid_type), intent(in) :: grid
  type(matwrk_elem_type), pointer :: matwrk(:)
  type(matpro_elem_type), intent(inout) :: matpro(grid%nelem)
  integer, intent(in) :: ndof
  double precision, intent(in) :: dt
  integer, pointer :: elist(:)
  integer :: e,e1
  logical :: flat_grid

  if (echo_init) write(iout,fmt1,advance='no') 'Defining material work arrays'

  allocate(matwrk(grid%nelem))

 ! for memory optimization in box grids with homogeneous materials:
 ! get the list of the first element of each material
  flat_grid = SE_isFlat(grid)
  elist => SE_firstElementTagged(grid)

  do e=1,grid%nelem

   ! KV is the only non-exclusive material (it can be combined with other materials)
   ! so we put it on a separate 'if'
    if (MAT_isKelvinVoigt(matpro(e))) then
      allocate(matwrk(e)%kv)
      call MAT_KV_init_elem_work(matwrk(e)%kv,matpro(e),grid%ngll,dt)
    endif

    if (MAT_isELastic(matpro(e))) then
      e1 = elist(grid%tag(e))
     ! memory optimization in box grids with homogeneous materials:
      if ( flat_grid .and. MAT_isElasticHomogeneous(matpro(e)) .and. e>e1 ) then
        matwrk(e)%elast => matwrk(e1)%elast
      else
        allocate(matwrk(e)%elast)
        call MAT_ELAST_init_elem_work(matwrk(e)%elast,matpro(e),grid,e,ndof,flat_grid)
      endif 

    elseif (MAT_isPlastic(matpro(e))) then
      if (ndof/=2) call IO_abort('MAT_init_work: plasticity requires ndof=2 (P-SV) ')
      allocate(matwrk(e)%derint)
      allocate(matwrk(e)%plast)
      call MAT_set_derint(matwrk(e)%derint,grid,e)
      call MAT_PLAST_init_elem_work(matwrk(e)%plast,matpro(e),grid%ngll,dt)
      !call MAT_PLAST_init_elem_work(matwrk(e)%plast,matpro(e),grid%ngll,dt, SE_elem_coord(grid,e))

    elseif (MAT_isDamage(matpro(e))) then
      allocate(matwrk(e)%derint)
      allocate(matwrk(e)%dmg)
      call MAT_set_derint(matwrk(e)%derint,grid,e)
      call MAT_DMG_init_elem_work(matwrk(e)%dmg,matpro(e),grid%ngll)
   
    elseif (MAT_isVisco(matpro(e))) then
      if (ndof/=2) call IO_abort('MAT_init_work: visco-elasticity requires ndof=2 (P-SV) ')
      allocate(matwrk(e)%derint)
      allocate(matwrk(e)%visco)
      call MAT_set_derint(matwrk(e)%derint,grid,e)
      call MAT_VISCO_init_elem_work(matwrk(e)%visco,matpro(e),grid%ngll)


  !!  elseif (MAT_isUser(matpro(e)))
  !!    allocate(matwrk(e)%user)
  !!    call MAT_USER_init_work(matwrk(e)%user,matpro,grid, ... )
  !! if required, pass additional arguments to MAT_init_work 
    endif

  enddo

  deallocate(elist)

  if (echo_init) write(iout,fmtok)

 ! report memory allocations
  call storearray('matwrk', size(transfer(matwrk(1),(/0/)))*grid%nelem,iinteg)
  call storearray('matwrk%elast',MAT_ELAST_memwrk,idouble)
  call storearray('matwrk%kv',MAT_KV_memwrk,idouble)
  call storearray('matwrk%derint',MAT_DERINT_memwrk,idouble)
  call storearray('matwrk%dmg',MAT_DMG_memwrk,idouble)
  call storearray('matwrk%plast',MAT_PLAST_memwrk,idouble)
  call storearray('matwrk%visco',MAT_VISCO_memwrk,idouble)
  !!call storearray('matwrk%user',MAT_USER_memwrk,idouble)

end subroutine MAT_init_work


!=======================================================================
! Compute the internal forces on a single element 
! and update internal variables
! Called by the solver
!
subroutine MAT_Fint(f,d,v,matpro,matwrk,ngll,ndof,dt,grid, E_ep,E_el,sg,sgp)

  use spec_grid, only : sem_grid_type

  integer, intent(in) :: ngll,ndof
  double precision, dimension(ngll,ngll,ndof), intent(out) :: f
  double precision, dimension(ngll,ngll,ndof), intent(inout) :: d,v
  type(matpro_elem_type), intent(in) :: matpro
  type(matwrk_elem_type), intent(inout) :: matwrk
  double precision, intent(in) :: dt
  type(sem_grid_type), intent(in) :: grid
  double precision, intent(out) :: E_ep !increment of plastic energy
  double precision, intent(out) :: E_el !total elastic energy
  double precision, intent(out) :: sg(3),sgp(3)   !stress glut

  double precision, dimension(ngll,ngll,ndof+1) :: e,s

  if (MAT_isKelvinVoigt(matpro)) call MAT_KV_add_etav(d,v,matwrk%kv,ngll,ndof)
  if (MAT_isElastic(matpro)) then
   ! elastic material has a specialized scheme
   ! that does not require intermediate computation of strain and stress
    call MAT_ELAST_f(f,d,matwrk%elast,grid%hprime,grid%hTprime,ngll,ndof)
    if (grid%W < huge(1d0)) call MAT_ELAST_add_25D_f(f,d,matwrk%elast,ngll,ndof)

    E_ep = 0d0
    E_el = 0d0
    sg = 0d0
    sgp = 0d0
  else
    e  = MAT_strain(d,matwrk,ngll,ndof)
! DEVEL: the Kelvin-Voigt term should involve only the elastic strain rate (total - plastic)
    !e = e + eta*( MAT_strain(v,matwrk,ngll,ndof) - ep_rate )
    call MAT_stress(s,e,matwrk,matpro,ngll,ndof,.true.,dt,E_ep,E_el,sg,sgp)
    f = MAT_forces(s,matwrk%derint,ngll,ndof)
  endif

end subroutine MAT_Fint

!=======================================================================
 subroutine MAT_stress(s,e,matwrk,matpro,ngll,ndof,update,dt,E_ep,E_el,sg,sgp)

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: e(ngll,ngll,ndof+1)
  double precision, intent(out) :: s(ngll,ngll,ndof+1)
  type (matwrk_elem_type) :: matwrk
  type(matpro_elem_type), intent(in) :: matpro
  logical, intent(in) :: update
  double precision, optional, intent(in) :: dt
  double precision, optional, intent(out) :: E_ep, E_el, sg(3),sgp(3)

  double precision, dimension(ngll,ngll) :: E_ep_local, E_el_local
  double precision, dimension(ngll,ngll,3) :: sg_local,sgp_local
  integer :: k

  sg_local = 0d0
  sgp_local = 0d0

  if (MAT_isElastic(matpro)) call MAT_ELAST_stress(s,e,matpro,ngll,ndof) 
  if (MAT_isPlastic(matpro)) call MAT_PLAST_stress(s,e,matwrk%plast,ngll,update, &
                                                   E_ep_local,E_el_local,sgp_local) 
  if (MAT_isDamage(matpro)) call MAT_DMG_stress(s,e,matwrk%dmg,ngll,update,dt, &
                                                E_ep_local,E_el_local,sg_local,sgp_local)
  if (MAT_isVisco(matpro)) call MAT_VISCO_stress(s,e,matwrk%visco,ngll,dt)
!!  if (MAT_isUser(matpro)) call MAT_USER_stress(s,e,matwrk,ngll,update,...)

  if (present(E_ep)) E_ep = sum( matwrk%derint%weights * E_ep_local )
  if (present(E_el)) E_el = sum( matwrk%derint%weights * E_el_local )
  if (present(sg)) then
    do k =1,3
      sg(k) = sum( matwrk%derint%weights * sg_local(:,:,k) )
      sgp(k) = sum( matwrk%derint%weights * sgp_local(:,:,k) )
    enddo
  endif

end subroutine MAT_stress

!=======================================================================
! write a snapshot of internal variables
subroutine MAT_write(matwrk,matpro,tag)

  type(matwrk_elem_type), intent(in) :: matwrk(:)
  type(matpro_elem_type), intent(in) :: matpro(:)
  character(*), intent(in) :: tag

  call MAT_PLAST_write(matwrk,matpro,tag)
  call MAT_DMG_write(matwrk,matpro,tag)
  !!call MAT_USER_write(matwrk,matpro,tag)

end subroutine MAT_write

!-----------------------------------------------------------------------
! Export snapshots of damage variables to binary file pla_xxx_sem2d.dat
! (where xxx is the snapshot index):
!   
! The ASCII file pla_elems_sem2d.tab contains the list of elements concerned,
! in two columns:
!   col 1 = element index in the pla_xxx_sem2d.dat file
!   col 2 = element index in the global mesh
!
  subroutine MAT_PLAST_write(mw,mp,tag)

  use stdio, only : IO_new_unit

  type (matwrk_elem_type), intent(in) :: mw(:)
  type (matpro_elem_type), intent(in) :: mp(:)
  character(*), intent(in) :: tag
  character(30) :: fname

  integer :: e,k,nelem, iout,iol=0
  logical :: initialized = .false.

  if (.not. MAT_anyPlastic() ) return

  nelem = size(mp)

 ! initialize
  if (.not.initialized) then 
    iout = IO_new_unit()
    open(unit=iout,file='pla_elems_sem2d.tab')
    k=0
    do e = 1,nelem
      if ( .not. MAT_isPlastic(mp(e)) ) cycle
      k=k+1
      write(iout,*) k,e
      if (k==1) iol = e
    enddo
    close(iout)
    inquire( IOLENGTH=iol )  MAT_PLAST_export(mw(iol)%plast)
    initialized = .true.
  endif

 ! export snapshots
  write(fname,'("pla_",A,"_sem2d.dat")') trim(tag)
  iout = IO_new_unit()
  open(unit=iout,file=trim(fname),status='replace',access='direct',recl=iol)
  k=0
  do e = 1,nelem
    if ( .not. MAT_isPlastic(mp(e)) ) cycle
    k=k+1
    write(iout,rec=k) MAT_PLAST_export(mw(e)%plast)
  enddo
  close(iout)

  end subroutine MAT_PLAST_write

!-----------------------------------------------------------------------
! Export snapshots of damage variables to binary file dmg_xxx_sem2d.dat
! (where xxx is the snapshot index):
!   
! The ASCII file dmg_elems_sem2d.tab contains the list of elements concerned,
! in two columns:
!   col 1 = element index in the dmg_xxx_sem2d.dat file
!   col 2 = element index in the global mesh
!
  subroutine MAT_DMG_write(mw,mp,tag)

  use stdio, only : IO_new_unit

  type (matwrk_elem_type), intent(in) :: mw(:)
  type (matpro_elem_type), intent(in) :: mp(:)
  character(*), intent(in) :: tag
  character(30) :: fname

  integer :: e,k,nelem, iout,iol=0
  logical :: initialized = .false.

  if (.not. MAT_anyDamage() ) return

  nelem = size(mp)

 ! initialize
  if (.not.initialized) then 
    iout = IO_new_unit()
    open(unit=iout,file='dmg_elems_sem2d.tab')
    k=0
    do e = 1,nelem
      if ( .not. MAT_isDamage(mp(e)) ) cycle
      k=k+1
      write(iout,*) k,e
      if (k==1) iol = e
    enddo
    close(iout)
    inquire( IOLENGTH=iol ) MAT_DMG_export(mw(iol)%dmg)
    initialized = .true.
  endif

 ! export snapshots
  write(fname,'("dmg_",A,"_sem2d.dat")') trim(tag)
  iout = IO_new_unit()
  open(unit=iout,file=trim(fname),status='replace',access='direct',recl=iol)
  k=0
  do e = 1,nelem
    if ( .not. MAT_isDamage(mp(e)) ) cycle
    k=k+1
    write(iout,rec=k) MAT_DMG_export(mw(e)%dmg)
  enddo
  close(iout)

  end subroutine MAT_DMG_write


!=======================================================================
subroutine MAT_stress_dv(s,d,v,matwrk,matpro,grid,e,ngll,ndof)

  use spec_grid, only : sem_grid_type

  integer, intent(in) :: ngll,ndof,e
  double precision, dimension(ngll,ngll,ndof+1), intent(out) :: s
  double precision, dimension(ngll,ngll,ndof), intent(inout) :: d,v
  type(matwrk_elem_type) :: matwrk
  type(matpro_elem_type), intent(in) :: matpro
  type(sem_grid_type), intent(in) :: grid
  
  if (MAT_isKelvinVoigt(matpro)) call MAT_KV_add_etav(d,v,matwrk%kv,ngll,ndof)
  call MAT_stress(s, MAT_strain(d,matwrk,grid,e,ngll,ndof) &
                 ,matwrk,matpro,ngll,ndof,update=.false.)
  
end subroutine MAT_stress_dv

!=======================================================================
! coefficients for derivation and integration
  subroutine MAT_set_derint(d,g,e)

  use spec_grid, only : sem_grid_type, SE_InverseJacobian, SE_VolumeWeights

  type(derint_type), intent(inout) :: d
  type(sem_grid_type), intent(in) :: g
  integer, intent(in) :: e

  integer :: n
  double precision :: xjaci(2,2,g%ngll,g%ngll)
  
  n = g%ngll

  allocate(d%dxi_dx(n,n))
  allocate(d%dxi_dy(n,n))
  allocate(d%deta_dx(n,n))
  allocate(d%deta_dy(n,n))
  allocate(d%weights(n,n))

 ! coefficients to compute derivatives
  d%H => g%hprime
  d%Ht => g%hTprime
  xjaci = SE_InverseJacobian(g,e)
  d%dxi_dx  = xjaci(1,1,:,:)
  d%dxi_dy  = xjaci(1,2,:,:)
  d%deta_dx = xjaci(2,1,:,:)
  d%deta_dy = xjaci(2,2,:,:)

 ! integration weights (quadrature weights and jacobian)
  d%weights = SE_VolumeWeights(g,e)
  
 ! for memory usage report
  MAT_DERINT_memwrk = MAT_DERINT_memwrk &
                   + size( transfer(d, (/ 0d0 /) )) &
                   + size(d%weights)*5

  end subroutine MAT_set_derint

  
!=======================================================================
  function MAT_strain_gen(d,m,g,enum,ngll,ndof) result(e)

  use spec_grid, only : sem_grid_type
  use fields_class, only : FIELD_strain_elem

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: d(ngll,ngll,ndof)
  type(matwrk_elem_type), intent(in) :: m
  type(sem_grid_type), intent(in) :: g
  integer, intent(in) :: enum
  double precision :: e(ngll,ngll,ndof+1)
  
  if (associated(m%derint)) then
    e = MAT_strain_m(d,m,ngll,ndof)
  else
    e = FIELD_strain_elem(d,ngll,ndof,g,enum)
  endif
  
  end function MAT_strain_gen

!-----------------------------------------------------------------------
  function MAT_strain_m(d,m,ngll,ndof) result(e)
  
  use stdio, only : IO_abort

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: d(ngll,ngll,ndof)
  type(matwrk_elem_type), intent(in) :: m
  double precision :: e(ngll,ngll,ndof+1)

  if (.not.associated(m%derint)) &
    call IO_abort('mat_common:MAT_strain: undefined jacobian in argument m')

  if (ndof==1) then
    e = MAT_strain_SH(d,m%derint,ngll)
  elseif(ndof==2) then
    e = MAT_strain_PSV(d,m%derint,ngll)
  else
    call IO_abort('mat_common:MAT_strain: illegal value for argument ndof')
  endif

  end function MAT_strain_m

!-----------------------------------------------------------------------
  function MAT_strain_SH(d,m,ngll) result(e)

  use mxmlib

  integer, intent(in) :: ngll
  double precision, intent(in) :: d(ngll,ngll,1)
  type(derint_type), intent(in) :: m
  double precision :: e(ngll,ngll,2)

  double precision, dimension(ngll,ngll) :: dUz_dxi,dUz_deta

 ! local gradient
  dUz_dxi  = mxm( m%Ht, d(:,:,1), ngll )
  dUz_deta = mxm( d(:,:,1), m%H, ngll )

 ! strain 
  e(:,:,1) = 0.5d0*( dUz_dxi*m%dxi_dx + dUz_deta*m%deta_dx ) ! e13 = dUz_dx
  e(:,:,2) = 0.5d0*( dUz_dxi*m%dxi_dy + dUz_deta*m%deta_dy ) ! e23 = dUz_dy

  end function MAT_strain_SH

!-----------------------------------------------------------------------
  function MAT_strain_PSV(d,m,ngll) result(e)

  use mxmlib

  integer, intent(in) :: ngll
  double precision, intent(in) :: d(ngll,ngll,2)
  type(derint_type), intent(in) :: m
  double precision :: e(ngll,ngll,3)

  double precision, dimension(ngll,ngll) :: dUx_dxi,dUy_dxi,dUx_deta,dUy_deta

 ! local gradient
  dUx_dxi  = mxm( m%Ht, d(:,:,1), ngll )
  dUy_dxi  = mxm( m%Ht, d(:,:,2), ngll )
  dUx_deta = mxm( d(:,:,1), m%H, ngll )
  dUy_deta = mxm( d(:,:,2), m%H, ngll )

 ! strain 
  e(:,:,1) = dUx_dxi * m%dxi_dx + dUx_deta * m%deta_dx
  e(:,:,2) = dUy_dxi * m%dxi_dy + dUy_deta * m%deta_dy
  e(:,:,3) = 0.5d0*( dUx_dxi * m%dxi_dy + dUx_deta * m%deta_dy  &
                   + dUy_dxi * m%dxi_dx + dUy_deta * m%deta_dx  )

  end function MAT_strain_PSV

  
!=======================================================================
  function MAT_divcurl(d,m,g,enum,ngll,ndof) result(e)

  use spec_grid, only : sem_grid_type
  use fields_class, only : FIELD_divcurl_elem
  use stdio, only : IO_abort

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: d(ngll,ngll,ndof)
  type(matwrk_elem_type), intent(in) :: m
  type(sem_grid_type), intent(in) :: g
  integer, intent(in) :: enum
  double precision :: e(ngll,ngll,2)
  
  if (ndof/=2) call IO_abort('mat_gen:MAT_divcurl_gen: only for P-SV')

  if (associated(m%derint)) then
    e = MAT_divcurl_PSV(d,m%derint,ngll)
  else
    e = FIELD_divcurl_elem(d,ngll,g,enum)
  endif
  
  end function MAT_divcurl

!-----------------------------------------------------------------------
  function MAT_divcurl_PSV(d,m,ngll) result(e)

  use mxmlib

  integer, intent(in) :: ngll
  double precision, intent(in) :: d(ngll,ngll,2)
  type(derint_type), intent(in) :: m
  double precision :: e(ngll,ngll,2)

  double precision, dimension(ngll,ngll) :: dUx_dxi,dUy_dxi,dUx_deta,dUy_deta

 ! local gradient
  dUx_dxi  = mxm( m%Ht, d(:,:,1), ngll )
  dUy_dxi  = mxm( m%Ht, d(:,:,2), ngll )
  dUx_deta = mxm( d(:,:,1), m%H, ngll )
  dUy_deta = mxm( d(:,:,2), m%H, ngll )

 ! div = dUx/dx + dUy/dy
  e(:,:,1) = dUx_dxi * m%dxi_dx + dUx_deta * m%deta_dx &
           + dUy_dxi * m%dxi_dy + dUy_deta * m%deta_dy

 ! curl = dUx/dy - dUy/dx
  e(:,:,2) = dUx_dxi * m%dxi_dy + dUx_deta * m%deta_dy &
           - dUy_dxi * m%dxi_dx - dUy_deta * m%deta_dx  

  end function MAT_divcurl_PSV

!=======================================================================
! computes element forces from element stresses

  function MAT_forces(s,m,ngll,ndof) result(f)

  use mxmlib
  use stdio, only : IO_abort

  integer, intent(in) :: ngll,ndof
  double precision, intent(in) :: s(ngll,ngll,ndof+1)
  type(derint_type), intent(in) :: m
  double precision :: f(ngll,ngll,ndof)

  double precision, dimension(ngll,ngll) :: tmp1,tmp2 

  if (ndof==1) then
    tmp1 = -m%weights * ( m%dxi_dx  * s(:,:,1) + m%dxi_dy  * s(:,:,2) )
    tmp2 = -m%weights * ( m%deta_dx * s(:,:,1) + m%deta_dy * s(:,:,2) )
    f(:,:,1) = mxm(m%H, tmp1, ngll) + mxm(tmp2, m%Ht, ngll)

  elseif(ndof==2) then
    tmp1 = -m%weights * ( m%dxi_dx  * s(:,:,1) + m%dxi_dy  * s(:,:,3) )
    tmp2 = -m%weights * ( m%deta_dx * s(:,:,1) + m%deta_dy * s(:,:,3) )
    f(:,:,1) = mxm(m%H, tmp1, ngll) + mxm(tmp2, m%Ht, ngll)

    tmp1 = -m%weights * ( m%dxi_dx  * s(:,:,3) + m%dxi_dy  * s(:,:,2) )
    tmp2 = -m%weights * ( m%deta_dx * s(:,:,3) + m%deta_dy * s(:,:,2) )
    f(:,:,2) = mxm(m%H, tmp1, ngll) + mxm(tmp2, m%Ht, ngll)

  else
    call IO_abort('mat_common:MAT_forces: illegal value for argument ndof')
  endif

  end function MAT_forces

!=======================================================================
! finds the diagonal of the stiffness matrix, and stores its inverse

  subroutine MAT_diag_stiffness_init(invKDiag,g,f,matwrk)

  use spec_grid, only : sem_grid_type
  use fields_class, only : fields_type

  double precision, dimension(:,:), pointer :: invKDiag
  type(sem_grid_type), intent(in) :: g
  type(fields_type), intent(in) :: f
  type(matwrk_elem_type), intent(in) :: matwrk(:)

  double precision, dimension(g%ngll,g%ngll,f%ndof) :: dloc,floc
  integer :: e, i, j, k, dof
    
  allocate(invKDiag(f%npoin,f%ndof))
  call storearray('invKDiag',size(invKDiag),idouble)
  invKDiag = 0.d0
  dloc = 0.0d0

  do e=1,g%nelem
  do j=1,g%ngll
  do i=1,g%ngll
    k = g%ibool(i,j,e)
    do dof=1,f%ndof
      dloc(i,j,dof) = 1.0d0
      call MAT_ELAST_f(floc,dloc,matwrk(e)%elast,g%hprime,g%hTprime,g%ngll,f%ndof)
      invKDiag(k,dof) = invKDiag(k,dof) + floc(i,j,dof)
      dloc(i,j,dof) = 0.0d0
    enddo
  enddo
  enddo
  enddo

  invKDiag = 1d0/invKDiag

  end subroutine MAT_diag_stiffness_init

end module mat_gen
