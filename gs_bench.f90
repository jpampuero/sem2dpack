program gs_bench
! Microbenchmark for the SEM2DPACK gather (FIELD_get_elem_sub) and scatter
! (FIELD_add_elem_2) inner loops, using the real ibool from 2.5D_inplane.
! Compares the current global-field layout (npoin,ndof) against (ndof,npoin).
  implicit none
  integer, parameter :: nelem=12800, ngll=5, npoin=205761, ndof=2
  integer, parameter :: nrep=300
  integer :: ibool(ngll,ngll,nelem)
  integer :: ib_flat(ngll*ngll*nelem)
  integer :: eord(nelem)
  double precision :: rr
  integer :: ee, itmp
  double precision, allocatable :: Fnd(:,:)   ! (npoin,ndof)  current layout
  double precision, allocatable :: Fdn(:,:)   ! (ndof,npoin)  proposed layout
  double precision :: floc(ngll,ngll,ndof)
  integer(8) :: c0,c1,cr
  integer :: e,i,j,k,rep,idx
  double precision :: chk, t_g_nd, t_g_dn, t_s_nd, t_s_dn

  ! read ibool (stream of int32, 25 per element)
  open(10,file='ibool_sem2d.dat',access='stream',form='unformatted',status='old')
  read(10) ib_flat
  close(10)
  idx=0
  do e=1,nelem
    do j=1,ngll
    do i=1,ngll
      idx=idx+1
      ibool(i,j,e)=ib_flat(idx)
    enddo
    enddo
  enddo

  allocate(Fnd(npoin,ndof), Fdn(ndof,npoin))
  call random_number(Fnd); Fdn = 0d0
  do k=1,npoin
    Fdn(1,k)=Fnd(k,1); Fdn(2,k)=Fnd(k,2)
  enddo

  ! random element-visit order, to defeat the RCM locality and probe the
  ! cold / memory-bound regime
  do e=1,nelem
    eord(e)=e
  enddo
  do e=nelem,2,-1
    call random_number(rr)
    ee=1+int(rr*e)
    if (ee<1) ee=1
    if (ee>e) ee=e
    itmp=eord(e); eord(e)=eord(ee); eord(ee)=itmp
  enddo

  !============ GATHER, current layout (npoin,ndof) ============
  chk=0d0
  call system_clock(c0,cr)
  do rep=1,nrep
    do e=1,nelem
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        floc(i,j,:)=Fnd(k,:)
      enddo
      enddo
      chk=chk+sum(floc)
    enddo
  enddo
  call system_clock(c1)
  t_g_nd=dble(c1-c0)/dble(cr)
  print '(A,F8.3,A,ES12.4)','gather (npoin,ndof): ',t_g_nd,' s   chk=',chk

  !============ GATHER, proposed layout (ndof,npoin) ============
  ! (chk kept running to avoid dead-code elimination)
  call system_clock(c0,cr)
  do rep=1,nrep
    do e=1,nelem
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        floc(i,j,:)=Fdn(:,k)
      enddo
      enddo
      chk=chk+sum(floc)
    enddo
  enddo
  call system_clock(c1)
  t_g_dn=dble(c1-c0)/dble(cr)
  print '(A,F8.3,A,ES12.4)','gather (ndof,npoin): ',t_g_dn,' s   chk=',chk

  !============ SCATTER, current layout (npoin,ndof) ============
  floc=1d-6
  call system_clock(c0,cr)
  do rep=1,nrep
    do e=1,nelem
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        Fnd(k,:)=Fnd(k,:)+floc(i,j,:)
      enddo
      enddo
    enddo
  enddo
  call system_clock(c1)
  t_s_nd=dble(c1-c0)/dble(cr)
  print '(A,F8.3,A,ES12.4)','scatter (npoin,ndof):',t_s_nd,' s   chk=',sum(Fnd(1:10,1))

  !============ SCATTER, proposed layout (ndof,npoin) ============
  floc=1d-6
  call system_clock(c0,cr)
  do rep=1,nrep
    do e=1,nelem
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        Fdn(:,k)=Fdn(:,k)+floc(i,j,:)
      enddo
      enddo
    enddo
  enddo
  call system_clock(c1)
  t_s_dn=dble(c1-c0)/dble(cr)
  print '(A,F8.3,A,ES12.4)','scatter (ndof,npoin):',t_s_dn,' s   chk=',sum(Fdn(1,1:10))

  print '(/A)','--- RCM order (warm, representative of within-call) ---'
  print '(A,F6.2,A,F6.2,A,F5.2,A)','gather : ',t_g_nd/dble(nrep*nelem*ngll*ngll)*1d9, &
        ' -> ',t_g_dn/dble(nrep*nelem*ngll*ngll)*1d9,' ns   (',t_g_nd/t_g_dn,'x)'
  print '(A,F6.2,A,F6.2,A,F5.2,A)','scatter: ',t_s_nd/dble(nrep*nelem*ngll*ngll)*1d9, &
        ' -> ',t_s_dn/dble(nrep*nelem*ngll*ngll)*1d9,' ns   (',t_s_nd/t_s_dn,'x)'

  !============ GATHER random order, both layouts ============
  ! (chk kept running to avoid dead-code elimination)
  call system_clock(c0,cr)
  do rep=1,nrep
    do ee=1,nelem
      e=eord(ee)
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        floc(i,j,:)=Fnd(k,:)
      enddo
      enddo
      chk=chk+sum(floc)
    enddo
  enddo
  call system_clock(c1); t_g_nd=dble(c1-c0)/dble(cr)
  ! (chk kept running to avoid dead-code elimination)
  call system_clock(c0,cr)
  do rep=1,nrep
    do ee=1,nelem
      e=eord(ee)
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        floc(i,j,:)=Fdn(:,k)
      enddo
      enddo
      chk=chk+sum(floc)
    enddo
  enddo
  call system_clock(c1); t_g_dn=dble(c1-c0)/dble(cr)

  !============ SCATTER random order, both layouts ============
  floc=1d-6
  call system_clock(c0,cr)
  do rep=1,nrep
    do ee=1,nelem
      e=eord(ee)
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        Fnd(k,:)=Fnd(k,:)+floc(i,j,:)
      enddo
      enddo
    enddo
  enddo
  call system_clock(c1); t_s_nd=dble(c1-c0)/dble(cr)
  floc=1d-6
  call system_clock(c0,cr)
  do rep=1,nrep
    do ee=1,nelem
      e=eord(ee)
      do j=1,ngll
      do i=1,ngll
        k=ibool(i,j,e)
        Fdn(:,k)=Fdn(:,k)+floc(i,j,:)
      enddo
      enddo
    enddo
  enddo
  call system_clock(c1); t_s_dn=dble(c1-c0)/dble(cr)

  print '(/A,ES10.2)','--- random element order (poor locality) --- chk=',chk+sum(Fnd(1:10,1))+sum(Fdn(1,1:10))
  print '(A,F6.2,A,F6.2,A,F5.2,A)','gather : ',t_g_nd/dble(nrep*nelem*ngll*ngll)*1d9, &
        ' -> ',t_g_dn/dble(nrep*nelem*ngll*ngll)*1d9,' ns   (',t_g_nd/t_g_dn,'x)'
  print '(A,F6.2,A,F6.2,A,F5.2,A)','scatter: ',t_s_nd/dble(nrep*nelem*ngll*ngll)*1d9, &
        ' -> ',t_s_dn/dble(nrep*nelem*ngll*ngll)*1d9,' ns   (',t_s_nd/t_s_dn,'x)'
end program gs_bench
