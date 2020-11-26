module memory_info

!=======================================================================
!
!     "memory_info" : for directory of dynamically allocated arrays
!      ----------
!
!=======================================================================

  implicit none
  private

  integer, parameter, public :: iinteg = 1, isngl = 2, idouble = 3

! This is the array directory, it is GLOBAL
! WARNING: a better implementation would use linked lists
  integer, parameter :: maxnbarrays = 250
  integer, save :: nbarrays =0
  integer*8, dimension(maxnbarrays),save :: arraysizes=0
  integer, dimension(maxnbarrays),save :: arraytypes=0
  integer, parameter :: NAME_LEN = 16
  character(len=NAME_LEN), dimension(maxnbarrays), save :: arraynames  = ' '

! ici codage en dur des tailles des variables en octets
! mettre iratio = 1 pour le Cray, iratio = 2 pour les autres machines
  integer, parameter :: iratio = 2

  public :: MEMO_echo,storearray

contains

subroutine MEMO_echo
!=======================================================================
!
!     Dynamic storage allocation :
!     --------------------------
!
!     Print a directory listing of all dynamically allocated arrays
!       and their properties
!
!=======================================================================

  use stdio, only : IO_new_unit

  integer :: iout
  ! devel: should be long integers, these overflow for very large simulations
  integer*8 :: itotsize,iarray,isize
  character(len=7) :: label(3)
  integer :: isizevars(3)

  isizevars(1) = 8/iratio  ! integer, in bytes = bit_size(0)/8
  isizevars(2) = 8/iratio  ! single precision
  isizevars(3) = 8         ! double precision

  label(1) = 'Integer'
  label(2) = 'Real   '
  label(3) = 'Double '


  iout = IO_new_unit()
  open(iout,file='MemoryInfo_sem2d.txt')
  write(iout,90) ''

  itotsize = 0
  do iarray = 1,nbarrays
   ! size in bytes
    isize = arraysizes(iarray)*isizevars(arraytypes(iarray))
   ! total size in bytes
    itotsize = itotsize + isize
    write(iout,110) iarray,arraynames(iarray),isize,itotsize
    !arraysizes(iarray),label(arraytypes(iarray))
  enddo

  write(iout,100) nbarrays,dble(itotsize)/dble(1024*1024),itotsize

  close(iout)

  90   format(/1x,51('=')/ &
  ' =   S E M 2 D P A C K   m e m o r y   u s a g e   ='/1x,51('=')// &
  '       Name            Size (bytes)     Cumulative'/1x,51('-'),a/)
  100   format(//,' Total number of allocated arrays. . . . .',i11/ &
  ' Total size of arrays in megabytes . . . .',f11.3/ &
  '                      in bytes . . . . . .',i11///) 
  110   format(i4,3x,a16,2x,i10,5x,i10) !WARNING: "a" must be NAME_LEN long

end subroutine MEMO_echo

!=====================================================================

! or other types try:
! size_of_foo = size( transfer(foo, (/ 0d0 /) ))

subroutine storearray(name,isize,itype)
!
!=======================================================================
!
!     Dynamic storage : store the array properties
!     ----------------
!
!=======================================================================

  use stdio, only : IO_abort

  character*(*) name
  integer, intent(in) :: isize,itype

  if (itype /= iinteg .and. itype /= isngl .and. itype /= idouble) &
    call IO_abort('Wrong array type in dynamic allocation')

  if (isize < 0) call IO_abort('Incoherent array size in dynamic allocation')
  if (isize ==0) return

  nbarrays = nbarrays + 1
  if(nbarrays > maxnbarrays) call IO_abort('Maximum number of arrays reached')

  arraysizes(nbarrays) = isize
  arraytypes(nbarrays) = itype
  arraynames(nbarrays) = name

  call MEMO_echo()

end subroutine storearray

end module memory_info
