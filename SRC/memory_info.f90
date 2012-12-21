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
! WARNING: a better implementation would  linked lists
  integer, parameter :: maxnbarrays = 250
  integer, save :: nbarrays =0
  integer, dimension(maxnbarrays),save :: arraysizes=0,arraytypes=0
  character(len=12), dimension(maxnbarrays),save :: &
    arraynames  = '                    '

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
  integer :: itotsize,iarray
  character(len=7) :: label(3)
  integer :: isizevars(3)

  isizevars(1) = 8/iratio  ! integer
  isizevars(2) = 8/iratio  ! single precision
  isizevars(3) = 8         ! double precision

  label(1) = 'Integer'
  label(2) = 'Real   '
  label(3) = 'Double '

! compute total size in bytes
  itotsize = 0
  do iarray = 1,nbarrays
    itotsize = itotsize + arraysizes(iarray)*isizevars(arraytypes(iarray))
  enddo

  iout = IO_new_unit()
  open(iout,file='MemoryInfo_sem2d.txt')

  write(iout,100) nbarrays,dble(itotsize)/dble(1024*1024),itotsize, &
                      itotsize/isizevars(3)

  do iarray = 1,nbarrays
    write(iout,110) iarray,arraysizes(iarray),arraynames(iarray), &
        label(arraytypes(iarray))
  enddo

  close(iout)

  100   format(//1x,41('=')/ &
  ' =  D i r e c t o r y     l i s t i n g  ='/1x,41('=')// &
  ' Total number of allocated arrays. . . . . . . . . .',i11/ &
  ' Total size of arrays in megabytes . . . . . . . . .',f11.3/ &
  ' Total size of arrays in bytes . . . . . . . . . . .',i11/ &
  ' Total size of arrays in double precision words. . .',i11/// &
  '  Array nb    Size         Name        Type'/1x,47('=')/)
  110   format(i6,3x,i10,5x,a12,2x,a7)

end subroutine MEMO_echo

!=====================================================================

subroutine storearray(name,isize,itype)
!
!=======================================================================
!
!     Dynamic storage : store the array properties
!     ----------------
!
!=======================================================================

  character*(*) name
  integer, intent(in) :: isize,itype

  if(itype /= iinteg .and. itype /= isngl .and. itype /= idouble) &
    stop 'Wrong array type in dynamic allocation'

  if(isize <= 0) stop 'Incoherent array size in dynamic allocation'

  nbarrays = nbarrays + 1
  if(nbarrays > maxnbarrays) stop 'Maximum number of arrays reached'

  arraysizes(nbarrays) = isize
  arraytypes(nbarrays) = itype
  arraynames(nbarrays) = name

end subroutine storearray

end module memory_info
