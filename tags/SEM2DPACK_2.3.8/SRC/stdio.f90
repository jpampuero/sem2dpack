! SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
!                            with emphasis on computational seismology and earthquake source dynamics.
! 
! Copyright (C) 2003-2007 Jean-Paul Ampuero
! All Rights Reserved
! 
! Jean-Paul Ampuero
! 
! California Institute of Technology
! Seismological Laboratory
! 1200 E. California Blvd., MC 252-21 
! Pasadena, CA 91125-2100, USA
! 
! ampuero@gps.caltech.edu
! Phone: (626) 395-6958
! Fax  : (626) 564-0715
! 
! http://web.gps.caltech.edu/~ampuero/
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
module stdio

  implicit none
  private
  
  interface IO_rw_field
    module procedure IO_field_node_single, IO_field_node_multi, &
                     IO_field_elem_single, IO_field_elem_multi
  end interface IO_rw_field

!! Error output unit, default: standard output                              
  integer, save :: ierror = 6 

  public :: IO_read_skip, IO_new_unit, IO_test_msg, IO_rw_field, IO_abort, &
    IO_file_length, IO_file_columns

contains

!===================================================================
! Jumps in a input ASCII file
  subroutine IO_read_skip (nlines,iounit)

    integer, intent(in) :: nlines,iounit
    integer :: i

    do i = 1, nlines ; read(iounit,*) ; end do

  end subroutine IO_read_skip

!===================================================================
!! Creates a new io unit, avoiding conflicts
  integer function IO_new_unit()

    integer :: i 
    logical :: not_available 

    i = 100
    not_available = .true.
    do while (not_available)
      i = i+1
      inquire(i,OPENED=not_available)
    enddo
    IO_new_unit = i

  end function IO_new_unit


!=====================================================================
!! File length
  function IO_file_length(filename) result(n)
  
  character(*), intent(in) :: filename
  integer :: n, io
  
  io = IO_new_unit()
  open(io,file=filename,status='old')
  n=0
  do
    read(io,*,END=100)
    n=n+1
  enddo
100 continue
  close(io)
  
  end function IO_file_length

!=====================================================================
!! Count the number of columns in a text file
  function IO_file_columns(filename) result(ncol)

  character(*) :: filename
  character(500) :: line !WARNING: assumed max line length
  integer  :: ncol,i,iin
  double precision :: xread !WARNING: assume double precision inputs

  iin = IO_new_unit()
  open(iin,file=filename)
  read(iin,'(A)') line
  close(iin)

  do ncol=1,huge(ncol)
    read(line,*,end=200) (xread, i=1,ncol)
  enddo
  200 continue
  ncol = ncol-1

  end function IO_file_columns

!=====================================================================
!! Tests a section name in input file
  logical function IO_test_msg(iin,msg)
    
  integer, intent(in) :: iin
  character(*), intent(in) :: msg
  character(LEN(msg)) :: msg_read

  read(iin,'(A)',advance='no') msg_read
  IO_test_msg =   msg_read == msg

  end function IO_test_msg


!===================================================================
!!--- Reads or writes a field from binary file

! a time-tagged input name can be useful:
!    write(filename,222) fieldname,it
!  222 format(a,i5.5,'.dat') 

!-------------------------------------------------------------------
! field in node-wise storage, single component
  subroutine IO_field_node_single(field,filename,action)

    double precision, intent(inout) :: field(:)
    character(*), intent(in) :: filename
    character, intent(in) :: action

    integer :: iunit,iol,k

    iunit = IO_new_unit()
    if (action=='r') then
      INQUIRE( IOLENGTH=iol ) field(1)
      open(unit=iunit,file=trim(filename)//'_sem2d.dat',status='old',access='direct',recl=iol)
      do k=1,size(field)
        read(iunit,rec=k) field(k)
      enddo
    elseif (action=='w') then 
      INQUIRE( IOLENGTH=iol ) real(field(1))
      open(unit=iunit,file=trim(filename)//'_sem2d.dat',status='replace',access='direct',recl=iol)
      do k=1,size(field)
        write(iunit,rec=k) real(field(k))
      enddo
    else
      call IO_abort('IO_rw_field: illegal action')
    endif
    close(iunit)

  end subroutine IO_field_node_single

!-------------------------------------------------------------------
! field in node-wise storage, multiple components
  subroutine IO_field_node_multi(fields,filenames,action)

    double precision, intent(inout) :: fields(:,:)
    character(*), intent(in) :: filenames(:)
    character, intent(in) :: action

    integer :: i

    do i=1,size(fields,2)
      call IO_field_node_single(fields(:,i),filenames(i),action)
    enddo

  end subroutine IO_field_node_multi

!-------------------------------------------------------------------
! field in element-wise storage, single component
  subroutine IO_field_elem_single(field,filename,action)

    double precision, intent(inout) :: field(:,:,:)
    character(*), intent(in) :: filename
    character, intent(in) :: action

    integer :: iunit,iol,k

    iunit = IO_new_unit()
    if (action=='r') then
      INQUIRE( IOLENGTH=iol ) field(:,:,1)
      open(unit=iunit,file=trim(filename)//'_sem2d.dat',status='old',access='direct',recl=iol)
      do k=1,size(field,3)
        read(iunit,rec=k) field(:,:,k)
      enddo
    elseif (action=='w') then 
      INQUIRE( IOLENGTH=iol ) real(field(:,:,1))
      open(unit=iunit,file=trim(filename)//'_sem2d.dat',status='replace',access='direct',recl=iol)
      do k=1,size(field,3)
        write(iunit,rec=k) real(field(:,:,k))
      enddo
    else
      call IO_abort('IO_rw_field: illegal action')
    endif
    close(iunit)

  end subroutine IO_field_elem_single

!-------------------------------------------------------------------
! field in element-wise storage, multiple components
  subroutine IO_field_elem_multi(fields,filenames,action)

    double precision, intent(inout) :: fields(:,:,:,:)
    character(*), intent(in) :: filenames(:)
    character, intent(in) :: action

    integer :: i

    do i=1,size(fields,4)
      call IO_field_elem_single(fields(:,:,:,i),filenames(i),action)
    enddo

  end subroutine IO_field_elem_multi

!===================================================================
!! Gives an abort message and leaves program 
  subroutine IO_abort(message)

    character(*),intent(in) :: message

    write(ierror,*) message
    write(ierror,*) 'FATAL ERROR, aborting.'
    stop
    
  end subroutine IO_abort


end module stdio
 