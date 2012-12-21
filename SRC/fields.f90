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
module fields_class

  implicit none

  type fields_type
    integer :: ndof,npoin
    double precision, dimension(:,:), pointer ::  displ,veloc,accel,tmp
  end type fields_type

contains

!===================================================================================
!
  subroutine FIELDS_read(fields,file)

  use echo, only : echo_input,iout
  use stdio, only : IO_new_unit

  type(fields_type), intent(inout) :: fields
  character(50), intent(in) :: file

  integer :: n,inump,npoin,iunit
  
  if (echo_input) write(iout,'(/,"Reading initial fields from external file",/)')

  iunit = IO_new_unit()
  open(iunit,file=trim(file),status='old')
  do n = 1,fields%npoin
    read(iunit,*) inump, fields%displ(inump,:), fields%veloc(inump,:),fields%accel(inump,:)
  enddo
  close(iunit)
  
  end subroutine FIELDS_read

!===================================================================================
!
  double precision function FIELDS_get_max(field)

  double precision, intent(in) :: field(:,:)
  
!  FIELDS_get_max = sqrt(maxval(field(:,1)**2 + field(:,2)**2))
  FIELDS_get_max = max( abs(maxval(field)), abs(minval(field)) )

  end function FIELDS_get_max


!===================================================================================
!
  subroutine FIELDS_init(fields,npoin)

  type(fields_type), intent(out) :: fields
  integer, intent(in) :: npoin

  fields%npoin = npoin
  call FIELD_init(fields%displ,'displ')
  call FIELD_init(fields%veloc,'veloc')
  call FIELD_init(fields%accel,'accel')
  call FIELD_init(fields%tmp,'tmp')

  contains

!-- single field:
  subroutine FIELD_init(field,name)

  use memory_info

  double precision, pointer :: field(:,:)
  character(*), intent(in) :: name

  allocate(field(npoin,fields%ndof))
  call storearray(name,size(field),idouble)
  field = 0.d0

  end subroutine FIELD_init

  end subroutine FIELDS_init



end module fields_class
