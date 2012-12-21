! SEM2DPACK version 2.2.12c -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module bnd_grid
    
  use stdio, only : IO_abort
  use constants

  implicit none
  private

!-----------------------------------------------------------------------
!-- boundary mesh :
!
!     tag    = unique id
!     nelem  = total number of bundary elements
!     npoin  = total number of boundary nodes (non redundant)
!     ngnod  = number of nodes per element
!     elem   = #bnd_element --> #bulk_element
!     edge   = #bnd_element --> #edge in bulk element
!     node   = #bnd_node    --> #bulk_node
!     ibool  = (#bnd_node_local_index,#bnd_element) --> #bnd_node


  type bnd_grid_type
    integer :: nelem=0,npoin=0,ngnod=0,tag=0
    integer, pointer :: ibool(:,:)=>null()
    integer, dimension(:), pointer :: elem=>null() &
                                     ,edge=>null() &
                                     ,node=>null()
  end type bnd_grid_type

  public :: bnd_grid_type

!  contains

end module bnd_grid
