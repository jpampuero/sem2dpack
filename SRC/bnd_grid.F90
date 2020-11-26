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
