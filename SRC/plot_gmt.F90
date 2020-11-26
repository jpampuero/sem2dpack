module plot_gmt

  implicit none
  private

  public :: PLOT_GMT_init

contains
  
subroutine PLOT_GMT_init(ibool)

  use stdio, only : IO_new_unit
  use echo, only : echo_init,iout,fmt1,fmtok

  integer, intent(in) :: ibool(:,:,:)

  integer :: ngll,nelem,i,j,e,fid

  ngll = size(ibool,1)
  nelem = size(ibool,3)

  fid = IO_new_unit()
  open(unit=fid,file='grid_sem2d.gmt',status='unknown')
  if (echo_init) write(iout,fmt1,advance='no') "Write GMT triangulation file grid_sem2d.gmt"

  do e =1,nelem
  do j =1,ngll-1
  do i =1,ngll-1
    write(fid,*) ibool(i,j,e),ibool(i+1,j,e),ibool(i,j+1,e)
    write(fid,*) ibool(i+1,j,e),ibool(i+1,j+1,e),ibool(i,j+1,e)
  enddo
  enddo
  enddo
  
  close(fid)

  if (echo_init) write(iout,fmtok)

end subroutine PLOT_GMT_init

end module plot_gmt
