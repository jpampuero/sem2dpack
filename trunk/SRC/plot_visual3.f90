module plot_visual3

  implicit none

contains

! routine sauvegarde fichier Visual3
  subroutine plotvisual3(field,grid,mat,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit
  use prop_mat, only : matpro_elem_type, MAT_getProp
  use echo, only : echo_run,iout,fmt1,fmtok

  type (sem_grid_type), intent(in) :: grid
  type(matpro_elem_type), intent(in) :: mat(:)
  double precision    , intent(in) :: field(:,:)
  integer             , intent(in) :: it

  double precision rmaxnorm,zoffset,xdiff,zdiff
  double precision :: celem(grid%ngll,grid%ngll)
  integer :: icell,i,j,e,fid,ip
  character(50) :: name

  print *,'Entering Visual3 file generation...'

  fid = IO_new_unit()
  write(name,"('Visual3_',i5.5,'_sem2d.inp')") it
  if (echo_run) write(iout,fmt1,advance='no') "Dump Visual3 file "//trim(name)
  open(unit=fid,file=trim(name),status='unknown')

  ! dummy thickness = size of first element
  xdiff = grid%coord(1,grid%ibool(grid%ngll,1,1)) - grid%coord(1,grid%ibool(1,1,1))
  zdiff = grid%coord(2,grid%ibool(grid%ngll,1,1)) - grid%coord(2,grid%ibool(1,1,1))
  zoffset = sqrt(xdiff**2 + zdiff**2)

  ! number of nodes, cells, data per cell
  !write(fid,180) 2*grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),zoffset
  write(fid,180) size(field,2)*grid%npoin &
                 ,grid%nelem*(grid%ngll-1)*(grid%ngll-1),zoffset

  ! index and coordinates of nodes
  do ip=1,grid%npoin
    write(fid,200) ip,grid%coord(:,ip)
  enddo

  ! index and topology of the cells
  icell = 0
  do e=1,grid%nelem
    call MAT_getProp(celem,mat(e),'cp')

    do i=1,grid%ngll-1
    do j=1,grid%ngll-1

      icell = icell + 1

      write(fid,210) icell,celem(i,j), &
        grid%ibool(i,j,e),grid%ibool(i+1,j,e), &
        grid%ibool(i+1,j+1,e),grid%ibool(i,j+1,e)

    enddo
    enddo
  enddo

  ! normalized node data
  rmaxnorm = maxval(abs(field))
  do ip=1,grid%npoin
    write(fid,205) ip,field(ip,:)/rmaxnorm
  enddo

  close(fid)

  if (echo_run) write(iout,fmtok)

180 format(I0,1x,I0,1x,e12.5)
200 format(I0,1x,e12.5,1x,e12.5)
205 format(I0,1x,e12.5,1x,e12.5)
210 format(I0,1x,e12.5,1x,I0,1x,I0,1x,I0,1x,I0)

  end subroutine plotvisual3

end module plot_visual3
