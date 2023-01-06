module plot_avs

  implicit none

contains

! routine sauvegarde fichier AVS
  subroutine plotavs(field,grid,it)

  use spec_grid, only : sem_grid_type
  use stdio, only : IO_new_unit
  use echo, only : echo_run,iout,fmt1,fmtok

  type (sem_grid_type), intent(in) :: grid
  
  integer, intent(in) :: it
  double precision, intent(in) :: field(:,:)

  integer :: icell,i,j,e,fid,ip
  double precision :: rmaxnorm
  character(50) :: name

  fid = IO_new_unit()
  write(name,"('AVS_',I5.5,'_sem2d.inp')") it
  if (echo_run) write(iout,fmt1,advance='no') "Dump AVS file "//trim(name)
  open(unit=fid,file=trim(name),status='unknown')

  ! number of nodes, cells, data per cell
  write(fid,180) grid%npoin,grid%nelem*(grid%ngll-1)*(grid%ngll-1),1,0,0

  ! index and coordinates of nodes (3D coord = 0)
  do ip=1,grid%npoin
    write(fid,200) ip,grid%coord(:,ip)
  enddo

  ! index and topology of the cells
  icell = 0
  do e=1,grid%nelem
  do i=1,grid%ngll-1
  do j=1,grid%ngll-1

    icell = icell + 1
    write(fid,210) icell,grid%tag(e),grid%ibool(i,j+1,e), &
      grid%ibool(i,j,e),grid%ibool(i+1,j,e),grid%ibool(i+1,j+1,e)

  enddo
  enddo
  enddo

  ! dummy structure data vector and labels 
  write(fid,*) ' 1 1'
  write(fid,*) ' Label1, Label2'

  ! normalized node data (scalar)
  rmaxnorm = maxval(abs(field))
  if ( size(field,2)==1 ) then 
    do ip=1,grid%npoin
      write(fid,205) ip,sqrt(field(ip,1)**2 +field(ip,2)**2)/rmaxnorm
    enddo
  else
    do ip=1,grid%npoin
      write(fid,205) ip,field(ip,1)/rmaxnorm
    enddo
  endif

  close(fid)

  if (echo_run) write(iout,fmtok)

180 format(5(I0,1x))
200 format(I0,1x,e12.5,' 0. ',e12.5)
205 format(I0,1x,e12.5)
210 format(I0,1x,I0,' quad ',4(I0,1x))

  end subroutine plotavs

end module plot_avs
