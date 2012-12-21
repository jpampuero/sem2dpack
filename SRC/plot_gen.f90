! SEM2DPACK version 2.2.12beta -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
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
module plot_gen

  use plot_visual3
  use plot_avs
  use plot_postscript
  use plot_gmt

  implicit none
  private

  integer, parameter :: nb_fields = 5
  logical, save :: selected_fields(nb_fields)
  character(nb_fields), parameter :: field_names='DVAES'

  integer, parameter :: nb_comps = 4
  logical, save :: selected_comps(nb_comps)
  character(nb_comps), parameter :: comp_names='xyza'

  logical, save :: visual3,avs,postscript,bin,gmt

  public :: read_plot_gen,PLOT_FIELD, PLOT_init

contains

!=======================================================================
!
! BEGIN INPUT BLOCK
!
! NAME   : PLOTS
! PURPOSE: Selects a format to export snapshots
! SYNTAX : &PLOTS fields,components, bin,visual3,avs,postscript,gmt /
!
! ARG: fields   [char*] ['V'] fields to export in snapshots
!                 (begining of output file names given in parenthesis)
!                 'D'     displacements (dx,dy,dz,da)
!                 'V'     velocity (vx,vy,vz,va)
!                 'A'     acceleration (ax,ay,az,aa)
!                 'E'     strain (e11,e22,e12,e23,e13)
!                 'S'     stress (s11,s22,s12,s33,e13,e23)
! ARG: components [char*] ['ya'] components for PostScript outputs
!			'x','y','z' and/or 'a' (amplitude) 
!			(in SH only 'y' is considered)
! ARG: postscript [log] [T] PostScript
! ARG: gmt      [log] [F] output triangulation file grid_sem2d.gmt
!                 to be used in "pscontour -T" of the General Mapping Tool (GMT)
! ARG: avs      [log] [F] AVS
! ARG: visual3  [log] [F] Visual3
! ARG: bin      [log] [T] binary
!               
! NOTE   : If you choose PostScript you may need also a $POSTSCRIPT input block.
!          Other formats apply only to 'DVA' fields, 'ES' are exported as binary.
!
! END INPUT BLOCK

  subroutine read_plot_gen(iin,ndof)
  
  use echo, only : iout,echo_input

  integer, intent(in) :: iin,ndof

  character(nb_fields) :: fields
  character(nb_comps) :: components
  integer :: i

  NAMELIST / PLOTS / visual3,avs,gmt,postscript,bin,fields,components

  visual3    = .false.
  avs        = .false.
  gmt        = .false.
  postscript = .true.
  bin        = .true.
  fields     = 'V'
  components = 'ya'
  
  rewind(iin)
  read(iin,PLOTS,END=100) 
100 continue

  do i=1,nb_fields
    selected_fields(i) = scan(fields,field_names(i:i)) > 0
  enddo
  
  do i=1,nb_comps
    selected_comps(i) = scan(components,comp_names(i:i)) > 0
  enddo
  
  if (ndof==1) then
    selected_comps(1)=.false.
    selected_comps(2)=.true.
    selected_comps(3)=.false.
    selected_comps(4)=.false.
  else
    selected_comps(2)=.false.
  endif

  if (echo_input) write(iout,200) postscript,gmt,avs,visual3,bin, &
                                  selected_fields(:),selected_comps(:)

  if (postscript) call POST_PS_read(iin)

  return

  200 format(//1x,'S n a p s h o t   O u t p u t s',/1x,31('='),//5x,&
  'Save results in PS file or not . . . . .(postscript) = ',L1/5x, &
  'Save grid triangulation for GMT. . . . . . . . (gmt) = ',L1/5x, &
  'Save results in AVS file or not. . . . . . . . (avs) = ',L1/5x, &
  'Save results in Visual3 file or not. . . . (visual3) = ',L1/5x, &
  'Save results in binary file or not . . . . . . (bin) = ',L1/5x, &
  'Selected fields :',/5x, &
  '  Displacement . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Velocity . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Acceleration . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Strain . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Stress . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  'Selected components for PostScript snapshots :',/5x, &
  '  X  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Y  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Z  . . . . . . . . . . . . . . . . . . . . . . . . = ',L1/5x, &
  '  Amplitude  . . . . . . . . . . . . . . . . . . . . = ',L1)

  end subroutine read_plot_gen


!=======================================================================

  subroutine PLOT_FIELD(pb,tag,it,stitle)

  use problem_class
  use stdio, only : IO_rw_field
  use mat_elastic, only : ELAST_strain_stress

  type(problem_type), intent(in) :: pb
  integer, intent(in) :: it
  character(*), intent(in) :: stitle,tag

  double precision, pointer :: f(:,:)
  double precision, allocatable :: ssf(:,:,:,:)
  character(30) :: fnames(10)
  character(50) :: ps_file_name
  character :: fchar
  integer :: i,nssf,k
  integer :: comp_indx(nb_comps)

  comp_indx = (/ 1,1,2,0 /)

  do i = 1,3

    if (.not.selected_fields(i)) cycle

   ! Select the displayed field
    select case(i)
      case(1)
        f => pb%fields%displ
        fchar = 'd'
      case(2)
        f => pb%fields%veloc
        fchar = 'v'
      case(3)
        f => pb%fields%accel
        fchar = 'a'
    end select
  
    if (pb%fields%ndof==1) then
      write(fnames(1),'(A,"y_",A)') fchar,trim(tag)
    else
      write(fnames(1),'(A,"x_",A)') fchar,trim(tag) 
      write(fnames(2),'(A,"z_",A)') fchar,trim(tag)
    endif

    if (postscript) then
      do k=1,nb_comps
      if (selected_comps(k)) then
        write(ps_file_name,'( A,A,"_",A,"_sem2d.ps" )') fchar,comp_names(k:k),trim(tag)
        call PLOT_PS(file=trim(ps_file_name), vfield=f,comp=comp_indx(k) &
                  ,it_in=it,time_in=pb%time%time &
                  ,grid=pb%grid,mat=pb%matpro,stitle=stitle &
                  ,src=pb%src,rec=pb%rec)
      endif
      enddo
    endif
    if (avs) call plotavs(f,pb%grid,it)
    if (visual3) call plotvisual3(f,pb%grid,pb%matpro,it)
    if (bin) call IO_rw_field(f,fnames(1:2),'w')

  enddo

  nssf = 0

  if (pb%fields%ndof==1) then

    if (selected_fields(4)) then
      nssf=2
      write(fnames(1),'("e13_",A)') trim(tag) 
      write(fnames(2),'("e23_",A)') trim(tag) 
    endif
    if (selected_fields(5)) then
      write(fnames(nssf+1),'("s13_",A)') trim(tag) 
      write(fnames(nssf+2),'("s23_",A)') trim(tag) 
      nssf=nssf+2
    endif

  else

    if (selected_fields(4)) then
      nssf=3
      write(fnames(1),'("e11_",A)') trim(tag) 
      write(fnames(2),'("e22_",A)') trim(tag) 
      write(fnames(3),'("e12_",A)') trim(tag) 
    endif
    if (selected_fields(5)) then
      write(fnames(nssf+1),'("s11_",A)') trim(tag) 
      write(fnames(nssf+2),'("s22_",A)') trim(tag) 
      write(fnames(nssf+3),'("s33_",A)') trim(tag) 
      write(fnames(nssf+4),'("s12_",A)') trim(tag) 
      nssf=nssf+4
    endif

  endif

  if (nssf>0) then
    allocate(ssf(pb%grid%ngll,pb%grid%ngll,pb%grid%nelem,nssf))
    call ELAST_strain_stress(pb%matpro,pb%grid,pb%fields%displ, &
                             selected_fields(4),selected_fields(5),ssf)
    call IO_rw_field(ssf,fnames,'w')
    deallocate(ssf)
  endif
  
  end subroutine PLOT_FIELD
  
!=======================================================================
  subroutine PLOT_init(grid)

  use spec_grid, only : sem_grid_type

  type (sem_grid_type), intent(inout) :: grid

  if (postscript) call POST_PS_init(grid)
  if (gmt) call PLOT_GMT_init(grid%ibool)
  
  end subroutine PLOT_init

end module plot_gen
