!=======================================================================
! Module: color_elem
! Description: Element-conflict graph construction and greedy coloring
!              for conflict-free parallel assembly in SEM2DPACK.
!=======================================================================
module color_elem

  use stdio, only : IO_abort
  use echo, only : echo_init, iout, fmt1, fmtok

  implicit none
  private

  type color_group_type
    integer :: nelem = 0
    integer, pointer :: elem(:) => null()      ! (nelem)
  end type color_group_type

  type elem_coloring_type
    integer :: ncolors = 0
    integer :: nelem = 0
    type(color_group_type), pointer :: colors(:) => null()
  end type elem_coloring_type

  public :: color_group_type, elem_coloring_type, &
            COLOR_build_and_validate, COLOR_free

contains

!=======================================================================
! Build element conflict graph from control nodes (knods), execute
! greedy coloring, group each color into a flat element list, validate
! conflict-freedom, and report color statistics.
!
! Using knods (4 or 9 control nodes per element) instead of ibool
! (ngll*ngll GLL nodes per element) produces the same conflict graph
! because two elements share a GLL node if and only if they share a
! control node, but with much less work.
!=======================================================================
  subroutine COLOR_build_and_validate(knods, nelem, ngnod, npoin, coloring)

    integer, intent(in) :: nelem, ngnod, npoin
    integer, intent(in) :: knods(ngnod, nelem)
    type(elem_coloring_type), intent(inout) :: coloring

    integer, allocatable :: node_count(:), node_ptr(:), node_cur(:), node_elem_list(:)
    integer, allocatable :: elem_color(:), node_owner(:)
    logical, allocatable :: color_used(:)
    integer :: total_entries, e, i, k, p, nbr, c, max_color, icol
    integer :: n_c, idx
    integer :: max_possible_colors

    if (echo_init) then
      write(iout, *)
      write(iout, '(a)') ' E l e m e n t   c o l o r i n g'
      write(iout, '(a)') ' ==============================='
      write(iout, fmt1, advance='no') 'Building element-conflict graph & coloring'
    endif

    ! Step 1: Count elements per control node (CSR structure)
    allocate(node_count(npoin))
    node_count = 0

    do e = 1, nelem
      do i = 1, ngnod
        p = knods(i, e)
        node_count(p) = node_count(p) + 1
      enddo
    enddo

    allocate(node_ptr(npoin + 1))
    node_ptr(1) = 1
    do p = 1, npoin
      node_ptr(p + 1) = node_ptr(p) + node_count(p)
    enddo

    total_entries = node_ptr(npoin + 1) - 1
    allocate(node_elem_list(total_entries))
    allocate(node_cur(npoin))
    node_cur = node_ptr(1:npoin)

    do e = 1, nelem
      do i = 1, ngnod
        p = knods(i, e)
        node_elem_list(node_cur(p)) = e
        node_cur(p) = node_cur(p) + 1
      enddo
    enddo

    deallocate(node_cur)
    deallocate(node_count)

    ! Step 2: Greedy coloring in element order e = 1 .. nelem
    max_possible_colors = 100
    allocate(color_used(max_possible_colors))
    color_used = .false.

    allocate(elem_color(nelem))
    elem_color = 0

    do e = 1, nelem
      do i = 1, ngnod
        p = knods(i, e)
        do k = node_ptr(p), node_ptr(p + 1) - 1
          nbr = node_elem_list(k)
          if (nbr < e) then
            c = elem_color(nbr)
            if (c > 0) then
              if (c > max_possible_colors) then
                call IO_abort('COLOR_build: max_possible_colors exceeded')
              endif
              color_used(c) = .true.
            endif
          endif
        enddo
      enddo

      c = 1
      do while (color_used(c))
        c = c + 1
        if (c > max_possible_colors) then
          call IO_abort('COLOR_build: max_possible_colors exceeded')
        endif
      enddo
      elem_color(e) = c

      do i = 1, ngnod
        p = knods(i, e)
        do k = node_ptr(p), node_ptr(p + 1) - 1
          nbr = node_elem_list(k)
          if (nbr < e) then
            c = elem_color(nbr)
            if (c > 0) color_used(c) = .false.
          endif
        enddo
      enddo
    enddo

    deallocate(color_used)
    deallocate(node_elem_list)
    deallocate(node_ptr)

    max_color = maxval(elem_color)
    coloring%ncolors = max_color
    coloring%nelem = nelem

    allocate(coloring%colors(max_color))

    ! Step 3: Group each color into a flat element list
    do icol = 1, max_color
      n_c = count(elem_color == icol)
      coloring%colors(icol)%nelem = n_c

      allocate(coloring%colors(icol)%elem(n_c))

      idx = 0
      do e = 1, nelem
        if (elem_color(e) == icol) then
          idx = idx + 1
          coloring%colors(icol)%elem(idx) = e
        endif
      enddo
    enddo

    if (echo_init) write(iout, fmtok)

    ! Step 4: Validate that within every color no two elements share a control node
    if (echo_init) write(iout, fmt1, advance='no') 'Validating conflict-free coloring'

    allocate(node_owner(npoin))
    node_owner = 0

    do icol = 1, max_color
      do idx = 1, coloring%colors(icol)%nelem
        e = coloring%colors(icol)%elem(idx)
        do i = 1, ngnod
          p = knods(i, e)
          if (node_owner(p) /= 0 .and. node_owner(p) /= e) then
            write(iout, '(A,I0,A,I0,A,I0,A,I0)') &
              'ERROR: Conflict in color ', icol, ' at node ', p, &
              ' between elements ', node_owner(p), ' and ', e
            call IO_abort('COLORING VALIDATION FAILED: Elements share node within color')
          endif
          node_owner(p) = e
        enddo
      enddo

      node_owner = 0
    enddo

    deallocate(node_owner)
    deallocate(elem_color)

    if (echo_init) then
      write(iout, fmtok)
      write(iout, 100) 'Total number of colors . . . . . . . . . = ', coloring%ncolors
      do icol = 1, coloring%ncolors
        write(iout, 110) 'Color ', icol, ': ', coloring%colors(icol)%nelem, ' elements'
      enddo
      write(iout, 120) 'Coloring validation . . . . . . . . . . = PASSED (no conflicts)'
      write(iout, *)
    endif

100 format(5X,A,I0)
110 format(7X,A,I0,A,I0,A)
120 format(5X,A)

  end subroutine COLOR_build_and_validate

!=======================================================================
! Deallocate dynamic structures in elem_coloring_type
!=======================================================================
  subroutine COLOR_free(coloring)
    type(elem_coloring_type), intent(inout) :: coloring
    integer :: icol

    if (associated(coloring%colors)) then
      do icol = 1, coloring%ncolors
        if (associated(coloring%colors(icol)%elem)) deallocate(coloring%colors(icol)%elem)
      enddo
      deallocate(coloring%colors)
      coloring%colors => null()
    endif
    coloring%ncolors = 0
    coloring%nelem = 0
  end subroutine COLOR_free

end module color_elem
