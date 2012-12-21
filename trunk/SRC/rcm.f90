module rcmlib
! Reverse Cuthill-McKee reordering
! Subroutines by John Burkardt
! http://www.csit.fsu.edu/~burkardt/f_src/rcm/rcm.html

  implicit none
  private

  public :: genrcm,perm_inverse 

contains


subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
  node_num )

!*******************************************************************************
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Modified:
!
!   05 January 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), is nonzero for those nodes which are
!    to be considered.
!
!    Output, integer DEG(NODE_NUM), contains, for each  node in the connected
!    component, its degree.
!
!    Output, integer ICCSIZE, the number of nodes in the connected component.
!
!    Output, integer LS(NODE_NUM), stores in entries 1 through ICCSIZE the nodes
!    in the connected component, starting with ROOT, and proceeding 
!    by levels.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer deg(node_num)
  integer i
  integer iccsze
  integer ideg
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer ls(node_num)
  integer lvlend
  integer lvsize
  integer mask(node_num)
  integer nbr
  integer node
  integer root
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            iccsze = iccsze + 1
            ls(iccsze) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  return
end subroutine degree 


subroutine genrcm ( node_num, adj_num, adj_row, adj, perm )

!*******************************************************************************
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Modified:
!
!    04 January 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused 
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer iccsze
  integer mask(node_num)
  integer level_num
  integer level_row(node_num+1)
  integer num
  integer perm(node_num)
  integer root

  mask(1:node_num) = 1

  num = 1

  do i = 1, node_num
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!
      call root_find ( root, adj_num, adj_row, adj, mask, level_num, &
        level_row, perm(num), node_num )
!
!  RCM orders the component using ROOT as the starting node.
!
      call rcm ( root, adj_num, adj_row, adj, mask, perm(num), iccsze, &
        node_num )

      num = num + iccsze
!
!  We can stop once every node is in one of the connected components.
!
      if ( node_num < num ) then
        return
      end if

    end if

  end do

  return
end subroutine genrcm 


subroutine i_swap ( i, j )

!*******************************************************************************
!
!! I_SWAP swaps two integer values.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end subroutine i_swap 



subroutine ivec_reverse ( n, a )

!*******************************************************************************
!
!! IVEC_REVERSE reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n/2
    call i_swap ( a(i), a(n+1-i) )
  end do

  return
end subroutine ivec_reverse 


subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*******************************************************************************
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!    The root node chosen by the user is assigned level 1, and masked.
!    All (unmasked) nodes reachable from a node in level 1 are
!    assigned level 2 and masked.  The process continues until there
!    are no unmasked nodes adjacent to any node in the current level.
!    The number of levels may vary between 2 and NODE_NUM.
!
!  Modified:
!
!    28 October 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM).  On input, only nodes with nonzero
!    MASK are to be processed.  On output, those nodes which were included
!    in the level set have MASK set to 1.
!
!    Output, integer LEVEL_NUM, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the rooted 
!    level structure.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer level_num
  integer level_row(node_num+1)
  integer level(node_num)
  integer lvlend
  integer lvsize
  integer mask(node_num)
  integer nbr
  integer node
  integer root

  mask(root) = 0
  level(1) = root
  level_num = 0
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
    level_num = level_num + 1
    level_row(level_num) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
    do i = lbegin, lvlend

      node = level(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1)-1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          iccsze = iccsze + 1
          level(iccsze) = nbr
          mask(nbr) = 0
        end if

      end do

    end do
!
!  Compute the current level width (the number of nodes encountered.)
!  If it is positive, generate the next level.
!
    lvsize = iccsze - lvlend

    if ( lvsize <= 0 ) then
      exit
    end if

  end do

  level_row(level_num+1) = lvlend + 1
!
!  Reset MASK to 1 for the nodes in the level structure.
!
  mask(level(1:iccsze)) = 1

  return
end subroutine level_set 



subroutine perm_inverse ( n, perm, perm_inv )

!*******************************************************************************
!
!! PERM_INVERSE produces the inverse of a given permutation.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items permuted.
!
!    Input, integer PERM(N), a permutation.
!
!    Output, integer PERM_INV(N), the inverse permutation.
!
  implicit none

  integer n

  integer i
  integer perm(n)
  integer perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end subroutine perm_inverse 

subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

!*******************************************************************************
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!    It is used as the starting point for the RCM ordering.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM), a mask for the nodes.  Only 
!    those nodes with nonzero input mask values are considered by the 
!    routine.  The nodes numbered by RCM will have their mask values 
!    set to zero.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!  Local parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold 
!    the degree of the nodes in the section graph specified by mask and root.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer deg(node_num)
  integer fnbr
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer k
  integer l
  integer lbegin
  integer lnbr
  integer lperm
  integer lvlend
  integer mask(node_num)
  integer nbr
  integer node
  integer perm(node_num)
  integer root
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

  mask(root) = 0

  if ( iccsze <= 1 ) then
    return
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  If no neighbors, skip to next node in this level.
!
      if ( lnbr <= fnbr ) then
        cycle
      end if
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      k = fnbr

      do while ( k < lnbr )

        l = k
        k = k + 1
        nbr = perm(k)

        do while ( fnbr < l )

          lperm = perm(l)

          if ( deg(lperm) <= deg(nbr) ) then
            exit
          end if

          perm(l+1) = lperm
          l = l-1

        end do

        perm(l+1) = nbr

      end do

    end do

  end do
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call ivec_reverse ( iccsze, perm )

  return
end subroutine rcm 

subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*******************************************************************************
!
!! ROOT_FIND finds a pseudo-peripheral node.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    NODE_NUM-1 from the pseudo-peripheral node.
!
!  Modified:
!
!    28 October 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!    Gibbs, Poole, Stockmeyer,
!    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, pages 236-250, 1976.
!
!    Gibbs,
!    Algorithm 509: A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, pages 378-387, 1976.
!
!  Parameters:
!
!    Input/output, integer ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), specifies a section subgraph.  Nodes 
!    for which MASK is zero are ignored by FNROOT.
!
!    Output, integer LEVEL_NUM, is the number of levels in the level structure
!    rooted at the node ROOT.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), integer LEVEL(NODE_NUM), the 
!    level structure array pair containing the level structure found.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer iccsze
  integer j
  integer jstrt
  integer k
  integer kstop
  integer kstrt
  integer level(node_num)
  integer level_num
  integer level_num2
  integer level_row(node_num+1)
  integer mask(node_num)
  integer mindeg
  integer nabor
  integer ndeg
  integer node
  integer root
!
!  Determine the level structure rooted at ROOT.
!
  call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
    level_row, level, node_num )
!
!  Count the number of nodes in this level structure.
!
  iccsze = level_row(level_num+1) - 1
!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
  if ( level_num == 1 ) then
    return
  end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
  if ( level_num == iccsze ) then
    return
  end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
  do

    mindeg = iccsze

    jstrt = level_row(level_num)
    root = level(jstrt)

    if ( jstrt < iccsze ) then

      do j = jstrt, iccsze

        node = level(j)
        ndeg = 0
        kstrt = adj_row(node)
        kstop = adj_row(node+1)-1

        do k = kstrt, kstop
          nabor = adj(k)
          if ( 0 < mask(nabor) ) then
            ndeg = ndeg+1
          end if
        end do

        if ( ndeg < mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do

    end if
!
!  Generate the rooted level structure associated with this node.
!
    call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
      level_row, level, node_num )
!
!  If the number of levels did not increase, accept the new ROOT.
!
    if ( level_num2 <= level_num ) then
      exit
    end if

    level_num = level_num2
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
    if ( iccsze <= level_num ) then
      exit
    end if

  end do

  return
end subroutine root_find 


end module rcmlib
