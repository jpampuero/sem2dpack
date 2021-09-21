! This is a set of subroutines useful for the material "dynamic Damage" 
! (cf mat_dyn_damage.f90)

module utils_dyn_damage

  !modules used
  use constants, only : PI
  use stdio, only : IO_abort
  use prop_mat

  implicit none
  private

public :: stress_invariant, strain_invariant, sort_desc, kronecker &
        , vec2mat, mat2vec
!public :: KI_parameters, compute_KI2, stress_invariant, strain_invariant &
!		, stress_damage, kronecker, rupture_velocity, ini_dyndmg_cst, &
!		, eigs, sort_desc, damage_var, derivs

contains

!=======================================================================
! Compute the invariants of the stress tensor matrix.
!=======================================================================
  subroutine vec2mat(mat,vec,ndof)
!=======================================================================
  !Input variable
  !----------------
  !vec	: vector
  !ndof : integer defining antiplane/ plane strain case
  
  !output variable
  !----------------
  !mat  : matrice
!=======================================================================
  implicit none
  
  !input/ouput variables
  integer, intent(in) :: ndof
  double precision, intent(in) :: vec(ndof+1)
  double precision, intent(out) :: mat(3,3)
 
   mat(:,:) = 0.0d0
  	if (ndof==2) then !PSV
		mat(1,1)=vec(1)
		mat(2,2)=vec(2)
		mat(1,2)=vec(3)
		mat(2,1)=vec(3)
	else
		mat(1,3)=vec(1)
		mat(3,1)=vec(1)
		mat(2,3)=vec(2)
		mat(3,2)=vec(2)
    endif   	

  end subroutine vec2mat


!=======================================================================
! Compute the invariants of the stress tensor matrix.
!=======================================================================
  subroutine mat2vec(mat,vec,ndof)
!=======================================================================
  !Input variable
  !----------------
  !mat  : matrice
  !ndof : integer defining antiplane/ plane strain case
  
  !output variable
  !----------------
  !vec	: vector
!=======================================================================
  implicit none
  
  !input/ouput variables
  integer, intent(in) :: ndof
  double precision, intent(out) :: vec(ndof+1)
  double precision, intent(in) :: mat(3,3)

  	if (ndof==2) then !PSV
		vec(1)=mat(1,1)
		vec(2)=mat(2,2)
		vec(3)=mat(1,2)
	else
		vec(1)=mat(1,3)
		vec(2)=mat(2,3)
    endif   	
  end subroutine mat2vec
  

  
!===========================================================

!=======================================================================
! Compute the invariants of the stress tensor matrix.
!=======================================================================
  subroutine stress_invariant(stress_tensor,sigma,tau)
!=======================================================================
  !Input variable
  !----------------
  !stress	: stress tensor
  
  !output variable
  !----------------
  !sigma	: first invariant of the stress tnesor
  !tau		: second invariant of the stress tensor
!=======================================================================
  !modules used
  use constants, only : PI
  
  implicit none
  
  !input/ouput variables
  integer, parameter :: n=3
  double precision, intent(in) :: stress_tensor(n,n)
  double precision, intent(out) :: sigma,tau
  
  !local variables
  double precision :: abserr2=1.0e-09
  double precision :: delta(n,n), dev_stress(n,n), dev(n,n), eige(n), eigv(n,n)
  integer i, j
  
  !compute the eigenvalues
  !call eigs(stress_tensor,eige,eigv,abserr2,n)
  
  !compute 1st invariant
  sigma=0.0d0
  do i=1,n
    sigma=sigma+stress_tensor(i,i)
  end do
  sigma=sigma/dble(n)
  
 !compute 2nd invariant
  call kronecker(delta,n)
  dev_stress=stress_tensor-(sigma*delta)
  do i=1,n
    do j=1,n
    	dev(i,j)=dev_stress(i,j)*dev_stress(j,i)
    end do
  end do
  tau = sqrt(0.5*sum(dev)) 
 
  end subroutine stress_invariant
  
!=======================================================================
! Compute the invariants of the strain tensor matrix.
!=======================================================================
  subroutine strain_invariant(strain,epsi,gamma)
!=======================================================================
  !Input variable
  !----------------
  !strain	: strain tensor
   
  !output variable
  !----------------
  !epsi 	: first invariant of the strain tensor
  !gamma	: second invariant of the strain tensor
!=======================================================================
  !modules used
  use constants, only : PI
  
  implicit none
  
  !input/ouput variables
  integer, parameter :: n=3
  double precision, intent(in) :: strain(n,n)
  double precision, intent(out) :: epsi, gamma
  
  !local variables
  double precision :: abserr2=1.0e-09
  double precision :: delta(n,n), dev_strain(n,n), dev(n,n), eige(n), eigv(n,n)
  integer i, j
  
  !compute the eigenvalues
  !call eigs(strain,eige,eigv,abserr2,n)
  
  !compute 1st invariant
  epsi = 0.0d0
  do i=1,n
    epsi=epsi+strain(i,i)
  end do

 !compute 2nd invariant
  call kronecker(delta,n)
  dev_strain=strain-((epsi/3)*delta)
  do i=1,n
    do j=1,n
    	dev(i,j)=dev_strain(i,j)*dev_strain(j,i)
    end do
  end do
  gamma = sqrt(2*sum(dev)) 
 
  end subroutine strain_invariant


!=======================================================================
! Evaluate eigenvalues and eigenvectors of a real symmetric matrix
! modified from 
! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch07/Jacobi.f90
!=======================================================================
!  subroutine eigs(a,eige,eigvs,abserr2,n)
!!subroutine Jacobi(a,x,abserr,n)
!!=======================================================================
!! Evaluate eigenvalues and eigenvectors
!! of a real symmetric matrix a(n,n): a*x = lambda*x 
!! method: Jacoby method for symmetric matrices 
!! Alex G. (December 2009)
!!-----------------------------------------------------------
!! input ...
!! a(n,n) - array of coefficients for matrix A
!! n      - number of equations
!! abserr2 - abs tolerance [sum of (off-diagonal elements)^2]
!! output ...
!! eige(i,i) - eigenvalues
!! eigv(i,j) - eigenvectors
!! comments ...
!!===========================================================
!  implicit none
!
!  !input/ouput variables
!  integer, intent(in) :: n
!  double precision, intent(in) :: a(n,n), abserr2
!  double precision, intent(out) :: eige(n), eigvs(n,n)
!
!  !local variables
!  integer i, j, k, ind(n)
!  double precision :: b2, bar
!  double precision :: beta, coeff, c, s, cs, sc
!  double precision :: eigemat(n,n), eigv(n,n)
!  
!
!! special case if e12=e13=e23=0
!if (a(1,2)==0 .and. a(1,3)==0 .and. a(2,3)==0) then
!	
!	!eigenvalues are the trace of the matrix
!	do i = 1,n
!  		eige(i) = a(i,i)
!	end do
!	!eigenvactors are the Koenecker matrix 
!	call kronecker(eigv,n)
!
!!-------------------JACOBI code ----------------------------------------
!else
!! initialize eigv(i,j)=0, eigv(i,i)=1
!! *** the array operation eigv=0.0 is specific for Fortran 90/95
!eigv = 0.0
!do i=1,n
!  eigv(i,i) = 1.0
!end do
!
!! find the sum of all off-diagonal elements (squared)
!b2 = 0.0
!do i=1,n
!  do j=1,n
!    if (i.ne.j) b2 = b2 + a(i,j)**2
!  end do
!end do
!
!if (b2 <= abserr2) return
!
!! average for off-diagonal elements /2
!bar = 0.5*b2/float(n*n)
!
!!Find eigenvalues and eigenvectors
!eigemat = a
!do while (b2.gt.abserr2)
!  do i=1,n-1
!    do j=i+1,n
!      if (eigemat(j,i)**2 <= bar) cycle  ! do not touch small elements
!      b2 = b2 - 2.0*eigemat(j,i)**2
!      bar = 0.5*b2/float(n*n)
!! calculate coefficient c and s for Givens matrix
!      beta = (eigemat(j,j)-eigemat(i,i))/(2.0*eigemat(j,i))
!      coeff = 0.5*beta/sqrt(1.0+beta**2)
!      s = sqrt(max(0.5+coeff,0.0))
!      c = sqrt(max(0.5-coeff,0.0))
!! recalculate rows i and j
!      do k=1,n
!        cs =  c*eigemat(i,k)+s*eigemat(j,k)
!        sc = -s*eigemat(i,k)+c*eigemat(j,k)
!        eigemat(i,k) = cs
!        eigemat(j,k) = sc
!      end do
!! new matrix a_{k+1} from a_{k}, and eigenvectors 
!      do k=1,n
!        cs =  c*eigemat(k,i)+s*eigemat(k,j)
!        sc = -s*eigemat(k,i)+c*eigemat(k,j)
!        eigemat(k,i) = cs
!        eigemat(k,j) = sc
!        cs =  c*eigv(k,i)+s*eigv(k,j)
!        sc = -s*eigv(k,i)+c*eigv(k,j)
!        eigv(k,i) = cs
!        eigv(k,j) = sc
!      end do
!    end do
!  end do
!end do
!do i = 1,n
!  eige(i) = eigemat(i,i)
!end do
!end if
!!-------------------JACOBI code ----------------------------------------
!
!!sort out (decrease) the eigenvalues
!call sort_desc(eige,ind,n)
!!accordingly sort the eigenvectors
!do i = 1,n
!  eigvs(1:,i)=eigv(1:,ind(i))
!end do
!
!
!end subroutine eigs
  
!=======================================================================
! sorts array x() into descending order.
!=======================================================================
  subroutine sort_desc(x,ind,size)
!=======================================================================
  !Input variable
  !----------------
  !x		: array to be sort
  !size		: size of array
  
  !output variable
  !----------------
  !x		: arrays sorted
  !ind 		: indices of array's element before sorting
!=======================================================================
  implicit none

  !input/ouput variables
  double precision, dimension(1:), intent(inout) :: x
  integer, intent(in) :: size
  integer, intent(out) :: ind(size) 

  !local variables
  integer :: i,j,location
  double precision :: xtemp(size), temp
  INTEGER, DIMENSION (1)  :: Tloc
  INTEGER, DIMENSION (1)  :: Tloc2
  
 !sort the matrix
  xtemp = x
  do i = 1, Size
 	temp = x(i)
  	tloc = maxloc(x(i:size))+i-1
  	x(i) = maxval(x(i:size))
  	x(tloc(1)) = temp
  	tloc2 = minloc(abs(xtemp-x(i)))
  	ind(i)=tloc2(1)
  end do
  
end subroutine sort_desc

!=======================================================================
! Kronecker delta function
!=======================================================================
  subroutine kronecker(x,n)
!=======================================================================
  !Input variable
  !----------------
  !n : size of array to be created
  
  !output variable
  !----------------
  !x : Kronecker delta function
!=======================================================================
  implicit none

  !input/ouput variables
  integer, intent(in) :: n
  double precision, intent(out) :: x(n,n)

  !local variables
  integer :: i,j
  
  !create the matrix
  do i = 1, n
	do j= 1, n
	  if (i==j) then
		x(i,j) = 1
	  else
		x(i,j) = 0
	  end if
	end do
  end do

end subroutine kronecker



end module utils_dyn_damage
