
!-------------------------------------------------------------------
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 26.3.2018 TV: Cleaning
!
!NOTE: Whose code is this?
!-----------------------------
module possu
use common
implicit none 


contains

!____________________________________________________________________________
!
! Inverse of a real matrix (lapack)
!___________________________________________________________________
function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(:,:), allocatable :: Ainv

  real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
  integer, dimension(:), allocatable :: ipiv   ! pivot indices
  integer :: n, info

 
  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK

  n = size(A,1)
  allocate(Ainv(n,n))
  allocate(work(n))
  allocate(ipiv(n))
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

!____________________________________________________________________________
!
! Inverse of a complex matrix (lapack)
!___________________________________________________________________

function Cinv(A) result(Ainv)
  complex(dp), dimension(:,:), intent(in) :: A
  complex(dp), dimension(:,:), allocatable :: Ainv

  complex(dp), dimension(:), allocatable :: work  ! work array for LAPACK
  integer, dimension(:), allocatable :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  n = size(A,1)
  allocate(Ainv(n,n))
  allocate(work(n))
  allocate(ipiv(n))
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function Cinv


!___________________________________________________________________
!
!       SVD interface for (lapack)
!____________________________________________________________

subroutine svd(A,U,S,VT,M,N)
  complex(dp), dimension(:,:), intent(in) :: A

  complex(dp) :: U(M,M), VT(N,N)
  real(dp) :: S(N), RWORK(5*N)
      
  complex(dp), dimension(:), allocatable :: WORK
  INTEGER           :: LDA,LDU,M,N,LWORK,LDVT,INFO
  CHARACTER         ::  JOBU, JOBVT


  JOBU='A'
  JOBVT='A'
  LDA=M
  LDU=M
  LDVT=N

  LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))

  ALLOCATE(work(lwork))

  CALL ZGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK,RWORK, INFO )

end subroutine svd

!___________________________________________________________________________
!
!          Pseudoinverse
!___________________________________________________________________________

function pinv(A) result(Ainv)
  complex(dp), dimension(:,:), intent(in) :: A
  complex(dp), dimension(:,:), allocatable :: Ainv

  complex(dp), dimension(:,:), allocatable :: U, V, SU
  real(dp), dimension(:), allocatable :: S
  integer :: M, N, i1

  M= size(A,1)
  N= size(A,2)

  allocate(U(M,M), V(N,N), S(N))
  allocate(SU(N,M))
  allocate(Ainv(N,M))

  call svd(A,U,S,V,M,N)

  do i1 = 1,N
    if(s(i1) > 1e-6) then 
        SU(i1,:) = 1/s(i1) * conjg(U(:,i1))
    else 
        SU(i1,:) = 0.0 * conjg(U(:,i1)) 
    end if 
    
  end do

  Ainv = matmul(transpose(conjg(V)),SU)
end function pinv

end module possu
