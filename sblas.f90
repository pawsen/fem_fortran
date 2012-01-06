module sblas
! module containing drivers to NIST sparse blas routines, downloaded from
! http://math.nist.gov/~KRemington/fspblas/
! This version only have fortran bindings

! Notice that
! http://math.nist.gov/spblas/original.html
! is also an implementation of sparse NIST with fortran AND c binding. This version, however, seems to support fewer types of storage of the sparse matrix.

  implicit none

  public :: sparse_dense_multiply


contains

  subroutine sparse_dense_multiply(iK_in,jK_in,sK_in, m_in, k_in, b, c)

    ! Multiply a sparse matrix(A) by a dense(B), eg C= A*B
    ! using sparse blas routine from NIST

    !A(m,k)
    !B((k,n))
    !C(m,n)

    integer, intent(in):: iK_in(:), jK_in(:)
    real(8), intent(in) :: sK_in(:), B(:,:)
    integer, intent(in) :: m_in, k_in

    real(8), intent(out) :: C(:,:)


    real(8), allocatable :: work(:,:)
    integer, allocatable :: indx(:), jndx(:)

    integer :: m,k,n,nnz,ldb,ldc,lwork, i,j

    integer :: descra(9), transa
    real(8) :: alpha, beta

    external :: dcoomm

    ! TEST structure
    nnz = size(ik_in,1)

    alpha = 1d0
    beta = 0d0

    m= m_in !numbers of row in sparse matrix(A)
    n= size(C,2) !number of colums in output-matrix(C)
    k= k_in !numbers of colums in sparse matrix(A)
    ldb=k; ldc=m; lwork=22

    allocate(work(m,n))

    !         ------------ begin interface description ------------
    !   Toolkit interface:
    !   dcoomm -- compressed sparse row format matrix-matrix multiply
    !  
    !   C <- alpha A B + beta C
    !  
    !   Arguments:
    !  
    !   int transa	Indicates how to operate with the sparse matrix
    !  		0 : operate with matrix
    !  		1 : operate with transpose matrix
    !  
    !   int m	Number of rows in matrix A
    !  
    !   int n	Number of columns in matrix c
    !  
    !   int k	Number of columns in matrix A
    !  
    !   double alpha Scalar parameter
    !  
    !   double beta  Scalar parameter
    !  
    !   int descra()	Descriptor argument.  Nine element integer array
    !  		descra(1) matrix structure
    !  			0 : general
    !  			1 : symmetric
    !  			2 : Hermitian
    !  			3 : Triangular
    !  			4 : Skew(Anti)-Symmetric
    !  			5 : Diagonal
    !  		descra(2) upper/lower triangular indicator
    !  			1 : lower
    !  			2 : upper
    !  		descra(3) main diagonal type
    !  			0 : non-unit
    !  			1 : unit
    !  		descra(4) Array base 
    !  			0 : C/C++ compatible
    !  			1 : Fortran compatible
    !  		descra(5) repeated indices?
    !  			0 : unknown
    !  			1 : no repeated indices
    !  
    !  
    !
    !   double val()  scalar array of length nnz containing matrix entries.
    !  
    !   int indx()    integer array of length nnz containing row indices.
    !
    !   int jndx()    integer array of length nnz containing column indices.
    !
    !   int nnz       number of non-zero elements in A.
    !
    !   double b()    rectangular array with first dimension ldb.
    !  
    !   double c()    rectangular array with first dimension ldc.
    !  
    !   double work() scratch array of length lwork.  lwork should be at least
    !                 max(m,n) 


    ! Settings for this specific case
    descra = 0 
    descra(1) = 0 ! generel matrix stucture
    descra(4) = 1 !fortran based arrays
    descra(5) = 0 !indicier (kan vÃ¦re) gentaget
    transa = 0 ! not transposed


    call dcoomm( transa, m, n, k, alpha, descra, &
         sK_in, iK_in, jK_in, nnz, &
         B, ldb, beta, c, ldc, work, lwork)



    !     !########## TEST structure ##########
    !     nnz = 9

    !     alpha = 1d0
    !     beta = 0d0

    !     m=4 !numbers of row in sparse matrix(A)
    !     n=4 !number of colum sin output-matrix(C)
    !     k=4 !numbers of colums in sparse matrix(A)
    !     ldb=m; ldc=m; lwork=22

    !     allocate(a(nnz),indx(nnz),jndx(nnz))
    !     allocate(b(k,n),c(m,n), work(m,n))

    !     indx = (/1,1,2,3,3,3,4,4,1/)
    !     jndx = (/1,4,2,1,1,3,2,4,1/)
    !     a = (/0.5,3.0,2.,2.,2.,3.,6.,5.,0.5/)

    !     b(1,:) = [3, 2, 3, 3 ]
    !     b(2,:) = [4, 4, 1, 1 ]
    !     b(3,:) = [4, 3, 4, 3 ]
    !     b(4,:) = [1, 2, 2, 1 ]

    !     ! RESULT:
    !     c(1,:) = [6, 8, 9, 6 ]
    !     c(2,:) = [8, 8, 2, 2 ]
    !     c(3,:) = [24, 7, 4, 21 ]
    !     c(4,:) = [29, 4, 6, 11 ]

    !     Write(*,*)
    !     WRITE(*,*) 'C:'
    !     DO I = 1, size(c,1)
    !        !print*,(b(i,j),j=1,size(b,2))
    !        WRITE(*,9998) ( c( I, J ), J = 1, size(c,2) )
    !     END DO
    ! 9998 FORMAT( 11(:,1X,f6.2) )
    !     !##############################

  end subroutine sparse_dense_multiply


!   subroutine sparse_dense_multiply(iK_in,jK_in,sK_in, m_in, k_in, b, c) 
    
!     ! Multiply a sparse matrix(A) by a dense(B), eg C= A*B
!     ! using sparse blas routine from NIST

!     !A(m,k)
!     !B((k,n))
!     !C(m,n)

!     integer, intent(in):: iK_in(:), jK_in(:)
!     real(8), intent(in) :: sK_in(:), B(:,:)
!     integer, intent(in) :: m_in, k_in

!     real(8), intent(out) :: C(:,:)
    

!     real(8), allocatable :: work(:)
!     integer, allocatable :: indx(:), jndx(:)

!     integer :: m,k,n,nnz,ldb,ldc,lwork, i,j

!     integer :: descra(9), transa
!     real(8) :: alpha, beta

!     external :: dcoomm

!     ! TEST structure
!     nnz = size(ik_in,1)

!     alpha = 1d0
!     beta = 0d0

!     m= m_in !numbers of row in sparse matrix(A)
!     n= size(C,2) !number of colums in output-matrix(C)
!     k= k_in !numbers of colums in sparse matrix(A)
!     ldb=k; ldc=m; lwork=n*m

!     allocate(work(n*m))
!     !allocate(work(MAX(m,n)))

!     ! Settings for this specific case
!     descra = 0 
!     descra(1) = 0 ! generel matrix stucture
!     descra(3) = 0 ! diagonal -non-unit
!     descra(4) = 1 !fortran based arrays
!     descra(5) = 0 !indicier (kan vÃ¦re) gentaget
!     transa = 0 ! not transposed


!     call dcoomm( transa, m, n, k, alpha, descra, &
!          sK_in, iK_in, jK_in, nnz, &
!          B, ldb, beta, c, ldc, work, lwork)


!     !         ------------ begin interface description ------------
!     !   Toolkit interface:
!     !   dcoomm -- compressed sparse row format matrix-matrix multiply
!     !  
!     !   C <- alpha A B + beta C
!     !  
!     !   Arguments:
!     !  
!     !   int transa	Indicates how to operate with the sparse matrix
!     !  		0 : operate with matrix
!     !  		1 : operate with transpose matrix
!     !  
!     !   int m	Number of rows in matrix A
!     !  
!     !   int n	Number of columns in matrix c
!     !  
!     !   int k	Number of columns in matrix A
!     !  
!     !   double alpha Scalar parameter
!     !  
!     !   double beta  Scalar parameter
!     !  
!     !   int descra()	Descriptor argument.  Nine element integer array
!     !  		descra(1) matrix structure
!     !  			0 : general
!     !  			1 : symmetric
!     !  			2 : Hermitian
!     !  			3 : Triangular
!     !  			4 : Skew(Anti)-Symmetric
!     !  			5 : Diagonal
!     !  		descra(2) upper/lower triangular indicator
!     !  			1 : lower
!     !  			2 : upper
!     !  		descra(3) main diagonal type
!     !  			0 : non-unit
!     !  			1 : unit
!     !  		descra(4) Array base 
!     !  			0 : C/C++ compatible
!     !  			1 : Fortran compatible
!     !  		descra(5) repeated indices?
!     !  			0 : unknown
!     !  			1 : no repeated indices
!     !  
!     !  
!     !
!     !   double val()  scalar array of length nnz containing matrix entries.
!     !  
!     !   int indx()    integer array of length nnz containing row indices.
!     !
!     !   int jndx()    integer array of length nnz containing column indices.
!     !
!     !   int nnz       number of non-zero elements in A.
!     !
!     !   double b()    rectangular array with first dimension ldb.
!     !  
!     !   double c()    rectangular array with first dimension ldc.
!     !  
!     !   double work() scratch array of length lwork.  lwork should be at least
!     !                 max(m,n) 
    

! !     !########## TEST structure ##########
! !     nnz = 9

! !     alpha = 1d0
! !     beta = 0d0

! !     m=4 !numbers of row in sparse matrix(A)
! !     n=4 !number of colum sin output-matrix(C)
! !     k=4 !numbers of colums in sparse matrix(A)
! !     ldb=m; ldc=m; lwork=22

! !     allocate(a(nnz),indx(nnz),jndx(nnz))
! !     allocate(b(k,n),c(m,n), work(m,n))

! !     indx = (/1,1,2,3,3,3,4,4,1/)
! !     jndx = (/1,4,2,1,1,3,2,4,1/)
! !     a = (/0.5,3.0,2.,2.,2.,3.,6.,5.,0.5/)

! !     b(1,:) = [3, 2, 3, 3 ]
! !     b(2,:) = [4, 4, 1, 1 ]
! !     b(3,:) = [4, 3, 4, 3 ]
! !     b(4,:) = [1, 2, 2, 1 ]

! !     ! RESULT:
! !     c(1,:) = [6, 8, 9, 6 ]
! !     c(2,:) = [8, 8, 2, 2 ]
! !     c(3,:) = [24, 7, 4, 21 ]
! !     c(4,:) = [29, 4, 6, 11 ]

! !     Write(*,*)
! !     WRITE(*,*) 'C:'
! !     DO I = 1, size(c,1)
! !        !print*,(b(i,j),j=1,size(b,2))
! !        WRITE(*,9998) ( c( I, J ), J = 1, size(c,2) )
! !     END DO
! ! 9998 FORMAT( 11(:,1X,f6.2) )
! !     !##############################

!   end subroutine sparse_dense_multiply


end module sblas
