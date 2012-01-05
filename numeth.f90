module numeth

   implicit none

   integer, dimension(:), allocatable :: work
   save

   private
   public :: factor, solve, bfactor, bsolve, vector_mul, gauss_points
   public :: index_dups, rotation_matrix, sparse_multiply
   public :: correct_bound_vector, get_edof
   PUBLIC :: lapack_eigenval, lapack_solve_general, lapack_inverse_matrix
   public :: blas_matmul

   public :: make_dir

   !print to screen
   PUBLIC :: display_matrix

contains

   subroutine factor(a)

   ! This subroutine factors a general matrix [a].
 
   integer :: i, j, k, neqn
   real(8), dimension(:, :), intent(inout) :: a

   neqn = size(a, 1)

   do j = 2, neqn
      do i = 1, j-1
         a(j, i) = a(j, i) - dot_product(a(j, 1:i-1), a(1:i-1, i))
         a(i, j) = a(i, j) - dot_product(a(i, 1:i-1), a(1:i-1, j))
         a(j, i) = a(j, i)/a(i, i)
      end do
      a(j, j) = a(j, j) - dot_product(a(j, 1:j-1), a(1:j-1, j))
   end do

end subroutine factor

subroutine solve(a, b)

   ! This subroutine solves [a]{x} = {b} using the previously factored
   ! coefficient matrix [a] from subroutine 'factor'.
   ! The subroutine returns the solution {x} by overwriting b.

   integer :: i, j, k, neqn
   real(8), dimension(:, :), intent(in) :: a
   real(8), dimension(:), intent(inout) :: b

   neqn = size(a, 1)

   ! Forward substitution
   do i = 2, neqn
      b(i) = b(i) - dot_product(a(i, 1:i-1), b(1:i-1))
   end do
 
   ! Backward substitution
   b(neqn) = b(neqn)/a(neqn, neqn)
   do i = neqn-1, 1, -1
      b(i) = (b(i) - dot_product(a(i, i+1:neqn), b(i+1:neqn)))/a(i, i)
   end do

end subroutine solve
 
subroutine bfactor(a)

   ! This subroutine factors a matrix [a] stored in banded form.

   integer :: i, j, n, l, k, bw, neqn
   real(8) :: c
   real(8), dimension(:, :), intent(inout) :: a

   bw = size(a, 1)
   neqn = size(a, 2)

   do n = 1, neqn
      do l = 2, bw
         if (a(l, n) == 0.) cycle
         i = n+l-1
         c = a(l, n)/a(1, n)
         j = 0
      do k = l, bw
         j = j+1
         a(j, i) = a(j, i) - c*a(k, n)
      end do
      a(l, n) = c
      end do
   end do 

end subroutine bfactor

subroutine bsolve(a, b)

   ! This subroutine solves [a]{x} = {b} using the previously factored
   ! coefficient matrix [a] from subroutine 'bfactor'.
   ! The subroutine returns the solution {x} by overwriting b.

   integer :: i, n, l, m, k, bw, neqn
   real(8), dimension(:, :), intent(inout) :: a
   real(8), dimension(:) :: b

   bw = size(a, 1)
   neqn = size(a, 2)

   do n = 1,neqn
      do l = 2, bw
         if (a(l, n) == 0.) cycle 
         i = n+l-1
         b(i) = b(i)-a(l, n)*b(n)
      end do
      b(n) = b(n)/a(1, n)
   end do

   do m = 2, neqn
      n = neqn+1-m
      do l = 2, bw
         if (a(l, n) == 0.) cycle 
         k = n+l-1
         b(n) = b(n)-a(l, n)*b(k)
      end do
   end do

end subroutine bsolve

subroutine vector_mul(A,B,C)
  ! ganger to vektorer med hinanden, således resultatet bilver en matrix

  real(8),dimension(:), intent(IN) :: A, B
  real(8),dimension(:,:), intent(OUT) :: C
  integer :: i,j

  do i=1,size(A,1)
     do j=1,size(B,1)
        C(i,j) = A(i)*B(j)
     end do
  end do

end subroutine vector_mul

SUBROUTINE gauss_points(w,xi,eta,ng)

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), DIMENSION(:) , INTENT(out):: w, xi, eta

  SELECT CASE ( ng )
  case (1)
     w = 2.0
     xi = 0.0
     eta =xi
  case (2)
     w = 1.0
     xi(1) = -3.0**(-0.5)
     xi(2) = 3.0**(-0.5)
     eta = xi
  case (3)
     w(1) = 5.0/9.0
     w(2) = 8.0/9.0
     w(3) = 5.0/9.0
     xi(1) = -0.6**(0.5)
     xi(2) = 0.0
     xi(3) = 0.6**(0.5)
     eta = xi
  case (4)  
     w(1) = (18.0+dsqrt(30d0))/36.0
     w(2) = (18.0+dsqrt(30d0))/36.0
     w(3) = (18.0-dsqrt(30d0))/36.0
     w(4) = (18.0-dsqrt(30d0))/36.0
     xi(1) = dsqrt((3d0-2d0*dsqrt(6d0/5d0))/7d0)
     xi(2) = -dsqrt((3d0-2d0*dsqrt(6d0/5d0))/7d0)
     xi(3) = dsqrt((3d0+2d0*dsqrt(6d0/5d0))/7d0)
     xi(4) = -dsqrt((3d0+2d0*dsqrt(6d0/5d0))/7d0)
     eta = xi
  case default
     print *, 'Error in number of gauss points'
     print *, 'ng = [1;4]'
     stop
  END SELECT


END SUBROUTINE gauss_points

  subroutine index_dups(in_vec,out_vec,k)
    !return the index of non-dups indicies(nodes)

    implicit none
    integer, INTENT(IN) :: in_vec(:,:)
    integer, INTENT(OUT) :: out_vec(:),k
    integer, allocatable :: buffer(:,:)
    !integer :: k                   ! The number of unique elements
    integer :: i, j

    allocate(buffer(size(in_vec,1),2))

    k = 1
    buffer(1,:) = in_vec(1,:)
    out_vec(1) = 1
    outer: do i=2,size(in_vec,1)
       do j=1,k
          if ((buffer(j,1) == in_vec(i,1)) .and. &
               (buffer(j,2) == in_vec(i,2)) ) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       buffer(k,:) = in_vec(i,:)
       out_vec(k) = i
    end do outer

!!$    write(*,advance='no',fmt='(a,i0,a)') 'Unique list has ',k,' elements: '
!!$    write(*,*) buffer(1:k,1)
  end subroutine index_dups

  subroutine rotation_matrix(phi,T_p,T_t,lambda)
    !this subroutine calculates the rotation matrix

    real(8), intent(IN) :: phi
    real(8), intent(out) :: T_p(3,3), T_t(2,2), lambda(2,2)
    real(8) :: T(6,6), c, s

    c = DCOS(phi)
    s = DSIN(phi)

    T = 0d0
    T(1,1) = c**2
    T(1,2) = s**2
    T(1,3) = c*s
    T(2,1) = s**2
    T(2,2) = c**2
    T(2,3) = -c*s
    T(3,1) = -2d0*c*s
    T(3,2) = 2d0*c*s
    T(3,3) = c**2 - s**2
    T(4,4) = c
    T(4,5) = s
    T(5,4) = -s
    T(5,5) = c
    T(6,6) = 1d0

    lambda = 0d0
    lambda(1,1) = c
    lambda(1,2) = s
    lambda(2,1) = -s
    lambda(2,2) = c
    !lambda(3,3) = 1

    T_p = T(1:3,1:3)
    T_t = T(4:5,4:5)

  end subroutine rotation_matrix

  subroutine make_dir

    use fedata

    CHARACTER :: delimiter
    logical :: exist

    !CALL get_environment_variable('DELIMITER',delimiter) ! find ud af om '/' eller '\' skal bruges ved mapper
    ! kræver at variablen DELIMITER eksistere. Fx på linux: export DELIMITER='/'
    !find out if folder exist by checking for . in the folder. Works on windows/unix/mac
    !http://www.rhinocerus.net/forum/lang-fortran/354299-how-determine-whether-not-there-exists-directory.html

    !INQUIRE (DIRECTORY=trim(filename)//'dir',EXIST=exist) ! virker med ifort, men ikke gfortran
    INQUIRE(file='./'//trim(filename)//'dir/'//'.',EXIST=exist)
    !INQUIRE(file=trim(filename)//'dir'//delimiter//'.',EXIST=exist)

    ! print*
    ! print*,delimiter
    ! print*

    if (.not. exist) then
       CALL system('mkdir '//trim(filename)//'dir') ! make new directory
    end if

  end subroutine make_dir

  subroutine sparse_multiply(row,col,val,xx,b,correct_bound)!,nb,bound)
    !multiply a sparse matrix, A, with a dense vector,x,, eg
    ! A*x = b

    real(8), dimension(:), intent(IN) :: val,xx
    integer, dimension(:), intent(IN) :: row,col
    real(8), dimension(:), intent(OUT) :: b
    logical, optional, intent(IN) :: correct_bound

    integer :: i,ii

    b = 0d0
    do ii=1,size(row,1)
       i = row(ii)
       b(i) = b(i) + val(ii)*xx(col(ii))
    end do

    ! Arggh ikke pænt. HVAD laver nedenstående i min fine matrix-rutine. ARGHHH!
    if (present(correct_bound)) then
       if ( correct_bound .eqv. .true.) then
          call correct_bound_vector(b)
       end if
    end if

!!$    !! test structure
!!$    real(8) :: val(9),xx(4),b(4)
!!$    integer :: row(9),col(9)
!!$    
!!$    row = (\1,1,2,3,3,3,4,4,1\)
!!$    col = (\1,4,2,1,1,3,2,4,1\)
!!$    val = (\0.5,3.0,2.,2.,2.,3.,6.,5.,0.5\)
!!$    xx = (\2.,4.,6.,8.\)
!!$
!!$    Call sparse_multiply(row,col,val,xx,b)
!!$    
!!$    write(*,'("b = ")')
!!$    do i=1,4
!!$       write (*, '(i6, 1x, f25.15)') i, b(i)
!!$    end do
!!$
!!$    !!A.x=b
!!$    !!   | 1     0     0     3| |2| |26|
!!$    !!   | 0     2     0     0| |4| |8 |
!!$    !!A: | 4     0     3     0|*|6|=|26|
!!$    !!   | 0     6     0     5| |8| |64|



  end subroutine sparse_multiply

  subroutine correct_bound_vector(vec)
  
    use fedata
    real(8), intent(INOUT) :: vec(:)
    integer :: i,idof, kk
    
    select case( element(1)%id )
    case(2)
       kk = 2
    case default
       kk = 3
    end select
    do i = 1, nb ! OUTVEC corrected for BC
       if  (bound(i,2) < 6) then! forskydning
          idof = kk*(bound(i,1)-1) + bound(i,2)
          vec(idof) = 0
       end if
    end do


  end subroutine correct_bound_vector


  Function lapack_inverse_matrix(a) result(ainv) 
    implicit none

    ! Bruger LAPACK
    !Fra http://www.megasolutions.net/fortran/Inverse-of-a-square-matrix-35510.aspx
    !En anden kode der ikke bruger lapack er http://www.dreamincode.net/code/snippet1272.htm

    integer, parameter :: dp = kind(1.0d0) 
    Real(kind=dp), Intent(in):: a(:,:) 
    real(kind=dp) :: ainv(size(a,1),size(a,2)) 
    Real(kind=dp), Allocatable:: WORK(:) 
    Integer, Allocatable:: IPIV(:) 
    Integer:: m, N, Info, LWORK 

    !########## test-struktur ##########
    !mat(1,:) = [1. , 5. , 2.]
    !mat(2,:) = [1. , 1. , 7.]
    !mat(3,:) = [0. , -3. , 4.]
    ! skal resultere i
    !inv_mat(1,:) = [-25. , 26. , -33.]
    !inv_mat(2,:) = [ -4. , -4. ,   5.]
    !inv_mat(3,:) = [  3. , -3. ,   4.]
    !print*,'inv_mat'
    !do i=1,3
    !   print*,(inv_mat(i,j),j=1,3)
    !end do
    !###################################

    m=Size(a,Dim=1) 
    N=Size(a,Dim=2) 
    If (m/=N) Then 
       Ainv=0.0_dp 
    Else 
       Allocate (IPIV(N)) 
       Ainv=A 
       Call DGETRF(N,N,Ainv,N,IPIV,Info) 
       If(Info == 0) Then 
          Allocate (WORK(N)) 
          LWORK=-1 
          Call DGETRI(N,Ainv,N,IPIV,WORK,LWORK,Info) 
          LWORK=WORK(1) 
          DeAllocate (WORK) 
          Allocate (WORK(LWORK)) 
          Call DGETRI(N,Ainv,N,IPIV,WORK,LWORK,Info) 
          If (Info /= 0) Ainv=0 
          DeAllocate (WORK) 
       Elseif(Info >0) then !if INFO = i:
          Ainv=0 
          print*
          print*,'Error in inv_matrix'
          print*,'U(i,i) is exactly zero. The factorization'
          print*,'has been completed, but the factor U is exactly'
          print*,'singular, and division by zero will occur if it is used'
          print*,'to solve a system of equations., i: ', info
          print*
          print*,'På dansk: Inputmatricen er singulær, numeth.f90!'
          error stop
       else
          Ainv=0 
          print*
          print*,'Error in inv_matrix'
          print*,'the i-th argument had an illegal value, -i: ', info
          print*,'numeth.f90'
          error stop
       EndIf
       DeAllocate (IPIV) 
    EndIf
  End Function lapack_inverse_matrix


  subroutine get_edof(e,nen,xe,edof41,edof)

    use fedata

    integer, intent(in) ::  e
    integer, intent(out) :: edof(:), nen
    integer, intent(out), optional :: edof41(:)
    real(8), intent(out) :: xe(:)
    integer :: i

    nen = element(e)%numnode
    do i = 1, nen
       xe(2*i-1) = x(element(e)%ix(i),1)
       xe(2*i  ) = x(element(e)%ix(i),2)
       if (present(edof41)) then
          edof41(i) = element(e)%ix(i)
       end if

       if ((element(1)%id == 2 )) then ! plane42 element
          edof(2*i-1) = 2 * element(e)%ix(i) - 1
          edof(2*i)   = 2 * element(e)%ix(i)
       else ! plate
          edof(3*i-2) = 3 * element(e)%ix(i) - 2 
          edof(3*i-1) = 3 * element(e)%ix(i) - 1
          edof(3*i)   = 3 * element(e)%ix(i)
       end if
    end do

  end subroutine get_edof


  subroutine lapack_eigenval(mat,eigenval)
    ! this subroutine calculates eigenvalue and vectors of a symmetric matrix (standard eigenvalue problem). It uses dsyev.f from LAPACK
    
    !
    !  DSYEV Example.
    !  ==============
    ! From http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/dsyev_ex.f.htm
    !
    !  Program computes all eigenvalues and eigenvectors of a real symmetric
    !  matrix A:
    !
    !  mat(1,:) = [ 1.96,-6.49, -0.47, -7.20, -0.65]
    !  mat(2,:) = [-6.49, 3.80, -6.39,  1.50, -6.34]
    !  mat(3,:) = [-0.47,-6.39,  4.17, -1.51,  2.67]
    !  mat(4,:) = [-7.20, 1.50, -1.51,  5.70,  1.80]
    !  mat(5,:) = [-0.65,-6.34,  2.67,  1.80, -7.10]
    !
    !  Description.
    !  ============
    !
    !  The routine computes all eigenvalues and, optionally, eigenvectors of an
    !  n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
    !
    !  A*v(j) = lambda(j)*v(j)
    !
    !  where lambda(j) is its eigenvalue. The computed eigenvectors are
    !  orthonormal.
    !
    !  Example Program Results.
    !  ========================
    !
    ! DSYEV Example Program Results
    !
    ! Eigenvalues
    ! -11.07  -6.23   0.86   8.87  16.09
    !
    ! Eigenvectors (stored columnwise)
    !  -0.30  -0.61   0.40  -0.37   0.49
    !  -0.51  -0.29  -0.41  -0.36  -0.61
    !  -0.08  -0.38  -0.66   0.50   0.40
    !   0.00  -0.45   0.46   0.62  -0.46
    !  -0.80   0.45   0.17   0.31   0.16

    real(8), intent(inout) :: mat(:,:)
    real(8), intent(out) :: eigenval(:)
    !real(8) :: mat(5,5), eigenval(5), eigenvec(5,5)

    external :: DSYEV
    integer :: n, lda, lwork,info,i,j
    character :: jobz*1, uplo*1
    integer,parameter :: LWMAX = 1000
    real(8) :: work(LWMAX)
    real(8), allocatable :: A(:,:), W(:)


    jobz = 'V' !compute both eigenvalues and vectors
    uplo = 'U' !upper triangle of A is stored
    N    = size(mat,2)   !dimension of A 
    LDA  = N   !The leading dimension of the array A.  LDA >= max(1,N)
    LWORK= -1
    
    allocate(A(N,N),W(N))
    A = mat
    
    ! If mat i not completely symmetric(due to numerical errors), we can do:
    !A = 0.5d0 * (A + TRANSPOSE(A))

    ! Make A as an upper diagonal matrix
    do i=1,N
       do j=1,N
          if (i>j) then
             A(i,j) = 0d0
          end if
       end do
    end do

    !call DISPLAY_MATRIX( 'eigenvec', A )

    !Query the optimal workspace.
    call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

    !Solve eigenproblem.
    call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

    if (info < 0) then
       print*,'if INFO = -i, the i-th argument had an illegal value'
       error stop
    elseif (info > 0) then
       print*,' if INFO = i, the algorithm failed to converge; i'
       print*,'         off-diagonal elements of an intermediate tridiagonal'
       print*,'         form did not converge to zero.'
       error stop
    end if

    !call DISPLAY_MATRIX( 'eigenvec', A )

    eigenval = W
    mat = A
    !eigenvec = A
    
  end subroutine lapack_eigenval

  SUBROUTINE DISPLAY_MATRIX( title, A )

    character(len=*), intent(in) :: title
    REAL(8), intent(in) ::A(:,:)
    INTEGER:: i,j,ii

    ii = size(A,1)

    WRITE(*,*)
    WRITE(*,*) title
    DO I = 1, ii
       print*,(A(i,j),j=1,size(A,2))
       !WRITE(*,9998) ( A( I, J ), J = 1, size(A,1) )
    END DO

9998 FORMAT( 11(:,1X,f6.2) )
    RETURN

    ! !Another way to print the matrix
    ! do i=1,12
    !    print*,(ke_piezo(i,j),j=1,4)
    !    !write (*,'(8(D10.4 2x))') (real(ke_stiff(i,j)),j=1,8)
    ! end do

  END SUBROUTINE DISPLAY_MATRIX

  subroutine lapack_solve_general(A,B)

    !solves Ax = b
    ! Can solve for multiple rhs. Change nrhs

    real(8), intent(in) :: A(:,:)
    real(8), intent(inout) :: B(:)

    integer :: n, i, j, info
    integer :: lda, ldb
    integer, allocatable :: ipiv(:)
    integer, parameter :: nrhs = 1

    external :: dgesv

    n = size(A,1)
    lda = n
    ldb = n
    allocate(ipiv(n))
    
    CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )


    ! Check for the exact singularity.
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The diagonal element of the triangular factor of A,'
       WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
       WRITE(*,*)'A is singular; the solution could not be computed.'
       error STOP
    elseif (info < 0) then
       write(*,*)'Fejl i lapack_solve_general'
       write(*,*)'the i-th argument had an illegal value'
       error stop
    END IF

    !  Purpose
    !  =======
    !
    !  DGESV computes the solution to a real system of linear equations
    !     A * X = B,
    !  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    !
    !  The LU decomposition with partial pivoting and row interchanges is
    !  used to factor A as
    !     A = P * L * U,
    !  where P is a permutation matrix, L is unit lower triangular, and U is
    !  upper triangular.  The factored form of A is then used to solve the
    !  system of equations A * X = B.
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The number of linear equations, i.e., the order of the
    !          matrix A.  N >= 0.
    !
    !  NRHS    (input) INTEGER
    !          The number of right hand sides, i.e., the number of columns
    !          of the matrix B.  NRHS >= 0.
    !
    !  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    !          On entry, the N-by-N coefficient matrix A.
    !          On exit, the factors L and U from the factorization
    !          A = P*L*U; the unit diagonal elements of L are not stored.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,N).
    !
    !  IPIV    (output) INTEGER array, dimension (N)
    !          The pivot indices that define the permutation matrix P;
    !          row i of the matrix was interchanged with row IPIV(i).
    !
    !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
    !          On entry, the N-by-NRHS matrix of right hand side matrix B.
    !          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
    !
    !  LDB     (input) INTEGER
    !          The leading dimension of the array B.  LDB >= max(1,N).
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    !                has been completed, but the factor U is exactly
    !                singular, so the solution could not be computed.

  end subroutine lapack_solve_general

  subroutine blas_matmul(A,B,C)
    ! multiply two matrices to get C=A*B by using blas/atlas routine. Should be better than matmul for large arrays

    ! http://www.netlib.org/blas/dgemm.f
    ! http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/cpp/win/mkl/refman/bla/functn_gemm.html
    
    real(8), intent(in) :: A(:,:), B(:,:)
    real(8), intent(inout) :: C(:,:)

    real(8) :: alpha, beta
    integer :: m, n, k, LDA, LDB, LDC
    character:: TRANSA*1, TRANSB*1

    external :: DGEMM

    
    TRANSA = 'N'
    TRANSB = 'N'
    m = size(A,1)
    n = size(B,2)
    k = size(A,2)
    lda = m
    ldb = k
    ldc = m
    alpha = 1d0
    beta = 0d0


    call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

    !  Arguments
    !  ==========
    !
    !  TRANSA - CHARACTER*1.
    !           On entry, TRANSA specifies the form of op( A ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSA = 'N' or 'n',  op( A ) = A.
    !
    !              TRANSA = 'T' or 't',  op( A ) = A**T.
    !
    !              TRANSA = 'C' or 'c',  op( A ) = A**T.
    !
    !           Unchanged on exit.
    !
    !  TRANSB - CHARACTER*1.
    !           On entry, TRANSB specifies the form of op( B ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSB = 'N' or 'n',  op( B ) = B.
    !
    !              TRANSB = 'T' or 't',  op( B ) = B**T.
    !
    !              TRANSB = 'C' or 'c',  op( B ) = B**T.
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry,  M  specifies  the number  of rows  of the  matrix
    !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry,  N  specifies the number  of columns of the matrix
    !           op( B ) and the number of columns of the matrix C. N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  K      - INTEGER.
    !           On entry,  K  specifies  the number of columns of the matrix
    !           op( A ) and the number of rows of the matrix op( B ). K must
    !           be at least  zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
    !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !           part of the array  A  must contain the matrix  A,  otherwise
    !           the leading  k by m  part of the array  A  must contain  the
    !           matrix A.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !           least  max( 1, k ).
    !           Unchanged on exit.
    !
    !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
    !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !           part of the array  B  must contain the matrix  B,  otherwise
    !           the leading  n by k  part of the array  B  must contain  the
    !           matrix B.
    !           Unchanged on exit.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !           least  max( 1, n ).
    !           Unchanged on exit.
    !
    !  BETA   - DOUBLE PRECISION.
    !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !           supplied as zero then C need not be set on input.
    !           Unchanged on exit.
    !
    !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
    !           Before entry, the leading  m by n  part of the array  C must
    !           contain the matrix  C,  except when  beta  is zero, in which
    !           case C need not be set on entry.
    !           On exit, the array  C  is overwritten by the  m by n  matrix
    !           ( alpha*op( A )*op( B ) + beta*C ).
    !
    !  LDC    - INTEGER.
    !           On entry, LDC specifies the first dimension of C as declared
    !           in  the  calling  (sub)  program.   LDC  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.


  end subroutine blas_matmul
  

end module numeth

