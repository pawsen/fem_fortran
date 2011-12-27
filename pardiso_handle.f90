! this has to come before the module
INCLUDE 'mkl_pardiso.f90'

module pardiso_handle


  implicit NONE

  PRIVATE
  PUBLIC :: pardiso_format, pardiso_solve

CONTAINS


subroutine mumps

  INCLUDE 'mpif.h'
  INCLUDE 'dmumps_struc.h'
  TYPE (DMUMPS_STRUC) id
  INTEGER IERR, I

  integer :: n, nz, ia(12),ja(12)
  real(8) :: A(12), b(5)
  
  
  ! input info
  n = 5
  nz = 12
  ia = \(1,2,4,5,2,1,5,3,2,3,1,3\)
  ja = \(2,3,3,5,1,1,2,4,5,2,3,3\)
  A = \(3,-3,2,1,3,2,4,2,6,-1,4,1\)
  b = \(20,24,9,6,13\)

!!$1 2 3.0
!!$2 3 -3.0
!!$4 3 2.0
!!$5 5 1.0
!!$2 1 3.0
!!$1 1 2.0
!!$5 2 4.0
!!$3 4 2.0
!!$2 5 6.0
!!$3 2 -1.0
!!$1 3 4.0
!!$3 3 1.0        :values
!!$20.0
!!$24.0
!!$9.0
!!$6.0
!!$13.0           :RHS

  CALL MPI_INIT(IERR)
  !Define a communicator for the package
  id%COMM = MPI_COMM_WORLD
  !Ask for unsymmetric code
  id%SYM = 0
  !Host working
  id%PAR = 1
  !Initialize an instance of the package
  id%JOB = -1
  CALL DMUMPS(id)
  !Define problem on the host (processor 0)
  IF ( id%MYID .eq. 0 ) THEN
     READ(5,*) id%N
     READ(5,*) id%NZ
     ALLOCATE( id%IRN ( id%NZ ) )
     ALLOCATE( id%JCN ( id%NZ ) )
     ALLOCATE( id%A( id%NZ ) )
     ALLOCATE( id%RHS ( id%N ) )
     READ(5,*) ( id%IRN(I) ,I=1, id%NZ )
     READ(5,*) ( id%JCN(I) ,I=1, id%NZ )
     READ(5,*) ( id%A(I),I=1, id%NZ )
     READ(5,*) ( id%RHS(I) ,I=1, id%N )
  END IF
  Call package for solution
  id%JOB = 6
  CALL DMUMPS(id)
  !Solution has been assembled on the host
  IF ( id%MYID .eq. 0 ) THEN
     WRITE( 6, * ) ' Solution is ',(id%RHS(I),I=1,id%N)
  END IF
  Deallocate user data
  IF ( id%MYID .eq. 0 )THEN
     DEALLOCATE( id%IRN )
     DEALLOCATE( id%JCN )
     DEALLOCATE( id%A
     )
     DEALLOCATE( id%RHS )
  END IF
  !Destroy the instance (deallocate internal data structures)
  id%JOB = -2
  CALL DMUMPS(id)
  CALL MPI_FINALIZE(IERR)
  STOP
END


end subroutine mumps


!  subroutine pardiso_format(iK_out,jK_out,sK_out,iK_ind,jK_ind,sK_ind)
  subroutine pardiso_format(iK_in,jK_in,sK_in)

    
    USE mkl_pardiso
    use fedata

    ! parametre for pardiso interface
    INTEGER :: job(8), n, nnz, info,i
!    integer, allocatable, dimension(:) , intent(out) :: iK_out, jK_out
!    real(8), allocatable, dimension(:) , intent(out) :: sK_out
    integer,dimension(:), intent(in) :: iK_in,jK_in
    real(8),dimension(:), intent(in) :: sK_in

    integer, allocatable, dimension(:) :: iK_out, jK_out
    real(8), allocatable, dimension(:) :: sK_out

    ! n = dimension of A = neqn, nnz = number non-zero,
    n = neqn
    allocate(iK(n+1))
    allocate(iK_out(n+1),sK_out(n),jK_out(n))
 
    data(job(1) = 1) ! This setting indicates that we want to convert a coordinate/triplet based sparse to 
    data(job(2:3) = 1)  ! Both input and output matrices are 1-based (not 0-bases as in c).  
 !   job(5) = neqn
    data(job(6)=1)

do i=1,n
   print*,'sk',sK_out(i)
end do


    ! mkl_dcsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
!    call mkl_dcsrcoo (job, n, sK_out, jK_out,iK_out,nnz, sK_ind, iK_ind,jK_ind,info)
!    call mkl_dcsrcoo (job, n, sK_in, jK_in,iK_in,nnz, sK, iK,jK,info)
   call mkl_dcsrcoo (job, n, sK_out, jK_out,iK_out,nnz, sK_in, iK_in,jK_in,info)


    print*,'info =',info

  end subroutine pardiso_format

  subroutine pardiso_solve(iK,jK,sK)
!    INCLUDE 'mkl_pardiso.f90'
    USE mkl_pardiso
    use fedata

    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    !.. Internal solver memory pointer 
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
    !.. All other variables

    INTEGER :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
    INTEGER :: error1
    INTEGER, ALLOCATABLE :: iparm( : )
    INTEGER, dimension(:), intent(in) :: iK,jK
    real(8), dimension(:), intent(in) :: sK !,p
    !real(8), dimension(:), intent(in) :: d

    
    INTEGER i, idum(1)
    REAL(8) ddum(1)
    !.. Fill all arrays containing matrix data.
    n = size(jK,1)!neqn
    nnz = size(sK,1)

    nrhs = 1 
    maxfct = 1 
    mnum = 1

    !..
    !.. Set up PARDISO control parameter
    !..
    ALLOCATE( iparm ( 64 ) )
    do i = 1, 64
       iparm(i) = 0
    end do

    iparm(1) = 1 ! no solver default
    iparm(2) = 2 ! fill-in reordering from METIS
    iparm(4) = 0 ! no iterative-direct algorithm
    iparm(5) = 0 ! no user fill-in reducing permutation
    iparm(6) = 0 ! =0 solution on the first n compoments of x
    iparm(8) = 9 ! numbers of iterative refinement steps
    iparm(10) = 13 ! perturbe the pivot elements with 1E-13
    iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
    iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
    iparm(14) = 0 ! Output: number of perturbed pivots
    iparm(18) = -1 ! Output: number of nonzeros in the factor LU
    iparm(19) = -1 ! Output: Mflops for LU factorization
    iparm(20) = 0 ! Output: Numbers of CG Iterations

    error  = 0 ! initialize error flag
    msglvl = 0 ! print statistical information
    mtype  = 11! -2 ! symmetric, indefinite

    !.. Initiliaze the internal solver memory pointer. This is only
    ! necessary for the FIRST call of the PARDISO solver.

    ALLOCATE ( pt ( 64 ) )
    do i = 1, 64
       pt( i )%DUMMY =  0 
    end do

    !.. Reordering and Symbolic Factorization, This step also allocates
    ! all memory that is necessary for the factorization

    phase = 11 ! only reordering and symbolic factorization

    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, sK, iK, jK, &
         idum, nrhs, iparm, msglvl, ddum, ddum, error)

    WRITE(*,*) 'Reordering completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       GOTO 1000
    END IF
    WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
    WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

    !.. Factorization.
    phase = 22 ! only factorization
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, sK, iK, jK, &
         idum, nrhs, iparm, msglvl, ddum, ddum, error)
    WRITE(*,*) 'Factorization completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       GOTO 1000
    ENDIF

    !.. Back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only factorization
    
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, sK, iK, jK, &
         idum, nrhs, iparm, msglvl, p, d, error)
    WRITE(*,*) 'Solve completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       GOTO 1000
    ENDIF
    WRITE(*,*) 'The solution of the system is '
    DO i = 1, n
       WRITE(*,*) ' d(',i,') = ', d(i)
    END DO

1000 CONTINUE
    !.. Termination and release of memory
    phase = -1 ! release internal memory
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
         idum, nrhs, iparm, msglvl, ddum, ddum, error1)
!!$
!!$    IF ( ALLOCATED( ia ) )      DEALLOCATE( ia )
!!$    IF ( ALLOCATED( ja ) )      DEALLOCATE( ja )
!!$    IF ( ALLOCATED( a ) )       DEALLOCATE( a )
!!$    IF ( ALLOCATED( b ) )       DEALLOCATE( b )
!!$    IF ( ALLOCATED( x ) )       DEALLOCATE( x )
!!$    IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )

    IF (error1 /= 0) THEN
       WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
       STOP 1
    ENDIF

    IF ( error /= 0 ) STOP 1
    STOP 0


end subroutine pardiso_solve

end module pardiso_handle
