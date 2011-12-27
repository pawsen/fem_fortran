module numeth

   implicit none

   integer, dimension(:), allocatable :: work
   save

   private
   public :: factor, solve, bfactor, bsolve, vector_mul, gauss_points
   public :: index_dups, rotation_matrix, sparse_multiply, inv_matrix
   public :: correct_bound_vector, get_edof

   public :: make_dir

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

    CALL get_environment_variable('DELIMITER',delimiter) ! find ud af om '/' eller '\' skal bruges ved mapper
    !find out if folder exist by checking for . in the folder. Works on windows/unix/mac
    !INQUIRE(file=trim(filename)//'dir/'//'.',EXIST=exist)
    INQUIRE(file=trim(filename)//'dir'//delimiter//'.',EXIST=exist)

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
    if (present(correct_bound) .and. correct_bound .eqv. .true.) then
       call correct_bound_vector(b)
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


  Function Inv_matrix(a) result(ainv) 
    implicit none 
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
  End Function Inv_matrix


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

end module numeth

