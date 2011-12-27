MODULE arpack


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: arpack_init, arpack_plane

CONTAINS

  subroutine arpack_init

    use fedata
    use fea
    use piezo
    use solve_handle_real
    use exodus

    real(8), allocatable :: tmp_mat(:,:)
    integer :: nzz,i, n_eigen, j
    real(8) :: sigma
    logical :: shift
    real(8), allocatable :: eigenval(:,:), eigenvec(:,:)

    sigma = 1000d0
    eigenvalue%sigma = sigma
    shift = eigenvalue%shift    
    
    n_eigen = 10

    allocate(eigenval(n_eigen,3), eigenvec(neqn,n_eigen))

    call exodus_write(eigenval,eigenvec)

    if (.not.shift) then
       call build_mvec !build mass-vector
       
       select case( element(1)%id )
       case(2)
          call buildstiff_fea(0) ! build sparse stiffness-matrix
       case default
          call buildstiff_eigenvalue_piezo
       end select

!!$       call buildmass_fea(0)
!!$       call mumps_init_real
!!$       call mumps_solve_real(4)
    else
       call build_mvec
       call buildstiff_fea(0)
       call mumps_init_real
       call mumps_solve_real(4)
    end if


    call arpack_plane(n_eigen,neqn,shift,sigma,iK,jK,sK,eigenval,eigenvec)
    call exodus_write(eigenval,eigenvec)

  end subroutine arpack_init

  subroutine arpack_plane(n_eigenvalue,neqn,shift,sigma,iK,jK,sK,eigenval,eigenvec)
    
    use fea
    use solve_handle_real
    use numeth
    use sort_array
   
    integer, intent(in) :: n_eigenvalue, neqn, iK(:), jK(:)
    real(8), intent(out) :: eigenval(:,:), eigenvec(:,:)
    real(8), intent(in), optional :: sK(:)
    real(8), intent(IN) :: sigma
    logical, intent(IN) :: shift

    integer, parameter :: maxn=1000000, maxnev=20, maxncv=25,ldv=maxn
    integer :: i, index_list(n_eigenvalue), print_info
    !%--------------%
    !| Local Arrays |
    !%--------------%

    integer ::         iparam(11), ipntr(14)
    logical ::         select(maxncv)
    real(8) ::         ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),&
         v(ldv,maxncv), workd(3*maxn), workev(3*maxncv), &
         workl(3*maxncv*maxncv+6*maxncv)
    !v: indeholder egenvektorer
    !d: indeholder egenværdier
    
    !%---------------%
    !| Local Scalars |
    !%---------------%

    character ::       bmat*1, which*2
    integer ::         ido, n, nev, ncv, lworkl, info, ierr, j,&
         nconv, maxitr, ishfts, mode
    real(8) ::         tol, sigmar, sigmai
    logical ::         first, rvec


    !%-----------------------------%
    !| BLAS & LAPACK routines used |
    !%-----------------------------%
    real(8) ::          dnrm2, dlapy2
    !external ::        daxpy, dnrm2, dpttrf, dpttrs, dlapy2

    !%----------------------------------------------------%
    !| The number N is the dimension of the matrix.  A    |
    !| generalized eigenvalue problem is solved (BMAT =   |
    !| 'G').  NEV is the number of eigenvalues to be      |
    !| approximated.  The user can modify NEV, NCV, WHICH |
    !| to solve problems of different sizes, and to get   |
    !| different parts of the spectrum.  However, The     |
    !| following conditions must be satisfied:            |
    !|                    N <= MAXN,                      |
    !|                  NEV <= MAXNEV,                    |
    !|              NEV + 2 <= NCV <= MAXNCV              |
    !%----------------------------------------------------%

    n     = neqn
    nev   = n_eigenvalue
    ncv   = 20 
    if ( n .gt. maxn ) then
       print *, ' ERROR with _NDRV3: N is greater than MAXN '
       stop
    else if ( nev .gt. maxnev ) then
       print *, ' ERROR with _NDRV3: NEV is greater than MAXNEV '
       stop
    else if ( ncv .gt. maxncv ) then
       print *, ' ERROR with _NDRV3: NCV is greater than MAXNCV '
       stop
    end if
    bmat  = 'G'
    which = 'SM' !smallest magnitude
    sigmar = sigma ! shift freqvency
    sigmai = 0d0

    !%-----------------------------------------------------%
    !| The work array WORKL is used in DNAUPD as           |
    !| workspace.  Its dimension LWORKL is set as          |
    !| illustrated below.  The parameter TOL determines    |
    !| the stopping criterion. If TOL<=0, machine          |
    !| precision is used.  The variable IDO is used for    |
    !| reverse communication, and is initially set to 0.   |
    !| Setting INFO=0 indicates that a random vector is    |
    !| generated in DNAUPD to start the Arnoldi iteration. |
    !%-----------------------------------------------------%

    lworkl = 3*ncv**2+6*ncv 
    tol    = 0d0!1E-!0.0 
    ido    = 0
    info   = 0

    !%---------------------------------------------------%
    !| This program uses exact shifts with respect to    |
    !| the current Hessenberg matrix (IPARAM(1) = 1).    |
    !| IPARAM(3) specifies the maximum number of Arnoldi |
    !| iterations allowed.  Mode 2 of DNAUPD is used     |
    !| (IPARAM(7) = 2).  All these options can be        |
    !| changed by the user. For details, see the         |
    !| documentation in DNAUPD.                          |
    !%---------------------------------------------------%

    ishfts = 1
    maxitr = 300
    if (shift) then
       mode   = 3
    else
       mode   = 2
    end if

    iparam(1) = ishfts 
    iparam(3) = maxitr  
    iparam(7) = mode 


print*,'mode ', mode
    !%-------------------------------------------%
    !| M A I N   L O O P (Reverse communication) |
    !%-------------------------------------------%

    do

       !%---------------------------------------------%
       !| Repeatedly call the routine DNAUPD and take | 
       !| actions indicated by parameter IDO until    |
       !| either convergence is indicated or maxitr   |
       !| has been exceeded.                          |
       !%---------------------------------------------%    

       call dnaupd ( ido, bmat, n, which, nev, tol, resid, & 
            ncv, v, ldv, iparam, ipntr, workd, &
            workl, lworkl, info )

       if ((ido .eq. -1 .or. ido .eq. 1) .and. mode == 2 ) then

          !%----------------------------------------%
          !| Perform  y <--- OP*x = inv[M]*A*x      |
          !| The user should supply his/her own     |
          !| matrix vector routine and a linear     |
          !| system solver.  The matrix-vector      |
          !| subroutine should take workd(ipntr(1)) |
          !| as input, and the final result should  |
          !| be returned to workd(ipntr(2)).        |
          !%----------------------------------------%

          ! these two are user supplied
          ! computes y=A.v
          call sparse_multiply(iK,jK,sK,workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1),.true.) ! multiply K.v

          ! solve M.w=y for w
          ! M is lumped an saved as a vector. eg the system is solved by division
           call mvec_mul((/0d0/), workd(ipntr(2):ipntr(2)+n-1),2)

          !call mumps_solve_arpack_real(workd(ipntr(2):ipntr(2)+n-1),3) 

          !%-----------------------------------------%
          !| L O O P   B A C K to call DNAUPD again. |
          !%-----------------------------------------%
          cycle
       else if ( ido .eq. -1 .and. mode == 3 ) then
          
          !%-------------------------------------------%
          !| Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
          !| to force starting vector into the range   |
          !| of OP.   The user should supply his/her   |
          !| own matrix vector multiplication routine  |
          !| and a linear system solver.  The matrix   |
          !| vector multiplication routine should take |
          !| workd(ipntr(1)) as the input. The final   |
          !| result should be returned to              |
          !| workd(ipntr(2)).                          |
          !%-------------------------------------------%

          ! Mv = y
          call mvec_mul(workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1),1)

          ! solve (A-sigma*M)w = y
          call mumps_solve_arpack_real(workd(ipntr(2):ipntr(2)+n-1),3)
          cycle
       else if ( ido .eq. 1 .and. mode == 3 ) then

          !%-----------------------------------------%
          !| Perform y <-- OP*x = inv[A-sigma*M]*M*x |
          !| M*x has been saved in workd(ipntr(3)).  |
          !| The user only need the linear system    |
          !| solver here that takes workd(ipntr(3))  |
          !| as input, and returns the result to     |
          !| workd(ipntr(2)).                        |
          !%-----------------------------------------%
          
          ! solve (A-sigma*M)w = v
          call dcopy( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1) ! copy vector
          call mumps_solve_arpack_real(workd(ipntr(2):ipntr(2)+n-1),3)       
          cycle
       else if ( ido .eq. 2) then

          !%-------------------------------------%
          !|        Perform  y <--- M*x          |
          !| The matrix vector multiplication    |
          !| routine should take workd(ipntr(1)) |
          !| as input and return the result to   |
          !| workd(ipntr(2)).                    |
          !%-------------------------------------%

          call mvec_mul(workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1),1)
          cycle

       end if

       !%-----------------------------------------%
       !| Either we have convergence, or there is |
       !| an error.                               |
       !%-----------------------------------------% 

       if ( info .lt. 0 ) then

          !%--------------------------%
          !| Error message. Check the |
          !| documentation in DNAUPD. |
          !%--------------------------%

          print *, ' '
          print *, ' Error with _naupd, info = ', info
          print *, ' Check the documentation of _naupd.'
          print *, ' ' 
       else 

          !%-------------------------------------------%
          !| No fatal errors occurred.                 |
          !| Post-Process using DNEUPD.                |
          !|                                           |
          !| Computed eigenvalues may be extracted.    |  
          !|                                           |
          !| Eigenvectors may also be computed now if  |
          !| desired.  (indicated by rvec = .true.)    | 
          !%-------------------------------------------%

          rvec = .true.
          call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, &
               sigmar, sigmai, workev, bmat, n, which, nev, tol, &
               resid, ncv, v, ldv, iparam, ipntr, workd, &
               workl, lworkl, ierr )

          !%-----------------------------------------------%
          !| The real part of the eigenvalue is returned   |
          !| in the first column of the two dimensional    |
          !| array D, and the IMAGINARY part is returned   |
          !| in the second column of D.  The corresponding |
          !| eigenvectors are returned in the first NEV    |
          !| columns of the two dimensional array V if     |
          !| requested.  Otherwise, an orthogonal basis    |
          !| for the invariant subspace corresponding to   |
          !| the eigenvalues in D is returned in V.        |
          !%-----------------------------------------------%

          if ( ierr .ne. 0 ) then

             !%------------------------------------%
             !| Error condition:                   |
             !| Check the documentation of DNEUPD. |
             !%------------------------------------%

             print *, ' ' 
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd'
             print *, ' '

          else 

             first = .true.
             nconv = iparam(5)
             do j=1, iparam(5)

                !%---------------------------%
                !| Compute the residual norm |
                !|                           |
                !|  ||  A*x - lambda*M*x ||  |
                !|                           |
                !| for the NCONV accurately  |
                !| computed eigenvalues and  |
                !| eigenvectors.  (iparam(5) |
                !| indicates how many are    |
                !| accurate to the requested |
                !| tolerance)                |
                !%---------------------------%

                if (d(j,2) .eq. 0)  then

                   !%--------------------%
                   !| Ritz value is real |
                   !%--------------------%

                   call mmul_fea(v(1:n,j),ax(1:n),1) !av
                   call mmul_fea(v(1:n,j),mx(1:n),2) !mv
                   !call av(n, v(1,j), ax)
                   !call mv(n, v(1,j), mx)
                   call daxpy(n, -d(j,1), mx, 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)
                   d(j,3) = d(j,3) / abs(d(j,1))

                else if (first) then

                   !%------------------------%
                   !| Ritz value is complex  |
                   !| Residual of one Ritz   |
                   !| value of the conjugate |
                   !| pair is computed.      |
                   !%------------------------%

                   !call av(n, v(1,j), ax)
                   !call mv(n, v(1,j), mx)
                   call mmul_fea(v(1:n,j),ax(1:n),1) !av
                   call mmul_fea(v(1:n,j),mx(1:n),2) !mv
                   call daxpy(n, -d(j,1), mx, 1, ax, 1)
                   call mmul_fea(v(1:n,j+1),mx(1:n),2) !mv      
                   !call mv(n, v(1,j+1), mx)
                   call daxpy(n, d(j,2), mx, 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)**2
                   call mmul_fea(v(1:n,j+1),ax(1:n),1) !av
                   call mmul_fea(v(1:n,j+1),mx(1:n),2) !mv
                   !call av(n, v(1,j+1), ax)
                   !call mv(n, v(1,j+1), mx)
                   call daxpy(n, -d(j,1), mx, 1, ax, 1)
                   call mmul_fea(v(1:n,j),mx(1:n),2) !mv  
                   !call mv(n, v(1,j), mx)
                   call daxpy(n, -d(j,2), mx, 1, ax, 1)
                   d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                   d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if

             end do

             !%-----------------------------%
             !| Display computed residuals. |
             !%-----------------------------%

             print_info = 0
             if (print_info == 1) then
                call dmout(6, nconv, 3, d, maxncv, -6,'Ritz values (Real,Imag) and relative residuals')
             end if

             do i=1,nev
                print*,'eigenfrekvency ', SQRT(d(i,1)) /(2.0*3.1415927 )
             end do
             print*

          end if

          !%------------------------------------------%
          !| Print additional convergence information |
          !%------------------------------------------%


          if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
          else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',&
                  ' Arnoldi update, try increasing NCV.'
             print *, ' '
          end if

          if (print_info == 1) then
             print *, ' '
             print *, ' _NDRV3 '
             print *, ' ====== '
             print *, ' '
             print *, ' Size of the matrix is ', n
             print *, ' The number of Ritz values requested is ', nev
             print *, ' The number of Arnoldi vectors generated',&
                  ' (NCV) is ', ncv
             print *, ' What portion of the spectrum: ', which
             print *, ' The number of converged Ritz values is ', &
                  nconv 
             print *, ' The number of Implicit Arnoldi update',&
                  ' iterations taken is ', iparam(3)
             print *, ' The number of OP*x is ', iparam(9)
             print *, ' The convergence criterion is ', tol
             print *, ' '
          end if


          exit ! quit loop

       end if


       !%---------------------------%
       !| Done with program dndrv3. |
       !%---------------------------%
    end do

    if (nconv >2) then
       eigenval(1:nconv,:) = d(1:nconv,1:3)
       eigenvec(:,1:nconv) = v(1:neqn, 1:nconv)

    else
       stop
    end if

!!$    eigenval = d(1:n_eigenvalue,1:3)
!!$    call Shell_Sort_real(eigenval(:,1),index_list)
!!$    eigenval(:,2:3) = d(index_list,2:3)
!!$    eigenvec = v(1:neqn, index_list)

!!$    do i=1,nev
!!$       print*,'eigenfrekvency ', SQRT(eigenval(i,1)) /(2.0*3.1415927 )
!!$    end do
!!$    print*
!!$
!!$    do i=1,nev
!!$       print*,'eigenfrekvency ', index_list(i)
!!$    end do
!!$    print*


  end subroutine arpack_plane



end MODULE arpack
