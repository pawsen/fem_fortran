MODULE arpack


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: arpack_plane

CONTAINS

  subroutine arpack_plane

    integer, parameter :: maxn=256, maxnev=10, maxncv=25,ldv=maxn

    !%--------------%
    !| Local Arrays |
    !%--------------%

    integer ::         iparam(11), ipntr(14)
    logical ::         select(maxncv)
    real(8) ::         ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),&
         v(ldv,maxncv), workd(3*maxn), workev(3*maxncv), &
         workl(3*maxncv*maxncv+6*maxncv), md(maxn), me(maxn-1)

    !%---------------%
    !| Local Scalars |
    !%---------------%

    character ::       bmat*1, which*2
    integer ::         ido, n, nev, ncv, lworkl, info, ierr, j,&
         nconv, maxitr, ishfts, mode
    real(8) ::         tol, sigmar, sigmai, h
    logical ::         first, rvec

    !%------------%
    !| Parameters |
    !%------------%

    real(8), parameter:: zero = 0.0D+0, one = 1.0D+0

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

    n     = 100 
    nev   = 4 
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
    which = 'LM'

    !%------------------------------------------------%
    !| M is the mass matrix formed by using piecewise |
    !| linear elements on [0,1].                      |

    !%------------------------------------------------%

    h = one / dble(n+1) !DBLE(A): Converts A to double precision real type. 
    do j = 1, n-1
       md(j) = 4.0D+0*h 
       me(j) = one*h 
    end do
    md(n) = 4.0D+0*h

    ! from dpttrf, md and me is:
    !  md      (input/output) DOUBLE PRECISION array, dimension (N)
    !          On entry, the n diagonal elements of the tridiagonal matrix
    !          A.  On exit, the n diagonal elements of the diagonal matrix
    !          D from the L*D*L**T factorization of A

    !  me      (input/output) DOUBLE PRECISION array, dimension (N-1)
    !          On entry, the (n-1) off-diagonal elements of the tridiagonal
    !          matrix A.
    !          On exit, the (n-1) off-diagonal elements of the unit
    !          bidiagonal factor L or U from the factorization of A.

    call dpttrf(n, md, me, ierr)
    if ( ierr .ne. 0 ) then
       print*, ' '
       print*, ' ERROR with _pttrf. ' 
       print*, ' ' 
       stop
    end if

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
    tol    = 0.0 
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
    mode   = 2

    iparam(1) = ishfts 
    iparam(3) = maxitr  
    iparam(7) = mode 

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

       if (ido .eq. -1 .or. ido .eq. 1) then

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
          call av (n, workd(ipntr(1)), workd(ipntr(2))) ! computes y=A.v
          call dpttrs(n, 1, md, me, workd(ipntr(2)), n,ierr) ! solve M.w=y for w

          if ( ierr .ne. 0 ) then
             print*, ' ' 
             print*, ' ERROR with _pttrs. ' 
             print*, ' '
             exit !quit loop
          end if

          !%-----------------------------------------%
          !| L O O P   B A C K to call DNAUPD again. |
          !%-----------------------------------------%
          cycle

       else if ( ido .eq. 2) then

          !%-------------------------------------%
          !|        Perform  y <--- M*x          |
          !| The matrix vector multiplication    |
          !| routine should take workd(ipntr(1)) |
          !| as input and return the result to   |
          !| workd(ipntr(2)).                    |
          !%-------------------------------------%

          call mv (n, workd(ipntr(1)), workd(ipntr(2)))

          !%-----------------------------------------%
          !| L O O P   B A C K to call DNAUPD again. |
          !%-----------------------------------------%
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

                if (d(j,2) .eq. zero)  then

                   !%--------------------%
                   !| Ritz value is real |
                   !%--------------------%

                   call av(n, v(1,j), ax)
                   call mv(n, v(1,j), mx)
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

                   call av(n, v(1,j), ax)
                   call mv(n, v(1,j), mx)
                   call daxpy(n, -d(j,1), mx, 1, ax, 1)
                   call mv(n, v(1,j+1), mx)
                   call daxpy(n, d(j,2), mx, 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)**2
                   call av(n, v(1,j+1), ax)
                   call mv(n, v(1,j+1), mx)
                   call daxpy(n, -d(j,1), mx, 1, ax, 1)
                   call mv(n, v(1,j), mx)
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

             call dmout(6, nconv, 3, d, maxncv, -6,'Ritz values (Real,Imag) and relative residuals')

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
          exit ! quit loop

       end if


       !%---------------------------%
       !| Done with program dndrv3. |
       !%---------------------------%
    end do

  end subroutine arpack_plane


  !==========================================================================
  !matrix vector multiplication subroutine

  subroutine av (n, v, w)
    integer           n, j
    real(8) :: v(n), w(n), dd, dl, du, s, h
    real(8), parameter :: rho = 1.0D+1, one = 1.0D+0,two = 2.0D+0

    ! Compute the matrix vector multiplication y<---A*x
    ! where A is stiffness matrix obtained from the finite element
    ! discretization of the 1-dimensional convection diffusion operator
    !                 d^2u/dx^2 + rho*(du/dx)
    ! on the interval [0,1] with zero Dirichlet boundary condition using
    ! linear elements.

    h = one / dble(n+1)
    s = rho / two
    dd = two / h
    dl = -one/h - s
    du = -one/h + s

    w(1) =  dd*v(1) + du*v(2)
    do j = 2,n-1
       w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1) 
    end do
    w(n) =  dl*v(n-1) + dd*v(n) 
    return
  end subroutine av


  !------------------------------------------------------------------------
  subroutine mv (n, v, w)
    integer           n, j
    real(8) ::  v(n), w(n), h
    real(8),  parameter :: one = 1.0D+0, four = 4.0D+0

    ! Compute the matrix vector multiplication y<---M*x
    ! where M is the mass matrix formed by using piecewise linear 
    ! elements on [0,1].

    !w = y, v = x

    w(1) =  four*v(1) + one*v(2)
    do j = 2,n-1
       w(j) = one*v(j-1) + four*v(j) + one*v(j+1) 
    end do
    w(n) =  one*v(n-1) + four*v(n) 

    h = one / dble(n+1)
    call dscal(n, h, w, 1)
    return

  end subroutine mv

end MODULE arpack
