module eigen


  implicit none

  private
  public :: 

contains



  subroutine mmul(Xp,Yp)
    ! This subroutine calculates Y = M*X

    use fedata
    use link1
    use plane42
    use mindlin41

    real(8), dimension(:), intent(in) :: Xp
    real(8), dimension(:), intent(out) :: Yp  

    integer :: e, i, nen
    integer, parameter :: mdim = 12
    integer, dimension(mdim) :: edof
    real(8), dimension(8) :: xe
    real(8), dimension(mdim, mdim) :: ke, me
    real(8) :: young, nu, area, dens, thk, bmat(3,12), N(3,12), xi, eta,B_M(5,12),jac(2,2),detjac

    Yp = 0.

    do e = 1, ne

       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          edof41(i) = element(e)%ix(i)
          ! Remember that the dimension of xe is only dependent on problem type(2d/3d), and not DOF pr node.
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)

          if ((element(1)%id == 2 )) then ! plane42 element
             edof(2*i-1) = 2 * element(e)%ix(i) - 1
             edof(2*i)   = 2 * element(e)%ix(i)
          else ! plate
             edof(3*i-2) = 3 * element(e)%ix(i) - 2 
             edof(3*i-1) = 3 * element(e)%ix(i) - 1
             edof(3*i)   = 3 * element(e)%ix(i)
          end if
       end do


       select case( element(e)%id )
       case( 1 )
          young = mprop(element(e)%mat)%young
          area  = mprop(element(e)%mat)%area
          call link1_ke(xe, young, area, ke)
       case( 2 )
          ! calculate plane42 element mass matrix
          young = mprop(element(e)%mat)%young
          nu  = mprop(element(e)%mat)%nu
          dens  = mprop(element(e)%mat)%dens
          thk = mprop(element(e)%mat)%thk
          call plane42_ke(xe, young, nu, dens, thk, ke, bmat, me)         
       case( 3 )
          ! calculate mindlin41 element mass matrix
          dens  = mprop(element(e)%mat)%dens
          thk = mprop(element(e)%mat)%thk
          xi = 0.
          eta = 0.
          call shape(xe,xi,eta,B_M,jac,detjac,N)
          call mindlin41_me(xe,dens,thk,N,me)
       end select

       ! Calculate Yp
       Yp(edof) = Yp(edof) + matmul(me,Xp(edof))

    end do

  end subroutine mmul


  subroutine eigen

    ! This subroutine calculates eigenvectors and eigenvalues

    use numeth
    use processor
    use fedata

    integer :: p_it, pmax
    real(8) :: eps, X_old(neqn),XP(neqn), YP(neqn), Y_old(neqn), X0(neqn), Y0(neqn), rp, ww(neqn), lambda, omega, error
    !allocate()

    ! Build load-vector
    call buildload

    ! Build stiffness and mass matrix
    call buildstiff

    ! Remove rigid body modes
    call enforce
    if (.not. banded) then
       call factor(k)
    else    
       call bfactor(k)
    endif

    X0 = 1.0
    call mmul(X0,Y0)
    Yp = Y0
    Xp = X0

    pmax = 100        
    do p_it = 1,pmax
       Y_old = Yp
       X_old = Xp
       if (.not. banded) then
          call solve(k, Yp)
       else    
          call bsolve(k,Yp)
       endif
       Xp = Yp
       Yp = Y_old

       call mmul(Xp,Yp)

       rp = SQRT(DOT_PRODUCT(Xp,Yp))
       Yp = Yp/rp

       !   eps = 10**(-12)d0
       !           123456789012
       eps = 0.000000000001d0
       if (p_it .ne. 1) then
          error = abs(MAXVAL(Xp-X_old)/MAXVAL(Xp))
          if(error .lt. eps) then
             print*, "Break, p_it =", p_it    
             exit
          end if
       end if
    end do

    ! Eigenvector D and eigenvalue lambda: "(32(f5.2,tr1))"
    ww = Xp/rp
    lambda = DOT_PRODUCT(Xp,Y_old)/rp**2
	omega = SQRT(lambda)

    print*, "omega =",omega

    ! Plot vibration modes
    call plot('eigenmode', 'xwin', 'color', ' ',(/lambda/),ww(1:neqn),(/1000.d0,1.d0/))
    call pgend


  end subroutine eigen

  subroutine multi_eigen

    ! This subroutine calculates eigenvectors and eigenvalues

    use numeth
    use processor
    use fedata

    integer :: p_it, pmax, i, j, idof, b
    integer, parameter :: neig = 1
    real(8) :: eps, X_old(neqn),XP(neqn), YP(neqn), Y_old(neqn), X0(neqn), Y0(neqn), rp, ww(neqn,neig), lambda(neig), omega(neig)
    real(8) :: error, Z(neqn,neig), cj
    !allocate()

    ! Build load-vector
    call buildload

    ! Build stiffness and mass matrix
    call buildstiff

    ! Remove rigid body modes
    call enforce

    if (.not. banded) then
       call factor(k)
    else    
       call bfactor(k)
    endif

    ! Calculation of multiple eigenvalues and eigenvectors

    do i = 1,neig
       X0 = 1.0
       call mmul(X0,Y0)
       Yp = Y0
       Xp = X0
       Z = 0.0
       do j=1,i-1
          call mmul(ww(1:neqn,j),Z(1:neqn,j))
       end do


       pmax = 1000        
       do p_it = 1,pmax
          Y_old = Yp
          X_old = Xp

          do b = 1, nb
             idof = 3*(bound(b,1)-1) + bound(b,2)
             Yp(1:neqn) = Yp(1:neqn) - k(1:neqn, idof) * bound(b, 3)
             Yp(idof) = bound(b, 3)
             !$             k(1:neqn, idof) = 0.
             !$             k(idof, 1:neqn) = 0.
             !$             k(idof, idof) = 1.
          end do

          if (.not. banded) then
             call solve(k, Yp)
          else    
             call bsolve(k,Yp)
          endif
          Xp = Yp
          Yp = Y_old
          do j = 1,i-1
             cj = dot_product(Xp,Z(1:neqn,j))
             Xp = Xp - cj*ww(1:neqn,j)
          end do

          call mmul(Xp,Yp)

          rp = SQRT(DOT_PRODUCT(Xp,Yp))
          Yp = Yp/rp

          !   eps = 10**(-12)d0
          !           123456789012
          eps = 0.000000000001d0
          if (p_it .ne. 1) then
             error = abs(MAXVAL(Xp-X_old)/MAXVAL(Xp))
             if(error .lt. eps) then
                print*, "Break, p_it =", p_it    
                exit
             end if
          end if
       end do

       ! Eigenvector D and eigenvalue lambda:
       ww(:,i) = Xp/rp
       lambda(i) = DOT_PRODUCT(Xp,Y_old)/rp**2
       omega(i) = SQRT(lambda(i))

       print*, "omega",i,"=",omega(i)

    end do

    ! Plot vibration modes
    call plot('eigenmode', 'matlab', 'color', ' ',(/lambda(neig)/),ww(1:neqn,neig),(/200.d0,1.d0/),ww(1:neqn,1:neig),neig)
    call pgend

  end subroutine multi_eigen

end module eigen
