MODULE plane42nonlin

 ! This module contains subroutines specific to non-linear calculations

 IMPLICIT NONE
 
 PRIVATE :: gauss_points
 PUBLIC :: plane42nonlin_ke, plane42nonlin_ss, plane42nonlin_shape,plane42nonlin_residual

CONTAINS

 SUBROUTINE plane42nonlin_ke(xe,de, young, nu, thk,ng, ke,eresidual)! Tangent stiffness, non-lin deformation. WITHOUT element stress stifness

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng
  REAL(8), INTENT(IN) :: young, nu, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe,de
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke
  REAL(8), DIMENSION(:), INTENT(OUT) :: eresidual
	
  INTEGER :: i, j
  REAL(8) :: fact, Cmat(3,3), detjac, jac(2,2), Bmat(3,8), Nmat(2,8), helpproduct1(8,3),helpproduct2(8,4)
  REAL(8) :: estress(3), estrain(3), B0mat(3,8), BLmat(3,8), Gmat(4,8), Mmat(4,4)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)
 
  ! Build constitutive matrix (plane stress)
  
  Cmat = 0.0
  fact = young/(1.-nu**2)
  Cmat(1, 1) = fact
  Cmat(1, 2) = fact*nu
  Cmat(2, 1) = fact*nu
  Cmat(2, 2) = fact
  Cmat(3, 3) = fact*(1.-nu)/2.
  
  ke=0.0
  eresidual =0.0
	! Måske skal Mmat flyttes ind i løkken. Jeg tror det ikke, da spændingen altid evalueres i midten!
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42nonlin_shape(xe,de, xi(i), eta(j), Nmat,Bmat,B0mat,BLmat,Gmat, jac, detjac)
          	call plane42nonlin_ss(xe, de, young, nu,B0mat,BLmat, estress, estrain)
  			Mmat = 0.0
  			Mmat(1,1) = estress(1)
  			Mmat(2,2) = estress(1)
  			Mmat(3,1) = estress(3)
  			Mmat(4,2) = estress(3)
  
  			Mmat(1,3) = estress(3)
  			Mmat(2,4) = estress(3)
  			Mmat(3,3) = estress(2)
  			Mmat(4,4) = estress(2)

        	helpproduct1 = MATMUL(TRANSPOSE(Bmat),Cmat)
            helpproduct2 = MATMUL(TRANSPOSE(Gmat),Mmat)
            ke =ke+ W(i)*W(j)*thk*detjac*(MATMUL(helpproduct1,Bmat)+MATMUL(helpproduct2,Gmat))

            eresidual = eresidual+ W(i)*W(j)*thk*MATMUL(transpose(Bmat),estress)*detjac
         end do
	end do
	
 END SUBROUTINE plane42nonlin_ke


 SUBROUTINE plane42nonlin_ss(xe, de, young, nu,B0mat, BLmat, estress, estrain)

  ! This subrotuine constructs the element stress and strain
 
  REAL(8), INTENT(IN) :: young, nu
  REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
  REAL(8), DIMENSION(:,:), INTENT(IN)  :: B0mat, BLmat
  REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain

  REAL(8) :: fact, helpproduct(3,8), Cmat(3, 3)

 ! Build strain-displacement matrix 
  estrain = 0.0
  ! Compute non-lin element strain
  helpproduct = (B0mat+ 0.5*BLmat)
  estrain = MATMUL(helpproduct,de)

  ! Build constitutive matrix (plane stress)
  Cmat = 0.0
  fact = young/(1.-nu**2)
  Cmat(1, 1) = fact
  Cmat(1, 2) = fact*nu
  Cmat(2, 1) = fact*nu
  Cmat(2, 2) = fact
  Cmat(3, 3) = fact*(1.-nu)/2.

  ! Compute element stress
  estress = MATMUL(Cmat, estrain)

  ! Compute principal stress and direction

 END SUBROUTINE plane42nonlin_ss


SUBROUTINE plane42nonlin_shape(xe,de, xi, eta, Nmat,Bmat,B0mat,BLmat,Gmat, jac, detjac)
 ! This subrotuine implement isoparametric elements
 
    REAL(8), INTENT(IN) :: xi, eta
  REAL(8), DIMENSION(:), INTENT(IN)  :: xe,de
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: bmat, jac
  REAL(8), INTENT(OUT) :: Nmat(2,8),B0mat(3, 8),BLmat(3, 8),Gmat(4,8)
  REAL(8), INTENT(OUT) :: detjac

  REAL(8) :: Ntilde(4,8), GAMMA(2,2), GAMMAtilde(4,4), N1,N2,N3,N4 , L(3,4), produkthelp(3,4)
  REAL(8) :: JACtilde(4,4),Amat(3,4),theta(4)
    Nmat = 0.
    N1 = 1.0/4.0*(1-xi)*(1-eta)!(8.1) i noterne
    N2 = 1.0/4.0*(1+xi)*(1-eta)
    N3 = 1.0/4.0*(1+xi)*(1+eta)
    N4 = 1.0/4.0*(1-xi)*(1+eta)

    Nmat(1,1) = N1
    Nmat(2,2) = N1
    Nmat(1,3) = N2
    Nmat(2,4) = N2
    Nmat(1,5) = N3
    Nmat(2,6) = N3
    Nmat(1,7) = N4
    Nmat(2,8) = N4
    
    L(1,1) = 1. !(8.6) i noterne
    L(2,1) = 0.
    L(3,1) = 0.
    
    L(1,2) = 0.
    L(2,2) = 0.
    L(3,2) = 1.
    
    L(1,3) = 0.
    L(2,3) = 0.
    L(3,3) = 1.
    
    L(1,4) = 0.
    L(2,4) = 1.
    L(3,4) = 0.

    jac(1,1) = 1.0/4.0*(-(1-eta)*xe(1)+(1-eta)*xe(3)+(1+eta)*xe(5)-(1+eta)*xe(7)) ! Jacobi-matrix s.207 i cook
    jac(1,2) = 1.0/4.0*(-(1-eta)*xe(2)+(1-eta)*xe(4)+(1+eta)*xe(6)-(1+eta)*xe(8))
    jac(2,1) = 1.0/4.0*(-(1-xi)*xe(1)-(1+xi)*xe(3)+(1+xi)*xe(5)+(1-xi)*xe(7))
    jac(2,2) = 1.0/4.0*(-(1-xi)*xe(2)-(1+xi)*xe(4)+(1+xi)*xe(6)+(1-xi)*xe(8))

    detjac = jac(1,1)*jac(2,2)-jac(2,1)*jac(1,2) ! determinant af Jacobi

    GAMMA(1,1) = 1.0/detjac*jac(2,2) ! invers af Jacobi-matrix
    GAMMA(2,1) = -1.0/detjac*jac(2,1)
    GAMMA(1,2) = -1.0/detjac*jac(1,2)
    GAMMA(2,2) = 1.0/detjac*jac(1,1)

    GAMMAtilde(1,1) = GAMMA(1,1)
    GAMMAtilde(2,1) = GAMMA(2,1)
    GAMMAtilde(1,2) = GAMMA(1,2)
    GAMMAtilde(2,2) = GAMMA(2,2)

    GAMMAtilde(3,1) = 0.0
    GAMMAtilde(4,1) = 0.0
    GAMMAtilde(3,2) = 0.
    GAMMAtilde(4,2) = 0.

    GAMMAtilde(1,3) = 0.
    GAMMAtilde(2,3) = 0.
    GAMMAtilde(1,4) = 0.
    GAMMAtilde(2,4) = 0.  

    GAMMAtilde(3,3) = GAMMA(1,1)
    GAMMAtilde(4,3) = GAMMA(2,1)
    GAMMAtilde(3,4) = GAMMA(1,2)
    GAMMAtilde(4,4) = GAMMA(2,2)
    
    Ntilde = 0.0
    Ntilde(1,1) = -1.0/4.0*(1-eta) !N1,xi
    Ntilde(2,1) = -1.0/4.0*(1-xi) !N1,eta

    Ntilde(3,2) = -1.0/4.0*(1-eta) !N1,xi
    Ntilde(4,2) = -1.0/4.0*(1-xi) !N1,eta 
    
    Ntilde(1,3) = 1.0/4.0*(1-eta) !N2
    Ntilde(2,3) = -1.0/4.0*(1+xi)

    Ntilde(3,4) = 1.0/4.0*(1-eta)
    Ntilde(4,4) = -1.0/4.0*(1+xi)

    Ntilde(1,5) = 1.0/4.0*(1+eta) !N3
    Ntilde(2,5) = 1.0/4.0*(1+xi)

    Ntilde(3,6) = 1.0/4.0*(1+eta)
    Ntilde(4,6) = 1.0/4.0*(1+xi)
    
    Ntilde(1,7) = -1.0/4.0*(1+eta) !N4
    Ntilde(2,7) = 1.0/4.0*(1-xi)

    Ntilde(3,8) = -1.0/4.0*(1+eta)
    Ntilde(4,8) = 1.0/4.0*(1-xi)

    JACtilde = 0.0 ! For calculating G in non-lin deformation
    JACtilde(1,1) = GAMMA(1,1)
    JACtilde(3,1) = GAMMA(2,1)
    JACtilde(1,2) = GAMMA(1,2)
    JACtilde(3,2) = GAMMA(2,2)
    
    JACtilde(2,3) = GAMMA(1,1)
    JACtilde(4,3) = GAMMA(2,1)
    JACtilde(2,4) = GAMMA(1,2)
    JACtilde(4,4) = GAMMA(2,2)

    Gmat = MATMUL(JACtilde,Ntilde)

    theta=MATMUL(Gmat,de)
    Amat = 0.0

    Amat(1,1)=theta(1)
    Amat(1,2)=theta(2)

    Amat(2,3)=theta(3)
    Amat(2,4)=theta(4)
    
    Amat(3,1)=theta(3)
    Amat(3,2)=theta(4)
    Amat(3,3)=theta(1)
    Amat(3,4)=theta(2)
    
    produkthelp = MATMUL(L, GAMMAtilde) !hjælpeprodukt ved B=[L]*[GAMMAtilde]*[Ntilde]
    b0mat = MATMUL(produkthelp, Ntilde)! Bmat for linear deformation
    bLmat = MATMUL(Amat,Gmat)! Bmat non-lin deformation. (8) i extra note
    bmat = b0mat + bLmat !Bmat used for calculating tangent stiffness
END SUBROUTINE plane42nonlin_shape

      
 SUBROUTINE plane42nonlin_residual(xe,de, young,nu,thk,ng, eresidual)! Tangent stiffness, non-lin deformation. WITHOUT element stress stifness

  IMPLICIT NONE
  ! This subrotuine calculate residualet
 
  INTEGER, INTENT(IN) :: ng
  REAL(8), INTENT(IN) :: thk,young,nu
  REAL(8), DIMENSION(:), INTENT(IN) :: xe,de
  REAL(8), DIMENSION(:), INTENT(OUT) :: eresidual

  integer :: i,j
  REAL(8) :: Bmat(3, 8), Nmat(2,8), jac(2,2),estress(3),estrain(3)
  REAL(8) :: B0mat(3, 8),BLmat(3, 8), Gmat(4,8),detjac
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

 eresidual =0.0
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42nonlin_shape(xe,de, xi(i), eta(j), Nmat,Bmat,B0mat,BLmat,Gmat, jac, detjac)
   		 	call plane42nonlin_ss(xe, de, young, nu,B0mat, BLmat, estress, estrain)
        	eresidual = eresidual+ W(i)*W(j)*thk*MATMUL(transpose(Bmat),estress)*detjac
        end do
	end do
	
 END SUBROUTINE plane42nonlin_residual

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

END MODULE plane42nonlin