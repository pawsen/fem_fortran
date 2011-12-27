MODULE plane42transient

 ! This module contains subroutines specific to the PLANE42 element.

 use numeth ! guss points ligger i numeth
 use plane42 ! til shape
 IMPLICIT NONE

 PRIVATE
 PUBLIC :: plane42transient_ke,plane42transient_me,plane42transient_ce,plane42transient_KRe,&
 		plane42transient_CZe, PLANE42transient_DKE, PLANE42transient_DME, PLANE42TRANSIENT_CZE_PUNKT

CONTAINS

 SUBROUTINE plane42transient_ke(xe, young1,young2, nu1, nu2, thk,ng, ke1,ke2)

  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), INTENT(IN) :: young1, young2, nu1, nu2, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke1,ke2
	
  INTEGER :: i, j
  ! TopOpt
  REAL(8) :: fact, Cmat1(3,3),Cmat2(3,3), helpproduct1(8,3), helpproduct2(8,3)
  
  real(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta


  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  ! Build constitutive matrix (plane stress)
  
  Cmat1 = 0.0
  fact = young1/(1.0d0-nu1**2)
  Cmat1(1, 1) = fact
  Cmat1(1, 2) = fact*nu1
  Cmat1(2, 1) = fact*nu1
  Cmat1(2, 2) = fact
  Cmat1(3, 3) = fact*(1.0d0-nu1)/2.0d0

  Cmat2 = 0.0d0
  fact = young2/(1.0d0-nu2**2)
  Cmat2(1, 1) = fact
  Cmat2(1, 2) = fact*nu2
  Cmat2(2, 1) = fact*nu2
  Cmat2(2, 2) = fact
  Cmat2(3, 3) = fact*(1.0d0-nu2)/2.0d0

  ke1=0.0d0
  ke2=0.0d0

	
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
        	helpproduct1 = MATMUL(TRANSPOSE(bmat),Cmat1)
            helpproduct2 = MATMUL(TRANSPOSE(bmat),Cmat2)
        	ke1 =ke1+ W(i)*W(j)*thk*MATMUL(helpproduct1,bmat)*detjac
            ke2 =ke2+ W(i)*W(j)*thk*MATMUL(helpproduct2,bmat)*detjac
         end do
	end do
 
	
 END SUBROUTINE plane42transient_ke


 
 SUBROUTINE plane42transient_me(xe,thk,ng, me)

  ! This subroutine constructs the damping matrix for
  ! a rectangular 4-noded quad element.


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), INTENT(IN) ::  thk
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: me
	
  INTEGER :: i, j
  REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta


  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

    me = 0.0d0
	
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
            me =me+ W(i)*W(j)*thk*MATMUL(TRANSPOSE(Nmat),Nmat)*detjac 
		end do
	end do

	
 END SUBROUTINE plane42transient_me

  SUBROUTINE plane42transient_dke(xe, young1,young2, nu1, nu2, thk,ng, dke)

  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), INTENT(IN) :: young1, young2, nu1, nu2, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: dke
	
  INTEGER :: i, j
  ! TopOpt
  REAL(8) :: fact, Cmat1(3,3),Cmat2(3,3), helpproduct1(8,3), helpproduct2(8,3)
  REAL(8) :: ke1(8,8), ke2(8,8)
  
  real(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta


  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  ! Build constitutive matrix (plane stress)
  
  Cmat1 = 0.0
  fact = young1/(1.0d0-nu1**2)
  Cmat1(1, 1) = fact
  Cmat1(1, 2) = fact*nu1
  Cmat1(2, 1) = fact*nu1
  Cmat1(2, 2) = fact
  Cmat1(3, 3) = fact*(1.0-nu1)/2.0

  Cmat2 = 0.0
  fact = young2/(1.0d0-nu2**2)
  Cmat2(1, 1) = fact
  Cmat2(1, 2) = fact*nu2
  Cmat2(2, 1) = fact*nu2
  Cmat2(2, 2) = fact
  Cmat2(3, 3) = fact*(1.0-nu2)/2.0

  ke1=0.0d0
  ke2=0.0d0

	
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
        	helpproduct1 = MATMUL(TRANSPOSE(bmat),Cmat1)
            helpproduct2 = MATMUL(TRANSPOSE(bmat),Cmat2)
        	ke1 =ke1+ W(i)*W(j)*thk*MATMUL(helpproduct1,bmat)*detjac
            ke2 =ke2+ W(i)*W(j)*thk*MATMUL(helpproduct2,bmat)*detjac
         end do
	end do
    dke = ke1 - ke2
	
 END SUBROUTINE plane42transient_dke


 
 SUBROUTINE plane42transient_dme(xe,thk,dens1, dens2, ng, dme)

  ! This subroutine constructs the damping matrix for
  ! a rectangular 4-noded quad element.


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), INTENT(IN) :: dens1, dens2, thk
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: dme
	
  INTEGER :: i, j
  REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta


  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

    dme = 0.0
	
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
            dme =dme+ W(i)*W(j)*thk*MATMUL(TRANSPOSE(Nmat),Nmat)*detjac
		end do
	end do
    dme = ( dens1 - dens2 ) * dme
	
 END SUBROUTINE plane42transient_dme


 SUBROUTINE plane42transient_ce(xe,c_damp, ng, ce)

  ! This subroutine constructs the damping matrix for
  ! a rectangular 4-noded quad element.


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng ! number of gauss-points
  REAL(8), INTENT(IN) :: c_damp
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ce
	
  INTEGER :: i, j
  REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta


  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  ce=0.0

	
	do i= 1,ng !number of gauss points
		do j = 1,ng
        	call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
            ce =ce+ W(i)*W(j)*c_damp*MATMUL(TRANSPOSE(Nmat),Nmat)*detjac
		end do
	end do
	
 END SUBROUTINE plane42transient_ce


SUBROUTINE plane42transient_CZe(xe,eface, young1,young2,dens1,dens2, nu1, nu2, thk, length,ng,rho, CZe)


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng,eface ! number of gauss-points
  REAL(8), INTENT(IN) :: young1,young2,dens1,dens2, nu1, nu2, thk, length, rho
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: CZe
	
  INTEGER :: i, j, numerisk
  
  REAL(8) :: jac(2,2), detjac, bmat(3,8), Nmat(2,8), Zmat(2,2)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

   real(8) :: dens, young, nu, lambda, mu, cL, cT, rr1, rr2
  

  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  dens = rho*dens1+(1.0-rho)*dens2
  young = rho*young1+(1.0-rho)*young2
  nu = rho*nu1+(1.0-rho)*nu2
  !Lame parameters  
  lambda = young*nu/ ((1.0+nu)*(1.0-2.0*nu)) 
  mu = young/(2.0*(1.0+nu))! shear modulus
  cL = dsqrt( (lambda+2.0*mu) / dens)
  cT = dsqrt( mu / dens)

  
! F�lgende er kun for b�lgeretning vinkelret p� fladen
!
	numerisk = 0
  
! Numerisk integration - Virker, men langsommere.
 select case (numerisk)
 case(1)

  ! rr1 er b�lgens(og dermed flades) 1. komponemt ganget sammen
  ! rr2 er b�lgens(og dermed flades) 2. komponemt ganget sammen
  select case(eface)
  case(1)
  	rr1 = 0.0d0
    rr2 = -1.0d0
  	eta = -1.0d0
  case(2)
  	rr1 = 1.0d0
    rr2 = 0.0d0
  	xi = 1.0d0
  case(3)
  	rr1 = 0.0d0
    rr2 = 1.0d0
  	eta = 1.0d0
  case(4)
  	rr1 = -1.0d0
    rr2 = 0.0d0
  	xi = -1.0d0
  end select 

  Zmat = 0.0d0
  Zmat(1,1) = dens*cL*rr1+dens*cT*(1.0d0-rr1)
  Zmat(2,2) = dens*cL*rr2+dens*cT*(1.0d0-rr2)

  CZe = 0
  do i=1,ng
    call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
     CZe = Cze+length*w(i)*thk*MATMUL(MATMUL(TRANSPOSE(Nmat),Zmat),Nmat)
  end do

 case(0)! Analytisk udtryk fundet ved integration.
  select case(eface)	
  case(1)
    
    CZe=0.0
    Cze(1,1) = 2.D0/3.D0*length*dens*cT
    Cze(1,3) = length*dens*cT/3
    Cze(3,1) = Cze(1,3)
    Cze(2,2) = 2.D0/3.D0*length*dens*cL
    Cze(2,4) = length*dens*cL/3
    Cze(4,2) = Cze(2,4)
    Cze(3,3) = Cze(1,1)
    Cze(4,4) = Cze(2,2)

  case(2)
  	Cze = 0.0d0
    Cze(3,3) = 2.D0/3.D0*dens*cL*length
    Cze(5,5) = Cze(3,3)
    Cze(3,5) = dens*cL*length/3
    CZe(5,3) = CZe(3,5)
    Cze(4,4) = 2.D0/3.D0*dens*cT*length
    CZe(6,6) = CZe(4,4)
    Cze(4,6) = dens*cT*length/3
    Cze(6,4) = Cze(4,6)
    
  case(3)
  	Cze = 0.0d0
    Cze(5,5) = 2.D0/3.D0*length*dens*cT
    Cze(7,7) = Cze(5,5)
    Cze(5,7) = length*dens*cT/3
    Cze(7,5) = Cze(5,7)
    Cze(6,6) = 2.D0/3.D0*length*dens*cL
    Cze(8,8) = Cze(6,6)
    Cze(6,8) = length*dens*cL/3
    Cze(8,6) = Cze(6,8)

  case(4)
  	Cze = 0.0d0
    Cze(1,1) = 2.D0/3.D0*dens*cL*length
    Cze(7,7) = Cze(1,1)
    Cze(1,7) = dens*cL*length/3
    Cze(7,1) = Cze(1,7) 
    Cze(2,2) = 2.D0/3.D0*dens*cT*length
    Cze(8,8) = Cze(2,2)
    Cze(2,8) = dens*cT*length/3
    Cze(8,2) = Cze(2,8)

  end select
  
    CZe = CZe*thk
 end select

	
 END SUBROUTINE plane42transient_CZe


SUBROUTINE plane42transient_CZe_punkt(xe, parameters,eface,ng,rho, CZe)


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ng,eface ! number of gauss-points
  REAL(8), INTENT(IN) :: rho, parameters(:)
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: CZe
	
  INTEGER :: i, j
  
  REAL(8) :: jac(2,2), detjac, bmat(3,8), Nmat(2,8), Zmat(2,2), helpmat1(2,2),helpmat2(2,2)
  REAL(8) :: young1,young2,dens1,dens2, nu1, nu2, thk, length, r, rvec(2), nvec(2)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
   real(8) :: dens, young, nu, lambda, mu, cL, cT
  

  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  young1 = parameters(1)
  young2 = parameters(2)
  dens1 = parameters(3)
  dens2 = parameters(4)
  nu1 = parameters(5)
  nu2 = parameters(6)
  thk = parameters(7)
  length = parameters(8)
  r = parameters(9)
  rvec(1) = parameters(10)
  rvec(2) = parameters(11)
  nvec(1) = parameters(12)
  nvec(2) = parameters(13)


  dens = rho*dens1+(1.0-rho)*dens2
  young = rho*young1+(1.0-rho)*young2
  nu = rho*nu1+(1.0-rho)*nu2
  !Lame parameters  
  lambda = young*nu/ ((1.0+nu)*(1.0-2.0*nu)) 
  mu = young/(2.0*(1.0+nu))! shear modulus
  cL = dsqrt( (lambda+2.0*mu) / dens)
  cT = dsqrt( mu / dens)


  helpmat1 = 0.0d0
  do i=1,2
  	do j=1,2
    	helpmat1(i,j) = rvec(i)*rvec(j)
    end do
  end do
  Zmat = 0.0d0
  Zmat(1,1) = 1.0d0
  Zmat(2,2) = 1.0d0
  Zmat = mu/cT*dot_product(nvec,rvec) * Zmat - 2*mu*(1/cT-1/cL)*dot_product(nvec,rvec)*helpmat1
  do i=1,2
  	do j=1,2
    	helpmat1(i,j) = nvec(i)*rvec(j)
    end do
  end do
  Zmat = Zmat + lambda/cL * helpmat1
  
  do i=1,2
  	do j=1,2
    	helpmat1(i,j) = rvec(i)*nvec(j)
    end do
  end do
  Zmat = Zmat + mu/cT * helpmat1
  
  select case(eface)
  case(1)
  	eta = -1.0d0
  case(2)
  	xi = 1.0d0
  case(3)
  	eta = 1.0d0
  case(4)
  	xi = -1.0d0
  end select

  CZe = 0
  do i=1,ng
    call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
     CZe = Cze+length*w(i)*thk*MATMUL(MATMUL(TRANSPOSE(Nmat),Zmat),Nmat)
  end do

 END SUBROUTINE plane42transient_CZe_punkt

  
 SUBROUTINE plane42transient_KRe(xe, parameters,eface,ng,rho, KRe)

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ng, eface
  REAL(8), INTENT(IN) :: rho, parameters(:)
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: KRe

	
  INTEGER :: i, j
   
  REAL(8) :: jac(2,2), detjac, bmat(3,8), Nmat(2,8), Rmat(2,2), helpmat1(2,2),helpmat2(2,2)
  REAL(8) :: young1,young2,dens1,dens2, nu1, nu2, thk, length, r, rvec(2), nvec(2)
  REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
   real(8) :: dens, young, nu, lambda, mu, cL, cT
  

  allocate(w(ng),xi(ng),eta(ng))
  call gauss_points(w,xi,eta,ng)

  young1 = parameters(1)
  young2 = parameters(2)
  dens1 = parameters(3)
  dens2 = parameters(4)
  nu1 = parameters(5)
  nu2 = parameters(6)
  thk = parameters(7)
  length = parameters(8)
  r = parameters(9)
  rvec(1) = parameters(10)
  rvec(2) = parameters(11)
  nvec(1) = parameters(12)
  nvec(2) = parameters(13)

  dens = rho*dens1+(1.0-rho)*dens2
  young = rho*young1+(1.0-rho)*young2
  nu = rho*nu1+(1.0-rho)*nu2
  !Lame parameters  
  lambda = young*nu/ ((1.0+nu)*(1.0-2.0*nu)) 
  mu = young/(2.0*(1.0+nu))! shear modulus
  cL = dsqrt( (lambda+2.0*mu) / dens)
  cT = dsqrt( mu / dens)
 
  helpmat1 = 0.0d0
  do i=1,2
  	do j=1,2
    	helpmat1(i,j) = nvec(i)*rvec(j)
    end do
  end do
  
  Rmat = 0.0d0
  Rmat(1,1) = 1.0d0
  Rmat(2,2) = 1.0d0
  Rmat = 2.0d0*mu/r*dot_product(nvec,rvec) * Rmat - (lambda+ 2.0d0*mu)/r*helpmat1
  do i=1,2
  	do j=1,2
    	helpmat1(i,j) = rvec(i)*nvec(j)
    end do
  end do
  Rmat = Rmat + 2.0d0*mu/r* helpmat1
  
  select case(eface)
  case(1)
  	eta = -1.0d0
  case(2)
  	xi = 1.0d0
  case(3)
  	eta = 1.0d0
  case(4)
  	xi = -1.0d0
  end select
  
  KRe = 0.0
  do i=1,ng
    call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
    KRe = KRe+length*w(i)*thk*MATMUL(MATMUL(TRANSPOSE(Nmat),Rmat),Nmat)
  end do
	
 END SUBROUTINE plane42transient_KRe


END MODULE plane42transient
