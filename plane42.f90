MODULE plane42
  ! This module contains subroutines specific to the PLANE42 element.

  use numeth ! gauss points ligger i numeth
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: plane42_ke,plane42_me, plane42_re, plane42_vol, plane42_ss
  PUBLIC :: plane42_thermLoad, plane42_thermKobling
  PUBLIC :: plane42_shape, plane42_area
  PUBLIC :: plane42_piezo, plane42_ke_piezo

CONTAINS

  SUBROUTINE plane42_piezo(xe, mat_vec, thk,ng, ke,phi)

    ! This subroutine constructs the stiffness matrix for
    ! a rectangular 4-noded quad element.
    use plane41

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke(8,4)
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j
    REAL(8) :: Emat(2,3), helpproduct(8,2)
    real(8) :: detjac, jac(2,2), Nmat(2,8), bmat(3,8) ! output fra plane42_shape
    REAL(8) :: Bmat_t(2,4),  Nvec(4), jac_t(2,2), detjac_t ! output fra plane41_shape
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp(2,3) ! rotation
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    !$$$$$$   REAL(8), DIMENSION(:,:) , allocatable:: helpproduct

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    Emat = 0d0
    Emat(2,1) = mat_vec(7) - mat_vec(2)*mat_vec(7)/mat_vec(1)
    Emat(2,2) = mat_vec(8) - mat_vec(3)*mat_vec(7)/mat_vec(1)
    Emat(1,3) = mat_vec(9)

    if (present(phi)) then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp = MATMUL(TRANSPOSE(Lambda),Emat)
       Emat = MATMUL(helpp,T_p)
    end if

    ke=0d0
    do i= 1,ng !number of gauss points
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat_t, detjac_t, jac_t)
          helpproduct = MATMUL(TRANSPOSE(Bmat),TRANSPOSE(Emat))
          ke =ke+ W(i)*W(j)*thk*MATMUL(helpproduct,bmat_t)*detjac
       end do
    end do

  END SUBROUTINE plane42_piezo

  SUBROUTINE plane42_ke_piezo(xe,mat_vec, young, nu, thk,ng, ke,phi)

    ! This subroutine constructs the stiffness matrix for
    ! a rectangular 4-noded quad element.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: young, nu, thk, mat_vec(:)
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j
    REAL(8) :: fact, Cmat(3,3), detjac, jac(2,2), bmat(3,8), Nmat(2,8), helpproduct(8,3)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp1(3,3) ! rotation
    !$$$$$$   REAL(8), DIMENSION(:,:) , allocatable:: helpproduct

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    ! Build constitutive matrix (plane stress)

    Cmat = 0d0
    Cmat(1, 1) = mat_vec(1) - mat_vec(2)**2/mat_vec(1)
    Cmat(1, 2) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
    Cmat(2, 1) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
    Cmat(2, 2) = mat_vec(4) - mat_vec(3)**2/mat_vec(1)
    Cmat(3, 3) = mat_vec(5)

    if (present(phi)) then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp1 = MATMUL(TRANSPOSE(T_p),Cmat)
       Cmat = MATMUL(helpp1,T_p)  
    end if

!!$    Cmat = 0d0
!!$    fact = young/(1d0-nu**2)
!!$    Cmat(1, 1) = fact
!!$    Cmat(1, 2) = fact*nu
!!$    Cmat(2, 1) = fact*nu
!!$    Cmat(2, 2) = fact
!!$    Cmat(3, 3) = fact*(1d0-nu)/2d0


    ke=0d0

    do i= 1,ng !number of gauss points
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          helpproduct = MATMUL(TRANSPOSE(bmat),Cmat)
          ke =ke+ W(i)*W(j)*thk*MATMUL(helpproduct,bmat)*detjac
       end do
    end do



  END SUBROUTINE plane42_ke_piezo

  SUBROUTINE plane42_ke(xe, young, nu, thk,ng, ke)

    ! This subroutine constructs the stiffness matrix for
    ! a rectangular 4-noded quad element.


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: young, nu, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke

    INTEGER :: i, j
    REAL(8) :: fact, Cmat(3,3), detjac, jac(2,2), bmat(3,8), Nmat(2,8), helpproduct(8,3)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    !$$$$$$   REAL(8), DIMENSION(:,:) , allocatable:: helpproduct

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    ! Build constitutive matrix (plane stress)

    Cmat = 0d0
    fact = young/(1d0-nu**2)
    Cmat(1, 1) = fact
    Cmat(1, 2) = fact*nu
    Cmat(2, 1) = fact*nu
    Cmat(2, 2) = fact
    Cmat(3, 3) = fact*(1d0-nu)/2d0

    ke=0d0

    do i= 1,ng !number of gauss points
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          helpproduct = MATMUL(TRANSPOSE(bmat),Cmat)
          ke =ke+ W(i)*W(j)*thk*MATMUL(helpproduct,bmat)*detjac
       end do
    end do

  END SUBROUTINE plane42_ke

  SUBROUTINE plane42_me(xe,dens,thk, ng, me)

    ! This subroutine constructs the mass matrix for
    ! a rectangular 4-noded quad element.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng
    REAL(8), INTENT(IN) :: dens, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: me

    INTEGER :: i, j
    REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    me = 0.0

	do i= 1,ng !number of gauss points
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          me =me+ W(i)*W(j)*dens*MATMUL(TRANSPOSE(Nmat),Nmat)*detjac*thk
       end do
	end do

  END SUBROUTINE plane42_me

  SUBROUTINE plane42_re(xe, eface, fe, thk, ng, re)
    ! This subroutine assembles the element surface loads.

    INTEGER, INTENT(IN) :: eface, ng
    REAL(8), INTENT(IN) :: fe, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: re(8)

    REAL(8) :: Nmat(2,8), bmat(3,8), jac(2,2), detjac, helpmatrix(2)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    INTEGER:: i

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)
    !      3
    !  l ______k
    !   |     |
    ! 4 |     |  2
    !   |_____|
    !  i       j
    !      1
    !
    ! node numbers: i, j, k, l
    ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)

    re =0.


    IF (eface == 1) THEN
       do i=1,ng
          eta = -1.0
          call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
          helpmatrix(1) = -jac(1,2)
          helpmatrix(2) = jac(1,1)
          re =re+ W(i)*thk*(-fe)*MATMUL(TRANSPOSE(Nmat),helpmatrix)
       end do

    ELSEIF (eface == 2) THEN
       do i=1,ng
          xi = 1.0
          call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
          helpmatrix(1) = -jac(2,2)
          helpmatrix(2) = jac(2,1)
          re =re+W(i)*thk*(fe)*MATMUL(TRANSPOSE(Nmat),helpmatrix)
       end do
    ELSEIF (eface == 3) THEN
       do i=1,ng
          eta = 1.0
          call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
          helpmatrix(1) = jac(1,2)
          helpmatrix(2) = -jac(1,1)

          re =re+ W(i)*thk*(fe)*MATMUL(TRANSPOSE(Nmat),helpmatrix)
       end do
    ELSEIF (eface == 4) THEN
       do i=1,ng
          xi = -1.0
          call plane42_shape(xe, xi(i), eta(i),Nmat, bmat, jac, detjac)
          helpmatrix(1) = jac(2,2)
          helpmatrix(2) = -jac(2,1)
          re =re+ W(i)*thk*(-fe)*MATMUL(TRANSPOSE(Nmat),helpmatrix)
       end do
    ENDIF

  END SUBROUTINE plane42_re

  SUBROUTINE plane42_vol(xe, accel, dens, thk,ng , rvol)

    ! This subroutine assembles the element surface loads.
    ! accel er acceleration af element, dens er densitet


    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk, accel(2), dens
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: rvol(8)

    INTEGER :: i, j
    REAL(8) :: phi(2), detjac, jac(2,2), bmat(3,8), Nmat(2,8)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    !      3
    !  l ______k
    !   |     |
    ! 4 |     |  2
    !   |_____|
    !  i       j
    !      1
    !
    ! node numbers: i, j, k, l
    ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)


    ! Volumenkrafts-intensitet
    phi(1) = dens*accel(1)
    phi(2) = dens*accel(2)
    rvol = 0.0
    do i= 1,ng !number of gauss points
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          rvol =rvol+ W(i)*W(j)*thk*MATMUL(TRANSPOSE(Nmat),phi)*detjac
       end do
    end do

  END SUBROUTINE plane42_vol


  SUBROUTINE plane42_ss(xe, de, young, nu, estress, estrain)

    ! This subrotuine constructs the element stress and strain

    REAL(8), INTENT(IN) :: young, nu!, de(8,1)
    REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
    REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain

    REAL(8) :: fact, bmat(3, 8), Cmat(3, 3), Nmat(2,8), jac(2,2)
    REAL(8) :: xi, eta, detjac


    ! Build strain-displacement matrix ! for controid(midten) , dvs xi=eta=0
    xi = 0.0 ! Tøjning evalueres i centrum
    eta = 0.0
    call plane42_shape(xe, xi, eta, Nmat,bmat, jac, detjac)

    estrain = 0.0
    ! Compute element strain
    estrain = MATMUL(bmat, de)

    ! Build constitutive matrix (plane stress)
    Cmat = 0d0
    fact = young/(1d0-nu**2)
    Cmat(1, 1) = fact
    Cmat(1, 2) = fact*nu
    Cmat(2, 1) = fact*nu
    Cmat(2, 2) = fact
    Cmat(3, 3) = fact*(1d0-nu)/2d0

    ! Compute element stress
    estress = MATMUL(Cmat, estrain)

    ! Compute principal stress and direction

  END SUBROUTINE plane42_ss


  SUBROUTINE plane42_thermLoad(xe, ng, young, nu, alpha, thk, deltaT, rThermal)

    ! This subroutine calcutale the nodal loads from element temperatures calculated with thermal routine.
    ! Used in coupled problems
    ! rThermal = int( B^T C epsilon deltaT )dV

    INTEGER, INTENT(IN) :: ng
    REAL(8), INTENT(IN) :: young, thk, alpha, nu, deltaT
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: rThermal(8)

    INTEGER :: i, j
    REAL(8) :: detjac, jac(2,2), Bmat(3,8), Nmat(2,8)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    REAL(8) :: Cmat(3,3), tstart, tend, epsilon(3), fact, helpproduct(8,3)

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)


    ! Build constitutive matrix (plane stress)
    Cmat = 0d0
    fact = young/(1d0-nu**2)
    Cmat(1, 1) = fact
    Cmat(1, 2) = fact*nu
    Cmat(2, 1) = fact*nu
    Cmat(2, 2) = fact
    Cmat(3, 3) = fact*(1d0-nu)/2d0

    epsilon = 0.0d0
    ! Build thermal strain vector
    epsilon(1:2) = alpha

    rThermal = 0.0d0

    do i = 1 , ng
       do j = 1 , ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          rThermal = rThermal + w(i)*w(j)*thk*MATMUL(MATMUL(TRANSPOSE(Bmat),Cmat), epsilon)*deltaT*detjac
       end do
    end do

    !$$$$$$     do i= 1,ng
    !$$$$$$         do j = 1,ng
    !$$$$$$             call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
    !$$$$$$             helpproduct = MATMUL(TRANSPOSE(bmat),Cmat)
    !$$$$$$             rThermal =rThermal+ W(i)*W(j)*thk*deltaT*MATMUL(helpproduct,epsilon)*detjac
    !$$$$$$         end do
    !$$$$$$     end do


  END SUBROUTINE plane42_thermLoad


  SUBROUTINE plane42_thermKobling(xe, ng, young, nu, alpha, thk,Ae)

    ! This subroutine calculates koblingsmatricen 
    ! Ae = int( B^T C epsilon Nt )dV
    use plane41
    use numeth

    INTEGER, INTENT(IN) :: ng
    REAL(8), INTENT(IN) :: young, thk, alpha, nu
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    !$$$$$$   REAL(8), INTENT(OUT) :: Ae(8,4)
    REAL(8), INTENT(OUT) :: Ae(4,8)

    INTEGER :: i, j
    REAL(8) :: detjac, jac(2,2), Bmat(3,8), Nmat(2,8)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    REAL(8) :: Cmat(3,3), tstart, tend, deltaT, epsilon(3), fact, helpproduct(8)!, helpproduct2(8,4)
    REAL(8) :: Bmat_t(2,4),  Nvec(4), jac_t(2,2), detjac_t ! output fra plane41_shape
    real(8) :: helpproduct2(4,8)

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)


    ! Build constitutive matrix (plane stress)
    Cmat = 0d0
    fact = young/(1d0-nu**2)
    Cmat(1, 1) = fact
    Cmat(1, 2) = fact*nu
    Cmat(2, 1) = fact*nu
    Cmat(2, 2) = fact
    Cmat(3, 3) = fact*(1d0-nu)/2d0


    ! Build thermal strain vector
    epsilon = 0.0d0
    epsilon(1:2) = alpha

    Ae = 0.0d0

    do i = 1 , ng
       do j = 1 , ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat_t, detjac_t, jac_t)
          helpproduct = MATMUL(MATMUL(TRANSPOSE(Bmat),Cmat), epsilon)
          !$$$$$$       call vector_mul(helpproduct,Nvec,helpproduct2) ! dimension af rThermal afhænger af rækkefølgen af input-vektorerne til vector_mul
          !$$$$$$       Ae = Ae + W(i)*W(j)*thk* helpproduct2 *detjac

          call vector_mul(Nvec,helpproduct,helpproduct2) ! dimension af rThermal afhænger af rækkefølgen af input-vektorerne til vector_mul
          Ae = Ae + W(i)*W(j)*thk* helpproduct2 *detjac
       end do
    end do


  END SUBROUTINE plane42_thermKobling


  SUBROUTINE plane42_shape(xe, xi, eta, Nmat,bmat, jac, detjac)
    ! This subrotuine implement isoparametric elements

    REAL(8), INTENT(IN) :: xi, eta							! New coordinates of the gauss point
    REAL(8), DIMENSION(:), INTENT(IN)  :: xe					! Nodal coordinates of element
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: bmat				! B0
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: jac				! Jacobian matrix for plane integration
    REAL(8), INTENT(OUT) :: Nmat(2,8)							! Shape function matrix
    REAL(8), INTENT(OUT) :: detjac							! Determinant of Jacobian matrix

    REAL(8) :: Ntilde(4,8), GAMMA(2,2), GAMMAtilde(4,4), N1,N2,N3,N4 , L(3,4), produkthelp(3,4)
    Nmat = 0
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

    L(1,1) = 1 !(8.6) i noterne
    L(2,1) = 0
    L(3,1) = 0

    L(1,2) = 0
    L(2,2) = 0
    L(3,2) = 1

    L(1,3) = 0
    L(2,3) = 0
    L(3,3) = 1

    L(1,4) = 0
    L(2,4) = 1
    L(3,4) = 0

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

    GAMMAtilde(3,1) = 0
    GAMMAtilde(4,1) = 0
    GAMMAtilde(3,2) = 0
    GAMMAtilde(4,2) = 0

    GAMMAtilde(1,3) = 0
    GAMMAtilde(2,3) = 0
    GAMMAtilde(1,4) = 0
    GAMMAtilde(2,4) = 0    

    GAMMAtilde(3,3) = GAMMA(1,1)
    GAMMAtilde(4,3) = GAMMA(2,1)
    GAMMAtilde(3,4) = GAMMA(1,2)
    GAMMAtilde(4,4) = GAMMA(2,2)

    Ntilde = 0
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

    produkthelp = MATMUL(L, GAMMAtilde) !hjælpeprodukt ved B=[L]*[GAMMAtilde]*[Ntilde]
    bmat = MATMUL(produkthelp, Ntilde)

  END SUBROUTINE plane42_shape


  SUBROUTINE plane42_area(xe,ng , area)


    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: area

    INTEGER :: i, j
    REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)


    area = 0.0
    do i= 1,ng
       do j = 1,ng
          call plane42_shape(xe, xi(i), eta(j),Nmat, bmat, jac, detjac)
          area = area + w(i)*w(j)*detjac
       end do
    end do

  END SUBROUTINE plane42_area


END MODULE plane42
