MODULE plane41

  ! This module contains subroutines specific to the PLANE41 element.
  ! that is, isoparametric Q4 element with 1 DOF pr. node. Use for thermal problem.

  use numeth

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: plane41_ke_t, plane41_he_t, plane41_rbe_t, plane41_rhe_t, plane41_rqe_t, &
       plane41_temp, plane41_shape, plane41_dielectric

CONTAINS


  SUBROUTINE plane41_dielectric(xe, mat_vec, thk,ng, ke,phi)

    ! This subroutine constructs the element conductivity matrix [k_e].

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(4,4), INTENT(OUT) :: ke
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j
    REAL(8) :: epsilon(2,2)
    REAL(8) :: Bmat(2,4),  Nvec(4), jac(2,2), detjac
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp(2,2) ! rotation
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    epsilon = 0d0
    epsilon(1,1) = mat_vec(10)
    epsilon(2,2) = mat_vec(11)+mat_vec(7)**2/mat_vec(1)

   if (present(phi)) then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp = MATMUL(TRANSPOSE(Lambda),epsilon)
       epsilon = MATMUL(helpp,Lambda)
    end if

    ke = 0.0d0
    do i = 1 , ng
       do j = 1 , ng
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat, detjac, jac)
          ke = ke - w(i)*w(j) *thk* MATMUL(MATMUL(TRANSPOSE(bmat),epsilon), bmat) * detjac ! minus is intended!
       end do
    end do


  END SUBROUTINE plane41_dielectric



  SUBROUTINE plane41_ke_t(xe,ng, kcond, thk, ke_t)

    ! This subroutine constructs the element conductivity matrix [k_e].

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: kcond, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(4,4), INTENT(OUT) :: ke_t

    INTEGER :: i, j
    REAL(8) :: kappa(2,2)
    REAL(8) :: Bmat(2,4),  Nvec(4), jac(2,2), detjac
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    kappa = 0.0d0
    kappa(1,1) = kcond
    kappa(2,2) = kappa(1,1)

    ke_t = 0.0d0
    do i = 1 , ng
       do j = 1 , ng
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat, detjac, jac)
          ! Element thermal stiffness matrix
          ke_t = ke_t + w(i)*w(j) *thk* MATMUL(MATMUL(TRANSPOSE(bmat),kappa), bmat) * detjac
       end do
    end do


  END SUBROUTINE plane41_ke_t

  SUBROUTINE plane41_he_t(xe, eface, hconv, thk, he_t)

    ! This subroutine constructs the element boundary convection matrix [h_e].

    INTEGER, INTENT(IN) :: eface
    REAL(8), INTENT(IN) :: hconv, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: he_t
    INTEGER :: i, j
    REAL(8) :: hmat(4,4), length


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

    IF (eface == 1) THEN
       i = 1
       j = 2
    ELSEIF (eface == 2) THEN
       i = 2
       j = 3
    ELSEIF (eface == 3) THEN
       i = 3
       j = 4
    ELSEIF (eface == 4) THEN
       i = 4
       j = 1
    end if

    ! he_t fås ved eksakt integration, se s. 461 i COOK, af het = \int_{-1}^{1} N^T N h ,dS

    he_t = 0.0d0
    hmat = 0.0d0
    length = dsqrt( ( xe(2*j-1) - xe(2*i-1) )**2 + ( xe(2*j) - xe(2*i) )**2 )
    hmat(i,i) = 2.0d0
    hmat(j,j) = 2.0d0
    hmat(i,j) = 1.0d0
    hmat(j,i) = 1.0d0
    he_t = thk*hconv*length*hmat/6.0d0

  END SUBROUTINE plane41_he_t

  SUBROUTINE plane41_rbe_t(xe, eface, fb, thk, rbe_t)

    ! This subroutine constructs the element heat flux vector {r_B_e}.

    INTEGER, INTENT(IN) :: eface
    REAL(8), INTENT(IN) :: fb, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: rbe_t(4)
    INTEGER :: i, j
    REAL(8) :: rbevec(4), length

    IF (eface == 1) THEN
       i = 1
       j = 2
    ELSEIF (eface == 2) THEN
       i = 2
       j = 3
    ELSEIF (eface == 3) THEN
       i = 3
       j = 4
    ELSEIF (eface == 4) THEN
       i = 4
       j = 1
    end if

    ! rbe_t fås ved eksakt integration, se s. 462 i COOK, af rbe_t = \int_{-1}^{1} N^T fb ,dS

    rbe_t = 0.0d0
    rbevec = 0.0d0

    length = dsqrt( ( xe(2*j-1) - xe(2*i-1) )**2 + ( xe(2*j) - xe(2*i) )**2)
    rbevec(i) = 1.0d0
    rbevec(j) = 1.0d0
    rbe_t = fb*thk*length*rbevec/2.0d0

  END SUBROUTINE plane41_rbe_t

  SUBROUTINE plane41_rhe_t(xe, eface, t_inf, thk, hconv, rhe_t)

    ! This subroutine constructs the element boundary convection vector {r_h_e}.

    INTEGER, INTENT(IN) :: eface
    REAL(8), INTENT(IN) :: t_inf, thk, hconv
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: rhe_t(4)
    INTEGER :: i, j
    REAL(8) :: rhevec(4), length


    IF (eface == 1) THEN
       i = 1
       j = 2
    ELSEIF (eface == 2) THEN
       i = 2
       j = 3
    ELSEIF (eface == 3) THEN
       i = 3
       j = 4
    ELSEIF (eface == 4) THEN
       i = 4
       j = 1
    end if

    rhe_t = 0.0d0
    rhevec = 0.0d0

    length = sqrt( ( xe(2*j-1) - xe(2*i-1) )**2 + ( xe(2*j) - xe(2*i) )**2)
    rhevec(i) = 1.0d0
    rhevec(j) = 1.0d0
    rhe_t = t_inf*thk*hconv*length*rhevec/2.0d0

  END SUBROUTINE plane41_rhe_t

  SUBROUTINE plane41_rqe_t(xe,ng, qint, thk, rqe_t)

    ! This subroutine constructs the element internal heat generation vector {r_Q_e}

    INTEGER, INTENT(IN) :: ng
    REAL(8), INTENT(IN) :: thk, qint
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: rqe_t(4)
    INTEGER :: i, j
    REAL(8) :: Bmat(2,4),  Nvec(4), jac(2,2), detjac
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    rqe_t = 0.0d0
    do i = 1 , ng
       do j = 1 , ng
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat, detjac, jac)
          rqe_t = rqe_t + w(i)*w(j)*thk*Nvec*qint*detjac
       end do
    end do

  END SUBROUTINE plane41_rqe_t


  SUBROUTINE plane41_temp(xe, tne, te)

    ! This subrotuine constructs the element temperature, te, from the nodal temperature, tne.
    ! Temp evalueres i midten af elementet, jf. xi = 0 og eta = 0.

    REAL(8), DIMENSION(:), INTENT(IN)  :: xe, tne
    REAL(8), INTENT(OUT) :: te

    REAL(8) :: Bmat(2,4),  Nvec(4), jac(2,2), detjac
    REAL(8) :: xi, eta

    xi = 0.0d0 ! Temp evalueres i centrum
    eta = 0.0d0
    call plane41_shape(xe, xi, eta, Nvec, bmat, detjac, jac)

    ! Compute element temperature
    te = DOT_PRODUCT(Nvec,tne) ! (12.2-1) i COOK


  END SUBROUTINE plane41_temp


  SUBROUTINE plane41_shape(xe, xi, eta, Nvec, bmat, detjac, jac)

    ! This subroutine returns several outputs related to the shape functions for elements with
    ! 1 DOF pr node. 

    REAL(8), INTENT(IN) :: xi, eta                    ! New coordinates of the gauss point
    REAL(8), DIMENSION(:), INTENT(IN) :: xe           ! Nodal coordinates of element
    REAL(8), INTENT(OUT) :: detjac                    ! Det of Jacobian matrix
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: Bmat    	! B0
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: jac       ! Jacobian matrix
    REAL(8), DIMENSION(:), INTENT(OUT) :: Nvec		! shapefunktioner, vektor
    REAL(8), DIMENSION(2,2) :: GAMMA					! Invers jacobi-matrix, J^(-1)
    REAL(8)  :: N_xi(4), N_eta(4), Ntilde(2,4)	! Differentiated Shape function matrix Ntilde
    INTEGER :: i

    ! shapefunktion differentieret mht xi og eta
    N_xi(1) = (-1+eta)/4.0d0	 	!N1,xi
    N_xi(2) = (1-eta) /4.00		!N2,xi
    N_xi(3) = (1+eta) /4.0d0		!N3,xi
    N_xi(4) = (-1-eta)/4.0d0		!N4,xi
    N_eta(1)= (-1+xi) /4.0d0		!N1,eta
    N_eta(2)= (-1-xi) /4.0d0		!N2,eta
    N_eta(3)= (1+xi)  /4.0d0		!N3,eta
    N_eta(4)= (1-xi)  /4.0d0		!N4,eta

    ! Shapefunktion for isoparametrisk Q4 element. 
    Nvec(1) = (1-xi)*(1-eta)/4.0d0
    Nvec(2) = (1+xi)*(1-eta)/4.0d0
    Nvec(3) = (1+xi)*(1+eta)/4.0d0
    Nvec(4) = (1-xi)*(1+eta)/4.0d0

    ! Bemærk at ved 1 DOF pr knude er er Ntilde[2x4] i modsætning til Ntilde[4x8] ved 2 DOF pr knude., se evt. s. 57 i uge 8 noter
    Ntilde = 0.0d0
    do i = 1, 4
       Ntilde(1,i) = N_xi(i)
       Ntilde(2,i) = N_eta(i)
    end do

    !$$$$$$   Nmat=0.0d0 ! Skal ikke bruges, så derfor er den udkommenteret
    !$$$$$$   do i = 1, 4
    !$$$$$$     Nmat(1,2*i-1) = Nvec(i)
    !$$$$$$     Nmat(2,2*i)   = Nvec(i)
    !$$$$$$   end do 

    jac = 0.0d0
    do i = 1, 4
       jac(1,1) = jac(1,1) + N_xi(i)*xe(2*i-1)
       jac(1,2) = jac(1,2) + N_xi(i)*xe(2*i)
       jac(2,1) = jac(2,1) + N_eta(i)*xe(2*i-1)
       jac(2,2) = jac(2,2) + N_eta(i)*xe(2*i)
    end do

    detjac = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

    GAMMA = 0.0d0
    GAMMA(1,1) = jac(2,2)
    GAMMA(2,2) = jac(1,1)
    GAMMA(1,2) = -jac(1,2)
    GAMMA(2,1) = -jac(2,1)
    GAMMA = GAMMA*(1.0d0/detjac)

    ! Bemærk at Bmat er defineret anderledes her, end for plane42, se s. 462 (12.2-7) og s 208 (6.2-13)
    Bmat = MATMUL(GAMMA, Ntilde)

  END SUBROUTINE plane41_shape

END MODULE plane41
