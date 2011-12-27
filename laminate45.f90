MODULE laminate45

  ! This module contains subroutines specific to the MINDLIN41 element.

  use numeth
  IMPLICIT NONE

  PRIVATE :: shape, constitutive_matrix,piezo_matrix,dielectric_matrix
  PUBLIC :: laminate45_ke, laminate45_re
  PUBLIC :: laminate45_piezo

CONTAINS

  SUBROUTINE laminate45_piezo(xe,int_parameters,real_parameters,mat_vec, thk,ng, ke)

    ! This subroutine constructs the stiffness matrix for
    ! a rectangular 4-noded quad element.
    use plane41

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke !(20,4)
    integer, intent(IN) :: int_parameters(:)
    real(8), intent(in) :: real_parameters(:)

    INTEGER :: i, j
    REAL(8) :: helpproduct(20,2), ke_m(20,4),ke_b(20,4)

    REAL(8) :: detjac,jac(2,2), Nmat(3,12),Bm(3,20),Bb(3,20),Bs(2,20) ! output fra laminate_shape
    REAL(8) :: Bmat_t(2,4),  Nvec(4), jac_t(2,2), detjac_t ! output fra plane41_shape
    REAL(8) :: H1mat(2,3), H2mat(2,3)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    !$$$$$$   REAL(8), DIMENSION(:,:) , allocatable:: helpproduct

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    call piezo_matrix(int_parameters,real_parameters,mat_vec,H1mat,H2mat)
    
    ke_m=0d0
    ke_b=0d0
    do i = 1, ng
       do j = 1, ng
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat_t, detjac_t, jac_t)
          call shape(xe,xi(i),eta(j),Bm,Bb,Bs,jac,detjac,Nmat)

          helpproduct = MATMUL(TRANSPOSE(Bm),TRANSPOSE(H1mat))
          ke_m = ke_m + W(i)*W(j)*thk*MATMUL(helpproduct,bmat_t)*detjac

          helpproduct = MATMUL(TRANSPOSE(Bb),TRANSPOSE(H2mat))
          ke_b = ke_b + W(i)*W(j)*thk*MATMUL(helpproduct,bmat_t)*detjac

          ! Notice Bs is not used since H3mat is zero
       end do
    end do

    ke = ke_m + ke_b

!!$    do i= 1,ng !number of gauss points
!!$       do j = 1,ng
!!$          call shape(xe,xi(i),eta(j),Bm,Bf,Bs,jac,detjac,Nmat)
!!$          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat_t, detjac_t, jac_t)
!!$          helpproduct = MATMUL(TRANSPOSE(Bmat),TRANSPOSE(Emat))
!!$          ke =ke+ W(i)*W(j)*thk*MATMUL(helpproduct,bmat_t)*detjac
!!$       end do
!!$    end 
  END SUBROUTINE laminate45_piezo


  SUBROUTINE laminate45_ke(xe, young, nu, dens, thk, shear,int_parameters,&
       real_parameters,ng, ke,type,mat_vec)


    ! This subroutine constructs the stiffness matrix for
    ! a quadrilateral 4-noded element.

    INTEGER, INTENT(IN) :: ng,type
    REAL(8), INTENT(IN) :: young, nu, dens, thk
    real(8), intent(INOUT) :: shear
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:), OPTIONAL, INTENT(IN) ::  mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke!ke(20,20)
    integer, intent(IN) :: int_parameters(:)
    real(8), intent(in) :: real_parameters(:)

    INTEGER :: i, j
    REAL(8) :: ke_m(20,20),ke_mb(20,20), ke_b(20,20), ke_s(20,20), xi_1, &
         eta_1, w_1
    REAL(8) :: Cmat(3,3), Dmat(3,3), k, Ds_mat(2,2)
    REAL(8) :: detjac,jac(2,2), Nmat(3,12),Bm(3,20),Bb(3,20),Bs(2,20)
    
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)



    ke = 0d0
    ke_s = 0d0
    ke_b = 0d0
    ke_mb = 0d0
 
    do i = 1, ng
       do j = 1, ng
          call shape(xe,xi(i),eta(j),Bm,Bb,Bs,jac,detjac,Nmat)
          ke_m = ke_m + w(i)*w(j)*MATMUL(MATMUL(TRANSPOSE(Bm),Dmat),Bm)*detjac
          ke_b = ke_b + w(i)*w(j)*MATMUL(MATMUL(TRANSPOSE(Bb),Dmat),Bb)*detjac
          ke_mb = ke_mb + w(i)*w(j)*MATMUL(MATMUL(TRANSPOSE(Bm),Dmat),Bb)*detjac
       end do
    end do

    ! For at undgå shear-locking, skal k_s integreres med "single Gauss kvadrature"(dvs ng = 1), jf
    ! "Element Behavior", s 543 i cook
    xi_1 = 0d0
    eta_1 = 0d0
    w_1 = 2d0
    call shape(xe,xi_1,eta_1,Bm,Bb,Bs,jac,detjac,Nmat)
    ke_s = w_1*w_1*MATMUL(MATMUL(TRANSPOSE(Bs),Ds_mat),Bs)*detjac
    ke = ke_m + ke_b + ke_mb + TRANSPOSE(ke_mb) + ke_s

  END SUBROUTINE laminate45_ke

  subroutine constitutive_matrix(int_parameters,real_parameters, &
       mat_vec,Amat,Bmat,Dmat,G3mat)
    ! This sunroutine calculates the element constitutive maxtrix for laminates
    
    real(8), intent(OUT), DIMENSION(3,3) :: Amat,Bmat,Dmat
    real(8), intent(OUT) :: G3mat(2,2)
    integer, intent(IN) :: int_parameters(:)
    real(8), intent(in) :: real_parameters(:), mat_vec(:)
    
    real(8) :: Cmat(3,3), Cs_mat(2,2), thk,nu,young,shear,k
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i,type
    ! nc: number of layers, z: layer-coordinate
    nc = int_parameters(1)
    type = int_parameters(2)
    nu = real_parameters(1)
    young = real_parameters(2)
    shear = real_parameters(3)

    allocate(z(nc+1))
    do i=1,nc+1
       z(i) = real_parameters(i+3)
    end do

    ! værdien af faktoren k, står på s 535 i cook
    k = 5d0/6d0
    
    Cmat = 0d0
    Cs_mat = 0d0
    Amat = 0d0
    Bmat = 0d0
    Dmat = 0d0
    G3mat = 0d0
  
    select case (type)
    case (1) ! normalt isotropisk materiale
       ! Build constitutive matrix (plane stress)
       Cmat(1,1) = 1d0
       Cmat(1,2) = nu
       Cmat(2,1) = nu
       Cmat(2,2) = 1d0
       Cmat(3,3) = (1d0-nu)/2d0
       Cmat = young/(1d0-nu**2)*Cmat

       if (shear < 1E-6 ) then! in case it's forgotten to give shear a value in input-file
          shear = young/(2*(1+nu))
          print*,'shear er ikke angivet i inputfil. Værdi er udregnet, mindlin41.f90'
       end if

       ! shear stiffness
       Cs_mat(1,1) = k*shear
       Cs_mat(2,2) = k*shear

    case (2) ! piezoelectric material
       Cmat(1, 1) = mat_vec(1) - mat_vec(2)**2/mat_vec(1)
       Cmat(1, 2) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
       Cmat(2, 1) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
       Cmat(2, 2) = mat_vec(4) - mat_vec(3)**2/mat_vec(1)
       Cmat(3, 3) = mat_vec(5)

!!$       !shear stiffness
!!$       D_M(4,4) = mat_vec(5)*thk*k
!!$       D_M(5,5) = mat_vec(6)*thk*k
    end select

    do i=1,nc
       thk = z(i+1) - z(i)
       ! z_bar = (1d0/2d0 * (z(i+1) - z(i)) ) ! = position of the middle of layer i.

       Amat = Amat + Cmat*thk
       Bmat = Bmat + Cmat*thk* (1d0/2d0 * (z(i+1) - z(i)) )
       Dmat = Dmat + Cmat*1d0/3d0* (z(i+1)**3 - z(i)**3)
       G3mat = G3mat + Cs_mat*thk
    end do


  end subroutine constitutive_matrix

  subroutine piezo_matrix(int_parameters,real_parameters,&
       mat_vec,H1mat,H2mat)
    ! This sunroutine calculates the element electric  constitutive maxtrix for laminates

    real(8), intent(OUT) :: H1mat(2,3), H2mat(2,3)
    integer, intent(IN) :: int_parameters(:)
    real(8), intent(in) :: real_parameters(:), mat_vec(:)
    
    real(8) :: Emat(2,3), thk
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i
    ! nc: number of layers, z: layer-coordinate
    nc = int_parameters(1)
    allocate(z(nc+1))

    do i=1,nc+1
       z(i) = real_parameters(i)
    end do

    H1mat = 0d0
    H2mat = 0d0

    Emat = 0d0
    Emat(2,1) = mat_vec(7) - mat_vec(2)*mat_vec(7)/mat_vec(1)
    Emat(2,2) = mat_vec(8) - mat_vec(3)*mat_vec(7)/mat_vec(1)
    Emat(1,3) = mat_vec(9)

    do i=1,nc
       thk = z(i+1) - z(i)
       
       H1mat = H1mat + Emat*thk
       H2mat = H2mat + Emat*thk*(1d0/2d0 * (z(i+1) - z(i)) )
    end do

  end subroutine piezo_matrix

  subroutine dielectric_matrix(int_parameters,real_parameters,&
       mat_vec,Jmat)
    ! This sunroutine calculates the element electric  constitutive maxtrix for laminates

    real(8), intent(OUT) :: Jmat(2,2)
    integer, intent(IN) :: int_parameters(:)
    real(8), intent(in) :: real_parameters(:), mat_vec(:)
    
    real(8) :: thk, epsilon(2,2)
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i
    ! nc: number of layers, z: layer-coordinate
    nc = int_parameters(1)
    allocate(z(nc+1))

    do i=1,nc+1
       z(i) = real_parameters(i)
    end do

    Jmat = 0d0

    epsilon = 0d0
    epsilon(1,1) = mat_vec(10)
    epsilon(2,2) = mat_vec(11)+mat_vec(7)**2/mat_vec(1)
    
    do i=1,nc
       thk = z(i+1) - z(i)
       
       Jmat = Jmat + thk
    end do

  end subroutine dielectric_matrix

  SUBROUTINE laminate45_re(xe, eface, fe, thk, re)

    ! This subroutine assembles the element surface loads.

    INTEGER, INTENT(IN) :: eface
    REAL(8), INTENT(IN) :: fe, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(OUT) :: re(12)

    REAL(8) :: aa, bb

    !      
    !  l ______k
    !   |     |
    !   |  1  |  
    !   |_____|
    !  i       j
    !      
    !
    ! node numbers: i, j, k, l
    ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)

    aa = (xe(4)-xe(1))/2d0
    bb = (xe(11)-xe(2))/2d0
    re = 0d0

    !Påført moment i x-retning
    IF (eface == 1) THEN
       re(2) = fe
       re(5) = fe

       !Påført moment i y-retning
    ELSEIF (eface == 2) THEN
       re(6) = fe
       re(9) = fe

    ELSEIF (eface == 3) THEN
       re(8) = fe
       re(11)= fe

    ELSEIF (eface == 4) THEN
       re(12)= fe
       re(3) = fe

       !fordelt last på latteral side
    ELSEIF (eface == 5) THEN
       !Forces in z-dir.
       re(1) = aa*bb*fe !Node 1
       re(4) = aa*bb*fe !Node 2
       re(7) = aa*bb*fe !Node 3
       re(10)= aa*bb*fe !Node 4

       !Moments x-dir.
       re(2) = aa*bb*fe*bb !Node 1
       re(5) = aa*bb*fe*bb !Node 2
       re(8) = -aa*bb*fe*bb !Node 3
       re(11)= -aa*bb*fe*bb !Node 4

       !Moments y-dir.
       re(3) = -aa*bb*fe*aa !Node 1
       re(6) = aa*bb*fe*aa !Node 2
       re(9) = aa*bb*fe*aa !Node 3
       re(12)= -aa*bb*fe*aa !Node 4

    ELSE
       write (*, *) 'Error: Wrong face number!!'
       stop

    ENDIF

!!!!
!!!!! GANG MED TYKKELSE!!!
!!!!!
    re = re*thk

  END SUBROUTINE laminate45_re


  SUBROUTINE shape(xe,xi,eta,B_m,B_b,B_s,jac,detjac,N)

    ! This subroutine calculates the shape function matrix, the Jacobian
    ! matrix and the determinant of the Jacobian matrix.

    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(IN) :: xi, eta
    REAL(8), DIMENSION(:,:), INTENT(OUT) ::  B_m(3,20),B_b(3,20),B_s(2,20) , jac(2,2), N
    REAL(8), INTENT(OUT) :: detjac

    REAL(8) :: N1, N2, N3, N4, Nvec(4)
    REAL(8), DIMENSION(3,4) :: L
    REAL(8), DIMENSION(2,4) :: dN, dNxy
    REAL(8), DIMENSION(2,2) :: Gamma
    REAL(8), DIMENSION(4,4) :: Gtilde
    REAL(8), DIMENSION(3,4) :: B1
    INTEGER :: i


    ! Calculate shape functions
    N1 = 1d0/4d0*(1d0-xi)*(1d0-eta)
    N2 = 1d0/4d0*(1d0+xi)*(1d0-eta)
    N3 = 1d0/4d0*(1d0+xi)*(1d0+eta)
    N4 = 1d0/4d0*(1d0-xi)*(1d0+eta)

    !N_i,eta
    dN(1,1) = -1d0/4d0*(1d0-eta)
    dN(1,2) = 1d0/4d0*(1d0-eta)
    dN(1,3) = 1d0/4d0*(1d0+eta)
    dN(1,4) = -1d0/4d0*(1d0+eta)
    dN(2,1) = -1d0/4d0*(1d0-xi)
    dN(2,2) = -1d0/4d0*(1d0+xi)
    dN(2,3) = 1d0/4d0*(1d0+xi)
    dN(2,4) = 1d0/4d0*(1d0-xi)

    Nvec = [N1,N2,N3,N4]

    ! NB SIKKERT FORKERT: TIL MINDLIN
    n = 0.d0
    n(1,1) = N1
    n(1,4) = N2
    n(1,7) = N3
    n(1,10)= N4
    n(2,2) = N1
    n(2,5) = N2
    n(2,8) = N3
    n(2,11)= N4
    n(3,3) = N1
    n(3,6) = N2
    n(3,9) = N3
    n(3,12)= N4
    ! Calculate Jacobian
    jac = 0.d0
    do i=1,4
       jac(1,1) = jac(1,1) + dN(1,i)*xe(1+(i-1)*2)! 1,4,7...
       jac(1,2) = jac(1,2) + dN(1,i)*xe(2+(i-1)*2)! 2,5,8
       jac(2,1) = jac(2,1) + dN(2,i)*xe(1+(i-1)*2)
       jac(2,2) = jac(2,2) + dN(2,i)*xe(2+(i-1)*2)
    end do

    ! Calculate determinant of Jacobian
    detjac = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

    Gamma(1,1) = jac(2,2)/detjac
    Gamma(1,2) = -jac(1,2)/detjac
    Gamma(2,1) = -jac(2,1)/detjac
    Gamma(2,2) = jac(1,1)/detjac

    ! omregn fra isoparametrisk koordinat til fysisk, jf (15.3-5) i cook
    dNxy(1,1) = Gamma(1,1)*dN(1,1)+Gamma(1,2)*dN(2,1)
    dNxy(1,2) = Gamma(1,1)*dN(1,2)+Gamma(1,2)*dN(2,2)
    dNxy(1,3) = Gamma(1,1)*dN(1,3)+Gamma(1,2)*dN(2,3)
    dNxy(1,4) = Gamma(1,1)*dN(1,4)+Gamma(1,2)*dN(2,4)
    dNxy(2,1) = Gamma(2,1)*dN(1,1)+Gamma(2,2)*dN(2,1)
    dNxy(2,2) = Gamma(2,1)*dN(1,2)+Gamma(2,2)*dN(2,2)
    dNxy(2,3) = Gamma(2,1)*dN(1,3)+Gamma(2,2)*dN(2,3)
    dNxy(2,4) = Gamma(2,1)*dN(1,4)+Gamma(2,2)*dN(2,4)

    ! B splittes op i membran, bending og shear komponent
    
    B_m = 0d0
    B_b = 0d0
    B_s = 0d0
    do i=1,4
       ! membran
       B_m(1,1+(i-1)*5) = dNxy(1,i)
       B_m(2,2+(i-1)*5) = dNxy(2,i)
       B_m(3,1+(i-1)*5) = dNxy(2,i)
       B_m(3,1+(i-1)*5) = dNxy(1,i)

       ! bending
       B_b(1,4+(i-1)*5) = dNxy(1,i)
       B_b(2,5+(i-1)*5) = dNxy(2,i)
       B_b(3,4+(i-1)*5) = dNxy(2,i)
       B_b(3,5+(i-1)*5) = dNxy(1,i)

       ! shear
       B_s(1,3+(i-1)*5) = -dNxy(1,i)
       B_s(1,4+(i-1)*5) = Nvec(i)
       B_s(2,3+(i-1)*5) = -dNxy(2,i)
       B_s(2,5+(i-1)*5) = Nvec(i)
    end do

  END SUBROUTINE shape

END MODULE laminate45


