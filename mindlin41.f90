MODULE mindlin42

  ! This module contains subroutines specific to the MINDLIN41 element.

  use numeth
  IMPLICIT NONE

  PRIVATE :: mindlin42_shape, constitutive_matrix,piezo_matrix,dielectric_matrix
  PRIVATE :: isotrop_constitutive_matrix

  PUBLIC :: mindlin42_ke, mindlin42_me, mindlin42_re, mindlin42_ss
  PUBLIC :: mindlin42_piezo, mindlin42_dielectric

CONTAINS

  SUBROUTINE mindlin42_dielectric(xe, mat_vec,ep,thk,ng, ke,type,phi)

    ! This subroutine constructs the element conductivity matrix [k_e].
    use plane41, only: plane41_shape

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: type,ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk, ep
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(4,4), INTENT(OUT) :: ke
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j
    REAL(8) :: Jmat(2,2), helpmat(4,2) !epsilon(2,2)
    REAL(8) :: Bmat(2,4),  Nvec(4), jac(2,2), detjac
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp(2,2) ! rotation

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    call dielectric_matrix(type,ep,thk,mat_vec,Jmat)

    if (present(phi)) then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp = MATMUL(TRANSPOSE(Lambda),Jmat)
       Jmat = MATMUL(helpp,Lambda)
    end if

    ke = 0.0d0
    do i = 1 , ng
       do j = 1 , ng
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat, detjac, jac)
          helpmat = MATMUL(TRANSPOSE(bmat),Jmat)
          ke = ke - w(i)*w(j) * MATMUL(helpmat, bmat) * detjac ! minus is intended!
       end do
    end do

  END SUBROUTINE mindlin42_dielectric

  SUBROUTINE mindlin42_piezo(xe, mat_vec, thk,ng, ke,phi)
    ! This subroutine constructs the stiffness matrix for
    ! a rectangular 4-noded quad element.
    use plane41, only : plane41_shape

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng ! number of gauss-points
    REAL(8), INTENT(IN) :: thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke!ke(12,4)
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j, ii, jj
    REAL(8) :: helpproduct(12,2), H1mat(2,3), H2mat(2,3)!Emat(2,3)
    real(8) :: detjac, jac(2,2), Nmat(3,12), Bb(3,12), Bs(2,12) ! output fra mindlin42_shape
    REAL(8) :: Bmat_t(2,4),  Nvec(4), jac_t(2,2), detjac_t ! output fra plane41_shape
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp(2,3) ! rotation
    !$$$$$$   REAL(8), DIMENSION(:,:) , allocatable:: helpproduct

    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    call piezo_matrix(mat_vec,thk,H1mat,H2mat)

    if (present(phi)) then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp = MATMUL(TRANSPOSE(Lambda),H2mat)
       H2mat = MATMUL(helpp,T_p)
    end if
!!$  

    ke=0d0
    do i= 1,ng !number of gauss points
       do j = 1,ng
          call mindlin42_shape(xe,xi(i),eta(j),Bb,Bs,jac,detjac,Nmat)
          call plane41_shape(xe, xi(i), eta(j), Nvec, bmat_t, detjac_t, jac_t)

          helpproduct = MATMUL(TRANSPOSE(Bb),TRANSPOSE(H2mat))
          ke = ke + W(i)*W(j)*MATMUL(helpproduct,bmat_t)*detjac
       end do
    end do

  END SUBROUTINE mindlin42_piezo


  SUBROUTINE mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke,type,mat_vec,phi)


    ! This subroutine constructs the stiffness matrix for
    ! a quadrilateral 4-noded element.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ng,type
    REAL(8), INTENT(IN) :: young, nu, dens, thk
    real(8), intent(INOUT) :: shear
    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), DIMENSION(:), OPTIONAL, INTENT(IN) ::  mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke
    REAL(8),optional, INTENT(in) :: phi

    INTEGER :: i, j
    REAL(8) :: ke_b(12,12), ke_s(12,12),detjac, xi_1, &
         eta_1, w_1
    REAL(8) :: jac(2,2), Nmat(3,12), k, B_b(3,12), B_s(2,12)
    real(8), dimension(3,3) :: Dmat, Amat, Bmat
    real(8) :: Ds_mat(2,2), real_parameters(5)
    REAL(8) :: T_p(3,3), T_t(2,2), Lambda(2,2), helpp1(3,3), helpp2(2,2) ! rotation
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    ! værdien af faktoren k, står på s 535 i cook
    k = 5d0/6d0
    !k = 5d0/8d0

    real_parameters(1) = nu
    real_parameters(2) = young
    real_parameters(3) = shear
    real_parameters(4) = thk
    real_parameters(5) = k

    call constitutive_matrix(type,real_parameters, &
       mat_vec,Amat,Bmat,Dmat,Ds_mat)

    if ((present(phi)) .and. (type == 2) )then
       call rotation_matrix(phi,T_p,T_t,lambda)
       helpp1 = MATMUL(TRANSPOSE(T_p),Dmat)
       Dmat = MATMUL(helpp1,T_p)
       
       helpp2 = MATMUL(TRANSPOSE(T_t),Ds_mat)
       Ds_mat = MATMUL(helpp2,T_t)      
    end if


    ke = 0d0
    ke_b = 0d0
    ke_s = 0d0
    do i = 1, ng
       do j = 1, ng
          call mindlin42_shape(xe,xi(i),eta(j),B_b,B_s,jac,detjac,Nmat)
          ke_b = ke_b + w(i)*w(j)*MATMUL(MATMUL(TRANSPOSE(B_b),Dmat),B_b)*detjac
       end do
    end do

    ! For at undgå shear-locking, skal k_s integreres med "single Gauss kvadrature"(dvs ng = 1), jf
    ! "Element Behavior", s 543 i cook
    xi_1 = 0d0
    eta_1 = 0d0
    w_1 = 2d0
    call mindlin42_shape(xe,xi_1,eta_1,B_b,B_s,jac,detjac,Nmat)
    ke_s = w_1*w_1*MATMUL(MATMUL(TRANSPOSE(B_s),Ds_mat),B_s)*detjac

    ke = ke_b + ke_s

  END SUBROUTINE MINDLIN42_ke


  subroutine constitutive_matrix(type,real_parameters, &
       mat_vec,Amat,Bmat,Dmat,Ds_mat)
    ! This sunroutine calculates the element constitutive maxtrix for laminates
    
    real(8), intent(OUT), DIMENSION(3,3) :: Amat,Bmat,Dmat
    real(8), intent(OUT) :: Ds_mat(2,2)
    real(8), intent(in) :: real_parameters(:), mat_vec(:)
    integer, intent(in) :: type
    
    real(8) :: Cmat(3,3), Cs_mat(2,2), thk2,thk_ptz,thk_pcb,nu,young,shear,k, thk
    real(8) :: Cmat_pcb(3,3), Cs_mat_pcb(2,2)
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i
    real(8) :: Dmat_test(3,3), Ds_mat_test(2,2)
    ! nc: number of layers, z: layer-coordinate
    thk = real_parameters(4)
    k = real_parameters(5)
    
    ! nc: number of layers, z: layer-coordinate
     nc = int(mat_vec(13))
    if (nc> 2) then
       allocate(z(nc+1))
       call calculate_z(mat_vec,z)
       thk_pcb = z(3)-z(2) ! for 3-layer laminate, hvor det midterste er pcb
    else 
       nc = 1
       allocate(z(nc+1))
       z(1) = -0.5d0*thk
       z(2) = 0.5d0*thk
       thk_pcb = thk
    end if
    Amat = 0d0
    Bmat = 0d0
    Dmat = 0d0
    Ds_mat = 0d0
  
    select case(type)
    case(1)! normalt isotropisk materiale
       call isotrop_constitutive_matrix(real_parameters,Cmat_pcb,Cs_mat_pcb)
       Dmat = thk_pcb**3/12d0*Cmat_pcb
       Ds_mat = Cs_mat_pcb*thk_pcb
    case(2) ! piezoelectric material

       Cmat(1, 1) = mat_vec(1) - mat_vec(2)**2/mat_vec(1)
       Cmat(1, 2) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
       Cmat(2, 1) = mat_vec(3) - mat_vec(2)*mat_vec(3)/mat_vec(1)
       Cmat(2, 2) = mat_vec(4) - mat_vec(3)**2/mat_vec(1)
       Cmat(3, 3) = mat_vec(5)

       Cs_mat(1,1) = mat_vec(6)*k!det er korrekt at mat_vec(6) kommer før mat_vec(5)
       Cs_mat(2,2) = mat_vec(5)*k

       if (nc == 3) then
          thk_ptz = z(4)-z(3)
          thk_pcb = z(3)-z(2) ! thk pcb
          call isotrop_constitutive_matrix(real_parameters,Cmat_pcb,Cs_mat_pcb)
          Dmat = 2d0/3d0*Cmat*(z(4)**3 - z(3)**3) + 2d0/3d0*Cmat_pcb*(z(3)**3 - 0d0)
          Ds_mat = 2*Cs_mat*thk_ptz + Cs_mat_pcb*thk_pcb
       else
          Dmat = thk**3/12d0*Cmat
          Ds_mat = Cs_mat*thk
       endif

    end select

    Dmat_test = 0d0
    Ds_mat_test = 0d0
    do i=1,nc
       thk2 = z(i+1) - z(i)
       ! z_bar = (1d0/2d0 * (z(i+1) - z(i)) ) ! = position of the middle of layer i.
       Amat = Amat + Cmat*thk2
       Bmat = Bmat + Cmat*thk2* (1d0/2d0 * (z(i+1) - z(i)) )
       Dmat_test = Dmat_test + Cmat*1d0/3d0* (z(i+1)**3 - z(i)**3)
       Ds_mat_test = Ds_mat_test + Cs_mat*thk2
    end do

  end subroutine constitutive_matrix
  

  subroutine isotrop_constitutive_matrix(real_parameters,Cmat,Cs_mat)
    real(8), intent(OUT) :: Cmat(3,3), Cs_mat(2,2)
    real(8), intent(in) :: real_parameters(:)
    real(8) :: nu, young, shear, k

    nu = real_parameters(1)
    young = real_parameters(2)
    shear = real_parameters(3)
    k = real_parameters(5)

    Cmat = 0d0
    Cs_mat = 0d0

    ! Build constitutive matrix (plane stress)
    Cmat(1,1) = 1d0
    Cmat(1,2) = nu
    Cmat(2,1) = nu
    Cmat(2,2) = 1d0
    Cmat(3,3) = (1d0-nu)/2d0
    Cmat = young/(1d0-nu**2)*Cmat

    if (shear < 1E-3 ) then! in case it's forgotten to give shear a value in input-file
       shear = young/(2*(1+nu))
       !print*,'shear er ikke angivet i inputfil. Værdi er udregnet, mindlin41.f90'
    end if

    ! shear stiffness
    Cs_mat(1,1) = k*shear
    Cs_mat(2,2) = k*shear

  end subroutine isotrop_constitutive_matrix

  subroutine calculate_z(mat_vec,z)
    real(8), INTENT(in) :: mat_vec(:)
    real(8), INTENT(out) :: z(:)
    integer :: nc

    nc= int(mat_vec(13))

    select case(nc)
    case(2)
       z(1) = -mat_vec(14)
       z(2) = 0d0
       z(3) = mat_vec(15)
    case(3)
       z(1) = -mat_vec(14) - 0.5d0*mat_vec(15)
       z(2) = - 0.5d0*mat_vec(15)
       z(3) = + 0.5d0*mat_vec(15)
       z(4) = mat_vec(14) + 0.5d0*mat_vec(15)
    end select

  end subroutine calculate_z

  subroutine piezo_matrix(mat_vec,thk,H1mat,H2mat)
    ! This sunroutine calculates the element electric  constitutive maxtrix for laminates

    real(8), intent(OUT) :: H1mat(2,3), H2mat(2,3)
    real(8), intent(in) :: mat_vec(:), thk
    
    real(8) :: Emat(2,3), thk1, z_mean_ptz, thk_ptz, thk_pcb
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i
    real(8):: H1mat_test(2,3), H2mat_test(2,3)

    ! nc: number of layers, z: layer-coordinate
    nc = int(mat_vec(13))
    allocate(z(nc+1))

    call calculate_z(mat_vec,z)

    H1mat = 0d0
    H2mat = 0d0

    Emat = 0d0
    Emat(2,1) = mat_vec(7) - mat_vec(2)*mat_vec(7)/mat_vec(1)
    Emat(2,2) = mat_vec(8) - mat_vec(3)*mat_vec(7)/mat_vec(1)
    Emat(1,3) = mat_vec(9)

    select case(nc)
    case(1:2)
       ! integration of two laminate plates, each with thickness 1/2*thk
       ! Emat ==  H2mat in this case
       H2mat = 1d0/4d0 * thk**2 * Emat
    case(3)!
       thk_pcb = z(3)-z(2)
       thk_ptz = z(4)-z(3)!thk_ptz
       z_mean_ptz = 0.5d0*thk_ptz+0.5d0*thk_pcb
       H2mat = 2*Emat*thk_ptz*z_mean_ptz
    end select

    do i=1,nc
       thk1 = z(i+1) - z(i)

       H1mat_test = H1mat_test + Emat*thk1
       H2mat_test = H2mat_test + Emat*thk1*(1d0/2d0 * (z(i+1) - z(i)) )
    end do


  end subroutine piezo_matrix

  subroutine dielectric_matrix(type,ep,thk,mat_vec,Jmat)
    ! This sunroutine calculates the element electric  constitutive maxtrix for laminates

    real(8), intent(OUT) :: Jmat(2,2)
    real(8), intent(in) :: ep,thk, mat_vec(:)
    integer, intent(in) :: type

    real(8) :: thk_pcb,thk_ptz, epsilon(2,2), epsilon_pcb(2,2)
    real(8), allocatable, DIMENSION(:) :: z
    integer :: nc,i
    ! nc: number of layers, z: layer-coordinate
    nc = int(mat_vec(13))
    allocate(z(nc+1))
    call calculate_z(mat_vec,z)

    Jmat = 0d0
    epsilon = 0d0

    select case(type)
       case(1) ! normalt isotropisk materiale
          epsilon_pcb(1,1) = ep
          epsilon_pcb(2,2) = ep
          if (nc == 1) then
             thk_pcb = thk
          else
             thk_pcb = (z(3)-z(2))
          end if
          Jmat = epsilon_pcb*thk_pcb
       case(2)! piezo-material
          epsilon(1,1) = mat_vec(10)
          epsilon(2,2) = mat_vec(11)+mat_vec(7)**2/mat_vec(1)

          if (nc == 3) then ! include dielectric from pcb-middle-layer
             epsilon_pcb(1,1) = ep
             epsilon_pcb(2,2) = ep

             thk_ptz = z(4)-z(3)
             thk_pcb = z(3)-z(2)
             Jmat = 2d0*epsilon*thk_ptz + epsilon_pcb* thk_pcb

             !print *,ep, mat_vec(10),  mat_vec(11)+mat_vec(7)**2/mat_vec(1)
          else
             ! integration of two laminate plates, each with thickness 1/2*thk
             Jmat = epsilon*thk
          end if
          
       end select

  end subroutine dielectric_matrix


  SUBROUTINE mindlin42_me(xe,dens_pcb,mat_vec,thk,ng,me,type)

    REAL(8), INTENT(IN) :: dens_pcb, thk
    REAL(8), DIMENSION(:), INTENT(IN) :: xe, mat_vec
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: me
    integer, INTENT(IN) :: ng, type
    
    INTEGER :: i, j, nc
    real(8) :: detjac, jac(2,2), Nmat(3,12), Bb(3,12), Bs(2,12) ! output fra mindlin42_shape
    real(8) :: dens1,dens2, dens_ptz,z_mean_ptz, thk_ptz, thk_pcb, rho(3,3)
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    real(8), allocatable, DIMENSION(:) :: z
    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    
    dens_ptz = mat_vec(12)
    ! nc: number of layers
    nc = int(mat_vec(13))
    allocate(z(nc+1))
    call calculate_z(mat_vec,z)

    if (nc == 3) then! for 3-layer laminate, hvor det midterste er pcb
       thk_pcb = mat_vec(15)
       thk_ptz = mat_vec(14)
    else 
       thk_pcb = thk
       thk_ptz = thk
    end if
    
    select case(type)
    case(1) ! normalt 
       dens1 = thk_pcb*dens_pcb
       dens2 = dens_pcb*thk_pcb**3 / 12d0
    case(2)
       select case(nc)
       case(3)
          z_mean_ptz = 0.5d0*thk_ptz+0.5d0*thk_pcb

          dens1 = thk_pcb*dens_pcb + 2*thk_ptz*dens_ptz
          dens2 =2d0/3d0* dens_ptz *(z(4)**3 - z(3)**3) + &
               2d0/3d0* dens_pcb *(z(3)**3 - 0d0)

!!$          dens1 = thk_pcb*dens_pcb + 2*thk_ptz*dens_ptz
!!$          dens2 = 2*(thk_ptz**3/12d0*dens_ptz + z_mean_ptz**2*dens_ptz) + &
!!$               thk_pcb**3/12d0*dens_pcb
       case default
          dens1 = dens_ptz*thk_ptz
          dens2 = dens_ptz*thk_ptz**3 / 12d0
       end select
    end select

    rho = 0d0
    rho(1,1) = dens1
    rho(2,2) = dens2
    rho(3,3) = dens2

    me =  0d0
	do i= 1,ng !number of gauss points
       do j = 1,ng
          call mindlin42_shape(xe,xi(i),eta(j),Bb,Bs,jac,detjac,Nmat)
          me =me+ W(i)*W(j)*MATMUL(MATMUL(TRANSPOSE(Nmat),rho) ,Nmat)*detjac!thk er inkluderet i dens
       end do
    end do

  END SUBROUTINE MINDLIN42_me

  SUBROUTINE MINDLIN42_re(xe, eface, fe, thk, re)

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

    aa = (xe(3)-xe(1))/2d0
    bb = (xe(8)-xe(2))/2d0
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
       ! Ved mindlin-plader skal en fordelt latteral kraft kun ligges til i Uz-frihederne, dvs intet moment-bidrag, modsat kirchoff hvor der er kobling

       !Forces in z-dir.
       re(1) = aa*bb*fe !Node 1
       re(4) = aa*bb*fe !Node 2
       re(7) = aa*bb*fe !Node 3
       re(10)= aa*bb*fe !Node 4
    ELSE
       write (*, *) 'Error: Wrong face number!!'
       stop

    ENDIF

!!!!
!!!!! GANG MED TYKKELSE!!!
!!!!!
    re = re*thk

  END SUBROUTINE MINDLIN42_re

  SUBROUTINE MINDLIN42_ss(xe, de, young, nu, thk, shear, estress, estrain)

    ! This subrotuine constructs the element stress and strain

    REAL(8), INTENT(IN) :: young, nu, thk, shear! G = shear
    REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
    REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain

    REAL(8) :: Cmat(3, 3) ,k, Dmat(3,3), D_M(5,5), kappa(5), MCR(5),xi,eta
    REAL(8) :: jac(2,2), N(3,12), detjac,B_M(5,12), B_b(3,12), B_s(2,12)
    REAL(8) :: estress_2(3)

    ! værdien af faktoren k, står på s 535 i cook
    k = 5d0/6d0

    xi = 0
    eta = 0
    call mindlin42_shape(xe,xi,eta,B_b,B_s,jac,detjac,N)

    B_M(1:3,:) = B_b
    B_M(4:5,:) = B_s

    ! Compute element strain
    ! kappa er givet som, jf (15.3-3), s. 542
    kappa = MATMUL(B_M, de)

    ! strain er givet ved (15.1-1), s 532 i cook
    estrain = -kappa
    estrain(1:3) = estrain(1:3)*thk/2


    ! Build constitutive matrix (plane stress)
    Cmat = 0.0d0
    Cmat(1,1) = 1d0
    Cmat(1,2) = nu
    Cmat(2,1) = nu
    Cmat(2,2) = 1d0
    Cmat(3,3) = (1.d0-nu)/2

    Cmat = young/(1.d0-nu*nu)*Cmat

    Dmat = thk**3/12*Cmat
    D_M = 0.d0
    D_M(1:3,1:3) = Dmat
    D_M(4,4) = k*shear*thk
    D_M(5,5) = k*shear*thk

    !Momenter og kræfter(en kraft, flere kræfter!), jf (15.1-5), s 534 i cook
    MCR = -MATMUL(D_M, kappa)
    ! spændinger angives i teksten øverst s 533 i cook
    estress(1:3) = 6*MCR(1:3)/thk**2
    estress(4:5) = 1.5*MCR(4:5)/thk

    !alternativt kan spændingerne gives ud fra plane-stress udtrykket, 15.1-3
    estress_2 =   MATMUL(Cmat, estrain(1:3))


  END SUBROUTINE MINDLIN42_ss

  SUBROUTINE mindlin42_shape(xe,xi,eta,B_b,B_s,jac,detjac,N)

    ! This subroutine calculates the shape function matrix, the Jacobian
    ! matrix and the determinant of the Jacobian matrix.

    IMPLICIT NONE

    REAL(8), DIMENSION(:), INTENT(IN) :: xe
    REAL(8), INTENT(IN) :: xi, eta
    REAL(8), DIMENSION(:,:), INTENT(OUT) ::  B_b(3,12),B_s(2,12) , jac(2,2), N
    REAL(8), INTENT(OUT) :: detjac

    REAL(8) :: N1, N2, N3, N4, Nvec(4)
    REAL(8), DIMENSION(3,4) :: L
    REAL(8), DIMENSION(2,4) :: dN, dNxy
    REAL(8), DIMENSION(2,2) :: Gamma
    REAL(8), DIMENSION(4,4) :: Gtilde
    INTEGER :: i


    ! Calculate shape functions
    N1 = 1d0/4d0*(1d0-xi)*(1d0-eta)
    N2 = 1d0/4d0*(1d0+xi)*(1d0-eta)
    N3 = 1d0/4d0*(1d0+xi)*(1d0+eta)
    N4 = 1d0/4d0*(1d0-xi)*(1d0+eta)

    !N_i,eta
    dN(1,1) = -(1d0-eta)
    dN(1,2) = (1d0-eta)
    dN(1,3) = (1d0+eta)
    dN(1,4) = -(1d0+eta)
    dN(2,1) = -(1d0-xi)
    dN(2,2) = -(1d0+xi)
    dN(2,3) = (1d0+xi)
    dN(2,4) = (1d0-xi)
    dN = dN*1d0/4d0

    Nvec = [N1,N2,N3,N4]
    
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
       jac(1,1) = jac(1,1) + dN(1,i)*xe(2*i-1)! 1,3,5...
       jac(1,2) = jac(1,2) + dN(1,i)*xe(2*i)! 2,4,6
       jac(2,1) = jac(2,1) + dN(2,i)*xe(2*i-1)
       jac(2,2) = jac(2,2) + dN(2,i)*xe(2*i)
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

    ! B splittes op i bending og shear komponent
    
    B_b = 0d0
    B_s = 0d0
    do i=1,4
       ! bending
       B_b(1,2+(i-1)*3) = dNxy(1,i)
       B_b(2,3+(i-1)*3) = dNxy(2,i)
       B_b(3,2+(i-1)*3) = dNxy(2,i)
       B_b(3,3+(i-1)*3) = dNxy(1,i)

       ! shear
       B_s(1,1+(i-1)*3) = -dNxy(1,i)
       B_s(1,2+(i-1)*3) = Nvec(i)
       B_s(2,1+(i-1)*3) = -dNxy(2,i)
       B_s(2,3+(i-1)*3) = Nvec(i)
    end do


  END SUBROUTINE mindlin42_shape

END MODULE MINDLIN42
