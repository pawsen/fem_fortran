module nonlin

  implicit none

  private
  public ::  non_lin, cauchy, residual ! non-linear routines

contains


   ! Ikke-line�re rutiner

  subroutine non_lin(flag,rho,rho_min)

    ! This subroutine calculates displacements using non-linear teory

    use numeth
    use fea
    use processor
    use plot_routiner
    use fedata
    use file_init

	integer, intent(IN) :: flag
    real(8), dimension(:), optional, intent(IN) :: rho
    real(8), optional, intent(in) :: rho_min
    integer :: n,i,e,j,kk,dev, n_incr,n_inner , i1, i2,a, animation,idof
    real(8), dimension(neqn)::deltaP, deltaD, R, PP
    real(8):: normR,NormP, tol, normK
    real(8) :: plotval(ne)
    real(8) :: b

    ! data for NR algorithm
    tol=1e-6
    n_incr=5
    n_inner = 100
    animation = 0 !	0: ingen animation

	! indl�s datafil
    !subroutine inputfile(file,filter_type,solver_type,problem_type,vol_type,rmin_type,info,save_rho,n_iter,&
	!				stop_krit,animation,penal,damp_fact,max_vol,rmin,rho_min,tol,movelimit)
    !if ((present(rho) .and. (.not. present(rho_min)) ) then
    !if ((present(rho)) then
   	!	call inputfile(a,a,a,a,a,a,a,a,n_incr,a,animation,penal,b,b,b,rho_min,tol,b)
    !end if

    print*,'n_incr',n_incr
    allocate( t_elem(ne) )
    t_elem = 0.000100d0 ! temp-stigning
	if (present(rho)) then
       antype = 'TOPTHERMSTRUCT_hard'
       call buildload(rho,rho_min)
    else
       call buildload
    end if

    do i = 1, nk
       ! Add spring contribution
       if (springs(i, 5) == 1) then ! input/output => 0/1
          idof = 2*(springs(i, 2)-1)+springs(i, 3)
          ! sat minus for at maksimere forskydning
       end if
    end do


    deltaP = P/real(n_incr) ! load incr.
	PP = P
    P = 0.0! reset load vector
    d = 0.0 !reset displacement

	kk = 0
	n = 0 ! t�llevariabel til animation
    do i=1,n_incr
       P = P+deltaP

       !$$$$$$     print*,'i',i

       ! Start of load increment loop   
       do j=1,n_inner
          n = n+1

          if (present(rho)) then
             ! calculate K and factorize. Rho skal med ind her
             call buildstiff_nonlin(flag,R,rho,rho_min)
          else
             call buildstiff_nonlin(flag,R)
          end if

          ! calculate residual
          R=R-P

          ! break criterion
          normR = dsqrt(DOT_PRODUCT(R,R))
          normP = dsqrt(DOT_PRODUCT(PP,PP))

          !$$$$$$         if (present(rho)) then
          !$$$$$$             ! calculate K and factorize. Rho skal med ind her
          !$$$$$$             call buildstiff(flag,rho,rho_min)
          !$$$$$$         else
          !$$$$$$             call buildstiff(flag)
          !$$$$$$         end if

          !$$$$$$         normK = 0.0
          !$$$$$$         do i1=1,bw
          !$$$$$$             do i2 =1,neqn
          !$$$$$$                 normK = normK + k(i1,i2)**2
          !$$$$$$             end do
          !$$$$$$         end do
          !$$$$$$         normK = dsqrt(normK)
          !$$$$$$         print*,'normK', normK

          call enforce_fea
          k=-1*k
          call bfactor(k)

          call bsolve(k, R)
          deltaD(1:neqn) = R(1:neqn)
          d=d+deltaD

          ! break skal f�rst ske efter k er blevet faktoriseret, 
          if (normR/normP<tol) then
             !$$$$$$           print*,j 
             exit
          end if
          if (j == n_inner) then
             print*,'residual: ',normR/normP
          end if


          ! update displacements
          d=d+deltaD

          select case(animation)
          case(1)
             call  plotanim(n, 0, 1, .true., .true., .false., .true.,.true., 'Udboej', &
                  0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1d0,rho)
          case(2)! udskriv til fil
             dev = 5
             if (mod(n,dev) == 0) then ! if the remainder after deviding by 10 is zero, then call output. Eg. output kaldes ved hvert 10'ende n.
                kk = kk+1
                if (kk == 1) then; call write_structure('_plotdeformed_topology.m.m') 
                end if
                call output_anim(D,kk) !,plotval) udskriver forskydning for hvert skridt
             end if
          end select
       end do
       print*,'forskydning',D(idof)

    end do
    ! Luk animationsvonduet ved at kalde med iter = -1
    select case( animation )
    case(1)
       call  plotanim(-1, 0, 1, .true., .true., .false., .true.,.true., 'Udboej', &
             0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1d0,rho)
    end select

    call cauchy ! does the same af recover - Just cauchy stresses instead.
    !call buildload

    !$$$$$$ compliance=dot_product(d,p)
    do e=1,ne
       plotval(e) = (stress(e,1)**2+stress(e,2)**2-stress(e,1)*stress(e,2)+3.0*stress(e,3)**2)**(0.50d0) !von mises sp�nding
    end do
    !$$$$$$     
    !$$$$$$      print*,compliance
    !$$$$$$     call output  
    !$$$$$$ 
    !$$$$$$     ! Plot deformed shape
    !$$$$$$     call output_deformed('deformeret',rho)
    !$$$$$$    call plot('un/deformed', 'xwin', 'color', ' ')
    !$$$$$$    call plot('elements', 'matlab', 'color', 'spaendingsintensitet', plotval)

  end subroutine non_lin

  subroutine buildstiff_nonlin(flag,outvector,rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata
    use plane42nonlin

	integer, INTENT(IN) :: flag
	real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), dimension(neqn), INTENT(OUT) :: outvector
   	real(8), optional, INTENT(IN) :: rho_min

    integer :: e, i, j
    integer :: nen,idof
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe,residual
    real(8), dimension(mdim, mdim) :: ke
    real(8) :: young, nu, area, dens, thk
    bw = size(k,1)  
    ! Reset stiffness matrix
    outvector = 0d0
    if (banded /=2) then
       K = 0d0
    end if

    do e = 1, ne

       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
       end do

       ! da ke laves p� baggrund af forskydningen for det enkelte element, er det ikke muligt at bruge "flag" option
       young = mprop(element(e)%mat)%young
       thk = mprop(element(e)%mat)%thk
       nu = mprop(element(e)%mat)%nu
       call plane42nonlin_ke(xe,d(edof), young, nu, thk,ng, ke,residual)

       do i = 1, 2*nen
          do j =1, 2*nen
             if (edof(i)>=edof(j)) then ! brug kun elementer der st�r under eller p� diagonalen i global K+
                !if ((antype == 'STATIC') .or. (antype == 'COUPLED') .or. (antype == 'TRANSIENT')) then
                if (.not. present(rho)) then
                   k((edof(i)-edof(j)+1),edof(j)) = k((edof(i)-edof(j)+1),edof(j)) + ke(i,j)
                else
                   k((edof(i)-edof(j)+1),edof(j)) = k((edof(i)-edof(j)+1),edof(j)) + &
                        (rho_min+(1d0-rho_min)*rho(e)**penal)*ke(i,j)! Modificeret SIMP-v�gtning af elasticitetsmodul
                end if
             end if
          end do
       end do

       if (present(rho)) then
          outvector(edof) = outvector(edof)+(rho_min+(1d0-rho_min)*rho(e)**penal)*residual 
       else
          outvector(edof) = outvector(edof)+residual    
       end if

    end do

    do i = 1, nb
       idof = 2*(bound(i,1)-1) + bound(i,2)
       outvector(idof) = 0d0 ! Correct for BC
    end do


    !�������������������������������������������������������������
    ! springs
    if (nk > 0) then
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 1) == 1) then
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             k(1,idof) = k(1,idof) + springs(i, 4)
          end if
       end do
    end if
    !�������������������������������������������������������������

  end subroutine buildstiff_nonlin



  subroutine residual(outvector,rho,rho_min)
    ! This subroutine calculates the outvector(residual)

    use fedata
    use plane42nonlin

    real(8), dimension(neqn), INTENT(OUT) :: outvector
	real(8), optional, dimension(:), INTENT(IN) :: rho
   	real(8), optional, INTENT(IN) :: rho_min

    integer :: e, i ,idof
    integer :: nen
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe,de, residualet
    real(8) :: thk, young, nu

	outvector = 0.0

    do e = 1, ne
       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
       end do

       thk = mprop(element(e)%mat)%thk
       young = mprop(element(e)%mat)%young
       nu = mprop(element(e)%mat)%nu

       call plane42nonlin_residual(xe,d(edof), young,nu,thk,ng,residualet)
       if (present(rho)) then
          outvector(edof) = outvector(edof)+(rho_min+(1d0-rho_min)*rho(e)**penal)*residualet 
       else
          outvector(edof) = outvector(edof)+residualet    
       end if

    end do

    do i = 1, nb
       idof = 2*(bound(i,1)-1) + bound(i,2)
       outvector(idof) = 0 ! Correct for BC
    end do


  end subroutine residual


  subroutine cauchy
	use fedata
    use plane42nonlin

    integer :: e,i, nen
    integer, parameter :: mdim = 8
    integer :: edof(mdim)
    real(8), dimension(mdim) :: xe, de
    real(8) :: young, nu, thk, dens, xi, eta
    REAL(8) :: Bmat(3, 8), B0mat(3, 8),BLmat(3, 8), Nmat(2,8),Gmat(4,8), jac(2,2), detjac !Output from plane42nonlin_shape
    real(8) :: helpproduct(2,2), THETA(4), F(2,2), cauchy_stress(2,2), detF, PKstress(2,2), Identity(2,2), DD(2,2)
    real(8), dimension(3) :: estrain, estress

	! Calculate stresses i centroid
	xi=0.0
    eta=0.0

	Identity=0.0d0 ! Identity matrix
    Identity(1,1)=1.0d0
    Identity(2,2)=1.0d0

  	do e = 1, ne
      ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i)   = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
          de(2*i-1) = d(edof(2*i-1))
          de(2*i)   = d(edof(2*i))
       end do

       young = mprop(element(e)%mat)%young
       nu = mprop(element(e)%mat)%nu

       call plane42nonlin_shape(xe,de,xi,eta, Nmat,Bmat,B0mat,BLmat,Gmat, jac, detjac)
       call plane42nonlin_ss(xe, de,young,nu,B0mat, BLmat, estress, estrain)

       THETA=MATMUL(Gmat,de) !(5) i ekstra noter
       DD(1:2,1)=THETA(1:2) !(3.10) krenk noter
       DD(1:2,2)=THETA(3:4)

       F = Identity+DD!(3.11) Krenk noter
       detF=F(1,1)*F(2,2)-F(2,1)*F(1,2)

       ! Piola-Kirchhoff stress tensor
       PKstress(1,1)=estress(1)
       PKstress(1,2)=estress(3)
       PKstress(2,1)=estress(3)
       PKstress(2,2)=estress(2)

       helpproduct = MATMUL(F,PKstress)
       cauchy_stress=1/detF*MATMUL(helpproduct,transpose(F))
       stress(e, 1) = cauchy_stress(1,1)
       stress(e, 2) = cauchy_stress(2,2)
       stress(e, 3) = cauchy_stress(2,1)
       !$$$$$$         plotval(e) = (stress(e,1)**2+stress(e,2)**2-stress(e,1)*stress(e,2)+3.0*stress(e,3)**2)**(0.50d0) !von mises sp�nding

       strain(e, 1:3) = estrain
	end do

  end subroutine cauchy

  !�������������������������������������������������

end module nonlin
