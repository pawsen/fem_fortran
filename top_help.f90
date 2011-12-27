MODULE top_help

  ! This module contains subroutines specific to the PLANE42 element.

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mma_handle, oc, neighbourMatrix,filter, density_filter, gradient, &
       finite_check, lambda_calc, force_rho,force_elements_init

CONTAINS

  SUBROUTINE mma_handle(iteration,inak,low,upp,max_vol,objekt,vol,rho, rho_old, rho_old2, &
       df0dx, dg,rho_tilde,nc,movelimit_in,constraint,dconstraint)
    
	use fedata

	REAL(8), DIMENSION(:), intent(INOUT) :: low, upp, rho_old, rho_old2, rho
    INTEGER, INTENT(IN) :: iteration, nc! nc = number of constraint
    REAL(8), DIMENSION(:),optional, intent(IN) :: rho_tilde
    REAL(8), DIMENSION(:), intent(IN) :: df0dx,vol,dg

    REAL(8), DIMENSION(:,:), optional, intent(IN) :: dconstraint ! gradienter for constraint (udover volumen)

    real(8), INTENT(IN) :: objekt,max_vol
    real(8),optional, INTENT(INOUT) :: movelimit_in, constraint(:)
    type(INAKTIV_def), intent(in) :: inak

    INTEGER:: e, j, ne_aktiv

    ! MMA 
   	! Initialisering af hjælpevariable
   	INTEGER,dimension(:), allocatable :: IYFREE
    REAL(8),dimension(:), allocatable :: ALFA, BETA, P0, Q0 ! size=ne
	! Output
    REAL(8) :: zmma
    REAL(8),DIMENSION(:), allocatable :: rho_new ! size=ne
    ! Input
    REAL(8) :: geps, top_vol, movelimit
    REAL(8), DIMENSION(:), allocatable :: amma,cmma, fval, fmax ! st parameter
    REAL(8), DIMENSION(:), allocatable :: xmin, xmax !size=ne
	REAL(8), DIMENSION(:), allocatable :: dfdx
    REAL(8), DIMENSION(:), allocatable :: UU,GRADF, DSRCH,PMMA, QMMA,BMMA, ymma, ulam,HESSF

    if (inak%bool) then
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if
    j = ne_aktiv ! number of elements that is to be optimized
    allocate(xmin(j),xmax(j))
    allocate(ALFA(j),BETA(j),P0(j),Q0(j))
    allocate(rho_new(j))
	allocate( UU(nc),GRADF(nc), DSRCH(nc),BMMA(nc),IYFREE(nc),ymma(nc), ulam(nc),HESSF(nc*(nc+1)/2))
	allocate(amma(nc),cmma(nc),fval(nc),fmax(nc)) 
	allocate( dfdx(nc*j),PMMA(nc*j), QMMA(nc*j)  ) ! alle gradienter til contraints      

 ! Initialisering af MMA-parametre
    amma = 0.0
    cmma = 1000.0
    geps = 1.0E-6
    fmax = 0.0 ! rhs. af constraint, sum(rho)/V_max - 1 <=0
    top_vol = 0.0
    IYFREE = 0


    if (.not. present(movelimit_in)) then
       movelimit = 0.2d0
    else
       movelimit = movelimit_in
    end if

    if (antype == 'TRANSIENT') then
       ! Da der ikke er nogen begrænsning(constraint), skal fval sættes så den altid er opfyldt. Dvs fval < fmax =0
       fval = -1.0d0
    else
       top_vol = dot_product(rho_tilde,vol) ! det skal være den fysiske volumen og ikke matematisk
       fval(1) = top_vol/max_vol - 1 ! skaleret
       print *,'fval', fval(1)
    end if

	if (nc == 1) then ! kun volumen
       dfdx= dg
       amma = 0d0
    else
       amma = 1d0
       amma(nc) = 0D0
       fval(nc) = fval(1)
       do j=1,nc-1
          fval(j) = constraint(j)
       end do
       ! k = (j-1)*M + i
       do e = 1,ne_aktiv
          do j=1,nc
             dfdx((e-1)*nc+j) = dconstraint(e,j)
          end do
       end do
    end if

   	do e=1,ne_aktiv
       xmin(e) = DMAX1(0.001, rho(e) - movelimit) ! Set move limits on design variables
       xmax(e) = DMIN1(1.000, rho(e) + movelimit)
    end do

	! nc = number of constraint.
    call mmasub(iteration,nc,ne_aktiv,geps,IYFREE,rho,rho_new,xmin,xmax,rho_old,rho_old2,low,upp,&
         ALFA,BETA,amma,BMMA,cmma,ymma,zmma,ulam,objekt,fval,fmax,df0dx,dfdx,PMMA,QMMA,P0,Q0,UU, &
         GRADF,DSRCH,HESSF)

    if (iteration > 2) then
       rho_old2 = rho_old
    end if
    if (iteration > 1) then
       rho_old = rho
    end if
    rho = rho_new


  END SUBROUTINE mma_handle


  SUBROUTINE oc(filter_type,max_vol,vol,rho_new,rho_old,dc,dg)
	! Bisection, bestemmelse af lambda

    use fedata

    REAL(8), DIMENSION(:), INTENT(OUT) :: rho_new
    REAL(8), DIMENSION(:), INTENT(IN)  :: vol, rho_old, dc, dg
    integer, INTENT(IN) :: filter_type
    REAL(8), INTENT(IN) :: max_vol

    INTEGER :: e
    REAL(8) :: lambda1, lambda2, lambda_mid, epsilon, g, rho_min

    lambda1 = 1.0E-10 !nedre startgï¿½t
    lambda2 = 1.0E10 !ï¿½vre startgï¿½t
    !epsilon = 1.0E-8 ! stopkoefficient
    rho_min = 1.0E-6 ! Mindste vï¿½rdi densiteten kan antage

    do while ( (lambda2-lambda1)/(lambda1+lambda2) > 1.0E-8 )

       lambda_mid = (lambda1 + lambda2)/2.0;

       do e=1,ne
          rho_new(e) = rho_old(e) * ( -dc(e)/(lambda_mid*dg(e)) )**damp_fact
          if (rho_new(e)<= rho_min) then
             rho_new(e) = rho_min
          else if (rho_new(e) >= 1) then
             rho_new(e) = 1.0
          end if
       end do
       select case( filter_type )
       case( 3 )   ! For Density filter: The physical densities rho_hat should satisfy g.
          ! Computation of filteret densities rho_hat
          !call density_filter( vol, rho_new)
       end select
       g = 0.0
       g = DOT_PRODUCT(dg,rho_new) - max_vol ! Finder rod for "subject to" funktion, se (4.1) i noter
       if  (g > 0) then
          lambda1 = lambda_mid
       else
          lambda2 = lambda_mid
       end if
    end do

  END SUBROUTINE oc

  SUBROUTINE neighbourMatrix(inak,rmin_type,rmin)

    use fedata

    type(INAKTIV_def), intent(in) :: inak
    integer, INTENT(IN) :: rmin_type
    REAL(8), INTENT(INOUT) :: rmin

    integer :: e, f, i, nen, max_neigh, count, ii, jj, ne_aktiv
    real(8) :: dist_ei
    real(8), allocatable :: n_neigh(:), center(:,:)
    real(8) ::  xLength, yLength


    if (inak%bool) then
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if
    i = ne_aktiv
    allocate(n_neigh(i), center(i,2))

	max_neigh = 0
    center = 0.0d0

    do ii = 1, ne_aktiv

       if (inak%bool) then
          e = inak%act(ii)
       else
          e = ii
       end if

       nen = element(e)%numnode
       do i = 1,nen
          center(ii,1) =center(ii,1)+ 1.0/4.0 * x(element(e)%ix(i), 1) ! x-koordinat, Centrum af element e
          center(ii,2) =center(ii,2)+ 1.0/4.0 * x(element(e)%ix(i), 2)
       end do
   	end do

    select case ( rmin_type )
   	case( 2)! omregn r_min til abs- afstand
       xLength = real(x(element(1)%ix(2),1)) - real(x(element(1)%ix(1),1)) 
       yLength = real(x(element(1)%ix(2),2)) - real(x(element(1)%ix(1),2))
       if (yLength > xLength) then
          rmin = rmin*yLength
       else
          rmin = rmin*xLength
       end if
    end select

	n_neigh = 0.0d0 ! Antallet af elementer der ligger inden for rmin
 
   do e=1,ne_aktiv
       do i = 1,ne_aktiv

          dist_ei = dsqrt( (center(e,1)-center(i,1))**2 + (center(e,2)-center(i,2))**2 ) ! Afstand mellem centre af e og f
          if (dist_ei <= rmin .and. e /= i) then
             n_neigh(e) = n_neigh(e) + 1
          end if
       end do
       if (n_neigh(e) > max_neigh) then
          max_neigh = n_neigh(e)
       end if
    end do

    allocate( neigh(ne_aktiv,2*max_neigh+3) ) ! Matrix med naboskabs-info. Initialiseret i fedata
    neigh = 0.0d0

	print*,'Maks Naboelementer, excl elementet selv', max_neigh

    do e = 1, ne_aktiv
       count = 0
       neigh(e,1) = n_neigh(e)+1 ! antal af naboelementer - incl elementet selv
       neigh(e,2) = e !elementet selv tï¿½lles med
       neigh(e,3) = rmin ! det har vï¿½gtningen rmin
       do i = 1, ne_aktiv
          dist_ei = dsqrt( (center(e,1)-center(i,1))**2 + (center(e,2)-center(i,2))**2 )
          if (dist_ei <= rmin .and. e /= i .and. count+1 < n_neigh(e)) then
             count = count + 1
             neigh(e,2*(count+1)) = i ! Naboelementnummer
             neigh(e,2*(count+1)+1) = rmin - dist_ei ! vï¿½gtning
          elseif (dist_ei <= rmin .and. e /= i .and. count+1 == n_neigh(e)) then ! sidste element
             count = count + 1
             neigh(e,2*(count+1)) = i
             neigh(e,2*(count+1)+1) = rmin- dist_ei
             EXIT
          end if
       end do
	end do

 ! Selvom elementnumrene bliver gemt som real(dvs floating point) og senere konverteret til int(i routinerne hvor de bruges), synes det ikke at være et problem.

 ! print neighbourhood
 !$$$$$$              print*,'Weight ='
 !$$$$$$              do i = 1,10 ! elementer der skal printes for
 !$$$$$$                 write (*,'(10(f6.2,1x))') (neigh(i,e*2+1), e=1,max_neigh)
 !$$$$$$                 print*,'---------------------------------------------------------'
 !$$$$$$              end do
 !$$$$$$              print*,'Naboskab ='
 !$$$$$$              do i = 1,10
 !$$$$$$                write (*,'(10(f4.0,1x))') (neigh(i,e*2), e=1,max_neigh)
 !$$$$$$                print*,'----------------------------------------------------------'
 !$$$$$$              end do


  END SUBROUTINE neighbourMatrix

  subroutine filter(filter_type,inak,vol,rho,dc,dg,beta,eta)
    ! This subroutine select and apply filters:

    use fedata
    use plane42
    use fea

    integer, intent(IN) :: filter_type
    type(INAKTIV_def), intent(in) :: inak
    real(8), dimension(:), intent(IN) :: vol, rho
    real(8), dimension(:), intent(INOUT) :: dc, dg
    real(8), optional, intent(IN) :: beta,eta

    integer :: e, e2, e3, f, i, j, neigh_num, ne_aktiv
    real(8) :: numer_dc, numer_dg, dd, denum, weigth, numer, numer2
    real(8), dimension(:), allocatable :: dcn, dgn
    !    real(8), dimension(ne) :: dcn, dgn

    if (inak%bool) then
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if
    j = ne_aktiv
    allocate(dcn(j),dgn(j))

    select case( filter_type )
    case( 0 )

    case( 1 ) ! SENSITIVITY 
       ! This subroutine filters the sensitivities using a linearly decaying weighting function
       ! and accounts for non-regular meshes with varying element volumes, vol(e)

       dcn = 0.0d0
       do e = 1, ne_aktiv
          denum = 0.0d0
          do f = 1,int(Neigh(e,1)) !naboelementer incl elementet selv 
             neigh_num = Neigh(e,2*f) ! naboelementets nummer
             !print*,neigh_num
             if (neigh_num == 0) then
                exit
             else
                weigth = Neigh(e,2*f+1)
                dcn(e) = dcn(e) +  weigth* rho(neigh_num)*dc(neigh_num)/vol(neigh_num)
                denum = denum + weigth !vægtning
             end if
             !print*,weigth
          end do
          dcn(e) = dcn(e)/( rho(e) * denum / vol(e) ) 
          !$$$$$$         if (e == 5) then
          !$$$$$$           stop
          !$$$$$$         end if
       end do
       dc = dcn
       !print*,'max-element',maxval(Neigh(1:ne,1))

    case( 2 ) ! DENSITY
       dcn = 0.0d0
       dgn = 0.0d0

       do e = 1, ne_aktiv
          numer_dc = 0.0d0
          numer_dg = 0.0d0
          do e2 = 1,int( Neigh(e,1))
             dd = 0.0d0
             i = Neigh(e,2*e2)
             denum = 0.0d0
             do e3 = 1,int(Neigh(i,1))
                j = Neigh(i,2*e3)
                denum = denum + Neigh(i,2*e3+1)*vol(j)
             end do
             !$$$$$$             dd = Neigh(e,2*e2+1)*vol(i)/denum
             dd = Neigh(e,3)*vol(e)/denum
             dcn(e) = dcn(e) + dd*dc(i)
             dgn(e) = dgn(e) + dd*dg(i)
             !$$$$$$             numer_dc = numer_dc + dd*dc(i)
             !$$$$$$             numer_dg = numer_dg + dd*dg(i)
          end do
          !$$$$$$         dcn(e) = numer_dc
          !$$$$$$         dgn(e) = numer_dg
       end do
       dc = dcn
       dg = dgn

	case ( 3 )
       do e = 1, ne_aktiv
          numer = 0.0d0
          numer2 = 0.0d0
          do e2 = 1, int(Neigh(e,1))
           	 dd = 0.0d0
             i = Neigh(e,2*e2)
             denum = 0.0d0
             do e3 = 1, int(Neigh(i,1))
                j = Neigh(i,2*e3)
                denum = denum + Neigh(i,2*e3+1)*vol(j)
             end do
             dd = Neigh(e,2*e2+1)*vol(i)/denum

             dd = (beta*EXP(-beta*(1.0-rho(i)))+EXP(-beta))*dd ! kædereglen: rho er rho_tilde, hvor den kaldes  

             numer = numer + dd*dc(i)
             numer2 = numer2 + dd*dg(i)
          end do
          dcn(e) = numer
          dgn(e) = numer2
       end do
       dc = dcn
       dg = dgn

    case ( 4 )
       do e = 1, ne_aktiv
          numer = 0.0d0
          numer2 = 0.0d0
          do e2 = 1,int( Neigh(e,1))
           	 dd = 0.0d0
             i = Neigh(e,2*e2)
             denum = 0.0d0
             do e3 = 1, int(Neigh(i,1))
                j = Neigh(i,2*e3)
                denum = denum + Neigh(i,2*e3+1)*vol(j)
             end do
             dd = Neigh(e,2*e2+1)*vol(i)/denum

             dd = ((1-dtanh(beta*(rho(i)-eta))**2)*beta/&
                  (dtanh(beta*eta)+dtanh(beta*(1-eta)))) *dd ! kædereglen: rho er rho_tilde, hvor den kaldes  

             numer = numer + dd*dc(i)
             numer2 = numer2 + dd*dg(i)
          end do
          dcn(e) = numer
          dgn(e) = numer2
       end do
       dc = dcn
       dg = dgn
    end select



  end subroutine filter

 
  subroutine density_filter(inak,vol, rho)

    ! This subroutine apply density filter on rho to obtain rho_dens:
    use fedata

    !$$$$$$        integer, intent(IN) :: filter_type
    real(8), dimension(:), intent(IN) :: vol
    real(8), dimension(:), intent(INOUT) :: rho
    type(INAKTIV_def), intent(in) :: inak
    
    integer :: i, e, neigh_num, ne_aktiv
    real(8) :: numer, denum,weigth, rho_filter(ne)

    if (inak%bool) then
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if

    do e = 1, ne_aktiv ! #5 i iterationsskema s. 414

       numer = 0.0d0
       denum = 0.0d0
       do i = 1,int(neigh(e,1))! naboelementer
          neigh_num = neigh(e,2*i)! Naboelementets nummer
          weigth = Neigh(e,2*i+1)
          numer = numer +  weigth* rho(neigh_num)*vol(neigh_num)! Se (12) i filter-paper by sigmund
          denum = denum + weigth*vol(neigh_num)
       end do
       rho_filter(e) = numer / denum
    end do
    rho = rho_filter


  END SUBROUTINE density_filter


  subroutine gradient(flag,problem_type,inak ,D_in,vol,max_vol,rho, rho_min, dc, dg, &
       lambda1,lambda2,dc2, eigenval,eigenvec)

	! This subroutine calculates the gradients:

    use fedata
    use numeth
    use fea
    use plane42
    use plane41
    use plane42nonlin
    use mindlin42

    integer, INTENT(IN) :: flag, problem_type
    real(8), dimension(:), intent(IN) :: vol, rho, D_in
    real(8), intent(IN) :: rho_min, max_vol
    real(8), dimension(:), optional, intent(IN) :: lambda1! lambda1 er til force_inverter
    real(8), dimension(:), optional, intent(IN) :: lambda2
    real(8), dimension(:), intent(OUT) :: dc, dg
    real(8), dimension(:), optional, intent(OUT) :: dc2
    real(8), optional, intent(IN) :: eigenval, eigenvec(:)
    type(INAKTIV_def), intent(in) :: inak

    integer :: e, i, nen, mdim, kk, ne_aktiv
    integer, parameter :: ndim = 4
    integer:: idof(ndim)
    real(8), dimension(ndim) :: helpproduct1,lamb2
    real(8) :: ke_t(ndim,ndim), tn_e(ndim),res_vek(neqn)
    real(8), allocatable, dimension(:) :: xe, de,  helpproduct2, residual
    real(8), allocatable, dimension(:,:) :: ke, me, helpmat
    integer, allocatable :: edof(:)
    real(8) :: lamb1(8), rtemp(8), rtemp2(8)
    real(8) :: young, nu, thk, alpha, tstart, tend, deltaT, kcond, Ae(8,4), fact
    real(8) :: shear, dens

    select case( element(1)%id )
    case(2)
       mdim = 8
       kk = 2
    case default
       mdim = 12
       kk = 3
    end select
    allocate(xe(mdim), de(mdim), helpproduct2(mdim), residual(mdim))
    allocate(ke(mdim,mdim), me(mdim,mdim), helpmat(mdim,mdim), edof(mdim))

    if (problem_type /= 6) then

       if (inak%bool) then
          ne_aktiv = size(inak%act,1)
       else
          ne_aktiv = ne
       end if

       do i = 1, ne_aktiv!#7

          if (inak%bool) then
             e = inak%act(i)
          else
             e = i
          end if

          call get_edof(e,nen,xe,idof,edof)
          
          if (problem_type /= 12) then
             de = D_in(edof)
          end if

          if ( ((flag == 0) .or. (e == 1))) then
             young = mprop(element(e)%mat)%young
             nu    = mprop(element(e)%mat)%nu
             thk   = mprop(element(e)%mat)%thk
             select case( element(e)%id )
             case( 2 )!plane
                call plane42_ke(xe, young, nu, thk, ng, ke)
                select case (problem_type)
                case(2)
                   alpha = mprop(element(e)%mat)%alpha
                   tstart = mprop(element(e)%mat)%tstart
                   tend = t_elem(e)
                   ! Temperature change of element
                   deltaT = tend - tstart  
                   call plane42_thermLoad(xe, ng, young, nu, alpha, thk, deltaT, rtemp)
                case(3)
                   alpha = mprop(element(e)%mat)%alpha

                   call plane42_thermKobling(xe, ng, young, nu, alpha, thk,Ae)
                   kcond = mprop(element(e)%mat)%kcond
                   call plane41_ke_t(xe,ng, kcond, thk, ke_t)
                end select
             case(3:6)!plate
                shear = mprop(element(e)%mat)%shear
                dens  = mprop(element(e)%mat)%dens
                call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke,1,mat_vec)
                call mindlin42_me(xe,dens,mat_vec,thk,ng,me,1)
             end select
          end if

          select case (problem_type)
          case(1)! Force inverter
             lamb1 = lambda1(edof)
             dc(i) = -penal*(1d0-rho_min)*rho(i)**(penal-1d0) *DOT_PRODUCT(lamb1,matmul(ke,de))
          case(2)! mekanisme med "hardcoded" temp-stigning, givet ved t_elem
             lamb1 = lambda1(edof)
             !$$$$$$             dc(e) = penal*(rho(e)**(penal-1)) * DOT_PRODUCT(lamb1, rtemp - MATMUL(ke,de)) 
             dc(i) = penal*(1.0d0-rho_min)*rho(i)**(penal-1d0) * DOT_PRODUCT(lamb1, rtemp - MATMUL(ke,de)) 
             !print*,'rtemp',rtemp(1)
          case(3)! koblet mekanisme
             tn_e = tn(idof) ! knude_temp
             rtemp = MATMUL(Ae,tn_e) ! themisk bidrag til lastvektoren

             tstart = mprop(element(i)%mat)%tstart
             tend = t_elem(i)
             ! Temperature change of element
             deltaT = tend - tstart  
             call plane42_thermLoad(xe, ng, young, nu, alpha, thk, deltaT, rtemp2)

             lamb1 = lambda1(edof) ! lambda2 i papir
             lamb2 = lambda2(idof)

             fact = penal*(1.0d0-rho_min)*rho(i)**(penal-1)
             helpproduct1 = -fact * MATMUL(ke_t , tn_e)
             helpproduct2 = fact * (rtemp2 - MATMUL(ke , de))
             dc(i) = DOT_PRODUCT( helpproduct1, lambda2) + DOT_PRODUCT( helpproduct2, lambda1)
          case(4)
             print *, 'Error in top_help->gradient; Jeg skal ikke kaldes, da det er et termisk topOpt problem'
             error stop
          case(7)! Force inverter med begrï¿½nsning pï¿½ krydskobling
             lamb1 = lambda1(edof)
             dc(i) = -penal*(1d0-rho_min)*rho(i)**(penal-1d0) *DOT_PRODUCT(lamb1,matmul(ke,de))

             lamb1 = lambda2(edof)
             dc2(i) = -penal*(1d0-rho_min)*rho(i)**(penal-1d0) *DOT_PRODUCT(lamb1,matmul(ke,de))

          case(12) !eigenfrequency
             de = eigenvec(edof)
             !normalization of eigenvector wrt. mass matrix
             de = de/SQRT(DOT_PRODUCT(de,matmul(me,de)))

             if (rho(i)<0.1) then
                helpmat =  penal*(1d0-rho_min)*rho(i)**(penal-1d0)*ke - eigenval * 6d0*rho(i)**5 *me
             else
                helpmat =  penal*(1d0-rho_min)*rho(i)**(penal-1d0)*ke - eigenval *me
             end if
             dc(i) = DOT_PRODUCT(de,matmul(helpmat,de))

          case default ! compliance = gradient på "normal" vis.
             dc(i) = - penal*(1d0-rho_min)*rho(i)**(penal-1d0) * DOT_PRODUCT(de,matmul(ke,de))! de=eigenvector
          end select

          dg(i) = vol(i)/max_vol ! skaleret!
       end do

    elseif (problem_type == 6) then !ikke-lineï¿½r Force inverter
       res_vek = 0d0

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

          call plane42nonlin_residual(xe,d(edof), young,nu,thk,ng,residual)
          res_vek(edof) = res_vek(edof) - (penal*(1d0-rho_min)*rho(e)**(penal-1d0)) *residual 
       end do

       ! correct for bounds
       do i = 1, nb
          idof = 2*(bound(i,1)-1) + bound(i,2)
          res_vek(idof) = 0 ! Correct for BC
       end do

       do e = 1, ne
          ! Find coordinates and degrees of freedom
          nen = element(e)%numnode
          do i = 1, nen
             edof(2*i-1) = 2 * element(e)%ix(i) - 1  
             edof(2*i)   = 2 * element(e)%ix(i)
          end do

          lamb1 = lambda1(edof)
          residual = res_vek(edof)

          dc(e) =  DOT_PRODUCT(lamb1,residual)

          dg(e) = vol(e)/max_vol ! skaleret!
       end do
    end if

  end subroutine gradient


  subroutine lambda_calc(flag,problem_type,rho,rho_min,lambda1,lambda2)

	! This subroutine calculates the gradients:

    use fedata
    use numeth
    use fea
    use plane42
    use plane41

    integer, INTENT(IN) :: flag, problem_type
    real(8), dimension(:), intent(IN) :: rho, lambda1
    real(8), intent(IN) :: rho_min
    real(8), dimension(:), intent(OUT) :: lambda2


    integer :: e, i, nen
    integer, parameter :: mdim = 8, ndim = 4
    integer:: edof(mdim), idof(ndim)
    real(8), dimension(mdim) :: xe, de, lamb1
    real(8) :: young, nu, thk, alpha, Ae(8,4), Ae_opr(8,4)

    do e = 1, ne

       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
          idof(i)	= element(e)%ix(i)
       end do
       de = d(edof)
       if ( (flag == 0) .or. (e == 1) ) then
          young = mprop(element(e)%mat)%young
          nu    = mprop(element(e)%mat)%nu
          thk   = mprop(element(e)%mat)%thk
          alpha = mprop(element(e)%mat)%alpha
          call plane42_thermKobling(xe, ng, young, nu, alpha, thk,Ae_opr)
       end if
       select case (problem_type)
       case(3)! koblet mekanisme
          lamb1 = lambda1(edof) ! lambda2 i papir
          Ae = (rho_min+(1.0d0-rho_min)*rho(e)**penal) * Ae_opr ! koblingsmatrice. Dvs d f / d t
          lambda2(idof) = MATMUL(lamb1,Ae) ! lambda1 i papir
       case default 
          print *, 'Error in top_help->gradient; Jeg skal ikke kaldes, da det er et termisk topOpt problem'
          stop
       end select
    end do
    call bsolve(k_t,lambda2)

  end subroutine lambda_calc


  subroutine finite_check(flag,problem_type,inak, filter_type,vol,max_vol, rho,rho_min2, L, lambda, dc,&
       df_approx,f_scale,beta,eta,parameters,deltaT,U0,dotU0,dof_x, Lmulti)

    use numeth
    use processor
    use fedata
    use plane42
    use fea
    use thermal
    use transient
    use nonlin

  	! This subroutine calculates the finite difference check of the gradients:

	integer, intent(IN) :: flag,  filter_type, problem_type
    type(INAKTIV_def), intent(in) :: inak
  	REAL(8), DIMENSION(:), INTENT(IN) :: dc, rho, vol, L, lambda
   !$$$$$$     REAL(8), DIMENSION(:,:) L, lambda
    real(8), dimension(:,:), intent(OUT) :: df_approx
    real(8), intent(in) ::  rho_min2, max_vol, f_scale
    real(8), optional, intent(in) :: beta, eta
    real(8), optional, dimension(:,:), intent(in) :: Lmulti

    integer :: e
	real(8) :: fa, fb, delta_pert
    real(8), dimension(ne) :: rho_pert, fd_dc, dg,dc_analytisk, rho_old, rho_pert_tilde
    real(8), dimension(ne) :: dc_dummy

    ! Ekstra input tilfï¿½jet til transient
    integer, optional, intent(inout) :: parameters(:)
    integer, optional, intent(in) :: dof_x(:)
    real(8), optional, intent (in) :: deltaT, U0(:), dotU0(:)
    real(8), allocatable :: saved_U(:,:)
    real(8) :: objekt, dummy_mat(1,1) = 0.0d0 ! til kald af transient
    integer :: i,n,nmax

    ! tilfï¿½jet pga. thermal - skal "optimeres"
    real(8) :: compliance, rho_min
    !linux_fejl : rho_min i input er omdøbt til rho_min2, fordi der var en fejl ved kompilering

    ! til transient
    select case (problem_type) ! dï¿½rlig lï¿½sning
    case (10)
    case default
       if (present(parameters)) then
          parameters(3) = 0 ! animation
          nmax = parameters(6)
          allocate( saved_U(neqn,nmax+1) )
       end if
    end select


	delta_pert = 0.0000001d0
    rho_pert = rho

    select case( filter_type )
    case(0:1)

    case(2)
       dg = vol
       ! filter rho to obtain rho_hat
       call density_filter(inak,vol, rho_pert)

       select case (problem_type)

       case default
          !call gradient(flag,problem_type,0,D,vol,max_vol,rho_pert,rho_min, dc_analytisk, dg, lambda)

       case(4)!sensitiviteterne beregnes i en anden routine for termiske problemer
          call object_sens_t(flag,rho_pert,rho_min,compliance,dc_analytisk)

       case(10)
          dc_analytisk = dc

       end select

       call filter(filter_type,inak,vol,rho_pert,dc_analytisk,dg) ! filtrer dc; finder dc analytisk-

    case(3)
       dg = vol

       ! filter rho to obtain rho_tilde
       call density_filter(inak,vol, rho_pert)

       do e=1,ne ! rho_tilde -> rho_bar
          rho_pert_tilde(e) = rho_pert(e) ! Oprindelig densitetsvï¿½gtet designvar. Skal bruges i filter ved kï¿½deregel
          rho_pert(e) = EXP(-beta*(1.0-rho_pert(e))) - (1.0-rho_pert(e))*EXP(-beta) ! Modificeret heavyside funktion (29) -> rho_bar/rho_tilde
       end do

       select case (problem_type)

       case default
          !call gradient(flag,problem_type,0,D,vol,max_vol,rho_pert,rho_min, dc_analytisk, dg, lambda)

       case(4)!sensitiviteterne beregnes i en anden routine for termiske problemer
          call object_sens_t(flag,rho_pert,rho_min,compliance,dc_analytisk)              

       case(10)
          dc_analytisk = dc

       end select

       call filter(filter_type,inak,vol,rho_pert_tilde,dc_analytisk,dg,beta) ! filtrer dc; finder dc analytisk -> heavyside bruger rho_hat/rho__tilde

    case(4) ! ROBUST virker pt. kun med problem (9)
       dg = vol
       call density_filter(inak,vol, rho_pert) ! filter rho to obtain rho_tilde
       do e=1,ne
          rho_pert_tilde(e) = rho_pert(e)
          rho_pert(e) = (dtanh(beta*eta)+dtanh(beta*(rho_pert_tilde(e)-eta)))/&
               (dtanh(beta*eta)+dtanh(beta*(1-eta)))
       end do

       select case (problem_type)

       case default
          !call gradient(flag,problem_type,0,D,vol,max_vol,rho_pert,rho_min, dc_analytisk, dg, lambda)

       case(4)!sensitiviteterne beregnes i en anden routine for termiske problemer
          call object_sens_t(flag,rho_pert,rho_min,compliance,dc_analytisk)

       case(10)
          dc_analytisk = dc

       end select

       call filter(filter_type,inak,vol,rho_pert_tilde,dc_analytisk,dg,beta,eta) ! filtrer dc; finder dc analytisk -> heavyside bruger rho_hat/rho__tilde

    end select

    rho_old = rho_pert

    do e = 1,ne ! Der benyttes central-differens og IKKE forward-differens som i noterne fra FEM-kurset
       ! NB kun ï¿½t element ï¿½ndres ad gangen. Resten holdes "oprindelig". Derfor nulsï¿½ttes rho_pert i bunden af lï¿½kken
       ! For rhoa:
       rho_pert(e) = rho_pert(e) - delta_pert
       if ( rho_pert(e) < 0) then
          print*,'ERROR, rho < 0 i fd_check. element nr: ', e
          stop
       end if



       ! Objektfunktion
       select case (problem_type)
       case(0)
          !buildload kaldes ikke displ, sï¿½ hvis den er design-afhï¿½ngig skal den kaldes nu
          !buildstiff kaldes fra displ
          call displ(flag,rho_pert,rho_min)
          fa = DOT_PRODUCT(d,p)
       case(1)! Force inverter
          call displ(flag,rho_pert,rho_min)
          fa = dot_product(d,L)
       case(2)! mekanisme med "hardcoded" temp-stigning, givet ved t_elem
          call buildload(rho_pert,rho_min)
          call displ(flag,rho_pert,rho_min)
          fa = DOT_PRODUCT(d,L)
       case(3)
          call buildtermload
          call therm_displ(flag,rho_min,rho_pert)
          call buildload(rho_pert,rho_min)
          call displ(flag,rho_pert,rho_min)
          fa = DOT_PRODUCT(d,L)
       case(4) !rent thermisk
          call therm_displ(flag,rho_min,rho_pert)
          call object_sens_t(flag,rho_pert,rho_min,compliance,dc_dummy) ! object and gradient
          fa = compliance
       case(5) ! transient
          !call inputfile(...) laves senere
          saved_U = 0.0d0
          call half_step_CD(parameters,(/0.0d0/),flag,deltaT,rho_pert,U0,dotU0,saved_U) ! STEP 1
          ! Objekt-funktionen er integralet af "forskydning af endeknuden", dvs d(neqn-1)^2, over tiden
          objekt = 0.0d0
          do n=1,nmax ! Objekt findes ved trapez-integration
             if ( (n > 1) .and. (n < nmax) ) then 
                objekt = objekt + 2.0d0* dot_product( saved_U(:,n) , L* saved_U(:,n))
             else
                objekt = objekt + dot_product( saved_U(:,n) , L* saved_U(:,n))
             end if
          end do
          objekt = (deltaT)/(2.0d0) * objekt
          fa = objekt
       case(6)
          call non_lin(flag,rho_pert,rho_min)
          fa = dot_product(d,L)
       case(9)! Force inverter
          call displ(flag,rho_pert,rho_min)
          fa = dot_product(d,L)
       case(10)
          call displ(flag,rho,rho_min) ! #6, solve system
          fa = DOT_PRODUCT(d,Lmulti(:,1)) ! Compliance for forskydning 1
          objekt = DOT_PRODUCT(d,Lmulti(:,2)) ! objekt er blot hjï¿½lpestï¿½rrelse
          fa = fa+objekt+0.1*(fa-objekt)**2
       case default ! gradient pï¿½ "normal" vis.
          fa = DOT_PRODUCT(d,p)
       end select

       ! For rhob:
       rho_pert(e) = rho_pert(e) + 2.0d0*delta_pert

       ! Objektfunktion
       select case (problem_type)
       case(0)
          !buildload kaldes ikke displ, sï¿½ hvis den er design-afhï¿½ngig skal den kaldes nu
          !buildstiff kaldes fra displ
          call displ(flag,rho_pert,rho_min)
          fb = DOT_PRODUCT(d,p)
       case(1)! Force inverter
          call displ(flag,rho_pert,rho_min)
          fb = dot_product(d,L)
       case(2)! mekanisme med "hardcoded" temp-stigning, givet ved t_elem
          call buildload(rho_pert,rho_min)
          call displ(flag,rho_pert,rho_min)
          fb = DOT_PRODUCT(d,L)
       case(3)
          call buildtermload
          call therm_displ(flag,rho_min,rho_pert)
          call buildload(rho_pert,rho_min) 
          call displ(flag,rho_pert,rho_min)
          fb = DOT_PRODUCT(d,L)
       case(4) !rent thermisk
          call therm_displ(flag,rho_min,rho_pert)
          call object_sens_t(flag,rho_pert,rho_min,compliance,dc_dummy) ! object and gradient
          fb = compliance
       case(5) ! transient
          saved_U = 0.0d0
          call half_step_CD(parameters,(/0.0d0/),flag,deltaT,rho_pert,U0,dotU0,saved_U)
          objekt = 0.0d0 ! (7) i paper
          do n=1,nmax ! Objekt findes ved trapez-integration
             if ( (n > 1) .and. (n < nmax) ) then 
                objekt = objekt + 2.0d0* dot_product( saved_U(:,n) , L* saved_U(:,n))
             else
                objekt = objekt + dot_product( saved_U(:,n) , L* saved_U(:,n))
             end if
          end do
          objekt = (deltaT)/(2.0d0) * objekt
          fb = objekt
       case(6)
          call non_lin(flag,rho_pert,rho_min)
          fb = dot_product(d,L)
       case(9)! Force inverter
          call displ(flag,rho_pert,rho_min)
          fb = dot_product(d,L)
       case(10)
          call displ(flag,rho,rho_min) ! #6, solve system
          fb = DOT_PRODUCT(d,Lmulti(:,1)) ! Compliance for forskydning 1
          objekt = DOT_PRODUCT(d,Lmulti(:,2)) ! objekt er blot hjï¿½lpestï¿½rrelse
          fb = fb+objekt+0.1*(fb-objekt)**2
       case default ! gradient pï¿½ "normal" vis.
          fb = DOT_PRODUCT(d,p)
       end select

       ! Finite Difference
       df_approx(e,1) = (fb - fa)/(2.0d0*delta_pert) ! "beregnede" sensitivities
       df_approx(e,1) = df_approx(e,1)/f_scale       

       rho_pert(e) = rho_pert(e) - delta_pert ! nulstil det ï¿½ndrede element

 	end do

    select case( filter_type )
    case( 0:1 )
       df_approx(1:ne,2) = dc/f_scale ! analytiske sensitivities, skaleret FILTRERET
    case(2)
       !de beregnede sensitiviteter skal filtreres
       fd_dc(1:ne) = df_approx(1:ne,1)
       call filter(filter_type,inak,vol,rho_pert,fd_dc,dg) ! filtering af beregnede sensitivities

       df_approx(1:ne,1) = fd_dc !filtrerede beregnede sensitivities
       df_approx(1:ne,2) = dc_analytisk(1:ne)/f_scale  ! analytiske sensitivities
    case(3)
       ! de beregnede sensitiviteter skal filtreres
       fd_dc(1:ne) = df_approx(1:ne,1)
       call filter(filter_type,inak,vol,rho_pert_tilde,fd_dc,dg,beta) ! filtering af beregnede sensitivities

       df_approx(1:ne,1) = fd_dc !filtrerede beregnede sensitivities
       df_approx(1:ne,2) = dc_analytisk(1:ne)/f_scale  ! analytiske sensitivities
    case(4)
       ! de beregnede sensitiviteter skal filtreres
       fd_dc(1:ne) = df_approx(1:ne,1)
       call filter(filter_type,inak,vol,rho_pert_tilde,fd_dc,dg,beta,eta) ! filtering af beregnede sensitivities

       df_approx(1:ne,1) = fd_dc !filtrerede beregnede sensitivities
       df_approx(1:ne,2) = dc_analytisk(1:ne)/f_scale  ! analytiske sensitivities


  	end select

  end subroutine finite_check


  subroutine force_elements_init ! (Tag evt. bï¿½de 'force_elements'& 'force_rho(rho)' rutinerne herunder...)
    ! Her findes elementer, hvor krï¿½fter virker: 
    ! Alle elementer gennemgï¿½s og de for hvilke en knudepunkt svarende til et knudepunkt i loads registreres.

    use fedata

    integer :: e, n, nen, i,elem
    real(8) :: center(ne,2), xLength,yLength,radius,dist_ei

    i = 0
    do e=1,ne ! Jeg tï¿½ller fï¿½rst...
       nen = element(e)%numnode
       do n=1,nen
          if ( element(e)%ix(n) == springs(1,2) ) then
             i = i + 1
             elem = e
          end if
       end do
    end do



    center = 0.0d0
	do e=1,ne
       nen = element(e)%numnode
       do n = 1,nen
          center(e,1) =center(e,1)+ 1.0/4.0 * x(element(e)%ix(n), 1) ! x-koordinat, Centrum af element e
          center(e,2) =center(e,2)+ 1.0/4.0 * x(element(e)%ix(n), 2)
       end do
   	end do

    xLength = real(x(element(1)%ix(2),1)) - real(x(element(1)%ix(1),1)) 
    yLength = real(x(element(1)%ix(2),2)) - real(x(element(1)%ix(1),2))
    radius = dsqrt( xLength**2 + yLength**2) * 2

    do e = 1,ne
       dist_ei = dsqrt( (center(e,1)-center(elem,1))**2 + (center(e,2)-center(elem,2))**2 ) ! Afstand mellem centre af e og f
       if (dist_ei <= radius .and. e /= elem) then
          i = i + 1
       end if
    end do

    allocate(force_elements(i))

    i = 0
    do e=1,ne ! Nu noterer jeg ogsï¿½...
       nen = element(e)%numnode
       do n=1,nen
          if ( element(e)%ix(n) == springs(1,2) ) then
             i = i + 1
             force_elements(i) = e
          end if
       end do
    end do

    do e = 1,ne
       dist_ei = dsqrt( (center(e,1)-center(elem,1))**2 + (center(e,2)-center(elem,2))**2 ) ! Afstand mellem centre af e og f
       if (dist_ei <= radius .and. e /= elem) then
          i = i + 1
          force_elements(i) = e
       end if
    end do

    print*,'e'



  end subroutine force_elements_init

  subroutine force_rho(rho)

    ! Gennemtvinger at rho lig 1 i elementer, hvor krï¿½fter virker.

    use fedata

	real(8), dimension(:), intent(INOUT) :: rho
	integer :: i

	do i = 1,size(force_elements,1)
       rho(force_elements(i))=1.0d0
  	end do

  end subroutine force_rho



END MODULE top_help
