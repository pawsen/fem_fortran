MODULE topology

  ! 

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: topopt, topopt_trans, inactive_elements, elements_equal

CONTAINS

  subroutine topopt(flag)

    use fea
    use thermal
    use numeth
    use processor
    use fedata
    use top_help
    use file_init
    use plot_routiner
    use nonlin
    use arpack
    use piezo
    use solve_handle_real
    use exodus
    use plate, only : buildstiff_plate

    integer, intent(IN) :: flag
    ! generelt fil
    integer :: i,j, e, iter,n_iter, itend, idof, idof2, nc, ne_aktiv
    integer :: n, na

    ! TopOpt - stuff
    real(8), dimension(:), allocatable :: rho,rho_tilde, vol, dc, dg, dc_filter
    real(8), dimension(:), allocatable :: plotval
    real(8), dimension(:), allocatable :: xval, xval_mat, xval_old, xval_old2
    real(8), dimension(:), allocatable :: dobject, constraint
    real(8) :: rho_min

    ! Parametre
    integer :: file, filter_type, solver_type,vol_type, rmin_type, fd_check, animation, info
    integer :: problem_type, save_rho, stop_krit
    real(8) :: max_vol,total_vol, rmin,element_areal, tol, movelimit ! parametre til hver fil
    real(8) :: beta, change, Mnd, parameters(20)
    real(8), dimension(neqn) :: L, lambda1 ! FORCE inverter

    ! fd-check
    REAL(8), DIMENSION(:,:), allocatable :: df_approx
    real(8), dimension(10) :: para

    ! skalering
    real(8) :: f_scale

    ! min(max(compliance))
    real(8), dimension(:), allocatable :: L2, lambda2
    real(8), dimension(:,:), allocatable :: dconstraint, P_mat, D_mat,compliance_mat

    ! Robust design
    real(8) :: eta(3), max_vol_start, beta1, beta2, zeta
    integer :: iter1, iter2
    real(8), dimension(:,:), allocatable :: rho_bar,D_robust
    real(8), dimension(:,:), allocatable :: lambda_mat

    ! Inaktive elementer
    type(INAKTIV_def), allocatable :: inak

    ! Mekanisme med flere outputs
    real(8), dimension(:,:), allocatable :: Lmulti
    real(8) :: rpenal

    !eigenfrequency
    real(8) :: sigma, start_guess
    integer :: n_eigen, n_conv
    logical :: shift
    real(8), allocatable :: eigenval(:,:), eigenvec(:,:)
    integer, allocatable ::list(:)
    logical, parameter :: scale = .true.


    ! Mn_iter = 100MA 
    REAL(8), DIMENSION(:), allocatable:: low, upp
    REAL(8):: compliance

    real(8) :: a1,a2,a3,a4,a5,a6,a7,a8

    start_guess = 233101.
    !start_guess = 50d0

    ! TOPOPT & FILTER INITIAL:
    file 		= 5		!  1:bridge1 			2 = MBB5, 3: T01, 4 = MBB fra filter-artikel
    filter_type = 2		! 0: No filter 			1: Sensivity filter. 2: Density
    solver_type = 2		! 1: Bisection			2: MMA
    vol_type 	= 2		! 1: absolut volumen	2: Volfrac
    rmin_type 	= 2		! 1: absolut afstand	2: afstand i forhold til elementets stï¿½rrelse
    fd_check   	= 1		! 0: Fra				1: Finite difference check
    animation 	= 0		! 0: Fra				1: Animation		2: ï¿½ndring af rho
    info 		= 1		! information
    problem_type = 1	! 0: compliance			1: forceinverter	2: mekanisme thermal	3: mekanisme koblet termisk		4: thermal
    !flag				! sættes i main
    save_rho	= 0		! 0: fra				1: gem rho som filename_rho
    stop_krit	= 1		! 0: fra				1: brug stopkriterie

    n_iter = 1000
    tol =1.0E-3
    penal = 3.0
    rho_min = 0.001 ! til thermal
    max_vol = 0.25d0 ! forceinverter
    rmin = 1.3d0! The filter radius percent of max element side length
    damp_fact = 0.6 ! til OC

    ! hvis en fil med "filename_para.txt" eksisterer, bliver parametrene indlæst her.
    call inputfile(file,filter_type,solver_type,problem_type,vol_type,rmin_type,info,save_rho,n_iter,&
         stop_krit,animation,penal,damp_fact,max_vol,rmin,rho_min,tol,movelimit)

    beta = 1.0
    iter1 = 50
    iter2= 2000
    beta1=1.0d0
    beta2=100.0d0
    zeta=1e-2
    
    allocate(inak)
    inak%bool = .true.
    !inak%bool = .false.

    ! Inaktive elementer
    if (inak%bool) then
       j = 0
       do e=1,ne
          if (element(e)%mat /= 2) then ! aktivt element
             j = j+1
          end if
       end do
       allocate(inak%act(j),inak%inact(ne-j))
       j = 0
       i = 0
       do e=1,ne
          if (element(e)%mat /= 2) then
             j = j+1
             inak%act(j) = e
          else ! inaktive elements
             i = i+1
             inak%inact(i) = e
          end if
       end do
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if
    
    select case(problem_type)
    case(12)! bound formulation: we introduce a new designvariabel, 
       n = ne_aktiv+1 ! total number of design variabel
       j = ne_aktiv ! number of active elements, eg to be filtered
    case default
       n = ne_aktiv
       j = ne_aktiv
    end select

    allocate(rho(j),rho_tilde(j),vol(j), dc(j), dg(j), dc_filter(j))
    allocate(xval(n), xval_old(n), xval_old2(n), low(n), upp(n))
    allocate(dobject(n), df_approx(j,2))
    allocate( xval_mat(n) ) ! skal bruges til mma

    allocate ( compliance_out(n_iter), plotval(ne) )
    compliance_out = 0.0

    if (inak%bool) then
       call volume(flag,inak,vol)
    else
       call volume(flag,inak,vol)
    end if

    total_vol = SUM(vol)
    select case( vol_type )
    case( 2 ) ! volfrac
       max_vol = total_vol*max_vol
    end select
    rho = max_vol/(total_vol*1.2) ! Hvert element(uanset størrelse) tildeles samme "densitet"/"designvariabel-værdi".

    if(rho(1)>1) then ! Tjek at tildelt densitet ikke er større end 1.
       write (*, *) 'Error: max_vol/vol_frac er for stor, torsk!'
       error stop
    end if
    
    na = ne_aktiv
    xval(1:na) = rho
    select case(problem_type)
    case(12)
       xval(na+1) = start_guess! gæt på beta
    end select

    select case( solver_type)
    case(1)
       print*, 'OC'
    case(2)
       xval_old = xval
       xval_old2 = xval
       low = 0.0d0
       upp = 1.0d0
       print*, 'MMA-solver'
    end select

    ! INIT. AF NEIGHBOURHOOD MATRIX SAMT FILTERAFHÆNGIGE STØRRELSER:
    ! skal ændres så "flag" også tages som input
    select case( filter_type )
    case( 1:2 )
       call neighbourMatrix(inak,rmin_type,rmin)
    case( 3 )
       call neighbourMatrix(inak,rmin_type,rmin)
    case( 4 )
       max_vol_start = max_vol
       call neighbourMatrix(inak,rmin_type,rmin)
       allocate(D_robust(neqn,3),rho_bar(ne,3),compliance_mat(n_iter,3))
       eta(1) = 0.2d0
       eta(2) = 0.5d0
       eta(3) = 0.8d0
       !eta = 1d0
    end select

    select case( animation )
    case(1)
       call exodus_init
    end select

    ! Test af indlæsning/skrivning af rho. Det virker
    !$$$$$$     call output_vector(rho,'rho')
    !$$$$$$     call rho_input(rho_vec)
    !$$$$$$     call output_vector(rho_vec,'rho_ny')

    print*,'max_vol',max_vol
    print*,'rmin',rmin


    !##############################
    ! Initialisering -choose problem-type:
    L = 0.0d0
    select case (problem_type)
    case(0) ! statisk problem
       call buildload !Buildload kaldes her da den ikke er design-afhængig. Dvs den skal kun kaldes een gang
       if (banded == 2) then
           call buildstiff_fea(flag,rho,rho_min)
           call mumps_init_real
           call mumps_solve_real(1)
       end if
    case(1) ! (1): force-inverter
       call buildload
       solver_type = 2
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             L(idof) = 1.0d0
          end if
       end do

       ! NB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! find elemeter der sidder hvor fjederen er monteret
       call force_elements_init
    case(2)!(2)mekanisme med "hardcoded" termisk last
       allocate( t_elem(ne_aktiv) )
       t_elem = 100.0 ! temp-stigning
       solver_type = 2 ! MMA
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             ! sat minus for at maksimere forskydning
             L(idof) = -1.0d0
          end if
       end do
       ! find elemeter der sidder hvor fjederen er monteret
       call force_elements_init
    case(3)! mekanisme med koblet termisk last
       allocate(lambda2(nn))! Bemï¿½rk der er byttet om  pï¿½ lambda1 og lambda2 i forhold til nummerering i bogen
       allocate( t_elem(ne) )
       solver_type = 2
       do i = 1, nk
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             L(idof) = 1.0d0
          end if
       end do
    case(4) ! Thermal
       call buildtermload

    case(6) ! (1): geometrisk ikke-lineær
       call buildload
       solver_type = 2
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             L(idof) = 1.0d0
          end if
       end do
    case(7) ! Mekanisme med flere constraint
       allocate(L2(neqn),lambda2(neqn))
       allocate(constraint(1),dconstraint(ne,1))
       call buildload
       nc = 2
       solver_type = 2
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             L = 0d0
             L(idof) = 1.0d0
             if (springs(i, 3) == 1) then
                idof2 = 2*(springs(i, 2)-1) + 2 ! krydsretningen er y-retningen
             else
                idof2 = 2*(springs(i, 2)-1) + 1 ! krydsretningen er x-retningen
             end if
             L2 = 0d0
             L2(idof2) = 1d0
          end if
       end do

    case(8) ! min( max( compliance))
       call buildload
       solver_type = 2
       nc = nk+1
       allocate(constraint(nk),dconstraint(ne,nc))
       allocate(P_mat(neqn,nk),D_mat(neqn,nk),compliance_mat(n_iter,nk))
       P_mat = 0d0
       D_mat = 0d0
       compliance_mat = 0d0
       dconstraint = 0d0
       ! opdel de påførte krafter i P1/P2/P3... ud fra "fjeder-inputtet"
       do i = 1, nk
          idof = 2*(springs(i, 2)-1)+springs(i, 3)
          P_mat(idof,i) = P(idof)
       end do
    case(9) ! Robust force-inverter
       call buildload
       allocate(lambda_mat(neqn,3))
       nc = 4 ! kun tre udover volumen!
       allocate(constraint(nc-1),dconstraint(ne,nc))
       dconstraint = 0d0
       solver_type = 2
       L = 0.0d0
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             L(idof) = 1.0d0
          end if
       end do
    case(10) ! Mekanisme med parrallelforskydning - maximerer forskydning
       call buildload
       solver_type = 2
       rpenal = 1d0
       allocate(Lmulti(neqn,2),lambda_mat(neqn,2),dconstraint(ne,2),constraint(2)) ! her bruges constraint m.m. blot som hjï¿½lpestï¿½rrelser
       Lmulti = 0d0
       j=0
       do i = 1, nk
          ! Add spring contribution
          if (springs(i, 5) == 1) then ! input/output => 0/1
             j=j+1
             idof = 2*(springs(i, 2)-1)+springs(i, 3)
             Lmulti(idof,j) = -1.0d0
          end if
       end do
       nc = 1 ! kun volumen!
    case(11) ! Thermal
       call buildtermload
       allocate(lambda_mat(neqn,3))
       nc = 4 ! kun tre udover volumen!
       allocate(constraint(nc-1),dconstraint(ne,nc))
       dconstraint = 0d0

    case(12) !eigenfrequency
       call mumps_init_real

       eigenvalue%calc = .true.
       shift = eigenvalue%shift
       sigma = eigenvalue%sigma ! shift value
       n_eigen = 4
       allocate(eigenval(n_eigen,3), eigenvec(neqn,n_eigen))

       nc = n_eigen+1 !n_eigen + volumen!
       n = ne_aktiv +1
       allocate(constraint(nc-1),dconstraint(n,nc))
       allocate(list(2))!contains id of double eigenvalue
    end select


    !##############################
    ! TOPOLOGY OPTIMIZATION LOOP
    iter = 0
    f_scale = 1d0
    do i = 1,n_iter
       ! sætter rho hvor fjederen er monteret til 1.
       select case( problem_type)
       case(1)
          !$$$$$$                 call force_rho(rho)
       end select

       iter = i
       
       !rho_mat = rho ! MMA skal have den "oprindelige" rho ind og ikke rho_hat. Derfor gemmes den her
       ! Vi bruger nu xval. Og da den ikke filtreres, er det ikke nødvendigt at kopiere xval over mere
       select case (filter_type) ! overvej at ligge denne ned efter mma, routinen. Her gï¿½r den ingen godt!
       case( 2)
          ! Computation of filteret densities rho_hat
          call density_filter(inak,vol, rho) ! #5 i paper
       case ( 3 )
          call density_filter(inak,vol, rho)
          rho_tilde = rho ! skal bruges i "filter" pga. kï¿½deregel
          do e=1,ne_aktiv
             rho(e) = EXP(-beta*(1.0-rho(e))) - (1.0-rho(e))*EXP(-beta) ! Modificeret heavyside (29): rho -> rho_bar
             if (rho(e) < 0.0d0) then
                rho(e) = 0.0d0
             end if
          end do
       case(4) ! densitetsfilter med robustprojektion
          call density_filter(inak,vol, rho)
          rho_tilde = rho ! skal bruges i "filter" pga. kï¿½deregel
          do j = 1,3
             do e=1,ne_aktiv
                rho_bar(e,j) = (dtanh(beta*eta(j))+dtanh(beta*(rho_tilde(e)-eta(j))))/&
                     (dtanh(beta*eta(j))+dtanh(beta*(1-eta(j))))
             end do
          end do
          rho = rho_bar(:,1)
       end select
       
       ! Animation:
       
       select case( animation )
       case(1)
          if (inak%bool) then
             plotval(inak%act) = rho
             plotval(inak%inact) = 1d0
          else
             plotval = rho
          end if

          call exodus_write_elem(i,plotval)
          call exodus_write_time(i,real(i,8))

          call  plotanim(iter, 0, 1, .true., .false., .false., .true.,.true., 'anim', &
               0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1.0d0,plotval)
       case(2) ! change in rho
       end select

       !##############################
       ! COMPLIANCE. Solve system and find the value of the objective function
       select case (problem_type)
       case(0) ! static
          call displ(flag,rho,rho_min) ! #6, solve system
          compliance_out(i) = DOT_PRODUCT(d,p) ! Compliance
       case(1) ! force_inverter
          call displ(flag,rho,rho_min) ! #6, solve system
          compliance_out(i) = DOT_PRODUCT(d,L) ! Compliance for force inverter
          ! SOLVE
          lambda1 = L
          call bsolve(k,lambda1)
       case(2) ! mekanisme med "hardcoded" temp-stigning, givet ved t_elem
          antype = 'TOPTHERMSTRUCT_hard'
          call buildload(rho,rho_min) ! den termiske last er designafhï¿½ngig=> Pt = [B]*[C(rho)]*[alpha]*dT
          call displ(flag,rho,rho_min)
          compliance_out(i) = DOT_PRODUCT(d,L) ! Compliance for force inverter
          lambda1 = L
          call bsolve(k,lambda1)
       case(3) ! mekanisme med koblet termisk last
          call buildtermload(rho)
          call therm_displ(flag,rho_min,rho)
          call buildload(rho,rho_min) 
          call displ(flag,rho,rho_min)
          compliance_out(i) = DOT_PRODUCT(d,L)
          lambda1 = L
          call bsolve(k,lambda1)
          call lambda_calc(flag,problem_type,rho,rho_min,lambda1,lambda2)
       case(4) ! thermal
          call therm_displ(flag,rho_min,rho)
          call object_sens_t(flag,rho,rho_min,compliance,dc) ! object and gradient
          compliance_out(i) = compliance
          dg = vol/max_vol
       case(6) ! geometrisk ikke-lineï¿½r
          call non_lin(flag,rho,rho_min) ! #6, solve system
          compliance_out(i) = DOT_PRODUCT(d,L) ! Compliance for force inverter
          ! SOLVE
          lambda1 = L
          ! bemï¿½rk at k er tangent stiffness
          call bsolve(k,lambda1)
       case(7)! mekanisme med begrænsning på krydsbevægelse
          call displ(flag,rho,rho_min) ! #6, solve system
          compliance_out(i) = DOT_PRODUCT(d,L) ! Compliance for force inverter
          ! SOLVE
          lambda1 = L
          call bsolve(k,lambda1)
          lambda2 = L2
          call bsolve(k,lambda2)
          constraint = (dot_product(D,L2)/dot_product(D,L))**2
       case(8) ! min(max(compliance))P = P1
          ! call displ kan ikke bruges da displ laver stivhedsmatricen hver gang
          call buildstiff_fea(flag,rho,rho_min)
          call enforce_fea
          call bfactor(k)                
          do j=1,nk
             D_mat(:,j) = P_mat(:,j)
             call bsolve(k, D_mat(:,j))
             constraint(j) = DOT_PRODUCT(D_mat(:,j),P_mat(:,j)) ! Compliance
             compliance_mat(i,j) = constraint(j)
          end do
          compliance_out(i) = maxval(dabs(constraint)) !
       case(9) ! Robust design for force_inverter
          do j=1,3
             call displ(flag,rho_bar(:,j),rho_min) ! #6, solve system
             D_robust(:,j) = D
             constraint(j) = DOT_PRODUCT(D,L)
             lambda_mat(:,j) = L
             call bsolve(k,lambda_mat(:,j))
             compliance_mat(i,j) = constraint(j)
          end do
          constraint = constraint + 10d0 ! pga mma
          compliance_out(i) = maxval(constraint) ! Compliance er maxvï¿½rdien af de tre compliance
       case(10) ! Elevator
          call displ(flag,rho,rho_min) ! #6, solve system
          do j=1,2
             constraint(j) = DOT_PRODUCT(d,Lmulti(:,j)) ! Compliance for forskydning j
             ! SOLVE
             lambda_mat(:,j) = Lmulti(:,j)
             call bsolve(k,lambda_mat(:,j))
          end do
          compliance_out(i) = SUM(constraint) + rpenal*(constraint(1)-constraint(2))**2
       case(11) ! thermal med robust
          do j=1,3
             call therm_displ(flag,rho_min,rho_bar(:,J))
             call object_sens_t(flag,rho_bar(:,J),rho_min,compliance,dc) ! object and gradient
             dconstraint(:,j) = dc
             constraint(j) = compliance
             compliance_mat(i,j) = constraint(j)
          end do
          do e=1,ne
             dg(e) = vol(e)/max_vol
          end do
          compliance_out(i) = maxval(constraint) ! Compliance er maxværdien af de tre compliance
       case(12)! eigenfrequency
          na = ne_aktiv
          n = na +1
          call build_mvec(rho) !build mass-vector
          if (elem_type == 'PLATE_GMSH') then
             call buildstiff_eigenvalue_piezo(shift,rho,rho_min)
             !call buildstiff_plate(rho,rho_min)
          else
             call buildstiff_fea(0,rho,rho_min)
          end if

          if (iter == 1) then
             call mumps_solve_real(1)! symbolic factorizing
          end if 
          call mumps_solve_real(2)! numerical factorizing
          call arpack_plane(n_eigen,neqn,shift,sigma,iK,jK,sK,nnz_ub, &
               n_conv,eigenval,eigenvec)

          ! scale eigenvector. The derivations of constraints require that \phi^T * [M] * \phi = 1
          do j = 1,n_eigen
             eigenvec(:,j) = eigenvec(:,j)/ &
                  SQRT( DOT_PRODUCT(eigenvec(:,j),mvec*eigenvec(:,j)) )
          end do
          

          if (iter == 1) then
             if (scale) then
                xval(n) = log(1.0 / eigenval(1,1)) +100
             else
                xval(n) = eigenval(1,1)
             end if
             xval_old(n) = xval(n)! beta from MMA
             xval_old2(n) = xval(n)
          end if
          if ( scale ) then
             !Compliance_out(i) = log( 1.0/xval(n) )
             compliance_out(i) = xval(n)
             Compliance_out(i) =Compliance_out(i) +100 ! Make sure that compliance is negative
          else
             compliance_out(i) = - xval(n)
          end if

          print*,'xval(n): ',xval(n)
          do j=1,nc-1 ! for all eigenvalue constraint
             ! scale constraint
             constraint(j) = xval(n)/eigenval(j,1)-1
             print*,'constraint: ',constraint(j)
             ! Constraint skaleres når gradienten er udregnet
          end do
          do j=1,nc-1
            ! print*,'eigenval: ', eigenval(j,1)
          end do
          

       end select

       if (( MOD(iter,50) == 0 ) .OR. ( iter == 1 )) then
          f_scale = 1d0!dabs(compliance_out(i))
       end if


       ! CALCULATION OF THE GRADIENTS:
       select case( problem_type )
       case(0) ! statik
          call gradient(flag,problem_type,inak,D,rho,rho_min, dc,vol,max_vol, dg) ! #7
          dobject = dc
       case(1:2,6) !force_inverter
          call gradient(flag,problem_type,inak,D,rho,rho_min, dc,vol,max_vol, dg, lambda1)
          dobject = dc
       case(3) ! koblet mekanisme
          call gradient(flag,problem_type,inak,D,rho,rho_min, dc,vol,max_vol, dg, lambda1,lambda2)
          dobject = dc
       case(7)!mekanisme med begrænsning på krydsbevægelsen
          call gradient(flag,problem_type,inak,D,rho,rho_min, dc,vol,max_vol, dg, lambda1,lambda2,dc_filter)! dc_filter er dc_hat
          dconstraint(:,1) = 2 * (DOT_PRODUCT(L2,D)/DOT_PRODUCT(L,D)) * (dc_filter*DOT_PRODUCT(L,D) - &
               DOT_PRODUCT(L2,D)*dc)/DOT_PRODUCT(L,D)**2
          dobject = dc
       case(8)! min(max(compliance))
          do j=1,nk
             call gradient(flag,problem_type,inak,D_mat(:,j),rho,rho_min, dc,vol,max_vol, dg)
             dconstraint(:,j) = dc
          end do
       case(9)! mekanisme - Robust
          problem_type = 1
          do j = 1,3
             call gradient(flag,problem_type,inak,D_robust(:,j),rho_bar(:,j),rho_min,&
                  dc,vol,max_vol, dg,lambda_mat(:,j))
             dconstraint(:,j) = dc
          end do
          problem_type = 9
       case(10) ! elevator
          problem_type = 1
          do j=1,2
             call gradient(flag,problem_type,inak,D,rho,rho_min, dc,vol,max_vol, dg, lambda_mat(:,j))
             dconstraint(:,j) = dc
          end do
          do e=1,ne
             dc(e) = dconstraint(e,1) + dconstraint(e,2) + 2 * rpenal*(constraint(1)-constraint(2))&
                  *(dconstraint(e,1)-dconstraint(e,2))
          end do
          problem_type = 10
          dobject = dc
       case(12) ! eigenfrequency
          call compare_eigenval(eigenval(:,1),list)
          if(list(1) /= 0) then ! double eigenfrequency
             call gradient(flag,problem_type,inak,(/0d0/),rho,rho_min,&
                  dc,dc2= dconstraint(:,list(2)), &
                  eigenval=eigenval(list(1),1), eigenvec_d=eigenvec(:,list), double_eigen=.true.)
             dconstraint(:,list(1)) = dc
             
          end if
          na = ne_aktiv

          do j=1,nc-1
             if ((j == list(1) ).or. (j==list(2))) cycle

             call gradient(flag,problem_type,inak,(/0d0/),rho,rho_min, dc,&
                  eigenval=eigenval(j,1), eigenvec=eigenvec(:,j),double_eigen=.false.)

             if (scale) then
                dconstraint(1:na,j) = -1.0/constraint(j)*dc
                dconstraint(na+1,j) = 1.0/constraint(j)*1d0
             else
                dconstraint(1:na,j) = -dc
                dconstraint(na+1,j) = 1d0
             end if

          end do

          dconstraint(1:na,nc) = vol/max_vol !dg, skaleret volumen
          dconstraint(na+1,nc) = 0d0

          if (scale) then
             dobject = 0d0
             dobject(na+1) = -1.0/xval(na+1) !afledte mht beta
             constraint = log(constraint)
          else
             dobject = 0d0
             dobject(na+1) = -1d0 !afledte mht beta
          end if
          


          do j=1,nc-1
             !dconstraint(:,j) = dconstraint(:,j) 
          end do
       end select


       ! Finite difference check for sensitivity
       select case( fd_check )
       case( 1 )
          if (iter == 1) then
          ! kun for eigen
             call finite_check_eigen(problem_type,inak,filter_type,vol,max_vol, rho,rho_min, &
                  f_scale,dconstraint(1:na,1))

          end if

          ! if  (iter == 2) then!( MOD(iter,5) == 1 ) then
          !    select case (filter_type)
          !    case (4)! Robust
          !       problem_type = 1 ! kun for problem (9) med samme gradienter som (1)
          !       call finite_check(flag, problem_type,inak, filter_type,vol,max_vol,& 
          !            xval(1:n),rho_min,L,lambda_mat(:,3),dobject,df_approx,f_scale,beta,eta(3))
          !       problem_type = 9
          !    case default
          !       select case (problem_type)
          !       case(10)
          !          call finite_check(flag, problem_type,inak, filter_type,vol,&
          !               max_vol,xval(1:n),rho_min,L,lambda1,dobject, df_approx,&
          !               f_scale,beta,0d0, Lmulti=Lmulti)
          !       case default
          !          call finite_check(flag, problem_type,inak, filter_type,vol,max_vol,&
          !               xval(1:n),rho_min,L,lambda1,dobject, df_approx,f_scale,beta)
          !       end select
          !    end select
          !    call output_matrix(df_approx,'finite_diff_tjeck')
          ! end if
       end select


       select case( problem_type )! #8
       case default
          select case( filter_type )! #8
          case (1:2)
             call filter(filter_type,inak,vol,rho,dc,dg) 
          case (3)
             call filter(filter_type,inak,vol,rho_tilde,dc,dg,beta) ! #8
          end select
       case (7:9,11)
          select case( filter_type )! #8
          case (1:2)
             do j=1,size(dconstraint,2)-1
                dc_filter = dconstraint(:,j)
                call filter(filter_type,inak,vol,rho,dc_filter,dg)
                dconstraint(:,j) = dc_filter
             end do
          case (3:4)
             do j=1,size(dconstraint,2)-1
                dc_filter = dconstraint(:,j)
                call filter(filter_type,inak,vol,rho_tilde,dc_filter,dg,beta,eta(j)) ! #8
                dconstraint(:,j) = dc_filter
             end do
          end select
       end select


       ! Finding rho by bi-section or MMA method:
       select case( solver_type ) ! #9
       case(1)
          xval_old = xval
          call oc(filter_type,max_vol,vol,xval(1:n),xval_old(1:n),dc,dg)
       case(2)
          ! rho er de fysiske densiteter(dvs de der kommer fra densitetsfilteret).
          !rho_tilde = rho ! fysiske densiteter. Skal bruges til at beregne den aktuelle volumen i MMA-handle
          !rho = rho_mat ! MMA skal have den "oprindelige" rho ind og ikke rho_hat.
          ! da det er rho og ikke xval der filtreres, skal der ikke gøres noget ved xval


          compliance = compliance_out(i)/f_scale
          dobject = dobject/f_scale

          select case(problem_type)
          case default
             dobject = dc
             call mma_handle(iter,inak,rho_min,low,upp,max_vol,compliance,vol,xval, xval_old, xval_old2,dobject,dg,&
                  rho,1,movelimit)
          case(7,12)! mekanisme med begrænsning på krydsbevægelse
             call mma_handle(iter,inak,rho_min,low,upp,max_vol,compliance,vol,xval, xval_old, xval_old2, dobject, (/0d0/), &
                  rho,nc ,movelimit, constraint, dconstraint)
          case(8)! min(max(compliance))P = P1
             compliance = 0d0
             dc = 0d0
             dconstraint = dconstraint / f_scale
             constraint = constraint / f_scale ! skalering af compliance
             dconstraint(:,nc) = dg
             call mma_handle(iter,inak,rho_min,low,upp,max_vol,compliance,vol,xval, xval_old, xval_old2,dc,dg,&
                  rho,nc,movelimit,constraint,dconstraint) ! korrigeret, da g ,dconstraintensitetsfilteret
          case(9,11)! Robust design
             compliance = 0d0
             dc = 0d0

             dconstraint = dconstraint / f_scale
             constraint = constraint / f_scale ! skalering af compliance
             dconstraint(:,nc) = dg
             call mma_handle(iter,inak,rho_min,low,upp,max_vol,compliance,vol,xval, xval_old, xval_old2,dc,dg,&
                  rho_bar(:,1), nc, movelimit, constraint, dconstraint) ! korrigeret, da g ,dconstraintensitetsfilteret
          end select
       end select

       rho = xval(1:na)
       ! do e=1,ne
       !    if (rho(e) < 0) then
       !       print*,'efter mma'
       !       print*,'e',e
       !    end if
       ! end do

       change = maxval(abs(xval_old(1:na) - rho))

       select case (filter_type)
       case( 1:2 )
          ! STOP CRITERIA (På sædvanlig vis)
          select case(stop_krit)
          case(1)
             if ( ( change < tol*maxval(abs(rho)) ) .and.  (i > 40) ) then
                exit
             end if
          end select
       case( 3:4 )
          select case (filter_type)
          case(4) ! opdater volumen-begrï¿½nsning
             if (  MOD(iter,30) == 1 ) then
                max_vol = dot_product(rho_bar(:,1),vol) * max_vol_start/dot_product(rho_bar(:,2),vol)
             end if
          end select
          ! Opdater beta
          !if ( ( ( MOD(iter,50) == 0 )  .or. ( change < 0.01  ) ) .and. ( beta < 200 ) .and. (iter > 40) ) then ! overvej at bruge tolerance
          if ( ( ( MOD(iter,70) == 0 ) ) .and. ( beta < 200 ) ) then ! overvej at bruge tolerance
             beta = 1.3d0*beta
             change = 0.5 ! Lader iterationsprocessen lï¿½be lidt lï¿½ngere med ny beta.
          end if
          ! STOP CRITERIA (Med Heavyside)
          select case(stop_krit)
          case(1)
             if ( ( change < tol*maxval(abs(rho)) ) .and.  (i > 100) ) then
                exit
             end if
          end select
          if (iter == 1500) then
             exit
          end if

          !$$$$$$                 ! Opdater beta Yuriy-metode
          !$$$$$$                 if ( MOD(iter,30) == 0 ) then
          !$$$$$$                     if (iter < iter1) then
          !$$$$$$                         beta = beta1
          !$$$$$$                     elseif (iter > iter2) then
          !$$$$$$                         beta = beta2
          !$$$$$$                     else
          !$$$$$$                         beta = beta1 + (beta2-beta1)*((exp(zeta*i)-exp(zeta*iter1))/&
          !$$$$$$                                 (exp(zeta*iter2)-exp(zeta*iter1)))
          !$$$$$$                     end if
          !$$$$$$ !$$$$$$                 end if
          !$$$$$$                 ! STOP CRITERIA (heavyside + robust)
          !$$$$$$                 select case(stop_krit)
          !$$$$$$                     case(1)
          !$$$$$$                         if ( ( change < tol*maxval(abs(rho)) ) .and.  (i > iter2) ) then
          !$$$$$$                             exit
          !$$$$$$                         end if
          !$$$$$$                 end select
       end select


       select case( info )
       case(1)
          ! Printing relevant information:
          print*,'----------------------------------'
          print*,'compliance = ',compliance_out(i)
          print*,'mma_compliance = ',compliance
          print*,'iter = ',iter
          !print*,'Volumen = ',dot_product(rho,vol) ! volumen udregnes med "fysisk" densitet
          !print*,'change = ',maxval(abs(xval_old(1:na) - rho))
          print*,'tolerance = ',tol*maxval(abs(rho))
          print*,'max dc',maxval(dc)
          print*,'max dg',maxval(dg)
          print*,'beta',beta
       end select
       if (i == 10) then
          !$$$$$$           stop
       end if
    end do


    ! Computation of filteret densities rho_hat
    ! overvej om ikke denne skal flyttes efter MMA, sï¿½ledes kaldet i begyndelsen kan fjernes

    select case (filter_type) ! En sidste densitetsfiltrering: "matematisk" rho fra MMA -> "fysisk" rho
    case( 2 )
       call density_filter(inak,vol, rho)
    case ( 3 )
       call density_filter(inak,vol, rho)           
       do e=1,ne
          rho(e) = EXP(-beta*(1.0-rho(e))) - (1.0-rho(e))*EXP(-beta) ! Modificeret heavyside (29): rho -> rho_bar
       end do
    case(4) ! densitetsfilter med robustprojektion
       call density_filter(inak,vol, rho)
       do j=1,3
          do e=1,ne
             rho_bar(e,j) = (dtanh(beta*eta(j))+dtanh(beta*(rho(e)-eta(j))))/&
                  (dtanh(beta*eta(j))+dtanh(beta*(1-eta(j))))
          end do
       end do
       rho = rho_bar(:,1)
    end select

    ! Recover stress
    select case (problem_type)
    case(0:1) ! static
       call recover
    case(8)
       call output_matrix(compliance_mat,'compliance_mat')
    end select

    select case (filter_type)
    case(4)
       call output_matrix(compliance_mat,'compliance_mat')
       call output_matrix(rho_bar,'rho_mat')
    end select

    call output_vector(compliance_out,'objekt')	
    call output_elements('Densitet',rho)

    select case( save_rho)
    case(1) ! udskriv rho
       call output_vector(rho,trim(filename)//'_rho')
       ! Indlï¿½s rho med nedenstï¿½ende
       !$$$$$$     call rho_input(rho_vec,trim(filename)//'dir\'//trim(filename)//'_rho'//'.m')
    end select

    ! procent gråskala. Formel fra sigmund_filterartikel
    Mnd = 0.0d0
    do e=1,ne_aktiv
       Mnd =Mnd+ 4*rho(e)*(1d0 - rho(e))
    end do
    Mnd = Mnd/real(ne_aktiv)*100
    print*,'Mnd',Mnd

    ! udskriv parametre til fil
    parameters = 0d0
    parameters(1) = filter_type
    parameters(2) = solver_type
    parameters(3) = problem_type
    parameters(4) = vol_type
    parameters(5) = rmin_type
    parameters(6) = info
    parameters(7) = save_rho
    parameters(8) = n_iter

    parameters(9) = penal
    parameters(10) = damp_fact
    parameters(11) = max_vol
    parameters(12) = rmin
    parameters(13) = rho_min
    parameters(14) = tol
    parameters(15) = Mnd
    parameters(16) = total_vol
    parameters(17) = element_areal
    parameters(18) = compliance
    parameters(19) = i-1 ! antal iterationer
    parameters(20) = beta

    call output_vector(parameters,'parametre')

    ! Luk animationsvonduet ved at kalde med iter = -1
    select case( animation )
    case(1:2)
       call exodus_finalize

       call  plotanim(-1, 0, 1, .true., .false., .false., .true.,.true., 'anim', &
            0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1.0d0,rho)
    end select

    call output_deformed('undeformed','denstitet',rho)
    call output_deformed('deformed','densitet',rho)

  end subroutine topopt


  SUBROUTINE TopOpt_trans(flag,rho_in)

    use fea
    use transient
    use file_init
    use fedata
    use plane42transient
    use plot_routiner
    use processor
    use top_help

    integer, intent(inout) :: flag
    real(8), dimension(:),optional, intent(in) :: rho_in
    integer :: nmax, n_iter
    integer, dimension(:) , allocatable :: dof_x
    real(8) :: deltaT , fc, T_0, delta, alpha1, beta1
    !fc : center frequency, T_0 : center of wave packet in time domain, delta = bandwidth

    real(8),dimension(neqn) :: U0,dotU0, L
    real(8), dimension(:,:),allocatable :: saved_U, saved_lambda

    integer :: file, problem_type, fd_check, filter_type, info, save_rho
    integer :: animation, output_type, loes_type, parameters(10)
    real(8) :: parameters_real(8), rho_min

    real(8) :: lmin, cL ! til bestemmelse af deltaT
    real(8) :: objekt, tids_int, change, thk, deltaT_max
    integer :: e,i,n,iter, nen, idof

    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe, vector,vector2, dKU, dCU, dMU! Tidsintegration (7)
    real(8), dimension(mdim,mdim) :: dme,dke,dce

    real(8) :: df_approx(ne,2), value, dummy_mat(1,1) = 0.0d0 ! til kald af transient
    real(8), dimension(:), allocatable :: objekt_out

    ! Inaktive elementer
    type(INAKTIV_def), allocatable :: inak

    ! til fd_chech

    ! MMA 
    REAL(8), DIMENSION(ne):: low, upp, rho, rho_old, rho_old2, dg, dc
    ! bemï¿½rk at dc = d_xi

    real(8) :: f_scale

    allocate(inak)
    inak%bool = .false.

    ! ##############################
    ! Her starter du
    file = 1			! 1 : longsving2		2: disk
    output_type = 1		! 1 : udskriv forskydning for hvert tidsskridt, !	2 : udskriv forskydning for udvalgt DOF, !	3 : udskriv forskydning for alle centerknuder, ! 4 : udskriv energi i bjï¿½lken
    loes_type= 2		! 1 : Transient			2: TopOpt
    animation = 0 		! 1 : Til !				0: fra
    !$$$$$$     flag = 1            ! Benyttes kun ved transient. 0 : fra   1 : alle elementer er ens
    save_rho = 1		! 1 : gem rho

    filter_type = 0		! 0: No filter 			1: Sensivity filter. 2: Density
    fd_check   	= 0		! 0: Fra				1: Finite difference check
    info 		= 1		! information
    problem_type = 5	! 0: Fra				3: bandgap
    rand = 1			! 0: ingen dï¿½mpning		1: dï¿½mp hï¿½jre side	2: dï¿½mp venstre og hï¿½jre	4: dï¿½mp alle sider
    normal = 1 ! er bï¿½lgeretningen normal pï¿½ fladen?	0: nej	1: ja

    !call initial_transient
    if (present(rho_in)) then
       rho = rho_in
    else
       rho = 0.50d0
    end if
    call inputfile_trans(file,rand,normal,loes_type,output_type &
         ,rho,n_iter,nmax,deltaT,cL,young1,dens1,nu1,young2,dens2,nu2)

    ! Initialise external force
    call buildload

    ! knuder der skal minimeres
    L = 0d0
    if (nd> 0) then
       do i = 1, nd
          ! find dof og indsï¿½t i L
          idof = 2*(nodes(i, 1)-1)+nodes(i, 2)
          L(idof) = 1.0d0
       end do
    else
       L(neqn-1) = 1.0d0
    end if

    call elements_init(file,L) ! find elementer i enderne

    ! sï¿½tter inaktive elementer
    ! Nb udkommenteret. SKal slï¿½s til igen
    !$$$$$$     value = 0.5d0
    !$$$$$$     call inactive_elements(rho,value)


    ! MMA init
    low = 0.0
    upp = 1.0
    rho_old = rho
    rho_old2 = rho 
    dg = 0.0d0 ! ingen constrain ts => dg = 0

    U0 = 0.0 ! initial displacement
    dotU0 = 0.0
    dc = 0.0d0
    rho_min = 1E-4


    !$$$$$$     call plot_gauss(nmax,deltaT, fc, T_0, delta)      

    !deltaT = 0.050d0 ! til longsving5
    !$$$$$$     nmax = 500 ! til longsving5
    !nmax = 2000
    alpha1 = 0.00
    beta1 = 0.0000
    !$$$$$$     n_iter = 50

    parameters(1) = file
    parameters(2) = flag
    parameters(3) = animation
    parameters(4) = loes_type
    parameters(5) = output_type
    parameters(6) = nmax
    parameters(7) = rand
    parameters(8) = normal
    parameters(9) = fd_check
    parameters_real(1) = cL

    select case ( loes_type )
    case ( 2 ) ! initialiser matricer til TopOpt
       allocate( saved_U(neqn,nmax+1) )
       allocate( saved_lambda(neqn,nmax+1) )
       allocate( objekt_out(n_iter) )
       objekt_out = 0.0d0
    end select

    iter = 0
    change =1.0d0
    do while( (change > 1E-6) .AND. (iter < n_iter) )
       iter = iter +1
       print*,'iter',iter
       parameters(10) = iter

       select case ( loes_type )
       case ( 1 ) !Kun transient
          call half_step_CD(parameters,parameters_real,flag,deltaT,rho,U0,dotU0)
          EXIT
       end select

       !$$$$$$     if (iter == 2) then
       !$$$$$$       pause
       !$$$$$$     end if

       ! Solve the transient problem
       parameters(1) = file ! skal vï¿½re her, da den nulstilles senere i lï¿½kken
       call half_step_CD(parameters,parameters_real,flag,deltaT,rho,U0,dotU0,saved_U)
!!!!! Efter time-integration !!!!
       ! STEP 2
       ! Objekt-funktionen er integralet af "forskydning af endeknuden", dvs d(neqn-1)^2, over tiden
       objekt = 0.0d0 ! (7) i paper
       do n=1,nmax ! Objekt ligges sammen. Evt kan trapez-reglen benyttes
          if ( (n > 1) .and. (n < nmax) ) then 
             objekt = objekt + 2.0d0* dot_product( saved_U(:,n) ,L* saved_U(:,n))
          else
             objekt = objekt + dot_product( saved_U(:,n) , L* saved_U(:,n))
          end if

       end do
       objekt = (deltaT)/(2.0d0) * objekt
       ! Bemï¿½rk at trapezreglen hedder (b-a)/(2*N) * (f(x1)+2*f(x2)+...+2*f(x_N-1)+f(x_N))
       ! Her er (b-a)/(2N) = deltaT*real(nmax)/(2*real(nmax))=deltaT/2 
       objekt_out( iter ) = objekt ! opr. objekt til plot
       print*,'objekt',objekt
       print*,'d*d',dot_product(saved_U(:,1),saved_U(:,1))

       if (( MOD(iter,20) == 0 ) .OR. ( iter == 1 )) then
          f_scale = dabs(objekt)
       end if


       ! Jf ovenstï¿½ende er G = d(neqn-1)^2, dvs dG/du = 2*d(neqn-1).
       ! Pga variabel-transformationen tau = T(slut-tid)-t, integres fra tau = T gï¿½ende mod tau = 0

       ! Prï¿½v at flytte nedenstï¿½ende ind i objekt-lï¿½kken. Mï¿½ske er det hurtigere.
       ! Dette er dG. Af "plads-hensyn" kaldes den saved_lambda
       saved_lambda = 0.0d0
       do n=0,nmax-1 ! bemï¿½rk "kï¿½rer" baglï¿½ns
          !$$$$$$       do i=1,size(dof_x,1)
          !$$$$$$         saved_lambda(dof_x(i),n+1) = (2.0d0*saved_U(dof_x(i),nmax-n)) ! Objektfunktion forholder sig kun til sidste knude.
          saved_lambda(:,n+1) = (2.0d0* L *saved_U(:,nmax-n)) ! Objektfunktion forholder sig kun til sidste knude.
          !print*,'Gmat',dGmat(dof_x,n)
          !$$$$$$       end do
       end do

       ! Begyndelsesbetingelser for lambda.
       U0 = 0.0d0
       dotU0 = 0.0d0
       parameters(1) = 0 !=> file = 0, da lasten er giver ved dGmat
       call half_step_CD(parameters,parameters_real,flag,deltaT,rho,U0,dotU0,saved_lambda ) ! STEP 3

       dMU = 0.0d0
       dCU = 0.0d0
       dKU = 0.0d0
       dc = 0.0d0
       do e=1,ne
          ! Find coordinates and degrees of freedom
          nen = element(e)%numnode
          do i = 1, nen
             xe(2*i-1) = x(element(e)%ix(i),1)
             xe(2*i  ) = x(element(e)%ix(i),2)
             edof(2*i-1) = 2 * element(e)%ix(i) - 1  
             edof(2*i)   = 2 * element(e)%ix(i)
          end do

          if (e == 1 .and. iter == 1) then
             thk = mprop(element(e)%mat)%thk
             call plane42transient_dke(xe, young1,young2, nu1, nu2, thk,ng, dke)! calculate element derivative of stifness-matrix
             call plane42transient_dme(xe,thk,dens1, dens2, ng, dme)! calculate element derivative of mass-matrix
             print*,'tidsint kald dme'
          end if

          tids_int = 0.0d0
          do n=1,nmax! STEP 4 - Integration (17)

             vector = saved_U(edof,n)
             dKU = MATMUL(dke,vector)

             if ( n == 1 ) then
            	!vector2 = ( saved_U(edof,n+1)  - U0(edof) ) / ( 2.0d0*deltaT )! (11.12-1a)
                vector = ( saved_U(edof,n+1)  - 2.0d0*saved_U(edof,n)+ U0(edof) )  / deltaT**2 ! (11.12-1b)
             elseif ( n > 1 ) then
            	!vector2 = ( saved_U(edof,n+1)  - saved_U(edof,n-1) ) / ( 2.0d0*deltaT )! (11.12-1a)
                vector = ( saved_U(edof,n+1)  - 2.0d0*saved_U(edof,n)+ saved_U(edof,n-1) )  / deltaT**2 ! (11.12-1b)
             end if
             dMU = MATMUL(dme,vector)

             if (alpha1>0 .or. beta1>0) then
                dce = alpha1*dme+beta1*dke
                dCU = MATMUL(dce,vector2)
             end if

             do i=1,8
              	vector2(i) = ( -dMU(i) -dCU(i) -dKU(i) )
             end do
             ! Husk at lambda er fundet ved "bagvendt" integration. Derfor skal den "vendes" om igen her.
             vector = saved_lambda(edof,nmax+1-n)
             tids_int = DOT_PRODUCT(vector , vector2 )

             if ( (n > 1) .and. (n < nmax) ) then ! Trapez-integration
                dc(e) = dc(e) + 2.0d0*tids_int
             else
                dc(e) = dc(e) + tids_int
             end if
          end do
          dc(e) = (deltaT)/(2.0d0) * dc(e)
          dc(e) = dc(e)/f_scale ! skalering
       end do
       print*,'dc*dc',DOT_PRODUCT(dc,dc)

       select case( fd_check )
       case (1)
          if  (iter == 2) then!( MOD(iter,5) == 1 ) then
             parameters(1) = file
             !call finite_check(flag,problem_type, filter_type,(/0.0d0/),0d0,rho,rho_min,L,(/0.0d0/),& 
              !    dc, df_approx,f_scale,0d0,0d0,parameters,deltaT,U0,dotU0,dof_x)
             call output_matrix(df_approx,'finite_diff_tjeck')
          end if
       end select

       objekt=objekt/f_scale ! skaleret objekt til mma
       dc = dc/f_scale
       call mma_handle(iter,inak,rho_min,low,upp,0.0d0,objekt,(/0.0d0/),rho,rho_old,rho_old2,dc,dg,(/0.0d0/),1) ! STEP 5
       call inactive_elements(rho,value)
       call  plotanim(iter, 0, 1, .true., .false., .false., .true.,.true., 'Udboej', &
               0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1.0d0,rho)

       change = MAXVAL(DABS( rho - rho_old )) ! # 10
       print*,'change in rho',change

    end do

    select case ( loes_type )
    case(2)
       call output_vector(objekt_out,'objekt')
       select case( save_rho)
       case(1) ! udskriv rho
          call output_vector(rho,trim(filename)//'_rho')
          ! Indlï¿½s rho med nedenstï¿½ende
          !$$$$$$     call rho_input(rho_vec,trim(filename)//'dir\'//trim(filename)//'_rho'//'.m')
       end select
    end select

    call output_elements('Densitet',rho)



  END SUBROUTINE TopOpt_trans

  subroutine inactive_elements(rho,value)
    use fedata

    real(8), intent(in) :: value
    real(8), dimension(:), intent(INOUT) :: rho
    integer :: e1, e2,i,e, nen, kk, n

    select case(rand)
    case(0)
       return
    case default
       e1 = 0
       e2 = 0
       kk = 0
       n = 2 ! "lag" af inaktive elementer = n+1, dvs n=2 giver 3 inaktive elementer
       !n = 39 ! til stï¿½nger, dvs aflevering specialkursus
       !n = 59! til plots
       do e=1,ne
          do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
             if (e == element_rand(i,1) ) then

                select case(element_rand(i,2))!face
                case(1) ! face 1(bottom of structure)
                   rho(e:e+n) = value
                case(2) ! hï¿½jre side
                   if (e1 ==0) then ! finder fï¿½rste element i hï¿½re side. Dvs det nederste
                      e1 = e
                   end if
                   kk = kk+1! tï¿½ller antal elementer i hï¿½jden
                case(3) ! top
                   rho(e-n:e) = value
                case(4) ! finder sidste element i venstre side. Dvs det ï¿½verste
                   if (e>e2) then
                      e2 = e
                   end if
                end select
             end if
          end do
       end do

       if (e1>0) then
          rho(e1-n*kk:e1+kk) = value
       end if
       if (e2>0) then
          rho(1:(n+1)*e2) = value
       end if
    end select


  end subroutine inactive_elements

  subroutine elements_equal(flag)

    use fedata
    use fea    

    integer, intent(OUT) :: flag

    integer, parameter :: mdim = 8
    integer :: e,i,nen
    real(8), dimension(mdim) :: xe
    real(8), dimension(ne,2) :: l

    do e =1,ne
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
       end do
       l(e,1) = xe(5) - xe(1) ! Element side length(x)
       l(e,2) = xe(6) - xe(2) ! Element side length(y)
       !$$$$$$         print*,'l(e,1)',l(e,1)
       !$$$$$$         print*,'l(e,2)',l(e,2)
    end do

    ! Test if the elements are equal:
    if ( (maxval(l(:,1)) .eq. minval(l(:,1))) .and. (maxval(l(:,2)) .eq. minval(l(:,2))) ) then
       flag = 1
    else
       flag = 0
    end if
    !$$$$$$      print*,'minval1',minval(l(:,1))
    !$$$$$$      print*,'minval2',minval(l(:,2))
    !$$$$$$ 
    !$$$$$$      print*,'maxval1',maxval(l(:,1))
    !$$$$$$      print*,'maxval2',maxval(l(:,2))

  end subroutine elements_equal


END MODULE topology

