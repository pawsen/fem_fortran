module transient

  implicit none

  private
  public :: half_step_CD, buildmass, mmul_init, mmul, elements_init, deltaT_init, explosive_source
  public :: abs_init,abs_mmul, initial_transient

contains

  subroutine initial_transient(flag)

    ! This subroutine is mainly used to allocate vectors and matrices
    use fedata
    use link1
    use plane42

    integer, intent(IN) :: flag
    integer :: n, e, nen, min_edof, max_edof, diff_edof

    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof

    ! This subroutine computes the number of global equation,
    ! half bandwidth, etc and allocates global arrays.

    ! Calculate number of equations
    bw=0 ! giver bw en vï¿½rdi til if-sï¿½tningen
    neqn = 2*nn
    select case(flag)
    case(0)  
       if (banded == 0) then
          allocate (k(neqn, neqn))
       elseif(banded == 1)then !allocate k=banded
          do e=1,ne
             nen = element(e)%numnode ! Antallet af knudepunkter i ï¿½t element
             do n = 1, nen !nen = number of element nodes. Just for Yuriy
                edof(2*n-1) = 2 * element(e)%ix(n) - 1  !dof for element e, knude n. i x-retning
                edof(2*n)   = 2 * element(e)%ix(n)  !dof for element e, knude n. i x-retning
             end do
             min_edof = minval(edof(1:2*nen)) ! For truss, er edof 4 lang(nen=2). For kontinuum er edof 8(nen=4). Derfor benyttes edof(1:2*nen)
             max_edof = maxval(edof(1:2*nen))
             diff_edof =max_edof-min_edof+1
             if (diff_edof > bw) then
                bw = (max_edof-min_edof+1)
             end if
          end do

          allocate (k(bw,neqn))
          print*,'bw',bw
       end if
    end select

    allocate (p(neqn), d(neqn))
    allocate (strain(ne, 3), stress(ne, 3))


    ! MODIFIED FOR WINDOWS COMPILER
    do e=1,ne
       strain(e,1:3) = 0
       stress(e,1:3) = 0
    enddo

  end subroutine initial_transient


  subroutine half_step_CD(parameters,parameters_real,flag,deltaT,rho,U0,dotU0,saved_U)
    ! This subroutine does direct integration by explicit half step method

    use plot_routiner
    use processor
    use fedata
    use fea
    use WriteEnsight

    real(8), dimension(:,:),optional, INTENT(INOUT) :: saved_U
    real(8), dimension(neqn), INTENT(IN) :: U0,dotU0
    real(8), dimension(ne), INTENT(IN) :: rho
    real(8), INTENT(IN) :: deltaT
    integer, dimension(:) , intent(IN):: parameters ! indeholder: nmax, file, flag,animation,loes_type,output_type
    real(8), dimension(:),intent(IN) :: parameters_real
    integer, intent(IN) :: flag

    integer :: i, n,e, kk, file,animation,loes_type,output_type, nmax, nen, dev, fd_check, iter
    real(8), dimension(neqn) ::  Dprev, Dprev2, Dcurr, vec_rhs, mlumped, r_int, MDcurr, CDcurr, MDprev, CDprev, r_ext
    real(8), dimension(neqn) :: dotD, ddotD, r_int_Z, M
    real(8), dimension(parameters(6)) :: output_vec, energy2,output_kin,output_pot ! dim = nmax
    real(8), dimension(parameters(6)) :: output_d,output_dotD ! dim = nmax
    real(8), dimension(parameters(6)+1) :: kraftvec
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8) :: time, t_shift
    real(8) , allocatable, dimension(:,:) ::  ke1 , ke2 ,me
    real(8) :: plotval(ne), radius, center(2), t0, f0, ft, fak, cL, thk
    real(8), allocatable, dimension(:) :: me_lumped, g! stedvektor til eksplosiv kilde
    character(len=40):: buffer, buffer2


    mlumped = 0.0
    output_vec = 0.0

    file = parameters(1)
    !flag = parameters(2)
    animation = parameters(3)
    loes_type = parameters(4)
    output_type = parameters(5)
    nmax = parameters(6)   
    cL = parameters_real(1)

    kk = 0 ! Til plot af forskydning

    ddotD = 0 !initail acc. eq. 11.12-12
    dotD = dotU0 !initial velocity
    Dcurr = U0 ! initial displacement

    Dprev = Dcurr-deltaT*dotD+(deltaT**2/2.0)*ddotD
    CDcurr = 0.0d0
    CDprev = 0.0d0

    r_ext = 0.0d0
    select case(flag)
    case(1)!alle elementer er ens. Derfor bruge samme ke,me, me_lumped
       allocate( ke1(8,8), ke2(8,8), me(8,8), me_lumped(8) )
       call mmul_init(rho,ke1,ke2,me,me_lumped)
       call buildmass(flag,rho,M) ! bruges til beregning af energi
    case(0)
       call buildstiff_fea(flag,rho)
       call buildmass(flag,rho,M)! 0 fortï¿½ller at elementerne ikke er ens
    end select

    ! HUSK RB pï¿½ M!!!!!!!! VIGTIGT !!!!!!
    select case(file)
    case(2)
       !call explosive_source(radius,center,g)
       center(1) = 15.0d0
       center(2) = 15.0d0
       thk = mprop(element(1)%mat)%thk
       radius = 5d0*thk
       f0 = cL/(thk*10)
       t0 = 1.0d0/f0
       print*,'t0',t0
       if (.not. allocated(g)) then
          allocate(g(neqn))
       end if
       call explosive_source(Dcurr,radius,center,g)
       call abs_init(Dcurr,center)
       do n=1,size(g,1)
          if ( g(n) /= 0) then
             kk = kk +1
          end if
       end do
       print*,'antal DOF i g',kk
    case(3:4)
       !call explosive_source(radius,center,g)
       center(1) = 15.0d0
       center(2) = 15.0d0
       radius = 6*mprop(element(1)%mat)%thk
       if (.not. allocated(g)) then
          allocate(g(neqn))
       end if
       call explosive_source(Dcurr,radius,center,g)
       do n=1,size(g,1)
          if ( g(n) /= 0) then
             kk = kk +1
          end if
       end do
       print*,'antal DOF i g',kk
    end select

    kk = 0
    do n=1,nmax+1
       time = deltaT*real(n)
       select case(file)
       case(1)!Jacobs fil
          !$$$$$$             ft = DSIN(0.5*time)
          ft = DEXP(-0.005d0*(time-25d0)*(time-25d0))*DSIN(0.5*time) ! Jakob kraft
          !ft = DEXP(-0.005d0*(time*7-25d0)**2)*DSIN(0.5*time*7) ! Jakob kraft
          !ft = DEXP(-0.5d0*(time*7-25d0)**2)*DSIN(0.5*time) ! Jakob kraft
          r_ext = 0.010d0*p(1:neqn)*ft
          kraftvec(n) = ft
       case(2)! disk med eksploderende kilde
          if (time <= 2*t0) then
             !call explosive_source(Dcurr,radius,center,g)
             ft = 0.10d0*(-2)*pi**2*f0**2*(time-t0)*DEXP(-pi**2*f0**2*(time-t0)**2)
             r_ext =ft*g ! Eksploderende kilde
          else
             ft = 0.0d0
             r_ext = ft
             !print*,'time',time
          end if
          !$$$$$$             fak = 10.0d0
          !$$$$$$             if (time < fak*deltaT*30*2.5) then
          !$$$$$$             ft = DEXP(-0.005d0*(time*fak-25d0)**2)*DSIN(0.5*time*fak) ! Jakob kraft
          !$$$$$$             else
          !$$$$$$               ft = 0
          !$$$$$$             end if
          !$$$$$$             r_ext = ft*g

          kraftvec(n) = ft
          !$$$$$$             call abs_init(Dcurr,center)! beregner enhedsvektor og radius r 
       case(3)! disk med gauss-profil
          fak = 25.0d0
          ft = DEXP(-0.005d0*(time*fak-25d0)**2)*DSIN(0.5*time*fak) ! Jakob kraft
          r_ext =0.1d0*ft*g

          kraftvec(n) = ft
       case(4)! disk med harmonisk kraft
          ft = DSIN(1*2*pi*time)
          r_ext =0.1d0*ft*g

          kraftvec(n) = ft
       case(0)! kraft er specificeret i matrix
          if ( n /= nmax+1 ) then
             r_ext = saved_U(:,n)
          end if
       case(5) ! gauss-puls til at generere struktur
          t_shift = 10d0
          T0 = 5.9d0
          f0 = 0.5d0
          ft = sin(2*pi*f0*(time-t_shift))*dexp(-(time-t_shift)**2/T0**2)
          r_ext = 1.0d0*p(1:neqn)*ft
          kraftvec(n) = ft
       case(6) ! smal gauss-puls til at lave transmisstans plot
          t_shift = 3d0
          !$$$$$$             T0 = 1.3d0 ! virker nogenlunde
          T0 = 0.3d0
          f0 = 0.5d0
          ft = dsin(2*pi*f0*(time-t_shift))*dexp(-(time-t_shift)**2/T0**2)
          !$$$$$$             ft = dcos(2*pi*f0*(time-t_shift))*dexp(-(time-t_shift)**2/T0**2)
          r_ext = 10.0d0*p(1:neqn)*ft
          kraftvec(n) = ft
       end select
       !$$$$$$     if (n==nmax) then
       !$$$$$$       pause
       !$$$$$$     end if


       select case(flag)
       case(1)! Alle elementer er ens
          call mmul(flag,Dcurr,r_int,1,rho,(/0.0d0/),ke1,ke2) ! internal forces
          !$$$$$$         call mmul(flag,Dcurr,MDcurr,3,rho,me_lumped) 
          !$$$$$$         call mmul(flag,Dprev,MDprev,3,rho,me_lumped)
          MDcurr =  Dcurr*M
          MDprev = Dprev*M
       case(0)
          call band_mul(r_int,K,Dcurr) ! internal forces
          MDcurr =  Dcurr*M
          MDprev = Dprev*M
          ! Husk RB IKKE er "monteret"
       end select

       MDcurr = (2.0d0/deltaT**2)*MDcurr
       MDprev = (1.0d0/deltaT**2)*MDprev

       select case(rand) 
       case(0) ! inegn dï¿½mpning
       case default
          select case(normal)
          case(0) ! punktkilde, dvs ekstra stivhed ligges til
             !$$$$$$                 call abs_init(Dcurr,center)
             call abs_mmul(0,0,rho,Dcurr,CDcurr,r_int_Z) ! dï¿½mpning fra abs. rand
             call abs_mmul(0,0,rho,Dcurr,CDcurr) ! dï¿½mpning fra abs. rand
             CDcurr = CDcurr/deltaT
             call abs_mmul(0,0,rho,Dprev,CDprev)	
             CDprev = CDprev/deltaT

             r_int = r_int + r_int_Z
          case(1)
             call abs_mmul(0,1,rho,Dcurr,CDcurr) ! dï¿½mpning fra abs. rand
             CDcurr = CDcurr/deltaT
             call abs_mmul(0,1,rho,Dprev,CDprev)	
             CDprev = CDprev/deltaT
          case(2) ! punktkilde, dvs ekstra stivhed ligges til
             !$$$$$$                 call abs_init(Dcurr,center)
             call abs_mmul(0,0,rho,Dcurr,CDcurr) ! dï¿½mpning fra abs. rand
             call abs_mmul(0,0,rho,Dcurr,CDcurr) ! dï¿½mpning fra abs. rand
             CDcurr = CDcurr/deltaT
             call abs_mmul(0,0,rho,Dprev,CDprev)	
             CDprev = CDprev/deltaT
          end select
       end select

       vec_rhs = r_ext - r_int + MDcurr-CDcurr - MDprev + CDprev ! rhs of (11.12-8)
       mlumped = M/(deltaT**2)

       Dprev2 = Dprev
       Dprev = Dcurr

       Dcurr = vec_rhs/mlumped   
       dotD = ( Dcurr - Dprev2 ) / ( 2.0d0*deltaT )! (11.12-1a)

       ! beregn plot-vï¿½rdi
       D = Dcurr ! plotanim skal have forskydninger for D
       do e = 1, ne
          nen = element(e)%numnode
          do i = 1, nen
             edof(2*i-1) = 2 * element(e)%ix(i) - 1  
             edof(2*i)   = 2 * element(e)%ix(i)
          end do
          plotval(e) = dsqrt(dot_product(Dcurr(edof),Dcurr(edof)))
          !plotval(e) = 0!(stress(e,1)**2+stress(e,2)**2-stress(e,1)*stress(e,2)+3.0*stress(e,3)**2)**(0.50d0) !von mises spï¿½nding
       end do

       select case ( animation )
       case ( 1 )
          D = Dcurr ! plotanim skal have forskydninger for D
          do e = 1, ne
             nen = element(e)%numnode
             do i = 1, nen
                edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                edof(2*i)   = 2 * element(e)%ix(i)
             end do
             plotval(e) = dsqrt(dot_product(Dcurr(edof),Dcurr(edof)))
             !plotval(e) = 0!(stress(e,1)**2+stress(e,2)**2-stress(e,1)*stress(e,2)+3.0*stress(e,3)**2)**(0.50d0) !von mises spï¿½nding
          end do
          call  plotanim(n, 0, 1, .true., .true., .false., .true.,.true., 'Udboej', &
              0.0d0, (/0.0d0/), (/0.0d0/), (/0.0d0/),0.0d0,1d0,rho)
       end select

       select case ( loes_type )
       case ( 1 ) ! Transient lï¿½sning
          select case ( output_type )
          case ( 1 )
             dev = 5
             if (mod(n,dev) == 0) then ! if the remainder after deviding by 10 is zero, then call output. Eg. output kaldes ved hvert 10'ende n.
                kk = kk+1
                if (kk == 1) then;call write_structure('_plotdeformed_topology.m.m'); call output_data('data_anim',rho,deltaT,dev)
                end if
                   !call output_anim(Dcurr,kk) !,plotval) udskriver forskydning for hvert skridt


                   WRITE (buffer2,"('000000',TL4,i4)") kk
                   write (buffer, '(a,a,a)') trim(filename),'_',buffer2
                   print*,buffer
                   call WriteEnsightVar(plotval,buffer)
                end if
                if (n>1) then !udskriv energien i strukturen
                   output_vec(n-1) = 0.50d0*DOT_PRODUCT(dotD,M*dotD) + 0.50d0*DOT_PRODUCT(Dprev,r_int)
                end if
             case ( 2 )
                ! D(Midter_dof_x) el evt neqn-1
                !call output_1DOF(D(neqn-1),n,nmax,deltaT,'DOFdisp') ! Udskriver forskydning for bestemt knude
                ! knuden er givet ved: "F,	"node",	hflu,	1,	1" i input.
                do i = 1, np !np = number loads
                   if (loads(i, 1) == 5) then
                      Midter_dof_x = 2*(loads(i,2)-1) + 1
                   end if
                end do
                !$$$$$$                     print*,'Midter_dof_x',Midter_dof_x
                output_vec(n) = Dcurr(Midter_dof_x)
                !print*,'midterknude-endelig',midterknude
             case ( 3 ) ! udskriver forskydning for alle center DOF i x-retning
                call output_center(Dcurr(center_dofs_x),n,nmax,deltaT,rho)
             case ( 4 ) ! udskriver "energien" i bjï¿½lken
                ! HUSK r_int = K*Dprev, Derfor skal det ogsï¿½ vï¿½re Dprev herunder
                if (n>1) then
                   output_vec(n-1) = 0.50d0*DOT_PRODUCT(dotD,M*dotD) + 0.50d0*DOT_PRODUCT(Dprev,r_int) !
                   output_kin(n-1) = 0.50d0*DOT_PRODUCT(dotD,M*dotD) !
                   output_pot(n-1) = 0.50d0*DOT_PRODUCT(Dprev,r_int) !
                   output_d(n-1) = DOT_PRODUCT(Dprev,Dprev) !
                   output_dotD(n-1) = DOT_PRODUCT(dotD,dotD) !
                end if
                ! output skives efter tidsintegrationen - Hurtigere.
                ! call output_1DOF(energy,n,nmax,deltaT,'Energi') ! Udskriver forskydning for bestemt knude
             end select
          end select

          if (present(saved_U)) then
             saved_U(:,n) = Dcurr
             ! Da dotD og ddotD ikke lï¿½ngere gemmes og sendes tilbage til topology er nedenstï¿½ende ikke lï¿½ngere nï¿½dvendigt og derfor udkommenteret

             ! dotD og ddotD beregnes begge pï¿½ baggrund af bl.a. D_n+1. Derfor kï¿½res lï¿½kken nmax+1 gange og derfor fï¿½lgende:
             ! bemï¿½rk at for n haves: D_n, dotD_n-1 og ddotD_n-1. Derfor gemmes n=nmax+1 ikke for D.

             !$$$$$$         if (n == nmax +1) then
             !$$$$$$         ! dette er D_n+1.
             !$$$$$$             dotD = ( Dcurr  - saved_U(:,n-2) ) / ( 2.0d0*deltaT )! (11.12-1a)
             !$$$$$$             ddotD = ( Dcurr  - 2.0d0*saved_U(:,n-1)+ saved_U(:,n-2) )  / deltaT**2 ! (11.12-1b)
             !$$$$$$ 
             !$$$$$$         else
             !$$$$$$             
             !$$$$$$             if ( n == 2 ) then
             !$$$$$$                 dotD = ( saved_U(:,n)  - U0 ) / ( 2.0d0*deltaT )! (11.12-1a)
             !$$$$$$                 ddotD = ( saved_U(:,n)  - 2.0d0*saved_U(:,n-1)+ U0 )  / deltaT**2 ! (11.12-1b)
             !$$$$$$             elseif ( n > 2 ) then
             !$$$$$$                 dotD = ( saved_U(:,n)  - saved_U(:,n-2) ) / ( 2.0d0*deltaT )! (11.12-1a)
             !$$$$$$                ! print*,'n ',n
             !$$$$$$                 ddotD = ( saved_U(:,n)  - 2.0d0*saved_U(:,n-1)+ saved_U(:,n-2) )  / deltaT**2 ! (11.12-1b)
             !$$$$$$             end if
             !$$$$$$         end if
             !$$$$$$ 
             !$$$$$$         if (n > 1) then
             !$$$$$$             !if (present(saved_dotU)) then; saved_dotU(:,n-1) = dotD; endif
             !$$$$$$             if (present(saved_ddotU)) then; saved_ddotU(:,n-1) = ddotD; endif
             !$$$$$$         end if

          end if

       end do 



       select case ( loes_type )
       case ( 1 ) ! Transient lï¿½sning
          select case ( output_type )
          case(1) ! udskriver deformation yil hvert tilskridt
             call WriteEnsightGeo(filename)
             call EnsightCase('VonMises',1,kk,1)
             call output_vector(kraftvec,'kraftvec')
             call output_vector(output_vec,'energy')
          case( 2)
             call output_vector(kraftvec,'kraftvec')
             call output_vector(output_vec,'DOFdisp')
             call output_data('data',rho,deltaT)
          case( 4) ! udskriver "energien" i bjï¿½lken
             call output_vector(output_vec,'energy')
             call output_vector(output_kin,'kin_energy')
             call output_vector(output_pot,'pot_energy')
             call output_vector(output_d,'d')
             call output_vector(output_dotd,'dotd')
             call output_data('data',rho,deltaT)
             call output_vector(kraftvec,'kraftvec')
          end select
       case(2) ! top-opt
          fd_check = parameters(9)
          iter = parameters(10)
          if ((file/=0) .and. (iter == 1) .and. (fd_check == 0) ) then
             call output_vector(kraftvec,'kraftvec')
          end if
       end select


     end subroutine half_step_CD


     subroutine buildmass(flag,rho,M)
       ! This subroutine calculate the global lumped mass-matrix (diagonal),
       ! eq. its in vector format, from element masses

       use fedata
       use plane42transient

       real(8), dimension(neqn), INTENT(OUT) :: M
       real(8), dimension(:), INTENT(in) :: rho
       integer, intent(IN) :: flag
       integer :: e, i,j, nen, idof
       integer, parameter :: mdim = 8
       integer, dimension(mdim) :: edof
       real(8), dimension(mdim) :: xe, me_lumped
       real(8), dimension(mdim, mdim) :: me, me1
       real(8) ::thk

       M = 0.0d0

       do e = 1, ne

          ! Find coordinates and degrees of freedom
          nen = element(e)%numnode
          do i = 1, nen
             xe(2*i-1) = x(element(e)%ix(i),1)
             xe(2*i  ) = x(element(e)%ix(i),2)
             edof(2*i-1) = 2 * element(e)%ix(i) - 1  
             edof(2*i)   = 2 * element(e)%ix(i)
          end do

          if ((flag ==0) .or. (e==1)) then
             thk = mprop(element(e)%mat)%thk
             call plane42transient_me(xe,thk, ng, me1)
             me = ( rho(e)*dens1+(1.0d0-rho(e))*dens2 ) * me1
             call lumped_mass(lumped_type,me,me_lumped)
          else ! elementerne er ens, men pga. forskelle i rho skal me beregnes pï¿½ ny hver gang
             me = ( rho(e)*dens1+(1.0d0-rho(e))*dens2 ) * me1
             call lumped_mass(lumped_type,me,me_lumped)
          end if

          M(edof) = M(edof)+me_lumped
       end do

       do i = 1, nb
          idof = 2*(bound(i,1)-1) + bound(i,2)
          M(idof) = 1 ! Her skal diagonalelementer ved RB sï¿½ttes lig 1
       end do

     end subroutine buildmass

     subroutine band_mul(vector,A,u)
       ! inspiration fra http://maldun.lima-city.de/introduction_to_python/Example.html

       ! This subroutine multiply a symmetric bandmatrix with a vector
       ! A = bandmatrix, u = input-vector, vector = result

       use fedata


       real(8), dimension(neqn), INTENT(OUT) :: vector
       real(8), INTENT(IN) :: A(:,:), u(:)
       integer :: i,j

       bw = size(k,1)
       do i =1,neqn
          vector(i) = A(1,i)*u(i);
       end do

       do j = 2,bw
          do i =1,neqn
             if ((i+j-1) <= neqn) then
                vector(i) = vector(i) + A(j,i)*u(i+j-1)
                vector(i+j-1) = vector(i+j-1) + A(j,i)*u(i)
             end if
          end do
       end do

     end subroutine band_mul


     subroutine mmul_init(rho,ke1,ke2,me,me_lumped)
       ! This subroutine calculate element stiffnes and mass matrix and me_lumped vector,
       ! which is the mass matric lumped. Only to be used when all elements are equal.

       use fedata
       use plane42transient

       integer, parameter :: mdim = 8
       real(8), dimension(mdim,mdim), INTENT(OUT) :: ke1,ke2, me
       real(8), dimension(mdim), INTENT(OUT) :: me_lumped
       real(8), dimension(:), intent(IN) :: rho
       integer :: e, i ,j, nen
       integer, dimension(mdim) :: edof
       real(8), dimension(mdim) :: xe
       !real(8), dimension(mdim,mdim) :: ke1,ke2
       real(8) :: me_diag_sum, me_tot, thk


       ! alle elementer er ens => brug samme elementmatrice
       e = 1
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
       end do
       thk = mprop(element(e)%mat)%thk
       call plane42transient_ke(xe, young1,young2, nu1, nu2, thk,ng,ke1,ke2)
       call plane42transient_me(xe,thk, ng, me)

       call lumped_mass(lumped_type,me,me_lumped)

     end subroutine mmul_init

     subroutine lumped_mass(lumped_type,me,me_lumped)
       ! This subroutine lump element mass into a vector.
       ! Two types of lumping.

       integer, parameter :: mdim = 8
       integer, INTENT(IN) :: lumped_type
       real(8), dimension(mdim,mdim), INTENT(IN) ::  me
       real(8), dimension(mdim), INTENT(OUT) :: me_lumped
       integer :: i ,j, nen
       real(8) :: me_diag_sum, me_tot

       me_lumped = 0.0
       me_diag_sum  = 0.0
       me_tot = 0.0
       nen = 4 !number of nodes
       select case( lumped_type)
       case(1)! Yuriy lumping, Add all row element to corressponding diagonal
          do j=1,2*nen
             do i=1,2*nen
                me_lumped(j) = me_lumped(j)+ me(j,i)
             end do
          end do
       case(2)! HRZ - mass lumping, in COOK, s. 380
          do i=1,2*nen
             me_lumped(i) = me(i,i)
          end do
          me_tot = sum(me)! m in cook
          me_diag_sum = sum(me_lumped) ! S in COOK, s. 380
          do i = 1, 2*nen
             me_lumped(i) = (me_tot/me_diag_sum)*me_lumped(i)
          end do
       end select

     end subroutine lumped_mass


     subroutine mmul(flag,invector,outvector,mtype,rho,me_lumped_init,ke1,ke2)
       ! This subroutine multiply element matrix with the inputvector and correct it for bc.

       ! Only to be used when all elements are equal. If this is not the case, its faster
       ! to assemble the global matrices. Or maybe not with many dof.

       ! element_init: is either element stiffnes- og mass matrix. Only used when all
       ! elements are equal


       use fedata
       use plane42transient

       real(8), dimension(neqn), INTENT(IN) :: invector
       real(8), dimension(ne), INTENT(IN) :: rho
       integer, INTENT(IN) :: mtype, flag
       real(8), dimension(neqn), INTENT(OUT) :: outvector
       real(8) ,optional, INTENT(IN) :: ke1(:,:),ke2(:,:), me_lumped_init(:)
       integer :: e, i ,j,ii, eface
       integer :: nen, idof
       integer, parameter :: mdim = 8
       integer, dimension(mdim) :: edof
       real(8), dimension(mdim) :: xe , me_lumped
       real(8), dimension(mdim, mdim) :: ke, me,me1, ce, CZe, KRe
       real(8) :: young, nu, area, dens, thk, c_damp, length, x_koor, y_koor, p1, p2
       real(8) :: me_diag_sum, me_tot, lambda, mu, b, lambda2, mu2, lambda1, mu1, height
       real(8) :: parameters(13),r, r1,r2, nvec(2)

       ! Lame parametre angives i mmul 4.
       ! alpha1, beta1 angives i topology

       outvector = 0.0
       c_damp = 0.0

       select case( flag)
       case(1) ! alle elementer er ens
          do e = 1,ne
             nen = element(e)%numnode
             do i = 1, nen
                edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                edof(2*i)   = 2 * element(e)%ix(i)
             end do
             select case( mtype)
             case(1)!product of stiffness and invec
                ke = rho(e)*ke1 + (1.0d0-rho(e))*ke2
                outvector(edof) = outvector(edof)+MATMUL(ke,invector(edof))
             case(2)!product of mass and invec
                me = ( rho(e)*dens1+(1.0d0-rho(e))*dens2 ) * me1
                outvector(edof) = outvector(edof)+MATMUL(me,invector(edof))
             case(3)!product of lumped mass and invec
                me_lumped = ( rho(e)*dens1+(1.0d0-rho(e))*dens2 ) * me_lumped_init
                outvector(edof) = outvector(edof)+me_lumped*invector(edof)
             case(4)!product of damping and invec

                !$$$$$$     ! Forsï¿½g pï¿½ at gï¿½re routinen hurtigere. Men indtil videre afhï¿½nger CZe af rho, sï¿½ det er ikke lykkedes.
                !$$$$$$     do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
                !$$$$$$         if (e == element_rand(i,1) ) then
                !$$$$$$         select case(element_rand(i,2))!face
                !$$$$$$             case(1) ! face 1(bottom of structure)
                !$$$$$$                 outvector(edof) = outvector(edof)+MATMUL(CZe1,invector(edof))
                !$$$$$$             case(2) ! right side
                !$$$$$$                 outvector(edof) = outvector(edof)+MATMUL(CZe2,invector(edof))
                !$$$$$$             case(3) ! top
                !$$$$$$                 outvector(edof) = outvector(edof)+MATMUL(CZe3,invector(edof))
                !$$$$$$             case(4) ! left side
                !$$$$$$                 outvector(edof) = outvector(edof)+MATMUL(CZe4,invector(edof))
                !$$$$$$         end select
                !$$$$$$     
                !$$$$$$     outvector(edof) = outvector(edof)+MATMUL(CZe,invector(edof)) 
                !$$$$$$ 
                !$$$$$$         end if
                !$$$$$$ 
                !$$$$$$     end do

             end select
          end do
          do i = 1, nb ! outvector corrected for BC
             idof = 2*(bound(i,1)-1) + bound(i,2)
             outvector(idof) = 0.0
          end do

       case default ! alle elementer er ikke ens. Brug i stedet global matrix. Hurtigere.

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

             if (mtype == 1) then 
        	!call plane42transient_ke(xe, young1,young2, nu1, nu2, thk,ng,rho(e), ke)

                !outvector(edof) = outvector(edof)+MATMUL(ke,invector(edof))
                print*,'error i mmul. Flag = 0 ikke implementeret'
                stop

             elseif (mtype == 2) then!product of mass and invec
                !call plane42transient_me(xe,thk,dens1, dens2, ng,rho(e), me)! calculate element mass-matrix
                !outvector(edof) = outvector(edof)+MATMUL(me,invector(edof))
                print*,'error i mmul. Flag = 0 ikke implementeret'
                stop

             elseif (mtype == 3) then !product of lumped mass and invec
                !call lumped_mass(lumped_type,me,me_lumped)
                !outvector(edof) = outvector(edof)+me_lumped*invector(edof)
                print*,'error i mmul. Flag = 0 ikke implementeret'
                stop
             elseif (mtype == 4) then !product of damping and invec
                print*,'error i mmul. Flag = 0 ikke implementeret'
                stop	
             end if
          end do
          do i = 1, nb ! Y corrected for BC
             idof = 2*(bound(i,1)-1) + bound(i,2)
             outvector(idof) = 0!bound(i, 3) ! Her skal diagonalelementer ved RB sï¿½ttes lig 0
          end do
       end select

     end subroutine mmul


     subroutine abs_mmul(flag,mtype,rho,invector,outvector,outvector2)
       ! This subroutine multiply element matrix with the inputvector and correct it for bc.

       ! Only to be used when all elements are equal. If this is not the case, its faster
       ! to assemble the global matrices. Or maybe not with many dof.

       ! element_init: is either element stiffnes- og mass matrix. Only used when all
       ! elements are equal


       use fedata
       use plane42transient

       real(8), dimension(neqn), INTENT(IN) :: invector
       real(8), dimension(ne), INTENT(IN) :: rho
       integer, INTENT(IN) :: mtype, flag
       real(8), dimension(neqn), INTENT(OUT) :: outvector
       real(8),optional, dimension(neqn), INTENT(OUT) :: outvector2

       integer :: e, i ,j,ii, eface
       integer :: nen, idof
       integer, parameter :: mdim = 8
       integer, dimension(mdim) :: edof
       real(8), dimension(mdim) :: xe
       real(8), dimension(mdim, mdim) ::  CZe, KRe
       real(8) :: thk, length, p1, p2
       real(8) :: parameters(13),r, r1,r2, nvec(2)

       ! Lame parametre angives i mmul 4.
       ! alpha1, beta1 angives i topology

       outvector = 0.0d0
       if (present(outvector2)) then
          outvector2 = 0.0d0
       end if

       select case(mtype)

       case(1) ! er bï¿½lgeretningen normal pï¿½ fladen? ja
          do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
             e = element_rand(i,1)

             nen = element(e)%numnode
             do j = 1, nen
                xe(2*j-1) = x(element(e)%ix(j),1)
                xe(2*j  ) = x(element(e)%ix(j),2)
                edof(2*j-1) = 2 * element(e)%ix(j) - 1  
                edof(2*j)   = 2 * element(e)%ix(j)
             end do

             thk = mprop(element(e)%mat)%thk

             select case(element_rand(i,2))!face
             case(1) ! face 1(bottom of structure)
                p1 = (x(element(e)%ix(2),1)-x(element(e)%ix(1),1))
                p2 = (x(element(e)%ix(2),2)-x(element(e)%ix(1),2))
                nvec(1) = 0.0d0
                nvec(2) = -1.0d0
             case(2) ! right side
                p1 = (x(element(e)%ix(3),1)-x(element(e)%ix(2),1))
                p2 = (x(element(e)%ix(3),2)-x(element(e)%ix(2),2))
                nvec(1) = 1.0d0
                nvec(2) = 0.0d0
             case(3) ! top
                p1 = (x(element(e)%ix(4),1)-x(element(e)%ix(3),1))
                p2 = (x(element(e)%ix(4),2)-x(element(e)%ix(3),2))
                nvec(1) = 0.0d0
                nvec(2) = 1.0d0
             case(4) ! left side
                p1 = (x(element(e)%ix(1),1)-x(element(e)%ix(4),1))
                p2 = (x(element(e)%ix(1),2)-x(element(e)%ix(4),2))
                nvec(1) = -1.0d0
                nvec(2) = 0.0d0
             end select

             length = 0.50d0*dsqrt(p1**2+p2**2)! length er randens halve lï¿½ngde.
             eface = element_rand(i,2)

             ! er bï¿½lgeretningen normal pï¿½ fladen? ja
             call plane42transient_CZe(xe,eface, young1,young2,dens1,dens2, nu1, nu2, thk, length,ng,rho(e), CZe)
             outvector(edof) = outvector(edof)+MATMUL(CZe,invector(edof)) 
          end do

       case (0)! er bï¿½lgeretningen normal pï¿½ fladen? Nej
          do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
             e = element_rand(i,1)

             nen = element(e)%numnode
             do j = 1, nen
                xe(2*j-1) = x(element(e)%ix(j),1)
                xe(2*j  ) = x(element(e)%ix(j),2)
                edof(2*j-1) = 2 * element(e)%ix(j) - 1  
                edof(2*j)   = 2 * element(e)%ix(j)
             end do

             thk = mprop(element(e)%mat)%thk

             select case(element_rand(i,2))!face
             case(1) ! face 1(bottom of structure)
                p1 = (x(element(e)%ix(2),1)-x(element(e)%ix(1),1))
                p2 = (x(element(e)%ix(2),2)-x(element(e)%ix(1),2))
                nvec(1) = 0.0d0
                nvec(2) = -1.0d0
             case(2) ! right side
                p1 = (x(element(e)%ix(3),1)-x(element(e)%ix(2),1))
                p2 = (x(element(e)%ix(3),2)-x(element(e)%ix(2),2))
                nvec(1) = 1.0d0
                nvec(2) = 0.0d0
             case(3) ! top
                p1 = (x(element(e)%ix(4),1)-x(element(e)%ix(3),1))
                p2 = (x(element(e)%ix(4),2)-x(element(e)%ix(3),2))
                nvec(1) = 0.0d0
                nvec(2) = 1.0d0
             case(4) ! left side
                p1 = (x(element(e)%ix(1),1)-x(element(e)%ix(4),1))
                p2 = (x(element(e)%ix(1),2)-x(element(e)%ix(4),2))
                nvec(1) = -1.0d0
                nvec(2) = 0.0d0
             end select

             length = 0.50d0*dsqrt(p1**2+p2**2)! length er randens halve lï¿½ngde.
             eface = element_rand(i,2)

             ! er bï¿½lgeretningen normal pï¿½ fladen? Nej
             r = abs_rand(i)%r
             r1 = abs_rand(i)%r_enhed(1)
             r2 = abs_rand(i)%r_enhed(2)

             parameters(1) = young1
             parameters(2) = young2
             parameters(3) = dens1
             parameters(4) = dens2
             parameters(5) = nu1
             parameters(6) = nu2
             parameters(7) = thk
             parameters(8) = length
             parameters(9) = r
             parameters(10) = r1
             parameters(11) = r2
             parameters(12) = nvec(1)
             parameters(13) = nvec(2)

             ! Nedenstï¿½ende er for at teste den numeriske integration i forhold til de eksakte indtastede matricer i plane42transient_CZe
             !$$$$$$         r = 0.0d0
             !$$$$$$         parameters(1) = young1
             !$$$$$$         parameters(2) = young2
             !$$$$$$         parameters(3) = dens1
             !$$$$$$         parameters(4) = dens2
             !$$$$$$         parameters(5) = nu1
             !$$$$$$         parameters(6) = nu2
             !$$$$$$         parameters(7) = thk
             !$$$$$$         parameters(8) = length
             !$$$$$$         parameters(9) = r
             !$$$$$$         parameters(10) = nvec(1)
             !$$$$$$         parameters(11) = nvec(2)
             !$$$$$$         parameters(12) = nvec(1)
             !$$$$$$         parameters(13) = nvec(2)


             ! Dï¿½mpningsbidraget ligges til outvector.
             ! Stivhedsbidraget ligges til den globale stivhedsmatrix.
             call plane42transient_CZe_punkt(xe,parameters,eface,ng,rho(e), CZe)
             outvector(edof) = outvector(edof)+MATMUL(CZe,invector(edof))

             if (present(outvector2)) then
        	call plane42transient_KRe(xe, parameters,eface,ng,rho(e), KRe)
        	outvector2(edof) = outvector2(edof)+MATMUL(KRe,invector(edof))
             end if


          end do
       end select

       do i = 1, nb ! Y corrected for BC
          idof = 2*(bound(i,1)-1) + bound(i,2)
          outvector(idof) = 0!bound(i, 3) ! Her skal diagonalelementer ved RB sï¿½ttes lig 0
          if (present(outvector2)) then
             outvector2 = 0.0d0
          end if
       end do

     end subroutine abs_mmul

     subroutine abs_stivhed(rho)

       use fedata
       use plane42transient

       real(8), dimension(:), intent(in) :: rho
       integer :: e, i ,j,ii, eface
       integer :: nen
       integer, parameter :: mdim = 8
       integer, dimension(mdim) :: edof
       real(8), dimension(mdim) :: xe
       real(8), dimension(mdim, mdim) :: KRe
       real(8) :: thk, length, p1, p2
       real(8) :: parameters(13),r, r1,r2, nvec(2)

       do e = 1, ne

          do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
             if (e == element_rand(i,1) ) then

                ! Find coordinates and degrees of freedom
                nen = element(e)%numnode
                do ii = 1, nen
                   xe(2*ii-1) = x(element(e)%ix(ii),1)
                   xe(2*ii  ) = x(element(e)%ix(ii),2)
                   edof(2*ii-1) = 2 * element(e)%ix(ii) - 1  
                   edof(2*ii)   = 2 * element(e)%ix(ii)
                end do
                thk = mprop(element(e)%mat)%thk

                select case(element_rand(i,2))!face
                case(1) ! face 1(bottom of structure)
                   p1 = (x(element(e)%ix(2),1)-x(element(e)%ix(1),1))
                   p2 = (x(element(e)%ix(2),2)-x(element(e)%ix(1),2))
                   nvec(1) = 0.0d0
                   nvec(2) = -1.0d0
                case(2) ! right side
                   p1 = (x(element(e)%ix(3),1)-x(element(e)%ix(2),1))
                   p2 = (x(element(e)%ix(3),2)-x(element(e)%ix(2),2))
                   nvec(1) = 1.0d0
                   nvec(2) = 0.0d0
                case(3) ! top
                   p1 = (x(element(e)%ix(4),1)-x(element(e)%ix(3),1))
                   p2 = (x(element(e)%ix(4),2)-x(element(e)%ix(3),2))
                   nvec(1) = 0.0d0
                   nvec(2) = 1.0d0
                case(4) ! left side
                   p1 = (x(element(e)%ix(1),1)-x(element(e)%ix(4),1))
                   p2 = (x(element(e)%ix(1),2)-x(element(e)%ix(4),2))
                   nvec(1) = -1.0d0
                   nvec(2) = 0.0d0
                end select

                length = 0.50d0*dsqrt(p1**2+p2**2)! length er randens halve lï¿½ngde.
                eface = element_rand(i,2)

                r = abs_rand(i)%r
                r1 = abs_rand(i)%r_enhed(1)
                r2 = abs_rand(i)%r_enhed(2)

                parameters(1) = young1
                parameters(2) = young2
                parameters(3) = dens1
                parameters(4) = dens2
                parameters(5) = nu1
                parameters(6) = nu2
                parameters(7) = thk
                parameters(8) = length
                parameters(9) = r
                parameters(10) = r1
                parameters(11) = r2
                parameters(12) = nvec(1)
                parameters(13) = nvec(2)

                ! Stivhedsbidraget ligges til den globale stivhedsmatrix.
                call plane42transient_KRe(xe, parameters,eface,ng,rho(e), KRe)

                do ii = 1, 2*nen
                   do j =1, 2*nen
                      if (edof(ii)>=edof(j)) then ! brug kun elementer der stï¿½r under eller pï¿½ diagonalen i global K+
                         k((edof(ii)-edof(j)+1),edof(j)) = k((edof(ii)-edof(j)+1),edof(j)) + KRe(ii,j)
                      end if
                   end do
                end do
             end if
          end do
       end do

     end subroutine abs_stivhed


     subroutine elements_init(file,L)
       !  This subroutine find the nodes in the end and the beginning af a beam.

       use fedata
       use plot_routiner

       integer, optional, intent(IN) :: file
       real(8), optional, intent(INOUT) :: L(:)
       integer :: i, e, nen , number_elements_x, number_elements_y, number_nodes_y, number_nodes_x, &
            center_node_y, center_node, node, index1, index2, ny,nx
       real(8) :: xmin, xmax, ymin, ymax, lx, ly, length_first_element, hight_first_element
       integer(8), dimension(:),allocatable :: element_top, element_bottom


       ! NB husk at neqn = 2*nn

       ! Find maximum size of undeformed structure - ONLY valid for rectangular elements.
       xmin = 1d10
       xmax = -1d10
       ymin = 1d10
       ymax = -1d10
       do e = 1,ne 
          nen = element(e)%numnode
          xmin = min(real(minval(x(element(e)%ix(1:nen),1))),xmin)
          xmax = max(real(maxval(x(element(e)%ix(1:nen),1))),xmax)  
          ymin = min(real(minval(x(element(e)%ix(1:nen),2))),ymin)
          ymax = max(real(maxval(x(element(e)%ix(1:nen),2))),ymax)
       end do
       lx = xmax - xmin
       ly = ymax - ymin

       length_first_element= x(element(1)%ix(2),1)-x(element(1)%ix(1),1)
       hight_first_element= x(element(1)%ix(4),2)-x(element(1)%ix(1),2)
       number_elements_x = NINT(lx/length_first_element) ! Antal elemnter i x-retningen ! runder op til nï¿½rmeste integer. Hvis nint ikke medtages, rundes der forkert af
       number_nodes_x = number_elements_x+1 ! Alternativt nn/number_nodes_y
       number_elements_y = NINT(ly/hight_first_element) ! Antal elemnter i y-retningen
       number_nodes_y = number_elements_y+1 !Antal knuder i y-retningen

       ! Getting all dofs i y-direction in the beginning of beam
       allocate(element_beginning(number_elements_y),element_end(number_elements_y))
       allocate( center_dofs_x(number_nodes_x), dof_x_end(number_nodes_y) )
       allocate(element_top(number_elements_x),element_bottom(number_elements_x))

       ! Elements i venstre og hï¿½jre side
       do i=1,number_elements_y
          element_beginning(i) = i
          element_end(i) = (ne-number_elements_y+i)
          !print*,'elemet_end',element_end(i)
       end do

       ! elements top and bottom. Bemï¿½rk at hjï¿½rneelementer gï¿½r igen, da bï¿½de venstre/hï¿½jre og ï¿½vre/nedre rand skal med
       do i=1,number_elements_x
          element_top(i) = i*number_elements_y
          element_bottom(i) = 1+(i-1)*number_elements_y
       end do

       !$$$$$$     rand = 2 ! hvilke sider skal dï¿½mpes.
       select case(rand)
       case(4) ! alle fire sider skal dï¿½mpes
          allocate(element_rand(2*number_elements_y+2*(number_elements_x) , 2))
          element_rand = 0
          index1 = 1
          index2 = number_elements_y
          element_rand(index1:index2,1) = element_beginning
          element_rand(index1:index2,2) = 4 !face

          index1 = number_elements_y+1
          index2 = 2*number_elements_y
          element_rand(index1:index2,1) = element_end
          element_rand(index1:index2,2) = 2

          index1 = 2*number_elements_y+1
          index2 = 2*number_elements_y+(number_elements_x)
          element_rand(index1:index2,1) = element_bottom
          element_rand(index1:index2,2) = 1

          index1 = 2*number_elements_y+(number_elements_x)+1
          index2 = 2*number_elements_y+2*(number_elements_x)
          element_rand(index1:index2,1) = element_top
          element_rand(index1:index2,2) = 3

       case(1) ! kun hï¿½jre side skal dï¿½mpes
          allocate(element_rand(number_elements_y, 2))
          element_rand = 0
          index1 = 1
          index2 = number_elements_y
          element_rand(index1:index2,1) = element_end
          element_rand(index1:index2,2) = 2

       case(2) ! Venstre og hï¿½jre side skal dï¿½mpes
          allocate(element_rand(2*number_elements_y, 2))
          element_rand = 0
          index1 = 1
          index2 = number_elements_y
          element_rand(index1:index2,1) = element_beginning
          element_rand(index1:index2,2) = 4 !face

          index1 = number_elements_y+1
          index2 = 2*number_elements_y
          element_rand(index1:index2,1) = element_end
          element_rand(index1:index2,2) = 2
       end select

       !call output_matrix(real(element_rand),'element_rand') ! Hvis element_rand ikke kan udskrives pga integer: ï¿½ndre element_rand til real(8) i fedata
       
       ! undersï¿½g om frihedsgrader i hjï¿½rneknuder indgï¿½r i L sï¿½ de tï¿½ller med i objektfunktionen
       ny = number_nodes_y
       nx = number_nodes_x
       if (present(file)) then
          select case(file) ! 2:3 angiver disk
          case(2:3)
             if (L(1) == 0) then! nederste venstre hjï¿½rne
        	L(1) = 1d0; end if
             if (L(2) == 0) then
                L(2) = 1d0; end if
             if (L(2*ny-1) == 0) then! ï¿½verste venstre hjï¿½rne
                L(2*ny-1) = 1d0; end if
             if (L(2*ny) == 0) then
                L(2*ny) = 1d0; end if
             if (L(2*ny*(nx-1)+1) == 0) then! nederste hï¿½jre hjï¿½rne
                L(2*ny*(nx-1)+1) = 1d0; end if
             if (L(2*ny*(nx-1)+2) == 0) then
                L(2*ny*(nx-1)+2) = 1d0; end if
             if (L(2*ny*nx-1) == 0) then! ï¿½verste hï¿½jre hjï¿½rne
                L(2*ny*nx-1) = 1d0; end if
             if (L(2*ny*nx) == 0) then
                L(2*ny*nx) = 1d0
             end if
          end select
       end if
       
       
       if (MOD(number_nodes_x,2) == 0) then ! lige antal knuder i x-retning(dvs. der er ikke en knude prï¿½cis i x_retning midten)
          center_node = NINT(number_nodes_x/2.0)*number_nodes_y + NINT(number_nodes_y/2.0) ! runder op
       else ! der er en knude prï¿½cis i midten
          center_node = NINT(number_nodes_x/2.0)*number_nodes_y - NINT(number_nodes_y/2.0)+1
       end if
       midter_dof_x = center_node*2-1
       midter_dof_y = center_node*2
       print*,'Midter_dof_x', Midter_dof_x  

       if (MOD(number_nodes_y,2) == 0) then ! lige antal knuder i y-retning(dvs. der er ikke en knude prï¿½cis i midten)
          center_node_y = NINT(number_nodes_y/2.0)
       else ! der er en knude prï¿½cis i midten
          center_node_y = NINT(number_nodes_y/2.0)
          print*,'der er en knude prï¿½cis i midten i hï¿½jde-retning(y)'
       end if

       do i=1,(number_nodes_x )
          node = center_node_y + (i-1)*number_nodes_y
          center_dofs_x(i) = 2 * node - 1 ! Det er DOF i x-retning der skal bruges og ikke selve knudenummeret
          !$$$$$$          print*,'center_nodes', center_nodes_x(i)
       end do

       ! DOF_x for enden af elementet. Det forudsï¿½tter at element_end(1) er det fï¿½rste element af endeelementerne, samt at
       ! endeelementerne er sammenhï¿½ngende, sï¿½ nummereringen er kontinuert. Ved en normal bjï¿½lkestruktur er det opfyldt.
       e = element_end(1)
       node = (element(e)%ix(2) )*2 -1
       do i=1,number_nodes_y
          dof_x_end(i) = node + (i-1)*2
       end do

     end subroutine elements_init

     subroutine deltaT_init(rho,deltaT)
       ! calculate deltaT from the CLF-condition, page 413 i COOK.

       use fedata

       real(8), INTENT(OUT) :: deltaT
       real(8), INTENT(in) :: rho(:)
       real(8) :: dens, young, nu, lambda, mu, cL1,cL2,cT1, cT2, lengthX, lengthY, c


       lengthX= x(element(1)%ix(2),1)-x(element(1)%ix(1),1)
       lengthY= x(element(1)%ix(4),2)-x(element(1)%ix(1),2)

       ! Lame parameters  
       lambda = young1*nu1/ ((1d0+nu1)*(1d0-2d0*nu1)) 
       mu = young1/(2d0*(1d0+nu1))


       cL1 = dsqrt( (lambda+2.0*mu) / dens1)
       cT1 = dsqrt( mu / dens1)

       ! Lame parameters  
       lambda = young2*nu2/ ((1.0+nu2)*(1.0-2.0*nu2)) 
       mu = young2/(2.0*(1.0+nu2))! shear modulus

       cL2 = dsqrt( (lambda+2.0*mu) / dens1)
       cT2 = dsqrt( mu / dens1)

       c = max(cL1,cT1,cL2,cT2)

       deltaT = MIN(lengthX,lengthY)/c
       print*,'deltaT_max',deltaT

     end subroutine deltaT_init

     subroutine explosive_source(Dcurr,radius,center,g)
       ! beregner vï¿½gtnings-vektoren (g(r)) ved en eksplosiv source med radius a og centrum c(x,y).

       use fedata

       real(8), intent(in) :: radius, center(:), Dcurr(:)!Fin(:)
       real(8), dimension(:), intent(out) :: g

       integer :: i
       real(8) :: p1, p2, r, c1, c2,a

       c1 = center(1)
       c2 = center(2) 
       a = radius
       g = 0.0d0
       do i=1,nn

          p1 = x(i,1)+Dcurr(2*i-1)! x-coordinate of node
          p2 = x(i,2)+Dcurr(2*i)! y-coordinate of node

          r = dsqrt( (p1-c1)**2 + (p2-c2)**2 )

          if(r<= a) then! punktet ligger pï¿½ disken
             if (r==0) then
                cycle
                !r = r+1E-12
             end if
             !$$$$$$             print*,'hej'
             g(2*i-1) = ( 1.0d0-r**2/a**2 ) * ( p1 - c1 )/r
             g(2*i) = ( 1.0d0-r**2/a**2 ) * ( p2 - c2 )/r
             !print*,'x-ret',( 1.0d0-r**2/a**2 ) * ( p1 - c1 )/r
             !print*,'y-ret',( 1.0d0-r**2/a**2 ) * ( p2 - c2 )/r
          end if
       end do

     end subroutine explosive_source


     subroutine abs_init(Dcurr,center)

       use fedata
       use plot_routiner
       real(8), intent(in) :: center(:), Dcurr(:)
       !    real(8), dimension(:), intent(out) :: g
       real(8) , allocatable :: vector(:)

       integer :: i,e, kk, edof(8), j, nen
       real(8) :: p1, p2, r,r1,r2, c1, c2, xe(8)

       select case(rand)
       case(0) ! ingen dï¿½mpning
          return
       case default

          if (.not. allocated(abs_rand)) then
             allocate(abs_rand(SIZE(element_rand,1)))
             !$$$$$$         allocate(vector(SIZE(element_rand,1)))
          end if

          c1 = center(1)
          c2 = center(2) 
          kk = 0
          do i = 1, SIZE(element_rand,1) ! Arbsorbing BC's
             e = element_rand(i,1)
             !$$$$$$ 
             nen = element(e)%numnode
             do j = 1, nen
                xe(2*j-1) = x(element(e)%ix(j),1)
                xe(2*j  ) = x(element(e)%ix(j),2)
                edof(2*j-1) = 2 * element(e)%ix(j) - 1  
                edof(2*j)   = 2 * element(e)%ix(j)
             end do

             select case(element_rand(i,2))!face
             case(1) ! face 1(bottom of element)
                ! koordinater for midten af randen mellem knude 1 og 2 i elementet
                p1 = 0.50d0*(xe(1)+xe(3)) + & ! udeformeret
                     0.50d0*(Dcurr(edof(1)) + Dcurr(edof(3)) ) ! deformation

                p2 = 0.50d0*(xe(2)+xe(4)) + &
                     0.50d0*(Dcurr(edof(2)) + Dcurr(edof(4)) )
             case(2) ! right side
                p1 = 0.50d0*(xe(3)+xe(5)) + &
                     0.50d0*(Dcurr(edof(3)) + Dcurr(edof(5)) )

                p2 = 0.50d0*(xe(4)+xe(6)) + &
                     0.50d0*(Dcurr(edof(4)) + Dcurr(edof(6)) )
             case(3) ! top
                p1 = 0.50d0*(xe(5)+xe(7)) + &
                     0.50d0*(Dcurr(edof(5)) + Dcurr(edof(7)) )

                p2 = 0.50d0*(xe(6)+xe(8)) + &
                     0.50d0*(Dcurr(edof(6)) + Dcurr(edof(8)) )
             case(4) ! left side
                p1 = 0.50d0*(xe(7)+xe(1)) + &
                     0.50d0*(Dcurr(edof(7)) + Dcurr(edof(1)) )

                p2 = 0.50d0*(xe(8)+xe(2)) + &
                     0.50d0*(Dcurr(edof(8)) + Dcurr(edof(2)) )
             end select

             ! enhedsvektor fra centrum mod midten af randen
             r1 = p1-c1
             r2 = p2-c2
             r = dsqrt( (r1)**2 + (r2)**2 )
             r1 = r1 /r 
             r2 = r2 /r

             abs_rand(i)%r = r
             abs_rand(i)%r_enhed(1) = r1
             abs_rand(i)%r_enhed(2) = r2
             !kk = kk+1
             !vector(kk) = i

          end do
          !call output_vector(vector,'index_rand')

       end select


     end subroutine abs_init

   end module transient
