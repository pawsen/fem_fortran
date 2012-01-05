module fea

  implicit none

  private
  public :: displ, initial_fea, buildload, buildstiff_fea, enforce_fea, recover, volume, areal
  public :: buildmass_fea, mmul_fea, build_mvec, mvec_mul, eigen

  public :: assemble_sparse, sparse_add_lagrangian


contains

  subroutine initial_fea

    ! This subroutine is mainly used to allocate vectors and matrices
    use fedata
    use link1
    use plane42

    integer :: n, e, nen, min_edof, max_edof, diff_edof, nzz, nzz2

    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof

    ! This subroutine computes the number of global equation,
    ! half bandwidth, etc and allocates global arrays.

    ! Calculate number of equations
    bw=0 ! giver bw en v�rdi til if-s�tningen
    neqn = 2*nn
    if (banded == 0) then ! dense format
       allocate (k(neqn, neqn))
    elseif(banded == 1) then !allocate k=banded
       do e=1,ne
          nen = element(e)%numnode ! Antallet af knudepunkter i et element
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


       !$$$$$$         bw=0
       !$$$$$$         do e=1,ne
       !$$$$$$             max_edof=MAXVAL(element(e)%ix(1:4))
       !$$$$$$             min_edof=MINVAL(element(e)%ix(1:4))
       !$$$$$$             diff_edof=(max_edof-min_edof+1)*2
       !$$$$$$             if (diff_edof > bw) then
       !$$$$$$                 bw=diff_edof
       !$$$$$$             end if
       !$$$$$$         end do


       allocate (k(bw,neqn))
       print*,'bw',bw

    end if

    if (banded ==  2) then ! sparse format with lagrangian multiplier

       if ( ((antype == 'EIGEN') .and.  (eigenvalue%shift .eqv. .false.)) ) then
          neqn_nb = neqn
          nzz = 8*8*ne
       else
          neqn_nb = neqn+nb
          nzz = 8*8*ne+2*nb 
          !8*8*ne is the number of times the local stiffness matrix adds a value to the global(including dublicates, eg. is't not nz(number of non-zeros)). Because of lagrangien multiplier there is added xtra two non-zero entries for each RB
       end if

       if (antype /= 'EIGEN') then
          allocate (p(neqn+nb), d(neqn+nb))
       end if

       allocate (iK(nzz),jK(nzz),sK(nzz))
    else! normalt
       allocate (p(neqn), d(neqn))
    endif

    if (antype /= 'EIGEN') then
       allocate (strain(ne, 3), stress(ne, 3))
       strain = 0d0
       stress = 0d0
    end if

  end subroutine initial_fea

  subroutine displ(flag,rho,rho_min)

    ! This subroutine calculates displacements
    use numeth
    use fedata
    use processor
    use plot_routiner
    use solve_handle_real
    use exodus

    integer, INTENT(IN) :: flag
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(INOUT) :: rho_min

    integer :: e
    real(8) :: plotval(ne)

    !linux_test
    !call input

    ! det er ikke nødvendigt at kalde buildload ved topology-opt. Det sker i topology
    if ((antype == 'STATIC') .or. (antype == 'COUPLED')) then
       call buildload
    end if

    if (present(rho)) then
       call buildstiff_fea(flag,rho,rho_min)
    else! Build stiffness matrix
       call buildstiff_fea(flag)
    endif

    ! Remove rigid body modes
    call enforce_fea

    d = p
    if (banded /=2) then ! not sparse
       if (banded == 0) then
          if (sum(k) > nb) then
             call factor(k)
             ! solve for displacement vector
             call solve(k, d)
          else
             print*,'WARNING in DISPL: Stiffness equal to zero'
             stop
          end if
       elseif (banded == 1) then !banded
          if (sum(k) > nb) then
             call bfactor(k)
             call bsolve(k, d) !notice that d=p when bsolve is called
          else
             print*,'WARNING in DISPL: Stiffness equal to zero'
             stop
          end if
       end if
    else ! sparse format
       if (present(rho)) then
          call mumps_solve_real(5)
       else
          call mumps_init_real
          call mumps_solve_real(6)
          call mumps_finalize_real
       end if
    endif


    print*,'asd',antype
    ! Ikke topOpt problem, dvs print solution
    if ( ((antype == 'STATIC') .or. (antype == 'COUPLED')) ) then
       call recover ! p bliver beregnet igen
       do e = 1, ne
          plotval(e) = (stress(e,1)**2+stress(e,2)**2-stress(e,1)*stress(e,2)+3.0*stress(e,3)**2)**(0.50d0) !von mises stresses
       end do

       call exodus_init
       call exodus_write_node(1,d)
       call exodus_write_time(1,1.0d0)
       call exodus_finalize

       call output !linux_fejl
      
       if (present(rho)) then
          !call plot('elements', 'xwin', 'color', 'spaendingsintensitet', plotval)
          call output_deformed('deformeret','densitet',rho)
          call output_elements('', plotval)
          call plot('elements', 'xwin', 'color', 'spaendingsintensitet', plotval)
          call plot('deformed', 'xwin', 'color', 'deformeret',rho)

       else

          !call plot('elements', 'xwin', 'color', 'spaendingsintensitet', plotval) !linux_fejl
          !call plot('un/deformed', 'xwin', 'color', ' ')
          !call output_elements('', plotval)
       end if
    endif


  end subroutine displ

  subroutine buildload(rho,rho_min)

    ! This subroutine builds the global load vector.

    use fedata
    use plane42

    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(IN) :: rho_min
    integer :: i, idof, e, nen, eface , n ,j, idof2(4)
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(mdim) :: re ,rvol, rtemp
    real(8) :: fe , thk, dens, alpha, nu, young, tstart, tend, deltaT,Ae(8,4), tn_e(4)


    ! Build load vector from static force
    p = 0.0d0
    do i = 1, np !np = number loads
       ! Build nodal load contribution
       if (loads(i, 1) == 1) then
          idof = 2*(loads(i,2)-1) + loads(i,3)
          p(idof) = loads(i, 4)
          ! Build uniformly distributed surface (pressure) load contribution
       else if (loads(i, 1) == 2) then
          e = loads(i,2) ! element SFE virker p�
          nen = element(e)%numnode ! Antallet af knudepunkter i �t element
          do n = 1, nen
             xe(2*n-1) = x(element(e)%ix(n),1)
             xe(2*n  ) = x(element(e)%ix(n),2)
             edof(2*n-1) = 2 * element(e)%ix(n) - 1  
             edof(2*n)   = 2 * element(e)%ix(n)
          end do
          eface = loads(i,3) !face id
          fe = loads(i,4)!v�rdi af load
          thk = mprop(element(e)%mat)%thk ! tykkelse af strukturen
          call plane42_re(xe, eface, fe, thk, ng, re)

          do j = 1, 2*nen
             p(edof(j)) = p(edof(j)) + re(j)
          end do


       end if
    end do


    do e=1,ne ! adding volumenforces
       if (element(e)%id == 2) then ! kontinuum!
          nen = element(e)%numnode ! Antallet af knudepunkter i eet element
          do n = 1, nen
             xe(2*n-1) = x(element(e)%ix(n),1)
             xe(2*n  ) = x(element(e)%ix(n),2)
             edof(2*n-1) = 2 * element(e)%ix(n) - 1  
             edof(2*n)   = 2 * element(e)%ix(n)

          end do
          dens = mprop(element(e)%mat)%dens ! densitet af element
          thk = mprop(element(e)%mat)%thk ! tykkelse af element
          call plane42_vol(xe, accel, dens, thk, ng, rvol)

          do j = 1, 2*nen
             p(edof(j)) = p(edof(j)) + rvol(j)
          end do
       end if
    end do

    ! Add calculated temp to the forcevector
    !rho_min = 0.001 ! til thermal
    if ((antype == 'COUPLED') .or. (antype == 'TOPTHERMSTRUCT') .or. (antype == 'TOPTHERMSTRUCT_hard')) then
       do e = 1, ne
          nen = element(e)%numnode
          do j = 1, nen	! Do 2
             xe(2*j-1) = x(element(e)%ix(j),1)
             xe(2*j  ) = x(element(e)%ix(j),2)
             edof(2*j-1) = 2 * element(e)%ix(j) - 1  
             edof(2*j)   = 2 * element(e)%ix(j)
             idof2(j)	= element(e)%ix(j)
          end do
          thk = mprop(element(e)%mat)%thk
          young = mprop(element(e)%mat)%young
          alpha = mprop(element(e)%mat)%alpha
          nu = mprop(element(e)%mat)%nu
          tstart = mprop(element(e)%mat)%tstart
          tend = t_elem(e)
          ! Temperature change of element
          deltaT = tend - tstart
          call plane42_thermLoad(xe, ng, young, nu, alpha, thk, deltaT, rtemp)
          !$$$$$$         call plane42_thermKobling(xe, ng, young, nu, alpha, thk,Ae)
          !$$$$$$         tn_e = tn(idof2) ! knude_temp
          !$$$$$$         rtemp = MATMUL(Ae,tn_e) ! themisk bidrag til lastvektoren

          if (present(rho)) then
             do j = 1, 2*nen
                p(edof(j)) = p(edof(j)) + (rho_min+(1.0d0-rho_min)*rho(e)**penal)*rtemp(j)
             end do
          else
             p(edof) = p(edof) + rtemp
          end if
       end do
       !print*,'p*p',dot_product(p,p)
    end if


  end subroutine buildload


  subroutine buildstiff_fea(flag,rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata
    use link1
    use plane42
    use plane42nonlin
    use plane42transient
 
    integer, INTENT(IN) :: flag
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(IN) :: rho_min


    integer :: e, i, j, ii, jj
    integer :: nen, irow, icol, idof
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(mdim, mdim) :: ke, ke1,ke2
    real(8) :: young, nu, area, dens, thk
    real(8), allocatable :: me_lumped(:), me(:,:),sigma


    ! If we use shift, the stiffness should be modified as specified by ARPACK
    if (eigenvalue%shift) then
       allocate(me(mdim,mdim),me_lumped(mdim),sigma)
    end if

    bw = size(k,1)
    if (banded /= 2) then
       ! Reset stiffness matrix
       K = 0d0
    else ! sparse format
       ii = 0
       jj = 0
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
       ! hvis alle elementer er ens, skal nedenst�ende kun k�res �n gang
       if ( ((flag == 0) .or. (e == 1)) .and. (antype /= 'NONLIN') .and. (antype /= 'TRANSIENT')) then
          ! Gather material properties and find element stiffness matrix
          select case( element(e)%id )
          case( 1 )
             young = mprop(element(e)%mat)%young
             area  = mprop(element(e)%mat)%area
             call link1_ke(xe, young, area, ke)
          case( 2 )
             ! calculate plane42rect element stiffness matrix here
             young = mprop(element(e)%mat)%young
             thk = mprop(element(e)%mat)%thk
             nu = mprop(element(e)%mat)%nu
             call plane42_ke(xe, young, nu, thk, ng, ke)
          end select

          ! da ke laves p� baggrund af forskydningen for det enkelte element, er det ikke muligt at bruge "flag" option
       elseif ( antype == 'NONLIN' ) then
          young = mprop(element(e)%mat)%young
          thk = mprop(element(e)%mat)%thk
          nu = mprop(element(e)%mat)%nu
          !$$$$$$         call plane42nonlin_ke(xe,d(edof), young, nu, thk,ng, ke)
          print*,'error her skal jeg ikke være - fea->buildstiff->non_lin'
          stop
       elseif ( ((flag == 0) .or. (e == 1)) .and. (antype == 'TRANSIENT') ) then
          thk = mprop(element(e)%mat)%thk
          call plane42transient_ke(xe, young1,young2, nu1, nu2, thk, ng, ke1,ke2)
       end if

       if ((antype == 'TRANSIENT') ) then
          ke = rho(e)*ke1 + (1.0d0-rho(e))*ke2
       end if

       ! if(eigenvalue%shift) then
       !    dens = mprop(element(e)%mat)%dens
       !    sigma = eigenvalue%sigma
       !    call plane42_me(xe,dens,thk, ng,me)
       !    call lump_mass(nen,2,me,me_lumped)
       !    !ke = ke-sigma*me
       !    do i=1,2*nen
       !       ke(i,i) = ke(i,i)-sigma*me_lumped(i)
       !    end do
       ! end if

       ! Assemble into global matrix
       if (banded == 0) then
          do i = 1, 2*nen
             do j = 1, 2*nen
              	if (.not. present(rho)) then! STATIC, COUPLED, TRANSIENT
                   k(edof(i), edof(j)) = k(edof(i), edof(j)) + ke(i, j)
                else
                   k(edof(i), edof(j)) = k(edof(i), edof(j)) + (rho_min+(1d0-rho_min)*rho(e)**penal)*ke(i,j)
                end if
             end do
          end do
       elseif (banded == 2) then ! sparse, Benytter ikke symmetrien
          if (present(rho) .and. ((element(e)%mat /= elem_id))) then
             call assemble_sparse(nen,2,ii,edof,ke,jj,rho,rho_min)
          else
             call assemble_sparse(nen,2,ii,edof,ke)
          end if
       else ! banded
          do i = 1, 2*nen
             do j =1, 2*nen
                if (edof(i)>=edof(j)) then ! brug kun elementer der står under eller på diagonalen i global K
                   if (.not. present(rho)) then! STATIC, COUPLED, TRANSIENT
                      k((edof(i)-edof(j)+1),edof(j)) = k((edof(i)-edof(j)+1),edof(j)) + ke(i,j)
                   else
                      k((edof(i)-edof(j)+1),edof(j)) = k((edof(i)-edof(j)+1),edof(j)) + &
                           (rho_min+(1d0-rho_min)*rho(e)**penal)*ke(i,j)! Modificeret SIMP-v�gtning af elasticitetsmodul
                   end if
                end if
             end do
          end do
       end if

    end do

    !add values from lagrangian multipliers. BUT not if it's a eigenvalue problem vithout shift
    if ( banded == 2 .and. ((antype /= 'EIGEN') .or. eigenvalue%shift)) then
       call sparse_add_lagrangian(2,ii)
    end if

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
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
    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

  end subroutine buildstiff_fea

  subroutine assemble_sparse(nen,kk,ii,edof,ke,jj,rho,rho_min)
    
    use fedata
    integer, intent(in) :: nen,kk, edof(:)
    integer, intent(inout) :: ii
    real(8), intent(in) :: ke(:,:)
    integer, intent(inout), optional :: jj
    real(8), intent(in), optional :: rho(:), rho_min
    integer :: i,j

    if (present(rho)) then
       jj = jj +1
       do i = 1, kk*nen
          do j =1, kk*nen
             ii = ii+1
             iK(ii) = edof(i)
             jK(ii) = edof(j)
             sK(ii) = (rho_min+(1d0-rho_min)* rho(jj)**penal) * ke(i,j)
             ! Modificeret SIMP-vægtning af elasticitetsmodul
          end do
       end do
    else
       do i = 1, kk*nen
          do j =1, kk*nen
             ii = ii+1
             iK(ii) = edof(i)
             jK(ii) = edof(j)
             sK(ii) = ke(i,j)
          end do
       end do
    end if

  end subroutine assemble_sparse

  subroutine sparse_add_lagrangian(kk,ii)

    use fedata
    
    integer, intent(in) :: kk
    integer, intent(inout) :: ii
    integer :: idof, i, jj
    
    ! UPS: bound(i,3) is used both for temp and moments for plates. Stupid!

    jj = 0
    do i=1,nb
       if ((bound(i,2) < 6) .and. &! forskydning
            ((kk /= 2) .or. (bound(i,2) /= 3) )  ) then ! bound(i,2) = 3 => temp på randen
          jj = jj +1
          idof = 3*(bound(i,1)-1) + bound(i,2)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+jj
          jK(ii) = idof
          sK(ii ) = 1d0
          ! matrix to the right
          ii = ii+1
          iK(ii) = idof
          jK(ii) = neqn+jj
          sK(ii ) = 1d0
       end if
    end do

  end subroutine sparse_add_lagrangian


  subroutine buildmass_fea(flag,rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata
    use plane42

    integer, INTENT(IN) :: flag
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(IN) :: rho_min

    integer :: e, i, j, ii
    integer :: nen, idof
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(mdim, mdim) :: me, me_lumped
    real(8) :: dens, thk

    ii = 0
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
       dens = mprop(element(e)%mat)%dens
       call plane42_me(xe,dens,thk, ng,me)

       me_lumped = 0d0
       do j=1,2*nen
          do i=1,2*nen
             me_lumped(j,j) =me_lumped(j,j)+ me(j,i)
          end do
       end do

       ! Assemble into global matrix
       do i = 1, 2*nen
          do j =1, 2*nen
             ii = ii+1
             ! iK2(ii) = edof(i)
             ! jK2(ii) = edof(j)
             ! sK2(ii) = me(i,j)

             !sK2(ii) = me_lumped(i,j)
          end do
       end do
    end do

    do i=1,nb
       if (bound(i,2) /= 3) then! bound(i,2) = 3 => temp p� randen
          idof = 2*(bound(i,1)-1) + bound(i,2)
          ! nederste matrix
          ii = ii+1
          ! iK2(ii) = neqn+i
          ! jK2(ii) = idof
          ! sK2(ii ) = 1
          ! ! matrix to the right
          ! ii = ii+1
          ! iK2(ii) = idof
          ! jK2(ii) = neqn+i
          ! sK2(ii ) = 1
       end IF
    end do

  end subroutine buildmass_fea
  
  subroutine enforce_fea

    ! This subroutine enforces the support boundary conditions.
    use fedata

    integer :: i, j, idof
    real(8), dimension(neqn) :: kl
    real(8):: penalty_factor

    bw = size(k,1) !half bandwitdh   

    ! Correct for supports
    if ( banded == 0) then
       if (.not. penalty) then
          do i = 1, nb
             idof = 2*(bound(i,1)-1) + bound(i,2)
             p(1:neqn) = p(1:neqn) - k(1:neqn, idof) * bound(i, 3)
             p(idof) = bound(i, 3)
             k(1:neqn, idof) = 0.
             k(idof, 1:neqn) = 0.
             k(idof, idof) = 1.
          end do
       else
          ! Penalty method
          penalty_factor = 10E+9*maxval(K)
          do i = 1, nb
             if (bound(i,2) /= 3) then
                idof = 2*(bound(i,1)-1) + bound(i,2)
                k(idof, idof) = k(idof, idof) + penalty_factor
                p(idof) = penalty_factor * bound(i, 3)
             end if
          end do
       end if
    elseif(banded == 1) then !banded
       if (.not. penalty) then
          ! NB!!!! AT HAVE FORSKYDNINGER FORSKELLIG FRA NUL, VIRKER IKKE FOR BANDED. DET VIRKER FOR DENSE OG SPARSE!!
!!! HELLER IKKE SELVOM K SÆTTES TIL NUL EFTER DE RELEVANTE ELEMENTER I KL ER FUNDET
          do i = 1, nb
             if (bound(i,2) /= 3) then! bound(i,2) = 3 => temp på randen
                idof = 2*(bound(i,1)-1) + bound(i,2)

                k(1:bw, idof) = 0.0d0
                k(1, idof) = 1.0d0
                do j = 1 , bw-1
                   if (idof-j > 0) then
                      k(1+j,idof-j) = 0.0d0 ! s. fig 6.2 i noter(s. 42). Række i 'full matrix' bliver til skrå linje i 'banded matrix'.
                   end if
                end do

                ! nedenstående(kl) finder den værdi der skal lægges til load-vektoren ved imposing af BC, jf (2.7-6) i COOK
                kl(1:neqn)=0.0d0

                do j = 1, idof-1
                   if (idof+1-j<bw) then
                      kl(j) = k(idof+1-j , j) 
                   end if
                end do

                do j = idof , bw
                   kl(j) = k(j-idof+1 , idof)
                end do

                p(1:neqn) = p(1:neqn) - kl(1:neqn) * bound(i, 3)
                p(idof) = bound(i, 3)

             end if
          end do
       else
          ! Penalty method
          penalty_factor = 10E+9 * maxval(K)
          do i = 1, nb
             if (bound(i,2) /= 3) then
                idof = 2*(bound(i,1)-1) + bound(i,2)
                k(1, idof) = k(1, idof) + penalty_factor
                p(idof) = penalty_factor * bound(i, 3)
             end if
          end do
       end if
    elseif(banded == 2) then !sparse. Add constraint to the load vector
       do i = 1,nb
          p(neqn+i) = bound(i, 3)
       end do
    endif

  end subroutine enforce_fea

  subroutine recover

    ! This subroutine recovers the element stress, element strain, 
    ! and nodal reaction forces

    use fedata
    use link1
    use plane42
    use plane42nonlin

    integer :: e, i, nen
    integer, parameter :: mdim = 8
    integer :: edof(mdim)
    real(8), dimension(mdim) :: xe, de
    real(8), dimension(mdim, mdim) :: ke
    real(8) :: young, nu, area, thk, xi ,eta
    real(8), dimension(3) :: estrain, estress

    !Output from plane42nonlin_shape
    REAL(8) :: Bmat(3, 8), B0mat(3, 8),BLmat(3, 8), Nmat(2,8),Gmat(4,8), jac(2,2), detjac

    p = 0.

    do e = 1, ne

       ! Find coordinates etc...
       nen = element(e)%numnode
       do i = 1,nen
          xe(2*i-1) = x(element(e)%ix(i), 1)
          xe(2*i)   = x(element(e)%ix(i), 2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1
          edof(2*i)   = 2 * element(e)%ix(i)
          de(2*i-1) = d(edof(2*i-1))
          de(2*i)   = d(edof(2*i))
       end do

       ! Find stress and strain
       select case( element(e)%id )
       case( 1 )
          young = mprop(element(e)%mat)%young
          area  = mprop(element(e)%mat)%area
          call link1_ke(xe, young, area, ke)
          p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))!beregner loads igen, da p bliver overskrevet ved l�sning af KD = P.
          call link1_ss(xe, de, young, estress, estrain)
          stress(e, 1:3) = estress
          strain(e, 1:3) = estrain
       case( 2 )
          young = mprop(element(e)%mat)%young
          thk = mprop(element(e)%mat)%thk
          nu = mprop(element(e)%mat)%nu
          if ( (antype == 'NONLIN') ) then ! bem�rk dette er ikke de "fysiske" sp�ndinger. Se cauchy-rutinen i stedet. Denne b�r ikke bruges.
             xi = 0.0d0
             eta = 0.0d0
             call plane42nonlin_shape(xe,d(edof),xi,eta, Nmat,Bmat,B0mat,BLmat,Gmat, jac, detjac)
             call plane42nonlin_ss(xe, d(edof),young,nu,B0mat, BLmat, estress, estrain)
          else
             call plane42_ke(xe, young, nu, thk,ng, ke)
             p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
             call plane42_ss(xe, de, young, nu, estress, estrain)
          end if
          stress(e, 1:3) = estress
          strain(e, 1:3) = estrain
       end select
    end do

  end subroutine recover

  subroutine areal(flag,arealet)
    use plane42
    use fedata

    integer, INTENT(IN) :: flag
    real(8), INTENT(OUT) :: arealet
    integer :: e, i, nen, ii
    integer, parameter :: mdim = 8
    real(8), dimension(mdim) :: xe
    real(8) :: element_areal

    arealet =0.0d0
    do e = 1, ne
       if ( (flag == 0) .or. (e == 1) ) then
          nen = element(e)%numnode
          do i = 1,nen
             xe(2*i-1) = x(element(e)%ix(i), 1)
             xe(2*i)   = x(element(e)%ix(i), 2)
          end do
          call plane42_area(xe,ng,element_areal)
       end if
       arealet = arealet + element_areal
    end do

  end subroutine areal

  subroutine volume(flag,inak,vol)
    ! Beregning af areal

    use fedata
    use plane42
    use numeth ! guss points ligger i numeth

    integer, INTENT(IN) :: flag
    type(INAKTIV_def), intent(in) :: inak
    !real(8), optional, INTENT(OUT) :: element_areal
    REAL(8), DIMENSION(:), INTENT(OUT) :: vol
    integer :: e, i, j, ii, nen, ne_aktiv
    integer, parameter :: mdim = 8
    real(8), dimension(mdim) :: xe
    !real(8) :: arealet
    REAL(8) :: detjac, jac(2,2), bmat(3,8), Nmat(2,8),thk, element_areal
    REAL(8), DIMENSION(:) , allocatable:: w, xi, eta
    allocate(w(ng),xi(ng),eta(ng))
    call gauss_points(w,xi,eta,ng)

    vol = 0.0d0
    negative_detjac = 0

    if (inak%bool) then
       ne_aktiv = size(inak%act,1)
    else
       ne_aktiv = ne
    end if

    ! only include volume of active elements
    do ii = 1, ne_aktiv
       
       if ( (flag == 0) .or. (ii == 1) ) then

          if (inak%bool) then
             e = inak%act(ii)
          else
             e = ii
          end if

          nen = element(e)%numnode
          do i=1,nen
             xe(2*i-1) = x(element(e)%ix(i),1)
             xe(2*i  ) = x(element(e)%ix(i),2)
          end do

          element_areal =0.0
          do i=1,ng
             do j=1,ng
                call plane42_shape(xe, xi(i), eta(i), Nmat,bmat, jac, detjac)
            	element_areal = element_areal+W(i)*W(j)*detjac
              	if (detjac<0) then ! Hvis elementerne er for 'forvredne' bliver detjac negativ => artificial ekstra stivhed
                   negative_detjac = negative_detjac+1
            	end if
             end do
          end do
          thk = mprop(element(e)%mat)%thk
       end if
       vol(ii) = element_areal*thk
    end do

    if (negative_detjac /= 0) then
       print*,'Antal negative Jacobi-determinanter', negative_detjac
       error stop
    end if
    !print*,'Areal af et element',element_areal

  end subroutine volume

  subroutine sort_sparse(flag)
    ! sort sparse Matrix given i coordinate format with dublicates to CRS format
    use fedata
    use sort_array
    use plane42

    integer, intent(IN) :: flag

    TYPE :: sparse_index
       integer, DIMENSION(:), POINTER :: colum
       integer :: k
    END TYPE sparse_index
    TYPE(sparse_index), ALLOCATABLE :: row(:)
    integer, POINTER :: tmp(:)
    integer :: e,i,j, m, nz, jj, jj_old, index, nen

    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(mdim, mdim) :: ke
    real(8) :: young, nu, area, dens, thk


    m = 10
    
    ALLOCATE(row(neqn))
    DO i = 1,neqn
       ALLOCATE(row(i)%colum(m))
    END DO


    row(:)%k = 0 ! keeps track of elements in corressponding vector
    do e = 1, ne
       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
       end do
       do i = 1, 2*nen
          do j =1, 2*nen
             row(i)%k = row(i)%k +1
             row(edof(i))%colum(row(i)%k) = edof(j)
          end do
       end do
    end do


    do i=1,neqn
       ! Remove dublicates and sort afterwards
       call remove_dups(row(i)%colum(:),tmp,row(i)%k)
       DEALLOCATE(row(i)%colum)
       ALLOCATE(row(i)%colum(row(i)%k))
       row(i)%colum = tmp
       ! sorting
       call heapsort(row(i)%colum(:))
    end do

    ! Global colum vector, containing colum-index for each row
    nz = SUM(row(:)%k)!non-zero
    allocate(jK(nz),sK(nz),iK(neqn+1))

    
    jj_old = 1
    do i=1,neqn
       jj = sum(row(1:i)%k)
       jK(jj_old:jj) = row(i)%colum(:)
       jj_old = SUM(row(1:i)%k)
    end do

    ! Global row_index vector
    iK(1) = 1
    do i=1,neqn
       ik(i+1) = iK(i) + row(i)%k
    end do


    ! assemble value vector
    !finde den tilhørende frihedsgrad i jK.
       ! Ved først at slå op i iK, fås intervallet i hvilket frihedsgraden ligger i jK.
       ! I dette interval søges så for frihedsgraden. Index for frihedsgraden fortæller
       ! hvor ke(i,j) skal placeres i sK
       do i = 1, 2*nen
          do j =1, 2*nen
             jj_old = iK(edof(i))
             jj = iK(edof(i)+1)
             index = binarySearch(jK(jj_old:jj),edof(j))
             index = index + jj_old
             sK(index) = ke(i,j)
          end do
       end do
   ! end do


  end subroutine sort_sparse

  subroutine mmul_fea(invector,outvector,mtype)
    ! This subroutine calculates the outvector and corrects it for BC

    use fedata
    use plane42
    use numeth

	real(8), dimension(neqn), INTENT(IN) :: invector
    real(8), dimension(neqn), INTENT(OUT) :: outvector
    integer, INTENT(IN) :: mtype
    integer :: e, i, j
    integer :: nen, idof
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe, me_lumped
    real(8), dimension(mdim, mdim) :: ke, me, ce
    real(8) :: young, nu, area, dens, thk
    real(8) :: me_diag_sum, me_tot

	outvector = 0d0
    do e = 1, ne

       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i  ) = x(element(e)%ix(i),2)
          edof(2*i-1) = 2 * element(e)%ix(i) - 1  
          edof(2*i)   = 2 * element(e)%ix(i)
       end do

       ! calculate element mass-matrix
       young = mprop(element(e)%mat)%young
       thk = mprop(element(e)%mat)%thk
       nu = mprop(element(e)%mat)%nu
       dens = mprop(element(e)%mat)%dens

       select case(mtype)
       case(1)!product of stiffness and invec
          call plane42_ke(xe, young, nu, thk,ng, ke)
          outvector(edof) = outvector(edof) +MATMUL(ke,invector(edof))
       case(2)!product of mass and invec
          call plane42_me(xe,dens,thk, ng,me)
          outvector(edof) = outvector(edof)+MATMUL(me,invector(edof))
       case(3)!product of lumped mass and invec
          call plane42_me(xe,dens,thk, ng, me)! calculate element mass-matrix
          me_lumped = 0d0
          do j=1,2*nen
             do i=1,2*nen
                me_lumped(j) =me_lumped(j)+ me(j,i)
             end do
          end do
          outvector(edof) = outvector(edof)+me_lumped*invector(edof)
       end select
       
    end do

    ! OUTVEC corrected for BC. NB sætter boud til 0.
    call correct_bound_vector(outvector)

  end subroutine mmul_fea

  subroutine build_mvec(rho)
    ! This subroutine calculates the outvector and corrects it for BC

    use fedata
    use plane42
    use mindlin42
    use numeth

    real(8),optional, INTENT(IN) :: rho(:)

    integer :: e, i, j
    integer :: nen, idof, mdim, kk
   
    real(8), dimension(8) :: xe
    real(8), allocatable :: me_lumped(:), me(:,:)
    integer, allocatable :: edof(:)
    real(8) ::  dens, thk

    select case( element(1)%id )
    case(2)
       mdim = 8
       kk = 2
    case default
       mdim = 12
       kk = 3
    end select

    allocate(me_lumped(mdim), me(mdim,mdim), edof(mdim))
    if (.not. allocated(mvec)) then
       allocate(mvec(neqn))
    end if
    mvec = 0d0

    i = 0
    do e=1,ne
       call get_edof(e,nen,xe,edof=edof)

       ! calculate element mass-matrix
       thk = mprop(element(e)%mat)%thk
       dens = mprop(element(e)%mat)%dens

       select case( element(e)%id )
       case( 2 ) ! plane42
          call plane42_me(xe,dens,thk, ng, me)
       case(3:5) ! laminate
          call mindlin42_me(xe,dens,mat_vec,thk,ng,me,2)
       case(6) ! plate, no piezo
          call mindlin42_me(xe,dens,mat_vec,thk,ng,me,1)
       end select
       call lump_mass(nen,kk,me,me_lumped)

       if (present(rho) .and. (element(e)%mat /=2)) then
          i = i +1
          if (rho(i)<0.1) then
             mvec(edof) = mvec(edof)+rho(i)**6 * me_lumped
          else
             mvec(edof) = mvec(edof)+rho(i)*me_lumped
          end if
       else
          mvec(edof) = mvec(edof)+me_lumped
       end if

    end do

    do i = 1, nb ! OUTVEC corrected for BC
       if  (bound(i,2) < 6) then! forskydning
          idof = kk*(bound(i,1)-1) + bound(i,2)
          mvec(idof) = 1d0
       end if
    end do

    do i=1,neqn
       if (mvec(i) == 0) then
          print*,'arggh__mvec(i)==0. i= ',i
          error stop
       end if
    end do

  end subroutine build_mvec


  subroutine mvec_mul(invec,outvec,mtype)
    ! This subroutine calculates the outvector. Since mvec is corrected, we dont need to correct outvec.

    use fedata
    use numeth

    real(8), intent(IN) :: invec(:)
    real(8), intent(INOUT) :: outvec(:)
    integer, intent(in) :: mtype
    integer :: e, i, j

    select case(mtype)
    case(1)!product of mass and invec
       outvec = mvec*invec
    case(2)!Solve equation. For ARPACK
       outvec = outvec/mvec
    case(3)!product of lumped mass and invec
       print*,'errro mvec_mul, fea'
       stop
    end select

    call correct_bound_vector(outvec)
    
  end subroutine mvec_mul


  subroutine lump_mass(nen,kk,me,me_lumped)

    use fedata

    real(8), intent(IN):: me(:,:)
    real(8), intent(OUT):: me_lumped(:)
    integer, intent(IN) :: nen, kk
    integer :: e, i,j! lumped_type
    real(8) :: me_diag_sum, me_tot

    !lumped_type = 1

    me_lumped = 0d0
    me_diag_sum  = 0d0
    me_tot = 0d0

    if (lumped_type ==1) then! Yuriy lumping, Add all row element to corressponding diagonal
       do j=1,kk*nen
          do i=1,kk*nen
             me_lumped(j) =me_lumped(j)+ me(j,i)
          end do
       end do
    elseif (lumped_type == 2) then ! HRZ - mass lumping, in COOK, s. 380
       do i=1,kk*nen
          me_lumped(i) = me(i,i)
       end do
       me_tot = sum(me)! m in cook
       me_diag_sum = sum(me_lumped) ! S in COOK, s. 380
       do i = 1, kk*nen
          me_lumped(i) = (me_tot/me_diag_sum)*me_lumped(i)
       end do
    end if


  end subroutine lump_mass


  subroutine eigen
    ! This subroutine calculates eigenvalues and eigenvectors

    !NB. VIRKER SIKKERT IKKE. BLIVER IKKE BRUGT 

    use numeth
    use fedata

    integer :: e ,i, j, p_ite, pmax, n, rneig
    real(8), dimension(neqn) ::  X_current, X_previous ,Y_current, Y_previous, differens
    integer, parameter :: neig = 10 ! Number of eigenvalues you are trying to find
    integer, dimension(neig) :: konv

    real(8), dimension(neqn,neig) :: D_eigen
    real(8), dimension(neig) :: eigenvalues, lambda
    real(8) :: norm, normdifferens, r, Z(neqn,neig-1), c, lambda0
    integer, parameter :: mdim = 8
    integer, dimension(mdim) :: idof

    pmax = 100
    D_eigen=0 
    Z=0
    lambda = 0
    lambda0 = 0 ! Corresponding to the 2th eigenfrequency

    call buildload ! Build load-vector

    call buildstiff_fea(0) ! Build stiffness matrix

    if (lambda0 /= 0) then ! only use SHIFT when lambda0 is given a nonzero value
       !call shift(lambda0)
    end if

    call enforce_fea ! Remove rigid body modes

    if (banded == 0) then ! factorization of K
       call factor(k)
    else !banded
       call bfactor(k)
    endif

	do i=1,neig
       X_previous = 1.! Initial guess for eigenvector

       do j = 1, nb ! correct X for BC
          idof = 2*(bound(i,1)-1) + bound(i,2)
          X_previous(idof) = bound(i, 3)
       end do

       call mmul_fea(X_previous,Y_current,3) ! Y corrected for BC!
       do j=1,i-1
          call mmul_fea(D_eigen(1:neqn,j),Z(1:neqn,j),3)
       end do
       do p_ite=1,pmax

          if (banded == 0) then
             Y_previous = Y_current ! Y_previous updated
             call solve(k, Y_previous) ! The solution X is written to Y_previous
             X_current = Y_previous   
             Y_previous = Y_current ! Y_previous is restored (new Y_current not jet calculated)
          else !banded
             Y_previous = Y_current
             call bsolve(k, Y_previous)  ! The solution X is written to Y_previous
             X_current = Y_previous
             Y_previous = Y_current ! Y_previous is restored
          end if

          c = 0
          do j=1,i-1
             do n=1, neqn ! Dot-product of vectors X and Z
                c = c+X_current(n)*Z(n,j)
             end do
             X_current = X_current - c*D_eigen(1:neqn,j)
          end do
          call mmul_fea(X_current,Y_current,3)
          r = 0
          do n=1, neqn ! Dot-product of vectors X and Y
             r = r+X_current(n)*Y_current(n)
          end do
          r = dsqrt(r)
          Y_current = real(Y_current/real(r))
          differens = X_current-X_previous
          normdifferens = 0
          norm =0
          do n=1, neqn !
             normdifferens = normdifferens+differens(n)*differens(n) ! Sums squared differens-components
             norm =norm + X_current(n)*X_current(n) ! Sums squared X-components
          end do
          if (dsqrt(normdifferens)/dsqrt(norm)<1e-12) then ! stop criteria
             exit
          end if
          X_previous = X_current ! update X_previous

       end do
       if( p_ite<pmax) then ! Checks convergence
          konv(i) = 1 ! Convergence
       else
          konv(i) = 0 ! No convergence
       end if

       D_eigen(1:neqn,i) = X_current/r
       do n=1, neqn !
          lambda(i) = lambda(i)+X_current(n)*Y_previous(n)/(r*r)
       end do
    end do

    do i = 1,neig
       lambda(i) = lambda(i)+lambda0 ! when SHIFT with lambda0=0 as default
       eigenvalues(i) = dsqrt(lambda(i))
       print*,'eigenval',eigenvalues(i)
	end do

    !call sturm(lambda(neig),rneig) ! real number of eigenvalues supposed to be found

	!call output_modal(D_eigen,eigenvalues,neig,rneig,konv,D_eigen)

 !$$$$$$     do i=1,neig
 !$$$$$$             D = D_eigen(1:neqn,i) ! Saves mode shapes for plot in 
 !$$$$$$             D=D/abs(maxval(D)) ! D is normlized with maxval(D) to get a visibly mode shape      
 !$$$$$$             call plot('deformed', 'gif', 'color', 'mode shapes')
 !$$$$$$             call plot('deformed', 'xwin', 'color', 'mode shapes')
 !$$$$$$             PAUSE ! enables the user to save .gif before overwritten
 !$$$$$$     end do

  end subroutine eigen

end module fea
