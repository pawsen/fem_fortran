MODULE thermal

  implicit none

  private
  public :: initial_t, therm_displ, recover_temp, enforce_term, buildtermK, buildtermload, &
       object_sens_t

contains

  !������������������������������������������������������������������
  ! Termiske rutiner. 1 frigedsgrad pr. knude., dvs 4 DOF pr element

  subroutine buildtermload(rho,rho_min)

    ! This subroutine builds the global thermal load vector.

    use fedata
    use plane41

    real(8), optional, dimension(ne), intent(in) :: rho
    real(8), optional, intent(in) :: rho_min
    integer :: i,n , e, nen, eface, nt, idof
    integer, parameter :: mdim = 8, ndim = 4
    integer, dimension(ndim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(ndim) :: rt
    real(8) :: fb, thk, t_inf, hconv
    integer :: ce(3*ne) ! ce skal v�re st�rre end ne, da hvert element jo i princippet kan have konvektion p� alle flader.
    ! og hver flade t�ller 1 i ce.

    ! Build global termisk load vector, se (12.2-5), s. 460 i COOK for mere
    r_t = 0.0d0
    nt = 0 ! number termisk-convection, t�llevariabel
    do i = 1, np	
       ! Build convection contribution [r_h]
       if (loads(i, 1) == 3) then !thermisk load specificeret i form af konvektion. Se linje ~540 i prosessor for v�rdier for load
          nt = nt + 1
          e = loads(i, 2)	! element SFE virker p�
          eface = loads(i, 3) !face id
          t_inf = loads(i, 4) ! temp af omgivende fluid, t_inf
          hconv = loads(i, 5) ! h for element
          ce(nt) = i			! element numbers with boundary convection
          !print*,'nt',nt

          nen = element(e)%numnode
          do n = 1, nen
             edof(n) = element(e)%ix(n)
             xe(2*n-1) = x(edof(n),1)
             xe(2*n  ) = x(edof(n),2)
          end do

          thk = mprop(element(e)%mat)%thk
          call plane41_rhe_t(xe, eface, t_inf, thk, hconv, rt)

          if (present(rho)) then
             r_t(edof) = r_t(edof) + rt!(rho_min+(1.0d0-rho_min)*rho(e)**penal)*rt
          else
             r_t(edof) = r_t(edof) + rt
          end if

       else if (loads(i, 1) == 5) then
          ! nodal flux
          idof  = loads(i,2)
          r_t(idof) = loads(i,4)

          ! Build heat flux contribution [r_B]
       else if (loads(i, 1) == 4) then
          e = loads(i, 2)	
          eface = loads(i, 3)
          fb = loads(i, 4) ! heat flux, [W/m^2]

          nen = element(e)%numnode
          do n = 1, nen
             edof(n) = element(e)%ix(n)
             xe(2*n-1) = x(edof(n),1)
             xe(2*n  ) = x(edof(n),2)
          end do
          thk = mprop(element(e)%mat)%thk
          call plane41_rbe_t(xe, eface, fb, thk, rt)

          if (present(rho)) then
             r_t(edof) = r_t(edof) + rt!(rho_min+(1.0d0-rho_min)*rho(e)**penal)*rt
          else
             r_t(edof) = r_t(edof) + rt
          end if

       end if



    end do


    ! Element faces with convection boundary, conv_elem = [element face hconv]
    ! Til brug ved samling af global convection matrix
    if (.not. allocated(conv_elem)) allocate(conv_elem(nt))
    conv_elem = ce(1:nt)

    ! Build volumetric heat generation contribution
    ! - Kan sammenlignes med accelerationskrafter der virker "lige" meget p� alle elementer
    if (qint /= 0.0d0) then	
       do e = 1, ne
          nen = element(e)%numnode
          do n = 1, nen
             edof(n) = element(e)%ix(n)
             xe(2*n-1) = x(edof(n),1)
             xe(2*n  ) = x(edof(n),2)
          end do

          thk = mprop(element(e)%mat)%thk
          call plane41_rqe_t(xe,ng, qint, thk, rt)

          if (present(rho)) then
             r_t(edof) = r_t(edof) + rt!(rho_min+(1.0d0-rho_min)*rho(e)**penal)*rt
          else
             r_t(edof) = r_t(edof) + rt
          end if

       end do
    end if


  end subroutine buildtermload


  subroutine buildtermK(flag,rho_min, rho)

    ! This subroutine builds the combination of the globalconductivity matrix [k] and the global boundary
    ! convection matrix [h] from the local element matrices, jf (12.2-5) i cook

    use fedata
    use plane41

    integer, INTENT(IN) :: flag
    real(8), optional, INTENT(IN) :: rho_min, rho(ne)
    integer :: e, i, j,n, eface, g, nt
    integer :: nen
    integer, parameter :: mdim = 8 , ndim = 4
    integer, dimension(ndim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(ndim, ndim) :: ke_t, he_t
    real(8) :: kcond, hconv, thk

    ! Reset thermal "stiffness" matrix
    if ( banded /= 2) then
       k_t = 0.0d0
    end if

    ! Finds the global conductivity matrix [k]
    do e = 1, ne
       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do i = 1, nen
          edof(i) = element(e)%ix(i)
          xe(2*i-1) = x(edof(i),1)
          xe(2*i  ) = x(edof(i),2) 
       end do

       ! hvis alle elementer er ens, skal nedenst�ende kun k�res �n gang
       if ( (flag == 0) .or. (e == 1) ) then
          ! Gather material properties and find element matrices
          select case( element(e)%id )
          case( 1 )
             print*,'ERROR in buildtermK: Thermal conductivity not implemented for truss elements.'
          case( 2 )
             ! calculates the element conductivity matrix [k_e]
             thk = mprop(element(e)%mat)%thk
             kcond = mprop(element(e)%mat)%kcond
             call plane41_ke_t(xe,ng, kcond, thk, ke_t)
          end select
       end if

       ! Assemble into global matrix

       if ( banded == 0) then
          do i = 1, nen
             do j = 1, nen
               	if (present(rho)) then
                   k_t(edof(i), edof(j)) = k_t(edof(i), edof(j)) + (rho_min+(1.0d0-rho_min)*rho(e)**penal)*ke_t(i,j)
                else
                   k_t(edof(i), edof(j)) = k_t(edof(i), edof(j)) + ke_t(i, j)
               	end if
             end do
          end do
       elseif(banded == 1) then
          do i = 1, nen
             do j = 1, nen
                g=edof(i)-edof(j)
                if (g >= 0) then ! brug kun elementer der st�r under eller p� diagonalen i global K+
                   if (present(rho)) then !topOpt
                      k_t(g+1,edof(j)) = k_t(g+1,edof(j)) + (rho_min+(1.0d0-rho_min)*rho(e)**penal)*ke_t(i,j)
                      !print*,'hej'
                   else
                      k_t(g+1,edof(j)) = k_t(g+1,edof(j)) + ke_t(i,j)
                   end if
                end if
             end do
          end do
       end if
    end do


    ! Finds the global boundary convection matrix [h]
    nt = size(conv_elem,1) ! conv_elem laves i buildthermload
    do n = 1, nt
       i = conv_elem(n)
       e = loads(i, 2)	! element SFE virker p�
       thk = mprop(element(e)%mat)%thk
       eface = loads(i, 3) !face id
       hconv = loads(i, 5) 	! h for element

       ! Find coordinates and degrees of freedom
       nen = element(e)%numnode
       do j = 1, nen
          edof(j) = element(e)%ix(j)
          xe(2*j-1) = x(edof(j),1)
          xe(2*j  ) = x(edof(j),2)
       end do


       ! calculates the element boundary convection matrix [h_e]
       call plane41_he_t(xe, eface, hconv, thk, he_t)
       ! Assemble into global matrix
       if (banded == 0) then
          do i = 1, nen
             do j = 1, nen
                if (present(rho)) then
                   k_t(edof(i), edof(j)) = k_t(edof(i), edof(j)) + (rho_min+(1.0d0-rho_min)*rho(e)**penal)*he_t(i,j)
                else
                   k_t(edof(i), edof(j)) = k_t(edof(i), edof(j)) + he_t(i,j)
               	end if
             end do
          end do
       elseif(banded == 1) then
          do i = 1, nen
             do j = 1, nen
                g=edof(i)-edof(j)
                if (g >= 0) then
                   if (present(rho)) then !topOpt
                      k_t(g+1,edof(j)) = k_t(g+1,edof(j)) + (rho_min+(1.0d0-rho_min)*rho(e)**penal)*he_t(i,j)
                   else
                      k_t(g+1,edof(j)) = k_t(g+1,edof(j)) + he_t(i,j)
                   end if
                end if
             end do
          end do
       end if

    end do


  end subroutine buildtermK

  subroutine enforce_term!(r_t)

    ! This subroutine enforces the thermal boundary conditions.
    ! Se kap 2.7 i COOK

    use fedata

    !   real(8), dimension(neqn), intent(INOUT) :: r_t
    integer :: i, j, idof
    real(8), dimension(nn) :: kl_t
    integer :: bw_t
    bw_t = size(k_t,1)
    kl_t=0.0d0

    if (banded == 0) then
       if (.not. penalty) then
          do i = 1, nb
             if (bound(i,2) == 3) then ! temp prescribed
                idof = bound(i,1) ! DOF = node_id, da der kun er �n dof pr knude
                r_t = r_t - k_t(1:nn, idof) * bound(i, 3)
                r_t(idof) = bound(i,3)
                k_t(1:nn, idof) = 0.0d0
                k_t(idof, 1:nn) = 0.0d0
                k_t(idof, idof) = 1.0d0
             end if
          end do
       else
          print*,'ERROR in enforce_term: Penalty method not implemented for thermal loads. - Men s� g�r det dog!!'
       endif
    elseif(banded == 1) then ! banded
       if (.not. penalty) then
          do i = 1, nb
             if (bound(i,2) == 3) then ! temp p� rand
                idof = bound(i,1)
                ! nedenst�ende(kl_t) finder den v�rdi der skal l�gges til load-vektoren ved imposing af BC, jf (2.7-6) i COOK
                kl_t=0.0d0
                do j = 1, idof-1
                   if (idof-j<bw_t) then
                      kl_t(j) = k_t(idof-j+1 , j)
                   end if
                end do
                do j = idof , nn
                   if (j-idof+1<=bw_t) then
                      kl_t(j) = k_t(j-idof+1 , idof)
                   end if
                end do

                r_t = r_t - kl_t * bound(i, 3)
                r_t(idof) = bound(i, 3)

                k_t(1, idof) = 1.0d0
                k_t(2:bw_t,idof) = 0.0d0
                do j=2,bw_t
                   if ((idof-(j-1)) > 0) then
                      k_t(j,idof-(j-1)) = 0.0d0 ! s. fig 6.2 i noter(s. 42). R�kke i 'full matrix' bliver til sk� linje i 'banded matrix'.
                   end if
                end do

             end if
          end do
       else
          print*,'ERROR in enforce_term: Penalty method not implemented for thermal loads.'; stop 
       end if


    end if


  end subroutine enforce_term

  subroutine recover_temp!(t)!,tn)

    ! This subroutine recovers the temperature for each element 

    use fedata
    use plane41

    !   real(8), dimension(ne), intent(OUT) :: t
    !   real(8), dimension(neqn), intent(IN) ::  tn
    integer :: e, i, nen
    integer, parameter :: mdim = 8, ndim = 4
    integer :: edof(ndim)
    real(8) :: xe(mdim), tne(ndim), te


    do e = 1, ne
       ! Find coordinates etc...
       nen = element(e)%numnode
       do i = 1,nen
          edof(i) = element(e)%ix(i)
          xe(2*i-1) = x(edof(i),1)
          xe(2*i  ) = x(edof(i),2)
       end do

       ! Nodal temperatures for element
       tne = tn(edof)

       ! Find temperature
       select case( element(e)%id )
       case( 1 )
          print*,'ERROR in recover_temp: Thermal conductivity not implemented for truss elements.'
       case( 2 )
          call plane41_temp(xe, tne, te)
          t_elem(e) = te
       end select
    end do

  end subroutine recover_temp

  subroutine initial_t

    ! This subroutine find the bandwidth and allocate thermal stiffness matrix
    use fedata
    use plane42


    integer :: i ,j, e, nen, bwe
    neqn = nn

    if (banded == 0) then
       allocate (k_t(neqn, neqn))
    elseif (banded == 1) then !allocate k=banded
       bw = 0
       do e = 1 , ne
          nen = element(e)%numnode
          do i = 1 , (nen-1)
             do j = i+1 , nen
                bwe = ABS(element(e)%ix(i) - element(e)%ix(j))+1
                if (bwe>bw) then
                   bw = bwe
                end if
             end do
          end do
       end do
       allocate (k_t(bw,neqn))
    end if

    allocate (tn(neqn), r_t(neqn))

  end subroutine initial_t


  subroutine therm_displ(flag, rho_min,rho)

    ! This subroutine calculates displacements

    use numeth
    use processor
    use plot_routiner
    use fedata

    integer, INTENT(IN) :: flag
    real(8), optional, INTENT(IN) :: rho_min, rho(ne)
    real(8) :: tstart, plotval(ne)

    ! r_t = 'termisk load' vektor
    ! t = element temperatur
    ! tn = node temperatur

    if (antype == 'THERMAL' .or. antype == 'COUPLED') then
       call buildtermload
    end if
    if (present(rho) ) then
       call buildtermK(flag,rho_min,rho)
    else! Build thermal stiffness matrix
       call buildtermK(flag)
    endif

    ! Enforce thermal boundary conditions
    call enforce_term

    tn = r_t
    if (banded == 1) then
       call factor(k_t)
       call solve(k_t, tn)
    elseif (banded == 0) then    
       call bfactor(k_t)
       call bsolve(k_t, tn)
    endif

    ! lig start-temperatur for alle elementer til resultatet. IKKE smart!
    tstart = mprop(element(1)%mat)%tstart
    ! tn is nodal temp
    tn = tstart + tn

    if (antype == 'TOPTHERMSTRUCT') then

       call recover_temp!(t) ! Calculate element center temperatures from node temperatures
    elseif (antype == 'THERMAL' .or. antype == 'COUPLED') then
       allocate(t_elem(ne))
       call recover_temp!(t) ! Calculate element center temperatures from node temperatures
       call output_term(t_elem,tn) ! skriver temp til fil
       ! Plot element values
       plotval = t_elem
       call plot('elements', 'xwin', 'color', ' ', plotval) ! liux_fejl
       !call plot('elements', 'matlab', 'color', 'spaendingsintensitet', plotval)
    endif

  end subroutine therm_displ


  subroutine object_sens_t(flag,rho,rho_min,f,df)
    ! d.03_05_2011
    ! This subroutine calculate the objectfunction and sensitivities for pure thermal
    ! optimization

    use fedata
    use plane41

    integer, INTENT(IN) :: flag
    real(8), dimension(:), intent(in) :: rho
    real(8), intent(out) :: f
    real(8), dimension(:), intent(out) :: df
    real(8), intent(IN) :: rho_min

    integer :: e, i, j, nen, n_conv
    integer, parameter :: mdim = 8, ndim = 4
    integer :: edof(ndim), eface
    real(8) :: xe(mdim),fac, hconv, kcond, thk
    real(8), dimension(ndim, ndim) :: ke_t, he_t



    f = 0.0d0       ! object-funktion
    he_t = 0.0d0    ! element convection-matrix
    ke_t = 0.0d0    ! element conductivity-matrix


    do e = 1, ne
       nen = element(e)%numnode
       do i = 1,nen
          edof(i) = element(e)%ix(i)
          xe(2*i-1) = x(edof(i),1)
          xe(2*i  ) = x(edof(i),2)
       end do

       if ( (flag == 0) .or. (e == 1) ) then
          ! Local thermal stiffness matrix
          thk = mprop(element(e)%mat)%thk
          kcond = mprop(element(e)%mat)%kcond
          call plane41_ke_t(xe,ng, kcond, thk, ke_t)
       end if

       he_t = 0.0d0  
       n_conv = size(conv_elem,1) !conv_elem indeholder index(i forhold til load-vektor) p� de elementer med konvektion
       if (n_conv > 1) then ! m�ske overfl�digt?
          do j = 1, n_conv
             i = conv_elem(j)
             if (e == loads(i, 2) ) then ! loads(i, 2) = element convection virker p�
                eface = loads(i, 3) !face id
                hconv = loads(i, 5) ! h for element

                ! calculates the element boundary convection matrix [h_e]
                if (eface == 5) then
                   ! Konvection fra top/bottem
                else
                   ! calculates the element boundary convection matrix [h_e]
                   call plane41_he_t(xe, eface, hconv, thk, he_t)
                   exit ! stopper loop og fors�tter loppet over e
                end if
             endif
          end do
       end if

       ! Compute compliance
       f = f + DOT_PRODUCT( MATMUL(( (rho_min+(1.0d0-rho_min)*rho(e)**penal)* ke_t + he_t), tn(edof)), tn(edof) )

       !Sensitivity 
       df(e) = -(1.0d0-rho_min)*penal*(rho(e)**(penal-1))*DOT_PRODUCT(MATMUL(ke_t + he_t,tn(edof)),tn(edof)) ! NB FJERNET MINUS FORAN UDTRYK d_04_05_11
    end do

  end subroutine object_sens_t

END MODULE thermal
