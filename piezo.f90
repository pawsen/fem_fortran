module piezo

  implicit none

  private :: buildstiff_piezo, rb_init_piezo, recover_potential, enforce_piezo
  public :: initial_piezo, displ_piezo

  public :: buildstiff_eigenvalue_piezo

  private :: buildS_piezo,get_element_matrix_piezo

contains

  subroutine initial_piezo

    ! This subroutine is mainly used to allocate vectors and matrices
    use fedata
    use file_init

    !Yesterday I spent (+ 2234.34 3423.4 (* 12.0 12.2))

    !    logical, INTENT(IN) :: harmonic
    integer :: nzz,dof, rand_int(4), int_dummy(11)
    real(8) ::  rand_val(4), real_dummy(7)

    if ((elem_type /= 'PLATE_GMSH') .and. &
         (elem_type /= 'PLANE_GMSH')) then
       ! add bound from potential
       call parameter_input(real_dummy,int_dummy,rand_int,rand_val)
       call rb_init_piezo(rand_int,rand_val)
       print*
       print*,'automatic added potential & RB'
    end if

    if (elem_type == 'PLATE') then
       neqn = 3*nn
       dof = 12 ! mechanical dof for one element
       nzz = (dof*dof+4*4)*ne+2*4*dof*ne+2*nb
       !nzz = (dof*dof)*ne+2*nb

    else if ((elem_type == 'PLATE_GMSH' .and. eigenvalue%calc) .or. &
        ( antype == 'TOPSTRUCT_EIGEN') ) then
       neqn = 3*nn
       dof = 12   
       nzz = 12*12*ne
    else if (elem_type == 'PLATE_GMSH' .and. eigenvalue%calc .eqv. .false.)then
       neqn = 3*nn
       dof = 12 
       nzz = (dof*dof+4*4)*ne_ptz+2*4*dof*ne_ptz+& ! ptz-material
            (12*12+4*4)*(ne-ne_ptz) + & ! non-ptz material
            2*nb !RB
       !ptz = 8192
       !nb =  114
       !pcb = 158976
       !tot = 167282
    else
       dof = 8
       neqn = 2*nn
       nzz = (8*8+4*4)*ne+2*4*8*ne+2*nb

    end if

    if (harmonic) then
       allocate(sKZ(nzz), dZ(neqn+nn+nb),pZ(neqn+nn+nb))
       sKZ = CMPLX(0d0,0d0)
       dZ = CMPLX(0d0,0d0)
       pZ = CMPLX(0d0,0d0)
    else
       allocate(sK(nzz), d(neqn+nn+nb),p(neqn+nn+nb))
       sK = 0d0
       d = 0d0
       p = 0d0
    end if

    allocate(iK(nzz),jK(nzz))
    allocate (strain(ne, 5), stress(ne, 5))
    iK = 0
    jk = 0


  end subroutine initial_piezo

  subroutine displ_piezo(flag,rho,rho_min)

    ! This subroutine calculates displacements
    use numeth
    use fedata
    use processor
    use plot_routiner
    use solve_handle_real
    use solve_handle_complex
    use file_init
    use plate
    use fea ! to build mechanical load for beam
    use exodus


    integer, INTENT(IN) :: flag
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(INOUT) :: rho_min

    integer :: e, rand_int(4), int_dummy(11), flagger
    real(8) :: plotval(ne), rand_val(4), real_dummy(7), omega,frekvens
    !    logical :: harmonic = .false.!.true.
    !    logical :: harmonic = .true.


    flagger = 0
    if (harmonic) then
       frekvens = mat_vec(20)
       !frekvens = 4.143E+4!5.94E+4!4.624E+4!
       omega = 2*pi*frekvens
       print*,'frekvens= ', frekvens
       call buildS_piezo(flagger,omega)
       call enforce_piezo
       call mumps_solve_complex

       call output_deformed('Deformeret','potentiale',(/0d0/),dZ(neqn+1:neqn+nn))
    else

       if (element(1)%id  == 2) then ! plane42
          call buildload
       else
          call buildload_plate
       end if
       call buildstiff_piezo(flagger)
       ! Remove rigid body modes
       call enforce_piezo

       call mumps_init_real
       call mumps_solve_real(6)
       call mumps_finalize_real

       call output_deformed('Deformeret','potentiale',d(neqn+1:neqn+nn))
    end if

    call exodus_write

    call output(2)


    !    call output
    if (present(rho)) then
       !call plot('elements', 'xwin', 'color', 'spaendingsintensitet', plotval)
       call output_deformed('deformeret','densitet',rho)
       call output_elements('', plotval)
       call plot('elements', 'xwin', 'color', 'spaendingsintensitet', plotval)
       call plot('deformed', 'xwin', 'color', 'deformeret',rho)
    endif

  end subroutine displ_piezo

  subroutine buildstiff_piezo(flag,rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata

    integer, INTENT(INOUT) :: flag !! NB!!!!
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(IN) :: rho_min

    integer :: e,i,j,idof,mdim,kk,nen,ii
    integer, parameter :: ndim = 4

    integer, dimension(ndim) :: edof41
    integer,dimension(:), allocatable :: edof
    real(8), dimension(8) :: xe
    real(8), dimension(:,:), allocatable :: ke_stiff, me
    real(8), dimension(:,:), allocatable :: ke_piezo, ke_dielectric, ke_piezo_transpose

    select case( element(1)%id )
    case(2)
       mdim = 8
       kk = 2
    case default
       mdim = 12
       kk = 3
    end select
    allocate(ke_stiff(mdim,mdim), me(mdim,mdim), edof(mdim))
    allocate(ke_piezo(mdim,ndim), ke_dielectric(ndim,ndim), ke_piezo_transpose(ndim,mdim) )

    ke_stiff = 0d0
    ke_piezo = 0d0
    ke_dielectric = 0d0
    ke_piezo_transpose = 0d0

    flag = 0
    ii = 0
    do e = 1, ne

       call get_element_matrix_piezo(flag,e,edof,edof41,nen,&
            ke_stiff,ke_piezo,ke_piezo_transpose,ke_dielectric)

       ! Assemble into global matrix
       ! Stiffness
       do i = 1, kk*nen
          do j =1, kk*nen
             ii = ii+1
             iK(ii) = edof(i)
             jK(ii) = edof(j)
             sK(ii) = ke_stiff(i,j)
          end do
       end do

       ! Dielectric
       do i = 1, nen
          do j = 1, nen
             ii = ii+1
             iK(ii) = neqn+ edof41(i)
             jK(ii) = neqn+ edof41(j)
             sK(ii) = ke_dielectric(i,j)
          end do
       end do

       ! KOBLINGS-MATRICER. ONLY for ptz-materral
       if ( element(e)%id /= 6 ) then
!!$       if (e==1) then
!!$          print*,'test'
!!$          do i=1,12
!!$             print*,(ke_piezo(i,j),j=1,4)
!!$             !write (*,'(8(D10.4 2x))') (real(ke_stiff(i,j)),j=1,8)
!!$          end do
!!$       end if

          ! piezo coupling effect
          do i = 1, kk*nen
             do j = 1, nen
                !matrix to the right            
                ii = ii+1
                iK(ii) = edof(i)
                jK(ii) = neqn+ edof41(j)
                sK(ii) = ke_piezo(i,j)
             end do
          end do

          do i = 1, nen
             do j = 1, kk*nen
                !nederste matrix
                ii = ii+1
                iK(ii) = neqn+ edof41(i)
                jK(ii) = edof(j)
                sK(ii) = ke_piezo_transpose(i,j)
             end do
          end do

       end if
    end do

    ! add lagrangian multipliers
    !stiffness:
    do i=1,nb
       if (bound(i,2) < 6) then! forskydning
          idof = kk*(bound(i,1)-1) + bound(i,2)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+nn+i
          jK(ii) = idof
          sK(ii ) = 1d0
          ! matrix to the right
          ii = ii+1
          iK(ii) = idof
          jK(ii) = neqn+nn+i
          sK(ii ) = 1d0
       else !potentiale
          idof = bound(i,1)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+nn+i
          jK(ii) = neqn + idof
          sK(ii) = 1d0
          ! matrix to the right
          ii = ii + 1
          iK(ii) = neqn + idof
          jK(ii) = neqn+nn + i
          sK(ii) = 1d0
          !print*,'potential RB'
       end if
    end do

  end subroutine buildstiff_piezo


  subroutine get_element_matrix_piezo(flag,e,edof,edof41,nen,&
       ke_stiff,ke_piezo,ke_piezo_transpose,ke_dielectric, &
       me)

    use fedata
    use mindlin42
    use plane42
    use plane41
    use numeth

    integer, intent(in) :: flag, e
    integer, intent(out) :: edof(:), edof41(:), nen
    !real(8), intent(out) :: xe(:)
    real(8), intent(out), dimension(:,:) :: ke_stiff,ke_piezo,ke_piezo_transpose,ke_dielectric
    real(8), intent(out), dimension(:,:), optional :: me

    real(8), dimension(8) :: xe
    real(8) :: young, nu, area, dens, thk, shear, phi, ep
    integer :: i

    ! Find coordinates and degrees of freedom
    call get_edof(e,nen,xe,edof41,edof)

    ! hvis alle elementer er ens, skal nedenst�ende kun k�res �n gang
    if ( ((flag == 0) .or. (e == 1))) then
       ! Gather material properties and find element stiffness matrix
       young = mprop(element(e)%mat)%young
       thk = mprop(element(e)%mat)%thk
       nu = mprop(element(e)%mat)%nu
       shear = mprop(element(e)%mat)%shear
       dens  = mprop(element(e)%mat)%dens
       phi = 0!

       select case( element(e)%id )
       case( 2 )
          ! calculate plane42 element stiffness matrix here
          !call plane42_ke(xe,young, nu, thk, ng, ke_stiff)
          phi = 0 ! pi/4d0
          call plane42_ke_piezo(xe,mat_vec, young, nu, thk, ng, ke_stiff,phi)
          call plane42_piezo(xe, mat_vec, thk,ng, ke_piezo,phi)
          call plane41_dielectric(xe,mat_vec, thk,ng, ke_dielectric,phi)
          ke_piezo_transpose = TRANSPOSE(ke_piezo)
          !call plane42_me(xe, dens, ng, me, thk)
          !             ce = alpha1*me+beta1*ke_stiff
       case(3) ! plate
          phi =  0! pi!pi/4d0
          ep = mprop(element(e)%mat)%ep

          call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke_stiff,2,mat_vec,phi)
          call mindlin42_piezo(xe, mat_vec, thk,ng,ke_piezo,phi)
          call mindlin42_dielectric(xe,mat_vec,ep, thk,ng, ke_dielectric,2,phi)
          ke_piezo_transpose = TRANSPOSE(ke_piezo)

          if (present(me)) then
             call mindlin42_me(xe,dens,mat_vec,thk,ng,me,2)
          end if
       case(4) ! polariseret piezo, ringmotor
          ep = mprop(element(e)%mat)%ep
          phi = element(e)%ptz
          call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke_stiff,2,mat_vec,phi)
          call mindlin42_piezo(xe, mat_vec, thk,ng,ke_piezo,phi)
          call mindlin42_dielectric(xe,mat_vec,ep, thk,ng, ke_dielectric,2,phi)
          ke_piezo_transpose = TRANSPOSE(ke_piezo)

          if (present(me)) then
             call mindlin42_me(xe,dens,mat_vec,thk,ng,me,2)
          end if
       case(5) ! polariseret piezo, ringmotor
          ep = mprop(element(e)%mat)%ep
          phi = element(e)%ptz + pi ! polariseret modsat
          call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke_stiff,2,mat_vec,phi)
          call mindlin42_piezo(xe, mat_vec, thk,ng,ke_piezo,phi)
          call mindlin42_dielectric(xe,mat_vec,ep, thk,ng, ke_dielectric,2,phi)
          ke_piezo_transpose = TRANSPOSE(ke_piezo)

          if (present(me)) then
             call mindlin42_me(xe,dens,mat_vec,thk,ng,me,2)
          end if
       case(6) !rent pcb, ringmotor
          ep = mprop(element(e)%mat)%ep
          call mindlin42_dielectric(xe,mat_vec,ep, thk,ng, ke_dielectric,1)
          call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke_stiff,1,mat_vec) !type = 1
          if (present(me)) then
             call mindlin42_me(xe,dens,mat_vec,thk,ng,me,1)
          end if
       end select
    end if

  end subroutine get_element_matrix_piezo

  subroutine buildstiff_eigenvalue_piezo(rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata
    use numeth

    real(8), optional, intent(in) :: rho(:), rho_min

    integer :: e,i,j,idof,nen,ii,jj, flag
    integer, parameter ::kk=3, mdim = 12 , ndim = 4 ! mindlin
    !    integer, parameter ::kk = 2, mdim = 8 , ndim = 4 ! bj�lke

    integer, dimension(ndim) :: edof41
    integer, dimension(mdim) :: edof
    real(8), dimension(8) :: xe
    real(8), dimension(mdim, mdim) :: ke_stiff, me
    real(8) :: ke_piezo(mdim,ndim), ke_dielectric(ndim,ndim), ke_piezo_transpose(ndim,mdim)
    real(8) :: helpmat0(4,4),helpmat1(12,4),helpmat2(12,12)


    ke_stiff = 0d0
    ke_piezo = 0d0
    ke_dielectric = 0d0
    ke_piezo_transpose = 0d0

    helpmat2 = 0d0
    
    flag = 0
    ii = 0
    jj = 0
    do e = 1, ne
       
       call get_element_matrix_piezo(flag,e,edof,edof41,nen,&
            ke_stiff,ke_piezo,ke_piezo_transpose,ke_dielectric)

       if ( element(e)%id /= 6 ) then
          ! helpmat0 = inv_matrix(ke_dielectric)
          ! helpmat1 = MATMUL(ke_piezo, helpmat0)
          ! helpmat2 = MATMUL(helpmat1,ke_piezo_transpose)
       else
          helpmat2 = 0d0
       end if

       ! Assemble into global matrix
       ! Stiffness
       if (present(rho) .and. (element(e)%mat /=2) ) then
          jj = jj +1
          do i = 1, kk*nen
             do j =1, kk*nen
                ii = ii+1
                iK(ii) = edof(i)
                jK(ii) = edof(j)
                sK(ii) = (rho_min+(1d0-rho_min)*rho(jj)**penal)* ke_stiff(i,j) - &
                     0*helpmat2(i,j)
             end do
          end do
       else
          do i = 1, kk*nen
             do j =1, kk*nen
                ii = ii+1
                iK(ii) = edof(i)
                jK(ii) = edof(j)
                sK(ii) = ke_stiff(i,j) -0*helpmat2(i,j)
             end do
          end do
       end if

       !sinve we do sparse matrix-vector multiplication, wo dont need lagrangian RB.
    end do

  end subroutine buildstiff_eigenvalue_piezo


  subroutine recover_potential(node_vec,elem_vec)

    ! This subroutine calculates the potential in each element from the node values
    ! Only for element having 1 dof pr node.
    use fedata
    use plane41

    real(8), dimension(ne), intent(OUT) :: elem_vec
    real(8), dimension(nn), intent(IN) ::  node_vec
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

       ! Nodal values for element
       tne = node_vec(edof)

       ! Find element value
       select case( element(e)%id )
       case( 1 )
          print*,'ERROR in recover_potential: Thermal conductivity not implemented for truss elements.'
          error stop
       case( 2 )
          call plane41_temp(xe, tne, te)
          elem_vec(e) = te
       end select
    end do

  end subroutine recover_potential


  subroutine rb_init_piezo(rand_int,rand_val)
    !  This subroutine find the nodes in the end and the beginning af a beam.

    use fedata
    use plot_routiner

    integer, intent(IN) :: rand_int(:) ! m� ikke hedde rand, da rand er initialiseret i fedata, DOH!!!
    real(8), intent(IN) :: rand_val(:)
    integer :: i, e, nen, nb_temp
    integer :: number_elements_x, number_elements_y, number_nodes_y, number_nodes_x, ny,nx
    real(8) :: xmin, xmax, ymin, ymax, lx, ly, length_first_element, hight_first_element
    integer(8), dimension(:,:),allocatable :: bound_temp

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
    number_elements_x = NINT(lx/length_first_element) ! Antal elemnter i x-retningen ! runder op til n�rmeste integer. Hvis nint ikke medtages, rundes der forkert af
    number_nodes_x = number_elements_x+1 ! Alternativt nn/number_nodes_y
    number_elements_y = NINT(ly/hight_first_element) ! Antal elemnter i y-retningen
    number_nodes_y = number_elements_y+1 !Antal knuder i y-retningen

    ny = number_nodes_y
    nx = number_nodes_x

    allocate (bound_temp(nb, 3))
    bound_temp = bound
    DEALLOCATE(bound)

    nb_temp = nb
    do e=1,size(rand_int)! find number of bounds including potential
       ! remember this is potential, so there's only 1 DOF pr node, eg. we need only node numbers
       select case(rand_int(e))
       case(1) ! bottom
          nb_temp = nb_temp + nx
       case(2)
          nb_temp = nb_temp + ny
       case(3)
          nb_temp = nb_temp + nx
       case(4)
          nb_temp = nb_temp + ny
       end select
    end do
    allocate (bound(nb_temp,3))
    bound(1:nb,:) = bound_temp

    do e=1,size(rand_int) ! set RB for potential
       select case(rand_int(e))
       case(1) ! bottom
          do i=1,nx
             nb = nb+1
             bound(nb,2) = 10 !identifyer for potential
             bound(nb,1) = 1+(i-1)*ny ! node number
             bound(nb,3) = rand_val(e)! Value of rb, eg V=0 or V=10 or..
          end do
       case(2) ! right side of structure
          do i=1,ny
             nb = nb+1
             bound(nb,2) = 10
             bound(nb,1) = nn-ny+i
             bound(nb,3) = rand_val(e)
          end do
       case(3) ! top
          do i=1,nx
             nb = nb+1
             bound(nb,2) = 10
             bound(nb,1) = i*ny
             bound(nb,3) = rand_val(e)
          end do
       case(4) ! left side of structure
          do i=1,ny
             nb = nb+1
             bound(nb,2) = 10
             bound(nb,1) = i
             bound(nb,3) = rand_val(e)
          end do
       end select
    end do


  end subroutine rb_init_piezo

  subroutine buildS_piezo(flag,omega,rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata

    integer, INTENT(IN) :: flag
    real(8), INTENT(IN) :: omega ! exciationsfrekvens
    real(8), optional, dimension(:), INTENT(IN) :: rho
    real(8), optional, INTENT(IN) :: rho_min

    integer :: e, i, j, ii, type
    integer :: nen, idof
    integer, parameter ::kk=3, mdim = 12 , ndim = 4 ! mindlin
    !    integer, parameter ::kk = 2, mdim = 8 , ndim = 4 ! bj�lke

    integer, dimension(ndim) :: edof41
    integer, dimension(mdim) :: edof
    real(8), dimension(mdim) :: xe
    real(8), dimension(mdim, mdim) :: ke_stiff, me, ce
    real(8) :: ke_piezo(mdim,ndim), ke_dielectric(ndim,ndim), ke_piezo_transpose(ndim,mdim)

    real(8) :: alpha1, beta1

    ke_stiff = 0d0
    ke_piezo = 0d0
    ke_dielectric = 0d0
    ke_piezo_transpose = 0d0
    me = 0d0
    ce = 0.0

    ii = 0
    do e = 1, ne

       call get_element_matrix_piezo(flag,e,edof,edof41,nen,&
            ke_stiff,ke_piezo,ke_piezo_transpose,ke_dielectric, &
            me)
       !ce = alpha1*me+beta1*ke_stiff

       ! Assemble into global matrix
       ! Stiffness
       do i = 1, kk*nen
          do j =1, kk*nen
             ii = ii+1
             iK(ii) = edof(i)
             jK(ii) = edof(j)
             sKZ(ii) = CMPLX( ke_stiff(i,j) - omega**2 * me(i,j) , omega*ce(i,j) )
          end do
       end do
       ! Dielectric
       do i = 1, nen
          do j = 1, nen
             ii = ii+1
             iK(ii) = neqn+ edof41(i)
             jK(ii) = neqn+ edof41(j)
             sKZ(ii) =CMPLX( ke_dielectric(i,j) , 0.0)
          end do
       end do

       ! KOBLINGS-MATRICER. ONLY for ptz-materral
       if ( element(e)%id /= 6 ) then
          ! piezo coupling effect
          do i = 1, kk*nen
             do j = 1, nen
                !matrix to the right            
                ii = ii+1
                iK(ii) = edof(i)
                jK(ii) = neqn+ edof41(j)
                sKZ(ii) = CMPLX( ke_piezo(i,j) , 0.0 )
             end do
          end do
          do i = 1, nen
             do j = 1, kk*nen
                !nederste matrix
                ii = ii+1
                iK(ii) = neqn+ edof41(i)
                jK(ii) = edof(j)
                sKZ(ii) = CMPLX( ke_piezo_transpose(i,j) ,0.0 )
             end do
          end do
       end if
    end do


    ! add lagrangian multipliers
    !stiffness:
    do i=1,nb
       if (bound(i,2) < 6) then! forskydning
          idof = kk*(bound(i,1)-1) + bound(i,2)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+nn+i
          jK(ii) = idof
          sKZ(ii ) = CMPLX( 1d0, 0d0)
          ! matrix to the right
          ii = ii+1
          iK(ii) = idof
          jK(ii) = neqn+nn+i
          sKZ(ii ) = CMPLX( 1d0, 0d0)
       else
          idof = bound(i,1)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+nn+i
          jK(ii) = neqn + idof
          select case(INT(bound(i,2)))
          case(10) ! gnd
             sKZ(ii) = CMPLX( 1d0, 0d0)
          case(11) ! cos
             sKZ(ii) = CMPLX( 1d0, 0d0)
          case(12) ! sin
             sKZ(ii) = CMPLX( 1d0, 0d0)
          end select

          ! matrix to the right
          ii = ii + 1
          iK(ii) = neqn + idof
          jK(ii) = neqn+nn + i
          select case(INT(bound(i,2)))
          case(10) ! gnd
             sKZ(ii) = CMPLX( 1d0, 0d0)
          case(11) ! cos
             sKZ(ii) = CMPLX( 1d0, 0d0)
          case(12) ! sin
             sKZ(ii) = CMPLX( 1d0, 0d0)
          end select
       end if
    end do


  end subroutine buildS_piezo

  subroutine enforce_piezo

    ! This subroutine enforces the support boundary conditions.
    use fedata

    integer :: i, j, idof, nn_buf

    !nn_buf = nn
    !nn = 0 

    !sparse. Add constraint to the load vector
    if (harmonic) then
       do i = 1,nb
          select case(INT(bound(i,2)))
          case(10) ! gnd
             pZ(neqn+nn+i) = CMPLX(bound(i, 3),0d0)
          case(11) ! cos
             pZ(neqn+nn+i) = CMPLX(bound(i, 3),0d0)
          case(12) ! sin
             !(a+ib)*(c+id) = (a*c-b*d) + i(a*d + b*c)
             ! S� for at f� en komplex p�virkning s�ttes a=1,b=0 & c=0,d=-pot
             pZ(neqn+nn+i) = CMPLX(0d0,-bound(i, 3))

          case default ! resten, incl geometrisk constraint
             pZ(neqn+nn+i) = CMPLX(bound(i, 3),0d0)
          end select
       end do
    else
       do i = 1,nb
          p(neqn+nn+i) = bound(i, 3)
       end do
    end if
    !nn = nn_buf
  end subroutine enforce_piezo


  SUBROUTINE sweep_piezo(flag,rho,sweep) ! Frecuency Sweep

	use fedata
    use plot_routiner
	use solve_handle_complex

    integer, intent(in) :: flag!, points
    !    character(len=*), intent(in) :: title
	integer :: i
 !    real(8),pointer, intent(in) :: sweep(:)
    real(8), intent(in) :: rho(ne) , sweep(:)
    real(8), dimension(:,:), allocatable :: sweepout ! Matrix med frekvens (:,1) og eftergivenhed (:,2)
    real(8) :: lower, upper, omega
    integer :: points, flagger
    flagger = 0

    lower = sweep(1)!300!0.0E+4
    upper = sweep(2)!7.0E+4
    points = (upper-lower)/50d0

    allocate(sweepout(points,4))
    sweepout = 0.0d0

	print*, 'sweeping from omega ', lower, ' to ', upper

	do i=1,points-1
       sweepout(i,1) = lower+(real(i)-1.0)*(upper-lower)/(real(points)-1.0) !Line�r interpolation af frekvens
	end do
    sweepout(points,1) = upper

    do i=1,points

       omega = 2*pi*sweepout(i,1) ! vinkelfrekvens
       call buildS_piezo(flagger,omega)
       call enforce_piezo
       if (i == 1) then ! initialise mumps
          call mumps_init_complex
       end if
       call mumps_solve_complex_sweep
       !$$$$$$       sweepout(i,2) = ABS( SUM( dz * p) ) ! Dynamisk eftergivenhed I

       sweepout(i,2) = SQRT(DOT_PRODUCT(dz(1:3:neqn),dz(1:3:neqn)))! NOTICE: dot_product is SUM(CONJG(X)*Y) for complex vectors
       sweepout(i,3) = SQRT( DOT_PRODUCT(REAL(dz(1:3:neqn)),REAL(dz(1:3:neqn))) )
       sweepout(i,4) = SUM(DABS(REAL(dz(1:3:neqn))))/(neqn/3d0)
       print*, 'Frekvens ',sweepout(i,1) ,'Point ',i,'of ', points
       print*, 'Eftergivenhed', sweepout(i,2)

       pZ = CMPLX(0d0,0d0)
    end do

    ! destry mumps
    call mumps_finilize_complex

    call output_matrix(sweepout,'sweep_'//'ring')

  END SUBROUTINE sweep_piezo


  subroutine pade_approx_piezo




  end subroutine pade_approx_piezo



end module piezo