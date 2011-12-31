module plate

  implicit none

  private :: enforce, recover, rb_init_plate, load_init_plate
  private :: size_of_stucture
  public :: displ_plate, initial_plate, buildload_plate, buildstiff_plate

contains

  subroutine initial_plate!(plate_type)

    ! This subroutine is mainly used to allocate vectors and matrices
    use fedata
    use link1
    use plane42
    use mindlin42
    use file_init
    
    !integer, intent(in) :: plate_type
    integer :: e, nen, bw_e, nzz
    integer:: rand_int(4), int_dummy(11)
    real(8) ::  rand_val(4), real_dummy(7)

    integer, parameter :: mdim = 12
    integer:: plate_type = 2

    ! add rb along one of the sides

    call parameter_input(real_dummy,int_dummy,rand_int,rand_val)
    call rb_init_plate(rand_int,rand_val)
    !call load_init_plate ! inspænder plade langs randen og giver jævnt fordelt last på fladen
    
    
    ! This subroutine computes the number of global equation,
    ! half bandwidth, etc and allocates global arrays.
    
    ! Calculate number of equations
    neqn = 3*nn

    if ( banded == 0) then
       allocate (k(neqn, neqn))
    elseif (banded == 1) then
       bw = 0
       do e=1,ne
          nen = element(e)%numnode
          bw_e = (maxval(element(e)%ix(1:nen))-(minval(element(e)%ix(1:nen)))+1)*3
          if (bw_e > bw) then
             bw = bw_e
          end if
       end do
       allocate (k(bw, neqn))
       print*, 'Bandwidth =', bw

    end if

    if(banded == 2) then ! sparse
       !12*12*ne is the number of times the local stiffness matrix adds a value to the global(including dublicates, eg. is't not nz(number of non-zeros)). Because of lagrangien multiplier there is added xtra two non-zero entries for each RB

       if ( (eigenvalue%calc .and.  (eigenvalue%shift .eqv. .false.)) ) then
          nzz = 12*12*ne
          neqn_nb = neqn
       else
          nzz = 12*12*ne+2*nb
          neqn_nb = neqn+nb
       end if

       allocate(iK(nzz),jK(nzz),sK(nzz))

       if (.not. eigenvalue%calc) then
          allocate (p(neqn+nb), d(neqn+nb))
       end if
    else
       allocate (p(neqn), d(neqn))
    end if

    if (.not. eigenvalue%calc) then
       select case(plate_type)
       case (1)!kirchhoff
          allocate (strain(ne, 3), stress(ne, 3))
       case (2) !mindlin
          allocate (strain(ne, 5), stress(ne, 5))
       end select

       ! MODIFIED FOR WINDOWS COMPILER
       strain = 0d0
       stress = 0d0
    end if

    print*,'total DOF',neqn
  end subroutine initial_plate

  subroutine displ_plate

    ! This subroutine calculates displacements

    use numeth
    use processor
    use fedata
    use plot_routiner
    use solve_handle_real

    integer :: e, plate_type !1 = kirchhoff, 2 = mindlin
    real(8), dimension(:), allocatable :: plotval, p1,p2,angle,vonMises

    !1 = kirchhoff, 2 = mindlin
    plate_type = 2

    call initial_plate

    ! Build load-vector
    call buildload_plate

    ! Build stiffness matrix
    call buildstiff_plate

    ! Remove rigid body modes
    call enforce
    
    
    if (banded /= 2) then
       d = p
       if ( banded == 0) then
          call factor(k)
          call solve(k, d)
       elseif(banded == 1) then
          call bfactor(k)
          call bsolve(k, d)
       end if
       ! Recover stress
       !call recover(plate_type)

    else
       call mumps_init_real
       call mumps_solve_real(6)
       call mumps_finalize_real
    end if

    allocate(plotval(ne),p1(ne),p2(ne),angle(ne),vonMises(ne))
    call recover(plate_type,p1,p2,angle,vonMises)

    ! Output results
    call output(plate_type)
    plotval = 1d0
    call output_deformed('mindlin','vonMises', vonMises)  
    call plot_vector('principal stresses',p1,p2,angle)


  end subroutine displ_plate
  
  subroutine buildload_plate

    ! This subroutine builds the global load vector.

    use fedata
    use plate41rect
    use mindlin42
    use numeth, only : get_edof
    
    integer :: i, j, e, nen, pdof, eface
    integer, parameter :: mdim = 12
    integer, dimension(mdim) :: edof
    real(8), dimension(8) :: xe
    real(8), dimension(mdim) :: re
    real(8) :: fe, thk


    p = 0.d0
    do i = 1, np

       ! Build nodal load contribution
       if (loads(i, 1) == 1) then 
          pdof = loads(i, 2)*3-3+loads(i, 3)
          p(pdof) = loads(i, 4)

          ! Build uniformly distributed surface (pressure) load contribution
       else if (loads(i, 1) == 2) then
          ! Find coordinates and degrees of freedom
          e   = loads(i,2)  
          nen = element(e)%numnode
          fe    = loads(i,4)
          eface = loads(i,3)
          thk   = mprop(element(e)%mat)%thk

          call get_edof(e,nen,xe,edof=edof)

          ! NB der er forskel på hvordan lastvektoren laves, alt efter om det er mindlin eller kirchoff. Derfor virker nedenstående kun for mindlin!
          call mindlin42_re(xe,eface,fe,thk,re)

          p(edof) = p(edof) + re    
       else
          print *, 'Error in plate/buildload' 
          print *, 'Load type not known'
          error stop
       end if
    end do

  end subroutine buildload_plate


  subroutine buildstiff_plate(rho,rho_min)

    ! This subroutine builds the global stiffness matrix from
    ! the local element stiffness matrices.

    use fedata
    use plate41rect
    use mindlin42
    use fea, only : assemble_sparse
    use numeth, only : get_edof

    !integer, intent(IN):: plate_type
    real(8), optional, INTENT(IN) :: rho_min, rho(:)
    integer :: e, i, j, ii, jj
    integer :: nen, idof
    integer, parameter :: mdim = 12, kk = 3
    integer, dimension(mdim) :: edof
    ! Remember that dimension of xe is only dependent on problem
    !  type(2d/3d), and not DOF pr node.
    real(8), dimension(8) :: xe
    real(8), dimension(mdim, mdim) :: ke, me
    real(8) :: young, nu, dens, thk, shear
    integer, parameter :: plate_type = 2

    ! Reset stiffness matrix
    if ( banded == 0) then
       k(1:neqn, 1:neqn) = 0.
    elseif(banded == 1) then
       k(1:bw, 1:neqn) = 0.
    end if

    ii = 0! sparse
    jj = 0
    
    do e = 1, ne

       call get_edof(e,nen,xe,edof=edof)

       ! calculate plate41 element stiffness matrix here
       young = mprop(element(e)%mat)%young
       nu    = mprop(element(e)%mat)%nu
       thk   = mprop(element(e)%mat)%thk
       dens  = mprop(element(e)%mat)%dens

       select case(plate_type)
       case (1)
          call plate41rect_ke(xe, young, nu, thk, ke)
       case (2)
          shear = mprop(element(e)%mat)%shear
          call mindlin42_ke(xe, young, nu, dens, thk, shear,ng, ke,1,mat_vec) !type = 1
       end select

       ! Building global stiffness matrix
       if ( banded == 0) then
          k(edof, edof) = k(edof, edof) + ke
       elseif (banded == 1) then
          do i = 1, 3*nen
             do j = 1, 3*nen
                if (edof(i) .ge. edof(j)) then
                   k (1+edof(i)-edof(j),edof(j)) = k(1+edof(i)&
                        &-edof(j),edof(j)) + ke(i,j)
                end if
             end do
          end do

       elseif (banded == 2) then
          if (present(rho) .and. ((element(e)%mat /= elem_id))) then
             call assemble_sparse(nen,kk,ii,edof,ke,jj,rho,rho_min)
          else
             call assemble_sparse(nen,kk,ii,edof,ke)
          end if
       end if

    end do

    !add values from lagrangian multipliers
    if ( banded == 2 .and. ((antype /= 'EIGEN') .or. eigenvalue%shift)) then
       do i=1,nb
          ! UPS: bound(i,3) is used both for temp and moments for
          !  plates. Stupid!
          idof = 3*(bound(i,1)-1) + bound(i,2)
          ! nederste matrix
          ii = ii+1
          iK(ii) = neqn+i
          jK(ii) = idof
          sK(ii ) = 1d0
          ! matrix to the right
          ii = ii+1
          iK(ii) = idof
          jK(ii) = neqn+i
          sK(ii ) = 1d0
       end do

    end if

  end subroutine buildstiff_plate

  subroutine enforce

    ! This subroutine enforces the support boundary conditions.
    use fedata

    integer :: i,s, idof
    real(8) :: penal_fak
    real(8), dimension(neqn) :: k_vector

    ! Correct for supports
    if ( banded == 0) then
       if (.not. penalty) then
          do i = 1, nb
             idof = 3*(bound(i,1)-1) + bound(i,2)
             p(1:neqn) = p(1:neqn) - k(1:neqn, idof) * bound(i, 3)
             p(idof) = bound(i, 3)
             k(1:neqn, idof) = 0.
             k(idof, 1:neqn) = 0.
             k(idof, idof) = 1.
          end do
       else
          penal_fak = 10e9*maxval(K)
          do i = 1, nb
             idof = 2*(bound(i,1)-1) + bound(i,2)
             k(idof, idof) = k(idof, idof) + penal
             p(idof) = penal * bound(i, 3)  
          end do
       endif
    elseif (banded == 1) then
       do i = 1, nb
          idof = 3*(bound(i,1)-1) + bound(i,2)
          k_vector(1:neqn) = 0 
          if (idof<bw) then
             k_vector(idof:idof+bw-1) = k(1:bw, idof)
             do s = 1,(idof-1)
                k_vector(s) = k(idof-(s-1),s)
             end do
          elseif (idof<neqn-bw) then
             k_vector(idof:idof+bw-1) = k(1:bw, idof)
             do s = 1,bw
                k_vector(s) = k(bw-(s-1),idof-(bw-(s-1)))
             end do
          else
             k_vector(idof:neqn) = k(1:(neqn-idof+1),idof)
             do s = 1,bw
                k_vector(s) = k(bw-(s-1),idof-(bw-(s-1)))
             end do
          end if
          p(1:neqn) = p(1:neqn) - k_vector * bound(i, 3)
          p(idof) = bound(i, 3)
          k(1:bw, idof) = 0.
          if (idof<bw) then
             do s = 1,idof
                k(idof-(s-1),s) = 0d0
             end do
          else
             do s = 1,bw
                k(bw-(s-1),idof-(bw-(s))) = 0.
             end do
          end if
          k(1, idof) = 1.
       end do

    elseif(banded == 2) then !sparse. Add constraint to the load
       ! vector
       do i = 1,nb
          p(neqn+i) = bound(i, 3)
       end do
    endif


  end subroutine enforce


  subroutine recover(plate_type,p1,p2,angle,vonMises)

    ! This subroutine recovers the element stress, element strain, 
    ! and nodal reaction forces

    use fedata
    use plate41rect
    use mindlin42

    integer, intent(IN):: plate_type
    real(8),dimension(:), intent(out):: p1,p2,angle,vonMises

    integer :: e, i, nen
    integer, parameter :: mdim = 12
    integer :: edof(mdim)
    real(8) :: xe(8), de(mdim)

    real(8), dimension(mdim, mdim) :: ke, me
    real(8) :: young, nu, area, dens, thk, shear
    real(8) :: c_angle, s_angle
    real(8), dimension(:), allocatable :: estrain, estress
    select case(plate_type)
    case (1)
       allocate(estrain(3),estress(3))
    case(2)
       allocate(estrain(5),estress(5))
    end select

    p = 0d0

    do e = 1, ne

       ! Find coordinates etc...
       nen = element(e)%numnode
       do i = 1,nen
          xe(2*i-1) = x(element(e)%ix(i),1)
          xe(2*i) = x(element(e)%ix(i),2)

          edof(3*i-2) = 3 * element(e)%ix(i) - 2 
          edof(3*i-1) = 3 * element(e)%ix(i) - 1
          edof(3*i)   = 3 * element(e)%ix(i)

          de(3*i-2) = d(edof(3*i-2))
          de(3*i-1) = d(edof(3*i-1))
          de(3*i)   = d(edof(3*i))
       end do
       
       ! Find stress and strain
       select case( element(e)%id )
       case( 1 )
          print*,'This is only for shells'
       case( 2 )
          print*,'Error in recovery_plate, forgot to change element_id to shell'
       case( 3 )
          young = mprop(element(e)%mat)%young
          nu    = mprop(element(e)%mat)%nu
          thk   = mprop(element(e)%mat)%thk
          dens  = mprop(element(e)%mat)%dens
          select case(plate_type)
          case (1)
             call plate41rect_ke(xe, young, nu, thk, ke)
             call plate41rect_ss(xe, de, young, nu, thk, estress, estrain)
          case (2)
             shear = mprop(element(e)%mat)%shear
             call mindlin42_ke(xe,young, nu, dens, thk, shear,ng, ke,1,mat_vec)
             call MINDLIN42_ss(xe, de, young, nu, thk, shear, estress, estrain)
          end select
       end select

          p(edof) = p(edof) + matmul(ke, de)
          stress(e, :) = estress
          strain(e, :) = estrain

          vonMises(e) = (estress(1)**2+estress(2)**2-estress(1)*estress(2)+3.0*estress(3)**2)**(0.50d0) !von mises spænding
          p1(e) = (estress(1)+estress(2))*0.50d0+dsqrt(((estress(1)-estress(2))*0.50d0)**2+estress(3)**2) !hovedspænding 1
          p2(e) = (estress(1)+estress(2))*0.50d0-dsqrt(((estress(1)-estress(2))*0.50d0)**2+estress(3)**2) !hovedspænding 2
          c_angle = (estress(1)-estress(2))/(p1(e)-p2(e)) ! cosins til vinlen for hovedspænding 1
          s_angle = (-2.0d0*estress(3))/(p1(e)-p2(e))! sinus til vinlen for hovedspænding 2
          angle(e) = datan2(s_angle,c_angle)/2.0d0! vinklen for hovedspænding 1. positiv i retning mod uret. [rad]

       end do

     end subroutine recover

  subroutine rb_init_plate(rand_int,rand_val)
    ! This subroutine adds a distributed load along an edge given by rand_int.

    use fedata
    use plot_routiner

    integer, intent(IN) :: rand_int(:) ! må ikke hedde rand, da rand er initialiseret i fedata, DOH!!!
    real(8), intent(IN) :: rand_val(:)
    integer :: i, e, nen, np_temp
    integer :: ny,nx
    integer(8), dimension(:,:),allocatable :: loads_temp
    integer, dimension(:) ,pointer :: int_parameters

    call size_of_stucture(int_parameters)
    nx = int_parameters(1)
    ny = int_parameters(2)


    allocate (loads_temp(np, 5))
    loads_temp = loads
    DEALLOCATE(loads)

    np_temp = np
    do e=1,size(rand_int)! find number of bounds including the new force
       select case(rand_int(e))
       case(1) ! bottom
          np_temp = np_temp + nx
       case(2)
          np_temp = np_temp + ny
       case(3)
          np_temp = np_temp + nx
       case(4)
          np_temp = np_temp + ny
       end select
    end do
    allocate (loads(np_temp,5))
    loads(1:np,:) = loads_temp

    ! NB hjørne-knuder vil indgå to gange hvis, hvis to tilstødende rande påskrives.
    ! DET ER ET PROBLEM!!!
    do e=1,size(rand_int) ! set RB for force
       select case(rand_int(e))
       case(1) ! bottom
          do i=1,nx
             np = np+1

             loads(np,1) = 1
             loads(np,2) = 1+(i-1)*ny 			! node number
             loads(np,3) = 1 					! direction of load, fz.
             if (i == 1 .or. i == nx) then
                loads(np,4) = 0.5*rand_val(e)	! Value of rb, eg Fx = 10 N ..
             else
                loads(nb,4) = rand_val(e)
             end if
          end do
       case(2) ! right side of structure
          do i=1,ny
             np = np+1

             loads(np,1) = 1
             loads(np,2) = nn-ny+i 
             loads(np,3) = 1 
             if (i == 1 .or. i == ny) then
                loads(np,4) = 0.5*rand_val(e)
             else
                loads(np,4) = rand_val(e)
             end if

          end do
       case(3) ! top
          do i=1,nx
             np = np+1

             loads(np,1) = 1
             loads(np,2) = i*ny 
             loads(np,3) = 1 
             if (i == 1 .or. i == nx) then
                loads(np,4) = 0.5*rand_val(e)
             else
                loads(np,4) = rand_val(e)
             end if
          end do
       case(4) ! left side of structure
          do i=1,ny
             np = np+1

             loads(np,1) = 1
             loads(np,2) = i 
             loads(np,3) = 1 
             if (i == 1 .or. i == ny) then
                loads(np,4) = 0.5*rand_val(e)
             else
                loads(np,4) = rand_val(e)
             end if
          end do
       end select
    end do


  end subroutine rb_init_plate

  subroutine load_init_plate
    ! tilføjer fordelt latteral last på alle knuder
    ! og simpelt understøttet langs randen

    use fedata
    use numeth

    integer :: e, nx, ny, rand_int(4),i
    integer, dimension(:) ,pointer :: int_parameters
    integer, allocatable :: index_nb(:)
    real(8), allocatable :: bound_buf(:,:)

    ! FIRST: add force
    np = ne
    DEALLOCATE(loads)
    allocate(loads(np,5))

    do e=1,ne
       loads(e, 1) = 2
       loads(e, 2) = e          ! Element number
       loads(e, 3) = 5          ! Surface number
       loads(e, 4) = -0.00001d0         ! Size of load (pressure)
    end do 

    ! SECOND: add boundary
    call size_of_stucture(int_parameters)
    nx = int_parameters(1)
    ny = int_parameters(2)


    allocate(bound_buf(2*(nx+ny),3))
    rand_int = (/1,2,3,4/)
    nb = 0
   do e=1,size(rand_int) ! set RB for force
       select case(rand_int(e))
       case(1) ! bottom
          do i=1,nx
             nb = nb+1
             bound_buf(nb,1) = 1+(i-1)*ny 	! node number
             bound_buf(nb,2) = 1 			! direction
             bound_buf(nb,3) = 0			! prescribed
          end do
       case(2) ! right side of structure
          do i=1,ny
             nb = nb+1
             bound_buf(nb,1) = nn-ny+i 
             bound_buf(nb,2) = 1 			
             bound_buf(nb,3) = 0			
          end do
       case(3) ! top
          do i=1,nx
             nb = nb+1
             bound_buf(nb,1) = i*ny
             bound_buf(nb,2) = 1 			
             bound_buf(nb,3) = 0			

          end do
       case(4) ! left side of structure
          do i=1,ny
             nb = nb+1
             bound_buf(nb,1) = i
             bound_buf(nb,2) = 1 			
             bound_buf(nb,3) = 0
          end do
       end select
    end do

    !NB hjørne-knuder indgår to gange, hvor to tilstødende rande mødes.
    ! De fjernes
    allocate(index_nb(nb))
    call index_dups(INT(bound_buf(:,1:2)),index_nb,i)
    DEALLOCATE(bound)
    allocate(bound(i,3))
    bound = bound_buf( index_nb(1:i),: )
    nb = i

  end subroutine load_init_plate

  subroutine size_of_stucture(int_parameters)

    use fedata

    integer, pointer, intent(out) :: int_parameters(:)
!    real(8), intent(IN) ::
    integer :: i, e, nen
    integer :: number_elements_x, number_elements_y, number_nodes_y, number_nodes_x, ny,nx
    real(8) :: xmin, xmax, ymin, ymax, lx, ly, length_first_element, hight_first_element

    ! tjekker at int_parameter ikke er allokeret
    if(associated(int_parameters))deallocate(int_parameters) 
    allocate(int_parameters(4)) 

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
    number_elements_x = NINT(lx/length_first_element) ! Antal elemnter i x-retningen ! runder op til nærmeste integer. Hvis nint ikke medtages, rundes der forkert af
    number_nodes_x = number_elements_x+1 ! Alternativt nn/number_nodes_y
    number_elements_y = NINT(ly/hight_first_element) ! Antal elemnter i y-retningen
    number_nodes_y = number_elements_y+1 !Antal knuder i y-retningen

    ny = number_nodes_y
    nx = number_nodes_x

    int_parameters(1) = nx
    int_parameters(2) = ny
    int_parameters(3) = number_elements_x
    int_parameters(4) = number_elements_y


  end subroutine size_of_stucture


end module plate
