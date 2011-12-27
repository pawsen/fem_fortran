MODULE laminate
  ! This module contains subroutines specific to the MINDLIN41 element.
  IMPLICIT NONE

  private :: initial_laminate
  public :: displ_laminate


contains
  
  subroutine initial_laminate

    ! This subroutine is mainly used to allocate vectors and matrices
    use fedata
    use link1
    use plane42
    use mindlin42

    integer :: e, nen, bw_e
    
    integer, parameter :: mdim = 12

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
    allocate (p(neqn), d(neqn))

    allocate (strain(ne, 5), stress(ne, 5))

    ! MODIFIED FOR WINDOWS COMPILER
    strain = 0d0
    stress = 0d0

  end subroutine initial_laminate

  subroutine displ_laminate

    ! This subroutine calculates displacements

    use numeth
    use processor
    use fedata
    use plot_routiner

    integer :: e
    real(8), dimension(:), allocatable :: plotval

    call initial_laminate

    ! Build load-vector
!    call buildload_plate

    ! Build stiffness matrix
    !call buildstiff(plate_type)

    ! Remove rigid body modes
    !call enforce

    if ( banded == 0) then
       call factor(k)
       call solve(k, p)
    elseif(banded == 1) then
       call bfactor(k)
       call bsolve(k, p)
    end if
    ! Transfer results
    d(1:neqn) = p(1:neqn)
    ! Recover stress
    !call recover(plate_type)
    ! Output results
    call output
    call plot('un/deformed', 'matlab', 'color', ' ')
    allocate(plotval(ne))
    plotval = 1d0
    !call output_deformed_plate('mindlin',plotval )


  end subroutine displ_laminate
END MODULE laminate


