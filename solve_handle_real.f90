module solve_handle_real

  implicit none

  private
  public :: mumps_init_real, mumps_finalize_real
  public :: mumps_solve_real, mumps_solve_arpack_real

  INCLUDE 'mpif.h'
  INCLUDE 'dmumps_struc.h'

  TYPE (DMUMPS_STRUC):: id
  INTEGER:: IERR

contains

  subroutine mumps_init_real!(eigenvalue)

    use fedata
    
    !logical, intent(in) :: eigenvalue
    !integer, intent(IN) :: type
    integer :: i

    CALL MPI_INIT(IERR)
    !Define a communicator for the package
    id%COMM = MPI_COMM_WORLD
    !Ask for unsymmetric code
    id%SYM = 0
    !Host working
    id%PAR = 1
    !Initialize an instance of the package
    id%JOB = -1

    CALL DMUMPS(id)
    IF ( id%MYID .eq. 0 ) THEN

       id%N = neqn_nb
       id%NZ = size(iK,1)
       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sK

       if ((antype == 'EIGEN') .or. eigenvalue%calc) then
          ALLOCATE( id%RHS ( id%N ) )
       else
          ! id%RHS => p
          ! pointer virker ikke for p.
          ALLOCATE( id%RHS ( id%N ) )
       end if

    ELSE
       print*,'error, JOB not initialized-> solve_handle_real.f90'
       stop 
    END IF

    ! message 
    !id%ICNTL(1) = 0 ! error message
    id%ICNTL(3) = 0 !
    id%ICNTL(4) = 1

  end subroutine mumps_init_real


  subroutine mumps_finalize_real

    IF (ASSOCIATED(id%IRN)) then
       ! deallocate pointers
       NULLIFY (id%IRN)
       NULLIFY (id%JCN)
       NULLIFY (id%A)
       NULLIFY (id%RHS)
    end if

    !Destroy the instance (deallocate internal data structures)
    id%JOB = -2
    CALL DMUMPS(id)
    CALL MPI_FINALIZE(IERR)

  end subroutine mumps_finalize_real


  subroutine mumps_solve_real(type)

    use fedata
    use plot_routiner
    use numeth

    integer, intent(IN) :: type

    INTEGER::  i ! , n ,nz,e, n2, display_print

    !real(8), allocatable, dimension(:) ::
    
    if (.not. eigenvalue%calc) then
       id%RHS = d
    end if

    !Call package for solution
    select case(type)
    case(1)! perform analysis. Ignores numericel values, eg symbolic factorising
       id%JOB = 1
    case(2)! performs the factorization. Uses numerical values and information from JOB=1
       id%JOB = 2
    case(3)!performs the solution. eg JOB=1,2 must be run before this
       id%JOB = 3
    case(4)!combines the actions of JOB=1 with those of JOB=2.
       id%JOB = 4
    case(5)! combines the actions of JOB=2 and JOB=3.
       id%JOB = 5
    case(6)!combines the actions of calls with JOB=1, 2, and 3., eg only JOB=-1 must be run
       id%JOB = 6
    end select

    CALL DMUMPS(id)
    !Solution has been assembled on the host
    IF ( id%MYID /= 0 ) THEN
       print*,'error, Solution has not been assembled in mumps -> solve_handle_real.f90'
       print*,'JOB ID: ',id%JOB
       stop
    END IF

    select case(type)
    case(3,5:6)
       d = id%RHS
    end select

  end subroutine mumps_solve_real


  subroutine mumps_solve_arpack_real(rhs,type)

    use fedata
    use plot_routiner
    use numeth

    integer, intent(IN) :: type
    real(8), intent(inout) :: rhs(:)
    integer :: i, n, nz

    n = size(rhs,1)
    nz = id%N

    id%RHS(1:n) = rhs
    id%RHS(n+1:nz) = 0d0

    !Call package for solution
    select case(type)
    case(3)!performs the solution. eg JOB=1,2 must be run before this
       id%JOB = 3
    end select

    CALL DMUMPS(id)
    !Solution has been assembled on the host
    IF ( id%MYID /= 0 ) THEN
       print*,'error, Solution has not been assembled in mumps -> solve_handle_arpack_real.f90'
       print*,'JOB ID: ',id%JOB
       stop
    END IF

    rhs = id%RHS(1:n)

  end subroutine mumps_solve_arpack_real


  
!!$! Calculate forces
!!$if (antype == 'PIEZO' .and. &
!!$     ((elem_type == 'PLANE42').or.(elem_type == 'PLANE_GMSH')) ) then
!!$   n = (8*8+4*4)*ne+2*4*8*ne
!!$   n2 = neqn+nn
!!$
!!$elseif (antype == 'PIEZO' .and. (elem_type == 'PLATE')) then
!!$   n = (12*12+4*4)*ne+2*4*12*ne
!!$   n2 = neqn+nn
!!$   print*,'n2',n2
!!$  n = (12*12)*ne
!!$  n2 = neqn
!!$elseif (antype == 'PIEZO' .and. (elem_type == 'PLATE_GMSH') ) then
!!$   n = (12*12+4*4)*ne_ptz+2*4*12*ne_ptz + 12*12*(ne-ne_ptz)
!!$   n2 = neqn+nn
!!$   print*,'n2',n2
!!$  n = (12*12)*ne
!!$  n2 = neqn
!!$else
!!$   select case(antype)
!!$   case ("PLATE")
!!$      n = 12*12*ne
!!$   case default
!!$      n = 8*8*ne
!!$   end select
!!$   n2 = neqn
!!$end if
!!$call sparse_multiply(iK(1:n),jK(1:n),sK(1:n),id%RHS(1:n2),p(1:n2))
!!$
!!$display_print = 0
!!$if( display_print == 1) then
!!$
!!$   if ((antype == 'PIEZO') .and. (elem_type ==  'PLANE42')) then
!!$
!!$      write (*, '(" Node            Displacements")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9)') i, id%RHS(2*i-1), id%RHS(2*i)
!!$      end do
!!$
!!$      write (*, '(" Node             Applied/reaction forces")')
!!$      write (*, '("number     1-direction     2-direction")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9)') i, p(2*i-1), p(2*i)
!!$      end do
!!$
!!$      write (*, '(" Node            Potential")')
!!$      do i =1,nn !neqn+1,neqn+nn
!!$         write (*, '(i6, 1x, f25.15)') i, id%RHS(neqn+i)
!!$      end do
!!$
!!$      write (*, '(" dof            lambda")')
!!$      do i = neqn+nn+1, neqn+nn + nb
!!$         write (*, '(i6, 1x, f15.9)') i, id%RHS(i)
!!$      end do
!!$
!!$   else if ((antype == 'PIEZO') .and. (elem_type ==  'PLATE')) then
!!$
!!$      write (*, '(" Node      Displacements              Slope            ")')
!!$      write (*, '("number     w-direction     x-direction      y-direction")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i,id%RHS(3*i-2), id%RHS(3*i-1),id%RHS(3*i)
!!$      end do
!!$
!!$      write (*, '(" Node            Potential")')
!!$      do i =1,nn !neqn+1,neqn+nn
!!$         write (*, '(i6, 1x, f25.15)') i, id%RHS(neqn+i)
!!$      end do
!!$
!!$      write (*, '(" Node                Applied/reaction forces")')
!!$      write (*, '("number     1-direction     2-direction     3-direction")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, p(3*i-2), p(3*i-1), p(3*i)
!!$      end do
!!$
!!$      write (*, '(" dof            lambda")')
!!$      do i = neqn+nn+1, neqn+nn + nb
!!$         write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
!!$      end do
!!$
!!$
!!$   else if (elem_type ==  'PLATE') then
!!$
!!$      write (*, '(" Node      Displacements              Slope            ")')
!!$      write (*, '("number     w-direction     x-direction      y-direction")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i,id%RHS(3*i-2), id%RHS(3*i-1),id%RHS(3*i)
!!$      end do
!!$
!!$      write (*, '(" Node                Applied/reaction forces")')
!!$      write (*, '("number     1-direction     2-direction     3-direction")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, p(3*i-2), p(3*i-1), p(3*i)
!!$      end do
!!$
!!$      write (*, '(" dof            lambda")')
!!$      do i = neqn+1, neqn + nb
!!$         write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
!!$      end do
!!$
!!$   else
!!$      write (*, '(" Node            Displacements")')
!!$      do i = 1, nn
!!$         write (*, '(i6, 1x, f25.15, 1x, f25.15)') i, id%RHS(2*i-1), id%RHS(2*i)
!!$      end do
!!$
!!$      write (*, '(" dof            lambda")')
!!$      do i = neqn+1, neqn + nb
!!$         write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
!!$      end do
!!$   end if

!!$    write (*, '("#####FORCE#####")')
!!$    write (*, '(" Node            P")')
!!$
!!$    do i = 1, nn
!!$       write (*, '(i6, 1x, f25.15, 1x, f25.15)') i, pp(2*i-1), pp(2*i)
!!$    end do
!!$
!!$    if (antype == 'PIEZO') then
!!$       write (*, '(" dof            potentiale")')
!!$       do i = 1,nn!neqn+1, neqn + nn
!!$          write (*, '(i6, 1x, f25.15)') i, pp(neqn+i)
!!$       end do
!!$    else
!!$       write (*, '(" dof            lambda")')
!!$       do i = neqn+1, neqn + nb
!!$          write (*, '(i6, 1x, f25.15)') i, pp(i)
!!$       end do
!!$    end if

end module solve_handle_real
