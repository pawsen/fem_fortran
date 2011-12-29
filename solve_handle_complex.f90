module solve_handle_complex

  implicit none
  
  private
  public :: mumps_init_complex, mumps_finilize_complex, mumps_solve_complex_sweep
  public :: mumps_solve_complex


    INCLUDE 'mpif.h'
    INCLUDE 'zmumps_struc.h'

    TYPE (ZMUMPS_STRUC):: id
    INTEGER:: IERR

contains

  subroutine mumps_solve_complex

    use fedata
    use plot_routiner


    TYPE (ZMUMPS_STRUC):: id
!    TYPE (DMUMPS_STRUC):: id
    INTEGER:: IERR, i , n ,nz,e, n2, display_print
    !      integer, allocatable, dimension(:) :: iK, jK
    real(8), allocatable, dimension(:) :: pp!sK, p


    CALL MPI_INIT(IERR)
    id%COMM = MPI_COMM_WORLD
    id%SYM = 0
    id%PAR = 1
    id%JOB = -1
    
    CALL ZMUMPS(id)
    IF ( id%MYID .eq. 0 ) THEN
       id%N = neqn_nb
       id%NZ = size(iK,1)

       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sKZ
       id%RHS => pZ

    END IF

    !set level for message output
    id%ICNTL(3) = 0 !
    id%ICNTL(4) = 1

    !Call package for solution    
    id%JOB = 6
    CALL ZMUMPS(id)

    !Solution has been assembled on the host
    IF ( id%MYID .eq. 0 ) THEN
       dZ = id%RHS
    
       ! Calculate forces
       if (antype == 'PIEZO' .and. &
            ((elem_type == 'PLANE42').or.(elem_type == 'PLANE_GMSH')) ) then
          n = (8*8+4*4)*ne+2*4*8*ne
          n2 = neqn+nn

       elseif (antype == 'PIEZO' .and. (elem_type == 'PLATE')) then
          n = (12*12+4*4)*ne+2*4*12*ne
          n2 = neqn+nn
          print*,'n2',n2
!!$          n = (12*12)*ne
!!$          n2 = neqn
       elseif (antype == 'PIEZO' .and. (elem_type == 'PLATE_GMSH') ) then
          n = (12*12+4*4)*ne_ptz+2*4*12*ne_ptz + 12*12*(ne-ne_ptz)
          n2 = neqn+nn
          print*,'n2',n2
!!$          n = (12*12)*ne
!!$          n2 = neqn
       else
          select case(antype)
          case ("PLATE")
             n = 12*12*ne
          case default
             n = 8*8*ne
          end select
          n2 = neqn
       end if   
    end if

    !call sparse_multiply(iK(1:n),jK(1:n),sK(1:n),id%RHS(1:n2),p(1:n2))

    print*,'nb',nb

    ! deallocate pointers
    NULLIFY (id%IRN)
    NULLIFY (id%JCN)
    NULLIFY (id%A)
    NULLIFY (id%RHS)
    !Destroy the instance (deallocate internal data structures)
    id%JOB = -2
    CALL ZMUMPS(id)
!    CALL DMUMPS(id)

    CALL MPI_FINALIZE(IERR)

  end subroutine mumps_solve_complex

  subroutine mumps_solve_complex_sweep

    use fedata

!!$    ! point to stiffness vectors
!!$    id%IRN => iK
!!$    id%JCN => jK
!!$    id%A => sKZ
!!$    id%RHS => pZ

    ! factorization and solving
    id%JOB = 5
    CALL ZMUMPS(id)
    
    !Solution has been assembled on the host
    IF ( id%MYID /= 0 ) THEN
       print*,'error, Solution has not been assembled in mumps -> solve_handle_complex.f90'
       stop
    end if
    dZ = id%RHS

!!$    ! deallocate pointers
!!$    NULLIFY (id%IRN)
!!$    NULLIFY (id%JCN)
!!$    NULLIFY (id%A)
!!$    NULLIFY (id%RHS)

    !status af pointer kan tjekkes med, se !http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/lref_for/source_files/rfpoin.htm
    !IF (.NOT.ASSOCIATED(arrow))

  end subroutine mumps_solve_complex_sweep


  subroutine mumps_init_complex

    use fedata

    CALL MPI_INIT(IERR)
    id%COMM = MPI_COMM_WORLD
    id%SYM = 0
    id%PAR = 1
    id%JOB = -1


    CALL ZMUMPS(id)
    IF ( id%MYID .eq. 0 ) THEN
       id%N = neqn_nb
       id%NZ = size(iK,1)
       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sKZ
       id%RHS => pZ
    END IF

    ! message 
    !id%ICNTL(1) = 0 ! error message
    id%ICNTL(3) = 0 !
    id%ICNTL(4) = 1

    ! perform analysis. Ignores numericel values, eg symbolic factorising
    id%JOB = 1
    CALL ZMUMPS(id)
   
  end subroutine mumps_init_complex
  
  subroutine mumps_finilize_complex

    IF (ASSOCIATED(id%IRN)) then
       ! deallocate pointers
       NULLIFY (id%IRN)
       NULLIFY (id%JCN)
       NULLIFY (id%A)
       NULLIFY (id%RHS)
    end if

    !Destroy the instance (deallocate internal data structures)
    id%JOB = -2
    CALL ZMUMPS(id)
    CALL MPI_FINALIZE(IERR)

  end subroutine mumps_finilize_complex


end module solve_handle_complex
