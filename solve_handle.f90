module solve_handle

  implicit none
  
  private
  public :: mumps_solve_real, mumps_solve_complex
  public :: sparse_multiply


contains

  subroutine mumps_solve_real

    use fedata
    use plot_routiner

    INCLUDE 'mpif.h'
    INCLUDE 'dmumps_struc.h'

    TYPE (DMUMPS_STRUC):: id
    INTEGER:: IERR, i , n ,nz,e, n2, display_print
    !      integer, allocatable, dimension(:) :: iK, jK
    real(8), allocatable, dimension(:) :: pp!sK, p

!    n = neqn + nb
    n = size(p,1)
    nz = size(iK,1)
    print*,'n: ',n


!!$    do e = 1,nz
!!$       print*,e,iK(e),jK(e),sK(e)
!!$    end do


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
    !Define problem on the host (processor 0)
    ! id%ICNTL(14) = 20

    IF ( id%MYID .eq. 0 ) THEN
       id%N = n
       id%NZ = nz
       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sK
       id%RHS => p
    END IF
    !Call package for solution
    id%JOB = 6
    CALL DMUMPS(id)
    !Solution has been assembled on the host
    IF ( id%MYID .eq. 0 ) THEN

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
       call sparse_multiply(iK(1:n),jK(1:n),sK(1:n),id%RHS(1:n2),p(1:n2))
       
       display_print = 0
       if( display_print == 1) then

          if ((antype == 'PIEZO') .and. (elem_type ==  'PLANE42')) then

             write (*, '(" Node            Displacements")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9)') i, id%RHS(2*i-1), id%RHS(2*i)
             end do

             write (*, '(" Node             Applied/reaction forces")')
             write (*, '("number     1-direction     2-direction")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9)') i, p(2*i-1), p(2*i)
             end do

             write (*, '(" Node            Potential")')
             do i =1,nn !neqn+1,neqn+nn
                write (*, '(i6, 1x, f25.15)') i, id%RHS(neqn+i)
             end do

             write (*, '(" dof            lambda")')
             do i = neqn+nn+1, neqn+nn + nb
                write (*, '(i6, 1x, f15.9)') i, id%RHS(i)
             end do

          else if ((antype == 'PIEZO') .and. (elem_type ==  'PLATE')) then

             write (*, '(" Node      Displacements              Slope            ")')
             write (*, '("number     w-direction     x-direction      y-direction")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i,id%RHS(3*i-2), id%RHS(3*i-1),id%RHS(3*i)
             end do

             write (*, '(" Node            Potential")')
             do i =1,nn !neqn+1,neqn+nn
                write (*, '(i6, 1x, f25.15)') i, id%RHS(neqn+i)
             end do

             write (*, '(" Node                Applied/reaction forces")')
             write (*, '("number     1-direction     2-direction     3-direction")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, p(3*i-2), p(3*i-1), p(3*i)
             end do

             write (*, '(" dof            lambda")')
             do i = neqn+nn+1, neqn+nn + nb
                write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
             end do


          else if (elem_type ==  'PLATE') then

             write (*, '(" Node      Displacements              Slope            ")')
             write (*, '("number     w-direction     x-direction      y-direction")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i,id%RHS(3*i-2), id%RHS(3*i-1),id%RHS(3*i)
             end do

             write (*, '(" Node                Applied/reaction forces")')
             write (*, '("number     1-direction     2-direction     3-direction")')
             do i = 1, nn
                write (*, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, p(3*i-2), p(3*i-1), p(3*i)
             end do

             write (*, '(" dof            lambda")')
             do i = neqn+1, neqn + nb
                write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
             end do

          else
             write (*, '(" Node            Displacements")')
             do i = 1, nn
                write (*, '(i6, 1x, f25.15, 1x, f25.15)') i, id%RHS(2*i-1), id%RHS(2*i)
             end do

             write (*, '(" dof            lambda")')
             do i = neqn+1, neqn + nb
                write (*, '(i6, 1x, f25.15)') i, id%RHS(i)
             end do
          end if
       end if
    END IF

    d = id%RHS

    !call output_forskydninger(id%RHS,'forskydning')

    print*,'nb',nb


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

    ! deallocate pointers
    NULLIFY (id%IRN)
    NULLIFY (id%JCN)
    NULLIFY (id%A)
    NULLIFY (id%RHS)


    !Destroy the instance (deallocate internal data structures)
    id%JOB = -2
    CALL DMUMPS(id)
    CALL MPI_FINALIZE(IERR)

  end subroutine mumps_solve_real

  subroutine mumps_solve_complex

    use fedata
    use plot_routiner

    INCLUDE 'mpif.h'
    INCLUDE 'zmumps_struc.h'
!    INCLUDE 'dmumps_struc.h'

    TYPE (ZMUMPS_STRUC):: id
!    TYPE (DMUMPS_STRUC):: id
    INTEGER:: IERR, i , n ,nz,e, n2, display_print
    !      integer, allocatable, dimension(:) :: iK, jK
    real(8), allocatable, dimension(:) :: pp!sK, p

    n = size(pZ,1)
    nz = size(iK,1)
    print*,'n: ',n


    CALL MPI_INIT(IERR)
    id%COMM = MPI_COMM_WORLD
    id%SYM = 0
    id%PAR = 1
    id%JOB = -1
    
    CALL ZMUMPS(id)
    IF ( id%MYID .eq. 0 ) THEN
       id%N = n
       id%NZ = nz

       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sKZ
       id%RHS => pZ

    END IF
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
    use plot_routiner

    INCLUDE 'mpif.h'
    INCLUDE 'zmumps_struc.h'

    TYPE (ZMUMPS_STRUC):: id
    INTEGER:: IERR, i , n ,nz,e, n2, display_print


    n = size(pZ,1)
    nz = size(iK,1)

    CALL MPI_INIT(IERR)
    id%COMM = MPI_COMM_WORLD
    id%SYM = 0
    id%PAR = 1
    id%JOB = -1

    CALL ZMUMPS(id)
    IF ( id%MYID .eq. 0 ) THEN
       id%N = n
       id%NZ = nz
       ! point to stiffness vectors
       id%IRN => iK
       id%JCN => jK
       id%A => sKZ
       id%RHS => pZ
    END IF

    ! perform analysis. Ignores numericel values.
    id%JOB = 1
    CALL ZMUMPS(id)
   

    ! factorization and solving
    id%JOB = 5
    CALL ZMUMPS(id)
    
    !Solution has been assembled on the host
    IF ( id%MYID .eq. 0 ) THEN

    end if
    dZ = id%RHS
    print*,'nb',nb

    ! deallocate pointers
    NULLIFY (id%IRN)
    NULLIFY (id%JCN)
    NULLIFY (id%A)
    NULLIFY (id%RHS)

    !status af pointer kan tjekkes med, se !http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/lref_for/source_files/rfpoin.htm
    !IF (.NOT.ASSOCIATED(arrow))


    !Destroy the instance (deallocate internal data structures)
    id%JOB = -2
    CALL ZMUMPS(id)

    CALL MPI_FINALIZE(IERR)

  end subroutine mumps_solve_complex_sweep


  subroutine sparse_multiply(row,col,val,x,b)
    !multiply a sparse matrix, A, with a dense vector,x,, eg
    ! A*x = b
    
    real(8), dimension(:), intent(IN) :: val,x
    integer, dimension(:), intent(IN) :: row,col
    real(8), dimension(:), intent(OUT) :: b
    integer :: i,ii

    b = 0d0
    do ii=1,size(row,1)
       i = row(ii)
       b(i) = b(i) + val(ii)*x(col(ii))
    end do


!!$    !! test structure
!!$    real(8) :: val(9),xx(4),b(4)
!!$    integer :: row(9),col(9)
!!$    
!!$    row = (\1,1,2,3,3,3,4,4,1\)
!!$    col = (\1,4,2,1,1,3,2,4,1\)
!!$    val = (\0.5,3.0,2.,2.,2.,3.,6.,5.,0.5\)
!!$    xx = (\2.,4.,6.,8.\)
!!$
!!$    Call sparse_multiply(row,col,val,xx,b)
!!$    
!!$    write(*,'("b = ")')
!!$    do i=1,4
!!$       write (*, '(i6, 1x, f25.15)') i, b(i)
!!$    end do
!!$
!!$    !!A.x=b
!!$    !!   | 1     0     0     3| |2| |26|
!!$    !!   | 0     2     0     0| |4| |8 |
!!$    !!A: | 4     0     3     0|*|6|=|26|
!!$    !!   | 0     6     0     5| |8| |64|



  end subroutine sparse_multiply

end module solve_handle
