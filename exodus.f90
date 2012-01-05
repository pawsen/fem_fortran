module exodus

  implicit none


  public :: exodus_write
  public :: exodus_init, exodus_write_node, exodus_write_elem
  public :: exodus_write_time, exodus_finalize
  private
  include 'exodusII.inc'

  ! initialisering
  integer:: exoid, num_dim, num_nodes, num_elem, num_elem_blk
  integer:: num_elem_in_block(10), num_nodes_per_elem(10), numattr(10)
  integer:: ebids(10) !unique id
  integer:: ierr

  ! data output, Initialisering
  integer:: num_glo_vars, num_nod_vars, num_elem_vars

CONTAINS

  subroutine exodus_write(eigenval,eigenvec)
    !
    ! This is a test program for the Fortran binding of the EXODUS II
    ! database write routines.
    !
    use fedata
    use numeth

    real(8), intent(IN), optional ::eigenval(:,:), eigenvec(:,:)


    integer:: num_node_sets, num_side_sets
    integer:: cpu_word_size, io_word_size

    ! data output, Initialisering
    integer:: whole_time_step, num_time_steps

    real(4):: elem_var_vals(ne), nodal_var_vals(nn)
    real(4):: elem_val, time_value

    real(4):: glob_var_vals(100)
    real(4):: x1(nn), y1(nn), z1(nn)
    real(4):: frekvens, omega, period, lower, upper, time

    character*(MXSTLN):: coord_names(3)
    character*(MXSTLN):: cname
    character*(MXSTLN):: var_names(4)
    character(len = 40):: filename_exodus !denne må åbenbart godt være længere end MXSTLN

    integer :: ierr, e
    integer :: kk, i, j,  m

    integer, allocatable :: connect(:,:), e_num(:)



    cpu_word_size = 0
    io_word_size = 0

    ! create EXODUS II files
    call make_dir
    filename_exodus = trim(filename_out)//'.exo'
    exoid = excre (filename_exodus,EXCLOB, cpu_word_size, io_word_size, ierr)

    ! initialize file with parameters
    num_dim = 3
    num_nodes = nn
    num_elem = ne
    if ((ne == ne_ptz) .or. (ne_ptz == 0)) then
       num_elem_blk = 1 ! different blok set(eg. groups of elements)
    else
       num_elem_blk = 2
    end if
    num_node_sets =0
    num_side_sets = 0
    call expini (exoid, "FORTRAN PIEZO_RING ", num_dim, num_nodes,& 
         num_elem, num_elem_blk, num_node_sets,&
         num_side_sets, ierr)


    x1(1:nn) = real(x(:,1),4)
    y1(1:nn) = real(x(:,2),4)
    z1(1:nn) = 0.0
    ! write nodal coordinates values and names to database
    call expcor (exoid, x1(1:num_nodes), y1(1:num_nodes), z1(1:num_nodes) , ierr)

    coord_names(1) = "xcoor"
    coord_names(2) = "ycoor"
    coord_names(3) = "zcoor"
    call expcon (exoid, coord_names, ierr)

    !write element block parameters
    ! pcb-elements
    if ((ne == ne_ptz) .or. (ne_ptz == 0)) then
       num_elem_in_block(1) = ne
    else
       num_elem_in_block(1) = ne-ne_ptz
    end if
    num_nodes_per_elem(1) = 4

    ! piezo-elements
    num_elem_in_block(2) = ne_ptz
    num_nodes_per_elem(2) = 4

    do i=1,num_elem_blk
       ebids(i) = i ! unique ID
       numattr(i) = 0 ! number of attributes ! NOT working
    end do

    cname = "shell"
    do i=1,num_elem_blk
       call expelb (exoid,ebids(i),cname,num_elem_in_block(i),&
            num_nodes_per_elem(i),numattr(i),ierr)
    end do

    ! write element connectivity, s16
    allocate(e_num(num_elem_blk))
    allocate(connect(MAX(ne-ne_ptz,ne_ptz)*4,num_elem_blk) )
    e_num(:) = 0
    do e=1,ne
       select case( element(e)%id )
       case(2) !quad
          e_num(1) = e_num(1)+1
          do i = 1,4
             connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
          end do
       case(3:5)
          if (num_elem_blk == 1) then !der er kun ptz materiale
             e_num(1) = e_num(1)+1
          else
             e_num(2) = e_num(2)+1
          end if
          do i = 1,4
             if (num_elem_blk == 1) then
                connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
             else
                connect(4*(e_num(2)-1) + i,2) = (element(e)%ix(i))
             end if
          end do
       case(6)
          e_num(1) = e_num(1)+1
          do i = 1,4
             connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
          end do
       end select
    end do
    do i=1,num_elem_blk
       call expelc (exoid, ebids(i), connect(1:e_num(i)*4,i), ierr)
    end do

    !write results variables parameters and names
    !nodes
    if ( (antype == 'PIEZO') .and. (eigenvalue%calc .eqv. .false.)) then
       num_nod_vars = 4
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"
       var_names(4) = "POT"

    elseif((antype == 'EIGEN')) then
       num_nod_vars = 3
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"
    end if

    call expvp (exoid, "n", num_nod_vars, ierr)!s143, write parameters
    call expvan (exoid, "n", num_nod_vars, var_names, ierr) !s147, write names

    if (.false.) then
       num_elem_vars = 1
       var_names(1) = "rho"
       call expvp (exoid, "e", num_nod_vars, ierr)
       call expvan (exoid, "e", num_nod_vars, var_names, ierr)
    end if



    if ( harmonic) then

       ! bestem frekvens og lav lineær intepolation for tiden i een periode.
       frekvens = mat_vec(20)
       if (int(frekvens) == 0) then
          omega = 0d0
          period = 0d0
          num_time_steps = 2
          print*
          print*,'Frekvensen er 0, så jeg håber DU har sat potentialet på som cos.'
          print*,'Hvis ikke, er både forskydninger og potential nul i .exo filen'
          print*
       else
          omega = 2*pi* frekvens
          period = 1./frekvens
          num_time_steps = 100
       end if
       lower = 0
       upper = period

       
       do i = 1, num_time_steps
          time = lower+(real(i)-1.0)*(upper-lower)/(real(num_time_steps)-1.0)! lineær intepolation

          ! write time value
          call exptim (exoid, i, time, ierr)!s152

          ! write node values
          do kk = 1, num_nod_vars
             select case(kk)
             case(3)! displacement z-direction
                nodal_var_vals =  real(dZ(1:neqn:3),4)*cos(omega*time)-AIMAG(dZ(1:neqn:3))*sin(omega*time)
             case(4)! potential
                nodal_var_vals =  real(dZ(neqn+1:neqn+nn),4)*cos(omega*time)-AIMAG(dZ(neqn+1:neqn+nn))*sin(omega*time)
             case default
                nodal_var_vals = 0.0
             end select
             call expnv (exoid, i, kk, num_nodes,nodal_var_vals, ierr)!s187

          end do
          
          ! updata exodus-file. This also flush netCDF I/O buffer., s.3 as3
          call exupda (exoid, ierr)
       end do
    else if((antype == 'EIGEN')) then
       num_time_steps = size(eigenval,1) ! number of converged eigenvalues
       do i = 1, num_time_steps
          time =  SQRT(eigenval(i,1))/(2*pi)!egenfrekvens. Kun real-part
          ! write time value
          call exptim (exoid, i, time, ierr)!s152
          ! write node values
          do kk = 1, num_nod_vars
             select case(kk)
             case(1)! x-direction
                if (element(1)%id == 2) then
                   nodal_var_vals = real(eigenvec(1:neqn-1:2,i),4)
                else
                   nodal_var_vals = 0.
                end if
             case(2)! y-direction 
                if (element(1)%id == 2) then
                   nodal_var_vals = real(eigenvec(2:neqn:2,i),4)
                else
                   nodal_var_vals = 0.
                end if
             case(3)
                if (element(1)%id /= 2) then
                   nodal_var_vals = real(eigenvec(1:neqn:3,i),4)
                else
                   nodal_var_vals = 0.
                end if
             case(4)! potential)
             case default
                nodal_var_vals = 0.0
             end select
             call expnv (exoid, i, kk, num_nodes,nodal_var_vals, ierr)!s187
          end do
          ! updata exodus-file. This also flush netCDF I/O buffer., s.3 as3
          call exupda (exoid, ierr)
       end do
    else

       num_time_steps = 1
       do i = 1, num_time_steps
          time = 1
          ! write time value
          call exptim (exoid, i, time, ierr)!s152
          
          ! write node values
          do kk = 1, num_nod_vars
             select case(kk)
             case(1)! x-direction
                if (element(1)%id == 2) then
                   nodal_var_vals = real(d(1:neqn-1:2),4)
                else
                   nodal_var_vals = 0.
                end if
             case(2)! y-direction 
                if (element(1)%id == 2) then
                   nodal_var_vals = real(d(2:neqn:2),4)
                else
                   nodal_var_vals = 0.
                end if
             case(3)! displacement z-direction
                if (element(1)%id == 2) then
                   nodal_var_vals = 0.
                else
                   nodal_var_vals = real(d(1:neqn:3),4)
                end if
                
             case(4)! potential
                nodal_var_vals =  real(d(neqn+1:neqn+nn),4)
             case default
                nodal_var_vals = 0.0
             end select

             call expnv (exoid, i, kk, num_nodes,nodal_var_vals, ierr)
          end do

          ! write time value
          call exptim (exoid, i, i, ierr)!s152

          ! updata exodus-file. This also flush netCDF I/O buffer., s.3 as3
          call exupda (exoid, ierr)
       end do
    end if

    !close the EXODUS files
    call exclos (exoid, ierr)


  end subroutine exodus_write

  subroutine exodus_init

    use numeth
    use fedata

    integer:: num_node_sets, num_side_sets
    integer:: cpu_word_size, io_word_size

    real(4):: x1(nn), y1(nn), z1(nn)

    character*(MXSTLN):: coord_names(3)
    character*(MXSTLN):: cname
    character*(MXSTLN):: var_names(4), elem_var_names(4)
    character(len = 40):: filename_exodus !denne må åbenbart godt være længere end MXSTLN


    integer :: e
    integer :: kk, i, j,  m

    integer, allocatable :: connect(:,:), e_num(:)

    cpu_word_size = 0
    io_word_size = 0

    ! create EXODUS II files
    call make_dir
    filename_exodus = trim(filename_out)//'.exo'
    exoid = excre (filename_exodus,EXCLOB, cpu_word_size, io_word_size, ierr)


    ! initialize file with parameters
    num_dim = 3
    num_nodes = nn
    num_elem = ne
    if (ne == ne_ptz .or. ne_ptz == 0) then
       num_elem_blk = 1 ! different blok set(eg. groups of elements)
    else
       num_elem_blk = 2
    end if
    num_node_sets =0
    num_side_sets = 0
    call expini (exoid, "FORTRAN PIEZO_RING ", num_dim, num_nodes,& 
         num_elem, num_elem_blk, num_node_sets,&
         num_side_sets, ierr)


    x1(1:nn) = real(x(:,1),4)
    y1(1:nn) = real(x(:,2),4)
    z1(1:nn) = 0.0
    ! write nodal coordinates values and names to database
    call expcor (exoid, x1(1:num_nodes), y1(1:num_nodes), z1(1:num_nodes) , ierr)

    coord_names(1) = "xcoor"
    coord_names(2) = "ycoor"
    coord_names(3) = "zcoor"
    call expcon (exoid, coord_names, ierr)

    !write element block parameters
    ! pcb-elements
    if (ne == ne_ptz .or. ne_ptz == 0) then
       num_elem_in_block(1) = ne
    else
       num_elem_in_block(1) = ne-ne_ptz
    end if
    num_nodes_per_elem(1) = 4

    ! piezo-elements
    num_elem_in_block(2) = ne_ptz
    num_nodes_per_elem(2) = 4

    do i=1,num_elem_blk
       ebids(i) = i ! unique ID
       numattr(i) = 0 ! number of attributes ! NOT working
    end do

    cname = "shell"
    do i=1,num_elem_blk
       call expelb (exoid,ebids(i),cname,num_elem_in_block(i),&
            num_nodes_per_elem(i),numattr(i),ierr)
    end do

    ! write element connectivity, s16
    allocate(e_num(num_elem_blk))
    allocate(connect(MAX(ne-ne_ptz,ne_ptz)*4,num_elem_blk) )
    e_num(:) = 0
    do e=1,ne
       select case( element(e)%id )
       case(2) !quad
          e_num(1) = e_num(1)+1
          do i = 1,4
             connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
          end do
       case(3:5)
          if (num_elem_blk == 1) then !der er kun ptz materiale
             e_num(1) = e_num(1)+1
          else
             e_num(2) = e_num(2)+1
          end if
          do i = 1,4
             if (num_elem_blk == 1) then
                connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
             else
                connect(4*(e_num(2)-1) + i,2) = (element(e)%ix(i))
             end if
          end do
       case(6)
          e_num(1) = e_num(1)+1
          do i = 1,4
             connect(4*(e_num(1)-1) + i,1) = (element(e)%ix(i))
          end do
       end select
    end do
    do i=1,num_elem_blk
       call expelc (exoid, ebids(i), connect(1:e_num(i)*4,i), ierr)
    end do

    !write results variables parameters and names
    !nodes
    if ( (antype == 'PIEZO') .and. ( eigenvalue%calc .eqv. .false.)) then
       num_nod_vars = 4
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"
       var_names(4) = "POT"

    elseif(antype == 'EIGEN') then
       num_nod_vars = 3
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"

    else if ( antype == 'STATIC') then
       num_nod_vars = 3
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"
    else if( (antype == 'TOPSTRUCT') .or. (antype == 'TOPSTRUCT_EIGEN') ) then
       num_nod_vars = 3
       var_names(1) = "DISPLX"
       var_names(2) = "DISPLY"
       var_names(3) = "DISPLZ"
       
       num_elem_vars = 1
       elem_var_names(1) = "rho"
       call expvp (exoid, "e", num_elem_vars, ierr)
       call expvan (exoid, "e", num_elem_vars, elem_var_names, ierr)
    end if
    

    call expvp (exoid, "n", num_nod_vars, ierr)!s143, write parameters
    call expvan (exoid, "n", num_nod_vars, var_names, ierr) !s147, write names

  end subroutine exodus_init

  subroutine exodus_write_time(i,time)
    !write time as the last thing before next timestep
    !This is due to the flus-buffer command
    
    real(8), intent(in) :: time
    integer, intent(in) :: i ! step number

    ! write time value
    call exptim (exoid, i, real(time,4), ierr)!s152
    ! updata exodus-file. This also flush netCDF I/O buffer., s.3 as3
    call exupda (exoid, ierr)

  end subroutine exodus_write_time


  subroutine exodus_write_node(i, node_val)

    use fedata

    integer, intent(in) :: i ! step number
    real(8), intent(in) :: node_val(:)
    integer :: kk, j, m
    real(4):: node_var_vals(num_nodes)


    ! write node values
    do kk = 1, num_nod_vars
       select case(kk)
       case(1)! x-direction
          if (element(1)%id == 2) then
             node_var_vals = real(node_val(1:neqn-1:2),4)
          else
             node_var_vals = 0.
          end if
       case(2)! y-direction 
          if (element(1)%id == 2) then
             node_var_vals = real(node_val(2:neqn:2),4)
          else
             node_var_vals = 0.
          end if
       case(3)! displacement z-direction
          if (element(1)%id == 2) then
             node_var_vals = 0.0
          else
             node_var_vals =  real(node_val(1:neqn:3),4)
          end if
       case(4)! potential
          node_var_vals =  real(node_val(neqn+1:neqn+nn),4)
       case default
          node_var_vals = 0.0
       end select

       call expnv (exoid, i, kk, num_nodes,node_var_vals, ierr)!s187

    end do

  end subroutine exodus_write_node

  subroutine exodus_write_elem(i, elem_val)

    use fedata

    integer, intent(in) :: i ! step number
    real(8), intent(in) :: elem_val(:)
    integer :: kk, j, m
    real(4):: elem_var_vals(ne)
    real(4):: elem_value, time_value

    !write element values
    do  kk = 1, num_elem_vars
       do j = 1, num_elem_blk
          elem_var_vals = elem_val
          ! do m = 1, num_elem_in_block(j)
          !    select case(j)
          !    case(1)!PCB
          !       elem_value = 10*(-1)**(i-1)!-real((i-1)*999.,4)
          !    case(2)!PIEZO
          !       elem_value = 10*(-1)**i
          !    end select
          !    elem_var_vals(m) = elem_value
          ! end do

          call expev (exoid, i, kk, ebids(j),num_elem_in_block(j), elem_var_vals, ierr)!s 163
       end do
    end do

  end subroutine exodus_write_elem

  subroutine exodus_finalize

    !close the EXODUS files
    call exclos (exoid, ierr)

  end subroutine exodus_finalize

end module exodus
