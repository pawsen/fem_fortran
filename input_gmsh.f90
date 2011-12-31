module input_gmsh
  implicit none

  PRIVATE ::  read_material_data
  PUBLIC :: read_gmsh, open_file

CONTAINS

  subroutine read_gmsh

    use fedata
    use numeth

    integer, parameter :: IOgmsh= 24                              !Gmsh mesh file fid
    integer:: i, j, ie, icel, ibnd, iloop

    !
    ! there are 19 gmsh element types: \
    !
    !    1 : 2-node line
    !    2 : 3-node triangle (face)
    !    3 : 4-node quadrangle (face)
    !    4 : 4-node tetrahedron
    !    5 : 8-node hexahedron (eg cube)
    !    6 : 6-node triangular-prism
    !    7 : 5-node pyramid
    !
    !  8-14: 'second-order' elements.  Ref Gmsh manual.
    !   15 : 1-node point
    ! 16-19: more second-order FEM elements
    !
    ! the nodes/vertices for each element are read into the
    ! v array.
    !
    ! each element can have several tags.
    ! the first tag gives the physical_group number.
    ! all elements on the same boundary have the same physical_group number.
    !

    integer, parameter :: element_type(19) = &
         (/ 2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13 /)

    character(len=72)  :: c_input1, c_input2, c_input3 ! til sammenlignig af $-strings i msh filen
    character(len = 12) :: command, cvalue, cvalue2
    integer :: n_elements, ielement, ielement_type, n_tags, n_nodes, ivalue
    integer ::  n_ptz_buf, kk
    integer, allocatable :: index_nb(:)
    integer :: tags(64), v(27)
    real(8) :: rvalue, gnd, potential
    real(8), allocatable :: bound_buf(:,:)

    kk = 0

    ! da read_antype er kaldt i main, er 'filename' kendt. 
    allocate (mprop(2))
    call read_material_data(gnd,potential)
    mprop(2) = mprop(1)

    open(IOgmsh, file = trim(filename)//'.msh')

    nn = 0
    ne = 0
    nb = 0
    np = 0
    ne_ptz = 0
    do
       read (IOgmsh,*) command
       if (command == '$Elements') then
          read (IOgmsh,*) n_elements
          do ie=1,n_elements
             read(IOgmsh,*) ielement, ielement_type, &
                  n_tags, (tags(i),i=1,n_tags)

             select case(ielement_type)
             case(3) ! 4-node quad
                ne = ne+1
                ! only for plates
                if (( (tags(1) >= 100001) .and. (tags(1) <= 100070)) .or. &
                     tags(1) == 10015) then
                   n_ptz_buf = MOD(tags(1),100000)
                   if (n_ptz_buf > n_ptz) then
                      n_ptz =n_ptz_buf
                   end if
                   ne_ptz = ne_ptz +1
                end if
             case(1) ! line
                nb = nb +2
             case(15) !single node, point
                select case(tags(1))
                case(1000:1001)!point load
                   np = np+1
                case(10000:10009)! RB
                   nb = nb+1
                end select
             end select

          end do
          exit
       else
          cycle
       end if
    end do

    close (IOgmsh)

    open(IOgmsh, file = trim(filename)//'.msh')
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$MeshFormat')
    read(IOgmsh,*) c_input1,c_input2,c_input3
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$EndMeshFormat')


    ! read the nodes from the .msh file
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$Nodes')
    read(IOgmsh,*) n_nodes

    nn = n_nodes
    allocate (x(nn, 2))
    allocate (element(ne))
    allocate (bound_buf(nb, 3), loads(np, 5))

    do iloop=1,n_nodes
       read(IOgmsh,*) ivalue,x(ivalue, 1), x(ivalue, 2), rvalue
    enddo
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$EndNodes')

    ! read the elements from the .msh file
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$Elements')
    read(IOgmsh,*) n_elements

    ivalue = 0
    nb = 0
    np = 0
    do ie=1,n_elements

       !
       read(IOgmsh,*) ielement, ielement_type, &
            n_tags, (tags(i),i=1,n_tags),&
            (v(i),i=1,element_type(ielement_type))

       select case(ielement_type)
       case(1) ! 2-node line
          nb = nb +2
          bound_buf(nb-1:nb, 1) = v(1:2) ! node id
          select case(tags(1)) ! find boundary type
          case(10000) ! gnd
             bound_buf(nb-1:nb,3) = gnd ! V =0
             bound_buf(nb-1:nb,2) = 10 ! id
          case(10001) ! cos
             !potential is read from para-file
             bound_buf(nb-1:nb,3) = potential
             bound_buf(nb-1:nb,2) = 11
          case(10002) ! sin
             bound_buf(nb-1:nb,3) = potential
             bound_buf(nb-1:nb,2) = 12
          case(10005) ! RB, UX/UZ=0
             bound_buf(nb-1:nb,3) = 0 ! prescribed value
             bound_buf(nb-1:nb,2) = 1 !UZ/UX
          case(10006) ! RB, UY/PSIX=0
             bound_buf(nb-1:nb,3) = 0
             bound_buf(nb-1:nb,2) = 2 !PSIX/UY
          case(10007) ! RB, PSIY=0
             bound_buf(nb-1:nb,3) = 0
             bound_buf(nb-1:nb,2) = 3 !PSIY
          end select
       case(15) ! 1-node point
          select case(tags(1))
          case(1000:1001)!point load
             np = np + 1
             loads(np, 1) = 1
             loads(np, 2) = v(1) ! node
             loads(np, 4) = -1d0 ! magnitude and direction of load
             select case(tags(1)) ! find point load type
             case(1000) !f_x
                loads(np, 3) = 1
             case(1001) !f_y
                loads(np, 3) = 2
             case default
                print*,'point load type not known, gmsh'
                error stop
             end select
          case(10000:10009) !point RB
             nb = nb +1
             bound_buf(nb,1) = v(1) ! node id
             select case(tags(1))
             case(10005) ! RB, UX/UZ=0
                bound_buf(nb,3) = 0 ! prescribed value
                bound_buf(nb,2) = 1 !UZ/UX
             case(10006)
                bound_buf(nb,3) = 0 ! prescribed value
                bound_buf(nb,2) = 2 !UZ/UX
             end select
          end select
       case(3) ! 4-node quad
          ivalue = ivalue +1
          do j = 1,4
             element(ivalue)%ix(j) =  v(j)
          end do
          element(ivalue)%mat = 1! vi giver alle elementer mat == 1. Det betyder ikke noget
          element(ivalue)%numnode = 4

          select case(tags(1)) ! find material type
          case(10010) !pcb domain
             element(ivalue)%id  = 6
          case(10011) !pcb domain !skal ikke optimeres
             element(ivalue)%mat = 2
             if (elem_type == 'PLANE_GMSH') then
                element(ivalue)%id  = 2
             else
                element(ivalue)%id  = 6
             end if
          case(10015) !piezo domain, Plate
             if (elem_type == 'PLANE_GMSH') then
                element(ivalue)%id  = 2
             else
                element(ivalue)%id  = 3
             end if
          case(10020) !piezo domain, beam
             element(ivalue)%id  = 2
          case(100001 :) ! find piezo-nummeret
             n_ptz_buf = MOD(tags(1),100000)

             element(ivalue)%ptz = 2*Pi/n_ptz*(n_ptz_buf-1)!-1 fordi ring.geo er skrevet med c-syntax, dvs første piezo-element har nummer 0
             !print*,'phi= ',  element(ivalue)%ptz

             ! Polarisering
             select case (MOD(n_ptz_buf,8))
             case(1)
                element(ivalue)%id =  5 ! right
                kk = kk+1
             case(2)
                element(ivalue)%id =  4 ! left
                kk = kk+1
             case(3)
                element(ivalue)%id =  4
                kk = kk+1
             case(4)
                element(ivalue)%id =  5
                kk = kk+1
             case(5)
                element(ivalue)%id =  4
                kk = kk+1
             case(6)
                element(ivalue)%id =  5
                kk = kk+1
             case(7)
                element(ivalue)%id =  5
                kk = kk+1
             case(0)
                element(ivalue)%id =  4
                kk = kk+1
             end select

!!$          case(10021) ! left 'polarized'
!!$             element(ivalue)%mat =  3
!!$             n_ptz = npz +1
!!$          case(10022) ! right 'polarized'
!!$             element(ivalue)%mat =  4
!!$             n_ptz = npz +1
          end select

       case default

          print*,'Forkert elementtype, input_gmsh.f90, type=',ielement_type
          stop
       end select
    end do

 print*,'kk = ', kk

    ! Remove dublicate nodes. Notice that if duplicate nodes have different id, it's not removed
    allocate(index_nb(nb))
    call index_dups(INT(bound_buf(:,1:2)),index_nb,i)

    nb = i
    allocate(bound(i, 3))
    bound = bound_buf( index_nb(1:i),: )

    print*,'n_ptz= ',n_ptz
    print*,'ne_ptz= ',ne_ptz
    print*,'ne= ',ne
    print*,'nb= ',nb


    if (.not. (ne == ivalue)) then
       print*,'error ne /= ivalue'
    end if
    read(IOgmsh,*) c_input1
    call check_input_character(c_input1,'$EndElements')


    close(IOgmsh)

  end subroutine read_gmsh



  subroutine read_material_data(gnd,pot)

    use fedata

    real(8), INTENT(OUT) :: pot, gnd
    integer, parameter :: IOgmsh= 24 
    character(len = 12) :: command, cvalue, cvalue2
    integer :: ivalue,i
    real(8) :: rvalue, rvec(10)

    rvec = 0d0

    gnd = 0d0
    pot = 0d0

    !##############################
    ! read material data from seperate file

    ! Make thinkness equal to one as default
    mprop%thk = 1.
    ! And reset all other mprop parameters ! Sigmund 2007 addition
    mprop%young = 0.
    mprop%nu = 0.
    mprop%dens = 0.
    mprop%youngy = 0.
    mprop%shear = 0.
    mprop%alpha = 0.
    mprop%kcond = 0.
    mprop%tstart = 0.
    mprop%ep = 0d0
    !##############################
    mat_vec = 0d0
    mat_vec(12) = 1 ! standard antal lag

    !call open_file(IOgmsh,'_para.txt')
    open(IOgmsh, file = trim(filename)//'_para.txt')


    do

       read (IOgmsh, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit
       elseif (command == 'MP' .or. command == 'mp') then
          backspace (IOgmsh)
          read (IOgmsh, *) command, cvalue, ivalue, rvalue
          if (cvalue == 'EX' .or. cvalue == 'ex') then
             mprop(ivalue)%young = rvalue
          elseif (cvalue == 'PRXY' .or. cvalue == 'prxy') then
             mprop(ivalue)%nu = rvalue
          elseif (cvalue == 'DENS' .or. cvalue == 'dens') then
             mprop(ivalue)%dens = rvalue
          elseif (cvalue == 'EY' .OR. cvalue == 'ey') THEN
             mprop(ivalue)%youngy = rvalue
          elseif (cvalue == 'GXY' .OR. cvalue == 'gxy') THEN
             mprop(ivalue)%shear = rvalue
          elseif (cvalue == 'EP' .OR. cvalue == 'ep') THEN
             mprop(ivalue)%ep = rvalue*eps0
          elseif (cvalue == 'GND' .OR. cvalue == 'gnd') THEN
             gnd = rvalue
          elseif (cvalue == 'POT' .OR. cvalue == 'pot') THEN
             pot = rvalue
             !##############################  Heat transfer material properties  !added 28_4_11 for thermal
          elseif (cvalue == 'ALPX' .OR. cvalue == 'alpx') THEN       ! alpha = Thermal expansion coefficient
             mprop(ivalue)%alpha = rvalue
          elseif (cvalue == 'KCOND' .OR. cvalue == 'kcond') THEN     ! k = thermal conductivity
             mprop(ivalue)%kcond = rvalue
          elseif (cvalue == 'TSTART' .OR. cvalue == 'tstart') THEN   ! T0 = start temperature
             mprop(ivalue)%tstart = rvalue
             !##############################  piezo-elektrisk element
          elseif (cvalue == 'C11' .OR. cvalue == 'c11') THEN   
             mat_vec(1) = rvalue
          elseif (cvalue == 'C12' .OR. cvalue == 'c12') THEN
             mat_vec(2) = rvalue
          elseif (cvalue == 'C13' .OR. cvalue == 'c13') THEN
             mat_vec(3) = rvalue
          elseif (cvalue == 'C33' .OR. cvalue == 'c33') THEN
             mat_vec(4) = rvalue
          elseif (cvalue == 'C44' .OR. cvalue == 'c44') THEN
             mat_vec(5) = rvalue
          elseif (cvalue == 'C66' .OR. cvalue == 'c66') THEN
             mat_vec(6) = rvalue
          elseif (cvalue == 'e31' .OR. cvalue == 'E31') THEN
             mat_vec(7) = rvalue
          elseif (cvalue == 'e33' .OR. cvalue == 'E33') THEN
             mat_vec(8) = rvalue
          elseif (cvalue == 'e15' .OR. cvalue == 'E15') THEN
             mat_vec(9) = rvalue
          elseif (cvalue == 'ep11' .OR. cvalue == 'EP11') THEN
             mat_vec(10) = rvalue*eps0 ! ganger med vacuumpermativity
          elseif (cvalue == 'ep33' .OR. cvalue == 'EP33') THEN
             mat_vec(11) = rvalue*eps0 ! ganger med vacuumpermativity
          elseif (cvalue == 'dens_ptz' .OR. cvalue == 'DENS_PTZ') THEN
             mat_vec(12) = rvalue
          else
             write (*, *) 'Warning: Undefined material property: ',cvalue
          endif
       elseif (command == 'R' .or. command == 'r') then 
          backspace (IOgmsh)
          read (IOgmsh, *) command, ivalue, rvalue
          mprop(ivalue)%thk = rvalue
       elseif (command == 'LAYER' .OR. command == 'layer') THEN
          backspace (IOgmsh)
          read (IOgmsh, *) command, ivalue
          mat_vec(13) = ivalue
       elseif (command == 'LAYER_THK' .OR. command == 'layer_thk') THEN
          backspace (IOgmsh)
          read (IOgmsh, *) command, (rvec(i),i=1,int(mat_vec(13)))
          do i=1,int(mat_vec(13))
             mat_vec(13+i) = rvec(i)
          end do
!!$             mat_vec(13) = rvalue
!!$             mat_vec(14) = rvalue1
!!$             mat_vec(15) = rvalue2
       elseif (command == 'frekvens' .OR. command == 'FREKVENS') THEN
          backspace (IOgmsh)
          read (IOgmsh, *) command, rvalue
          mat_vec(20) = rvalue
       endif
    end do

    close(IOgmsh)

  end subroutine read_material_data

  !------------------------------------------------------------------------------------

  subroutine check_input_character(c1,c2)

    implicit none

    character (len=*) :: c1, c2

    if( c1(1:len(c2)) /= c2 )then
       write(*,*)  'error reading Gmsh input file: ',&
            'the following two characters should be the ',&
            'same but differ ',c1(1:len(c2)),c2
       stop
    endif

  end subroutine check_input_character

  subroutine open_file(fid,ext)

    use fedata

    implicit none

    integer, INTENT(IN) :: fid
    character(len=*), intent(in) :: ext
    logical :: fnexist


   
    !write (*, '("File(s) in current directory: ")')
    !call system('dir')
    !write (*, *)
    !write (*, '("Enter fem input file (<filename> or enter for previous <filename>).")')
    !write (*, '("Filename?")')
    !$$$$$$    read (*, '(a)') filename
    ! Dont ask for filename: Comment out the line above and uncomment the line below
    filename = ' '
    if (filename == ' ') then
       inquire(file = '.fem_filename', exist = fnexist)
       if (.not. fnexist) then
          write (*, *)
          write (*, '("Error: previous filename does not exist.")') 
          stop
       endif
       open (10, file = '.fem_filename')
       do
          read (10, *,end = 11) filename! goto 11 when end of file is reached
          if (SCAN( trim(filename), '!', .false. ) == 1 ) then ! if '!' is the first character in filename, then take next line in .fem_filename
             cycle
          else
             exit
          end if
       end do
11     close (10)
    end if

    write (*, *)
    ! try firt to open filename.ext
    ! If this file doesn't exist, then try open
    inquire(file = trim(filename)//trim(ext), exist = fnexist)
    if (fnexist) then
       write (*, 20) "Opening file: ", trim(filename),trim(ext)
       open (fid, file = trim(filename)//trim(ext))
    else
       inquire(file = trim(filename), exist = fnexist)
       if (fnexist) then
          write (*, 20) "Opening file: ", trim(filename)
          open (fid, file = trim(filename))
       else
          write (*, *)
          write (*, '("Error in input_gmsh.f90: file ", a, " does not exist.")') trim(filename) 
          stop
       end if
    endif
20  format(a,a,a)
    ! DONT overwrite .fem_filename
!!$    open (10, file = '.fem_filename')
!!$    write (10, *) trim(filename)
!!$    close (10)


  end subroutine open_file


end module input_gmsh

