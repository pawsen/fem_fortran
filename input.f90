module input_test

  implicit NONE

  PRIVATE
  PUBLIC :: input, rho_input, read_antype
  PUBLIC :: matrix_input

CONTAINS

  subroutine rho_input(rho_vec,title)

    use fedata

    !  This subroutine reads in the design-variabel (rho). It also tjek that the number of rho in the inputfile
    !  equals ne(number of elements). Notice that the file containing rho, should have a line only saying 
    !  'rho[' and endline '];'(eg disse skal stï¿½ seperat og vï¿½rdierne for rho stï¿½r imellem disse tegn.

    real(8), dimension(:), intent(OUT) :: rho_vec
    character(len=*), intent(in), optional :: title
    integer :: nr
    character(len = 40) :: command, rho_filename
    real(8) :: rvalue
    logical :: fnexist

    if (present(title)) then ! indlï¿½s fil givet fra 'title'-input
       rho_filename = trim(title)
    else ! spï¿½rg om filnavn.

       write (*, *)
       write (*, '("Enter rho input file (<filename> or enter for previous <filename>).")')
       write (*, '("Filename?")')

       !$$$$$$    read (*, '(a)') rho_filename
       !! Dont ask for filename: Comment out the lines above and uncomment the line below
       rho_filename = ' '

       if (rho_filename == ' ') then
          inquire(file = '.rho_filename', exist = fnexist)
          if (.not. fnexist) then
             write (*, *)
             write (*, '("Error: previous filename(rho) does not exist.")') 
             !stop ! i stedet for at stoppe programmet tildeles rho vï¿½rdien 1, hvorfor den ikke fï¿½r nogen betydning.
             write (*, '("Rho has been given the default value 1: eq no influence on deformation")') 
             write (*, *)
             rho_vec = 1.0d0
             return
          endif
          open (10, file = '.rho_filename')
          read (10, *) rho_filename
          close (10)
       end if

    end if

    ! Tjek at fil eksisterer
    inquire(file = trim(rho_filename), exist = fnexist)
    if (.not. fnexist) then 
       write (*, *)
       write (*, '("Error: file(rho) ", a, " does not exist.")') trim(rho_filename) 
       !stop ! i stedet for at stoppe programmet tildeles rho vï¿½rdien 1, hvorfor den ikke fï¿½r nogen betydning.
       write (*, '("Rho has been given the default value 1: eq no influence on deformation")')
       write (*, *)
       rho_vec = 1.0d0
       return
    endif

    open (10, file = '.rho_filename')
    write (10, *) trim(rho_filename)
    close (10)

    write (*, *)
    write (*, '("Reading rho_input file ", a)') trim(rho_filename) 


    open (10, file = trim(rho_filename))

    nr = 0 ! antal rho i filen - dvs length(rho)
    do
       read (10, *) command ! bemï¿½rk *. Den fjerner 'blank spaces', sï¿½ den er MEGET vigtig
       if (trim(command) == 'vec=[') then !Trim: Removes trailing blank characters of a string. Dvs ogsï¿½ MEGET vigtig ved sammenligning af strings
          do
             read (10, *) command
             if (trim(command) == '];') then
            	exit ! slutningen af rho. Gï¿½r ud af fï¿½rste do-lï¿½kke
             else
            	nr = nr + 1
             endif
          enddo
          exit ! Gï¿½r ud af anden do-lï¿½kke
       endif
    end do

    close (10) ! lukker filen, sï¿½ der lï¿½ses fra toppen nï¿½ste gang den ï¿½bnes
    open (10, file = trim(rho_filename))

    if (nr == ne) then
       nr = 0
       do
          read (10, *) command
          if (trim(command) == 'vec=[') then
             do while (nr <= ne)
                read (10, *) command
            	if (trim(command) == '];') then
                   exit
            	else
                   read (command,'(g30.20)') rvalue ! konverterer string til tal
                   nr = nr +1
                   rho_vec(nr) = rvalue
                endif
             enddo
             exit
          endif
       end do
    else
       print*,'Error in rho_input - Number of rhos in file is different than ne(number of elements'
       stop
    endif

    close (10)

  end subroutine rho_input

  subroutine input

    use fedata
    use input_gmsh
    !---------------------------------------------------------------------
    !  This subroutine reads in the input in ansys format.


    integer :: i, j, e, ixsave, eixsave(4), stat
    integer :: nncheck, necheck
    integer :: nen, net
    integer :: ivalue, ivalue2, ivalue3, ifield(19)
    integer :: currentet, currentmp, currentr
    character(len = 12), dimension(10) :: ename
    character(len = 12) :: command, cvalue, cvalue2
    real(8) :: rvalue, rvalue2, rvalue3
    logical :: fnexist
    integer, dimension(:,:), allocatable :: eix

    nn = 0
    nncheck = 0
    net = 0
    ne = 0
    necheck = 0
    nm = 0

    nb = 0
    np = 0
    nk = 0 ! spring
    nd = 0 ! knuder der skal minimeres

    call open_file(10,'')

    ! read in /PREP7 data
    do
       read (10, '(a)', end = 100) command
       if (command == '/PREP7') exit
       if (command == '/' .or. command == ' ') cycle
    end do
    do
       read (10, '(a1)', end = 200) command
       if (command == '/' .or. command == ' ') cycle

       backspace (10)
       read (10, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit
       elseif (command == 'N' .or. command == 'n') then
          nn = nn + 1
          backspace (10)
          read (10, *) command, ivalue
          if (ivalue > nncheck) nncheck = ivalue
       elseif (command == 'ET' .or. command == 'et') then
          net = net + 1
       elseif (command == 'EN' .or. command == 'en') then
          ne = ne + 1
          backspace (10)
          read (10, *) command, ivalue
          if (ivalue > necheck) necheck = ivalue
       elseif (command == 'MP' .or. command == 'mp') then
          backspace (10)
          read (10, *) command, cvalue, ivalue
          if (ivalue > nm) nm = ivalue
       elseif (command == 'D' .or. command == 'd') then
          nb = nb + 1
       elseif (command == 'F' .or. command == 'f') then
          np = np + 1
       elseif (command == 'SFE' .or. command == 'sfe') then
          np = np + 1
          !##############################
       elseif (command == 'K' .or. command == 'k') then
          nk = nk + 1
       elseif (command == 'nodes' .or. command == 'NODES') then
          nd = nd + 1
          !##############################
       endif
    end do


    if (nn == 0) then
       write (*, *) 'Error: No nodes defined.'
       stop
    elseif (nn /= nncheck) then
       write (*, *) 'Error: Node number(s) skipped.'
       stop
    elseif (net == 0) then
       write (*, *) 'Error: No element types defined.'
       stop
    elseif (ne == 0) then
       write (*, *) 'Error: No elements defined.'
       stop
    elseif (ne /= necheck) then
       write (*, *) 'Error: Element number(s) skipped.'
       stop
    elseif (nm == 0) then
       write (*, *) 'Error: No material properties defined.'
       stop
    elseif (nb == 0) then
       write (*, *) 'Error: No supports defined.'
    elseif (np == 0) then
       write (*, *) 'Warning: No loads defined.'
    elseif (nk == 0) then
       write (*, *) 'Message: No springs defined.'
    endif
    write (*, *)
    close (10)

    allocate (x(nn, 2))
    allocate (element(ne))
    allocate (eix(ne,4))
    allocate (mprop(nm))
    allocate (bound(nb, 3), loads(np, 5))
    !############################## added for springs and thermal
    if (nk > 0) then
       print*,'Number of springs',nk
       allocate (springs(nk,5))
    endif
    if (nd > 0) then
       allocate (nodes(nd,2))
       nodes = 0.
    endif
    !allocate (hconv(np)) ! thermal ! Bruges ikke mere, ligger i 5'te sï¿½jle i load i stedet
    ! Make thinkness equal to one as default
    mprop%thk = 1.
    ! And reset all other mprop parameters ! Sigmund 2007 addition
    mprop%young = 0.
    mprop%nu = 0.
    mprop%dens = 0.
    mprop%youngy = 0.
    mprop%shear = 0.
    !############################## !added 28_4_11 for thermal
    mprop%alpha = 0.
    mprop%kcond = 0.
    mprop%tstart = 0.
    !##############################

    ! MODIFIED FOR WINDOWS COMPILER
    do e = 1,ne
       element(e)%ix(1:4) = 0
    enddo

    open (10, file = trim(filename))

    nb = 0
    np = 0
    nk = 0 ! springs
    nd = 0 ! knuder der skal minimeres ved topopt
    currentet = 0
    currentmp = 0
    currentr = 0
    accel=0.0
    !############################## !added 28_4_11 for thermal
    qint=0.0d0       ! skalar, til brug ved internal heat generation
    !hconv = 0.0d0   ! vektor, til brug ved konvektion ! bruges ikke mere
    !##############################


    do
       read (10, '(a1)', end = 200) command
       if (command == '/' .or. command == ' ') cycle

       backspace (10)
       read (10, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit
       elseif (command == 'N' .or. command == 'n') then
          backspace (10)
          read (10, *) command, ivalue, x(ivalue, 1), x(ivalue, 2), rvalue
       elseif (command == 'ET' .or. command == 'et') then
          backspace (10)
          read (10, *) command, ivalue, cvalue
          if (cvalue == 'LINK1' .or. cvalue == 'link1') then
             ename(ivalue) = 'link1'
          elseif (cvalue == 'PLANE42' .or. cvalue == 'plane42') then
             ename(ivalue) = 'plane42'
          elseif (cvalue == 'PLANE42RECT' .or. cvalue == 'plane42rect') then
             ename(ivalue) = 'plane42'
             !############################## added d. 31.8.11 for plade
          elseif (cvalue == 'PLATE' .or. cvalue == 'plate41') then
             ename(ivalue) = 'plate'
             !##############################
          else
             write (*, *) 'Error: Undefined element. ', cvalue
             stop
          endif
          currentet = ivalue
          if (elem_type == 'PLATE') then ! added because I tend to forget to set ET,1,PLATE in structure file.
             ename(ivalue) = 'plate'
          end if
       elseif (command == 'TYPE' .or. command == 'type') then
          backspace (10)
          read (10, *) command, ivalue
          currentet = ivalue
       elseif (command == 'EN' .or. command == 'en') then
          if (currentet == 0) then
             write (*, *) 'Error: No previous element type pointer defined.'
             stop
          endif
          if (currentmp == 0) then
             write (*, *) 'Error: No previous material property pointer defined.'
             stop
          endif
          backspace (10)
          select case( ename(currentet) )
          case( 'link1' )
             read (10, *) command, ivalue, (element(ivalue)%ix(j), j = 1, 2)
             element(ivalue)%mat = currentmp
             element(ivalue)%id  = 1
             element(ivalue)%numnode = 2
          case( 'plane42' )
             read (10, *) command, ivalue, (element(ivalue)%ix(j), j = 1, 4)
             element(ivalue)%mat = currentmp
             element(ivalue)%id  = 2
             element(ivalue)%numnode = 4
             !############################## d. 31.8.11 for plade
          case( 'plate' )
             read (10, *) command, ivalue, (element(ivalue)%ix(j), j = 1, 4)
             element(ivalue)%mat = currentmp
             element(ivalue)%id  = 3
             element(ivalue)%numnode = 4
             !##############################
          case default
             write(*, *) 'Error: Unknown element type'
             stop
          end select
       elseif (command == 'MP' .or. command == 'mp') then
          backspace (10)
          read (10, *) command, cvalue, ivalue, rvalue
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
             mat_vec(10) = rvalue
          elseif (cvalue == 'ep33' .OR. cvalue == 'EP33') THEN
             mat_vec(11) = rvalue
          elseif (cvalue == 'dens_ptz' .OR. cvalue == 'DENS_PTZ') THEN
             mat_vec(12) = rvalue
             !##############################
          else
             write (*, *) 'Warning: Undefined material property: ',cvalue
          endif
          currentmp = ivalue
       elseif (command == 'MAT' .or. command == 'mat') then
          backspace (10)
          read (10, *) command, ivalue
          currentmp = ivalue
       elseif (command == 'R' .or. command == 'r') then 
          backspace (10)
          read (10, *) command, ivalue, rvalue
          if (ename(currentet) == 'link1') then
             ! define truss area
             mprop(ivalue)%area = rvalue
          elseif (ename(currentet) == 'plane42') then
             ! define 4-noded quad element thickness
             mprop(ivalue)%thk = rvalue
             !############################## added  d. 31.8.11 for plade
          elseif (ename(currentet) == 'plate') then
             ! define 4-noded quad plate element thickness
             mprop(ivalue)%thk = rvalue
             !##############################
          else
             write (*, *)
             write (*, '("Error: Undefined real constant (r) card.")')
             stop
          endif
          currentr = ivalue
       elseif (command == 'REAL' .or. command == 'real') then
          backspace (10)
          read (10, *) command, ivalue
          currentr = ivalue
       elseif (command == 'D' .or. command == 'd') then !boundary condition
          backspace (10)
          nb = nb + 1 
          read (10, *) command, ivalue, cvalue, rvalue
          bound(nb, 1) = ivalue ! node id
          bound(nb, 3) = rvalue ! RB, eg V=100 or D=0
          if (cvalue == 'UX' .or. cvalue == 'ux') then
             bound(nb, 2) = 1
          elseif (cvalue == 'UY' .or. cvalue == 'uy') then
             bound(nb, 2) = 2
             !############################## Boundary temperature !added 28_4_11 for thermal
          elseif (cvalue == 'TEMP' .or. cvalue == 'temp') then
             bound(nb, 2) = 3
             !############################## d. 31.8.11 for plade
          elseif (cvalue == 'UZ' .or. cvalue == 'uz') then
             bound(nb, 2) = 1
          elseif (cvalue == 'PSIX' .or. cvalue == 'psix') then
             bound(nb, 2) = 2
          elseif (cvalue == 'PSIY' .or. cvalue == 'psiy') then
             bound(nb, 2) = 3
             !############################## Boundary potential, added 17.10 for piezo
          elseif (cvalue == 'V' .or. cvalue == 'v') then
             bound(nb,2) = 4
          else
             write (*, *) 'Error: Nodal dof not "ux", "uy" or "temp"'
             stop
          endif

       elseif (command == 'F' .or. command == 'f') then ! load
          backspace (10)
          np = np + 1
          read (10, *) command, ivalue, cvalue, rvalue
          loads(np, 1) = 1
          loads(np, 2) = ivalue
          if (cvalue == 'FX' .or. cvalue == 'fx') then
             loads(np, 3) = 1
          elseif (cvalue == 'FY' .or. cvalue == 'fy') then
             loads(np, 3) = 2
             !############################## heat flux in nodes !added 28_4_11 for thermal
          elseif (cvalue == 'HFLU' .or. cvalue == 'hflu') then
             loads(np, 1) = 5 ! "genkendelsesparameter" til buildload
             loads(np, 3) = 3 ! bruges ikke
             !############################## added d. 31.8.11 for plade
             ! NOTICE: nummereringen skal svare til DOF
          elseif (cvalue == 'FZ' .or. cvalue == 'fz') then
             loads(np, 3) = 1
          elseif (cvalue == 'MX' .or. cvalue == 'mx') then
             loads(np, 3) = 2
          elseif (cvalue == 'MY' .or. cvalue == 'my') then
             loads(np, 3) = 3
             !##############################
          else
             write (*, *) 'Error: command not fx, fy, fz or hflu for F in input-routine(processor)'
             stop
          endif
          loads(np, 4) = rvalue
       elseif (command == 'SFE' .or. command == 'sfe') then
          backspace (10)
          np = np + 1
          read (10, *) command, ivalue, ivalue2, cvalue, rvalue2, rvalue
          !############################## Check type of surface load !added 28_4_11 for thermal
          ! note that cvalue2 er ï¿½ndret til rvalue 2, da cvalue er string og rvalue er real
          if (cvalue == 'PRESS' .or. cvalue == 'press' .or. cvalue == 'PRES' .or. cvalue == 'pres') then
             loads(np, 1) = 2                ! Type of load = pressure
             loads(np, 2) = ivalue           ! Element number
             loads(np, 3) = ivalue2          ! Surface number
             loads(np, 4) = rvalue           ! Size of load (pressure)
          elseif (cvalue == 'CONV' .or. cvalue == 'conv') then
             loads(np, 1) = 3                ! Type of load = convection
             loads(np, 2) = ivalue           ! Element number
             loads(np, 3) = ivalue2          ! Surface number
             loads(np, 4) = rvalue           ! Size of load (temperature of surrounding fluid = T_inf)
             loads(np, 5) = rvalue2          ! convective heat transfer coefficient [W/m^2 * K]
             !hconv(np) = rvalue2                ! convective heat transfer coefficient [W/m^2 * K] ! Bruges ikke mere
             ! grunden til at hconv ikke indgår i loads, som loads(np,5), er er loads bliver allokeret som loads(np,4).
             ! Og jeg ved ikke hvad det vil betyde, hvis det ændres til loads(np,4). Så er det gjort alligevel...!
          elseif (cvalue == 'HFLU' .or. cvalue == 'hflu') then
             loads(np, 1) = 4                ! Type of load = heat flux
             loads(np, 2) = ivalue           ! Element number
             loads(np, 3) = ivalue2          ! Surface number
             loads(np, 4) = rvalue           ! Size of load (heat flux) [W/m^2]
          else
             print*,'Error: command not press/pres, conv or hflu for SFE in input-routine(processor)'
             stop
          end if
          ! Internal heat generation = Q. Samme vï¿½rdi for hele strukturen
       elseif (command == 'INHG' .or. command == 'inhg') then
          backspace (10)
          read (10, *) command, rvalue
          qint = rvalue
          !##############################
       elseif (command == 'ACEL' .or. command == 'acel') then
          backspace (10)
          read (10, *) command, rvalue, rvalue2
          accel(1) = rvalue
          accel(2) = rvalue2
          !############################## !added 28_3_11 for springs
       elseif (command == 'K' .or. command == 'k') then
          backspace (10)
          nk = nk + 1
          read (10, *) command, ivalue, cvalue, rvalue, rvalue2
          springs(nk, 1) = 1
          springs(nk, 2) = ivalue ! node the spring is fixed at
          if (cvalue == 'FX' .or. cvalue == 'fx') then
             springs(nk, 3) = 1
          elseif (cvalue == 'FY' .or. cvalue == 'fy') then
             springs(nk, 3) = 2
          else
             write (*, *) 'Error: Nodal dof not "fx" or "fy" for spring'
             stop
          endif
          springs(nk, 4) = rvalue  ! fjeder-konstant
          springs(nk, 5) = rvalue2 ! input/output => 0/1
          !############################## !added 17_6_11. Giver de knuder hvis forskydning skal minimeres
       elseif (command == 'nodes' .or. command == 'NODES') then
          backspace (10)
          nd = nd + 1
          read (10, *) command, ivalue, cvalue
          nodes(nd, 1) = ivalue ! node the spring is fixed at
          if (cvalue == 'FX' .or. cvalue == 'fx') then
             nodes(nd, 2) = 1
          elseif (cvalue == 'FY' .or. cvalue == 'fy') then
             nodes(nd, 2) = 2
          else
             write (*, *) 'Error: Nodal dof not "fx" or "fy" for minimering af knuder'
             stop
          endif

       endif
    end do

    close (10)

    open (10, file = trim(filename))

    ! read in /SOLU data (if it exists)
    DO
       READ (10, '(a)', END = 400) command
       IF (command == '/SOLU') EXIT
       IF (command == '/' .OR. command == ' ') CYCLE
    END DO
    DO
       READ (10, '(a1)', END = 300) command
       IF (command == '/' .OR. command == ' ') CYCLE

       BACKSPACE (10)
       READ (10, *) command
       IF (command == 'FINISH' .OR. command == 'finish') THEN
          EXIT
       ELSEIF (command == 'ANTYPE' .OR. command == 'antype') THEN
          BACKSPACE (10)
          READ (10, *) command, cvalue
          IF (cvalue == 'STATIC' .OR. cvalue == 'static') THEN
             antype = 'STATIC'
          ELSEIF (cvalue == 'STATIC_NL' .OR. cvalue == 'static_nl') THEN
             antype = 'STATIC_NL'
          ELSEIF (cvalue == 'MODAL' .OR. cvalue == 'modal') THEN
             antype = 'MODAL'
          ELSEIF (cvalue == 'ANGLE' .OR. cvalue == 'angle') THEN
             antype = 'ANGLE'
          ELSEIF (cvalue == 'TRANS' .OR. cvalue == 'trans') THEN
             antype = 'TRANS'
             !##############################
          ELSEIF (cvalue == 'THERMAL' .OR. cvalue == 'thermal') THEN
             antype = 'THERMAL'
          ELSEIF (cvalue == 'COUPLED' .OR. cvalue == 'coupled') THEN
             antype = 'COUPLED'
             !##############################
          ELSE
             WRITE (*, *)
             WRITE (*, '("ERROR: only STATIC, STATIC_NL, MODAL, ANGLE and TRANS analyses are implemented.")')
             STOP
          ENDIF
          write(*,*)
          write(*,'("ANTYPE IS GIVEN In model_file. BE sure that it doenst overwrite the value given in _para.txt file!!!")' )
          WRITE (*, *)
       ENDIF
    END DO

    !$$$$$$    antype = 'STATIC' ! default

400 continue

    ! check and reorder element numbering for pressure loads
    do e = 1,ne
       eix(e,:) = element(e)%ix
    end do
    do i = 1, np
       if (loads(i, 1) == 2) then
          e = loads(i, 2)
          if (element(e)%id /= 2) cycle
          if (x(element(e)%ix(1), 1) <= x(element(e)%ix(2), 1) .and. &
               x(element(e)%ix(1), 1) <= x(element(e)%ix(3), 1) .and. &
               x(element(e)%ix(1), 2) <= x(element(e)%ix(3), 2) .and. &
               x(element(e)%ix(1), 2) <= x(element(e)%ix(4), 2)) cycle
          write (*, '("Warning: Reordering element connectivity for element #", i6)') e
          eixsave = element(e)%ix
          do j = 1, 4
             ixsave = eixsave(1)
             eixsave(1:3) = eixsave(2:4)
             eixsave(4) = ixsave
             loads(i, 3) = loads(i, 3) - 1
             if (loads(i, 3) < 1) loads(i, 3) = 4
             if (x(eixsave(1), 1) <= x(eixsave(2), 1) .and. &
                  x(eixsave(1), 1) <= x(eixsave(3), 1) .and. &
                  x(eixsave(1), 2) <= x(eixsave(3), 2) .and. &
                  x(eixsave(1), 2) <= x(eixsave(4), 2)) then
                eix(e,:) = eixsave
                exit
             end if
             if (j == 4) then
                write (*, *) 
                write (*, '("Error: Something wrong with element connectivity.")') 
                stop
             endif
          end do
       endif
    end do
    do e = 1,ne
       element(e)%ix = eix(e,:)
    end do

    return

100 write (*, *)
    write (*, '("ERROR: no /PREP7 in input file")')
    stop

200 write (*, *)
    write (*, '("ERROR: no /PREP7 FINISH in input file")')
    stop

300 write (*, *)
    write (*, '("ERROR: no /SOLU FINISH in input file")')
    stop
    
  end subroutine input

  subroutine matrix_input(title,mat)

    use fedata

    real(8), intent(OUT) :: mat(:,:)
    character(len=*), intent(in), optional :: title
    integer :: i,j, jj, jj_max
    logical :: fnexist
    character(len=40) :: file


    file = trim(dir_out)//trim(title)
    inquire(file = trim(file), exist = fnexist)
    if (.not. fnexist) then
       file = trim(title)
       inquire(file = trim(file), exist = fnexist)
       if (.not. fnexist) then
          print*,'Error fil: ', trim(title), 'kan ikke findes'
          error stop
       end if
    end if
    write (*, *) "Opening file: ", trim(file)
    open (10, file =trim(file) )

    j = size(mat,2)
    jj_max = size(mat,1)
    jj = 0

    do jj = 1,jj_max

       !jj = jj+1
       if (jj> jj_max) then
          print*,'Error: Allokeret matrix er for lille i forhold til filen der indlæses'
          print*,'fil: ', trim(title)
          error stop
       end if
       read (10,*) (mat(jj,i),i=1,j)  ! end = 10 => går til 10 stop

    end do

  close (10)

  return

  end subroutine matrix_input


  subroutine read_antype(sweep)
    use fedata
    use input_gmsh

    real(8), pointer, intent(out), optional :: sweep(:)

    character(len = 40) :: command, rho_filename
    logical :: fnexist
    character(len = 40) :: cvalue
    real(8) :: rvalue
    integer :: ivalue, ivalue2(2), i

    allocate(sweep(3))
    ! DEFAULT values
    sweep(3) = 0 ! dont do sweep by default
    sweep(1) = 2E+4 ! lower bound
    sweep(2) = 6E+4 ! upper bound
    harmonic = .false.

    allocate(eigenvalue)
    eigenvalue%calc = .false.
    eigenvalue%shift = .true.! we use shift as standard. Much better for finding smallest eigenvalues
    eigenvalue%n_eigen = 10 ! default number of eigenvalues to calculate

    elem_type = ' ' ! for checking if elem_type is read from para-file
    print*,'antype_read'

    ! bestem filename
    call open_file(10,'.msh')
    close(10)

    rho_filename = trim(filename)//'_para.txt'


    ! Tjek at fil eksisterer
    inquire(file = trim(rho_filename), exist = fnexist)
    if (.not. fnexist) then 
       write (*, *)
       write (*, '("Error in read_antype: file(parametre) ", a, " does not exist.")') trim(rho_filename) 
       !stop ! i stedet for at stoppe programmet tildeles rho vÃ¦rdien 1, hvorfor den ikke fÃ¥r nogen betydning.
       write (*, '("Antype er sat til static")')
       write (*, *)
       antype = 'STATIC'
       return
    endif

    !clear output file if exist. Needed because we use append(write to bottom of file) for writing to the file in all other writing commands
    dir_out = trim(filename)//'dir/'
    filename_out = trim(dir_out)//trim(filename)
    inquire(file = trim(filename_out)//'.out', exist = fnexist)
    if (fnexist) then
       open (10, file = trim(filename_out)//'.out')
       write(10,'()') ! clear file
       close (10)
    end if
    

    open (10, file = trim(rho_filename))
    do
       read (10, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit ! quit do loop

       ELSEIF (command == 'ANTYPE' .OR. command == 'antype') THEN
          BACKSPACE (10)
          READ (10, *) command, cvalue
          IF (cvalue == 'STATIC' .OR. cvalue == 'static') THEN
             antype = 'STATIC'
          ELSEIF (cvalue == 'PIEZO' .OR. cvalue == 'piezo') THEN
             antype = 'PIEZO'
          ELSEIF (cvalue == 'RATIO_LENGTH_THICKNESS' .OR. &
               (cvalue == 'ratio_length_thickness')) THEN
             antype = 'RATIO_LENGTH_THICKNESS'
          ELSEIF (cvalue == 'EIGEN' .OR. cvalue == 'eigen') THEN
             antype = 'EIGEN'
             eigenvalue%calc = .true.
          ELSEIF (cvalue == 'NONLIN' .OR. cvalue == 'nonlin') THEN
             antype = 'NONLIN'
          ELSEIF (cvalue == 'MODAL' .OR. cvalue == 'modal') THEN
             antype = 'MODAL'
          ELSEIF (cvalue == 'PLATE' .OR. cvalue == 'plate') THEN
             antype = 'PLATE'
          ELSEIF (cvalue == 'KIRCHHOFF' .OR. cvalue == 'kirchhoff') THEN
             antype = 'KIRCHHOFF'
          ELSEIF (cvalue == 'MINDLIN' .OR. cvalue == 'mindlin') THEN
             antype = 'MINDLIN'
          ELSEIF (cvalue == 'THERMAL' .OR. cvalue == 'thermal') THEN
             antype = 'THERMAL'
          ELSEIF (cvalue == 'COUPLED' .OR. cvalue == 'coupled') THEN
             antype = 'COUPLED'
          ELSEIF (cvalue == 'TRANSIENT' .OR. cvalue == 'transient') THEN
             antype = 'TRANSIENT'
          ELSEIF (cvalue == 'TOPTHERM' .OR. cvalue == 'toptherm') THEN
             antype = 'TOPTHERM'
          ELSEIF (cvalue == 'TOPSTRUCT_EIGEN' .OR. cvalue == 'topstruct_eigen') THEN
             antype = 'TOPSTRUCT_EIGEN'
             eigenvalue%calc = .true.
          ELSEIF (cvalue == 'TOPSTRUCT' .OR. cvalue == 'topstruct') THEN
             antype = 'TOPSTRUCT'
          ELSEIF (cvalue == 'TOPTHERMSTRUCT' .OR. cvalue == 'topthermstruct') THEN
             antype = 'TOPTHERMSTRUCT'
          ELSEIF (cvalue == 'TOPTHERMSTRUCT_hard' .OR. cvalue == 'topthermstruct_hard') THEN
             antype = 'TOPTHERMSTRUCT_hard'
          ELSE
             print*,'ERROR - ANTYPE ikke genkendt - sÃ¦tter antype = static'
             antype = 'STATIC'
          END IF
       ELSEIF (command == 'ELEM_TYPE' .OR. command == 'elem_type') THEN
          BACKSPACE (10)
          READ (10, *) command, cvalue
          IF (cvalue == 'PLATE' .OR. cvalue == 'plate') THEN
             elem_type = 'PLATE'
          ELSEIF (cvalue == 'PLANE42' .OR. cvalue == 'plane42') THEN
             elem_type = 'PLANE42'
          ELSEIF (cvalue == 'PLATE_GMSH' .OR. cvalue == 'plate_gmsh') THEN
             elem_type = 'PLATE_GMSH'
          ELSEIF (cvalue == 'PLANE_GMSH' .OR. cvalue == 'plane_gmsh') THEN
             elem_type = 'PLANE_GMSH'
          ELSE
             print*,'ERROR - ELEM_TYPE ikke genkendt - sets elem_type = PLANE42'
             elem_type = 'PLANE42'
          endif
       ELSEIF (command == 'HARMONIC' .OR. command == 'harmonic') THEN
          BACKSPACE (10)
          READ (10, *) command, ivalue
          if (ivalue == 1) then
             harmonic = .true.
          else
             harmonic = .false.
          end if
       ELSEIF (command == 'SWEEP' .OR. command == 'sweep') THEN
          BACKSPACE (10)
          if (present(sweep)) then
             READ (10, *) command, sweep(1), sweep(2)
          end if
       ELSEIF (command == 'SWEEP_BOOL' .OR. command == 'sweep_bool') THEN
          BACKSPACE (10)    
          READ (10, *) command, sweep(3)
          if (sweep(3) == 1) then
             harmonic = .true. !has to be harmonic because the harmonic bool is used to initialize complex system vectors
          end if
       ELSEIF (command == 'EIGENVALUE' .OR. command == 'eigenvalue') THEN
          BACKSPACE (10)
          READ (10, *) command, ivalue
          !if (ivalue == 1) then
          !eigenvalue%calc = .true.
          eigenvalue%n_eigen = ivalue
          !end if
       ELSEIF (command == 'SHIFT' .OR. command == 'shift') THEN
          BACKSPACE (10)
          READ (10, *) command, ivalue
          ! We use shift as standard.
          if (ivalue == 0) then
             eigenvalue%shift = .false.
          else !read shift value
             BACKSPACE (10)
             READ (10, *) command, ivalue, rvalue
             eigenvalue%sigma = AINT(rvalue,8) ! round the shift value towards zero
          end if
       end IF
    end do

    !check if elem_type is assigned an element type. If not, then give it a default value.
    if (trim(elem_type)=='') then 
       elem_type = 'PLANE42'
       write(*,*) "elem_type is not read from file. sets elem_type = PLANE42" 
    end if


    close (10) ! lukker filen, sÃ¥ der lÃ¦ses fra toppen nÃ¦ste gang den Ã¥bnes

  end subroutine read_antype

end module input_test
