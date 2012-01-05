MODULE file_init

  ! 

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: inputfile_trans, inputfile, parameter_input
CONTAINS

  subroutine inputfile_trans(file,rand,normal,loes_type,output_type &
       ,rho,n_iter,nmax,deltaT,cL,young1,dens1,nu1,young2,dens2,nu2)

   use transient !udkommenteret d. 21/10. Uden at vide hvad konsekvensen er
   use plot_routiner, only : output_vector

    integer, INTENT(INOUT) :: file,rand,normal, loes_type,output_type
    real(8), intent(IN) :: rho(:)
    integer, INTENT(OUT) :: n_iter
    integer, optional, INTENT(OUT) :: nmax
    real(8), optional, INTENT(OUT) ::deltaT,cL,young1,dens1,nu1,young2,dens2,nu2
    real(8) :: cT, dens,young,nu, lambda,mu, deltaT_max, real_parameters(10), output_vec(20),deltaT_faktor
    integer :: nr, read_para, cL_set, int_parameters(10)

    !indlæs	materialedata fra fil
    call parameter_input_transient(real_parameters,int_parameters)
    deltaT = 0.0d0

    write(*,*)
    read_para = int_parameters(4)
    if (read_para == 1) then ! parametre er sat i 'filename_para.txt'
       rand = int_parameters(1)
       normal = int_parameters(2)
       cL_set = int_parameters(3)
       loes_type = int_parameters(5)
       output_type = int_parameters(6)
       file = int_parameters(7)
       nmax = int_parameters(8)
       n_iter = int_parameters(9)

       deltaT_faktor = real_parameters(9)
       if (cL_set == 1) then!materialeparametre bestemmes ud fra cL og cT
          cL = real_parameters(1)
          cT = real_parameters(2)
          dens1 = real_parameters(7)
          young1 = cT**2*dens1*(3*cL**2-4*cT**2)/(cL**2-cT**2)
          nu1 = (cL**2-2*cT**2)/(cL**2-cT**2)/2
          young2 = young1
          dens2 = dens1
          nu2 = nu1
       else
          young1 = real_parameters(3)
          young2 = real_parameters(4)
          nu1 = real_parameters(5)
          nu2 = real_parameters(6)
          dens1 = real_parameters(7)
          dens2 = real_parameters(8)
       end if
       write(*,*)'materialedata indlæst fra fil!'
       write(*,*)
    else ! angiv materialedata her
       write(*,*)'materialedata angivet manuelt i fortran kode!'
       write(*,*)
       select case( file )
       case(1)! Jacobs fil
          !$$$$$$             deltaT = 0.10d0
          deltaT = 0.050d0 ! til longsving5
          nmax = 2000

          ! Materiale-parametre
          young1 = 1.0!2.0/3.0d0
          dens1 = 1.0d0
          nu1 = 0
          young2 = 1.5!2.0d0
          dens2 = 1.3!2.0d0
          nu2 = 0.0

          !$$$$$$         cL = 1
          !$$$$$$         cT = cL /2
          !$$$$$$         ! Materiale-parametre
          !$$$$$$         dens1 = 1.0d0
          !$$$$$$         young1 = cT**2*dens1*(3*cL**2-4*cT**2)/(cL**2-cT**2)
          !$$$$$$         nu1 = (cL**2-2*cT**2)/(cL**2-cT**2)/2
          !$$$$$$         !nu1 = 0.0d0!
          !$$$$$$         young2 = young1
          !$$$$$$         dens2 = dens1
          !$$$$$$         nu2 = nu1

       case(2)! eksploderende kilde
          !deltaT = 0.00670d0
          deltaT = 0.0030d0
          nmax = 2000

          young1 = 1.0!2.0/3.0d0
          dens1 = 1.0d0
          nu1 = 0
          young2 = 1.5!2.0d0
          dens2 = 1.3!2.0d0
          nu2 = 0.0

          cL = 1d0
          cT = cL/2d0
          ! Materiale-parametre
          dens1 = 1.0d0
          young1 = cT**2*dens1*(3*cL**2-4*cT**2)/(cL**2-cT**2)
          nu1 = (cL**2-2*cT**2)/(cL**2-cT**2)/2
          !$$$$$$         nu1 = 0.0d0!
          young2 = young1
          dens2 = dens1
          nu2 = nu1
       case(3)! disk med harmonisk kraft
          deltaT = 0.10d0
          nmax = 2000

          ! Materiale-parametre
          young1 = 1.0!2.0/3.0d0
          dens1 = 1.0d0
          nu1 = 0
          young2 = 1.5!2.0d0
          dens2 = 1.3!2.0d0
          nu2 = 0.0

       end select
       n_iter = 500
       deltaT_faktor = 0.7d0
    end if

    call deltaT_init(rho,deltaT_max)
    if ((deltaT > deltaT_max) .and. deltaT /= 0 )then
       print*,'Error, deltaT er for stor. Den må maks være ',deltaT_max
       stop
    end if

    deltaT = deltaT_faktor*deltaT_max ! sætter deltaT til at være 70% af den masimale værdi

    ! print materialeoplysninger
    nr = 10 ! elementnummeret skal ikke være i enden hvor de inaktive elementer sidder. Ellers er det uden betydning da alle elementer er ens
    dens = rho(nr)*dens1+(1.0-rho(nr))*dens2
    young = rho(nr)*young1+(1.0-rho(nr))*young2
    nu = rho(nr)*nu1+(1.0-rho(nr))*nu2

    lambda = young*nu/ ((1.0+nu)*(1.0-2.0*nu)) 
    mu = young/(2.0*(1.0+nu))! shear modulus (G)

    cL = dsqrt( (lambda+2.0*mu) / dens)
    cT = dsqrt( mu / dens)

    write(*,*)
    write(*,30) "deltaT: ",deltaT, "nmax: ",nmax,"n_iter: ",n_iter
    write(*,10) "cL: ",cL,"cT: ",cT, "young1: ",young1, "young2: ",young2
    write(*,20) "nu1: ",nu1, "nu2: ",nu2, "dens1: ",dens1, "dens2: ",dens2
    write(*,40) "rand: ",rand, "normal: ",normal, "file: ",file
10  format(a,f5.2,3x,a,f5.2,3x,a,f5.2,3x,a,f5.2)
20  format(a,f5.2,3x,a,f5.2,3x,a,f5.2,3x,a,f5.2)
30  format(a,e10.3,3x,a,i5,3x,a,i5)
40  format(a,i2,3x,a,i2,3x,a,i2)
    write(*,*)
    ! for info om format se: http://www.math.hawaii.edu/lab/197/fortran/fort3.htm

    real_parameters = 0d0
    int_parameters = 0
    ! udskriver parametre til fil
    real_parameters(1) = cL
    real_parameters(2) = cT
    real_parameters(3) = young1
    real_parameters(4) = young2
    real_parameters(5) = nu1
    real_parameters(6) = nu2
    real_parameters(7) = dens1
    real_parameters(8) = dens2
    real_parameters(9) = deltaT_faktor
    real_parameters(10) = deltaT


    if (read_para == 1) then ! parametre er sat i 'filename_para.txt'
       int_parameters(1) = cL_set
       int_parameters(2) = read_para
    end if

    int_parameters(3) = rand
    int_parameters(4) = normal
    int_parameters(5) = loes_type
    int_parameters(6) = output_type
    int_parameters(7) = file
    int_parameters(8) = nmax
    int_parameters(9) = n_iter

    output_vec(1:10) = real_parameters
    output_vec(11:20) = real(int_parameters)
    call output_vector(output_vec,'parametre')


  end subroutine inputfile_trans

  !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

  subroutine parameter_input_transient(real_parameters,int_parameters)
    use fedata

    !  This subroutine reads in the design-variabel (rho). It also tjek that the number of rho in the inputfile
    !  equals ne(number of elements). Notice that the file containing rho, should have a line only saying 
    !  'rho[' and endline '];'(eg disse skal stå seperat og værdierne for rho står imellem disse tegn.

    real(8), dimension(:), intent(OUT) :: real_parameters
    integer, dimension(:), intent(OUT) :: int_parameters
    character(len = 40) :: command, rho_filename
    real(8) :: rvalue
    integer :: ivalue
    logical :: fnexist
    real(8) :: cL_r, cT_r, young1_r, young2_r,nu1_r, nu2_r, dens1_r,dens2_r,deltaT_faktor
    integer :: cL_set, rand_r, normal_r, read_para, filetype, loes_type, output_type, n_iter, nmax

    real_parameters = 0
    int_parameters = 0


    rho_filename = trim(filename)//'_para.txt'

    ! Tjek at fil eksisterer
    inquire(file = trim(rho_filename), exist = fnexist)
    if (.not. fnexist) then 
       write (*, *)
       write (*, '("Error: file(parametre) ", a, " does not exist.")') trim(rho_filename) 
       !stop ! i stedet for at stoppe programmet tildeles rho værdien 1, hvorfor den ikke får nogen betydning.
       write (*, '("parametre er sat som i file-init")')
       write (*, *)
       read_para = 0
       return
    endif

    open (10, file = trim(rho_filename))
    do
       read (10, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit ! quit do løkke
       elseif (command == 'cL') then
          backspace (10)
          read (10, *) command, rvalue
          cL_r = rvalue
       elseif (command == 'cT') then
          backspace (10)
          read (10, *) command, rvalue
          cT_r = rvalue
       elseif (command == 'young1') then
          backspace (10)
          read (10, *) command, rvalue
          young1_r = rvalue
       elseif (command == 'young2') then
          backspace (10)
          read (10, *) command, rvalue
          young2_r = rvalue
       elseif (command == 'nu1') then
          backspace (10)
          read (10, *) command, rvalue
          nu1_r = rvalue
       elseif (command == 'nu2') then
          backspace (10)
          read (10, *) command, rvalue
          nu2_r = rvalue
       elseif (command == 'dens1') then
          backspace (10)
          read (10, *) command, rvalue
          dens1_r = rvalue
       elseif (command == 'dens2') then
          backspace (10)
          read (10, *) command, rvalue
          dens2_r = rvalue
       elseif (command == 'deltaT_faktor') then
          backspace (10)
          read (10, *) command, rvalue
          deltaT_faktor = rvalue
       elseif (command == 'rand') then
          backspace (10)
          read (10, *) command, ivalue
          rand_r = ivalue
       elseif (command == 'normal') then
          backspace (10)
          read (10, *) command, ivalue
          normal_r = ivalue
       elseif (command == 'cL_set') then
          backspace (10)
          read (10, *) command, ivalue
          cL_set = ivalue
       elseif (command == 'read_para') then
          backspace (10)
          read (10, *) command, ivalue
          read_para = ivalue
       elseif (command == 'file') then
          backspace (10)
          read (10, *) command, ivalue
          filetype = ivalue
       elseif (command == 'loes_type') then
          backspace (10)
          read (10, *) command, ivalue
          loes_type = ivalue
       elseif (command == 'output_type') then
          backspace (10)
          read (10, *) command, ivalue
          output_type = ivalue
       elseif (command == 'nmax') then
          backspace (10)
          read (10, *) command, ivalue
          nmax = ivalue
       elseif (command == 'n_iter') then
          backspace (10)
          read (10, *) command, ivalue
          n_iter = ivalue
       endif
    end do

    int_parameters(1) = rand_r
    int_parameters(2) = normal_r
    int_parameters(3) = cL_set
    int_parameters(4) = read_para
    int_parameters(5) = loes_type
    int_parameters(6) = output_type
    int_parameters(7) = filetype
    int_parameters(8) = nmax
    int_parameters(9) = n_iter

    real_parameters(1) = cL_r
    real_parameters(2) = cT_r
    real_parameters(3) = young1_r
    real_parameters(4) = young2_r
    real_parameters(5) = nu1_r
    real_parameters(6) = nu2_r
    real_parameters(7) = dens1_r
    real_parameters(8) = dens2_r
    real_parameters(9) = deltaT_faktor

    close (10) ! lukker filen, så der læses fra toppen næste gang den åbnes


  end subroutine parameter_input_transient


  subroutine inputfile(file,filter_type,solver_type,problem_type,vol_type,rmin_type,info,save_rho,n_iter,&
       stop_krit,animation,penal,damp_fact,max_vol,rmin,rho_min,tol,movelimit)

    integer, INTENT(IN) :: file
    integer, intent(INOUT) :: filter_type,solver_type,problem_type,vol_type,rmin_type,info
    integer, intent(INOUT) :: save_rho,n_iter,stop_krit,animation
    real(8), intent(inout) :: penal,damp_fact,max_vol,rmin,rho_min,tol,movelimit

    real(8) :: vol_frac

    real(8) :: real_parameters(10)
    integer :: int_parameters(15), read_para

    !indlæs	materialedata fra fil
    call parameter_input(real_parameters,int_parameters)

    read_para = int_parameters(9)
    select case(read_para)
    case(1)! parametre indlæses fra fil
       !file = int_parameters(1)
       filter_type = int_parameters(1)
       solver_type = int_parameters(2)
       problem_type = int_parameters(3)
       vol_type = int_parameters(4)
       rmin_type = int_parameters(5)
       info = int_parameters(6)
       save_rho = int_parameters(7)
       n_iter = int_parameters(8)
       stop_krit = int_parameters(10)
       animation = int_parameters(11)

       penal = real_parameters(1)
       damp_fact = real_parameters(2)
       max_vol = real_parameters(3)
       rmin = real_parameters(4)
       rho_min = real_parameters(5)
       tol = real_parameters(6)
       movelimit = real_parameters(7)

    case(0)! angives manuelt her
       select case( file )
       case(1)! Brigde1
          penal = 3
          max_vol = 2000.d0
          vol_frac = 0.5
          rmin  = 2.2d0	! Filter Radius. [m]
          damp_fact   = 0.4d0	! Parameter to make solution stable. (eta = 0.65)

       case(2)! MBB5
          penal = 3.0
          !$$$$$$             max_vol = 100.0
          vol_frac = 0.5
          damp_fact = 0.5
          rmin = 0.8
       case(3)! TO1
          penal = 3
          max_vol = 30.d0
          vol_frac = 0.5
          rmin  = 0.6d0	! Filter Radius. [m]
          damp_fact   = 0.4d0	! Parameter to make solution stable. (eta = 0.65)
       case(4)! MBB fra filter-artikel
          penal = 3
          vol_frac = 0.50d0
          rmin  = 3.0	! Filter Radius. [m]
          damp_fact   = 0.4d0	! Parameter to make solution stable. (eta = 0.65)
       end select
       select case( vol_type )
       case( 2 ) ! så er vol_frac specificeret
          !	max_vol = vol_frac
       end select
       n_iter = 1000
       select case(problem_type)
       case(1)! movelimit sat ned pga mekanisme-design
          movelimit = 0.10d0
       case default
          movelimit = 0.20d0
       end select
    end select


    write(*,*)
    write(*,10) "penal: ",penal, "max_vol: ",max_vol,"rmin: ",rmin,"movelimit: ",movelimit
    write(*,20) "tol: ",tol,"rho_min: ",rho_min, "damp_fact: ",damp_fact, "n_iter: ",n_iter,"stop_krit: ",stop_krit
    write(*,30) "filter: ",filter_type, "solver: ",solver_type, "problem: ",problem_type, "rmin_type: ",rmin_type
10  format(a,f5.2,3x,a,f5.2,3x,a,f5.2,3x,a,f5.2)
20  format(a,e10.3,3x,a,e10.2,3x,a,f5.2,3x,a,i5,3x,a,i3)
30  format(a,i3,3x,a,i3,3x,a,i3,3x,a,i3)
    write(*,*)
    ! for info om format se: http://www.math.hawaii.edu/lab/197/fortran/fort3.htm


  end subroutine inputfile


  subroutine parameter_input(real_parameters,int_parameters,rand_int,rand_val)
    use fedata

    !  This subroutine reads in the design-variabel (rho). It also tjek that the number of rho in the inputfile
    !  equals ne(number of elements). Notice that the file containing rho, should have a line only saying 
    !  'rho[' and endline '];'(eg disse skal stå seperat og værdierne for rho står imellem disse tegn.

    real(8), dimension(:), intent(OUT) :: real_parameters
    integer, dimension(:), intent(OUT) :: int_parameters
    real(8) , optional, intent(out) :: rand_val(4)
    integer , optional, intent(out) :: rand_int(4)
    character(len = 40) :: command, rho_filename
    real(8) :: rvalue, rvalue2, rvalue3, rvalue4
    integer :: ivalue, ivalue2, ivalue3, ivalue4
    logical :: fnexist
    real(8) :: p_factor,damp_factor, max_vol, rmin, rho_min, tol, movelimit
    integer :: file, filter_type, solver_type, problem_type, vol_type, rmin_type, info
    integer :: save_rho, n_iter, read_para, stop_krit, anim

    real_parameters = 0d0
    int_parameters = 0
    damp_factor = 0
    command = ''

    rho_filename = trim(filename)//'_para.txt'

    ! Tjek at fil eksisterer
    inquire(file = trim(rho_filename), exist = fnexist)
    if (.not. fnexist) then 
       write (*, *)
       write (*, '("Error: file(parametre) ", a, " does not exist.")') trim(rho_filename) 
       !stop ! i stedet for at stoppe programmet tildeles rho værdien 1, hvorfor den ikke får nogen betydning.
       write (*, '("parametre er sat som i file-init")')
       write (*, *)
       read_para = 0
       return
    endif

    open (10, file = trim(rho_filename))
    do
       read (10, *) command
       if (command == 'FINISH' .or. command == 'finish') then
          exit ! quit do løkke
          ! real(8) - dvs data mv.
       elseif (command == 'penal') then
          backspace (10)
          read (10, *) command, rvalue
          p_factor = rvalue
       elseif (command == 'rho_min') then
          backspace (10)
          read (10, *) command, rvalue
          rho_min = rvalue
       elseif (command == 'rmin') then
          backspace (10)
          read (10, *) command, rvalue
          rmin = rvalue
       elseif (command == 'damp_fact') then
          backspace (10)
          read (10, *) command, rvalue
          damp_factor = rvalue
       elseif (command == 'tol') then
          backspace (10)
          read (10, *) command, rvalue
          tol = rvalue
       elseif (command == 'max_vol') then
          backspace (10)
          read (10, *) command, rvalue
          max_vol = rvalue
       elseif (command == 'movelimit') then
          backspace (10)
          read (10, *) command, rvalue
          movelimit = rvalue

          ! integer, dvs parametre mv.
       elseif (command == 'filter') then
          backspace (10)
          read (10, *) command, ivalue
          filter_type = ivalue
       elseif (command == 'solver') then
          backspace (10)
          read (10, *) command, ivalue
          solver_type = ivalue
       elseif (command == 'problem') then
          backspace (10)
          read (10, *) command, ivalue
          problem_type = ivalue
       elseif (command == 'vol_type') then
          backspace (10)
          read (10, *) command, ivalue
          vol_type = ivalue
       elseif (command == 'rmin_type') then
          backspace (10)
          read (10, *) command, ivalue
          rmin_type = ivalue
       elseif (command == 'info') then
          backspace (10)
          read (10, *) command, ivalue
          info = ivalue
       elseif (command == 'n_iter') then
          backspace (10)
          read (10, *) command, ivalue
          n_iter = ivalue
       elseif (command == 'read_para') then
          backspace (10)
          read (10, *) command, ivalue
          read_para = ivalue
       elseif (command == 'save_rho') then
          backspace (10)
          read (10, *) command, ivalue
          save_rho = ivalue
       elseif (command == 'stop_krit') then
          backspace (10)
          read (10, *) command, ivalue
          stop_krit = ivalue
       elseif (command == 'anim') then
          backspace (10)
          read (10, *) command, ivalue
          anim = ivalue
       elseif (command == 'rand') then
          backspace (10)
          if (present(rand_int)) then
             read (10, *) command, ivalue, ivalue2, ivalue3, ivalue4
             rand_int = [ivalue, ivalue2, ivalue3, ivalue4]
          else ! skal være her, ellers bliver det en uendelig løkke
             read (10, *) command
          end if
       elseif (command == 'rand_val') then
          backspace (10)
          if (present(rand_val)) then
             read (10, *) command, rvalue, rvalue2, rvalue3, rvalue4
             rand_val = [rvalue, rvalue2, rvalue3, rvalue4]
          else
             read (10, *) command
          end if
       endif
    end do

    close (10) ! lukker filen, så der læses fra toppen næste gang den åbnes

    int_parameters(1) = filter_type
    int_parameters(2) = solver_type
    int_parameters(3) = problem_type
    int_parameters(4) = vol_type
    int_parameters(5) = rmin_type
    int_parameters(6) = info
    int_parameters(7) = save_rho
    int_parameters(8) = n_iter
    int_parameters(9) = read_para
    int_parameters(10) = stop_krit
    int_parameters(11) = anim

    real_parameters(1) = p_factor
    real_parameters(2) = damp_factor
    real_parameters(3) = max_vol
    real_parameters(4) = rmin
    real_parameters(5) = rho_min
    real_parameters(6) = tol
    real_parameters(7) = movelimit




  end subroutine parameter_input


END MODULE file_init
