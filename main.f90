program main

  use fedata
  !   use processor
  use fea
  use input_test
  use nonlin
  use thermal
  use topology
  use transient
  use plot_routiner
  use input_test
  use plate
  use file_init
  use piezo
  use input_gmsh
  use arpack

  implicit none

  integer :: flag, read_rho
  real(8) :: rho_min
  real(8), dimension(:) ,pointer :: sweep
  real(8), allocatable, dimension(:) :: rho

  mat_vec(13) = 1! eet lag som standard
  call read_antype(sweep)
  if ( ((elem_type == 'PLATE_GMSH').or.(elem_type == 'PLANE_GMSH'))) then
  ! &.and.( antype == 'PIEZO')) then
     call read_gmsh
  else
     ! Read model data, Also read antype fra filename_para.txt
     call input
  end if

  read_rho = 0 ! 0: fra	1: indl�s rho fra fil ifm beregning af displ/varmefordeling
  if (read_rho == 1) then
     allocate(rho(ne))
     call rho_input(rho,trim(filename)//'dir/'//trim(filename)//'_rho'//'.m') !linux_fejl
  end if

  write(*,10) "antype: ",antype
  write(*,10) "elem_type: ",elem_type
10 format(a,a)
  write(*,*)


  ! unders�g om elementer er ens
  !  call elements_equal(flag)
  !  write(*,10) "flag: ",flag, "antype: ",antype
  !  10 format(a,i3,3x,a,a)
!!$  call initial_fea
!!$  call eigen
!!$
!!$  stop

  flag = 0
  if (antype == 'STATIC') then
     ! Calculate linear displacements
     call initial_fea
     if (allocated(rho)) then
        penal = 3d0
        rho_min = 0.001
        call displ(flag,rho,rho_min)
     else
        call displ(flag)
     end if
  else if (eigenvalue%calc) then
     ! calculate eigenvalues and modeshapes
     if (elem_type == 'PLATE_GMSH') then
        call initial_piezo
        call arpack_init
     else
        call initial_fea
        call arpack_init
     end if
  else if (antype == 'PIEZO') then
     !harmonic = .true.
     print*,'Harmonic= ', harmonic

     call initial_piezo
     if (harmonic) then
        call sweep_piezo(0,(/0d0/),sweep)
     else
        call displ_piezo(flag)
     end if
  else if (antype == 'TOPSTRUCT_EIGEN') then
     call initial_piezo
     call TopOpt(flag)
  else if (antype == 'NONLIN') then
     ! Displacement for geometric nonlinear problems.
     flag = 0 ! da ke beregnes p� baggrund af forkydninger for det enkelte element, kan flag = 1 ikke bruges
     call initial_fea
     if (allocated(rho)) then
        penal = 3d0
        rho_min = 0.001
        call non_lin(flag,rho,rho_min)
     else
        call non_lin(flag)
     end if
  else if (antype == 'MODAL') then
     ! snydt:D
  else if (antype == 'KIRCHHOFF' .or. antype =='PLATE') then
     !linux_fejl
     call displ_plate
     !$$$$$$    else if (antype == 'MINDLIN') then
     !$$$$$$         call initial_plate
     !$$$$$$         call displ_plate
  else if (antype == 'THERMAL') then
     ! Calculate temperatures
     call initial_t
     if (allocated(rho)) then
        rho_min = 0.001d0
        penal = 3d0
        call therm_displ(flag,rho_min,rho)
     else
        call therm_displ(flag)
     end if
  else if (antype == 'COUPLED') then
     ! Calculate temperatures, stresses and strains  
     call initial_t
     !call initial
     if (allocated(rho)) then
        rho_min = 0.001d0
        penal = 3d0
        call therm_displ(flag,rho_min,rho)
        call initial_fea
        call displ(flag,rho,rho_min)
     else
        call therm_displ(flag)
        call initial_fea
        call displ(flag)
     end if
  else if (antype == 'TRANSIENT') then
     ! TopOpt for thermal problems
     flag = 1
     write(*,*)
     print*,'flag: ',flag
     call initial_transient(flag)
     if (allocated(rho)) then
        call TopOpt_trans(flag,rho)
     else
        call TopOpt_trans(flag)
     end if
  else if (antype == 'TOPTHERM') then
     ! TopOpt for thermal problems
     call initial_t
     call TopOpt(flag)
  else if (antype == 'TOPTHERMSTRUCT') then
     ! TopOpt for koblet mekanisme med termisk stress
     allocate(t_elem(ne))
     call initial_t
     call initial_fea
     call TopOpt(flag)
  else if (antype == 'TOPTHERMSTRUCT_hard') then
     ! TopOpt for mekanisme med "hardcoded" temp-stigning givet ved t_elem
     call initial_fea
     call TopOpt(flag)
  else if (antype == 'TOPSTRUCT') then
     ! TopOpt for structurel problems (eq force-inverter, compliance,... )
     call initial_fea
     call TopOpt(flag)
  else
     print*,'ERROR: Wrong analysis type. Only STATIC, MODAL, THERMAL, COUPLED, TOPSTRUCT, TOPTHERM are implemented.'
  end if

end program main
