module fedata
! This module is used to define global values and matrices

   ! Variables
   integer :: ne, nn, nb, nd, np, nm, nk, neqn, bw, neqn_nb
   integer :: n_ptz, ne_ptz
   ! nk = springs, nd = knuder der skal minimeres

   ! Coordinates:
   real(8), dimension(:,:), allocatable :: x

   ! Elements:
   type element_def
      integer, dimension(4) :: ix
      integer :: id, mat, numnode
      real(8) :: ptz
   end type element_def
   type(element_def), dimension(:), allocatable :: element

   ! Material properties:
   type matprop
      real(8) :: young, nu, thk, youngy, shear, dens, area, ep!permativity
      real(8) :: tconv, kcond, alpha, tstart ! added d_28_4_11 for thermal
   end type matprop
   type(matprop), dimension(:), allocatable :: mprop
   
   ! Boundary conditions
   real(8), dimension(:,:), allocatable :: bound, loads, springs, nodes
   real(8) :: accel(2)

   ! Internal heat generation
   real(8) :: Qint
   
   ! Working arrays:
   real(8), dimension(:,:), allocatable :: k, strain, stress
   real(8), dimension(:), allocatable, target :: d
   real(8), dimension(:), allocatable, target :: p

   ! Thermal added d_28_4_11 for thermal
   real(8), dimension(:,:), allocatable :: k_t !thermal stifness matrix -> k+h
   real(8), dimension(:),   allocatable :: t_elem, tn, r_t ! t_elem = element temp, tn = temp in nodes, r_t = thermal load vector
   integer, allocatable :: conv_elem(:) ! bruges til at finde index for de elementer der har konvektion i loads => size(conv_elem,1) = antal elementer med convection

   ! i/o
   character(len = 20) :: filename
   character(len = 40) :: filename_out, dir_out
   logical, parameter :: plot2screen = .true.

   character(len = 40) :: antype, elem_type

   ! Constants
   real(8), parameter ::  pi = 3.1415927 
   real(8), parameter ::  eps0 =  8.854187E-12 !Vacuumpermativity

   ! Parameters
   real(8), parameter :: scale_def = 1.0
   real(8), parameter :: scale_vec = 1.0
!   logical, parameter :: banded = .true.
   integer, parameter :: banded = 2!2
   logical, parameter :: penalty = .false.
   integer, parameter :: ng = 2 !Number of gauss points
   integer, parameter :: lumped_type = 1 ! 1=Yuriy lumped, 2=HRZ mass lumping
   ! Timing
   real(8) start_time, finish_time

   integer :: negative_detjac  
    
   ! TopOpt
   real(8), dimension(:), allocatable ::  compliance_out
   real(8) :: penal, damp_fact

   ! Naboskabsmatrice - Til filtering
   REAL(8), DIMENSION(:,:), allocatable :: neigh

   ! ¤¤¤¤¤¤¤¤¤¤¤ Transient ¤¤¤¤¤¤¤¤¤¤¤¤
   ! TopOpt
   real(8):: nu1, nu2, young1, young2, dens1, dens2

   integer, dimension(:), allocatable :: element_beginning, element_end, center_dofs_x, dof_x_end
   integer, dimension(:,:), allocatable :: element_rand
   INTEGER :: midter_dof_x, midter_dof_y,rand, normal

   ! ABS_rand hvor der ikke er normal-indicies
   type abs_rand_def
      real(8), dimension(2) :: r_enhed
      real(8) :: r
   end type abs_rand_def
   type(abs_rand_def), dimension(:), allocatable :: abs_rand

   ! Mekanisme. De elementer der ligger omkring knudepunktet hvor fjederen virker
   integer, dimension(:), allocatable :: force_elements

   ! sparse format
   integer, dimension(:), allocatable,target :: iK,jK
   real(8), dimension(:), allocatable, target :: sK
   COMPLEX(8), DIMENSION(:), allocatable, target :: sKZ, dZ, Pz
   logical :: harmonic
   real(8) :: mat_vec(20) ! vektor der indeholder materialedata for piezo-elektrisk
   ! mat_vec = [C11 C12 C13 C33 C44 C66 e31 e33 e15 ep11 ep33 layer layer_thk]
   ! integer, parameter :: plate_type = 1

   integer, dimension(:), allocatable,target :: iK2,jK2
   real(8), dimension(:), allocatable, target :: sK2  

   real(8), dimension(:), allocatable :: mvec

   type eigenvalue_def
      logical :: calc, shift
      real(8) :: sigma
      integer :: n_eigen
   end type eigenvalue_def
   type(eigenvalue_def),  allocatable :: eigenvalue

   type INAKTIV_def
      logical :: bool
      integer, allocatable :: act(:), inact(:)
   end type INAKTIV_def
   integer, parameter :: elem_id = 2

end module fedata
