!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMOUN MODULES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  in this file, there are the main and most important parameters and variables used in this !!!!!!!
!                                            program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Copyright (C) 2018 Quentin Vanhaelen


MODULE variables


logical, dimension(1) :: mol_added 
type binary_tree_structure
  integer :: childgauche
  integer :: childdroit
  integer ::  parents
  integer :: label
  real, pointer :: psum => NULL() 
 end type binary_tree_structure
 type (binary_tree_structure), allocatable, dimension(:)  :: tree_structure
 type stochioupdate
   integer, pointer :: species => NULL()  
 end type stochioupdate
type (stochioupdate), allocatable, dimension(:,:)  :: number_species,number_species_product,number_species_tilde
integer,  allocatable, dimension(:) ::  stochio_number,stochio_number_product,stochio_number_tilde
 real    :: t !  temporal  variable
 real ::  sum
 real, allocatable,  dimension(:)  :: d_mu ! absolute delay for each reaction
 real,  allocatable, dimension(:) :: react_t_delay !  vector containing the varying delay for each reaction
 integer, allocatable, dimension(:) :: mu_delayed ! vector for labels of the reaction having a delay
 integer,  allocatable, dimension(:) ::  binary_tree
integer,  allocatable, dimension(:) ::  leftchild,  rightchild
!  assuming that N  chemical species could react through M reaction channels:
real, allocatable, dimension(:)    :: c_nu !  probability that the reaction nu occur
real, allocatable, dimension(:)    :: a_nu !  a_nu() = c_nu()*h_nu()
 real, allocatable, dimension(:) ,target   ::  p_sum
real, allocatable, dimension(:)  :: h_nu !  stochiometric coefficient
integer, allocatable, dimension(:) ,target :: Y !  vector of the N chemical  species
real, allocatable, dimension(:,:) ::  matrix_result, matrix_result_global
integer, allocatable, dimension(:,:) :: matrix_nbre_transition
real, allocatable, dimension(:,:) :: matrix_A_transition, matrix_Adt_transition
real, allocatable, dimension(:,:) :: matrix_nbre_transition_global, matrix_A_transition_global, matrix_Adt_transition_global
integer,  allocatable, dimension(:,:) ::  update_matrix,  update_matrix_tilde, update_matrix_product  !  matrix for the update of the number of molecules
integer,  allocatable, dimension(:,:) ::  update,update_tilde, update_product
 real    ::  a_0  
 integer :: mu  ! global channel  label  
 integer :: mu_yes ! channel  label delayed reaction only
integer :: mu_no ! channel  label non delayed reaction only
logical, allocatable, dimension(:) :: consuming 
real, allocatable, dimension(:) :: kinetic_rate,vector_mu_tilde
integer :: etat_initial, etat_final

real, allocatable, dimension(:)  :: Y_concentration !  vector of the N chemical  species before scaling
END MODULE variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE parameters

 integer :: nrun, nruntotal !  number of runs of set of chemical  reactions to  be done
 integer :: freq_run_saving  !frequence de sauvegarde des fichiers every X run
 integer ::  ntotal !  number of chemical  reactions to  be simulated
 integer :: N_species ! number of chemical  species in the system (see below for model specification)
 parameter ( N_species = 20)
 integer :: M_react_channel ! number of reaction channels   (see below for model specification)
 parameter (M_react_channel = 36)
 integer :: M_delayed !  number of reaction including a delay
 parameter (M_delayed = 100000)

! frequency of saving for chemical  events
integer  :: saving_frequency


! size of the output file
integer :: number_events_saved
parameter (number_events_saved= 50000)
 real ::  scaling 
integer :: last_node_tree
integer :: deep_tree
 parameter (deep_tree = 5)
END MODULE parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE statistics
real :: SSA_duration
real :: SSA_physical_time_av 
END MODULE statistics






















































































































