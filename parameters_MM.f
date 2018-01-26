module MM_input

use type_m         
use constants_m
use parameters_m    , only : driver , Pop_size , N_generations , Top_Selection , Pop_range , Mutation_rate , Mutate_Cross
use MM_types        , only : MM_molecular , MM_system 

type(MM_system)                  :: MM
type(MM_molecular) , allocatable :: species(:) 

real*8                 :: temperature, pressure, cutoff_radius, thermal_relaxation_time, pressure_relaxation_time, damping_Wolf
integer                :: MM_log_step, MM_frame_step , N_of_CGSteps 
logical                :: read_velocities, gmx_input_format , Selective_Dynamics
character (5)          :: Units_MM
character (len=6)      :: OPT_driver
character (len=11)     :: driver_MM 
character (len=14)     :: thermostat
type(integer_interval) :: nmd_window

logical , parameter :: T_ = .true. , F_ = .false. 

contains

!================================
 subroutine Define_MM_Environment
!================================
implicit none

!------------------------------------------------------------------------------
! SYSTEM  INFO
!
  MM % N_of_molecules = 1                   ! <== total number of molecules
  MM % N_of_species   = 1                   ! <== total number of species

  CALL allocate_species( MM % N_of_species )

!------------------------------------------------------------------------------
! repeat the following information filling for all the different species ...
!
  species(1) % residue         = "PSB"      ! <== Residue label for species i ; character(len3)
  species(1) % N_of_molecules  = 1          ! <== Number of molecules of species i
  species(1) % N_of_atoms      = 66         ! <== Number of atoms comprosing a single molecule of species i
  species(1) % flex            = T_         ! <== Flexible : T_ , F_
  
  Selective_Dynamics = F_                   ! <== ad_hoc_MM_tuning sets MegaMass to selected atoms

!------------------------------------------------------------------------------
! ENVIRONMENT parameters ...
!

  thermostat                = "Microcanonical"  ! <== Berendsen, Nose_Hoover, Microcanonical

  temperature               = 300.d0            ! <== Bath Temperature (K)
  pressure                  = 1.d0              ! <== Pressure

  thermal_relaxation_time   = 1.d-1             ! <== Temperature coupling term with the bath
                                                ! <== SMALL = STRONG ; use "= infty" to decouple

  pressure_relaxation_time  = infty             ! <== Pressure coupling term 
                                                ! <== SMALL = STRONG ; use "= infty" to decouple

  cutoff_radius             = 15.d0             ! <== Cut off radius (Angs.) for electrostatic and LJ interactions
  damping_Wolf              = 0.02d0            ! <== damping parameter (Angs.^-1) ; reasonable values: R_c*Wolf ~ ....
                                                ! <== Wolf's method damping parameter (length^{-1}) ; (J. Chem. Phys. 1999; 110(17):8254)
!------------------------------------------------------------------------------
! GENERAL INFO ...
!

  driver_MM              = "MM_Dynamics"       ! <== MM_Dynamics , MM_Optimize , NormalModes , Parametrize

  read_velocities        = T_                   ! <== reads the initial velocities : T_ , F_
  gmx_input_format       = T_                   ! <== reads FF parameters from gmx input files : T_ , F_  

  MM_log_step            =  100                 ! <== step for saving MM results & parameters
  MM_frame_step          =  20                  ! <== step for saving MM results & parameters

  Units_MM               = "eV"                 ! <== choose OUTPUT energy units: "eV" or "kj-mol" 
!--------------------------------------------------------------------
!           Genetic_Alg and CG OPTIMIZATION parameters
!

  OPT_driver     = "GACGRc"                     ! <== justCG , GACG , GACGAd , GACGRc

  Pop_Size       =  100    
  N_generations  =  15    
  Top_Selection  =  10                          ! <== top selection < Pop_Size
  Pop_range      =  0.15d0                      ! <== range of variation of parameters
  Mutation_rate  =  0.4           
  Mutate_Cross   =  F_                          ! <== false -> pure Genetic Algorithm 

  N_of_CGSteps   =  4                           ! <== number of CG steps for OPT_drivers: GACGAd or GACGRc
  nmd_window     = integer_interval( 0 , 00 )   ! <== only modes within the window entre cost evaluation ; 0 means no cutoff 

! =====================================================================================

end subroutine Define_MM_Environment

!
!
!
!================================
 subroutine allocate_species( N )
!================================
implicit none
integer , intent(in)    :: N

! local variables ...
integer :: i

allocate( species ( N ) )

do i = 1 , N
    species(i) % my_species     = 0
    species(i) % N_of_atoms     = 0
    species(i) % N_of_molecules = 0
    species(i) % cm(3)          = 0.0d0
    species(i) % mass           = 0.0d0
    species(i) % flex           = .false.
    species(i) % residue        = "XXX"
    species(i) % nr             = 0
    species(i) % Nbonds         = 0
    species(i) % Nangs          = 0
    species(i) % Ndiheds        = 0
    species(i) % Nharm          = 0
    species(i) % Nbonds14       = 0
    species(i) % NIntraLJ       = 0
end do

end subroutine allocate_species
!
!
!
end module MM_input

