module MM_input

use constants_m
use MM_parms_module 
use EH_parms_module , only : driver , Pop_size , N_generations , Top_Selection , Pop_range , Mutation_rate , Mutate_Cross

contains

!================================
 subroutine Define_MM_Environment
!================================
implicit none

!local parameters ...
logical, parameter :: T_ = .true. , F_ = .false. 

!------------------------------------------------------------------------------
! SYSTEM  INFO
!
  MM % N_of_molecules = 1                   ! <== total number of molecules
  MM % N_of_species   = 1                   ! <== total number of species

  CALL allocate_species( MM % N_of_species )

!------------------------------------------------------------------------------
! repeat the following information filling for all the different species ...
!
  species(1) % residue         = "AZP"      ! <== Residue label for species i ; character(len3)
  species(1) % N_of_molecules  = 1          ! <== Number of molecules of species i
  species(1) % N_of_atoms      = 35         ! <== Number of atoms comprosing a single molecule of species i
  species(1) % flex            = T_         ! <== Flexible : T_ , F_

  Selective_Dynamics = F_                   ! <== ad_hoc_MM_tuning sets MegaMass to selected atoms

!------------------------------------------------------------------------------
! ENVIRONMENT parameters ...
!

!  thermostat                = "Nose_Hoover"     ! <== Berendsen, Nose_Hoover, Microcanonical
!  thermostat                = "Berendsen"       ! <== Berendsen, Nose_Hoover, Microcanonical
  thermostat                = "Microcanonical"  ! <== Berendsen, Nose_Hoover, Microcanonical

  temperature               = 200.d0            ! <== Bath Temperature (K)
  pressure                  = 1.d0              ! <== Pressure

  thermal_relaxation_time   = 5.d-1             ! <== Temperature coupling term with the bath
                                                ! <== SMALL = STRONG ; use "= infty" to decouple

  pressure_relaxation_time  = infty             ! <== Pressure coupling term 
                                                ! <== SMALL = STRONG ; use "= infty" to decouple

  cutoff_radius             = 15.d0             ! <== Cut off radius (Angs.) for electrostatic and LJ interactions
  damping_Wolf              = 0.022d0           ! <== damping parameter (Angs.^-1) ; reasonable values: R_c*Wolf ~ ....
                                                ! <== Wolf's method damping parameter (length^{-1}) ; (J. Chem. Phys. 1999; 110(17):8254)
!------------------------------------------------------------------------------
! GENERAL INFO ...
!
  driver_MM              = "MM_Dynamics"       ! <== MM_Dynamics , MM_Optimize , NormalModes , Parametrize

  read_velocities        = T_                   ! <== reads the initial velocities : T_ , F_

  MM_input_format        = "GMX"                ! <== GMX, NAMD, GAFF

  MM_log_step            =  50                  ! <== step for saving MM results & parameters

  MM_frame_step          =  50                 ! <== step for saving MM results & parameters

  Units_MM               = "eV"                 ! <== choose OUTPUT energy units: "eV" or "kj-mol" 
!--------------------------------------------------------------------
!           Genetic_Alg and CG OPTIMIZATION parameters
!
  Pop_Size       =  50       
  N_generations  =  30    
  Top_Selection  =  5                           ! <== top selection < Pop_Size ; Top selection will optimized by CG
  Pop_range      =  0.22                        ! <== range of variation of parameters
  Mutation_rate  =  0.1           
  Mutate_Cross   =  F_                          ! <== false -> pure Genetic Algorithm 

  nmd_window     = integer_interval( 0 , 1500 )    ! <== only modes within the window are considered for cost evaluation 
                                                   ! (0,0) means no cutoff 

  OPT_driver     = "use_weights"                ! "just_Erg"       ; just optmize with respect to Energy
                                                ! "use_weights"    ; weight = weight(:)
                                                ! "use_no_weights" : weight = D_one
                                                ! "use_overweight" : weight = weight + |W - W_REF|/W_REF
                                                ! "LineUpCost"     : weight = weight + extra cost assigned for normal_modes that are out of place 
!=====================================================================================
If( .NOT.  any(["just_Erg","use_weights","use_no_weights","use_overweight","LineUpCost"] == OPT_driver) ) Stop ">>> check your OPT_driver options <<<"

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

