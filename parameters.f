MODULE parameters_m 

use EH_parms_module

contains
!
!
!=============================
 subroutine Define_Environment
!=============================
implicit none

! local variables ...
character (len=7) :: argument
logical           :: dynamic

! local parameter ...
logical, parameter :: T_ = .true. , F_ = .false. 

!--------------------------------------------------------------------
! ACTION	flags
!
  DRIVER         = "MM_Dynamics"              ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO, FSSH, CSDM] , MM_Dynamics
!			
  nuclear_matter = "MDynamics"               ! <== solvated_sys , extended_sys , MDynamics
!			
!			
  Survival       = T_                       
  DP_Moment      = F_                       
  QMMM           = F_
  OPT_parms      = F_                        ! <== read OPT_basis parameters from "opt_eht_parms.input"
  ad_hoc         = T_                        ! <== ad hoc tuning of parameters
  Band_structure = F_
!----------------------------------------------------------------------------------------
!           MOLECULAR MECHANICS parameters are defined separately @ parameters_MM.f 
!----------------------------------------------------------------------------------------

!--------------------------------------------------------------------
!           READING FILE FORMAT
!
  file_type    =  "structure"                 ! <== structure or trajectory
  file_format  =  "pdb"                       ! <== xyz , pdb or vasp
!--------------------------------------------------------------------
!           DIAGNOSTIC & DATA-ANALYSIS & VISUALIZATION flags
!
  HFP_Forces        = F_                      ! <== Hellman-Feynman-Pulay forces; MUST be T_ for QMMM calcs and F_ otherwise
                                              
  SPECTRUM          = F_                          
  Alpha_Tensor      = F_                      ! <== Embeded Finite Field Polarizability 

  GaussianCube      = F_                       
  GaussianCube_step = 5000000                 ! <== time step for saving Gaussian Cube files

  NetCharge         = F_                      ! <== pdb format charge Occupancy 
  CH_and_DP_step    = 1000000                 ! <== time step for saving charge and Induced DP values
                                              ! <== pdb format: charge --> Occupancy ; DP --> next to occupancy

  DensityMatrix     = F_                      ! <== generates data for postprocessing 
  AutoCorrelation   = F_             
  VDOS_             = F_
!--------------------------------------------------------------------
!           POTENTIALS
!
  EnvField_    =  F_                          ! <== Potential produced by Environment
  Environ_Type =  "Ch_MM"                     ! <== point charges: Ch_MM ; dipoles: { DP_QM , DP_MM } ...
  Environ_step =  5                           ! <== step for updating EnvField

  Coulomb_     =  F_                          ! <== use dipole potential for solvent molecules

  Induced_     =  F_                          ! <== use induced dipole potential 
!--------------------------------------------------------------------
!           SAMPLING parameters
!
  frame_step   =  1                           ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           SECURITY COPY
!
  restart       = F_                          ! <== TRUE for restarting dynamics
  step_security = 1000                       ! <== step for saving backup files
                                              ! <== default = 100 (QMMM) ; 1000 (MM) 
!--------------------------------------------------------------------
!           QDynamics parameters
!
  t_i  =  0.d0                              
  t_f  =  5.d1                               ! <== final time in PICOseconds
  n_t  =  50000                              ! <== number of time steps

  n_part = 2                                  ! <== # of particles to be propagated: default is e=1 , e+h=2 

  hole_state    = 50                          ! <== GROUND STATE calcs     = 0 (ZERO)
                                              ! <== case STATIC & DP_calcs = hole state of special FMO
                                              ! <== case DYNAMIC           = intial MO for < HOLE > wavepacket in DONOR fragment

  electron_state = 52                          ! <== case STATIC & DP_calcs = excited state of special FMO
                                              ! <== case DYNAMIC           = intial MO for < ELECTRON > wavepacket in DONOR fragment

  LCMO = F_                                   ! <== initial wavepackets as Linear Combination of Molecular Orbitals (LCMO)
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
  nnx = 0  ; nny = 0                          ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

  PBC = [ 1 , 1 , 0 ]                         ! <== PBC replicas : 1 = yes , 0 = no

!--------------------------------------------------------------------
!           DOS parameters
!
  sigma     =  0.040d0                                     !

  DOS_range = real_interval( -15.d0 , 0.d0 )            

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
  occupied  =  real_interval( -15.50d0 , -9.501d0 )       

  empty     =  real_interval( -9.500d0 , -4.00d0 )        

!--------------------------------------------------------------------
!           Genetic_Alg and CG OPTIMIZATION parameters
!

  Pop_Size       =  200  
  N_generations  =  200    
  Pop_range      =  0.92          ! <== range of variation of parameters [0:1]
  selection_by   =  'roullete'     ! option={roullete,ranking,sorting}
  Mutation_rate  =  0.9     

  Adaptive_      =  T_            ! <== true  -> Adaptive GA method
  Mutate_Cross   =  T_            ! <== false -> pure Genetic Algorithm ; prefer false for fine tunning !

  CG_            =  F_            ! <== use conjugate gradient method after genetic algorithm
  Top_Selection  =  5             ! <== top selection to undergo CG_
  profiling      =  T_            ! <== for tuning the optimization parameters of the code

!--------------------------------------------------------------------
!  hereafter only CHECKLIST and  WARNINGS !!!

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" , "slice_FSSH" , "slice_CSDM" )
        
        dynamic = T_ 

    case( "avrg_confgs" , "Genetic_Alg" , "diagnostic" )

        dynamic = ( F_ .OR. Survival )

        If( Top_Selection > Pop_size ) stop ">> Top_Selection > Pop_size; execution aborted"

    case( "MM_Dynamics" )

        QMMM = F_
        dynamic = F_
        
    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

static = .not. dynamic

! verbose is T_ only if ...
verbose = (DRIVER /= "Genetic_Alg") .AND. (DRIVER /= "slice_AO") .AND. (DRIVER /= "slice_Cheb") .AND. (DRIVER /= "slice_FSSH") .AND. (DRIVER /= "slice_CSDM")

If ( DRIVER(1:5)=="slice" .AND. nuclear_matter=="extended_sys" .AND. file_type=="structure" ) then
    Print*," >>> halting: " 
    Print*,"     for fixed nuclei use DRIVER = q_dynamics; " 
    Print*,"     for slice use either file_type=trajectory or nuclear_matter=MDynamics <<<" 
    stop
End If    

If ( QMMM == T_ .AND. HFP_Forces == F_ ) then
    stop ">>> conflict between QMMM and HFP_Forces; execution halted, check parameters.f <<<"
elseif ( QMMM == F_ .AND. HFP_Forces == T_ .AND. driver /= "diagnostic" ) then
    CALL system("sed '11i>>> MUST turn off HFP_Forces; execution halted, check parameters.f <<<' .warning.signal |cat")
    stop 
end if

If ( nuclear_matter == "MDynamics" ) NetCharge = T_

!-------------------------------------------------------------
! get command line argument to preview input data, then stop  ...
CALL GET_COMMAND_ARGUMENT( 1 , argument )
preview = F_
if( COMMAND_ARGUMENT_COUNT() /= 0 ) then
    select case ( argument )

        case( "preview" )
        preview = .true.

        case( "resume" )
        resume = .true.

    end select
end if
!--------------------------------------------------------------

include 'formats.h'

end subroutine Define_Environment
!
!
!
end module parameters_m

