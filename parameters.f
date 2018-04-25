MODULE parameters_m 

use type_m

integer                 :: nnx , nny , n_t , step_security , PBC(3)
integer                 :: n_part , initial_state , hole_state , frame_step , GaussianCube_step , CH_and_DP_step
integer                 :: Pop_Size , N_generations , Top_Selection , file_size 
real*8                  :: t_i , t_f , sigma
real*8                  :: Pop_range , Mutation_rate  
type (real_interval)    :: occupied , empty , DOS_range 
type (integer_interval) :: holes , electrons , rho_range
character (len=4)       :: file_format
character (len=11)      :: DRIVER , file_type 
character (len=12)      :: nuclear_matter
character (len=7)       :: argument
logical                 :: DensityMatrix , AutoCorrelation , VDOS_ , Mutate_Cross , QMMM , LCMO , exist , preview
logical                 :: GaussianCube , Survival , SPECTRUM , DP_Moment , Alpha_Tensor , OPT_parms , ad_hoc , restart
logical                 :: verbose , static , DP_field_ , Coulomb_ , CG_ , profiling , Induced_ , NetCharge , EHM_Forces 
logical , parameter     :: T_ = .true. , F_ = .false. 

contains
!
!
!=============================
 subroutine Define_Environment
!=============================
implicit none

! local variables ...
logical :: dynamic

!--------------------------------------------------------------------
! ACTION	flags
!
  DRIVER         = "MM_Dynamics"             ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO, ElHl ] , MM_Dynamics

!			
  nuclear_matter = "MDynamics"               ! <== solvated_sys , extended_sys , MDynamics
!			
!			
  Survival       = F_                       
  DP_Moment      = F_                       
  QMMM           = F_
  OPT_parms      = F_                        ! <== read OPT_basis parameters from "OPT_eht_parameters.input.dat"
  ad_hoc         = T_                        ! <== ad hoc tuning of parameters

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
  
  SPECTRUM          = F_                          
  Alpha_Tensor      = F_                      ! <== Embeded Finite Field Polarizability 
  EHM_Forces        = F_                      ! <== for diagnostic only: Hellman-Feynman-Pulay forces for Ext. Huckel 

  GaussianCube      = F_                       
  GaussianCube_step = 500                     ! <== time step for saving Gaussian Cube files

  NetCharge         = F_                      ! <== pdb format charge Occupancy 
  CH_and_DP_step    = 1                       ! <== time step for saving charge and Induced DP values
                                              ! <== pdb format: charge --> Occupancy ; DP --> next to occupancy

  DensityMatrix     = F_                      ! <== generates data for postprocessing 
  AutoCorrelation   = F_             
  VDOS_             = F_
!--------------------------------------------------------------------
!           SECURITY COPY
!
  restart       = F_                          ! <== TRUE for restarting dynamics
  step_security = 1                           ! <== step for saving backup files
!--------------------------------------------------------------------
!           POTENTIALS
!
  DP_field_    =  F_                          ! <== use dipole potential for solvent molecules

  Coulomb_     =  F_                          ! <== use dipole potential for solvent molecules

  Induced_     =  F_                          ! <== use induced dipole potential 
!--------------------------------------------------------------------
!           SAMPLING parameters
!
  frame_step   =  1                           ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           QDynamics parameters
!
  t_i  =  0.d0                              
  t_f  =  1.0d-2                              ! <== final time in PICOseconds
  n_t  =  40                                  ! <== number of time steps

  n_part = 1                                  ! <== # of particles to be propagated: default is e=1 , e+h=2 

  hole_state    = 120                         ! <== GROUND STATE calcs     = 0 (ZERO)
                                              ! <== case STATIC & DP_calcs = hole state of special FMO
                                              ! <== case DYNAMIC           = intial MO for < HOLE >     wavepacket in DONOR fragment

  initial_state = 121                         ! <== case STATIC & DP_calcs = excited state of special FMO
                                              ! <== case DYNAMIC           = intial MO for < ELECTRON > wavepacket in DONOR fragment

  LCMO = F_                                   ! <== initial wavepackets as Linear Combination of Molecular Orbitals (LCMO)
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
  nnx = 0  ; nny = 0                          ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

  PBC = [ 0 , 0 , 0 ]                         ! <== PBC replicas : 1 = yes , 0 = no

!--------------------------------------------------------------------
!           DOS parameters
!
  sigma     =  0.080d0                                     !

  DOS_range = real_interval( -16.d0 , -4.d0 )            

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
  occupied  =  real_interval( -15.50d0 , -9.501d0 )       

  empty     =  real_interval( -9.500d0 , -4.00d0 )        

!--------------------------------------------------------------------
!           Genetic_Alg and CG OPTIMIZATION parameters
!

  Pop_Size       =  50  
  N_generations  =  10    
  Top_Selection  =  20                     ! <== top selection < Pop_Size
  Pop_range      =  0.05d0                 ! <== range of variation of parameters
  Mutation_rate  =  0.5           
  Mutate_Cross   =  F_                     ! <== false -> pure Genetic Algorithm ; prefer false for fine tunning !

  CG_            =  F_                     ! <== use conjugate gradient method after genetic algorithm
  profiling      =  T_                     ! <== for tuning the optimization parameters of the code

!--------------------------------------------------------------------
!  WARNINGS !!!

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" , "slice_ElHl" , "slice_M0" , "slice_MOt" )
        
        dynamic = T_ .OR. Survival 

    case( "avrg_confgs" , "Genetic_Alg" , "diagnostic" )

        dynamic = F_ .OR. Survival

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
verbose = (DRIVER /= "Genetic_Alg") .AND. (DRIVER /= "slice_AO") .AND. (DRIVER /= "slice_Cheb") 

If ( DRIVER(1:5)=="slice" .AND. nuclear_matter=="extended_sys" .AND. file_type=="structure" ) then
    Print*," >>> halting: " 
    Print*,"     for fixed nuclei use DRIVER = q_dynamics; " 
    Print*,"     for slice use either file_type=trajectory or nuclear_matter=MDynamics <<<" 
    stop
End If    

If ( nuclear_matter == "MDynamics" ) NetCharge = T_

If ( (frame_step /= 1) .AND. (file_type /= "trajectory") ) then
    Print*, " >>> halting: frame_step /= 1, only for avrg_confgs or time-slice dynamics <<<" 
    stop
End If    

!-------------------------------------------------------------
! get command line argument to preview input data, then stop  ...
CALL GET_COMMAND_ARGUMENT( 1 , argument )
preview = F_
if( COMMAND_ARGUMENT_COUNT() /= 0 ) then
    select case ( argument )

        case( "preview" )
        preview = .true.

    end select
end if
!--------------------------------------------------------------

include 'formats.h'

end subroutine Define_Environment
!
!
!
end module parameters_m
