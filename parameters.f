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
  DRIVER         = "avrg_confgs"             ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO, FSSH] , MM_Dynamics
!			
  nuclear_matter = "extended_sys"            ! <== solvated_sys , extended_sys , MDynamics
!			
!			
  Survival       = T_                       
  DP_Moment      = F_                       
  QMMM           = F_                        ! <== Hellman-Feynman-Pulay ; HFP_Forces MUST be T_ for QMMM calcs 
  OPT_parms      = T_                        ! <== read OPT_basis parameters from "opt_eht_parameters.input.dat"
  ad_hoc         = T_                        ! <== ad hoc tuning of parameters

!----------------------------------------------------------------------------------------
!           MOLECULAR MECHANICS parameters are defined separately @ parameters_MM.f 
!----------------------------------------------------------------------------------------

!--------------------------------------------------------------------
!           READING FILE FORMAT
!
!  file_type    =  "structure"                 ! <== structure or trajectory
  file_type    =  "trajectory"                ! <== structure or trajectory
  file_format  =  "pdb"                       ! <== xyz , pdb or vasp
!--------------------------------------------------------------------
!           DIAGNOSTIC & DATA-ANALYSIS & VISUALIZATION flags
!
  HFP_Forces        = F_                      ! <== Hellman-Feynman-Pulay forces; MUST be T_ for QMMM calcs and F_ otherwise
                                              
  SPECTRUM          = F_                          
  Alpha_Tensor      = F_                      ! <== Embeded Finite Field Polarizability 

  GaussianCube      = F_                       
  GaussianCube_step = 5000000                 ! <== time step for saving Gaussian Cube files

  DensityMatrix     = F_                      ! <== generates data for postprocessing 
  AutoCorrelation   = F_             
  VDOS_             = F_
!--------------------------------------------------------------------
!           POTENTIALS
!
  EnvField_    =  F_                          ! <== Potential produced by Environment
  Environ_Type =  "Ch_MM"                     ! <== point charges: Ch_MM ; dipoles: { DP_QM , DP_MM } ...
  Environ_step =  10                          ! <== step for updating EnvField

  Coulomb_     =  F_                          ! <== use dipole potential for solvent molecules

  Induced_     =  F_                          ! <== use induced dipole potential 
!--------------------------------------------------------------------
!           SAMPLING parameters
!
  frame_step   =  5                           ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           SECURITY COPY
!
  restart       = F_                          ! <== TRUE for restarting dynamics
  step_security = 10000                       ! <== step for saving backup files
                                              ! <== default = 100 (QMMM) ; 1000 (MM) 
!--------------------------------------------------------------------
!           QDynamics parameters
!
  t_i  =  0.d0                              
  t_f  =  2.0d-1                              ! <== final time in PICOseconds
  n_t  =  1000                                ! <== number of time steps

  CT_dump_step = 1                            ! <== step for saving El&Hl survival charge density  

  n_part = 2                                  ! <== # of particles to be propagated: default is e=1 , e+h=2 

  hole_state     = 65                         ! <== GROUND STATE calcs     = 0 (ZERO)
                                              ! <== case STATIC & DP_calcs = hole state of special FMO
                                              ! <== case DYNAMIC           = intial MO for < HOLE > wavepacket in DONOR fragment

  electron_state = 66                         ! <== case STATIC & DP_calcs = excited state of special FMO
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

  Pop_Size       =  100    
  N_generations  =  1000
  Pop_range      =  0.3     ! <== range of variation of parameters [0:1]
  Mutation_rate  =  0.7     

  Adaptive_      =  T_       ! <== true  -> Adaptive GA method
  Mutate_Cross   =  T_       ! <== false -> pure Genetic Algorithm ; prefer false for fine tunning !

  CG_            =  F_       ! <== use conjugate gradient method after genetic algorithm
  Top_Selection  =  5        ! <== top selection to undergo CG_
  profiling      =  T_       ! <== for tuning the optimization parameters of the code

!--------------------------------------------------------------------
!  hereafter only CHECKLIST and  WARNINGS !!!

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" )
        
        dynamic = T_ 
        NetCharge = T_

    case( "avrg_confgs" , "Genetic_Alg" , "diagnostic" )

        dynamic = ( F_ .OR. Survival )
        NetCharge = F_

        If( Top_Selection > Pop_size ) stop ">> Top_Selection > Pop_size; execution aborted"

    case( "MM_Dynamics" )

        QMMM = F_
        dynamic = T_                                                                                                                                                                                            
        NetCharge = F_
        nuclear_matter = "MDynamics"
        
    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

static = .not. dynamic

! verbose is T_ only if ...
verbose = merge( T_ , F_ , (DRIVER /= "Genetic_Alg") .AND. (DRIVER /= "slice_AO") .AND. (DRIVER /= "slice_Cheb") .AND. (DRIVER /= "avrg_confgs") )



If ( DRIVER(1:5)=="slice" .AND. nuclear_matter=="extended_sys" .AND. file_type=="structure" ) then
    Print*," >>> halting: " 
    Print*,"     for fixed nuclei use DRIVER = q_dynamics; " 
    Print*,"     for slice use either file_type=trajectory or nuclear_matter=MDynamics <<<" 
    stop
End If    

If ( QMMM == T_ .AND. HFP_Forces == F_ ) then
    stop ">>> conflict between QMMM and HFP_Forces; execution halted, check parameters.f <<<"
elseif ( QMMM == F_ .AND. HFP_Forces == T_ .AND. driver /= "diagnostic" ) then
    CALL warning("MUST turn off HFP_Forces; execution halted, check parameters.f")
    stop 
end if

If ( driver == "slice_Cheb" .AND. electron_state*hole_state == 0 )  &
stop ">>> execution halted: driver=[slice_Cheb] only propagates ElHl pairs; for individual wvpckts use other drivers <<<"

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
