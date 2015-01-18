Program qdynamo

use type_m
use constants_m
use parameters_m            , only : Define_Environment , driver , nuclear_matter              
use MM_input                , only : driver_MM
use Semi_Empirical_Parms    , only : read_EHT_parameters
use Structure_Builder       , only : Read_Structure
use qdynamics_m             , only : qdynamics
use Sampling_m              , only : Avrg_Confgs
use GA_driver_m             , only : GA_driver
use diagnostic_m            , only : diagnostic
use Chebyshev_driver_m      , only : Chebyshev_driver
use Eigen_driver_m          , only : Eigen_driver
use MMechanics_m            , only : MMechanics
use MD_read_m               , only : Build_MM_Environment
use good_vibrations_m       , only : Optimize_Structure , normal_modes , Optimize_Parameters_Driver

! local variables ...

! Initialize GPU if necessary 
call GPU_Init(0,1)

!========================================================
!                   DRIVER ROUTINE
!========================================================

CAll Define_Environment
   
If( driver /= "avrg_confgs") CALL system( "./env.sh" )

CALL read_EHT_parameters

CALL Read_Structure

If( driver == "MM_Dynamics" .OR. nuclear_matter == "MDynamics" ) CALL Build_MM_Environment

select case ( driver )

    case ( "q_dynamics" )
        CALL qdynamics

    case ( "slice_AO" , "slice_ElHl" , "slice_MO0" , "slice_MOt" )
        CALL Eigen_driver

    case ( "slice_Cheb" )
        CALL Chebyshev_driver

    case ( "avrg_confgs" )
        CALL Avrg_Confgs

    case ( "Genetic_Alg" )
        CALL GA_driver

    case ( "diagnostic" )
        CALL diagnostic

    case ( "MM_Dynamics" )

        select case ( driver_MM )

            case ( "MM_Dynamics" )
                CALL MMechanics        

            case ( "MM_Optimize" )
                CALL Optimize_Structure

            case ( "NormalModes" )
                CALL normal_modes

            case ( "Parametrize" )
                CALL Optimize_Parameters_Driver

            case default
                Print*, " >>> Check your driver options <<< :" , driver
                stop

        end select

            case default
                Print*, " >>> Check your driver options <<< :" , driver
                stop

        end select                

    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

include 'formats.h'

! Finalize GPU if necessary 
call GPU_Finalize

END
