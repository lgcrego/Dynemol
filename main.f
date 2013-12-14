Program qdynamo

use type_m
use constants_m
use parameters_m            , only : driver , Define_Environment               
use Semi_Empirical_Parms    , only : read_EHT_parameters
use Structure_Builder       , only : Read_Structure
use qdynamics_m             , only : qdynamics
use Sampling_m              , only : Avrg_Confgs
use GA_driver_m             , only : GA_driver
use diagnostic_m            , only : diagnostic
use Chebyshev_driver_m      , only : Chebyshev_driver
use Eigen_driver_m          , only : Eigen_driver
use MMechanics_m            , only : MMechanics

! local variables ...
 

!========================================================
!                   DRIVER ROUTINE
!========================================================

CAll Define_Environment

If( driver /= "avrg_confgs") CALL system( "./env.sh" )

CALL read_EHT_parameters

CALL Read_Structure

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
        CALL MMechanics

    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

include 'formats.h'

END
