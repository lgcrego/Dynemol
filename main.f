Program qdynamo

use type_m
use constants_m
use Semi_Empirical_Parms    , only : read_EHT_parameters
use Structure_Builder       , only : Read_Structure
use qdynamics_m             , only : qdynamics
use Sampling_m              , only : Solvated_M

! local variables ...
 

!========================================================
!                   DRIVER ROUTINE
!========================================================
 
CALL read_EHT_parameters

CALL Read_Structure

select case ( driver )

    case ( "q_dynamics" )
        CALL qdynamics

    case ( "solvated_M" )
        CALL Solvated_M

    case ( "AlgGenetic" )

end select

include 'formats.h'

END
