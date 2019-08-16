Program qdynamo

use MPI
use type_m
use constants_m
use MPI_definitions_m       , only : launch_MPI , master , world , myid
use parameters_m            , only : Define_Environment , driver , nuclear_matter              
use MM_input                , only : driver_MM
use Semi_Empirical_Parms    , only : read_EHT_parameters
use Structure_Builder       , only : Read_Structure
use qdynamics_m             , only : qdynamics
use Sampling_m              , only : Avrg_Confgs
use GA_driver_m             , only : GA_driver
use diagnostic_m            , only : diagnostic
use Chebyshev_driver_m      , only : Chebyshev_driver
use QMDynamicSlice_driver_m , only : QMDynamicSlice_driver
use MMechanics_m            , only : MMechanics
use MD_read_m               , only : Build_MM_Environment
use good_vibrations_m       , only : Optimize_Structure , normal_modes , Optimize_Parameters_Driver
use Ehrenfest_Builder       , only : EhrenfestForce

! local variables ...
integer :: err

CALL Define_Environment
   
! Initialize MPI if necessary ...
call launch_MPI

! Initialize GPU if necessary ...
CALL GPU_Init(myid,1)

!========================================================
!                   DRIVER ROUTINE
!========================================================

If( master .and. driver /= "avrg_confgs") CALL system( "./env.sh" )

CALL read_EHT_parameters

CALL Read_Structure

If( driver == "MM_Dynamics" .OR. nuclear_matter == "MDynamics" ) CALL Build_MM_Environment

select case ( driver )

    case ( "q_dynamics" )
        CALL qdynamics

    case ( "slice_AO" , "slice_Cheb" )
        CALL QMDynamicSlice_driver

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

include 'formats.h'

! Finalize MPI if necessary 
Call mpi_barrier( world , err )
call MPI_FINALIZE(err)

! Finalize GPU if necessary 
call GPU_Finalize

END
