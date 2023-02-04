Program qdynamo

use MPI
use type_m
use constants_m
use setup_checklist      
use card_reading            , only : ReadInputCard
use MPI_definitions_m       , only : launch_MPI , setup_MPI_labor_force , master , myid
use parameters_m            , only : Define_Environment , driver , nuclear_matter , restart
use MM_input                , only : Define_MM_Environment , driver_MM
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

! local variables ...
integer :: err
logical :: go_without_card


! Initialize MPI if necessary ...
call launch_MPI

! Initialize GPU if necessary ...
CALL GPU_Init(myid,1)

!========================================================
!               THE TRUTH IS OUT THERE
!========================================================

CALL get_environment_vars

If( master ) then 
    If( .not. restart ) CALL system( dynemoldir//"env.sh" )
end If

inquire( file=dynemolworkdir//"makefile" , EXIST = go_without_card )

if( go_without_card ) then
     CAll Define_Environment
     CALL Define_MM_Environment
else
     call ReadInputCard
    If( master ) then
        ! cloning the psf files into trunk directories ...
        call system("cp card.inpt log.trunk/.")
        CALL system("echo dyn.trunk/ dos.trunk/ opt.trunk/ ancillary.trunk/ | xargs -n 1 cp log.trunk/card.inpt ")
        end if
end if

CALL setup_MPI_labor_force

If( master ) CALL checklist

CALL read_EHT_parameters

CALL Read_Structure

If( driver == "MM_Dynamics" .OR. nuclear_matter == "MDynamics" ) CALL Build_MM_Environment

If( master ) CALL dump_driver_parameters_and_tuning

select case ( driver )

    case ( "q_dynamics" )
        CALL qdynamics

    case ( "slice_AO" , "slice_FSSH" , "slice_Cheb" , "slice_CSDM" )
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

if( master ) then
    CALL system("echo dyn.trunk/ dos.trunk/ opt.trunk/ ancillary.trunk/ | xargs -n 1 cp log.trunk/driver_parms_and_tuning.log ")
    end if

include 'formats.h'

! Finalize MPI 
call MPI_FINALIZE(err)

! Finalize GPU 
call GPU_Finalize

END
