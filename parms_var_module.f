module EH_parms_module

use type_m ,       only :  real_interval , warning

integer                 :: nnx , nny , n_t , step_security , PBC(3)
integer                 :: n_part , electron_state , hole_state , frame_step , GaussianCube_step , CH_and_DP_step
integer                 :: Pop_Size , N_generations , Top_Selection , CT_dump_step , Environ_step
real*8                  :: t_i , t_f , sigma
real*8                  :: Pop_range , Mutation_rate
type (real_interval)    :: occupied , empty , DOS_range
character (len=5)       :: Environ_type
character (len=4)       :: file_format
character (len=11)      :: DRIVER , file_type        
character (len=12)      :: nuclear_matter
character (len=8)       :: selection_by
logical                 :: DensityMatrix , AutoCorrelation , VDOS_ , Mutate_Cross , QMMM , LCMO , preview , Adaptive_
logical                 :: GaussianCube , Survival , SPECTRUM , DP_Moment , Alpha_Tensor , OPT_parms , ad_hoc , restart
logical                 :: verbose , static , EnvField_ , Coulomb_ , CG_ , profiling , Induced_ , NetCharge , HFP_Forces , Band_structure
logical                 :: resume

end module EH_parms_module


module MM_parms_module

use type_m          , only : integer_interval
use MM_types        , only : MM_molecular , MM_system , MM_atomic
use EH_parms_module , only : Pop_size , N_generations , Top_Selection , Pop_range , Mutation_rate , Mutate_Cross

real*8                 :: temperature, pressure, cutoff_radius, thermal_relaxation_time, pressure_relaxation_time, damping_Wolf
integer                :: MM_log_step, MM_frame_step, spawn_step
logical                :: read_velocities, Selective_Dynamics, spawn
character (4)          :: MM_input_format
character (5)          :: Units_MM
character (len=14)     :: OPT_driver
character (len=11)     :: driver_MM 
character (len=14)     :: thermostat
type(integer_interval) :: nmd_window
type(MM_system)        :: MM
type(MM_molecular) , allocatable :: species(:) 

end module MM_parms_module
