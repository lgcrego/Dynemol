module EH_parms_module

use type_m ,       only :  real_interval , warning

integer                 :: nnx , nny , n_t , step_security , PBC(3)
integer                 :: n_part , electron_state , hole_state , frame_step , GaussianCube_step 
integer                 :: Pop_Size , N_generations , Top_Selection , CT_dump_step , Environ_step
real*8                  :: t_i , t_f , sigma
real*8                  :: Pop_range , Mutation_rate
type (real_interval)    :: occupied , empty , DOS_range
character (len=5)       :: Environ_type
character (len=6)       :: DK_of_mixing
character (len=4)       :: file_format
character (len=11)      :: DRIVER , file_type        
character (len=12)      :: nuclear_matter
character (len=8)       :: selection_by
logical                 :: DensityMatrix , AutoCorrelation , VDOS_ , Mutate_Cross , QMMM , LCMO , preview , Adaptive_
logical                 :: GaussianCube , Survival , SPECTRUM , DP_Moment , Alpha_Tensor , OPT_parms , ad_hoc , restart
logical                 :: verbose , static , EnvField_ , Coulomb_ , CG_ , profiling , Induced_ , NetCharge , HFP_Forces , Band_structure
logical                 :: resume , rnd_seed , ad_hoc_droplet

end module EH_parms_module
!
!
!
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
!
!
!
module VDOS_tuning                                                                                                                                            
integer :: Nsteps_per_sample , projection_rule
character(2) , allocatable ::my_symbols(:)
character(3) , allocatable ::my_residues(:)
character(1) , allocatable ::my_fragments(:)
character(2) , allocatable ::my_MMSymbols(:)
end module VDOS_tuning
!
!
!
module atomicmass
implicit none
    real*8 , parameter  , dimension(1:107)  :: Atomic_mass = (/                                 &
    1.00795d0,   4.00260d0,   6.94122d0,   9.01218d0,  10.81172d0,  12.01078d0,  14.00672d0,    &
   15.99943d0,  18.99840d0,  20.17976d0,  22.98970d0,  24.30506d0,  26.98153d0,  28.08553d0,    &
   30.97376d0,  32.06552d0,  35.45322d0,  39.94812d0,  39.09830d0,  40.07842d0,  44.95591d0,    &
   47.86710d0,  50.94151d0,  51.99616d0,  54.93805d0,  55.84500d0,  58.93320d0,  58.69344d0,    &
   63.54600d0,  65.38200d0,  69.72310d0,  72.64100d0,  74.92160d0,  78.96340d0,  79.90410d0,    &
   83.79822d0,  85.46783d0,  87.62120d0,  88.90585d0,  91.22422d0,  92.90638d0,  95.96220d0,    &
   98.94000d0, 101.07220d0, 102.90550d0, 106.42120d0, 107.86820d0, 112.41182d0, 114.81830d0,    &
  118.71000d0, 121.76000d0, 127.60320d0, 126.90477d0, 131.29362d0, 132.90540d0, 137.32770d0,    &
  138.90548d0, 140.11612d0, 140.90765d0, 144.24232d0, 146.91510d0, 150.36220d0, 151.96412d0,    &
  157.25320d0, 158.92535d0, 162.50012d0, 164.93032d0, 167.25932d0, 168.93421d0, 173.05452d0,    &
  174.96681d0, 178.49200d0, 180.94788d0, 183.84000d0, 186.20710d0, 190.23320d0, 192.21730d0,    &
  195.08490d0, 196.96654d0, 200.59000d0, 204.38332d0, 207.20000d0, 208.98040d0, 208.98400d0,    &
  209.98710d0, 222.01760d0, 223.01970d0, 226.02540d0, 227.02780d0, 232.03806d0, 231.03588d0,    &
  238.02891d0, 237.00000d0, 244.00000d0, 243.00000d0, 247.00000d0, 247.00000d0, 251.00000d0,    &
  252.00000d0, 257.00000d0, 256.00000d0, 254.00000d0, 257.00000d0,  13.01900d0,  14.02700d0,    &
   15.03500d0,  15.03500d0  /)
end module atomicmass

