MODULE parameters_m 

use type_m

integer                 :: nnx , nny , mmx , mmy , mmz , n_t 
integer                 :: initial_state , hole_state , frame_step , GaussianCube_step
real*8                  :: t_i , t_f , sigma
type (real_interval)    :: occupied , empty , DOS_range 
type (integer_interval) :: holes , electrons , rho_range
character (len=4)       :: file_format
character (len=11)      :: DRIVER , file_type 
character (len=12)      :: state_of_matter
logical                 :: GaussianCube , Survival , SPECTRUM , DP_Moment , OPT_basis , ad_hoc
logical                 :: verbose , static , DP_field_
logical , parameter     :: T_ = .true. , F_ = .false. 

integer , parameter     :: n_part = 1

contains
!
!
!=============================
 subroutine Define_Environment
!=============================
implicit none

!--------------------------------------------------------------------
! ACTION	flags
!
  DRIVER          = "slice_AO"                ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO , MO0 , MOt] 
!			
  state_of_matter = "extended_sys"            ! <== solvated_sys , extended_sys 
!			
  GaussianCube    = F_                       
  Survival        = T_                       
  SPECTRUM        = F_                          
  DP_Moment       = F_                       
  OPT_basis       = T_                        ! <== read OPT_basis parameters from "OPT_eht_parameters.input.dat"
  ad_hoc          = T_                        ! <== ad hoc tuning of parameters
!--------------------------------------------------------------------
!           READING FILE FORMAT
!
  file_type    =  "trajectory"                ! <= structure or trajectory
  file_format  =  "pdb"                       ! <= xyz , pdb or vasp
!--------------------------------------------------------------------
!           POTENTIALS
!
  DP_field_    =  T_                          ! <== use dipole potential for solvent molecules
!--------------------------------------------------------------------
!           SAMPLING parameters
!
  frame_step   =  1                           ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           QDynamics parameters
!
  t_i  =  0.d0                               
  t_f  =  5.d-1                               ! <== final time in PICOseconds
  n_t  =  1001                                ! <== number of time steps

  GaussianCube_step = 100                     ! <== time step for saving Gaussian Cube files

  hole_state    =  29                         ! <== 0 for ground state (GS) calculation
  initial_state =  30                         ! <== intial MO of DONOR fragment
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
  nnx = 0  ; nny = 0                          ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

  mmx = 1  ; mmy = 1   ; mmz = 1              ! <== PBC replicas : 1 = yes , 0 = no

!--------------------------------------------------------------------
!           DOS parameters
!
  sigma     =  0.080d0                                     !

  DOS_range = real_interval( -20.d0 , -1.d0 )            

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
  occupied  =  real_interval( -14.50d0 , -11.01d0 )       

  empty     =  real_interval( -11.00d0 , -6.00d0 )        

!--------------------------------------------------------------------

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" , "slice_MO0" , "slice_MOt" )
        
        static   = F_ 

    case( "avrg_confgs" , "Genetic_Alg" , "diagnostic" )

        static = T_ .OR. survival

    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

verbose = ( (DRIVER /= "Genetic_Alg")  .AND. (DRIVER /= "slice_AO" .OR. (.NOT. DP_Moment)) )   ! verbose OFF for those DRIVERS

end subroutine Define_Environment
!
!
!
end module parameters_m

