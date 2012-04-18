MODULE parameters_m 

use type_m

integer                 :: nnx , nny , mmx , mmy , mmz , n_t , step_security
integer                 :: n_part , initial_state , hole_state , frame_step , GaussianCube_step
real*8                  :: t_i , t_f , sigma
type (real_interval)    :: occupied , empty , DOS_range 
type (integer_interval) :: holes , electrons , rho_range
character (len=4)       :: file_format
character (len=11)      :: DRIVER , file_type 
character (len=12)      :: state_of_matter
logical                 :: GaussianCube , Survival , SPECTRUM , DP_Moment , OPT_basis , ad_hoc , restart
logical                 :: verbose , static , DP_field_ , Coulomb_
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
  DRIVER          = "slice_AO"                ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO, ElHl ] 
!			
  state_of_matter = "extended_sys"            ! <== solvated_sys , extended_sys 
!			
  GaussianCube    = F_                       
  Survival        = T_                       
  SPECTRUM        = F_                          
  DP_Moment       = T_                       
  OPT_basis       = T_                        ! <== read OPT_basis parameters from "OPT_eht_parameters.input.dat"
  ad_hoc          = T_                        ! <== ad hoc tuning of parameters

  GaussianCube_step = 100                     ! <== time step for saving Gaussian Cube files

!--------------------------------------------------------------------
!           READING FILE FORMAT
!
  file_type    =  "trajectory"                ! <= structure or trajectory
  file_format  =  "pdb"                       ! <= xyz , pdb or vasp
!--------------------------------------------------------------------
!           SECURITY COPY
!
  restart       = F_                          ! TRUE for restarting dynamics
  step_security = 100                         ! not yet implemented
!--------------------------------------------------------------------
!           POTENTIALS
!
  DP_field_    =  T_                          ! <== use dipole potential for solvent molecules

  Coulomb_     =  F_                          ! <== use dipole potential for solvent molecules
!--------------------------------------------------------------------
!           SAMPLING parameters
!
  frame_step   =  1                           ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           QDynamics parameters
!
  t_i  =  0.d0                               
  t_f  =  5.0d0                                 ! <== final time in PICOseconds
  n_t  =  1000                                  ! <== number of time steps

  n_part = 2                                    ! <== # of particles to be propagated: default is e=1 , e+h=2 

  hole_state    =  5                            ! <== GROUND STATE calcs     = 0 (ZERO)
                                                ! <== case STATIC & DP_calcs = hole state of special FMO
                                                ! <== case DYNAMIC           = intial MO for < HOLE >     wavepacket in DONOR fragment

  initial_state =  30                           ! <== case STATIC & DP_calcs = excited state of special FMO
                                                ! <== case DYNAMIC           = intial MO for < ELECTRON > wavepacket in DONOR fragment
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
  nnx = 0  ; nny = 0                          ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

  mmx = 0  ; mmy = 0   ; mmz = 1              ! <== PBC replicas : 1 = yes , 0 = no

!--------------------------------------------------------------------
!           DOS parameters
!
  sigma     =  0.080d0                                     !

  DOS_range = real_interval( -20.d0 , 20.d0 )            

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
  occupied  =  real_interval( -15.50d0 , -9.501d0 )       

  empty     =  real_interval( -9.500d0 , -4.00d0 )        

!--------------------------------------------------------------------

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" , "slice_ElHl" , "slice_MO0" , "slice_MOt" , "Coulomb" )
        
        dynamic = T_ .OR. Survival 

    case( "avrg_confgs" , "Genetic_Alg" , "diagnostic" )

        dynamic = F_ .OR. Survival

    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

static = .not. dynamic

! verbose is T_ only if ...
verbose = .true. !(DRIVER /= "Genetic_Alg") .AND. (DRIVER /= "slice_AO") 

end subroutine Define_Environment
!
!
!
end module parameters_m

