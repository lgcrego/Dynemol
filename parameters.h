 integer 					:: initial_state , HOMO_mol , frame_step , GaussianCube_step
 logical 					:: GaussianCube , Survival , SPECTRUM , DP_Moment , OPT_basis , ad_hoc
 logical 					:: verbose , DP_field_
 type (real_interval) 		:: occupied , empty , DOS_range 
 type (integer_interval) 	:: holes , electrons , rho_range
 character (len=4)			:: file_format
 character (len=11)			:: DRIVER , file_type 
 character (len=12)			:: state_of_matter
 logical , parameter 		:: T_ = .true. , F_ = .false. 

 parameter (& 
!--------------------------------------------------------------------
!           ACTION	flags
!
            DRIVER		 	= "avrg_confgs"  ,	   & ! <== q_dynamics , avrg_confgs , Genetic_Alg , diagnostic , slice_[Cheb, AO , MO0 , MOt] 
!			
            state_of_matter = "solvated_sys" ,     & ! <== solvated_sys , extended_sys 
!			
            GaussianCube 	= F_ ,                 &
			Survival     	= F_ ,                 &
            SPECTRUM     	= F_ ,                 & 
			DP_Moment    	= T_ ,                 &
			OPT_basis    	= T_ ,                 & ! <== read OPT_basis parameters from "OPT_eht_parameters.input.dat"
			ad_hoc       	= T_ ,                 & ! <== ad hoc tuning of parameters
!--------------------------------------------------------------------
!           READING FILE FORMAT
!
            file_type	 =  "trajectory",       & ! <= structure or trajectory
            file_format  =  "pdb"  ,            & ! <= xyz , pdb or vasp
!--------------------------------------------------------------------
!           POTENTIALS
!
			DP_field_    =  T_  ,               & ! <== use dipole potential for solvent molecules
!--------------------------------------------------------------------
!           SAMPLING parameters
!
			frame_step    =  1000 ,             & ! <== step for avrg_confgs and time-slice dynamics ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           QDynamics parameters
!
            t_i           =  0.d0 ,             &
            t_f           =  1.0d0 ,            & ! <== final time in PICOseconds
            n_t           =  1001  ,            & ! <== number of time steps

			GaussianCube_step = 100,   			& ! <== time step for saving Gaussian Cube files

            initial_state =  30  ,              & ! <== intial MO
			HOMO_mol      =  29  ,              & ! <== HOMO of the molecule 
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
            nnx = 0    , nny = 0 ,              		& ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

            mmx = 0    , mmy = 0	, mmz = 0   ,		& ! <== PBC replicas : 1 = yes , 0 = no

!--------------------------------------------------------------------
!           SLATER  parameters
!
            n_part     =   1 ,                  & ! <== # of particals in Slater determinant
            nSL        =   1 ,                  & ! <== nSL = (n_part)! = size of Slater wave-function
!--------------------------------------------------------------------
!           DOS parameters
!
            sigma       =   0.080d0 ,                            &  !

			DOS_range   =  real_interval( -20.d0 , -1.d0 ) ,     &

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
            occupied    =  real_interval( -14.50d0 , -11.01d0 ) , & 

            empty       =  real_interval( -11.00d0 , -6.00d0 )  , & 

!--------------------------------------------------------------------
!           Readfield  parameters
!

			hole_state	=	864									,	&

            holes     	= 	integer_interval( 864 , ABOVE ) 	,   & 

            electrons 	= 	integer_interval( 864 , ABOVE )  	,   & 

			rho_range 	= 	integer_interval( 864 , 1300  ) 	,   & 

!--------------------------------------------------------------------
!           2-Pi  COHERENT  CONTROL
!
            start_2P  = 17.000d0 ,              & ! <= start of the train 
            end_2P    = 250.00d0 ,              & ! <= end of the train 
            delta_2P  = 34.00d0  ,              & ! <= interval between 2Pi pulses
!--------------------------------------------------------------------
!			Environment
!
			verbose = ( DRIVER /= "Genetic_Alg" ) 					&
)

