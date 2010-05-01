 integer 					:: initial_state , HOMO_mol , frame_step 
 real*8						:: MD_dt
 logical 					:: GaussianCube , Survival , SPECTRUM , DP_Moment , OPT_basis , ad_hoc
 logical 					:: verbose , DP_field_
 type (real_interval) 		:: occupied , empty , DOS_range 
 type (integer_interval) 	:: holes , electrons , rho_range
 character (len=4)			:: file_format
 character (len=11)			:: DRIVER , file_type
 logical , parameter 		:: T_ = .true. , F_ = .false. 

 parameter (& 
!--------------------------------------------------------------------
!           ACTION	flags
!
			DRIVER		 = "eigen_slice"   ,	& ! <== q_dynamics , solvated_M , Genetic_Alg , solid_sys , diagnostic , eigen_slice
!			
            GaussianCube = F_ ,                 &
			Survival     = T_ ,                 &
            SPECTRUM     = F_ ,                 & 
			DP_Moment    = F_ ,                 &
			OPT_basis    = F_ ,                 & ! <== read OPT_basis parameters from "OPT_eht_parameters.input.dat"
			ad_hoc       = F_ ,                 & ! <== ad hoc tuning of parameters
!--------------------------------------------------------------------
!           READING FILE FORMAT
!
            file_type	 =  "trajectory" ,      & ! <= structure or trajectory
            file_format  =  "vasp" ,            & ! <= xyz , pdb or vasp
!--------------------------------------------------------------------
!           POTENTIALS
!
			DP_field_    =  F_  ,               & ! <== use dipole potential for solvent molecules
!--------------------------------------------------------------------
!           SAMPLING parameters
!
			frame_step    =  25   ,             & ! <== step for avrg_confgs ; frame_step =< size(trj)
!--------------------------------------------------------------------
!           QDynamics parameters
!
            t_i           =  0.d0 ,             &
            t_f           =  1.d-1 ,            & ! <== final time in PICOseconds
            n_t           =  500  ,             & ! <== number of time steps
			MD_dt		  =  8.d-4 ,			& ! <== time step of MD simulation in PICOseconds

            initial_state =  22  ,              & ! <== intial MO
			HOMO_mol      =  21  ,              & ! <== HOMO of the molecule 
!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
            nnx = 0    , nny = 0 ,              & ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

            mmx = 0    , mmy = 0 ,              & ! <== PBC replicas : 1 = yes , 0 = no
            n_unit     =  (2*mmx+1)*(2*mmy+1) , & ! <== # unit cells repeated periodically 
!--------------------------------------------------------------------
!           SLATER  parameters
!
            n_part     =   1 ,                  & ! <== # of particals in Slater determinant
            nSL        =   1 ,                  & ! <== nSL = (n_part)! = size of Slater wave-function
!--------------------------------------------------------------------
!           DOS parameters
!
            sigma       =   0.080d0 ,                            &  !

			DOS_range   =  real_interval( -17.d0 , -6.d0 ) ,     &

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

