 integer 					:: initial_state , HOMO_mol  
 logical 					:: GaussianCube , Survival , SPECTRUM , DP_Moment
 type (real_interval) 		:: occupied , empty , DOS_range 
 type (integer_interval) 	:: holes , electrons , rho_range
 character (len=4)			:: file_format
 character (len=10)			:: DRIVER , file_type

 parameter (& 
!--------------------------------------------------------------------
!           ACTIONS
!
			DRIVER		 = "q_dynamics",	 	& ! <== q_dynamics , solvated_M , Genetic_Ag
!			
            GaussianCube = .false. ,            &
			Survival     = .false. ,            &
            SPECTRUM     = .false. ,            & 
			DP_Moment    = .false. ,            &
!--------------------------------------------------------------------
!           READING FILE FORMAT
!
            file_type	 =  "structure"  ,      & ! <= structure or trajectory
            file_format  =  "pdb"  ,            & ! <= xyz , pdb or vasp

!--------------------------------------------------------------------
!           INITIAL  CONDITIONS
!
            t_i           =  0.d0 ,             &
            t_f           =  1.d1 ,             & ! <== final time in PICOseconds
            n_t           =  500   ,            & ! <== number of time steps
            initial_state =  99  ,              & ! <== intial MO
			HOMO_mol      =  96  ,              & ! <== HOMO of the molecule 
!--------------------------------------------------------------------
!           STRUCTURAL  PARAMETERS
!
            nnx = 0    , nny = 0 ,              & ! <==  (nnx,nny) = (extended) REAL copies on each side
!
!           Periodic Boundary Conditions 

            mmx = 1    , mmy = 1 ,              & ! <== PBC replicas : 1 = yes , 0 = no
            n_unit     =  (2*mmx+1)*(2*mmy+1) , & ! <== # unit cells repeated periodically 
!--------------------------------------------------------------------
!           SLATER  PARAMETERS
!
            n_part     =   1 ,                  & ! <== # of particals in Slater determinant
            nSL        =   1 ,                  & ! <== nSL = (n_part)! = size of Slater wave-function
!--------------------------------------------------------------------
!           DOS PARAMETERS
!
            sigma       =   0.080d0 ,                            &  !

			DOS_range   =  real_interval( -17.d0 , -6.d0 ) ,     &

!--------------------------------------------------------------------
!           SPECTRUM  PARAMETERS
!
            occupied    =  real_interval( -13.30d0 , -11.01d0 ) , & 

            empty       =  real_interval( -11.00d0 , -6.00d0 )  , & 

!--------------------------------------------------------------------
!           Readfield  PARAMETERS
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
            delta_2P  = 34.00d0                 & ! <= interval between 2Pi pulses
!--------------------------------------------------------------------
)

