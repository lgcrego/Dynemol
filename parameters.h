 integer :: initial_state , HOMO_state
 logical :: CC , GaussianCube , Survival , POSCAR , SPECTRUM , DP_Moment
 type(real_interval) :: occupied , empty , DOS_range 
 type(integer_interval) :: holes , electrons

 parameter (& 
!--------------------------------------------------------------------
!           ACTIONS
!
            GaussianCube = .false. ,            &
			Survival     = .false. ,            &
            CC           = .false. ,            & 
            POSCAR       = .false. ,            & 
            SPECTRUM     = .true.  ,            & 
			DP_Moment    = .true.  ,            &
!--------------------------------------------------------------------
!           INITIAL  CONDITIONS
!
            t_i           =  0.d0 ,             &
            t_f           =  100.00 ,            & ! <== final time in PICOseconds
            n_t           =  500   ,             & ! <== number of time steps
            initial_state =  96  ,             & ! <== intial MO
			HOMO_state    =  95  ,             & ! <== HOMO of the molecule 
!--------------------------------------------------------------------
!           MOLECULAR DYNAMICS
!
            delta_t   =  8.d-4 ,                & ! <= time step of MD
!--------------------------------------------------------------------
!           STRUCTURAL  PARAMETERS
!
            nnx = 0    , nny = 0 ,              & ! <==  (nnx,nny) = (extended) REAL copies on each side
!
            T_x        =  20.44568 ,            & ! <== translation parameter
			T_y        =  22.69200 ,            & ! <== translation parameter
			T_z        =  50.00000 ,            & ! <== translation parameter

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
            sigma       =   0.160d0 ,                            &  !0.01d0

			DOS_range   =  real_interval( -17.d0 , -6.d0 ) ,     &

!--------------------------------------------------------------------
!           SPECTRUM  PARAMETERS
!
            occupied    =  real_interval( -17.30d0 , -11.01d0 ) , & 

            empty       =  real_interval( -11.00d0 , -5.00d0 )  , & 

!--------------------------------------------------------------------
!           Readfield  PARAMETERS
!
            holes     = integer_interval( 120 , ABOVE ) ,    & 

            electrons = integer_interval( 120 , ABOVE )  , 	  & 

!--------------------------------------------------------------------
!           2-Pi  COHERENT  CONTROL
!
            start_2P  = 17.000d0 ,              & ! <= start of the train 
            end_2P    = 250.00d0 ,              & ! <= end of the train 
            delta_2P  = 34.00d0                 & ! <= interval between 2Pi pulses
!--------------------------------------------------------------------
)

