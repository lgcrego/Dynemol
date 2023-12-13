module Chebyshev_driver_m

    use MPI
    use type_m
    use constants_m
    use MPI_definitions_m           , only : master , world , myid,         &
                                             ForceComm , ForceCrew ,        &
                                             ChebyCrew, KernelComm ,        &
                                             myid , world 
    use parameters_m                , only : t_i , n_t , t_f , n_part ,     &
                                             frame_step , nuclear_matter ,  &
                                             EnvField_ , Induced_ , QMMM ,  &
                                             GaussianCube , static ,        &
                                             GaussianCube_step ,            &
                                             hole_state , restart ,         &
                                             step_security, HFP_Forces ,    &
                                             preview , Environ_step
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj , MD_dt
    use Allocation_m                , only : DeAllocate_UnitCell ,          &
                                             DeAllocate_Structures ,        &
                                             Allocate_Brackets 
    use Structure_Builder           , only : Unit_Cell ,                    &
                                             Extended_Cell ,                &
                                             Generate_Structure ,           &
                                             Basis_Builder ,                &
                                             ExCell_basis
    use tuning_m                    , only : eh_tag
    use DP_main_m                   , only : Dipole_Matrix                 
    use Solvated_M                  , only : Prepare_Solvated_System 
    use TD_Dipole_m                 , only : wavepacket_DP                                        
    use Dielectric_Potential        , only : Environment_SetUp                                              
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Polarizability_m            , only : Build_Induced_DP
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format
    use Data_Output                 , only : Populations 
    use Backup_m                    , only : Security_Copy ,                &
                                             Restart_state ,                &
                                             Restart_Sys
    use MM_dynamics_m               , only : MolecularMechanics ,           &
                                             preprocess_MM , MoveToBoxCM 
    use DiabaticEhrenfest_Builder   , only : EhrenfestForce 
#ifdef USE_GPU
    use gpu_ElHl_Chebyshev_m        , only : ElHl_Chebyshev => gpu_ElHl_Chebyshev,    &
                                             preprocess_ElHl_Chebyshev => gpu_preprocess_ElHl_Chebyshev
#warning "Compiling Chebyshev for GPU. Subroutine calls redirected to GPU ones"
#else
    use ElHl_Chebyshev_m            , only : ElHl_Chebyshev  ,              &
                                             preprocess_ElHl_Chebyshev
#endif

    public :: Chebyshev_driver

    private

    ! module variables ...
    Complex*16      , allocatable , dimension(:,:)  :: AO_bra , AO_ket , DUAL_ket , DUAL_bra 
    real*8          , allocatable , dimension(:)    :: Net_Charge_MM
    real*8                                          :: t
    integer                                         :: N , it

contains
!
!
!
!==============================================
 subroutine Chebyshev_driver( Qdyn , final_it )
!==============================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(out)   :: final_it

! local variables ...
integer         :: mpi_D_R = mpi_double_precision , err
integer         :: frame , frame_init , frame_final , frame_restart
real*8          :: t_rate 
type(universe)  :: Solvated_System

it = 1
t  = t_i

!--------------------------------------------------------------------------------
! Quantum Dynamics & All that Jazz ...

! time is PICOseconds in EHT & seconds in MM ...
If( nuclear_matter == "MDynamics" ) then
    t_rate      = t_f / float(n_t)
    frame_final = n_t
else
    t_rate      = merge( t_f / float(n_t) , MD_dt * frame_step , MD_dt == epsilon(1.d0) )
    frame_final = size(trj)
end If

If( restart ) then
    CALL Restart_stuff( QDyn , frame_restart )  
else
    CALL Preprocess( QDyn )
end If

frame_init = merge( frame_restart+1 , frame_step+1 , restart )

do frame = frame_init , frame_final , frame_step

    If( (it >= n_t) .OR. (t >= t_f) ) exit    

    it = it + 1

    CALL ElHl_Chebyshev( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , t , t_rate , it )

    ! calculate, for using in MM ...
    If( QMMM ) then
        CALL EhrenfestForce( Extended_Cell , ExCell_basis )
    end If

    If( GaussianCube .AND. mod(frame,GaussianCube_step) < frame_step ) CALL  Send_to_GaussianCube( frame )

    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    ! build new UNI(t + t_rate) ...
    !============================================================================

    select case ( nuclear_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL DeAllocate_UnitCell ( Unit_Cell )
            CALL Coords_from_Universe( Unit_Cell , Solvated_System )

        case( "extended_sys" )

            CALL DeAllocate_UnitCell ( Unit_Cell )
            CALL Coords_from_Universe( Unit_Cell , trj(frame) )

        case( "MDynamics" )

            ! MM preprocess ...
            if( frame == frame_step+1 ) CALL preprocess_MM()   
            ! MM precedes QM ; notice calling with frame -1 ...
            CALL MolecularMechanics( t_rate , frame - 1 )   ! <== MM precedes QM ...

            ! IF QM_erg < 0 => turn off QMMM ; IF QM_erg > 0 => turn on QMMM ...
            QMMM = (.NOT. (Unit_Cell% QM_erg < D_zero)) .AND. (HFP_Forces == .true.)

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
            stop

    end select

    CALL Generate_Structure ( frame )

    ! export new coordinates for ForceCrew, only if QMMM = true, to avoid halting ...
    If( QMMM ) CALL MPI_BCAST( Extended_Cell%coord , Extended_Cell%atoms*3 , mpi_D_R , 0 , ForceComm, err )

    CALL Basis_Builder ( Extended_Cell , ExCell_basis )

    If( EnvField_ ) CALL DP_stuff ( "EnvField" )

    If( Induced_ ) CALL DP_stuff ( "Induced_DP" )

    if( mod(frame,step_security) == 0 ) CALL Security_Copy( Dual_bra , Dual_ket , AO_bra , AO_ket , t , it , frame )

    print*, frame 

end do

call GPU_unpin( AO_bra )
call GPU_unpin( AO_ket )
deallocate( AO_bra , AO_ket , DUAL_bra , DUAL_ket )

final_it = it

include 'formats.h'

end subroutine Chebyshev_driver
!
!
!
!=============================
 subroutine Preprocess( QDyn )
!=============================
implicit none
type(f_time)    , intent(out)   :: QDyn

! local variables
integer         :: hole_save , err
integer         :: mpi_D_R = mpi_double_precision
type(universe)  :: Solvated_System

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

select case ( nuclear_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , 1 )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(1) )

    case( "MDynamics" )

        CALL MoveToBoxCM

    case default

        Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
        stop

end select

CALL Generate_Structure ( 1 )

CALL Basis_Builder ( Extended_Cell , ExCell_basis )

If( Induced_ ) CALL Build_Induced_DP( basis = ExCell_basis , instance = "allocate" )

If( EnvField_ ) then

    hole_save  = hole_state
    hole_state = 0
    static     = .true. 

    ! DP potential in the static GS configuration ...
    CALL Environment_SetUp( Extended_Cell )

    hole_state = hole_save
    static     = .false.

    CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

end If

N = size(ExCell_basis)
CALL Allocate_Brackets( N , AO_bra , AO_ket , DUAL_bra , DUAL_ket )

CALL preprocess_ElHl_Chebyshev( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )

! done for ForceCrew ; ForceCrew dwells in EhrenfestForce ...
If( ForceCrew ) CALL EhrenfestForce( Extended_Cell , ExCell_basis , AO_bra , AO_ket )

If( QMMM ) allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

If( GaussianCube ) CALL Send_to_GaussianCube  ( it )

If( Induced_ ) CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

! export new coordinates for ForceCrew on stand-by ...
If( master ) CALL MPI_BCAST( Extended_Cell%coord , Extended_Cell%atoms*3 , mpi_D_R , 0 , ForceComm, err )
!..........................................................................

include 'formats.h'

end subroutine Preprocess
!
!
! 
!========================================
 subroutine Send_to_GaussianCube( frame )
!========================================
implicit none
integer , intent(in) :: frame

! local variables ...
integer :: n

! LOCAL representation for film STO production ...
do n = 1 , n_part
    if( eh_tag(n) == "XX" ) cycle
    CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , frame ,t , eh_tag(n) )
end do

!----------------------------------------------------------

end subroutine Send_to_GaussianCube
!
!
!
!
!===============================
 subroutine DP_stuff( instance )
!===============================
implicit none
character(*) , intent(in) :: instance

!local variables ...

!----------------------------------------------------------
!       LOCAL representation for DP calculation ...

select case( instance )

    case( "DP_matrix" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

    case( "EnvField" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        ! wavepacket component of the dipole vector ...
        ! decide what to do with this ############ 
        !CALL wavepacket_DP( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_ket )

        If( mod(it-1,Environ_step) == 0 ) CALL Environment_SetUp( Extended_Cell )

    case( "Induced_DP" ) 

        If( .NOT. EnvField_ ) CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

end select

!----------------------------------------------------------

end subroutine DP_stuff
!
!
!
!================================================
 subroutine Restart_stuff( QDyn , frame_restart )
!================================================
implicit none
type(f_time)    , intent(out) :: QDyn
integer         , intent(out) :: frame_restart

!local variables ...
integer :: err
integer :: MPI_D_R = mpi_double_precision

CALL DeAllocate_QDyn ( QDyn , flag="alloc" )

CALL Restart_State ( DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame_restart )

CALL Restart_Sys ( Extended_Cell , ExCell_basis , Unit_Cell , DUAL_ket , AO_bra , AO_ket , frame_restart )

CALL Preprocess_ElHl_Chebyshev( Extended_Cell , ExCell_basis , DUAL_ket , AO_bra , AO_ket )

! done for ForceCrew ; ForceCrew dwells in EhrenfestForce ...
If( ForceCrew ) CALL EhrenfestForce( Extended_Cell , ExCell_basis , AO_bra , AO_ket )

If( Induced_ ) then
      CALL Build_Induced_DP( basis = ExCell_basis , instance = "allocate" )
      CALL DP_stuff ( "Induced_DP" )
end If

allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

! ForceCrew is on stand-by for this ...
If( master ) CALL MPI_BCAST( Extended_Cell%coord , Extended_Cell%atoms*3 , mpi_D_R , 0 , ForceComm , err )

end subroutine Restart_stuff
!
!
!
end module Chebyshev_driver_m
