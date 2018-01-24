module Chebyshev_driver_m

    use type_m
    use constants_m
    use parameters_m                , only : t_i , n_t , t_f , n_part ,     &
                                             frame_step , nuclear_matter ,  &
                                             DP_Field_ , Induced_ , QMMM ,  &
                                             GaussianCube , static ,        &
                                             GaussianCube_step ,            &
                                             hole_state , restart ,         &
                                             step_security
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
    use FMO_m                       , only : eh_tag
    use DP_main_m                   , only : Dipole_Matrix                 
    use Solvated_M                  , only : Prepare_Solvated_System 
    use TD_Dipole_m                 , only : wavepacket_DP                                        
    use DP_potential_m              , only : Molecular_DPs                                              
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Polarizability_m            , only : Build_Induced_DP
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format
    use Data_Output                 , only : Populations ,                  &
                                             Net_Charge
    use Backup_m                    , only : Security_Copy ,                &
                                             Restart_state ,                &
                                             Restart_Sys
    use MM_dynamics_m               , only : MolecularMechanics ,           &
                                             preprocess_MM , MoveToBoxCM 
    use Ehrenfest_Builder           , only : EhrenfestForce 
    use Chebyshev_m                 , only : Chebyshev ,                    &
                                             preprocess_Chebyshev
    use ElHl_Chebyshev_m            , only : ElHl_Chebyshev  ,              &
                                             preprocess_ElHl_Chebyshev

    public :: Chebyshev_driver

    private

    ! module variables ...
    Complex*16      , allocatable , dimension(:,:)  :: AO_bra , AO_ket , DUAL_ket , DUAL_bra , past_AO_bra , past_AO_ket
    real*8          , allocatable , dimension(:)    :: Net_Charge_MM
    integer                                         :: nn

contains
!
!
!
!========================================
 subroutine Chebyshev_driver( Qdyn , it )
!========================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(out)   :: it

! local variables ...
real*8          :: t , t_rate 
integer         :: frame , frame_init , frame_final , frame_restart
type(universe)  :: Solvated_System

it = 1
t  = t_i
nn = n_part

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
    CALL Restart_stuff( QDyn , t , it , frame_restart )  
else
    CALL Preprocess( QDyn , it )
end If

frame_init = merge( frame_restart+1 , frame_step+1 , restart )

do frame = frame_init , frame_final , frame_step

    If( (it >= n_t) .OR. (t >= t_f) ) exit    

    it = it + 1

    ! for using in Ehrenfest; Chebyshev also delivers data to Ehrenfest ...
    If( QMMM ) then
        past_AO_ket = AO_ket
        past_AO_bra = AO_bra
    end If    

    If( nn == 1) then
        CALL Chebyshev( Extended_Cell , ExCell_basis , AO_bra(:,1) , AO_ket(:,1) , Dual_bra(:,1) , Dual_ket(:,1) , QDyn , t , t_rate , it )
    else
        CALL ElHl_Chebyshev( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , t , t_rate , it )
    end If    

    ! calculate, for using in MM ...
    If( QMMM ) then
        Net_Charge_MM = Net_Charge
        CALL EhrenfestForce( Extended_Cell , ExCell_basis , past_AO_bra , past_AO_ket )
    end If

    If( GaussianCube .AND. mod(frame,GaussianCube_step) < frame_step ) CALL  Send_to_GaussianCube( frame , t )

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
            if( frame == frame_step+1 ) CALL preprocess_MM( Net_Charge = Net_Charge_MM )   
            ! MM precedes QM ; notice calling with frame -1 ...
            CALL MolecularMechanics( t_rate , frame - 1 , Net_Charge = Net_Charge_MM )   ! <== MM precedes QM ...

            ! IF QM_erg < 0 => turn off QMMM ; IF QM_erg > 0 => turn on QMMM ...
            QMMM = .NOT. (Unit_Cell% QM_erg < D_zero)

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
            stop

    end select

    CALL Generate_Structure ( frame )

    CALL Basis_Builder ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) CALL DP_stuff ( "DP_field" )

    If( Induced_ .OR. QMMM ) CALL DP_stuff ( "Induced_DP" )

    if( mod(frame,step_security) == 0 ) CALL Security_Copy( Dual_bra , Dual_ket , AO_bra , AO_ket , t , it , frame )

    print*, frame 

end do

call GPU_unpin( AO_bra )
call GPU_unpin( AO_ket )
deallocate( AO_bra , AO_ket , DUAL_bra , DUAL_ket )

include 'formats.h'

end subroutine Chebyshev_driver
!
!
!
!==================================
 subroutine Preprocess( QDyn , it )
!==================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(in)    :: it

! local variables
integer         :: hole_save 
logical         :: el_hl_
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

el_hl_ = any( Unit_Cell%Hl )  

CALL Generate_Structure ( 1 )

CALL Basis_Builder ( Extended_Cell , ExCell_basis )

If( Induced_ .OR. QMMM ) CALL Build_Induced_DP( basis = ExCell_basis , instance = "allocate" )

If( DP_field_ ) then

    hole_save  = hole_state
    hole_state = 0
    static     = .true. 

    ! DP potential in the static GS configuration ...
    CALL Molecular_DPs  ( Extended_Cell )

    hole_state = hole_save
    static     = .false.

end If

CALL Allocate_Brackets( size(ExCell_basis) , AO_bra , AO_ket , DUAL_bra , DUAL_ket , past_AO_bra , past_AO_ket )

If( nn == 1) then
    CALL preprocess_Chebyshev( Extended_Cell , ExCell_basis , AO_bra(:,1) , AO_ket(:,1) , Dual_bra(:,1) , Dual_ket(:,1) , QDyn , it )
else
    CALL preprocess_ElHl_Chebyshev( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )
end If

If( QMMM ) allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

If( GaussianCube ) CALL Send_to_GaussianCube  ( it , t_i )

If( Induced_ .OR. QMMM ) CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

!..........................................................................

include 'formats.h'

end subroutine Preprocess
!
!
! 
!=========================================
 subroutine Send_to_GaussianCube( it , t )
!=========================================
implicit none
integer     , intent(in)    :: it
real*8      , intent(in)    :: t

! local variables ...
integer :: n

! LOCAL representation for film STO production ...
do n = 1 , n_part
    CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
end do

!----------------------------------------------------------

end subroutine Send_to_GaussianCube
!
!
!
!
!===================================
 subroutine DP_stuff( instance )
!===================================
implicit none
character(*)  , intent(in)    :: instance

!local variables ...

!----------------------------------------------------------
!       LOCAL representation for DP calculation ...

select case( instance )

    case( "DP_matrix" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

    case( "DP_field" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        ! wavepacket component of the dipole vector ...
        CALL wavepacket_DP( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_ket )

        CALL Molecular_DPs( Extended_Cell )

    case( "Induced_DP" ) 

        If( .NOT. DP_field_ ) CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

end select

!----------------------------------------------------------

end subroutine DP_stuff
!
!
!
!=========================================================
 subroutine Restart_stuff( QDyn , t , it , frame_restart )
!=========================================================
implicit none
type(f_time)    , intent(out) :: QDyn
real*8          , intent(out) :: t
integer         , intent(out) :: it
integer         , intent(out) :: frame_restart

CALL DeAllocate_QDyn ( QDyn , flag="alloc" )

CALL Restart_State ( DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame_restart )

CALL Restart_Sys ( Extended_Cell , ExCell_basis , Unit_Cell , DUAL_ket , AO_bra , AO_ket , frame_restart )

If( QMMM ) then 

    allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

    CALL Build_Induced_DP( basis = ExCell_basis , instance = "allocate" )

    CALL DP_stuff ( "Induced_DP" )

end If

end subroutine Restart_stuff
!
!
!
end module Chebyshev_driver_m
