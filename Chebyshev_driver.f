! Subroutine for computing time evolution through Chebyshev expansion ...
! Propagation of ONE WAVE PACKET ONLY : (n_part = 1) ...
module Chebyshev_driver_m

    use type_m
    use constants_m
    use parameters_m                , only : t_f , n_t , nuclear_matter ,   &
                                             frame_step , DP_Field_ ,       &
                                             n_part
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj
    use Allocation_m                , only : Allocate_UnitCell ,            &
                                             DeAllocate_UnitCell ,          &
                                             DeAllocate_Structures ,        &
                                             Allocate_Brackets
    use Structure_Builder           , only : Unit_Cell ,                    &
                                             Extended_Cell ,                &
                                             Generate_Structure ,           &
                                             Basis_Builder ,                &
                                             ExCell_basis
    use Solvated_M                  , only : Prepare_Solvated_System                                              
    use Schroedinger_m              , only : DeAllocate_QDyn
    use DP_potential_m              , only : Molecular_DPs     
    use Data_Output                 , only : Dump_stuff 
    use ElHl_Chebyshev_m            , only : ElHl_Chebyshev  ,              &
                                             preprocess_ElHl_Chebyshev
    use Chebyshev_m                 , only : Chebyshev  ,                   &
                                             preprocess_Chebyshev

    public :: Chebyshev_driver

    private

contains
!
!
!
!===========================
 subroutine Chebyshev_driver
!===========================

select case( n_part )

    case( 1 )
        CALL El_Chebyshev_driver

    case( 2 )
        CALL ElHl_Chebyshev_driver

end select

end subroutine Chebyshev_driver
!
!
!
!================================
 subroutine ElHl_Chebyshev_driver
!================================
implicit none

! local variables ...
integer                     :: it , frame 
real*8                      :: t 
real*8       , allocatable  :: QDyn_temp(:,:,:)
complex*16   , allocatable  :: Psi_bra(:,:) , Psi_ket(:,:)
type(f_time)                :: QDyn
type(universe)              :: Solvated_System

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

it = 1

select case ( nuclear_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , 1 )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System , 1 )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(1) , 1 )

    case default

        Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter

end select

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) &
CALL Molecular_DPs          ( Extended_Cell )

CALL preprocess_ElHl_Chebyshev( Extended_Cell , ExCell_basis , Psi_bra , Psi_ket , QDyn , it=1 )

do frame = 2 , size(trj) , frame_step

    if( (it >= n_t) .OR. (t >= t_f) ) exit

    it = it + 1

    CALL ElHl_Chebyshev( Extended_Cell , ExCell_basis , Psi_bra , Psi_ket , QDyn , t , it )

    CALL DeAllocate_UnitCell   ( Unit_Cell     )
    CALL DeAllocate_Structures ( Extended_Cell )
    DeAllocate                 ( ExCell_basis  )

    select case ( nuclear_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter

    end select

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Molecular_DPs          ( Extended_Cell )

    print*, frame

end do

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 , n_part ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 , : ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn )   

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

deallocate( Psi_bra , Psi_ket )

end subroutine ElHl_Chebyshev_driver
!
!
!
!==============================
 subroutine El_Chebyshev_driver
!==============================
implicit none

! local variables ...
integer                     :: it , frame 
real*8                      :: t 
real*8       , allocatable  :: QDyn_temp(:,:,:)
complex*16   , allocatable  :: Psi_bra(:) , Psi_ket(:)
type(f_time)                :: QDyn
type(universe)              :: Solvated_System

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

it = 1

select case ( nuclear_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , 1 )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System , 1 )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(1) , 1 )

    case default

        Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter

end select

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) &
CALL Molecular_DPs          ( Extended_Cell )

CALL preprocess_Chebyshev( Extended_Cell , ExCell_basis , Psi_bra , Psi_ket , QDyn , it )

do frame = 2 , size(trj) , frame_step

    if( (it >= n_t) .OR. (t >= t_f) ) exit

    it = it + 1

    CALL Chebyshev( Extended_Cell , ExCell_basis , Psi_bra , Psi_ket , QDyn , t , it )

    CALL DeAllocate_UnitCell   ( Unit_Cell     )
    CALL DeAllocate_Structures ( Extended_Cell )
    DeAllocate                 ( ExCell_basis  )

    select case ( nuclear_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter

    end select

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Molecular_DPs          ( Extended_Cell )

    print*, frame

end do

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 , n_part ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 , : ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn )   

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

deallocate( Psi_bra , Psi_ket )

end subroutine El_Chebyshev_driver
!
!
!
end module Chebyshev_driver_m
