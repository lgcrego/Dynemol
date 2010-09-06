! Subroutine for computing time evolution through time slices
module Chebyshev_driver_m

    use type_m
    use constants_m
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
    use Dipole_potential_m          , only : Solvent_Molecule_DP     
    use Data_Output                 , only : Dump_stuff 
    use Chebyshev_m                 , only : Chebyshev  ,                   &
                                             preprocess_Chebyshev

    public :: Chebyshev_driver

    private

contains
!
!
!===========================
 subroutine Chebyshev_driver
!===========================
implicit none

! local variables ...
integer                     :: it , frame 
real*8                      :: t 
real*8       , allocatable  :: QDyn_temp(:,:)
complex*16   , allocatable  :: Psi(:)
logical                     :: done = .false.
type(f_time)                :: QDyn
type(universe)              :: Solvated_System

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

it = 0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

do frame = 1 , size(trj) , frame_step

    if( (it >= n_t) .OR. (t >= t_f) ) exit

    it = it + 1

    select case ( state_of_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

        case default

            Print*, " >>> Check your state_of_matter options <<< :" , state_of_matter

    end select

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Solvent_Molecule_DP    ( Extended_Cell )

    If( .NOT. done ) then

        CALL preprocess_Chebyshev( Extended_Cell , ExCell_basis , Psi , QDyn , it )
        done = .true.

    else

        CALL Chebyshev( Extended_Cell , ExCell_basis , Psi , QDyn , t , it )

    end if

    CALL DeAllocate_UnitCell   ( Unit_Cell     )
    CALL DeAllocate_Structures ( Extended_Cell )
    DeAllocate                 ( ExCell_basis  )

    print*, frame

end do

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn )

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

end subroutine Chebyshev_driver
!
!
!
end module Chebyshev_driver_m
