! Subroutine for computing average properties over nuclear configurations ...

module Sampling_m

    use type_m
    use constants_m
    use Babel_m             , only : System_Characteristics ,       &
                                     Coords_from_Universe ,         &
                                     trj
    use Allocation_m        , only : Allocate_UnitCell ,            &
                                     DeAllocate_UnitCell ,          &
                                     DeAllocate_Structures 
    use Solvated_M          , only : Prepare_Solvated_System ,      &
                                     DeAllocate_TDOS ,              &
                                     DeAllocate_PDOS ,              &
                                     DeAllocate_SPEC 
    use QCModel_Huckel      , only : EigenSystem
    use FMO_m               , only : FMO_analysis
    use Structure_Builder   , only : Unit_Cell ,                    &
                                     Extended_Cell ,                &
                                     Generate_Structure ,           &
                                     Basis_Builder 
    use DOS_m
    use Oscillator_m        , only : Optical_Transitions
    use Multipole_Core      , only : Dipole_Matrix 
    use Schroedinger_m      , only : Huckel_dynamics ,              &
                                     DeAllocate_QDyn
    use Data_Output         , only : Dump_stuff

    use dipole_potential_m  , only : Solvent_Molecule_DP

    public :: Avrg_Confgs 

    private

contains
!
!
!=======================
 subroutine Avrg_Confgs
!=======================
 implicit none

! local variables ...
integer                         :: i , frame , nr , N_of_residues
real*8                          :: internal_sigma
character(3)                    :: residue
logical                         :: FMO_ , DIPOLE_
type(C_eigen)                   :: UNI , FMO
type(f_grid)                    :: TDOS , SPEC
type(f_grid)    , allocatable   :: PDOS(:) 
type(f_time)                    :: QDyn
type(universe)                  :: Solvated_System

! preprocessing stuff .....................................................

FMO_    = ( spectrum .OR. survival  )
DIPOLE_ = ( spectrum .OR. DP_Moment )

internal_sigma = sigma / ( float(size(trj))/float(frame_step) )      

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!..........................................................................

! Average over samples : Quantum Dynamics & All that Jazz ...

do frame = 1 , size(trj) , frame_step

    select case ( state_of_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

        case default

            Print*, " >>> Check your state_of_matter options <<< :" , state_of_matter
            stop

    end select

    CALL Generate_Structure( frame )

    CALL Basis_Builder( Extended_Cell , ExCell_basis )

    If( DP_field_ ) CALL Solvent_Molecule_DP( Extended_Cell )

    CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

    CALL Total_DOS( UNI%erg , TDOS , internal_sigma )                                             

    do nr = 1 , N_of_residues
!        CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr , internal_sigma )            
    end do

    If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

    If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )

    If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC , internal_sigma )

    If( survival ) CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO , QDyn )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
         DeAllocate             ( ExCell_basis )

    print*, frame

end do

! average over configurations ...
If( file_type == "trajectory" ) QDyn%dyn = QDyn%dyn / ( float(size(trj))/float(frame_step) )

CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

end subroutine Avrg_Confgs

end module Sampling_m
