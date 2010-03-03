module qdynamics_m

 use type_m
 use constants_m
 use Solvated_M             , only : DeAllocate_TDOS ,      &
                                     DeAllocate_PDOS ,      &
                                     DeAllocate_SPEC 
 use FMO_m                  , only : FMO_analysis
 use QCModel_Huckel         , only : EigenSystem
 use DOS_m
 use Structure_Builder      , only : Generate_Structure ,   &
                                     Basis_Builder
 use Multipole_Core         , only : Dipole_Matrix
 use Oscillator_m           , only : Optical_Transitions
 use Schroedinger_m         , only : Huckel_dynamics ,      &
                                     DeAllocate_QDyn
 use RK_m
 use Data_Output            , only : Dump_stuff

 contains
!
!
!
!====================
subroutine qdynamics
!====================
implicit none 

! local variables ...
 integer                        :: nr , N_of_residues
 real*8          , allocatable  :: QDyn(:,:)
 character(3)                   :: residue
 character(3)    , allocatable  :: list_of_residues(:) 
 character(1)    , allocatable  :: list_of_fragments(:)
 logical                        :: FMO_ , DIPOLE_
 type(eigen)                    :: UNI
 type(eigen)                    :: FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 

 
! preprocessing stuff ...................................

FMO_    = ( spectrum .AND. survival  )
DIPOLE_ = ( FMO_     .OR.  DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

allocate( list_of_residues ( size(Unit_Cell     % list_of_residues ) ) , source = Unit_Cell     % list_of_residues  )
allocate( list_of_fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    residue = Unit_Cell % list_of_residues(nr)
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , residue , nr )            
 end do

 If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

 If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

! If( survival ) CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO , QDyn )

! CALL RK4_dynamics( Extended_Cell, ExCell_basis, UNI, FMO )

CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn , list_of_residues , list_of_fragments )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

include 'formats.h'

end subroutine qdynamics
!
!
!
end module qdynamics_m
