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

 public :: qdynamics

 private

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
 character(3)                   :: residue
 logical                        :: FMO_ , DIPOLE_
 type(C_eigen)                  :: UNI
 type(C_eigen)                  :: FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(f_time)                   :: QDyn

 
! preprocessing stuff ...................................

FMO_    = ( spectrum .AND. survival  )
DIPOLE_ = ( FMO_     .OR.  DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

 If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( survival ) CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO , QDyn )

! CALL RK4_dynamics( Extended_Cell, ExCell_basis, UNI, FMO )

CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn )

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
