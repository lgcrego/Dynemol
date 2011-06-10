module qdynamics_m

 use type_m
 use constants_m
 use parameters_m           , only : spectrum , DP_Moment , &
                                     survival , DP_Field_ 
 use Solvated_M             , only : DeAllocate_TDOS ,      &
                                     DeAllocate_PDOS ,      &
                                     DeAllocate_SPEC 
 use FMO_m                  , only : FMO_analysis ,         &
                                     eh_tag 
 use QCModel_Huckel         , only : EigenSystem
 use DOS_m                  , only : Total_DOS ,            &
                                     Partial_DOS
 use Structure_Builder      , only : Generate_Structure ,   &
                                     Basis_Builder ,        &
                                     Unit_Cell ,            &
                                     Extended_Cell ,        &
                                     ExCell_Basis
 use DP_main_m              , only : Dipole_Matrix
 use DP_potential_m         , only : Molecular_DPs 
 use Oscillator_m           , only : Optical_Transitions
 use Schroedinger_m         , only : Huckel_dynamics ,      &
                                     DeAllocate_QDyn
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
 logical                        :: FMO_ , DIPOLE_ , el_hl_
 type(R_eigen)                  :: UNI
 type(R_eigen)                  :: el_FMO , hl_FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(f_time)                   :: QDyn

 
! preprocessing stuff ...................................

el_hl_  = any( Unit_Cell%fragment == "H")
FMO_    = ( spectrum .OR. survival  )
DIPOLE_ = ( spectrum .OR. DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 If( DP_field_ )CALL Molecular_DPs( Extended_Cell )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, el_FMO , instance="D")

 If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( survival ) then
    
    select case( el_hl_ )

        case( .false. )

            CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, el_FMO , QDyn=QDyn )

        case( .true. )

            CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, hl_FMO , instance="H")

            CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, el_FMO , hl_FMO , QDyn )

        end select

 end If

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
