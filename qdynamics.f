module qdynamics_m

 use type_m
 use constants_m
 use parameters_m         , only : spectrum , DP_Moment , &
                                   survival , EnvField_ , &
                                   NetCharge
 use Solvated_M           , only : DeAllocate_TDOS ,      &
                                   DeAllocate_PDOS ,      &
                                   DeAllocate_SPEC 
 use QCModel_Huckel       , only : EigenSystem
 use DOS_m                , only : Total_DOS ,            &
                                   Partial_DOS
 use Structure_Builder    , only : Generate_Structure ,   &
                                   Basis_Builder ,        &
                                   Unit_Cell ,            &
                                   Extended_Cell 
 use DP_main_m            , only : Dipole_Matrix
 use Dielectric_Potential , only : Environment_SetUp 
 use Oscillator_m         , only : Optical_Transitions
 use Schroedinger_m       , only : Simple_dynamics ,      &
                                   DeAllocate_QDyn
 use Data_Output          , only : Dump_stuff , Net_Charge

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
 type(R_eigen)                  :: UNI
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(f_time)                   :: QDyn
 type(STO_basis) , allocatable  :: ExCell_basis(:)

 
! preprocessing stuff ...................................

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 If( NetCharge .AND. (.NOT. allocated(Net_Charge)) ) allocate( Net_Charge(Extended_Cell%atoms) , source = D_zero )

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 If( any([DP_Moment,Spectrum,EnvField_]) ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( EnvField_ )CALL Environment_SetUp( Extended_Cell )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( survival ) CALL Simple_dynamics( Extended_Cell , ExCell_basis , UNI , QDyn )

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
