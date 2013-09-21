module diagnostic_m

 use type_m
 use omp_lib
 use constants_m
 use parameters_m               , only : spectrum , DP_Moment , &
                                         survival , DP_Field_ , &
                                         Alpha_Tensor ,         &
                                         GaussianCube
 use Solvated_M                 , only : DeAllocate_TDOS ,      &
                                         DeAllocate_PDOS ,      &
                                         DeAllocate_SPEC 
 use QCModel_Huckel             , only : EigenSystem
 use DOS_m
 use Structure_Builder          , only : Unit_Cell ,            &
                                         Extended_Cell ,        &
                                         Generate_Structure ,   &
                                         Basis_Builder ,        &
                                         ExCell_basis
 use GA_QCModel_m               , only : Mulliken
 use DP_main_m                  , only : Dipole_Matrix
 use DP_potential_m             , only : Molecular_DPs
 use Oscillator_m               , only : Optical_Transitions
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format
 use Data_Output                , only : Dump_stuff
 use Embedded_FF_Alpha          , only : AlphaPolar

 public :: diagnostic

 private

 contains
!
!
!
!====================
subroutine diagnostic
!====================
implicit none 

! local variables ...
 integer                        :: i , nr , N_of_residues
 character(3)                   :: residue
 real*8                         :: DP(3)
 type(R_eigen)                  :: UNI
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 

 
! preprocessing stuff ...................................

IF ( survival ) pause " >>> quit: diagnostic driver does not carry q_dynamics calculations <<< "

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

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

 If( DP_Moment .OR. Spectrum ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R , DP )  

 If( Alpha_Tensor .AND. DP_Moment ) CALL AlphaPolar( Extended_Cell, ExCell_basis ) 

 If( Spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(4,:) , UNI%R(:,4) , 4 , 0.d0 )
 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(5,:) , UNI%R(:,5) , 5 , 0.d0 )

 CALL Dump_stuff( TDOS , PDOS , SPEC )

 CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
 CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
 CALL DeAllocate_SPEC( SPEC , flag="dealloc" )

 include 'formats.h'

end subroutine diagnostic
!
!
!
end module diagnostic_m
