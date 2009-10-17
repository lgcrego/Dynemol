 Program qdynamo

 use type_m
 use constants_m
 use Sampling_m             , only : Solvated_M
 use FMO_m
 use QCModel_Huckel
 use projectors
 use DOS_m
 use Structure_Builder
 use Multipole_Core
 use Oscillator_m
 use Dynamics_m
 use RK_m

 type(eigen)    :: UNI
 type(eigen)    :: FMO
 type(f_grid)   :: TDOS , PDOS , SPEC
 
 complex*16 , parameter :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

!========================================================
!     starting up 
!========================================================
 
 CALL read_EHT_parameters

 CALL Read_Structure
 
 CALL Solvated_M

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

! CALL Partial_DOS( Extended_Cell, UNI , PDOS )

 CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

 CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

! CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO )

! CALL RK4_dynamics( Extended_Cell, ExCell_basis, UNI, FMO )

 include 'formats.h'

 END
