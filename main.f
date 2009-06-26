 Program qdynamo

 use type_m
 use constants_m
 use FMO_m
 use QCModel_Huckel
 use projectors
 use DOS_m
 use Structure_Builder
 use Multipole_Core
 use Oscillator_m
 use Dynamics_m
 use QOptics_m

 type(eigen) :: UNI
 type(eigen) :: FMO
 
 complex*16 , parameter :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

!========================================================
!     starting up 
!========================================================
 
 CALL read_EHT_parameters

 CALL Read_Structure

 CALL Generate_Structure( t_i )

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL TDOS( UNI%erg )

 CALL PDOS( Extended_Cell, UNI )

 CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

 CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI )

 CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO )

! CALL Redfield_Equations( Extended_Cell, ExCell_basis, UNI, FMO )

 include 'formats.h'

 END
