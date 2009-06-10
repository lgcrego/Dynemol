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

 complex*16 , allocatable :: FMO_L(:,:) , FMO_R(:,:) 
 complex*16 , allocatable :: zL(:,:)    , zR(:,:)
 real*8     , allocatable :: erg(:)     , erg_FMO(:) 
 
 complex*16 , parameter   :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

!========================================================
!     starting up 
!========================================================
 
 CALL read_EHT_parameters

 CALL Read_Structure

 CALL Generate_Structure(t_i)

 CALL Basis_Builder(Extended_Cell,ExCell_basis)

 CALL eigen(Extended_Cell,ExCell_basis,zL,zR,erg)

 CALL TDOS(erg)

 CALL PDOS(extended_cell,zL,zR,erg)

 CALL FMO_analysis(Extended_Cell,zR,FMO_L,FMO_R,erg_FMO)

 CALL Dipole_Matrix(Extended_Cell,ExCell_basis,zL,zR)  

 CALL Optical_Transitions(Extended_Cell, ExCell_basis, DP_matrix_AO, zL, zR, erg)

 CALL Huckel_dynamics(Extended_Cell, ExCell_basis, zL, zR, FMO_L, FMO_R, erg)

! CALL Redfield_Equations(Extended_Cell, ExCell_basis, DP_matrix_AO, zL, zR, erg)

 include 'formats.h'

 END
