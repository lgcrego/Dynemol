module GA_driver_m

 use type_m
 use constants_m
 use parameters_m               , only : spectrum , DP_Moment , GaussianCube
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use QCModel_Huckel             , only : EigenSystem
 use GA_m                       , only : Genetic_Algorithm 
 use GA_QCModel_m               , only : Mulliken
 use DOS_m 
 use Multipole_Routines_m       , only : Util_multipoles
 use Structure_Builder          , only : Generate_Structure , Extended_Cell , Unit_Cell , Basis_Builder , ExCell_basis
 use DP_main_m                  , only : Dipole_Matrix
 use Oscillator_m               , only : Optical_Transitions
 use Data_Output                , only : Dump_stuff
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format

 public :: GA_driver

 private

 contains
!
!
!
!====================
 subroutine GA_driver
!====================
implicit none 

! local variables ...
 integer                        :: i , nr , N_of_residues
 real*8                         :: DP(3)
 character(3)                   :: residue
 logical                        :: DIPOLE_
 type(R_eigen)                  :: UNI
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(STO_basis) , allocatable  :: OPT_basis(:)

 
! preprocessing stuff ...................................

DIPOLE_ = ( spectrum .OR. DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! setting up the system ...

CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell, ExCell_basis )

! setting up constraints ...
CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

If( DIPOLE_ ) CALL Util_multipoles

! Optimization of Huckel parameters ... 
CALL Genetic_Algorithm( Extended_Cell, ExCell_basis, OPT_basis )

! calculations with new parameters ...
CALL EigenSystem( Extended_Cell, OPT_basis, UNI )

CALL Total_DOS( UNI%erg, TDOS )

do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
end do

If( DIPOLE_ ) CALL Dipole_Matrix( Extended_Cell, OPT_basis, UNI%L, UNI%R, DP )  

If( spectrum ) CALL Optical_Transitions( Extended_Cell, OPT_basis, UNI , SPEC )

!----------------------------------------------
! print zone ...
Print 154, DP, sqrt( dot_product(DP,DP) )
Print*, " " 
Print*, UNI%erg(16) - UNI%erg(15)
Print*, UNI%erg(17) - UNI%erg(15)
Print*, UNI%erg(17) - UNI%erg(16)
Print*, UNI%erg(15) - UNI%erg(14)
Print*, " "

! Population analysis ...

If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(14,:) , UNI%R(:,14) , 14 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(15,:) , UNI%R(:,15) , 15 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(16,:) , UNI%R(:,16) , 16 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(17,:) , UNI%R(:,17) , 17 , 0.d0 )
!----------------------------------------------

CALL Dump_stuff( TDOS , PDOS , SPEC )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )

include 'formats.h'

end subroutine GA_driver
!
!
!
end module GA_driver_m
