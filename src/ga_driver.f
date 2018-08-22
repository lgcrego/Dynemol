module GA_driver_m

 use type_m
 use constants_m
 use parameters_m               , only : spectrum , DP_Moment , GaussianCube , Alpha_Tensor
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use GA_m                       , only : Genetic_Algorithm 
 use GA_QCModel_m               , only : GA_eigen , Mulliken , GA_DP_Analysis , AlphaPolar
 use DOS_m 
 use Multipole_Routines_m       , only : Util_multipoles
 use Structure_Builder          , only : Generate_Structure , Extended_Cell , Unit_Cell , Basis_Builder , ExCell_basis
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
 integer                        :: nr , N_of_residues
 real*8                         :: DP(3) , Alpha_ii(3)
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
CALL GA_eigen( Extended_Cell, ExCell_basis, UNI )

If( DIPOLE_ ) CALL Util_multipoles

! Optimization of Huckel parameters ... 
CALL Genetic_Algorithm( ExCell_basis, OPT_basis )

! calculations with new parameters ...
CALL GA_eigen( Extended_Cell, OPT_basis, UNI )

CALL Total_DOS( UNI%erg, TDOS )

do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
end do

If( DIPOLE_ ) CALL GA_DP_Analysis( Extended_Cell, OPT_basis, UNI%L, UNI%R, DP )  

If( Alpha_Tensor .AND. DP_Moment ) CALL AlphaPolar( Extended_Cell, OPT_basis , Alpha_ii )

If( spectrum ) CALL Optical_Transitions( Extended_Cell, OPT_basis, UNI , SPEC )

!----------------------------------------------
! print zone ...

Print*, " " 
Print 154, DP, sqrt( dot_product(DP,DP) )
Print 189 , Alpha_ii , sum( Alpha_ii ) / three 

Print*, " " 
Print*, "dE1 = ",UNI%erg(56) - UNI%erg(55) , 3.49
Print*, "dE2 = ",UNI%erg(55) - UNI%erg(54) , 1.16
Print*, "dE3 = ",UNI%erg(57) - UNI%erg(56) , 1.7
Print*, "dE4 = ",UNI%erg(58) - UNI%erg(56) , 2.62
Print*, "dE5 = ",UNI%erg(58) - UNI%erg(57) , 0.92
Print*, " "
Print*, "dE1 = ",UNI%erg(166) - UNI%erg(165) , 2.16
Print*, "dE2 = ",UNI%erg(165) - UNI%erg(164) , 1.49
Print*, "dE3 = ",UNI%erg(167) - UNI%erg(166) , 1.94
Print*, "dE4 = ",UNI%erg(168) - UNI%erg(166) , 3.22
Print*, "dE5 = ",UNI%erg(168) - UNI%erg(167) , 1.28


! Population analysis ...
! print*,  Mulliken(UNI,ExCell_basis,MO=87,residue="TRI") 

!If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(6,:) , UNI%R(:,6) , 6 , 0.d0 )

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
