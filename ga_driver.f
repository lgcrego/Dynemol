module GA_driver_m

 use type_m
 use constants_m
 use MPI_definitions_m          , only : slave
 use parameters_m               , only : spectrum , DP_Moment , GaussianCube , Alpha_Tensor
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use GA_m                       , only : Genetic_Algorithm 
 use GA_QCModel_m               , only : GA_eigen , Mulliken , GA_DP_Analysis , AlphaPolar
 use cost_EH                    , only : REF_DP , REF_Alpha
 use Multipole_Routines_m       , only : Util_multipoles
 use Structure_Builder          , only : Generate_Structure , Extended_Cell , Unit_Cell , Basis_Builder , ExCell_basis
 use Oscillator_m               , only : Optical_Transitions
 use Data_Output                , only : Dump_stuff
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format
 use DOS_m 

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
if( slave ) return

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
Print 152, REF_DP, sqrt( dot_product(REF_DP,REF_DP) )
!Print 189 , Alpha_ii , sum( Alpha_ii ) / three 

Print*, " "
Print*, "dE1 = ",UNI%erg(121) - UNI%erg(120) , 2.016d0
Print*, "dE2 = ",UNI%erg(120) - UNI%erg(119) , 0.842d0
Print*, "dE3 = ",UNI%erg(122) - UNI%erg(121) , 1.625d0
Print*, "dE4 = ",UNI%erg(123) - UNI%erg(121) , 2.359d0
Print*, "dE5 = ",UNI%erg(123) - UNI%erg(122) , 0.734d0
Print*, " "
Print*, "dE1 = ",UNI%erg(345) - UNI%erg(344) , 2.000d0
Print*, "dE2 = ",UNI%erg(344) - UNI%erg(343) , 0.102d0
Print*, "dE3 = ",UNI%erg(346) - UNI%erg(345) , 1.730d0
Print*, "dE4 = ",UNI%erg(347) - UNI%erg(345) , 3.429d0
Print*, "dE5 = ",UNI%erg(347) - UNI%erg(346) , 1.700d0

! Population analysis ...
! print*,  Mulliken(UNI,ExCell_basis,MO=87,residue="TRI") 

!If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(2,:) , UNI%R(:,2) , 2 , 0.d0 )

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
