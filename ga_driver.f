module GA_driver_m

 use type_m
 use constants_m
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use QCModel_Huckel             , only : EigenSystem
 use GA_m                       , only : Genetic_Algorithm
 use DOS_m
 use Structure_Builder          , only : Generate_Structure , Basis_Builder
 use Multipole_Core             , only : Dipole_Matrix
 use Oscillator_m               , only : Optical_Transitions
 use Data_Output                , only : Dump_stuff
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format

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
 real*8                         :: DP(3)
 character(3)                   :: residue
 character(3)    , allocatable  :: list_of_residues(:)
 logical                        :: DIPOLE_
 type(eigen)                    :: UNI
 type(eigen)                    :: FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(OPT)                      :: REF
 type(STO_basis) , allocatable  :: GA_basis(:)

 
! preprocessing stuff ...................................

DIPOLE_ = ( spectrum .OR. DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

allocate( list_of_residues ( size(Unit_Cell%list_of_residues) ) , source=Unit_Cell%list_of_residues )
!.........................................................

! setting up the system ...
CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell, ExCell_basis )

CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

! setting up constraints ...
allocate( REF%erg(size(UNI%erg)) )
REF%erg = UNI%erg
REF%DP  = [ 3.84 , 0.0 , 0.0 ]

! Optimization of Huckel parameters ... 
CALL Genetic_Algorithm( Extended_Cell, ExCell_basis, REF , GA_basis )

! calculations with new parameters ...
CALL EigenSystem( Extended_Cell, GA_basis, UNI )

CALL Total_DOS( UNI%erg, TDOS )

do nr = 1 , N_of_residues
    residue = Unit_Cell % list_of_residues(nr)
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , residue , nr )            
end do

If( DIPOLE_ ) CALL Dipole_Matrix( Extended_Cell, GA_basis, UNI%L, UNI%R, DP )  

Print 154, DP, sqrt( dot_product(DP,DP) )
print*, UNI%erg(9) - UNI%erg(8)
print*, UNI%erg(9) - UNI%erg(7)
print*, UNI%erg(9) - UNI%erg(6)

If( spectrum ) CALL Optical_Transitions( Extended_Cell, GA_basis, UNI , SPEC )

If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(6,:) , UNI%R(:,6) , 6 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(7,:) , UNI%R(:,7) , 7 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(8,:) , UNI%R(:,8) , 8 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(9,:) , UNI%R(:,9) , 9 , 0.d0 )

CALL Dump_stuff( TDOS , PDOS , SPEC , list_of_residues=list_of_residues )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )

include 'formats.h'

end subroutine GA_driver
!
!
!
end module GA_driver_m
