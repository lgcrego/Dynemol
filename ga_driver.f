module GA_driver_m

 use type_m
 use constants_m
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use QCModel_Huckel             , only : EigenSystem
 use GA_m                       , only : Genetic_Algorithm 
 use GA_QCModel_m               , only : Mulliken
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
 integer                        :: i , nr , N_of_residues
 real*8                         :: DP(3)
 character(3)                   :: residue
 logical                        :: DIPOLE_
 type(C_eigen)                  :: UNI
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

!.........................................................

! setting up the system ...

CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell, ExCell_basis )

! setting up constraints ...
CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

allocate( REF%erg(size(UNI%erg)) )
REF%erg = UNI%erg
REF%DP  = [ 0.00 , 0.0 , 0.0 ]

! Optimization of Huckel parameters ... 
CALL Genetic_Algorithm( Extended_Cell, ExCell_basis, REF , GA_basis )

! calculations with new parameters ...
CALL EigenSystem( Extended_Cell, GA_basis, UNI )

CALL Total_DOS( UNI%erg, TDOS )

do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
end do

If( DIPOLE_ ) CALL Dipole_Matrix( Extended_Cell, GA_basis, UNI%L, UNI%R, DP )  

If( spectrum ) CALL Optical_Transitions( Extended_Cell, GA_basis, UNI , SPEC )

!----------------------------------------------
! print zone ...
Print 154, DP, sqrt( dot_product(DP,DP) )
Print*, " " 
Print*, UNI%erg(12) - UNI%erg(11)
Print*, UNI%erg(12) - UNI%erg(10)
Print*, UNI%erg(11) - UNI%erg(10)
Print*, UNI%erg(11) - UNI%erg(09)
Print*, " " 
! Population analysis ...
print*, "09-Ix-l=1", Mulliken(UNI,GA_basis,MO=09,atom=1,AO_ang=1) 
print*, "09-Ic-l=0", Mulliken(UNI,GA_basis,MO=09,atom=2,AO_ang=0)
print*, "09-Ix-l=1", Mulliken(UNI,GA_basis,MO=09,atom=3,AO_ang=1)
Print*, " " 
print*, "10-Ix-l=1", Mulliken(UNI,GA_basis,MO=10,atom=1,AO_ang=1) 
print*, "10-Ic-l=1", Mulliken(UNI,GA_basis,MO=10,atom=2,AO_ang=1)
print*, "10-Ix-l=1", Mulliken(UNI,Ga_basis,MO=10,atom=3,AO_ang=1)
 Print*, " " 
print*, "11-Ix-l=1", Mulliken(UNI,GA_basis,MO=11,atom=1,AO_ang=1)
print*, "11-Ix-l=1", Mulliken(UNI,GA_basis,MO=11,atom=2,AO_ang=1)
print*, "11-Ix-l=1", Mulliken(UNI,GA_basis,MO=11,atom=3,AO_ang=1)
Print*, " " 
print*, "12-Ix-l=0", Mulliken(UNI,GA_basis,MO=12,atom=1,AO_ang=0)
print*, "12-Ix-l=1", Mulliken(UNI,GA_basis,MO=12,atom=1,AO_ang=1)
print*, "12-Ic-l=1", Mulliken(UNI,GA_basis,MO=12,atom=2,AO_ang=1)
print*, "12-Ix-l=0", Mulliken(UNI,GA_basis,MO=12,atom=3,AO_ang=0)
print*, "12-Ix-l=1", Mulliken(UNI,GA_basis,MO=12,atom=3,AO_ang=1)

If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(09,:) , UNI%R(:,09) , 09 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(10,:) , UNI%R(:,10) , 10 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(11,:) , UNI%R(:,11) , 11 , 0.d0 )
If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(12,:) , UNI%R(:,12) , 12 , 0.d0 )
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
