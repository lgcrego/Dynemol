module GA_driver_m

 use type_m
 use constants_m
 use parameters_m               , only : spectrum , DP_Moment , GaussianCube , Alpha_Tensor , OPT_parms
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use GA_m                       , only : Genetic_Algorithm 
 use GA_QCModel_m               , only : GA_eigen , Mulliken , GA_DP_Analysis , AlphaPolar
 use DOS_m 
 use Semi_Empirical_Parms       , only : EH_atom
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
 integer                        :: i , nr , N_of_residues , MO_total
 integer         , allocatable  :: MOnum(:)
 real*8                         :: DP(3) , Alpha_ii(3)
 character(6)                   :: MOstr
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

! reading command line arguments for plotting MO cube files ...
MO_total  = COMMAND_ARGUMENT_COUNT()
allocate( MOnum(MO_total) )
do i = 1 , MO_total
    CALL GET_COMMAND_ARGUMENT(i, MOstr)
    read( MOstr,*) MOnum(i)
end do
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


! Population analysis ...
print*,  "45  = " , Mulliken(UNI,OPT_basis,MO=29,atom=[4,5]   ) 

If( GaussianCube ) then
    do i = 1 , MO_total
        CALL Gaussian_Cube_Format( UNI%L(MOnum(i),:) , UNI%R(:,MOnum(i)) , MOnum(i) , 0.d0 )
    end do
    Print*, '>> Gaussian Cube done <<',MOnum(:)
end if

If( OPT_parms ) then    
     Print 445
     Print 45 , EH_atom%EHSymbol
else 
     Print*, ">> OPT_parms were not used <<"
end if

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
