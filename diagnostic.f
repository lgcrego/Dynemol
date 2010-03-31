module diagnostic_m

 use type_m
 use constants_m
 use Solvated_M                 , only : DeAllocate_TDOS ,      &
                                         DeAllocate_PDOS ,      &
                                         DeAllocate_SPEC 
 use FMO_m                      , only : FMO_analysis
 use QCModel_Huckel             , only : EigenSystem
 use DOS_m
 use Structure_Builder          , only : Generate_Structure ,   &
                                         Basis_Builder
 use Multipole_Core             , only : Dipole_Matrix
 use Oscillator_m               , only : Optical_Transitions
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format
 use Data_Output                , only : Dump_stuff

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
 logical                        :: FMO_ , DIPOLE_
 type(eigen)                    :: UNI
 type(eigen)                    :: FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 

 
! preprocessing stuff ...................................

FMO_    = ( spectrum .AND. survival  )
DIPOLE_ = ( FMO_     .OR.  DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO )

 If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 do i = 1 , 4
     print*, UNI%erg(504+i) + 12.d0
 enddo

 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(505,:) , UNI%R(:,505) , 505 , 0.d0 )
 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(506,:) , UNI%R(:,506) , 506 , 0.d0 )
 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(507,:) , UNI%R(:,507) , 507 , 0.d0 )
 If( GaussianCube ) CALL Gaussian_Cube_Format( UNI%L(508,:) , UNI%R(:,508) , 508 , 0.d0 )

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
