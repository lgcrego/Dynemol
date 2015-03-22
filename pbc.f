module PBC_m

    use type_m
    use parameters_m            , only  : PBC
    use Babel_m
    use constants_m
    use Allocation_m
    use Structure_Builder 
    use Semi_Empirical_Parms

    public ::  Generate_Periodic_Structure , Generate_Periodic_DPs

    private

    interface Generate_Periodic_Structure 
        module procedure default_Periodic_Structure 
        module procedure GACG_Periodic_Structure 
    end interface

contains
!
!
!
!====================================================================
 subroutine default_Periodic_Structure( cell , pbc_cell , pbc_basis )
!====================================================================
 implicit none
 type(structure)               , intent(in)  :: cell
 type(structure)               , intent(out) :: pbc_cell
 type(STO_basis) , allocatable , intent(out) :: pbc_basis(:)

! local variables ... 
 integer :: ix , iy , iz , k , n , copy

!----------------------------------------------------------
! (VIRTUAL) REPLICAS for Period Boundary Conditions

 CALL Allocate_Structures( product(2*PBC(:)+1)*cell%atoms , pbc_cell )

 pbc_cell % coord              (1:cell%atoms,1:3)  =  cell % coord
 pbc_cell % nr                 (1:cell%atoms)      =  cell % nr 
 pbc_cell % BasisPointer       (1:cell%atoms)      =  cell % BasisPointer
 pbc_cell % fragment           (1:cell%atoms)      =  cell % fragment
 pbc_cell % symbol             (1:cell%atoms)      =  cell % symbol
 pbc_cell % residue            (1:cell%atoms)      =  cell % residue
 pbc_cell % MMSymbol           (1:cell%atoms)      =  cell % MMSymbol
 pbc_cell % AtNo               (1:cell%atoms)      =  cell % AtNo
 pbc_cell % Nvalen             (1:cell%atoms)      =  cell % Nvalen
 pbc_cell % k_WH               (1:cell%atoms)      =  cell % k_WH
 pbc_cell % solvation_hardcore (1:cell%atoms)      =  cell % solvation_hardcore
 pbc_cell % hardcore           (1:cell%atoms)      =  cell % hardcore
 pbc_cell % polar              (1:cell%atoms)      =  cell % polar
 pbc_cell % solute             (1:cell%atoms)      =  cell % solute
 pbc_cell % DPF                (1:cell%atoms)      =  cell % DPF   
 pbc_cell % El                 (1:cell%atoms)      =  cell % El    
 pbc_cell % Hl                 (1:cell%atoms)      =  cell % Hl    
 pbc_cell % copy_No  (1:cell%atoms)                =  0


! include the replicas        

 k = cell%atoms
 copy = 0

 DO iz = -PBC(3) , PBC(3)
 DO iy = -PBC(2) , PBC(2)
 DO ix = -PBC(1) , PBC(1)

    If( (ix /= 0) .OR. (iy /= 0) .OR. (iz /= 0) ) THEN 

        copy = copy + 1
        DO n = 1 , cell%atoms

            k = k + 1

            pbc_cell % coord              (k,1) =  cell % coord(n,1) + ix * cell%T_xyz(1)
            pbc_cell % coord              (k,2) =  cell % coord(n,2) + iy * cell%T_xyz(2)
            pbc_cell % coord              (k,3) =  cell % coord(n,3) + iz * cell%T_xyz(3)

            pbc_cell % nr                 (k)   =  cell % nr                 (n) 
            pbc_cell % BasisPointer       (k)   =  cell % BasisPointer       (n) 
            pbc_cell % fragment           (k)   =  cell % fragment           (n)
            pbc_cell % Symbol             (k)   =  cell % Symbol             (n)         
            pbc_cell % residue            (k)   =  cell % residue            (n)
            pbc_cell % MMSymbol           (k)   =  cell % MMSymbol           (n)
            pbc_cell % AtNo               (k)   =  cell % AtNo               (n)
            pbc_cell % Nvalen             (k)   =  cell % Nvalen             (n)
            pbc_cell % k_WH               (k)   =  cell % k_WH               (n)
            pbc_cell % solvation_hardcore (k)   =  cell % solvation_hardcore (n)
            pbc_cell % hardcore           (k)   =  cell % hardcore           (n)
            pbc_cell % polar              (k)   =  cell % polar              (n)
            pbc_cell % solute             (k)   =  cell % solute             (n) 
            pbc_cell % DPF                (k)   =  cell % DPF                (n) 
            pbc_cell % El                 (k)   =  cell % El                 (n) 
            pbc_cell % Hl                 (k)   =  cell % Hl                 (n) 
            pbc_cell % copy_No  (k)             =  copy

        END DO
    END IF

 END DO
 END DO
 END DO

 pbc_cell % atoms = product(2*PBC(:)+1)*cell%atoms 

!------------------------------------------------------------

 CALL Basis_Builder( pbc_cell , pbc_basis )

 end subroutine default_Periodic_Structure 
! 
!
!
!
!=======================================================================================================
 subroutine Generate_Periodic_DPs( a , CC_Mol , DP_Mol , nr_Mol , pbc_CC_Mol , pbc_DP_MOL , pbc_nr_MOL )
!=======================================================================================================
 implicit none
 type(structure)                  , intent(in)  :: a
 real*8                           , intent(in)  :: CC_Mol(:,:)
 real*8                           , intent(in)  :: DP_Mol(:,:)
 integer                          , intent(in)  :: nr_Mol(:)
 real*8           , allocatable   , intent(out) :: pbc_CC_Mol(:,:)
 real*8           , allocatable   , intent(out) :: pbc_DP_MOL(:,:)
 integer          , allocatable   , intent(out) :: pbc_nr_MOL(:)

! local variables ... 
integer :: ix , iy , iz , j , k , n , N_of_DP_MOLs , N_of_pbc_mols

! number of DPs in the original cell ...
N_of_DP_MOLs = size( nr_MOL )

! (VIRTUAL) REPLICAS for Period Boundary Conditions ...
N_of_pbc_mols = product(2*PBC(:)+1) * N_of_DP_MOLs

allocate( pbc_CC_MOL ( N_of_pbc_mols , 3 ) )
allocate( pbc_DP_MOL ( N_of_pbc_mols , 3 ) )
allocate( pbc_nr_MOL ( N_of_pbc_mols )     )

! original DPs cell ...
do j = 1 , 3  
    pbc_CC_MOL(:,j) = CC_MOL(:,j)
    pbc_DP_MOL(:,j) = DP_Mol(:,j)
end do

pbc_nr_MOL = nr_Mol 

! include the replicas        
k = N_of_DP_MOLs

DO iz = -PBC(3) , PBC(3)
DO iy = -PBC(2) , PBC(2)
DO ix = -PBC(1) , PBC(1)

    If( (ix /= 0) .OR. (iy /= 0) .OR. (iz /= 0) ) THEN 

        DO n = 1 , N_of_DP_MOLs

            k = k + 1

            pbc_CC_MOL (k,1) =  CC_MOL(n,1) + ix * a % T_xyz(1)
            pbc_CC_MOL (k,2) =  CC_MOL(n,2) + iy * a % T_xyz(2)
            pbc_CC_MOL (k,3) =  CC_MOL(n,3) + iz * a % T_xyz(3)

            forall( j =1:3 ) pbc_DP_MOL(k,j) = DP_Mol(n,j)

            pbc_nr_MOL(k) = nr_Mol(n)

        END DO
    END IF

END DO
END DO
END DO

end subroutine Generate_Periodic_DPs
!
!
!
!
!=========================================================================
 subroutine GACG_Periodic_Structure( cell , basis , pbc_cell , pbc_basis )
!=========================================================================
 implicit none
 type(structure)               , intent(in)  :: cell
 type(STO_basis)               , intent(in)  :: basis(:)
 type(structure)               , intent(out) :: pbc_cell
 type(STO_basis) , allocatable , intent(out) :: pbc_basis(:)

! local variables ... 
 integer :: ix , iy , iz , k , n , copy

!----------------------------------------------------------
! (VIRTUAL) REPLICAS for Period Boundary Conditions

 CALL Allocate_Structures( product(2*PBC(:)+1)*cell%atoms , pbc_cell )

 pbc_cell % coord              (1:cell%atoms,1:3)  =  cell % coord
 pbc_cell % nr                 (1:cell%atoms)      =  cell % nr 
 pbc_cell % BasisPointer       (1:cell%atoms)      =  cell % BasisPointer
 pbc_cell % fragment           (1:cell%atoms)      =  cell % fragment
 pbc_cell % symbol             (1:cell%atoms)      =  cell % symbol
 pbc_cell % residue            (1:cell%atoms)      =  cell % residue
 pbc_cell % MMSymbol           (1:cell%atoms)      =  cell % MMSymbol
 pbc_cell % AtNo               (1:cell%atoms)      =  cell % AtNo
 pbc_cell % Nvalen             (1:cell%atoms)      =  cell % Nvalen
 pbc_cell % k_WH               (1:cell%atoms)      =  cell % k_WH
 pbc_cell % solvation_hardcore (1:cell%atoms)      =  cell % solvation_hardcore
 pbc_cell % hardcore           (1:cell%atoms)      =  cell % hardcore
 pbc_cell % polar              (1:cell%atoms)      =  cell % polar
 pbc_cell % solute             (1:cell%atoms)      =  cell % solute
 pbc_cell % DPF                (1:cell%atoms)      =  cell % DPF   
 pbc_cell % El                 (1:cell%atoms)      =  cell % El    
 pbc_cell % Hl                 (1:cell%atoms)      =  cell % Hl    
 pbc_cell % copy_No  (1:cell%atoms)                =  0


! include the replicas        

 k = cell%atoms
 copy = 0

 DO iz = -PBC(3) , PBC(3)
 DO iy = -PBC(2) , PBC(2)
 DO ix = -PBC(1) , PBC(1)

    If( (ix /= 0) .OR. (iy /= 0) .OR. (iz /= 0) ) THEN 

        copy = copy + 1
        DO n = 1 , cell%atoms

            k = k + 1

            pbc_cell % coord              (k,1) =  cell % coord(n,1) + ix * cell%T_xyz(1)
            pbc_cell % coord              (k,2) =  cell % coord(n,2) + iy * cell%T_xyz(2)
            pbc_cell % coord              (k,3) =  cell % coord(n,3) + iz * cell%T_xyz(3)

            pbc_cell % nr                 (k)   =  cell % nr                 (n) 
            pbc_cell % BasisPointer       (k)   =  cell % BasisPointer       (n) 
            pbc_cell % fragment           (k)   =  cell % fragment           (n)
            pbc_cell % Symbol             (k)   =  cell % Symbol             (n)         
            pbc_cell % residue            (k)   =  cell % residue            (n)
            pbc_cell % MMSymbol           (k)   =  cell % MMSymbol           (n)
            pbc_cell % AtNo               (k)   =  cell % AtNo               (n)
            pbc_cell % Nvalen             (k)   =  cell % Nvalen             (n)
            pbc_cell % k_WH               (k)   =  cell % k_WH               (n)
            pbc_cell % solvation_hardcore (k)   =  cell % solvation_hardcore (n)
            pbc_cell % hardcore           (k)   =  cell % hardcore           (n)
            pbc_cell % polar              (k)   =  cell % polar              (n)
            pbc_cell % solute             (k)   =  cell % solute             (n) 
            pbc_cell % DPF                (k)   =  cell % DPF                (n) 
            pbc_cell % El                 (k)   =  cell % El                 (n) 
            pbc_cell % Hl                 (k)   =  cell % Hl                 (n) 
            pbc_cell % copy_No  (k)             =  copy

        END DO
    END IF

 END DO
 END DO
 END DO

 pbc_cell % atoms = product(2*PBC(:)+1)*cell%atoms 

!------------------------------------------------------------

 CALL Basis_Builder( pbc_cell , pbc_basis , GACG_flag = .true. )

 do k = 1 , size(basis)

    where( (pbc_basis%EHSymbol == basis(k)%EHSymbol) .AND. (pbc_basis%l == basis(k)%l) )

        pbc_basis % IP       =  basis(k) % IP    
        pbc_basis % Nzeta    =  basis(k) % Nzeta 
        pbc_basis % coef(1)  =  basis(k) % coef  (1)
        pbc_basis % coef(2)  =  basis(k) % coef  (2)
        pbc_basis % zeta(1)  =  basis(k) % zeta  (1)
        pbc_basis % zeta(2)  =  basis(k) % zeta  (2)
        pbc_basis % k_WH     =  basis(k) % k_WH

    end where

 end do

 end subroutine GACG_Periodic_Structure 
! 
!
!
!
!
end module PBC_m
