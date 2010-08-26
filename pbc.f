module PBC_m

    use type_m
    use Babel_m
    use constants_m
    use Allocation_m
    use Structure_Builder 
    use Semi_Empirical_Parms

contains
!
!
!
!=====================================================================
 subroutine Generate_Periodic_Structure( cell , pbc_cell , pbc_basis )
!=====================================================================
 implicit none
 type(structure)               , intent(in)  :: cell
 type(structure)               , intent(out) :: pbc_cell
 type(STO_basis) , allocatable , intent(out) :: pbc_basis(:)

! local variables ... 
 integer :: ix , iy , iz , k , n , copy

!----------------------------------------------------------
! (VIRTUAL) REPLICAS for Period Boundary Conditions

 CALL Allocate_Structures( (2*mmx+1)*(2*mmy+1)*(2*mmz+1)*cell%atoms , pbc_cell )

 pbc_cell % coord    (1:cell%atoms,1:3)  =  cell % coord
 pbc_cell % symbol   (1:cell%atoms)      =  cell % symbol
 pbc_cell % AtNo     (1:cell%atoms)      =  cell % AtNo
 pbc_cell % k_WH     (1:cell%atoms)      =  cell % k_WH
 pbc_cell % fragment (1:cell%atoms)      =  cell % fragment
 pbc_cell % Symbol   (1:cell%atoms)      =  cell % Symbol
 pbc_cell % MMSymbol (1:cell%atoms)      =  cell % MMSymbol
 pbc_cell % residue  (1:cell%atoms)      =  cell % residue
 pbc_cell % nr       (1:cell%atoms)      =  cell % nr 
 pbc_cell % copy_No  (1:cell%atoms)      =  0


! include the replicas        

 k = cell%atoms
 copy = 0

 DO iz = -mmz , mmz
 DO iy = -mmy , mmy
 DO ix = -mmx , mmx

    If( (ix /= 0) .OR. (iy /= 0) .OR. (iz /= 0) ) THEN 

        copy = copy + 1
        DO n = 1 , cell%atoms

            k = k + 1

            pbc_cell % coord    (k,1) =  cell % coord(n,1) + ix * cell%T_xyz(1)
            pbc_cell % coord    (k,2) =  cell % coord(n,2) + iy * cell%T_xyz(2)
            pbc_cell % coord    (k,3) =  cell % coord(n,3) + iz * cell%T_xyz(3)
            pbc_cell % AtNo     (k)   =  cell % AtNo     (n)
            pbc_cell % k_WH     (k)   =  cell % k_WH     (n)
            pbc_cell % fragment (k)   =  cell % fragment (n)
            pbc_cell % Symbol   (k)   =  cell % Symbol   (n)         
            pbc_cell % MMSymbol (k)   =  cell % MMSymbol (n)
            pbc_cell % residue  (k)   =  cell % residue  (n)
            pbc_cell % nr       (k)   =  cell % nr       (n) 
            pbc_cell % copy_No  (k)   =  copy

        END DO
    END IF

 END DO
 END DO
 END DO

!------------------------------------------------------------

 CALL Basis_Builder(pbc_cell,pbc_basis)

 end subroutine Generate_Periodic_Structure 
! 
!
!
!
!
!
end module PBC_m
