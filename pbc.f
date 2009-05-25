module PBC_m

    use type_m
    use Babel_m
    use constants_m
    use Allocation_m
    use EHT_parameters
    use Structure_Builder 

contains
!
!
!
!
 subroutine Generate_Periodic_Structure(cell,pbc_cell,pbc_basis)

 type(structure)               , intent(in)  :: cell
 type(structure)               , intent(out) :: pbc_cell
 type(STO_basis) , allocatable , intent(out) :: pbc_basis(:)

 integer :: ix , iy

!----------------------------------------------------------
! (VIRTUAL) REPLICAS for Period Boundary Conditions

 CALL Allocate_Structures( (2*mmx+1)*(2*mmy+1)*cell%atoms , pbc_cell )

 pbc_cell%coord   (1:cell%atoms,1:3)  =  cell%coord
 pbc_cell%symbol  (1:cell%atoms)      =  cell%symbol
 pbc_cell%AtNo    (1:cell%atoms)      =  cell%AtNo
 pbc_cell%k_WH    (1:cell%atoms)      =  cell%k_WH
 pbc_cell%copy_No (1:cell%atoms)      =  0

! include the replicas        

 k = cell%atoms
 copy = 0

 DO ix = -mmx , mmx
 DO iy = -mmy , mmy

    If( (ix /= 0) .OR. (iy /= 0) ) THEN 

        copy = copy + 1
        DO n = 1 , cell%atoms

            k = k + 1

            pbc_cell%coord   (k,1) =  cell%coord(n,1) + ix * cell%T_xyz(1)
            pbc_cell%coord   (k,2) =  cell%coord(n,2) + iy * cell%T_xyz(2)
            pbc_cell%coord   (k,3) =  cell%coord(n,3)
            pbc_cell%AtNo    (k)   =  cell%AtNo(n)
            pbc_cell%k_WH    (k)   =  cell%k_WH(n)
            pbc_cell%copy_No (k)   =  copy
            pbc_cell%symbol  (k)   =  atom( extended_cell%AtNo(n) )%symbol

        END DO
    END IF

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
