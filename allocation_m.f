 module Allocation_m

    use type_m

 contains
!
!
!
!
 subroutine Allocate_UnitCell( unit_cell , n_residues )
 implicit none
 type(structure) , intent(inout) :: unit_cell
 integer         , intent(in)    :: n_residues

    allocate( unit_cell % Symbol            (unit_cell%atoms)   )
    allocate( unit_cell % AtNo              (unit_cell%atoms)   )
    allocate( unit_cell % Nvalen            (unit_cell%atoms)   )
    allocate( unit_cell % coord             (unit_cell%atoms,3) )
    allocate( unit_cell % k_WH              (unit_cell%atoms)   )
    allocate( unit_cell % fragment          (unit_cell%atoms)   )
    allocate( unit_cell % nr                (unit_cell%atoms)   )
    allocate( unit_cell % residue           (unit_cell%atoms)   )
    allocate( unit_cell % MMSymbol          (unit_cell%atoms)   )
    allocate( unit_cell % solute            (unit_cell%atoms)   )
    allocate( unit_cell % FMO               (unit_cell%atoms)   )
    allocate( unit_cell % list_of_residues  (n_residues)        )
    allocate( unit_cell % list_of_fragments (n_residues)        )

    unit_cell%N_of_Solvent_Molecules = 0
 
 end subroutine Allocate_UnitCell
! 
!
!
 subroutine Allocate_Structures(System_Size,System)

    integer         , intent(in)    :: System_Size
    type(structure) , intent(inout) :: System

    System%atoms = System_Size

    allocate( System % BasisPointer (System_Size)   )
    allocate( System % Symbol       (System_Size)   )
    allocate( System % AtNo         (System_size)   )
    allocate( System % Nvalen       (System_size)   )
    allocate( System % coord        (System_Size,3) )
    allocate( System % k_WH         (System_size)   )
    allocate( System % copy_No      (System_size)   )
    allocate( System % fragment     (System_size)   )
    allocate( System % nr           (System_size)   )
    allocate( System % residue      (System_size)   )
    allocate( System % MMSymbol     (System_size)   )
    allocate( System % solute       (System_size)   )
    allocate( System % FMO          (System_size)   )
 
    System%N_of_Solvent_Molecules = 0

 end subroutine Allocate_Structures
!
!
!
 subroutine Allocate_BracKets(Basis_Size, MO_bra, MO_ket, AO_bra, AO_ket, DUAL_bra, DUAL_ket, bra, ket, phase)

    integer                  , intent(in)  :: Basis_Size
    complex*16 , ALLOCATABLE , intent(out) :: MO_bra   (:,:) , MO_ket   (:,:)
    complex*16 , ALLOCATABLE , intent(out) :: AO_bra   (:,:) , AO_ket   (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: DUAL_ket (:,:) , DUAL_bra (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: phase(:) , bra(:) , ket(:)

    allocate( MO_bra   (Basis_Size,n_part) , MO_ket   (Basis_Size,n_part) )
    allocate( AO_bra   (Basis_Size,n_part) , AO_ket   (Basis_Size,n_part) )
    allocate( DUAL_bra (Basis_Size,n_part) , DUAL_ket (Basis_Size,n_part) )
    allocate( bra      (Basis_Size)        , ket      (Basis_Size)        )
    allocate( phase    (Basis_Size) )

 end subroutine Allocate_BracKets
!
!
!
!
 subroutine DeAllocate_UnitCell(unit_cell)

    type(structure) , intent(inout) :: unit_cell

    deallocate( unit_cell % symbol            )
    deallocate( unit_cell % AtNo              )
    deallocate( unit_cell % Nvalen            )
    deallocate( unit_cell % coord             )
    deallocate( unit_cell % k_WH              )
    deallocate( unit_cell % fragment          )
    deallocate( unit_cell % nr                )
    deallocate( unit_cell % residue           )
    deallocate( unit_cell % MMSymbol          )
    deallocate( unit_cell % solute            )
    deallocate( unit_cell % FMO               )
    deallocate( unit_cell % list_of_residues  )
    deallocate( unit_cell % list_of_fragments )
 
 end subroutine DeAllocate_UnitCell
!
!
!
!
 subroutine DeAllocate_Structures(System)

    type(structure) , intent(inout) :: System

    deallocate( System % BasisPointer )
    deallocate( System % symbol       )
    deallocate( System % AtNo         )
    deallocate( System % Nvalen       )
    deallocate( System % coord        )
    deallocate( System % k_WH         )
    deallocate( System % copy_No      )
    deallocate( System % fragment     )
    deallocate( System % nr           )
    deallocate( System % residue      )
    deallocate( System % MMSymbol     )
    deallocate( System % solute       )
    deallocate( System % FMO          )

    If( allocated(System%list_of_residues ) ) deallocate( System%list_of_residues  )
    If( allocated(System%list_of_fragments) ) deallocate( System%list_of_fragments )
 
 end subroutine DeAllocate_Structures
!
!
!
 end module Allocation_m
