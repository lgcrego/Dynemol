 module Allocation_m

    use type_m
    use parameters_m    , only : n_part

    interface Allocate_BracKets
        module procedure Allocate_BracKets
        module procedure Allocate_BracKets_ElHl
        module procedure Allocate_BracKets_Chebyshev
    end interface Allocate_BracKets

 contains
!
!
!
!
 subroutine Allocate_UnitCell( unit_cell , n_residues )
 implicit none
 type(structure)            , intent(inout) :: unit_cell
 integer         , optional , intent(in)    :: n_residues

    allocate( unit_cell % Symbol             (unit_cell%atoms)   )
    allocate( unit_cell % AtNo               (unit_cell%atoms)   )
    allocate( unit_cell % Nvalen             (unit_cell%atoms)   )
    allocate( unit_cell % coord              (unit_cell%atoms,3) )
    allocate( unit_cell % k_WH               (unit_cell%atoms)   )
    allocate( unit_cell % fragment           (unit_cell%atoms)   )
    allocate( unit_cell % nr                 (unit_cell%atoms)   )
    allocate( unit_cell % residue            (unit_cell%atoms)   )
    allocate( unit_cell % MMSymbol           (unit_cell%atoms)   )
    allocate( unit_cell % QMMM               (unit_cell%atoms)   )
    allocate( unit_cell % solute             (unit_cell%atoms)   )
    allocate( unit_cell % DPF                (unit_cell%atoms)   )
    allocate( unit_cell % El                 (unit_cell%atoms)   )
    allocate( unit_cell % Hl                 (unit_cell%atoms)   )
    allocate( unit_cell % flex               (unit_cell%atoms)   )
    allocate( unit_cell % solvation_hardcore (unit_cell%atoms)   )
    allocate( unit_cell % hardcore           (unit_cell%atoms)   )
    allocate( unit_cell % V_shift            (unit_cell%atoms)   )
    allocate( unit_cell % polar              (unit_cell%atoms)   )

    If( present(n_residues) ) allocate( unit_cell % list_of_residues   (n_residues) )
    If( present(n_residues) ) allocate( unit_cell % list_of_fragments  (n_residues) )

    unit_cell%N_of_Solvent_Molecules = 0
 
 end subroutine Allocate_UnitCell
! 
!
!
!
 subroutine Allocate_Structures(System_Size,System)
    implicit none
    integer         , intent(in)    :: System_Size
    type(structure) , intent(inout) :: System

    System%atoms = System_Size

    allocate( System % BasisPointer       (System_Size)   )
    allocate( System % Symbol             (System_Size)   )
    allocate( System % AtNo               (System_size)   )
    allocate( System % Nvalen             (System_size)   )
    allocate( System % coord              (System_Size,3) )
    allocate( System % k_WH               (System_size)   )
    allocate( System % copy_No            (System_size)   )
    allocate( System % fragment           (System_size)   )
    allocate( System % nr                 (System_size)   )
    allocate( System % residue            (System_size)   )
    allocate( System % MMSymbol           (System_size)   )
    allocate( System % QMMM               (System_size)   )
    allocate( System % solute             (System_size)   )
    allocate( System % DPF                (System_size)   )
    allocate( System % El                 (System_size)   )
    allocate( System % Hl                 (System_size)   )
    allocate( System % flex               (System_size)   )
    allocate( System % solvation_hardcore (System_size)   )
    allocate( System % hardcore           (System_size)   )
    allocate( System % V_shift            (System_size)   )
    allocate( System % polar              (System_size)   )

    System%N_of_Solvent_Molecules = 0

 end subroutine Allocate_Structures
!
!
!
!
 subroutine Allocate_BracKets(Basis_Size, MO_bra, MO_ket, AO_bra, AO_ket, DUAL_bra, DUAL_ket, phase)
    implicit none
    integer                  , intent(in)  :: Basis_Size
    complex*16 , ALLOCATABLE , intent(out) :: MO_bra   (:,:) , MO_ket   (:,:)
    complex*16 , ALLOCATABLE , intent(out) :: AO_bra   (:,:) , AO_ket   (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: DUAL_ket (:,:) , DUAL_bra (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: phase(:) 

    allocate( MO_bra   (Basis_Size,n_part) , MO_ket   (Basis_Size,n_part) )
    allocate( AO_bra   (Basis_Size,n_part) , AO_ket   (Basis_Size,n_part) )
    allocate( DUAL_bra (Basis_Size,n_part) , DUAL_ket (Basis_Size,n_part) )
    allocate( phase    (Basis_Size) )

 end subroutine Allocate_BracKets
!
!
!
!
 subroutine Allocate_BracKets_ElHl(Basis_Size, MO_bra, MO_ket, AO_bra, AO_ket, DUAL_bra, DUAL_ket, phase)
    implicit none
    integer                  , intent(in)  :: Basis_Size
    complex*16 , ALLOCATABLE , intent(out) :: MO_bra   (:,:) , MO_ket   (:,:)
    complex*16 , ALLOCATABLE , intent(out) :: AO_bra   (:,:) , AO_ket   (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: DUAL_ket (:,:) , DUAL_bra (:,:) 
    complex*16 , ALLOCATABLE , intent(out) :: phase(:,:) 

    allocate( MO_bra   (Basis_Size,n_part) , MO_ket   (Basis_Size,n_part) )
    allocate( AO_bra   (Basis_Size,n_part) , AO_ket   (Basis_Size,n_part) )
    allocate( DUAL_bra (Basis_Size,n_part) , DUAL_ket (Basis_Size,n_part) )
    allocate( phase    (Basis_Size,n_part)                                )

 end subroutine Allocate_BracKets_ElHl
!
!
!
!
 subroutine Allocate_BracKets_Chebyshev(Basis_Size, AO_bra, AO_ket, DUAL_bra, DUAL_ket)
    implicit none
    integer                             , intent(in)  :: Basis_Size
    complex*16            , ALLOCATABLE , intent(out) :: AO_bra      (:,:) , AO_ket      (:,:) 
    complex*16            , ALLOCATABLE , intent(out) :: DUAL_ket    (:,:) , DUAL_bra    (:,:) 

    allocate( AO_bra      (Basis_Size,n_part) , AO_ket      (Basis_Size,n_part) )
    allocate( DUAL_bra    (Basis_Size,n_part) , DUAL_ket    (Basis_Size,n_part) )

    call GPU_Pin( AO_bra, Basis_Size*n_part*16 )
    call GPU_Pin( AO_ket, Basis_Size*n_part*16 )

 end subroutine Allocate_BracKets_Chebyshev
!
!
!
!
 subroutine DeAllocate_UnitCell(unit_cell)
    implicit none
    type(structure) , intent(inout) :: unit_cell

    deallocate( unit_cell % symbol             )
    deallocate( unit_cell % AtNo               )
    deallocate( unit_cell % Nvalen             )
    deallocate( unit_cell % coord              ) 
    deallocate( unit_cell % k_WH               )
    deallocate( unit_cell % fragment           )
    deallocate( unit_cell % nr                 )
    deallocate( unit_cell % residue            )
    deallocate( unit_cell % MMSymbol           )
    deallocate( unit_cell % QMMM               )
    deallocate( unit_cell % solute             )
    deallocate( unit_cell % DPF                )
    deallocate( unit_cell % El                 )
    deallocate( unit_cell % Hl                 )
    deallocate( unit_cell % flex               )
    deallocate( unit_cell % solvation_hardcore )
    deallocate( unit_cell % hardcore           )
    deallocate( unit_cell % V_shift            )
    deallocate( unit_cell % polar              )

    If( allocated(unit_cell%list_of_residues ) ) deallocate( unit_cell%list_of_residues  )
    If( allocated(unit_cell%list_of_fragments) ) deallocate( unit_cell%list_of_fragments )

 end subroutine DeAllocate_UnitCell
!
!
!
!
 subroutine DeAllocate_Structures(System)
    implicit none
    type(structure) , intent(inout) :: System

    deallocate( System % BasisPointer       )
    deallocate( System % symbol             )
    deallocate( System % AtNo               )
    deallocate( System % Nvalen             )
    deallocate( System % coord              )
    deallocate( System % k_WH               )
    deallocate( System % copy_No            )
    deallocate( System % fragment           )
    deallocate( System % nr                 )
    deallocate( System % residue            )
    deallocate( System % MMSymbol           )
    deallocate( System % QMMM               )
    deallocate( System % solute             )
    deallocate( System % DPF                )
    deallocate( System % El                 )
    deallocate( System % Hl                 )
    deallocate( System % flex               )
    deallocate( System % solvation_hardcore )
    deallocate( System % hardcore           )
    deallocate( System % V_shift            )
    deallocate( System % polar              )

    If( allocated(System%list_of_residues ) ) deallocate( System%list_of_residues  )
    If( allocated(System%list_of_fragments) ) deallocate( System%list_of_fragments )

 end subroutine DeAllocate_Structures
!
!
!
 end module Allocation_m
