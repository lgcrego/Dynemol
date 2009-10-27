module type_m

    use constants_m 

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

    type structure
        integer                    :: atoms 
        integer                    :: N_of_electrons
        integer                    :: molecule
        integer      , allocatable :: copy_No(:)
        integer      , allocatable :: BasisPointer(:)
        character(1) , allocatable :: fragment(:)
        character(2) , allocatable :: symbol(:)
        character(3) , allocatable :: residue(:)
        character(1) , allocatable :: list_of_fragments(:)
        character(3) , allocatable :: list_of_residues(:)
        character(3) , allocatable :: MMSymbol(:)
        integer      , allocatable :: AtNo(:)
        real*8       , allocatable :: coord(:,:)
        real*8       , allocatable :: k_WH(:)
        real*8                     :: BoxSides(3)
        real*8                     :: Center_of_Mass(3)
        real*8                     :: Center_of_Charge(3)
        real*8                     :: T_xyz(3)
    end type structure

!-----------------------------------------------------------------------------
    type atomic
        real*8                        :: xyz(3)
        real*8                        :: mass
        real*8                        :: charge
        integer                       :: AtNo
        integer                       :: nresid 
        integer                       :: copy_No
        character(3)                  :: residue
        character(3)                  :: Symbol
        character(3)                  :: MMSymbol
        character(1)                  :: TorF(3)
        character(1)                  :: fragment
    end type atomic

    type molecular
        type(atomic)    , allocatable :: atom(:) 
        real*8                        :: CG(3)
        real*8                        :: radius
        integer                       :: N_of_Atoms 
        integer                       :: nresid   
        integer                       :: copy_No
        character(3)                  :: residue 
        character(72)                 :: Solvent_Characteristics
    end type molecular

    type universe
        type(atomic)    , allocatable :: atom(:)
        type(molecular) , allocatable :: solvent(:)
        type(molecular)               :: dye
        real*8                        :: box(3)
        real*8                        :: Surface_Boundary
        integer                       :: N_of_Atoms
        integer                       :: N_of_Surface_Atoms
        integer                       :: N_of_Solvent_Atoms
        integer                       :: N_of_Solvent_Molecules
        character(1)    , allocatable :: list_of_fragments(:)
        character(3)    , allocatable :: list_of_residues(:)
        character(72)                 :: System_Characteristics
    end type universe
!-----------------------------------------------------------------------------


    type EHT
        character(len=2) :: symbol
        integer          :: AtNo
        integer          :: Nvalen
        integer          :: Nzeta(0:3)
        integer          :: Nquant(0:3)
        integer          :: AngMax
        real*8           :: IP(0:3)
        real*8           :: zeta(0:3,2)
        real*8           :: coef(0:3,2)
        integer          :: DOS 
    end type EHT     


    type STO_basis
        integer          :: n
        integer          :: l
        integer          :: m
        integer          :: atom
        integer          :: copy_No
        integer          :: AtNo
        integer          :: Nzeta
        real*8           :: IP
        real*8           :: k_WH
        real*8           :: coef(1:2)
        real*8           :: zeta(1:2)
        real*8           :: x
        real*8           :: y
        real*8           :: z
        character(len=2) :: symbol
        character(len=2) :: fragment
        character(len=3) :: MMSymbol
        character(len=3) :: residue 
    end type STO_basis


    type eigen
        complex*16 , allocatable :: R(:,:)
        complex*16 , allocatable :: L(:,:)
        real*8     , allocatable :: erg(:)
    end type eigen


    type spin_orbital
        integer :: orbital
        integer :: spin
    end type spin_orbital


    type R3_vector
        real*8 , dimension(3) :: dp
    end type R3_vector
 

    type C3_vector
        complex*16 , dimension(3) :: dp
    end type C3_vector
 

    type real_interval
        real*8 :: inicio 
        real*8 :: fim
    end type real_interval


    type integer_interval
        integer :: inicio
        integer :: fim
    end type integer_interval


    type int_pointer
        integer , pointer :: PTR => null()
    end type


    type transition
        type(R3_vector)        , allocatable :: matrix(:,:)
        type(real_interval)                  :: bra_range
        type(real_interval)                  :: ket_range
        type(integer_interval)               :: bra_indx_range
        type(integer_interval)               :: ket_indx_range
        integer                , allocatable :: bra_PTR(:)
        integer                , allocatable :: ket_PTR(:)
        character(len=8)                     :: flag
    end type transition


    type f_grid
        real*8  , allocatable   :: grid (:)
        real*8  , allocatable   :: func (:)
        real*8  , allocatable   :: peaks(:)
        real*8  , allocatable   :: average(:)
    end type f_grid


    include 'parameters.h'    

end module type_m
