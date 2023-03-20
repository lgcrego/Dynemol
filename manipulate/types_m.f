module types_m

    type atomic
        real*8              :: xyz(3)
        real*8              :: mass
        real*8              :: charge
        integer             :: indx
        integer             :: AtNo
        integer             :: my_intra_id
        integer             :: nrcg
        integer             :: nresid
        character(len=4)    :: resid
        character(len=3)    :: Symbol
        character(len=4)    :: MMSymbol
        character(len=4)    :: AASymbol
        character(len=4)    :: FFSymbol
        character(len=1)    :: TorF(3)
        character(len=1)    :: fragment
        logical             :: copy
        logical             :: delete
        logical             :: translate
        logical             :: rotate
        logical             :: group
    end type atomic


    type EHT
        character(len=3)    :: symbol
        integer             :: AtNo
        integer             :: Nvalen
        integer             :: Nzeta(0:3)
        integer             :: Nquant(0:3)
        integer             :: AngMax
        real*8              :: IP(0:3)
        real*8              :: zeta(0:3,2)
        real*8              :: coef(0:3,2)
        integer             :: DOS 
    end type EHT     

    type molecular
        type(atomic) , allocatable :: atom(:)
        real*8                     :: xyz(3)
        real*8                     :: radius
        integer                    :: indx
        integer                    :: N_of_atoms 
        character(len=4)           :: resid
        character(len=4)           :: residAA
        character(len=72)          :: Solvent_Characteristics
    end type molecular

    type universe
        type(atomic)    , allocatable :: atom(:)
        type(molecular)               :: solvent
        type(molecular)               :: dye
        type(molecular) , allocatable :: amino(:)
        real*8                        :: box(3)
        real*8                        :: Surface_Boundary
        real*8                        :: time
        integer                       :: N_of_atoms
        integer                       :: Total_conect
        integer                       :: N_of_aminos
        integer                       :: N_of_Surface_Atoms
        integer                       :: N_of_Solvent_Atoms
        integer                       :: N_of_Solvent_Molecules
        logical         , allocatable :: topol(:,:)
        character(1)    , allocatable :: list_of_fragments(:)
        character(3)    , allocatable :: list_of_residues(:)
        character(len=72)             :: Surface_Characteristics
        character(len=72)             :: System_Characteristics
    end type universe

    type real_interval
        real*8 :: inicio
        real*8 :: fim
    end type real_interval

    type integer_interval
        integer :: inicio
        integer :: fim
    end type integer_interval

    type real_pointer
        real*8 , pointer :: PTR => null()
    end type real_pointer

    type R3_vector
        real*8 , dimension(3) :: xyz
    end type R3_vector

end module types_m
