module project

    type MM_molecular
        integer                             :: my_species
        integer                             :: N_of_atoms
        integer                             :: N_of_molecules
        real*8                              :: cm(3)
        real*8                              :: mass
        logical                             :: flex
        character(3)                        :: resid
        integer                             :: nresid
        integer                             :: Nbonds
        integer             , allocatable   :: bonds(:,:)
        real*8              , allocatable   :: kbond0(:,:)
        integer                             :: Nangs
        integer             , allocatable   :: angs(:,:)
        real*8              , allocatable   :: kang0(:,:)
        integer                             :: Ndiheds
        integer             , allocatable   :: diheds(:,:)
        real*8              , allocatable   :: kdihed0(:,:)
        integer                             :: Nharm
        integer                             :: Nbonds14
        integer             , allocatable   :: bonds14(:,:)
        real*8              , allocatable   :: fact14(:)
    end type MM_molecular

    type MM_atomic
        integer                             :: AtNo
        integer                             :: my_id
        integer                             :: my_species
        integer                             :: nresid
        character(3)                        :: resid
        character(2)                        :: Symbol
        character(3)                        :: MMSymbol
        real*8                              :: xyz(3)
        real*8                              :: vel(3)
        real*8                              :: fbond(3)
        real*8                              :: fang(3)
        real*8                              :: fdihed(3)
        real*8                              :: fnonbd(3)
        real*8                              :: fnonch(3)
        real*8                              :: fgeraldih(3)
        real*8                              :: ftotal(3)
        real*8                              :: fch(3)
        real*8                              :: fsr(3)
        real*8                              :: mass
        real*8                              :: charge
        real*8                              :: MM_charge
        real*8                              :: eps
        real*8                              :: sig
        logical                             :: free
    end type MM_atomic

    type MM_system
        type(MM_atomic)     , allocatable   :: atom(:)
        type(MM_molecular)  , allocatable   :: molecule(:)
        real*8                              :: box(3)
        real*8                              :: ibox(3)
        integer                             :: N_of_atoms
        integer                             :: N_of_species
        integer                             :: N_of_molecules
        character(72)                       :: System_Characteristics
    end type MM_system

end module project
!
