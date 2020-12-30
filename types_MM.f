module MM_types

use constants_m

public :: MM_atomic , MM_molecular , MM_system , DefineBonds , DefineAngles , DefinePairs , DefineMorse, debug_MM
public :: MMOPT_Control, Logicalkey

    type MM_atomic
        integer                             :: AtNo
        integer                             :: my_id
        integer                             :: my_intra_id
        integer                             :: my_species
        integer                             :: nr
        character(3)                        :: residue
        character(2)                        :: Symbol
        character(4)                        :: MMSymbol
        character(4)                        :: EHSymbol
        real*8                              :: xyz(3)
        real*8                              :: vel(3)
        real*8                              :: fbond(3)
        real*8                              :: fang(3)
        real*8                              :: fdihed(3)
        real*8                              :: fnonbd14(3)
        real*8                              :: fnonch14(3)
        real*8                              :: fnonbd(3)
        real*8                              :: fnonch(3)
        real*8                              :: fcoupling(3)
        real*8                              :: Ehrenfest(3)
        real*8                              :: f_MM(3)
        real*8                              :: ftotal(3)
        real*8                              :: fch(3)
        real*8                              :: fsr(3)
        real*8                              :: fMorse(3)
        real*8                              :: mass
        real*8                              :: charge
        real*8                              :: MM_charge
        real*8                              :: eps
        real*8                              :: eps14
        real*8                              :: sig
        real*8                              :: sig14
        logical                             :: flex
    end type MM_atomic

    type MM_molecular
        type(MM_atomic)     , allocatable   :: atom(:)
        integer                             :: my_species
        integer                             :: N_of_atoms
        integer                             :: N_of_molecules
        real*8                              :: cm(3)
        real*8                              :: mass
        logical                             :: flex
        character(3)                        :: residue
        integer                             :: nr
        integer                             :: Nbonds
        integer             , allocatable   :: bonds(:,:)
        real*8              , allocatable   :: kbond0(:,:)
        integer                             :: Nangs
        integer             , allocatable   :: angs(:,:)
        real*8              , allocatable   :: kang0(:,:)
        integer                             :: Ndiheds
        integer                             :: NTorsions
        integer                             :: NImpropers
        integer             , allocatable   :: diheds(:,:)
        character(3)        , allocatable   :: funct_bond(:)
        character(3)        , allocatable   :: funct_angle(:)
        integer             , allocatable   :: funct_dihed(:)
        real*8              , allocatable   :: kdihed0(:,:)
        character(4)        , allocatable   :: bond_type(:)
        character(4)        , allocatable   :: angle_type(:)
        character(4)        , allocatable   :: dihedral_type(:)
        integer                             :: Nharm
        integer                             :: Nbonds14
        integer             , allocatable   :: bonds14(:,:)
        integer                             :: NintraLJ
        integer             , allocatable   :: IntraLJ(:,:)
    end type MM_molecular

    type MM_system
        type(MM_atomic)     , allocatable   :: atom(:)
        type(MM_molecular)  , allocatable   :: molecule(:)        
        real*8                              :: box(3)
        real*8                              :: ibox(3)
        real*8                              :: fudgeLJ
        real*8                              :: fudgeQQ
        integer                             :: N_of_atoms
        integer                             :: N_of_species
        integer                             :: N_of_molecules
        integer                             :: N_of_AtomTypes        
        integer                             :: CombinationRule
    end type MM_system

    type DefineBonds
        character(15)                       :: label
        real*8                              :: kbond0(2)
    end type DefineBonds

    type DefineAngles
        character(15)                       :: label
        real*8                              :: kang0(2)
    end type DefineAngles

    type DefinePairs
        character(4)                        :: MMSymbols(2)
        real*8                              :: Parms(2)
    end type DefinePairs

    type DefineMorse
        character(3)                        :: MMSymbols(2)
        real*8                              :: Parms(3)
    end type DefineMorse

    type LogicalKey
        logical       :: bonds(3)    = .false.
        logical       :: angs(4)     = .false.
        logical       :: diheds(15)  = .false.
        logical       :: adiabatic   = .false.
        integer       :: dihedtype   
        character(20) :: comment
    end type LogicalKey

    type MMOPT_Control
        logical       :: adiabatic_OPT  = .false.
        logical       :: preprocess     = .true.
        logical       :: use_no_weights = .false.
        logical       :: new_adiabat    = .false.
        logical       :: LineUpCost     = .false.
        logical       :: use_overweight = .false.
    end type MMOPT_Control

    interface debug_MM
        module procedure debug_MM_atomic 
        module procedure debug_MM_molecular
        module procedure debug_MM_system
    end interface debug_MM
 
    private

contains
!
!
!
!===============================
 subroutine debug_MM_atomic( a )
!===============================
implicit none
type(MM_atomic) , intent(in) :: a(:)

! local variables ...
integer :: option

do

    print'("")'      
    write(*,*) ' (0) QUIT        '
    print'("")'      
    write(*,*) ' (1) Symbol      '
    write(*,*) ' (2) MMSymbol    '
    write(*,*) ' (3) EHSymbol    '
    write(*,*) ' (4) AtNo        '
    write(*,*) ' (5) my_id       '
    write(*,*) ' (6) my_intra_id '
    write(*,*) ' (7) my_species  '
    write(*,*) ' (8) nr          '
    write(*,*) ' (9) residue     '
    write(*,*) ' (10) charge     '
    write(*,*) ' (11) MM_charge  '
    write(*,*) ' (12) eps        '
    write(*,*) ' (13) eps14      '
    write(*,*) ' (14) sig        '
    write(*,*) ' (15) sig14      '
    write(*,*) ' (16) flex       '
    write(*,*) ' (17) mass       '

    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,10) a(:) % Symbol

        case(2)
            write(*,20) a(:) % MMSymbol

        case(3)
            write(*,20) a(:) % EHSymbol

        case(4)
            write(*,30) a(:) % AtNo

        case(5)
            write(*,40) a(:) % my_id

        case(6)
            write(*,40) a(:) % my_intra_id

        case(7)
            write(*,30) a(:) % my_species

        case(8)
            write(*,50) a(:) % nr

        case(9)
            write(*,20) a(:) % residue

        case(10)
            write(*,60) a(:) % charge

        case(11)
            write(*,60) a(:) % MM_charge

        case(12)
            write(*,60) a(:) % eps

        case(13)
            write(*,60) a(:) % eps14

        case(14)
            write(*,60) a(:) % sig
            
        case(15)
            write(*,60) a(:) % sig14

        case(16)
            write(*,70) a(:) % flex

        case(17)
            write(*,80) a(:) % mass

        case default
            exit

    end select

end do

10 Format(50a3)
20 Format(37a4)
30 Format(50i3)
40 Format(30i5)
50 Format(37i4)
60 Format(18F8.4)
70 Format(50L3)
80 Format(14F12.3)

end subroutine debug_MM_atomic
!
!
!
!====================================
 subroutine debug_MM_molecular( mol )
!====================================
implicit none
type(MM_molecular) , intent(in) :: mol(:)

! local variables ...
integer :: option

!local parameters ...
real*8  , parameter :: kmol = 6.02214129d26  ! mol X 1000

do

    print'("")'      
    write(*,*) ' (0) QUIT            '
    print'("")'      
    write(*,*) ' (1)  my_species     '
    write(*,*) ' (2)  N_of_atoms     '
    write(*,*) ' (3)  N_of_molecules '
    write(*,*) ' (4)  nr             '
    write(*,*) ' (5)  residue        '
    write(*,*) ' (6)  flex           '
    write(*,*) ' (7)  Nbonds         '
    write(*,*) ' (8)  Nangs          '
    write(*,*) ' (9)  Ndiheds        '
    write(*,*) ' (10) Nharm          '
    write(*,*) ' (11) Nbonds14       '
    write(*,*) ' (12) mass (in u.a.)   '

    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,10) mol(:) % my_species

        case(2)
            write(*,20) mol(:) % N_of_atoms

        case(3)
            write(*,20) mol(:) % N_of_molecules

        case(4)
            write(*,20) mol(:) % nr

        case(5)
            write(*,30) mol(:) % residue

        case(6)
            write(*,40) mol(:) % flex

        case(7)
            write(*,20) mol(:) % Nbonds

        case(8)
            write(*,20) mol(:) % Nangs

        case(9)
            write(*,20) mol(:) % Ndiheds

        case(10)
            write(*,20) mol(:) % Nharm

        case(11)
            write(*,20) mol(:) % Nbonds14

        case(12)
            write(*,50) mol(:) % mass * kmol

        case default
            exit

    end select

end do

10 Format(50i3)
20 Format(33i5)
30 Format(37a4)
40 Format(50L3)
50 Format(14F12.3)

end subroutine debug_MM_molecular
!
!
!
!==================================
 subroutine debug_MM_system ( sys )
!==================================
implicit none
type(MM_system) , intent(in) :: sys

! local variables ...
integer :: option

do

    print'("")'      
    write(*,*) ' (0) QUIT            '
    print'("")'      
    write(*,*) ' (1) N_of_atoms      '
    write(*,*) ' (2) N_of_species    '
    write(*,*) ' (3) N_of_molecules  '
    write(*,*) ' (4) N_of_AtomTypes  '

    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,'(a13,i5)') "N_of_atoms = ", sys % N_of_atoms 

        case(2)
            write(*,'(a15,i3)') "N_of_species = ", sys % N_of_species   

        case(3)
            write(*,'(a17,i4)') "N_of_molecules = ", sys  % N_of_molecules 

        case(4)
            write(*,'(a17,i3)') "N_of_AtomTypes = ", sys  % N_of_AtomTypes 

        case default
            exit

    end select

end do

end  subroutine debug_MM_system 
!
!
!
end module MM_types
