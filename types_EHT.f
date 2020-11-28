module type_m

    use constants_m 
    use execution_time_m

    type structure
        integer                    :: atoms 
        integer                    :: N_of_electrons
        integer                    :: N_of_Solvent_Molecules
        integer                    :: N_of_Solute_Molecules
        integer      , allocatable :: nr(:)
        integer      , allocatable :: copy_No(:)
        integer      , allocatable :: BasisPointer(:)
        character(1) , allocatable :: fragment(:)
        character(2) , allocatable :: symbol(:)
        character(3) , allocatable :: residue(:)
        character(1) , allocatable :: list_of_fragments(:)
        character(3) , allocatable :: list_of_residues(:)
        character(4) , allocatable :: MMSymbol(:)
        character(2) , allocatable :: QMMM(:)
        integer      , allocatable :: AtNo(:)
        integer      , allocatable :: Nvalen(:)
        real*8       , allocatable :: coord(:,:)
        real*8       , allocatable :: k_WH(:)
        real*8       , allocatable :: solvation_hardcore(:)
        real*8       , allocatable :: hardcore(:)
        real*8       , allocatable :: V_shift(:)
        real*8       , allocatable :: polar(:)
        real*8                     :: BoxSides(3)
        real*8                     :: Center_of_Mass(3)
        real*8                     :: Center_of_Charge(3)
        real*8                     :: T_xyz(3)
        real*8                     :: MD_Kin
        real*8                     :: MD_Pot
        real*8                     :: QM_erg
        complex*16                 :: QM_wp_erg(2)
        real*8                     :: Total_erg
        logical     , allocatable  :: solute(:)
        logical     , allocatable  :: DPF(:)
        logical     , allocatable  :: El(:)
        logical     , allocatable  :: Hl(:)
        logical     , allocatable  :: flex(:)
    end type structure

!-----------------------------------------------------------------------------
    type atomic
        real*8                        :: xyz(3)
        real*8                        :: mass
        real*8                        :: charge
        real*8                        :: solvation_hardcore
        real*8                        :: hardcore
        real*8                        :: polar
        real*8                        :: V_shift
        integer                       :: my_id
        integer                       :: AtNo
        integer                       :: Nvalen
        integer                       :: nr
        integer                       :: copy_No
        character(3)                  :: residue
        character(2)                  :: Symbol
        character(4)                  :: MMSymbol
        character(2)                  :: QMMM
        character(1)                  :: TorF(3)
        character(1)                  :: fragment
        logical                       :: solute
        logical                       :: DPF
        logical                       :: El
        logical                       :: Hl
        logical                       :: flex
    end type atomic

    type Point_Charges
        integer ,  allocatable  :: nr (:)
        real*8  ,  allocatable  :: Q  (:)
        real*8  ,  allocatable  :: xyz(:,:)
    end type Point_Charges


    type molecular
        type(Point_Charges)           :: PC
        type(atomic)    , allocatable :: atom(:) 
        real*8                        :: radius
        real*8                        :: CG(3)
        real*8                        :: CC(3)
        real*8                        :: DP(3)
        integer                       :: N_of_Atoms 
        integer                       :: nr
        integer                       :: copy_No
        character(3)                  :: residue 
        character(72)                 :: Solvent_Characteristics
        logical                       :: solute
        logical                       :: DPF
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
        character(len=2) :: Symbol
        character(len=4) :: EHSymbol
        character(len=3) :: residue 
        integer          :: AtNo
        integer          :: Nvalen
        integer          :: Nzeta(0:3)
        integer          :: Nquant(0:3)
        integer          :: AngMax
        integer          :: n
        integer          :: l
        real*8           :: IP(0:3)
        real*8           :: zeta(0:3,2)
        real*8           :: coef(0:3,2)
        real*8           :: k_WH(0:3)
        integer          :: DOS 
        real*8           :: polar
    end type EHT     


    type STO_basis
        integer          :: indx
        integer          :: atom
        integer          :: nr
        integer          :: copy_No
        integer          :: AtNo
        integer          :: Nzeta
        integer          :: n
        integer          :: l
        integer          :: m
        integer          :: s
        real*8           :: j
        real*8           :: IP
        real*8           :: k_WH
        real*8           :: coef(1:2)
        real*8           :: zeta(1:2)
        real*8           :: x
        real*8           :: y
        real*8           :: z
        real*8           :: solvation_hardcore
        real*8           :: hardcore
        real*8           :: V_shift
        character(len=2) :: symbol
        character(len=1) :: fragment
        character(len=4) :: EHSymbol
        character(len=3) :: residue 
        logical          :: solute
        logical          :: DPF
        logical          :: El
        logical          :: Hl
        logical          :: flex
        logical          :: modified = .false.
    end type STO_basis


    type C_eigen
        complex*16 , allocatable :: R(:,:)
        complex*16 , allocatable :: L(:,:)
        real*8     , allocatable :: erg(:)
        integer                  :: Fermi_state
    end type C_eigen


    type R_eigen
        real*8     , allocatable :: R(:,:)
        real*8     , allocatable :: L(:,:)
        real*8     , allocatable :: erg(:)
        real*8                   :: time
        integer                  :: Fermi_state
    end type R_eigen


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


    type real_pointer
        real*8 , pointer :: PTR => null()
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
        real*8       , allocatable   :: grid       (:)
        real*8       , allocatable   :: func       (:)
        real*8       , allocatable   :: peaks      (:)
        real*8       , allocatable   :: occupation (:)
        real*8       , allocatable   :: average    (:)
        character(1)                 :: fragment
        character(3)                 :: residue
    end type f_grid


    type f_time
        real*8       , allocatable   :: dyn(:,:,:)     ! <== dyn( time , fragments , el/hl )
        real*8       , allocatable   :: std(:,:,:)     ! <== std( time , fragments , el/hl )
        character(1) , allocatable   :: fragments(:)
        character(3) , allocatable   :: residues(:)
    end type f_time

    type on_the_fly 
        logical :: mode = .false.  ! <== must turn on before use ...
        integer :: gen
        integer :: Ngen
    end type on_the_fly

    type dipoles
        integer ,  allocatable  :: nr   (:)
        real*8  ,  allocatable  :: CC   (:,:)
        real*8  ,  allocatable  :: DP   (:,:)
        real*8  ,  allocatable  :: el_DP(:,:)
        real*8  ,  allocatable  :: hl_DP(:,:)
    end type dipoles


    interface debug_EH
        module procedure debug_EH_structure
        module procedure debug_EH_atomic 
        module procedure debug_EH_STO_basis
    end interface debug_EH
 
contains
!
!
!
!==================================
 subroutine debug_EH_structure( a )
!==================================
implicit none
type(structure) , intent(in) :: a

! local variables ...
integer :: option

do

    print'("")'      
    write(*,*) ' (0) QUIT                       '
    print'("")'      
    write(*,*) ' (1)  Symbol                    '
    write(*,*) ' (2)  MMSymbol                  '
    write(*,*) ' (3)  AtNo                      '
    write(*,*) ' (4)  fragment                  '
    write(*,*) ' (5)  nr                        '
    write(*,*) ' (6)  residue                   '
    write(*,*) ' (7)  copy_No                   '
    write(*,*) ' (8)  BasisPointer              '
    write(*,*) ' (9)  list_of_fragments         '
    write(*,*) ' (10) list_of_residues          ' 
    write(*,*) ' (11) polar                     '
    write(*,*) ' (12) k_WH                      '
    write(*,*) ' (13) DPF                       '
    write(*,*) ' (14) El                        '
    write(*,*) ' (15) Hl                        '
    write(*,*) ' (16) flex                      '
    write(*,*) ' (17) solute                    '
    write(*,*) ' (18) QM-MM                     '
    write(*,*) ' (19) atoms                     '
    write(*,*) ' (20) N_of_electrons            '
    write(*,*) ' (21) N_of_solvent_Molecules    '
    write(*,*) ' (22) N_of_solute_Molecules     '
    write(*,*) ' any other number cotinues      '


    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,10) a % Symbol(:)

        case(2)
            write(*,20) a % MMSymbol(:)

        case(3)
            write(*,30) a % AtNo(:)

        case(4)
            write(*,10) a % fragment(:)

        case(5)
            write(*,50) a % nr(:)

        case(6)
            write(*,20) a % residue(:)

        case(7)
            write(*,30) a % copy_No(:)

        case(8)
            write(*,80) a % BasisPointer(:)

        case(9)
            write(*,10) a % list_of_fragments(:)

        case(10)
            write(*,40) a % list_of_residues(:)

        case(11)
            write(*,60) a % polar(:)

        case(12)
            write(*,60) a % k_WH(:)

        case(13)
            write(*,70) a % DPF(:)

        case(14)
            write(*,70) a % El(:)

        case(15)
            write(*,70) a % Hl(:)

        case(16)
            write(*,70) a % flex(:)

        case(17)
            write(*,70) a % solute(:)

        case(18)
            write(*,20) a % QMMM(:)

        case(19)
            write(*,'(a14,i5)') "N. of atoms = ", a % atoms 

        case(20)
            write(*,'(a18,i5)') "N. of electrons = ", a % N_of_electrons

        case(21)
            write(*,'(a26,i5)') "N. of Solvent Molecules = ", a % N_of_solvent_Molecules

        case(22)
            write(*,'(a25,i5)') "N. of Solute Molecules = ", a % N_of_solute_Molecules

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
80 Format(25i6)

end subroutine debug_EH_structure
!
!
!
!===============================
 subroutine debug_EH_atomic( a )
!===============================
implicit none
type(atomic) , intent(in) :: a(:)

! local variables ...
integer :: option

do

    print'("")'      
    write(*,*) ' (0) QUIT        '
    print'("")'      
    write(*,*) ' (1) Symbol      '
    write(*,*) ' (2) MMSymbol    '
    write(*,*) ' (3) AtNo        '
    write(*,*) ' (4) my_id       '
    write(*,*) ' (5) copy_No     '
    write(*,*) ' (6) fragment    '
    write(*,*) ' (7) nr          '
    write(*,*) ' (8) residue     '
    write(*,*) ' (9) charge      '
    write(*,*) ' (10) DPF        '
    write(*,*) ' (11) El         '
    write(*,*) ' (12) Hl         '
    write(*,*) ' (13) flex       '
    write(*,*) ' (14) QM-MM      '
    write(*,*) ' any other number continues     '

    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,10) a(:) % Symbol

        case(2)
            write(*,20) a(:) % MMSymbol

        case(3)
            write(*,30) a(:) % AtNo

        case(4)
            write(*,80) a(:) % my_id

        case(5)
            write(*,30) a(:) % copy_No

        case(6)
            write(*,10) a(:) % fragment

        case(7)
            write(*,50) a(:) % nr

        case(8)
            write(*,20) a(:) % residue

        case(9)
            write(*,90) a(:) % charge

        case(10)
            write(*,70) a(:) % DPF

        case(11)
            write(*,70) a(:) % El

        case(12)
            write(*,70) a(:) % Hl
            
        case(13)
            write(*,70) a(:) % flex

        case(14)
            write(*,20) a(:) % QMMM
            
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
80 Format(25i6)
90 Format(15F10.3)

end subroutine debug_EH_atomic
!
!
!
!==================================
 subroutine debug_EH_STO_basis( a )
!==================================
implicit none
type(STO_basis) , intent(in) :: a(:)

! local variables ...
integer :: option

do

    print'("")'      
    write(*,*) ' (0) QUIT            '
    print'("")'      
    write(*,*) ' (1)  Symbol         '
    write(*,*) ' (2)  EHSymbol       '
    write(*,*) ' (3)  AtNo           '
    write(*,*) ' (4)  fragment       '
    write(*,*) ' (5)  nr             '
    write(*,*) ' (6)  residue        '
    write(*,*) ' (7)  copy_No        '
    write(*,*) ' (8)  indx           '
    write(*,*) ' (9)  atom           '
    write(*,*) ' (10) Nzeta          '
    write(*,*) ' (11) IP             '
    write(*,*) ' (12) k_WH           '
    write(*,*) ' (13) solute         '
    write(*,*) ' (14) DPF            '
    write(*,*) ' (15) El             '
    write(*,*) ' (16) Hl             '
    write(*,*) ' (17) flex           '
    write(*,*) ' any other number continues     '


    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,10) a(:) % Symbol

        case(2)
            write(*,20) a(:) % EHSymbol

        case(3)
            write(*,30) a(:) % AtNo

        case(4)
            write(*,10) a(:) % fragment

        case(5)
            write(*,50) a(:) % nr

        case(6)
            write(*,20) a(:) % residue

        case(7)
            write(*,30) a(:) % Copy_No

        case(8)
            write(*,80) a(:) % indx

        case(9)
            write(*,80) a(:) % atom

        case(10)
            write(*,30) a(:) % Nzeta

        case(11)
            write(*,90) a(:) % IP

        case(12)
            write(*,60) a(:) % k_WH

        case(13)
            write(*,70) a(:) % solute

        case(14)
            write(*,70) a(:) % DPF

        case(15)
            write(*,70) a(:) % El

        case(16)
            write(*,70) a(:) % Hl

        case(17)
            write(*,70) a(:) % flex

        case default
            exit

    end select

end do

10 Format(50a3)
20 Format(37a4)
30 Format(50i3)
40 Format(30i5)
50 Format(37i4)
60 Format(18F8.3)
70 Format(50L3)
80 Format(25i6)
90 Format(15F10.3)

end subroutine debug_EH_STO_basis

end module type_m
