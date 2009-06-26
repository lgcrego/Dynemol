module type_m

    use constants_m 

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

    type structure
        integer                                         :: atoms 
        integer                                         :: N_of_electrons
        integer                                         :: molecule
        integer          , dimension(:)   , allocatable :: copy_No
        integer          , dimension(:)   , allocatable :: BasisPointer
        character(len=2) , dimension(:)   , allocatable :: fragment
        character(len=2) , dimension(:)   , allocatable :: symbol
        integer          , dimension(:)   , allocatable :: AtNo
        real*8           , dimension(:,:) , allocatable :: coord
        real*8           , dimension(:)   , allocatable :: k_WH
        real*8           , dimension(3)                 :: BoxSides
        real*8           , dimension(3)                 :: Center_of_Mass
        real*8           , dimension(3)                 :: Center_of_Charge
        real*8           , dimension(3)                 :: T_xyz
    end type structure

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
 

    type real_interval
        real*8 :: inicio 
        real*8 :: fim
    end type real_interval


    type integer_interval
        integer :: inicio
        integer :: fim
    end type integer_interval
   

    type transition
        type(R3_vector)        , allocatable :: matrix(:,:)
        type(real_interval)                  :: bra_range
        type(real_interval)                  :: ket_range
        type(integer_interval)               :: bra_indx_range
        type(integer_interval)               :: ket_indx_range
        integer                , allocatable :: bra_POINTER(:)
        integer                , allocatable :: ket_POINTER(:)
        character(len=8)                     :: flag
    end type transition

    include 'parameters.h'    

end module type_m
