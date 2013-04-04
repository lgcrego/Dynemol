module tuning_m

    use type_m
    use constants_m
    use parameters_m    , only  : T_ , F_

    public :: Setting_Fragments , ad_hoc_tuning

    private

    logical , save  :: ad_hoc_verbose_ = T_

    contains
!
!
!
!================================
 subroutine ad_hoc_tuning( univ )
!================================
implicit none
type(universe) , intent(inout) :: univ

! local variables ...
integer :: i , ioerr

! edit structure  ...

!-----------------------------------
!      define %atom
!-----------------------------------

 where( univ % atom % AtNo == 1 ) univ % atom % hardcore = 2.04d0
 where( univ % atom % AtNo == 6 ) univ % atom % hardcore = 2.45d0
 where( univ % atom % AtNo == 7 ) univ % atom % hardcore = 2.04d0
 where( univ % atom % AtNo == 8 ) univ % atom % hardcore = 1.75d0

!-----------------------------------
!      define %residue
!-----------------------------------

!-----------------------------------
!      define %nr
!-----------------------------------
 
!------------------------------------
!      define %DPF (Dipole Fragment) 
!------------------------------------

!default: %DPF = F_

!use default: %DPF = %solute  
! where( univ % atom % DPF ) univ % atom % solute = .true.

!-----------------------------------
!      define %El   : mandatory !!
!-----------------------------------

 where( univ % atom % residue == "BE2" ) univ % atom % El = .true.

!---------------------------------------------------
!      define %Hl   : must be T_ for El/Hl calcs ...
!---------------------------------------------------

 where( univ % atom % residue == "BE1" ) univ % atom % Hl = .true.

!------------------------------------------------
!      define %fragments   : Donor fragment ...
!------------------------------------------------

!default: %El => DONOR
 where( univ % atom % El ) univ % atom % fragment = "D"

!......................................................................

If( ad_hoc_verbose_ ) then
    Print 46
    ad_hoc_verbose_ = F_
end If

include 'formats.h'

end subroutine ad_hoc_tuning
!
!
!
!=================================
 subroutine Setting_Fragments( a )
!=================================
implicit none
type(universe)  , intent(inout) :: a

! local variables ...
integer  :: i 

! --------- Table of STANDARD fragments ----------------
!
! STANDARD fragments are set based on RESIDUE names ...
! 
!   Acceptor    =   A       
!   Bridge      =   B
!   Donor       =   D  (defined ONLY in ad_hoc)
!   Electron    =   E  (defined ONLY in ad_hoc)
!   Hole        =   H 
!   Molecule    =   M
!   Solvent     =   S
!   Cluster     =   C 
!   System      =   #
!
! some typical cases are used below ...
!--------------------------------------------------------

 DO i = 1 , size(a%atom)
 
    select case(a%atom(i)%residue)

        case( 'BE1') 
            a%atom(i)%fragment = '1' 

        case( 'ETH') 
            a%atom(i)%fragment = 'B' 

        case( 'H2O' , 'SOL' ) 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 2.0d0
        
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 3.d0

        case( 'PPH') 
            a%atom(i)%solvation_hardcore = 2.d0

        case( 'C60') 
            a%atom(i)%fragment = '0' 
            a%atom(i)%solvation_hardcore = 2.d0

        case( 'LIG') 
            a%atom(i)%fragment = '1' 
            a%atom(i)%solvation_hardcore = 2.d0

    end select

 END DO

end subroutine Setting_Fragments

end module tuning_m 
!
!
!
!
module MM_tuning_routines

use MM_types , only : MM_atomic, DefineBonds, DefineAngles

private

public :: ad_hoc_MM_tuning , SpecialBonds, SpecialAngs

! module variables ...
type(DefineBonds) , allocatable :: SpecialBonds(:)
type(DefineAngles), allocatable :: SpecialAngs(:)

contains

!================================================
 subroutine ad_hoc_MM_tuning( system , instance )
!================================================
implicit none
type(MM_atomic) , optional , intent(inout) :: system(:)
character(*)               , intent(in)    :: instance


select case ( instance ) 

!=================================================
    case ("General")

!----------------------------------
!      define SPECIAL atoms 
!----------------------------------

!----------------------------------
!      define MM atom types 
!----------------------------------

!----------------------------------
!      define my_species
!----------------------------------

!----------------------------------
!      define resid's
!----------------------------------

!----------------------------------
!      define nresid's
!----------------------------------


!----------------------------------
!     Selective_Dynamics
!----------------------------------

where(system % my_id == 50 ) system % free = .false.
where(system % my_id == 52 ) system % free = .false.
where(system % my_id ==  2 ) system % free = .false.
where(system % my_id == 51 ) system % free = .false.

!----------------------------------
!       charge of the atoms 
!----------------------------------


!=================================================

    case( 'SpecialBonds' )

!----------------------------------
!      define SPECIAL bonds
!----------------------------------
 allocate(SpecialBonds(2))
 SpecialBonds(1) % nome      = 'bond_gb16'
 SpecialBonds(1) % kbond0(1) = 392459.2
 SpecialBonds(1) % kbond0(2) = 0.14010

 SpecialBonds(2) % nome      = 'bond_gb53'
 SpecialBonds(2) % kbond0(1) = 392459.2
 SpecialBonds(2) % kbond0(2) = 0.14580
!----------------------------------
!      define SPECIAL angles 
!----------------------------------
 allocate(SpecialAngs(2))
 SpecialAngs(1) % nome     = 'angle_ga07'
 SpecialAngs(1) % kang0(1) = 527.184
 SpecialAngs(1) % kang0(2) = 108.000

 SpecialAngs(2) % nome     = 'angle_ga27'
 SpecialAngs(2) % kang0(1) = 527.184
 SpecialAngs(2) % kang0(2) = 120.000

!=================================================

end select

end subroutine ad_hoc_MM_tuning

end module MM_tuning_routines
!
!
!
module syst
 real*8                                 :: temper, press, talt, talp, Initial_density
 real*8                                 :: Ekin     = 0.d0
 real*8                                 :: DensTot  = 0.d0
 real*8                                 :: TempTot  = 0.d0
 real*8                                 :: PressTot = 0.d0
end module syst
!
!
!
module for_force
 integer                               :: forcefield
 real*8                                :: rcut, vrecut, frecut, rcutsq, pot, ecoul, eintra, evdw, bdpot, angpot, dihpot
 real*8, dimension(:,:)  , allocatable :: vscut, fscut
 real*8, dimension(:,:,:), allocatable :: tmp_fsr, tmp_fch
 real*8, dimension(:,:)  , allocatable :: erfkr
 real*8                                :: KAPPA, vself, lj14pot, coul14pot, pot2
 character(4)                          :: Dihedral_Potential_Type
end module for_force
!
!
!
module atomicmass
 real*8, dimension(1:107) :: atmas = (/                               &
    1.008,   4.003,   6.939,   9.012,  10.811,  12.011,  14.007,      &
   15.999,  18.998,  20.183,  22.989,  24.312,  26.982,  28.086,      &
   30.974,  32.064,  35.453,  39.948,  39.102,  40.080,  44.956,      &
   47.900,  50.942,  51.996,  54.938,  55.847,  58.933,  58.710,      &
   63.540,  65.370,  69.720,  72.590,  74.922,  78.960,  79.909,      &
   83.800,  85.470,  87.620,  88.905,  91.220,  92.906,  95.940,      &
   98.000, 101.070, 102.905, 106.400, 107.870, 112.400, 114.820,      &
  118.690, 121.750, 127.600, 126.904, 131.300, 132.905, 137.340,      &
  138.910, 140.120, 140.907, 144.240, 147.000, 150.350, 151.960,      &
  157.250, 158.924, 162.500, 164.930, 167.260, 168.934, 173.040,      &
  174.970, 178.490, 180.948, 183.850, 186.200, 190.200, 192.200,      &
  195.090, 196.967, 200.590, 204.370, 207.190, 208.980, 210.000,      &
  210.000, 222.000, 223.000, 226.000, 227.000, 232.038, 231.000,      &
  238.030, 237.000, 242.000, 243.000, 247.000, 247.000, 249.000,      &
  254.000, 253.000, 256.000, 254.000, 257.000,  13.019,  14.027,      &
   15.035,  15.035  /)

 character*2, dimension(1:107) :: aicon = (/                          &
     ' H',    'HE',    'LI',    'BE',    ' B',    ' C',    ' N',      &
     ' O',    ' F',    'NE',    'NA',    'MG',    'AL',    'SI',      &
     ' P',    ' S',    'CL',    'AR',    ' K',    'CA',    'SC',      &
     'TI',    ' V',    'CR',    'MN',    'FE',    'CO',    'NI',      &
     'CU',    'ZN',    'GA',    'GE',    'AS',    'SE',    'BR',      &
     'KR',    'RB',    'SR',    ' Y',    'ZR',    'NB',    'MO',      &
     'TC',    'RU',    'RH',    'PD',    'AG',    'CD',    'IN',      &
     'SN',    'SB',    'TE',    ' I',    'XE',    'CS',    'BA',      &
     'LA',    'CE',    'PR',    'ND',    'PM',    'SM',    'EU',      &
     'GD',    'TB',    'DY',    'HO',    'ER',    'TM',    'YB',      &
     'LU',    'HF',    'TA',    ' W',    'RE',    'OS',    'IR',      &
     'PT',    'AU',    'HG',    'TL',    'PB',    'BI',    'PO',      &
     'AT',    'RN',    'FR',    'RA',    'AC',    'TH',    'PA',      &
     ' U',    'NP',    'PU',    'AM',    'CM',    'BK',    'CF',      &
     'ES',    'FM',    'MD',    'NO',    'LW',    ' C',    ' C',      &
     ' C',    ' C'  /)
end module atomicmass
