module tuning_m

    use type_m
    use constants_m
    use parameters_m    , only  : T_ , F_ , static

    public :: ad_hoc_tuning

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

! edit structure  ...

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      univ % atom(:) % etc
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------
!      define %atom
!-----------------------------------

!-----------------------------------
!      define %residue
!-----------------------------------

!-----------------------------------
!      define %nr
!-----------------------------------

!---------------------------------------------------
!      define %DPF     (Dipole Fragment) 
!      define %V_shift (FMO offset shift)
!---------------------------------------------------
where( univ% atom% residue == "ADN" .OR. univ% atom% residue == "THY") univ% atom% V_shift = -0.3d0
where( univ% atom% residue == "GUA" .OR. univ% atom% residue == "CYT") univ% atom% V_shift = +0.3d0

!---------------------------------------------------
!      define %QMMM  
!      default is QMMM = QM; set QMMM = MM for classical atoms ... 
!---------------------------------------------------

!---------------------------------------------------
!      define %El   : mandatory !!
!---------------------------------------------------
univ % atom (138:151) % El = .true.
!---------------------------------------------------
!      define %Hl   : must be T_ for El/Hl calcs ...
!---------------------------------------------------
univ % atom (138:151) % Hl = .true.
!----------------------------------------------------
!      define %fragment 
!----------------------------------------------------
where( univ% atom% residue == "BKB" ) univ% atom% fragment = "Y"
univ % atom (359:372) % fragment = "1"
univ % atom (170:183) % fragment = "2"
univ % atom (327:340) % fragment = "3"
univ % atom (391:404) % fragment = "4"
univ % atom (106:119) % fragment = "5"
univ % atom (423:436) % fragment = "6"
univ % atom ( 74: 87) % fragment = "7"

!......................................................................
!default: %El => DONOR
If( any(univ % atom%El) ) then
    where( univ % atom % El ) univ % atom % fragment = "D"
else
    if(.NOT. static) stop ">> execution stopped, must define eletron ...%El in ad_hoc_tuning; is ad_hoc = T_? <<"
end If

If( ad_hoc_verbose_ ) then
    Print 46
    ad_hoc_verbose_ = F_
end If

! just touching univ ...
univ = univ

include 'formats.h'

end subroutine ad_hoc_tuning

end module tuning_m
!
!
!
!
module MM_tuning_routines

    use constants_m     , only: large
    use parameters_m    , only: T_ , F_ , static
    use MM_types        , only: MM_atomic, DefineBonds, DefineAngles

    private

    public :: ad_hoc_MM_tuning , SpecialBonds, SpecialAngs

    ! module variables ...
    type(DefineBonds) , allocatable :: SpecialBonds(:)
    type(DefineAngles), allocatable :: SpecialAngs(:)

    contains

!================================================
 subroutine ad_hoc_MM_tuning( atom , instance )
!================================================
implicit none
type(MM_atomic) , optional , intent(inout) :: atom(:)
character(*)               , intent(in)    :: instance

select case ( instance ) 

!==================================
    case ("General")

!----------------------------------
!      define SPECIAL atoms 
!----------------------------------

!----------------------------------
!      define flex
!----------------------------------

!----------------------------------
!      define MM atom types 
!----------------------------------

!----------------------------------
!      define my_species
!----------------------------------

!----------------------------------
!      define residues
!----------------------------------

!----------------------------------
!      define nr
!----------------------------------

!----------------------------------
!       charge of the atoms 
!----------------------------------


    case ("MegaMass")
!----------------------------------
!     Selective_Dynamics 
!----------------------------------


    case( 'SpecialBonds' )
!----------------------------------
!      define SPECIAL bonds
!----------------------------------

!----------------------------------
!      define SPECIAL angles 
!----------------------------------

!=====================================

end select

end subroutine ad_hoc_MM_tuning

end module MM_tuning_routines
!
!
!
module syst
real*8 :: bath_T, press, talt, talp, Initial_density
end module syst
!
!
!
module for_force
 integer                               :: forcefield
 real*8, dimension(:,:)  , allocatable :: vscut, fscut
 real*8                                :: rcut, vrecut, frecut, rcutsq, KAPPA
 real*8                                :: ecoul, eintra, evdw
 real*8                                :: bdpot, harm_bond, morse_bond, Morspot, angpot
 real*8                                :: dihpot, proper_dih, ryck_dih, harm_dih, imp_dih
 real*8                                :: LJ_14, LJ_intra, Coul_14, Coul_intra
 real*8                                :: pot_INTER, pot_total
 character(4)                          :: Dihedral_Potential_Type
end module for_force
!
!
!
module atomicmass
    real*8 , parameter  , dimension(1:107)  :: Atomic_mass = (/                                 &
    1.00795d0,   4.00260d0,   6.94122d0,   9.01218d0,  10.81172d0,  12.01078d0,  14.00672d0,    &
   15.99943d0,  18.99840d0,  20.17976d0,  22.98970d0,  24.30506d0,  26.98153d0,  28.08553d0,    &
   30.97376d0,  32.06552d0,  35.45322d0,  39.94812d0,  39.09830d0,  40.07842d0,  44.95591d0,    &
   47.86710d0,  50.94151d0,  51.99616d0,  54.93805d0,  55.84500d0,  58.93320d0,  58.69344d0,    &
   63.54600d0,  65.38200d0,  69.72310d0,  72.64100d0,  74.92160d0,  78.96340d0,  79.90410d0,    &
   83.79822d0,  85.46783d0,  87.62120d0,  88.90585d0,  91.22422d0,  92.90638d0,  95.96220d0,    &
   98.94000d0, 101.07220d0, 102.90550d0, 106.42120d0, 107.86820d0, 112.41182d0, 114.81830d0,    &
  118.71000d0, 121.76000d0, 127.60320d0, 126.90477d0, 131.29362d0, 132.90540d0, 137.32770d0,    &
  138.90548d0, 140.11612d0, 140.90765d0, 144.24232d0, 146.91510d0, 150.36220d0, 151.96412d0,    &
  157.25320d0, 158.92535d0, 162.50012d0, 164.93032d0, 167.25932d0, 168.93421d0, 173.05452d0,    &
  174.96681d0, 178.49200d0, 180.94788d0, 183.84000d0, 186.20710d0, 190.23320d0, 192.21730d0,    &
  195.08490d0, 196.96654d0, 200.59000d0, 204.38332d0, 207.20000d0, 208.98040d0, 208.98400d0,    &
  209.98710d0, 222.01760d0, 223.01970d0, 226.02540d0, 227.02780d0, 232.03806d0, 231.03588d0,    &
  238.02891d0, 237.00000d0, 244.00000d0, 243.00000d0, 247.00000d0, 247.00000d0, 251.00000d0,    &
  252.00000d0, 257.00000d0, 256.00000d0, 254.00000d0, 257.00000d0,  13.01900d0,  14.02700d0,    &
   15.03500d0,  15.03500d0  /)

    character*2 , dimension(1:107)  :: aicon = (/                     &
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
