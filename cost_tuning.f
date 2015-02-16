module cost_tuning_m

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken


    public :: evaluate_cost , list_of_modes

    private 

    interface evaluate_cost
        module procedure evaluate_cost_EHT
        module procedure evaluate_cost_nmd
    end interface

    ! module variables ...
    integer  , allocatable :: list_of_modes(:)
contains
!
!
!
!=============================================================
 function evaluate_cost_EHT( OPT_UNI , basis , DP , Alpha_ii )
!=============================================================
implicit none
type(R_eigen)               , intent(in)  :: OPT_UNI
type(STO_basis)             , intent(in)  :: basis(:)
real*8          , optional  , intent(in)  :: DP(3)
real*8          , optional  , intent(in)  :: Alpha_ii(3)
real*8                                    :: evaluate_cost_EHT

! local variables ...
real*8   :: chi(20) , weight(20)
real*8   :: REF_DP(3) , REF_Alpha(3)

! general definitions ...
chi(:) = 0.d0   ;   weight(:) = 0.d0

!--------------------
! HOMO-LUMO gaps ...     
!--------------------

chi(1) = ( OPT_UNI%erg(5) - OPT_UNI%erg(4) ) - 7.6d0                             ; weight(1) = 1.0d0

!--------------------------------------------------------------------
! Population analysis ...
! Mulliken( GA , basis , MO , atom , AO_ang , EHSymbol , residue )
!--------------------------------------------------------------------

!chi(6) =  Mulliken(OPT_UNI,basis,MO=75,residue="TRI") - 0.88d0               ; weight(6) =  2.0d0
!chi(7) =  Mulliken(OPT_UNI,basis,MO=75,residue="TPH") - 0.15d0               ; weight(7) =  2.0d0
!chi(8) =  Mulliken(OPT_UNI,basis,MO=75,residue="CBX") - 0.15d0               ; weight(8) =  2.0d0

!chi(9) =  Mulliken(OPT_UNI,basis,MO=76,residue="TRI") - 0.15d0               ; weight(9) =  2.0d0
!chi(10)=  Mulliken(OPT_UNI,basis,MO=76,residue="TPH") - 0.40d0               ; weight(10)=  2.0d0
!chi(11)=  Mulliken(OPT_UNI,basis,MO=76,residue="CBX") - 0.55d0               ; weight(11)=  2.0d0

!chi(12)=  Mulliken(OPT_UNI,basis,MO=77,residue="TRI") - 0.45d0               ; weight(12)=  12.0d0
!chi(13)=  Mulliken(OPT_UNI,basis,MO=77,residue="TPH") - 0.30d0               ; weight(13)=  12.0d0
!chi(14)=  Mulliken(OPT_UNI,basis,MO=77,residue="CBX") - 0.35d0               ; weight(14)=  12.0d0

!-------------------------
! Total DIPOLE moment ...
!-------------------------

REF_DP = [ 0.d-4 , 1.85d0 , 0.0000d0 ]

chi(6)  = DP(1) - REF_DP(1)     ; weight(6) = 1.d0
chi(7)  = DP(2) - REF_DP(2)     ; weight(7) = 2.d0
chi(8)  = DP(3) - REF_DP(3)     ; weight(8) = 1.d0

!chi(15)  = dot_product( DP , DP ) - dot_product( REF_DP , REF_DP )     ; weight(15) = 2.d0

!-----------------------------------------------------
! Polarizability: Alpha tensor diagonal elements  ...
!-----------------------------------------------------

REF_Alpha = [ 9.2d0 , 8.5d0 , 7.8d0 ]

chi(9)  = Alpha_ii(1) - REF_Alpha(1)     ; weight(9)  = 1.4d0
chi(10) = Alpha_ii(2) - REF_Alpha(2)     ; weight(10) = 1.d0
chi(11) = Alpha_ii(3) - REF_Alpha(3)     ; weight(11) = 1.4d0

!......................................................................
! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost_EHT = sqrt( dot_product(chi,chi) )

end function evaluate_cost_EHT
!
!
!
!
!=============================================================
 function evaluate_cost_nmd( Hesse_erg , nmd_list , instance )
!=============================================================
implicit none
real*8           , optional , intent(in)  :: Hesse_erg(:)
integer          , optional , intent(in)  :: nmd_list(:)
character(*)     , optional , intent(in)  :: instance
real*8           :: evaluate_cost_nmd

! local variables ...
real*8           :: chi(20)    = D_zero
real*8           :: weight(20) = D_zero

select case (instance)
    case("preprocess")
            
        allocate( list_of_modes , source = [10,13,16,17] )

    case default

!------------------------
! NMD frequencies (cm-1})
!------------------------

chi(1) = Hesse_erg(nmd_list(1))  - 113.d0                             ; weight(1) = 1.0d0

chi(2) = Hesse_erg(nmd_list(2))  - 255.d0                             ; weight(2) = 1.0d0

chi(3) = Hesse_erg(nmd_list(3))  - 289.d0                             ; weight(3) = 0.9d0

chi(4) = Hesse_erg(nmd_list(4))  - 528.d0                             ; weight(4) = 0.7d0

!--------------------------------------------------------------------

! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost_nmd = sqrt( dot_product(chi,chi) )

end select

end function evaluate_cost_nmd
!
!
!
end module cost_tuning_m
