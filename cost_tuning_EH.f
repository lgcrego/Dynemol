module cost_EH

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken


    public :: evaluate_cost , REF_DP , REF_Alpha

    ! module variables ...
    real*8 :: REF_DP(3) , REF_Alpha(3)

    private 

contains
!
!
!
!=========================================================
 function evaluate_cost( OPT_UNI , basis , DP , Alpha_ii )
!=========================================================
implicit none
type(R_eigen)               , intent(in)  :: OPT_UNI
type(STO_basis)             , intent(in)  :: basis(:)
real*8          , optional  , intent(in)  :: DP(3)
real*8          , optional  , intent(in)  :: Alpha_ii(3)
real*8                                    :: evaluate_cost

! local variables ...
integer  :: dumb
real*8   :: chi(20) , weight(20)

! general definitions ...
chi(:) = 0.d0   ;   weight(:) = 0.d0

!--------------------
! HOMO-LUMO gaps ...     
!--------------------

chi(1) = ( OPT_UNI%erg(121) - OPT_UNI%erg(120) )  - 2.016d0                           ; weight(1) = 2.0d0
chi(2) = ( OPT_UNI%erg(120) - OPT_UNI%erg(119) )  - 0.842d0                           ; weight(2) = 2.0d0
chi(3) = ( OPT_UNI%erg(122) - OPT_UNI%erg(121) )  - 1.625d0                           ; weight(3) = 1.5d0
chi(4) = ( OPT_UNI%erg(123) - OPT_UNI%erg(121) )  - 2.359d0                           ; weight(4) = 1.0d0
chi(5) = ( OPT_UNI%erg(123) - OPT_UNI%erg(122) )  - 0.734d0                           ; weight(5) = 1.0d0

chi(6)  = ( OPT_UNI%erg(345) - OPT_UNI%erg(344) )  - 2.000d0                          ; weight(6) = 2.0d0
chi(7)  = ( OPT_UNI%erg(344) - OPT_UNI%erg(343) )  - 0.102d0                          ; weight(7) = 2.0d0
chi(8)  = ( OPT_UNI%erg(346) - OPT_UNI%erg(345) )  - 1.730d0                          ; weight(8) = 1.5d0
chi(9)  = ( OPT_UNI%erg(347) - OPT_UNI%erg(345) )  - 3.429d0                          ; weight(9) = 1.0d0
chi(10) = ( OPT_UNI%erg(347) - OPT_UNI%erg(346) )  - 1.700d0                          ; weight(10) = 1.0d0

!--------------------------------------------------------------------
! Population analysis ...
! Mulliken( GA , basis , MO , atom , AO_ang , EHSymbol , residue )
!--------------------------------------------------------------------

!chi(6) =  Mulliken(OPT_UNI,basis,MO=75,residue="TRI") - 0.88d0               ; weight(6) =  2.0d0
!chi(7) =  Mulliken(OPT_UNI,basis,MO=75,residue="TPH") - 0.15d0               ; weight(7) =  2.0d0
!chi(8) =  Mulliken(OPT_UNI,basis,MO=75,residue="CBX") - 0.15d0               ; weight(8) =  2.0d0

!-------------------------
! Total DIPOLE moment ...
!-------------------------

REF_DP = [ -8.45 , -2.92 , -5.325 ]

!chi(11)  = DP(1) - REF_DP(1)     ; weight(11) = 1.d-1
!chi(12)  = DP(2) - REF_DP(2)     ; weight(12) = 1.d-1
!chi(13)  = DP(3) - REF_DP(3)     ; weight(13) = 1.d-1

!-----------------------------------------------------
! Polarizability: Alpha tensor diagonal elements  ...
!-----------------------------------------------------

REF_Alpha = [ 9.2d0 , 8.5d0 , 7.8d0 ]

!chi(9)  = Alpha_ii(1) - REF_Alpha(1)     ; weight(9)  = 1.4d0
!chi(10) = Alpha_ii(2) - REF_Alpha(2)     ; weight(10) = 1.d0
!chi(11) = Alpha_ii(3) - REF_Alpha(3)     ; weight(11) = 1.4d0

!......................................................................
! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost = sqrt( dot_product(chi,chi) )

! just touching basis ...
dumb = basis(1)%indx

end function evaluate_cost
!
!
!
end module cost_EH
