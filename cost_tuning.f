module cost_tuning_m

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken


    public :: evaluate_cost

    private 

contains
!
!
!
!==============================================
 function evaluate_cost( OPT_UNI , basis , DP )
!==============================================
implicit none
type(R_eigen)   , intent(in)  :: OPT_UNI
type(STO_basis) , intent(in)  :: basis(:)
real*8          , intent(in)  :: DP(3)
real*8                        :: evaluate_cost

! local variables ...
real*8   :: chi(20) , weight(20)
real*8   :: REF_DP(3)

! general definitions ...
chi(:) = 0.d0   ;   weight(:) = 0.d0

!-----------------------------------
! HOMO-LUMO gaps ...     
!-----------------------------------

chi(1) = ( OPT_UNI%erg(16) - OPT_UNI%erg(15) ) - 8.8425d0                             ; weight(1) = 1.0d0
chi(2) = ( OPT_UNI%erg(17) - OPT_UNI%erg(15) ) - 9.5315d0                             ; weight(2) = 1.0d0
chi(3) = ( OPT_UNI%erg(17) - OPT_UNI%erg(16) ) - 0.6890d0                             ; weight(3) = 1.0d0
chi(4) = ( OPT_UNI%erg(15) - OPT_UNI%erg(14) ) - 1.8800d0                             ; weight(4) = 1.0d0
chi(5) = ( OPT_UNI%erg(16) - OPT_UNI%erg(14) ) - 10.7300d0                             ; weight(5) = 1.0d0

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

!-----------------------------------
! Total DIPOLE moment ...
!-----------------------------------

REF_DP = [ 3.d-4 , 1.63d0 , 0.0015d0 ]

chi(6)  = DP(1) - REF_DP(1)     ; weight(6) = 2.d0
chi(7)  = DP(2) - REF_DP(2)     ; weight(7) = 2.d0
chi(8)  = DP(3) - REF_DP(3)     ; weight(8) = 2.d0

!chi(15)  = dot_product( DP , DP ) - dot_product( REF_DP , REF_DP )     ; weight(15) = 2.d0

!......................................................................
! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost = sqrt( dot_product(chi,chi) )

end function evaluate_cost
!
!
!
end module cost_tuning_m
