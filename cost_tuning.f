module cost_tuning_m

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken


    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG

    private 

    interface evaluate_cost
        module procedure evaluate_cost_EHT
        module procedure evaluate_cost_nmd
    end interface

    ! module variables ...
    real*8 , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:)
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
 function evaluate_cost_nmd( Hesse_erg , nmd_indx , instance )
!=============================================================
implicit none
real*8           , optional               , intent(in)     :: Hesse_erg(:)
integer          , optional , allocatable , intent(inout)  :: nmd_indx(:)
character(*)     , optional               , intent(in)     :: instance
real*8           :: evaluate_cost_nmd 

! local variables ...
real*8           :: chi(30)    = D_zero
real*8           :: weight(30) = D_zero

select case (instance)
    case("preprocess")
    evaluate_cost_nmd = real_large
            
allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,24,15,16,17,18,23,21,22,25,30,26,27,29,28] )
!allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,24,15,16,17,18,23,21,22,25,30,26,27,29,28,31] )

    case default

!------------------------
! NMD frequencies (cm-1})
!------------------------

chi(1) = Hesse_erg(nmd_indx(1))  - 395.d0                             ; weight(1) = 1.0d0

chi(2) = Hesse_erg(nmd_indx(2))  - 395.d0                             ; weight(2) = 1.0d0

chi(3) = Hesse_erg(nmd_indx(3))  - 599.d0                             ; weight(3) = 2.0d0

chi(4) = Hesse_erg(nmd_indx(4))  - 599.d0                             ; weight(4) = 1.0d0

chi(5) = Hesse_erg(nmd_indx(5))  - 677.d0                             ; weight(5) = 1.0d0

chi(6) = Hesse_erg(nmd_indx(6))  - 700.d0                             ; weight(6) = 1.0d0

chi(7) = Hesse_erg(nmd_indx(7))  - 838.d0                             ; weight(7) = 1.0d0

chi(8) = Hesse_erg(nmd_indx(8))  - 838.d0                             ; weight(8) = 1.0d0

chi(9) = Hesse_erg(nmd_indx(9))  - 954.d0                             ; weight(9) = 1.0d0

chi(10)= Hesse_erg(nmd_indx(10))  - 954.d0                            ; weight(10)= 1.0d0

chi(11)= Hesse_erg(nmd_indx(11))  - 981.d0                            ; weight(11)= 1.0d0

chi(12)= Hesse_erg(nmd_indx(12))  - 994.d0                            ; weight(12)= 1.0d0

chi(13)= Hesse_erg(nmd_indx(13))  -  1014.d0                          ; weight(13)= 2.0d0

chi(14)= Hesse_erg(nmd_indx(14))  -  1048.d0                          ; weight(14)= 2.0d0

chi(15)= Hesse_erg(nmd_indx(15))  -  1048.d0                          ; weight(15)= 2.0d0

chi(16)= Hesse_erg(nmd_indx(16))  - 1137.d0                           ; weight(16)= 1.0d0

chi(17)= Hesse_erg(nmd_indx(17))  - 1162.d0                           ; weight(17)= 1.0d0

chi(18)= Hesse_erg(nmd_indx(18))  - 1162.d0                           ; weight(18)= 1.0d0

chi(19)= Hesse_erg(nmd_indx(19))  - 1332.d0                           ; weight(19)= 1.0d0

chi(20)= Hesse_erg(nmd_indx(20))  - 1376.d0                           ; weight(20)= 1.0d0

chi(21)= Hesse_erg(nmd_indx(21))  - 1477.d0                           ; weight(21)= 1.0d0

chi(22)= Hesse_erg(nmd_indx(22))  - 1477.d0                           ; weight(22)= 1.0d0

chi(23)= Hesse_erg(nmd_indx(23))  - 1612.d0                           ; weight(23)= 1.0d0

chi(24)= Hesse_erg(nmd_indx(24))  - 1612.d0                           ; weight(24)= 1.0d0

!chi(25)= Hesse_erg(nmd_indx(25))  - 3080.d0                           ; weight(25)= 1.0d0

!--------------------------------------------------------------------

If( .NOT. allocated(nmd_REF_erg ) ) allocate( nmd_REF_erg  , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx)) )
If( .NOT. allocated(nmd_NOPT_erg) ) allocate( nmd_NOPT_erg , source = Hesse_erg(nmd_indx) )

! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost_nmd = sqrt( dot_product(chi,chi) )

end select

end function evaluate_cost_nmd
!
!
!
end module cost_tuning_m
