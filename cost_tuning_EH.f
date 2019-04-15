module cost_EH

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken , Bond_Type


    public :: evaluate_cost , REF_DP , REF_Alpha

    ! module variables ...
    real*8 :: REF_DP(3) , REF_Alpha(3)

    private 

contains
!
!
!
!==================================================================
 function evaluate_cost( system , OPT_UNI , basis , DP , Alpha_ii )
!==================================================================
implicit none
type(structure)             , intent(in) :: system
type(R_eigen)               , intent(in) :: OPT_UNI
type(STO_basis)             , intent(in) :: basis(:)
real*8          , optional  , intent(in) :: DP(3)
real*8          , optional  , intent(in) :: Alpha_ii(3)
real*8                                   :: evaluate_cost

! local variables ...
integer  :: dumb
real*8   :: chi(40)    = D_zero
real*8   :: weight(40) = D_one
real*8   :: REF_DP(3) , REF_Alpha(3)

!--------------------
! HOMO-LUMO gaps ...     
!--------------------
chi(1) = ( OPT_UNI%erg(29) - OPT_UNI%erg(28) )  - 5.5190d0                         ; weight(1) = 1.0d0
chi(2) = ( OPT_UNI%erg(30) - OPT_UNI%erg(28) )  - 5.7242d0                         ; weight(2) = 1.0d0
chi(3) = ( OPT_UNI%erg(30) - OPT_UNI%erg(29) )  - 0.2050d0                         ; weight(3) = 1.0d0
chi(4) = ( OPT_UNI%erg(30) - OPT_UNI%erg(27) )  - 6.9960d0                         ; weight(4) = 1.0d0

!-------------------------------------------------------------------------
! Population analysis ...
! Mulliken( GA , basis , MO , atom=[.,.,.] , AO_ang , EHSymbol , residue )
!-------------------------------------------------------------------------
! NO charge in these atoms ...


! missing charge on these atoms ...
!chi(7)  =  Mulliken(OPT_UNI,basis,MO=29,residue="NH2")   - 1.00d0               ; weight(7)  =  5.0d0
!chi(8)  =  Mulliken(OPT_UNI,basis,MO=29,atom=[8]     )   - 1.00d0               ; weight(8)  =  3.0d0

!-------------------------------------------------------------------------
! Bond Type analysis ...
! Bond_Type( system , GA , MO , atom1 , atom2 , AO , "+" or "-" )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!-------------------------------------------------------------------------
chi(31) =  Bond_Type(system, OPT_UNI, 29, 4 , 5 , 'Pz', '+')                       ! ; weight(31) = 1.0d0          
chi(32) =  Bond_Type(system, OPT_UNI, 29, 10, 11, 'Pz', '+')                       ! ; weight(32) = 1.0d0          
chi(33) =  Bond_Type(system, OPT_UNI, 29, 4 , 11, 'Pz', '-')                       ! ; weight(33) = 1.0d0          
                                                                                               
chi(35) =  Bond_Type(system, OPT_UNI, 30, 8 , 7 , 'Pz', '+')                       ! ; weight(35) = 1.0d0         
chi(36) =  Bond_Type(system, OPT_UNI, 30, 1 , 11, 'Pz', '+')                       ! ; weight(36) = 1.0d0         
                                                                                    
!-------------------------                                                         
! Total DIPOLE moment ...
!-------------------------

REF_DP = [ 0.d-4 , 1.85d0 , 0.0000d0 ]

!chi(6)  = DP(1) - REF_DP(1)     ; weight(6) = 1.d0
!chi(7)  = DP(2) - REF_DP(2)     ; weight(7) = 2.d0
!chi(8)  = DP(3) - REF_DP(3)     ; weight(8) = 1.d0

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
