module cost_EH

    use type_m
    use constants_m
    use GA_QCModel_m            , only : Mulliken , Bond_Type , MO_character


    public :: evaluate_cost , REF_DP , REF_Alpha

    ! module variables ...
    real*8 :: REF_DP(3) , REF_Alpha(3)

    private 

contains
!
!
!
!==================================================================
 function evaluate_cost( system , OPT_UNI , basis , DP , Alpha_ii , ShowCost)
!==================================================================
implicit none
type(structure)             , intent(in) :: system
type(R_eigen)               , intent(in) :: OPT_UNI
type(STO_basis)             , intent(in) :: basis(:)
real*8          , optional  , intent(in) :: DP(3)
real*8          , optional  , intent(in) :: Alpha_ii(3)
logical         , optional  , intent(in) :: ShowCost
real*8                                   :: evaluate_cost

! local variables ...
integer  :: i , dumb
real*8   :: chi(55)    = D_zero
real*8   :: weight(55) = D_one
real*8   :: REF_DP(3) , REF_Alpha(3)

!--------------------
! HOMO-LUMO gaps ...     
!--------------------
chi(1) = ( OPT_UNI%erg(115) - OPT_UNI%erg(114) )  - 2.6470d0           ; weight(1) = 1.0d0 
chi(2) = ( OPT_UNI%erg(114) - OPT_UNI%erg(113) )  - 0.3040d0           ; weight(2) = 1.0d0
chi(3) = ( OPT_UNI%erg(115) - OPT_UNI%erg(113) )  - 2.9510d0           ; weight(3) = 1.0d0
chi(4) = ( OPT_UNI%erg(113) - OPT_UNI%erg(112) )  - 0.8950d0           ; weight(4) = 1.0d0
chi(5) = ( OPT_UNI%erg(112) - OPT_UNI%erg(111) )  - 0.4360d0           ; weight(5) = 1.0d0
chi(6) = ( OPT_UNI%erg(117) - OPT_UNI%erg(116) )  - 1.6000d0           ; weight(6) = 1.0d0

!-------------------------------------------------------------------------
! Population analysis ...
! Mulliken( GA , basis , MO , atom=[.,.,.] , AO_ang , EHSymbol , residue )
!-------------------------------------------------------------------------
! NO charge in these atoms ...

chi(11) =  Mulliken(OPT_UNI, basis, MO=112, atom=[33:40])
chi(12) =  Mulliken(OPT_UNI, basis, MO=112, atom=[44:51])
chi(13) =  Mulliken(OPT_UNI, basis, MO=112, atom=[66:73])
chi(14) =  Mulliken(OPT_UNI, basis, MO=112, atom=[55:62])

!-------------------------------------------------------------------------
! MO character ...
! MO_character( system , GA , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!-------------------------------------------------------------------------

!chi(50) =  MO_character(OPT_UNI, basis, MO=20, AO='Py') 

!-------------------------------------------------------------------------
! Bond Type analysis ...
! Bond_Type( system , GA , MO , atom1 , atom2 , AO , "+" or "-" )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!-------------------------------------------------------------------------

chi(55) =  Bond_Type(system, OPT_UNI, 113,  23,  28, 'pz', '+')                                 
chi(54) =  Bond_Type(system, OPT_UNI, 113,  22,  27, 'pz', '+')                                 
chi(53) =  Bond_Type(system, OPT_UNI, 113,  22,  19, 'pz', '+')                                 
chi(52) =  Bond_Type(system, OPT_UNI, 113,   7,  12, 'pz', '+')                                 
chi(51) =  Bond_Type(system, OPT_UNI, 113,  23,  19, 'pz', '-')                                 
chi(50) =  Bond_Type(system, OPT_UNI, 113,  27,  28, 'pz', '-')                                 
chi(49) =  Bond_Type(system, OPT_UNI, 113,  22,   6, 'pz', '-')                                 
chi(48) =  Bond_Type(system, OPT_UNI, 113,  22,  20, 'pz', '-')                                 


chi(47) =  Bond_Type(system, OPT_UNI, 112,  22,   6, 'pz', '+')                                 
chi(46) =  Bond_Type(system, OPT_UNI, 112,  22,  20, 'pz', '+')                                 
chi(41) =  Bond_Type(system, OPT_UNI, 112,  30,  27, 'pz', '+')                                 
chi(42) =  Bond_Type(system, OPT_UNI, 112,  16,  19, 'pz', '+')                                 
chi(45) =  Bond_Type(system, OPT_UNI, 112,  22,  19, 'pz', '+')                                 
chi(44) =  Bond_Type(system, OPT_UNI, 112,  22,  27, 'pz', '+')                                 



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

!show the cost ...
If( present(ShowCost) ) then

   open( unit=33 , file='opt_trunk/view_cost.dat' , status='unknown' )

   do i = 1 , size(chi)
      write(33,*) i , sqrt(chi(i)*chi(i)) , weight(i)
   end do 

end If

! just touching basis ...
dumb = basis(1)%indx

end function evaluate_cost
!
!
!
end module cost_EH
