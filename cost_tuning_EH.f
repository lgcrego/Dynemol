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
chi(1) = ( OPT_UNI%erg(22) - OPT_UNI%erg(21) )  - 5.3780d0           ; weight(1) = 5.0d-1 
chi(2) = ( OPT_UNI%erg(23) - OPT_UNI%erg(21) )  - 6.8730d0           ; weight(2) = 5.0d-1
chi(3) = ( OPT_UNI%erg(23) - OPT_UNI%erg(22) )  - 1.4950d0           ; weight(3) = 5.0d-1
chi(4) = ( OPT_UNI%erg(23) - OPT_UNI%erg(20) )  - 7.3805d0           ; weight(4) = 5.0d-1
chi(5) = ( OPT_UNI%erg(21) - OPT_UNI%erg(20) )  - 0.5075d0           ; weight(5) = 5.0d-1

!-------------------------------------------------------------------------
! Population analysis ...
! Mulliken( GA , basis , MO , atom=[.,.,.] , AO_ang , EHSymbol , residue )
!-------------------------------------------------------------------------
! NO charge in these atoms ...

!chi(49) =  Mulliken(OPT_UNI, basis, MO=20, atom=[ 7])    ; weight(49)= 3.0   

!chi(47) =  Mulliken(OPT_UNI, basis, MO=21, atom=[ 5])    ; weight(47)= 3.0 !2.5    
chi(8)  =  Mulliken(OPT_UNI, basis, MO=21, atom=[ 8])    
chi(9)  =  Mulliken(OPT_UNI, basis, MO=21, atom=[ 2])   - 0.3 
chi(45) =  Mulliken(OPT_UNI, basis, MO=21, atom=[11])   - 0.4
!chi(46) =  Mulliken(OPT_UNI, basis, MO=23, atom=[13]) -0.1

chi(10) =  Mulliken(OPT_UNI, basis, MO=21, atom=[13])    

! missing charge on these atoms ...
chi(13) =  Mulliken(OPT_UNI, basis, MO=22, atom=[11]) - 0.4   
chi(6)  =  Mulliken(OPT_UNI, basis, MO=22, atom=[ 2]) - 0.4   

chi(15) =  Mulliken(OPT_UNI, basis, MO=23, atom=[ 2]) - 0.2 
chi(16) =  Mulliken(OPT_UNI, basis, MO=23, atom=[12]) - 0.2 
chi(14) =  Mulliken(OPT_UNI, basis, MO=23, atom=[11]) - 0.4
chi(14) =  Mulliken(OPT_UNI, basis, MO=23, atom=[ 7]) - 0.4
chi(51) =  Mulliken(OPT_UNI, basis, MO=23, atom=[ 3]) - 0.4   
!-------------------------------------------------------------------------
! MO character ...
! MO_character( system , GA , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!-------------------------------------------------------------------------

chi(50) =  MO_character(OPT_UNI, basis, MO=20, AO='Py') 
chi(17) =  MO_character(OPT_UNI, basis, MO=21, AO='Pz') 
chi(18) =  MO_character(OPT_UNI, basis, MO=22, AO='Pz') 
chi(19) =  MO_character(OPT_UNI, basis, MO=23, AO='Pz') 

!-------------------------------------------------------------------------
! Bond Type analysis ...
! Bond_Type( system , GA , MO , atom1 , atom2 , AO , "+" or "-" )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!-------------------------------------------------------------------------

chi(48) =  Bond_Type(system, OPT_UNI, 21,  5,  7, 'pz', '-')                                 
chi(20) =  Bond_Type(system, OPT_UNI, 21,  5,  3, 'Pz', '+')                                 
chi(21) =  Bond_Type(system, OPT_UNI, 21, 11, 12, 'Pz', '-')                                 
chi(22) =  Bond_Type(system, OPT_UNI, 21, 12,  2, 'Pz', '-')                                 
chi(23) =  Bond_Type(system, OPT_UNI, 21,  3,  2, 'Pz', '-')                                 
chi(24) =  Bond_Type(system, OPT_UNI, 21,  2, 13, 'Pz', '-')                                 
chi(25) =  Bond_Type(system, OPT_UNI, 21, 13, 11, 'Pz', '-')                                 


chi(26) =  Bond_Type(system, OPT_UNI, 22,  5,  7, 'Pz', '+')                                 
chi(27) =  Bond_Type(system, OPT_UNI, 22,  3,  4, 'Pz', '+')                                 
chi(28) =  Bond_Type(system, OPT_UNI, 22, 11, 12, 'Pz', '+')                                 
chi(29) =  Bond_Type(system, OPT_UNI, 22,  5,  3, 'Pz', '-')                                 
chi(30) =  Bond_Type(system, OPT_UNI, 22,  3,  2, 'Pz', '-')                                 
chi(31) =  Bond_Type(system, OPT_UNI, 22,  2, 11, 'Pz', '-')                                 
chi(32) =  Bond_Type(system, OPT_UNI, 22,  7,  8, 'Pz', '-')                                 
chi(33) =  Bond_Type(system, OPT_UNI, 22,  7, 11, 'Pz', '-')                                 


chi(34) =  Bond_Type(system, OPT_UNI, 23,  3,  2, 'Pz', '+')                                 
chi(35) =  Bond_Type(system, OPT_UNI, 23,  7, 11, 'Pz', '+')                                 
chi(36) =  Bond_Type(system, OPT_UNI, 23,  5,  6, 'Pz', '+')                                 
chi(37) =  Bond_Type(system, OPT_UNI, 23, 12, 13, 'Pz', '-')                                 
chi(38) =  Bond_Type(system, OPT_UNI, 23,  2, 12, 'Pz', '-')                                 
chi(39) =  Bond_Type(system, OPT_UNI, 23,  7,  8, 'Pz', '-')                                 
chi(40) =  Bond_Type(system, OPT_UNI, 23,  7,  5, 'Pz', '-')                                
chi(41) =  Bond_Type(system, OPT_UNI, 23,  5,  3, 'Pz', '-')                                
chi(42) =  Bond_Type(system, OPT_UNI, 23, 12, 11, 'Pz', '-')                                

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
