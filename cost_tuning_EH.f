module cost_EH

    use type_m
    use constants_m
    use GA_QCModel_m , only : MO_erg_diff,        &  
                              Mulliken,           &
                              Bond_Type,          &
                              MO_character,       &
                              Localize,           &
                              Exclude,            &
                              Adaptive_GA,        &
                              me => i_       

    public :: evaluate_cost , REF_DP , REF_Alpha 

    ! module variables ...
    real*8 :: REF_DP(3) , REF_Alpha(3)

    ! module parameters ...
    logical :: lock = .false.

    private 

contains
!
!===========================================================================
! 
!                  THE TRUE PARAMETERS ARE OUT THERE
!
!                          TRUST NO OTHERS
!
!                       DENY EVERYTHING ELSE
!
!===========================================================================
!
!
!==========================================================================
 function evaluate_cost( sys , OPT_UNI , basis , DP , Alpha_ii , ShowCost )
!==========================================================================
implicit none
type(structure)             , intent(in) :: sys
type(R_eigen)               , intent(in) :: OPT_UNI
type(STO_basis)             , intent(in) :: basis(:)
real*8          , optional  , intent(in) :: DP(3)
real*8          , optional  , intent(in) :: Alpha_ii(3)
logical         , optional  , intent(in) :: ShowCost
real*8                                   :: evaluate_cost

! local variables ...
integer  :: i , dumb
real*8   :: eval(200) = D_zero
logical  :: input_mode

input_mode = Adaptive_GA% mode

!-------------------------------------------------------------------------
! Energy gaps ...     
! MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , {weight} )
! {...} terms are optional 
!-------------------------------------------------------------------------
eval(me) = MO_erg_diff( OPT_UNI, 41 , 40 , 3.28d0 )
eval(me) = MO_erg_diff( OPT_UNI, 42 , 40 , 5.10d0 )
eval(me) = MO_erg_diff( OPT_UNI, 43 , 40 , 5.10d0 )
eval(me) = MO_erg_diff( OPT_UNI, 41 , 39 , 5.58d0 )
eval(me) = MO_erg_diff( OPT_UNI, 41 , 38 , 5.58d0 )
eval(me) = MO_erg_diff( OPT_UNI, 41 , 37 , 6.03d0 )
!----------------------------------------------------------------------------------------------
! ==> MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!
! ==> Bond_Type( sys , OPT_UNI , MO , atom1 , AO1 , atom2 , AO2 , "+" or "-" )
! Bond Topolgy analysis ...
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
! weight < 0  ==> does not update "me" when Mulliken in called
!
! ==> Exclude( OPT_UNI , basis , MO , {atom}=[:] , {AO} , {EHSymbol} , {residue} , {reference} , {from_to} , {adaptive} )
! NO charge on these atoms ...
! {...} terms are optional  
! default reference < 0.001 
! from_to = real_interval( begin , end ) : no need to use {reference} if {from_to} is used
! adaptive = {input_mode,lock} : logical flag to enable adpative GA method, lock sets reference = end
!
! ==> Localize( OPT_UNI , basis , MO , {atom}=[:] , {AO} , {EHSymbol} , {residue} , {reference} , {from_to} , {adaptive} )
! {...} terms are optional
! default criterium (reference=0.85): localized > 85% of total population
! from_to = real_interval( begin , end ) : no need to use {reference} if {from_to} is used
! adaptive = {input_mode,lock} : logical flag to enable adpative GA method , lock sets reference = end
!----------------------------------------------------------------------------------------------

!39 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NH") - 0.6
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NB") - 0.35
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "CP") - 0.05
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "CA") - 0.05

!38 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "NH") - 0.6	
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "NB") - 0.35
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "CP") - 0.05
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "CA") - 0.05
!
!37 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=37, AO="Py",  EHSymbol = "NB") - 0.4
!eval(me) =  Mulliken(OPT_UNI, basis, MO=37, AO="Px",  EHSymbol = "NB") - 0.4

!36 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Px",  EHSymbol = "NB") - 0.4
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Py",  EHSymbol = "NB") - 0.2
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Px",  EHSymbol = "NA") - 0.2

!35 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Px",  EHSymbol = "NB") - 0.2
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Py",  EHSymbol = "NB") - 0.4
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Py",  EHSymbol = "NA") - 0.2


!!39 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NH", from_to = real_interval( 0.0 , 0.60 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.33 ), adaptive  = input_mode) 
!
!!38 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "NH", from_to = real_interval( 0.0 , 0.60 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.33 ), adaptive  = input_mode) 
!!
!!37 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=37, AO="Py",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.40 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=37, AO="Px",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.40 ), adaptive  = input_mode) 
!
!!36 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Px",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.40 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Py",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.20 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=36, AO="Px",  EHSymbol = "NA", from_to = real_interval( 0.0 , 0.20 ), adaptive  = input_mode) 
!
!!35 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Px",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.20 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Py",  EHSymbol = "NB", from_to = real_interval( 0.0 , 0.40 ), adaptive  = input_mode) 
!eval(me) =  Mulliken(OPT_UNI, basis, MO=35, AO="Py",  EHSymbol = "NA", from_to = real_interval( 0.0 , 0.20 ), adaptive  = input_mode) 

!-------------------------                                                         
! Total DIPOLE moment ...
!-------------------------
!REF_DP = [ 0.0d0 , 0.0d0 , 2.2d0 ]
!eval(me+1) = DP(1) - REF_DP(1)     
!eval(me+2) = DP(2) - REF_DP(2)    
!eval(me+3) = DP(3) - REF_DP(3) 
!me = me + 3

!-----------------------------------------------------
! Polarizability: Alpha tensor diagonal elements  ...
!-----------------------------------------------------
!REF_Alpha = [ 9.2d0 , 8.5d0 , 7.8d0 ]
!eval() = Alpha_ii(1) - REF_Alpha(1)   
!eval() = Alpha_ii(2) - REF_Alpha(2)  
!eval() = Alpha_ii(3) - REF_Alpha(3) 

!......................................................................
! at last, show the cost ...
If( present(ShowCost) ) then

   open( unit=33 , file='opt.trunk/view_cost.dat' , status='unknown' )

   do i = 1 , me
      write(33,*) i , dabs(eval(i)) 
   end do 

   CALL system( dynemoldir//"env.sh save_cost_statement " )

   Print 218

end If
!......................................................................

! evaluate total cost ...
evaluate_cost = sum( abs(eval) )

! just touching variables ...
dumb = basis(1)%atom

!reset index for next round ...
me = 0

include 'formats.h'

end function evaluate_cost
!
!
!
end module cost_EH
