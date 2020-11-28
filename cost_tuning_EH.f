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

    private 

contains
!
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
real*8   :: REF_DP(3) , REF_Alpha(3)
logical  :: mode

mode = Adaptive_GA% mode

!-------------------------------------------------------------------------
! Energy gaps ...     
! MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , {weight} )
! {...} terms are optional 
!-------------------------------------------------------------------------
eval(me) = MO_erg_diff( OPT_UNI, 115, 114,  3.2d0 )

!----------------------------------------------------------------------------------------------
! ==> MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!
! ==> Localize( OPT_UNI , basis , MO , {atom}=[:] , {EHSymbol} , {residue} , {threshold} , {slide} , {adaptive} )
! {...} terms are optional
! default criterium (threshold=0.85): localized > 85% of total population
! slide = real_interval( begin , end ) : no need to use {threshold} if {slide} is used
! adaptive = mode : to enable adpative GA method
!
! ==> Bond_Type( sys , OPT_UNI , MO , atom1 , AO1 , atom2 , AO2 , "+" or "-" )
! Bond Topolgy analysis ...
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!
! ==> Exclude( OPT_UNI , basis , MO , {atom}=[:] , {EHSymbol} , {residue} , {threshold} , {slide} , {adaptive} )
! NO charge on these atoms ...
! {...} terms are optional  
! default threshold < 0.001 
! slide = real_interval( begin , end ) : no need to use {threshold} if {slide} is used
! adaptive = mode : to enable adpative GA method
!
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
! weight < 0  ==> does not update "me" when Mulliken in called
!----------------------------------------------------------------------------------------------

!114 ===================

!115 ===================
eval(me) = Bond_Type( sys , OPT_UNI , 115 , 4 , 'Pz' , 8 , 'Pz' , '+' )
eval(me) = Bond_Type( sys , OPT_UNI , 115 , 7 , 'Pz' , 5 , 'Pz' , '+' )
eval(me) = Bond_Type( sys , OPT_UNI , 115 , 6 , 'Pz' , 3 , 'Pz' , '+' )

eval(me) = Localize( OPT_UNI , basis , 115 , EHSymbol='SA' , slide=real_interval( 0.20 , 0.26 ), adaptive=mode )
eval(me) = Exclude ( OPT_UNI , basis , 115 , atom = [3,6]  , slide=real_interval( 0.12 , 0.08 ), adaptive=mode )
!-------------------------                                                         
! Total DIPOLE moment ...
!-------------------------
REF_DP = [ 0.d-4 , 1.85d0 , 0.0000d0 ]
!eval()  = DP(1) - REF_DP(1)     
!eval()  = DP(2) - REF_DP(2)    
!eval()  = DP(3) - REF_DP(3)   

!-----------------------------------------------------
! Polarizability: Alpha tensor diagonal elements  ...
!-----------------------------------------------------
REF_Alpha = [ 9.2d0 , 8.5d0 , 7.8d0 ]
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

   CALL system( "./env.sh save_cost_statement " )

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
