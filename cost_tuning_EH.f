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
eval(me) = MO_erg_diff( OPT_UNI, 17, 16, 0.29d0 )
eval(me) = MO_erg_diff( OPT_UNI, 16, 15, 2.81d0 , 0.2 )
eval(me) = MO_erg_diff( OPT_UNI, 15, 14, 0.06d0 )
eval(me) = MO_erg_diff( OPT_UNI, 14, 13, 1.08d0 )
eval(me) = MO_erg_diff( OPT_UNI, 13, 12, 8.76d0 )
eval(me) = MO_erg_diff( OPT_UNI, 12, 11, 0.21d0 , 10.0)
eval(me) = MO_erg_diff( OPT_UNI, 11, 10, 0.28d0 , 10.0)
eval(me) = MO_erg_diff( OPT_UNI, 10, 9 , 1.22d0 )
eval(me) = MO_erg_diff( OPT_UNI, 9 , 8 , 0.83d0 )
eval(me) = MO_erg_diff( OPT_UNI, 8 , 7 , 0.13d0 )
eval(me) = MO_erg_diff( OPT_UNI, 7 , 6 , 3.16d0 )
eval(me) = MO_erg_diff( OPT_UNI, 6 , 5 , 0.10d0 )
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

!3 ===================
!eval(me) =  exclude (OPT_UNI, basis, MO=39, AO="S",  EHSymbol = "HW", from_to = real_interval( 0.91 , 0.30 ), adaptive  = input_mode) 

!5 LUMO ===================
!eval(me) =  MO_character( OPT_UNI , basis , MO=41 , AO='S')
eval(me) =  exclude (OPT_UNI, basis,  MO=41, AO="Py", EHSymbol = "$$", from_to = real_interval( 0.42 , 0.40 ), adaptive  = input_mode) 
eval(me) =  exclude(OPT_UNI, basis,   MO=41, AO="S",  EHSymbol = "$$", from_to = real_interval( 0.325 , 0.32 ), adaptive  = input_mode) 
!eval(me) =  exclude (OPT_UNI, basis,   MO=41, AO="S",  EHSymbol = "OW", from_to = real_interval( 0.93 , 0.12 ), adaptive  = input_mode) 
!eval(me) =  exclude (OPT_UNI, basis,   MO=41, AO="Pz", EHSymbol = "OW", from_to = real_interval( 0.97 , 0.36 ), adaptive  = input_mode) 
!eval(me) =  Bond_Type( sys , OPT_UNI , MO=41 , atom1=1 , AO1="S" , atom2=2 , AO2="S" , instance="-" ) 
!eval(me) =  Bond_Type( sys , OPT_UNI , MO=41 , atom1=1 , AO1="S" , atom2=3 , AO2="S" , instance="-" ) 
!
!6 LUMO+1 ===================
!eval(me) =  MO_character( OPT_UNI , basis , MO=42 , AO='Py')
!eval(me) =  Localize (OPT_UNI, basis, MO=42, AO="S", EHSymbol = "HW", from_to = real_interval( 0.1, 0.5), adaptive  = input_mode) 
!eval(me) =  Localize (OPT_UNI, basis, MO=42, AO="Py", EHSymbol = "OW", from_to = real_interval( 0.1, 0.4), adaptive = input_mode) 

!-------------------------                                                         
! Total DIPOLE moment ...
!-------------------------
REF_DP = [ 0.0d0 , 0.0d0 , 2.2d0 ]
eval(me+1) = DP(1) - REF_DP(1)     
eval(me+2) = DP(2) - REF_DP(2)    
eval(me+3) = DP(3) - REF_DP(3) 
me = me + 3

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
