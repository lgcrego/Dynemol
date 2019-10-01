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
eval(me) = MO_erg_diff( OPT_UNI, 123, 122,  2.8670, weight = 2.0 )
eval(me) = MO_erg_diff( OPT_UNI, 122, 121,  0.0930, weight = 2.0 )
eval(me) = MO_erg_diff( OPT_UNI, 123, 121,  2.9600 )
eval(me) = MO_erg_diff( OPT_UNI, 121, 120,  1.0970, weight = 2.0 )
eval(me) = MO_erg_diff( OPT_UNI, 120, 119,  0.2020 )
eval(me) = MO_erg_diff( OPT_UNI, 125, 124,  1.6310 )

!----------------------------------------------------------------------------------------------
! ==> MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!
! ==> Localize( OPT_UNI , basis , MO , {atom}=[:] , {residue} , {threshold} )
! {...} terms are optional 
! default criterium (threshold=0.85): localized > 85% of total population
!
! ==> Bond_Type( sys , OPT_UNI , MO , atom1 , atom2 , AO , "+" or "-" )
! Bond Topolgy analysis ...
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!
! ==> Exclude( OPT_UNI , basis , MO , {atom}=[:] , {residue} , {threshold} )
! NO charge on these atoms ...
! {...} terms are optional  
! default threshold < 0.001 
!
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
! weight < 0  ==> does not update "me" when Mulliken in called
!----------------------------------------------------------------------------------------------

!119 ===================
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 23, 'Pz', 28, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 19, 'Pz', 11, 'Pz', '-')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 20, 'Pz', 22, 'Pz', '-')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 27, 'Pz', 22, 'Pz', '+')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 19, 'Pz', 22, 'Pz', '+')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 19, 'Pz', 23, 'Pz', '-')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 20, 'Pz', 24, 'Pz', '+')  
eval(me) =  Bond_Type(sys, OPT_UNI, 119, 20, 'Pz', 18, 'Pz', '+')  

eval(me) =  Exclude (OPT_UNI, basis, MO=119, atom = [30], threshold =0.05 ) 
eval(me) =  Exclude (OPT_UNI, basis, MO=119, atom = [16], threshold =0.05 ) 
eval(me) =  Exclude (OPT_UNI, basis, MO=119, residue = "COO", threshold =0.13 ) 

!120 ===================
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 23, 'Pz', 28, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 12, 'Pz',  7, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 16, 'Pz', 19, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 16, 'Pz', 11, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 22, 'Pz', 19, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 22, 'Pz', 27, 'Pz', '+')                                

eval(me) =  Bond_Type(sys, OPT_UNI, 120, 79, 'Px', 80, 'Px', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 78, 'Px', 79, 'Px', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 79, 'Px', 81, 'S ', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 78, 'Py', 81, 'S ', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 120, 78, 'Pz', 81, 'S ', '+')                                

eval(me) =  Localize(OPT_UNI, basis, MO=120, atom = [81], vary=real_interval(0.1,0.15), adaptive=mode )
eval(me) =  Localize(OPT_UNI, basis, MO=120, residue = "COO", vary=real_interval(0.4,0.60), adaptive=mode )    

eval(me) =  Exclude (OPT_UNI, basis, MO=120, atom = [77], threshold = 0.15 ) 
eval(me) =  Exclude (OPT_UNI, basis, MO=120, atom = [5,7,9,12,21,23,25,28], threshold = 0.30 ) 

!121 ===================
eval(me) =  Bond_Type(sys, OPT_UNI, 121, 23, 'Pz', 28, 'Pz', '-')                                
eval(me) =  Exclude (OPT_UNI, basis, MO=121, atom = [30], threshold =0.025 ) 
eval(me) =  Exclude (OPT_UNI, basis, MO=121, atom = [6,22], threshold =0.05 ) 
eval(me) =  Localize(OPT_UNI, basis, MO=121, atom=[1:30], threshold =0.7 )    

!122 ===================
eval(me) =  Bond_Type(sys, OPT_UNI, 122, 23, 'Pz', 28, 'Pz', '+')                                
eval(me) =  Localize(OPT_UNI, basis, MO=122, atom=[1:30], threshold = 0.7 )    

!123 ===================
eval(me) =  Exclude(OPT_UNI, basis, MO=123, atom = [6 ], threshold =0.02 ) 
eval(me) =  Bond_Type(sys, OPT_UNI, 123, 25, 'Pz', 21, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 123, 23, 'Pz', 28, 'Pz', '+')                                

!124 ===================
eval(me) =  Exclude(OPT_UNI, basis, MO=124, atom = [22], threshold =0.02 ) 
eval(me) =  Bond_Type(sys, OPT_UNI, 124, 25, 'Pz', 21, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 124, 23, 'Pz', 28, 'Pz', '-')                                

!125 ===================

eval(me) =  Localize(OPT_UNI, basis, MO=125, EHSymbol = "CA", threshold = 0.60 )    

eval(me) =  Bond_Type(sys, OPT_UNI, 125, 25, 'Pz', 21, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 125, 23, 'Pz', 28, 'Pz', '-')                                

eval(me) =  Localize(OPT_UNI, basis, MO=125, atom=[33:40], threshold = 0.18 ) 
eval(me) =  Localize(OPT_UNI, basis, MO=125, atom=[44:51], threshold = 0.15 )
eval(me) =  Localize(OPT_UNI, basis, MO=125, atom=[55:62], threshold = 0.15 )
eval(me) =  Localize(OPT_UNI, basis, MO=125, atom=[66:73], threshold = 0.15 )

!126 ===================

eval(me) =  Exclude(OPT_UNI, basis, MO=126, atom=[ 1:30] , threshold = 0.2 ) 

eval(me) =  Localize(OPT_UNI, basis, MO=126, atom=[33:40], threshold = 0.20 ) 
eval(me) =  Localize(OPT_UNI, basis, MO=126, atom=[44:51], threshold = 0.20 )
eval(me) =  Localize(OPT_UNI, basis, MO=126, atom=[55:62], threshold = 0.20 )
eval(me) =  Localize(OPT_UNI, basis, MO=126, atom=[66:73], threshold = 0.20 )

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

   open( unit=33 , file='opt_trunk/view_cost.dat' , status='unknown' )

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
