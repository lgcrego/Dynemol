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
eval(me) = MO_erg_diff( OPT_UNI, 50, 49, 3.8153d0 )
eval(me) = MO_erg_diff( OPT_UNI, 51, 49, 5.1876d0 )
eval(me) = MO_erg_diff( OPT_UNI, 50, 48, 5.3631d0 )
eval(me) = MO_erg_diff( OPT_UNI, 49, 48, 1.5478d0 )
eval(me) = MO_erg_diff( OPT_UNI, 51, 50, 1.3722d0 )
eval(me) = MO_erg_diff( OPT_UNI, 48, 47, 0.2046d0 )
eval(me) = MO_erg_diff( OPT_UNI, 47, 46, 0.0734d0 )
eval(me) = MO_erg_diff( OPT_UNI, 52, 51, 0.2484d0 )
!===========
!  46    H-3
!  47    H-2
!  48    H-1
!  49    HOMO
!  50    LUMO
!  51    L+1
!  52    L+2
!  53    L+3
!----------------------------------------------------------------------------------------------
! ==> MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!
! ==> Localize( OPT_UNI , basis , MO , {atom}=[:] , {EHSymbol} , {residue} , {threshold} , {slide} , {adaptive} )
! {...} terms are optional
! default criterium (threshold=0.85): localized > 85% of total population
! slide = real_interval( begin , end ) : no need to use {threshold} if {slide} is used
! adaptive = {mode,lock} : logical flag to enable adpative GA method , lock sets up threshold = end
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
! adaptive = {mode,lock} : logical flag to enable adpative GA method, lock sets up threshold = end
!
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
! weight < 0  ==> does not update "me" when Mulliken in called
!----------------------------------------------------------------------------------------------

!46 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=46 , AO='Pz')

eval(me) =  Localize(OPT_UNI, basis, MO=46, residue = "CYT", slide = real_interval( 0.15 , 0.55), adaptive = lock)    
eval(me) =  Localize(OPT_UNI, basis, MO=46, residue = "GUA", slide = real_interval( 0.10 , 0.40), adaptive = lock)    

eval(me) =  Exclude (OPT_UNI, basis, MO=46, atom = [14:17], slide = real_interval( 0.21 , 0.01), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 46, 13, 'Pz', 22, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46,  2, 'Pz', 11, 'Pz', '-')                                

!47 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=47 , AO='Pz', y_or_n='n' )

eval(me) =  Localize(OPT_UNI, basis, MO=47, residue = "GUA", slide = real_interval( 0.31 , 0.93), adaptive = lock )    

!48 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=48 , AO='Pz')

eval(me) =  Localize(OPT_UNI, basis, MO=48, residue = "CYT", slide = real_interval( 0.45, 0.98 ), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 48,  3, 'Pz',  5, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48,  2, 'Pz', 11, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48,  2, 'Pz',  3, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48,  2, 'Pz', 13, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 11, 'Pz', 13, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48,  5, 'Pz',  8, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48,  5, 'Pz',  7, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 11, 'Pz', 20, 'Pz', '-')                                

eval(me) =  Exclude (OPT_UNI, basis, MO=48, atom = [ 7], slide = real_interval( 0.28, 0.15), adaptive = lock) 

!49 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=49, residue = "GUA", slide = real_interval( 0.43, 0.97 ), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 49, 18, 'Pz', 17, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 17, 'Pz', 24, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 24, 'Pz', 14, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 16, 'Pz', 15, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 21, 'Pz', 23, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 18, 'Pz', 19, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 16, 'Pz', 17, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 14, 'Pz', 15, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 23, 'Pz', 24, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 21, 'Pz', 22, 'Pz', '-')                                

eval(me) =  Exclude (OPT_UNI, basis, MO=49, atom=[20], slide = real_interval( 0.21 , 0.005 ), adaptive = lock)    

!50 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=50, residue = "CYT", slide = real_interval( 0.21 , 0.91), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 50,  5, 'Pz',  7, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50,  3, 'Pz',  5, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50,  2, 'Pz',  3, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50,  7, 'Pz', 11, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50,  7, 'Pz',  8, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 12, 'Pz',  2, 'Pz', '-')                                

eval(me) =  Exclude(OPT_UNI, basis, MO=50, atom = [12], slide = real_interval( 0.11, 0.07 ), adaptive = lock) 

!51 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=51, residue = "CYT", slide = real_interval( 0.30 , 0.95), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 51,  2, 'Pz',  3, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  7, 'Pz', 11, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  3, 'Pz',  5, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  2, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  3, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51, 11, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  5, 'Pz',  7, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  7, 'Pz',  8, 'Pz', '-')                                

!52 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=52, residue="CYT", slide = real_interval( 0.00 , 0.85), adaptive = lock )

!53 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=53, residue="GUA", slide = real_interval( 0.00, 0.85), adaptive = lock )

eval(me) =  Bond_Type(sys, OPT_UNI, 53, 18, 'Pz', 17, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 16, 'Pz', 24, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 19, 'Pz', 18, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 18, 'Pz', 20, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 17, 'Pz', 24, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 20, 'Pz', 12, 'Pz', '-')                                

eval(me) =  Exclude(OPT_UNI, basis, MO=53, atom=[21], slide = real_interval( 0.303, 0.003 ) , adaptive = lock ) 
eval(me) =  Exclude(OPT_UNI, basis, MO=53, atom=[22], slide = real_interval( 0.31, 0.007 ), adaptive = lock ) 

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
