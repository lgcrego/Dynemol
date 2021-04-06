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
real*8   :: eval(200) 
real*8   :: REF_DP(3) , REF_Alpha(3)
logical  :: mode

mode = Adaptive_GA% mode

eval(:) = d_zero
!-------------------------------------------------------------------------
! Energy gaps ...     
! MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , {weight} )
! {...} terms are optional 
!-------------------------------------------------------------------------
eval(me) = MO_erg_diff( OPT_UNI, 50, 49, 4.9013d0 )
eval(me) = MO_erg_diff( OPT_UNI, 51, 49, 5.4700d0 )
eval(me) = MO_erg_diff( OPT_UNI, 50, 48, 5.5321d0 )
eval(me) = MO_erg_diff( OPT_UNI, 49, 48, 0.6308d0 )
eval(me) = MO_erg_diff( OPT_UNI, 51, 50, 0.5687d0 )
eval(me) = MO_erg_diff( OPT_UNI, 48, 47, 0.2857d0 )
eval(me) = MO_erg_diff( OPT_UNI, 47, 46, 0.4522d0 )
eval(me) = MO_erg_diff( OPT_UNI, 52, 51, 0.6248d0 )
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
eval(me) =  Localize(OPT_UNI, basis, MO=46, residue = "ADN", slide = real_interval( 0.9 , 0.91), adaptive = lock)    

eval(me) =  MO_character( OPT_UNI , basis , MO=46 , AO='Pz')

eval(me) =  Bond_Type(sys, OPT_UNI, 46, 15, 'Pz', 14, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46, 15, 'Pz',  6, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46, 12, 'Pz', 11, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46,  2, 'Pz',  3, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46, 12, 'Pz', 14, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46, 15, 'Pz',  2, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46,  8, 'Pz', 11, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 46,  8, 'Pz',  6, 'Pz', '+')                                

eval(me) =  Exclude (OPT_UNI, basis, MO=46, atom = [ 7], slide = real_interval( 0.30, 0.05), adaptive = lock) 
eval(me) =  Exclude (OPT_UNI, basis, MO=46, atom = [ 5], slide = real_interval( 0.30, 0.08), adaptive = lock) 

!47 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=47 , AO='Pz', y_or_n='n' )

eval(me) =  Localize(OPT_UNI, basis, MO=47, residue = "ADN", slide = real_interval( 0.5 , 0.91), adaptive = lock )    

!48 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=48 , AO='Pz')

eval(me) =  Localize(OPT_UNI, basis, MO=48, residue = "THY", slide = real_interval( 0.3, 0.97 ), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 48, 20, 'Pz', 18, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 20, 'Pz', 25, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 29, 'Pz', 17, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 17, 'Pz', 18, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 21, 'Pz', 20, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 30, 'Pz', 17, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 20, 'Pz', 26, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 20, 'Pz', 27, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 48, 17, 'Pz', 27, 'Pz', '+')                                

eval(me) =  Exclude (OPT_UNI, basis, MO=48, atom = [27], slide = real_interval( 0.15, 0.05), adaptive = lock) 
eval(me) =  Exclude (OPT_UNI, basis, MO=48, atom = [27], slide = real_interval( 0.15, 0.08), adaptive = lock) 

!49 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=49, residue = "ADN", slide = real_interval( 0.3, 0.97 ), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 49, 14, 'Pz', 12, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  3, 'Pz',  5, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  2, 'Pz', 15, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 15, 'Pz',  6, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  6, 'Pz',  7, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  7, 'Pz', 11, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  8, 'Pz',  7, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 14, 'Pz', 15, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  3, 'Pz',  2, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49,  5, 'Pz',  6, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 49, 12, 'Pz', 11, 'Pz', '-')                                

!50 ===================
eval(me) =  MO_character(OPT_UNI, basis, MO=50, AO='Pz')  
eval(me) =  Localize(OPT_UNI, basis, MO=50, residue = "THY", slide = real_interval( 0.3, 0.95 ), adaptive = mode)    

eval(me) =  Bond_Type(sys, OPT_UNI, 50, 20, 'Pz', 25, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 29, 'Pz', 17, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 26, 'Pz', 25, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 20, 'Pz', 18, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 17, 'Pz', 18, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 50, 27, 'Pz', 25, 'Pz', '-')                                

!51 ===================
eval(me) =  Localize(OPT_UNI, basis, MO=51, residue = "ADN", slide = real_interval( 0.30 , 0.95), adaptive = lock)    

eval(me) =  Bond_Type(sys, OPT_UNI, 51,  5, 'Pz',  6, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  2, 'Pz',  6, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  3, 'Pz',  2, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  3, 'Pz',  5, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  6, 'Pz',  7, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  7, 'Pz',  8, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51, 14, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51,  7, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 51, 14, 'Pz',  2, 'Pz', '-')                                

eval(me) =  Exclude(OPT_UNI, basis, MO=51, atom=[15], slide = real_interval( 0.07 , 0.05 ), adaptive = lock ) 
eval(me) =  Exclude(OPT_UNI, basis, MO=51, atom=[11], slide = real_interval( 0.05 , 0.03 ), adaptive = lock ) 

!52 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=52 , AO='Pz')

eval(me) =  Localize(OPT_UNI, basis, MO=52, residue="ADN", slide = real_interval( 0.30 , 0.6), adaptive = mode )

eval(me) =  Bond_Type(sys, OPT_UNI, 52,  6, 'Pz',  7, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52, 14, 'Pz', 12, 'Pz', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52,  6, 'Pz',  5, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52,  7, 'Pz',  8, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52,  7, 'Pz', 11, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52, 11, 'Pz', 12, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52,  2, 'Pz', 15, 'Pz', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 52, 15, 'Pz',  6, 'Pz', '-')                                

eval(me) =  Exclude(OPT_UNI, basis, MO=52, atom=[ 3], slide = real_interval( 0.007, 0.003), adaptive = mode ) 

!53 ===================
eval(me) =  MO_character( OPT_UNI , basis , MO=53 , AO='Pz', y_or_n='n' )

eval(me) =  Localize(OPT_UNI, basis, MO=53, residue="THY", slide = real_interval( 0.9 , 0.95), adaptive = mode )

eval(me) =  Bond_Type(sys, OPT_UNI, 53, 16,  's', 19,  's', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 17, 'Px', 18, 'Px', '+')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 16,  'S', 17, 'Px', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 18, 'Px', 19,  'S', '-')                                
eval(me) =  Bond_Type(sys, OPT_UNI, 53, 20, 'Px', 18, 'Px', '+')                                

eval(me) =  Exclude (OPT_UNI, basis, MO=53, atom=[   29] , slide = real_interval( 0.05,0.008), adaptive = mode )
eval(me) =  Localize(OPT_UNI, basis, MO=53, atom=[16,19] , slide = real_interval( 0.45, 0.55), adaptive = mode )

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
