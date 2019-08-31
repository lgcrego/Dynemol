module cost_EH

    use type_m
    use constants_m
    use GA_QCModel_m , only : MO_erg_diff, Mulliken, Bond_Type, MO_character, Localize, Exclude, me => i_


    public :: evaluate_cost , REF_DP , REF_Alpha

    ! module variables ...
    real*8 :: REF_DP(3) , REF_Alpha(3)

    private 

contains
!
!
!
!============================================================================
 function evaluate_cost( system , OPT_UNI , basis , DP , Alpha_ii , ShowCost)
!============================================================================
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
real*8   :: eval(200) = D_zero
real*8   :: REF_DP(3) , REF_Alpha(3)

!-------------------------------------------------------------------------
! Energy gaps ...     
! MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , {weight} )
! {...} terms are optional 
!-------------------------------------------------------------------------
eval(me) = MO_erg_diff( OPT_UNI, 115, 114,  2.6470, weight = 3.0 )
eval(me) = MO_erg_diff( OPT_UNI, 114, 113,  0.3040 )
eval(me) = MO_erg_diff( OPT_UNI, 115, 113,  2.9510 )
eval(me) = MO_erg_diff( OPT_UNI, 113, 112,  0.8950 )
eval(me) = MO_erg_diff( OPT_UNI, 112, 111,  0.4360 )
eval(me) = MO_erg_diff( OPT_UNI, 117, 116,  1.6000 )

!-------------------------------------------------------------------------
! MO character ...
! MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!-------------------------------------------------------------------------
eval(me) =  MO_character(OPT_UNI, basis, MO=111, AO='Pz') 
eval(me) =  MO_character(OPT_UNI, basis, MO=112, AO='Pz') 
eval(me) =  MO_character(OPT_UNI, basis, MO=117, AO='Pz') 
eval(me) =  MO_character(OPT_UNI, basis, MO=118, AO='Pz') 

!-------------------------------------------------------------------------
! Localize( OPT_UNI , basis , MO , {atom}=[:] , {residue} , {threshold} )
! {...} terms are optional 
! default criterium (threshold=0.85): localized > 85% of total population
!-------------------------------------------------------------------------
eval(me) =  Localize(OPT_UNI, basis, MO=117, EHSymbol = "CA", threshold = 0.5 )    

eval(me) =  Exclude(OPT_UNI, basis, MO=115, atom   = [6 ], threshold =2.5d-2) 
eval(me) =  Exclude(OPT_UNI, basis, MO=116, atom   = [22], threshold =2.5d-2) 

eval(me) =  Mulliken(OPT_UNI, basis, MO=117, atom=[33:40]) - Mulliken(OPT_UNI, basis, MO=117, atom=[66:73])

if(1==2) then
!-------------------------------------------------------------------------
! Bond Type analysis ...
! Bond_Type( system , OPT_UNI , MO , atom1 , atom2 , AO , "+" or "-" )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!-------------------------------------------------------------------------
eval(me) =  Bond_Type(system, OPT_UNI, 112, 23, 28, 'Pz', '+')                                
eval(me) =  Bond_Type(system, OPT_UNI, 112, 16, 19, 'Pz', '+')                                
eval(me) =  Bond_Type(system, OPT_UNI, 112, 16, 11, 'Pz', '+')                                
eval(me) =  Bond_Type(system, OPT_UNI, 112, 23, 19, 'Pz', '-')                                
eval(me) =  Bond_Type(system, OPT_UNI, 112, 28, 27, 'Pz', '-')                                

eval(me) =  Bond_Type(system, OPT_UNI, 113, 23, 28, 'Pz', '-')                                

eval(me) =  Bond_Type(system, OPT_UNI, 114, 23, 28, 'Pz', '+')                                

eval(me) =  Bond_Type(system, OPT_UNI, 115, 25, 21, 'Pz', '-')                                
eval(me) =  Bond_Type(system, OPT_UNI, 115, 23, 28, 'Pz', '+')                                

eval(me) =  Bond_Type(system, OPT_UNI, 116, 25, 21, 'Pz', '+')                                
eval(me) =  Bond_Type(system, OPT_UNI, 116, 23, 28, 'Pz', '-')                                

eval(me) =  Bond_Type(system, OPT_UNI, 117, 25, 21, 'Pz', '-')                                
eval(me) =  Bond_Type(system, OPT_UNI, 117, 23, 28, 'Pz', '-')                                

!-------------------------------------------------------------------------
! NO charge in these atoms ...
! Exclude( OPT_UNI , basis , MO , {atom}=[:] , {residue} , {threshold} )
! {...} terms are optional  
! default threshold < 1.d-3 
!-------------------------------------------------------------------------
eval(me) =  Exclude(OPT_UNI, basis, MO=111, atom   = [16], threshold =2.5d-2) 
eval(me) =  Exclude(OPT_UNI, basis, MO=111, atom   = [30], threshold =2.5d-2) 

!-------------------------------------------------------------------------
! Population analysis ...
! Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!-------------------------------------------------------------------------
! too much charge on these atoms ...
eval(me) =  Mulliken(OPT_UNI, basis, MO=112, residue="BZN",  weight = 3.0 ) 
eval(me) =  Mulliken(OPT_UNI, basis, MO=118, residue="PPH",  weight = 3.0 )      

end if

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

end If
!......................................................................

! evaluate total cost ...
evaluate_cost = dot_product(eval,eval) 

! just touching variables ...
dumb = basis(1)%atom
dumb = DP(1)
dumb = Alpha_ii(1) 

!reset index for next round ...
me = 0

end function evaluate_cost
!
!
!
end module cost_EH
