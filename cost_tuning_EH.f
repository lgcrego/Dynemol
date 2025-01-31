module cost_EH

    use type_m
    use constants_m
    use util_m       , only : Lp_norm,            &
                              split_line,         &
                              count_lines,        &
                              parse_this,         &
                              truncate_array,     &
                              change_single_character_in_string
    use GA_QCModel_m , only : MO_erg_diff,        &  
                              Mulliken,           &
                              Bond_Type,          &
                              MO_character,       &
                              Localize,           &
                              Exclude,            &
                              Adaptive_GA,        &
                              me => i_       

    public :: evaluate_cost , REF_DP , REF_Alpha , parse_EH_cost_function

    private 

    ! module variables ...
    real*8  :: REF_DP(3) , REF_Alpha(3)
    type(GA_features) , allocatable :: MO_ERG_DIFF_parms(:) 
    type(GA_features) , allocatable :: MO_CHARACTER_parms(:)
    type(GA_features) , allocatable :: BOND_TYPE_parms(:)
    type(GA_features) , allocatable :: Mulliken_parms(:)
    type(GA_features) , allocatable :: Exclude_parms(:)
    type(GA_features) , allocatable :: Localize_parms(:)
    type(GA_features) :: append_parms

    ! module parameters ...

contains
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
logical  :: adaptive

integer          :: MO_up , MO_down , MO , atom1 , atom2
integer , allocatable :: atom(:)
real             :: de_ref , weight , ref
character(len=1) :: pm
character(len=2) :: EHSymbol, Symbol
character(len=3) :: residue
character(len=5) :: AO , AO1 , AO2
type(real4_interval) :: from_to

adaptive = Adaptive_GA% mode

!-------------------------------------------------------------------------
! Energy gaps ...     
! MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , {weight} )
! {...} terms are optional 
!-------------------------------------------------------------------------
if(allocated(MO_ERG_DIFF_parms)) then
   do i = 1 , MO_ERG_DIFF_parms(1)% entries
   
            MO_up    = MO_ERG_DIFF_parms(i)% MO_up
            MO_down  = MO_ERG_DIFF_parms(i)% MO_down
            dE_ref   = MO_ERG_DIFF_parms(i)% dE_ref
            weight   = MO_ERG_DIFF_parms(i)% weight
   
            eval(me) = MO_erg_diff( OPT_UNI , MO_up , MO_down , dE_ref , weight )
   print*, i, me, MO_up , MO_down , dE_ref , weight
   end do
end if
!----------------------------------------------------------------------------------------------
! ==> MO_character( OPT_UNI , basis , MO , AO )
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!----------------------------------------------------------------------------------------------
if(allocated(MO_CHARACTER_parms)) then
   do i = 1 , MO_CHARACTER_parms(1)% entries
   
            MO = MO_CHARACTER_parms(i)% MO
            AO = MO_CHARACTER_parms(i)% AO
   
            eval(me) = MO_character( OPT_UNI , basis , MO , AO )
   print*, i, me, MO , AO
   end do
end if
!----------------------------------------------------------------------------------------------
! ==> Bond_Type( sys , OPT_UNI , MO , atom1 , AO1 , atom2 , AO2 , "+" or "-" )
! Bond Topolgy analysis ...
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!  + = Bonding               &         - = Anti_Bonding
!----------------------------------------------------------------------------------------------
if(allocated(BOND_TYPE_parms)) then
   do i = 1 , BOND_TYPE_parms(1)% entries
   
            MO    = BOND_TYPE_parms(i)% MO
            atom1 = BOND_TYPE_parms(i)% atom_1
            AO1   = BOND_TYPE_parms(i)% AO1
            atom2 = BOND_TYPE_parms(i)% atom_2
            AO2   = BOND_TYPE_parms(i)% AO2
            pm    = BOND_TYPE_parms(i)% pm_sign
   
            eval(me) = Bond_Type( sys , OPT_UNI , MO , atom1 , AO1 , atom2 , AO2 , pm )
   
   print*, i, me, MO , atom1 , AO1 , atom2 , AO2 , pm
   end do
end if
!----------------------------------------------------------------------------------------------
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {weight} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
! weight < 0  ==> does not update "me" when Mulliken is called
!----------------------------------------------------------------------------------------------
if(allocated(Mulliken_parms)) then
   do i = 1 , Mulliken_parms(1)% entries
   
            MO       = Mulliken_parms(i)% MO
            AO       = Mulliken_parms(i)% AO
            EHSymbol = Mulliken_parms(i)% EHSymbol
            Symbol   = Mulliken_parms(i)% Symbol
            residue  = Mulliken_parms(i)% residue
            weight   = Mulliken_parms(i)% weight
            ref      = Mulliken_parms(i)% ref
   
            eval(me) = Mulliken( OPT_UNI , basis , MO , Mulliken_parms(i)%atom , AO , EHSymbol , Symbol , residue , weight ) - max(ref,0.0)
   
   print*, i, me, MO , size(Mulliken_parms(i)% atom) ,  AO , EHSymbol , Symbol , residue , weight , ref
   end do
end if
!----------------------------------------------------------------------------------------------
! ==> Exclude( OPT_UNI , basis , MO , {atom}=[:] , {AO} , {EHSymbol} , {residue} , {reference} , {from_to} , {adaptive} )
! NO charge on these atoms ...
! {...} terms are optional  
! default reference < 0.001 
! from_to = real_interval( begin , end ) : no need to use {reference} if {from_to} is used
! adaptive = {input_mode,lock} : logical flag to enable adpative GA method, lock sets reference = end
!----------------------------------------------------------------------------------------------
if(allocated(Exclude_parms)) then
   do i = 1 , Exclude_parms(1)% entries
   
            MO       = Exclude_parms(i)% MO
            AO       = Exclude_parms(i)% AO
            EHSymbol = Exclude_parms(i)% EHSymbol
            residue  = Exclude_parms(i)% residue
            ref      = Exclude_parms(i)% ref
            from_to  = Exclude_parms(i)% from_to
            adaptive = Exclude_parms(i)% adaptive
   
            eval(me) = Exclude( OPT_UNI, basis, MO, Exclude_parms(i)%atom, AO, EHSymbol, residue, ref, from_to, adaptive )
   
   print*, MO , size(Exclude_parms(i)% atom) ,  AO , EHSymbol , residue , ref , from_to, adaptive
   end do
end if
!----------------------------------------------------------------------------------------------
! ==> Localize( OPT_UNI , basis , MO , {atom}=[:] , {AO} , {EHSymbol} , {residue} , {reference} , {from_to} , {adaptive} )
! {...} terms are optional
! default criterium (reference=0.85): localized > 85% of total population
! from_to = real_interval( begin , end ) : no need to use {reference} if {from_to} is used
! adaptive = {input_mode,lock} : logical flag to enable adpative GA method , lock sets reference = end
!----------------------------------------------------------------------------------------------
if(allocated(Localize_parms)) then
   do i = 1 , Localize_parms(1)% entries
   
            MO       = Localize_parms(i)% MO
            AO       = Localize_parms(i)% AO
            EHSymbol = Localize_parms(i)% EHSymbol
            residue  = Localize_parms(i)% residue
            ref      = Localize_parms(i)% ref
            from_to  = Localize_parms(i)% from_to
            adaptive = Localize_parms(i)% adaptive
   
            eval(me) = Localize( OPT_UNI, basis, MO, Localize_parms(i)%atom, AO, EHSymbol, residue, ref, from_to, adaptive )
   
   print*, MO , size(Localize_parms(i)% atom) ,  AO , EHSymbol , residue , ref , from_to, adaptive
   end do
end if
!39 ===================
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NH") - 0.6
!eval(me) =  Mulliken(OPT_UNI, basis, MO=39, AO="Pz",  EHSymbol = "NB") - 0.35
!eval(me) =  Mulliken(OPT_UNI, basis, MO=38, AO="Pz",  EHSymbol = "CP") - 0.05
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
evaluate_cost = Lp_norm(eval,p=1)

! just touching variables ...
dumb = basis(1)%atom

!reset index for next round ...
me = 0

include 'formats.h'

end function evaluate_cost
!
!
!
!
!==========================================================
 subroutine parse_EH_cost_function
!==========================================================
use util_m
implicit none

! local variables ...
 integer :: f_unit , ioerr 
 logical :: exist
 character(len=30) :: keyword
 character(len=80) :: line

! file error msg (more reliable than iostat) ...
inquire( file=dynemolworkdir//"cost_tuning.inpt", EXIST=exist )
if( .not. exist ) then
     CALL warning('file  "cost_tuning.inpt"  not found; terminating execution')
     stop
end if

open(file='cost_tuning.inpt', status='old', newunit=f_unit, iostat=ioerr )

!=====================================================================================
!  reading  the input CARD, one line at a time ...

read_loop: do 

      read(f_unit,'(A)',iostat=ioerr) line
      if ( ioerr /= 0 ) exit read_loop

      read(line,*,iostat=ioerr) keyword   ! <== keyword = first contiguous string from line 

      ! commented line ...
      if( index(keyword,"!") /= 0 ) cycle read_loop 

      keyword = TO_UPPER_CASE( keyword )

      select case ( keyword(1:7) )
                  case( "MO_ERG_" )
                         allocate( MO_ERG_DIFF_parms(1) )
                         call Build_GA_Parms( this_parms = MO_ERG_DIFF_parms, in=f_unit )

                  case( "MO_CHAR" )
                         allocate( MO_CHARACTER_parms(1) )
                         call Build_GA_Parms( this_parms = MO_CHARACTER_parms, in=f_unit )

                  case( "BOND_TY" )
                         allocate( BOND_TYPE_parms(1) )
                         call Build_GA_Parms( this_parms = BOND_TYPE_parms, in=f_unit )

                  case( "MULLIKE" )
                         allocate( Mulliken_parms(1) )
                         call Build_GA_Parms( this_parms = Mulliken_parms, in=f_unit )

                  case( "EXCLUDE" )
                         allocate( Exclude_parms(1) )
                         call Build_GA_Parms( this_parms = Exclude_parms, in=f_unit )

                  case( "LOCALIZ" )
                         allocate( Localize_parms(1) )
                         call Build_GA_Parms( this_parms = Localize_parms, in=f_unit )
      end select
      ! reset keyword ...
      keyword = "XXXXXXX"

end do read_loop

close(f_unit)

end subroutine parse_EH_cost_function
!
!
!
!
!===============================================
 subroutine Build_GA_Parms( this_parms , in )
!===============================================
use util_m
implicit none
type(GA_features) , allocatable , intent(inout) :: this_parms(:)
integer                         , intent(in)    :: in

! local variables ...
 integer :: i, j, row, ioerr, detecting_field
 integer :: bra, ket, field_positions(18), token_positions(18)
 real    :: inicio, fim
 character(len=:) , allocatable :: tokens(:)
 character(len= 6) :: flag
 character(len=80) :: line
 character(len=17) :: fields(18)=["MO=","MO_UP","MO_DOWN","DE_REF","ATOM1","AO1","ATOM2","AO2","PM_SIGN", &
                      "ATOM","AO","EHSYMBOL","SYMBOL","RESIDUE","WEIGHT","FROM_TO","ADAPTIVE","REF"]

row = 0 
do 
      read(in,'(A)',iostat=ioerr) line
      line = TO_UPPER_CASE( line )
      allocate( tokens , source=split_line( line , token_length=30 ) )
      if (index(tokens(1), "!") /= 0) then
         deallocate(tokens)
         cycle
      elseif (tokens(1) == "END") then
         deallocate(tokens)
         exit
      end if

      row = row + 1

      call initialize(append_parms)

      ! default value ...
      append_parms% weight = 1. 

      ! parsing the command line ...
      !----------------------------------------------------
      field_positions = 0
      token_positions = 0
      do i = 1 , size(fields)
      do j = 1 , size(tokens)
           detecting_field =  index( tokens(j) , trim(adjustL(fields(i))) ) 
           if( detecting_field == 1 ) &
           then
               field_positions(i) = i
               token_positions(i) = j
               exit
           end if
      end do
      end do

      do i = 1 , size(fields)

         if( field_positions(i) == 0 ) cycle

         select case (field_positions(i))

                case(1) ! <== MO
                     read(tokens(token_positions(i))(4:),*) append_parms% MO

                case(2) ! <== MO_up
                     read(tokens(token_positions(i))(7:),*) append_parms% MO_up

                case(3) ! <== MO_down
                     read(tokens(token_positions(i))(9:),*) append_parms% MO_down

                case(4) ! <== dE_ref
                     read(tokens(token_positions(i))(8:),*) append_parms% dE_ref

                case(5) ! <== atom1
                     read(tokens(token_positions(i))(7:),*) append_parms% atom_1

                case(6) ! <== AO1
                     append_parms% AO1 = tokens(token_positions(i))(5:)

                case(7) ! <== atom2
                     read(tokens(token_positions(i))(7:),*) append_parms% atom_2

                case(8) ! <== AO2
                     append_parms% AO2 = tokens(token_positions(i))(5:)

                case(9) ! <== pm_sign
                     append_parms% pm_sign = tokens(token_positions(i))(9:)

                case(10) ! <== Atom
                     bra = 6
                     ket= index( tokens(token_positions(i)) , "]" ) 

                     append_parms% atom = parse_this( tokens(token_positions(i)) (bra+1:ket-1) )

                case(11) ! <== AO
                     append_parms% AO = tokens(token_positions(i))(4:)
                     
                case(12) ! <== EHSymbol
                     append_parms% EHSymbol = tokens(token_positions(i))(10:)

                case(13) ! <== Symbol
                     append_parms% EHSymbol = tokens(token_positions(i))(8:)

                case(14) ! <== residue
                     append_parms% residue = tokens(token_positions(i))(9:)

                case(15) ! <== weight
                     read(tokens(token_positions(i))(8:),*) append_parms% weight

                case(16) ! <== from_to
                     bra = 9
                     ket= index( tokens(token_positions(i)) , ")" ) 
                     call change_single_character_in_string( tokens(token_positions(i)) , remove=":" , insert=" " )

                     read( tokens(token_positions(i)) (bra+1:ket-1) , * )  inicio , fim
                     append_parms% from_to = real4_interval(inicio,fim)

                case(17) ! <== adaptive
                     read(tokens(token_positions(i))(10:),*) flag
                     append_parms% adaptive = merge( .true. , .false. , any( [".TRUE.","TRUE","T","T_"] == flag ) ) 

                case(18) ! <== reference
                     read(tokens(token_positions(i))(5:),*) append_parms% ref
         end select

      end do
      !----------------------------------------------------
      deallocate(tokens)

      if( row > 1) then 
          this_parms = [ this_parms , append_parms ] 
      else
          this_parms = append_parms
      end if
end do

this_parms(1)%entries = row

end subroutine Build_GA_Parms
!
!
!
!=========================
subroutine initialize( a )
!=========================
implicit none
type(GA_features) , intent(inout) :: a

a% AO       = "XX"
a% AO1      = "XX"
a% AO2      = "XX"
a% residue  = "XXX"
a% EHSymbol = "XX"
a% Symbol   = "XX"
a% weight   = 1.0
a% ref      = -1.0
a% from_to  = real4_interval(-1.0,-1.0)
a% adaptive = .false.
if(allocated(a%atom)) deallocate(a%atom)

end subroutine initialize
!
!
!
!
end module cost_EH
