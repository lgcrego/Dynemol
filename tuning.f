module tuning_m

    use type_m
    use constants_m

    public :: Setting_Fragments

    contains
!
!
!
!=================================
 subroutine Setting_Fragments( a )
!=================================
implicit none
type(universe)  , intent(inout) :: a

! local variables ...
integer  :: i 

! ---------- Table of fragments -------------
!   Acceptor    =   A       
!   Donor       =   D 
!   Molecule    =   M
!   Solvent     =   S
!   Cluster     =   C 
!   Passivator  =   P 
!   ghost       =   #
!--------------------------------------------

 DO i = 1 , size(a%atom)

    select case(a%atom(i)%residue)
        case( 'CCC') 
            a%atom(i)%fragment = 'C' 
        case( 'ALQ') 
            a%atom(i)%fragment = 'M' 
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
        case( 'PYR') 
            a%atom(i)%fragment = 'P' 
        case default
            a%atom(i)%fragment = '#' 
    end select

 END DO

end subroutine Setting_Fragments

end module tuning_m 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TiO2 ...
!where( system % atom % MMSymbol == "TiB" ) system % atom % residue = "BLK"
!where( system % atom % MMSymbol == "OB"  ) system % atom % residue = "BLK"
!where( system % atom % MMSymbol == "TiD" ) system % atom % residue = "TiD"
!where( system % atom % MMSymbol == "O2c" ) system % atom % residue = "O2c"
!where( system % atom % MMSymbol == "O3c" ) system % atom % residue = "O3c"
