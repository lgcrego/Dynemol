module tuning_m

    use type_m
    use constants_m
    use parameters_m    , only  : T_ , F_

    public :: Setting_Fragments , ad_hoc_tuning

    private

    logical , save  :: ad_hoc_verbose_ = T_

    contains
!
!
!
!=========================================
 subroutine ad_hoc_tuning( struct , univ )
!=========================================
implicit none
type(structure) , optional  , intent(inout) :: struct
type(universe)  , optional  , intent(inout) :: univ

! local variables ...
integer :: i , ioerr

! edit structure  .....................................................

If( present(struct) ) then

    where( struct % FMO       ) struct % nr = 1
    where( .not. struct % FMO ) struct % nr = struct % nr - 3

end If


! edit structure  .....................................................

If( present(univ) ) then


    where( univ % atom % residue == "ION" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP1" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP2" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP3" ) univ % atom % FMO = .true.

    where( univ % atom % FMO       ) univ % atom % nr = 1
    where( .not. univ % atom % FMO ) univ % atom % nr = univ % atom % nr - 3

    univ % atom % solute = univ % atom % FMO    

end if

!......................................................................

If( ad_hoc_verbose_ ) then
    Print 46
    ad_hoc_verbose_ = F_
end If

include 'formats.h'

end subroutine ad_hoc_tuning
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
!   Solute      =   R
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

        case( 'PYR') 
            a%atom(i)%fragment = 'P' 

        case( 'H2O' , 'SOL' ) 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 2.0d0
        
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 3.d0

        case( 'ION') 
            a%atom(i)%fragment = 'I' 
            a%atom(i)%solvation_hardcore = 7.d0

        case( 'BP1') 
            a%atom(i)%fragment = 'D' 
            a%atom(i)%solvation_hardcore = 7.d0

        case( 'BP2') 
            a%atom(i)%fragment = '1' 
            a%atom(i)%solvation_hardcore = 7.d0

        case( 'BP3') 
            a%atom(i)%fragment = '2' 
            a%atom(i)%solvation_hardcore = 7.d0

        case default
            a%atom(i)%fragment = '#' 

    end select

 END DO

end subroutine Setting_Fragments

end module tuning_m 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
