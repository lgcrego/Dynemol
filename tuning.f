module tuning_m

    use type_m
    use constants_m

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

    !===================================
    !      define LIGAND_FMO atoms
    !===================================


end If


! edit structure  .....................................................

If( present(univ) ) then


    !===================================
    !      define LIGAND_FMO atoms
    !===================================

    where( univ % atom % residue == "ION" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP1" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP2" ) univ % atom % FMO = .true.
    where( univ % atom % residue == "BP3" ) univ % atom % FMO = .true.

    where( univ % atom % FMO ) univ % atom % solute = .true.

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
        case( 'ION') 
            a%atom(i)%fragment = 'I' 
        case( 'PYR') 
            a%atom(i)%fragment = 'P' 
        case( 'BP1') 
            a%atom(i)%fragment = 'D' 
        case( 'BP2') 
            a%atom(i)%fragment = '1' 
        case( 'BP3') 
            a%atom(i)%fragment = '2' 
        case default
            a%atom(i)%fragment = '#' 
    end select

 END DO

end subroutine Setting_Fragments

end module tuning_m 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
