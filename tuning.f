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

    where( struct % DPF       ) struct % nr = 1
    where( .not. struct % DPF ) struct % nr = struct % nr - 3

end If


! edit structure  .....................................................

If( present(univ) ) then

    where( univ % atom % residue == "ION" ) univ % atom % DPF = .true.
    where( univ % atom % residue == "BP1" ) univ % atom % DPF = .true.
    where( univ % atom % residue == "BP2" ) univ % atom % DPF = .true.
    where( univ % atom % residue == "BP3" ) univ % atom % DPF = .true.

    where( univ % atom % DPF       ) univ % atom % nr = 1
    where( .not. univ % atom % DPF ) univ % atom % nr = univ % atom % nr - 3

    univ % atom % solute = univ % atom % DPF        

    where( univ % atom % residue == "BP1" ) univ % atom % El = .true.
    where( univ % atom % residue == "BP1" .OR.  &
           univ % atom % residue == "BP2" .OR.  &
           univ % atom % residue == "BP3" .OR.  &
           univ % atom % residue == "ION" ) univ % atom % Hl = .true.

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
!
! fragments are set based on RESIDUE names ...
! 
!   Acceptor    =   A       
!   Donor       =   D 
!   Exciton     =   E 
!   Hole        =   H 
!   Molecule    =   M
!   Solvent     =   S
!   Solute      =   R
!   Cluster     =   C 
!   ghost       =   #
!
! some typical cases are used below ...
!--------------------------------------------

 DO i = 1 , size(a%atom)
 
    select case(a%atom(i)%residue)
        case( 'LFT') 
            a%atom(i)%fragment = 'L' 

        case( 'DON') 
            a%atom(i)%fragment = 'D' 

        case( 'RGT') 
            a%atom(i)%fragment = 'R' 

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
