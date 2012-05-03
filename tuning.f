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
!================================
 subroutine ad_hoc_tuning( univ )
!================================
implicit none
type(universe) , intent(inout) :: univ

! local variables ...
integer :: i , ioerr

! edit structure  ...

!-----------------------------------
!      define %residue
!-----------------------------------

!-----------------------------------
!      define %nr
!-----------------------------------

!------------------------------------
!      define %DPF (Dipole Fragment) 
!------------------------------------

 where( univ % atom % residue == "ION" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP1" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP2" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP3" ) univ % atom % DPF = .true.

!-----------------------------------
!      define %solute 
!-----------------------------------

!-----------------------------------
!      define %El 
!-----------------------------------

!-----------------------------------
!      define %Hl
!-----------------------------------

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
!   Bridge      =   B
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
