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

 where( univ % atom % nr <= 4 ) univ % atom % nr = 1
 where( univ % atom % nr >  4 ) univ % atom % nr = univ % atom % nr - 3
 
!------------------------------------
!      define %DPF (Dipole Fragment) 
!------------------------------------

!default: %DPF = F_
 where( univ % atom % residue == "ION" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP1" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP2" ) univ % atom % DPF = .true.
 where( univ % atom % residue == "BP3" ) univ % atom % DPF = .true.

!use default: %DPF = %solute  
 where( univ % atom % DPF ) univ % atom % solute = .true.

!-----------------------------------
!      define %El   : mandatory !!
!-----------------------------------

 where( univ % atom % DPF ) univ % atom % El = .true.

!---------------------------------------------------
!      define %Hl   : must be T_ for El/Hl calcs ...
!---------------------------------------------------

 where( univ % atom % residue == "ION" ) univ % atom % Hl = .true.

!------------------------------------------------
!      define %fragments   : Donor fragment ...
!------------------------------------------------

!default: %El => DONOR
 where( univ % atom % El ) univ % atom % fragment = "D"

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

! --------- Table of STANDARD fragments ----------------
!
! STANDARD fragments are set based on RESIDUE names ...
! 
!   Acceptor    =   A       
!   Bridge      =   B
!   Donor       =   D  (defined in ad_hoc)
!   Electron    =   E  (defined in ad_hoc)
!   Hole        =   H 
!   Molecule    =   M
!   Solvent     =   S
!   Cluster     =   C 
!   System      =   #
!
! some typical cases are used below ...
!--------------------------------------------------------

 DO i = 1 , size(a%atom)
 
    select case(a%atom(i)%residue)

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
