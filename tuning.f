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

 univ % atom(49:50) % residue = "CCC"

!-----------------------------------
!      define %nr
!-----------------------------------

!------------------------------------
!      define %DPF (Dipole Fragment) 
!------------------------------------

!default: %DPF = F_
!use default: %DPF = %solute  
! where( univ % atom % DPF ) univ % atom % solute = .true.

!-----------------------------------
!      define %El   : mandatory !!
!-----------------------------------

 where( univ % atom % residue == "TA1" ) univ % atom % El = .true.

!---------------------------------------------------
!      define %Hl   : must be T_ for El/Hl calcs ...
!---------------------------------------------------

 where( univ % atom % residue == "TA1" ) univ % atom % Hl = .true.

!------------------------------------------------
!      define %fragments   : Donor fragment define here ...
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

        case( 'CCC') 
            a%atom(i)%fragment = 'A' 

        case( 'H2O' , 'SOL' ) 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 2.0d0
        
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
            a%atom(i)%solvation_hardcore = 3.d0

    end select

 END DO

end subroutine Setting_Fragments

end module tuning_m 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
