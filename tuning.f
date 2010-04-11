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

    where( struct % residue == "ALQ" ) struct % k_WH = 2.d0

    where( struct % residue == "ION" ) struct % fragment = "S"

end If


! edit structure  .....................................................

If( present(univ) ) then

    !===================================
    !       charge of the atoms ...
    !===================================

    OPEN(unit=33,file='charge.dat',status='old',iostat=ioerr,err=11)
    do i = 1 , count(univ % atom % residue =="CCC")
        read(33,*,iostat=ioerr) univ%atom(i)%charge
    end do
    close(33)

11  if( ioerr > 0 ) stop 'charge.dat file not found; terminating execution'

    !===================================
    !      define MM atom types ...
    !===================================

    where( univ % atom % charge ==  2.1960d0 ) univ % atom % MMSymbol = "TiB"
    where( univ % atom % charge ==  1.6470d0 ) univ % atom % MMSymbol = "TiD"
    where( univ % atom % charge == -1.098d0  ) univ % atom % MMSymbol = "OB" 
    where( univ % atom % charge == -0.8235d0 ) univ % atom % MMSymbol = "O2c"
    where( univ % atom % charge == -1.647d0  ) univ % atom % MMSymbol = "O3c" 

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
