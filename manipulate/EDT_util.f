module EDT_util_m

use types_m
use util_m        , only: split_line , seq_in_range
use Read_Parms
use Constants_m
use GMX_routines

public :: on_the_fly_tuning , parse_this

private

        !module variables ...
        integer               :: nresid
        integer , allocatable :: displace(:) !displace is USED IN THE CALLING ROUTINE to accomodate the residue numbers after delete operations
        character(len=4) :: resid
        character(len=4) :: MMSymbol
        character(len=1) :: fragment

contains

!
!=====================================
subroutine on_the_fly_tuning( system )
!=====================================
implicit none
type(universe), intent(inout) :: system
    
! local variables
integer :: option

CALL display_menu(option)

select case ( option )

    case (1)
        Write(*,*) "not implemented"

    case (2)
        write(*,'(/a)',advance='no') ">>> enter new MMSymbol: "
        read (*,'(a)') MMSymbol
        CALL change_this( system , instance="MMSymbol" )

    case (3)
        write(*,'(/a)',advance='no') ">>> enter new fragment: "
        read (*,'(a)') fragment
        CALL change_this( system , instance="fragment" )

    case (4)
        write(*,'(/a)',advance='no') ">>> enter new residue name: "
        read (*,'(a)') resid
        CALL change_this( system , instance="resid" )

    case (5)
        write(*,'(/a)',advance='no') ">>> enter new residue number: "
        read (*,*) nresid
        CALL change_this( system , instance="nresid" )

    case (6)
        CALL copy_mask( system )

    case (7)
        CALL delete_mask( system )

    case (8)
        CALL translation_mask( system )

    case (9)
        CALL rotation_mask( system )

    case default
        return

end select

end subroutine on_the_fly_tuning
!
!
!
!====================================
subroutine copy_mask( system )
!====================================
implicit none 
type(universe)              , intent(inout) :: system

! local variables ...
integer                        :: i
integer          , allocatable :: nr(:) , indx(:)
character(len=1)               :: choice
character(len=3) , allocatable :: residue(:)
character(len=80)              :: line

! reset varible ...
system% atom(:)% translate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose keyword of stuff to TRANSLATE : '
write(*,'(/a)') ' (1) = tuning already done'
write(*,'(/a)') ' (2) = residue number '
write(*,'(/a)') ' (3) = residue name '
write(*,'(/a)') ' (4) = atom index '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )
    case( '1' ) 
        ! do nothing, proceed ...

    case( '2' )
        write(*,'(1x,3/a)') "enter the residue numbers to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        nr =  parse_this(line)
        do i = 1 , size(nr)
           where( system% atom(:)% nresid == nr(i) ) system% atom(:)% copy = .true.
        end do

    case( '3' )
        write(*,'(1x,3/a)') "enter the residue names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        residue = split_line(line)
        do i = 1 , size(residue)
           where( system% atom(:)% resid == residue(i) ) system% atom(:)% copy = .true.
        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        indx =  parse_this(line)

        system% atom(indx)% copy = .true.

end select

end subroutine copy_mask
!
!
!
!==========================================
subroutine change_this( system , instance )
!==========================================
implicit none 
type(universe)   , intent(inout) :: system
character(len=*) , intent(in)    :: instance

! local variables ...
integer                        :: i 
integer          , allocatable :: nr(:) , indx(:)
character(len=1)               :: choice
character(len=3) , allocatable :: residue(:) , MM_name(:)
character(len=80)              :: line

! reset varible ...
system% atom(:)% translate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose the feature that controls the operation : '
write(*,'(/a)') ' (1) = tuning already done'
write(*,'(/a)') ' (2) = residue number '
write(*,'(/a)') ' (3) = residue name '
write(*,'(/a)') ' (4) = atom index '
write(*,'(/a)') ' (5) = MMSymbol '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )
    case( '1' ) 
        ! do nothing, proceed ...

    case( '2' )
        write(*,'(1x,3/a)') "enter the residue numbers to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        nr =  parse_this(line)
        do i = 1 , size(nr)
           CALL change_via_nr( system , nr(i) , instance )
        end do

    case( '3' )
        write(*,'(1x,3/a)') "enter the residue names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        residue = split_line(line)
        do i = 1 , size(residue)
           CALL change_via_residue( system , residue(i) , instance )
        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        indx =  parse_this(line)
        CALL change_via_index( system , indx , instance )

    case( '5' )
        write(*,'(1x,3/a)') "enter the MMSymbol names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        MM_name = split_line(line)
        do i = 1 , size(MM_name)
           CALL change_via_MM_name( system , MM_name(i) , instance )
        end do
 
end select

CALL systemQQ( "clear" )

end subroutine change_this
!
!
!
!====================================
subroutine translation_mask( system )
!====================================
implicit none 
type(universe)              , intent(inout) :: system

! local variables ...
integer                        :: i
integer          , allocatable :: nr(:) , indx(:)
character(len=1)               :: choice
character(len=3) , allocatable :: residue(:)
character(len=80)              :: line

! reset varible ...
system% atom(:)% translate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose keyword of stuff to TRANSLATE : '
write(*,'(/a)') ' (1) = tuning already done'
write(*,'(/a)') ' (2) = residue number '
write(*,'(/a)') ' (3) = residue name '
write(*,'(/a)') ' (4) = atom index '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )
    case( '1' ) 
        ! do nothing, proceed ...

    case( '2' )
        write(*,'(1x,3/a)') "enter the residue numbers to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        nr =  parse_this(line)
        do i = 1 , size(nr)
           where( system% atom(:)% nresid == nr(i) ) system% atom(:)% translate = .true.
        end do

    case( '3' )
        write(*,'(1x,3/a)') "enter the residue names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        residue = split_line(line)
        do i = 1 , size(residue)
           where( system% atom(:)% resid == residue(i) ) system% atom(:)% translate = .true.
        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        indx =  parse_this(line)

        system% atom(indx)% translate = .true.

end select

end subroutine translation_mask
!
!
!
!=================================
subroutine rotation_mask( system )
!=================================
implicit none 
type(universe) , intent(inout) :: system

!	local variables
integer           :: i
integer          , allocatable :: nr(:) , indx(:)
character(1)      :: choice
character(len=3) , allocatable :: residue(:)
character(len=80) :: line

! reset varible ...
system% atom(:)% rotate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose keyword of stuff to Rotate : '
write(*,'(/a)') ' (1) = tuning already done'
write(*,'(/a)') ' (2) = residue number '
write(*,'(/a)') ' (3) = residue name '
write(*,'(/a)') ' (4) = atom index '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )
    case( '1' ) 
        ! do nothing, proceed ...

    case( '2' )
        write(*,'(1x,3/a)') "enter the residue numbers to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        nr =  parse_this(line)
        do i = 1 , size(nr)
           where( system% atom(:)% nresid == nr(i) ) system% atom(:)% rotate = .true.
        end do

    case( '3' )
        write(*,'(1x,3/a)') "enter the residue names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        residue = split_line(line)
        do i = 1 , size(residue)
           where( system% atom(:)% resid == residue(i) ) system% atom(:)% rotate = .true.
        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        indx =  parse_this(line)

        system% atom(indx)% rotate = .true.

end select

end subroutine rotation_mask
!
!
!
!===============================
subroutine delete_mask( system )
!===============================
implicit none
type(universe) , intent(inout) :: system

!local variables
integer           :: i
character(len=1)  :: option
character(len=80) :: line
integer           , allocatable :: nr(:), indx(:)
character(len=1)  , allocatable :: fragment(:)
character(len=2)  , allocatable :: MMSymbol(:)
character(len=3)  , allocatable :: residue(:)

allocate( displace(system%N_of_atoms) , source = 0 ) 

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose keyword of stuff to Delete : '
write(*,'(/a)') ' (1) = MMSymbol'
write(*,'(/a)') ' (2) = fragment'
write(*,'(/a)') ' (3) = residue number '
write(*,'(/a)') ' (4) = residue name '
write(*,'(/a)') ' (5) = atom indices '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') option

select case( option )
    case( '1' ) 
        write(*,'(1x,3/a)') "enter the MMSymbol names to be DELETED, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        MMSymbol = split_line(line)
        do i = 1 , size(MMSymbol)
             where( system% atom(:)% MMSymbol == MMSymbol(i) ) system% atom(:)% delete = .true.
        end do

    case( '2' ) 
        write(*,'(1x,3/a)') "enter the fragment names to be DELETED, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        fragment = split_line(line)
        do i = 1 , size(fragment)
             where( system% atom(:)% fragment == fragment(i) ) system% atom(:)% delete = .true.
        end do

    case( '3' )
        write(*,'(1x,3/a)') "enter the residue numbers to be DELETED: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        nr =  parse_this(line)
        do i = 1 , size(nr)
             where( system% atom(:)% nresid == nr(i) ) system% atom(:)% delete = .true.
        end do

    case( '4' )
        write(*,'(1x,3/a)') "enter the residue names to be changed, separated by spaces (press ENTER to send) : "
        read (*,'(a)') line
        residue = split_line(line)
        do i = 1 , size(residue)
           where( system% atom(:)% resid == residue(i) ) system% atom(:)% delete = .true.
        end do

    case( '5' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed: separated by spaces, or in the format 'first:last' (press ENTER to send) : "
        read (*,'(a)') line
        indx =  parse_this(line)

        system% atom(indx)% delete = .true.

end select

CALL systemQQ( "clear" )

end subroutine delete_mask
!
!
!
!
!================================================
subroutine change_via_nr( system , nr , instance)
!================================================
implicit none 
type(universe)   , intent(inout) :: system
integer          , intent(in)    :: nr
character(len=*) , intent(in)    :: instance

! local variables ...

select case( instance )
      
       case("MMSymbol")
           where( system% atom(:)% nresid == nr ) system% atom(:)% MMSymbol = MMSymbol
       case("Fragment")                                                                 
           where( system% atom(:)% nresid == nr ) system% atom(:)% fragment = fragment
       case("resid")                                                                    
           where( system% atom(:)% nresid == nr ) system% atom(:)% resid = resid
       case("nresid")                                                                   
           where( system% atom(:)% nresid == nr ) system% atom(:)% nresid = nresid

end select

end subroutine change_via_nr
!
!
!
!
!==========================================================
subroutine change_via_residue( system , residue , instance)
!==========================================================
implicit none 
type(universe)   , intent(inout) :: system
character(len=*) , intent(in)    :: residue
character(len=*) , intent(in)    :: instance

! local variables ...

select case( instance )
      
       case("MMSymbol")
           where( system% atom(:)% resid == residue ) system% atom(:)% MMSymbol = MMSymbol
       case("Fragment")                                                                   
           where( system% atom(:)% resid == residue ) system% atom(:)% fragment = fragment
       case("resid")                                                                      
           where( system% atom(:)% resid == residue ) system% atom(:)% resid = resid
       case("nresid")                                                                     
           where( system% atom(:)% resid == residue ) system% atom(:)% nresid = nresid

end select

end subroutine change_via_residue
!
!
!
!
!============================================================
subroutine change_via_MM_name( system , MM_name , instance)
!============================================================
implicit none 
type(universe)   , intent(inout) :: system
character(len=*) , intent(in)    :: MM_name
character(len=*) , intent(in)    :: instance

! local variables ...

select case( instance )
      
       case("MMSymbol")
           where( system% atom(:)% MMSymbol == MM_name ) system% atom(:)% MMSymbol = MMSymbol
       case("Fragment")                                                                   
           where( system% atom(:)% MMSymbol == MM_name ) system% atom(:)% fragment = fragment
       case("resid")                                                                      
           where( system% atom(:)% MMSymbol == MM_name ) system% atom(:)% resid = resid
       case("nresid")                                                                     
           where( system% atom(:)% MMSymbol == MM_name ) system% atom(:)% nresid = nresid

end select

end subroutine change_via_MM_name
!
!
!
!
!=================================================
subroutine change_via_index( system , indx , instance)
!=================================================
implicit none 
type(universe)   , intent(inout) :: system
integer          , intent(in)    :: indx(:)
character(len=*) , intent(in)    :: instance

! local variables ...

select case( instance )
      
       case("MMSymbol")
           system% atom(indx)% MMSymbol = MMSymbol
       case("Fragment")                             
           system% atom(indx)% fragment = fragment
       case("resid")                                
           system% atom(indx)% resid = resid
       case("nresid")                               
           system% atom(indx)% nresid = nresid

end select

end subroutine change_via_index
!
!
!
!=====================================
function parse_this(line) result(indx)
!=====================================
implicit none 
character(len=*)        , intent(in)  :: line

! local variables ...
integer                        :: i , num
character(len=9) , allocatable :: tokens(:)
integer          , allocatable :: indx(:)

allocate( tokens , source = split_line(line) ) 

allocate(indx(0))

do i = 1 , size(tokens)
	if( scan(tokens(i),":") /= 0 ) &
	then
        indx = [ indx , seq_in_range(tokens(i)) ]
    else
        read(tokens(i),*) num 
        indx = [ indx , num ] 
    end if
end do

deallocate( tokens )

end function parse_this
!
!
!
!
!===============================================
subroutine display_menu( option )
!===============================================
implicit none
integer , intent(out) :: option

! local variables


!character(len=1) :: yn

!CALL systemQQ( "clear" )
!
!--------------------------------------------------------------
! invoke xterm terminal to call visualization tool ...
!write(*,'(/a)',advance='no') ">>> preview structue ? (y/n)  "
!read(*,'(a)') yn
!if( yn == "y") &
!then
!    CALL systemQQ( "nohup xterm &" )
!end if
!--------------------------------------------------------------

CALL systemQQ( "clear" )

write(*,'(/a)',advance='no') ">>> choose keyword for Edition Operation (1..9): "
print*, ""
write(*,'(/a)') 'Modify Symbol (1)'
write(*,'(/a)') 'Modify MMSymbol  (2)'
write(*,'(/a)') 'Modify Fragment (3)'
write(*,'(/a)') 'Modify residue name (4)'
write(*,'(/a)') 'residue number (5)'
write(*,'(/a)') 'COPY stuff (6)'
write(*,'(/a)') 'DELETE stuff (7)'
write(*,'(/a)') 'TRANSLATE stuff (8)'
write(*,'(/a)') 'ROTATE stuff (9)'

write(*,'(/a)',advance='no') '>>>   '                                                                                                                     
read (*,'(I)') option 

end subroutine display_menu
!
!
!
!
end module EDT_util_m

