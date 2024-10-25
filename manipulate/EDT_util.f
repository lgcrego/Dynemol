module EDT_util_m

use types_m
use util_m
use Read_Parms
use Constants_m
use GMX_routines

public :: on_the_fly_tuning , parse_this

private

        !module variables ...
        integer          :: nresid
        character(len=4) :: resid
        character(len=4) :: MMSymbol
        character(len=1) :: fragment

contains

!
!===============================================
subroutine on_the_fly_tuning( system , displace)
!===============================================
implicit none
type(universe)                   , intent(inout) :: system
integer , allocatable , optional , intent(out)   :: displace(:)

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
        read (*,'(a)') MMSymbol
        CALL change_this( system , instance="fragment" )

    case (4)
        write(*,'(/a)',advance='no') ">>> enter new residue name: "
        read (*,'(a)') resid
        CALL change_this( system , instance="resid" )

    case (5)
        write(*,'(/a)',advance='no') ">>> enter new residue number: "
        read (*,'(a)') nresid
        CALL change_this( system , instance="nresid" )

    case (6)
        Write(*,*) "not implemented"

    case (7)
        CALL delete_mask( system , displace )

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
!==========================================
subroutine change_this( system , instance )
!==========================================
implicit none 
type(universe)   , intent(inout) :: system
character(len=*) , intent(in)    :: instance

! local variables ...
integer                        :: nr 
integer          , allocatable :: indx(:)
character(len=1)               :: choice
character(len=3)               :: residue
character(len=80)              :: line

! reset varible ...
system% atom(:)% translate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose the feature that controls the operation : '
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
        write(*,'(/a)') "choose the residue numbers to be changed (0 to finish) : "

        do
           read*, nr
           If( nr == 0 ) exit
           CALL change_nr( system , nr , instance )
        end do

    case( '3' )
        write(*,'(1x,3/a)') "choose the names of the residues to be changed (press ENTER after each residue; use @ to finish) : "

        do
           read*, residue
           If( residue == "@" ) exit
           CALL change_residue( system , residue , instance )
        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be changed separated by spaces (press ENTER) : "
        read (*,'(a)') line

        indx =  parse_this(line)

        CALL change_index( system , indx , instance )

end select

end subroutine change_this
!
!
!
!
!
!
!====================================
subroutine translation_mask( system )
!====================================
implicit none 
type(universe)              , intent(inout) :: system

! local variables ...
integer                        :: nr
integer          , allocatable :: indx(:)
character(len=1)               :: choice
character(len=3)               :: residue
character(len=80)              :: line

! reset varible ...
system% atom(:)% translate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose stuff to Copy : '
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
        write(*,'(/a)') "choose the residue numbers to be translated (0 to finish) : "

        do
           read*, nr
           If( nr == 0 ) exit
           
           where( system% atom(:)% nresid == nr ) system% atom(:)% translate = .true.

        end do

    case( '3' )
        write(*,'(1x,3/a)') "choose the names of the residues to be translated (press ENTER after each residue; use @ to finish) : "

        do
           read*, residue
           If( residue == "@" ) exit

           where( system% atom(:)% resid == residue ) system% atom(:)% translate = .true.

        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be translated separated by spaces (press ENTER) : "
        read (*,'(a)') line

        indx = parse_this(line)

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
integer           :: nr
integer           , allocatable :: indx(:)
character(1)      :: choice
character(len=3)  :: residue
character(len=80) :: line

! reset varible ...
system% atom(:)% rotate = .false.

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose stuff to Rotate : '
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
        write(*,'(/a)') "choose the residue numbers to be rotated (0 to finish) : "

        do
           read*, nr
           If( nr == 0 ) exit
           
           where( system% atom(:)% nresid == nr ) system% atom(:)% rotate = .true.

        end do

    case( '3' )
        write(*,'(1x,3/a)') "choose the names of the residues to be rotated (press ENTER after each residue; use @ to finish) : "

        do
           read*, residue
           If( residue == "@" ) exit

           if( .not. any(system% atom(:)% resid == residue) ) Print*, ">> no residue ", residue ," found"

           where( system% atom(:)% resid == residue ) system% atom(:)% rotate = .true.

        end do
 
    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be rotated separated by spaces (press ENTER) : "
        read (*,'(a)') line

        indx = parse_this(line)

        system% atom(indx)% rotate = .true.

end select

end subroutine rotation_mask
!
!
!
!==========================================
subroutine delete_mask( system , displace )
!==========================================
implicit none
type(universe)               , intent(inout) :: system
integer        , allocatable , intent(out)   ::  displace(:)

!local variables
integer           :: nr , begin_of_block
character(len=1)  :: option
character(len=3)  :: residue
character(len=80) :: line
integer           , allocatable :: nr_mask(:) ,indx(:)
character(len=3)  , allocatable :: residue_mask(:)

allocate( displace(system%N_of_atoms) , source = 0 ) 

CALL systemQQ( "clear" )

write(*,'(/a)') ' Choose stuff to Delete : '
write(*,'(/a)') ' (1) = tuning already done'
write(*,'(/a)') ' (2) = residue number '
write(*,'(/a)') ' (3) = residue name '
write(*,'(/a)') ' (4) = atom indices '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') option

select case( option )
    case( '1' ) 
        ! do nothing, proceed ...

    case( '2' )
        write(*,'(/a)') "choose the residue numbers to be deleted (0 to finish) : "

        allocate ( nr_mask(system%N_of_atoms) )
        nr_mask(:) = system%atom(:)%nresid 
        do
           read*, nr
           If( nr == 0 ) exit
           
           where( system% atom(:)% nresid == nr ) system% atom(:)% delete = .true.

           begin_of_block = minloc( nr_mask , mask = nr_mask==nr , dim=1 )
           displace(begin_of_block:) = displace(begin_of_block:) + 1
        end do
        deallocate( nr_mask )

    case( '3' )
        write(*,'(1x,3/a)') "choose the name of the residue to be deleted (@ to finish) : "

        allocate ( residue_mask(system%N_of_atoms) ) 
        residue_mask(:) = system%atom(:)%resid 
        do
           read*, residue
           If( residue == "@" ) exit

           where( system% atom(:)% resid == residue ) system% atom(:)% delete = .true.

           where( system% atom(:)% resid == residue ) system% atom(:)% nresid = -1

           begin_of_block = minloc( residue_mask , mask = residue_mask==residue , dim=1 )
           displace(begin_of_block:) = displace(begin_of_block:) + 1
        end do
        deallocate( residue_mask )

    case( '4' )
        write(*,'(1x,3/a)') "enter the indices of the atoms to be deleted separated by spaces (press ENTER) : "
        read (*,'(a)') line

        indx = parse_this(line)

        system% atom(indx)% delete = .true.

end select


end subroutine delete_mask
!
!
!
!
!============================================
subroutine change_nr( system , nr , instance)
!============================================
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

end subroutine change_nr
!
!
!
!
!======================================================
subroutine change_residue( system , residue , instance)
!======================================================
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

end subroutine change_residue
!
!
!
!
!=================================================
subroutine change_index( system , indx , instance)
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

end subroutine change_index
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
character(len=1) :: yn

CALL systemQQ( "clear" )

!--------------------------------------------------------------
! invoke xterm terminal to call visualization tool ...
write(*,'(/a)',advance='no') ">>> preview structue ? (y/n)  "
read(*,'(a)') yn
if( yn == "y") &
then
    CALL systemQQ( "nohup xterm &" )
end if
!--------------------------------------------------------------

CALL systemQQ( "clear" )

write(*,'(/a)',advance='no') ">>> choose field to EDIT (1..9): "
print*, ""
write(*,'(/a)') 'Symbol (1)'
write(*,'(/a)') 'MMSymbol  (2)'
write(*,'(/a)') 'Fragment (3)'
write(*,'(/a)') 'resid (4)'
write(*,'(/a)') 'nresid (5)'
write(*,'(/a)') 'copy (6)'
write(*,'(/a)') 'delete (7)'
write(*,'(/a)') 'translate (8)'
write(*,'(/a)') 'rotate (9)'

read (*,*) option

end subroutine display_menu
!
!
!
!
end module EDT_util_m

