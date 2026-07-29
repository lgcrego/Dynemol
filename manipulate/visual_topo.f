module visual_topology

use constants_m
use ansi_colors
use types_m            , only : dynemolworkdir, universe
use util_m             , only : read_file_name, TO_UPPER_CASE
use Topology_routines  , only : dump_topol

public :: bond_check

private
 
    ! module variables ...
    integer, allocatable :: bond_pairs(:,:)

contains
!
!
!
!=============================
 subroutine bond_check(system)
!=============================
implicit none
type(universe), intent(inout) :: system

! local variables ...
character(len=1)  :: option
character(len=30) :: f_name

CALL systemQQ( "clear" ) 
                                                                                                                                                              
write(*,'(/a)') bold // cyan // ' Select topology file format:' // reset
                                                                                                                                                              
write(*,'(a)') green // ' (1) ' // reset // '= psf'                                                                                                           
write(*,'(a)') green // ' (2) ' // reset // '= itp'                                                                                                          
                                                                                                                                                              
write(*,'(/a)', advance='no') bold // yellow // '>>> ' // reset                                                                                               
read (*,'(a)') option                                                                                                                                    
                                                                                                                                                              
select case( option )                                                                                                                                    
                                                                                                                                                              
    case( '1' )                                                                                                                                               
        CALL read_file_name( f_name , file_type="psf" )                                                                                                       
        CALL psf_file_reader (f_name )                                                                                                            
                                                                                                                                                              
    case( '2' )                                                                                                                                               
        CALL read_file_name( f_name , file_type="itp" )                                                                                                       
        CALL itp_file_reader (f_name )                                                                                                            

    case default
        write(*,'(/a)') red // 'Invalid option.' // reset

end select                

system% topol = build_topo_mtx(system)

call save_pdb_with_conect( system )

call systemQQ("vmd -e .load.tcl")

call systemQQ("rm .load.tcl")

end subroutine bond_check
!
!
!
subroutine psf_file_reader (f_name)
implicit none
character(len=*), intent(in):: f_name

!local variables
integer :: f_unit, k, j, n, Nbonds, ioerr
character(len=60) :: line

open(newunit=f_unit, file=dynemolworkdir//f_name, status='old',action='read', iostat=ioerr)

if (ioerr /= 0) then
    write(*,*) trim(f_name),' not found.'
    stop
end if

!-----------------------------------------
! Locate !NBOND section
!-----------------------------------------
do
    read(f_unit,'(A)',iostat=ioerr) line
    if (ioerr /= 0) stop 'Could not locate !NBOND section.'

    line = to_upper_case(line)
    if( verify( "!NBOND" , line ) == 0 ) exit
end do

backspace(f_unit)

read(f_unit,*) Nbonds

if( Nbonds > 0 ) then
    allocate( bond_pairs (Nbonds, 2) , source = 0 ) 
    
    do k = 1, ceiling(Nbonds/four)-1
        read(f_unit, *)  ( ( bond_pairs((k-1)*4+n,j) , j=1,2 ) , n=1,4 )
    end do 
    read(f_unit, *)  ( ( bond_pairs((k-1)*4+n,j) , j=1,2 ) , n=1,merge(4,mod(NBonds,4),mod(NBonds,4)==0) )
end if

close(f_unit)

end subroutine psf_file_reader
!
!
!
subroutine itp_file_reader (f_name)
implicit none
character(len=*), intent(in):: f_name

!local variables
integer :: f_unit, n, Nbonds, ioerr
character(len=60) :: line

open(newunit=f_unit, file=dynemolworkdir//f_name, status='old', action='read', iostat=ioerr)

if (ioerr /= 0) then
    write(*,*) trim(f_name),' not found.'
    stop
end if

!-----------------------------------------
! Locate [ bonds ] section
!-----------------------------------------
do
    read(f_unit,'(A)',iostat=ioerr) line
    if (ioerr /= 0) stop 'Could not locate [ bonds ] section.'

    if (index(to_upper_case(line),'[ BONDS ]') > 0) exit
end do

!-----------------------------------------
! Count bond records
!-----------------------------------------
Nbonds = 0
do
    read(f_unit,'(A)',iostat=ioerr) line
    if (ioerr /= 0) exit

    line = adjustl(line)

    if (line == "") cycle
    if (line(1:1) == ";") cycle
    if (line(1:1) == "[") exit

    Nbonds = Nbonds + 1
end do

!-----------------------------------------
! Return to [ bonds ]
!-----------------------------------------
if( Nbonds > 0 ) then

    rewind(f_unit)
    allocate(bond_pairs(Nbonds,2))

    do
        read(f_unit,'(A)') line
        if (index(to_upper_case(line),'[ BONDS ]') > 0) exit
    end do
    
    n = 0
    do
        read(f_unit,'(A)', iostat=ioerr) line
        if (ioerr /= 0) exit
    
        line = adjustl(line)
    
        if (line == "") cycle
        if (line(1:1) == ";") cycle
        if (line(1:1) == "[") exit
    
        n = n + 1
    
        read(line,*) bond_pairs(n,1), bond_pairs(n,2)
    end do
endif

close(f_unit)

do n = 1, Nbonds
    print*, n, bond_pairs(n,1), bond_pairs(n,2)
end do

end subroutine itp_file_reader
!
!
!============================================
function build_topo_mtx(system) result(topo)
!============================================
implicit none 
type(universe), intent(in) :: system

! local variables ...
integer :: i, j, k, n_atoms
logical, allocatable :: topo(:,:)

n_atoms = system% N_of_atoms

if( allocated(bond_pairs) ) then
    allocate( topo(n_atoms,n_atoms), source = .false. )
    do k = 1, size(bond_pairs, dim=1)
        i = bond_pairs(k,1)
        j = bond_pairs(k,2)
        topo(i,j) = .true.
        topo(j,i) = .true.
    end do
end if

end function build_topo_mtx
!
!
!============================================
subroutine save_pdb_with_conect(sys, string)
!============================================
implicit none 
type(universe)          , intent(inout) :: sys
character(*)  , optional, intent(in)    :: string

! local variables ...
integer :: i, k, f_unit
character(len=:), allocatable :: file_name

!----------------------------------------------
!    generate pdb file with conect record
!----------------------------------------------
file_name = "seed.pdb"

If( present(string) ) then
    file_name = string
end if

OPEN(newunit=f_unit, file=file_name, status='unknown')

write(f_unit,6) sys%System_Characteristics
write(f_unit,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms

            write(f_unit,2)  'ATOM  '                   ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        sys%atom(i)%MMSymbol            ,  &    ! <== atom type
                        ' '                             ,  &    ! <== alternate location indicator
                        sys%atom(i)%resid               ,  &    ! <== residue name
                        ' '                             ,  &    ! <== chain identifier
                        sys%atom(i)%nresid              ,  &    ! <== residue sequence number
                        ' '                             ,  &    ! <== code for insertion of residues
                        ( sys%atom(i)%xyz(k) , k=1,3 )  ,  &    ! <== xyz coordinates 
                        1.00                            ,  &    ! <== occupancy
                        0.00                            ,  &    ! <== temperature factor
                        ' '                             ,  &    ! <== segment identifier
                        ' '                             ,  &    ! <== here only for tabulation purposes
                        sys%atom(i)%symbol              ,  &    ! <== chemical element symbol
                        sys%atom(i)%charge                      ! <== charge on the atom
end do

! check and print topological conections ...
If( allocated(sys%topol) ) CALL dump_topol(sys,f_unit)

write(f_unit,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , sys%N_of_atoms , 0 , sys%N_of_atoms , 0
write(f_unit,*) 'END'

close(f_unit)

call write_vmd_load_script(file_name)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a72)

end subroutine save_pdb_with_conect
!
!
!
!=========================================
subroutine write_vmd_load_script(pdb_file)
!=========================================
implicit none

character(len=*), intent(in) :: pdb_file

integer :: unit

open(newunit=unit, file=".load.tcl", status="replace", action="write")

write(unit,'(a)') "# load.tcl"
write(unit,'(a)') "package require topotools"

write(unit,'(a)') "mol new " // trim(pdb_file) // " autobonds off waitfor all"

write(unit,'(a)') ""
write(unit,'(a)') "# Apply CPK representation"
write(unit,'(a)') "mol representation CPK 0.500000 0.300000 10.000000 8.000000"
write(unit,'(a)') "mol color Name"
write(unit,'(a)') "mol selection {all}"
write(unit,'(a)') "mol material Opaque"
write(unit,'(a)') "mol addrep top"

close(unit)

end subroutine write_vmd_load_script
!
!
!
end module visual_topology
