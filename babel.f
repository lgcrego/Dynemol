 module Babel_m

    use type_m
    use Allocation_m

    PUBLIC :: Read_from_XYZ , Read_from_Poscar 
    PUBLIC :: Read_PDB , Read_VASP

    character(len=72) , PUBLIC :: System_Characteristics

    type atomic
        real*8                        :: xyz(3)
        real*8                        :: mass
        real*8                        :: charge
        integer                       :: AtNo
        integer                       :: nresid 
        integer                       :: copy_No
        character(3)                  :: residue
        character(3)                  :: Symbol
        character(3)                  :: MMSymbol
        character(1)                  :: TorF(3)
        character(1)                  :: fragment
    end type atomic

    type molecular
        type(atomic)    , allocatable :: atom(:) 
        real*8                        :: CG(3)
        real*8                        :: radius
        integer                       :: N_of_Atoms 
        integer                       :: nresid   
        integer                       :: copy_No
        character(3)                  :: residue 
        character(72)                 :: Solvent_Characteristics
    end type molecular

    type universe
        type(atomic)    , allocatable :: atom(:)
        type(molecular) , allocatable :: solvent(:)
        type(molecular)               :: dye
        real*8                        :: box(3)
        real*8                        :: Surface_Boundary
        integer                       :: N_of_Atoms
        integer                       :: N_of_Surface_Atoms
        integer                       :: N_of_Solvent_Atoms
        integer                       :: N_of_Solvent_Molecules
        character(1)    , allocatable :: list_of_fragments(:)
        character(3)    , allocatable :: list_of_residues(:)
        character(72)                 :: System_Characteristics
    end type universe

 contains
!
!
!
!===================================
 subroutine Read_from_XYZ(unit_cell)
!===================================

 type(structure) , intent(out) :: unit_cell

 character(len=2) :: dumb_char 
 character(len=1) :: fragment
 real*8           :: k_WH
 integer          :: i , j , j1 , j2 , n_kWH , n_fragments , N_of_Configurations

 OPEN(unit=3,file='xyz.dat',status='old')   

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell

! allocating arrays
 CALL Allocate_UnitCell(unit_cell)

! Acceptor (A) , Donor (D) or Molecule (M)
! fragment atoms must be packed together
 read(3,*) n_fragments
 do i = 1 , n_fragments
    read(3,*) j1 , j2 , fragment
    forall(j=j1:j2) unit_cell%fragment(j) = fragment
 end do

! defining the k_WH parameter for the atom 
! No of k_WH to be used 
 read(3,*) n_kWH                    
 do i = 1 , n_kWH
    read(3,*) j1 , j2 , k_WH
    forall(j=j1:j2) unit_cell%k_WH(j) = k_WH
 end do

! reading coordinates
 read(3,*) dumb_char
 DO j = 1 , unit_cell%atoms
     read(3,*) unit_cell%symbol(j), (unit_cell%coord(j,i), i=1,3)
 END DO   

 CLOSE(3)

! unit_cell dimensions
 unit_cell%T_xyz(1) = T_x
 unit_cell%T_xyz(2) = T_y
 unit_cell%T_xyz(3) = T_z

 print 70, System_Characteristics

 include 'formats.h'

 end subroutine Read_from_XYZ
!
!
!
!======================================
 subroutine Read_from_Poscar(unit_cell)
!======================================

 type(structure) , intent(out) :: unit_cell

 real*8            :: x0 , y0 , z0
 real*8            :: a , b , c
 real*8            :: kWH
 integer           :: n_kWH , n_fragments , i , j , j1 , j2 , n , boundary_atom
 integer           :: N_of_atoms , N_of_elements , N_of_cluster_atoms , N_of_Configurations
 character(len=1)  :: TorF , fragment
 logical           :: flag

 real*8           , allocatable :: xyz(:,:)
 integer          , allocatable :: atom_No(:)
 character(len=2) , allocatable :: element(:) , symbol(:)

 OPEN(unit=3,file='poscar.in',status='old')

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell
 
! allocating arrays
 CALL Allocate_UnitCell(unit_cell)

! Acceptor (A) , Donor (D) or Molecule (M)
! fragment atoms must be packed together
 read(3,*) n_fragments
 do i = 1 , n_fragments
    read(3,*) j1 , j2 , fragment
    forall(j=j1:j2) unit_cell%fragment(j) = fragment
 end do

! defining the k_WH parameter for the atom 
! No of k_WH to be used 
 read(3,*) n_kWH                    
 do i = 1 , n_kWH
    read(3,*) j1 , j2 , k_WH
    forall(j=j1:j2) unit_cell%k_WH(j) = k_WH
 end do

! reads the unit cell vectors
 read(3,*) x0 , y0 , z0
 a = dsqrt(x0*x0 + y0*y0 + z0*z0)
 read(3,*) x0 , y0 , z0
 b = dsqrt(x0*x0 + y0*y0 + z0*z0)
 read(3,*) x0 , y0 , z0
 c = dsqrt(x0*x0 + y0*y0 + z0*z0)

! reads the chemical elements
 read(3,*) N_of_elements 
 allocate(element(N_of_elements),atom_No(N_of_elements))
 read(3,*) (element(i),i=1,N_of_elements)

! reads the number of atoms for each species
 read(3,*) (atom_No(i),i=1,N_of_elements)

 N_of_atoms = sum(atom_No(1:N_of_elements))
 if( N_of_atoms /= unit_cell%atoms) pause 'Number of atoms do not match !!!'

 allocate(xyz(N_of_atoms,3),symbol(N_of_atoms))

! reads the atom that sets the interface 
 read(3,*) boundary_atom

! reads the coordinates from vasp-POSCAR
 do i = 1 , N_of_atoms

    read(3,*) xyz(i,1) , xyz(i,2) , xyz(i,3) , TorF , TorF , TorF
 
    do j = 1 , N_of_elements
        if( (sum(atom_No(:j-1))+1 <= i) .AND. (i <= sum(atom_No(:j))) ) symbol(i) = element(j)
    end do

 end do
 close(3)

! rescale coordinates by the unit cell vectors
 xyz(:,1) = xyz(:,1)*a
 xyz(:,2) = xyz(:,2)*b
 xyz(:,3) = xyz(:,3)*c

! separate cluster from molecule atoms
 n = 1
 do i = 1 , N_of_atoms
    flag = xyz(i,3) > xyz(boundary_atom,3) ! <== cluster atoms first
    if (flag) then
        unit_cell%symbol(n)  = symbol(i)
        unit_cell%coord(n,:) = xyz(i,:)
        n = n + 1
    end if
 end do
 N_of_cluster_atoms = n - 1

 do i = 1 , N_of_atoms
    flag = xyz(i,3) <= xyz(boundary_atom,3) ! <== molecule atoms at the end
    if (flag) then
        unit_cell%symbol(n)  = symbol(i)
        unit_cell%coord(n,:) = xyz(i,:)
        n = n + 1
    end if
 end do

 deallocate(xyz,symbol,element,atom_No)

 include 'formats.h'

 end subroutine Read_from_Poscar
!
!
!
!========================
 subroutine Read_PDB(trj)
!========================
type(universe)  , allocatable   , intent(out)   :: trj(:)

! local variables ...
integer                         :: openstatus , inputstatus , i , j , k , model , number_of_atoms , n , m
real*8                          :: time_1 , time_2 , delta_t 
character(4)                    :: keyword
character(1)                    :: test

open(unit = 31, file = 'frames.pdb', status = 'old', action = 'read', iostat = openstatus)
if (openstatus > 0) stop " *** Cannot open the file *** "

! find the number of model frames ...
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( inputstatus /= 0 ) exit
        
    if ( keyword == 'info' ) then
        backspace 31
        read(unit = 31, fmt = 43, iostat = inputstatus) System_Characteristics
    end if
    if ( keyword == 'MODE' ) then
        backspace 31
        read(unit = 31, fmt = 36, iostat = inputstatus) model ! <==
    end if
end do

! return to the top of the file ...
rewind 31

! read number the atoms and time ...
read(unit = 31, fmt = 35, iostat = inputstatus) keyword
do
    if ( keyword == 'MODE' ) then
        exit
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_1 
                    exit    ! <== time_1 read
                end if
            end do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        else
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        end if
    end if
end do

do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'ATOM' ) then
        backspace 31
        read(unit = 31, fmt = 32, iostat = inputstatus) number_of_atoms ! <==
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_2 
                    exit    ! <== time_2 read
                end if
            end do
            exit    ! <== leave outer do loop
        end if
    end if
end do

delta_t = time_2 - time_1

! return to the top of the file ...
rewind 31

! allocate array of frames ...
allocate( trj(model) )

trj(:)%System_Characteristics = System_Characteristics

do j = 1 , model
    if( j == 1 ) then
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 33, iostat = inputstatus)     &
            trj(j)%atom(i)%MMSymbol , trj(j)%atom(i)%residue , trj(j)%atom(i)%nresid , ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
        CALL MMSymbol_2_Symbol( trj(j)%atom )
        CALL Symbol_2_AtNo( trj(j)%atom )
        CALL Initial_Setting_of_Fragment( trj(j) )
    else
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if ( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 37, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if
end do

close(31)

trj%N_of_atoms = number_of_atoms

! get list of residues in trj ...
CALL Identify_Residues( trj(1) )

! get list of fragments in trj ...
CALL Identify_Fragments( trj(1) )

! Copy information from trj(1) to trj(:) ...
forall(i = 2:model )
    trj(i)%atom%AtNo         =  trj(1)%atom%AtNo    
    trj(i)%atom%MMSymbol     =  trj(1)%atom%MMSymbol
    trj(i)%atom%residue      =  trj(1)%atom%residue
    trj(i)%atom%nresid       =  trj(1)%atom%nresid
    trj(i)%atom%Symbol       =  trj(1)%atom%Symbol
    trj(i)%atom%fragment     =  trj(1)%atom%fragment
end forall

! SORT and GROUP residues ...
do i = 1 , size(trj)
    CALL Sort_Residues( trj(i) )
end do

! Information about the Solvent System ...
trj%N_of_Solvent_Molecules =  maxval( trj(1) % atom % nresid, trj(1) % atom % fragment == 'S' )         &
                            - minval( trj(1) % atom % nresid, trj(1) % atom % fragment == 'S' ) + 1

do i = 1 , size(trj)
    allocate( trj(i)%solvent(trj(i)%N_of_Solvent_Molecules) )

    trj(i)%solvent%N_of_Atoms = count( trj(i)%atom%fragment == 'S' ) / trj(i)%N_of_Solvent_molecules

    CALL Center_of_Gravity( trj(i) )
end do

! Formats ...
32 format(5x, i6)
33 format(13x, a3, t18, a3, t24, i3, t33, f6.3, t41, f6.3, t49, f6.3)
34 format(2x, a1, t10, a1, t18, i2, t27, f10.6, t42, f10.6, t59, f10.6)
35 format(a4)
36 format(7x, i7)
37 format(32x, f6.3, t41, f6.3, t49, f6.3)
38 format(10x, a1)
39 format(81x, f7.0)
40 format(6x, 3f9.3)
41 format(f10.5)
42 format(a1)
43 format(10x,a72)

end subroutine Read_PDB
!
!
!=========================
 subroutine Read_VASP(trj)
!=========================
type(universe)  , allocatable   , intent(out) :: trj(:)

! local variables ....
character(1)                    :: idx
real*8          , allocatable   :: distance_ligation(:,:) , distance_T(:)
integer                         :: openstatus , inputstatus , atoms , i , j , k , ligation , indx1 , indx2

open(unit = 13, file = 'VASP.trj', status = 'old', action = 'read', iostat = openstatus)
if( openstatus > 0 ) stop '*** Cannot open the file ***'

! read the number of models ...
model = 0
do
    read(unit = 13, fmt = 21, iostat = inputstatus) idx
    if( inputstatus /= 0 ) exit
    if( idx == '$' ) then
        model = model + 1
    end if
end do

model = model - 1

rewind 13

! read the number the atoms ...
read(unit = 13, fmt = 21, iostat = inputstatus) idx
atoms = 0
do
    read(unit = 13, fmt = 21, iostat = inputstatus) idx
    if( idx /= '$' ) then
        atoms = atoms + 1
    else
        exit
    end if
end do

rewind 13

! read the system ...
allocate( trj(model) )

do j = 1 , model
    if( j == 1 ) then
        read(unit = 13, fmt = 21, iostat = inputstatus) idx

        allocate( trj(j)%atom(atoms) )

        do i = 1 , atoms
            read(unit = 13, fmt = 22, iostat = inputstatus) trj(j)%atom(i)%Atno, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
        CALL AtNo_2_Symbol(trj(j)%atom)
    else
        read(unit = 13, fmt = 21, iostat = inputstatus) idx

        allocate( trj(j)%atom(atoms) )

        do i = 1 , atoms
            read(unit = 13, fmt = 23, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if
end do

forall(i=2:model)
    trj(i)%atom%Atno = trj(1)%atom%Atno
    trj(i)%atom%Symbol = trj(1)%atom%Symbol
end forall

! formats ...
21 format(a1)
22 format(i2, t4, f8.5, t13, f8.5, t22, f8.5)
23 format(3x, f8.5, t13, f8.5, t22, f8.5)

end subroutine Read_VASP
!
!
!
!===========================
 subroutine Symbol_2_AtNo(a)
!===========================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(a(i)%symbol)
        case( 'H') 
            a(i)%AtNo = 1 
        case( 'C') 
            a(i)%AtNo = 6 
        case( 'N') 
            a(i)%AtNo = 7 
        case( 'O') 
            a(i)%AtNo = 8 
        case( 'F') 
            a(i)%AtNo = 9 
        case( 'AL','Al') 
            a(i)%AtNo = 13 
        case( 'S','s') 
            a(i)%AtNo = 16 
        case( 'CL','Cl') 
            a(i)%AtNo = 17 
        case( 'TI','Ti') 
            a(i)%AtNo = 22 
        case( 'MN','Mn') 
            a(i)%AtNo = 25 
        case default
            print*, ' >> unkown atom found ; execution terminated <<' 
            stop
    end select

 END DO

 end subroutine Symbol_2_AtNo
!
!
!
!===========================
 subroutine AtNo_2_Symbol(a)
!===========================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(a(i)%Atno)
        case( 1) 
            a(i)%Symbol = 'H'
        case( 6) 
            a(i)%Symbol = 'C'
        case( 7) 
            a(i)%Symbol = 'N'
        case( 8) 
            a(i)%Symbol = 'O'
        case( 9) 
            a(i)%Symbol = 'F'
        case( 13) 
            a(i)%Symbol = 'Al'
        case( 16) 
            a(i)%Symbol = 'S '
        case( 17) 
            a(i)%Symbol = 'Cl '
        case( 22) 
            a(i)%Symbol = 'Ti '
        case( 25) 
            a(i)%Symbol = 'Mn'
        case default
            print*, ' >> unkown atom found ; execution terminated <<' 
            stop
    end select

 END DO

 end subroutine AtNo_2_Symbol
!
!
!
!===============================
 subroutine MMSymbol_2_Symbol(a)
!===============================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer             :: i
character(len=1)    :: element

 DO i = 1 , size(a)

    write( element,'(A1)' ) adjustl( a(i)%MMSymbol )

    select case( element )
        case( 'C' ) 
            a(i)%Symbol = 'C' 
        case( 'N' ) 
            a(i)%Symbol = 'N' 
        case( 'O' ) 
            a(i)%Symbol = 'O' 
        case( 'H' ) 
            a(i)%Symbol = 'H' 
    end select

    select case( adjustl( a(i)%MMSymbol) )
        case( 'YC' ) 
            a(i)%Symbol = 'C' 
        case( 'YN' ) 
            a(i)%Symbol = 'N' 
        case( 'Al' ) 
            a(i)%Symbol = 'Al' 
        case( 'Ti' ) 
            a(i)%Symbol = 'Ti' 
    end select

 END DO

end subroutine MMSymbol_2_Symbol
!
!
!
!==========================================
subroutine Initial_Setting_of_Fragment( a )
!==========================================
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
!--------------------------------------------

 DO i = 1 , size(a%atom)

    select case(a%atom(i)%residue)
        case( 'CCC') 
            a%atom(i)%fragment = 'C' 
        case( 'Alq') 
            a%atom(i)%fragment = 'M' 
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
    end select

 END DO

end subroutine Initial_Setting_of_Fragment
!
!
!
!=============================================
subroutine Identify_Fragments( a )
!=============================================
implicit none
type(universe)  , intent(inout) :: a

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag

allocate( temp(a%N_of_Atoms) )

temp(1) = a % atom(1) % fragment
counter = 1

do i = 1 , a%N_of_Atoms

    flag = .true.
    do j = 1 , counter
        flag = flag .AND. (temp(j) /= a%atom(i)%fragment)
    end do

    if( flag ) then
        counter = counter + 1
        temp(counter) = a%atom(i)%fragment
    end if

end do

! build list of fragments in a ...
allocate( a%list_of_fragments(counter) )
a%list_of_fragments = temp(1:counter)
deallocate( temp )

end subroutine Identify_Fragments
!
!
!
!=============================================
subroutine Identify_Residues( a )
!=============================================
implicit none
type(universe)  , intent(inout) :: a

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag

allocate( temp(a%N_of_Atoms) )

temp(1) = a % atom(1) % residue
counter = 1

do i = 1 , a%N_of_Atoms

    flag = .true.
    do j = 1 , counter
        flag = flag .AND. (temp(j) /= a%atom(i)%residue)
    end do

    if( flag ) then
        counter = counter + 1
        temp(counter) = a%atom(i)%residue
    end if

end do

! build list of residues in a ...
allocate( a%list_of_residues(counter) )
a%list_of_residues = temp(1:counter)
deallocate( temp )

end subroutine Identify_Residues
!
!
!
!===============================
subroutine Sort_Residues(system)
!===============================
implicit none
type(universe) , intent(inout) :: system

! local variables ...
integer        :: i , j , N_of_elements , N_of_atoms ,  iptr
type(universe) :: temp

allocate( temp%atom(1) )

N_of_atoms = size(system%atom)

do i = 1 , N_of_atoms-1

    iptr = i

    do j = i+1 , N_of_atoms
        if( LGT( system % atom(j) % residue , system % atom(iptr) % residue ) ) then
            iptr = j
        end if
    end do

    if( i /= iptr ) then
        temp   % atom (1)    = system % atom (i)
        system % atom (i)    = system % atom (iptr)
        system % atom (iptr) = temp   % atom (1)
    end if

end do

deallocate( temp%atom )

end subroutine Sort_Residues
!
!
!
!=====================================
subroutine Center_of_Gravity( trj )
!=====================================
implicit none
type(universe) , intent(inout) :: trj

! local variables ...
integer :: i , j , n , mol_atoms

! initial position of S fragments in trj%atom array ...
n = minloc( trj%atom%fragment , 1 , trj%atom%fragment == "S" ) 

! center of gravity ...
do i = 1 , trj%N_of_Solvent_Molecules 

    mol_atoms = trj%solvent(i)%N_of_Atoms 

    forall( j=1:3 ) trj%solvent(i)%CG(j) = sum( trj%atom(n:n+mol_atoms)%xyz(j) ) / trj%solvent(i)%N_of_Atoms

    trj % solvent(i) % nresid  = trj % atom(n) % nresid
    trj % solvent(i) % residue = trj % atom(n) % residue

    n = n + mol_atoms

end do

end subroutine Center_of_Gravity

end module Babel_m

