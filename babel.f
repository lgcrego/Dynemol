 module Babel_m

    use type_m                  
    use Allocation_m            , only : Allocate_UnitCell
    use Semi_Empirical_Parms    , only : Define_EH_Parametrization

    PUBLIC :: Read_from_XYZ , Read_from_Poscar , Read_from_PDB
    PUBLIC :: Read_PDB , Read_VASP
    PUBLIC :: Coords_from_Universe

    type(universe)  , allocatable , public  :: trj(:)

    character(len=72) , PUBLIC :: System_Characteristics

    interface Symbol_2_AtNo
        module procedure Sym_2_AtNo_TRJ
        module procedure Sym_2_AtNo_XYZ
    end interface

    interface Identify_Fragments
        module procedure Identify_Fragments_Structure
        module procedure Identify_Fragments_Universe
    end interface

contains
!
!
!
!=====================================
 subroutine Read_from_XYZ( Unit_Cell )
!=====================================

 type(structure) , intent(out) :: unit_cell

 character(len=1) :: fragment
 character(len=2) :: dumb_char 
 character(len=3) :: residue
 real*8           :: k_WH
 integer          :: i , j , j1 , j2 , n_kWH 
 integer          :: n_residues , N_of_Configurations

 OPEN(unit=3,file='xyz.dat',status='old')   

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell
 read(3,*) n_residues

! allocating arrays ...
 CALL Allocate_UnitCell( Unit_Cell , n_residues )

! residues must be packed together ...
 do i = 1 , n_residues
    read(3,*) j1 , j2 , residue , fragment

    unit_cell % list_of_residues  (i) = residue
    unit_cell % list_of_fragments (i) = residue

    forall(j=j1:j2) unit_cell % residue  (j)  = residue
    forall(j=j1:j2) unit_cell % fragment (j) = fragment
 end do

! defining the k_WH parameter for the atom 
! No of k_WH to be used 
 read(3,*) n_kWH                    
 do i = 1 , n_kWH
    read(3,*) j1 , j2 , k_WH
    forall(j=j1:j2) unit_cell%k_WH(j) = k_WH
 end do

! read unit_cell dimensions ...
 read(3,*) unit_cell%T_xyz(1) 
 read(3,*) unit_cell%T_xyz(2)
 read(3,*) unit_cell%T_xyz(3)

! reading coordinates ...
 read(3,*) dumb_char
 DO j = 1 , unit_cell%atoms
     read(3,*) unit_cell%symbol(j), (unit_cell%coord(j,i), i=1,3)
 END DO   

 CLOSE(3)

 CALL Symbol_2_AtNo(Unit_Cell)

 include 'formats.h'

 end subroutine Read_from_XYZ
!
!
!
!=====================================
 subroutine Read_from_PDB( Unit_Cell )
!=====================================
 implicit none 
 type(structure)   , intent(out) :: Unit_Cell

! local variables ...
integer             :: i , j , ioerr , N_of_atoms , nresidue
character(len=5)    :: keyword
type(universe)      :: system        

OPEN(unit=3,file='input.pdb',status='old',iostat=ioerr,err=11)
read(3,'(a72)') System_Characteristics

! reading unit cell vectors ...
read(unit=3,fmt=105,iostat=ioerr) keyword
if ( keyword == "CRYST" ) then
    backspace 3
    read(3,fmt=100) system%box(1) , system%box(2) , system%box(3)
end if
    
! scan file for N_of_Atoms ...   
N_of_atoms = 0
do
    read(unit=3,fmt=105,iostat=ioerr) keyword
    if( keyword == "HETAT" ) then
        N_of_atoms = N_of_atoms + 1
    end if
    if ( keyword == "MASTE" ) exit
end do
system%N_of_atoms = N_of_atoms
        
allocate( system%atom(system%N_of_atoms) )

!read data ...    
rewind 3
do
    read(unit=3,fmt=105,iostat = ioerr) keyword
    if( keyword == "HETAT" ) then
        backspace 3
        do i = 1 , system%N_of_atoms
            read(3,115)  system%atom(i)%MMSymbol            ,  &    ! <== atom type
                         system%atom(i)%residue             ,  &    ! <== residue name
                         nresidue                           ,  &    ! <== residue sequence number
                         (system%atom(i)%xyz(j) , j=1,3)    ,  &    ! <== xyz coordinates 
                         system%atom(i)%Symbol                      ! <== chemical element symbol
        end do
    end if
    if ( keyword == "MASTE" ) exit
end do

close(3)

! convert residues to upper case ...
forall( i=1:system%N_of_atoms ) system%atom(i)%residue = TO_UPPER_CASE( system%atom(i)%residue )

! preprocessing the universe system ...
CALL Symbol_2_AtNo      ( system%atom )
CALL Identify_Residues  ( system      )
CALL Setting_Fragments  ( system      )
CALL Identify_Fragments ( system      )
CALL Sort_Residues      ( system      )

! transfer structure <-- universe 
CALL Coords_from_Universe( Unit_Cell , system )

deallocate( system%atom , system%list_of_fragments , system%list_of_residues )
11 if( ioerr > 0 ) stop "input.pdb file not found; terminating execution"

100 format(t10, f6.3, t19, f6.3, t28, f6.3)
105 format(a5)
115 FORMAT(t12,a5,t18,a3,t23,i4,t31,f8.3,t39,f8.3,t47,f8.3,t77,a2)

end subroutine Read_from_PDB
!
!
!
!=====================================================
 subroutine Coords_from_Universe( Unit_Cell , System )
!=====================================================
 implicit none
 type(structure) , intent(out)    :: Unit_Cell
 type(universe)  , intent(inout)  :: System

! local variables ... 
integer         :: j , n_residues
character(72)   :: Characteristics

Unit_Cell%atoms = System%N_of_Atoms
select case( file_type )
    case( "structure" )
        n_residues = size( System%list_of_residues )
    case( "trajectory" )
        n_residues = size( trj(1)%list_of_residues )
    end select

! allocating Unit_Cell structure ...
CALL Allocate_UnitCell( Unit_Cell , n_residues )

! coordinates and other info from input data ...
forall( j=1:unit_cell%atoms )
    unit_cell % coord    (j,:) =  System % atom(j) % xyz(:)
    unit_cell % AtNo     (j)   =  System % atom(j) % AtNo  
    unit_cell % fragment (j)   =  System % atom(j) % fragment
    unit_cell % Symbol   (j)   =  System % atom(j) % Symbol
    unit_cell % residue  (j)   =  System % atom(j) % residue
    unit_cell % MMSymbol (j)   =  System % atom(j) % MMSymbol
end forall

! get list of residues ...
select case( file_type )
    case( "structure" )
        unit_cell%list_of_residues  = System%list_of_residues
        unit_cell%list_of_fragments = System%list_of_fragments
    case( "trajectory" )
        unit_cell%list_of_residues  = trj(1)%list_of_residues
        unit_cell%list_of_fragments = trj(1)%list_of_fragments
    end select

! unit_cell dimensions ...
unit_cell % T_xyz =  System % box

! get EH Parms ...
CALL Define_EH_Parametrization( Unit_Cell , Characteristics )

! verify consistency ...
if( System_Characteristics /= Characteristics ) then
    print*, System_Characteristics , Characteristics 
    Pause ">>> Input Parameter Inconsistency <<<"
end if

end subroutine Coords_from_Universe
!
!
!
!========================================
 subroutine Read_from_Poscar( Unit_Cell )
!========================================

 type(structure) , intent(out) :: unit_cell

 real*8            :: x0 , y0 , z0
 real*8            :: a , b , c
 real*8            :: kWH
 integer           :: n_kWH , n_residues , i , j , j1 , j2 , n , boundary_atom
 integer           :: N_of_atoms , N_of_elements , N_of_cluster_atoms , N_of_Configurations
 character(len=1)  :: TorF , fragment
 character(len=3)  :: residue
 logical           :: flag

 real*8           , allocatable :: xyz(:,:)
 integer          , allocatable :: atom_No(:)
 character(len=2) , allocatable :: element(:) , symbol(:)

 OPEN(unit=3,file='poscar.in',status='old')

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell
 read(3,*) n_residues

! allocating arrays ...
 CALL Allocate_UnitCell( Unit_Cell , n_residues )

! residues must be packed together ...
 do i = 1 , n_residues
    read(3,*) j1 , j2 , residue , fragment

    unit_cell % list_of_residues  (i) = residue
    unit_cell % list_of_fragments (i) = residue

    forall(j=j1:j2) unit_cell % residue  (j)  = residue
    forall(j=j1:j2) unit_cell % fragment (j) = fragment
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
        CALL Setting_Fragments( trj(j) )
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

! convert residues to upper case ...
forall( i=1:trj(1)%N_of_atoms ) trj(1)%atom(i)%residue = TO_UPPER_CASE( trj(1)%atom(i)%residue )

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
!============================
 subroutine Sym_2_AtNo_TRJ(a)
!============================
implicit none
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

end subroutine Sym_2_AtNo_TRJ
!
!
!============================
 subroutine Sym_2_AtNo_XYZ(a)
!============================
implicit none
type(structure) , intent(inout) :: a

! local variables ...
integer :: i

DO i = 1 , a%atoms

    select case(a%symbol(i))
        case( 'H') 
            a%AtNo(i) = 1 
        case( 'C') 
            a%AtNo(i) = 6 
        case( 'N') 
            a%AtNo(i) = 7 
        case( 'O') 
            a%AtNo(i) = 8 
        case( 'F') 
            a%AtNo(i) = 9 
        case( 'AL','Al') 
            a%AtNo(i) = 13 
        case( 'S','s') 
            a%AtNo(i) = 16 
        case( 'CL','Cl') 
            a%AtNo(i) = 17 
        case( 'TI','Ti') 
            a%AtNo(i) = 22 
        case( 'MN','Mn') 
            a%AtNo(i) = 25 
        case default
            print*, ' >> unkown atom found ; execution terminated <<' 
            stop
    end select

END DO

end subroutine Sym_2_AtNo_XYZ

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
subroutine Setting_Fragments( a )
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
!   Passivator  =   P 
!--------------------------------------------

 DO i = 1 , size(a%atom)

    select case(a%atom(i)%residue)
        case( 'CCC') 
            a%atom(i)%fragment = 'C' 
        case( 'ALQ') 
            a%atom(i)%fragment = 'M' 
        case( 'ACN') 
            a%atom(i)%fragment = 'S' 
        case( 'PYR') 
            a%atom(i)%fragment = 'P' 
    end select

 END DO

end subroutine Setting_Fragments
!
!
!
!==========================================
subroutine Identify_Fragments_Universe( a )
!==========================================
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

end subroutine Identify_Fragments_Universe
!
!
!
!=============================================
 subroutine Identify_Fragments_Structure ( a )
!=============================================
implicit none
type(structure)  , intent(inout) :: a

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag

allocate( temp(a % atoms) )

temp(1) = a % fragment(1)
counter = 1

do i = 1 , a % atoms

    flag = .true.
    do j = 1 , counter
     flag = flag .AND. ( temp(j) /= a%fragment(i) )
    end do

    if( flag ) then
        counter = counter + 1
        temp(counter) = a%fragment(i)
    end if

end do

! build list of fragments in a ...
allocate( a%list_of_fragments(counter) )
a%list_of_fragments = temp(1:counter)
deallocate( temp )

end subroutine Identify_Fragments_Structure
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

allocate( temp%atom(1) , temp%solvent(1) )

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

deallocate( temp%atom , temp%solvent )

end subroutine Sort_Residues
!
!
!
!======================================
 pure FUNCTION TO_UPPER_CASE ( STRING )
!======================================
 implicit none
 CHARACTER ( LEN = * )              , INTENT(IN)    :: STRING
 CHARACTER ( LEN = LEN ( STRING ) )                 :: TO_UPPER_CASE

! Local parameters ...
INTEGER, PARAMETER :: BIG_A = ICHAR ( "A" ), LITTLE_A = ICHAR ( "a" ), LITTLE_Z = ICHAR ( "z" )

! Local scalars ...
INTEGER :: I, ICHR

! Loop over the characters in the string ...
DO I = 1,LEN ( STRING )

!   Get the ASCII order for the character to be converted ...
    ICHR = ICHAR ( STRING(I:I) )

!   Use the order to change the case of the character ...
    IF ( ( ICHR >= LITTLE_A ) .AND. ( ICHR <= LITTLE_Z ) ) THEN
        TO_UPPER_CASE(I:I) = CHAR ( ICHR + BIG_A - LITTLE_A )
    ELSE
        TO_UPPER_CASE(I:I) = STRING(I:I)
    END IF
END DO

END FUNCTION TO_UPPER_CASE
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

