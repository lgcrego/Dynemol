 module Babel_m

    use IFPORT 
    use type_m                  
    use MM_input                , only : MM_input_format
    use parameters_m            , only : file_type,             &
                                         ad_hoc ,               &
                                         driver
    use Allocation_m            , only : Allocate_UnitCell
    use tuning_m                , only : ad_hoc_tuning
    use Babel_routines_m        , only : Symbol_2_AtNo ,        &
                                         Identify_Fragments ,   &
                                         AtNo_2_Symbol ,        &
                                         MMSymbol_2_Symbol ,    &
                                         Identify_Residues ,    &
                                         Pack_residues ,        &
                                         Sort_nr ,              &
                                         TO_UPPER_CASE ,        &
                                         Center_of_Gravity ,    &
                                         Initialize_System
    use Semi_Empirical_Parms    , only : atom                                         

!    private

    PUBLIC :: Read_from_XYZ , Read_from_Poscar , Read_from_PDB
    PUBLIC :: Read_PDB , Read_XYZ
    PUBLIC :: Coords_from_Universe
    PUBLIC :: Resume_from_TRJ 

    type(universe)      , allocatable   , public  :: trj(:)
    real*8                              , public  :: MD_dt
    integer             , allocatable   , public  :: QMMM_key(:)
    character(len=72)                   , public  :: System_Characteristics

contains
!
!
!
!=====================================================
 subroutine Coords_from_Universe( Unit_Cell , System )
!=====================================================
 implicit none
 type(structure)            , intent(out)   :: Unit_Cell
 type(universe)             , intent(inout) :: System

! local variables ... 
integer         :: j , n_residues

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

    unit_cell % coord              (j,:) =  System % atom(j) % xyz(:)
    unit_cell % AtNo               (j)   =  System % atom(j) % AtNo  
    unit_cell % fragment           (j)   =  System % atom(j) % fragment
    unit_cell % Symbol             (j)   =  System % atom(j) % Symbol
    unit_cell % residue            (j)   =  System % atom(j) % residue
    unit_cell % nr                 (j)   =  System % atom(j) % nr
    unit_cell % MMSymbol           (j)   =  System % atom(j) % MMSymbol
    unit_cell % solute             (j)   =  System % atom(j) % solute  
    unit_cell % DPF                (j)   =  System % atom(j) % DPF     
    unit_cell % El                 (j)   =  System % atom(j) % El      
    unit_cell % Hl                 (j)   =  System % atom(j) % Hl      
    unit_cell % flex               (j)   =  System % atom(j) % flex
    unit_cell % hardcore           (j)   =  System % atom(j) % hardcore
    unit_cell % solvation_hardcore (j)   =  System % atom(j) % solvation_hardcore
    unit_cell % V_shift            (j)   =  System % atom(j) % V_shift
    unit_cell % QMMM               (j)   =  System % atom(j) % QMMM

    unit_cell % Nvalen   (j)   =  atom(unit_cell%AtNo(j)) % Nvalen
    unit_cell % polar    (j)   =  atom(unit_cell%AtNo(j)) % polar

end forall

unit_cell%N_of_Solvent_Molecules = System%N_of_Solvent_Molecules

! get list of residues ...
select case( file_type )
    case( "structure" )
        If( size(unit_cell%list_of_fragments) /= size(system%list_of_fragments) ) then
            deallocate(unit_cell%list_of_fragments)
            allocate  (unit_cell%list_of_fragments(size(system%list_of_fragments)))
        end If
        unit_cell%list_of_residues  = System%list_of_residues
        unit_cell%list_of_fragments = System%list_of_fragments

    case( "trajectory" )
        If( size(unit_cell%list_of_fragments) /= size(trj(1)%list_of_fragments) ) then
            deallocate(unit_cell%list_of_fragments)
            allocate  (unit_cell%list_of_fragments(size(trj(1)%list_of_fragments)))
        end If
        unit_cell%list_of_residues  = trj(1)%list_of_residues
        unit_cell%list_of_fragments = trj(1)%list_of_fragments
end select

! sort the nr indices ...
if( MM_input_format == "GMX" ) CALL Sort_nr( unit_cell )

! nr indices must start with 1 ...
if( any(unit_cell % nr == 0) ) stop ">> halted: check input.pdb, nr indice = 0 found <<"

! unit_cell dimensions ...
unit_cell % T_xyz =  System % box

! standard Wolfgang-Helmholtz parameter ...
Unit_Cell % k_WH = 1.75d0

end subroutine Coords_from_Universe
!
!
!
!=====================================
 subroutine Read_from_XYZ( Unit_Cell )
!=====================================
 implicit none
 type(structure) , intent(out) :: unit_cell

 character(len=1)  :: dumb_char 
 integer           :: i , j 

 OPEN(unit=3,file='xyz.dat',status='old')   

! start reading the structure characteristics
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell
 read(3,'(A72)') System_Characteristics

! allocating arrays ...
 CALL Allocate_UnitCell( Unit_Cell )
 CALL Initialize_System( Unit_Cell )

! read unit_cell dimensions ...
 read(3,*) ( unit_cell%T_xyz(i) , i=1,3 )

! testing for AtNo or Symbol ...
 read(3,*) dumb_char
 backspace 3

! reading coordinates ...
 select case ( ichar(dumb_char) )

    ! Atomic Numbers ...    
    case( 48:57 )
        DO j = 1 , unit_cell%atoms
            read(3,*) unit_cell%AtNo(j), (unit_cell%coord(j,i), i=1,3)
        END DO   
        CALL AtNo_2_Symbol(Unit_Cell)

    ! Chemical Symbols ...
    case( 65:122 )
        DO j = 1 , unit_cell%atoms
            read(3,*) unit_cell%Symbol(j), (unit_cell%coord(j,i), i=1,3)
        END DO   
        CALL Symbol_2_AtNo(Unit_Cell)

 end select

 CLOSE(3)

 Unit_cell % Nvalen   (:)   =  atom(unit_cell%AtNo(:)) % Nvalen
 Unit_cell % polar    (:)   =  atom(unit_cell%AtNo(:)) % polar 
 Unit_cell % MMSymbol       =  Unit_Cell % Symbol

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
integer             :: i , j , N_of_atoms 
integer             :: k = 0
integer             :: file_err , io_err
character(len=5)    :: MMSymbol_char
character(len=6)    :: keyword
type(universe)      :: system        
logical             :: TorF

OPEN(unit=3 , file='input.pdb',status='old' , iostat=file_err , err=11)

!-------------------------------------------------------------------------------------------
read(unit=3 , fmt = 43 , iostat=io_err , err=12) System_Characteristics

! reading unit cell vectors ...
read(unit=3 , fmt=105 , iostat=io_err , err=12) keyword
if ( keyword == "CRYST1" ) then
    backspace 3
    read(3 , fmt=100 , iostat=io_err , err=12) system%box(1) , system%box(2) , system%box(3)
end if

! scan file for N_of_Atoms ...   
N_of_atoms = 0
do
    read(unit=3 , fmt=105 , iostat=io_err , err=12) keyword
    if ( keyword == "MASTER" ) exit
    N_of_atoms = N_of_atoms + 1
    if ( N_of_atoms > 20000 ) then
       TorF = systemQQ("sed '11i *** reading of input.pdb halted forcibly ; ckeck input.pdb file format ***  ' warning.signal |cat")
       STOP 
    end if
end do
system%N_of_atoms = N_of_atoms

allocate( system%atom(system%N_of_atoms) )
CALL Initialize_System( system )

!read data ...
rewind 3
do
    read(unit=3 , fmt=105 , iostat=io_err , err=12) keyword
    if( keyword == "CRYST1" ) then
        do i = 1 , system%N_of_atoms
                         system%atom(i)%my_id = i
            read(3 , 115 , iostat=io_err , err=12)             &
                         MMSymbol_char                      ,  &    ! <== atom type
                         system%atom(i)%residue             ,  &    ! <== residue name
                         system%atom(i)%nr                  ,  &    ! <== residue sequence number
                        (system%atom(i)%xyz(j) , j=1,3)     ,  &    ! <== xyz coordinates 
                         system%atom(i)%Symbol                      ! <== chemical element symbol

                         system%atom(i)%MMSymbol = adjustl(MMSymbol_char)
        end do
    end if
    k = k + 1
    if ( keyword == "MASTER" ) exit
    if ( k > 20000 ) stop " *** reading of input.pdb halted ; ckeck file format *** " 
end do
!-------------------------------------------------------------------------------------------

close(3)

! convert residues to upper case ...
forall( i=1:system%N_of_atoms ) system%atom(i)%residue = TO_UPPER_CASE( system%atom(i)%residue )

! use ad hoc tuning of EHT parameters ...
If( ad_hoc .and. (driver /= "MM_Dynamics") ) CALL ad_hoc_tuning( system )

! for use if atomic Symbols are not included in input.pdb ...
!If( sum( len_trim(system%atom%Symbol) ) == 0 ) CALL MMSymbol_2_Symbol( system%atom )
CALL MMSymbol_2_Symbol( system%atom )

! preprocessing the universe system ...
CALL Symbol_2_AtNo      ( system%atom )
CALL Identify_Residues  ( system      )
CALL Identify_Fragments ( system      )

!CALL Pack_Residues      ( system%atom , system%list_of_residues )

! vector QMMM_key is the key to exchange QM and MM atomic properties ...
allocate( QMMM_key(system%N_of_atoms) , source=system%atom%my_id)

If( (verify(driver,'slice') == 6) .AND. (file_type == "trajectory") ) then

    ! replicate system throughout trj ...
    CALL input_2_frames( system , trj )

else

    ! transfer Unit_Cell(structure) <-- system(universe) 
    CALL Coords_from_Universe( Unit_Cell , system )

end if

deallocate( system%atom , system%list_of_fragments , system%list_of_residues )
11 if( file_err > 0 ) stop "input.pdb file not found; terminating execution"
12 if( io_err   > 0 ) stop "problems reading input.pdb; check IO_file_formats"

43  format(a72)
100 format(t10, f6.3, t19, f6.3, t28, f6.3)
105 format(a6)
115 FORMAT(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3,t77,a2)

end subroutine Read_from_PDB
!
!
!
!========================================
 subroutine Read_from_Poscar( Unit_Cell )
!========================================
 implicit none
 type(structure) , intent(out) :: unit_cell

! local variables ... 
 real*8            :: x0 , y0 , z0
 real*8            :: a , b , c
 real*8            :: k_WH
 integer           :: n_kWH , n_residues , i , j , j1 , j2 , n , boundary_atom , openstatus
 integer           :: N_of_atoms , N_of_elements , N_of_cluster_atoms , N_of_Configurations
 character(len=1)  :: TorF , fragment
 character(len=3)  :: residue
 logical           :: flag , TF_flag

 real*8           , allocatable :: xyz(:,:)
 integer          , allocatable :: atom_No(:)
 character(len=2) , allocatable :: element(:) , symbol(:)

 open(unit = 3, file = 'poscar.in', status = 'old', action = 'read', iostat = openstatus)
 if (openstatus > 0) then
    TF_flag = systemQQ("sed '11i *** Cannot open the file poscar.in *** ' warning.signal |cat")
    STOP 
 end if

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
implicit none
type(universe)  , allocatable   , intent(out)   :: trj(:)

! local variables ...
integer      :: i , j , k  , dumb_number , openstatus , inputstatus , model , number_of_atoms
real*8       :: time_1 , time_2 , delta_t 
character(1) :: test
character(4) :: keyword
character(5) :: MMSymbol_char
logical      :: TorF

open(unit = 31, file = 'frames.pdb', status = 'old', action = 'read', iostat = openstatus)
if (openstatus > 0) then
    TorF = systemQQ("sed '11i *** Cannot open the file frames.pdb *** ' warning.signal |cat")
    STOP 
end if

read(unit = 31, fmt = 43, iostat = inputstatus) System_Characteristics

! find the number of model frames ...
model = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( inputstatus /= 0  ) exit
    if ( keyword == 'MODE' ) model = model + 1
end do

! return to the top of the file ...
rewind 31

! read number the atoms and time ...
read(unit = 31, fmt = 35, iostat = inputstatus) keyword
do
    if ( keyword == 'MODE' ) then
        exit
    else
        ! reading the time ...
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

number_of_atoms = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'ATOM' ) then
        backspace 31
        read(unit = 31, fmt = 32, iostat = inputstatus) dumb_number
        number_of_atoms = number_of_atoms + 1
    else
        ! reading the time ...
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

! Molecular Dynamics time step (pico-sec) ...
MD_dt = delta_t + epsilon(1.d0)

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
        CALL Initialize_System( trj(j) )

        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 33, iostat = inputstatus)   MMSymbol_char               ,     &
                                                              trj(j) % atom(i) % residue ,      &
                                                              trj(j) % atom(i) % nr ,           &
                                                            ( trj(j) % atom(i) % xyz(k) , k=1,3 )

            trj(j) % atom(i) % MMSymbol = adjustl(MMSymbol_char)
        end do

        ! convert residues to upper case ...
        forall( i=1:number_of_atoms ) trj(j)%atom(i)%residue = TO_UPPER_CASE( trj(j)%atom(i)%residue )

        CALL MMSymbol_2_Symbol  ( trj(j)%atom )
        CALL Symbol_2_AtNo      ( trj(j)%atom )
        ! use ad hoc tuning of EHT parameters ...
        If( ad_hoc .and. (driver /= "MM_Dynamics") ) CALL ad_hoc_tuning( trj(j) )

        trj(j)%atom % Nvalen    =  atom(trj(j)%atom%AtNo) % Nvalen
        trj(j)%atom % polar     =  atom(trj(j)%atom%AtNo) % polar 

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
        CALL Initialize_System( trj(j) )
    
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

    trj(i) % atom % AtNo                =  trj(1) % atom % AtNo    
    trj(i) % atom % MMSymbol            =  trj(1) % atom % MMSymbol
    trj(i) % atom % residue             =  trj(1) % atom % residue
    trj(i) % atom % nr                  =  trj(1) % atom % nr    
    trj(i) % atom % Symbol              =  trj(1) % atom % Symbol
    trj(i) % atom % fragment            =  trj(1) % atom % fragment
    trj(i) % atom % solute              =  trj(1) % atom % solute
    trj(i) % atom % DPF                 =  trj(1) % atom % DPF   
    trj(i) % atom % El                  =  trj(1) % atom % El    
    trj(i) % atom % Hl                  =  trj(1) % atom % Hl    
    trj(i) % atom % hardcore            =  trj(1) % atom % hardcore   
    trj(i) % atom % solvation_hardcore  =  trj(1) % atom % solvation_hardcore   
    trj(i) % atom % V_shift             =  trj(1) % atom % V_shift
    trj(i) % atom % QMMM                =  trj(1) % atom % QMMM   

    trj(i) % atom % Nvalen              =  atom(trj(1)%atom%AtNo) % Nvalen
    trj(i) % atom % polar               =  atom(trj(1)%atom%AtNo) % polar 

end forall

! GROUP residues ...
do i = 1 , size(trj)
    CALL Pack_Residues( trj(i)%atom , trj(1)%list_of_residues ) 
end do

! Information about the Solvent System ...

! solvent molecules must be contiguous ...
trj%N_of_Solvent_Molecules =  maxval( trj(1) % atom % nr, trj(1) % atom % fragment == 'S' )         &
                            - minval( trj(1) % atom % nr, trj(1) % atom % fragment == 'S' ) + 1

do i = 1 , size(trj)

    allocate( trj(i) % solvent( trj(i)%N_of_Solvent_Molecules ) )

    trj(i)%solvent%N_of_Atoms = count( trj(i)%atom%fragment == 'S' ) / trj(i)%N_of_Solvent_molecules

    CALL Center_of_Gravity( trj(i) )

end do

! Formats ...
32 format(5x, i6)
33 format(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)
35 format(a4)
36 format(7x, i7)
37 format(t31,f8.3,t39,f8.3,t47,f8.3)
38 format(10x, a1)
39 format(81x, f7.0)
40 format(6x, 3f9.3)
41 format(f13.7)
42 format(a1)
43 format(a72)

end subroutine Read_PDB
!
!
!
!========================
 subroutine Read_XYZ(trj)
!========================
implicit none
type(universe)  , allocatable   , intent(out) :: trj(:)

! local variables ....
character(1)  :: idx
integer       :: openstatus , inputstatus , atoms , i , j , k , model 
integer       :: j1 , j2 , n_residues 
character(1)  :: fragment
character(3)  :: residue 
logical       :: TorF

open(unit = 13, file = 'XYZ.trj', status = 'old', action = 'read', iostat = openstatus)
if( openstatus > 0 ) then
    TorF = systemQQ("sed '11i *** Cannot open the file XYZ.trj *** ' warning.signal |cat")
    STOP 
end if

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

! read the number the atoms ...
rewind 13
do
    read(unit = 13, fmt = 21, iostat = inputstatus) idx
    if( idx == '$' ) exit
end do
atoms = 0
do
    read(unit = 13, fmt = 21, iostat = inputstatus) idx
    if( idx /= '$' ) then
        atoms = atoms + 1
    else
        exit
    end if
end do

! prepare to read the system ...
allocate( trj(model) )
allocate( trj(1)%atom(atoms) )
CALL Initialize_System( trj(1) )

! start reading the structure characteristics ...
 rewind 13
 read(13,*) System_Characteristics
 trj(1) % System_Characteristics = System_Characteristics

! residues must be packed together ...
 read(13,*) n_residues
 do i = 1 , n_residues
    read(13,*) j1 , j2 , residue , fragment
    forall(j=j1:j2) trj(1) % atom(j) % residue  = residue
    forall(j=j1:j2) trj(1) % atom(j) % fragment = fragment
 end do

! read PBC vectors ... 
read(13,*) trj(1)%box(1) , trj(1)%box(2) , trj(1)%box(3)

! read coordinate from the frames ...
do j = 1 , model

    if( j == 1 ) then
        read(unit = 13, fmt = 21, iostat = inputstatus) idx

        do i = 1 , atoms
            read(unit = 13, fmt = 22, iostat = inputstatus) trj(j)%atom(i)%Atno, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
        CALL AtNo_2_Symbol(trj(j)%atom)
    else
        read(unit = 13, fmt = 21, iostat = inputstatus) idx

        allocate( trj(j)%atom(atoms) )
        CALL Initialize_System( trj(j) )

        do i = 1 , atoms
            read(unit = 13, fmt = 23, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if

    ! use ad hoc tuning of EHT parameters ...
    If( ad_hoc .and. (driver /= "MM_Dynamics") ) CALL ad_hoc_tuning( trj(j) )

end do

trj(1)%N_of_atoms = atoms

! convert residues to upper case ...
forall( i=1:trj(1)%N_of_atoms ) trj(1)%atom(i)%residue = TO_UPPER_CASE( trj(1)%atom(i)%residue )

! get list of residues in trj ...
CALL Identify_Residues( trj(1) )

! get list of fragments in trj ...
CALL Identify_Fragments( trj(1) )

! Copy information from trj(1) to trj(:) ...
forall(i = 2:model )
    trj(i) % N_of_atoms      = trj(1) % N_of_atoms
    trj(i) % box             = trj(1) % box
    trj(i) % atom % Atno     = trj(1) % atom % Atno
    trj(i) % atom % Symbol   = trj(1) % atom % Symbol
    trj(i) % atom % fragment = trj(1) % atom % fragment
    trj(i) % atom % residue  = trj(1) % atom % residue
    trj(i) % atom % hardcore = trj(1) % atom % hardcore   
    trj(i) % atom % V_shift  = trj(1) % atom % V_shift
    trj(i) % atom % QMMM     = trj(1) % atom % QMMM   
end forall

! GROUP residues ...
do i = 1 , size(trj)
    CALL Pack_Residues( trj(i)%atom , trj(1)%list_of_residues )
end do

! formats ...
21 format(a1)
22 format(i2, t4, f8.5, t13, f8.5, t22, f8.5)
23 format(3x, f8.5, t13, f8.5, t22, f8.5)

end subroutine Read_XYZ
!
!
!
!========================================
 subroutine Resume_from_TRJ ( Unit_Cell )
!========================================
implicit none
type(structure) , intent(out) :: Unit_Cell

! local variables ...
integer        :: i , j , k  , openstatus , inputstatus , model , number_of_atoms
real*8         :: dumb_xyz(3) , dumb_box(3)
character(4)   :: keyword
character(5)   :: MMSymbol_char
type(universe) :: system
logical        :: TorF

open(unit = 31, file = 'frames.pdb', status = 'old', action = 'read', iostat = openstatus)
if (openstatus > 0) then
    TorF = systemQQ("sed '11i *** Cannot open the file frames.pdb *** ' warning.signal |cat")
    STOP 
end if

read(unit = 31, fmt = 43, iostat = inputstatus) system% System_Characteristics 

!----------------------------------------------------------------------------
! find the number of model frames ...
model = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( inputstatus /= 0  ) exit
    if ( keyword == 'MODE' ) model = model + 1
end do

! return to the top of the file ...
rewind 31

!----------------------------------------------------------------------------
! read number of atoms and time ...
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'MODE' ) exit
end do

number_of_atoms = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'ATOM' ) then
        number_of_atoms = number_of_atoms + 1
    else
        exit   
    end if
end do

allocate( system% atom(number_of_atoms) )
CALL Initialize_System( system )

system% N_of_atoms = number_of_atoms

!----------------------------------------------------------------------------
! return to the top of the file ...
rewind 31

do j = 1 , model

    if( j < model ) then

        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( dumb_box(i) , i=1,3 )
            end if
            if( keyword == 'MODE' ) exit
        end do

        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 37, iostat = inputstatus) ( dumb_xyz(k) , k=1,3 )
        end do

    elseif( j == model ) then ! <== last frame ...

        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( system% box(i) , i=1,3 )
            end if
            if ( keyword == 'MODE' ) exit
        end do

        do i = 1 , number_of_atoms
            system% atom(i)% my_id = i
            read(unit = 31, fmt = 33, iostat = inputstatus)   MMSymbol_char               ,     &
                                                              system % atom(i) % residue ,      &
                                                              system % atom(i) % nr ,           &
                                                            ( system % atom(i) % xyz(k) , k=1,3 )

            system % atom(i) % MMSymbol = adjustl(MMSymbol_char)
        end do

        ! convert residues to upper case ...
        forall( i=1:number_of_atoms ) system %atom(i)% residue = TO_UPPER_CASE( system% atom(i)% residue )

        CALL MMSymbol_2_Symbol  ( system% atom )
        CALL Symbol_2_AtNo      ( system% atom )
        ! use ad hoc tuning of EHT parameters ...
        If( ad_hoc .and. (driver /= "MM_Dynamics") ) CALL ad_hoc_tuning( system )

        system%atom % Nvalen    =  atom(system%atom%AtNo) % Nvalen
        system%atom % polar     =  atom(system%atom%AtNo) % polar 
    
    end if

end do

close(31)

! preprocessing the universe system ...
CALL Symbol_2_AtNo      ( system%atom )
CALL Identify_Residues  ( system )
CALL Identify_Fragments ( system )

! GROUP residues ...
!CALL Pack_Residues( system%atom , system%list_of_residues ) 

! vector QMMM_key is the key to exchange QM and MM atomic properties ...
allocate( QMMM_key(system%N_of_atoms) , source=system%atom%my_id)

CALL Coords_from_Universe( Unit_Cell , system )

! save copy of the just-generated initial configuration ...
CALL Dump_pdb( system )

! Formats ...
33 format(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)
35 format(a4)
37 format(t31,f8.3,t39,f8.3,t47,f8.3)
40 format(6x, 3f9.3)
43 format(a72)

end subroutine Resume_from_TRJ
!
!
!
!
!=========================================
 subroutine input_2_frames( system , trj )
!=========================================
use parameters_m  , only: t_i , t_f , n_t  
implicit none
type(universe)                    , intent(in)    :: system
type(universe)    , allocatable   , intent(out)   :: trj(:)

! local variables ...
integer :: i  
real*8  :: delta_t

delta_t = ( t_f - t_i ) / float(n_t)

! Molecular Dynamics time step (pico-sec) ...
MD_dt = delta_t 

! allocate array of frames ...
allocate( trj(n_t) )

do i = 1 , n_t

    allocate( trj(i)%atom(system%N_of_atoms) )
    CALL Initialize_System( trj(i) )

end do

! Copy information from input.pdb to trj(:) ...
forall( i = 1:n_t ) trj(i) = system

end subroutine input_2_frames
!
!
!
!=======================
subroutine Dump_pdb(sys)
!=======================
implicit none 
type(universe)      , intent(inout) ::  sys

! local variables ...
integer ::  i , k 

OPEN(unit=4,file='input.pdb',status='unknown')

write(4,6) "Generated by Resume_from_TRJ, ",sys%System_Characteristics
write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms

            write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        sys%atom(i)%MMSymbol            ,  &    ! <== atom type
                        ' '                             ,  &    ! <== alternate location indicator
                        sys%atom(i)%residue             ,  &    ! <== residue name
                        ' '                             ,  &    ! <== chain identifier
                        sys%atom(i)%nr                  ,  &    ! <== residue sequence number
                        ' '                             ,  &    ! <== code for insertion of residues
                        ( sys%atom(i)%xyz(k) , k=1,3 )  ,  &    ! <== xyz coordinates 
                        1.00                            ,  &    ! <== occupancy
                        0.00                            ,  &    ! <== temperature factor
                        ' '                             ,  &    ! <== segment identifier
                        ' '                             ,  &    ! <== here only for tabulation purposes
                        sys%atom(i)%Symbol              ,  &    ! <== chemical element symbol
                        ' '                                     ! <== charge on the atom
end do

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , sys%N_of_atoms , 0 , sys%N_of_atoms , 0
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,a2)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a32,a72)

end subroutine Dump_pdb
!
!
!
end module Babel_m

