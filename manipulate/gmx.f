module GMX_routines

use types_m
use Read_Parms          , only : MMSymbol_2_Symbol , Symbol_2_AtNo , Atomic_mass
use RW_routines         , only : Initialize_System
use diagnosis_m
use FUNCTION_routines   , only : res , Solvent_residue_groups 
use Topology_routines   , only : connect , dump_topol


! module variables ...
integer, allocatable, private :: InputIntegers(:,:) 
logical, allocatable, private :: bond_matrix(:,:) , angle_matrix(:,:,:) , dihedral_matrix(:,:,:,:)
logical             , private :: done = .false.

contains
!
!
!
!============================
 subroutine save_GROMACS(sys)
!============================
implicit none 
type(universe) , intent(inout) :: sys

! local variables ...
character(1)                        :: answer
logical             , parameter     :: BACK = .TRUE. 

!determine charge groups for TiO2 itp file ...
!CALL TiO2_charge_groups(sys)

!determine residue groups for solvent ...
if( (count(sys%atom%fragment == "S") /= 0) .AND. (sum(sys%atom%nresid) == 0) ) CALL Solvent_residue_groups(sys)

!determine residue groups for fragment ...
if( count(sys%atom%fragment == "F") /= 0 ) then
end if

! sorting by fragment ...
CALL Sort_Fragments(sys)

! where MMSymbol is not defined MMSymbol = symbol ...
where( sys % atom % MMSymbol == "XXX" ) sys % atom % MMSymbol = sys % atom % Symbol

CALL diagnosis(sys)

CALL Connect(sys)

CALL Dump_pdb(sys)

write(*,'(/a)') ">>>  Save itp file ? (y/n)"
read (*,'(a)') answer

If( answer == "y" ) CALL Dump_itp(sys)

end subroutine save_GROMACS
!
!
!
!=================================
subroutine Dump_pdb( sys , title )
!=================================
implicit none 
type(universe)                  , intent(inout) ::  sys
character(*)        , optional  , intent(in)    :: title

! local variables ...
integer ::  i , k

!----------------------------------------------
!     generate pdb file for GROMACS
!----------------------------------------------

If( present(title) ) then
    write(4,'(A96)') title
else
    OPEN(unit=4,file='seed.pdb',status='unknown')
    write(4,6) sys%Surface_Characteristics
end if

write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms

            write(4,2)  'ATOM  '                        ,  &    ! <== non-standard atom
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

! check and print topological connections ...
If( allocated(sys%topol) ) CALL dump_topol(sys,4)

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , sys%N_of_atoms , 0 , sys%N_of_atoms , 0
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a6,a72)

end subroutine Dump_pdb
!
!
!
!=========================
subroutine Dump_itp( sys )
!=========================
implicit none 
type(universe), intent(inout) :: sys

! local variables ...
integer                       :: i, j, n, total_bonds, total_angs, total_diheds 
integer                       :: mol_conect

!----------------------------------------------
!         generate topology files 
!----------------------------------------------

OPEN(unit=10,file='seed.itp',status='unknown')

do n = 1 , maxval(sys%atom%nresid)
   
       !----------------------------------------
       ! heading 
       !----------------------------------------
       write(10,101) "[ moleculetype ]"
       write(10,*) , "SPECIES", 3
       write(10,*) 
       
       write(10,102) "[ atoms ]"
       do i = 1 , sys%N_of_atoms

           write(10,5) i                    ,  &  ! <== serial number within the residue
                       sys%atom(i)%Symbol   ,  &  ! <== force field descriptor of atom type
                       sys%atom(i)%nresid   ,  &  ! <== residue identifier
                       sys%atom(i)%resid    ,  &  ! <== residue name
                       sys%atom(i)%MMSymbol ,  &  ! <== atom type
                       sys%atom(i)%nrcg     ,  &  ! <== charge group
                       sys%atom(i)%charge   ,  &  ! <== charge of atom type       
                       sys%atom(i)%mass           ! <== mass of chemical element 
       end do
end do
write(10,*)
write(10,*)

!----------------------------------------
! start topology connections
!----------------------------------------
mol_conect = 0
do i = 1 , sys%total_conect
    j = size( pack( InputIntegers(i,:) , InputIntegers(i,:) /= 0 ) )
    mol_conect = j + mol_conect
end do

if( mol_conect == 0 ) then
    Write(*,*)
    Write(*,*) "======================= W A R N I N G ======================="
    Write(*,*) "No CONECT in the pdb file; cannot generate bonds, angles, etc"
    Write(*,*) "============================================================="
    Write(*,*)
endif

!----------------------------------------
! Assign CONECT to a logical bond matrix ...
!----------------------------------------
CALL get_bond_matrix( sys , total_bonds ) 

!--------------------------------------------------------
! Assign bond matrix to a logical i,k,j angle matrix ...
!--------------------------------------------------------
if( mol_conect > 3 ) then 
    CALL get_angle_matrix( sys , total_angs ) 
end if

!-----------------------------------------------------------
! Assign angle matrix to a logical i--l dihedral matrix ...
!-----------------------------------------------------------
if( mol_conect > 4 ) then 
    CALL get_dihedral_matrix( sys , total_diheds ) 
end if

CALL write_seed_itp( sys%atom )

deallocate( bond_matrix, angle_matrix, dihedral_matrix )

5   FORMAT(i6,a8,i8,a7,a8,i8,F9.4,F9.4)
101 FORMAT(a16)
102 FORMAT(a9)

end subroutine Dump_itp
!
!
!
!=================================
subroutine TiO2_Charge_groups(sys)
!=================================
implicit none
type(universe)  , intent(inout) :: sys

! local variables ...
integer                         :: i , j , indx , Ti_size , O2_size , rest_size
real*8                          :: distance , total_charge_of_group
logical                         :: flag(4)
type(atomic)    , allocatable   :: pbc(:) , temp_Ti(:) , temp_O2(:) , temp_rest(:)

! sort TiO2 atoms ... 
Ti_size = count( sys%atom%Symbol == 'Ti' .AND. sys%atom%resid == 'CCC' )
O2_size = count( sys%atom%Symbol == 'O'  .AND. sys%atom%resid == 'CCC' )

rest_size = sys%N_of_atoms - Ti_size - O2_size

allocate( temp_Ti(Ti_size) , temp_O2(O2_size) , temp_rest(rest_size) )

! separate TiO2 from the rest ...
temp_Ti   =  pack( sys%atom , (sys%atom%Symbol == 'Ti') .AND. (sys%atom%resid == 'CCC') , temp_Ti   )
temp_O2   =  pack( sys%atom , (sys%atom%Symbol == 'O' ) .AND. (sys%atom%resid == 'CCC') , temp_O2   )
temp_rest =  pack( sys%atom ,                                  sys%atom%resid /= 'CCC'  , temp_rest )

sys%atom( 1                 : Ti_size         )  =  temp_ti
sys%atom( Ti_size+1         : Ti_size+O2_size )  =  temp_O2
sys%atom( Ti_size+O2_size+1 : sys%N_of_atoms  )  =  temp_rest

deallocate( temp_Ti , temp_O2 , temp_rest )

! PBC ...
if( product(sys%box) == 1.0 ) pause ">>> check dimensions of the cell <<<"
CALL Generate_PBC(sys,pbc)

! establish charge groups for the cluster (nrcg) ...
where( sys%atom%resid == 'CCC' ) sys%atom%nrcg = 0

do i = 1 , sys%N_of_atoms 

    if( sys%atom(i)%symbol == 'Ti' ) then

        sys%atom(i)%nrcg = i

        do j = 1 , size(pbc)  

            if(pbc(j)%symbol == 'O') then

                indx = j - int( (j-1)/sys%N_of_atoms ) * sys%N_of_atoms 

                distance = dsqrt( sum( (sys%atom(i)%xyz - pbc(j)%xyz)**2 ) )

                flag(1) = ( dabs( pbc(j)%xyz(2)-sys%atom(i)%xyz(2) ) <= 0.5 )
                flag(2) = ( dabs( pbc(j)%xyz(3)-sys%atom(i)%xyz(3) ) <= 1.3 )
                flag(3) = ( distance < 2.2 )
                flag(4) = ( sys%atom(indx)%nrcg == 0 )

                if( flag(1) .AND. flag(2) .AND. flag(3) .AND. flag(4) ) sys%atom(indx)%nrcg = i

            end if

        end do
    end if
end do

deallocate( pbc )

! checking nrcg distribution for TiO2 ...
do i = 1 , sys%N_of_atoms

    if( sys%atom(i)%resid == 'CCC' ) then

        total_charge_of_group = sum( sys%atom%charge , sys%atom%nrcg == i )

        if( total_charge_of_group /= 0.0 ) Print*,'>>> Charge Group ',i,'is not neutral !!', total_charge_of_group
        if( sys%atom(i)%nrcg      == 0   ) Print*,'>>> Atom ',i,'(',sys%atom(i)%symbol,') did not find a group !!'

    end if

end do

end subroutine TiO2_Charge_Groups
!
!
!
!===============================
subroutine Generate_PBC(sys,pbc)
!===============================
implicit none
type(universe)                  , intent(inout) :: sys
type(atomic)    , allocatable   , intent(out)   :: pbc(:)

! local variables ...
integer :: PBC_N_of_atoms , i , j , n , counter , N_of_cells

N_of_cells = 9
PBC_N_of_atoms = sys%N_of_atoms * N_of_cells
allocate( pbc( PBC_N_of_atoms ) )

pbc(1:sys%N_of_atoms) = sys%atom

! replicating around the original cell ...
counter = sys%N_of_atoms

do j = -1 , +1 
do i = -1 , +1 

    If( (i /= 0) .OR. (j /= 0) ) THEN

        do n = 1 , sys%N_of_atoms

            counter = counter + 1

            pbc(counter) % xyz(1)   = sys % atom(n) % xyz(1) + i * sys%box(1)
            pbc(counter) % xyz(2)   = sys % atom(n) % xyz(2) + j * sys%box(2)
            pbc(counter) % xyz(3)   = sys % atom(n) % xyz(3) 
            pbc(counter) % symbol   = sys % atom(n) % symbol

        end do
    end if

end do
end do    

end subroutine Generate_PBC
!
!
!
!=================================
subroutine Sort_Fragments(system)
!=================================
implicit none
type(universe) , intent(inout) :: system

!	local variables
integer        :: i , j , N_of_atoms ,  iptr
type(universe) :: temp

allocate( temp%atom(1) )

N_of_atoms = size(system%atom)

do i = 1 , N_of_atoms-1

    iptr = i

    do j = i+1 , N_of_atoms
        if( LLT( system%atom(j)%fragment , system%atom(iptr)%fragment ) ) then
            iptr = j
        end if
    end do

    if( i /= iptr ) then
        temp%atom(1)      = system%atom(i)
        system%atom(i)    = system%atom(iptr)
        system%atom(iptr) = temp%atom(1)
    end if

end do

deallocate( temp%atom )

end subroutine Sort_Fragments
!
!
!
!========================================================
subroutine read_GROMACS( system , file_name , file_type )
!========================================================
implicit none 
type(universe)              , intent(out) :: system
character(*)    , optional  , intent(in)  :: file_name
character(*)    , optional  , intent(in)  :: file_type

! local variables ...
integer                         :: i, j, indx, ioerr, useless, N_of_atoms
character(len=80)               :: line
character(len=5)                :: MMSymbol_char , FFSymbol_char
character(len=6)                :: keyword
character(len=3)    , parameter :: ACN(6) = ['YN','YC','CT','HC','HC','HC']

! finds out what file to read ...
If( present(file_name) ) then

    OPEN(unit=3,file=file_name,status='old',iostat=ioerr,err=10)

else 
    if( present(file_type) ) then

        select case(file_type)
            case('pdb')
            OPEN(unit=3,file='input.pdb',status='old',iostat=ioerr,err=11)

            case('gro')
            OPEN(unit=3,file='input.gro',status='old',iostat=ioerr,err=12)
        end select

        end if
end if

! whether .gro or .pdb files ...
If( (verify("gro",file_name)==0) .OR. (verify("gro",file_type)==0) ) then

!---------------
! .gro files 
!---------------

    read(3,*) system%Surface_Characteristics
    read(3,*) system%N_of_atoms

    allocate( system%atom(system%N_of_atoms) )
    CALL Initialize_System( system )

!   reads the data ...
    do i = 1 , system%N_of_atoms 

        read(3,30,iostat=ioerr) system%atom(i)%nresid   ,   &
                                system%atom(i)%resid    ,   &
                                system%atom(i)%MMSymbol ,   &
                                useless ,                   &
                                (system%atom(i)%xyz(j),j=1,3) 

    end do

!   reads the unit cell vectors for Direct coordinate mode
    read(3,*) system%box(1) , system%box(2) , system%box(3)

!   nm --> Angs  ...
    system%box = system%box * 10.d0
    do i = 1 , system%N_of_atoms  
        system%atom(i)%xyz(:) = system%atom(i)%xyz(:) * 10.d0
    end do

else
!---------------
! .pdb files 
!---------------

    read(3,99) system%Surface_Characteristics
    
!   reads the unit cell vectors for Direct coordinate mode ...
    read(unit=3,fmt=105,iostat=ioerr) keyword

    if ( keyword == "CRYST1" ) then
        backspace 3
        read(3,fmt=100) system%box(1) , system%box(2) , system%box(3)
    end if

!   scan file for N_of_Atoms ...   
    N_of_atoms = 0
    do
        read(unit=3,fmt=105,iostat=ioerr) keyword
        if ( keyword == "MASTER" .or. keyword == "CONECT" ) exit
        N_of_atoms = N_of_atoms + 1
        print*, N_of_atoms
    end do
    system%N_of_atoms = N_of_atoms
        
    allocate( system%atom(system%N_of_atoms) )
    CALL Initialize_System( system ) 

    rewind 3

!   read data ...    
    do
        read(unit=3,fmt=105,iostat = ioerr) keyword
        if( keyword == "CRYST1" .or. keyword == "AUTHOR" ) then
            do i = 1 , system%N_of_atoms
                read(3,115)  MMSymbol_char                      ,  &    ! <== atom type
                             system%atom(i)%resid               ,  &    ! <== residue name
                             system%atom(i)%nresid              ,  &    ! <== residue sequence number
                             (system%atom(i)%xyz(j) , j=1,3)    ,  &    ! <== xyz coordinates 
                             system%atom(i)%symbol              ,  &    ! <== chemical element symbol
                             system%atom(i)%charge              ,  &    ! <== atom MM charge 
                             FFSymbol_char                              ! <== FF atom type

                system%atom(i)%MMSymbol = adjustl(MMSymbol_char)
                system%atom(i)%FFSymbol = adjustl(FFSymbol_char)
                system%atom(i)%my_intra_id = i 
                print*, system%atom(i)%MMSymbol 
            end do
        end if
        if ( keyword == "MASTER" .or. keyword == "CONECT" .or. keyword == "END" ) exit
    end do
    backspace(3)

   ! Generating a bond matrix for topology generation ...
   if ( keyword == "CONECT" ) then 
        allocate( InputIntegers(2*N_of_atoms,5), source = 0 )
        i = 0
        do
          read(3,103) line
          if( trim(line(1:6)) == "MASTER" ) exit
          i = i + 1
          read(line(7:11) ,'(I5)') InputIntegers(i,1)
          read(line(12:16),'(I5)') InputIntegers(i,2)
          read(line(17:21),'(I5)') InputIntegers(i,3)
          read(line(22:26),'(I5)') InputIntegers(i,4)
          read(line(27:31),'(I5)') InputIntegers(i,5)
        end do
        system%total_conect = i
   end if 

    system% atom% resid = adjustl(system% atom% resid)
!----------------------
! finished reading ...
!----------------------
end if

close(3)

! get Chemical symbol ...
CALL MMSymbol_2_Symbol( system%atom )

! get Atomic Number (AtNo) ...
CALL Symbol_2_AtNo( system%atom )

! fix the MMSymbols ...
indx=0
do 
    indx=indx+1
    select case( system%atom(indx)%resid )

        case('ACN')
            do j = 1 , size(ACN)
                system%atom(indx+j-1)%MMSymbol = ACN(j)
            end do
            indx = indx + size(ACN)-1

    end select
    if( indx >= system%N_of_atoms ) EXIT 
end do

10 if( ioerr > 0 ) stop "file_name file not found; terminating execution"
11 if( ioerr > 0 ) stop "input.pdb file not found; terminating execution"
12 if( ioerr > 0 ) stop "input.gro file not found; terminating execution"

30  format(I5,A3,A7,I5,3F8.4)
99  format(a72)
100 format(t10, f6.3, t19, f6.3, t28, f6.3)
103 format(a80)
105 format(a6)
110 format(t8, i4)
115 FORMAT(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3,t77,a2,t80,f8.4,t90,a3)


end subroutine read_GROMACS
!
!
!
!========================
subroutine gro_2_pdb(sys)
!========================
implicit none 
type(universe) , intent(inout) ::  sys

! local variables ...
integer ::  i , k

!----------------------------------------------
!     generate pdb file from gro file
!----------------------------------------------

OPEN(unit=4,file="seed.pdb",status="unknown")

write(4,6) 'COMPND' , '"',sys%Surface_Characteristics,'"'
write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms
    write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
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
                ' '                                     ! <== charge on the atom
end do

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , sys%N_of_atoms , 0 , sys%N_of_atoms , 0
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,a2)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a6,3x,a,a72,a)

end subroutine gro_2_pdb
!
!
!
!=============================================
subroutine Identify_Residues( sys , residues )
!=============================================
implicit none 
type(universe)                  , intent(inout) :: sys
character(3)    , allocatable   , intent(out)   :: residues(:)

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag


allocate( temp(sys%N_of_atoms) )

temp(1) = sys % atom(1) % resid

counter = 1

do i = 1 , sys%N_of_atoms

    flag = .true.
    do j = 1 , counter
        flag = flag .AND. (temp(j) /= sys%atom(i)%resid)
    end do

    if( flag ) then

        counter = counter + 1
        temp(counter) = sys%atom(i)%resid

    end if

end do

allocate( residues(counter) )
residues = temp(1:counter)
deallocate( temp )

end subroutine Identify_Residues 
!
!
!
!==============================================
subroutine get_bond_matrix( sys , total_bonds )
!==============================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_bonds

! local variables ...
integer :: i , j , k , l , N_atoms 
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( bond_matrix(N_atoms,N_atoms), source=.false. )

do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, sys%Total_conect
      flag1 = ( i == InputIntegers(k,1) )
      if( flag1 .eqv. .true.) then
        do l = 2, 5
          flag2 = ( j == InputIntegers(k,l) )
          if( flag2 .eqv. .true.) bond_matrix(i,j) = .true.
        end do
      end if
    end do
    bond_matrix(j,i) = bond_matrix(i,j)
  end do
end do

total_bonds = (size( pack( bond_matrix(:,:), bond_matrix(:,:) .eqv. .true. )))/ 2

end subroutine get_bond_matrix
!
!
!
!==============================================
subroutine get_angle_matrix( sys , total_angs )
!==============================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_angs

! local variables ...
integer :: i , j , k , m , N_atoms 
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( angle_matrix(N_atoms,N_atoms,N_atoms), source=.false. )

m = 0
do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, N_atoms
      flag1 = ( bond_matrix(i,k) .or. bond_matrix(k,i) ) .eqv. .true.
      flag2 = ( bond_matrix(j,k) .or. bond_matrix(k,j) ) .eqv. .true.
      if( flag1 .and. flag2 ) then
        angle_matrix(i,k,j) = .true.
        angle_matrix(j,k,i) = angle_matrix(i,k,j)
        angle_matrix(i,k,i) = .false.
        m = m + 1
      end if
    end do
  end do
end do
total_angs = m 

end subroutine get_angle_matrix
!
!
!
!===================================================
subroutine get_dihedral_matrix( sys , total_diheds )
!===================================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_diheds

! local variables ...
integer :: i , j , k , l , N_atoms
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( dihedral_matrix(N_atoms,N_atoms,N_atoms,N_atoms),source=.false. )

total_diheds=0
do i = 1, N_atoms
  do j = 1, N_atoms
    do k = 1, N_atoms
      if( angle_matrix(i,j,k) .eqv. .true. ) then
        do l = 1, N_atoms
          flag1 = bond_matrix(l,k) .or. bond_matrix(k,l)
          flag2 = bond_matrix(l,i) .or. bond_matrix(i,l)
          if( flag1 ) then 
            dihedral_matrix(i,j,k,l) = .true.
            dihedral_matrix(i,j,k,j) = .false.
            dihedral_matrix(l,k,j,i) = dihedral_matrix(i,j,k,l) 
          end if
          if( flag2 ) then 
            dihedral_matrix(l,i,j,k) = .true.
            dihedral_matrix(j,i,j,k) = .false.
            dihedral_matrix(k,j,i,l) = dihedral_matrix(l,i,j,k)
          end if
          if(dihedral_matrix(i,j,k,l) .eqv. .true.) total_diheds=total_diheds+1
        end do
      end if
    end do
  end do
end do
total_diheds = total_diheds/2 

end subroutine get_dihedral_matrix
!
!
!
!==================================
subroutine write_seed_itp( aux )
!==================================
implicit none 
type(atomic) , intent(in) :: aux(:)

! local variables ...
logical :: TorF
integer :: i , j , k , l , m , N_atoms , ati , atj , atk , atl , max_bonds
integer , allocatable :: bond_list(:,:) , angle_list(:,:) , dihedral_list(:,:) , tmp_list(:,:)

N_atoms = size(aux)

max_bonds = N_atoms*(N_atoms-1)/2

!-----------------------------------------------------------
! Writing stuff ...  
!-----------------------------------------------------------

allocate(tmp_list( max_bonds , 2 ))
write(10,105) "[ bonds ]"
m = 0
do i = 1, N_atoms - 1
  do j = i, N_atoms
    ati = aux(i) % my_intra_id
    atj = aux(j) % my_intra_id
    if( bond_matrix(ati,atj)  .eqv. .true. ) then
             write(10,102) i, j , 1
             m = m + 1
             tmp_list(m,1) = i
             tmp_list(m,2) = j
             endif
  end do
end do
allocate(bond_list , source=tmp_list(1:m,:))
deallocate(tmp_list)
write(10,*)
write(10,*)

allocate(tmp_list( max_bonds , 3 ))
write(10,106) "[ angles ]"
m = 0
do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, N_atoms
      ati = aux(i) % my_intra_id
      atj = aux(j) % my_intra_id
      atk = aux(k) % my_intra_id
      if( angle_matrix(ati,atk,atj) .eqv. .true. ) then
                write(10,103) i, k, j , 1
                m = m + 1
                tmp_list(m,1) = i
                tmp_list(m,2) = k
                tmp_list(m,3) = j
                endif
    end do
  end do
end do
allocate(angle_list , source=tmp_list(1:m,:))
deallocate(tmp_list)
write(10,*)
write(10,*)

allocate(tmp_list( max_bonds , 4 ))
write(10,107) "[ dihedrals ]"
m = 0
do i = 1 , N_atoms
  do j = 1, N_atoms
    do k = 1, N_atoms
      do l = 1, N_atoms
        if( i < l ) then
          ati = aux(i) % my_intra_id
          atj = aux(j) % my_intra_id
          atk = aux(k) % my_intra_id
          atl = aux(l) % my_intra_id
          if( dihedral_matrix(ati,atj,atk,atl) .eqv. .true. ) then
                       write(10,104) i, j, k, l,3
                       m = m + 1
                       tmp_list(m,1) = i
                       tmp_list(m,2) = j
                       tmp_list(m,3) = k
                       tmp_list(m,4) = l
                       endif
        end if
      end do
    end do
  end do
end do
allocate(dihedral_list , source=tmp_list(1:m,:))
deallocate(tmp_list)
write(10,*)
write(10,*)

close(10)

TorF = Checking_Topology( bond_list , angle_list , dihedral_list )
If( TorF ) then
    Print*, "error detected in Topology , check Topology.log"
    stop
End If

102 format(3I4)
103 format(4I4)
104 format(5I4)
105 FORMAT(a9)
106 FORMAT(a10)
107 FORMAT(a13)

end subroutine write_seed_itp
!
!
!
!================================================================
 function Checking_Topology( bonds , angs , diheds ) result(TorF)
!================================================================
implicit none
integer , intent(in) :: bonds (:,:)
integer , intent(in) :: angs  (:,:)
integer , intent(in) :: diheds(:,:)
logical              :: TorF
 
! local variables ... 
integer               :: i , x , y , z
integer               :: Nbonds , Nangs , Ndiheds , KeyLeft , KeyRight
integer , allocatable :: BondKeys(:) , AngKeys(:)
logical               :: flag

Nbonds  =  size(bonds (:,1)) 
Nangs   =  size(angs  (:,1))
Ndiheds =  size(diheds(:,1))

! checking bonds topology ...
allocate( BondKeys(Nbonds) )
do i = 1 , Nbonds

     x = bonds(i,1)  ;  y = bonds(i,2) 
     BondKeys(i) = PairingFunction( x , y , verbose = .true. ) 

end do

! checking angs topology ...
do i = 1 , Nangs

     flag = .false.

     x = angs(i,1)  ;  y = angs(i,2) 
     KeyLeft = PairingFunction(x,y) 
     If( .not. any(KeyLeft == BondKeys) ) call error_message(i,angs,flag,instance="ang")

     x = angs(i,2)  ;  y = angs(i,3) 
     KeyRight = PairingFunction(x,y) 
     If( .not. any(KeyRight == BondKeys) ) call error_message(i,angs,flag,instance="ang")

     If( KeyLeft == KeyRight ) call error_message(i,angs,flag,instance="ang")

end do

! checking diheds topology ...
allocate( AngKeys(Nangs) )
do i = 1 , Nangs

     x = angs(i,1)  ;  y = angs(i,2)   ;  z = angs(i,3) 
     AngKeys(i) = CantorPairing( x , y , z ) 

end do

do i = 1 , Ndiheds

     flag = .false.

     x = diheds(i,1)  ;  y = diheds(i,2)   ;  z = diheds(i,3) 
     KeyLeft = CantorPairing( x , y , z ) 
     If( .not. any(KeyLeft == AngKeys) ) call error_message(i,diheds,flag,instance="dihed")

     x = diheds(i,2)  ;  y = diheds(i,3)   ;  z = diheds(i,4) 
     KeyRight = CantorPairing( x , y , z ) 
     If( .not. any(KeyRight == AngKeys) ) call error_message(i,diheds,flag,instance="dihed")

end do

! prepare to leave ...
if( done ) then  
    TorF = .true.     ! <==  error detected
    close(10)
else
    TorF = .false.    ! <==  NO error detected
end If

end function Checking_Topology
!
!
!
!
!========================================
 function CantorPairing(i,j,k) result(R)
! 3-tupling Cantor Function ...
! f(i,j,k) = f(k,j,i)
!========================================
implicit none
integer            , intent(in) :: i,j,k

! local variables ... 
integer :: R , L , a , b

! Symmetric Pairing for (i,k)-tuple ...
a = max(i,k)  ;  b = min(i,k)
L = a*(a+1)/2 + b 

! Cantor pairing with the center pairing ...
R = (L+j)*(L+j+1)/2 + L 

end function CantorPairing
!
!
!
!===============================================
 function PairingFunction(i,j,verbose) result(k)
!===============================================
implicit none
integer            , intent(in) :: i,j
logical , optional , intent(in) :: verbose

! local variables ... 
integer :: k , a , b

If( (i == j) .and. present(verbose) ) then
    Print 232, i , j
    stop
end If

! the symmetric pairing satisfies f(i,j)=f(j,i) ...

! Symmetric Cantor Pairing ...
!k = (i+j)*(i+j+1)/2 + (i*j) 

! Symmetric Pairing ...
a = max(i,j)  ;  b = min(i,j)
k = a*(a+1)/2 + b 

232 format(/,1x,'>>>  Degenerate Pairing Function in Topology file.....: ',I4,I4)

end function PairingFunction
!
!
!
!==================================================
 subroutine error_message(i , a , flag , instance ) 
!==================================================
implicit none
integer          , intent(in) :: i
integer          , intent(in) :: a(:,:)
logical          , intent(inout) :: flag
character(len=*) , intent(in) :: instance

If( .not. done ) open (10, file='Topology.log', status='unknown')
done = .true.

If( flag == .true. ) return
select case (instance)

       case("ang")
       write(10,231) a(i,1) , a(i,2) , a(i,3) 

       case("dihed")
       write(10,233) a(i,1) , a(i,2) , a(i,3)  , a(i,4) 

end select
flag = .true.

231 format(/,1x,'>>> Error detected in Toplogy file .....: Angle (',I4,',',I4,',',I4,')' )
233 format(/,1x,'>>> Error detected in Toplogy file .....: Dihedral (',I4,',',I4,',',I4,',',I4,')' )

end subroutine error_message
!
!
!
end module GMX_routines
