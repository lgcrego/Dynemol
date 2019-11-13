module GMX_routines

use types_m
use Read_Parms          , only : MMSymbol_2_Symbol , Symbol_2_AtNo
use RW_routines         , only : Initialize_System
use diagnosis_m
use FUNCTION_routines   , only : res , Solvent_residue_groups 
use Topology_routines   , only : connect , dump_topol

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
integer                             ::  i , j , fragment_nresid
character(1)                        :: answer
logical             , parameter     ::  BACK = .TRUE. 

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
integer ::  i , j , k , nr 

!----------------------------------------------
!     generate pdb file for GROMACS
!----------------------------------------------

If( present(title) ) then
    write(4,'(A96)') title
else
    OPEN(unit=4,file='seed.pdb',status='unknown')
    write(4,6) 'COMPND' , sys%Surface_Characteristics
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
type(universe)                      , intent(inout) :: sys

! local variables ...
integer                             :: i , n , N_of_resid_atoms , N_of_resid_groups , N_of_resid_elements
type(atomic)        , allocatable   :: helper(:)
character(len=3)    , allocatable   :: residues(:)

CALL Identify_Residues( sys , residues )

!----------------------------------------------
!   generate .ipt files for GROMACS
!----------------------------------------------

allocate( helper(sys%N_of_atoms) )

OPEN(unit=10,file='seed.itp',status='unknown')

do n = 1 , size(residues)

    N_of_resid_atoms    =  count ( sys%atom%resid == residues(n) )

    if( n > 1 ) then
        N_of_resid_groups = maxval( sys%atom%nresid , sys%atom%resid == residues(n) ) - maxval( sys%atom%nresid , sys%atom%resid == residues(n-1) )  
    else
        N_of_resid_groups = maxval( sys%atom%nresid , sys%atom%resid == residues(n) ) 
    end if

    N_of_resid_elements =  N_of_resid_atoms / N_of_resid_groups

    helper = pack( sys%atom , sys%atom%resid == residues(n) , helper )  

!    do i = 1 , N_of_resid_elements
    do i = 1 , N_of_resid_atoms

        write(10,5) i                               ,  &        ! <== serial number within the residue
                    helper(i)%Symbol                ,  &        ! <== force field descriptor of atom type
                    helper(i)%nresid                ,  &        ! <== residue identifier
                    helper(i)%resid                 ,  &        ! <== residue name
                    helper(i)%MMSymbol              ,  &        ! <== atom type
                    helper(i)%nrcg                  ,  &        ! <== charge group
                    helper(i)%charge                ,  &        ! <== charge of atom type       
                    helper(i)%mass                              ! <== mass of chemical element 
                    
    end do
    write(10,*) '$$$'
end do
close(10)

deallocate( helper )

5 FORMAT(i6,a8,i8,a7,a8,i8,F9.4,F9.4)


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
real*8                          :: distance , total_charge , total_charge_of_group
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
integer        :: i , j , N_of_elements , N_of_atoms ,  iptr
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
integer                         :: i , j , indx , ioerr , useless , N_of_atoms
character(len=50)               :: dumb
character(len=5)                :: MMSymbol_char
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
        if ( keyword == "MASTER" .or. keyword == "END" ) exit
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
        if( keyword == "CRYST1" ) then
            do i = 1 , system%N_of_atoms
                read(3,115)  MMSymbol_char                      ,  &    ! <== atom type
                             system%atom(i)%resid               ,  &    ! <== residue name
                             system%atom(i)%nresid              ,  &    ! <== residue sequence number
                             (system%atom(i)%xyz(j) , j=1,3)    ,  &    ! <== xyz coordinates 
                             system%atom(i)%symbol                      ! <== chemical element symbol

                system%atom(i)%MMSymbol = adjustl(MMSymbol_char)

            end do
        end if
        if ( keyword == "MASTER" .or. keyword == "END") exit
    end do

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
105 format(a6)
110 format(t8, i4)
115 FORMAT(t12,a5,t18,a4,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3,t77,a2)

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
integer ::  i , j , k , nr 

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
end module GMX_routines
