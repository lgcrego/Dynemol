 module Babel_routines_m

    use type_m                  
    use Allocation_m            , only : Allocate_UnitCell

    PUBLIC :: Symbol_2_AtNo ,       &
              Identify_Fragments ,  &
              AtNo_2_Symbol ,       &
              MMSymbol_2_Symbol ,   &
              Identify_Residues ,   &
              Pack_residues ,       &
              Sort_nr ,             &
              Center_of_Gravity ,   &
              Initialize_System ,   &
              TO_UPPER_CASE

    private

    interface Initialize_System
        module procedure Initialize_Universe
        module procedure Initialize_Structure
    end interface

    interface Symbol_2_AtNo
        module procedure Sym_2_AtNo_TRJ
        module procedure Sym_2_AtNo_XYZ
    end interface

    interface AtNo_2_Symbol
        module procedure AtNo_2_Symbol_Uni
        module procedure AtNo_2_Symbol_Struct
    end interface

    interface Identify_Fragments
        module procedure Identify_Fragments_Structure
        module procedure Identify_Fragments_Universe
        module procedure Identify_Fragments_Basis
    end interface

contains
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

    select case( adjustl(a(i)%symbol) )
        case( 'H' ) 
            a(i)%AtNo = 1 
        case( 'LI','Li' ) 
            a(i)%AtNo = 3 
        case( 'BE','Be' ) 
            a(i)%AtNo = 4 
        case( 'C' ) 
            a(i)%AtNo = 6 
        case( 'N' ) 
            a(i)%AtNo = 7 
        case( 'O' ) 
            a(i)%AtNo = 8 
        case( 'F' ) 
            a(i)%AtNo = 9 
        case( 'AL','Al' ) 
            a(i)%AtNo = 13 
        case( 'SI','Si' ) 
            a(i)%AtNo = 14 
        case( 'P','p' ) 
            a(i)%AtNo = 15 
        case( 'S','s' ) 
            a(i)%AtNo = 16 
        case( 'CL','Cl' ) 
            a(i)%AtNo = 17 
        case( 'TI','Ti' ) 
            a(i)%AtNo = 22 
        case( 'MN','Mn' ) 
            a(i)%AtNo = 25 
        case( 'CU','Cu' ) 
            a(i)%AtNo = 29 
        case( 'ZN','Zn' ) 
            a(i)%AtNo = 30 
        case( 'SE','Se' ) 
            a(i)%AtNo = 34 
        case( 'BR','Br' ) 
            a(i)%AtNo = 35 
        case( 'Ru' ) 
            a(i)%AtNo = 44 
        case( 'CD','Cd' )
            a(i)%AtNo = 48 
        case( 'I' ) 
            a(i)%AtNo = 53 
        case( 'Pb' ) 
            a(i)%AtNo = 82
        case( '$$' ) 
            a(i)%AtNo = 90      ! <== special atoms
        case default
            print*, ' >> unknown atom found (1); execution terminated  << : ', a(i)%symbol , i
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

    select case( adjustl(a%symbol(i)) )
        case( 'H' ) 
            a%AtNo(i) = 1 
        case( 'LI','Li' ) 
            a%AtNo(i) = 3 
        case( 'C' ) 
            a%AtNo(i) = 6 
        case( 'N' ) 
            a%AtNo(i) = 7 
        case( 'O' ) 
            a%AtNo(i) = 8 
        case( 'F' ) 
            a%AtNo(i) = 9 
        case( 'AL','Al' ) 
            a%AtNo(i) = 13 
        case( 'SI','Si' ) 
            a%AtNo(i) = 14 
        case( 'P','p' ) 
            a%AtNo(i) = 15 
        case( 'S','s' ) 
            a%AtNo(i) = 16 
        case( 'CL','Cl' ) 
            a%AtNo(i) = 17 
        case( 'TI','Ti' ) 
            a%AtNo(i) = 22 
        case( 'MN','Mn' ) 
            a%AtNo(i) = 25 
        case( 'CU','Cu' ) 
            a%AtNo(i) = 29 
        case( 'ZN','Zn' ) 
            a%AtNo(i) = 30 
        case( 'SE','Se' ) 
            a%AtNo(i) = 34 
        case( 'Ru' ) 
            a%AtNo(i) = 44 
        case( 'CD','Cd' ) 
            a%AtNo(i) = 48 
        case( 'I' ) 
            a%AtNo(i) = 53 
        case( 'Pb' ) 
            a%AtNo(i) = 82
        case( '$$' ) 
            a%AtNo(i) = 90          ! <== special atoms ... 
        case default
            print*, ' >> unknown atom found (2); execution terminated << : ', a%symbol(i) , i
            stop
    end select

END DO

end subroutine Sym_2_AtNo_XYZ
!
!
!
!===============================
 subroutine AtNo_2_Symbol_Uni(a)
!===============================
implicit none
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(a(i)%Atno)
        case( 1 ) 
            a(i)%Symbol = 'H'
        case( 3 ) 
            a(i)%Symbol = 'Li'
        case( 4 ) 
            a(i)%Symbol = 'Be'
        case( 6 ) 
            a(i)%Symbol = 'C'
        case( 7 ) 
            a(i)%Symbol = 'N'
        case( 8 ) 
            a(i)%Symbol = 'O'
        case( 9 ) 
            a(i)%Symbol = 'F'
        case( 13 ) 
            a(i)%Symbol = 'Al'
        case( 14 ) 
            a(i)%Symbol = 'Si'
        case( 15 ) 
            a(i)%Symbol = 'P '
        case( 16 ) 
            a(i)%Symbol = 'S '
        case( 17 ) 
            a(i)%Symbol = 'Cl '
        case( 22 ) 
            a(i)%Symbol = 'Ti '
        case( 25 ) 
            a(i)%Symbol = 'Mn'
        case( 29 ) 
            a(i)%Symbol = 'Cu'
        case( 30 ) 
            a(i)%Symbol = 'Zn'
        case( 34 ) 
            a(i)%Symbol = 'Se'
        case( 44 ) 
            a(i)%Symbol = 'Ru'
        case( 48 ) 
            a(i)%Symbol = 'Cd'
        case( 53 ) 
            a(i)%Symbol = 'I'
        case( 82 ) 
            a(i)%Symbol = 'Pb'
        case( 90 ) 
            a(i)%Symbol = '$$'         ! <== special atoms 
        case default
            print*, ' >> unknown atom found (3); execution terminated << : ', a(i)%AtNo , i
            stop
    end select

 END DO

 end subroutine AtNo_2_Symbol_Uni
!
!
!
!==================================
 subroutine AtNo_2_Symbol_Struct(a)
!==================================
implicit none
type(structure) , intent(inout) :: a

! local variables ...
integer :: i

 DO i = 1 , a%atoms

    select case( a%Atno(i) )
        case( 1 ) 
            a%Symbol(i) = 'H'
        case( 3 ) 
            a%Symbol(i) = 'Li'
        case( 4 ) 
            a%Symbol(i) = 'Be'
        case( 6 ) 
            a%Symbol(i) = 'C'
        case( 7 ) 
            a%Symbol(i) = 'N'
        case( 8 ) 
            a%Symbol(i) = 'O'
        case( 9 ) 
            a%Symbol(i) = 'F'
        case( 13 ) 
            a%Symbol(i) = 'Al'
        case( 14 ) 
            a%Symbol(i) = 'Si'
        case( 15 ) 
            a%Symbol(i) = 'P '
        case( 16 ) 
            a%Symbol(i) = 'S '
        case( 17 ) 
            a%Symbol(i) = 'Cl '
        case( 22 ) 
            a%Symbol(i) = 'Ti '
        case( 25 ) 
            a%Symbol(i) = 'Mn'
        case( 29 ) 
            a%Symbol(i) = 'Cu'
        case( 30 ) 
            a%Symbol(i) = 'Zn'
        case( 34 ) 
            a%Symbol(i) = 'Se'
        case( 44 ) 
            a%Symbol(i) = 'Ru'
        case( 48 ) 
            a%Symbol(i) = 'Cd'
        case( 53 ) 
            a%Symbol(i) = 'I'
        case( 82 ) 
            a%Symbol(i) = 'Pb'
        case( 90 ) 
            a%Symbol(i) = '$$'          ! <== special atoms ...
        case default
            print*, ' >> unknown atom found (4); execution terminated << : ', a%AtNo(i) , i
            stop
    end select

 END DO

 end subroutine AtNo_2_Symbol_Struct
!
!
!
!===============================
 subroutine MMSymbol_2_Symbol(a)
!===============================
implicit none
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer             :: i
character(len=1)    :: element1
character(len=2)    :: element2

 DO i = 1 , size(a)

    write( element1,'(A1)' ) adjustl( a(i)%MMSymbol )

    select case( element1 )
        case( 'C' ) 
            a(i)%Symbol = 'C' 
        case( 'N' ) 
            a(i)%Symbol = 'N' 
        case( 'O' ) 
            a(i)%Symbol = 'O' 
        case( 'F' ) 
            a(i)%Symbol = 'F' 
        case( 'H' ) 
            a(i)%Symbol = 'H' 
        case( 'I' ) 
            a(i)%Symbol = 'I' 
        case( 'P' ) 
            a(i)%Symbol = 'P' 
    end select

    write( element2,'(A2)' ) adjustl( a(i)%MMSymbol )

    select case( element2 )
        case( 'Ow','Ox','OxH','Om','O_3','OH' )
            a(i)%Symbol = 'O' 

        case( 'Hw','HB','HC', 'GH' )
            a(i)%Symbol = 'H' 

        case( 'Ix','Ic' )
            a(i)%Symbol = 'I' 

        case( 'CT','YC','CM','C=','CC','CS','CTr','CA' ) 
            a(i)%Symbol = 'C' 

        case( 'YN','NTr','Nx' ) 
            a(i)%Symbol = 'N' 

        case( 'Al' ) 
            a(i)%Symbol = 'Al' 

        case( 'Si','SI' ) 
            a(i)%Symbol = 'Si' 

        case( 'Ti','TI' ) 
            a(i)%Symbol = 'Ti' 

        case( 'Li' ) 
            a(i)%Symbol = 'Li' 

        case( 'Be' ) 
            a(i)%Symbol = 'Be' 

        case( 'Se','SE' ) 
            a(i)%Symbol = 'Se' 

        case( 'Ru' ) 
            a(i)%Symbol = 'Ru' 

!        case( 'Cd','CD' ) 
!            a(i)%Symbol = 'Cd' 

        case( 'SS','S', 'SG', 'SD','SA' ) 
            a(i)%Symbol = 'S' 

        case( 'Pb','PB' ) 
            a(i)%Symbol = 'Pb' 

        case( '$$') 
            a(i)%Symbol = '$$'           !  <==special atoms ...

    end select

 END DO

end subroutine MMSymbol_2_Symbol
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

If( allocated(a%list_of_fragments) ) deallocate( a%list_of_fragments )

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

If( allocated(a%list_of_fragments) ) deallocate( a%list_of_fragments )

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
!================================================
 subroutine Identify_Fragments_Basis ( a , basis)
!================================================
implicit none
type(structure)  , intent(inout) :: a
type(STO_basis)  , intent(in)    :: basis(:)

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag

If( allocated(a%list_of_fragments) ) deallocate( a%list_of_fragments )

allocate( temp(a % atoms) )

temp(1) = basis(1) % fragment
counter = 1

do i = 1 , size(basis)

    flag = .true.
    do j = 1 , counter
     flag = flag .AND. ( temp(j) /= basis(i)%fragment )
    end do

    if( flag ) then
        counter = counter + 1
        temp(counter) = basis(i)%fragment
    end if

end do

! build list of fragments in a ...
allocate( a%list_of_fragments(counter) )
a%list_of_fragments = temp(1:counter)
deallocate( temp )

end subroutine Identify_Fragments_Basis
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
!================================================
 subroutine Pack_Residues( a , list_of_residues )
!================================================
implicit none
type(atomic)  , allocatable  ,  intent(inout) :: a(:)
character(*)                 , intent(in)     :: list_of_residues(:)

! local variables ...
integer                     :: i , N , first , last
type(atomic) , allocatable  :: temp(:) , buffer(:)

allocate( temp(size(a)) )

first = 1
last  = 0

do i = 1 , size(list_of_residues)

    N = count( a%residue == list_of_residues(i) )
    allocate( buffer(N) )

    buffer = pack( a , a%residue == list_of_residues(i) , buffer )

    last = last + N

    temp(first:last) = buffer

    first = last + 1

    deallocate(buffer)

end do

CALL move_alloc( from=temp , to=a )

end subroutine Pack_Residues
!
!
!
!=======================
 subroutine Sort_nr( a )
!=======================
 implicit none
 type(structure)  , intent(inout)  :: a

! local variables ... 
integer                 :: i ,  j , size_nr_list , last_nr
integer , allocatable   :: nr_list(:) , indx(:)
logical , allocatable   :: flag(:)

! check whether a%nr is already sorted ...
allocate( flag(size(a%nr)) , source = .false. )
forall( i = 2 : size(a%nr) ) flag(i) =  (a%nr(i) < a%nr(i-1)) .OR. (a%nr(i)-a%nr(i-1)) > 1
If( .NOT. any(flag) ) return

! the a%nr are not ordered .OR. there's a gap...
last_nr = 0

! pack => sort => reset a%nr ...
do i = 1 , size(a%list_of_residues)

    size_nr_list = count( a%residue == a%list_of_residues(i) ) 
    allocate( nr_list(size_nr_list) )
    allocate( indx   (size_nr_list) )

    ! pack ...
    nr_list = pack( a%nr , a%residue == a%list_of_residues(i) , nr_list ) 

    ! sort ...
    indx(1) = last_nr + 1
    do j = 2 , size_nr_list
        If( nr_list(j) == nr_list(j-1) ) then
            indx(j) = indx(j-1)
        else
            indx(j) = indx(j-1) + 1
        end If
    end do

    ! reset in the order of indx ...
    a%nr = unpack( indx , a%residue == a%list_of_residues(i) , a%nr )

    last_nr = maxval( indx )
            
    deallocate( nr_list , indx )

end do

end subroutine Sort_nr
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
integer :: i , j , n , n1 , n2 , mol_atoms
real*8  :: total_valence

! initial position of S fragments in trj%atom array ...
n = minloc( trj%atom%fragment , 1 , trj%atom%fragment == "S" ) 

! center of gravity ; the solvent molecules must be contiguous ...
do i = 1 , trj%N_of_Solvent_Molecules 

    mol_atoms = trj%solvent(i)%N_of_Atoms 

    n1 = n
    n2 = n+mol_atoms-1

    total_valence = sum( trj % atom(n1:n2) % Nvalen )

    forall( j=1:3 ) 

        trj % solvent(i) % CG(j) = sum( trj % atom(n1:n2) % xyz(j) ) / trj % solvent(i) % N_of_Atoms

        trj % solvent(i) % CC(j) = sum( trj % atom(n1:n2) % xyz(j) * trj % atom(n1:n2) % Nvalen ) / total_valence

    end forall

    trj % solvent(i) % nr      = trj % atom(n) % nr    
    trj % solvent(i) % residue = trj % atom(n) % residue

    n = n + mol_atoms

end do

end subroutine Center_of_Gravity
!
!
!
!==================================
subroutine Initialize_Universe( a )
!==================================
implicit none
type(universe)  :: a

! local variables ...
integer :: i

forall( i=1:3 ) 
    a % atom % xyz(i)  = 0.d0
    a % atom % TorF(i) = "X"
end forall

a % atom % mass               = 0.d0
a % atom % charge             = 0.d0
a % atom % solvation_hardcore = 3.0d0
a % atom % hardcore           = 3.0d0
a % atom % V_shift            = 0.d0   
a % atom % AtNo               = 0
a % atom % nr                 = 0
a % atom % residue            = "XXX"
a % atom % Symbol             = "XX"
a % atom % MMsymbol           = "XXX"
a % atom % QMMM               = "QM"      ! <== default is Quantum Mechanics for quantum dynamics calcs. ...
a % atom % fragment           = "X"
a % atom % solute             = .false.
a % atom % DPF                = .false.
a % atom % El                 = .false.
a % atom % Hl                 = .false.
a % atom % flex               = .false.

a % N_of_Surface_Atoms        = 0
a % N_of_Solvent_Atoms        = 0
a % N_of_Solvent_Molecules    = 0

end subroutine Initialize_Universe
!
!
!
!===================================
subroutine Initialize_Structure( a )
!===================================
implicit none
type(structure)  :: a

! local variables ...
integer :: i

forall( i=1:3 ) 
    a % coord(:,i)  = 0.d0
end forall

a % solvation_hardcore        = 3.0d0
a % hardcore                  = 3.0d0
a % V_shift                   = 0.0d0
a % AtNo                      = 0
a % nr                        = 1
a % residue                   = "XXX"
a % Symbol                    = "XX"
a % MMsymbol                  = "XXX"
a % QMMM                      = "QM"       ! <== default is Quantum Mechanics for quantum dynamics calcs. ...
a % fragment                  = "X"
a % solute                    = .false.
a % DPF                       = .false.
a % El                        = .false.
a % Hl                        = .false.
a % flex                      = .false.
a % Nvalen                    = 0
a % polar                     = 0
a % N_of_electrons            = 0
a % N_of_Solvent_Molecules    = 0

end subroutine Initialize_Structure
!
!
!
end module Babel_routines_m
