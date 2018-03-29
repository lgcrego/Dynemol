! Convert gmx data to mdflex program :: verify gmx format in IO_FORMATS 

module namd2mdflex

use constants_m
use for_force
use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles, DefinePairs, DefineMorse, debug_MM
use MM_tuning_routines     , only : SpecialBonds, SpecialAngs
use NonBondPairs           , only : Identify_NonBondPairs
use Babel_routines_m       , only : TO_UPPER_CASE

private
 
public :: prm2mdflex, psf2mdflex, SpecialPairs, SpecialPairs14, SpecialMorse 

    ! module variables ...
    character(3)      , allocatable  , save  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:), ImproperSymbols(:,:)
    real*8            , allocatable  , save  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:), ImproperParameters(:,:)
    type(DefinePairs) , allocatable :: SpecialPairs(:)
    type(DefinePairs) , allocatable :: SpecialPairs14(:)
    type(DefineMorse) , allocatable :: SpecialMorse(:)

contains
!
!
!
!================================================
 subroutine psf2mdflex( MM , atom , species , FF)
!================================================
implicit none
type(MM_system)     , intent(in)    :: MM
type(MM_atomic)     , intent(inout) :: atom(:)
type(MM_atomic)     , intent(inout) :: FF(:)
type(MM_molecular)  , intent(inout) :: species(:)

! local variables ...
character(15)   , allocatable   :: InputChars(:,:)
real*8          , allocatable   :: InputReals(:,:)
integer         , allocatable   :: InputIntegers(:,:) 
character(1)                    :: dummy_char
character(18)                   :: keyword 
character(10)                   :: string
character(200)                  :: line 
integer                         :: i1 , i2 , i3 , sp , nr , n
integer                         :: i , j , k , a , ioerr , dummy_int , counter , N_of_atoms
integer                         :: Nbonds , Nangs , Ndiheds , Nimpropers , Nbonds14 

allocate( InputChars    ( 20000 , 10 )                   )
allocate( InputReals    ( 20000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 20000 , 10 ) , source = I_zero )

! Reading different '.itp' species files ...
counter = 0
do a = 1 , MM % N_of_species

    string = species(a) % residue // '.psf'

    open(33, file=string, status='old',iostat=ioerr,err=101)

        101 if( ioerr > 0 ) then
            print*, string,' file not found; terminating execution' ; stop
        end if

        write(*,'(/2a9)',advance='no') "Reading ", string

        ! start reading the molecular structure of species(a) ...
        do
            read(33,100) keyword
            if( verify( "!NATOM" , keyword ) == 0 ) exit
        end do

        allocate( species(a) % atom ( species(a) % N_of_atoms ) )

        do i = 1 , species(a)%  N_of_atoms
            read(33, '(A)', iostat=ioerr) line
            read(line,*,iostat=ioerr) species(a) % atom(i) % my_id ,      &
                                      dummy_char ,                        &
                                      species(a) % atom(i) % nr ,         &
                                      species(a) % atom(i) % residue ,    &
                                      species(a) % atom(i) % EHSymbol ,   &
                                      species(a) % atom(i) % MMSymbol ,   &
                                      species(a) % atom(i) % MM_charge ,  &
                                      species(a) % atom(i) % mass

            species(a) % atom(i) % MMSymbol   = adjustr(species(a) % atom(i) % MMSymbol)
            species(a) % atom(i) % my_species = a
            species(a) % my_species           = a
            species(a) % atom(i) % flex       = species(a) % flex

            ! this is the standard; atomic flexibity can also be defined @ ad_hoc_MM_tuning ...    
            where( atom % my_species == a ) atom % flex = species(a) % flex

            counter = counter + 1
            FF(counter) % my_species = a
            FF(counter) % my_id      = species(a) % atom(i) % my_id
            FF(counter) % residue    = species(a) % atom(i) % residue
            FF(counter) % nr         = species(a) % atom(i) % nr     
            FF(counter) % EHSymbol   = species(a) % atom(i) % EHSymbol
            FF(counter) % MMSymbol   = species(a) % atom(i) % MMSymbol
            FF(counter) % MM_charge  = species(a) % atom(i) % MM_charge 
            FF(counter) % mass       = species(a) % atom(i) % mass 
        end do 
        backspace(33)

        N_of_atoms = species(a) % N_of_atoms

        ! convert residues to upper case ...
        forall( i=1:N_of_atoms ) species(a)% atom(i)% residue = TO_UPPER_CASE( species(a)% atom(i)% residue )

        i = 1
        do

            if( trim(atom(i) % residue) == trim(species(a) % atom(1) % residue) ) then
                atom(i:i+N_of_atoms-1) % MM_charge = species(a) % atom(:N_of_atoms) % MM_charge
                i = i + N_of_atoms
            else
                i = i + 1
            end if

            if( i > size(atom) ) exit

        end do
        rewind 33
!==============================================================================================
        ! Bonding parameters :: reading ...
        do
          read(33,100) keyword
          if( verify( "!NBOND:" , keyword ) == 0 ) exit
        end do
        backspace(33)
        read(33,*) Nbonds

        do k = 1 , ceiling(Nbonds/four)-1
          read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*4+n,j) , j=1,2 ) , n=1,4 )
        end do 
        read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*4+n,j) , j=1,2 ) , n=1,merge(4,mod(NBonds,4),mod(NBonds,4)==0) )

        species(a) % Nbonds = Nbonds

        allocate( species(a) % bonds      ( Nbonds , 2 ) )
        allocate( species(a) % funct_bond ( Nbonds     ) )
        allocate( species(a) % bond_type  ( Nbonds     ) )

        forall(i=1:2) species(a) % bonds(:Nbonds,i) = InputIntegers(:Nbonds,i)

        species(a) % funct_bond(:Nbonds) = "1"
        
        do i = 1 , Nbonds 
           select case ( species(a) % funct_bond(i) )
            case( "1" )  
                species(a) % bond_type(i) = "harm" 
            case( "3" )
                species(a) % bond_type(i) = "Mors"     
           end select 
        end do
        rewind 33
!==============================================================================================
        ! Angle parameters :: reading ...
        do
          read(33,100) keyword
          if( verify( "!NTHETA:" , keyword ) == 0 ) exit
        end do
        backspace(33)
        read(33,*) Nangs

        InputIntegers = I_zero

        do k = 1 , ceiling(Nangs/three)-1
          read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*3+n,j) , j=1,3 ) , n=1,3 )
        end do 
        read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*3+n,j) , j=1,3 ) , n=1,merge(3,mod(NAngs,3),mod(NAngs,3)==0) )

        species(a) % Nangs = Nangs

        allocate( species(a) % angs        ( Nangs , 3 ) )
        allocate( species(a) % funct_angle ( Nangs     ) )
        allocate( species(a) % angle_type  ( Nangs     ) )

        forall(i=1:3) species(a) % angs(:Nangs,i) = InputIntegers(:Nangs,i)

        species(a) % funct_angle(:Nangs) = "1"

        do i = 1 , Nangs
           select case ( species(a) % funct_angle(i) )
            case( "1" )
                species(a) % angle_type(i) = "harm"
            case( "5" )
                species(a) % angle_type(i) = "urba"
           end select
        end do
        rewind 33
!==============================================================================================
        ! Dihedral parameters :: reading ...

        !=========================================================================
        !reading normal dihedrals ...
        do
          read(33,100,iostat=ioerr) keyword
          if( verify( "!NPHI:" , keyword ) == 0 ) exit
        end do
        backspace(33)
        read(33,*) Ndiheds

        InputIntegers = I_zero

        do k = 1 , ceiling(Ndiheds/two)-1
          read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*2+n,j) , j=1,4 ) , n=1,2 )
        end do 
        read(33 , * , iostat=ioerr )  ( ( InputIntegers((k-1)*2+n,j) , j=1,4 ) , n=1,merge(2,mod(Ndiheds,2),mod(Ndiheds,2)==0) )

        rewind 33
        !=========================================================================
        !reading improper dihedrals ...
        do
          read(33,100,iostat=ioerr) keyword
          if( verify( "!NIMPHI:" , keyword ) == 0 ) exit
        end do
        backspace(33)
        read(33,*) Nimpropers

        do k = 1 , ceiling(Nimpropers/two)-1
          read(33 , * , iostat=ioerr )  ( ( InputIntegers(Ndiheds+(k-1)*2+n,j) , j=1,4 ) , n=1,2 )
        end do 
        read(33 , * , iostat=ioerr )  ( ( InputIntegers(Ndiheds+(k-1)*2+n,j) , j=1,4 ) , n=1,merge(2,mod(Ndiheds,2),mod(Ndiheds,2)==0) )
        !=========================================================================

        species(a)% Ndiheds = Ndiheds + Nimpropers

        allocate( species(a)% diheds      ( species(a)% Ndiheds , 4 ) )
        allocate( species(a)% funct_dihed ( species(a)% Ndiheds     ) )

        ! store normal dihedrals ...
        forall(i=1:4) species(a)% diheds( :Ndiheds , i ) = InputIntegers( :Ndiheds , i )
        species(a)% funct_dihed( :Ndiheds ) = 1

        ! store improper dihedrals ...
        forall(i=1:4) species(a)% diheds( Ndiheds+1:Ndiheds+Nimpropers , i ) = InputIntegers( :Ndiheds+Nimpropers , i )
        ! Está dando erro: funct_dihed está alocado para apenas Ndiheds
        !species(a)% funct_dihed( :Ndiheds+Nimpropers ) = 4

        ! define species(a) % dihedral_type ...
        CALL define_DihedralType( species(a) , species(a)% Ndiheds )

        rewind 33
!==============================================================================================
            ! Pairs 1-4 parameters :: reading ...
        do
            read(33,100,iostat=ioerr) keyword
            if( trim(keyword) == "[ pairs ]" .OR. ioerr /= 0 ) exit
        end do

        if( trim(keyword) == "[ pairs ]" ) then

            InputIntegers = I_zero
            i = 1
            read_loop5: do
                read(33, '(A)', iostat=ioerr) line
                if ( ioerr /= 0 ) exit read_loop5
                read(line,*,iostat=ioerr) InputChars(i,1)
                if( index(InputChars(i,1),";") /= 0 ) cycle read_loop5
                if( trim(InputChars(i,1)) == "[  "  ) exit
                if( ioerr > 0  ) exit
                if( ioerr /= 0 ) cycle read_loop5
                read(line,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,2 ) , InputReals(i,1)

                i = i + 1
            end do read_loop5
            backspace(33)

            Nbonds14 = i - 1
            species(a) % Nbonds14 = Nbonds14

            allocate( species(a) % bonds14 ( Nbonds14 , 2 ) )

            forall(i=1:2) species(a) % bonds14(:Nbonds14,i) = InputIntegers(:Nbonds14,i)

        else

            Nbonds14 = 0    

        end if

!==============================================================================================

        If( species(a) % Nbonds /= 0 ) then 

            CALL Identify_NonBondPairs( species , a )

        else 

            ! Intermediate variable ... 
            allocate( species(a) % IntraLJ ( (species(a) % N_of_Atoms * (species(a) % N_of_Atoms-1))/2, 2 ) , source = I_zero )

            k = 1
            do i = 1 , species(a) % N_of_Atoms - 1

                do j = i + 1, species(a) % N_of_Atoms
                    species(a) % IntraLJ(k,1) = i
                    species(a) % IntraLJ(k,2) = j
                    k = k + 1
                end do

            end do
            species(a) % NintraLJ = size( species(a) % IntraLJ(:,2) )
             
        end if

!==============================================================================================

    close(33)

    write(*,'(a9)') " << done "

end do

FF % residue  = adjustl(FF % residue)
FF % Symbol   = adjustl(FF % Symbol)
FF % MMSymbol = adjustl(FF % MMSymbol)

! passing MMSymbol from FF to atom ...
i1 = 1
    do nr = 1 , atom(MM%N_of_atoms) % nr
    sp = atom(i1) % my_species
    i3 = count(atom%nr==nr)
    i2 = i1 + (i3-1)
    atom(i1:i2)%MMSymbol = pack(FF % MMSymbol , FF % my_species == sp)
    i1 = i2+1
end do

atom % MMSymbol = adjustl(atom % MMSymbol)

deallocate( InputChars , InputIntegers )

100 format(a18)

end subroutine psf2mdflex
!
!
!
!=================================================
 subroutine prm2mdflex( MM , atom , species , FF )
!=================================================
implicit none 
type(MM_molecular)                  , intent(inout) :: species(:)
type(MM_system)                     , intent(inout) :: MM
type(MM_atomic)                     , intent(inout) :: atom(:)
type(MM_atomic)     , allocatable   , intent(inout) :: FF(:)
 
! local variables ...
character(3)    , allocatable   :: InputChars(:,:) , Input2Chars(:,:)
character(4)    , allocatable   :: funct_bond(:) , funct_angle(:)
real*8          , allocatable   :: InputReals(:,:) , Input2Reals(:,:)
integer         , allocatable   :: InputIntegers(:,:)
integer         , allocatable   :: Dihed_Type(:) , Bond_Type(:) , Angle_Type(:)
integer                         :: a , n , i , j , k , ioerr , dummy_int , N_of_AtomTypes 
integer                         :: NbondsTypes , NangsTypes , NdihedTypes , NImproperTypes, NBondParms, NPairsParms , NMorseParms 
character(3)                    :: dummy_char
character(18)                   :: keyword
character(200)                  :: line
logical                         :: flag1 , flag2 , flag3 , flag4

allocate( InputChars    ( 10000 , 10 )                   )
allocate( Input2Chars   ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( Input2Reals   ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

MM % CombinationRule = 2
MM % fudgeQQ = 1.0
MM % fudgeLJ = 1.0
  
open(33, file='input.prm', status='old', iostat=ioerr, err=10)

!   file error msg ...
    10 if( ioerr > 0 ) stop '"input.prm" file not found; terminating execution'

    ! reading the number of ATOMS ...
    do
        read(33,100) keyword
        if( verify( "ATOMS" , keyword ) == 0 ) exit
    end do

    i=1
    read_loop: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop
        read(line,*,iostat=ioerr) InputChars(i,1)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop
        if( trim(InputChars(i,1)) == "BON"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop
        read(line,*,iostat=ioerr) InputChars(i,1:5) , dummy_int, InputChars(i,6:9) , InputReals(i,1) 

        i = i + 1
    end do read_loop
    InputChars = adjustl(InputChars) 

    backspace(33)
 
    N_of_AtomTypes = i - 1

    MM % N_of_AtomTypes = i - 1

!=====================================================================================
!  reads BONDS ...
    do
        read(33,100) keyword
        if( trim(keyword) == "BONDS" ) exit
    end do
        
    i = 1
    read_loop1: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop1
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop1
        if( trim(InputChars(i,1)) == "ANG"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop1
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,2) , (InputReals(i,j) , j=1,2)

        i = i + 1
    end do read_loop1
    backspace(33)

    NbondsTypes = i - 1

    allocate( BondPairsSymbols    ( NbondsTypes , 2 ) )
    allocate( BondPairsParameters ( NbondsTypes , 3 ) , source = D_zero )
    allocate( funct_bond    ( NbondsTypes ) )

    forall(i=1:2) BondPairsSymbols(:NbondsTypes,i) = InputChars(:NbondsTypes,i)
    forall(i=1:2) BondPairsParameters(:NbondsTypes,i) = InputReals(:NbondsTypes,i)
  
   funct_bond( :NbondsTypes ) = "harm" 
   !do i = 1 , NbondsTypes 
       !BondPairsParameters(i,1) = InputReals(i,2) * factor2 * imol
       !BondPairsParameters(i,2) = InputReals(i,1) * nano_2_angs
   !end do

!=====================================================================================
!  reads ANGLES ...
    do
        read(33,100) keyword
        if( trim(keyword) == "ANGLES" ) exit
    end do

    i = 1
    read_loop2: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop2
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop2
        if( trim(InputChars(i,1)) == "DIH"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop2
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,3) , (InputReals(i,j) , j=1,2 )
       !backspace(33) 

       !select case( InputIntegers(i,1) )
       !    case( 1 ) ! Harmonic potential ...
       !        read(33,*) (dummy_char,k=1,3), dummy_int, (InputReals(i,j), j=1,2)
       !    case( 5 ) ! Urey-Bradley potential ...
       !        read(33,*) (dummy_char,k=1,3), dummy_int, (InputReals(i,j), j=1,4)
       !end select

        i = i + 1
    end do read_loop2
    backspace(33)

    NangsTypes = i - 1

    allocate( AngleSymbols    ( NangsTypes , 3 ) )
    allocate( AngleParameters ( NangsTypes , 4 ) , source = D_zero ) 
    allocate( Angle_Type      ( NangsTypes     ) )
    allocate( funct_angle     ( NangsTypes     ) )
    
    forall(i=1:3) AngleSymbols(:NangsTypes,i)     = InputChars(:NangsTypes,i)
    forall(i=1:2) AngleParameters(:NangsTypes,i)  = InputReals(:NangsTypes,i)
    
   Angle_Type( :NangsTypes ) = 1 
   funct_angle( :NangsTypes ) = "harm"

   !do i = 1 , NangsTypes
   !    select case( Angle_Type(i) )
   !    ! conversion 
   !    ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
   !        case( 1 ) ! Harmonic potential ...
   !            AngleParameters(:NangsTypes,1) = InputReals(:NangsTypes,2) * factor1 * imol
   !            AngleParameters(:NangsTypes,2) = InputReals(:NangsTypes,1) * deg_2_rad
   !            funct_angle(i) = "harm"
   !        case( 5 ) ! Urey-Bradley potential ...
   !            AngleParameters(:NangsTypes,1) = InputReals(:NangsTypes,2) * factor1 * imol
   !            AngleParameters(:NangsTypes,2) = InputReals(:NangsTypes,1) * deg_2_rad
   !            AngleParameters(:NangsTypes,3) = InputReals(:NangsTypes,4) * factor2 * imol
   !            AngleParameters(:NangsTypes,4) = InputReals(:NangsTypes,3) * nano_2_angs
   !            funct_angle(i) = "urba"
   !    end select
   !end do

!=====================================================================================
!  reads DIHEDRALS ...
    InputReals    = D_zero
    InputIntegers = I_zero
    do
        read(33,100) keyword
        if( trim(keyword) == "DIHEDRALS" ) exit
    end do

    i = 1
    read_loop3: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop3
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop3
        if( trim(InputChars(i,1)) == "IMP"    ) exit
        if( trim(InputChars(i,1)) == "CMA"    ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop3
        read(line,*,iostat=ioerr) (InputChars(i,k) , k=1,4) , (InputReals(i,j) , j=1,3)

        i = i + 1
    end do read_loop3
    backspace(33)

    NdihedTypes = i - 1

    allocate( Dihed_Type         ( NdihedTypes     )                   )
    allocate( DihedSymbols       ( NdihedTypes , 4 )                   )
    allocate( DihedParameters    ( NdihedTypes , 3 ) , source = D_zero )

    forall(k=1:4) DihedSymbols(:NdihedTypes,k)    = InputChars(:NdihedTypes,k)
    forall(k=1:3) DihedParameters(:NdihedTypes,k) = InputReals(:NdihedTypes,k)

    Dihed_Type(:NdihedTypes) = 1 

!=====================================================================================
!  reads IMPROPER ... 
    InputReals    = D_zero
    InputIntegers = I_zero
    do
        read(33,100) keyword
        if( trim(keyword) == "IMPROPER" ) exit
    end do

    i = 1
    read_loop4: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop4
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop4
        if( trim(InputChars(i,1)) == "CMA"    ) exit
        if( trim(InputChars(i,1)) == "NON"    ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop4
        read(line,*,iostat=ioerr) (InputChars(i,k) , k=1,4) , (InputReals(i,j) , j=1,3)

        i = i + 1
    end do read_loop4
    backspace(33)

    NImproperTypes = i - 1

   !allocate( Dihed_Type         ( NdihedTypes     )                   )
    allocate( ImproperSymbols       ( NImproperTypes , 4 )                   )
    allocate( ImproperParameters    ( NImproperTypes , 3 ) , source = D_zero )

    forall(k=1:4) ImproperSymbols(:NImproperTypes,k)    = InputChars(:NImproperTypes,k)
    forall(k=1:3) ImproperParameters(:NImproperTypes,k) = InputReals(:NImproperTypes,k)

   !Dihed_Type(:NdihedTypes) = InputIntegers(:NdihedTypes,1)

!=====================================================================================
!  NonBonding parameters :: reading ...
    do
        read(33,100) keyword
        if( trim(keyword(1:9)) == "NONBONDED" ) exit
    end do
    read(33,100)

    i = 1
    read_loop5: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop5
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop5
        if( trim(InputChars(i,1)) == "HBO"  ) exit
        if( trim(InputChars(i,1)) == "END"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop5
        read(line,*, iostat=ioerr) InputChars(i,1) , (InputReals(i,j) , j=1,3)

        i = i + 1
    end do read_loop5
    backspace(33)

    NBondParms = i - 1

    do i = 1 , N_of_AtomTypes
        where( FF % MMSymbol == InputChars(i,1) ) 
            FF % sig = InputReals(i,3)
            FF % eps = abs(InputReals(i,2))
        end where
    end do
    backspace(33)

    ! conversion 
    ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
 !  FF % eps = sqrt( FF % eps * factor1 * imol )
 !  FF % sig = FF % sig * nano_2_angs

    select case( MM % CombinationRule )

        case (2)

            FF % sig = FF % sig / TWO

        case (3)

            FF % sig = sqrt( FF % sig )

    end select

    do i = 1 , size(FF)
        where( atom % MMSymbol == FF(i) % MMSymbol )
            atom % eps = FF(i) % eps
            atom % sig = FF(i) % sig
        end where
    end do

    NPairsParms = i - 1

   ! allocate( SpecialPairs14 ( NPairsParms ) )

   ! forall(i=1) SpecialPairs14(:NPairsParms) % MMSymbols(i) = InputChars(:NPairsParms,i)

   ! SpecialPairs14(:NPairsParms) % Parms(2) = InputReals(:NPairsParms,2)
   ! SpecialPairs14(:NPairsParms) % Parms(1) = InputReals(:NPairsParms,1)

    ! conversion 
    ! factor1 = 1.0d26      <== Factor used to correct the units readed from Gromacs
    !SpecialPairs14(:NPairsParms) % Parms(2) = sqrt( SpecialPairs14(:NPairsParms)%Parms(2) * factor1 * imol )
    !SpecialPairs14(:NPairsParms) % Parms(1) = sqrt( SpecialPairs14(:NPairsParms)%Parms(1) * nano_2_angs    )

!=====================================================================================
!

close(33)

deallocate( InputChars , InputReals , InputIntegers )
deallocate( Input2Chars , Input2Reals )
!=====================================================================================

do a = 1 , MM % N_of_species

!   Assigning to each specie the corresponding parameter ...

    ! Bond parameters ...
    allocate( species(a) % kbond0( species(a) % Nbonds , 3 ) , source = D_zero )

    do k = 1 , NbondsTypes
        do n = 1 , species(a) % Nbonds

            flag1 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) )
            flag2 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) ) .AND. & 
                    ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) )  
            flag3 = ( adjustl(species(a) % bond_type(n)) == adjustl(funct_bond(k)) )

            if ( ( flag1 .OR. flag2 ) .AND. flag3 ) then 
                species(a) % kbond0(n,1) = BondPairsParameters(k,1)
                species(a) % kbond0(n,2) = BondPairsParameters(k,2)
                species(a) % kbond0(n,3) = BondPairsParameters(k,3) 
            end if
        end do
    end do

    if( allocated(SpecialBonds) ) then
        do k = 1, size(SpecialBonds)
            where( species(a) % funct_bond == SpecialBonds(k) % label ) species(a) % kbond0(:,1) = SpecialBonds(k) % kbond0(1) * factor2 * imol
            where( species(a) % funct_bond == SpecialBonds(k) % label ) species(a) % kbond0(:,2) = SpecialBonds(k) % kbond0(2) * nano_2_angs
        end do
    end if

    !=============================================================================
    ! Angle parameters ...
    allocate( species(a) % kang0(species(a) % Nangs , 4 ) , source = D_zero )

    do k = 1 , NangsTypes
        do n = 1 , species(a) % Nangs

            flag1 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,1)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,3)) )

            flag2 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,3)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,1)) )
            flag3 = ( adjustl(species(a) % angle_type(n)) == adjustl(funct_angle(k)) )

            if ( ( flag1 .OR. flag2 ) .AND. flag3 ) then 
                species(a) % kang0(n,1) = AngleParameters(k,1)
                species(a) % kang0(n,2) = AngleParameters(k,2)
                species(a) % kang0(n,3) = AngleParameters(k,3)
                species(a) % kang0(n,4) = AngleParameters(k,4)
            end if

        end do
    end do

    if( allocated(SpecialAngs) ) then
        do k = 1 , size(SpecialAngs)
            ! conversion 
            ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
            where( species(a) % funct_angle == SpecialAngs(k) % label ) species(a) % kang0(:,1) = SpecialAngs(k) % kang0(1) * factor1 * imol
            where( species(a) % funct_angle == SpecialAngs(k) % label ) species(a) % kang0(:,2) = SpecialAngs(k) % kang0(2) * deg_2_rad
        end do
    end if

    !=============================================================================
    ! Dihedral parameters ...
    allocate( species(a) % kdihed0 ( species(a) % Ndiheds , 6 ) , source = D_zero )

    read_loop0: do n = 1 , species(a) % Ndiheds
        do k = 1 , NdihedTypes 

            ! if funct = 1 (cos)
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        
            ! Eq. 4.60 (GMX manual 5.0.5)

            if( species(a) % funct_dihed(n) == 1 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 1 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 1 )


                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 1 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. & 
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 1 )

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    species(a) % kdihed0(n,1:3) = DihedParameters(k,1:3)
                    cycle read_loop0
                end if

            end if

            ! if funct = 2 (harm)
            ! V = 1/2.k ( xi - xi_0 )²
            ! Eq. 4.59 (GMX manual 5.0.5)

            if( species(a) % funct_dihed(n) == 2 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 2 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 2 )

                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 2 ) 

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 2 )

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 ) then
                    ! kdihed0(:,1) = xi_0 (deg)
                    ! kdihed0(:,2) = k_xi [ kJ/(mol.rad²) ]
                    species(a) % kdihed0(n,1:2) = DihedParameters(k,1:2)
                    cycle read_loop0
                end if

            end if

            ! if funct = 3 (cos3)
            ! V = C0 + C1 * cos( phi - 180 ) + C2 * cos^2( phi - 180 ) + C3 * cos^3( phi - 180 ) + C4 * cos^4( phi - 180 ) + C5 * cos^5( phi - 180 ) 
            ! Eq. 4.61 (GMX manual 5.0.5)
            
            if( species(a) % funct_dihed(n) == 3 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 3 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 3 )

                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 3 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 3 )
   
                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 ) then
                    species(a) % kdihed0(n,1:6) = DihedParameters(k,1:6) 
                    cycle read_loop0
                end if

            end if

        end do
    end do read_loop0
    !=============================================================================

end do

if( allocated(BondPairsParameters) ) deallocate( BondPairsParameters )
if( allocated(BondPairsSymbols)    ) deallocate( BondPairsSymbols    )
if( allocated(SpecialBonds)        ) deallocate( SpecialBonds        )
if( allocated(Bond_Type)           ) deallocate( Bond_Type           )
if( allocated(AngleParameters)     ) deallocate( AngleParameters     )
if( allocated(AngleSymbols)        ) deallocate( AngleSymbols        )
if( allocated(SpecialAngs)         ) deallocate( SpecialAngs         )
if( allocated(Angle_Type)          ) deallocate( Angle_Type          )
if( allocated(DihedParameters)     ) deallocate( DihedParameters     )
if( allocated(DihedSymbols)        ) deallocate( DihedSymbols        )
if( allocated(Dihed_Type)          ) deallocate( Dihed_Type          )

100 format(a18)
120  format(4a5,t22,I2,t26,6f14.4)
 
end subroutine prm2mdflex
!
!
!
!=======================================
 subroutine define_DihedralType( a , N )
!=======================================
implicit none
type(MM_molecular)  , intent(inout) :: a
integer             , intent(in)    :: N

! local variables ...
integer :: i

allocate( a % Dihedral_Type ( N ) )

do i = 1 , N

    select case( a % funct_dihed(i) )

        case( 1 ) ! V = k[1 + cos(n.phi - theta)]

            a % dihedral_type(i) = "cos" 

        case( 2 ) ! V = 1/2.k( xi - xi_0 )²

            a % dihedral_type(i) = "harm"

        case( 3 ) ! v = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
            
            a % dihedral_type(i) = "cos3"

        case( 4 ) ! V = k[1 + cos(n.phi - theta)] (improper)

            a % dihedral_type(i) = "imp" 

        case( 9 ) ! V = k[1 + cos(n.phi - theta)] (multiple; charmm FF)

            a % dihedral_type(i) = "chrm"

    end select

    a % dihedral_type(i) = adjustl( a % dihedral_type(i) )

end do

end subroutine define_DihedralType
!
!
!
end  module namd2mdflex
