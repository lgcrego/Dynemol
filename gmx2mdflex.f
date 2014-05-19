! Convert gmx data to mdflex program :: verify gmx format in IO FORMATS 
module gmx2mdflex

use constants_m
use for_force
use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles
use MM_tuning_routines     , only : SpecialBonds, SpecialAngs

private
 
public :: top2mdflex, itp2mdflex

    ! module variables ...
    real*8                           , save  :: fact14
    character(3)     , allocatable   , save  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:)
    real*8           , allocatable   , save  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:)

contains
!
!
!
!=================================================
 subroutine top2mdflex( MM , atom , species , FF )
!=================================================
implicit none 
type(MM_molecular)                  , intent(inout) :: species(:)
type(MM_system)                     , intent(inout) :: MM
type(MM_atomic)                     , intent(inout) :: atom(:)
type(MM_atomic)     , allocatable   , intent(inout) :: FF(:)
 
! local variables ...
type(MM_atomic) , allocatable   :: FF_tmp(:)
character(3)    , allocatable   :: InputChars(:,:)
real*8          , allocatable   :: InputReals(:,:)
integer         , allocatable   :: InputIntegers(:,:)
integer         , allocatable   :: Dihed_Type(:)
real*8                          :: factQQ , dummy_real , theta0 , ktheta0 , fudgeLJ , fudgeQQ
integer                         :: a , n , i , j , k , ioerr , dummy_int , N_of_AtomTypes , NbondsTypes , NangsTypes , NdihedTypes , Nbonds14Types

character(1)                    :: keyword_1
character(3)                    :: dummy_char
character(9)                    :: keyword_9
character(18)                   :: keyword
logical                         :: flag1 , flag2 , flag3 , flag4

allocate( InputChars    ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

forcefield = 2
  
open(33, file='topol.top', status='old', iostat=ioerr, err=10)
    
!   file error msg ...
    10 if( ioerr > 0 ) stop '"topol.top" file not found; terminating execution'

!   reading defaults ...
    CALL skip_lines(33,2)
    read(33,*) dummy_int, MM % CombinationRule, dummy_char, fudgeLJ , fudgeQQ

!=====================================================================================
!   reading he number of [ atomtypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ atomtypes ]" ) exit
    end do
    CALL skip_lines(33,2)

    i=1
    do
        read(33,*,iostat=ioerr) InputChars(i,1) , (InputReals(i,j) , j=1,2) , InputChars(i,2) , (InputReals(i,j) , j=3,4)
        if( ioerr /= 0 ) exit
        i = i + 1
    end do
    InputChars = adjustl(InputChars)

    backspace(33)
 
    N_of_AtomTypes = i - 1

    MM % N_of_AtomTypes = i - 1

    do i = 1 , N_of_AtomTypes
        where( FF % MMSymbol == InputChars(i,1) )
            FF % sig = InputReals(i,3)
            FF % eps = InputReals(i,4)
        end where   
    end do

    ! conversion 
    ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
    FF % eps = sqrt( FF % eps * factor1 * imol )
    FF % sig = FF % sig * nano_2_angs

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

!=====================================================================================
!  NonBonding parameters :: reading ...

!  if there are [ nonbonding parameters ] read here...

!=====================================================================================
!  reads [ bondtypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ bondtypes ]" ) exit
    end do
    CALL skip_lines(33,2)
        
    i = 1
    do
        read(33,*, iostat=ioerr) (InputChars(i,j) , j=1,2) , InputIntegers(i,1) , (InputReals(i,j) , j=1,2)
        if( ioerr /= 0 ) exit
        i = i + 1
    end do
    backspace(33)

    NbondsTypes = i - 1

    allocate( BondPairsSymbols    ( NbondsTypes , 2 ) )
    allocate( BondPairsParameters ( NbondsTypes , 2 ) )
    forall(i=1:2) BondPairsSymbols(:NbondsTypes,i) = InputChars(:NbondsTypes,i)
    BondPairsParameters(:NbondsTypes,1) = InputReals(:NbondsTypes,2) * factor2 * imol
    BondPairsParameters(:NbondsTypes,2) = InputReals(:NbondsTypes,1) * nano_2_angs

!=====================================================================================
!  reads [ angletypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ angletypes ]" ) exit
    end do
    CALL skip_lines(33,1)

    i = 1
    do
        read(33,*, iostat=ioerr) (InputChars(i,j) , j=1,3) , InputIntegers(i,1) , (InputReals(i,j) , j=1,2)
        if( ioerr /= 0 ) exit
        i = i + 1
    end do
    backspace(33)
    backspace(33)

    NangsTypes = i - 1

    allocate( AngleSymbols    ( NangsTypes , 3 ) )
    allocate( AngleParameters ( NangsTypes , 2 ) )
    forall(i=1:3) AngleSymbols(:NangsTypes,i)  = InputChars(:NangsTypes,i)
    
    ! conversion 
    ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
    AngleParameters(:NangsTypes,1) = InputReals(:NangsTypes,2) * factor1 * imol
    AngleParameters(:NangsTypes,2) = InputReals(:NangsTypes,1) * deg_2_rad

!=====================================================================================
!  reads [ dihedraltypes ] ...
    InputReals    = D_zero
    InputIntegers = I_zero
    do
        read(33,100) keyword
        if( trim(keyword) == "[ dihedraltypes ]" ) exit
    end do

    i = 1
    do
        read( 33 , * , iostat=ioerr ) (InputChars(i,k) , k=1,4) , InputIntegers(i,1) 
        if( ioerr /= 0 ) exit

        backspace(33)

        select case ( InputIntegers(i,1) )

            case( 1 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,3)

                ! conversion 
                ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
                ! kdihed0(:,1) = C0 (kJ/mol) * factor1 * imol
                ! kdihed0(:,2) = C1 (kJ/mol) * factor1 * imol
                ! kdihed0(:,3) = C2 (kJ/mol) * factor1 * imol
                ! kdihed0(:,4) = C3 (kJ/mol) * factor1 * imol
                ! kdihed0(:,5) = C4 (kJ/mol) * factor1 * imol
                ! kdihed0(:,6) = C5 (kJ/mol) * factor1 * imol
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 3 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,6)

                ! conversion 
                ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
                ! kdihed0(:,1) = C0 (kJ/mol) * factor1 * imol
                ! kdihed0(:,2) = C1 (kJ/mol) * factor1 * imol
                ! kdihed0(:,3) = C2 (kJ/mol) * factor1 * imol
                ! kdihed0(:,4) = C3 (kJ/mol) * factor1 * imol
                ! kdihed0(:,5) = C4 (kJ/mol) * factor1 * imol
                ! kdihed0(:,6) = C5 (kJ/mol) * factor1 * imol
                InputReals(i,1:6) = InputReals(i,1:6) * factor1 * imol

        end select            

        i = i + 1

    end do

    NdihedTypes = i - 1

    allocate( Dihed_Type         ( NdihedTypes     )                   )
    allocate( DihedSymbols       ( NdihedTypes , 4 )                   )
    allocate( DihedParameters    ( NdihedTypes , 6 ) , source = D_zero )

    forall(k=1:4) DihedSymbols(:NdihedTypes,k)    = InputChars(:NdihedTypes,k)
    forall(k=1:6) DihedParameters(:NdihedTypes,k) = InputReals(:NdihedTypes,k)

    Dihed_Type(:NdihedTypes) = InputIntegers(:NdihedTypes,1)

close(33)

deallocate( InputChars , InputReals , InputIntegers )
!=====================================================================================

do a = 1 , MM % N_of_species

!   Assigning to each specie the corresponding parameter ...
!   Bond parameters ...
    allocate( species(a) % kbond0( species(a) % Nbonds , 2 ) )

    do k = 1 , NbondsTypes
        do n = 1 , species(a) % Nbonds

            flag1 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) )
            flag2 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) ) .AND. & 
                    ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) ) 

            if ( flag1 .OR. flag2 ) then 
                species(a) % kbond0(n,1) = BondPairsParameters(k,1)
                species(a) % kbond0(n,2) = BondPairsParameters(k,2)
            end if
        end do
    end do

    if( allocated(SpecialBonds) ) then
        do k = 1, size(SpecialBonds)
            where( species(a) % funct_bond == SpecialBonds(k) % label ) species(a) % kbond0(:,1) = SpecialBonds(k) % kbond0(1) * factor2 * imol
            where( species(a) % funct_bond == SpecialBonds(k) % label ) species(a) % kbond0(:,2) = SpecialBonds(k) % kbond0(2) * nano_2_angs
        end do
    end if

    if( allocated(species(a) % funct_bond) ) deallocate( species(a) % funct_bond )

!   Angle parameters ...
    allocate( species(a) % kang0(species(a) % Nangs , 2 ) )

    do k = 1 , NangsTypes
        do n = 1 , species(a) % Nangs

            flag1 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,1)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,3)) )

            flag2 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,3)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,1)) )

            if( flag1 .OR. flag2 ) then
                species(a) % kang0(n,1) = AngleParameters(k,1)
                species(a) % kang0(n,2) = AngleParameters(k,2)
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

    if( allocated(species(a) % funct_angle) ) deallocate( species(a) % funct_angle )

    ! Dihedral parameters ...
    allocate( species(a) % kdihed0 ( species(a) % Ndiheds , 6 ) , source = D_zero )
    allocate( species(a) % harm    ( species(a) % Ndiheds     ) , source = I_zero )

    do n = 1 , species(a) % Ndiheds
        do k = 1 , NdihedTypes

            ! if funct = 1 (cos)
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        <== Eq. 4.61 (GMX manual 4.0)

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
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 1 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. & 
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 1 )

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! harm(:)      = n
                    species(a) % kdihed0(n,1:2) = DihedParameters(k,1:2)
                    species(a) % harm(n)        = int(DihedParameters(k,3))
                end if

            end if

            ! if funct = 3 (cos3)
            ! V = C0 + C1 * cos( phi - 180 ) + C2 * cos^2( phi - 180 ) + C3 * cos^3( phi - 180 ) + C4 * cos^4( phi - 180 ) + C5 * cos^5( phi - 180 ) 
            ! Eq. 4.62 (GMX manual 4.6.3)
            
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
                end if

            end if
        end do
    end do

end do

if( allocated(BondPairsParameters) ) deallocate( BondPairsParameters )
if( allocated(BondPairsSymbols)    ) deallocate( BondPairsSymbols    )
if( allocated(SpecialBonds)        ) deallocate( SpecialBonds        )
if( allocated(AngleParameters)     ) deallocate( AngleParameters     )
if( allocated(AngleSymbols)        ) deallocate( AngleSymbols        )
if( allocated(specialAngs)         ) deallocate( specialAngs         )
if( allocated(DihedParameters)     ) deallocate( DihedParameters     )
if( allocated(Dihed_Type)          ) deallocate( Dihed_Type          )
if( allocated(DihedSymbols)        ) deallocate( DihedSymbols        )

100 format(a18)
120  format(4a5,t22,I2,t26,6f14.4)

end subroutine top2mdflex
!
!
!
!================================================
 subroutine itp2mdflex( MM , atom , species , FF)
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
character(18)                   :: keyword , keyword_tmp
character(10)                   :: string
character(3)                    :: dummy_char , angatm1 , angatm2 , angatm3
real*8                          :: dummy_real , factor
integer                         :: i , j , k , n , a , ioerr , ilines , dummy_int , counter , Nbonds , Nangs , Ndiheds , Nbonds14 , N_of_atoms

allocate( InputChars    ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

! Reading different '.itp' species files ...
counter = 0
do a = 1 , MM % N_of_species

    string = species(a) % residue // '.itp'

    open(33, file=string, status='old',iostat=ioerr,err=101)

        101 if( ioerr > 0 ) then
            print*, string,' file not found; terminating execution' ; stop
        end if

        ! start reading the molecular structure of species(a) ...
        do
            read(33,100) keyword
            if( trim(keyword) == "[ atoms ]" ) exit
        end do
        CALL skip_lines(33,1)

        allocate( species(a) % atom ( species(a) % N_of_atoms ) )

        do i = 1 , species(a) % N_of_atoms

            read(33,*) species(a) % atom(i) % my_id ,      &
                       species(a) % atom(i) % MMSymbol ,   &
                       dummy_int ,                         &
                       species(a) % atom(i) % residue ,    &
                       species(a) % atom(i) % EHSymbol ,   &
                       dummy_int ,                         &
                       species(a) % atom(i) % MM_charge ,  &
                       species(a) % atom(i) % mass

            species(a) % atom(i) % MMSymbol   = adjustr(species(a) % atom(i) % MMSymbol)
            species(a) % atom(i) % my_species = a
            species(a) % my_species           = a
            species(a) % atom(i) % flex       = species(a) % flex

            ! this is the standard; atomic flexibity can also be defined @ ad_hoc_MM_tuning ...    
            where( atom % my_species == a ) atom % flex = species(a) % flex

            counter = counter + 1
            FF(counter) % my_id     = species(a) % atom(i) % my_id
            FF(counter) % residue   = species(a) % atom(i) % residue
            FF(counter) % EHSymbol  = species(a) % atom(i) % EHSymbol
            FF(counter) % MMSymbol  = species(a) % atom(i) % MMSymbol
            FF(counter) % MM_charge = species(a) % atom(i) % MM_charge

        end do

        N_of_atoms = species(a) % N_of_atoms

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

!==============================================================================================
        ! Bonding parameters :: reading ...
        do
            read(33,100) keyword
            if( trim(keyword) == "[ bonds ]" ) exit
        end do
        CALL skip_lines(33,1)

        i = 1
        do
            read(33,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,2 ) , InputChars(i,1)
            if( ioerr /= 0 ) exit
            i = i + 1
        end do
        backspace(33)

        Nbonds = i - 1
        species(a) % Nbonds = Nbonds

        allocate( species(a) % bonds      ( Nbonds , 2 ) )
        allocate( species(a) % funct_bond ( Nbonds     ) )

        forall(i=1:2) species(a) % bonds(:Nbonds,i) = InputIntegers(:Nbonds,i)

        species(a) % funct_bond(:Nbonds) = InputChars(:Nbonds,1)

!==============================================================================================
        ! Angle parameters :: reading ...
        do
            read(33,100) keyword
            if ( trim(keyword) == "[ angles ]" ) exit
        end do
        CALL skip_lines(33,1)

        InputIntegers = I_zero
        i = 1
        do
            read(33,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,3 ) , InputChars(i,1)
            if( ioerr /= 0 ) exit
            i = i + 1
        end do
        backspace(33)
           
        Nangs = i - 1
        species(a) % Nangs = Nangs

        allocate( species(a) % angs        ( Nangs , 3 ) )
        allocate( species(a) % funct_angle ( Nangs     ) )

        forall(i=1:3) species(a) % angs(:Nangs,i) = InputIntegers(:Nangs,i)

        species(a) % funct_angle(:Nangs) = InputChars(:Nangs,1)

!==============================================================================================
        ! expecting for dihedrals OR special-pairs ...
        do
            read(33,100,iostat=ioerr) keyword
            if ( trim(keyword) == "[ dihedrals ]" .OR. trim(keyword) == "[ pairs ]" .OR. ioerr /= 0 ) exit
        end do

        if( trim(keyword) == "[ dihedrals ]" ) then

            ! Dihedrals interactions :: reading ...
            CALL skip_lines(33,1)

            InputIntegers = I_zero
            i = 1
            do
                read(33,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,5 )
                if( ioerr /= 0 ) exit
                i = i + 1
            end do

            Ndiheds = i - 1
            species(a) % Ndiheds = Ndiheds

            allocate( species(a) % diheds  ( Ndiheds , 4 ) )
            allocate( species(a) % funct_dihed ( Ndiheds     ) )

            forall(i=1:4) species(a) % diheds(:Ndiheds,i) = InputIntegers(:Ndiheds,i)

            species(a) % funct_dihed(:Ndiheds) = InputIntegers(:Ndiheds,5)

            ! define species(a) % dihedral_type ...
            CALL define_DihedralType( species(a) , Ndiheds )

            ! expecting for special-pairs ...
            do
                read(33,100,iostat=ioerr) keyword
                if( trim(keyword) == "[ pairs ]" .OR. ioerr /= 0 ) exit
            end do
            CALL skip_lines(33,1)

        end if

        if( trim(keyword) == "[ pairs ]" ) then

            ! Special-pairs interactions :: reading ...
            CALL skip_lines(33,1)

            InputIntegers = I_zero
            i = 1
            do
                read(33,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,2 ) , InputReals(i,1)
                if( ioerr /= 0 ) exit
                i = i + 1
            end do

            Nbonds14 = i - 1
            species(a) % Nbonds14 = Nbonds14

            allocate( species(a) % bonds14 ( Nbonds14 , 2 ) )
            allocate( species(a) % fact14  ( Nbonds14     ) )

            forall(i=1:2) species(a) % bonds14(:Nbonds14,i) = InputIntegers(:Nbonds14,i)

            species(a) % fact14(:Nbonds14) = InputReals(:Nbonds14,1)

        end if

!==============================================================================================

    close(33)

end do

FF % residue  = adjustl(FF % residue)
FF % Symbol   = adjustl(FF % Symbol)
FF % MMSymbol = adjustl(FF % MMSymbol)

deallocate( InputChars , InputIntegers )

100 format(a18)

end subroutine itp2mdflex
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

        case( 1 )

            a % dihedral_type(i) = "cos"

        case( 3 )
            
            a % dihedral_type(i) = "cos3"

    end select

    a % dihedral_type(i) = adjustl( a % dihedral_type(i) )

end do

end subroutine define_DihedralType
!
!
!
!=========================================
subroutine skip_lines( file_unit , lines )
!=========================================
implicit none
integer , intent(in) :: file_unit
integer , intent(in) :: lines

!local variables ...
integer :: i , ioerr

do i = 1 , lines
    read(file_unit,*,iostat=ioerr) 
    if(ioerr < 0) exit
end do

end subroutine skip_lines
!
!
!
end  module gmx2mdflex
