! Convert gmx data to mdflex program :: verify gmx format in IO_FORMATS 

module namd2mdflex

use iso_fortran_env
use MM_input               , only : MM_input_format
use constants_m
use for_force
use MPI_definitions_m      , only : master
use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles, DefinePairs, DefineMorse, debug_MM
use MM_tuning_routines     , only : SpecialBonds, SpecialAngs
use NonBondPairs           , only : Identify_NonBondPairs
use Babel_routines_m       , only : TO_UPPER_CASE
use gmx2mdflex             , only : SpecialPairs , SpecialPairs14

public :: prm2mdflex, psf2mdflex, convert_NAMD_velocities, SpecialPairs, SpecialPairs14

private
 
    ! module variables ...
    character(4)      , allocatable  , save  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:)
    real*8            , allocatable  , save  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:)
    integer                                  :: GAFF_order(4) = [3,2,1,4] 

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
integer                         :: i , j , k , a , n , ioerr , counter , N_of_atoms
integer                         :: Nbonds , Nangs , Ndiheds , Nimpropers 

allocate( InputChars    ( 20000 , 10 )                   )
allocate( InputReals    ( 20000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 30000 , 10 ) , source = I_zero )

! Reading different '.itp' species files ...
counter = 0

do a = 1 , MM % N_of_species

    string = species(a) % residue // '.psf'

    open(33, file=string, status='old',iostat=ioerr,err=101)

        101 if( ioerr > 0 ) then
            print*, string,' file not found; terminating execution' ; stop
        end if

        If( master ) write(*,'(/2a9)',advance='no') "Reading ", string

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
          if( verify( "!NBOND" , keyword ) == 0 ) exit
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
        
        species(a) % bond_type(:) = "harm" 
        
        rewind 33

!==============================================================================================
        ! Angle parameters :: reading ...
        do
          read(33,100) keyword
          if( verify( "!NTHETA" , keyword ) == 0 ) exit
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
        species(a) % angle_type (:Nangs) = "harm"

        rewind 33

!==============================================================================================
        ! Dihedral parameters :: reading ...

        !=========================================================================
        !reading normal dihedrals ...
        do
          read(33,100,iostat=ioerr) keyword
          if( verify( "!NPHI" , keyword ) == 0 ) exit
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
          if( verify( "!NIMPHI" , keyword ) == 0 ) exit
        end do
        backspace(33)
        read(33,*) Nimpropers
  
        ! reading full lines ...
        do k = 1 , ceiling(Nimpropers/two)-1
            select case ( MM_input_format )
                case( "GAFF" )
                read(33 , * , iostat=ioerr )  ( ( InputIntegers(Ndiheds+(k-1)*2+n,GAFF_order(j)) , j=1,4 ) , n=1,2 )
                case default
                read(33 , * , iostat=ioerr )  ( ( InputIntegers(Ndiheds+(k-1)*2+n,           j ) , j=1,4 ) , n=1,2 )
            end select
        end do 

        ! reading incomplete lines ...
        select case ( MM_input_format )
            case( "GAFF" )
            read(33 , * , iostat=ioerr )  &
            ( ( InputIntegers(Ndiheds+(k-1)*2+n,GAFF_order(j)) , j=1,4 ) , n=1,merge(2,mod(Nimpropers,2),mod(Nimpropers,2)==0) ) 
            case default
            read(33 , * , iostat=ioerr )  &
            ( ( InputIntegers(Ndiheds+(k-1)*2+n,           j ) , j=1,4 ) , n=1,merge(2,mod(Nimpropers,2),mod(Nimpropers,2)==0) ) 
        end select

        !=========================================================================

        species(a)% NTorsions  = Ndiheds
        species(a)% NImpropers = NImpropers

        species(a)% Ndiheds = Ndiheds + Nimpropers

        allocate( species(a)% diheds      ( species(a)% Ndiheds , 4 ) )
        allocate( species(a)% funct_dihed ( species(a)% Ndiheds     ) )

        ! store normal dihedrals ...
        forall(i=1:4) species(a)% diheds( :Ndiheds , i ) = InputIntegers( :Ndiheds , i )
        species(a)% funct_dihed( :Ndiheds ) = 9

        ! store improper dihedrals ...
        forall(i=1:4) species(a)% diheds( Ndiheds+1:Ndiheds+Nimpropers , i ) = InputIntegers( Ndiheds+1:Ndiheds+Nimpropers , i )

        select case ( MM_input_format )
            case( "GAFF" , "NAMD" )   
                species(a)% funct_dihed( Ndiheds+1:Ndiheds+Nimpropers ) = 9
            case( "CHMM" )
                species(a)% funct_dihed( Ndiheds+1:Ndiheds+Nimpropers ) = 2
        end select

        ! define species(a) % dihedral_type ...
        CALL define_DihedralType( species(a) , species(a)% Ndiheds )

        rewind 33

!==============================================================================================

        If( master ) then
            TorF = Checking_Topology( species(a)%bonds , species(a)%angs , species(a)%diheds(:Ndiheds,:) )
            If( TorF ) then
                CALL system("sed '11i >>> error detected in Topology , check log.trunk/Topology.test.log <<<' warning.signal |cat")
                stop
            End If
        End If

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

    If( master ) write(*,'(a9)') " << done "

end do

FF % residue  = adjustl(FF % residue)
FF % Symbol   = adjustl(FF % Symbol)
FF % MMSymbol = adjustl(FF % MMSymbol)

deallocate( InputChars , InputIntegers )

100 format(a18)

end subroutine psf2mdflex
!
!
!
!==========================================
 subroutine prm2mdflex( MM , species , FF )
!==========================================
implicit none 
type(MM_molecular)                  , intent(inout) :: species(:)
type(MM_system)                     , intent(inout) :: MM
type(MM_atomic)     , allocatable   , intent(inout) :: FF(:)
 
! local variables ...
character(4)    , allocatable   :: InputChars(:,:) , Input2Chars(:,:)
character(4)    , allocatable   :: funct_bond(:) , funct_angle(:)
real*8          , allocatable   :: InputReals(:,:) , Input2Reals(:,:)
integer         , allocatable   :: InputIntegers(:,:)
integer         , allocatable   :: Dihed_Type(:) , Bond_Type(:) , Angle_Type(:)
integer                         :: a , n , i , j , k , l , l1 , j1 , dummy_int , ioerr , N_of_AtomTypes 
integer                         :: NbondsTypes , NangsTypes , NdihedTypes , NTorsionTypes , NImproperTypes, NBondParms, SpecialNBParms
real*8                          :: SCEE , SCNB
character(3)                    :: dummy_char
character(18)                   :: keyword
character(200)                  :: line
logical                         :: flag1 , flag2 , flag3 , flag4 , flag5 , flag6 

allocate( InputChars    ( 10000 , 10 )                   )
allocate( Input2Chars   ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( Input2Reals   ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

! NAMD FF definitions ... 
forcefield = 2   ! <== 1 = Born-Mayer (not implemented); 2 = Lennard-Jones (OK)
  
open(33, file='input.prm', status='old', iostat=ioerr, err=10)

!   file error msg ...
    10 if( ioerr > 0 ) stop '"input.prm" file not found; terminating execution'

!=====================================================================================
!  reading  FF DEFAULTS ...
    do
        read(33,100) keyword
        if( verify( "FF-sets" , keyword ) == 0 ) exit
    end do

    read(33,*) dummy_char
    read(33,*) MM% CombinationRule, SCNB, SCEE 

    MM % fudgeQQ = 1.d0 / SCEE
    MM % fudgeLJ = 1.d0 / SCNB

! checking input data ...
If( (MM_input_format == "GAFF") .AND. (SCNB/=1.0) ) stop " >>> WARNING: supposedely SCNB = 1 for GAFF MM_input_format "

!=====================================================================================
!  reading  ATOMS ...

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
        if( trim(InputChars(i,1)) == "BOND" ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop
        read(line,*,iostat=ioerr) InputChars(i,1) , dummy_int, InputChars(i,2) 

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
        if( trim(InputChars(i,1)) == "ANGL" ) exit
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

    forall(i=1:2) BondPairsSymbols(:NbondsTypes,i)    = adjustl(InputChars(:NbondsTypes,i))
    forall(i=1:2) BondPairsParameters(:NbondsTypes,i) = InputReals(:NbondsTypes,i)
  
   funct_bond( :NbondsTypes ) = "harm" 
   do i = 1 , NbondsTypes 
        BondPairsParameters(i,1) = InputReals(i,1) * TWO * factor1 * imol * cal_2_J 
        BondPairsParameters(i,2) = InputReals(i,2) 
   end do

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
        if( trim(InputChars(i,1)) == "DIHE" ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop2
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,3) , (InputReals(i,j) , j=1,4 )

        i = i + 1
    end do read_loop2
    backspace(33)

    NangsTypes = i - 1

    allocate( AngleSymbols    ( NangsTypes , 3 ) )
    allocate( AngleParameters ( NangsTypes , 4 ) , source = D_zero ) 
    allocate( Angle_Type      ( NangsTypes     ) )
    allocate( funct_angle     ( NangsTypes     ) )
    
    forall(i=1:3) AngleSymbols(:NangsTypes,i)     = InputChars(:NangsTypes,i)
    forall(i=1:4) AngleParameters(:NangsTypes,i)  = InputReals(:NangsTypes,i)
   
    do i = 1 , NangsTypes
        ! Harmonic potential ...
        ! factor1 = 1.0d26      <== Factor used to correct the units 
        if( AngleParameters(i,3) == D_zero ) then 
           Angle_Type(i) = 1 
           funct_angle(i) = "harm"
           AngleParameters(i,1) = InputReals(i,1) * TWO * factor1 * imol * cal_2_J
           AngleParameters(i,2) = InputReals(i,2) * deg_2_rad
        ! Urey-Bradley potential ...
        else 
           Angle_Type(i) = 3 
           funct_angle(i) = "urba"
           AngleParameters(i,1) = InputReals(i,1) * TWO * factor1 * imol * cal_2_J
           AngleParameters(i,2) = InputReals(i,2) * deg_2_rad
           AngleParameters(i,3) = InputReals(i,3) * TWO * factor1 * imol * cal_2_J 
           AngleParameters(i,4) = InputReals(i,4) 
        end if        
    end do

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
        if( trim(InputChars(i,1)) == "IMPR" ) exit
        if( trim(InputChars(i,1)) == "CMAP" ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop3
        read(line,*,iostat=ioerr) (InputChars(i,k) , k=1,4) , (InputReals(i,j) , j=1,3)

        i = i + 1
    end do read_loop3
    backspace(33)

    NdihedTypes   = i - 1
    NTorsionTypes = NdihedTypes

!=====================================================================================
!  reads IMPROPERS ... 
    do
        read(33,100) keyword
        if( trim(keyword) == "IMPROPERS" ) exit
    end do

    i = 1
    read_loop4: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop4
        read(line,*,iostat=ioerr) Input2Chars(i,1), Input2Chars(i,2)
        if( index(Input2Chars(i,1),"!") /= 0 ) cycle read_loop4
        if( trim(Input2Chars(i,1)) == "CMAP" ) exit
        if( trim(Input2Chars(i,1)) == "NONB" ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop4

        select case (MM_input_format) 
            case( "GAFF" ) 
            read(line,*,iostat=ioerr) (InputChars(Ndihedtypes+i,GAFF_order(k)) , k=1,4) , (InputReals(Ndihedtypes+i,j) , j=1,3)
            case default
            read(line,*,iostat=ioerr) (InputChars(Ndihedtypes+i,k            ) , k=1,4) , (InputReals(Ndihedtypes+i,j) , j=1,3)
        end select 

        i = i + 1
    end do read_loop4
    backspace(33)

    NImproperTypes = i - 1

    allocate( Dihed_Type         ( NdihedTypes+NImproperTypes     )                   )
    allocate( DihedSymbols       ( NdihedTypes+NImproperTypes , 4 )                   )
    allocate( DihedParameters    ( NdihedTypes+NImproperTypes , 3 ) , source = D_zero )

    ! Torsion ...
    Dihed_Type(:NdihedTypes) = 9 
    ! Improper
    select case ( MM_input_format )
        case( "GAFF" , "NAMD" )   
            Dihed_Type(NdihedTypes+1:NdihedTypes+NImproperTypes) = 9         
        case( "CHMM" )
            Dihed_Type(NdihedTypes+1:NdihedTypes+NImproperTypes) = 2        
    end select

    NdihedTypes = NdihedTypes + NImproperTypes

    forall(k=1:4) DihedSymbols(:NdihedTypes,k)    = InputChars(:NdihedTypes,k)
    forall(k=1:3) DihedParameters(:NdihedTypes,k) = InputReals(:NdihedTypes,k)

    do i = 1 , NdihedTypes
        select case( Dihed_Type(i) )
           case( 9 )
                ! V = k[1 + cos(n.phi - theta)] (charmm; multiple)
                ! factor1 = 1.0d26      <== Factor used to correct units 
                ! kdihed0(:,1) = phi_s   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(phi) ==> force constant (kcal.mol⁻¹) * factor1 * imol * cal_2_J
                ! kdihed0(:,3) = n       ==> multiplicity 
                DihedParameters(i,1) = InputReals(i,3) * deg_2_rad
                DihedParameters(i,2) = InputReals(i,1) * factor1 * imol * cal_2_J
                DihedParameters(i,3) = InputReals(i,2) 

           case( 2 )
                ! V = 1/2.k[cos(phi) - cos(phi0)]²
                ! factor1 = 1.0d26      <== Factor used to correct units 
                ! kdihed0(:,1) = xi_0   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(xi) ==> force constant (kcal.mol⁻¹.rad⁻²) * factor1 * imol * cal_2_J
                DihedParameters(i,1) = InputReals(i,3) * deg_2_rad
                DihedParameters(i,2) = InputReals(i,1) * TWO * factor1 * imol * cal_2_J

        end select 
    end do

!=====================================================================================
!  NonBonding parameters :: reading ...
    do
        read(33,100) keyword
        if( trim(keyword(1:9)) == "NONBONDED" ) exit
    end do
    read(33,100)
    
    InputReals = D_zero
    i = 1
    read_loop5: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop5
        read(line,*,iostat=ioerr) InputChars(i,1)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loop5
        if( trim(InputChars(i,1)) == "SPEC" ) exit
        if( trim(InputChars(i,1)) == "HBON" ) exit
        if( trim(InputChars(i,1)) == "END " ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop5
        read(line,*, iostat=ioerr) InputChars(i,1) , (InputReals(i,j) , j=1,6)
        
        i = i + 1
    end do read_loop5
    backspace(33)

    NBondParms = i - 1

    do i = 1 , NBondParms
        where( ( FF % MMSymbol == InputChars(i,1) ) .OR. ( ( adjustR(FF % MMSymbol(1:2)) // "*" )  == InputChars(i,1) ) )
            FF % sig = InputReals(i,3)
            FF % eps = abs(InputReals(i,2))
        end where
    end do
    backspace(33)
    FF % eps14 = FF % eps
    FF % sig14 = FF % sig

    do i = 1 , NBondParms
        if( InputReals(i,5) /= D_zero ) then
            where( ( FF % MMSymbol == InputChars(i,1) ) .OR. ( ( adjustR(FF % MMSymbol(1:2)) // "*" ) == InputChars(i,1) ) )
                FF % sig14 = InputReals(i,6)
                FF % eps14 = abs(InputReals(i,5)) 
            end where 
        end if
    end do

    ! conversion 
    ! factor1 = 1.0d26  <== Factor used to correct units 
    ! GAFF  vs  GMX  LJ parameters:
    ! -> epsilon_GAFF = epsilon_GMX / (cal_2_J * 2) 
    ! -> sigma_GAFF = (sigma_GMX*10/2 ) * 2^(1/6)

    FF % eps   = sqrt( FF % eps   * factor1 * imol * cal_2_J )
    FF % eps14 = sqrt( FF % eps14 * factor1 * imol * cal_2_J )
    FF % sig   = ( FF % sig   * TWO ) / (2**(1.d0/6.d0)) ! amber_LJ
    FF % sig14 = ( FF % sig14 * TWO ) / (2**(1.d0/6.d0)) ! amber_LJ

    select case( MM % CombinationRule )

        case (2)

            FF % sig   = FF % sig   / TWO
            FF % sig14 = FF % sig14 / TWO

        case (3)

            FF % sig   = sqrt( FF % sig   )
            FF % sig14 = sqrt( FF % sig14 )

    end select

!=====================================================================================
!  SPECIALNonBonding parameters :: reading ...
    do
        read(33,100,iostat=ioerr) keyword
        if(ioerr == iostat_end) goto 99     ! <== end of file 
        if( trim(keyword(1:7)) == "SPECIAL" ) exit
    end do
    read(33,100)

    InputReals = D_zero
    i = 1
    read_loopS: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loopS
        read(line,*,iostat=ioerr) InputChars(i,1)
        if( index(InputChars(i,1),"!") /= 0 ) cycle read_loopS
        if( trim(InputChars(i,1)) == "HBON" ) exit
        if( trim(InputChars(i,1)) == "END " ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loopS
        read(line,*, iostat=ioerr) (InputChars(i,j) , j=1,2) , (InputReals(i,j) , j=1,6)
        
        i = i + 1
    end do read_loopS
    backspace(33)

    SpecialNBParms = i - 1

    If( SpecialNBParms /= 0 ) then

        allocate( SpecialPairs ( SpecialNBParms ) )

        forall(i=1:2) SpecialPairs(:SpecialNBParms) % MMSymbols(i) = InputChars(:SpecialNBParms,i)

        SpecialPairs(:SpecialNBParms) % Parms(1) = InputReals(:SpecialNBParms,3)
        SpecialPairs(:SpecialNBParms) % Parms(2) = abs(InputReals(:SpecialNBParms,2))

        ! conversion 
        ! factor1 = 1.0d26  <== Factor used to correct units 
        ! GAFF  vs  GMX  LJ parameters:
        ! -> epsilon_GAFF = epsilon_GMX / (cal_2_J * 2) 
        ! -> sigma_GAFF = (sigma_GMX*10/2 ) * 2^(1/6)

        SpecialPairs(:SpecialNBParms) % Parms(1) = SpecialPairs(:SpecialNBParms) % Parms(1) * 2**(5.d0/6.d0) 
        SpecialPairs(:SpecialNBParms) % Parms(2) = SpecialPairs(:SpecialNBParms) % Parms(2) * factor1 * imol * cal_2_J 

        allocate( SpecialPairs14 ( SpecialNBParms ) )

        forall(i=1:2) SpecialPairs14(:SpecialNBParms) % MMSymbols(i) = InputChars(:SpecialNBParms,i)

        SpecialPairs14(:SpecialNBParms) % Parms(1) = InputReals(:SpecialNBParms,6)
        SpecialPairs14(:SpecialNBParms) % Parms(2) = abs(InputReals(:SpecialNBParms,5))

        ! conversion 
        SpecialPairs14(:SpecialNBParms) % Parms(1) = SpecialPairs14(:SpecialNBParms) % Parms(1) * 2**(5.d0/6.d0) 
        SpecialPairs14(:SpecialNBParms) % Parms(2) = SpecialPairs14(:SpecialNBParms) % Parms(2) * factor1 * imol * cal_2_J 

    endIf

!=====================================================================================
!
99 close(33)

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
                species(a) % kbond0(n,1:3) = BondPairsParameters(k,1:3)
            end if
        end do
    end do

    !=============================================================================
    ! Angle parameters ...
    allocate( species(a) % kang0(species(a) % Nangs , 4 ) , source = D_zero )

    do n = 1 , species(a) % Nangs
        do k = 1 , NangsTypes

            flag1 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,1)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,3)) )

            flag2 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,3)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) .AND. &
                    ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,1)) )

            if ( flag1 .OR. flag2 ) then 
                    species(a) % kang0(n,1:4) = AngleParameters(k,1:4) 
                if( species(a) % kang0(n,3) /= D_zero ) then 
                    species(a) % angle_type(n)  = "urba" 
                    species(a) % funct_angle(n) = "3" 
                end if
            end if

        end do 
    end do 

    !=============================================================================
    ! Dihedral parameters ...
    allocate( species(a)% kdihed0 ( species(a) % Ndiheds , 15 ) , source = D_zero )

    read_loop0: do n = 1 , species(a)% NTorsions ! <== Mind that some impropers may be included here ...
        ! control variables to multiple dihs ...
        j = 0 ; l = 0 ; l1 = 0 ; j1 = 0 

        read_loop7: do k = 1 , NTorsionTypes 

            ! if funct = 9 (chrm) (multiple) 
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple)       
            ! Eq. 4.60 (GMX manual 5.0.5)

            if( species(a) % funct_dihed(n) == 9 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )


                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 9 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' )                                                             .AND. & 
                        ( adjustl(DihedSymbols(k,4)) == 'X' )                                                             .AND. &
                        ( Dihed_Type(k) == 9 )

                flag5 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        (                                                             'X' == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag6 = (                                                             'X' == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                if( flag1 .OR. flag2 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    if( j1 > 0 ) species(a) % kdihed0(n,:) = D_zero 
                    if( j1 > 0 ) j1 = 0 
                    species(a) % kdihed0(n,3*(j+j1)+1) = DihedParameters(k,1) 
                    species(a) % kdihed0(n,3*(j+j1)+2) = DihedParameters(k,2) 
                    species(a) % kdihed0(n,3*(j+j1)+3) = DihedParameters(k,3) 
                    j = j + 1
                    cycle read_loop7
                end if

                if( flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n 
                    if( j > 0 ) cycle read_loop7
                    species(a) % kdihed0(n,3*j1+1) = DihedParameters(k,1)
                    species(a) % kdihed0(n,3*j1+2) = DihedParameters(k,2)
                    species(a) % kdihed0(n,3*j1+3) = DihedParameters(k,3)
                    j1= j1+ 1 
                    cycle read_loop7
                end if

            end if

        end do read_loop7
    end do read_loop0

    !=============================================================================
    ! Improper parameters ...

    read_loop6: do n = species(a)% NTorsions + 1 , species(a)% Ndiheds

        ! control variables to multiple dihs ...
        j = 0 ; l = 0 ; l1 = 0 ; j1 = 0 

        read_loop8: do k = NTorsionTypes + 1 , NdihedTypes   

            ! if funct = 9 (chrm) (multiple) 
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple)       
            ! Eq. 4.60 (GMX manual 5.0.5)

            if( species(a) % funct_dihed(n) == 9 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( adjustl(DihedSymbols(k,2)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,3)) == 'X' ) .AND. & 
                        ( Dihed_Type(k) == 9 ) 
                
                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( adjustl(DihedSymbols(k,2)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,3)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag5 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        (                                                             'X' == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                if( flag1 .OR. flag2 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    if( j1 > 0 ) species(a) % kdihed0(n,:) = D_zero 
                    if( j1 > 0 ) j1 = 0 
                    species(a) % kdihed0(n,3*(j+j1)+1) = DihedParameters(k,1) 
                    species(a) % kdihed0(n,3*(j+j1)+2) = DihedParameters(k,2) 
                    species(a) % kdihed0(n,3*(j+j1)+3) = DihedParameters(k,3) 
                    j = j + 1
                    cycle read_loop8
                end if

                if( flag3 .OR. flag4 .OR. flag5 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n 
                    if( j > 0 ) cycle read_loop8
                    species(a) % kdihed0(n,3*j1+1) = DihedParameters(k,1)
                    species(a) % kdihed0(n,3*j1+2) = DihedParameters(k,2)
                    species(a) % kdihed0(n,3*j1+3) = DihedParameters(k,3)
                    j1= j1+ 1 
                    cycle read_loop8
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

                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( adjustl(DihedSymbols(k,2)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,3)) == 'X' ) .AND. & 
                        ( Dihed_Type(k) == 2 ) 

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( adjustl(DihedSymbols(k,2)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,3)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 2 )

                flag5 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 2 )

                flag6 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 2 )

                if( flag1 .OR. flag2 ) then
                    ! kdihed0(:,1) = xi_0 (deg)
                    ! kdihed0(:,2) = k_xi [ kJ/(mol.rad²) ] 
                    if( l1 > 0 ) species(a) % kdihed0(n,:) = D_zero 
                    if( l1 > 0 ) l1 = 0 
                    species(a) % kdihed0(n,2*l+1) = DihedParameters(k,1)
                    species(a) % kdihed0(n,2*l+2) = DihedParameters(k,2) 
                    l = l + 1
                    cycle read_loop8
                end if

             if( flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    ! kdihed0(:,1) = xi_0 (deg)
                    ! kdihed0(:,2) = k_xi [ kJ/(mol.rad²) ] 
                    if( l > 1 ) cycle read_loop8
                    species(a) % kdihed0(n,2*l1+1) = DihedParameters(k,1)
                    species(a) % kdihed0(n,2*l1+2) = DihedParameters(k,2)
                    l1 = l1 + 1
                    cycle read_loop8
                end if

            end if

        end do read_loop8
    end do read_loop6
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
120 format(4a5,t22,I2,t26,6f14.4)
 
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
!
!================================================
 subroutine convert_NAMD_velocities( N_of_atoms )
!================================================
implicit none
integer , intent(in) :: N_of_atoms

! local variables ...
integer                        :: i
character(len=1)               :: choice , dumb
real*8           , allocatable :: vel(:,:)
logical                        :: exist

allocate( vel(N_of_atoms,3) )

write(*,'(/a)',advance='no') ">>> Want to convert velocity_MM.pdb to DynEMol format/units ? (y/n)"
read (*,'(a)') choice

if( choice == 'y' ) then

    inquire(file="velocity_MM.pdb", EXIST=exist)
    if (exist) then

        ! read velocity file in pdb format, and convert units ...
        open(unit=33 , file='velocity_MM.pdb' , status='old' , action='read')
        read(33,*) dumb
        do i = 1 , N_of_atoms
             read(33,15) vel(i,1) , vel(i,2) , vel(i,3)
             ! convert velocities from = Angs/ps  to = m/s ...
             vel(i,:) = vel(i,:) * 1.d2  
        end do
        close(33)

        ! write velocity file in xyz format ...
        open(unit=33 , file='velocity_MM.inpt' , status='unknown')
        do i = 1 , N_of_atoms
             write(33,*) vel(i,1) , vel(i,2) , vel(i,3)
        end do
        close(33)
        write(*,'(/a)') ' >>> velocity_MM.inpt generated <<<'

        write(*,'(/a)',advance='no') ">>> Want to continue (n) ? (y/n)"
        read (*,'(a)') choice

        If( choice /= "y" ) stop

    else
        STOP ' >> file velocity_MM.pdb was not found ! << '
    endif

endif

15 Format(t31,f8.3,t39,f8.3,t47,f8.3)

end subroutine convert_NAMD_velocities
!
!
!
end  module namd2mdflex
