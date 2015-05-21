! Convert gmx data to mdflex program :: verify gmx format in IO FORMATS 
module gmx2mdflex

use constants_m
use for_force
use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles, DefinePairs, DefineMorse, debug_MM
use MM_tuning_routines     , only : SpecialBonds, SpecialAngs
use NonBondPairs           , only : Identify_NonBondPairs

private
 
public :: top2mdflex, itp2mdflex, SpecialPairs, SpecialMorse 

    ! module variables ...
    character(3)     , allocatable   , save  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:)
    real*8           , allocatable   , save  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:)
    type(DefinePairs) , allocatable :: SpecialPairs(:)
    type(DefineMorse) , allocatable :: SpecialMorse(:)

contains
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
character(18)                   :: keyword 
character(10)                   :: string
character(200)                  :: line 
integer                         :: i1 , i2 , i3 , sp , nr
integer                         :: i , j , k , a , ioerr , dummy_int , counter , Nbonds , Nangs , Ndiheds , Nbonds14 , N_of_atoms

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

        write(*,'(/2a9)',advance='no') "Reading ", string

        ! start reading the molecular structure of species(a) ...
        do
            read(33,100) keyword
            if( trim(keyword) == "[ atoms ]" ) exit
        end do

        allocate( species(a) % atom ( species(a) % N_of_atoms ) )

        i = 1
        read_loop1: do
            read(33, '(A)', iostat=ioerr) line
            if ( ioerr /= 0 ) exit read_loop1
            read(line,*,iostat=ioerr) InputChars(i,1) 
            if( index(InputChars(i,1),";") /= 0 ) cycle read_loop1
            if( trim(InputChars(i,1)) == "[  "  ) exit
            if( ioerr > 0  ) exit
            if( ioerr /=  0 ) cycle read_loop1
            read(line,*,iostat=ioerr) species(a) % atom(i) % my_id ,      &
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
            FF(counter) % my_species = a
            FF(counter) % my_id      = species(a) % atom(i) % my_id
            FF(counter) % residue    = species(a) % atom(i) % residue
            FF(counter) % EHSymbol   = species(a) % atom(i) % EHSymbol
            FF(counter) % MMSymbol   = species(a) % atom(i) % MMSymbol
            FF(counter) % MM_charge  = species(a) % atom(i) % MM_charge 

            i = i + 1
 
        end do read_loop1
        backspace(33)
 
        ! DANDO ERRO NO COMENTÁRIO ABAIXO
        !If( size(species(a)%atom) /= count(atom(:)%residue == species(a)%atom(1)%residue) )  &
        !stop "residue size of this species differs from atom%residue; check tuning.f"

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

        i = 1
        read_loop2: do
            read(33, '(A)', iostat=ioerr) line
            if ( ioerr /= 0 ) exit read_loop2
            read(line,*,iostat=ioerr) InputChars(i,1)
            if( index(InputChars(i,1),";") /= 0 ) cycle read_loop2
            if( trim(InputChars(i,1)) == "[  "  ) exit 
            if( ioerr > 0  ) exit
            if( ioerr /= 0 ) cycle read_loop2
            read(line,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,2 ) , InputChars(i,1)

            i = i + 1
        end do read_loop2
        backspace(33)

        Nbonds = i - 1
        species(a) % Nbonds = Nbonds

        allocate( species(a) % bonds      ( Nbonds , 2 ) )
        allocate( species(a) % funct_bond ( Nbonds     ) )
        allocate( species(a) % bond_type  ( Nbonds     ) )

        forall(i=1:2) species(a) % bonds(:Nbonds,i) = InputIntegers(:Nbonds,i)

        species(a) % funct_bond(:Nbonds) = adjustl( (InputChars(:Nbonds,1)) )
        
        do i = 1 , Nbonds 
           select case ( species(a) % funct_bond(i) )
            case( "1" )  
                species(a) % bond_type(i) = "harm" 
            case( "3" )
                species(a) % bond_type(i) = "Mors"     
           end select 
        end do

!==============================================================================================
        ! Angle parameters :: reading ...
        do
            read(33,100) keyword
            if ( trim(keyword) == "[ angles ]" ) exit
        end do

        InputIntegers = I_zero
        i = 1
        read_loop3: do
            read(33, '(A)', iostat=ioerr) line
            if ( ioerr /= 0 ) exit read_loop3
            read(line,*,iostat=ioerr) InputChars(i,1)
            if( index(InputChars(i,1),";") /= 0 ) cycle read_loop3
            if( trim(InputChars(i,1)) == "[  "  ) exit
            if( ioerr > 0  ) exit
            if( ioerr /= 0 ) cycle read_loop3
            read(line,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,3 ) , InputChars(i,1)

            i = i + 1
        end do read_loop3
        backspace(33)
           
        Nangs = i - 1
        species(a) % Nangs = Nangs

        allocate( species(a) % angs        ( Nangs , 3 ) )
        allocate( species(a) % funct_angle ( Nangs     ) )

        forall(i=1:3) species(a) % angs(:Nangs,i) = InputIntegers(:Nangs,i)

        species(a) % funct_angle(:Nangs) = InputChars(:Nangs,1)

!==============================================================================================
        ! Dihedral parameters :: reading ...
        do
            read(33,100,iostat=ioerr) keyword
            if ( trim(keyword) == "[ dihedrals ]" .OR. ioerr /= 0 ) exit
        end do

        if( trim(keyword) == "[ dihedrals ]" ) then

            InputIntegers = I_zero
            i = 1
            read_loop4: do
                read(33, '(A)', iostat=ioerr) line
                if ( ioerr /= 0 ) exit read_loop4
                read(line,*,iostat=ioerr) InputChars(i,1)
                if( index(InputChars(i,1),";") /= 0 ) cycle read_loop4
                if( trim(InputChars(i,1)) == "[  "  ) exit
                if( ioerr > 0  ) exit
                if( ioerr /= 0 ) cycle read_loop4
                read(line,*, iostat=ioerr) ( InputIntegers(i,j) , j=1,5 )

                i = i + 1
            end do read_loop4
            backspace(33)

            Ndiheds = i - 1
            species(a) % Ndiheds = Ndiheds

            allocate( species(a) % diheds  ( Ndiheds , 4 ) )
            allocate( species(a) % funct_dihed ( Ndiheds     ) )

            forall(i=1:4) species(a) % diheds(:Ndiheds,i) = InputIntegers(:Ndiheds,i)

            species(a) % funct_dihed(:Ndiheds) = InputIntegers(:Ndiheds,5)

            ! define species(a) % dihedral_type ...
            CALL define_DihedralType( species(a) , Ndiheds )

        end if

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

end subroutine itp2mdflex
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
character(3)    , allocatable   :: InputChars(:,:) , Input2Chars(:,:)
character(4)    , allocatable   :: funct_bond(:)
real*8          , allocatable   :: InputReals(:,:) , Input2Reals(:,:)
integer         , allocatable   :: InputIntegers(:,:)
integer         , allocatable   :: Dihed_Type(:) , Bond_Type(:)
real*8                          :: fudgeLJ , fudgeQQ
integer                         :: a , n , i , j , k , ioerr , dummy_int , N_of_AtomTypes 
integer                         :: NbondsTypes , NangsTypes , NdihedTypes , NBondParms, NMorseParms
character(3)                    :: dummy_char
character(18)                   :: keyword
character(200)                  :: line
logical                         :: flag1 , flag2 , flag3 , flag4

allocate( InputChars    ( 10000 , 10 )                   )
allocate( Input2Chars   ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( Input2Reals   ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

forcefield = 2
  
open(33, file='topol.top', status='old', iostat=ioerr, err=10)

!   file error msg ...
    10 if( ioerr > 0 ) stop '"topol.top" file not found; terminating execution'

    do
      read(33,100) keyword
      if( trim(keyword) == "[ defaults ]" ) exit
    end do

    i=1
    read_loop: do
        read(33,*,iostat=ioerr) dummy_char, MM % CombinationRule, dummy_char, fudgeLJ, fudgeQQ 
        if ( index( dummy_char, ";") /= 0 ) cycle read_loop
        if( ioerr /= 0 ) exit
        i = i + 1
    end do read_loop

    backspace(33)
 
!=====================================================================================
!   reading the number of [ atomtypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ atomtypes ]" ) exit
    end do

    i=1
    read_loop1: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop1
        read(line,*,iostat=ioerr) InputChars(i,1)
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop1
        if( trim(InputChars(i,1)) == "[  "  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop1
        read(line,*,iostat=ioerr) InputChars(i,1) , (InputReals(i,j) , j=1,2) , InputChars(i,2) , (InputReals(i,j) , j=3,4)

        i = i + 1
    end do read_loop1
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
    do
        read(33,100) keyword
        if( trim(keyword) == "[ nonbond_params ]" ) exit
    end do

    i = 1
    k = 1
    read_loop5: do
        read(33,*, iostat=ioerr) (InputChars(i,j) , j=1,2) , InputIntegers(i,1) , (InputReals(i,j) , j=1,2)
        if( ioerr /= 0 ) exit 
        if( InputIntegers(i,1) == 3 ) then 
            backspace(33)
            read(33,*) (Input2Chars(k,j) , j=1,2) , dummy_int, (Input2Reals(k,j) , j=1,3) 
            k = k + 1 
            cycle read_loop5
        end if
        i = i + 1
    end do read_loop5
    backspace(33)

    NBondParms = i - 1

    allocate( SpecialPairs ( NbondParms ) )

    forall(i=1:2) SpecialPairs(:NBondParms) % MMSymbols(i) = InputChars(:NbondParms,i)

    SpecialPairs(:NBondParms) % Parms(2) = InputReals(:NbondParms,2)
    SpecialPairs(:NBondParms) % Parms(1) = InputReals(:NbondParms,1) 

    ! conversion 
    ! factor1 = 1.0d26      <== Factor used to correct the units readed from Gromacs
    SpecialPairs(:NBondParms) % Parms(2) = sqrt( SpecialPairs(:NBondParms) % Parms(2) * factor1 * imol )
    SpecialPairs(:NBondParms) % Parms(1) = sqrt( SpecialPairs(:NBondParms) % Parms(1) * nano_2_angs    )

    ! SpecialMorse Potential :: Nothing special about it ... 
    NMorseParms = k - 1  
    allocate( SpecialMorse ( NMorseParms ) ) 

    forall(i=1:2) SpecialMorse(:NMorseParms) % MMSymbols(i) = Input2Chars(:NMorseParms,i) 

    SpecialMorse(:NMorseParms) % Parms(3) = Input2Reals(:NMorseParms,3)
    SpecialMorse(:NMorseParms) % Parms(2) = Input2Reals(:NMorseParms,1)
    SpecialMorse(:NMorseParms) % Parms(1) = Input2Reals(:NMorseParms,2) 
    
    ! conversion   
    SpecialMorse(:NMorseParms) % Parms(1) = SpecialMorse(:NMorseParms) % Parms(1) * factor1 * imol
    SpecialMorse(:NMorseParms) % Parms(2) = SpecialMorse(:NMorseParms) % Parms(2) * nano_2_angs
    SpecialMorse(:NMorseParms) % Parms(3) = SpecialMorse(:NMorseParms) % Parms(3) / nano_2_angs

!=====================================================================================
!  reads [ bondtypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ bondtypes ]" ) exit
    end do
        
    i = 1
    read_loop2: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop2
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop2
        if( trim(InputChars(i,1)) == "[  "  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop2
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,2) , InputIntegers(i,1), (InputReals(i,j) , j=1,2)

        backspace(33)

        select case ( InputIntegers(i,1) )

            case( 1 )  

                read(33,*) (dummy_char, k=1,2) , dummy_int , (InputReals(i,j) , j=1,2)
   
            case( 3 ) 
  
                read(33,*) (dummy_char, k=1,2) , dummy_int , (InputReals(i,j) , j=1,3)

            end select 

        i = i + 1
    end do read_loop2
    backspace(33)

    NbondsTypes = i - 1

    allocate( BondPairsSymbols    ( NbondsTypes , 2 ) )
    allocate( BondPairsParameters ( NbondsTypes , 3 ) , source = D_zero )
    allocate( Bond_Type           ( NbondsTypes     ) )
    allocate( funct_bond          ( NbondsTypes     ) )

    forall(i=1:2) BondPairsSymbols(:NbondsTypes,i) = InputChars(:NbondsTypes,i)
  
    Bond_Type( :NbondsTypes ) = InputIntegers(:NbondsTypes,1)  

    do i = 1 , NbondsTypes 
        select case( Bond_Type(i) )
            case( 1 ) ! Harmonic potential ...
                BondPairsParameters(i,1) = InputReals(i,2) * factor2 * imol
                BondPairsParameters(i,2) = InputReals(i,1) * nano_2_angs
                funct_bond(i) = "harm"
            case( 3 ) ! Morse potential ...
                BondPairsParameters(i,1) = InputReals(i,2) * factor1 * imol
                BondPairsParameters(i,2) = InputReals(i,1) * nano_2_angs
                BondPairsParameters(i,3) = InputReals(i,3) / nano_2_angs
                funct_bond(i) = "Mors"
        end select 
    end do

!=====================================================================================
!  reads [ angletypes ] ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ angletypes ]" ) exit
    end do

    i = 1
    read_loop3: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop3
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop3
        if( trim(InputChars(i,1)) == "[  "  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop3
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,3) , InputIntegers(i,1), (InputReals(i,j) , j=1,2 )


        i = i + 1
    end do read_loop3
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
    read_loop4: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop4
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop4
        if( trim(InputChars(i,1)) == "#in"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop4
        read(line,*,iostat=ioerr) (InputChars(i,k) , k=1,4) , InputIntegers(i,1)

        backspace(33)

        select case ( InputIntegers(i,1) )

            case( 1 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,3)

                ! V = k[1 + cos(n.phi - theta)]
                ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
                ! kdihed0(:,1) = phi_s   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(phi) ==> force constant (kJ.mol⁻¹) * factor1 * imol
                ! kdihed0(:,3) = n       ==> multiplicity (it will be) 
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 2 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,2)

                ! V = 1/2.k[cos(phi) - cos(phi0)]²
                ! factor1 = 1.0d26      <== Factor used to correct the unis readed fom Gromacs
                ! kdihed0(:,1) = xi_0   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(xi) ==> force constant (kJ.mol⁻¹.rad⁻²) * factor1 * imol
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 3 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,6)

                ! V = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
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

    end do read_loop4

    NdihedTypes = i - 1

    allocate( Dihed_Type         ( NdihedTypes     )                   )
    allocate( DihedSymbols       ( NdihedTypes , 4 )                   )
    allocate( DihedParameters    ( NdihedTypes , 6 ) , source = D_zero )

    forall(k=1:4) DihedSymbols(:NdihedTypes,k)    = InputChars(:NdihedTypes,k)
    forall(k=1:6) DihedParameters(:NdihedTypes,k) = InputReals(:NdihedTypes,k)

    Dihed_Type(:NdihedTypes) = InputIntegers(:NdihedTypes,1)

close(33)

deallocate( InputChars , InputReals , InputIntegers )
deallocate( Input2Chars , Input2Reals )
!=====================================================================================

do a = 1 , MM % N_of_species

!   Assigning to each specie the corresponding parameter ...
!   Bond parameters ...
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

    ! Dihedral parameters ...
    allocate( species(a) % kdihed0 ( species(a) % Ndiheds , 6 ) , source = D_zero )
    allocate( species(a) % harm    ( species(a) % Ndiheds     ) , source = I_zero )

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
                    ! harm(:)      = n
                    species(a) % kdihed0(n,1:2) = DihedParameters(k,1:2)
                    species(a) % harm(n)        = int(DihedParameters(k,3)) 
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
