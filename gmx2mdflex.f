! Convert gmx data to mdflex program :: verify gmx format in IO FORMATS 

module gmx2mdflex

use constants_m
use for_force
use type_m                 , only : dynemolworkdir , warning
use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles, DefinePairs, DefineMorse, debug_MM
use MM_tuning_routines     , only : SpecialBonds, SpecialAngs
use NonBondPairs           , only : Identify_NonBondPairs
use Babel_routines_m       , only : TO_UPPER_CASE
use setup_checklist        , only : Checking_Topology

private
 
public :: top2mdflex, itp2mdflex, SpecialPairs, SpecialPairs14, SpecialMorse 

    ! module variables ...
    character(3)     , allocatable   , save  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:)
    real*8           , allocatable   , save  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:)
    type(DefinePairs) , allocatable :: SpecialPairs(:)
    type(DefinePairs) , allocatable :: SpecialPairs14(:)
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
character(10)                   :: string , word(3)
character(200)                  :: line 
logical                         :: TorF
integer                         :: i , j , k , a , ioerr , dummy_int , counter , Nbonds , Nangs , Ndiheds , Ntorsion , Nbonds14 , N_of_atoms

allocate( InputChars    ( 20000 , 10 )                   )
allocate( InputReals    ( 20000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 20000 , 10 ) , source = I_zero )

! Reading different '.itp' species files ...
counter = 0
do a = 1 , MM % N_of_species

    string = species(a) % residue // '.itp'

    ! cloning the itp files into log.trunk ...
    call systemQQ("cp "//string//" log.trunk/.") 

    open(33, file=dynemolworkdir//string, status='old',iostat=ioerr,err=101)

        101 if( ioerr > 0 ) then
            print*, string,' file not found; terminating execution' ; stop
        end if

        write(*,'(/2a9)',advance='no') "Reading ", string

        ! start reading the molecular structure of species(a) ...
        do
            read(33,100) keyword
            if( trim(keyword) == "[ atoms ]" ) exit   ! <== looking for [ atoms ] in *.itp
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
 

        N_of_atoms = species(a) % N_of_atoms

        ! convert MMSymbol to upper case ...
        forall( i=1:N_of_atoms ) species(a)% atom(i)% MMSymbol = TO_UPPER_CASE( species(a)% atom(i)% MMSymbol )
        ! convert residues to upper case ...
        forall( i=1:N_of_atoms ) species(a)% atom(i)% residue  = TO_UPPER_CASE( species(a)% atom(i)% residue )

        i = 1
        do

            if( i > size(atom) ) exit

            if( trim(atom(i) % residue) == trim(species(a) % atom(1) % residue) ) then
                atom(i:i+N_of_atoms-1) % MM_charge = species(a) % atom(:N_of_atoms) % MM_charge
                i = i + N_of_atoms
            else
                i = i + 1
            end if

        end do
        rewind 33
!==============================================================================================
        ! Bonding parameters :: reading ...
        do
            read(33,100) keyword
            if( trim(keyword) == "[ bonds ]" ) exit   ! <== looking for [ bonds ] in *.itp
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
        rewind 33
!==============================================================================================
        ! Angle parameters :: reading ...
        do
            read(33,100) keyword
            if ( trim(keyword) == "[ angles ]" ) exit   ! <== looking for [ angles ] in *.itp
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
        allocate( species(a) % angle_type  ( Nangs     ) )

        forall(i=1:3) species(a) % angs(:Nangs,i) = InputIntegers(:Nangs,i)

        species(a) % funct_angle(:Nangs) = InputChars(:Nangs,1)

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
        do
            read(33,100,iostat=ioerr) keyword
            if ( trim(keyword) == "[ dihedrals ]" .OR. ioerr /= 0 ) exit   ! <== looking for [ dihedrals ] in *.itp
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
            allocate( species(a) % funct_dihed ( Ndiheds ) )

            forall(i=1:4) species(a) % diheds(:Ndiheds,i) = InputIntegers(:Ndiheds,i)

            species(a) % funct_dihed(:Ndiheds) = InputIntegers(:Ndiheds,5)

            ! define species(a) % dihedral_type ...
            CALL define_DihedralType( species(a) , Ndiheds )

        end if
        rewind 33

!----------------------------------------------------------------------------------------------
        ! the IMPROPER dihedrals must be at the END OF THE LIST ...
        Ntorsion = count( species(a)%dihedral_type /= "imp" )

        species(a)% NTorsions  = Ntorsion
        species(a)% NImpropers = Ndiheds - Ntorsion

        TorF = Checking_Topology( species(a)%bonds , species(a)%angs , species(a)%diheds(:Ntorsion,:) )
        If( TorF ) then
            CALL warning("error detected in Topology , check log.trunk/Topology.test.log")
            stop
        End If
!----------------------------------------------------------------------------------------------

!==============================================================================================
            ! Pairs 1-4 parameters :: reading ...
        do
            read(33,100,iostat=ioerr) keyword
            if( trim(keyword) == "[ pairs ]" .OR. ioerr /= 0 ) exit      ! <== looking for [ pairs ] in *.itp
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

        rewind(33)

!==============================================================================================
                         ! AD_HOC parameters :: reading ...
        do
            read(33,100,iostat=ioerr) keyword
            if ( trim(keyword) == "[ ad-hoc ]" .OR. ioerr /= 0 ) exit  ! <== looking for [ ad-hoc ] in *.itp
        end do

        if( trim(keyword) == "[ ad-hoc ]" ) then

           read(33,'(A)',iostat=ioerr) line
           read(line,*,iostat=ioerr) keyword
           keyword = to_upper_case(keyword)

           select case (keyword)

                  case( ";" )
                  ! Go on now, go, walk out the door ...

                  case( "FLEX" , "FLEX:") 

                      read_loop: do
                          read(33, '(A)', iostat=ioerr) line
                          if( ioerr /= 0 )          exit  read_loop
                          if( len_trim(line) == 0 ) cycle read_loop  ! <== empty line

                          read(line,*,iostat=ioerr) string
                          if( index(string,";") /= 0 ) cycle read_loop  ! <== comment line
                          if( trim(string) == "[  "  ) exit  read_loop  ! <== end of block

                          read(line,*, iostat=ioerr) ( word(j) , j=1,3 ) 

                          read(word(1),'(i)') k

                          keyword = to_upper_case(word(3))                                                                                                                                        
                          TorF = merge( .true. , .false. , any( [".TRUE.","TRUE","T","T_"] == keyword ) )

                          atom(k) % flex = TorF 

                      end do read_loop

                  case default
                      CALL warning("halting: check AD-HOC section of *.itp file")                                                    
                      stop

           end select

        end if
!==============================================================================================

    close(33)

    write(*,'(a9)') " << done "

end do

FF % residue  = adjustl(FF % residue)
FF % Symbol   = adjustl(FF % Symbol)
FF % MMSymbol = adjustl(FF % MMSymbol)

! convert MMSymbol to upper case ...
forall( i=1:size(FF) ) FF(i)% MMSymbol = TO_UPPER_CASE( FF(i)% MMSymbol )

deallocate( InputChars , InputIntegers )

100 format(a18)

end subroutine itp2mdflex
!
!
!
!==========================================
 subroutine top2mdflex( MM , species , FF )
!==========================================
implicit none 
type(MM_molecular)                  , intent(inout) :: species(:)
type(MM_system)                     , intent(inout) :: MM
type(MM_atomic)     , allocatable   , intent(inout) :: FF(:)
 
! local variables ...
character(3)    , allocatable   :: InputChars(:,:) , Input2Chars(:,:)
character(4)    , allocatable   :: funct_bond(:) , funct_angle(:)
real*8          , allocatable   :: InputReals(:,:) , Input2Reals(:,:)
integer         , allocatable   :: InputIntegers(:,:)
integer         , allocatable   :: Dihed_Type(:) , Bond_Type(:) , Angle_Type(:)
integer                         :: a , n , i , j , j1, k , ioerr , dummy_int , N_of_AtomTypes 
integer                         :: NbondsTypes , NangsTypes , NdihedTypes , NBondParms, NPairsParms , NMorseParms
character(3)                    :: dummy_char
character(18)                   :: keyword
character(200)                  :: line
logical                         :: flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8

allocate( InputChars    ( 10000 , 10 )                   )
allocate( Input2Chars   ( 10000 , 10 )                   )
allocate( InputReals    ( 10000 , 10 ) , source = D_zero )
allocate( Input2Reals   ( 10000 , 10 ) , source = D_zero )
allocate( InputIntegers ( 10000 , 10 ) , source = I_zero )

forcefield = 2           ! 1 = Born-Mayer (not implemented); 2 = Lennard-Jones (OK)
  
! cloning the topol.top file into log.trunk ...
call systemQQ("cp topol.top log.trunk/.") 

open(33, file=dynemolworkdir//'topol.top', status='old', iostat=ioerr, err=10)

!   file error msg ...
    10 if( ioerr > 0 ) stop '"topol.top" file not found; terminating execution'

    do
      read(33,100) keyword
      if( trim(keyword) == "[ defaults ]" ) exit
    end do

    i=1
    read_loop: do
        read(33,*,iostat=ioerr) dummy_char, MM% CombinationRule, dummy_char, MM% fudgeLJ, MM% fudgeQQ 
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

        ! convert MMSymbol to upper case ...
        InputChars(i,1) = TO_UPPER_CASE( InputChars(i,1) )

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
    ! factor1 = 1.0d26      <== Factor used to correct units read from Gromacs
    ! GAFF  vs  GMX  LJ parameters:
    ! -> epsilon_GAFF = epsilon_GMX / cal_2_J  
    ! -> sigma_GAFF = (sigma_GMX*10/2 ) * 2^(1/6)

    FF % eps = sqrt( FF % eps * factor1 * imol )
    FF % sig = FF % sig * nano_2_angs

    select case( MM % CombinationRule )

        case (2) 

            FF % sig = FF % sig / TWO

        case (3)

            FF % sig = sqrt( FF % sig )
            
    end select
    FF % sig14 = FF % sig
    FF % eps14 = FF % eps 

!=====================================================================================
!  NonBonding parameters :: reading ...
    do
        read(33,100) keyword
        if( trim(keyword) == "[ nonbond_params ]" ) exit
    end do

    i = 1
    k = 1
    read_loop5: do
        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop5
        read(line,*,iostat=ioerr) InputChars(i,1)        
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop5
        if( trim(InputChars(i,1)) == "[  "  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 )  cycle read_loop5 

        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,2) , InputIntegers(i,1) , (InputReals(i,j) , j=1,2)
        if( InputIntegers(i,1) == 3 ) then 
            backspace(33)
            read(line,*) (Input2Chars(k,j) , j=1,2) , dummy_int, (Input2Reals(k,j) , j=1,3) 
            k = k + 1 
            cycle read_loop5
        end if
        i = i + 1
    end do read_loop5
    backspace(33)

    NBondParms = i - 1

    If( NBondParms /= 0 ) then  

        allocate( SpecialPairs ( NbondParms ) )

        forall(i=1:2) SpecialPairs(:NBondParms) % MMSymbols(i) = InputChars(:NbondParms,i)

        SpecialPairs(:NBondParms) % Parms(1) = InputReals(:NbondParms,1) 
        SpecialPairs(:NBondParms) % Parms(2) = InputReals(:NbondParms,2)

        ! conversion 
        ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
        SpecialPairs(:NBondParms) % Parms(1) = SpecialPairs(:NBondParms) % Parms(1) * nano_2_angs    
        SpecialPairs(:NBondParms) % Parms(2) = SpecialPairs(:NBondParms) % Parms(2) * factor1 * imol 

    EndIf

    ! SpecialMorse Potential :: Nothing special about it ... 
    NMorseParms = k - 1  

    If( NMorseParms /= 0 ) then

        allocate( SpecialMorse ( NMorseParms ) ) 

        forall(i=1:2) SpecialMorse(:NMorseParms) % MMSymbols(i) = Input2Chars(:NMorseParms,i) 

        SpecialMorse(:NMorseParms) % Parms(3) = Input2Reals(:NMorseParms,3)
        SpecialMorse(:NMorseParms) % Parms(2) = Input2Reals(:NMorseParms,1)
        SpecialMorse(:NMorseParms) % Parms(1) = Input2Reals(:NMorseParms,2) 
        
        ! conversion   
        SpecialMorse(:NMorseParms) % Parms(1) = SpecialMorse(:NMorseParms) % Parms(1) * factor1 * imol
        SpecialMorse(:NMorseParms) % Parms(2) = SpecialMorse(:NMorseParms) % Parms(2) * nano_2_angs
        SpecialMorse(:NMorseParms) % Parms(3) = SpecialMorse(:NMorseParms) % Parms(3) / nano_2_angs

    EndIf

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

        ! convert MMSymbol to upper case ...
        forall( k=1:2 ) InputChars(i,k) = TO_UPPER_CASE( InputChars(i,k) )

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

        ! convert MMSymbol to upper case ...
        forall( k=1:3 ) InputChars(i,k) = TO_UPPER_CASE( InputChars(i,k) )

        backspace(33) 

        select case( InputIntegers(i,1) )
            case( 1 ) ! Harmonic potential ...
                read(33,*) (dummy_char,k=1,3), dummy_int, (InputReals(i,j), j=1,2)
            case( 5 ) ! Urey-Bradley potential ...
                read(33,*) (dummy_char,k=1,3), dummy_int, (InputReals(i,j), j=1,4)
        end select

        i = i + 1
    end do read_loop3
    backspace(33)
    backspace(33)

    NangsTypes = i - 1

    allocate( AngleSymbols    ( NangsTypes , 3 ) )
    allocate( AngleParameters ( NangsTypes , 4 ) , source = D_zero ) 
    allocate( Angle_Type      ( NangsTypes     ) )
    allocate( funct_angle     ( NangsTypes     ) )
    
    forall(i=1:3) AngleSymbols(:NangsTypes,i)  = InputChars(:NangsTypes,i)
    
    Angle_Type( :NangsTypes ) = InputIntegers(:NangsTypes,1)  

    do i = 1 , NangsTypes
        select case( Angle_Type(i) )
        ! conversion 
        ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
            case( 1 ) ! Harmonic potential ...
                AngleParameters(:NangsTypes,1) = InputReals(:NangsTypes,2) * factor1 * imol
                AngleParameters(:NangsTypes,2) = InputReals(:NangsTypes,1) * deg_2_rad
                funct_angle(i) = "harm"
            case( 5 ) ! Urey-Bradley potential ...
                AngleParameters(:NangsTypes,1) = InputReals(:NangsTypes,2) * factor1 * imol
                AngleParameters(:NangsTypes,2) = InputReals(:NangsTypes,1) * deg_2_rad
                AngleParameters(:NangsTypes,3) = InputReals(:NangsTypes,4) * factor2 * imol
                AngleParameters(:NangsTypes,4) = InputReals(:NangsTypes,3) * nano_2_angs
                funct_angle(i) = "urba"
        end select
    end do

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
        if( trim(InputChars(i,1)) == "[  "  ) exit
        if( trim(InputChars(i,1)) == "#in"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop4
        read(line,*,iostat=ioerr) (InputChars(i,k) , k=1,4) , InputIntegers(i,1)

        ! convert MMSymbol to upper case ...
        forall( k=1:4 ) InputChars(i,k) = TO_UPPER_CASE( InputChars(i,k) )

        backspace(33)

        select case ( InputIntegers(i,1) )

            case( 1 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,3)

                !============================================================================
                ! V = k[1 + cos(n.phi - theta)]
                ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
                ! kdihed0(:,1) = phi_s   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(phi) ==> force constant (kJ/mol) * factor1 * imol
                ! kdihed0(:,3) = n       ==> multiplicity (it will be) 
                !============================================================================
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 2 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,2)

                !============================================================================
                ! V = 1/2.k[cos(phi) - cos(phi0)]²
                ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
                ! kdihed0(:,1) = xi_0   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(xi) ==> force constant (kJ/(mol.rad^2)) * factor1 * imol
                !============================================================================
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 3 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,6)

                !============================================================================
                ! V = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
                ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
                ! kdihed0(:,1) = C0 (kJ/mol) * factor1 * imol
                ! kdihed0(:,2) = C1 (kJ/mol) * factor1 * imol
                ! kdihed0(:,3) = C2 (kJ/mol) * factor1 * imol
                ! kdihed0(:,4) = C3 (kJ/mol) * factor1 * imol
                ! kdihed0(:,5) = C4 (kJ/mol) * factor1 * imol
                ! kdihed0(:,6) = C5 (kJ/mol) * factor1 * imol
                !============================================================================
                InputReals(i,1:6) = InputReals(i,1:6) * factor1 * imol

            case( 4 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,3)

                !============================================================================
                ! V = k[1 + cos(n.phi - theta)] (improper; same as 1)
                ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
                ! kdihed0(:,1) = phi_s   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(phi) ==> force constant (kJ.mol⁻¹) * factor1 * imol
                ! kdihed0(:,3) = n       ==> multiplicity (it will be) 
                !============================================================================
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

            case( 9 )

                read( 33 , * ) (dummy_char, k=1,4) , dummy_int , (InputReals(i,k) , k=1,3)

                !============================================================================
                ! V = k[1 + cos(n.phi - theta)] (multiple; same as 1)
                ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
                ! kdihed0(:,1) = phi_s   ==> angle (deg) * deg_2_rad
                ! kdihed0(:,2) = K_(phi) ==> force constant (kJ.mol⁻¹) * factor1 * imol
                ! kdihed0(:,3) = n       ==> multiplicity (it will be) 
                !============================================================================
                InputReals(i,1) = InputReals(i,1) * deg_2_rad
                InputReals(i,2) = InputReals(i,2) * factor1 * imol

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

!=====================================================================================
!  [ pairtypes ] parameters :: reading ...
    do
        read(33,100,iostat=ioerr) keyword
        if( trim(keyword) == "[ pairtypes ]" .OR. ioerr /= 0 ) exit
    end do

    if( trim(keyword) == "[ pairtypes ]" ) then
    i = 1
    k = 1
    read_loop6: do

        read(33, '(A)', iostat=ioerr) line
        if ( ioerr /= 0 ) exit read_loop6
        read(line,*,iostat=ioerr) InputChars(i,1), InputChars(i,2)
        if( index(InputChars(i,1),";") /= 0 ) cycle read_loop6
        if( trim(InputChars(i,1)) == "#in"  ) exit
        if( ioerr > 0  ) exit
        if( ioerr /= 0 ) cycle read_loop6
        read(line,*,iostat=ioerr) (InputChars(i,j) , j=1,2) , dummy_int, (InputReals(i,j), j=1,2)

        i = i + 1
    
    end do read_loop6
    backspace(33)

    NPairsParms = i - 1

    If( NPairsParms /= 0 ) then

        allocate( SpecialPairs14 ( NPairsParms ) )

        forall(i=1:2) SpecialPairs14(:NPairsParms) % MMSymbols(i) = InputChars(:NPairsParms,i)

        SpecialPairs14(:NPairsParms) % Parms(1) = InputReals(:NPairsParms,1)
        SpecialPairs14(:NPairsParms) % Parms(2) = InputReals(:NPairsParms,2)

        ! conversion 
        ! factor1 = 1.0d26      <== Factor used to correct the units read from Gromacs
        SpecialPairs14(:NPairsParms) % Parms(1) = sqrt( SpecialPairs14(:NPairsParms)%Parms(1) * nano_2_angs    )
        SpecialPairs14(:NPairsParms) % Parms(2) = sqrt( SpecialPairs14(:NPairsParms)%Parms(2) * factor1 * imol )

    EndIf

  end if
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
            ! factor1 = 1.0d26      <== Factor used to correct the units read fom Gromacs
            where( species(a) % funct_angle == SpecialAngs(k) % label ) species(a) % kang0(:,1) = SpecialAngs(k) % kang0(1) * factor1 * imol
            where( species(a) % funct_angle == SpecialAngs(k) % label ) species(a) % kang0(:,2) = SpecialAngs(k) % kang0(2) * deg_2_rad
        end do
    end if

    !=============================================================================
    ! Dihedral parameters ...
    allocate( species(a) % kdihed0 ( species(a) % Ndiheds , 15 ) , source = D_zero )

    read_loop0: do n = 1 , species(a) % Ndiheds
        ! control variables to multiple dihs ...
        j = 0 ; j1 = 0

        read_loop7: do k = 1 , NdihedTypes 

            !============================================
            ! if funct = 1 (cos)
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        
            ! Eq. 4.60 (GMX manual 5.0.5)
            !============================================

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

                flag5 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 1 )

                flag6 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 1 )

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    !===============================
                    species(a) % kdihed0(n,1:3) = DihedParameters(k,1:3)
                    cycle read_loop0
                end if

            end if

            !===============================
            ! if funct = 2 (harm)
            ! V = 1/2.k ( xi - xi_0 )^2
            ! Eq. 4.59 (GMX manual 5.0.5)
            !===============================

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

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    !======================================
                    ! kdihed0(:,1) = xi_0 (deg)
                    ! kdihed0(:,2) = k_xi [ kJ/(mol.rad^2) ]
                    !======================================
                    species(a) % kdihed0(n,1:2) = DihedParameters(k,1:2)
                    cycle read_loop0
                end if

            end if

            !==============================
            ! if funct = 3 (cos3)
            ! V = C0 + C1 * cos( phi - 180 ) + C2 * cos^2( phi - 180 ) + C3 * cos^3( phi - 180 ) + C4 * cos^4( phi - 180 ) + C5 * cos^5( phi - 180 ) 
            ! Eq. 4.61 (GMX manual 5.0.5)
            !==============================
            
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

                flag5 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 3 )

                flag6 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 3 )
   
                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    species(a) % kdihed0(n,1:6) = DihedParameters(k,1:6) 
                    cycle read_loop0
                end if

            end if

            !===========================================
            ! if funct = 4 (imp) (improper)
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        
            ! Eq. 4.60 (GMX manual 5.0.5)
            !===========================================
            if( species(a) % funct_dihed(n) == 4 ) then

                flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag3 = ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag5 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag6 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag7 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,2)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                flag8 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. & 
                        ( adjustl(DihedSymbols(k,3)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 4 )

                if( flag1 .OR. flag2 .OR. flag3 .OR. flag4 .OR. flag5 .OR. flag6 .OR. flag7 .OR. flag8 ) then
                    !================================
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    !================================
                    species(a) % kdihed0(n,1:3) = DihedParameters(k,1:3)
                    cycle read_loop0
                end if

            end if
           
            !===========================================
            ! if funct = 9 (chrm) (multiple)
            ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        
            ! Eq. 4.60 (GMX manual 5.0.5) 
            !===========================================
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
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag4 = ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag5 = ( adjustl(DihedSymbols(k,1)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                flag6 = ( adjustl(DihedSymbols(k,4)) == 'X' ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) .AND. &
                        ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) ) .AND. &
                        ( Dihed_Type(k) == 9 )

                if( flag1 .OR. flag2 ) then
                    !================================
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n
                    !================================
                    if( j1 > 0 ) species(a) % kdihed0(n,:) = D_zero
                    if( j1 > 0 ) j1 = 0
                    species(a) % kdihed0(n,3*(j+j1)+1) = DihedParameters(k,1)
                    species(a) % kdihed0(n,3*(j+j1)+2) = DihedParameters(k,2)
                    species(a) % kdihed0(n,3*(j+j1)+3) = DihedParameters(k,3)
                    j = j + 1
                    cycle read_loop7
                end if

                if( flag3 .OR. flag4 .OR. flag5 .OR. flag6 ) then
                    !================================
                    ! kdihed0(:,1) = phi_s (deg)
                    ! kdihed0(:,2) = k_phi (kJ/mol)
                    ! kdihed0(:,3) = n 
                    !================================
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

        case( 2 ) ! V = 1/2.k( xi - xi_0 )^2

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
end  module gmx2mdflex
