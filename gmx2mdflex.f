!======================================================================
! Convert gmx data to mdflex program :: verify gmx format in IO FORMATS 
!======================================================================
module gmx2mdflex

 use constants_m
 use for_force
 use MM_types               , only : MM_atomic, MM_molecular, MM_system, DefineBonds, DefineAngles
 use MM_tuning_routines     , only : SpecialBonds, SpecialAngs

 private
 
 public :: top2mdflex, itp2mdflex

   ! module variables ...
   integer                       :: Nbonds, Nangs, Ndiheds
   real*8                        :: fact14
   real*8         , allocatable  :: BondPairsParameters(:,:), AngleParameters(:,:), DihedParameters(:,:)
   character(3)   , allocatable  :: BondPairsSymbols(:,:), AngleSymbols(:,:), DihedSymbols(:,:)
   character(15)  , allocatable  :: special_bonds(:), special_angles(:)

contains
!
!
!=======================================
 subroutine top2mdflex( MM , atom , FF )
!=======================================
 implicit none 
 type(MM_system), intent(inout) :: MM
 type(MM_atomic), intent(inout) :: atom(:)
 type(MM_atomic), allocatable, intent(inout) :: FF(:)
 
 ! local variables ...
 integer      :: i, j, ioerr, dummy_int, ilines, dihed_type, k
 real*8       :: dummy_real, r0, k0, theta0, ktheta0
 character(3) :: dummy_char, string, bondatm1, bondatm2, angatm1, angatm2, angatm3
 character(3) :: dihedatm1, dihedatm2, dihedatm3, dihedatm4

 forcefield = 2
  
 open (8, file='topol.top', status='old', iostat=ioerr, err=10)
    ! msg file error ...
    10 if( ioerr > 0 ) then
       stop ' "topol.top" file not found; terminating execution'
    end if
             
    ! finding the number of atomtypes ...
    do 50 ilines = 1, 2
         read(8,*)
    50 continue
  
    read(8,*) dummy_int, MM % CombinationRule, dummy_char, fact14, fact14
    
    do 51 ilines = 1, 4
         read(8,*)
    51 continue

    i = 1
    do  
       read(8,*,iostat=ioerr) dummy_char, dummy_real, dummy_real, dummy_char, dummy_real, dummy_real
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do

    ! AtomTypes stuffs :: reading ...
    MM % N_of_AtomTypes = i-1
    rewind(8)

    do 60 ilines = 1, 7
       read(8,*)
    60  continue
  
 !   allocate( FF(MM % N_of_AtomTypes) )
    do i = 1 , MM % N_of_AtomTypes
       read(8,*,iostat=ioerr) FF(i) % MMSymbol, FF(i) % mass, FF(i) % MM_charge, dummy_char, FF(i) % sig, FF(i) % eps
       FF(i) % MMSymbol = adjustr( FF(i) % MMSymbol )
       FF(i) % eps   = FF(i) % eps * 1.d26 * imol
       FF(i) % eps   = SQRT( FF(i) % eps )
       FF(i) % sig   = FF(i) % sig * 1.d1
       select case ( MM % CombinationRule )
            case (2) 
            FF(i) % sig = FF(i) % sig / TWO
            case (3)
            FF(i) % sig = sqrt( FF(i) % sig )
       end select
       FF(i) % my_id = i

       where( atom % MMSymbol == FF(i) % MMSymbol ) atom % eps    = FF(i) % eps
       where( atom % MMSymbol == FF(i) % MMSymbol ) atom % sig    = FF(i) % sig
    end do

    ! Bonding parameters :: reading ...
    do 70 ilines = 1, 5
       read(8,*)
    70 continue

    i = 1
    do
       read(8,*,iostat=ioerr) dummy_char, dummy_char, dummy_int, dummy_real, dummy_real
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    Nbonds = i-1

    do j = 1, Nbonds+3
       backspace(8)
    end do

    allocate( BondPairsSymbols    (Nbonds,2) )
    allocate( BondPairsParameters (Nbonds,2) )
    do i = 1 , Nbonds
       read(8,*,iostat=ioerr) BondPairsSymbols(i,1:2), dummy_int, r0, k0
       BondPairsParameters(i,1) = k0
       BondPairsParameters(i,2) = r0
    end do
    BondPairsParameters(:,2) = BondPairsParameters(:,2) * 1.d1
    BondPairsParameters(:,1) = BondPairsParameters(:,1) * 1.d24 * imol

    ! Angle parameters :: reading ...
    do 80 ilines = 1, 5
       read(8,*)
    80 continue

    i = 1
    do
       read(8,*,iostat=ioerr) dummy_char, dummy_char, dummy_char, dummy_int, dummy_real, dummy_real
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    Nangs = i-1

    do j = 1, Nangs+4
       backspace(8)
    end do

    allocate( AngleSymbols    (Nangs,3) )
    allocate( AngleParameters (Nangs,2) )
    do i = 1 , Nangs
       read(8,*,iostat=ioerr) AngleSymbols(i,1:3), dummy_int, theta0, ktheta0
       AngleParameters(i,1) = ktheta0
       AngleParameters(i,2) = theta0
    end do
    AngleParameters(:,2) = AngleParameters(:,2) * pi / 180.d0
    AngleParameters(:,1) = AngleParameters(:,1) * 1.d26 * imol
 
   ! Dihedral parameters :: reading ...
    do 90 ilines = 1, 5
       read(8,*)
    90 continue

    i = 1
    do
       read(8,*,iostat=ioerr) dummy_char, dummy_char, dummy_char, dummy_char, dummy_int, dummy_real, dummy_real, &
                            & dummy_real, dummy_real, dummy_real, dummy_real
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    Ndiheds = i-1

    do j = 1, Ndiheds + 5
       backspace(8)
    end do

    allocate( DihedSymbols    (Ndiheds,4) )
    allocate( DihedParameters (Ndiheds,6) )
    do i = 1 , Ndiheds
       read(8,*,iostat=ioerr) DihedSymbols(i,1:4), dihed_type, DihedParameters(i,1:6) 
    end do
    DihedParameters(:,1:6) = DihedParameters(:,1:6) * 1.d26 * imol
    if ( Ndiheds == 0 ) then 
    Dihedral_Potential_Type = 'none'
    else
    Dihedral_Potential_Type = 'ryck'
    end if    

 close(8)


end subroutine top2mdflex
!
!
!
!================================================
 subroutine itp2mdflex( MM , atom , species , FF)
!================================================
 implicit none
 type(MM_system)   , intent(in)    :: MM
 type(MM_atomic)   , intent(inout) :: atom(:)
 type(MM_atomic)   , intent(inout) :: FF(:)
 type(MM_molecular), intent(inout) :: species(:)

 ! local variables ...
 integer          :: i, j, k, n, a, ioerr, ilines, dummy_int
 real*8           :: dummy_real, factor
 character(3)     :: dummy_char, angatm1, angatm2, angatm3, bondatm1, bondatm2
 character(3)     :: dihedatm1, dihedatm2, dihedatm3, dihedatm4
 character(10)    :: string
 logical          :: flag1, flag2

do a = 1, MM % N_of_species
 ! Reading different '.itp' species files ...
 string = species(a) % residue // '.itp'
 
 open (9, file=string, status='old',iostat=ioerr,err=101)
 
 101 if( ioerr > 0 )then
     print*, string,' file not found; terminating execution'; stop
     end if

    ! finding the number of atoms in the molecule ...
    do 50 ilines = 1, 6
       read(9,*)
    50 continue

    i = 1
    do
       read(9,*,iostat=ioerr) dummy_int, dummy_char, dummy_int, dummy_char, dummy_char, dummy_int, dummy_real, dummy_real
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do

    ! Atom & AtomTypes stuffs :: reading ...
    species(a) % N_of_atoms = i - 1 
    rewind(9)

    do 60 ilines = 1, 6
       read(9,*)
    60  continue
    
    allocate( species(a) % atom(species(a) % N_of_atoms) )

    do i = 1 , species(a) % N_of_atoms
       read(9,*,iostat=ioerr) species(a) % atom(i) % my_id,         &
                              species(a) % atom(i) % MMSymbol,      &
                              dummy_int,                            &
                              species(a) % atom(i) % residue,       &
                              species(a) % atom(i) % Symbol,        &
                              dummy_int,                            &
                              species(a) % atom(i) % MM_charge,     &
                              species(a) % atom(i) % mass

       species(a) % atom(i) % MMSymbol = adjustr(species(a) % atom(i) % MMSymbol)
    end do
 
    do i = 1 , species(a) % N_of_atoms
       where( atom % MMSymbol == species(a) % atom(i) % MMSymbol ) atom % MM_charge = species(a) % atom(i) % MM_charge
       FF(i) % eps = atom(i) % eps
       FF(i) % sig = atom(i) % sig
       FF(i) % MM_charge = atom(i) % MM_charge
       FF(i) % MMSymbol = atom(i) % MMSymbol
    end do

    ! Bonding parameters :: reading ...
    do 70 ilines = 1, 5
       read(9,*)
    70 continue

    i = 1
    do
       read(9,*,iostat=ioerr) dummy_int, dummy_int
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    species(a) % Nbonds = i

    do j = 1, species(a) % Nbonds+3
       backspace(9)
    end do

    allocate( species(a) % bonds(species(a)%Nbonds,2) )
    allocate( special_bonds(species(a)%Nbonds)        )
    do i = 1 , species(a) % Nbonds
       read(9,*,iostat=ioerr) species(a) % bonds(i,1:2), special_bonds(i)
    end do
    
   ! Pairs interactions :: reading ...
    do 71 ilines = 1, 4
       read(9,*)
    71 continue

    i = 1
    do
       read(9,*,iostat=ioerr) dummy_int, dummy_int, dummy_int
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    species(a) % Nbonds14 = i-1

    do j = 1, species(a) % Nbonds14+3
       backspace(9)
    end do
    
    allocate( species(a) % bonds14 (species(a)%Nbonds14,2) )
    allocate( species(a) % fact14  (species(a)%Nbonds14)   )
    do i = 1 , species(a) % Nbonds14
       read(9,*,iostat=ioerr) species(a) % bonds14(i,1:2), dummy_int
       species(a) % fact14(i) = fact14
    end do

    ! Angle parameters :: reading ...
    do 72 ilines = 1, 4
       read(9,*)
    72 continue

    i = 1
    do
       read(9,*,iostat=ioerr) dummy_int, dummy_int, dummy_int
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    species(a) % Nangs = i-1

    do j = 1, species(a) % Nangs+3
       backspace(9)
    end do

    allocate( species(a) % angs(species(a)%Nangs,3) )
    allocate( special_angles(species(a)%Nangs)      )
    do i = 1 , species(a) % Nangs
       read(9,*,iostat=ioerr) species(a) % angs(i,1:3), special_angles(i)
    end do

    ! Dihedrals parameters :: reading ...
    do 73 ilines = 1, 4
       read(9,*)
    73 continue

    i = 1
    do
       read(9,*,iostat=ioerr) dummy_int, dummy_int, dummy_int, dummy_int, dummy_int
       if ( ioerr /= 0 ) exit
       i = i + 1
    end do
    species(a) % Ndiheds = i-1
    
    do j = 1, species(a) % Ndiheds+1
       backspace(9)
    end do

    allocate( species(a) % diheds(species(a)%Ndiheds,4) )
    do i = 1 , species(a) % Ndiheds
       read(9,*,iostat=ioerr) species(a) % diheds(i,1:4), dummy_int
    end do

!====================================================================
! Assigning to each specie the corresponding parameter ...
 ! Bond parameters ...
 allocate( species(a) % kbond0(species(a) % Nbonds,2) )
 do k = 1 , Nbonds
    do n = 1 , species(a) % Nbonds
       
       flag1 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) ) &
            .AND.                                                                                                   &
               ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) ) 

       flag2 = ( adjustl(species(a) % atom(species(a) % bonds(n,1)) % MMSymbol) == adjustl(BondPairsSymbols(k,2)) ) &
            .AND.                                                                                                   & 
               ( adjustl(species(a) % atom(species(a) % bonds(n,2)) % MMSymbol) == adjustl(BondPairsSymbols(k,1)) ) 

       if ( flag1 .OR. flag2 ) then 
          species(a) % kbond0(n,1) = BondPairsParameters(k,1)
          species(a) % kbond0(n,2) = BondPairsParameters(k,2)
       end if
     end do
 end do     

 if ( allocated(SpecialBonds) ) then
    do k = 1, size(SpecialBonds)
       where( special_bonds == SpecialBonds(k) % nome ) species(a) % kbond0(:,1) = SpecialBonds(k) % kbond0(1) * 1.d24 * imol
       where( special_bonds == SpecialBonds(k) % nome ) species(a) % kbond0(:,2) = SpecialBonds(k) % kbond0(2) * 1.d1
    end do
 end if

 deallocate( BondPairsParameters, BondPairsSymbols, special_bonds, SpecialBonds ) 
 ! Angle parameters ...
 allocate( species(a) % kang0(species(a) % Nangs,2) )
 do k = 1 , Nangs
    do n = 1 , species(a) % Nangs

       flag1 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,1)) ) &
            .AND.                                                                                              &
               ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) &
            .AND.                                                                                              &
               ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,3)) )


       flag2 = ( adjustl(species(a) % atom(species(a) % angs(n,1)) % MMSymbol) == adjustl(AngleSymbols(k,3)) ) &
            .AND.                                                                                              &
               ( adjustl(species(a) % atom(species(a) % angs(n,2)) % MMSymbol) == adjustl(AngleSymbols(k,2)) ) &
            .AND.                                                                                              &
               ( adjustl(species(a) % atom(species(a) % angs(n,3)) % MMSymbol) == adjustl(AngleSymbols(k,1)) )

       if ( flag1 .OR. flag2 ) then
          species(a) % kang0(n,1) = AngleParameters(k,1)
          species(a) % kang0(n,2) = AngleParameters(k,2)
       end if

     end do
 end do

 if ( allocated(SpecialAngs) ) then
    do k = 1, size(SpecialAngs)
       where( special_angles == SpecialAngs(k) % nome ) species(a) % kang0(:,1) = SpecialAngs(k) % kang0(1) * 1.d26 * imol
       where( special_angles == SpecialAngs(k) % nome ) species(a) % kang0(:,2) = SpecialAngs(k) % kang0(2) * pi / 180.d0
    end do
 end if

 deallocate( AngleParameters, AngleSymbols, special_angles, SpecialAngs )

 ! Dihedral parameters ...
 allocate( species(a) % kdihed0(species(a) % Ndiheds,6) )
 do k = 1 , Ndiheds
    do n = 1 , species(a) % Ndiheds

       flag1 = ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,4)) )


       flag2 = ( adjustl(species(a) % atom(species(a) % diheds(n,4)) % MMSymbol) == adjustl(DihedSymbols(k,1)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,3)) % MMSymbol) == adjustl(DihedSymbols(k,2)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,2)) % MMSymbol) == adjustl(DihedSymbols(k,3)) ) &
            .AND.                                                                                                &
               ( adjustl(species(a) % atom(species(a) % diheds(n,1)) % MMSymbol) == adjustl(DihedSymbols(k,4)) )

       if ( flag1 .OR. flag2 ) then
          species(a) % kdihed0(n,1) = DihedParameters(k,1)
          species(a) % kdihed0(n,2) = DihedParameters(k,2)
          species(a) % kdihed0(n,3) = DihedParameters(k,3)
          species(a) % kdihed0(n,4) = DihedParameters(k,4)
          species(a) % kdihed0(n,5) = DihedParameters(k,5)
          species(a) % kdihed0(n,6) = DihedParameters(k,6)
       end if

     end do
 end do

 deallocate( DihedParameters, DihedSymbols )

!
end do

close(9)


end subroutine itp2mdflex
!
!
!
end  module gmx2mdflex
