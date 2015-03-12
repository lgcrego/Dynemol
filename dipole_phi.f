module DP_potential_m

    use type_m                  
    use constants_m
    use blas95
    use f95_precision
    use parameters_m            , only : PBC
    use Structure_Builder       , only : Extended_Cell
    use Semi_Empirical_Parms    , only : atom
    use DP_FMO_m                , only : DP_FMO_analysis
    use Multipole_Routines_m    , only : Util_Multipoles
    use PBC_m                   , only : Generate_Periodic_DPs

    public :: Molecular_DPs , DP_phi

    private

! module variables ...       
    type(dipoles) , save  :: DP_mols_pbc  

! module parameters ...       
    real*8  , parameter   :: DP_potential_factor = 2.9979255d0  ! <== e*p(Debye)/(4*Pi*epsilon_0) : eV * Angs^2

contains
!
!
!
!=============================
 subroutine Molecular_DPs( a )
!=============================
implicit none
type(structure) , intent(inout) :: a

! local variables ...
integer       :: i , j , N_of_DP_mols 
type(dipoles) :: solvent , solute , DP_mols

! initialize multipole routines ...
CALL Util_Multipoles

! garbage collection before restart calculations ...
CALL DeAllocate_DPs( DP_mols_pbc , flag = "garbage" )

! solvent and solute must not overlap ...
IF( any( (a%fragment=="S") .AND. (a%solute) ) ) then
    Print*, " >>> solvent and solute overlap in Molecular_DPs <<< : execution stopped" 
    stop
end If

If( any(a%solute) ) then
  
    CALL Build_DP_mols( a , solvent , instance = "solvent" )
    CALL Build_DP_mols( a , solute  , instance = "solute"  )

    N_of_DP_mols = a%N_of_Solvent_Molecules + a%N_of_Solute_Molecules

    CALL DeAllocate_DPs( DP_mols , N_of_DP_mols , flag = "alloc" )

    DP_mols % CC( 1:a%N_of_Solvent_Molecules , :) = solvent % CC 
    DP_mols % DP( 1:a%N_of_Solvent_Molecules , :) = solvent % DP 
    DP_mols % nr( 1:a%N_of_Solvent_Molecules    ) = solvent % nr 

    DP_mols % CC( a%N_of_Solvent_Molecules+1:N_of_DP_mols , : ) = solute % CC
    DP_mols % DP( a%N_of_Solvent_Molecules+1:N_of_DP_mols , : ) = solute % DP
    DP_mols % nr( a%N_of_Solvent_Molecules+1:N_of_DP_mols     ) = solute % nr

    CALL DeAllocate_DPs( solvent , flag = "dealloc" )
    CALL DeAllocate_DPs( solute  , flag = "dealloc" )

else
 
    CALL Build_DP_mols( a , solvent , instance = "solvent" )

    CALL DeAllocate_DPs( DP_mols , a%N_of_Solvent_Molecules , flag = "alloc" )

    DP_mols = solvent 

    CALL DeAllocate_DPs( solvent , flag = "dealloc" )

end If

! generate periodic structure of solvent molecules ; if PBCx=PBCy=PBCz=0 ==> DP_mols_pbc = DP_mols ...
CALL Generate_Periodic_DPs( a ,  DP_mols%CC     , DP_mols%DP     , DP_mols%nr     ,& 
                                 DP_mols_pbc%CC , DP_mols_pbc%DP , DP_mols_pbc%nr )

CALL DeAllocate_DPs( DP_mols , flag = "dealloc" )

end subroutine Molecular_DPs
!
!
!
!===================================================
 subroutine Build_DP_mols( a , molecule , instance )
!===================================================
implicit none
type(structure) , intent(inout) :: a
type(dipoles)   , intent(out)   :: molecule
character(*)    , intent(in)    :: instance


! local variables ...
integer                :: i , j , i1 , i2 , nr , last_nr , first_nr , nr_atoms
real*8                 :: total_valence , norm , vector(3)
real*8 , allocatable   :: Q_center(:,:) , DP_FMO(:,:) , Qi_Ri(:,:) 

select case( instance )

    case( "solvent" )

        ! find positions of solvent atoms in structure ...
        first_nr = minval( a%nr , a%fragment == "S" )
        last_nr  = maxval( a%nr , a%fragment == "S" )

        a%N_of_Solvent_Molecules = (last_nr - first_nr + 1)

        ! consistency check ...
        Print 157 , a%N_of_Solvent_Molecules

    case( "solute" )

        ! find positions of solute atoms in structure ...
        first_nr = minval( a%nr , a%solute )
        last_nr  = maxval( a%nr , a%solute )

        a%N_of_Solute_Molecules = (last_nr - first_nr + 1)

        ! consistency check ...
        Print 158 , a%N_of_Solute_Molecules 

end select

allocate(Q_center (last_nr,3) , source = 0.d0)
allocate(DP_FMO   (last_nr,3) , source = 0.d0)    
allocate(Qi_Ri    (a%atoms,3) , source = 0.d0)     

do nr = first_nr , last_nr 

    ! No of atoms with tag nr ...
    nr_atoms = count( a%nr == nr ) 

    ! position of nr residue in structure ...
    i1 = minloc( a%nr , 1 , a%nr == nr ) 
    i2 = (i1-1) + nr_atoms 

    !-----------------------------------------------------------------------------------------------
    ! calculate Center_of_Charge for each molecule: sum_i = (q_i * vec{r}_i) / sum_i q_i ...

    forall( j=1:3 , i=i1:i2 ) Qi_Ri(i,j) = a%Nvalen(i) * a%coord(i,j)

    total_valence = sum( [ (a%Nvalen(i) , i=i1,i2) ] )

    forall(j=1:3) Q_center(nr,j) = sum( Qi_Ri(i1:i2,j) ) / total_valence

    !-----------------------------------------------------------------------------------------------
    ! calculate the dipole vector ...

    CALL DP_FMO_analysis( a , Q_center(nr,:) , DP_FMO(nr,:) , nr ) 

    !-----------------------------------------------------------------------------------------------

end do

CALL DeAllocate_DPs( molecule , last_nr-first_nr+1 , flag = "alloc" )

! NOTICE: dipole moment is multiplied by DP_potential_factor ...
do j = 1 , 3 
    molecule%CC( : , j ) = Q_center( first_nr:last_nr , j ) 
    molecule%DP( : , j ) = DP_FMO  ( first_nr:last_nr , j ) * DP_potential_factor
end do

! save list of nr indices for use in DP_phi ...
molecule%nr = [ ( i , i=first_nr,last_nr ) ]

deallocate( Q_center , DP_FMO , Qi_Ri )

include 'formats.h'

end subroutine Build_DP_mols
!
!
!
!=====================================
 pure function DP_phi( a , b , basis )
!=====================================
implicit none
integer         , intent(in) :: a , b
type(STO_basis) , intent(in) :: basis(:)

! local variables ...
integer                 :: i , j , N_of_DP 
integer , allocatable   :: nr_Mols(:)
real*8                  :: hard_core , cut_off_radius 
real*8                  :: midpoint_ab(3)
real*8                  :: DP_phi(4)
real*8  , allocatable   :: distance(:) , distance_ALL(:) , mol_phi(:) 
real*8  , allocatable   :: vector(:,:) , vector_ALL(:,:) , decay(:,:) , DP_Mols(:,:) , mol_phi2(:,:)
logical , allocatable   :: mask(:)

! combination rule for solvation hardcore shell ...
hard_core = ( basis(a)%solvation_hardcore + basis(b)%solvation_hardcore ) / two

! midpoint between atoms a & b ...
midpoint_ab(1) = ( basis(a)%x + basis(b)%x ) / two
midpoint_ab(2) = ( basis(a)%y + basis(b)%y ) / two
midpoint_ab(3) = ( basis(a)%z + basis(b)%z ) / two

! total number of dipoles in PBC ...
N_of_DP = size( DP_mols_pbc%DP(:,1) )

!------------------------------------------------------------------------
! choice of solvent droplet .vs. PBC box ...
!------------------------------------------------------------------------

If( sum(PBC) == 0) then

    allocate( vector  ( N_of_DP , 3 ) , source = D_zero )
    allocate( DP_Mols ( N_of_DP , 3 ) , source = D_zero )
    allocate( decay   ( N_of_DP , 3 ) , source = D_zero )
    allocate( distance( N_of_DP     ) , source = D_zero )
    allocate( mol_phi ( N_of_DP     ) , source = D_zero )
    allocate( nr_Mols ( N_of_DP     ) , source = I_zero )

    do i = 1 , N_of_DP

        vector(i,:) = midpoint_ab(:) - DP_mols_pbc%CC(i,:)
        distance(i) = sqrt( sum( vector(i,:)*vector(i,:) ) )

    end do

    DP_Mols = DP_mols_pbc%DP
    nr_Mols = DP_mols_pbc%nr

else

    ! maximum distance from midpoint a-b ...
    cut_off_radius = minval(Extended_Cell%T_xyz) * two / three

    allocate( vector_ALL   ( N_of_DP , 3 ) , source = D_zero  )
    allocate( distance_ALL ( N_of_DP     ) , source = D_zero  )
    allocate( mask         ( N_of_DP     ) , source = .false. )

    do i = 1 , N_of_DP

        vector_ALL(i,:) = midpoint_ab(:) - DP_mols_pbc%CC(i,:)
        distance_ALL(i) = sqrt( sum( vector_ALL(i,:)*vector_ALL(i,:) ) )

    end do

    mask = ( (distance_ALL > hard_core) .AND. (distance_ALL < cut_off_radius) )

    ! redefine range of interaction with molecules inside cut_off radius ...
    N_of_DP = count( mask ) 

    allocate( vector  ( N_of_DP , 3 ) , source = D_zero )
    allocate( DP_Mols ( N_of_DP , 3 ) , source = D_zero )
    allocate( decay   ( N_of_DP , 3 ) , source = D_zero )
    allocate( distance( N_of_DP     ) , source = D_zero )
    allocate( mol_phi ( N_of_DP     ) , source = D_zero )
    allocate( nr_Mols ( N_of_DP     ) , source = I_zero )

    do j = 1 , 3 
        vector(:,j)  = pack( vector_ALL(:,j)     , mask , vector(:,j)  )
        DP_Mols(:,j) = pack( DP_mols_pbc%DP(:,j) , mask , DP_Mols(:,j) )
    end do

    nr_Mols  = pack( DP_mols_pbc%nr , mask , nr_Mols  )
    distance = pack( distance_ALL   , mask , distance )

    deallocate( vector_ALL , distance_ALL , mask )

end if

!------------------------------------------------------------------------
! calculate dipole potential at a-b midpoint ...
!------------------------------------------------------------------------

allocate( mol_phi2( N_of_DP , 3 ) , source = D_zero )

forall( i=1:N_of_DP )  

!   first order ...
    decay(i,:)  =  vector(i,:) / (distance(i)*distance(i)*distance(i))

    ! DP_potential due to dipole i ...
    mol_phi(i)  = - dot_product( DP_Mols(i,:),decay(i,:) )

! second order ...
    mol_phi2(i,:) = 2.0d0 * DP_Mols(i,:) / (distance(i)*distance(i)*distance(i))

end forall

! excluding self-interaction ...
where( (nr_Mols == basis(a)%nr) .OR. (nr_Mols == basis(b)%nr) ) 
    mol_phi = 0.d0
    mol_phi2(:,1) = 0.d0
    mol_phi2(:,2) = 0.d0
    mol_phi2(:,3) = 0.d0
end where

! first order ...
DP_phi(1) = sum( mol_phi )

! second order ...
forall( j=1:3 ) DP_phi(j+1) = sum( mol_phi2(:,j) )

deallocate( vector , decay , distance , DP_Mols , nr_Mols , mol_phi , mol_phi2 )

end function DP_phi
!
!
!
!
!===========================================
 subroutine DeAllocate_DPs( DPs , n , flag )
!===========================================
implicit none
type(dipoles)            ,   intent(inout) :: DPs
integer       , optional ,   intent(in)    :: n
character(*)             ,   intent(in)    :: flag

select case( flag )

    case( "alloc" )

        allocate( DPs%CC ( n , 3 ) , source = D_zero )
        allocate( DPs%DP ( n , 3 ) , source = D_zero )
        allocate( DPs%nr ( n )     , source = I_zero )

    case( "dealloc" )

        deallocate( DPs%CC )
        deallocate( DPs%DP )
        deallocate( DPs%nr )

    case( "garbage" )

        ! garbage collection before restart calculations ...
        If( allocated( DPs%CC ) ) deallocate( DPs%CC )
        If( allocated( DPs%DP ) ) deallocate( DPs%DP )
        If( allocated( DPs%nr ) ) deallocate( DPs%nr )

end select

end subroutine DeAllocate_DPs
!
!
!
end module DP_potential_m
