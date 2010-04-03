module dipole_potential_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Semi_Empirical_Parms    , only : atom
    use DP_FMO_m                , only : DP_FMO_analysis
    use Multipole_Core          , only : Util_Multipoles

    public :: Solvent_Molecule_DP , DP_phi

    private

    integer , allocatable , save :: nr_Solvent_Mol(:)
    real*8  , allocatable , save :: DP_Solvent_Mol(:,:) 
    real*8  , allocatable , save :: Center_of_Charge(:,:) 

    real*8  , parameter   :: DP_potential_factor = 2.9979255d0 ! <== e*p(Debye)/(4*Pi*epsilon_0) : eV * Angs^2

contains
!
!
!
!===================================
 subroutine Solvent_Molecule_DP( a )
!===================================
implicit none
type(structure) , intent(inout) :: a

! local variables ...
integer              :: i , j , nr , last_nr , first_nr , i1 , i2 , nr_atoms
real*8               :: total_valence , norm , vector(3)
real*8 , allocatable :: Q_center(:,:) , DP_FMO(:,:) , Qi_Ri(:,:) 

! local parameter ...
real*8 , parameter   :: DP_value = 3.84d0   ! <== Debye

! initialize multipole routines ...
CALL Util_Multipoles

! garbage collection before restart calculations ...
If( allocated(Center_of_Charge) ) deallocate(Center_of_Charge)
If( allocated(DP_Solvent_Mol  ) ) deallocate(DP_Solvent_Mol  )
If( allocated(nr_Solvent_Mol  ) ) deallocate(nr_Solvent_Mol  )

! finding positions of solvent atoms in structure ...
first_nr = minval( a%nr , a%fragment == "S" )
last_nr  = maxval( a%nr , a%fragment == "S" )

! consistency check ...
If( a%N_of_Solvent_Molecules /= (last_nr - first_nr + 1) ) pause ">>> N_of_Solvent_Molecules /= 'S' count <<<"

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
    ! calculate Center_of_Charge for each solvent molecule: sum_i = (q_i * vec{r}_i) / sum_i q_i ...

    forall( j=1:3 , i=i1:i2 ) Qi_Ri(i,j) = atom(a%AtNo(i))%Nvalen * a%coord(i,j)

    total_valence = sum( [ (atom(a%AtNo(i))%Nvalen , i=i1,i2) ] )

    forall(j=1:3) Q_center(nr,j) = sum( Qi_Ri(i1:i2,j) ) / total_valence
    !-----------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------
    ! calculate the dipole vector ...

    CALL DP_FMO_analysis( a , Q_center(nr,:) , DP_FMO(nr,:) , nr ) 

! for mechanical calculation of DP_FMO ...    
!    forall( j=1:3 ) vector(j) = (a%coord(i1+1,j) - a%coord(i1,j))  
!    norm = sqrt( dot_product(vector,vector) )
!    DP_FMO(nr,:) = ( vector(:) / norm ) * DP_value

    !-----------------------------------------------------------------------------------------------

end do

allocate( Center_of_Charge  ( last_nr-first_nr+1 , 3 ) )
allocate( DP_Solvent_Mol    ( last_nr-first_nr+1 , 3 ) )
allocate( nr_Solvent_Mol    ( last_nr-first_nr+1 )     )

forall( j=1:3 ) 

    Center_of_Charge( : , j ) = Q_center( first_nr:last_nr , j ) 
    DP_Solvent_Mol  ( : , j ) = DP_FMO  ( first_nr:last_nr , j ) * DP_potential_factor

end forall

! save list of nr indices for use in DP_phi ...
nr_Solvent_Mol = [ ( i , i=first_nr,last_nr ) ]

deallocate( Q_center , DP_FMO , Qi_Ri )

end subroutine solvent_molecule_DP
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
real*8                  :: DP_phi
real*8                  :: midpoint_ab(3)
real*8  , allocatable   :: vector(:,:) , decay(:,:) , distance(:) , mol_phi(:)

! midpoint between atoms i & j ...
midpoint_ab(1) = ( basis(a)%x + basis(b)%x ) / two
midpoint_ab(2) = ( basis(a)%y + basis(b)%y ) / two
midpoint_ab(3) = ( basis(a)%z + basis(b)%z ) / two

N_of_DP = size( DP_Solvent_Mol(:,1) )

allocate( vector   ( N_of_DP ,3 ) )
allocate( decay    ( N_of_DP ,3 ) )
allocate( distance ( N_of_DP    ) )
allocate( mol_phi  ( N_of_DP    ) )

forall( j=1:3 ) vector(:,j) = midpoint_ab(j) - Center_of_Charge(:,j)

forall( i=1:N_of_DP ) 

    distance(i) =  sqrt( dot_product(vector(i,:),vector(i,:)) )
    decay(i,:)  =  vector(i,:) / (distance(i)*distance(i)*distance(i))

    ! DP_potential due to solvent molecule i ...
    mol_phi(i)  = - dot_product( DP_Solvent_Mol(i,:),decay(i,:) )

end forall

! excluding self-interaction ...
where( (nr_Solvent_Mol == basis(a)%nr) .AND. (nr_Solvent_Mol == basis(b)%nr) ) mol_phi = 0.d0

DP_phi = sum( mol_phi )

deallocate( vector , decay , distance , mol_phi )

end function DP_phi
!
!
end module dipole_potential_m

