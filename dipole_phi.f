module dipole_potential_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Structure_Builder       , only : Extended_Cell
    use Semi_Empirical_Parms    , only : atom
    use DP_FMO_m                , only : DP_FMO_analysis
    use Multipole_Core          , only : Util_Multipoles
    use PBC_m                   , only : Generate_Periodic_Solvent

    public :: Solvent_Molecule_DP , DP_phi

    private

! module variables ...       
    integer , allocatable , save :: pbc_nr_Solvent_Mol(:)
    real*8  , allocatable , save :: pbc_DP_Solvent_Mol(:,:) 
    real*8  , allocatable , save :: pbc_Center_of_Charge(:,:)  

! module parameters ...       
    real*8  , parameter   :: DP_potential_factor = 2.9979255d0  ! <== e*p(Debye)/(4*Pi*epsilon_0) : eV * Angs^2

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
integer                 :: i , j , nr , last_nr , first_nr , i1 , i2 , nr_atoms
integer , allocatable   :: nr_Solvent_Mol(:)
real*8                  :: total_valence , norm , vector(3)
real*8  , allocatable   :: Q_center(:,:) , DP_FMO(:,:) , Qi_Ri(:,:) 
real*8  , allocatable   :: Center_of_Charge(:,:) , DP_Solvent_Mol(:,:)   

! local parameter ...
real*8 , parameter :: DP_value = 3.84d0   ! <== Debye

! initialize multipole routines ...
CALL Util_Multipoles

! garbage collection before restart calculations ...
If( allocated(pbc_Center_of_Charge) ) deallocate(pbc_Center_of_Charge)
If( allocated(pbc_DP_Solvent_Mol  ) ) deallocate(pbc_DP_Solvent_Mol  )
If( allocated(pbc_nr_Solvent_Mol  ) ) deallocate(pbc_nr_Solvent_Mol  )

! finding positions of solvent atoms in structure ...
first_nr = minval( a%nr , a%fragment == "S" )
last_nr  = maxval( a%nr , a%fragment == "S" )

a%N_of_Solvent_Molecules = (last_nr - first_nr + 1)

! consistency check ...
Print 157 , (last_nr - first_nr + 1)

allocate(Q_center (last_nr,3) , source = 0.d0)
allocate(DP_FMO   (last_nr,3) , source = 0.d0)      ! <== dipole moments of solvent molecules
allocate(Qi_Ri    (a%atoms,3) , source = 0.d0)     

do nr = first_nr , last_nr 

    ! No of atoms with tag nr ...
    nr_atoms = count( a%nr == nr ) 

    ! position of nr residue in structure ...
    i1 = minloc( a%nr , 1 , a%nr == nr ) 
    i2 = (i1-1) + nr_atoms 

    !-----------------------------------------------------------------------------------------------
    ! calculate Center_of_Charge for each solvent molecule: sum_i = (q_i * vec{r}_i) / sum_i q_i ...

    forall( j=1:3 , i=i1:i2 ) Qi_Ri(i,j) = a%Nvalen(i) * a%coord(i,j)

    total_valence = sum( [ (a%Nvalen(i) , i=i1,i2) ] )

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

do j = 1 , 3 
    Center_of_Charge( : , j ) = Q_center( first_nr:last_nr , j ) 
    DP_Solvent_Mol  ( : , j ) = DP_FMO  ( first_nr:last_nr , j ) * DP_potential_factor
end do

! save list of nr indices for use in DP_phi ...
nr_Solvent_Mol = [ ( i , i=first_nr,last_nr ) ]

! generate period structure of solvent molecules ; if mmx=mmy=mmz=0 ==> (pbc_DP_Solvent_Mol = DP_Solvent_Mol) ...
CALL Generate_Periodic_Solvent( a                       , &
                                Center_of_Charge        , &
                                DP_Solvent_Mol          , &
                                nr_Solvent_Mol          , &
                                pbc_Center_of_Charge    , &
                                pbc_DP_Solvent_Mol      , &
                                pbc_nr_Solvent_Mol      )

deallocate( Q_center , Center_of_Charge , DP_FMO , DP_Solvent_Mol , Qi_Ri , nr_Solvent_Mol )

include 'formats.h'

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
integer , allocatable   :: nr_Solvent_Mol(:)
real*8                  :: DP_phi , cut_off_radius
real*8                  :: midpoint_ab(3)
real*8  , allocatable   :: vector(:,:) , vector_ALL(:,:) , decay(:,:) , distance(:) , distance_ALL(:) , mol_phi(:) 
real*8  , allocatable   :: DP_Solvent_Mol(:,:) 
logical , allocatable   :: mask(:)

! local parameters ...
real*8  , parameter     :: hard_core = 2.d0  ! <== (Angs)

! midpoint between atoms a & b ...
midpoint_ab(1) = ( basis(a)%x + basis(b)%x ) / two
midpoint_ab(2) = ( basis(a)%y + basis(b)%y ) / two
midpoint_ab(3) = ( basis(a)%z + basis(b)%z ) / two

! total number of dipoles in PBC solvent ...
N_of_DP = size( pbc_DP_Solvent_Mol(:,1) )

!------------------------------------------------------------------------
! choice of solvent droplet .vs. PBC box ...
!------------------------------------------------------------------------

If( mmx+mmy+mmz == 0) then

    allocate( vector         ( N_of_DP , 3 ) , source = D_zero )
    allocate( DP_Solvent_Mol ( N_of_DP , 3 ) , source = D_zero )
    allocate( decay          ( N_of_DP , 3 ) , source = D_zero )
    allocate( distance       ( N_of_DP     ) , source = D_zero )
    allocate( mol_phi        ( N_of_DP     ) , source = D_zero )
    allocate( nr_Solvent_Mol ( N_of_DP     ) , source = I_zero )

    do i = 1 , N_of_DP

        vector(i,:) = midpoint_ab(:) - pbc_Center_of_Charge(i,:)
        distance(i) = sqrt( sum( vector(i,:)*vector(i,:) ) )

    end do

    DP_Solvent_Mol = pbc_DP_Solvent_Mol
    nr_Solvent_Mol = pbc_nr_Solvent_Mol

else

    ! maximum distance from midpoint a-b ...
    cut_off_radius = minval(Extended_Cell%T_xyz) * two / three

    allocate( vector_ALL   ( N_of_DP , 3 ) , source = D_zero  )
    allocate( distance_ALL ( N_of_DP     ) , source = D_zero  )
    allocate( mask         ( N_of_DP     ) , source = .false. )

    do i = 1 , N_of_DP

        vector_ALL(i,:) = midpoint_ab(:) - pbc_Center_of_Charge(i,:)
        distance_ALL(i) = sqrt( sum( vector_ALL(i,:)*vector_ALL(i,:) ) )

    end do

    mask = ( (distance_ALL > hard_core) .AND. (distance_ALL < cut_off_radius) )

    ! redefine range of interaction with solvent molecules inside cut_off radius ...
    N_of_DP = count( mask ) 

    allocate( vector         ( N_of_DP , 3 ) , source = D_zero )
    allocate( DP_Solvent_Mol ( N_of_DP , 3 ) , source = D_zero )
    allocate( decay          ( N_of_DP , 3 ) , source = D_zero )
    allocate( distance       ( N_of_DP     ) , source = D_zero )
    allocate( mol_phi        ( N_of_DP     ) , source = D_zero )
    allocate( nr_Solvent_Mol ( N_of_DP     ) , source = I_zero )

    do j = 1 , 3 
        vector(:,j)         = pack( vector_ALL(:,j)         , mask , vector(:,j)         )
        DP_Solvent_Mol(:,j) = pack( pbc_DP_Solvent_Mol(:,j) , mask , DP_Solvent_Mol(:,j) )
    end do

    nr_Solvent_Mol = pack( pbc_nr_Solvent_Mol , mask , nr_Solvent_Mol )
    distance       = pack( distance_ALL       , mask , distance       )

    deallocate( vector_ALL , distance_ALL , mask )

end if

!------------------------------------------------------------------------
! calculate dipole potential at a-b midpoint ...
!------------------------------------------------------------------------

forall( i=1:N_of_DP )  

    decay(i,:)  =  vector(i,:) / (distance(i)*distance(i)*distance(i))

    ! DP_potential due to solvent molecule i ...
    mol_phi(i)  = - dot_product( DP_Solvent_Mol(i,:),decay(i,:) )

end forall

! excluding self-interaction ...
where( (nr_Solvent_Mol == basis(a)%nr) .OR. (nr_Solvent_Mol == basis(b)%nr) ) mol_phi = 0.d0

DP_phi = sum( mol_phi )

deallocate( vector , decay , distance , DP_Solvent_Mol , nr_Solvent_Mol , mol_phi )

end function DP_phi
!
!
end module dipole_potential_m
