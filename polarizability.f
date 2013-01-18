 module polarizability_m

    use type_m
    use constants_m
    use parameters_m            , only : file_type
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : system => Extended_Cell    

    public :: Build_Induced_DP , Induced_DP_phi

    
    ! module variables ...       
    real*8 , allocatable , save :: Induced_DP(:,:)

    
contains
!
!
!
!=====================================================================
 subroutine Build_Induced_DP( basis , Dual_bra , Dual_ket , instance )
!=====================================================================
implicit none
type(STO_basis)   , optional    , intent(in)    :: basis(:)
complex*16        , optional    , intent(in)    :: Dual_bra(:,:)
complex*16        , optional    , intent(in)    :: Dual_ket(:,:)
character(*)      , optional    , intent(in)    :: instance

! local variables ...
integer                 :: ati , atj , N_of_atoms
real*8                  :: distance , decay , charge_El , charge_Hl 
real*8                  :: vector(3)
real*8  , allocatable   :: net_charge(:) 
real*8  , allocatable   :: local_E_field(:,:)  

! local parameters ...
real*8  , parameter     :: DP_factor = 4.803204d0               ! <== alpha[10^(-24)cm^3] * e / (4*Pi*epsilon_0) * 10^(20)  : Debye
real*8  , parameter     :: DP_potential_factor = 2.9979255d0    ! <== e*p[Debye]/(4*Pi*epsilon_0)                           : eV * Angs^2

N_of_atoms = system%atoms

! setup of Induced_DP ... 
If( instance == "allocate" ) then
    allocate( Induced_DP (N_of_atoms,3) , source = D_zero )
    return
end If

! calculate atomic net charges due to wavepacket ...
allocate( net_charge(N_of_atoms) , source = D_zero )

do ati = 1 , N_of_atoms

    charge_El = abs( sum( Dual_bra(:,1)*Dual_ket(:,1) , basis(:)%atom == ati ) )
    charge_Hl = abs( sum( Dual_bra(:,2)*Dual_ket(:,2) , basis(:)%atom == ati ) )

    net_charge(ati) = charge_HL - charge_EL

end do

print*, " "
print*, " --> net_charge = " , sum(net_charge)

! calculate induced atomic DP moment ...
allocate( local_E_field (N_of_atoms,3) , source = D_zero )

    do ati = 1 , N_of_atoms

    do atj = 1 , N_of_atoms

        If( ati /= atj ) then

            vector(:) = system%coord(ati,:) - system%coord(atj,:)
        
            distance = sqrt( sum(vector(:)*vector(:)) )

            decay = D_one / (distance*distance*distance)

            If( distance < system%solvation_hardcore(ati) ) decay = 0.d0

            local_E_field(ati,:) = local_E_field(ati,:) + net_charge(atj) * vector(:) * decay

        end If

    end do
    
    Induced_DP(ati,:) = system%polar(ati) * local_E_field(ati,:) * DP_factor  ! <== Debye

end do

! NOTICE: dipole moment is multiplied by DP_potential_factor ...
Induced_DP = Induced_DP * DP_potential_factor * half

deallocate( net_charge , local_E_field )

end subroutine Build_Induced_DP
!
!
!
!=============================================
 pure function Induced_DP_phi( a , b , basis )
!=============================================
implicit none
integer         , intent(in) :: a , b
type(STO_basis) , intent(in) :: basis(:)

! local variables ...
integer                 :: i , j , N_of_DP 
integer , allocatable   :: id_atom(:)
real*8                  :: hard_core 
real*8                  :: midpoint_ab(3)
real*8                  :: Induced_DP_phi(4)
real*8  , allocatable   :: distance(:) , distance_ALL(:) , mol_phi(:) 
real*8  , allocatable   :: vector(:,:) , vector_ALL(:,:) , decay(:,:) , DP_Mols(:,:) , mol_phi2(:,:)
logical , allocatable   :: mask(:)

! combination rule for solvation hardcore shell ...
hard_core = ( basis(a)%solvation_hardcore + basis(b)%solvation_hardcore ) / two

! midpoint between atoms a & b ...
midpoint_ab(1) = ( basis(a)%x + basis(b)%x ) / two
midpoint_ab(2) = ( basis(a)%y + basis(b)%y ) / two
midpoint_ab(3) = ( basis(a)%z + basis(b)%z ) / two

!------------------------------------------------------------------------
!------------------------------------------------------------------------

N_of_DP = size( Induced_DP(:,1) )

allocate( vector_ALL   ( N_of_DP , 3 ) , source = D_zero  )
allocate( distance_ALL ( N_of_DP     ) , source = D_zero  )
allocate( mask         ( N_of_DP     ) , source = .false. )

do i = 1 , N_of_DP
    vector_ALL(i,:) = midpoint_ab(:) - system%coord(i,:)
    distance_ALL(i) = sqrt( sum( vector_ALL(i,:)*vector_ALL(i,:) ) )
end do

mask = ( distance_ALL > hard_core )

! keep only the DPs outside hardcore ...
N_of_DP = count( mask ) 

allocate( vector  ( N_of_DP , 3 ) , source = D_zero )
allocate( DP_Mols ( N_of_DP , 3 ) , source = D_zero )
allocate( decay   ( N_of_DP , 3 ) , source = D_zero )
allocate( distance( N_of_DP     ) , source = D_zero )
allocate( mol_phi ( N_of_DP     ) , source = D_zero )
allocate( id_atom ( N_of_DP     ) , source = I_zero )   

do j = 1 , 3 
    vector(:,j)  = pack( vector_ALL(:,j) , mask , vector(:,j)  )
    DP_Mols(:,j) = pack( Induced_DP(:,j) , mask , DP_Mols(:,j) )
end do

id_atom  = pack( [ ( i , i=1,size(Induced_DP(:,1)) ) ] , mask , id_atom  )  
distance = pack( distance_ALL   , mask , distance )

deallocate( vector_ALL , distance_ALL , mask )

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
where( (id_atom == basis(a)%atom) .OR. (id_atom == basis(b)%atom) ) 
    mol_phi = 0.d0
    mol_phi2(:,1) = 0.d0
    mol_phi2(:,2) = 0.d0
    mol_phi2(:,3) = 0.d0
end where

! first order ...
Induced_DP_phi(1) = sum( mol_phi )

! second order ...
forall( j=1:3 ) Induced_DP_phi(j+1) = sum( mol_phi2(:,j) )

deallocate( vector , decay , distance , DP_Mols , id_atom , mol_phi , mol_phi2 )

end function Induced_DP_phi
!
!
!
!
end module polarizability_m      
