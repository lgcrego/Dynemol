 module polarizability_m

    use type_m
    use constants_m
    use parameters_m            , only : file_type ,                &
                                         CH_and_DP_step
    use Babel_m                 , only : System_Characteristics
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : system => Extended_Cell    

    public :: Build_Induced_DP , Induced_DP_phi , Induced_DP

    private
    
    ! module variables ...       
    integer              , save :: counter = 0
    real*8 , allocatable , save :: Induced_DP(:,:) , Induced_DP_Dressed(:,:) , net_charge(:)

    ! module parameters ...
    real*8  , parameter :: c0 = 1.5d0               ,&   ! 3/2
                           c1 = 8.888888888888d-2   ,&   ! 4/45
                           c2 = 7.5d0               ,&   ! 15/2
                           c3 = 15.d0               ,&
                           c4 = 22.5d0              ,&   ! 45/2
                           c5 = 11.25d0             ,&   ! 45/4
                           c0_inv = 1.d0/c0

contains
!
!
!
!=========================================================================
 subroutine Build_Induced_DP( basis , Dual_bra , Dual_ket , t , instance )
!=========================================================================
implicit none
type(STO_basis)   , optional    , intent(in)    :: basis(:)
complex*16        , optional    , intent(in)    :: Dual_bra(:,:)
complex*16        , optional    , intent(in)    :: Dual_ket(:,:)
real*8            , optional    , intent(in)    :: t
character(*)      , optional    , intent(in)    :: instance

! local variables ...
integer                 :: ati , atj , N_of_atoms 
real*8                  :: distance , decay , charge_El , charge_Hl 
real*8                  :: vector(3)
real*8  , allocatable   :: local_E_field(:,:)  

! local parameters ...
real*8  , parameter     :: DP_factor = 4.803204d0               ! <== alpha[10^(-24)cm^3] * e / (4*Pi*epsilon_0) * 10^(20)  : Debye
real*8  , parameter     :: DP_potential_factor = 2.9979255d0    ! <== e*p[Debye]/(4*Pi*epsilon_0)                           : eV * Angs^2


N_of_atoms = system%atoms

! setup of Induced_DP ... 
If( instance == "allocate" ) then
    allocate( Induced_DP         (N_of_atoms,3) , source = D_zero )
    allocate( Induced_DP_Dressed (N_of_atoms,3) , source = D_zero )
    allocate( net_charge         (N_of_atoms)   , source = D_zero )
    return
end If

! calculate atomic net charges due to wavepacket ...
do ati = 1 , N_of_atoms

    charge_El = abs( sum( Dual_bra(:,1)*Dual_ket(:,1) , basis(:)%atom == ati ) )
    charge_Hl = abs( sum( Dual_bra(:,2)*Dual_ket(:,2) , basis(:)%atom == ati ) )

    net_charge(ati) = charge_HL - charge_EL

end do

! calculate induced atomic DP moment ...
allocate( local_E_field (N_of_atoms,3) , source = D_zero )

do ati = 1 , N_of_atoms

    do atj = 1 , N_of_atoms

        If( ati /= atj ) then

            vector(:) = system%coord(ati,:) - system%coord(atj,:)
        
            distance = sqrt( sum(vector(:)*vector(:)) )

            decay = D_one / (distance*distance*distance)

            local_E_field(ati,:) = local_E_field(ati,:) + net_charge(atj) * vector(:) * decay * exclusion(ati,distance)

        end If

    end do
    
    Induced_DP(ati,:) = system%polar(ati) * local_E_field(ati,:) * DP_factor  ! <== Debye

end do

If( mod(counter,CH_and_DP_step)==0 ) CALL visualize_Induced_DP ( t )

! NOTICE: dipole moment is multiplied by DP_potential_factor ...
Induced_DP_Dressed = Induced_DP * DP_potential_factor * half

counter = counter + 1

deallocate( local_E_field )

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
real*8                  :: hardcore 
real*8                  :: midpoint_ab(3)
real*8                  :: Induced_DP_phi(4)
real*8  , allocatable   :: distance(:) , distance_ALL(:) , mol_phi(:) 
real*8  , allocatable   :: vector(:,:) , vector_ALL(:,:) , decay(:,:) , DP_Mols(:,:) , mol_phi2(:,:)
logical , allocatable   :: mask(:)

! combination rule for solvation hardcore shell ...
hardcore = ( basis(a) % hardcore + basis(b) % hardcore ) / two

! midpoint between atoms a & b ...
midpoint_ab(1) = ( basis(a)%x + basis(b)%x ) / two
midpoint_ab(2) = ( basis(a)%y + basis(b)%y ) / two
midpoint_ab(3) = ( basis(a)%z + basis(b)%z ) / two

!------------------------------------------------------------------------
!------------------------------------------------------------------------

N_of_DP = size( Induced_DP_Dressed(:,1) )

allocate( vector_ALL   ( N_of_DP , 3 ) , source = D_zero  )
allocate( distance_ALL ( N_of_DP     ) , source = D_zero  )
allocate( mask         ( N_of_DP     ) , source = .false. )

do i = 1 , N_of_DP
    vector_ALL(i,:) = midpoint_ab(:) - system%coord(i,:)
    distance_ALL(i) = sqrt( sum( vector_ALL(i,:)*vector_ALL(i,:) ) )
end do

mask = ( distance_ALL > hardcore )

! keep only the DPs outside hardcore ...
N_of_DP = count( mask ) 

allocate( vector  ( N_of_DP , 3 ) , source = D_zero )
allocate( DP_Mols ( N_of_DP , 3 ) , source = D_zero )
allocate( decay   ( N_of_DP , 3 ) , source = D_zero )
allocate( distance( N_of_DP     ) , source = D_zero )
allocate( mol_phi ( N_of_DP     ) , source = D_zero )
allocate( id_atom ( N_of_DP     ) , source = I_zero )   

do j = 1 , 3 
    vector(:,j)  = pack( vector_ALL(:,j)         , mask , vector(:,j)  )
    DP_Mols(:,j) = pack( Induced_DP_Dressed(:,j) , mask , DP_Mols(:,j) )
end do

id_atom  = pack( [ ( i , i=1,size(Induced_DP_Dressed(:,1)) ) ] , mask , id_atom  )  
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
!=====================================
 subroutine visualize_Induced_DP ( t )
!=====================================
implicit none
real*8  , intent(in) :: t

! local variables ...
integer :: ati , i , j , N_of_DP
real*8  , allocatable   :: mod_p(:) 

N_of_DP = size(Induced_DP(:,1))

allocate( mod_p(N_of_DP) )

do ati = 1 , N_of_DP
    mod_p(ati) = sqrt( sum(Induced_DP(ati,:)*Induced_DP(ati,:)) )
end do

!saving Induced atomic dipole sizes ...
OPEN(unit=113 , file="tmp_data/DipoleSize.inpt" , status = "unknown", action = "write" , position = "append" )
do ati = 1 , N_of_DP
    write(113,'(F9.5)',advance='no') mod_p(ati) 
end do
close(113)

end subroutine visualize_Induced_DP
!
!
!
!=============================
 function exclusion( ati , r )
!=============================
implicit none
integer , intent(in) :: ati
real*8  , intent(in) :: r

! local variables ...
integer :: N_level
real*8  :: a, a2, a3, a4, a5, a6
real*8  :: r2, r3, r4, r5, r6
real*8  :: exclusion

! N quantum number of s orbital ...
N_level = atom( system%AtNo(ati) ) % NQuant(0)

! zeta parameter of s orbital ...
a = atom( system%AtNo(ati) ) % zeta(0,1)

r2 = r*r
a2 = a*a

select case ( N_level )

    case( 1 )

        exclusion = -two*a2*exp(-two*r*a)*( r2 + r/a + D_one/(two*a2) ) + D_one 

    case( 2 )

        r3 = r2*r
        r4 = r2*r2
        
        a3 = a2*a
        a4 = a2*a2

        exclusion = -c0_inv*a4*exp(-two*r*a) * (r4 + two*r3/a + three*r2/a2 + three*r/a3 + c0/a4) + D_one

    case( 3 )

        r3 = r2*r
        r4 = r2*r2
        r5 = r4*r
        r6 = r3*r3

        a3 = a2*a
        a4 = a2*a2
        a5 = a4*a
        a6 = a3*a3

        exclusion = -c1*a6*exp(-two*r*a) * (r6 + three*r5/a + c2*r4/a2 + c3*r3/a3 + c4*r2/a4 + c4*r/a5 + c5/a6) + D_one

end select 

end function exclusion
!
!
!
end module polarizability_m      
