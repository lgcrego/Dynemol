module Dielectric_Potential

    use type_m                  
    use constants_m
    use blas95
    use f95_precision
    use parameters_m            , only : PBC , EnvField_ , Environ_type , verbose
    use MD_read_m               , only : atom
    use DP_potential_m          , only : Dipole_Potentials

    public :: Environment_SetUp , Q_phi

    private

    ! module variables ...       
    type(molecular) , allocatable :: MolPBC(:)

    ! module parameters ...       
    real*8 , parameter :: units =  14.39965173d0          ! <== e^2/Angs = 14.399 eV
    real*8 , parameter :: refractive_index = 1.33d0       ! <== refractive index of the dielectric medium ...

contains
!
!
!
!===================================
 subroutine Environment_SetUp( sys )
!===================================
implicit none
type(structure) , intent(in)  :: sys

If( EnvField_ .and. (.not. any(sys%fragment=="S")) ) stop 'execution halted: did not define solvent fragment'

select case (Environ_Type)

    case( 'DP_QM' , 'DP_MM')  ! <== leave to modulo DP_potential_m

        CALL Dipole_Potentials( sys )

    case( 'Ch_MM' )           ! <== stay in the modulo

        CALL Classical_Point_Charges( sys )

    case default

        stop 'execution halted: check your option for Environ_Type in parameters.f'

end select

end subroutine Environment_SetUp
!
!
!
!=========================================
 subroutine Classical_Point_Charges( sys )
!=========================================
implicit none
type(structure) , intent(in)  :: sys

! local variables ...
integer                       :: i, j, I1, I2, nr, na, last_nr, first_nr, N_of_Q, N_of_Mols
real*8                        :: total_valence
real*8          , allocatable :: Qi_Ri(:,:) 
type(molecular) , allocatable :: Env_Mols(:)

! find positions of environment molecules ...
! pdb entries must be in a single block ...
first_nr = minval( sys%nr , sys%fragment == "S" )
last_nr  = maxval( sys%nr , sys%fragment == "S" )

! total number of molecules comprising the dielectric domain ...
N_of_Mols = last_nr - first_nr + 1

CALL Allocation( sys, Env_Mols , N_of_Mols )

! total number of point-charges in dielectric domain ...
N_of_Q = sum( Env_Mols(:)% N_of_Atoms )

! consistency checks ...
If( N_of_Q /= count( sys%nr >= first_nr .and. sys%nr <= last_nr ) ) stop '>>> something wrong in Environment_SetUp <<<'
If( verbose ) Print 157 , N_of_Mols

!======================================================
! Setting Up Env_Atoms and Env_Mols ...
!
allocate( Qi_Ri(sys%atoms ,3) , source=D_zero )

Env_Mols(:)% nr = [( i , i=first_nr,last_nr )]

do i = 1 , N_of_Mols

    nr = Env_Mols(i)% nr
    na = Env_Mols(i)% N_of_Atoms

    ! position of nr residue in variable sys[1:atoms] ...
    I1 = minloc( sys%nr , 1 , sys%nr == nr ) 
    I2 = (I1-1) + na

    !-----------------------------------------------------------------------------------------------
    ! Point Charges of environment molecules ...

                  Env_Mols(i)% PC% Q  (1:na) = atom(I1:I2)% MM_charge 

    forall(j=1:3) Env_Mols(i)% PC% xyz(1:na,j) = sys% coord(I1:I2,j) 

                  Env_Mols(i)% PC% nr (1:na)   = nr

    !-----------------------------------------------------------------------------------------------
    ! calculate Center_of_Charge for each environment molecule: sum_i = (q_i * vec{r}_i) / sum_i q_i ...

    forall( j=1:3 , i=I1:I2 ) Qi_Ri(i,j) = sys%Nvalen(i) * sys%coord(i,j) 

    total_valence = sum( sys%Nvalen(I1:I2) )

    forall(j=1:3) Env_Mols(i)% CC(j) = sum( Qi_Ri(I1:I2,j) ) / total_valence

    !-----------------------------------------------------------------------------------------------

end do

deallocate( Qi_Ri)
!======================================================

! generate periodic structure of dielectric domain ; if PBCx=PBCy=PBCz=0 ==> Q_atoms_pbc = Q_atoms ...
CALL give_me_PBC( sys, Env_Mols, MolPBC )

!do i = 1 , size(MolPBC)
!   do j = 1 , molpbc(i)% N_of_Atoms
!      if( MOLpbc(i) %pc% Q(j) <0. ) then
!          write(33,'(A4,3F9.4)')  "O" , molpbc(i)%pc%xyz(j,1) , molpbc(i)%pc%xyz(j,2) , molpbc(i)%pc%xyz(j,3)
!      else
!          write(33,'(A4,3F9.4)')  "H" , molpbc(i)%pc%xyz(j,1) , molpbc(i)%pc%xyz(j,2) , molpbc(i)%pc%xyz(j,3)
!      end if
!   end do
!end do
!
!do i = 1 , size(MolPBC)
!       write(34,'(A4,3F9.4)')  "I" , MolPBC(i)%CC(1) , MolPBC(i)%CC(2) , MolPBC(i)%CC(3)
!end do

include 'formats.h'

end subroutine Classical_Point_Charges
!
!
!
!=============================
 function Q_phi( sys , a , b )
!=============================
implicit none
type(structure) , intent(in) :: sys
integer         , intent(in) :: a , b

! local variables ...
integer               :: i , j , k , na , N_of_Q , N_of_M
real*8                :: hardcore , cut_off_radius, CC_distance, midpoint_ab(3) , Q_phi(4)
real*8  , allocatable :: distance(:), V_phi(:), V_phi2(:,:), Q(:), versor(:,:), vector_ALL(:,:) , AT_Q(:), AT_versor(:,:) , AT_distance(:) 
logical               :: inside

! combination rule for solvation hardcore shell ...
hardcore = ( sys%solvation_hardcore(a) + sys%solvation_hardcore(b) ) / TWO

! midpoint between atoms a & b ...
midpoint_ab(:) = ( sys% coord(a,:) + sys% coord(b,:) ) / TWO

! total number of point charges in PBC ...
N_of_Q = sum( MolPBC(:)% N_of_Atoms )

!------------------------------------------------------------------------
! choice of solvent droplet .vs. PBC box ...
!------------------------------------------------------------------------
N_of_M = size(MolPBC)

allocate( vector_ALL ( N_of_M , 3 ) , source = D_zero )
allocate( AT_versor  ( N_of_Q , 3 ) , source = D_zero )
allocate( AT_distance( N_of_Q     ) , source = D_zero )
allocate( AT_Q       ( N_of_Q     ) , source = D_zero )

If( sum(PBC) == 0) then

       do j = 1 , 3
           vector_ALL(:,j) = midpoint_ab(j) - MolPBC(:)% CC(j)
       end do
      
       k = 0 
       do i = 1 , size(MolPBC) 
       
           na = MolPBC(i)% N_of_Atoms
       
           do j = 1 , na
       
                k = k + 1
       
                AT_Q(k) = MolPBC(i)% PC% Q(j)
          
                AT_distance(k) = sqrt( sum((midpoint_ab(:) - MolPBC(i)% PC% xyz(j,:))**2) )
       
                AT_versor(k,:) = ( midpoint_ab(:) - MolPBC(i)% PC% xyz(j,:) ) / AT_distance(k)
       
           end do 

        end do

        N_of_Q = k

else

        ! maximum distance from midpoint a-b ...
        cut_off_radius = minval(sys% T_xyz) / TWO
        
        do j = 1 , 3
            vector_ALL(:,j) = midpoint_ab(j) - MolPBC(:)% CC(j)
        end do

        k = 0 
        do i = 1 , size(MolPBC)
        
            CC_distance = sqrt( dot_product( vector_ALL(i,:),vector_ALL(i,:) ) )

            inside = ( (CC_distance > hardcore) .AND. (CC_distance < cut_off_radius) )
        
            If( inside ) then

                na = MolPBC(i)% N_of_Atoms

                do j = 1 , na
        
                   AT_Q(k+j) = MolPBC(i)% PC% Q(j)
           
                   AT_distance(k+j) = sqrt( sum((midpoint_ab(:) - MolPBC(i)% PC% xyz(j,:))**2) )

                   AT_versor(k+j,:) = ( midpoint_ab(:) - MolPBC(i)% PC% xyz(j,:) ) / AT_distance(k+j)
        
                end do 
                k = k + na
            end If 
        end do

        N_of_Q = k

end If

allocate( Q       ( N_of_Q     ) , source = AT_Q       (1:N_of_Q  ) )
allocate( distance( N_of_Q     ) , source = AT_distance(1:N_of_Q  ) )
allocate( versor  ( N_of_Q , 3 ) , source = AT_versor  (1:N_of_Q,:) )

deallocate( AT_Q , AT_distance , AT_versor )

!------------------------------------------------------------------------
! calculate dipole potential at a-b midpoint ...
!------------------------------------------------------------------------

allocate( V_phi   ( N_of_Q     ) , source = D_zero )
allocate( V_phi2  ( N_of_Q , 3 ) , source = D_zero )

! zeroth order potential due to point charges i ...
V_phi(:) = Q(:)/distance(:) 

do j = 1 , 3

    ! first order ...
    V_phi2(:,j) = Q(:) * versor(:,j) / (distance(:)*distance(:))

end do

! eliminate self-interactions ...
!where( (indx == a) .OR. (indx == b) ) 
!    V_phi = 0.d0
!    V_phi2(:,1) = 0.d0
!    V_phi2(:,2) = 0.d0
!    V_phi2(:,3) = 0.d0
!end where

! first order ...
Q_phi(1) = sum( V_phi(:) )

! second order ...
forall( j=1:3 ) Q_phi(j+1) = sum( V_phi2(:,j) )

! applying optical dielectric screening ; fix sign problem ...
Q_phi = - Q_phi * units / (refractive_index)**2

deallocate( versor , distance , Q , V_phi , V_phi2 )

end function Q_phi
!
!
!
!============================================
 subroutine give_me_PBC( sys , Mol , MolPBC )
!============================================
 implicit none
 type(structure) , intent(in)                  :: sys
 type(molecular) , allocatable , intent(inout) :: Mol(:) 
 type(molecular) , allocatable , intent(out)   :: MolPBC(:) 

! local variables ... 
integer  :: i , ix, iy, iz, j, L, n, na, N_of_M, N_of_M_pbc, nr_max

N_of_M = size(Mol)

! (VIRTUAL) REPLICAS for Period Boundary Conditions ...
N_of_M_pbc = product(2*PBC(:)+1) * N_of_M

If( .not. allocated(MolPBC) ) allocate( MolPBC(N_of_M_pbc) )

! original cell ...
forall(j=1:3) MolPBC(1:N_of_M)% CC(j)      = Mol(:)% CC(j)
              MolPBC(1:N_of_M)% nr         = Mol(:)% nr
              MolPBC(1:N_of_M)% N_of_Atoms = Mol(:)% N_of_Atoms

do i = 1 , N_of_M
    na = Mol(i)% N_of_Atoms
    If( .not. allocated(MolPBC(i)% PC% Q) ) Then
        allocate( MolPBC(i)% PC% Q(na)     )
        allocate( MolPBC(i)% PC% nr(na)    )
        allocate( MolPBC(i)% PC% xyz(na,3) )
    end If
    MolPBC(i)% PC% Q   = Mol(i)% PC% Q  
    MolPBC(i)% PC% nr  = Mol(i)% PC% nr 
    MolPBC(i)% PC% xyz = Mol(i)% PC% xyz
    end do

nr_max = Mol(N_of_M)%nr 

! including the replicas        
L = N_of_M
J = 0

DO iz = -PBC(3) , PBC(3)
DO iy = -PBC(2) , PBC(2)
DO ix = -PBC(1) , PBC(1)

    If( (ix /= 0) .OR. (iy /= 0) .OR. (iz /= 0) ) THEN

        J = J + 1

        DO n = 1 , N_of_M

            L = L + 1

            na = MolPBC(n)% N_of_Atoms

            MolPBC(L)% N_of_Atoms = na

            MolPBC(L)% CC(1) = Mol(n)% CC(1) + ix * sys% T_xyz(1)
            MolPBC(L)% CC(2) = Mol(n)% CC(2) + iy * sys% T_xyz(2)
            MolPBC(L)% CC(3) = Mol(n)% CC(3) + iz * sys% T_xyz(3)
            MolPBC(L)% nr    = Mol(n)% nr + nr_max*J

            If( .not. allocated(MolPBC(L)% PC% Q) ) Then
                allocate( MolPBC(L)% PC% Q(na)     )
                allocate( MolPBC(L)% PC% nr(na)    )
                allocate( MolPBC(L)% PC% xyz(na,3) )
            end If

            MolPBC(L)% PC% xyz(1:na,1) = Mol(n)% PC% xyz(1:na,1) + ix * sys% T_xyz(1)
            MolPBC(L)% PC% xyz(1:na,2) = Mol(n)% PC% xyz(1:na,2) + iy * sys% T_xyz(2)
            MolPBC(L)% PC% xyz(1:na,3) = Mol(n)% PC% xyz(1:na,3) + iz * sys% T_xyz(3)
            MolPBC(L)% PC% Q  (1:na)   = Mol(n)% PC% Q  (1:na)
            MolPBC(L)% PC% nr (1:na)   = MolPBC(L)% nr

        END DO


    END IF

END DO
END DO
END DO

! don't need these anymore ...
deallocate( Mol )

end subroutine give_me_PBC
!
!
!
!===================================
 subroutine Allocation( sys, a , n )
!===================================
implicit none
type(structure)               , intent(in)    :: sys
type(molecular) , allocatable , intent(inout) :: a(:)
integer                       , intent(in)    :: n

! local variables ...
integer :: i , I1, I2, nr_atoms, nr, first_nr, last_nr

allocate( a(n) )

! find positions of environment molecules ...
first_nr = minval( sys%nr , sys%fragment == "S" )
last_nr  = maxval( sys%nr , sys%fragment == "S" )

i = 0
do nr = first_nr , last_nr 

    i = i + 1

    ! # of atoms with tag nr ...
    nr_atoms = count( sys%nr == nr ) 

    a(i)% N_of_atoms = nr_atoms

    ! position of nr residue in variable sys[1:atoms] ...
    I1 = minloc( sys%nr , 1 , sys%nr == nr ) 
    I2 = (I1-1) + nr_atoms 

    allocate( a(i)% PC% Q  (nr_atoms)   )
    allocate( a(i)% PC% nr (nr_atoms)   )
    allocate( a(i)% PC% xyz(nr_atoms,3) )

end do
   
end subroutine Allocation
!
!
!
end module Dielectric_Potential
