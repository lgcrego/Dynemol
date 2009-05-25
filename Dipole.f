module Dipole_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use EHT_parameters
    use Structure_Builder 

    public :: Dipole_Moment , Center_of_Charge

    private

contains
!
!
!
!
subroutine Dipole_Moment(system, basis, M_matrix, L_vec, R_vec)

type(structure) , intent(inout) :: system
type(STO_basis) , intent(in)    :: basis(:)
type(dipole)    , intent(in)    :: M_matrix(:,:)
complex*16      , intent(in)    :: L_vec(:,:) , R_vec(:,:)

integer                         :: i, j, states, xyz, n_basis, Fermi_state
real*8                          :: Nuclear_DP(3), Electronic_DP(3), Total_DP(3) 
real*8          , allocatable   :: R_vector(:,:)
complex*16      , allocatable   :: a(:,:), b(:,:)
complex*16      , allocatable   :: origin_Dependent_DP(:), origin_Independent_DP(:)

real*8          , parameter     :: Debye_unit = 4.803204d0
complex*16      , parameter     :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

! atomic positions measured from the Center of Charge
 allocate(R_vector(system%atoms,3))
 CALL Center_of_Charge(system)
 forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - system%Center_of_Charge(xyz)

! Nuclear dipole ; if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)
 forall(xyz=1:3) Nuclear_DP(xyz) = sum( atom( system%AtNo(:) )%Nvalen * R_vector(:,xyz) )

! Electronic dipole 
 n_basis      =  size(basis)
 Fermi_state  =  system%N_of_electrons / 2
 
 allocate( a(n_basis,n_basis) )
 allocate( b(n_basis,n_basis) )
 allocate( origin_Dependent_DP(Fermi_state) )
 allocate( origin_Independent_DP(Fermi_state) )

 do xyz = 1 , 3

!   origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

    forall(states=1:Fermi_state)

        forall(i=1:n_basis) a(states,i) = L_vec(states,i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent_DP(states) = 2.d0 * sum( a(states,:) * R_vec(:,states) )

    end forall    
 
!   origin independent DP = sum{C_dagger * vec{M_matrix(i,j)} * C}

    a = M_matrix(:,:)%dp(xyz)
       
    CALL gemm(L_vec,a,b,'N','N',one,zero)    

    a = transpose(L_vec)

    forall(states=1:Fermi_state) origin_Independent_DP(states) = 2.d0 * sum(b(states,:)*a(:,states))

    Electronic_DP(xyz) = sum(origin_Dependent_DP + origin_Independent_DP)
 
    Total_DP(xyz) = Nuclear_DP(xyz) - Electronic_DP(xyz) 

 end do

 Total_DP = Total_DP * Debye_unit

 Print 154, Total_DP, dsqrt(sum(Total_DP(1:3)*Total_DP(1:3)))

 deallocate(R_vector,a,b)
 deallocate(origin_Dependent_DP)
 deallocate(origin_Independent_DP)

 include 'formats.h'

end subroutine Dipole_Moment
!
!
!
!
subroutine Center_of_Charge(a)

type(structure) , intent(inout) :: a

real*8 , allocatable :: Qi_Ri(:,:) 
real*8               :: total_valence

 allocate(Qi_Ri(a%atoms,3))

 forall(j=1:3,i=1:a%atoms) Qi_Ri(i,j) = atom(a%AtNo(i))%Nvalen * a%coord(i,j)

 total_valence = sum(atom(a%AtNo)%Nvalen)

 forall(j=1:3) a%Center_of_Charge(j) = sum(Qi_Ri(:,j)) / total_valence

 deallocate(Qi_Ri)

end subroutine Center_of_Charge
!
!
end module Dipole_m
