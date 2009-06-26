module QOptics_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Oscillator_m
    use FMO_m           , only : orbital
    use Multipole_Core

    public :: Redfield_Equations

    private

 contains
!
!
!
!-------------------------------------------------------------
subroutine Redfield_Equations(system, basis, QM, FMO)
!-------------------------------------------------------------
 type(structure) , intent(in)               :: system
 type(STO_basis) , intent(in)               :: basis(:)
 type(eigen)     , intent(in)               :: QM
 type(eigen)     , intent(in)               :: FMO

! . local variables
 complex*16            :: trace
 real*8  , parameter   :: Gamma_ij = 1.89808d-6  ! <== (omega_ij^3 e^2)/( 6 pi epsilon_0 c^3 hbar) ; unit = 1 / (eV^3 * ps * Angs^2) 
 real*8  , allocatable :: rho(:,:)
 real*8                :: Sparcity
 type(transition)      :: Trans_DP
 integer               :: All_basis, Nonzero, dim_bra, dim_ket
 real*8 , allocatable  :: bra(:) , ket(:)





 All_basis = size(basis)


!-----------------------------------------------------------
 trans_DP%bra_indx_range = electrons
 trans_DP%ket_indx_range = holes
 trans_DP%flag = 'Redfield'

 CALL Transition_Dipole_Builder(system, basis, QM, Trans_DP)

! . sets up rho(0)[i,j] = <i|L> <L|j> = ket[i] X bra[j]
 dim_bra = size(trans_DP%bra_POINTER)
 dim_ket = size(trans_DP%ket_POINTER)
 Allocate( bra(dim_bra) , ket(dim_ket) )
 bra = FMO%L(trans_DP%bra_POINTER,orbital(1))
 ket = FMO%R(trans_DP%ket_POINTER,orbital(1))
 print*, sum(ket*bra)



 print*, size(bra), size(ket), size(trans_DP%matrix(:,1)%DP(1)) , size(trans_DP%matrix(1,:)%DP(1))

 do i = 1 , dim_bra
 do j = 1 , dim_ket
    write(13,*) trans_DP%matrix(i,j)%DP(1) , trans_DP%matrix(j,i)%DP(1)
 end do
 end do





end subroutine Redfield_Equations
!
!
!
!
end module QOptics_m
