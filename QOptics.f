module QOptics_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Oscillator_m

    public :: Redfield_Equations

    private

 contains
!
!
!
subroutine Redfield_Equations(system, basis, DP_matrix_AO, zL, zR, erg)

type(structure) , intent(in)               :: system
type(STO_basis) , intent(in)               :: basis(:)
type(dipole)    , intent(in)               :: DP_matrix_AO(:,:)
complex*16      , intent(in) , allocatable :: zL(:,:) , zR(:,:)
real*8          , intent(in) , allocatable :: erg(:)

! . local variables
real*8  , parameter  :: Gamma_ij = 1.89808d-6  ! <== (omega_ij^3 e^2)/( 6 pi epsilon_0 c^3 hbar) ; unit = 1 / (eV^3 * ps * Angs^2) 

type(transition) :: Trans_DP

trans_DP%bra_range = electrons
trans_DP%ket_range = holes

CALL Transition_Dipole_Builder(system, basis, DP_matrix_AO, zL, zR, erg, Trans_DP)







end subroutine Redfield_Equations
!
!
!
!
end module QOptics_m
