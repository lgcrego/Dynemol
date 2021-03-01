module decoherence_m

    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use Structure_Builder       , only: Unit_Cell

    public :: apply_decoherence , DecoherenceRate

    private

    !module parameters ...

    !module variables ...

    interface apply_decoherence
        module procedure FirstOrderDecoherence
    end interface apply_decoherence

contains
!
!
!
!========================================================================
 subroutine FirstOrderDecoherence( MO_bra , MO_ket , erg , PES , t_rate )
!========================================================================
implicit none
complex*16 , intent(inout) :: MO_bra(:,:)
complex*16 , intent(inout) :: MO_ket(:,:)
real*8     , intent(in)    :: erg(:)
integer    , intent(in)    :: PES(:)
real*8     , intent(in)    :: t_rate

! local variables ...
integer :: n , i 
real*8  :: coeff , summ(2)
real*8 , allocatable :: tau_inv(:,:)

! J. Chem. Phys. 126, 134114 (2007)
CALL  DecoherenceRate( erg , PES , tau_inv )
! for the wavefunction tau(wvpckt) = 2.0*tau(rho) ...
tau_inv = tau_inv * HALF

summ = d_zero
do n = 1 , n_part 
do i = 1 , size(erg) 
     if( i == PES(n) ) cycle
     MO_bra(i,n) = MO_bra(i,n) * exp(-t_rate*tau_inv(i,n))
     MO_ket(i,n) = MO_ket(i,n) * exp(-t_rate*tau_inv(i,n))
     summ(n) = summ(n) + MO_bra(i,n)*MO_ket(i,n)
     end do
     end do

do n = 1 , n_part
     coeff = MO_bra(PES(n),n) * MO_ket(PES(n),n)
     coeff = (d_one - summ(n)) / coeff
     coeff = sqrt(coeff)
     MO_bra(PES(n),n) = MO_bra(PES(n),n) * coeff
     MO_ket(PES(n),n) = MO_ket(PES(n),n) * coeff
     end do

deallocate( tau_inv )

end subroutine FirstOrderDecoherence
!
!
!
!=======================================================
 subroutine DecoherenceRate( ergMO , PES , tau_inv )
!=======================================================
implicit none
real*8                , intent(in)  :: ergMO(:)
integer               , intent(in)  :: PES(:)
real*8  , allocatable , intent(out) :: tau_inv(:,:)

!local parameters ...
real*8 , parameter :: C = 0.1 * Hartree_2_eV   ! <== eV units
integer :: PPP(2) = [16,15]

!local variables ...                                                                                                                                                    
integer :: i , j
real*8  :: Const , tau , dE

! kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  !   E_kin()

allocate( tau_inv(size(ErgMO),2) , source = d_zero )

do j = 1 , n_part 
do i = 1 , size(ErgMO)
     if( i == PPP(j) ) cycle 
        dE = abs(ErgMO(i) - ErgMO(PPP(j)))
        tau_inv(i,j) = h_bar * Const / dE
        tau_inv(i,j) = d_one/tau_inv(i,j)
end do
end do

end subroutine DecoherenceRate
!
!
!
!================================
 function E_kin() result(kinetic)
!================================
use MD_read_m, only: MM, atom
implicit none

! local variables ...
integer :: i
real*8  :: kinetic

! calculation of the kinetic energy ...                                                                                                                                 
kinetic = d_zero
do i = 1 , MM % N_of_atoms
    kinetic = kinetic + atom(i)% mass * sum( atom(i)% vel(:) * atom(i)% vel(:) ) * half   ! <== J/kmol
end do
kinetic = kinetic * micro * kJmol_2_eV   ! <== eV

end function E_kin
!
!
!
end module decoherence_m
