module decoherence_m

    use type_m
    use constants_m
    use parameters_m            , only  : n_part
    use Structure_Builder       , only  : Unit_Cell

    public :: apply_decoherence 

    private

    !module parameters ...

    !module variables ...

contains
!
!
!
!====================================================================
 subroutine Apply_decoherence( MO_bra , MO_ket , erg , PES , t_rate )
!====================================================================
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

tau_inv = energy_based_decoherence( erg , PES )

summ = D_zero
do concurrent ( n=1:n_part , i=1:size(erg) , i /= PES(n) )

   MO_bra(i,n) = MO_bra(i,n) * exp(-t_rate*tau_inv(i,n))
   MO_ket(i,n) = MO_ket(i,n) * exp(-t_rate*tau_inv(i,n))

   summ(n) = summ(n) + MO_bra(i,n)*MO_ket(i,n)
   end do

do n = 1 , n_part

   coeff = MO_bra(PES(n),n) * MO_ket(PES(n),n)
   coeff = (d_one - summ(n)) / coeff
   coeff = sqrt(coeff)

   MO_bra(PES(n),n) = MO_bra(PES(n),n) * coeff
   MO_ket(PES(n),n) = MO_ket(PES(n),n) * coeff
   end do

deallocate( tau_inv )

end subroutine apply_decoherence
!
!
!
!==============================================================
 function energy_based_decoherence( erg , PES ) result(tau_inv)
!==============================================================
implicit none
real*8  , intent(in) :: erg(:)
integer , intent(in) :: PES(:)

!local parameters ...
real*8 , parameter :: C = 0.1 * Hartree_2_eV   ! <== eV units

!local variables ...
integer :: i , j
real*8  :: E_kin , tau
real*8 , allocatable :: tau_inv(:,:)

! kinetic energy in eV units ...
E_kin = Unit_Cell% MD_Kin + mid_prec

allocate( tau_inv(size(erg),2) , source = D_zero )

do concurrent ( j = 1:2 , i=1:size(erg) , i /= PES(j) )
   tau = h_bar/(erg(i)-erg(PES(j))) * (D_one + C/E_kin)
   tau_inv(i,j) = D_one/abs(tau)
   end do

end function energy_based_decoherence
!
!
!
end module decoherence_m
