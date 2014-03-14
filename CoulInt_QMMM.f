module QMMM_m
   
    use constants_m
    use for_force   , only : KAPPA
    use MD_read_m   , only : atom , MM

    private

    public :: QMMM_FORCE

contains
!
!
!
!=================================
subroutine QMMM_FORCE( NetCharge )
!=================================
implicit none
real*8  , intent(in)    :: NetCharge(:)

! local variables ...
real*8  , dimension (3) :: rij 
real*8                  :: rklq , rklsq , KRIJ , freal , expar
integer                 :: i , j

do j = 1 , size(atom)
    atom(j) % fcoupling = D_zero
end do

do i = 1 , size(atom)
    do j = i + 1 , size(atom)
    
        rij(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
        rij(:) = rij(:) - MM % box(:) * DINT( rij(:) * MM % ibox(:) )

        rklq = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
        rklsq = sqrt( rklq )

        KRIJ  = KAPPA * rklsq
        expar = exp( - KRIJ * KRIJ )

        freal = coulomb * NetCharge(i) * NetCharge(j) * ( 1.0d0 / rklq ) * ( 1.0d0 / rklsq )
        freal = freal * ( erfc( KRIJ ) + two * rsqpi * KAPPA * rklsq * expar )
        
        atom(i) % fcoupling(1:3) = atom(i) % fcoupling(1:3) + freal * rij(1:3)
        atom(j) % fcoupling(1:3) = atom(j) % fcoupling(1:3) - freal * rij(1:3)

    end do
end do

! New Get total force ...
do i = 1 , MM % N_of_atoms
    atom(i) % ftotal(:) = atom(i) % ftotal(:) + atom(i) % fcoupling(:) * 1.d-10
end do

end subroutine QMMM_FORCE
!
!
!
!==================
 function erfc( X )
!==================
implicit none
real*8              :: erfc
real*8              :: T , X , XSQ , TP
real*8  , parameter :: A1 =   0.254829592d0
real*8  , parameter :: A2 = - 0.284496736d0
real*8  , parameter :: A3 =   1.421413741d0
real*8  , parameter :: A4 = - 1.453122027d0
real*8  , parameter :: A5 =   1.061405429d0
real*8  , parameter :: P  =   0.327591100d0

T    = 1.0 / ( 1.0 + P * X )
XSQ  = X * X
TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
erfc = TP * exp( - XSQ )

end function erfc
!
!
!
end module QMMM_m
