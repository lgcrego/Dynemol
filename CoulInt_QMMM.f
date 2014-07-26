module QMMM_m
   
    use constants_m
    use MD_read_m           , only : atom , MM
    use polarizability_m    , only : Induced_DP

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
real*8  , dimension (3) :: rij , a(3) , f(3)
real*8                  :: rijq , rijsq , freal 
real*8                  :: V_ChgChg , V_ChgIndDP
integer                 :: i , j

! local parameters ...
real*8  :: D_2_eAngs = 0.20819434d0
real*8  :: cutoff    = 4.0d0

do j = 1 , size(atom)
    atom(j) % fcoupling = D_zero
end do

V_ChgChg   = D_zero
V_ChgIndDP = D_zero

do i = 1 , size(atom)

    ! charge/charge interaction ...
    do j = i + 1 , size(atom)
    
        rij(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
        rij(:) = rij(:) - MM % box(:) * DINT( rij(:) * MM % ibox(:) )

        rijq = sum( rij(:) * rij(:) )
        rijsq = sqrt( rijq )

        freal = coulomb * NetCharge(i) * NetCharge(j) * ( D_ONE / rijq  ) * ( D_ONE / rijsq )
        
        atom(i) % fcoupling(1:3) = atom(i) % fcoupling(1:3) + freal * rij(1:3)
        atom(j) % fcoupling(1:3) = atom(j) % fcoupling(1:3) - freal * rij(1:3)

        V_ChgChg = V_ChgChg + coulomb * NetCharge(i) * NetCharge(j) * ( D_ONE / rijsq )

    end do

    ! charge/induced-dipole interaction ...
    do j = 1 , size(atom)
    
        rij(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
        rij(:) = rij(:) - MM % box(:) * DINT( rij(:) * MM % ibox(:) )

        rijq = sum( rij(:) * rij(:) )
        rijsq = sqrt( rijq )

        if( rijsq > cutoff ) then

            ! a = (p.r)r ...
            a(1:3) = dot_product(Induced_DP(j,1:3),rij(1:3)) * rij(1:3)
            ! a = 3*(p.r)*r / ( |r|^2 ) - p ...
            a(1:3) = THREE * a(1:3) * ( D_ONE / rijq ) - Induced_DP(j,1:3)

            f(1:3) = coulomb * NetCharge(i) * a(1:3) * ( D_ONE / rijq ) * ( D_ONE / rijsq ) 

            atom(i) % fcoupling(1:3) = atom(i) % fcoupling(1:3) + f(1:3)
            atom(j) % fcoupling(1:3) = atom(j) % fcoupling(1:3) - f(1:3)

            V_ChgIndDP = V_ChgIndDP &
                       + coulomb*D_2_eAngs*NetCharge(i)*dot_product(Induced_DP(j,1:3),rij(1:3))*( D_ONE / rijq )*( D_ONE / rijsq )

        end if
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

T    = D_ONE / ( D_ONE + P * X )
XSQ  = X * X
TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
erfc = TP * exp( - XSQ )

end function erfc
!
!
!
end module QMMM_m
