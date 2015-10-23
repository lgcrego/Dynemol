module QMMM_m
   
    use constants_m
    use for_force            , only : pot_total
    use MD_read_m            , only : atom , MM
    use polarizability_m     , only : Induced_DP
    use Semi_Empirical_Parms , only : chemical_element => atom

    public :: QMMM_FORCE

    private

    ! module parameters ...
    real*8  , parameter :: c0 = 1.5d0               ,&   ! 3/2
                           c1 = 8.888888888888d-2   ,&   ! 4/45
                           c2 = 7.5d0               ,&   ! 15/2
                           c3 = 15.d0               ,&
                           c4 = 22.5d0              ,&   ! 45/2
                           c5 = 11.25d0             ,&   ! 45/4
                           c6 = 6.349206d-3         ,&   ! 2/315
                           c7 = 2.821869d-4         ,&   ! 4/14175
                           c8 = 85511197d-6         ,&   ! 4/467775
                           c0_inv = 1.d0/c0
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
real*8                  :: V_ChgChg , V_ChgIndDP , NetChrg_i , NetChrg_j , TotalChrg_i, ChargeProduct_ij
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

        NetChrg_i = NetCharge(i) * exclusion( i , rijsq*HALF )
        NetChrg_j = NetCharge(j) * exclusion( j , rijsq*HALF )

        ChargeProduct_ij = NetChrg_i * NetChrg_j 
        ChargeProduct_ij = ChargeProduct_ij + atom(i)%charge*NetChrg_j + atom(j)%charge*NetChrg_i

        freal = coulomb * ChargeProduct_ij * ( D_ONE / rijq  ) * ( D_ONE / rijsq )

        atom(i) % fcoupling(1:3) = atom(i) % fcoupling(1:3) + freal * rij(1:3)
        atom(j) % fcoupling(1:3) = atom(j) % fcoupling(1:3) - freal * rij(1:3)

        V_ChgChg = V_ChgChg + coulomb * ChargeProduct_ij * ( D_ONE / rijsq )

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

            TotalChrg_i = atom(i)%charge + NetCharge(i)

            f(1:3) = coulomb * TotalChrg_i * a(1:3) * ( D_ONE / rijq ) * ( D_ONE / rijsq ) 

            atom(i)% fcoupling(1:3) = atom(i)% fcoupling(1:3) + f(1:3)
            atom(j)% fcoupling(1:3) = atom(j)% fcoupling(1:3) - f(1:3)

            V_ChgIndDP = V_ChgIndDP &
                       + coulomb*D_2_eAngs*TotalChrg_i*dot_product(Induced_DP(j,1:3),rij(1:3))*( D_ONE / rijq )*( D_ONE / rijsq )

        end if
    end do
end do

! Append total force with Excited State Coulombic terms; force units = J/mts = Newtons ...
do i = 1 , MM % N_of_atoms
    atom(i)% ftotal(:) = atom(i)% ftotal(:) + atom(i)% fcoupling(:) * Angs_2_mts
end do

pot_total = pot_total + (V_ChgChg + V_ChgIndDP ) * mol * micro * factor3 / MM%N_of_molecules 

end subroutine QMMM_FORCE
!
!
!
!===========================
 function exclusion( i , r )
!===========================
implicit none
integer , intent(in) :: i
real*8  , intent(in) :: r

! local variables ...
integer :: N_level
real*8  :: a, a2, a3, a4, a5, a6, a7, a9, a11
real*8  :: r2, r3, r4, r5, r6, r7, r9, r11
real*8  :: exclusion

! N quantum number of s orbital ...
N_level = chemical_element( atom(i)%AtNo ) % NQuant(0)

! zeta parameter of s orbital ...
a = chemical_element( atom(i)%AtNo ) % zeta(0,1)

r2 = r*r
a2 = a*a

select case ( N_level )

    case( 1 )

        exclusion = -two*a2*exp(-two*r*a)*( r2 + r/a + D_one/(two*a2) ) + D_one 

    case( 2 )

        r3 = r2*r ; r4 = r2*r2
        
        a3 = a2*a ; a4 = a2*a2

        exclusion = -c0_inv*a4*exp(-two*r*a) * (r4 + two*r3/a + three*r2/a2 + three*r/a3 + c0/a4) + D_one

    case( 3 )

        r3 = r2*r ; r4 = r2*r2 ; r5 = r4*r ; r6 = r3*r3

        a3 = a2*a ; a4 = a2*a2 ; a5 = a4*a ; a6 = a3*a3

        exclusion = -c1*a6*exp(-two*r*a) * (r6 + three*r5/a + c2*r4/a2 + c3*r3/a3 + c4*r2/a4 + c4*r/a5 + c5/a6) + D_one

    case( 4 )

        r3 = r2*r ; r4 = r2*r2 ; r5 = r4*r ; r6 = r3*r3 ; r7 = r3*r4 

        a3 = a2*a ; a4 = a2*a2 ; a5 = a4*a ; a6 = a3*a3 ; a7 = a3*a4 

        exclusion = -c1*a6*exp(-two*r*a) * (r6 + three*r5/a + c2*r4/a2 + c3*r3/a3 + c4*r2/a4 + c4*r/a5 + c5/a6) + D_one &
                    -c6*r7*a7*exp(-two*r*a) * (four + a*r)

    case( 5 ) 

        r3 = r2*r ; r4 = r2*r2 ; r5 = r4*r ; r6 = r3*r3 ; r7 = r3*r4 ; r9 = r4*r5

        a3 = a2*a ; a4 = a2*a2 ; a5 = a4*a ; a6 = a3*a3 ; a7 = a3*a4 ; a9 = a4*a5

        exclusion = -c1*a6*exp(-two*r*a) * (r6 + three*r5/a + c2*r4/a2 + c3*r3/a3 + c4*r2/a4 + c4*r/a5 + c5/a6) + D_one &
                    -( c6*r7*a7*(four + a*r) + c7*r9*a9*(five + r*a) ) * exp(-two*r*a)

    case( 6 ) 

        r3 = r2*r ; r4 = r2*r2 ; r5 = r4*r ; r6 = r3*r3 ; r7 = r3*r4 ; r9 = r4*r5 ; r11 = r9*r2

        a3 = a2*a ; a4 = a2*a2 ; a5 = a4*a ; a6 = a3*a3 ; a7 = a3*a4 ; a9 = a4*a5 ; a11 = a9*a2

        exclusion = -c1*a6*exp(-two*r*a) * (r6 + three*r5/a + c2*r4/a2 + c3*r3/a3 + c4*r2/a4 + c4*r/a5 + c5/a6) + D_one &
                    -( c6*r7*a7*(four + a*r) + c7*r9*a9*(five + r*a) + c8*r11*a11*(six + r*a) ) * exp(-two*r*a)

end select 

end function exclusion
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
