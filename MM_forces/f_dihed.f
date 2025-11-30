module FF_diheds

    use type_m   
    use omp_lib
    use constants_m
    use MM_input     , only: MM_input_format
    use parameters_m , only: PBC 
    use for_force    , only: dihpot, ryck_dih, proper_dih, harm_dih, imp_dih 
    use MD_read_m    , only: atom , molecule , MM 

    private

    public :: f_dihed

    ! module variables
    real*8  :: gamma, pterm, rsinphi, phi, rijkj, rjkkl

contains
!
!
!=================
subroutine f_dihed
!=================
implicit none

! local_variables ...
real*8 , dimension (3):: rij, rjk, rkl, rijk, rjkl, rijkl
real*8  :: rjkq, rjksq, rijkj2, rjkkl2, rijkl2, rjksq2, rijkll
real*8  :: f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z
real*8  :: cosphi, sinphi, coephi
integer :: i, j, k, l, n, ati, atj, atk, atl

do j = 1 , MM % N_of_atoms
   atom(j) % fdihed(:) = D_zero  
end do
dihpot     = D_zero
proper_dih = D_zero
harm_dih   = D_zero
ryck_dih   = D_zero
imp_dih    = D_zero
 
!====================================================================
! Dihedral Potential Angle ... 
! USING IUPAC first convention for dihedral definitions: trans = 180 deg ...
do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Ndiheds
        ati = molecule(i) % diheds(j,1)
        atj = molecule(i) % diheds(j,2)
        atk = molecule(i) % diheds(j,3)
        atl = molecule(i) % diheds(j,4)
        if ( atom(atj) % flex .OR. atom(ati) % flex .OR. atom(atk) % flex .OR. atom(atl) % flex ) then
            ! Definition of vector rij = ri - rj
            rij(:)  = atom(ati) % xyz(1:3) - atom(atj) % xyz(1:3)
            rij(:)  = rij(1:3) - MM % box(1:3) * DNINT( rij(1:3) * MM % ibox(1:3) ) * PBC(1:3)
            ! Definition of vector rjk = rj - rk
            rjk(:)  = atom(atj) % xyz(1:3) - atom(atk) % xyz(1:3)
            rjk(:)  = rjk(1:3) - MM % box(1:3) * DNINT( rjk(1:3) * MM % ibox(1:3) ) * PBC(1:3)
            rjkq    = rjk(1)*rjk(1) + rjk(2)*rjk(2) + rjk(3)*rjk(3)
            rjksq   = 1.d0 / SQRT(rjkq)
            rjksq2  = rjksq * rjksq
            ! Definition of vector rkl = rk - rl
            rkl(:)  = atom(atk) % xyz(1:3) - atom(atl) % xyz(1:3)
            rkl(:)  = rkl(1:3) - MM % box(1:3) * DNINT( rkl(1:3) * MM % ibox(1:3) ) * PBC(1:3)
            ! Cross Product M = | rij X rjk | :: First dihedral vector ...
            rijk(1) = rij(2) * rjk(3) - rij(3) * rjk(2)
            rijk(2) = rij(3) * rjk(1) - rij(1) * rjk(3)
            rijk(3) = rij(1) * rjk(2) - rij(2) * rjk(1)
            rijkj   = rijk(1)*rijk(1) + rijk(2)*rijk(2) + rijk(3)*rijk(3)
            rijkj   = 1.d0 / SQRT(rijkj)
            rijkj2  = rijkj * rijkj
            ! Cross Product N = | rjk X rkl | :: Second dihedral vector ...
            rjkl(1) = rjk(2) * rkl(3) - rjk(3) * rkl(2)
            rjkl(2) = rjk(3) * rkl(1) - rjk(1) * rkl(3)
            rjkl(3) = rjk(1) * rkl(2) - rjk(2) * rkl(1)
            rjkkl   = rjkl(1)*rjkl(1) + rjkl(2)*rjkl(2) + rjkl(3)*rjkl(3)
            rjkkl   = 1.d0 / SQRT(rjkkl)
            rjkkl2  = rjkkl * rjkkl
            ! Cross Product O = | rjk X rkl | X | rij X rjk | 
            rijkl(1) = rjkl(2) * rijk(3) - rjkl(3) * rijk(2)
            rijkl(2) = rjkl(3) * rijk(1) - rjkl(1) * rijk(3)
            rijkl(3) = rjkl(1) * rijk(2) - rjkl(2) * rijk(1)
            rijkl2   = rijkl(1)*rijkl(1)+rijkl(2)*rijkl(2)+rijkl(3)*rijkl(3)
            rijkll   = SQRT(rijkl2)

            ! PHI is the dihedral angle defined by PHI = ACOS(B), where  B(rij,rjk,rkl) = [ (rij X rjk).(rjk X rkl) / |rij X rjk||rjk X rkl| ]
            ! and the sign of PHI is positive if the vector O is in the same direction as the bond vector rjk
            coephi = sum( rijk(:)*rjkl(:) )
            cosphi = coephi * rijkj * rjkkl 

            If( Abs(cosphi) > 1.0d0 ) cosphi=Sign(1.0d0,cosphi)
            sinphi = sum( rjk(:)*rijkl(:) ) * rjksq * rijkj * rjkkl  
            phi    = ATAN2( sinphi,cosphi )

            ! Avoid singularity in sinphi
            sinphi  = sign( max(1.d-10 , abs(sinphi) ) , sinphi )
            rsinphi = 1.d0 / sinphi
      
            ! selection of potential energy function type
            if( MM_input_format == "GMX" ) then
                CALL gmx(i,j)
            else
                CALL not_gmx(i,j)
            end if

            ! Calculate atomic forces ...
            f1x = gamma * ( (-rjkl(2) * rjk(3) + rjkl(3) * rjk(2)) - (-rijk(2) * rjk(3) + rijk(3) * rjk(2)) * coephi * rijkj2)
            f1y = gamma * ( ( rjkl(1) * rjk(3) - rjkl(3) * rjk(1)) - ( rijk(1) * rjk(3) - rijk(3) * rjk(1)) * coephi * rijkj2)
            f1z = gamma * ( (-rjkl(1) * rjk(2) + rjkl(2) * rjk(1)) - (-rijk(1) * rjk(2) + rijk(2) * rjk(1)) * coephi * rijkj2)
            f3x = gamma * ( (-rjkl(2) * rij(3) + rjkl(3) * rij(2)) - (-rijk(2) * rij(3) + rijk(3) * rij(2)) * coephi * rijkj2)
            f3y = gamma * ( ( rjkl(1) * rij(3) - rjkl(3) * rij(1)) - ( rijk(1) * rij(3) - rijk(3) * rij(1)) * coephi * rijkj2)
            f3z = gamma * ( (-rjkl(1) * rij(2) + rjkl(2) * rij(1)) - (-rijk(1) * rij(2) + rijk(2) * rij(1)) * coephi * rijkj2)
            f2x = gamma * ( (-rijk(2) * rkl(3) + rijk(3) * rkl(2)) - (-rjkl(2) * rkl(3) + rjkl(3) * rkl(2)) * coephi * rjkkl2)
            f2y = gamma * ( ( rijk(1) * rkl(3) - rijk(3) * rkl(1)) - ( rjkl(1) * rkl(3) - rjkl(3) * rkl(1)) * coephi * rjkkl2)
            f2z = gamma * ( (-rijk(1) * rkl(2) + rijk(2) * rkl(1)) - (-rjkl(1) * rkl(2) + rjkl(2) * rkl(1)) * coephi * rjkkl2)
            f4x = gamma * ( (-rijk(2) * rjk(3) + rijk(3) * rjk(2)) - (-rjkl(2) * rjk(3) + rjkl(3) * rjk(2)) * coephi * rjkkl2)
            f4y = gamma * ( ( rijk(1) * rjk(3) - rijk(3) * rjk(1)) - ( rjkl(1) * rjk(3) - rjkl(3) * rjk(1)) * coephi * rjkkl2)
            f4z = gamma * ( (-rijk(1) * rjk(2) + rijk(2) * rjk(1)) - (-rjkl(1) * rjk(2) + rjkl(2) * rjk(1)) * coephi * rjkkl2)

            atom(ati) % fdihed(1) = atom(ati) % fdihed(1) + f1x
            atom(ati) % fdihed(2) = atom(ati) % fdihed(2) + f1y
            atom(ati) % fdihed(3) = atom(ati) % fdihed(3) + f1z
            atom(atj) % fdihed(1) = atom(atj) % fdihed(1) - f1x - f3x + f2x
            atom(atj) % fdihed(2) = atom(atj) % fdihed(2) - f1y - f3y + f2y
            atom(atj) % fdihed(3) = atom(atj) % fdihed(3) - f1z - f3z + f2z
            atom(atk) % fdihed(1) = atom(atk) % fdihed(1) - f2x - f4x + f3x
            atom(atk) % fdihed(2) = atom(atk) % fdihed(2) - f2y - f4y + f3y
            atom(atk) % fdihed(3) = atom(atk) % fdihed(3) - f2z - f4z + f3z
            atom(atl) % fdihed(1) = atom(atl) % fdihed(1) + f4x
            atom(atl) % fdihed(2) = atom(atl) % fdihed(2) + f4y
            atom(atl) % fdihed(3) = atom(atl) % fdihed(3) + f4z

            dihpot = dihpot + pterm
    
        end if
    end do
end do

end subroutine f_dihed
!
!
!
!===================
 subroutine gmx(i,j)
!===================
implicit none
integer , intent(in) :: i, j

! local variables ...
real*8  :: psi, cos_Psi, dtheta
real*8  :: term, term1, term2, term3, term4 

select case( adjustl(molecule(i) % Dihedral_Type(j)) )

    case ('cos')    ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] 
                    ! Eq. 4.60 (GMX 5.0.5 manual)
        
        term  = int(molecule(i) % kdihed0(j,3)) * phi - molecule(i) % kdihed0(j,1)
        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        proper_dih = proper_dih + pterm
        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl

    case ('harm')   ! V = 1/2.k( xi - xi_0 )²
                    ! Eq. 4.59 (GMX 5.0.5 manual)

           dtheta = ( phi - molecule(i) % kdihed0(j,1) )
           dtheta = dtheta - Dnint( dtheta * 1.d0/TWOPI ) * TWOPI

           term  = molecule(i) % kdihed0(j,2) * dtheta 
           pterm = 0.5d0 * term * dtheta
           harm_dih = harm_dih + pterm 
           gamma = term * rsinphi * rijkj * rjkkl

    case('cos3')    ! V = C0 + C1*cos(phi - 180) + C2*cos^2(phi - 180) + C3*cos^3(phi - 180) + C4*cos^4(phi - 180) + C5*cos(phi - 180)  
                    ! Eq. 4.61 (GMX 5.0.5 manual)

        psi     = phi - PI
        cos_Psi = cos(psi)

        pterm = molecule(i) % kdihed0(j,1)                                                        + &
                molecule(i) % kdihed0(j,2) * cos_Psi                                              + &
                molecule(i) % kdihed0(j,3) * cos_Psi  * cos_Psi                                   + &
                molecule(i) % kdihed0(j,4) * cos_Psi  * cos_Psi  * cos_Psi                        + &
                molecule(i) % kdihed0(j,5) * cos_Psi  * cos_Psi  * cos_Psi  * cos_Psi             + &
                molecule(i) % kdihed0(j,6) * cos_Psi  * cos_Psi  * cos_Psi  * cos_Psi  * cos_Psi 

        ryck_dih = ryck_dih + pterm

        gamma = - sin(psi) * ( molecule(i) % kdihed0(j,2) +                                                   &
                               2.0d0 * molecule(i) % kdihed0(j,3) * cos_Psi  +                                &
                               3.0d0 * molecule(i) % kdihed0(j,4) * cos_Psi  * cos_Psi  +                     &
                               4.0d0 * molecule(i) % kdihed0(j,5) * cos_Psi  * cos_Psi  * cos_Psi  +          &
                               5.0d0 * molecule(i) % kdihed0(j,6) * cos_Psi  * cos_Psi  * cos_Psi  * cos_PSi  ) * rsinphi * rijkj * rjkkl

   case ('imp')    ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (improper) 
                   ! Eq. 4.60 (GMX 5.0.5 manual)

        term  = int(molecule(i) % kdihed0(j,3)) * phi - molecule(i) % kdihed0(j,1)
        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        imp_dih = imp_dih + pterm
        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl

    case ('chrm')   ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple) 
                    ! Eq. 4.60 (GMX 5.0.5 manual)

        term  = int(molecule(i) % kdihed0(j,3))  * phi - molecule(i) % kdihed0(j,1)
        term1 = int(molecule(i) % kdihed0(j,6))  * phi - molecule(i) % kdihed0(j,4)
        term2 = int(molecule(i) % kdihed0(j,9))  * phi - molecule(i) % kdihed0(j,7)
        term3 = int(molecule(i) % kdihed0(j,12)) * phi - molecule(i) % kdihed0(j,10)
        term4 = int(molecule(i) % kdihed0(j,15)) * phi - molecule(i) % kdihed0(j,13)

        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        pterm = pterm + molecule(i) % kdihed0(j,5)  * ( 1.d0 + cos(term1) )
        pterm = pterm + molecule(i) % kdihed0(j,8)  * ( 1.d0 + cos(term2) )
        pterm = pterm + molecule(i) % kdihed0(j,11) * ( 1.d0 + cos(term3) )
        pterm = pterm + molecule(i) % kdihed0(j,14) * ( 1.d0 + cos(term4) )

        proper_dih = proper_dih + pterm

        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,5)  * int( molecule(i) % kdihed0(j,6)) * sin(term1) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,8)  * int( molecule(i) % kdihed0(j,9)) * sin(term2) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,11) * int(molecule(i) % kdihed0(j,12)) * sin(term3) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,14) * int(molecule(i) % kdihed0(j,15)) * sin(term4) * rsinphi * rijkj * rjkkl

end select

end subroutine gmx
!
!
!
!=======================
 subroutine not_gmx(i,j)
!=======================
implicit none
integer , intent(in) :: i, j

! local variables ...
real*8 :: eme, dphi, rtwopi
real*8 :: A0, A1, A2, A3, C0, C1, C2, C3, C4, C5 
real*8 :: term, term1, term2, term3, term4 !, dphi1, dphi2 

rtwopi = 1.d0/twopi

select case( adjustl(molecule(i) % Dihedral_Type(j)) )

    case ('cos') ! V = k[1 + cos(n.phi - theta)]
        term  = int(molecule(i) % kdihed0(j,3)) * phi - molecule(i) % kdihed0(j,1)
        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        proper_dih = proper_dih + pterm
        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl
         
    case('harm') ! V = 1/2.k(phi - phi0)²
         stop "dead end"
!        dphi  = phi - molecule(i) % kdihed0(j,1)
!        dphi  = dphi - DNINT(dphi * rtwopi) * twopi
!        dphi1 = phi - molecule(i) % kdihed0(j,3)
!        dphi1 = dphi1 - DNINT(dphi1 * rtwopi) * twopi
!        dphi2 = phi - molecule(i) % kdihed0(j,5)
!        dphi2 = dphi2 - DNINT(dphi2 * rtwopi) * twopi
!        term  = molecule(i) % kdihed0(j,2) * dphi
!        term  = term + molecule(i) % kdihed0(j,4) * dphi1
!        term  = term + molecule(i) % kdihed0(j,6) * dphi2
!        pterm = HALF * term * dphi
!        pterm = pterm + HALF * term1 * dphi1
!        pterm = pterm + HALF * term2 * dphi2
!        harm_dih = harm_dih + pterm
!        gamma = term * rsinphi * rijkj * rjkkl          

    case('hcos') ! V = 1/2.k[cos(phi) - cos(phi0)]²
        dphi  = cos(phi) - cos( molecule(i) % kdihed0(j,2) )
        term  = molecule(i) % kdihed0(j,1) * dphi
        pterm = 0.5d0 * term * dphi
        gamma = term * rijkj * rjkkl          

    case('cos3') ! V = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
        pterm = 0.5d0 * ( molecule(i) % kdihed0(j,1) * ( 1.d0 + cos(phi)      ) + &
                          molecule(i) % kdihed0(j,2) * ( 1.d0 - cos(2.d0*phi) ) + &
                          molecule(i) % kdihed0(j,3) * ( 1.d0 + cos(3.d0*phi) ) )
        gamma = 0.5d0 * ( molecule(i) % kdihed0(j,1) * sin(phi)      - &
                2.0d0 * molecule(i) % kdihed0(j,2) * sin(2.d0*phi)   + &
                3.0d0 * molecule(i) % kdihed0(j,3) * sin(3.d0*phi) ) * rsinphi * rijkj * rjkkl 

    case('ryck') ! V = sum_i^5 Ci.[cos(phi)]^i
        eme   = cos(phi)
        C0    = molecule(i) % kdihed0(j,1) ; C1 = molecule(i) % kdihed0(j,2) ; C2 = molecule(i) % kdihed0(j,3) 
        C3    = molecule(i) % kdihed0(j,4) ; C4 = molecule(i) % kdihed0(j,5) ; C5 = molecule(i) % kdihed0(j,6) 
        pterm = C0 - C1 * eme + C2 * eme**2 - C3 * eme**3 + C4 * eme**4 - C5 * eme**5 
        ryck_dih = ryck_dih + pterm
        gamma = - ( -C1 + 2.d0 * C2 * eme - 3.d0 * C3 * eme**2 + 4.d0 * C4 * eme**3 - 5.d0 * C5 * eme**4 ) * rijkj * rjkkl

    case('opls') ! V = A0 + 1/2{A1[1 + cos(phi)] + A2[1 - cos(2.phi)] + A3[1 + cos(3.phi)]}
        dphi  = phi - molecule(i) % kdihed0(j,5)
        A0    = molecule(i) % kdihed0(j,1) ; A1 = molecule(i) % kdihed0(j,2) ; A2 = molecule(i) % kdihed0(j,3) 
        A3    = molecule(i) % kdihed0(j,4)
        pterm = A0 + 0.5d0 * ( A1 * (1.d0 + cos(dphi)) + A2 * (1.d0 - cos(2.d0*dphi)) + A3 * (1.d0 + cos(3.d0*dphi)) )
        gamma = 0.5d0 * ( A1 * sin(dphi) - 2.d0 * A2 * sin(2.d0*dphi) + 3.d0 * A3 * sin(3.d0*dphi) ) * rsinphi * rijkj * rjkkl
        

    case ('chrm')   ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple) 
                    ! Eq. 4.60 (GMX 5.0.5 manual)

        term  = int(molecule(i) % kdihed0(j,3))  * phi - molecule(i) % kdihed0(j,1)
        term1 = int(molecule(i) % kdihed0(j,6))  * phi - molecule(i) % kdihed0(j,4)
        term2 = int(molecule(i) % kdihed0(j,9))  * phi - molecule(i) % kdihed0(j,7)
        term3 = int(molecule(i) % kdihed0(j,12)) * phi - molecule(i) % kdihed0(j,10)
        term4 = int(molecule(i) % kdihed0(j,15)) * phi - molecule(i) % kdihed0(j,13)

        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        pterm = pterm + molecule(i) % kdihed0(j,5)  * ( 1.d0 + cos(term1) )
        pterm = pterm + molecule(i) % kdihed0(j,8)  * ( 1.d0 + cos(term2) )
        pterm = pterm + molecule(i) % kdihed0(j,11) * ( 1.d0 + cos(term3) )
        pterm = pterm + molecule(i) % kdihed0(j,14) * ( 1.d0 + cos(term4) )

        proper_dih = proper_dih + pterm

        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,5)  * int( molecule(i) % kdihed0(j,6)) * sin(term1) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,8)  * int( molecule(i) % kdihed0(j,9)) * sin(term2) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,11) * int(molecule(i) % kdihed0(j,12)) * sin(term3) * rsinphi * rijkj * rjkkl
        gamma = gamma - molecule(i) % kdihed0(j,14) * int(molecule(i) % kdihed0(j,15)) * sin(term4) * rsinphi * rijkj * rjkkl  

    case('none')
        pterm = 0.d0
        gamma = 0.d0
  
end select

end subroutine not_gmx
!                         
!                            
!
end module FF_diheds
