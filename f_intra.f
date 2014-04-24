module F_intra_m
   
    use constants_m
    use for_force   , only : rcut, vrecut, frecut, pot, bdpot, angpot, dihpot   , &
                             vscut, fscut, KAPPA, lj14pot, coul14pot, pot2      , &    
                             Dihedral_Potential_Type                            
    use MD_read_m   , only : atom , molecule , MM , read_from_gmx
    use MM_types    , only : MM_system , MM_molecular , MM_atomic , debug_MM

    private

    public :: FORCEINTRA

    ! module variables ...
    real*8  , dimension (3) :: rij , rjk , rkl , rijk , rjkl , rijkl , f1 , f2 , f3 , f4
    real*8                  :: rijq , rjkq , rklq , rijsq , rjksq , rklsq , fxyz , riju , riku , rijkj , rijkj2 , rjkkl , rjkkl2 ,     &
                               rijkl2 , rjksq2 , rijkll , f1x , f1y , f1z , f2x , f2y , f2z , f3x , f3y , f3z , f4x , f4y , f4z ,      &
                               sr2 , sr6 , sr12 , fs , phi , cosphi , sinphi , rsinphi , coephi , gamma , KRIJ , expar , eme , dphi ,  &
                               term , chrgi , chrgj , freal , sig , eps , pterm , A0 , A1 , A2 , A3 , rtwopi , qterm , rterm , sterm , &
                               tterm , C0 , C1 , C2 , C3 , C4 , C5
    integer                 :: i , j , k , l , ati , atj , atk , atl , loop

contains
!
!
!
!====================
subroutine FORCEINTRA
!====================
implicit none

rtwopi = 1.d0/twopi

do j = 1 , MM % N_of_atoms
    atom(j) % fnonbd(:) = 0.d0           ! Non-bonded
    atom(j) % fbond(:)  = 0.d0           ! Stretching/Bonding 
    atom(j) % fang(:)   = 0.d0           ! Bending/Angular
    atom(j) % fdihed(:) = 0.d0           ! Dihedral
    atom(j) % fnonch(:) = 0.d0           ! Non-bonded coulomb 1-4
end do

! new stretch ...
do i = 1 , MM % N_of_molecules
    do j = 1 ,  molecule(i) % Nbonds
        ati = molecule(i) % bonds(j,1)
        atj = molecule(i) % bonds(j,2)
        if( atom(atj) % free .OR. atom(ati) % free ) then
            rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) )
            rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq  = SQRT(rijq)
            qterm  = 0.5d0 * molecule(i) % kbond0(j,1)*( rijsq - molecule(i) % kbond0(j,2) )**2
            coephi = molecule(i) % kbond0(j,1)*( rijsq - molecule(i) % kbond0(j,2) )/rijsq
            atom(atj) % fbond(:) = atom(atj) % fbond(:) - coephi*rij(:)
            atom(ati) % fbond(:) = atom(ati) % fbond(:) + coephi*rij(:)
            bdpot = qterm + bdpot
        end if 
    end do
end do

! new bend ...
do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Nangs
        atj = molecule(i) % angs(j,1)
        ati = molecule(i) % angs(j,2)
        atk = molecule(i) % angs(j,3)
        if( atom(atj) % free .OR. atom(ati) % free .OR. atom(atk) % free ) then
            rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) )
            rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq  = SQRT(rijq)
            rjk(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
            rjk(:) = rjk(:) - MM % box(:)*DNINT( rjk(:) * MM % ibox(:) )
            rjkq   = rjk(1)*rjk(1) + rjk(2)*rjk(2) + rjk(3)*rjk(3)
            rjksq  = SQRT(rjkq)

            phi = ACOS( (rij(1)*rjk(1) + rij(2)*rjk(2) + rij(3)*rjk(3) ) / ( rijsq * rjksq ) )

            coephi = 0.d0
            if (phi < 1.d-12 .OR. abs(pi - phi) < 1.d-12) then
                coephi = 0.d0
                rterm  = 0.0d0 
            else 
                coephi = ( phi - molecule(i) % kang0(j,2) ) / SIN(phi)
                rterm  = 0.5d0 * molecule(i) % kang0(j,1) * ( phi - molecule(i) % kang0(j,2) )**2
            end if
            angpot = rterm + angpot
         
            do l = 1, 3
                if (l == 1) atl = atj
                if (l == 2) atl = ati
                if (l == 3) atl = atk
                do loop = 1, 3                    !eixos X,Y,Z (n = 1, 2 ou 3)
                    fxyz = 0.d0
                    riju = rij(loop)
                    riku = rjk(loop)
                    fxyz = ( molecule(i) % kang0(j,1) * coephi ) *                     &
                           ( (DEL(atl,atj)-DEL(atl,ati))*riku/(rijsq*rjksq) +          &
                           (DEL(atl,atk)-DEL(atl,ati))*riju/(rijsq*rjksq) -            &
                           COS(phi)*( (DEL(atl,atj)-DEL(atl,ati))*riju/(rijsq*rijsq) + &
                           (DEL(atl,atk)-DEL(atl,ati))*riku/(rjksq*rjksq)) )
                    atom(atl) % fang(loop) = atom(atl) % fang(loop) + fxyz
                end do
            end do
        end if
    end do
end do

! ################################################################
! Dihedral Potential Angle ... 
do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Ndiheds
        ati = molecule(i) % diheds(j,1)
        atj = molecule(i) % diheds(j,2)
        atk = molecule(i) % diheds(j,3)
        atl = molecule(i) % diheds(j,4)
        if ( atom(atj) % free .OR. atom(ati) % free .OR. atom(atk) % free .OR. atom(atl) % free ) then
            ! Definition of vector rij = rj - ri
            rij(:)  = atom(ati) % xyz(1:3) - atom(atj) % xyz(1:3)
            rij(:)  = rij(1:3) - MM % box(1:3) * DNINT( rij(1:3) * MM % ibox(1:3) )
            ! Definition of vector rjk = rj - rk
            rjk(:)  = atom(atj) % xyz(1:3) - atom(atk) % xyz(1:3)
            rjk(:)  = rjk(1:3) - MM % box(1:3) * DNINT( rjk(1:3) * MM % ibox(1:3) )
            rjkq    = rjk(1)*rjk(1) + rjk(2)*rjk(2) + rjk(3)*rjk(3)
            rjksq   = 1.d0 / SQRT(rjkq)
            rjksq2  = rjksq * rjksq
            ! Definition of vector rkl = rl - rk
            rkl(:)  = atom(atk) % xyz(1:3) - atom(atl) % xyz(1:3)
            rkl(:)  = rkl(1:3) - MM % box(1:3) * DNINT( rkl(1:3) * MM % ibox(1:3) )
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
            if( read_from_gmx ) then
                CALL gmx
            else
                CALL not_gmx
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
 
! ####################################################################
! New non-bonded 1,4 intramolecular interactions ...
do i = 1 , MM % N_of_molecules
    do j   = 1 , molecule(i) % Nbonds14
        ati    = molecule(i) % bonds14(j,1)
        atj    = molecule(i) % bonds14(j,2)
        if ( atom(atj) % free .OR. atom(ati) % free ) then

            chrgi  = atom(ati) % charge
            chrgj  = atom(atj) % charge
            rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) )
            rklq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            ! Lennard Jones ...
            select case ( MM % CombinationRule )

                 case (2)

                     sig = atom(ati) % sig + atom(atj) % sig    ! AMBER FF :: GMX COMB-RULE 2

                 case (3)

                     sig = atom(ati) % sig * atom(atj) % sig    ! AMBER FF :: GMX COMB-RULE 3

            end select

            eps   =  atom(ati) % eps * atom(atj) % eps 
            rklsq = sqrt(rklq)
            sr2   = sig * sig / rklq
            sr6   = sr2 * sr2 * sr2
            sr12  = sr6 * sr6
            fs    = 24.d0 * eps * ( TWO * sr12 - sr6 )
            fs    = (fs / rklq) * molecule(i) % fact14(j)
            atom(ati) % fnonbd(1:3) = atom(ati) % fnonbd(1:3) + fs * rij(1:3)
            atom(atj) % fnonbd(1:3) = atom(atj) % fnonbd(1:3) - fs * rij(1:3)
            ! factor used to compensate the factor1 and factor2 factors ...
            ! factor3 = 1.0d-20
            sterm  = 4.d0 * eps * factor3 * ( sr12 - sr6 ) 
!           alternative cutoff formula ...
!           sterm  = sterm - vscut(ati,atj) + fscut(ati,atj) * ( rklsq - rcut ) 
            sterm  = sterm * molecule(i) % fact14(j)

            !  Real part (Numero de cargas igual ao numero de sitios)
            sr2   = 1.d0 / rklq
            KRIJ  = KAPPA * rklsq
            expar = EXP( -(KRIJ*KRIJ) )
            freal = coulomb * chrgi * chrgj * ( sr2/rklsq )
            freal = freal * ( ERFC(KRIJ) + TWO * rsqpi * KAPPA * rklsq * expar ) * molecule(i) % fact14(j)
            atom(ati) % fnonch(1:3) = atom(ati) % fnonch(1:3) + freal * rij(1:3)
            atom(atj) % fnonch(1:3) = atom(atj) % fnonch(1:3) - freal * rij(1:3)
            ! factor used to compensate the factor1 and factor2 factors ...
            ! factor3 = 1.0d-20
            tterm = coulomb*factor3 * chrgi * chrgj * ERFC(KRIJ)/rklsq * molecule(i) % fact14(j)
!           alternative cutoff formula ...
!           tterm = tterm - vrecut * chrgi * chrgj + frecut * chrgi * chrgj * ( rklsq-rcut ) * factor3
            lj14pot   = lj14pot + sterm
            coul14pot = coul14pot + tterm

        end if
    end do
end do

! factor used to compensate the factor1 and factor2 factors ...
! factor3 = 1.0d-20
pot2 = pot + ( bdpot + angpot + dihpot) * factor3 + lj14pot + coul14pot

! New Get total force ...

do i = 1 , MM % N_of_atoms
    atom(i) % ftotal(:) = atom(i) % ftotal(:) + ( atom(i) % fbond(:)  +  &
                                                  atom(i) % fang(:)   +  &
                                                  atom(i) % fdihed(:) +  &
                                                  atom(i) % fnonbd(:) +  & 
                                                  atom(i) % fnonch(:)    &
                                                  ) * 1.d-10
end do

end subroutine FORCEINTRA
!
!
!
!==============
 subroutine gmx
!==============
implicit none

! local variables ...
real*8  :: psi

select case( adjustl(molecule(i) % Dihedral_Type(j)) )
    case ('cos')    ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]        <== Eq. 4.61

        term  =   molecule(i) % harm(j) * phi - molecule(i) % kdihed0(j,1)
        pterm =   molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        gamma = - molecule(i) % kdihed0(j,2) * molecule(i) % harm(j) * sin(term) * rsinphi * rijkj * rjkkl

    case('cos3')    ! V = C0 + C1*cos(phi - 180) + C2*cos^2(phi - 180) + C3*cos^3(phi - 180) + C4*cos^4(phi - 180) + C5*cos(phi - 180)      <== Eq. 4.62

        psi = phi - PI

        pterm = molecule(i) % kdihed0(j,1)                                                        + &
                molecule(i) % kdihed0(j,2) * cos(psi)                                             + &
                molecule(i) % kdihed0(j,3) * cos(psi) * cos(psi)                                  + &
                molecule(i) % kdihed0(j,4) * cos(psi) * cos(psi) * cos(psi)                       + &
                molecule(i) % kdihed0(j,5) * cos(psi) * cos(psi) * cos(psi) * cos(psi)            + &
                molecule(i) % kdihed0(j,6) * cos(psi) * cos(psi) * cos(psi) * cos(psi) * cos(psi)

        gamma = - sin(psi) * ( molecule(i) % kdihed0(j,2) +                                                   &
                               2.0d0 * molecule(i) % kdihed0(j,3) * cos(psi) +                                &
                               3.0d0 * molecule(i) % kdihed0(j,4) * cos(psi) * cos(psi) +                     &
                               4.0d0 * molecule(i) % kdihed0(j,5) * cos(psi) * cos(psi) * cos(psi) +          &
                               5.0d0 * molecule(i) % kdihed0(j,6) * cos(psi) * cos(psi) * cos(psi) * cos(psi) ) * rsinphi * rijkj * rjkkl

end select

end subroutine gmx
!
!
!
!==================
 subroutine not_gmx
!==================
implicit none

!local variables ...

select case( adjustl(Dihedral_Potential_Type) )
    case ('cos') ! V = k[1 + cos(n.phi - theta)]
        term  = molecule(i) % Nharm * phi - molecule(i) % kdihed0(j,2)
        pterm = molecule(i) % kdihed0(j,1) * ( 1.d0 + cos(term) )
        gamma = - molecule(i) % kdihed0(j,1) * molecule(i) % Nharm * sin(term) * rsinphi * rijkj * rjkkl
          
    case('harm') ! V = 1/2.k(phi - phi0)²
        dphi  = phi - molecule(i) % kdihed0(j,2)
        dphi  = dphi - NINT(dphi*rtwopi) * rtwopi
        term  = molecule(i) % kdihed0(j,1) * dphi
        pterm = 0.5d0 * term * dphi
        gamma = term * rsinphi * rijkj * rjkkl          

    case('hcos') ! V = 1/2.k[cos(phi) - cos(phi0)]²
        dphi  = cos(phi) - cos( molecule(i) % kdihed0(j,2) )
        term  = molecule(i) % kdihed0(j,1) * dphi
        pterm = 0.5d0 * term * dphi
        gamma = term * rijkj * rjkkl          

    case('cos3') ! v = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
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
        gamma = - ( -C1 + 2.d0 * C2 * eme - 3.d0 * C3 * eme**2 + 4.d0 * C4 * eme**3 - 5.d0 * C5 * eme**4 ) * rijkj * rjkkl

    case('opls') ! V = A0 + 1/2{A1[1 + cos(phi)] + A2[1 - cos(2.phi)] + A3[1 + cos(3.phi)]}
        dphi  = phi - molecule(i) % kdihed0(j,5)
        A0    = molecule(i) % kdihed0(j,1) ; A1 = molecule(i) % kdihed0(j,2) ; A2 = molecule(i) % kdihed0(j,3) 
        A3    = molecule(i) % kdihed0(j,4)
        pterm = A0 + 0.5d0 * ( A1 * (1.d0 + cos(dphi)) + A2 * (1.d0 - cos(2.d0*dphi)) + A3 * (1.d0 + cos(3.d0*dphi)) )
        gamma = 0.5d0 * ( A1 * sin(dphi) - 2.d0 * A2 * sin(2.d0*dphi) + 3.d0 * A3 * sin(3.d0*dphi) ) * rsinphi * rijkj * rjkkl
          
    case('none')
        pterm = 0.d0
        gamma = 0.d0
  
end select

end subroutine not_gmx
!
!
!
!===================
 function ERFC ( X )
!===================
 real*8 :: ERFC
 real*8 :: A1, A2, A3, A4, A5, P, T, X, XSQ, TP
 parameter ( A1 = 0.254829592, A2 = -0.284496736 ) 
 parameter ( A3 = 1.421413741, A4 = -1.453122027 ) 
 parameter ( A5 = 1.061405429, P  =  0.3275911   ) 

 T    = 1.0 / ( 1.0 + P * X )
 XSQ  = X * X
 TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
 ERFC = TP * EXP ( -XSQ )

end function ERFC
!
!
!
!====================
 function DEL ( X,Y )
!====================
 real*8  :: DEL
 integer :: X, Y

 if (X == Y) then
   DEL = 1.
 else if (X /= Y) then
   DEL = 0.
 end if

end function DEL
!
!
end module F_intra_m



