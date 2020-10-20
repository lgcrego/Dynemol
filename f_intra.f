module F_intra_m
   
    use constants_m
    use MM_input     , only : MM_input_format
    use parameters_m , only : PBC , QMMM
    use setup_m      , only : offset
    use for_force    , only : rcut, vrecut, frecut, pot_INTER, bdpot, angpot, dihpot, Morspot,      &
                              vscut, fscut, KAPPA, LJ_14, LJ_intra, Coul_14, Coul_intra, pot_total, &    
                              Dihedral_Potential_Type, forcefield, rcutsq, ryck_dih, proper_dih,    &
                              harm_dih, imp_dih, harm_bond, morse_bond
    use MD_read_m    , only : atom , molecule , MM 
    use MM_types     , only : MM_system , MM_molecular , MM_atomic , debug_MM
    use gmx2mdflex   , only : SpecialPairs , SpecialPairs14 , SpecialMorse

    private

    public :: FORCEINTRA, pot_INTRA

    ! module variables ...
    real*8  , dimension (3) :: rij , rjk , rkl , rik , rijk , rjkl , rijkl , f1 , f2 , f3 , f4
    real*8                  :: rijq , rjkq , rklq , rijsq , rjksq , rklsq , fxyz , riju , riku , rijkj , rijkj2 , rjkkl , rjkkl2 ,     &
                               rijkl2 , rjksq2 , rijkll , f1x , f1y , f1z , f2x , f2y , f2z , f3x , f3y , f3z , f4x , f4y , f4z ,      &
                               sr2 , sr6 , sr12 , fs , phi , cosphi , sinphi , rsinphi , coephi , gamma , KRIJ , expar , eme , dphi ,  &
                               term , chrgi , chrgj , freal , sig , eps , pterm , A0 , A1 , A2 , A3 , rtwopi , qterm , qterm0 , rterm, &
                               sterm , tterm , C0 , C1 , C2 , C3 , C4 , C5 , coephi0 , rterm0 , rikq , riksq , term1 , term2 , term3 , &
                               term4 , dphi1 , dphi2
    real*8                  :: pot_INTRA
    integer                 :: i , j , k , l , m , n , ati , atj , atk , atl , loop , ati1 , atj1 
    logical                 :: flag1, flag2, flag3, flag4, flag5
    logical                  :: there_are_NB_SpecialPairs   = .false.
    logical                  :: there_are_NB_SpecialPairs14 = .false.

    integer , allocatable   :: species_offset(:)


contains
!
!
!
!====================
subroutine FORCEINTRA
!====================
implicit none

! local_variables ...

rtwopi = 1.d0/twopi

do j = 1 , MM % N_of_atoms
    atom(j) % fnonbd14(:) = D_zero         ! Non-bonded 1-4
    atom(j) % fnonbd(:)   = D_zero         ! Non-bonded Intramolecular
    atom(j) % fbond(:)    = D_zero         ! Stretching/Bonding 
    atom(j) % fang(:)     = D_zero         ! Bending/Angular
    atom(j) % fdihed(:)   = D_zero         ! Dihedral
    atom(j) % fnonch14(:) = D_zero         ! Non-bonded coulomb 1-4
    atom(j) % fnonch(:)   = D_zero         ! Non-bonded coulomb Intramolecular
    atom(j) % fMorse(:)   = D_zero         ! Non-bonded Morse
end do

bdpot      = D_zero
angpot     = D_zero
dihpot     = D_zero
Morspot    = D_zero
harm_bond  = D_zero
morse_bond = D_zero
proper_dih = D_zero
harm_dih   = D_zero
ryck_dih   = D_zero
imp_dih    = D_zero
LJ_14      = D_zero
LJ_intra   = D_zero
Coul_14    = D_zero
Coul_intra = D_zero

CALL offset( species_offset )

! Bonding - stretching potential ...
do i = 1 , MM % N_of_molecules
    do j = 1 ,  molecule(i) % Nbonds
        ati = molecule(i) % bonds(j,1)
        atj = molecule(i) % bonds(j,2)
        if( atom(atj) % flex .OR. atom(ati) % flex ) then 

            rij(:)  = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:)  = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rijq    = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq   = SQRT(rijq)

            select case ( molecule(i) % bond_type(j) )

                case ( "harm" )
                ! harmonic potential ...
                qterm   = HALF * molecule(i) % kbond0(j,1) * ( rijsq - molecule(i) % kbond0(j,2) ) * ( rijsq - molecule(i) % kbond0(j,2) ) 
                coephi  = molecule(i) % kbond0(j,1)*( rijsq - molecule(i) % kbond0(j,2) )/rijsq
                harm_bond = qterm + harm_bond

                case ( "Mors" )
                ! Morse potential ...
                qterm0 = exp( -molecule(i) % kbond0(j,3) * ( rijsq - molecule(i) % kbond0(j,2) ) ) 
                qterm  = molecule(i) % kbond0(j,1) * ( D_ONE - qterm0 )*( D_ONE - qterm0 )
                coephi = TWO * molecule(i) % kbond0(j,1) * molecule(i) % kbond0(j,3) * qterm0 * ( D_ONE - qterm0 ) / rijsq 
                morse_bond = qterm + morse_bond               
 
            end select

            atom(atj) % fbond(:) = atom(atj) % fbond(:) - coephi*rij(:)
            atom(ati) % fbond(:) = atom(ati) % fbond(:) + coephi*rij(:)  
            bdpot = qterm + bdpot

        end if 
    end do
end do 

!====================================================================
! Angle - bending potential ...
do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Nangs
        atj = molecule(i) % angs(j,1)
        ati = molecule(i) % angs(j,2)
        atk = molecule(i) % angs(j,3)
        if( atom(atj) % flex .OR. atom(ati) % flex .OR. atom(atk) % flex ) then
            rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq  = SQRT(rijq)
            rjk(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
            rjk(:) = rjk(:) - MM % box(:)*DNINT( rjk(:) * MM % ibox(:) ) * PBC(:)
            rjkq   = rjk(1)*rjk(1) + rjk(2)*rjk(2) + rjk(3)*rjk(3)
            rjksq  = SQRT(rjkq)

            phi = ACOS( (rij(1)*rjk(1) + rij(2)*rjk(2) + rij(3)*rjk(3) ) / ( rijsq * rjksq ) )

            select case ( molecule(i) % angle_type(j) )
           
                case( "harm" , "urba" )
                ! Harmonic and Urey-Bradley potentials ...
 
                coephi = 0.d0
                if (phi < 1.d-12 .OR. abs(pi - phi) < 1.d-12) then
                    coephi = 0.d0
                    rterm  = 0.0d0 
                else 
                    coephi = ( phi - molecule(i) % kang0(j,2) ) / SIN(phi)
                    rterm  = HALF * molecule(i) % kang0(j,1) * ( phi - molecule(i) % kang0(j,2) )**2
                end if
                angpot = rterm + angpot
         
                do l = 1, 3
                    if (l == 1) atl = atj
                    if (l == 2) atl = ati
                    if (l == 3) atl = atk
                    do loop = 1, 3                    ! X,Y,Z axis (n = 1, 2 or 3)
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

                ! Urey-Bradley bonding term ...
                if( molecule(i) % angle_type(j) == "urba" ) then
                    rik(:) = atom(atk) % xyz(:) - atom(atj) % xyz(:)
                    rik(:) = rik(:) - MM % box(:) * DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
                    rikq   = rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3)
                    riksq  = SQRT(rikq)

                    coephi0 = molecule(i) % kang0(j,3) * ( riksq - molecule(i) % kang0(j,4) )/riksq
                    rterm0  = HALF * molecule(i) % kang0(j,3) * & 
                       ( riksq - molecule(i) % kang0(j,4) ) * ( riksq - molecule(i) % kang0(j,4) ) 
                    atom(atk) % fang(:) = atom(atk) % fang(:) - coephi0*rik(:)
                    atom(atj) % fang(:) = atom(atj) % fang(:) + coephi0*rik(:)  
                    angpot = rterm0 + angpot 
                end if
              
            end select 
        end if
    end do
end do

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

!====================================================================
! Non-bonded 1,4 intramolecular interactions ...

If( allocated(SpecialPairs14) ) there_are_NB_SpecialPairs14 = .true.

do i = 1 , MM % N_of_molecules
    do j   = 1 , molecule(i) % Nbonds14
        ati    = molecule(i) % bonds14(j,1)
        atj    = molecule(i) % bonds14(j,2)

        if ( atom(atj) % flex .OR. atom(ati) % flex ) then

            chrgi  = atom(ati) % charge
            chrgj  = atom(atj) % charge
            rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rklq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

            ! Lennard Jones ...
            select case ( MM % CombinationRule )

                case (2)
                    ! AMBER FF :: GMX COMB-RULE 2
                    sr2 = ( ( atom(ati) % sig14 + atom(atj) % sig14  ) * ( atom(ati) % sig14 + atom(atj) % sig14 ) ) / rklq

                case (3)
                    ! OPLS  FF :: GMX COMB-RULE 3
                    sr2 = ( ( atom(ati) % sig14 * atom(atj) % sig14 ) * ( atom(ati) % sig14 * atom(atj) % sig14 ) ) / rklq

            end select
            eps   =  atom(ati) % eps14 * atom(atj) % eps14 


            If( there_are_NB_SpecialPairs14 ) then    ! <== check whether (I,J) is a SpecialPair ... 

               read_loop1: do  n = 1, size(SpecialPairs14)

                   flag1 = ( adjustl( SpecialPairs14(n) % MMSymbols(1) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
                           ( adjustl( SpecialPairs14(n) % MMSymbols(2) ) == adjustl( atom(atj) % MMSymbol ) )
                   flag2 = ( adjustl( SpecialPairs14(n) % MMSymbols(2) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
                           ( adjustl( SpecialPairs14(n) % MMSymbols(1) ) == adjustl( atom(atj) % MMSymbol ) )
 
                   if ( flag1 .OR. flag2 ) then       ! <== apply SpecialPair parms ... 
                       sr2 = ( SpecialPairs14(n)%Parms(1) * SpecialPairs14(n)%Parms(1) ) / rklq
                       eps = SpecialPairs14(n) % Parms(2) 
                       exit read_loop1
                   end if
                
               end do read_loop1

            end if

            rklsq = sqrt(rklq)
            sr6   = sr2 * sr2 * sr2
            sr12  = sr6 * sr6
            fs    = 24.d0 * eps * ( TWO * sr12 - sr6 )
            fs    = (fs / rklq) * MM % fudgeLJ
            atom(ati) % fnonbd14(1:3) = atom(ati) % fnonbd14(1:3) + fs * rij(1:3)
            atom(atj) % fnonbd14(1:3) = atom(atj) % fnonbd14(1:3) - fs * rij(1:3)
            ! factor used to compensate the factor1 and factor2 factors ...
            ! factor3 = 1.0d-20
            sterm  = 4.d0 * eps * factor3 * ( sr12 - sr6 ) 
!           alternative cutoff formula ...
!           sterm  = sterm - vscut(ati,atj) + fscut(ati,atj) * ( rklsq - rcut ) 
            sterm  = sterm * MM % fudgeLJ

            !  Real part (Numero de cargas igual ao numero de sitios)
            sr2   = 1.d0 / rklq
            KRIJ  = KAPPA * rklsq
            expar = EXP( -(KRIJ*KRIJ) )
            freal = coulomb * chrgi * chrgj * ( sr2/rklsq )
            freal = freal * ( ERFC(KRIJ) + TWO * rsqpi * KAPPA * rklsq * expar ) * MM % fudgeQQ
            atom(ati) % fnonch14(1:3) = atom(ati) % fnonch14(1:3) + freal * rij(1:3)
            atom(atj) % fnonch14(1:3) = atom(atj) % fnonch14(1:3) - freal * rij(1:3)
            ! factor used to compensate the factor1 and factor2 factors ...
            ! factor3 = 1.0d-20
            tterm = coulomb * factor3 * chrgi * chrgj * ERFC(KRIJ)/rklsq * MM % fudgeQQ
!           alternative cutoff formula ...
!           tterm = tterm - vrecut * chrgi * chrgj + frecut * chrgi * chrgj * ( rklsq-rcut ) * factor3
            LJ_14   = LJ_14   + sterm
            Coul_14 = Coul_14 + tterm

        end if
    end do
end do

!====================================================================
! Lennard-Jones intramolecular interactions ...

If( allocated(SpecialPairs) ) there_are_NB_SpecialPairs = .true.

do i = 1 , MM % N_of_molecules
    do j   = 1 , molecule(i) % NintraLJ
        ati  = molecule(i) % IntraLJ(j,1) 
        ati1 = atom(ati) % my_intra_id + species_offset( atom(ati)%my_species )
        atj  = molecule(i) % IntraLJ(j,2) 
        atj1 = atom(atj) % my_intra_id + species_offset( atom(atj)%my_species ) 

        if ( atom(atj) % flex .OR. atom(ati) % flex ) then
            chrgi  = atom(ati) % charge
            chrgj  = atom(atj) % charge
            rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rklq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            if ( rklq < rcutsq ) then

            ! Lennard Jones ...
            select case ( MM % CombinationRule )

                case (2)
                    ! AMBER FF :: GMX COMB-RULE 2
                    sr2 = ( ( atom(ati) % sig + atom(atj) % sig ) * ( atom(ati) % sig + atom(atj) % sig ) ) / rklq

                case (3)
                    ! OPLS  FF :: GMX COMB-RULE 3
                    sr2 = ( ( atom(ati) % sig * atom(atj) % sig ) * ( atom(ati) % sig * atom(atj) % sig ) ) / rklq

            end select
            eps   =  atom(ati) % eps * atom(atj) % eps

            If( there_are_NB_SpecialPairs ) then    ! <== check whether (I,J) is a SpecialPair ... 

               read_loop: do  n = 1, size(SpecialPairs)

                  flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
                          ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(atj) % MMSymbol ) )
                  flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
                          ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(atj) % MMSymbol ) )

                  if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
                     sr2 = ( SpecialPairs(n)%Parms(1) * SpecialPairs(n)%Parms(1) ) / rklq 
                     eps = SpecialPairs(n) % Parms(2) 
                     exit read_loop
                  end if

               end do read_loop

            end if

            rklsq = SQRT(rklq)
            sr6   = sr2 * sr2 * sr2
            sr12  = sr6 * sr6
            fs    = 24.d0 * eps * ( TWO * sr12 - sr6 )
            ! with force cut-off ...
            fs    = (fs / rklq) - fscut(ati1,atj1) / rklsq     
            atom(ati) % fnonbd(1:3) = atom(ati) % fnonbd(1:3) + fs * rij(1:3)
            atom(atj) % fnonbd(1:3) = atom(atj) % fnonbd(1:3) - fs * rij(1:3)
            ! factor used to compensate factor1 ...
            ! factor3 = 1.0d-20
            sterm  = 4.d0 * eps * factor3 * ( sr12 - sr6 )
            ! alternative formula with cutoff ...
            sterm  = sterm - vscut(ati1,atj1) + fscut(ati1,atj1) * ( rklsq - rcut ) 

            !  Real part (Number of charges equal to the number of sites)
            sr2   = 1.d0 / rklq
            KRIJ  = KAPPA * rklsq
            expar = EXP( -(KRIJ*KRIJ) )
            freal = coulomb * chrgi * chrgj * ( sr2/rklsq )
            freal = freal * ( ERFC(KRIJ) + TWO * rsqpi * KAPPA * rklsq * expar )
            ! with force cut-off ...
            freal = freal - frecut / rklsq * chrgi * chrgj   
            atom(ati) % fnonch(1:3) = atom(ati) % fnonch(1:3) + freal * rij(1:3)
            atom(atj) % fnonch(1:3) = atom(atj) % fnonch(1:3) - freal * rij(1:3)
            ! factor used to compensate factor1 ...
            ! factor3 = 1.0d-20
            tterm = coulomb * factor3 * chrgi * chrgj * ERFC(KRIJ)/rklsq
            ! alternative formula with cutoff ...
            tterm = tterm - vrecut * chrgi * chrgj + frecut * chrgi * chrgj * ( rklsq-rcut ) * factor3
            LJ_intra   = LJ_intra   + sterm
            Coul_intra = Coul_intra + tterm

            end if
        end if
    end do
end do
!====================================================================
! Morse Intra/Inter potential for H transfer ...

If( allocated(SpecialMorse) ) then

   do k = 1 , MM % N_of_atoms - 1
       do l = k , MM % N_of_atoms
       read_loop2: do  n = 1, size(SpecialMorse) 
           flag1 = ( adjustl( SpecialMorse(n) % MMSymbols(1) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                   ( adjustl( SpecialMorse(n) % MMSymbols(2) ) == adjustl( atom(l) % MMSymbol ) )
           flag2 = ( adjustl( SpecialMorse(n) % MMSymbols(2) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                   ( adjustl( SpecialMorse(n) % MMSymbols(1) ) == adjustl( atom(l) % MMSymbol ) ) 
           if ( flag1 .OR. flag2 ) then
               atk = atom(k) % my_id
               atl = atom(l) % my_id
               rkl(:)  = atom(atk) % xyz(:) - atom(atl) % xyz(:)
               rkl(:)  = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)
               rklq    = rkl(1)*rkl(1) + rkl(2)*rkl(2) + rkl(3)*rkl(3)
               rklsq   = SQRT(rklq)
   
               ! Morse potential ...
               qterm0 = exp( -SpecialMorse(n) % Parms(3) * ( rklsq - SpecialMorse(n) % Parms(2) ) )
               qterm  = SpecialMorse(n) % Parms(1) * ( D_ONE - qterm0 )*( D_ONE - qterm0 )
               coephi = TWO * SpecialMorse(n) % Parms(1) * SpecialMorse(n) % Parms(3) * qterm0 * ( 1.d0 - qterm0 ) / rklsq
               atom(atk) % fMorse(:) = atom(atk) % fMorse(:) - coephi*rkl(:)
               atom(atl) % fMorse(:) = atom(atl) % fMorse(:) + coephi*rkl(:)
               Morspot = qterm + Morspot
           end if
       end do read_loop2
       end do
   end do

end If

!
!====================================================================
! factor used to compensate the factor1 and factor2 factors ...
! factor3 = 1.0d-20
pot_INTRA = ( bdpot + angpot + dihpot )*factor3 + LJ_14 + LJ_intra + Coul_14 + Coul_intra
pot_total = pot_INTER + pot_INTRA
pot_total = pot_total * mol * micro / MM % N_of_molecules

! Get total force; force units = J/mts = Newtons ...
do i = 1 , MM % N_of_atoms
    
    atom(i)% f_MM(:) = atom(i)% f_MM(:) + (atom(i) % fbond(:)    +  &
                                           atom(i) % fang(:)     +  &
                                           atom(i) % fdihed(:)   +  &
                                           atom(i) % fnonbd14(:) +  & 
                                           atom(i) % fnonch14(:) +  &
                                           atom(i) % fnonbd(:)   +  & 
                                           atom(i) % fMorse(:)   +  & 
                                           atom(i) % fnonch(:)      &
                                          ) * Angs_2_mts

    If ( QMMM ) then
       atom(i)% ftotal(:) = atom(i)% f_MM(:) + atom(i)% Ehrenfest(:)
       else
       atom(i)% ftotal(:) = atom(i)% f_MM(:) 
       end If

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
real*8  :: psi , cos_Psi , dtheta

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
!==================
 subroutine not_gmx
!==================
implicit none

! local variables ...

select case( adjustl(molecule(i) % Dihedral_Type(j)) )
    case ('cos') ! V = k[1 + cos(n.phi - theta)]
        term  = int(molecule(i) % kdihed0(j,3)) * phi - molecule(i) % kdihed0(j,1)
        pterm = molecule(i) % kdihed0(j,2) * ( 1.d0 + cos(term) )
        proper_dih = proper_dih + pterm
        gamma = - molecule(i) % kdihed0(j,2) * int(molecule(i) % kdihed0(j,3)) * sin(term) * rsinphi * rijkj * rjkkl
         
    case('harm') ! V = 1/2.k(phi - phi0)²
        dphi  = phi - molecule(i) % kdihed0(j,1)
        dphi  = dphi - DNINT(dphi * rtwopi) * twopi
        dphi1 = phi - molecule(i) % kdihed0(j,3)
        dphi1 = dphi1 - DNINT(dphi1 * rtwopi) * twopi
        dphi2 = phi - molecule(i) % kdihed0(j,5)
        dphi2 = dphi2 - DNINT(dphi2 * rtwopi) * twopi
        term  = molecule(i) % kdihed0(j,2) * dphi
        term  = term + molecule(i) % kdihed0(j,4) * dphi1
        term  = term + molecule(i) % kdihed0(j,6) * dphi2
        pterm = HALF * term * dphi
        pterm = pterm + HALF * term1 * dphi1
        pterm = pterm + HALF * term2 * dphi2
        harm_dih = harm_dih + pterm
        gamma = term * rsinphi * rijkj * rjkkl          

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



