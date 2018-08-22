
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

    use omp_lib

    public :: FORCEINTRA, pot_INTRA

    ! module variables ...
    real*8 :: pot_INTRA

contains

! Forceintra
subroutine FORCEINTRA()
    implicit none

    ! local
    integer :: ati, atj, atk, atl, ati1, atj1

    logical :: flag

    ! loop indexes
    integer :: i, j

    real*8 :: time_start
    real*8 :: time_end

    time_start = omp_get_wtime()

    do i = 1, MM%N_of_atoms
        atom(i)%fnonbd14(:) = 0         ! Non-bonded 1-4
        atom(i)%fnonbd(:)   = 0         ! Non-bonded Intramolecular
        atom(i)%fbond(:)    = 0         ! Stretching/Bonding
        atom(i)%fang(:)     = 0         ! Bending/Angular
        atom(i)%fdihed(:)   = 0         ! Dihedral
        atom(i)%fnonch14(:) = 0         ! Non-bonded coulomb 1-4
        atom(i)%fnonch(:)   = 0         ! Non-bonded coulomb Intramolecular
        atom(i)%fMorse(:)   = 0         ! Non-bonded Morse
    end do

    bdpot      = 0
    angpot     = 0
    dihpot     = 0
    Morspot    = 0
    harm_bond  = 0
    morse_bond = 0
    proper_dih = 0
    harm_dih   = 0
    ryck_dih   = 0
    imp_dih    = 0
    LJ_14      = 0
    LJ_intra   = 0
    Coul_14    = 0
    Coul_intra = 0


    ! Bonding - stretching potential ...
    do i = 1, MM%N_of_molecules
        do j = 1, molecule(i)%Nbonds

            ati = molecule(i)%bonds(j, 1)
            atj = molecule(i)%bonds(j, 2)

            if(.NOT. (atom(ati)%flex .OR. atom(atj)%flex)) then
                cycle
            end if

            call STRETCHING_POTENTIAL(atom(ati), atom(atj), molecule(i)%kbond0(j,:), molecule(i)%bond_type(j))
        end do
    end do


    ! Angle - bending potential ...
    do i = 1, MM%N_of_molecules
        do j = 1, molecule(i)%Nangs
            atj = molecule(i)%angs(j, 1)
            ati = molecule(i)%angs(j, 2)
            atk = molecule(i)%angs(j, 3)

            if (.NOT. (atom(atj)%flex .OR. atom(ati)%flex .OR. atom(atk)%flex)) then
                cycle
            end if

            ! Harmonic and Urey-Bradley potentials ...
            if (.NOT. (molecule(i)%angle_type(j) == 'harm' .OR. molecule(i)%angle_type(j) == 'urba')) then
                cycle
            end if

            call BENDING_POTENTIAL(ati, atj, atk, molecule(i)%kang0(j,:))

            ! Urey-Bradley bonding term ...
            if( molecule(i)%angle_type(j) == "urba" ) then
                call UREY_BRADLEY_BONDING(atom(atk), atom(atj), molecule(i)%kang0(j,:))
            end if
        end do
    end do


    ! Dihedral Potential Angle ...
    do i = 1, MM%N_of_molecules
        do j = 1, molecule(i)%Ndiheds
            ati = molecule(i)%diheds(j,1)
            atj = molecule(i)%diheds(j,2)
            atk = molecule(i)%diheds(j,3)
            atl = molecule(i)%diheds(j,4)

            if (atom(atj)%flex .OR. atom(ati)%flex .OR. atom(atk)%flex .OR. atom(atl)%flex) then
                call DIHEDDRAL_POTENTIAL_ANGLE(atom(ati), atom(atj), atom(atk), atom(atl), molecule(i)%kdihed0(j,:), molecule(i)%Dihedral_Type(j))
            end if

        end do
    end do


    ! Non-bonded 1,4 intramolecular interactions ...
    do i = 1, MM%N_of_molecules
        do j = 1, molecule(i)%Nbonds14
            ati = molecule(i)%bonds14(j,1)
            atj = molecule(i)%bonds14(j,2)

            if (atom(atj)%flex .OR. atom(ati)%flex) then
                call NON_BONDED_INTRAMOLECULAR_INTERACTIONS(atom(atj), atom(ati))
            end if

        end do
    end do

    ! Lennard Jones approximation
    call LENNARD_JONES()


    ! Morse Intra/Inter potential for H transfer ...
    call MORSE_POTENTIAL()

    pot_INTRA = (bdpot + angpot + dihpot) * factor3 + LJ_14 + LJ_intra + Coul_14 + Coul_intra
    pot_total = pot_INTER + pot_INTRA
    pot_total = pot_total * mol * micro / MM%N_of_molecules

    ! Get total force; force units = J/mts = Newtons ...
    do i = 1, MM%N_of_atoms
        atom(i)%ftotal(:) = atom(i)%ftotal(:) + (atom(i)%fbond(:)    +  &
                                                 atom(i)%fang(:)     +  &
                                                 atom(i)%fdihed(:)   +  &
                                                 atom(i)%fnonbd14(:) +  &
                                                 atom(i)%fnonch14(:) +  &
                                                 atom(i)%fnonbd(:)   +  &
                                                 atom(i)%fMorse(:)   +  &
                                                 atom(i)%fnonch(:)      &
                                                 ) * Angs_2_mts

        ! Append total forces with Excited State Ehrenfest terms; nonzero only if QMMM = T_ ...
        if(QMMM) then
            atom(i)%ftotal(:) = atom(i)%ftotal(:) + atom(i)%Ehrenfest(:)
        end if
    end do

    time_end = omp_get_wtime()

    print*, 'time==>', time_end - time_start
end subroutine FORCEINTRA


! Error function.
!
! Used to approximate values in statistics.
! For more info, see:
!   https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
!
! @param X[in]
!
! @return ERFC
pure function ERFC(X)
    ! args
    real*8, intent(in) :: X

    ! result
    real*8 :: ERFC

    ! local
    real*8 :: A1, A2, A3, A4, A5, P, T, TP

    ! Values taken from the wikipedia article
    parameter (A1 =  0.254829592)
    parameter (A2 = -0.284496736)
    parameter (A3 =  1.421413741)
    parameter (A4 = -1.453122027)
    parameter (A5 =  1.061405429)
    parameter (P  =  0.3275911)

    T = 1 / (1 + P * X)
    TP = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
    ERFC = TP * EXP(-(X ** 2))
end function ERFC


! Checks if number are equal.
!
! @param[in] X
! @param[in] Y
!
! @return 1.0 if equal, 0.0 otherwise;
pure function DEL(X, Y)
    implicit none

    ! args
    integer, intent(in) :: X
    integer, intent(in) :: Y

    ! result
    real*8  :: DEL

    if (X == Y) then
        DEL = 1
    else
        DEL = 0
    end if
end function DEL


! Get the flags status of the atoms
!
! @param[in] pair_number
! @param[in] atom_i
! @param[in] atom_j
pure function GET_FLAGS(pair_number, atom_i, atom_j)
    implicit none

    ! args
    integer, intent(in) :: pair_number
    type(MM_atomic), intent(in) :: atom_i
    type(MM_atomic), intent(in) :: atom_j

    ! result
    logical :: GET_FLAGS

    ! local
    logical :: flag1
    logical :: flag2
    character(4) :: mm_symbols1
    character(4) :: mm_symbols2
    character(4) :: symbol_i
    character(4) :: symbol_j

    symbol_i = atom_i%MMSymbol
    symbol_j = atom_j%MMSymbol
    mm_symbols1 = SpecialPairs(pair_number)%MMSymbols(1)
    mm_symbols2 = SpecialPairs(pair_number)%MMSymbols(2)

    flag1 =             adjustl(mm_symbols1) == adjustl(symbol_i)
    flag1 = flag1 .AND. adjustl(mm_symbols2) == adjustl(symbol_j)
    flag2 =             adjustl(mm_symbols2) == adjustl(symbol_i)
    flag2 = flag2 .AND. adjustl(mm_symbols1) == adjustl(symbol_j)

    GET_FLAGS = flag1 .OR. flag2
end function GET_FLAGS


! Get morse flag status for the atoms passed.
!
! @param pair_number[in]
! @param atom_i[in]
! @param atom_j[in]
!
! @return boolean
pure function GET_FLAGS_MORSE(pair_number, atom_i, atom_j)
    implicit none

    ! args
    integer, intent(in) :: pair_number
    type(MM_atomic), intent(in) :: atom_i
    type(MM_atomic), intent(in) :: atom_j

    ! result
    logical :: GET_FLAGS_MORSE

    ! local
    logical :: flag1
    logical :: flag2
    character(4) :: mm_symbols1
    character(4) :: mm_symbols2
    character(4) :: symbol_i
    character(4) :: symbol_j

    symbol_i = atom_i%MMSymbol
    symbol_j = atom_j%MMSymbol
    mm_symbols1 = SpecialMorse(pair_number)%MMSymbols(1)
    mm_symbols2 = SpecialMorse(pair_number)%MMSymbols(2)

    flag1 =             adjustl(mm_symbols1) == adjustl(symbol_i)
    flag1 = flag1 .AND. adjustl(mm_symbols2) == adjustl(symbol_j)
    flag2 =             adjustl(mm_symbols2) == adjustl(symbol_i)
    flag2 = flag2 .AND. adjustl(mm_symbols1) == adjustl(symbol_j)

    GET_FLAGS_MORSE = flag1 .OR. flag2

    ! flag1 = ( adjustl( SpecialMorse(n)%MMSymbols(1) ) == adjustl( atom(k)%MMSymbol ) ) .AND. &
    !         ( adjustl( SpecialMorse(n)%MMSymbols(2) ) == adjustl( atom(l)%MMSymbol ) )
    ! flag2 = ( adjustl( SpecialMorse(n)%MMSymbols(2) ) == adjustl( atom(k)%MMSymbol ) ) .AND. &
    !         ( adjustl( SpecialMorse(n)%MMSymbols(1) ) == adjustl( atom(l)%MMSymbol ) )
end function GET_FLAGS_MORSE


! Get the flags realted to atoms i and j
!
! @param pair_number[in]
! @param atom_i[in]
! @param atom_j[in]
pure function GET_FLAGS_1_4(pair_number, atom_i, atom_j)
    implicit none

    ! args
    integer, intent(in) :: pair_number
    type(MM_atomic), intent(in) :: atom_i
    type(MM_atomic), intent(in) :: atom_j

    ! result
    logical :: GET_FLAGS_1_4

    ! local
    character(4) :: mm_symbols1
    character(4) :: mm_symbols2
    character(4) :: symbol_i
    character(4) :: symbol_j

    logical :: flag1
    logical :: flag2

    symbol_i = atom_i%MMSymbol
    symbol_j = atom_j%MMSymbol
    mm_symbols1 = SpecialPairs14(pair_number)%MMSymbols(1)
    mm_symbols2 = SpecialPairs14(pair_number)%MMSymbols(2)


    flag1 =             adjustl(mm_symbols1) == adjustl(symbol_i)
    flag1 = flag1 .AND. adjustl(mm_symbols2) == adjustl(symbol_j)
    flag2 =             adjustl(mm_symbols2) == adjustl(symbol_i)
    flag2 = flag2 .AND. adjustl(mm_symbols1) == adjustl(symbol_j)

    GET_FLAGS_1_4 = flag1 .OR. flag2
end function GET_FLAGS_1_4


subroutine MORSE_POTENTIAL()
    implicit none

    ! local
    real*8 :: local_fmorse(SIZE(atom), 3)

    integer :: i
    integer :: j
    integer :: k

    logical :: flag

    local_fmorse = 0

    !$OMP parallel do schedule(guided) &
    !$OMP   default(shared) &
    !$OMP   private(i, j, k, flag) &
    !$OMP   reduction(+ : Morspot, local_fmorse)
    do i = 1, MM%N_of_atoms - 1
        do j = i, MM%N_of_atoms
            do k = 1, SIZE(SpecialMorse)
                flag = GET_FLAGS_MORSE(k, atom(i), atom(j))

                if (flag) then
                    call MORSE_POTENTIAL_AUX(i, j, SpecialMorse(k)%parms(:), local_fmorse(:,:))
                end if
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do schedule(static) private(i) default(shared)
    do i = 1, MM%N_of_atoms
        atom(i)%fMorse(:) = atom(i)%fMorse(:) + local_fmorse(i, :)
    end do
    !$OMP end parallel do
end subroutine MORSE_POTENTIAL

! Calculate Morse potential between 2 atoms.
!
! @param i[in]
! @param j[in]
! @param parms[in]
!
subroutine MORSE_POTENTIAL_AUX(i, j, parms, local_fmorse)
    implicit none

    ! args
    integer, intent(in) :: i
    integer, intent(in) :: j
    real*8, intent(in) :: parms(:)
    real*8, intent(inout) :: local_fmorse(:,:)

    ! local
    real*8 :: rkl(3)
    real*8 :: rkl_sqrt

    real*8 :: qterm0
    real*8 :: qterm
    real*8 :: coephi

    rkl(:) = DIFFERENCE_BOX(atom(i), atom(j))
    rkl_sqrt = SQRT(SUM(rkl(:) ** 2))

    ! Morse potential ...
    qterm0 = EXP(-parms(3) * (rkl_sqrt - parms(2)))
    qterm  = parms(1) * (1 - qterm0) ** 2
    coephi = 2 * parms(1) * parms(3) * qterm0 * (1 - qterm0) / rkl_sqrt

    local_fmorse(i, :) = local_fmorse(i, :) - (coephi * rkl(:))
    local_fmorse(j, :) = local_fmorse(j, :) + (coephi * rkl(:))

    Morspot = Morspot + qterm
end subroutine MORSE_POTENTIAL_AUX


! Calculate the box edges for the atoms passed.
!
! @param atom_i[in]
! @param atom_j[in]
!
! @return diff
pure function DIFFERENCE_BOX(atom_i, atom_j) result(diff)
    implicit none
    type(MM_atomic), intent(in) :: atom_i
    type(MM_atomic), intent(in) :: atom_j

    real*8 :: diff(3)

    diff(:) = atom_i%xyz(:) - atom_j%xyz(:)
    diff(:) = diff(:) - MM%box(:) * DNINT(diff(:) * MM%ibox(:)) * PBC(:)
end function DIFFERENCE_BOX


! Calculates the cross product of two vectors
!
! @param vec_a[in]
! @param vec_b[in]
!
! @return Cross product of vec_a and vec_b
pure function CROSS_PRODUCT(vec_a, vec_b) result(cross)
    implicit none
    real*8, intent(in) :: vec_a(3)
    real*8, intent(in) :: vec_b(3)

    real*8 :: cross(3)

    cross(1) = vec_a(2) * vec_b(3) - vec_a(3) * vec_b(2)
    cross(2) = vec_a(3) * vec_b(1) - vec_a(1) * vec_b(3)
    cross(3) = vec_a(1) * vec_b(2) - vec_a(2) * vec_b(1)
end function CROSS_PRODUCT


! Lennard Jones approximation
subroutine LENNARD_JONES()
    implicit none

    ! local
    integer :: ati
    integer :: ati1
    integer :: atj
    integer :: atj1

    real*8 :: local_fnonbd(SIZE(atom), 3)
    real*8 :: local_fnonch(SIZE(atom), 3)

    integer , allocatable :: species_offset(:)

    ! loop index
    integer :: i
    integer :: j
    integer :: n

    CALL offset(species_offset)

    local_fnonbd = 0
    local_fnonch = 0
    !$OMP parallel &
    !$OMP   default(shared) &
    !$OMP   private(i) &
    !$OMP   reduction(+ : LJ_intra, Coul_intra, local_fnonbd, local_fnonch)
    do i = 1, MM%N_of_molecules
        !$OMP do &
        !$OMP   schedule(static) &
        !$OMP   private(j, ati, ati1, atj, atj1)
        do j = 1, molecule(i)%NintraLJ
            ati  = molecule(i)%IntraLJ(j,1)
            ati1 = atom(ati)%my_intra_id + species_offset(atom(ati)%my_species)
            atj  = molecule(i)%IntraLJ(j,2)
            atj1 = atom(atj)%my_intra_id + species_offset(atom(atj)%my_species)

            if (atom(atj)%flex .OR. atom(ati)%flex) then
                call LENNARD_JONES_AUX(ati, ati1, atj, atj1, local_fnonbd, local_fnonch, MM%CombinationRule)
            end if

        end do
        !$OMP end do
    end do
    !$OMP end parallel

    !$OMP parallel do schedule(static) private(n) default(shared)
    do n = 1, SIZE(atom)
        atom(n)%fnonbd(:) = atom(n)%fnonbd(:) + local_fnonbd(n,:)
        atom(n)%fnonch(:) = atom(n)%fnonch(:) + local_fnonch(n,:)
    end do
    !$OMP end parallel do
end subroutine LENNARD_JONES


! Auxiliar function for the Lennard Jones subrouinte
!
! This subroutine is made this way to avoid some problema with multi-threading.
!
! @param ati[in]
! @param ati1[in]
! @param atj[in]
! @param atj1[in]
! @param combination_rule[in]
! @param local_fnonbd[inout]
! @param local_fnonch[inout]
subroutine LENNARD_JONES_AUX(ati, ati1, atj, atj1, local_fnonbd, local_fnonch, combination_rule)
    implicit none

    ! args
    integer, intent(in) :: ati
    integer, intent(in) :: ati1
    integer, intent(in) :: atj
    integer, intent(in) :: atj1
    integer, intent(in) :: combination_rule
    real*8, intent(inout) :: local_fnonbd(:,:)
    real*8, intent(inout) :: local_fnonch(:,:)

    ! local
    logical :: flag
    real*8 :: sr2
    real*8 :: eps
    real*8 :: chrgi
    real*8 :: chrgj
    real*8 :: rklq
    real*8 :: rklq_sqrt
    real*8 :: fs
    real*8 :: sterm
    real*8 :: KRIJ
    real*8 :: expar
    real*8 :: freal
    real*8 :: tterm

    real*8 :: rij(3)

    ! loop index
    integer :: n

    rij(:) = DIFFERENCE_BOX(atom(ati), atom(atj))
    rklq = SUM(rij(:) ** 2)

    if (.NOT. (rklq < rcutsq)) then
        return
    end if

    ! Lennard Jones ...
    ! AMBER FF :: GMX COMB-RULE 2
    if (combination_rule == 2) then
        sr2 = (atom(ati)%sig + atom(atj)%sig) ** 2 / rklq
    ! OPLS  FF :: GMX COMB-RULE 3
    else if (combination_rule == 3) then
        sr2 = (atom(ati)%sig ** 2) * (atom(atj)%sig ** 2) / rklq
    end if

    eps = atom(ati)%eps * atom(atj)%eps

    ! Nbond_parms directive on ...
    do n = 1, SIZE(SpecialPairs)
        flag = GET_FLAGS(n, atom(ati), atom(atj))

        if (.NOT. (flag)) then
            cycle
        end if

        ! AMBER FF :: GMX COMB-RULE 2
        if (combination_rule == 2) then
            sr2 = (SpecialPairs(n)%Parms(1) * 2) ** 2 / rklq
        ! OPLS  FF :: GMX COMB-RULE 3
        else if (combination_rule == 3) then
            sr2 = SpecialPairs(n)%Parms(1) ** 4 / rklq
        end if

        eps = SpecialPairs(n)%Parms(2) ** 2
        exit
    end do

    rklq_sqrt = SQRT(rklq)
    fs = 24 * eps * (2* sr2 ** 6 - sr2 ** 3)
    ! with force cut-off ...
    fs = (fs / rklq) - fscut(ati1, atj1) / rklq_sqrt

    ! factor3 used to compensate factor1 ...
    sterm = 4 * eps * factor3 * (sr2 ** 6 - sr2 ** 3)
    ! alternative formula with cutoff ...
    sterm = sterm - vscut(ati1, atj1) + fscut(ati1, atj1) * (rklq_sqrt - rcut)

    chrgi = atom(ati)%charge
    chrgj = atom(atj)%charge

    !  Real part (Number of charges equal to the number of sites)
    sr2 = 1.0 / rklq
    KRIJ  = KAPPA * rklq_sqrt
    expar = EXP(-(KRIJ ** 2.0))
    freal = coulomb * chrgi * chrgj * (sr2 / rklq_sqrt)
    freal = freal * (ERFC(KRIJ) + 2.0 * rsqpi * KAPPA * rklq_sqrt * expar)
    ! with force cut-off ...
    freal = freal - frecut / rklq_sqrt * chrgi * chrgj

    ! factor used to compensate factor1 ...
    ! factor3 = 1.0d-20
    tterm = (coulomb * factor3 * chrgi * chrgj * ERFC(KRIJ)) / rklq_sqrt
    ! alternative formula with cutoff ...
    tterm = tterm - vrecut * chrgi * chrgj
    tterm = tterm + frecut * chrgi * chrgj * (rklq_sqrt - rcut) * factor3

    local_fnonbd(ati,:) = local_fnonbd(ati,:) + fs * rij(:)
    local_fnonbd(atj,:) = local_fnonbd(atj,:) - fs * rij(:)
    local_fnonch(ati,:) = local_fnonch(ati,:) + freal * rij(:)
    local_fnonch(atj,:) = local_fnonch(atj,:) - freal * rij(:)

    LJ_intra   = LJ_intra   + sterm
    Coul_intra = Coul_intra + tterm
end subroutine LENNARD_JONES_AUX


! Calculates the sctetching potential of 2 atoms
!
! @param atom_i[inout]
! @param atom_j[inout]
! @param kbond0[in]
subroutine STRETCHING_POTENTIAL(atom_i, atom_j, kbond0, bond_type)
    implicit none

    ! args
    real*8, intent(in) :: kbond0(:)
    character(4), intent(in) :: bond_type
    type(MM_atomic), intent(inout) :: atom_i
    type(MM_atomic), intent(inout) :: atom_j

    ! local
    real*8 :: rij(3)
    real*8 :: rij_sqrt
    real*8 :: qterm
    real*8 :: qterm0
    real*8 :: coephi

    rij(:) = DIFFERENCE_BOX(atom_j, atom_i)
    rij_sqrt = SQRT(SUM(rij(:) ** 2))

    select case (bond_type)
        ! harmonic potential ...
        case ( "harm" )
        qterm   = 0.5 * kbond0(1) * ((rij_sqrt - kbond0(2)) ** 2)
        coephi  = (kbond0(1) * (rij_sqrt - kbond0(2))) / rij_sqrt
        harm_bond = qterm + harm_bond

        ! Morse potential ...
        case ( "Mors" )
        qterm0 = EXP(-kbond0(3) * (rij_sqrt - kbond0(2)))
        qterm  = kbond0(1) * ((1 - qterm0) ** 2)
        coephi = (2 * kbond0(1) * kbond0(3) * qterm0 * (1 - qterm0)) / rij_sqrt
        morse_bond = qterm + morse_bond
    end select

    atom_j%fbond(:) = atom_j%fbond(:) - (coephi * rij(:))
    atom_i%fbond(:) = atom_i%fbond(:) + (coephi * rij(:))

    bdpot = qterm + bdpot
end subroutine STRETCHING_POTENTIAL


! Calculates bending pontential between 3 atoms
!
! @param ati[in]
! @param atj[in]
! @param atk[in]
subroutine BENDING_POTENTIAL(ati, atj, atk, kang0)
    implicit none

    ! args
    integer, intent(in) :: ati
    integer, intent(in) :: atj
    integer, intent(in) :: atk
    real*8, intent(in) :: kang0(:)

    ! local
    real*8 :: rik(3)
    real*8 :: rij(3)
    real*8 :: rjk(3)
    real*8 :: rij_sqrt
    real*8 :: rjk_sqrt
    real*8 :: rik_sqrt

    integer :: atl

    ! loop indexes
    integer :: l
    integer :: loop

    real*8 :: phi
    real*8 :: fxyz
    real*8 :: riku
    real*8 :: riju
    real*8 :: coephi
    real*8 :: rterm

    rij(:) = DIFFERENCE_BOX(atom(atj), atom(ati))
    rjk(:) = DIFFERENCE_BOX(atom(atk), atom(ati))
    rij_sqrt  = SQRT(SUM(rij(:) ** 2))
    rjk_sqrt  = SQRT(SUM(rjk(:) ** 2))

    phi = ACOS(SUM(rij(:) * rjk(:)) / (rij_sqrt * rjk_sqrt))

    coephi = 0
    rterm  = 0
    if (.NOT. (phi < 1.d-12 .OR. ABS(pi - phi) < 1.d-12)) then
        coephi = (phi - kang0(2)) / SIN(phi)
        rterm  = 0.5 * kang0(1) * ((phi - kang0(2)) ** 2)
    end if

    angpot = rterm + angpot

    ! I have no idea what is happening in these loops
    do l = 1, 3
        if (l == 1) atl = atj
        if (l == 2) atl = ati
        if (l == 3) atl = atk
        do loop = 1, 3                    ! X,Y,Z axis (n = 1, 2 or 3)
            riju = rij(loop)
            riku = rjk(loop)
            fxyz = (kang0(1) * coephi) *                             &
                ((DEL(atl, atj) - DEL(atl, ati)) * riku / (rij_sqrt * rjk_sqrt) +   &
                 (DEL(atl, atk) - DEL(atl, ati)) * riju / (rij_sqrt * rjk_sqrt) -   &
                ((DEL(atl, atj) - DEL(atl, ati)) * riju / (rij_sqrt * rij_sqrt) +   &
                 (DEL(atl, atk) - DEL(atl, ati)) * riku / (rjk_sqrt * rjk_sqrt)) * COS(phi))

            atom(atl)%fang(loop) = atom(atl)%fang(loop) + fxyz
        end do
    end do
end subroutine BENDING_POTENTIAL


! Calculates the Urey-Bradley bonding values.
!
! @param atom_i[inout]
! @param atom_j[inout]
! @param kang0[in]
subroutine UREY_BRADLEY_BONDING(atom_i, atom_j, kang0)
    implicit none

    ! args
    type(MM_atomic), intent(inout) :: atom_i
    type(MM_atomic), intent(inout) :: atom_j
    real*8, intent(in) :: kang0(:)

    ! local
    real*8 :: rij(3)
    real*8 :: rij_sqrt
    real*8 :: coephi
    real*8 :: rterm

    rij(:) = DIFFERENCE_BOX(atom_i, atom_j)
    rij_sqrt  = SQRT(SUM(rij(:) ** 2))

    coephi = kang0(3) * (rij_sqrt - kang0(4)) / rij_sqrt
    rterm  = 0.5 * kang0(3) * ((rij_sqrt - kang0(4)) ** 2)

    atom_i%fang(:) = atom_i%fang(:) - (coephi * rij(:))
    atom_j%fang(:) = atom_j%fang(:) + (coephi * rij(:))

    angpot = rterm + angpot
end subroutine UREY_BRADLEY_BONDING


! Calculates the dihedral potential angle
!
! @param atom_i[inout]
! @param atom_j[inout]
! @param atom_k[inout]
! @param atom_l[inout]
! @param kdihed0[in]
! @param dihedral_type[in]
subroutine DIHEDDRAL_POTENTIAL_ANGLE(atom_i, atom_j, atom_k, atom_l, kdihed0, dihedral_type)
    implicit none

    ! args
    type(MM_atomic), intent(inout) :: atom_i
    type(MM_atomic), intent(inout) :: atom_j
    type(MM_atomic), intent(inout) :: atom_k
    type(MM_atomic), intent(inout) :: atom_l
    character(4), intent(in) :: dihedral_type
    real*8, intent(in) :: kdihed0(:)

    ! local
    real*8 :: gamma
    real*8 :: rjksq
    real*8 :: rjksq2
    real*8 :: rijkj
    real*8 :: rijkj2
    real*8 :: rjkkl
    real*8 :: rjkkl2
    real*8 :: coephi
    real*8 :: sinphi
    real*8 :: cosphi
    real*8 :: phi
    real*8 :: rsinphi
    real*8 :: pterm
    real*8 :: f1x, f1y, f1z, f3x, f3y, f3z, f2x, f2y, f2z, f4x, f4y, f4z

    real*8 ::rij(3)
    real*8 ::rjk(3)
    real*8 ::rkl(3)
    real*8 ::rijk(3)
    real*8 ::rjkl(3)
    real*8 ::rijkl(3)

    ! Definition of vector rij = rj - ri
    rij(:)  = DIFFERENCE_BOX(atom_i, atom_j)
    ! rij(:)  = atom_i%xyz(:) - atom_j%xyz(:)
    ! rij(:)  = rij(:) - MM%box(:) * DNINT( rij(:) * MM%ibox(:) ) * PBC(:)

    ! Definition of vector rjk = rj - rk
    rjk(:)  = DIFFERENCE_BOX(atom_j, atom_k)
    rjksq   = 1 / SQRT(SUM(rjk(:) ** 2))
    rjksq2  = rjksq ** 2

    ! Is this comment right?
    ! Definition of vector rkl = rk - rl
    rkl(:)  = DIFFERENCE_BOX(atom_k, atom_l)
    ! rkl(:)  = atom_k%xyz(:) - atom_l%xyz(:)
    ! rkl(:)  = rkl(:) - MM%box(:) * DNINT( rkl(:) * MM%ibox(:) ) * PBC(:)

    ! Cross Product M = | rij X rjk | :: First dihedral vector ...
    rijk(:) = CROSS_PRODUCT(rij, rjk)
    rijkj   = SUM(rijk(:) ** 2)
    rijkj   = 1 / SQRT(rijkj)
    rijkj2  = rijkj ** 2

    ! Cross Product N = | rjk X rkl | :: Second dihedral vector ...
    rjkl(:) = CROSS_PRODUCT(rjk, rkl)
    rjkkl   = SUM(rjkl(:) ** 2)
    rjkkl   = 1 / SQRT(rjkkl)
    rjkkl2  = rjkkl ** 2

    ! Cross Product O = | rjk X rkl | X | rij X rjk |
    rijkl(:) = CROSS_PRODUCT(rjkl, rijk)

    ! PHI is the dihedral angle defined by PHI = ACOS(B), where  B(rij,rjk,rkl) = [ (rij X rjk).(rjk X rkl) / |rij X rjk||rjk X rkl| ]
    ! and the sign of PHI is positive if the vector O is in the same direction as the bond vector rjk
    coephi = SUM(rijk(:) * rjkl(:))
    cosphi = coephi * rijkj * rjkkl

    if(ABS(cosphi) > 1) then
        cosphi = SIGN(1.0, cosphi)
    end if
    sinphi = SUM(rjk(:) * rijkl(:)) * rjksq * rijkj * rjkkl
    phi    = ATAN2(sinphi, cosphi)

    ! Avoid singularity in sinphi
    sinphi  = SIGN(MAX(1.d-10, ABS(sinphi)), sinphi)
    rsinphi = 1 / sinphi

    ! selection of potential energy function type
    if(MM_input_format == "GMX") then
        gamma = GMX(kdihed0, dihedral_type, phi, rsinphi, pterm, rijkj, rjkkl)
    else
        gamma = NOT_GMX(kdihed0, dihedral_type, phi, rsinphi, pterm, rijkj, rjkkl)
    end if

    ! Calculate atomic forces ...
    f1x = gamma * ((-rjkl(2) * rjk(3) + rjkl(3) * rjk(2)) - (-rijk(2) * rjk(3) + rijk(3) * rjk(2)) * coephi * rijkj2)
    f1y = gamma * (( rjkl(1) * rjk(3) - rjkl(3) * rjk(1)) - ( rijk(1) * rjk(3) - rijk(3) * rjk(1)) * coephi * rijkj2)
    f1z = gamma * ((-rjkl(1) * rjk(2) + rjkl(2) * rjk(1)) - (-rijk(1) * rjk(2) + rijk(2) * rjk(1)) * coephi * rijkj2)
    f3x = gamma * ((-rjkl(2) * rij(3) + rjkl(3) * rij(2)) - (-rijk(2) * rij(3) + rijk(3) * rij(2)) * coephi * rijkj2)
    f3y = gamma * (( rjkl(1) * rij(3) - rjkl(3) * rij(1)) - ( rijk(1) * rij(3) - rijk(3) * rij(1)) * coephi * rijkj2)
    f3z = gamma * ((-rjkl(1) * rij(2) + rjkl(2) * rij(1)) - (-rijk(1) * rij(2) + rijk(2) * rij(1)) * coephi * rijkj2)
    f2x = gamma * ((-rijk(2) * rkl(3) + rijk(3) * rkl(2)) - (-rjkl(2) * rkl(3) + rjkl(3) * rkl(2)) * coephi * rjkkl2)
    f2y = gamma * (( rijk(1) * rkl(3) - rijk(3) * rkl(1)) - ( rjkl(1) * rkl(3) - rjkl(3) * rkl(1)) * coephi * rjkkl2)
    f2z = gamma * ((-rijk(1) * rkl(2) + rijk(2) * rkl(1)) - (-rjkl(1) * rkl(2) + rjkl(2) * rkl(1)) * coephi * rjkkl2)
    f4x = gamma * ((-rijk(2) * rjk(3) + rijk(3) * rjk(2)) - (-rjkl(2) * rjk(3) + rjkl(3) * rjk(2)) * coephi * rjkkl2)
    f4y = gamma * (( rijk(1) * rjk(3) - rijk(3) * rjk(1)) - ( rjkl(1) * rjk(3) - rjkl(3) * rjk(1)) * coephi * rjkkl2)
    f4z = gamma * ((-rijk(1) * rjk(2) + rijk(2) * rjk(1)) - (-rjkl(1) * rjk(2) + rjkl(2) * rjk(1)) * coephi * rjkkl2)

    atom_i%fdihed(1) = atom_i%fdihed(1) + f1x
    atom_i%fdihed(2) = atom_i%fdihed(2) + f1y
    atom_i%fdihed(3) = atom_i%fdihed(3) + f1z
    atom_j%fdihed(1) = atom_j%fdihed(1) - f1x - f3x + f2x
    atom_j%fdihed(2) = atom_j%fdihed(2) - f1y - f3y + f2y
    atom_j%fdihed(3) = atom_j%fdihed(3) - f1z - f3z + f2z
    atom_k%fdihed(1) = atom_k%fdihed(1) - f2x - f4x + f3x
    atom_k%fdihed(2) = atom_k%fdihed(2) - f2y - f4y + f3y
    atom_k%fdihed(3) = atom_k%fdihed(3) - f2z - f4z + f3z
    atom_l%fdihed(1) = atom_l%fdihed(1) + f4x
    atom_l%fdihed(2) = atom_l%fdihed(2) + f4y
    atom_l%fdihed(3) = atom_l%fdihed(3) + f4z

    dihpot = dihpot + pterm
end subroutine DIHEDDRAL_POTENTIAL_ANGLE


! Auxiliay funcion to the DIHEDDRAL_POTENTIAL_ANGLE subroutine
!
! @param kdihed0[in]
! @param dihedral_type[in]
! @param phi[in]
! @param rsinphi[in]
! @param pterm[inout]
! @param rijkj[in]
! @param rjkkl[in]
!
! @return gamma
function GMX(kdihed0, dihedral_type, phi, rsinphi, pterm, rijkj, rjkkl) result(gamma)
    implicit none

    ! args
    real*8, intent(in) :: kdihed0(:)
    real*8, intent(in) :: phi
    real*8, intent(in) :: rsinphi
    real*8, intent(in) :: rijkj
    real*8, intent(in) :: rjkkl
    real*8, intent(out) :: pterm
    character(4), intent(in) :: dihedral_type

    ! result
    real*8 :: gamma

    ! local
    real*8  :: psi
    real*8  :: cos_Psi
    real*8  :: dtheta
    real*8 :: term, term1, term2, term3, term4

    select case(adjustl(dihedral_type))
        case ('cos')    ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] Eq. 4.60 (GMX 5.0.5 manual)
            term  = INT(kdihed0(3)) * phi - kdihed0(1)
            pterm = kdihed0(2) * (1 + COS(term))
            proper_dih = proper_dih + pterm
            gamma = -kdihed0(2) * INT(kdihed0(3)) * SIN(term) * rsinphi * rijkj * rjkkl

        case ('harm')   ! V = 1/2.k( xi - xi_0 )² Eq. 4.59 (GMX 5.0.5 manual)
               dtheta = phi - kdihed0(1)
               dtheta = dtheta - DNINT( dtheta * 1 / TWOPI ) * TWOPI

               term  = kdihed0(2) * dtheta
               pterm = 0.5 * term * dtheta
               harm_dih = harm_dih + pterm
               gamma = term * rsinphi * rijkj * rjkkl

        case('cos3')    ! V = C0 + C1*cos(phi - 180) + C2*cos^2(phi - 180) + C3*cos^3(phi - 180) + C4*cos^4(phi - 180) + C5*cos(phi - 180) Eq. 4.61 (GMX 5.0.5 manual)
            psi     = phi - PI
            cos_Psi = cos(psi)

            pterm =         kdihed0(1) * cos_Psi ** 0
            pterm = pterm + kdihed0(2) * cos_Psi ** 1
            pterm = pterm + kdihed0(3) * cos_Psi ** 2
            pterm = pterm + kdihed0(4) * cos_Psi ** 3
            pterm = pterm + kdihed0(5) * cos_Psi ** 4
            pterm = pterm + kdihed0(6) * cos_Psi ** 5

            ryck_dih = ryck_dih + pterm

            gamma = -SIN(psi)
            gamma = gamma + 1 * kdihed0(2) * cos_Psi ** 0
            gamma = gamma + 2 * kdihed0(3) * cos_Psi ** 1
            gamma = gamma + 3 * kdihed0(4) * cos_Psi ** 2
            gamma = gamma + 4 * kdihed0(5) * cos_Psi ** 3
            gamma = gamma + 5 * kdihed0(6) * cos_Psi ** 4
            gamma = gamma * rsinphi * rijkj * rjkkl

        case ('imp')    ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (improper) Eq. 4.60 (GMX 5.0.5 manual)
            term  = INT(kdihed0(3)) * phi - kdihed0(1)
            pterm = kdihed0(2) * (1 + COS(term))
            imp_dih = imp_dih + pterm
            gamma = - kdihed0(2) * INT(kdihed0(3)) * SIN(term) * rsinphi * rijkj * rjkkl

        case ('chrm')   ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple) Eq. 4.60 (GMX 5.0.5 manual)
            term  = INT(kdihed0(3))  * phi - kdihed0(1)
            term1 = INT(kdihed0(6))  * phi - kdihed0(4)
            term2 = INT(kdihed0(9))  * phi - kdihed0(7)
            term3 = INT(kdihed0(12)) * phi - kdihed0(10)
            term4 = INT(kdihed0(15)) * phi - kdihed0(13)

            pterm =         kdihed0(2)  * (1 + COS(term))
            pterm = pterm + kdihed0(5)  * (1 + COS(term1))
            pterm = pterm + kdihed0(8)  * (1 + COS(term2))
            pterm = pterm + kdihed0(11) * (1 + COS(term3))
            pterm = pterm + kdihed0(14) * (1 + COS(term4))

            proper_dih = proper_dih + pterm

            gamma =       - kdihed0(2)  * INT(kdihed0(3))  * SIN(term)  * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(5)  * INT(kdihed0(6))  * SIN(term1) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(8)  * INT(kdihed0(9))  * SIN(term2) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(11) * INT(kdihed0(12)) * SIN(term3) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(14) * INT(kdihed0(15)) * SIN(term4) * rsinphi * rijkj * rjkkl
    end select
end function GMX

! Auxiliay funcion to the DIHEDDRAL_POTENTIAL_ANGLE subroutine
!
! @param kdihed0[in]
! @param dihedral_type[in]
! @param phi[in]
! @param rsinphi[in]
! @param pterm[inout]
! @param rijkj[in]
! @param rjkkl[in]
!
! @return gamma
function NOT_GMX(kdihed0, dihedral_type, phi, rsinphi, pterm, rijkj, rjkkl) result(gamma)
    implicit none

    ! args
    real*8, intent(in) :: kdihed0(:)
    real*8, intent(in) :: phi
    real*8, intent(in) :: rsinphi
    real*8, intent(in) :: rijkj
    real*8, intent(in) :: rjkkl
    real*8, intent(out) :: pterm
    character(4), intent(in) :: dihedral_type

    ! result
    real*8 :: gamma

    ! local
    real*8 :: eme
    real*8 :: term, term1, term2, term3, term4
    real*8 :: dphi, dphi1, dphi2
    real*8 :: C0, C1, C2, C3, C4, C5
    real*8 :: A0, A1, A2, A3
    real*8 :: rtwopi

    ! select case( adjustl(molecule(i)%Dihedral_Type(j)) )
    select case(adjustl(dihedral_type))
        case ('cos') ! V = k[1 + cos(n.phi - theta)]
            term  = INT(kdihed0(3)) * phi - kdihed0(1)
            pterm = kdihed0(2) * (1 + COS(term))
            proper_dih = proper_dih + pterm
            gamma = -kdihed0(2) * INT(kdihed0(3)) * SIN(term) * rsinphi * rijkj * rjkkl

        case('harm') ! V = 1/2.k(phi - phi0)²
            rtwopi = 1.d0 / twopi

            dphi  = phi   - kdihed0(1)
            dphi  = dphi  - DNINT(dphi * rtwopi) * twopi
            dphi1 = phi   - kdihed0(3)
            dphi1 = dphi1 - DNINT(dphi1 * rtwopi) * twopi
            dphi2 = phi   - kdihed0(5)
            dphi2 = dphi2 - DNINT(dphi2 * rtwopi) * twopi

            term  =        kdihed0(2) * dphi
            term  = term + kdihed0(4) * dphi1
            term  = term + kdihed0(6) * dphi2

            pterm =         0.5 * term  * dphi
            pterm = pterm + 0.5 * term1 * dphi1
            pterm = pterm + 0.5 * term2 * dphi2

            harm_dih = harm_dih + pterm
            gamma = term * rsinphi * rijkj * rjkkl

        case('hcos') ! V = 1/2.k[cos(phi) - cos(phi0)]²
            dphi  = COS(phi) - COS( kdihed0(2) )
            term  = kdihed0(1) * dphi
            pterm = 0.5d0 * term * dphi
            gamma = term * rijkj * rjkkl

        case('cos3') ! V = 1/2.A1[1 + cos(phi)] + 1/2.A2[1 - cos(2.phi)] + 1/2.A3[1 + cos(3.phi)]
            pterm = 0.5d0 * (kdihed0(1) * (1.d0 + COS(phi))        + &
                             kdihed0(2) * (1.d0 - COS(2.d0 * phi)) + &
                             kdihed0(3) * (1.d0 + COS(3.d0 * phi)))
            gamma = 0.5d0 * (kdihed0(1) * SIN(phi)          - &
                    2.0d0 *  kdihed0(2) * SIN(2.d0 * phi)   + &
                    3.0d0 *  kdihed0(3) * SIN(3.d0 * phi)) * rsinphi * rijkj * rjkkl

        case('ryck') ! V = sum_i^5 Ci.[cos(phi)]^i
            eme   = COS(phi)
            C0    = kdihed0(1) ; C1 = kdihed0(2) ; C2 = kdihed0(3)
            C3    = kdihed0(4) ; C4 = kdihed0(5) ; C5 = kdihed0(6)
            pterm = C0 - (C1 * eme) + (C2 * eme ** 2) - (C3 * eme ** 3) + (C4 * eme ** 4) - (C5 * eme ** 5)
            ryck_dih = ryck_dih + pterm
            gamma = -(-C1 + (2.d0 * C2 * eme) - (3.d0 * C3 * eme ** 2) + (4.d0 * C4 * eme ** 3) - (5.d0 * C5 * eme ** 4)) * rijkj * rjkkl

        case('opls') ! V = A0 + 1/2{A1[1 + cos(phi)] + A2[1 - cos(2.phi)] + A3[1 + cos(3.phi)]}
            dphi  = phi - kdihed0(5)
            A0    = kdihed0(1) ; A1 = kdihed0(2) ; A2 = kdihed0(3)
            A3    = kdihed0(4)
            pterm = A0 + 0.5d0 * (A1 * (1.d0 + COS(dphi)) + A2 * (1.d0 - COS(2.d0 * dphi)) + A3 * (1.d0 + COS(3.d0 * dphi)) )
            gamma = 0.5d0 * ((A1 * SIN(dphi)) - (2.d0 * A2 * SIN(2.d0 * dphi)) + (3.d0 * A3 * SIN(3.d0 * dphi)) ) * rsinphi * rijkj * rjkkl

        case ('chrm')   ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple) Eq. 4.60 (GMX 5.0.5 manual)
            term  = INT(kdihed0(3))  * phi - kdihed0(1)
            term1 = INT(kdihed0(6))  * phi - kdihed0(4)
            term2 = INT(kdihed0(9))  * phi - kdihed0(7)
            term3 = INT(kdihed0(12)) * phi - kdihed0(10)
            term4 = INT(kdihed0(15)) * phi - kdihed0(13)

            pterm =         kdihed0(2)  * (1 + COS(term))
            pterm = pterm + kdihed0(5)  * (1 + COS(term1))
            pterm = pterm + kdihed0(8)  * (1 + COS(term2))
            pterm = pterm + kdihed0(11) * (1 + COS(term3))
            pterm = pterm + kdihed0(14) * (1 + COS(term4))

            proper_dih = proper_dih + pterm

            gamma =       - kdihed0(2)  * INT(kdihed0(3))  * SIN(term)  * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(5)  * INT(kdihed0(6))  * SIN(term1) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(8)  * INT(kdihed0(9))  * SIN(term2) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(11) * INT(kdihed0(12)) * SIN(term3) * rsinphi * rijkj * rjkkl
            gamma = gamma - kdihed0(14) * INT(kdihed0(15)) * SIN(term4) * rsinphi * rijkj * rjkkl

        case('none')
            pterm = 0.d0
            gamma = 0.d0
    end select
end function NOT_GMX


! Calculate the non-bonded intramolecular interactions
!
! @param atom_i[inout]
! @param atom_j[inout]
subroutine NON_BONDED_INTRAMOLECULAR_INTERACTIONS(atom_i, atom_j)
    implicit none

    ! args
    type(MM_atomic), intent(inout) :: atom_i
    type(MM_atomic), intent(inout) :: atom_j

    ! local
    real*8 :: rij(3)

    ! loop index
    integer :: pair_number

    real*8 :: sr2
    real*8 :: sr6
    real*8 :: sr12
    real*8 :: eps
    real*8 :: fs
    real*8 :: rklq
    real*8 :: rklq_sqrt
    real*8 :: charge_i
    real*8 :: charge_j
    real*8 :: expar
    real*8 :: freal
    real*8 :: tterm
    real*8 :: sterm
    real*8 :: KRIJ

    logical :: flag


    rij(:) = DIFFERENCE_BOX(atom_i, atom_j)
    rklq   = SUM(rij(:) ** 2)

    ! Lennard Jones ...
    ! AMBER FF :: GMX COMB-RULE 2
    if (MM%CombinationRule == 2) then
        sr2 = ((atom_i%sig14 + atom_j%sig14) ** 2) / rklq
    ! OPLS  FF :: GMX COMB-RULE 3
    else if (MM%CombinationRule == 3) then
        sr2 = ((atom_i%sig14 * atom_j%sig14) ** 2) / rklq
    end if

    eps = atom_i%eps14 * atom_j%eps14

    ! Nbond_parms directive on ...
    do  pair_number = 1, SIZE(SpecialPairs14)
        flag = GET_FLAGS_1_4(pair_number, atom_i, atom_j)

        if (.NOT. flag) then
            cycle
        end if

        ! AMBER FF :: GMX COMB-RULE 2
        if (MM%CombinationRule == 2) then
            sr2 = ((SpecialPairs14(pair_number)%Parms(1) * 2) ** 2) / rklq

        ! OPLS  FF :: GMX COMB-RULE 3
        else if (MM%CombinationRule == 3) then
            sr2 = (SpecialPairs14(pair_number)%Parms(1) ** 4) / rklq
        end if

        eps = SpecialPairs14(pair_number)%Parms(2) ** 2
        exit
    end do

    rklq_sqrt = SQRT(rklq)
    sr6   = sr2 ** 3
    sr12  = sr6 ** 2
    fs    = 24 * eps * (2 * sr12 - sr6)
    fs    = (fs / rklq) * MM%fudgeLJ
    atom_i%fnonbd14(:) = atom_i%fnonbd14(:) + (fs * rij(:))
    atom_j%fnonbd14(:) = atom_j%fnonbd14(:) - (fs * rij(:))
    ! factor used to compensate the factor1 and factor2 factors ...
    ! factor3 = 1.0d-20
    sterm  = 4 * eps * factor3 * ( sr12 - sr6 )
    ! alternative cutoff formula ...
    ! sterm  = sterm - vscut(ati,atj) + fscut(ati,atj) * ( rklq_sqrt - rcut )
    sterm  = sterm * MM%fudgeLJ

    charge_i = atom_i%charge
    charge_j = atom_j%charge

    !  Real part (Numero de cargas igual ao numero de sitios)
    sr2   = 1 / rklq
    KRIJ  = KAPPA * rklq_sqrt
    expar = EXP(-(KRIJ ** 2))
    freal = coulomb * charge_i * charge_j * (sr2 / rklq_sqrt)
    freal = freal * (ERFC(KRIJ) + 2 * rsqpi * KAPPA * rklq_sqrt * expar) * MM%fudgeQQ

    atom_i%fnonch14(:) = atom_i%fnonch14(:) + (freal * rij(:))
    atom_j%fnonch14(:) = atom_j%fnonch14(:) - (freal * rij(:))

    ! alternative cutoff formula ...
    tterm = (coulomb * factor3 * charge_i * charge_j * ERFC(KRIJ))/ (rklq_sqrt * MM%fudgeQQ)

    LJ_14   = LJ_14   + sterm
    Coul_14 = Coul_14 + tterm
end subroutine NON_BONDED_INTRAMOLECULAR_INTERACTIONS

end module F_intra_m
