! Program for computing two-center overlap integrals with STO

module Overlap_Builder
    use type_m
    use constants_m
    use parameters_m,         only : verbose, PBC, static , SO_coupl , extmagfield
    use PBC_m,                only : Generate_Periodic_Structure
    use Semi_Empirical_Parms, only : atom
    use Allocation_m,         only : DeAllocate_Structures

    public :: Overlap_Matrix

    private

    interface Overlap_matrix
        module procedure Overlap_Matrix
        module procedure Preprocess_OverlapMatrix
    end interface

    ! module variables
    logical,              save :: done = .false.
    logical,              save :: ready = .false.
    logical,              save :: Pulay = .false.
    integer, allocatable, save :: ija(:)
    real*8,  allocatable, save :: a_xyz(:,:)
    real*8,  allocatable, save :: b_xyz(:,:)
    real*8,  allocatable, save :: S_pack(:)
contains


subroutine OVERLAP_MATRIX(system, basis, S_matrix, purpose, site)
    implicit none

    ! args
    type(structure),               intent(in)    :: system
    type(STO_basis),               intent(in)    :: basis(:)
    real*8,           allocatable, intent(inout) :: S_matrix(:,:)
    character(len=*), optional,    intent(in)    :: purpose
    integer,          optional,    intent(in)    :: site

    ! local
    type(structure)              :: pbc_system
    type(STO_basis), allocatable :: pbc_basis(:)
    real*8         , allocatable :: S_tmp(:,:)
    integer                      :: NonZero, S_size , n
    real*8                       :: Sparsity

    CALL util_overlap

    S_size = SUM(atom(system%AtNo)%DOS, system%QMMM == "QM")

    allocate( S_tmp( S_size , S_size ) )

    select case (purpose)
        case('FMO')
            CALL Generate_Periodic_Structure( system, pbc_system, pbc_basis )
            CALL Build_Overlap_Matrix(system, basis, pbc_system, pbc_basis, S_tmp)

        case('GA-CG')
            ! if no PBC pbc_system = system ; do NOT use OPT_parms
            CALL Generate_Periodic_Structure(system, basis, pbc_system, pbc_basis)
            CALL Build_Overlap_Matrix(system, basis, pbc_system, pbc_basis, S_tmp)

        case('Pulay') ! <== used by diagnostic through Hellman_Feynman_Pulay function
            ! if no PBC pbc_system = system
            CALL Generate_Periodic_Structure(system, pbc_system, pbc_basis)
            CALL Pulay_Overlap(system, basis, pbc_system, pbc_basis, S_tmp, site)

        case default
            ! Overlap Matrix S(a,b) of the system
            if (verbose) then
                Print 53

                if (product(2 * PBC(:) + 1) /= 1) then
                    Print 51, PBC(:)
                else
                    Print 54
                end if
            end if

            ! if no PBC pbc_system = system
            CALL Generate_Periodic_Structure( system, pbc_system, pbc_basis )
            CALL Build_Overlap_Matrix( system, basis(1:S_size), pbc_system, pbc_basis, S_tmp , recycle = .true. )

            NonZero  = count(S_tmp /= 0.d0)
            Sparsity = float(NonZero) / float(S_size ** 2)

            if (verbose) then
                Print 69, Sparsity
                Print 55
                print*, '>> Overlap done <<'
            end if
    end select

    CALL Deallocate_Structures(pbc_system)
    if (allocated(pbc_basis)) then
        deallocate(pbc_basis)
    end if

    if( not(SO_coupl) .AND. not(extmagfield) ) then

        if( .NOT. allocated(S_matrix) ) allocate( S_matrix( S_size , S_size ) )

        S_matrix = S_tmp

    else

        S_size = 2 * S_size

        if( .NOT. allocated(S_matrix) ) allocate( S_matrix( S_size , S_size ) , source = 0.0d0 )

        n = S_size / 2

        S_matrix(     1 :      n ,     1 :      n ) = S_tmp( 1 : n , 1 : n )
        S_matrix( n + 1 : S_size , n + 1 : S_size ) = S_tmp( 1 : n , 1 : n )

    end if

    deallocate( S_tmp )
    
    include 'formats.h'
end subroutine OVERLAP_MATRIX


subroutine PREPROCESS_OVERLAPMATRIX(system, basis)
    implicit none

    ! args
    type(structure), intent(in) :: system
    type(STO_basis), intent(in) :: basis(:)

    ! local
    real*8,          allocatable :: S_matrix(:,:)
    integer                      :: S_size
    type(structure)              :: pbc_system
    type(STO_basis), allocatable :: pbc_basis(:)

    ! this subroutine simply instantiate the ORIGINAL Overlap Matrix for future use
    ! called by HFP force routines
    CALL util_overlap

    S_size = SUM(atom(system%AtNo)%DOS , system%QMMM == "QM")
    allocate(S_matrix(S_size,S_size))

    Pulay = .true.
    done = .false.

    ! if no PBC pbc_system = system
    CALL Generate_Periodic_Structure(system, pbc_system, pbc_basis)
    CALL Build_Overlap_Matrix(system, basis, pbc_system, pbc_basis, S_matrix)
    CALL Deallocate_Structures(pbc_system)

    if (allocated(pbc_basis)) then
        deallocate(pbc_basis)
    end if
    deallocate(S_matrix)
end subroutine PREPROCESS_OVERLAPMATRIX


subroutine BUILD_OVERLAP_MATRIX(b_system, b_basis, a_system, a_basis, S_matrix, recycle)
    use util_m , factorial => fact

    implicit none

    ! args
    type(structure),   intent(in)    :: b_system
    type(STO_basis),   intent(in)    :: b_basis(:)
    type(structure),   intent(in)    :: a_system
    type(STO_basis),   intent(in)    :: a_basis(:)
    real*8,            intent(inout) :: S_matrix(:,:)
    logical, optional, intent(in)    :: recycle

    ! local
    real*8  :: expa, expb, Rab, aux, anor
    real*8  :: sux(0:10)
    integer :: a, b, ia, ib, ja, jb
    integer :: na, la, ma
    integer :: nb, lb, mb
    integer :: msup, i, j, k, m
    logical :: motion_detector_ready, atom_not_moved

    ! local parameters
    integer, parameter :: mxn = 15
    integer, parameter :: mxl = 5

    ! local arrays
    real*8, dimension(0:mxl)                   :: solvec
    real*8, dimension(0:mxl)                   :: solnorm, sol_partial
    real*8, dimension(-mxl:mxl,-mxl:mxl,0:mxl) :: rl, rl2

    motion_detector_ready = present(recycle) .AND. ready

    S_matrix = 0

    !$OMP parallel do &
    !$OMP   default(shared) &
    !$OMP   schedule(dynamic, 1) &
    !$OMP   private(ib, ia, atom_not_moved, Rab, jb, ja, b, a, k, nb, na, lb, &
    !$OMP           la, mb, ma, aux, msup, solnorm, j, i, expb, expa, solvec, &
    !$OMP           anor, m, sol_partial, sux, rl, rl2)
    do ib = 1, b_system%atoms
        do ia = ib + 1, a_system%atoms

            if ((b_system%QMMM(ib) /= "QM") .OR. (a_system%QMMM(ia) /= "QM")) then
                cycle
            end if

            ! ckecking: if atoms ia and ib remain fixed => recover S_matrix
            if (motion_detector_ready) then

                Rab = GET_RAB(a_system%coord(ia,:), b_system%coord(ib,:))
                if (Rab > cutoff_Angs) then
                    cycle
                end if

                atom_not_moved = ALL(a_system%coord(ia,:) == a_xyz(ia,:))
                atom_not_moved = atom_not_moved .AND. ALL(b_system%coord(ib,:) == b_xyz(ib,:))
                if (atom_not_moved) then

                    do jb = 1, atom(b_system%AtNo(ib))%DOS
                        do ja = 1, atom(a_system%AtNo(ia))%DOS
                            b = b_system%BasisPointer(ib) + jb
                            a = a_system%BasisPointer(ia) + ja

                            a = a - (a_basis(a)%copy_No * size(b_basis))

                            if (S_matrix(a,b) /= D_zero) then
                                cycle
                            end if

                            do k = ija(a), ija(a + 1) - 1
                                if (ija(k) == b) then
                                    S_matrix(a,b) = S_pack(k)
                                end if
                            end do
                        end do
                    end do

                    cycle
                end if
            end if

            ! if atoms ia and ib move => calculate rotation matrix for the highest l
            call RotationOverlap( b_system, a_system, ia, ib, Rab, rl, rl2 )

            if (Rab > cutoff_Angs) then
                cycle
            end if

            do jb = 1, atom(b_system%AtNo(ib))%DOS
                do ja = 1, atom(a_system%AtNo(ia))%DOS
                    b = b_system%BasisPointer(ib) + jb
                    a = a_system%BasisPointer(ia) + ja

                    nb = b_basis(b)%n
                    lb = b_basis(b)%l
                    mb = b_basis(b)%m

                    na = a_basis(a)%n
                    la = a_basis(a)%l
                    ma = a_basis(a)%m

                    aux = 1.d0 / (factorial(na + na) * factorial(nb + nb))
                    msup = MIN(la, lb)

                    solnorm(0:msup) = 0

                    do j = 1, b_basis(b)%Nzeta
                        do i = 1, a_basis(a)%Nzeta
                            expb = b_basis(b)%zeta(j)
                            expa = a_basis(a)%zeta(i)

                            ! OVERLAP SUBROUTINE: lined-up frame
                            call solap(na, la, expa, nb, lb, expb, Rab, solvec)

                            anor = GET_ANOR(expa, expb, na, nb, la, lb, aux)

                            ! Introduces normalization of the STO in the integrals
                            do m = 0, msup - 1
                                sol_partial(m) = solvec(m) * anor
                                anor = anor / GET_ANOR_DIVISOR(la, lb, m)
                                if (m == 0) then
                                    anor = anor + anor
                                end if
                            end do
                            sol_partial(msup) = solvec(msup) * anor

                            forall(k=0:msup) solnorm(k) = solnorm(k) + a_basis(a)%coef(i) * b_basis(b)%coef(j) * sol_partial(k)
                        end do
                    end do

                    ! Rotation of overlap integrals to the molecular frame
                    sux(0) = solnorm(0) * rl(ma,0,la) * rl(mb,0,lb)
                    forall(k=1:msup) sux(k) = solnorm(k) * (rl(ma,-k,la) * rl(mb,-k,lb) + rl(ma,k,la) * rl(mb,k,lb))

                    a = a - a_basis(a)%copy_No * size(b_basis)
                    S_matrix(a,b) = S_matrix(a,b) + SUM(sux(0:msup))
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    ! diagonal elements
    forall(i=1:size(b_basis)) S_matrix(i,i) = 1

    ! save overlap data for reuse
    if ((.NOT. static .OR. Pulay) .AND. (.NOT. done)) then

        if(allocated(a_xyz)) then
            a_xyz = a_system%coord
            b_xyz = b_system%coord
        else
            allocate(a_xyz, source = a_system%coord)
            allocate(b_xyz, source = b_system%coord)
        end if

        CALL sprsin(S_matrix)

        done = .true.
        ready = .true.
    end if

    ! symmetric overlap matrix
    forall(i=1:size(b_basis)) S_matrix(i,i:size(b_basis)) = S_matrix(i:size(b_basis),i)

end subroutine BUILD_OVERLAP_MATRIX


subroutine PULAY_OVERLAP(b_system, b_basis, a_system, a_basis, S_matrix, site)
    use util_m , factorial => fact

    implicit none

    ! args
    type(structure), intent(in)    :: b_system
    type(STO_basis), intent(in)    :: b_basis(:)
    type(structure), intent(in)    :: a_system
    type(STO_basis), intent(in)    :: a_basis(:)
    real*8,          intent(inout) :: S_matrix(:,:)
    integer,         intent(in)    :: site

    ! local
    real*8  :: expa, expb, Rab, aux, anor
    real*8  :: sux(0:10)
    integer :: a, b, ia, ib, ja, jb
    integer :: na, la, ma
    integer :: nb, lb, mb
    integer :: msup, i, j, k, m
    logical :: atom_not_moved

    ! local parameters
    integer , parameter :: mxn = 15
    integer , parameter :: mxl = 5

    ! local arrays
    real*8 :: solvec(0:mxl)
    real*8 :: solnorm(0:mxl)
    real*8 :: sol_partial(0:mxl)
    real*8 :: rl(-mxl:mxl,-mxl:mxl,0:mxl)
    real*8 :: rl2(-mxl:mxl,-mxl:mxl,0:mxl)

    ib = site

    !$OMP parallel do schedule(static) &
    !$OMP   default(shared) private(ia, jb, ja, a, b)
    do ia = ib , a_system%atoms
        do jb = 1, atom(b_system%AtNo(ib))%DOS
            do ja = 1, atom(a_system%AtNo(ia))%DOS
                b = b_system%BasisPointer(ib) + jb
                a = a_system%BasisPointer(ia) + ja
                a = a - a_basis(a)%copy_No * size(b_basis)
                S_matrix(a,b) = 0
            end do
        end do
    end do
    !$OMP end parallel do

    ! The parallel region below is a very complicated one. It takes many public
    ! matrices, for reading, and performs lots of writes in a single, shared,
    ! matrix based on the values read from the other matrices. The order of
    ! these writes is very confusing and hard to predict.
    ! This causes race conditions when not using dynamic(1) schedule.
    ! The use of ATOMIC operations fixes the problem. But it is not a solution
    ! for concurrency problems that we face. And it makes the code slower.

    !$OMP parallel do &
    !$OMP default(shared) schedule(dynamic,1) &
    !$OMP private(ia, atom_not_moved, Rab, jb, ja, b, a, k, nb, na, lb, la, &
    !$OMP         mb, ma, aux, msup, solnorm, j, i, expb, expa, solvec, anor, &
    !$OMP         m, sol_partial, sux, rl, rl2)
    do ia = ib + 1, a_system%atoms

        if (a_system%QMMM(ia) /= "QM") then
            cycle
        end if

        ! ckecking: if atoms ia and ib remain fixed
        Rab = GET_RAB(a_system%coord(ia,:), b_system%coord(ib,:))
        if (Rab > cutoff_Angs) then
            cycle
        end if

        atom_not_moved = ALL(a_system%coord(ia,:) == a_xyz(ia,:))
        atom_not_moved = atom_not_moved .AND. ALL(b_system%coord(ib,:) == b_xyz(ib,:))
        if (atom_not_moved) then
         print*, "if I ever stop here there is a bug because atom b has always moved before calling this routine (18/09/20)" ;stop
            cycle
        end if

        ! if atoms ia and ib move => calculate rotation matrix for the highest l
        call ROTATIONOVERLAP(b_system, a_system, ia, ib, Rab, rl, rl2)

        do jb = 1, atom(b_system%AtNo(ib))%DOS
            do ja = 1, atom(a_system%AtNo(ia))%DOS
                b = b_system%BasisPointer(ib) + jb
                a = a_system%BasisPointer(ia) + ja

                nb = b_basis(b)%n
                lb = b_basis(b)%l
                mb = b_basis(b)%m

                na = a_basis(a)%n
                la = a_basis(a)%l
                ma = a_basis(a)%m

                aux = 1.d0 / (factorial(na * 2) * factorial(nb * 2))
                msup = MIN(la, lb)

                solnorm(0:msup) = 0.d0

                do j = 1, b_basis(b)%Nzeta
                    do i = 1, a_basis(a)%Nzeta
                        expb = b_basis(b)%zeta(j)
                        expa = a_basis(a)%zeta(i)

                        ! OVERLAP SUBROUTINE: lined-up frame
                        call SOLAP(na, la, expa, nb, lb, expb, Rab, solvec)

                        anor = GET_ANOR(expa, expb, na, nb, la, lb, aux)

                        ! Introduces normalization of the STO in the integrals
                        do m = 0, msup - 1
                            sol_partial(m) = solvec(m) * anor
                            anor = anor / GET_ANOR_DIVISOR(la, lb, m)
                            if (m == 0) then
                                anor = anor * 2
                            end if
                        end do
                        sol_partial(msup) = solvec(msup) * anor

                        forall(k=0:msup) solnorm(k) = solnorm(k) + a_basis(a)%coef(i) * b_basis(b)%coef(j) * sol_partial(k)
                    end do
                end do

                ! Rotation of overlap integrals to the molecular frame
                sux(0) = solnorm(0) * rl(ma,0,la) * rl(mb,0,lb)
                forall(k=1:msup) sux(k) = solnorm(k) * (rl(ma,-k,la) * rl(mb,-k,lb) + rl(ma,k,la) * rl(mb,k,lb))

                a = a - a_basis(a)%copy_No * size(b_basis)
                ! symmety in the block of atoms a and b 
                !$OMP atomic
                S_matrix(a,b) = S_matrix(a,b) + SUM(sux(0:msup))
                !$OMP atomic
                S_matrix(b,a) = S_matrix(b,a) + SUM(sux(0:msup))
            end do
        end do
    end do
    !$OMP end parallel do

    ! diagonal elements ...
    b = size(b_basis)
    !$OMP parallel do private(i) shared(S_matrix, b)
    do i = 1, b
        S_matrix(i,i) = 1
    end do
    !$OMP end parallel do

end subroutine PULAY_OVERLAP


pure function GET_ANOR_DIVISOR(la, lb, m) result(divisor)
    implicit none

    ! args
    integer, intent(in) :: la, lb
    integer, intent(in) :: m

    ! result
    real*8 :: divisor

    divisor = (la + m + 1.d0) * (la - m) * (lb + m + 1.d0) * (lb - m)
    divisor = DSQRT(divisor)
end function GET_ANOR_DIVISOR


pure function GET_RAB(a_coord, b_coord) result(rab)
    implicit none

    ! args
    real*8, intent(in) :: a_coord(:)
    real*8, intent(in) :: b_coord(:)

    ! result
    real*8 :: rab

    rab = SUM((a_coord - b_coord) ** 2)
    rab = SQRT(rab)
end function GET_RAB


pure function GET_ANOR(expa, expb, na, nb, la, lb, aux) result(anor)
    implicit none

    ! args
    real*8,  intent(in) :: expa, expb
    real*8,  intent(in) :: aux
    integer, intent(in) :: na, nb
    integer, intent(in) :: la, lb

    ! result
    real*8 :: anor

    ! parameter
    real*8, parameter :: PI = 4 * DATAN(1.d0)

    anor =        (expa * 2) ** (na * 2 + 1)
    anor = anor * (expb * 2) ** (nb * 2 + 1)
    anor = anor * (la * 2 + 1.d0)
    anor = anor * (lb * 2 + 1.d0)
    anor = anor * aux
    anor = DSQRT(anor) / (4 * PI)
end function GET_ANOR


subroutine ROTATIONOVERLAP(b_system, a_system, ia, ib, Rab, rl, rl2)
    implicit none

    ! args
    type(structure), intent(in)  :: b_system
    type(structure), intent(in)  :: a_system
    integer,         intent(in)  :: ia, ib
    real*8,          intent(out) :: Rab
    real*8,          intent(out) :: rl(:,:,:)
    real*8,          intent(out) :: rl2(:,:,:)

    ! local
    real*8  :: xa, ya, za, xb, yb, zb, xab, yab, zab, xy
    real*8  :: sinal, cosal, sinbet, cosbet, singa, cosga
    integer :: AtNo_a, AtNo_b
    integer :: la_max, lb_max, lmax

    ! local parameters
    real*8,  parameter :: tol = 1.d-10
    integer, parameter :: mxn = 15
    integer, parameter :: mxl = 5

    ! INPUT DATA
    ! na, la, expa, xa, ya, za:  N, L, EXPONENT and cartesian coordinates of function on A
    ! nb, lb, expb, xb, yb, zb:  N', L', EXPONENT' and cartesian coordinates of function on B
    AtNo_a = a_system%AtNo(ia)
    AtNo_b = b_system%AtNo(ib)

    xa = a_system%coord(ia,1)
    ya = a_system%coord(ia,2)
    za = a_system%coord(ia,3)
    xb = b_system%coord(ib,1)
    yb = b_system%coord(ib,2)
    zb = b_system%coord(ib,3)

    ! Loads the matrices for the rotation of (normalized) spherical harmonics
    ! from the lined-up system to the molecular frame
    ! Rotation matrices are stored in rl(m,m',l) with  (-l <= m,m' <= l)
    la_max = atom(AtNo_a)%AngMax
    lb_max = atom(AtNo_b)%AngMax

    xab = xb - xa
    yab = yb - ya
    zab = zb - za

    xy  = DSQRT(xab ** 2 + yab ** 2)
    if (xy > tol) then
        sinal = yab / xy
        cosal = xab / xy
    else
        sinal = 0
        cosal = 1
    endif

    Rab = DSQRT(xab ** 2 + yab ** 2 + zab ** 2)
    sinbet = xy / Rab
    cosbet = zab / Rab

    singa  = 0
    cosga  = 1
    lmax = MAX(la_max, lb_max)

    call ROTAR(lmax, mxl, cosal, sinal, cosbet, sinbet, cosga, singa, rl2, rl)
end subroutine ROTATIONOVERLAP


subroutine SPRSIN(S_matrix)
    implicit none

    ! args
    real*8 , intent(in) :: S_matrix(:,:)

    !=====================================================================
    ! converts a square matrix a(1:n,1:n) with physical dimension np into
    ! row-indexed sparse storage mode. Only elements of a with magnitude >=
    ! thresh are retained. Output is in two linear arrays with physical
    ! dimension nmax (an input parameter): sa(1:) contains array lines,
    ! indexed by ija(1:). The logical sizes of sa and ija on output are
    ! both  ija(ija(1)-1)-1
    !=====================================================================

    ! local
    integer                      :: i , j , k , n , length, nnz
    real*8                       :: thresh = mid_prec
    integer , allocatable , save :: tmp_ij(:)
    real*8  , allocatable , save :: tmp_S(:)
    logical               , save :: NOT_been_here = .true.

    !local parameters
    integer  :: nmax

    n = size(S_matrix(:,1))
    nmax = (n * (n + 1)) / 2 + 1

    if (NOT_been_here) then
        allocate(tmp_ij(nmax), ija(nmax))
        allocate(tmp_S(nmax), S_pack(nmax))
    end if

    tmp_ij = I_Zero
    tmp_S = D_zero

    do j = 1, n
        tmp_S(j) = S_matrix(j,j)
    end do

    tmp_ij(1) = n + 2
    k = n + 1
    do i = 1, n
        do j = 1, i
            if (ABS(S_matrix(i,j)) >= thresh) then
                if (i /= j) then
                    k = k + 1

                    if (k > nmax) then
                        pause 'nmax too small in sprsin'
                    end if

                    tmp_S(k) = S_matrix(i,j)
                    tmp_ij(k) = j
                end if
            end if
        end do
        tmp_ij(i + 1) = k + 1
    end do

    length = tmp_ij(tmp_ij(1) - 1) - 1          ! <== physical size of S_pack and ija
    nnz = tmp_ij(tmp_ij(1) - 1) - 2             ! <== number of nonzero elements

    ! prepare S_pack and ija for new data
    deallocate( S_pack , ija )
    allocate(S_pack(length), source = tmp_S (1:length) )
    allocate(ija(length)   , source = tmp_ij(1:length) )

    NOT_been_here = .false.
end subroutine SPRSIN


!-----------------  END OF MAIN     -------------------------
!     Subr0utine for overlap integrals:
!        < N L M | N' L' M >    M >= 0
!     with STO functions in a lined-up system
!     Starts with the basic overlap integrals:
!        < N 0 0 | N' 0 0 >
!     and uses the following recurrence relations:
!
!     For increasing  M :
!
!     < N M+1 M+1 | N' M+1 M+1 > = ((2M+1)**2/(4R**2)) * (
!             2 < N+1 M M | N'+1 M M > + 2 R**2 < N+1 M M | N'-1 M M >
!             + 2 R**2 < N-1 M M | N'+1 M M > - < N+3 M M | N'-1 M M >
!             - < N-1 M M | N'+3 M M > - R**4 < N-1 M M | N'-1 M M > )
!
!
!     For increasing  L:
!
!     < N L+1 M | N' L' M > = ((2L+1)/(2 (L-M+1) R) ) * (
!             < N+1 L M | N' L' M >  + R**2 < N-1 L M | N' L' M >
!             - < N-1 L M | N'+2 L' M > )
!          - ( (L+M)/(L-M+1) ) < N L-1 M | N' L' M >
!
!
!     For increasing  L':
!
!     < N L M | N' L'+1 M > = -((2L'+1)/(2 (L'-M+1) R) ) * (
!             < N L M | N'+1 L' M >  + R**2 < N L M | N'-1 L' M >
!             - < N+2 L M | N'-1 L' M > )
!          - ( (L'+M)/(L'-M+1) ) < N L M | N' L'-1 M >
!
!     The starting integrals are computed with the equation
!             (Jaime Fernandez Rico):
!     < N 0 0 | N' 0 0 > = ( 4 pi Exp[-z*R] N! N'! / (z+z')**(N+N'+1) )
!         * Sum[ Sum[
!            (N+N'-j-j')! [(z+z')*R]**(j+j') * 1F1[j'+1,j+j'+2,(z-z')R]
!            / ( (N-j)! (N'-j')! (j+j'+1)!
!         , {j',0,N'}] , {j,0,N}]
!
!     where  1F1[a,b,z] is the confluent hypergeometric function
!
!
!**********************************************************************
subroutine SOLAP(na, la, za, nb, lb, zb, R, sol)
    implicit real*8  (a-h,o-z)

    real*8 , parameter :: pi4 = 12.566370614359172953850573533118012d0
    common /const/ re(0:1000), reali(0:1000), fact(0:150) , facti(0:150)
    integer, parameter :: mxn = 15, mxl = 5
    integer, parameter :: mxh = mxn+3*mxl
    integer            :: na, la, nb, lb
    integer            :: nmin, npmin, npmax, ngmax, namax, iz, ia, imax, i, j, kint, n, np, k, m, ll, nmax, mtop
    real*8             :: zzpr, zzpi, zzpin, zzpinp, auxk, fact, sumk, sumj, facti, smat, smataux, smatm
    real*8             :: Rinv, R2, R4, updltm0i, saux, sbux, f11a
    real*8             :: za, zb, R, sol, z, zp, expzR, arg, auxa, re, reali, auxb, f11b, h1f1, aux
    common / kintcom / kint
    dimension h1f1(0:mxh+1,0:2*(mxh+1)), sol(0:mxl), smat(0:mxh,0:mxh)
    dimension smataux(0:mxh,0:mxh), smatm(0:mxh,0:mxh)

    if (za >= zb) then
        nmin = na - la
        nmax = na + la + lb + lb
        npmin = nb - lb
        npmax = nb + lb + la + la
        z = za
        zp = zb
    else
        nmin = nb - lb
        nmax = nb + lb + la + la
        npmin = na - la
        npmax = na + la + lb + lb
        z = zb
        zp = za
    endif
    expzR = dexp(-z*R)

    ! Computes the Confluent Hypergeometric functions (h1f1)
    ngmax = nmax + npmax + 2   ! Top value of gamma parameter
    namax = npmax + 1        ! Top value of alpha parameter
    arg = (z - zp) * R       ! Argument of the hypergeometric
    iz = z
    ia = (ngmax - iz) / 2

    if (ia < 1) then
        ia = 1
    end if

    if (ia >= namax) then
        ia = namax-1
    end if

    ! Computes 1F1[ia,ngmax,z]  and  1F1[ia+1,ngmax,z] by their series
    auxa = arg * re(ia) * reali(ngmax)
    f11a = 1.d0
    auxb = arg * re(ia + 1) * reali(ngmax)
    f11b = 1.d0
    imax = 1000

    do i = 1, imax
        f11a = f11a + auxa
        f11b = f11b + auxb
        auxa = auxa * re(ia + i)     * arg * reali(ngmax + i) * reali(i + 1)
        auxb = auxb * re(ia + i + 1) * arg * reali(ngmax + i) * reali(i + 1)
        if ( auxb < 1.d-35 * f11b) then
           go to 100
        end if
    end do

    100 h1f1(ia,ngmax) = f11a
    h1f1(ia + 1,ngmax) = f11b

    !     Computes the column with gamma = ngmax by the recursion:
    !        a * 1F1[a+1,g,x] = (x-2a-g) * 1F1[a,g,x] + (g-a) * 1F1[a-1,g,x]
    !     which is stable upwards for  a > (g - x) / 2
    !               and downwards for a < (g - x) / 2

    ! Upwards recursion:
    aux = arg - re(ngmax)
    do i = ia + 1, namax - 1
        h1f1(i + 1,ngmax) = ((aux + re(i + i)) * h1f1(i,ngmax) + re(ngmax - i) * h1f1(i - 1,ngmax)) * reali(i)
    end do

    ! Downwards recursion:
    do i = ia, 2, -1
        h1f1(i - 1,ngmax) = (re(i) * h1f1(i + 1,ngmax) - (aux + re(i + i)) * h1f1(i,ngmax)) * reali(ngmax - i)
    end do



    !     Computes the remaining columns (gamma < ngmax) by the recursion:
    !        g * 1F1[a,g,x] = a * 1F1[a+1,g+1,x] + (g-a) * 1F1[a,g+1,x]
    !     for the rows: a = 1,  and
    !        1F1[a+1,g,x]= (x/g) * 1F1[a+1,g+1,x] + 1F1[a,g,x]
    !     for the rows: 2 <= a < ngmax
    !     which are stable as they are used
    do j = ngmax - 1, 2, -1
        h1f1(1,j) = (h1f1(2,j + 1) + re(j - 1) * h1f1(1,j + 1)) * reali(j)
        aux = arg * reali(j)
        do i = 1, MIN(namax - 1,j - 2)
            h1f1(i + 1,j) = aux * h1f1(i + 1,j + 1) + h1f1(i,j)
        end do
    end do

    ! Computes the basic integrals  < N 0 0 | N' 0 0 >
    kint = 0

    if (za .ge. zb) then
        zzpR = (z + zp) * R
        zzpi = 1.d0 / (z + zp)
        zzpin = zzpi ** (nmin + npmin + 1)

        do n = nmin, nmax
            zzpinp = zzpin
            do np = npmin, npmax - n + nmin
                auxk = fact(n + np) * zzpinp
                sumk = 0.d0
                do k = 0, n + np
                    sumj = 0.d0
                    do j = max(0,k - n), min(k,np)
                        sumj = sumj + h1f1(j + 1,k + 2) * facti(n + j - k) * facti(np - j)
                    end do
                    sumk = sumk + sumj * auxk
                    auxk = auxk * zzpR * reali(n + np - k) * reali(k + 2)
                end do
                smat(n,np) = pi4 * expzR * fact(n) * fact(np) * sumk
                kint = kint + 1
                smataux(n,np) = smat(n,np)
                smatm(n,np) = smat(n,np)
                zzpinp = zzpinp * zzpi
            end do
            zzpin = zzpin * zzpi
        end do
    else
        zzpR = (z + zp) * R
        zzpi = 1.d0 / (z + zp)
        zzpin = zzpi ** (nmin + npmin + 1)
        do n = nmin, nmax
            zzpinp = zzpin
            do np = npmin, npmax - n + nmin
                auxk = fact(n + np) * zzpinp
                sumk = 0.d0
                do k = 0, n + np
                    sumj = 0.d0
                    do j = max(0,k - n), min(k,np)
                        sumj = sumj + h1f1(j + 1,k + 2) * facti(n + j - k) * facti(np - j)
                    end do
                    sumk = sumk + sumj * auxk
                    auxk = auxk * zzpR * reali(n + np - k) * reali(k + 2)
                end do
                smat(np,n) = pi4 * expzR * fact(n) * fact(np) * sumk
                kint = kint + 1
                smataux(np,n) = smat(np,n)
                smatm(np,n) = smat(np,n)
                zzpinp = zzpinp * zzpi
            end do
            zzpin = zzpin * zzpi
        end do
    end if

    ! Redefines the limits of the indices taking into account that
    ! matrix smat is stored as  smat(na,nb)
    nmin = na - la
    nmax = na + la + lb + lb
    npmin = nb - lb
    npmax = nb + lb + la + la

    Rinv = 1.d0 / R
    R2 = R * R
    R4 = R2 * R2

    updltm0i = .5d0
    mtop = MIN(la, lb)
    do m = 0, mtop

        ! Recurrence on LA
        if (la > m) then
            do np = npmin, npmax - 2
                saux = smat(nmin,np)
                do n = nmin + 1, nmax - 1 - np + npmin
                    sbux = re(m + m + 1) * .5d0 * Rinv * (smat(n + 1,np) + R2 * saux - smat(n - 1,np + 2))
                    saux = smat(n,np)
                    smataux(n,np) = saux
                    smat(n,np) = sbux
                    kint = kint + 1
                end do
            end do

            do ll = m + 1, la - 1
                do np = npmin, npmax - 2 * (ll + 1 - m)
                    saux = smat(nmin + ll - m,np)
                    do n = nmin + 1 + ll - m, nmax - 1 - np + npmin - ll + m
                        sbux = re(ll + ll + 1) * .5d0 * reali(ll - m + 1) * Rinv  * (smat(n + 1,np) + R2 * saux - smat(n - 1,np + 2)) &
                        - re(ll + m) * reali(ll - m + 1) * smataux(n,np)
                        saux = smat(n,np)
                        smataux(n,np) = saux
                        smat(n,np) = sbux
                        kint = kint + 1
                    end do
                end do
            end do
        end if

        ! Recurrence on LB
        if (lb > m) then
            do n = nmin + la - m, nmax - la + m - 2
                saux = smat(n,npmin)
                do np = npmin + 1, npmax - 2 * (la - m) - n + nmin + la - m - 1
                    sbux = -re(m + m + 1) * .5d0 * Rinv * (smat(n,np + 1) + R2 * saux - smat(n + 2,np - 1))
                    saux = smat(n,np)
                    smataux(n,np) = saux
                    smat(n,np) = sbux
                    kint = kint + 1
                end do
            end do

            do ll = m + 1, lb - 1
                do n = nmin + la - m, nmax - la + m - 2 * (ll + 1 - m)
                    saux = smat(n,npmin + ll - m)
                    do np = npmin + 1 + ll - m,npmax - 2 * (la - m) - n + nmin + la - m - ll + m - 1
                        sbux = -re(ll + ll + 1) * .5d0 * reali(ll - m + 1) * Rinv * (smat(n,np + 1) + R2 * saux - smat(n + 2,np - 1)) &
                        - re(ll + m) * reali(ll - m + 1) * smataux(n,np)
                        saux = smat(n,np)
                        smataux(n,np) = saux
                        smat(n,np) = sbux
                        kint = kint + 1
                    end do
                end do
            end do
        end if

        sol(m) = smat(na,nb)

        ! Recurrence on M

        if (m .eq. mtop) then
            go to 200
        end if

        nmin = nmin + 1
        nmax = nmax - 3
        npmin = npmin + 1
        npmax = npmax - 3

        aux = updltm0i * .25d0 * re(m + m + 1) * re(m + m + 1) * Rinv * Rinv

        do n = nmin, nmax
            do np = npmin, npmax - n + nmin
                smataux(n,np) = aux * (2.d0 * (smatm(n + 1,np + 1) + R2 * (smatm(n + 1,np - 1) + smatm(n - 1,np + 1))) &
                - smatm(n + 3,np - 1) - smatm(n - 1,np + 3) - R4 * smatm(n - 1,np - 1))
            end do
        end do

        do n = nmin, nmax
            do np = npmin, npmax - n + nmin
                smat(n,np) = smataux(n,np)
                kint = kint + 1
                smatm(n,np) = smataux(n,np)
            end do
        end do
        updltm0i = 1.d0
    end do

    200 continue
    return
end subroutine SOLAP


!***********************************************************************
!                                                                      *
!   subroutine rotar                                                   *
!                                                                      *
!   this subroutine yields the rotation matrices r(m',m;l) that are    *
!   necessary to perform a coordinate transformation used to align     *
!   two sets of real spherical harmonics centered at different points  *
!   (a and b).                                                         *
!   this transformation converts each original real spherical harmonic *
!   in a linear combination of the real spherical harmonics with the   *
!   same l and different m.                                            *
!   the maximum value for the orbital quantum number l is 12, to extend*
!   this program to greater values of l it is necessary to extend the  *
!   common sqroot (needed in the subroutine dlmn) with the values of   *
!   the square roots of the first 2*ltot+1 integers and their          *
!   reciprocals.                                                       *
!                                                                      *
!***********************************************************************
subroutine ROTAR(lmax, ltot, cosal, sinal, cosbet, sinbet, cosga, singa, dl, rl)
    implicit real*8 (a-h,o-z)

    integer :: ltot, lmax, l, l1
    real*8  :: root2, zero, one, pt5, cosal, sinal, cosbet, sinbet, cosga, singa, dl, rl
    real*8  :: cosag, cosamg, sinag, sinamg, tgbet2

    dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
    dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
    data root2/1.414213562373095d0/,zero/0.0d0/,one/1.0d0/,pt5/0.5d0/

    !
    !     computation of the initial matrices d0, r0, d1 and r1
    !
    dl(0,0,0)  = one
    rl(0,0,0)  = one

    if (lmax == 0) then
        go to 201
    end if

    dl(1,1,1)  = (one+cosbet) * pt5
    dl(1,0,1)  =-sinbet/root2
    dl(1,-1,1) = (one-cosbet) * pt5
    dl(0,1,1)  =-dl(1,0,1)
    dl(0,0,1)  = dl(1,1,1)-dl(1,-1,1)
    dl(0,-1,1) = dl(1,0,1)
    dl(-1,1,1) = dl(1,-1,1)
    dl(-1,0,1) = dl(0,1,1)
    dl(-1,-1,1)= dl(1,1,1)
    cosag  = cosal * cosga - sinal * singa
    cosamg = cosal * cosga + sinal * singa
    sinag  = sinal * cosga + cosal * singa
    sinamg = sinal * cosga - cosal * singa
    rl(0,0,1)  = dl(0,0,1)
    rl(1,0,1)  = root2 * dl(0,1,1) * cosal
    rl(-1,0,1) = root2 * dl(0,1,1) * sinal
    rl(0,1,1)  = root2 * dl(1,0,1) * cosga
    rl(0,-1,1) =-root2 * dl(1,0,1) * singa
    rl(1,1,1)  = dl(1,1,1) * cosag - dl(1,-1,1) * cosamg
    rl(1,-1,1) =-dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,1,1) = dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,-1,1)= dl(1,1,1) * cosag + dl(1,-1,1) * cosamg

    !     the remaining matrices are calculated using symmetry and
    !     recurrence relations by means of the subroutine dlmn.

    if (DABS(sinbet) < 1.d-14) then
        tgbet2 = zero
        sinbet = zero
        cosbet = DSIGN(one, cosbet)
    else if (DABS(sinbet) < 1.d-10 ) then
        tgbet2 = zero
        print*, 'WARNING in ROTAR: sinbet = ', sinbet, ' takes  0'
        sinbet = zero
        cosbet = DSIGN(one,cosbet)
    else
        tgbet2 = (one - cosbet) / sinbet
    end if
        do 10 l = 2, lmax
            l1 = l
            call dlmn(l1, ltot, sinal, cosal, cosbet, tgbet2, singa, cosga, dl, rl)
            10 continue
            201 continue
            !      write(6,902) lmax
            ! 902 format(2x,'the rotation matrices for lmax up to ',i2,' are:')
            !      call m2cmatprt(lmax,ltot,rl)
    return
end subroutine ROTAR


!***********************************************************************
!                                                                      *
!   subroutine dlmn                                                    *
!                                                                      *
!   this subroutine generates the matrices dl(m',m;l) for a fixed value*
!   of the orbital quantum number l, and it needs the dl(l-2;m',m) and *
!   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence*
!   relations. the matrices dl(m',m;l) are the rotation matrices for   *
!   complex spherical harmonics.                                       *
!                                                                      *
!***********************************************************************
subroutine DLMN(l, ltot, sinal, cosal, cosbet, tgbet2, singa, cosga, dl, rl)
    implicit real*8 (a-h,o-z)

    integer , parameter :: mxn = 15, mxl = 5
    integer :: ltot, l, iinf, isup, m, mp, laux, lbux, lauz, lbuz
    real*8  :: root2, one, pt5, sinal, cosal, cosbet, tgbet2, singa, cosga, dl, rl
    real*8  :: root, rooti, al, al1, tal1, ali, cosaux, amp, aux, cux, am, auz, factor, term, cuz
    real*8  :: sign, cosmal, sinmal, cosmga, sinmga, d1, d2, cosag, cosagm, sinag, sinagm
    dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
    dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
    common /roots/ root(0:2*mxl), rooti(0:2*mxl)
    data root2/1.414213562373095d0/, one/1.0d0/, pt5/0.5d0/
    iinf=1-l
    isup=-iinf

    !     computation of the dl(m',m;l) matrix, mp is m' and m is m.
    !     first row by recurrence: see equations 19 and 20 of reference (6)
    dl(l,l,l) = dl(isup,isup,l - 1) * (one + cosbet) * pt5
    dl(l,-l,l) = dl(isup,-isup,l - 1) * (one - cosbet) * pt5

    do 10 m = isup, iinf, -1
        dl(l,m,l) = -tgbet2 * root(l + m + 1) * rooti(l - m) * dl(l,m + 1,l)
        10 continue

        !     the rows of the upper quarter triangle of the dl(m',m;l) matrix
        !     see equation 21 of reference (6)
        al = l
        al1 = al - one
        tal1 = al + al1
        ali =one / al1
        cosaux = cosbet * al * al1
        do 11 mp = l - 1, 0, -1
            amp = mp
            laux = l + mp
            lbux = l - mp
            aux = rooti(laux) * rooti(lbux) * ali
            cux = root(laux - 1) * root(lbux - 1) * al

            do 12 m = isup, iinf, -1
                am = m
                lauz = l + m
                lbuz = l - m
                auz = rooti(lauz) * rooti(lbuz)
                factor = aux * auz
                term = tal1 * (cosaux - am * amp) * dl(mp,m,l - 1)

                if(lbuz /= 1 .AND. lbux /= 1) then
                    cuz = root(lauz - 1) * root(lbuz - 1)
                    term = term - dl(mp,m,l - 2) * cux * cuz
                end if
                dl(mp,m,l) = factor * term
                12 continue
                iinf = iinf + 1
                isup = isup - 1
                11 continue
                !     the remaining elements of the dl(m',m;l) matrix are calculated
                !     using the corresponding symmetry relations:
                !     reflexion ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
                !     inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
                !
                !     reflexion
                sign = one
                iinf = -l
                isup = l - 1
                do 13 m = l, 1, -1

                    do 14 mp = iinf, isup
                        dl(mp,m,l) = sign * dl(m,mp,l)
                        sign = -sign
                        14 continue
                        iinf = iinf + 1
                        isup = isup - 1
                        13 continue

                        !     inversion
                        iinf = -l
                        isup = iinf
                        do 15 m = l - 1, -l, -1
                            sign = -one

                            do 16 mp = isup, iinf, -1
                                dl(mp,m,l) = sign * dl(-mp,-m,l)
                                sign = -sign
                                16 continue
                                isup = isup + 1
                                15 continue

                                !     computation of the rotation matrices rl(m',m;l) for real spherical
                                !     harmonics using the matrices dl(m',m;l) for complex spherical
                                !     harmonics: see equations 10 to 18 of reference (6)
                                rl(0,0,l) = dl(0,0,l)
                                cosmal = cosal
                                sinmal = sinal
                                sign = - one

                                do 17 mp = 1, l
                                    cosmga = cosga
                                    sinmga = singa
                                    aux = root2 * dl(0,mp,l)
                                    rl(mp,0,l) = aux * cosmal
                                    rl(-mp,0,l)= aux * sinmal

                                    do 18 m = 1, l
                                        aux = root2 * dl(m,0,l)
                                        rl(0,m,l) = aux * cosmga
                                        rl(0,-m,l) = -aux * sinmga
                                        d1 = dl(-mp,-m,l)
                                        d2 = sign * dl(mp,-m,l)
                                        cosag = cosmal * cosmga - sinmal * sinmga
                                        cosagm = cosmal * cosmga + sinmal * sinmga
                                        sinag = sinmal * cosmga + cosmal * sinmga
                                        sinagm = sinmal * cosmga - cosmal * sinmga
                                        rl(mp,m,l)  = d1 * cosag + d2 * cosagm
                                        rl(mp,-m,l) = -d1 * sinag + d2 * sinagm
                                        rl(-mp,m,l) = d1 * sinag + d2 * sinagm
                                        rl(-mp,-m,l) = d1 * cosag - d2 * cosagm
                                        aux = cosmga * cosga - sinmga * singa
                                        sinmga = sinmga * cosga + cosmga * singa
                                        cosmga = aux
                                        18 continue
                                        sign = -sign
                                        aux= cosmal * cosal - sinmal * sinal
                                        sinmal = sinmal * cosal + cosmal * sinal
                                        cosmal = aux
                                        17 continue
                                        return
end subroutine DLMN


subroutine UTIL_OVERLAP
    implicit real*8 (a-h,o-z)

    real*8  :: re, reali, fact, facti, root, rooti, rli
    integer :: i

    integer, parameter :: mxn = 15
    integer, parameter :: mxl = 5

    common /const/ re(0:1000), reali(0:1000), fact(0:150), facti(0:150)
    common /roots/ root(0:2*mxl), rooti(0:2*mxl)

    fact(0) = 1.d0
    facti(0) = 1.d0

    do i = 1, 150
        fact(i) = fact(i-1) * i
        facti(i) = 1.d0 / fact(i)
    end do

    re(0) = 0.d0
    reali(0) = 0.d0 ! Just a trick to avoid problems in some recursion

    do i = 1, 1000
        re(i) = re(i-1) + 1.d0
        reali(i) = 1.d0 / re(i)
    end do

    ! These arrrays are REAL*8
    root(0) = 0.d0
    rooti(0) = 1.d200     ! Just a trick. It shouldn't be used

    do i = 1, 2*mxl
        rli = re(i)
        root(i) = dsqrt(rli)
        rooti(i) = 1.d0 / root(i)
    end do
end subroutine UTIL_OVERLAP

end module Overlap_Builder
