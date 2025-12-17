module Dissociative
    use types_m
    use constants_m
    use Read_parms
    use EDIT_routines, only : ReGroup, Bring_into_PBCBox
    use EDT_util_m   , only : parse_this
    
    implicit none    
    private

    public :: DWFF
 
    ! Derived types
    type statistics
       real*8 :: mean
       real*8 :: sdv
    end type statistics 

    ! Arrays describing connectivity and species
    integer, allocatable :: O_ptr(:)             ! O-atom indices
    integer, allocatable :: H_ptr(:)             ! H-atom indices
    integer, allocatable :: OH_bond_order(:)     ! bond order per-water
    integer, allocatable :: OO_neighbors(:)      ! count of O neighbors
    integer, allocatable :: HOH_indices(:,:)     ! per-water H indices (2)
    integer, allocatable :: Zundel_counter(:,:)  ! life-time counter
    integer, allocatable :: dimer_counter (:,:)  ! life-time counter
    integer, allocatable :: dimer_list(:)        ! indices of atoms in dimers

    real(8), allocatable :: table_OH(:,:)        ! OH distance table
    real(8), allocatable :: table_OO(:,:)        ! OO distance table

    logical, allocatable :: HOH_modified(:)
    type(universe), allocatable :: new_trj(:)

    ! module variables
    integer :: n_frames, n_atoms, n_mols, nHX, nOX, unit1, unit2, unit3

    ! module parameters 
    real*8, parameter :: OH_covalent_len = 1.25d0
    real*8, parameter :: OO_dimer_len    = 2.90d0
    real*8, parameter :: OH_distance_ref = 0.96d0
    real*8, parameter :: HH_distance_ref = 1.52d0
    real*8, parameter :: HOH_angle_ref   = 104.4d0

contains
!
!
!
!======================
subroutine DWFF(trj)
!======================
    implicit none
    type(universe), intent(in) :: trj(:)

    ! Local variables
    integer :: n
    type(statistics) :: OH_bond_length, HOH_angle

    CALL system( dynemoldir//"env.sh manipulate" )

    !------------------------------------
    ! Open output files
    !------------------------------------
    call open_output_files()

    !------------------------------------
    ! Trajectory size and initialization
    !------------------------------------
    n_frames = size(trj)

    call preprocess(trj(1)%atom, n_frames)

    ! Copy trajectory for optional visualization
    new_trj = trj

    !------------------------------------
    ! Main analysis loop
    !------------------------------------
    do n = 1, n_frames

        call bond_topology(n, trj(n)%atom, trj(n)%box)

        OH_bond_length = get_OH_covalent_bond_length( &
                              trj(n)%atom, trj(n)%box, verbose=.false. )

        HOH_angle = get_HOH_angles( &
                        trj(n)%atom, trj(n)%box, verbose=.false. )

        call data_output(n, n_frames, OH_bond_length, HOH_angle)

    end do

    !------------------------------------
    ! Optional visualization
    !------------------------------------
    call ask_user_visualization(trj)

    stop
end subroutine DWFF
!
!
!
!============================================
subroutine bond_topology(frame, atom, Txyz)
!============================================
    implicit none
    integer     , intent(in) :: frame
    type(atomic), intent(in) :: atom(:)
    real(8)     , intent(in) :: Txyz(3)

    ! Local variables
    integer :: i, j
    real(8) :: rij(3), rij2, rijlen

    ! --------------------------------------------------------
    ! Compute O–H and O–O distances (minimum image convention)
    ! --------------------------------------------------------
    table_OH = 0.0d0
    table_OO = 0.0d0

    associate (OX => atom(O_ptr), HX => atom(H_ptr))

        ! O–H distances
        do i = 1, nOX
            do j = 1, nHX
                rij = OX(i)%xyz - HX(j)%xyz
                rij = rij - Txyz * dnint(rij / Txyz)

                rij2   = dot_product(rij, rij)
                rijlen = sqrt(rij2)

                table_OH(i, j) = rijlen
            end do
        end do

        ! O–O distances
        do i = 1, nOX
            do j = i + 1, nOX
                rij = OX(i)%xyz - OX(j)%xyz
                rij = rij - Txyz * dnint(rij / Txyz)

                rij2   = dot_product(rij, rij)
                rijlen = sqrt(rij2)

                table_OO(i, j) = rijlen
                table_OO(j, i) = rijlen
            end do
            table_OO(i, i) = huge(1.0d0)
        end do

    end associate

    ! --------------------------------------------------------
    ! Count neighbors for each oxygen
    ! --------------------------------------------------------
    OH_bond_order = 0
    OO_neighbors  = 0

    do concurrent (i = 1:nOX)
        OH_bond_order(i) = count(table_OH(i, :) <  OH_covalent_len)
        OO_neighbors(i)  = count(table_OO(i, :) <= OO_dimer_len)
    end do

    where (OH_bond_order /= 2)
        HOH_modified = .true.
    end where

    ! --------------------------------------------------------
    ! Bond-order sanity check
    ! --------------------------------------------------------
    do i = 1, nOX
        select case (OH_bond_order(i))
        case (0)
            write(*,*) "Error: Oxygen", i, "has no covalent hydrogens."
            stop
        case (4:)
            write(*,*) "Error: Oxygen", i, "has bond order > 3."
            stop
        end select
    end do

    ! --------------------------------------------------------
    ! Build HOH triplets (H–O–H)
    ! --------------------------------------------------------
    HOH_indices = 0
    do i = 1, nOX
        j = 2*(i - 1)
        HOH_indices(i, 1) = H_ptr(j + 1)
        HOH_indices(i, 2) = O_ptr(i)
        HOH_indices(i, 3) = H_ptr(j + 2)
    end do

    call identify_species(frame, atom)

end subroutine bond_topology
!
!
!
!=======================================================================
function get_OH_covalent_bond_length(atom, Txyz, verbose) result(this)
!=======================================================================
    implicit none
    type(atomic), intent(in) :: atom(:)
    real(8)     , intent(in) :: Txyz(3)
    logical     , intent(in) :: verbose

    ! Local variables
    real(8) :: rij(3), rik(3)
    real(8) :: rijlen, riklen
    real(8) :: OH_bond_mean, OH_bond_sdv
    integer :: ati, atj, atk, n, n_bond
    real(8), allocatable :: OH_bond_length(:)
    type(statistics) :: this

    !================================
    ! OH covalent bond length
    !          J     K
    !           \   /
    !            \ /
    !             I
    !================================

    allocate(OH_bond_length(nOX), source=0.0d0)

    do n = 1, nOX

        if (HOH_modified(n)) cycle

        atj = HOH_indices(n, 1)   ! first hydrogen
        ati = HOH_indices(n, 2)   ! oxygen
        atk = HOH_indices(n, 3)   ! second hydrogen

        ! rij = r_j - r_i (minimum image)
        rij = atom(atj)%xyz - atom(ati)%xyz
        rij = rij - Txyz * dnint(rij / Txyz)
        rijlen = sqrt(dot_product(rij, rij))

        ! rik = r_k - r_i (minimum image)
        rik = atom(atk)%xyz - atom(ati)%xyz
        rik = rik - Txyz * dnint(rik / Txyz)
        riklen = sqrt(dot_product(rik, rik))

        OH_bond_length(n) = 0.5d0 * (rijlen + riklen)

    end do

    !-----------------------------------------
    ! Statistics: ignore dissociated molecules
    !-----------------------------------------
    n_bond = count(.not. HOH_modified)

    if (n_bond > 0) then
        OH_bond_mean = sum(OH_bond_length, mask=.not. HOH_modified) / n_bond
        OH_bond_sdv  = sqrt( sum( (OH_bond_length - OH_bond_mean)**2, &
                                  mask=.not. HOH_modified ) / n_bond )
    else
        OH_bond_mean = 0.0d0
        OH_bond_sdv  = 0.0d0
    end if

    this%mean = OH_bond_mean
    this%sdv  = OH_bond_sdv

    if (verbose) then
        write(*,*) OH_bond_mean, OH_bond_sdv
        do n = 1, nOX
            write(26,*) n, OH_bond_length(n)
        end do
    end if

    deallocate(OH_bond_length)

end function get_OH_covalent_bond_length
!
!
!
!==========================================================
function get_HOH_angles(atom, Txyz, verbose) result(this)
!==========================================================
    implicit none
    type(atomic), intent(in) :: atom(:)
    real(8)     , intent(in) :: Txyz(3)
    logical     , intent(in) :: verbose

    ! Local variables
    real(8) :: rij(3), rik(3)
    real(8) :: rijlen, riklen
    real(8) :: cos_theta, theta
    real(8) :: HOH_ang_mean, HOH_ang_sdv
    integer :: ati, atj, atk, n, n_ang
    real(8), allocatable :: HOH_angles(:)
    type(statistics) :: this

    !================================
    ! HOH angle (bending)
    !          J     K
    !           \   /
    !            \ /
    !             I
    !================================

    allocate(HOH_angles(nOX), source=0.0d0)

    do n = 1, nOX

        if (HOH_modified(n)) cycle

        atj = HOH_indices(n, 1)   ! first hydrogen
        ati = HOH_indices(n, 2)   ! oxygen
        atk = HOH_indices(n, 3)   ! second hydrogen

        ! rij = r_j - r_i (minimum image)
        rij = atom(atj)%xyz - atom(ati)%xyz
        rij = rij - Txyz * dnint(rij / Txyz)
        rijlen = sqrt(dot_product(rij, rij))

        ! rik = r_k - r_i (minimum image)
        rik = atom(atk)%xyz - atom(ati)%xyz
        rik = rik - Txyz * dnint(rik / Txyz)
        riklen = sqrt(dot_product(rik, rik))

        !-----------------------------------------
        ! HOH angle at oxygen: angle(j-i-k)
        !-----------------------------------------
        cos_theta = dot_product(rij, rik) / (rijlen * riklen)

        ! Numerical safety
        cos_theta = max(-1.0d0, min(1.0d0, cos_theta))

        theta = acos(cos_theta)
        HOH_angles(n) = theta * rad_2_deg

    end do

    !-----------------------------------------
    ! Statistics: ignore dissociated molecules
    !-----------------------------------------
    n_ang = count(.not. HOH_modified)

    if (n_ang > 0) then
        HOH_ang_mean = sum(HOH_angles, mask=.not. HOH_modified) / n_ang
        HOH_ang_sdv  = sqrt( sum( (HOH_angles - HOH_ang_mean)**2, &
                                  mask=.not. HOH_modified ) / n_ang )
    else
        HOH_ang_mean = 0.0d0
        HOH_ang_sdv  = 0.0d0
    end if

    this%mean = HOH_ang_mean
    this%sdv  = HOH_ang_sdv

    if (verbose) then
        write(*,*) HOH_ang_mean, HOH_ang_sdv
        do n = 1, nOX
            write(25,*) n, HOH_angles(n)
        end do
    end if

    deallocate(HOH_angles)

end function get_HOH_angles
!
!
!
!=================================================================
subroutine data_output(frame, n_frames, OH_bond_length, HOH_angle)
!=================================================================
    implicit none
    integer,          intent(in) :: frame, n_frames
    type(statistics), intent(in) :: OH_bond_length, HOH_angle

    ! Local variables
    integer :: i, un
    integer :: idx(2)
    logical :: done
    integer, allocatable :: lifetime(:)
    type(statistics) :: Zundel_lifetime, dimer_lifetime

    !------------------------------------------------------
    ! Write instantaneous averages
    !------------------------------------------------------
    if (frame == 1) then
        write(unit1, '(a)') '# frame  <OH>  <HOH>'
    end if

    write(unit1, '(i8,2f16.8)') frame, OH_bond_length%mean, HOH_angle%mean

    !------------------------------------------------------
    ! Final processing at last frame
    !------------------------------------------------------
    if (frame /= n_frames) return

    !======================================================
    ! Zundel lifetimes
    !======================================================
    allocate(lifetime(0))
    done = .false.

    do while (.not. done)
        idx = maxloc(Zundel_counter)
        if (Zundel_counter(idx(1), idx(2)) == 0) then
            done = .true.
        else
            lifetime = [lifetime, Zundel_counter(idx(1), idx(2))]
            Zundel_counter(idx(1), idx(2)) = 0
        end if
    end do

    if (size(lifetime) > 0) then
        Zundel_lifetime%mean = sum(lifetime) / real(size(lifetime), 8)
        Zundel_lifetime%sdv  = sqrt( sum((lifetime - Zundel_lifetime%mean)**2) &
                                    / real(size(lifetime), 8) )
    else
        Zundel_lifetime%mean = 0.0d0
        Zundel_lifetime%sdv  = 0.0d0
    end if

    open(newunit=un, file='DWFF.trunk/Zundel_lifetime.dat', status='unknown', action='write')
        write(un,'(a,2f16.8)') '# mean  sdv:', Zundel_lifetime%mean, Zundel_lifetime%sdv
        do i = 1, size(lifetime)
            write(un,'(i8)') lifetime(i)
        end do
    close(un)

    deallocate(lifetime)

    !======================================================
    ! Dimer lifetimes
    !======================================================
    allocate(lifetime(0))
    done = .false.

    do while (.not. done)
        idx = maxloc(dimer_counter)
        if (dimer_counter(idx(1), idx(2)) == 0) then
            done = .true.
        else
            lifetime = [lifetime, dimer_counter(idx(1), idx(2))]
            dimer_counter(idx(1), idx(2)) = 0
        end if
    end do

    if (size(lifetime) > 0) then
        dimer_lifetime%mean = sum(lifetime) / real(size(lifetime), 8)
        dimer_lifetime%sdv  = sqrt( sum((lifetime - dimer_lifetime%mean)**2) &
                                   / real(size(lifetime), 8) )
    else
        dimer_lifetime%mean = 0.0d0
        dimer_lifetime%sdv  = 0.0d0
    end if

    open(newunit=un, file='DWFF.trunk/dimer_lifetime.dat', status='unknown', action='write')
        write(un,'(a,2f16.8)') '# mean  sdv:', dimer_lifetime%mean, dimer_lifetime%sdv
        do i = 1, size(lifetime)
            write(un,'(i8)') lifetime(i)
        end do
    close(un)

    deallocate(lifetime)

end subroutine data_output
!
!
!
!=============================================
 subroutine ask_user_visualization(trj)
    implicit none
    type(universe), intent(in) :: trj(:)
!=============================================

    ! local variables
    character(1)      :: YorN
    character(len=80) :: line
    integer, allocatable :: O_list(:)

    write(*,'(/a)') 'Visualize a particular structure (y/n)?'
    write(*,'(a)', advance='no') '>>> '
    read(*,'(a)') YorN

    if (YorN /= 'y') return

    write(*,'(/a)') 'Enter the Oxygen atom indices for visualization:'
    write(*,'(a)') 'List them separated by spaces:'
    read(*,'(a)') line

    O_list = parse_this(line)
    call Save_PDB_Trajectory(trj, "frames-DWFF.pdb", O_list)

end subroutine ask_user_visualization
!
!
!
!=============================================
 subroutine open_output_files()
!=============================================
    open(newunit=unit1, file="DWFF.trunk/DWFF_data", status="unknown", action="write")
    open(newunit=unit2, file="DWFF.trunk/dimer_list",status="unknown", action="write")
    open(newunit=unit3, file="DWFF.trunk/ions_list", status="unknown", action="write")
end subroutine open_output_files
!
!
!
!======================================================
 subroutine Save_PDB_Trajectory(trj, file_name, O_list)
!======================================================
    implicit none
    type(universe)        , intent(in) :: trj(:)
    character(*), optional, intent(in) :: file_name
    integer               , intent(in) :: O_list(:)
    
    ! local variables ...
    integer      :: i, j, k, frame_step
    character(1) :: YorN 
    integer      :: un               ! output file unit
    integer, allocatable :: res_number(:)
    
    ! Ask user for a frame step
    write(*,'(/a)',advance='no') ' Saving with Frame step : '
    read (*,'(i3)'             ) frame_step
    
    ! Open output file
    If( present(file_name) ) then
        OPEN(newunit=un, file='DWFF.trunk/'//file_name, status='unknown', action='write')
    else
        OPEN(newunit=un, file='DWFF.trunk/frames-output.pdb', status='unknown', action='write')
    end if
    
    !-------------------------------------------
    ! Prepare: mark residues belonging to O_list
    !-------------------------------------------
    allocate(res_number(size(O_list)))
    
    do i = 1 , size(trj)
       associate( atom => new_trj(i)% atom )
           res_number = atom(O_list)% nresid
       
           do j = 1, size(res_number)
              where( atom% nresid == res_number(j) )  atom% resid = "XXX"
           end do
       
           atom% group = .true.
           call ReGroup( new_trj(i) )                                                                                                                
           call Bring_into_PBCBox(new_trj(i), res_name = "XXX")
       end associate
    end do
    
    deallocate(res_number)
    
    !-------------------------------------------
    !             Write frames 
    !-------------------------------------------
    do j = 1 , size(trj) , frame_step
    
        write(un,5) 'TITLE'  , 'manipulated by edview     t= ',trj(j)%time
        write(un,1) 'CRYST1' , trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
        write(un,3) 'MODEL' , j
    
        do i = 1 , size(trj(1)%atom)
            write(un,2) 'ATOM  '                      ,  &    ! <== non-standard atom
                 i                                    ,  &    ! <== global number
                 new_trj(j)%atom(i)%MMSymbol          ,  &    ! <== atom type
                 ' '                                  ,  &    ! <== alternate location indicator
                 new_trj(j)%atom(i)%resid             ,  &    ! <== residue name
                 ' '                                  ,  &    ! <== chain identifier
                 new_trj(j)%atom(i)%nresid            ,  &    ! <== residue sequence number
                 ' '                                  ,  &    ! <== code for insertion of residues
                 ( new_trj(j)%atom(i)%xyz(k) , k=1,3 ),  &    ! <== xyz coordinates 
                 1.00                                 ,  &    ! <== occupancy
                 0.00                                 ,  &    ! <== temperature factor
                 ' '                                  ,  &    ! <== segment identifier
                 ' '                                  ,  &    ! <== here only for tabulation purposes
                 new_trj(j)%atom(i)%symbol            ,  &    ! <== chemical element symbol
                 new_trj(j)%atom(i)%charge                    ! <== charge on the atom
        end do
        write(un,'(a)') 'MASTER'
        write(un,'(a)') 'END'
    
    end do
    
    close(un)
    
    1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
    2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)
    3 FORMAT(a5,i8)
    4 FORMAT(a6,t15,a21)
    5 FORMAT(a5,t1,a35,f12.7)
    6 FORMAT(a72)
    
    write(*,'(/a)') ' >>> frames-DWFF.pdb : writing done, press any key <<<'
    write(*,'(/a)') "That's all ? (y/n)"
    read (*,'(a)') YorN
    if( YorN /= "n" ) stop

end subroutine Save_PDB_Trajectory
!
!
!
!=========================================
subroutine identify_species(frame, atom)
!=========================================
    implicit none
    integer     , intent(in) :: frame
    type(atomic), intent(in) :: atom(:)

    ! Local variables
    integer :: i, j
    integer :: O_acceptor, O_donor, O_hydrox
    integer :: total_charge, hydrox_charge, dimer_charge
    real*8  :: d_OO, d_hydrox
    logical :: more_than_one

    more_than_one = .false.

    !--------------------------------------------------
    ! Identify dimers (updates dimer_counter internally)
    !--------------------------------------------------
    call identify_dimer(frame, atom)

    !--------------------------------------------------
    ! Loop over oxygen sites
    !--------------------------------------------------
    do i = 1, nOX

        ! Only oxygens with excess protons
        if (OH_bond_order(i) <= 2) cycle

        O_acceptor = i

        ! Nearest oxygen to acceptor
        O_donor = minloc(table_OO(:, i), dim=1)
        d_OO    = table_OO(O_donor, O_acceptor)

        ! Mask this pair to avoid reuse
        table_OO(O_acceptor, O_donor) = huge(1.0d0)
        table_OO(O_donor, O_acceptor) = huge(1.0d0)

        total_charge = OH_bond_order(O_acceptor) + OH_bond_order(O_donor) - 4

        select case (total_charge)

        !==================================================
        ! Candidate for hydronium-related species
        !==================================================
        case (1)

            if (d_OO < OO_dimer_len) then

                ! Look for hydroxyl oxygen bound to donor
                O_hydrox = minloc(table_OO(:, O_donor), dim=1)
                d_hydrox = table_OO(O_donor, O_hydrox)

                table_OO(O_donor, O_hydrox) = huge(1.0d0)
                table_OO(O_hydrox, O_donor) = huge(1.0d0)

                hydrox_charge = OH_bond_order(O_hydrox) - 2

                !-------------------------------
                ! Zundel-like configuration
                !-------------------------------
                if (hydrox_charge < 0) then

                    Zundel_counter(O_hydrox, O_acceptor) = &
                        Zundel_counter(O_hydrox, O_acceptor) + 1

                    if (more_than_one) write(unit3,*) repeat(".", 48)

                    write(unit3,*) frame, "Zundel", &
                                   O_ptr(O_donor), O_ptr(O_acceptor), d_OO, total_charge
                    write(unit3,*) frame, "Hydrox", &
                                   O_ptr(O_hydrox), O_ptr(O_acceptor), d_hydrox, hydrox_charge

                    more_than_one = .true.

                !-------------------------------
                ! Neutral dimer
                !-------------------------------
                elseif (hydrox_charge == 0) then

                    O_hydrox = minloc(table_OO(:, O_acceptor), dim=1)
                    d_hydrox = table_OO(O_acceptor, O_hydrox)

                    table_OO(O_acceptor, O_hydrox) = huge(1.0d0)
                    table_OO(O_hydrox, O_acceptor) = huge(1.0d0)

                    if (OH_bond_order(O_acceptor) + OH_bond_order(O_hydrox) == 4) then
                        dimer_charge = 0
                        write(unit2,2) frame, O_ptr(O_hydrox), O_ptr(O_acceptor), d_hydrox
                        dimer_counter(O_donor, O_acceptor) = &
                            dimer_counter(O_donor, O_acceptor) + 1
                    end if

                end if

            else
                ! Isolated hydronium
                write(unit3,*) frame, "H3O+", O_ptr(O_donor), O_ptr(O_acceptor), d_OO
            end if

        !==================================================
        ! Candidate for hydroxide
        !==================================================
        case (-1)
            write(unit3,*) frame, "HO-", O_ptr(O_donor), O_ptr(O_acceptor), d_OO

        end select

    end do

    write(unit3,*) repeat("=", 48)

2 format(4x, i8, 6x, i8, 5x, i8, 7x, f8.4)

end subroutine identify_species
!
!
!
!=========================================
subroutine identify_dimer(frame, atom)
!=========================================
    implicit none
    integer     , intent(in) :: frame
    type(atomic), intent(in) :: atom(:)

    ! Local variables
    integer :: i, j
    integer :: O_acceptor, O_donor
    integer :: total_charge
    integer :: HOH(3)
    real(8) :: d_OO

    if (frame /= 1) write(unit2,*) repeat("=", 60)
    write(unit2,*) "........frame......O_donor.....O_acceptor.......d_OO"

    !--------------------------------------------------
    ! Loop over oxygen sites
    !--------------------------------------------------
    do i = 1, nOX

        ! Only oxygens with excess protons
        if (OH_bond_order(i) <= 2) cycle

        O_acceptor = i

        ! Nearest oxygen donor
        O_donor = minloc(table_OO(:, i), dim=1)
        d_OO    = table_OO(O_donor, O_acceptor)

        total_charge = OH_bond_order(O_acceptor) + &
                       OH_bond_order(O_donor)    - 4

        !--------------------------------------------------
        ! Neutral water dimer: total charge = 0
        !--------------------------------------------------
        if (total_charge == 0 .and. d_OO < OO_dimer_len) then

            ! -- First H2O: acceptor
            j   = 2 * (O_acceptor - 1)
            HOH = [ H_ptr(j+1), O_ptr(O_acceptor), H_ptr(j+2) ]

            ! -- Second H2O: donor
            j   = 2 * (O_donor - 1)
            HOH = [ H_ptr(j+1), O_ptr(O_donor), H_ptr(j+2) ]

            write(unit2, 2) frame, O_ptr(O_donor), O_ptr(O_acceptor), d_OO

            dimer_counter(O_donor, O_acceptor) = &
                dimer_counter(O_donor, O_acceptor) + 1

        end if

    end do

2 format(4x, i8, 6x, i8, 5x, i8, 7x, f8.4)

end subroutine identify_dimer
!
!
!
!======================================
subroutine preprocess(atom, n_frames)
!======================================
    implicit none
    type(atomic), intent(in) :: atom(:)
    integer     , intent(in) :: n_frames

    ! Local variables
    integer :: i, j, k

    !------------------------------------
    ! Basic system sizes
    !------------------------------------
    n_atoms = size(atom)
    nOX     = count(atom%symbol == "O")
    nHX     = count(atom%symbol == "H")

    !------------------------------------
    ! Build pointer lists for O and H
    !------------------------------------
    allocate(O_ptr(nOX))
    allocate(H_ptr(nHX))

    i = 0
    j = 0
    do k = 1, n_atoms
        select case (atom(k)%symbol)
        case ("O")
            i = i + 1
            O_ptr(i) = k
        case ("H")
            j = j + 1
            H_ptr(j) = k
        end select
    end do

    !------------------------------------
    ! Allocate module-level work arrays
    !------------------------------------
    allocate(table_OH(nOX, nHX), source=0.0d0)
    allocate(table_OO(nOX, nOX), source=0.0d0)

    allocate(OH_bond_order(nOX), source=0)
    allocate(OO_neighbors(nOX), source=0)

    allocate(HOH_indices(nOX, 3), source=0)
    allocate(HOH_modified(nOX), source=.false.)

    allocate(Zundel_counter(nOX, nOX), source=0)
    allocate(dimer_counter (nOX, nOX), source=0)

    !------------------------------------
    ! Trajectory copy for visualization
    !------------------------------------
    allocate(new_trj(n_frames))

end subroutine preprocess
!
!
!
end module Dissociative
