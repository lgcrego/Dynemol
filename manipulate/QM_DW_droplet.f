module test_droplet

    use constants_m
    use types_m
    use ansi_colors
    use color_funcs
    use util_m       , only: TO_UPPER_CASE
    
    implicit none    
    private

    public :: DWFF_QM_droplet
 
    ! Arrays describing conectivity and species
    integer, allocatable :: O_ptr(:)             ! O-atom indices
    integer, allocatable :: H_ptr(:)             ! H-atom indices
    integer, allocatable :: OH_bond_order(:)     ! bond order per-water

    real(8), allocatable :: OH_distance_table(:,:)
    real(8), allocatable :: OO_distance_table(:,:)

    type(universe)       :: work_sys

    ! module variables
    integer :: n_atoms, nHX, nOX, unit3
    integer :: n_HOH_resids, n_WAT_resids

contains
!
!
!
!==============================
subroutine DWFF_QM_droplet(sys)
!==============================
    implicit none
    type(universe), intent(in):: sys

    ! local variables
    character(len=1) :: YorN

    call system( dynemoldir//"env.sh manipulate" ) 

    call preprocess(sys%atom)

    ! Keep an independent copy for topology and visualization operations
    work_sys = sys

    call bond_topology(sys%atom, sys%box)

    call save_pdb_file()

    call deallocate_and_leave

    write(*,'(/a)') bold//orange//'>>> ./DWFF.trunk/seed-DWFF.pdb : writing done <<<'//reset
    write(*,'(/a)') yellow_("That's all? (y/n)")
    read (*,'(a)') YorN
    if( YorN /= "n" ) stop

end subroutine DWFF_QM_droplet
!
!
!
!===================================
subroutine post_processing_analysis()
!===================================
    implicit none

    !local variables
    integer :: n_QM_HOH, zero_H, one_H, two_H, three_H, total_MM_charge

    print*, ""
    print*, cyan_(">>> diagnosis <<< ")
    print*, ""

    n_QM_HOH =  count( work_sys%atom(O_ptr)%fragment == "Q")
    print*, green_("number of HOH molecules in the droplet = "), n_QM_HOH
    print*, ""

    zero_H  = count(OH_bond_order(:) == 0 .and. work_sys%atom(O_ptr(:))%fragment == "Q")
    one_H   = count(OH_bond_order(:) == 1 .and. work_sys%atom(O_ptr(:))%fragment == "Q")
    two_H   = count(OH_bond_order(:) == 2 .and. work_sys%atom(O_ptr(:))%fragment == "Q")
    three_H = count(OH_bond_order(:) == 3 .and. work_sys%atom(O_ptr(:))%fragment == "Q")

    print*, green_("number of   O^2-   ions in the droplet = "), zero_H 
    print*, green_("number of   OH-    ions in the droplet = "), one_H  
    print*, green_("number of   H2O    mols in the droplet = "), two_H  
    print*, green_("number of   H3O+   ions in the droplet = "), three_H
    print*, ""

    total_MM_charge = -2*zero_H + -1*one_H + 1*three_H
    print*, green_("total MM charge of the HOH solvent = "), total_MM_charge
    print*, green_("total HOH solvent residues = "), n_HOH_resids
    print*, green_("total WAT solvent residues = "), n_WAT_resids
    print*, green_("total solvent residues     = "), n_HOH_resids + n_WAT_resids

end subroutine post_processing_analysis
!
!
!
!===================================
subroutine bond_topology(atom, Txyz)
!===================================
    implicit none
    type(atomic), intent(in) :: atom(:)
    real(8)     , intent(in) :: Txyz(3)

    ! Local variables
    integer :: i, j
    real(8) :: rij(3), rijlen

    ! --------------------------------------------------------
    ! Compute O–H and O–O distances (minimum image convention)
    ! --------------------------------------------------------
    OH_distance_table = 0.0d0
    OO_distance_table = 0.0d0

    associate (OX => atom(O_ptr), HX => atom(H_ptr))
        ! O–H distances
        do i = 1, nOX
            do j = 1, nHX
                rij = OX(i)%xyz - HX(j)%xyz
                rij = rij - Txyz * dnint(rij / Txyz)

                rijlen = norm2(rij)
                  
                ! not a symmetric matrix
                OH_distance_table(i, j) = rijlen
            end do
        end do

        ! O–O distances
        do i = 1, nOX
            do j = i + 1, nOX
                rij = OX(i)%xyz - OX(j)%xyz
                rij = rij - Txyz * dnint(rij / Txyz)

                rijlen = norm2(rij)

                OO_distance_table(i, j) = rijlen
                OO_distance_table(j, i) = rijlen
            end do
            OO_distance_table(i, i) = huge(1.0d0)
        end do
    end associate

    ! --------------------------------------------------------
    ! Build adjacency matrix
    ! given HX, find the closest OX
    ! --------------------------------------------------------
    allocate( work_sys%OH_pair(nHX,2) )

    do j = 1, nHx
       i = minloc( OH_distance_table(:, j), dim=1 )
       work_sys% OH_pair(j, 1) = O_ptr(i)
       work_sys% OH_pair(j, 2) = H_ptr(j)
    end do

    ! --------------------------------------------------------
    ! Count neighbors for each oxygen
    ! --------------------------------------------------------
    OH_bond_order = 0
    do i = 1, nOX
        OH_bond_order(i) = count( work_sys% OH_pair(:, 1) == O_ptr(i) )
    end do

    ! --------------------------------------------------------
    ! Bond-order sanity check
    ! --------------------------------------------------------
    do i = 1, nOX
        select case (OH_bond_order(i))
        case (0)
            write(*,*) "Error: Oxygen", i, "has no covalent hydrogens."
        case (4:)
            write(*,*) "Error: Oxygen", i, "has bond order > 3."
        end select
    end do

    call identify_species()

end subroutine bond_topology
!
!
!
!==========================
 subroutine save_pdb_file()
!==========================
    implicit none
    
    ! local variables ...
    integer :: i, j, k, n, O_idx
    integer :: out1, out2
    real(8) :: cutoff_radius
    character(3) :: res_name = "None"
    character(1) :: YorN, choice
    
    !-----------------------------------------------------------
    ! Assign each hydrogen's nr to its current oxygen owner's nr
    !-----------------------------------------------------------
    do n = 1, nHx                                                                                                                                       
       O_idx = work_sys% OH_pair(n, 1)
       work_sys%atom( H_ptr(n) )% nresid = work_sys%atom( O_idx )% nresid
    end do                                                                                                                                              

    !-------------------------------------------------                                                                                                     
    ! Regroup molecules using the updated conectivity                                                                                                     
    !-------------------------------------------------                                                                                                     
    write(*,'(/a)', advance='no') yellow_("radius of the Quantum droplet = ")
    read (*,*) cutoff_radius

    write(*,'(/a)', advance='no') yellow_("Place solute in the center of the box? (y/n): ")
    read (*,'(a)') YorN

    if( YorN == "y" ) &
    then
        write(*,'(/a)',advance='no') yellow_("Enter residue name of the solute: ")
        read (*,'(a)') res_name
        res_name = TO_UPPER_CASE(res_name)
        call translate_to_centroid(work_sys, res_name=res_name)
    else 
        call translate_to_centroid(work_sys, res_name="None")
    end if
    
    call ReGroup(work_sys)

    call QM_droplet( work_sys, cutoff_radius, res_name )

    write(*,'(/a/)') bold // orange_(">>> Save as: <<<") // reset
    write(*,'(a)') green // ' 1 :' // reset // ' keep only the QM atoms '
    write(*,'(a)') green // ' 2 :' // reset // ' embed the QM-HOH droplet in WAT box'
    read (*,'(a)') choice

    select case (choice)
        case("1")
            call eliminate_classical_atoms( work_sys )
            call pack_HOH_atoms(work_sys%atom)
            call post_processing_analysis
        case("2") 
            call pack_solvent_atoms(work_sys%atom)
            call post_processing_analysis
    end select 

    !-------------------------------------------
    !             Writing seed.pdb 
    !-------------------------------------------
    ! Open output files
    OPEN(newunit=out1, file='DWFF.trunk/seed-DWFF.pdb', status='replace', action='write')
    
    write(out1,5) 'TITLE' , 'manipulated by DynEMol    t= ', work_sys%time
    write(out1,1) 'CRYST1', (work_sys%box(j), j=1,3), 90.0, 90.0, 90.0, 'P 1', '1'
    
    do i = 1, size(work_sys%atom)
        write(out1,2) 'ATOM  '                 ,  &    ! <== non-standard atom
             i                                 ,  &    ! <== global number
             work_sys%atom(i)%MMSymbol         ,  &    ! <== atom type
             ' '                               ,  &    ! <== alternate location indicator
             work_sys%atom(i)%resid            ,  &    ! <== residue name
             work_sys%atom(i)%fragment         ,  &    ! <== fragment
             work_sys%atom(i)%nresid           ,  &    ! <== residue sequence number
             ' '                               ,  &    ! <== code for insertion of residues
             ( work_sys%atom(i)%xyz(k), k=1,3 ),  &    ! <== xyz coordinates 
             1.00                              ,  &    ! <== occupancy
             0.00                              ,  &    ! <== temperature factor
             ' '                               ,  &    ! <== segment identifier
             ' '                               ,  &    ! <== here only for tabulation purposes
             work_sys%atom(i)%symbol           ,  &    ! <== chemical element symbol
             work_sys%atom(i)%charge                   ! <== charge on the atom
    end do

    write(out1,'(a)') 'MASTER'
    write(out1,'(a)') 'END'
    close(out1)
    !--------------------------------------------------------

    !--------------------------------------------------------
    !             Writing velocities 
    !--------------------------------------------------------
    if( .not. any(work_sys%atom%vel(1) > low_prec) ) then
        ! do nothing
    else
        OPEN(newunit=out2, file='DWFF.trunk/seed-velocity_MM', status='unknown', action='write')
            do i = 1, size(work_sys%atom)
                write(out1,*) work_sys%atom(i)%vel
            end do
        close(out2)
    end if
    !--------------------------------------------------------

    1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
    2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)
    5 FORMAT(a5,t1,a35,f12.7)

end subroutine save_pdb_file
!
!
!
!============================
subroutine identify_species()
!============================
    implicit none

    ! Local variables
    integer :: i
    logical :: more_than_one

    more_than_one = .false.

    open(newunit=unit3, file="DWFF.trunk/charged_species_list", status="unknown", action="write")

    write(unit3,*) repeat("=", 70)

    !--------------------------
    ! Loop over oxygen sites
    !--------------------------
    do i = 1, nOX
        !----------------------
        ! charged species  
        !----------------------
        if (OH_bond_order(i) /= 2) &
        then
            if (more_than_one) write(unit3,*) repeat(".", 70)

            call save_charged_species(i)

            more_than_one = .true.
        end if
    end do

    write(unit3,*) repeat("=", 70)

    close(unit3) 

end subroutine identify_species
!
!
!
!==========================
subroutine preprocess(atom)
!==========================
    implicit none
    type(atomic), intent(in) :: atom(:)

    ! Local variables
    integer :: i, j, k

    !------------------------------------
    ! Basic system sizes
    !------------------------------------
    n_atoms = size(atom)
    nOX     = count(atom%MMsymbol == "OX")
    nHX     = count(atom%MMsymbol == "HX")

    !------------------------------------
    ! Build pointer lists for O and H
    !------------------------------------
    allocate(O_ptr(nOX))
    allocate(H_ptr(nHX))

    i = 0
    j = 0
    do k = 1, n_atoms
        select case (atom(k)%MMsymbol)
        case ("OX")
            i = i + 1
            O_ptr(i) = k
        case ("HX")
            j = j + 1
            H_ptr(j) = k
        end select
    end do

    !------------------------------------
    ! Allocate module-level work arrays
    !------------------------------------
    allocate(OH_distance_table(nOX, nHX), source=0.0d0)
    allocate(OO_distance_table(nOX, nOX), source=0.0d0)
    allocate(OH_bond_order(nOX)         , source=0)

end subroutine preprocess
!
!
!
!========================================
subroutine save_charged_species(i)
!========================================
    implicit none
    integer, intent(in) :: i

    ! Local variables
    integer :: total_charge

    select case(OH_bond_order(i))
         case(:1)
              total_charge = (-1)*OH_bond_order(i) 
              write(unit3,10) "OH-", "Oxygen = ",O_ptr(i), "charge = ", total_charge

         case(3:)
              total_charge = OH_bond_order(i) 
              write(unit3,10) "H3O+", "Oxygen = ",O_ptr(i), "charge = ", total_charge 
    end select

10  format(t3,a4,t33,a9,i4,t54,a9,i4)

end subroutine save_charged_species
!
!
!
!=========================
subroutine ReGroup(system)
!=========================
implicit none
type(universe) , intent(inout) :: system

!local variables
integer :: nr, i, k, ref
integer :: nr_min, nr_max 
real*8  :: dxyz(3), Txyz(3), centroid(3)
integer, allocatable :: in_range(:)

if (product(system%box) == 0) then
  Print*, "ERROR: simulation box has zero length in at least one dimension", system%box; stop
end if

Txyz = system%box

associate( atom => system%atom )

    nr_min = minval(atom%nresid)
    nr_max = maxval(atom%nresid)

    do nr = nr_min, nr_max

        in_range = pack( [(k, k=1,system%N_of_atoms)], atom%nresid == nr )

        ! sanity check
        if( size(in_range) == 0 ) cycle
    
        ! find the OX of this residue and use it as the unwrap reference
        ref = 0
        do k = 1, size(in_range)
           if ( atom(in_range(k))%MMSymbol == "OX" ) then
              ref = in_range(k); exit
           end if
        end do
        if( ref == 0 ) ref = in_range(1)  ! no OX in this group → use first atom
        
        do k = 1, size(in_range)
           if ( in_range(k) == ref ) cycle
           dxyz = atom(in_range(k))%xyz - atom(ref)%xyz
           dxyz = dxyz - Txyz * dnint(dxyz / Txyz)
           atom(in_range(k))%xyz = atom(ref)%xyz + dxyz
        end do

        ! centroid of the molecule of residue number nr 
        do i = 1, 3 
           centroid(i) = sum(atom(in_range)%xyz(i)) / size(in_range)
        end do  

        ! move molecule inside the box, if originally outside
        dxyz = Txyz * dnint( centroid / Txyz )
        do k = 1, size(in_range)
           atom(in_range(k))%xyz = atom(in_range(k))%xyz - dxyz
        end do

    end do
    
end associate

end subroutine ReGroup
!
!
!
!=================================================
subroutine translate_to_centroid(system, res_name)
!=================================================
    implicit none
    type(universe), intent(inout)        :: system
    character(*)  , intent(in), optional :: res_name
    
    !local variables
    integer :: i, N_of_solute_atoms
    real*8  :: centroid(3)

if( .not. present(res_name) ) then 

    do i = 1, 3
       centroid(i) = sum(system%atom(:)%xyz(i)) / system%N_of_atoms
       ! translate coordinates to the centroid of the box ...
       system%atom(:)%xyz(i) = system%atom(:)%xyz(i) - centroid(i)
    end do

else

    N_of_solute_atoms = count(system%atom(:)%resid==res_name)                                                                                             
    if( N_of_solute_atoms == 0 ) stop "No solute with this residue name"   

    ! place solute in the center of the PBC box
    do i=1,3
       centroid(i) = sum(system%atom(:)%xyz(i) , system%atom(:)%resid==res_name) / N_of_solute_atoms
       system%atom(:)%xyz(i) = system%atom(:)%xyz(i) - centroid(i)
    end do 

end if

end  subroutine translate_to_centroid
!
!
!
subroutine deallocate_and_leave

    deallocate(O_ptr,             &
               H_ptr,             &
               OH_bond_order,     &
               OH_distance_table, &
               OO_distance_table)

end subroutine deallocate_and_leave
!
!
!
!=========================================================
 subroutine QM_droplet( sys, QM_droplet_radius, res_name )
!=========================================================
implicit none
type(universe), intent(inout)        :: sys
real(8)       , intent(in)           :: QM_droplet_radius
character(*)  , intent(in), optional :: res_name

!local variables ...
integer :: i, nr , nr_max, N_of_solute_atoms, N_of_atoms_in_nr
real*8  :: distance
real*8  :: solvent_CG(3) , solute_CG(3)

sys % atom % fragment = "X"  !> default fragment
if( present(res_name) ) &
then
    where( sys % atom % resid == res_name ) sys % atom % fragment = "D"  !> donor fragment
else
    where( sys % atom % nresid == 1 ) sys % atom % fragment = "D" 
end if

! identify the centroid of the solute ...
N_of_solute_atoms = count(sys%atom%fragment == "D")
if( N_of_solute_atoms == 0 ) then
    stop "solute was not found"
end if

forall( i=1:3 ) solute_CG(i) = sum( sys%atom%xyz(i) , sys%atom%fragment == "D" ) / N_of_solute_atoms

nr_max =  maxval(sys%atom(:)%nresid)

do nr = 1 , nr_max

      N_of_atoms_in_nr = count(sys%atom%nresid==nr)
      if( N_of_atoms_in_nr == 0 ) cycle
 
      forall( i=1:3 ) solvent_CG(i) = sum( sys%atom%xyz(i) , sys%atom%nresid == nr ) / N_of_atoms_in_nr
      
      distance = norm2(solute_CG - solvent_CG)

      if( distance <= QM_droplet_radius ) then
          where( sys%atom%nresid == nr ) 
               sys%atom%fragment = "Q"
               end where
      end if
end do

end subroutine QM_droplet
!
!
!
!==============================================
subroutine eliminate_classical_atoms( system )
!==============================================
implicit none
type(universe), intent(inout) :: system

!local variables
type(universe)       :: temp
integer              :: New_No_of_atoms
logical, allocatable :: quantum_atoms(:)

! mask
quantum_atoms = (system% atom% fragment == "Q" )

New_No_of_atoms = count(quantum_atoms)
allocate( temp%atom( New_No_of_atoms ) )

temp%atom = pack( system%atom, quantum_atoms )

CALL move_alloc(from=temp%atom,to=system%atom)
system%N_of_atoms = New_No_of_atoms

end subroutine eliminate_classical_atoms
!
!
!
!==============================
subroutine pack_HOH_atoms(atom)
!==============================
    implicit none
    type(atomic), intent(inout) :: atom(:)

    ! local variables ...
    integer :: k, nr, ref, HOH_resid
    integer :: droplet_size
    integer :: lowest_nr, highest_nr, offset
    integer, allocatable :: in_range(:)
    type(atomic), allocatable :: aux_atom(:)

    droplet_size = size(atom)
    aux_atom     = atom

    associate( ref_nr => aux_atom%nresid, ref_name => aux_atom%resid )

        lowest_nr  = minval(ref_nr, mask=(ref_name=="HOH"))
        highest_nr = maxval(ref_nr, mask=(ref_name=="HOH"))

        offset = findloc(ref_name, value="HOH", dim=1) - 1

        call get_n_of_S_residues(offset, atom)

        ! highest residue number before the water block; 
        if ( offset==0 ) then
            HOH_resid = 1
        else
            HOH_resid = maxval(ref_nr, mask=(ref_nr < lowest_nr)) + 1
        end if

        do nr = lowest_nr, highest_nr

            ! atoms of THIS water molecule only
            in_range = pack( [(k, k=1,droplet_size)], (ref_nr==nr) .and. (ref_name=="HOH") )

            if ( size(in_range) == 0 ) cycle   ! nr is not a water molecule here; skip it

            ! find the OX of this residue and use it as the unwrap reference ...
            ref = 0
            do k = 1, size(in_range)
                if ( aux_atom(in_range(k))%MMSymbol == "OX" ) then
                    ref         = in_range(1)
                    in_range(1) = in_range(k)
                    in_range(k) = ref
                    exit
                end if
            end do

            if ( ref == 0 ) then
                write(*,'(a,i0)') "ERROR: lone proton found in residue = ", nr
            end if

            ! copy this molecule's atoms into contiguous slots, renumbering nresid in sequence ...
            do k = 1, size(in_range)
                atom(offset+k)        = aux_atom(in_range(k))
                atom(offset+k)%nresid = HOH_resid
            end do

            offset    = offset + size(in_range)
            HOH_resid = HOH_resid + 1

        end do

    end associate

end subroutine pack_HOH_atoms
!
!
!
!==================================
subroutine pack_solvent_atoms(atom)
!==================================
    implicit none
    type(atomic), intent(inout) :: atom(:)

    ! local variables ...
    integer :: k, nr, ref
    integer :: HOH_nr, WAT_nr
    integer :: lowest_nr, highest_nr
    integer :: droplet_size, n_Q_atoms
    integer :: offset, HOH_offset, WAT_offset
    integer     , allocatable :: in_range(:)
    type(atomic), allocatable :: aux_atom(:)

    droplet_size = size(atom)
    aux_atom     = atom

    associate( ref_nr   => aux_atom%nresid , &
               ref_name => aux_atom%resid  )

        lowest_nr  = minval(ref_nr, mask=(ref_name=="HOH"))
        highest_nr = maxval(ref_nr, mask=(ref_name=="HOH"))

        offset = findloc(ref_name, value="HOH", dim=1) - 1
        HOH_offset = offset
        WAT_offset = offset + count( atom%resid=="HOH" .and. atom%fragment=="Q")

        call get_n_of_S_residues(offset, atom)
        if( (n_HOH_resids + n_WAT_resids) /= (highest_nr - lowest_nr + 1) ) then
            stop "ERROR: (n_HOH_resids + n_WAT_resids) /= total number of solvent residues "
        end if

        ! highest residue number before the water block; 
        if ( HOH_offset==0 ) then
            HOH_nr = 1
        else
            HOH_nr = maxval( ref_nr(1:offset) ) + 1
        end if
        WAT_nr = HOH_nr + n_HOH_resids

        do nr = lowest_nr, highest_nr

            ! atoms of THIS water molecule only
            in_range = pack( [(k, k=1,droplet_size)], (ref_nr==nr) .and. (ref_name=="HOH") )

            ! find the OX of this residue and use it as the unwrap reference ...
            ref = 0
            do k = 1, size(in_range)
                if ( aux_atom(in_range(k))%MMSymbol == "OX" ) then
                    ref         = in_range(1)
                    in_range(1) = in_range(k)
                    in_range(k) = ref
                    exit
                end if
            end do

            if ( ref == 0 ) then
                Print*, red_bg("ERROR: lone proton found in residue = "), nr
                stop
            end if

            ! copy this molecule's atoms into contiguous slots, renumbering nresid in sequence ...
            if( aux_atom(in_range(1))%fragment == "Q" ) then
                do k = 1, size(in_range)
                    atom(HOH_offset+k)        = aux_atom(in_range(k))
                    atom(HOH_offset+k)%nresid = HOH_nr
                end do
                HOH_offset = HOH_offset + size(in_range)
                HOH_nr  = HOH_nr + 1
            else
                do k = 1, size(in_range)
                    atom(WAT_offset+k)          = aux_atom(in_range(k))
                    atom(WAT_offset+k)%nresid   = WAT_nr
                    atom(WAT_offset+k)%resid    = "WAT"
                    atom(WAT_offset+k)%fragment = "S"
                    if(atom(WAT_offset+k)%MMSymbol == "OX") atom(WAT_offset+k)%MMSymbol = "OW"
                    if(atom(WAT_offset+k)%MMSymbol == "HX") atom(WAT_offset+k)%MMSymbol = "HW"
                end do
                WAT_offset = WAT_offset + size(in_range)
                WAT_nr = WAT_nr + 1
            end if

        end do

    end associate

    n_Q_atoms = count( atom(offset+1:)%fragment == "Q" )
    ! consistency check: every solvent atom was copied exactly once
    if( HOH_offset /= offset + n_Q_atoms .or. WAT_offset /= droplet_size ) then
        print*, red_bg("ERROR in pack_solvent_atoms: solvent atoms were lost or duplicated")
        error stop 
    end if

end subroutine pack_solvent_atoms
!
!
!
!===========================================
subroutine get_n_of_S_residues(offset, atom)
!===========================================
! Counts the distinct solvent residues in atom(offset+1:), by fragment:
!     fragment "Q"  ->  HOH residues
!     fragment "X"  ->  WAT residues
!-------------------------------------------
    implicit none

    integer     , intent(in)  :: offset
    type(atomic), intent(in)  :: atom(:)

    ! Local variables
    integer :: i, first

    n_HOH_resids = 0
    n_WAT_resids = 0

    first = offset + 1

    if (first > size(atom)) return

    do i = first, size(atom)

        ! Skip residue if its nresid has already been encountered
        ! within the region being analyzed.
        if (i > first) then
            if (any(atom(first:i-1)%nresid == atom(i)%nresid)) cycle
        end if

        ! This is the first occurrence of this residue.
        select case (atom(i)%fragment)

        case ("Q")
            n_HOH_resids = n_HOH_resids + 1

        case ("X")
            n_WAT_resids = n_WAT_resids + 1

        case default
            write(*,'(a,a,a,i0)') &
            "ERROR: unrecognized fragment '", atom(i)%fragment, "' in residue ", atom(i)%nresid
            error stop

        end select

    end do

end subroutine get_n_of_S_residues
!
!
!
end module test_droplet
