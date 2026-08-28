module test_droplet

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

    call post_processing_analysis

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
    integer :: out1  
    real(8) :: cutoff_radius
    character(3) :: res_name = "None"
    character(1) :: YorN
    
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

    call sort_HOH_residue_numbers(work_sys, res_name)

    !-------------------------------------------
    !             Write seed 
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
!====================================================
subroutine sort_HOH_residue_numbers(system, res_name)
!====================================================
!
! Assigns a consistent, sequential nresid to each water molecule, using
! the current OH connectivity (OH_pair) rather than file position, so it
! is correct even when proton transfer has scrambled the OX,HX,HX order.
! Every oxygen and the hydrogens currently bonded to it share one nresid.
!
    implicit none
    type(universe), intent(inout) :: system
    character(*)  , intent(in), optional :: res_name

    integer :: io, max_donor, OX_offset

    ! start water numbering past the highest solute residue number
    if ( any(system%atom(:)%resid /= "HOH") ) then
        max_donor = findloc(system%atom(:)%resid == res_name, value = .true., dim = 1, back  = .true.)
    else
        max_donor = 0
    end if

    OX_offset = O_ptr(1) - max_donor

    select case( OX_offset)
         case(1)
              do io = 1, nOX
                  system%atom( O_ptr(io)+1 )%nresid = system%atom(O_ptr(io))%nresid
                  system%atom( O_ptr(io)+2 )%nresid = system%atom(O_ptr(io))%nresid
              end do

         case(2)
              do io = 1, nOX
                  system%atom( O_ptr(io)-1 )%nresid = system%atom(O_ptr(io))%nresid
                  system%atom( O_ptr(io)+1 )%nresid = system%atom(O_ptr(io))%nresid
              end do

         case(3)
              do io = 1, nOX
                  system%atom( O_ptr(io)-2 )%nresid = system%atom(O_ptr(io))%nresid
                  system%atom( O_ptr(io)-1 )%nresid = system%atom(O_ptr(io))%nresid
              end do
    end select   

end subroutine sort_HOH_residue_numbers
!
!
!
end module test_droplet
