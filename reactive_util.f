module reactive

    use type_m
    use MM_types
    use color_funcs
    use tuning_m     , only : solvent_residues
    use MD_read_m    , only : atom, MM
    use card_reading , only : solvent_QM_droplet_radius
    
    implicit none    
    private

    public :: DWFF_atom_indices, deal_with_proton_transfer
 
    ! Arrays describing conectivity and species
    integer, allocatable :: O_ptr(:)             ! O-atom indices
    integer, allocatable :: H_ptr(:)             ! H-atom indices
    real(8), allocatable :: OH_distance_table(:,:)

    ! Working copy of the atomic structure
    type(MM_atomic), allocatable :: atom_wrk(:)

    ! module variables
    integer :: n_atoms, nHX, nOX

contains
!
!
!
!==============================================
subroutine DWFF_atom_indices(sys, solvent_mols)
!==============================================
    implicit none
    type(structure)               , intent(in)  :: sys
    type(molecular) , allocatable , intent(out) :: solvent_mols(:)
    
    ! local variables ...
    integer :: i, k, nr
    integer :: nr_atoms
    integer :: lowest_nr, highest_nr
    integer :: N_of_S_Mols

    if( size(sys%nr) /= n_atoms ) &
        error stop "DWFF_atom_indices: sys and atom_wrk have inconsistent sizes"
 
    ! find positions of environment molecules ...
    lowest_nr  = minval(sys%nr, mask=is_solvent(sys%residue))
    highest_nr = maxval(sys%nr, mask=is_solvent(sys%residue))
    
    ! total number of molecules comprising the dielectric domain ...
    N_of_S_Mols = highest_nr - lowest_nr + 1
    
    allocate( solvent_mols(N_of_S_Mols) )
    
    i = 0
    do nr = lowest_nr, highest_nr
        i = i + 1

        solvent_mols(i)% nr = nr

        nr_atoms = count( atom_wrk% nr == nr )
        solvent_mols(i)% N_of_atoms = nr_atoms

        solvent_mols(i)% sys_id = pack( [(k, k = 1, n_atoms)], atom_wrk% nr == nr )
    
        allocate( solvent_mols(i)% PC% Q  (nr_atoms)   )
        allocate( solvent_mols(i)% PC% nr (nr_atoms)   )
        allocate( solvent_mols(i)% PC% xyz(nr_atoms,3) )
    end do

end subroutine DWFF_atom_indices
!
!
!
!===================================
subroutine deal_with_proton_transfer
!===================================
    implicit none

    ! local variables
    logical, save :: done = .false.

    if( .not. done ) then
        call setup
        done = .true. 
    end if

    ! Keep an independent copy for topology operations
    atom_wrk = atom

    call bond_topology( MM%box )

end subroutine deal_with_proton_transfer
!
!
!
!===============================
subroutine bond_topology( Txyz )
!===============================
    ! Operates on module-level atom_wrk (mutates it in place) ...
    implicit none
    real(8) , intent(in) :: Txyz(3)

    ! Local variables
    integer :: i, j
    real(8) :: rij(3), rijlen
    integer, allocatable :: OH_pair(:,:)
    integer, allocatable :: OH_bond_order(:)   ! bond order per-water

    ! ------------------------------------------------
    ! Compute O–H distances (minimum image convention)
    ! ------------------------------------------------
    OH_distance_table = 0.0d0

    associate (OX => atom_wrk(O_ptr), HX => atom_wrk(H_ptr))
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
    end associate

    ! --------------------------------------------------------
    ! Build OH_pair adjacency matrix
    ! given HX, find the closest OX
    ! --------------------------------------------------------
    allocate( OH_pair(nHX,2) )

    do j = 1, nHX
       i = minloc( OH_distance_table(:, j), dim=1 )
       OH_pair(j, 1) = O_ptr(i)
       OH_pair(j, 2) = H_ptr(j)

       ! Assign each hydrogen's nr to its current oxygen owner's nr
       atom_wrk( H_ptr(j) )% nr = atom_wrk( O_ptr(i) )% nr
    end do

    ! --------------------------------------------------------
    ! Count neighbors for each oxygen
    ! --------------------------------------------------------
    allocate(OH_bond_order(nOX), source=0)
    do i = 1, nOX
        OH_bond_order(i) = count( OH_pair(:, 1) == O_ptr(i) )
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

end subroutine bond_topology
!
!
!
!===============
subroutine setup
!===============
    implicit none

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

    ! Allocate module-level work array
    allocate(OH_distance_table(nOX, nHX), source=0.0d0)

end subroutine setup
!
!
!
elemental logical function is_solvent(resname) result(solvent)
   character(len=*), intent(in) :: resname
   solvent = any(resname == solvent_residues)
end function is_solvent
!
!
!
!==============================================
subroutine ReGroup_Reactive_Molecules(atom_wrk)
!==============================================
implicit none
type(MM_atomic), intent(inout) :: atom_wrk(:)

!local variables
integer :: nr, i, k, ref
integer :: nr_min, nr_max 
real*8  :: dxyz(3), Txyz(3), centroid(3)
integer, allocatable :: in_range(:)

Txyz = MM%box

associate( atom => atom_wrk )
    nr_min = minval(atom%nr)
    nr_max = maxval(atom%nr)

    do nr = nr_min, nr_max

        in_range = pack( [(k, k=1,MM% N_of_atoms)], atom% nr == nr )

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

end subroutine ReGroup_Reactive_Molecules
!
!
!
end module reactive
