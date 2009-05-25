 module orbitals_m

    use type_m
    use EHT_parameters
    use Structure_Builder

 contains
!
!
!
 subroutine orbitals( system, basis )


 molecule_basis = sum( atom( system%AtNo(system%molecule:system%atoms) )%DOS )

 Molecule_Basis_Size = count( basis%atom >= system%molecule )


 If( Orbital_Select ) then

      FMO(0)%orbital = HOMO_state    ; FMO(0)%spin = + 1 

      FMO(1)%orbital = initial_state ; FMO(1)%spin = + 1 

  else

      do i = 1 , 






!
!
 end module orbitals_m
