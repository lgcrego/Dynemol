module F_inter_m

    use constants_m
    use omp_lib
    use syst               , only : using_barostat
    use md_read_m          , only : MM, molecule
    use F_inter_nonbond    , only : f_inter_nonbonding
    use F_inter_DWFF       , only : f_DWFF_inter
    use Berendsen_Barostat , only : virial_tensor

    public :: FORCEINTER

    private

contains
!
!
!
!=====================
 subroutine FORCEINTER
!=====================
    implicit none
    
    if( using_barostat% inter ) call InitializeStressMatrix
    
    if ( MM % N_of_molecules > 1 ) &
    then
        call f_inter_nonbonding()
    
        if (any(molecule%DWFF)) then
           call f_DWFF_inter()
        endif
        
    end if
    
    if( using_barostat% inter ) call ConcludeStressMatrix

end subroutine FORCEINTER
!
!
!
!
!=================================
 subroutine InitializeStressMatrix
!=================================
    implicit none
    
    ! initializing variables for this integration step ...
    virial_tensor(:,:) = D_zero
    
end subroutine InitializeStressMatrix
!
!===============================
 subroutine ConcludeStressMatrix
!===============================
    implicit none
    
    ! local variables
    integer :: i,j
    
    ! symmetrizing the tensors ...
    virial_tensor = virial_tensor * factor3
    
    do concurrent (i = 1:2, j = 1:3, j>i)
      virial_tensor(j,i) = virial_tensor(i,j) 
    end do

end subroutine ConcludeStressMatrix
!
!
end module F_inter_m
