module F_inter_m

    use constants_m
    use omp_lib
    use syst            , only : using_barostat
    use md_read_m       , only : MM
    use F_inter_nonbond , only : f_inter_nonbonding, virial_tensor 

    public :: FORCEINTER

    private

    ! module variables ...
    real*8 , save :: virial_tensor(3,3)
    
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
 virial_tensor(:,:)   = D_zero

end subroutine InitializeStressMatrix
!
!===============================
 subroutine ConcludeStressMatrix
!===============================
 implicit none

 ! local variables
 integer :: i,j

 ! symmetrizing the tensors ...
 virial_tensor   = virial_tensor * factor3

 do concurrent (i = 1:2, j = 1:3, j>i)
   virial_tensor(j,i)   = virial_tensor(i,j) 
 end do

end subroutine ConcludeStressMatrix
!
!
end module F_inter_m
