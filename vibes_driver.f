module vibrational_modes_m

    use type_m
    use constants_m
    use MM_ERG_class_m          , only : CG_OPT
    use NonlinearMM_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: Optimize_Structure

    private 

    ! module variables ...
    type(CG_OPT) :: MM

contains
!
!
!
!=============================================================
 subroutine Optimize_Structure( )
!============================================================= 	
implicit none

! local variables ...
integer :: i , GlobalMinimum
real*8  :: local_minimum , InitialERG
real*8  ,  allocatable :: forces(:)

! instantiating MM ...
MM = CG_OPT( )

InitialERG = MM%cost() 

allocate(forces(MM%N_OF_Freedom) , source=D_zero)

call MM%cost_variation( forces )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM , MM%N_of_Freedom , local_minimum )

end subroutine Optimize_Structure
!
!
end module vibrational_modes_m
