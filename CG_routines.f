module CG_driver_m

    use type_m
    use constants_m
    use CG_class_m              , only : CG_OPT
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: CG_driver

    private 

    ! module variables ...
    type(CG_OPT) :: CG

contains
!
!
!
!=========================================================
 subroutine CG_driver( system , GA , GA_basis , CG_basis )
!========================================================= 	
implicit none
type(structure)                 , intent(in)  :: system
type(OPT)                       , intent(in)  :: GA
type(STO_basis)                 , intent(in)  :: GA_basis(:)
type(STO_basis) , allocatable   , intent(out) :: CG_basis(:)

! local variables ...
real*8 :: custo , local_minimum


! instantiating CG ...
CG = CG_OPT( GA_basis , GA )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( CG , GA%GeneSize , local_minimum )

allocate( CG_basis (size(CG%basis)) )

CG_basis = CG%basis
















end subroutine CG_driver
!
!
!
!
!
!
end module CG_driver_m
