module CG_driver_m

    use type_m
    use constants_m
    use CG_class_m              , only : CG_OPT
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: CG_driver

    private 

    ! module variables ...
    integer      :: Top_Selection = 15
    type(CG_OPT) :: CG

contains
!
!
!
!=============================================================
 subroutine CG_driver( system , GA , GA_Selection , CG_basis )
!============================================================= 	
implicit none
type(structure)                 , intent(in)    :: system
type(OPT)                       , intent(in)    :: GA
type(STO_basis)                 , intent(inout) :: GA_Selection(:,:)
type(STO_basis) , allocatable   , intent(out)   :: CG_basis(:)

! local variables ...
integer                 :: i , GlobalMinimum
real*8  , allocatable   :: local_minimum(:) , InitialCost(:)

Top_Selection = min( Top_Selection , size(GA_Selection(1,:)) )

allocate( local_minimum(Top_Selection) , source = real_large)
allocate( InitialCost(0:Top_Selection) )

InitialCost(0) = 0.d0

do i = 1 , Top_Selection

    ! instantiating CG ...
    CG = CG_OPT( GA_Selection(:,i) , GA )

    InitialCost(i) = CG%cost() 

    If( CG%cost() /= InitialCost(i-1) ) CALL Fletcher_Reeves_Polak_Ribiere_minimization( CG , GA%GeneSize , local_minimum(i) )

    ! temporarily stores CG optimized basis here ...
    GA_Selection(:,i) = CG%basis

end do

GlobalMinimum = minloc(local_minimum,dim=1)

do i = 1 , Top_Selection
    print*, InitialCost(i) , local_minimum(i)
end do
print*, GlobalMinimum

allocate( CG_basis (size(CG%basis)) )
CG_basis = GA_Selection(:,GlobalMinimum)

end subroutine CG_driver
!
!
end module CG_driver_m
