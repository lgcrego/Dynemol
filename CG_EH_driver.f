module EH_CG_driver_m

    use type_m
    use constants_m
    use parameters_m            , only : profiling
    use GA_QCModel_m            , only : eval_CG_cost
    use OPT_Parent_class_m      , only : GA_OPT
    use CG_class_m              , only : EH_OPT
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: CG_driver

    private 

    ! module variables ...

contains
!
!
!
!====================================================
 subroutine CG_driver( GA , GA_Selection , CG_basis )
!==================================================== 	
implicit none
type(GA_OPT)                    , intent(in)    :: GA
type(STO_basis)                 , intent(inout) :: GA_Selection(:,:)
type(STO_basis) , allocatable   , intent(out)   :: CG_basis(:)

! local variables ...
integer                 :: i , GlobalMinimum
integer                 :: Top_Selection 
real*8  , allocatable   :: local_minimum(:) , InitialCost(:)
character(len=2)        :: string
character(len=31)       :: f_name
type(EH_OPT)            :: CG

eval_CG_cost = .true.

If( profiling ) OPEN( unit=32 , file='opt_trunk/CG.log.dat' , status='unknown' )

Top_Selection = size(GA_Selection(1,:)) 

allocate( local_minimum(Top_Selection) , source = real_large)
allocate( InitialCost(0:Top_Selection) )

InitialCost(0) = 0.d0

do i = 1 , Top_Selection

    If( profiling ) then
        write(string,'(i2.2)') i
        f_name = 'opt_trunk/CG_OPT_parms'//string//'.dat'
        OPEN( unit = 42 , file=f_name , status='replace' )
    end If

    ! instantiating CG ...
    CG = EH_OPT( GA_Selection(:,i) , GA )

    InitialCost(i) = CG%cost() 

    Print 162 , i , Top_Selection 
    If( CG%cost() /= InitialCost(i-1) ) CALL Fletcher_Reeves_Polak_Ribiere_minimization( CG , GA%GeneSize , local_minimum(i) )

    ! temporarily stores CG optimized basis here ...
    If( local_minimum(i) < InitialCost(i) ) then

        GA_Selection(:,i) = CG%basis

    else ! <== CG minimization failed ...

        local_minimum(i) = InitialCost(i)

    End If

    If( profiling ) close(42)

end do

GlobalMinimum = minloc(local_minimum,dim=1)

do i = 1 , Top_Selection
    print*, InitialCost(i) , local_minimum(i)
end do
print*, GlobalMinimum

allocate( CG_basis (size(CG%basis)) )
CG_basis = GA_Selection(:,GlobalMinimum)

If( profiling ) close(32)

include 'formats.h'

end subroutine CG_driver
!
!
end module EH_CG_driver_m
