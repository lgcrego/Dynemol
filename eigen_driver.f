! Subroutine for computing time evolution through time slices
module QMDynamicSlice_driver_m

    use type_m
    use constants_m
    use parameters_m                , only : survival , driver , n_part , SOC
    use Data_Output                 , only : Dump_stuff 
    use Schroedinger_m              , only : DeAllocate_QDyn
    use AO_adiabatic_m              , only : AO_adiabatic
    use Chebyshev_driver_m          , only : Chebyshev_driver

    public :: QMDynamicSlice_driver

    private

contains
!
!
!================================
 subroutine QMDynamicSlice_driver
!================================
implicit none

! local variables ...
integer                :: it , n_spin
real*8  , allocatable  :: QDyn_temp(:,:,:,:)
type(f_time)           :: QDyn

If( .NOT. survival ) then
    Print*, ">>> Survival = .FALSE. <<<"
    Stop
end If

select case ( DRIVER )

    case( "slice_AO" , "slice_FSSH" )

        CALL AO_adiabatic ( QDyn , it )

    case( "slice_Cheb" )

        CALL Chebyshev_driver ( QDyn , it )

    case default
                Print*, " >>> Check your QMDynamicSlice_driver options <<< :" , driver

end select

n_spin = merge(2,1,SOC)

! this maneuver is used for "restart-jobs", because it/=nt at the end in this case ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 , n_part , n_spin ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 , : , : ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn ) 

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

include 'formats.h'

end subroutine QMDynamicSlice_driver
!
!
!
end module QMDynamicSlice_driver_m
