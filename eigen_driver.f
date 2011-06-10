! Subroutine for computing time evolution through time slices
module Eigen_driver_m

    use type_m
    use constants_m
    use parameters_m                , only : survival , driver , n_part
    use FMO_m                       , only : eh_tag 
    use Data_Output                 , only : Dump_stuff 
    use Schroedinger_m              , only : DeAllocate_QDyn
    use AO_adiabatic_m              , only : AO_adiabatic
    use ElHl_adiabatic_m            , only : ElHl_adiabatic

    public :: Eigen_driver

    private

contains
!
!
!=======================
 subroutine Eigen_driver
!=======================
implicit none

! local variables ...
integer                :: it 
real*8  , allocatable  :: QDyn_temp(:,:,:)
type(f_time)           :: QDyn

If( .NOT. survival ) then
    Print*, ">>> Survival = .FALSE. <<<"
    Stop
end If

select case ( DRIVER )

    case( "slice_MO0" )


    case( "slice_MOt" )


    case( "slice_ElHl" )

        CALL ElHl_adiabatic ( QDyn , it )

    case( "slice_AO" )

        CALL AO_adiabatic ( QDyn , it )

    case default
                Print*, " >>> Check your Eigen_driver options <<< :" , driver

end select

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 , n_part ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 , : ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn ) 

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

include 'formats.h'

end subroutine eigen_driver
!
!
!
end module Eigen_driver_m
