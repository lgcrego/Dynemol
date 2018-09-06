module MM_CG_driver_m

    use type_m                  
    use constants_m
    use parameters_m            , only : 
    use MM_input                , only : OPT_driver 
    use MD_read_m               , only : atom 
    use MM_types                , only : MM_atomic , LogicalKey
    use cost_MM                 , only : nmd_REF_erg , nmd_NOPT_erg , KeyHolder , overweight , chi
    use FF_OPT_class_m          , only : FF_OPT , atom0
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              
    use GA_m                    , only : Genetic_Algorithm

    public :: justCG 

    private 

    ! module variables ...
    type(FF_OPT) :: MM_parms

contains
!
!
!
!===========================================================================================
 subroutine justCG( MM_parms , GA_Selection , local_minimum , InitialCost , new_derictives )
!===========================================================================================
implicit none
type(FF_OPT)               , intent(inout) :: MM_parms
real*8       , allocatable , intent(inout) :: GA_Selection(:,:)
real*8                     , intent(inout) :: local_minimum(:)
real*8                     , intent(inout) :: InitialCost(:)
character(*)               , intent(in)    :: new_derictives

! local variables ...
integer          :: i , k 
integer          :: Top_Selection 
real*8           :: this_minimum 
type(LogicalKey) :: key 

Top_Selection = size(local_minimum)

do i = 1 , Top_Selection

    atom = atom0
    ! different KeyHolders(:) allow for sequential optimization ...
    do k = 1 , size(KeyHolder)

        key      =  KeyHolder(k)
        MM_parms =  FF_OPT( key , kernel = "NormalModes" , directives = new_derictives )

        write(*,190) i , KeyHolder(k)% comment , MM_parms% directives

        If( k == 1 ) then

            MM_parms% p    =  GA_Selection(:,i)
            InitialCost(i) =  MM_parms% cost()

        end If

        CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

        If( this_minimum == real_large ) exit

        ! temporarily stores CG optimized FF parameters here ...
        CALL save_temporary_results( GA_Selection(:,i) , MM_parms )
        local_minimum(i) = this_minimum

    end do
    print*, "  ==> done"

end do

include 'formats.h'

end subroutine justCG
!
!
!
!=================================================
 subroutine save_temporary_results( a , MM_parms )
!=================================================
implicit none
real*8       , intent(inout) :: a(:)
type(FF_OPT) , intent(inout) :: MM_parms

! local variables ...
type(LogicalKey) :: key

key = KeyHolder( size(KeyHolder) )

MM_parms = FF_OPT( key , kernel = "JustKey" )

a = MM_parms % p

end subroutine save_temporary_results
!
!
end module MM_CG_driver_m
