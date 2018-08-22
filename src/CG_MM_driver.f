module MM_CG_driver_m

    use type_m                  
    use constants_m
    use parameters_m            , only : 
    use MM_input                , only : OPT_driver , N_of_CGSteps
    use MD_read_m               , only : atom 
    use MM_types                , only : MM_atomic , LogicalKey
    use cost_MM                 , only : nmd_REF_erg , nmd_NOPT_erg , KeyHolder , overweight , chi
    use FF_OPT_class_m          , only : FF_OPT , atom0
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              
    use GA_m                    , only : Genetic_Algorithm

    public :: justCG , CGAd , CGRc 

    private 

    ! module variables ...
    type(FF_OPT) :: MM_parms

contains
!
!
!
!==========================================================================
 subroutine justCG( MM_parms , GA_Selection , local_minimum , InitialCost )
!==========================================================================
implicit none
type(FF_OPT)               , intent(inout) :: MM_parms
real*8       , allocatable , intent(inout) :: GA_Selection(:,:)
real*8                     , intent(inout) :: local_minimum(:)
real*8                     , intent(inout) :: InitialCost(:)

! local variables ...
integer          :: i , k 
integer          :: Top_Selection 
real*8           :: this_minimum 
type(LogicalKey) :: key 

Top_Selection = size(local_minimum)

do i = 1 , Top_Selection

    atom = atom0
    do k = 1 , size(KeyHolder)

        key      =  KeyHolder(k)
        MM_parms =  FF_OPT( key , kernel = "NormalModes" , directives = "use_no_weights" )

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
!============================================
 subroutine CGRc ( MM_parms , GA_Selection )
!============================================
implicit none
type(FF_OPT)               , intent(inout) :: MM_parms
real*8       , allocatable , intent(inout) :: GA_Selection(:,:)

! local variables ...
integer               :: i , k 
integer               :: Top_Selection 
real*8                :: this_minimum , setting_cost , signal = D_one
real*8  , allocatable :: previous_p(:)
type(LogicalKey)      :: key 

! local parameters ...
real*8 :: MaxOverweight = FIVE

Top_Selection = size(GA_Selection(1,:))

key  = KeyHolder(1)

MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "use_overweight" )

MM_parms% p = GA_Selection(:,1)

allocate( previous_p , source = GA_Selection(:,1) )

do k = 1 , N_of_CGSteps

    atom = atom0

    write(*,190) k , KeyHolder(1)% comment , MM_parms% directives

    CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

    If( this_minimum == real_large ) Then 
        Print*, ">>> CG Optimization failed , moving one step back and resuming OPT " 
        MM_parms% p  = previous_p 
        setting_cost = MM_parms% cost()
    EndIf

    forall(i=1:size(nmd_REF_erg)) overweight(i) = overweight(i) + signal*sqrt(abs(chi(i)/nmd_REF_erg(i)))

    overweight = merge(D_zero     , overweight    , overweight < D_zero       )
    overweight = merge(overweight , MaxOverweight , overweight < MaxOverweight)

    CALL Genetic_Algorithm( MM_parms , GA_Selection , directives = "use_overweigth" )

    If( k == N_of_CGSteps ) Exit

    previous_p = GA_Selection(:,1)

    signal = merge( D_one , -D_one , k < N_of_CGSteps/2 )

end do
Print 194  ! ==> done with Recursive steps

include 'formats.h'

end subroutine CGRc 
!
!
!
!============================================
 subroutine CGAd ( MM_parms , GA_Selection )
!============================================
implicit none
type(FF_OPT)               , intent(inout) :: MM_parms
real*8       , allocatable , intent(inout) :: GA_Selection(:,:)

! local variables ...
integer          :: k 
integer          :: Top_Selection 
real*8           :: this_minimum 
type(LogicalKey) :: key 

Top_Selection = size(GA_Selection(1,:))

key  = KeyHolder(1)

do k = 1 , N_of_CGSteps

    atom = atom0

    MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "proprocess_adiabatic"         )

    MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "use_overweights_adiabatic_OPT" )

    write(*,190) k , KeyHolder(1)% comment , MM_parms% directives

    CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

    CALL Genetic_Algorithm( MM_parms , GA_Selection , directives = "use_overweigths" )

    this_minimum = MM_parms% cost()
    If( this_minimum == real_large ) exit

end do
Print*, "  ==> done"

include 'formats.h'

end subroutine CGAd 
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
