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
        MM_parms =  FF_OPT( key , kernel = "NormalModes" , directives = "use_overweights" )

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
integer          :: i , k 
integer          :: Top_Selection 
real*8           :: this_minimum , signal = D_one
type(LogicalKey) :: key 

! local parameters ...
real*8           :: MaxOverweight = 5.d0

Top_Selection = size(GA_Selection(1,:))

key  = KeyHolder(1)

MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "use_overweight" )

MM_parms% p = GA_Selection(:,1)

do k = 1 , N_of_CGSteps

    atom =  atom0

    write(*,190) k , KeyHolder(1)% comment , MM_parms% directives

    CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

    forall(i=1:size(nmd_REF_erg)) overweight(i) = overweight(i) + signal*abs(chi(i)/nmd_REF_erg(i))

    overweight = merge(overweight , MaxOverweight , overweight < MaxOverweight)

    CALL Genetic_Algorithm( MM_parms , GA_Selection , directives = "use_overweigth" )
    this_minimum = MM_parms% cost()

    If( this_minimum == real_large ) Then ; Print*, ">>> Recursive Optimization failed " ; stop ; EndIf

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

    atom =  atom0

    MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "proprocess_adiabatic"         )

    MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = "use_no_weights_adiabatic_OPT" )

    write(*,190) k , KeyHolder(1)% comment , MM_parms% directives

    CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

    CALL Genetic_Algorithm( MM_parms , GA_Selection , directives = "use_no_weigths_LineUp" )

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
!
!
!
!
!=========================================================
!
!=========================================================
!
!
!
!
module EH_CG_driver_m

    use type_m
    use constants_m
    use parameters_m            , only : profiling
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
character(len=21)       :: f_name
type(EH_OPT)            :: CG

If( profiling ) OPEN( unit=32 , file='CG-log.dat' , status='unknown' )

Top_Selection = size(GA_Selection(1,:)) 

allocate( local_minimum(Top_Selection) , source = real_large)
allocate( InitialCost(0:Top_Selection) )

InitialCost(0) = 0.d0

do i = 1 , Top_Selection

    If( profiling ) then
        write(string,'(i2.2)') i
        f_name = 'CG_OPT_parms'//string//'.dat'
        OPEN( unit = 42 , file=f_name , status='replace' )
    end If

    ! instantiating CG ...
    CG = EH_OPT( GA_Selection(:,i) , GA )

    InitialCost(i) = CG%cost() 

    If( CG%cost() /= InitialCost(i-1) ) CALL Fletcher_Reeves_Polak_Ribiere_minimization( CG , GA%GeneSize , local_minimum(i) )

    ! temporarily stores CG optimized basis here ...
    GA_Selection(:,i) = CG%basis

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

end subroutine CG_driver
!
!
end module EH_CG_driver_m
