module OPT_Parent_class_m

    implicit none
    private

    type, public :: OPT_Parent
        integer                 :: N_of_Freedom 
        integer                 :: ITMAX 
        real*8                  :: BracketSize 
        real*8                  :: accuracy
        real*8  , allocatable   :: p(:)
        character (len=11)      :: driver
        character (len=72)      :: message
        logical                 :: profiling 
    contains
        procedure :: cost 
        procedure :: cost_variation
        procedure :: output 
    end type 


    type, extends(OPT_Parent), public :: GA_OPT
        real*8                          :: DP(3)
        real*8          , allocatable   :: erg(:)
        integer         , allocatable   :: key(:,:)
        integer                         :: GeneSize
        character(3)    , allocatable   :: EHSymbol(:)               
    end type GA_OPT


    contains

        function cost( me ) result(my_cost)
            class(OPT_Parent) , intent(inout) :: me
            real*8                            :: my_cost 
        end function

        subroutine cost_variation( me , vector ) 
            class(OPT_Parent) , intent(in)    :: me
            real*8            , intent(inout) :: vector(:) 
        end subroutine
        
        subroutine output( me , iter) 
            class(OPT_Parent) ,            intent(in) :: me
            Integer           , optional , intent(in) :: iter
        end subroutine


end module OPT_Parent_class_m
