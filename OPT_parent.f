module OPT_Parent_class_m

    implicit none
    private

    type , public :: OPT_Parent
        integer                 :: N_of_Freedom 
        integer                 :: ITMAX 
        real*8                  :: BracketSize 
        real*8  , allocatable   :: p(:)
        character (len=11)      :: driver
        logical                 :: profiling 
    contains
        procedure :: cost 
        procedure :: cost_variation
        procedure :: output 
    end type 


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
