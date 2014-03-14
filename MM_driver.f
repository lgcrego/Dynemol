module MMechanics_m

    use type_m
    use constants_m
    use parameters_m            , only : t_i , n_t , t_f , frame_step
    use MM_dynamics_m           , only : MolecularMechanics , preprocess_MM

    public :: MMechanics

    private

contains
!
!
!
!=====================
 subroutine MMechanics
!=====================
implicit none

! local variables ...
real*8  :: t , t_rate 
integer :: it , frame , frame_init , frame_final

it = 0
t  = t_i

CALL preprocess_MM

frame_init = 2

!-------------------------------------------------------

! time is PICOseconds in EHT & seconds in MM ... 
t_rate      = t_f / float(n_t)
frame_final = n_t

do frame = frame_init , frame_final , frame_step

    t = t + t_rate 

    If( it >= n_t .OR. t >= t_f + t_rate ) exit

    it = it + 1

    CALL MolecularMechanics( t_rate , frame )

end do

end subroutine MMechanics
!
!
!
end module MMechanics_m
