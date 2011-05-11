module execution_time_m

public :: start_clock , stop_clock

private

! module variables ...
real , save  :: start_time

contains
!
!
!=======================
 subroutine start_clock
!=======================
implicit none

! local variables ...
integer :: time_array_0(8)

call date_and_time(values=time_array_0)

start_time = time_array_0(5)*3600 + time_array_0(6)*60 + time_array_0(7) + 0.001*time_array_0(8)

end subroutine start_clock
!
!
!=====================
 subroutine stop_clock
!=====================
implicit none

! local variables ...
integer :: time_array_0(8)
real    :: finish_time

call date_and_time(values=time_array_0)

finish_time = time_array_0(5)*3600 + time_array_0(6)*60 + time_array_0(7) + 0.001*time_array_0(8)

print*, " >> Execution Time :", finish_time - start_time 

end subroutine stop_clock

end module execution_time_m
