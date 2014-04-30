module Backup_MM_m

    use parameters_m    , only : restart
    use MM_types        , only : MM_system , MM_atomic

    public  :: Security_Copy_MM , Restart_MM

contains
!
!                                                                    
!
!===============================================
subroutine Security_Copy_MM( MM , atom , frame )
!===============================================
implicit none
type(MM_system) , intent(in)    :: MM
type(MM_atomic) , intent(in)    :: atom(:)
integer         , intent(in)    :: frame

! local variables ...
integer         :: i , j , k , file_err
logical , save  :: first_time = .true. 

! check whether restart is properly set ...
if( (.NOT. restart) .AND. first_time ) then
    open(unit=3, file="Restart_copy_MM.dat", status="new", iostat=file_err , err=11 )
    11 if( file_err > 0 ) stop " <Restart_copy_MM.dat> exists; check restart parameter "
    close(3)
    CALL system( "rm Restart_copy_MM.dat" )
    first_time = .false.
end if

open(unit=33, file="Security_copy_MM.dat", status="unknown", form="unformatted", action="write")

write(33) frame

write(33) MM % N_of_atoms

do i = 1 , 3
    write(33) MM % box(i)
end do

do i = 1 , size(atom)
   write(33) atom(i) % charge
   write(33) atom(i) % MM_charge
   do j = 1 , 3
        write(33) atom(i) % xyz(j)
        write(33) atom(i) % vel(j)
        write(33) atom(i) % ftotal(j)
    end do
end do

close( 33 )

! erase restart_copy file ...
CALL system( "rm Restart_copy_MM.dat 2> qdynamo.err" )

end subroutine Security_Copy_MM
!
!
!
!=========================================
subroutine Restart_MM( MM , atom , frame )
!=========================================
implicit none
type(MM_system) , intent(inout) :: MM
type(MM_atomic) , intent(inout) :: atom(:)
integer         , intent(out)   :: frame

! local variables ...
integer :: i , j , file_err

open(unit=33, file="Restart_copy_MM.dat", form="unformatted", status="old", action="read" , iostat=file_err , err=11 )

read(33) frame

read(33) MM % N_of_atoms

do i = 1 , 3
    read(33) MM % box(i)
end do

do i = 1 , size(atom)
    read(33) atom(i) % charge
    read(33) atom(i) % MM_charge
    do j = 1 , 3
        read(33) atom(i) % xyz(j)
        read(33) atom(i) % vel(j)
        read(33) atom(i) % ftotal(j)
    end do
end do

close( 33 )

11 if( file_err > 0 ) stop " <Restart_copy_MM.dat> file not found; terminating execution; mv Security_copy_MM.dat Restart_copy_MM.dat "

end subroutine Restart_MM
!
!
!
end module Backup_MM_m
