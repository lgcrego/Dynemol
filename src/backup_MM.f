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
integer         :: i , j 
logical , save  :: first_time = .true. 
logical         :: exist

! check whether restart is properly set ...
If( first_time ) then

    If( restart ) then
        inquire( file="Security_copy_MM.dat", EXIST=exist )   
        If( exist ) stop " <Security_copy_MM.dat> exists; check restart parameter or move Security_copy_MM.dat to Restart_copy_MM.dat"
    else
        inquire( file="Restart_copy_MM.dat", EXIST=exist )
        If( exist ) stop " <Restart_copy_MM.dat> exists; check restart parameter or delete Restart_copy_MM.dat"
    end If

    ! get ride of Restart_copy_MM.dat for new Security_copy_MM.dat ...
    inquire( file="Restart_copy_MM.dat", EXIST=exist )
    If( exist ) CALL system( "rm Restart_copy_MM.dat" )

    first_time = .false.

end If

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
