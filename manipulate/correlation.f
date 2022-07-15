module Correlation_m

public :: Correlation_driver

private

contains
!
!
!
!============================
subroutine Correlation_driver
!============================
implicit none
character(1) :: operation

write(*,'(/a)') ' (a) = Auto-correlation of population '
write(*,'(/a)') ' (b) = Auto-correlation of dipole moment'
write(*,'(/a)') ' (c) = Correlation between populations '
write(*,'(/a)') ' (d) = Correlation between dipoles '
write(*,'(/a)') ' (e) = Vibration between populations '
write(*,'(/a)') ' (f) = Vibration between dipoles '
write(*,'(/a)') ' (g) = Generate random survival '
write(*,'(/a)') ' (h) = Auto-correlation of Energy Configurations'
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') operation

select case( operation )

    case( 'a' )
        CALL  Auto_Correlation_Population

    case( 'b' )
        CALL Auto_Correlation_Dipole

    case( 'c' )
        CALL Correlation_Population

    case( 'A' )
        CALL Correlation_Dipole

    case( 'e' )
        CALL Vibration_Population

    case( 'f' )
        CALL Vibration_Dipole

    case( 'g' )
        CALL Generate_Random

    case( 'h' )
        CALL AutoCorrl_ErgConfig

end select

end subroutine Correlation_driver
!
!
!
!================================
subroutine Correlation_Population
!================================
implicit none
real*8  , allocatable   :: survival_1(:) , Survival_2(:) , Corrl(:) , time(:) , time_1(:)
real*8                  :: weight , med_1 , med_2 , med_3 , med_4 , med_5 , med_6 , med_7
integer                 :: i , n
character               :: input_1 , input_2

print'("")'
write(*, '(1x,a)' ) "Input file = survival.dat"
write(*, '(1x,a)' ) "Output file = Correlation_population.dat"
print'("")'
write(*, '(1x,a)' ) "Enter with data for correlation of population"
write(*, '(1x,a)' ) " 1 ==> LIG0 = PPH "
write(*, '(1x,a)' ) " 2 ==> LIG1 = C60 "
write(*, '(1x,a)' ) " 3 ==> LIG2 = CAR "    
write(*, '(1x,a)' ) " 4 ==> dynamic variable"    
print'("")'
write(*, '(a17)', advance = 'no') ' Input 1 = '
read*, input_1
write(*, '(a17)', advance = 'no') ' Input 2 = '
read*, input_2

select case( input_1 )

    case( '1' )
        CALL Read_DONOR( time , survival_1 )

    case( '2' )
        CALL Read_LIG1( time , survival_1 )

    case( '3' )
        CALL Read_LIG2( time , survival_1 )

    case( '4' )
        CALL Read_Dyn( survival_1 )

end select

select case( input_2 )

    case( '1' )
        CALL Read_DONOR( time_1 , survival_2 )

    case( '2' )
        CALL Read_LIG1( time_1 , survival_2 )

    case( '3' )
        CALL Read_LIG2( time_1 , survival_2 )

    case( '4' )
        CALL Read_Dyn( survival_2 )

end select

n = size(survival_1)

allocate( Corrl ( n ) )

! calcule the correlation between the populations ...
do i = 1 , n

    med_1 = sum( survival_1(1:i) * survival_2(1:i) ) / float(i)
    
    med_2 = sum( survival_1(1:i) ) / float(i)
    med_3 = sum( survival_2(1:i) ) / float(i)

    med_4 = sum( survival_1(1:i) * survival_1(1:i) ) / float(i)
    med_5 = med_2 * med_2

    med_6 = sum( survival_2(1:i) * survival_2(1:i) ) / float(i)
    med_7 = med_3 * med_3

    weight = sqrt( ( med_4 - med_5 ) * ( med_6 - med_7 ) )

    Corrl(i) = ( med_1 - med_2 * med_3 ) / weight

end do

! save the results ...
open(unit=113, file="Correlation_population.dat", status="unknown", action="write")

do i = 2 , n

   write( unit=113 , fmt=120 ) time(i) , Corrl(i)
   
end do

close( 113 )

deallocate( survival_1 , survival_2 , Corrl , time )

!---------------------------------------------------------------------

120 format(f10.4 , t15 , f10.4)

end subroutine Correlation_Population
!
!
!
!============================
subroutine Correlation_Dipole
!============================
implicit none
real*8  , allocatable   :: dip_1(:) , dip_2(:) , Corrl(:) , time_1(:) , time_2(:) , time(:)
real*8                  :: weight , med_1 , med_2 , med_3 , med_4 , med_5 , med_6 , med_7
integer                 :: i , size_1 , size_2 , inputstatus , n
character(6)            :: word

print'("")'
write(*, '(1x,a)' ) "Input file = Dipole_1.dat and Dipole_2.dat"
write(*, '(1x,a)' ) "Output file = Correlation_dipole.dat"
print'("")'

! verificate size of system 1 ...
open(unit=33, file="Dipole_1.dat", status="old", action="read")

size_1 = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_1 = size_1 + 1

end do

size_1 = size_1 - 1

rewind 33

! verificate size of system 2 ...
open(unit=35, file="Dipole_2.dat", status="old", action="read")

inputstatus = 0
size_2      = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=35 , fmt=100 , iostat = inputstatus ) word
    size_2 = size_2 + 1

end do

size_2 = size_2 - 1

rewind 35

allocate( dip_1  ( size_1 ) )
allocate( dip_2  ( size_2 ) )
allocate( time_1 ( size_1 ) )
allocate( time_2 ( size_2 ) )

! read file 1 ...
do i = 1 , size_1

    read( unit=33 , fmt=330 ) time_1(i) , dip_1(i)

end do

close(33)

! read file 2 ...
do i = 1 , size_2

    read( unit=35 , fmt=330 ) time_2(i) , dip_2(i)

end do

close(35)

n = size_1

if( n > size_2 ) n = size_2

allocate( Corrl ( n ) )
allocate( time  ( n ) )

time(:) = time_1(1:size_1)

deallocate( time_1 , time_2 )

! calcule the correlation between the dipoles ...
do i = 1 , n

    med_1 = sum( dip_1(1:i) * dip_2(1:i) ) / float(i)
    
    med_2 = sum( dip_1(1:i) ) / float(i)
    med_3 = sum( dip_2(1:i) ) / float(i)

    med_4 = sum( dip_1(1:i) * dip_1(1:i) ) / float(i)
    med_5 = med_2 * med_2

    med_6 = sum( dip_2 (1:i) * dip_2(1:i) ) / float(i)
    med_7 = med_3 * med_3

    weight = sqrt( ( med_4 - med_5 ) * ( med_6 - med_7 ) )

    Corrl(i) = ( med_1 - med_2 * med_3 ) / weight

end do

! save the results ...
open(unit=113, file="Correlation_dipole.dat", status="unknown", action="write")

do i = 2 , n

   write( unit=113 , fmt=120 ) time(i) , Corrl(i)

end do

close( 113 )

deallocate( dip_1 , dip_2 , Corrl , time )

!---------------------------------------------------------------------

100 format(a6)
120 format(f9.4 , t19 , f10.4)
330 format(f9.4 , t40 , f10.5)

end subroutine Correlation_Dipole
!
!
!
!=============================
subroutine AutoCorrl_ErgConfig
!=============================
implicit none
real*8  , allocatable   :: Erg_Config(:) , Corrl(:) , Erg_tmp(:) , time_tmp(:) , time(:)
real*8                  :: sum1 , sum2 , sum3
integer                 :: i , inputstatus , size_sys
character(15)           :: word

print'("")'
write(*, '(1x,a)' ) "Input file = Erg_config.dat"
write(*, '(1x,a)' ) "Output file = AutoCorrl_ErgConfig.dat"
print'("")'

! read file ...
open(unit=33, file="Erg_config.dat", status="old", action="read")

allocate( Erg_tmp  ( 2500000 ) )
allocate( time_tmp ( 2500000 ) )

size_sys = 0
do

    if( inputstatus /= 0 ) exit

    read( unit=33 , fmt=100 , iostat = inputstatus ) word

    if( word == "           Step") then
        size_sys = size_sys + 1
        read( unit=33 , fmt=140 ) Time_tmp(size_sys)
    end if

    if( word == "   Coul. recip.") read( unit=33 , fmt=130 ) Erg_tmp(size_sys)

end do

close(33)

size_sys = size_sys - 2

allocate( Erg_config ( size_sys ) )
allocate( time       ( size_sys ) )
allocate( Corrl      ( size_sys ) )

Erg_config(:) = Erg_tmp(1:size_sys)
time(:)       = time_tmp(1:size_sys)

deallocate( Erg_tmp , time_tmp )

! calcule correlation energy ...
Corrl = 0.0d0

open(unit=20, file="AutoCorrl_ErgConfig.dat", status="unknown", action="write")

do i = 0 , size_sys-1

    sum1 = sum( Erg_config(1:size_sys-i) * Erg_config(i+1:size_sys) ) / ( size_sys - i )
    sum2 = sum( Erg_config(1:size_sys-i) ) / ( size_sys - i )
    sum3 = sum( Erg_config(i+1:size_sys) ) / ( size_sys - i )

    Corrl(i+1) = sum1 - sum2 * sum3

    write( unit=20 , fmt=120 ) time(i+1) , Corrl(i+1) / Corrl(1)
   
end do

close(20)

deallocate( Erg_config , Corrl , time )

!------------------------------------------------------------------------------------

100 format(a15)
120 format(f15.5 , t20 , f10.4)
130 format(t46 , E15.5E2 )
140 format(t16 , f15.5 )

end subroutine AutoCorrl_ErgConfig
!
!
!
!=====================================
subroutine Auto_Correlation_Population
!=====================================
implicit none
real*8  , allocatable   :: survival(:) , Corrl(:) , time(:)
real*8                  :: weight , sum1 , sum2 , sum3
integer                 :: i , counter , tau , n
character               :: input

print'("")'
write(*, '(1x,a)' ) "Input file = survival.dat"
write(*, '(1x,a)' ) "Output file = Auto_Correlation_population.dat"
print'("")'
write(*, '(1x,a)' ) "Enter with data for correlation of population"
write(*, '(1x,a)' ) " 1 ==> ACCEPTOR; PPH >> ACCEPTOR"
write(*, '(1x,a)' ) " 2 ==> LIG2; C60 >> LIG2"
write(*, '(1x,a)' ) " 3 ==> LIG3; CAR >> LIG3"        
print'("")'
write(*, '(a17)', advance = 'no') ' Input = '
read*, input

!   PPH >> ACCEPTOR 
!   C60 >> LIG2
!   CAR >> LIG3

select case( input )

    case( '1' )
        CALL Read_DONOR( time , survival )

    case( '2' )
        CALL Read_LIG1( time , survival )

    case( '3' )
        CALL Read_LIG2( time , survival )

end select

counter = size(survival)

allocate( Corrl ( counter ) )

! calcule correlation energy ...
Corrl  = 0.0d0

do i = 0 , counter-1

    sum1 = sum( survival(1:counter-i) * survival(i+1:counter) ) / ( counter - i )
    sum2 = sum( survival(1:counter-i) ) / ( counter - i )
    sum3 = sum( survival(i+1:counter) ) / ( counter - i )

    Corrl(i+1) = sum1 - sum2 * sum3
    
end do

weight = corrl(1)

! save the results ...
open(unit=20, file="Auto_Correlation_population.dat", status="unknown", action="write")

do i = 1 , counter

    write( unit=20 , fmt=120 ), time(i) , Corrl(i) / weight
    
end do   

! correlation time ( if the behavior is exponetial ) ...
do i = 1 , counter
    if( Corrl(i) < 0.0d0 ) then
        n = i
        exit
    end if
end do

tau = int(0.5d0 + sum(Corrl(1:n)) / Corrl(1))

print*, " Correlation time  ~ ", tau + tau

close(20)

deallocate( survival , Corrl , time )

!------------------------------------------------------------------------------------

120 format(f10.4 , t15 , f10.4)

end subroutine Auto_Correlation_Population
!
!
!
!=================================
subroutine Auto_Correlation_Dipole
!=================================
implicit none
real*8  , allocatable   :: dipole(:) , Corrl(:) , time(:)
real*8                  :: weight , sum1 , sum2 , sum3
integer                 :: i , counter , inputstatus , tau , n
character               :: word

print'("")'
write(*, '(1x,a)' ) "Input file = Dipole.dat"
write(*, '(1x,a)' ) "Output file = Auto_Correlation_Dipole.dat"
print'("")'

! verificate system size ...
open(unit=33, file="Dipole.dat", status="old", action="read")

counter = 0

do
    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    counter = counter + 1
end do

counter = counter - 1

rewind 33

allocate(dipole ( counter ) )
allocate(time   ( counter ) )

! read data file ...
do i = 1 , counter

    read( unit=33 , fmt=110 ) , time(i) , dipole(i)
    
end do

close(33)
 
allocate( Corrl ( counter ) )

! calcule correlation energy ...
Corrl  = 0.0d0

do i = 0 , counter-1

    sum1 = sum( dipole(1:counter-i) * dipole(i+1:counter) ) / ( counter - i )
    sum2 = sum( dipole(1:counter-i) ) / ( counter - i )
    sum3 = sum( dipole(i+1:counter) ) / ( counter - i )

    Corrl(i+1) = sum1 - sum2 * sum3
    
end do

weight = Corrl(1)

! save the results ...
open(unit=20, file="Auto_Correlation_dipole.dat", status="unknown", action="write")

do i = 1 , counter

    write( unit=20 , fmt=120 ), time(i) , Corrl(i) / weight
    
end do   

! correlation time ( if the behavior is exponetial ) ...
do i = 1 , counter
    if( Corrl(i) < 0.0d0 ) then
        n = i
        exit
    end if
end do

tau = int(0.5d0 + sum(Corrl(1:n)) / Corrl(1))

print*, " Correlation time  ~ ", tau + tau

close(20)

deallocate( dipole , Corrl , time )

!------------------------------------------------------------------------------------

100 format(a10)
110 format(f9.4 , t40 , f10.4)
120 format(f10.4 , t15 , f10.4)

end subroutine Auto_Correlation_Dipole
!
!
!
!==============================
subroutine Vibration_Population
!==============================
implicit none
real*8  , allocatable   :: pop_1(:) , pop_2(:) , Corrl(:) , time_1(:) , time(:)
real*8                  :: weight , med_1 , med_2 , med_3 , med_4 , med_5 , suport_1 , suport_2
integer                 :: i , n
character               :: input_1 , input_2

print'("")'
write(*, '(1x,a)' ) "Input file = survival.dat"
write(*, '(1x,a)' ) "Output file = Vibration_population.dat"
print'("")'
write(*, '(1x,a)' ) "Enter with data for correlation of population"
write(*, '(1x,a)' ) " D ==> ACCEPTOR;"
write(*, '(1x,a)' ) " I ==> ION;"
write(*, '(1x,a)' ) " 2 ==> LIG2;"
write(*, '(1x,a)' ) " 3 ==> LIG3;"
print'("")'
write(*, '(a17)', advance = 'no') ' Input 1 = '
read*, input_1
write(*, '(a17)', advance = 'no') ' Input 2 = '
read*, input_2

select case( input_1 )

    case( 'A' )
        CALL Read_DONOR( time , pop_1 )

    case( 'I' )
        CALL Read_ION( time , pop_1 )

    case( '2' )
        CALL Read_LIG1( time , pop_1 )

    case( '3' )
        CALL Read_LIG2( time , pop_1 )

end select

select case( input_2 )

    case( 'A' )
        CALL Read_DONOR( time_1 , pop_2 )

    case( 'I' )
        CALL Read_ION( time_1 , pop_2 )

    case( '2' )
        CALL Read_LIG1( time_1 , pop_2 )

    case( '3' )
        CALL Read_LIG2( time_1 , pop_2 )

end select

deallocate( time_1 )

n = size(pop_1)

allocate( Corrl ( n ) )

! calcule the correlation between f_1 and f_2 ...
do i = 1 , n

    suport_1 = sum( pop_1(1:i) ) / dfloat(i)
    suport_2 = sum( pop_2(1:i) ) / dfloat(i)

    med_1 = sum( (pop_1(1:i)-suport_1)*(pop_1(1:i)-suport_1) * (pop_2(1:i)-suport_2)*(pop_2(1:i)-suport_2) ) / dfloat(i)

    med_2 = sum( (pop_1(1:i)-suport_1) * (pop_1(1:i)-suport_1) ) / dfloat(i)
    med_3 = sum( (pop_2(1:i)-suport_2) * (pop_2(1:i)-suport_2) ) / dfloat(i)

    med_4 = sum( (pop_1(1:i)-suport_1) * (pop_1(1:i)-suport_1)* (pop_1(1:i)-suport_1) * (pop_1(1:i)-suport_1) ) / dfloat(i)
    med_5 = sum( (pop_2(1:i)-suport_2) * (pop_2(1:i)-suport_2)* (pop_2(1:i)-suport_2) * (pop_2(1:i)-suport_2) ) / dfloat(i)
    
    weight = sqrt( med_4 * med_5 )

    Corrl(i) = ( med_1 - med_2 * med_3 ) / weight

end do

! save the results ...
open(unit=113, file="Vibration_population.dat", status="unknown", action="write")

do i = 2 , n

   write( unit=113 , fmt=120 ) time(i) , Corrl(i)
   
end do

close( 113 )

deallocate( pop_1 , pop_2 , Corrl , time )

!---------------------------------------------------------------------

120 format(f9.4 , t19 , f9.4)

end subroutine Vibration_Population
!
!
!
!==========================
subroutine Vibration_Dipole
!==========================
implicit none
real*8  , allocatable   :: dipole_1(:) , dipole_2(:) , Corrl(:) , time_1(:) , time_2(:) , time(:)
real*8                  :: weight , med_1 , med_2 , med_3 , med_4 , med_5 , suport_1 , suport_2
integer                 :: i , size_1 , size_2 , inputstatus , n
character(6)            :: word

print'("")'
write(*, '(1x,a)' ) "Input file = Dipole_1.dat and Dipole_2.dat"
write(*, '(1x,a)' ) "Output file = Vibration_dipole.dat"
print'("")'

! verificate size of system 1 ...
open(unit=33, file="Dipole_1.dat", status="old", action="read")

size_1 = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_1 = size_1 + 1

end do

size_1 = size_1 - 1

rewind 33

! verificate size of system 2 ...
open(unit=35, file="Dipole_2.dat", status="old", action="read")

inputstatus = 0
size_2      = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=35 , fmt=100 , iostat = inputstatus ) word
    size_2 = size_2 + 1

end do

size_2 = size_2 - 1

rewind 35

allocate( dipole_1 ( size_1 ) )
allocate( dipole_2 ( size_2 ) )
allocate( time_1   ( size_1 ) )
allocate( time_2   ( size_2 ) )

! read file 1 ...
do i = 1 , size_1

    read( unit=33 , fmt=330 ) time_1(i) , dipole_1(i)

end do

close(33)

! read file 2 ...
do i = 1 , size_2

    read( unit=35 , fmt=330 ) time_2(i) , dipole_2(i)

end do

close(35)

n = size_1

if( n > size_2 ) n = size_2

allocate( Corrl ( n ) )
allocate( time  ( n ) )

time(:) = time_1(1:n)

deallocate( time_1 , time_2 )

! calcule the correlation between f_1 and f_2 ...
do i = 1 , n

    suport_1 = sum( dipole_1(1:i) ) / dfloat(i)
    suport_2 = sum( dipole_2(1:i) ) / dfloat(i)

    med_1 = sum( (dipole_1(1:i)-suport_1)*(dipole_1(1:i)-suport_1) * (dipole_2(1:i)-suport_2)*(dipole_2(1:i)-suport_2) ) / dfloat(i)

    med_2 = sum( (dipole_1(1:i)-suport_1) * (dipole_1(1:i)-suport_1) ) / dfloat(i)
    med_3 = sum( (dipole_2(1:i)-suport_2) * (dipole_2(1:i)-suport_2) ) / dfloat(i)

    med_4 = sum( (dipole_1(1:i)-suport_1) * (dipole_1(1:i)-suport_1)* (dipole_1(1:i)-suport_1) * (dipole_1(1:i)-suport_1) ) / dfloat(i)
    med_5 = sum( (dipole_2(1:i)-suport_2) * (dipole_2(1:i)-suport_2)* (dipole_2(1:i)-suport_2) * (dipole_2(1:i)-suport_2) ) / dfloat(i)
    
    weight = sqrt( med_4 * med_5 )

    Corrl(i) = ( med_1 - med_2 * med_3 ) / weight

end do

! save the results ...
open(unit=113, file="Vibration_dipole.dat", status="unknown", action="write")

do i = 2 , n

   write( unit=113 , fmt=120 ) time(i) , Corrl(i)

end do

close( 113 )

deallocate( dipole_1 , dipole_2 , Corrl , time )

!---------------------------------------------------------------------

100 format(a6)
120 format(f9.4 , t19 , f9.4)
330 format(f9.4 , t40 , f10.5)

end subroutine Vibration_Dipole
!
!
!
!=========================
subroutine Generate_Random
!=========================
implicit none
real*8  , allocatable   :: survival(:,:)
real*8  , allocatable   :: vec(:) , lim(:)
integer                 :: i , j , size_sys , n_level

write(*, '(1x,a)' ) "Enter with size of system"
print'("")'
write(*, '(a17)', advance = 'no') ' Size system = '
read*, size_sys

write(*, '(1x,a)' ) "Enter with the number of level in system (3 or 4)"
print'("")'
write(*, '(a17)', advance = 'no') ' Number of level = '
read*, n_level

allocate( vec      ( size_sys           ) )
allocate( survival ( size_sys , n_level ) )

CALL random_seed

do i = 1 , size_sys
    CALL random_number( vec(i) )
end do

allocate( lim ( 0 : n_level ) )

do i = 0 , n_level
    lim(i) = dfloat(i) / dfloat(n_level)
end do

survival = 0.0d0
do i = 1 , size_sys

    do j = 1 , n_level
    
        if( vec(i) > lim(j-1) .AND. vec(i) <= lim(j) ) then
            survival(i,j) = 1.0d0
            exit
        end if

    end do

end do

deallocate( vec )

! uotput of the results ...
open(unit=101, file="Random_Survival.dat", status="unknown", action="write")

write( unit=101 , fmt=120 ) "Number of ligands = " , n_level

if( n_level == 3 ) then
    do i = 1 , size_sys

       write( unit=101 , fmt=320 ) survival(i,1) , survival(i,2) , survival(i,3)

    end do
end if

if( n_level == 4 ) then
    do i = 1 , size_sys

       write( unit=101 , fmt=330 ) survival(i,1) , survival(i,2) , survival(i,3) , survival(i,4)

    end do
end if

close( 101 )

deallocate( survival )

!---------------------------------------------------------------------

120 format(a20 , t21 , i1)
320 format(t11 , f10.5 , t21 , f10.5 , t31 , f10.5)
330 format(t11 , f10.5 , t21 , f10.5 , t31 , f10.5 , t41 , f10.5)

end subroutine Generate_Random
!
!
!
!=================================
subroutine Read_DONOR( t , DONOR )
!=================================
implicit none
real*8  , allocatable   , intent(out)   :: t(:)
real*8  , allocatable   , intent(out)   :: DONOR(:)

! local variables ...
character(14)   :: word
integer         :: size_sys , inputstatus , i

! verificate size of system ...
open(unit=33, file="survival.dat", status="old", action="read")

size_sys = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_sys = size_sys + 1

end do

size_sys = size_sys - 2

rewind 33

allocate( t     ( size_sys ) )
allocate( DONOR ( size_sys ) )

! read file ...
read( unit = 33 , fmt = 100 ) word

do i = 1 , size_sys

    read( unit=33 , fmt=130 ) t(i) , DONOR(i)

end do

close(33)

!-------------------------------------------------

100 format(a14)
130 format(f10.5 , t11 , f10.5)

end subroutine Read_DONOR
!
!
!
!=============================
subroutine Read_ION( t , ION )
!=============================
implicit none
real*8  , allocatable   , intent(out)   :: t(:)
real*8  , allocatable   , intent(out)   :: ION(:)

! local variables ...
character(14)   :: word
integer         :: size_sys , inputstatus , i

! verificate size of system ...
open(unit=33, file="survival.dat", status="old", action="read")

size_sys = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_sys = size_sys + 1

end do

size_sys = size_sys - 2

rewind 33

allocate( t   ( size_sys ) )
allocate( ION ( size_sys ) )

! read file ...
read( unit=33 , fmt=100 , iostat = inputstatus ) word

do i = 1 , size_sys

    read( unit=33 , fmt=130 ),  t(i) , ION(i)

end do

close(33)

!-------------------------------------------------

100 format(a14)
130 format(f10.5 , t21 , f10.5)

end subroutine Read_ION
!
!
!
!===============================
subroutine Read_LIG1( t , LIG1 )
!===============================
implicit none
real*8  , allocatable   , intent(out)   :: t(:)
real*8  , allocatable   , intent(out)   :: LIG1(:)

! local variables ...
character(14)   :: word
integer         :: size_sys , inputstatus , i

! verificate size of system ...
open(unit=33, file="survival.dat", status="old", action="read")

size_sys = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_sys = size_sys + 1

end do

size_sys = size_sys - 2

rewind 33

allocate( t    ( size_sys ) )
allocate( LIG1 ( size_sys ) )

! read file ...
read( unit=33 , fmt=100 , iostat = inputstatus ) word

do i = 1 , size_sys

    read( unit=33 , fmt=130 ),  t(i) , LIG1(i)

end do

close(33)

!-------------------------------------------------

100 format(a14)
130 format(f10.5 , t21 , f10.5)

end subroutine Read_LIG1
!
!
!
!===============================
subroutine Read_LIG2( t , LIG2 )
!===============================
implicit none
real*8  , allocatable   , intent(out)   :: t(:)
real*8  , allocatable   , intent(out)   :: LIG2(:)

! local variables ...
character(14)   :: word
integer         :: size_sys , inputstatus , i

! verificate size of system ...
open(unit=33, file="survival.dat", status="old", action="read")

size_sys = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_sys = size_sys + 1

end do

size_sys = size_sys - 2

rewind 33

allocate( t    ( size_sys ) )
allocate( LIG2 ( size_sys ) )

! read file ...
 read( unit=33 , fmt=100 , iostat = inputstatus ) word

do i = 1 , size_sys

    read( unit=33 , fmt=130 ),  t(i) , LIG2(i)

end do

close(33)

!-------------------------------------------------

100 format(a14)
130 format(f10.5 , t31 , f10.5)

end subroutine Read_LIG2
!
!
!
!===============================
subroutine Read_Dyn( Dyn )
!===============================
implicit none
real*8  , allocatable   , intent(out)   :: Dyn(:)

! local variables ...
character(14)   :: word
integer         :: size_sys , inputstatus , i

! verificate size of system ...
open(unit=33, file="dynamic.dat", status="old", action="read")

size_sys = 0
do

    if( inputstatus /= 0 ) exit
    read( unit=33 , fmt=100 , iostat = inputstatus ) word
    size_sys = size_sys + 1

end do

size_sys = size_sys - 1

rewind 33

allocate( Dyn ( size_sys ) )

! read file ...
do i = 1 , size_sys

    read(33,*) Dyn(i)

end do

close(33)

!-------------------------------------------------

100 format(a14)

end subroutine Read_Dyn
!
!
!
end module Correlation_m
