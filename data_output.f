 module Data_Output

    use type_m
    use projectors
    use Structure_Builder

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

 contains
!
!
!
!===================================================
 subroutine Dump_Populations(system,basis,bra,ket,t)
!===================================================
 implicit none
 type(structure) , intent(in) :: system
 type(STO_basis) , intent(in) :: basis(:)
 complex*16      , intent(in) :: bra(:,:) , ket(:,:)
 real*8          , intent(in) :: t

! local variables
 real*8                :: pop_1 , pop_2 , pop_3 , pop_4 , pop_5 , pop_6 , pop_total
 character(len=2)      :: fragment

!----------------------------------------------------------c
!                SOME DATA EXTRACTION
!----------------------------------------------------------c

! ==> pop_1 = population of catechol + dye 

 fragment = 'D'

 pop_1 = pop_Slater(basis,bra,ket,fragment)

! ==> pop_2 = population of acceptor

 fragment = 'A'

 pop_2 = pop_Slater(basis,bra,ket,fragment)

! ==> pop_total = total population

 pop_total = pop_Slater(basis,bra,ket)        

 write(26,1001) t , pop_1 , pop_2 , pop_total

1001  format(8(1x,F8.4))
! ---------------------------------------------------- 
     
 include 'formats.h'

 end subroutine Dump_Populations 
!
!
!
!=========================================
 subroutine Dump_DOS( TDOS , PDOS , SPEC ) 
!=========================================
implicit none
type(f_grid)  , intent(in) :: TDOS
type(f_grid)  , intent(in) :: PDOS(:)
type(f_grid)  , intent(in) :: SPEC

! local variables ...
integer         :: i , nr
character(12)   :: string

! . save the TDOS ...
OPEN( unit=3 , file='TDOS.dat' , status='unknown' )
    do i = 1 , size(TDOS%func)
        write(3,10) TDOS%grid(i) , TDOS%average(i)
    end do
CLOSE(3)

! . save the PDOS ...
do nr = 1 , size(trj(1)%list_of_residues)
    string = "PDOS-"//trj(1)%list_of_residues(nr)//".dat"
    OPEN( unit=3 , file=string , status='unknown' )
        do i = 1 , size(PDOS(nr)%func)
            write(3,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i)
        end do
    CLOSE(3)
end do

! . save the peak and broadened specs ...
OPEN( unit=3 , file='spectrum.dat' , status='unknown' )
    do i = 1 , size(SPEC%func)
        write(3,11) SPEC%grid(i) , SPEC%average(i) 
    end do
CLOSE(3)

10   FORMAT(2F12.5)
11   FORMAT(3F13.9)

 end subroutine Dump_DOS
!
!
!
end module Data_Output
