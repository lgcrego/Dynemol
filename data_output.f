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
!---------------------------------------------------
 subroutine Dump_Populations(system,basis,bra,ket,t)
!---------------------------------------------------
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
 end module Data_Output
