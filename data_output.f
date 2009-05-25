 module Data_Output

    use type_m
    use projectors
    use Structure_Builder

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

 contains
!
!
 subroutine Dump_Populations(system,bra,ket,t)

 type(structure) , intent(in) :: system
 complex*16      , intent(in) :: bra(:,:) , ket(:,:)
 real*8          , intent(in) :: t

 integer :: list_of_atoms(system%atoms) , i , n0
 real*8  :: pop_1 , pop_2 , pop_3 , pop_4 , pop_5 , pop_6 , pop_total

!----------------------------------------------------------c
!                SOME DATA EXTRACTION
!----------------------------------------------------------c

! --- preparing the list of atoms for population analysis --- c

 n0 = ((2*nnx+1)*(2*nny+1) - 1) * unit_cell%atoms

! ==> pop_1 = population of catechol + dye 
     
 list_of_atoms(:) = 0

 forall(i=n0+313:system%atoms) list_of_atoms(i) = i

 pop_1 = pop_Slater(bra,ket,list_of_atoms,system)

! ==> pop_2 = population of dye
     
 list_of_atoms(:) = 0

 forall(i=n0+324:system%atoms) list_of_atoms(i) = i

 pop_2 = pop_Slater(bra,ket,list_of_atoms,system)

! ==> pop_3 = population of catechol

 list_of_atoms(:) = 0

 forall(i=n0+313:n0+323) list_of_atoms(i) = i

 pop_3 = pop_Slater(bra,ket,list_of_atoms,system)  

! ==> pop_4 = population of Q(a)+catechol
     
 list_of_atoms(:) = 0

 forall(i=n0+313:n0+335) list_of_atoms(i) = i

 pop_4 = pop_Slater(bra,ket,list_of_atoms,system)

! ==> pop_5 = population of Q(b)
     
 list_of_atoms(:) = 0

 forall(i=n0+336:n0+352) list_of_atoms(i) = i

 pop_5 = pop_Slater(bra,ket,list_of_atoms,system)

! ==> pop_6 = population of Q(c)
     
 list_of_atoms(:) = 0

 forall(i=n0+353:n0+369) list_of_atoms(i) = i

 pop_6 = pop_Slater(bra,ket,list_of_atoms,system)

! ==> pop_total = total population

 list_of_atoms(:) = 0

 forall(i=1:system%atoms) list_of_atoms(i) = i

 pop_total  = pop_Slater(bra,ket,list_of_atoms,system)        


 write(26,1001) t , pop_1 , pop_2 , pop_3 , pop_4 , pop_5 , pop_6 , pop_total

1001  format(8(1x,F8.4))
! ---------------------------------------------------- 
     
 include 'formats.h'

 end subroutine Dump_Populations 
!
!
 end module Data_Output
