 module projectors

    use type_m
    use constants_m
    use Structure_Builder
    use Semi_Empirical_Parms

 contains
!
!
!
!-----------------------------------------
 function pop_Slater(basis,za,zb,fragment)
!-----------------------------------------
 type(STO_basis)  , intent(in) :: basis(:)
 complex*16       , intent(in) :: za(:) , zb(:)
 character(len=1) , optional   :: fragment

! local variables
 real*8       :: pop_Slater
 complex*16   :: pop 
 integer      :: n , i 

! projection of | k(t) >  onto the atom k_atom ...

 pop = C_zero

 if( present(fragment) ) then
    
    pop = sum( za*zb , DIM = 1 , mask = basis%fragment == fragment )

 else

    pop = sum( za(:) * zb(:) )

 end if

 pop_Slater = pop

 end function

 end module projectors
