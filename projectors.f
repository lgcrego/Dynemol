 module projectors

    use type_m
    use constants_m
    use Structure_Builder
    use Semi_Empirical_Parms

 contains
!
!
!
!-----------------------------------------------
 function pop_Slater(basis,za,zb,fragment)
!-----------------------------------------------
 type(STO_basis)  , intent(in) :: basis(:)
 complex*16       , intent(in) :: za(:,:) , zb(:,:)
 character(len=1) , optional   :: fragment

! local variables
 real*8                                :: pop_Slater
 complex*16 , dimension(size(za(1,:))) :: pop
 integer                               :: n , i 

! . projection of | k(t) >  onto the atom k_atom 

 pop(:) = (0.d0,0.d0)

 if( present(fragment) ) then
    forall(n=1:n_part , i=1:size(basis) , basis(i)%fragment == fragment)  pop(n) = pop(n) + za(i,n)*zb(i,n)
 else
    forall(n=1:n_part , i=1:size(basis))  pop(n) = pop(n) + za(i,n)*zb(i,n)
 end if

 pop_Slater = real(sum(pop))

 end function

 end module projectors
