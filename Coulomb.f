module Coulomb_m
    use constants_m         , only : a_Bohr
    use util_m              , only : fact
    use Structure_Builder   , only : Extended_Cell              , &
                                     Generate_Structure         , &
                                     Basis_Builder              , &
                                     ExCell_basis


    public  :: Int_Coulomb

contains
!
!
!
!=====================
subroutine Int_Coulomb
!=====================
implicit none
real*8  , allocatable   :: coul(:,:,:,:)
real*8                  :: x1 , x2 , x3 , x4 , rn1 , rn2 , rn3 , rn4 , r_ab
integer                 :: i , j , k , l , n1 , l1 , n2 , l2 , n3 , l3 , n4 , l4
integer , parameter     :: mxl = 3

integer , parameter     :: atoms_Coul(2)=[1,2]  

!-------------------------------------------------
CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell , ExCell_basis )

!print*, ExCell_basis%symbol ; stop

! define n's ...
n1 = Excell_basis( 10 ) % n
n3 = Excell_basis( 15 ) % n
n2 = n1
n4 = n3

print*, n1 , n2 , n3 , n4

! define l's ...
l1 = Excell_basis( 10 ) % l
l3 = Excell_basis( 15 ) % l
l2 = l1
l4 = l3

print*, l1 , l2 , l3 , l4

! define x's ...
x1 = Excell_basis( 10 ) % zeta(1)
x3 = Excell_basis( 15 ) % zeta(1)
x2 = 9.d0
x4 = 1.d0

print*, x1 , x2 , x3 , x4

! define rn's ...
rn1 = dsqrt( (x1+x1)**(n1+n1+1)/fact(n1+n1) )
rn2 = dsqrt( (x2+x2)**(n2+n2+1)/fact(n2+n2) )
rn3 = dsqrt( (x3+x3)**(n3+n3+1)/fact(n3+n3) )
rn4 = dsqrt( (x4+x4)**(n4+n4+1)/fact(n4+n4) )

! calcule r_ab ...
r_ab = dsqrt( (Excell_basis(10)%x - Excell_basis(15)%x)** 2 + &
              (Excell_basis(10)%y - Excell_basis(15)%y)** 2 + &
              (Excell_basis(10)%z - Excell_basis(15)%z)** 2 )

! unit adjustment ...
!r_ab = a_Bohr * r_ab

print*, r_ab

CALL consta

allocate( coul(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl) , source=0.d0 )

! calcula Coulomb potential ...
CALL coul0sim( n1 , l1 , x1 , n2 , l2 , x2 , n3 , l3 , x3 , n4 , l4 , x4 , rn1 , rn2 , rn3 , rn4 , r_ab , coul )

do i = -mxl , mxl
    do j = -mxl , mxl
        do k = -mxl , mxl
            do l = -mxl , mxl
                if (coul(i,j,k,l) /= 0.0d0 ) print*, coul(i,j,k,l) , i , j , k , l
            end do
        end do
    end do
end do

!print*, coul( numbers_m(1) , numbers_m(1) , numbers_m(2) , numbers_m(2) )

deallocate( coul )

end subroutine Int_Coulomb

end module Coulomb_m
