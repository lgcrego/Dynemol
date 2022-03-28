module QCModel_Reciprocal_m

    use f95_precision
    use blas95
    use lapack95
    use omp_lib
    use type_m
    use constants_m
    use Babel_m                 , only : System_Characteristics
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Unit_Cell, Extended_Cell
    use Hamiltonians            , only : X_ij

    public :: EigenSystem_Reciprocal

    private

contains   
!
!
!
!===================================================
 subroutine  EigenSystem_Reciprocal( basis , H , S )
!===================================================
! The routine assumes that:
!   1 - The unit cells are ordered sequentially in the input file ...
!   2 - The atoms of all unit cells are ordered equally in the input file ...
!   3 - The input file has only one structure ==> Unit cell * Number of unit cells ...
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: H(:,:)
real*8          , intent(in)    :: S(:,:)

! local variables ...
complex*16  , allocatable   :: H_rec(:,:) , S_rec(:,:)
complex*16                  :: z
real*8      , allocatable   :: erg(:) , Ki(:,:) , Kf(:,:)
real*8                      :: K(3) , Tn(3) , a1(3) , a2(3) , a3(3) , b1(3) , b2(3) , b3(3) , dir(3)
real*8                      :: a , b , f , ang , step
integer                     :: i , j , l , l0 , p , n , k_ind , N_rec , N_uc , info
character*6                 :: jchar

integer :: ktot = 1000
integer :: pnts = 3

! Primitive vectors of real lattice ...
a1 = [ 2.494d0 , 0.000d0 , 0.000d0 ]
a2 = [ 1.247d0 , 2.160d0 , 0.000d0 ]
a3 = [ 0.000d0 , 0.000d0 , 1.000d0 ]

! Primitive vectors of reciproval lattice ...
CALL compute_bs( a1 , a2 , a3 , b1 , b2 , b3 )

N_rec = count( basis(:) % residue == "000" , 1 )
N_uc  = size(basis) / N_rec
l0    = ( N_uc - 1 ) * N_rec

do j = 1 , N_rec
    write(jchar,'(i6)') 1000+j
    open( unit=1000+j , name="orbital_"//trim(adjustl(jchar))//".dat" , action="write" , status="unknown" ) 
end do

allocate( H_rec ( N_rec , N_rec ) )
allocate( S_rec ( N_rec , N_rec ) )
allocate( erg   ( N_rec )         )

! Initial and final points ...
allocate( Ki ( 3 , pnts ) , source = D_zero )
allocate( Kf ( 3 , pnts ) , source = D_zero )
! First pair ...
Ki(:,1) = 0.5d0 * b1
Kf(:,1) = D_zero
! Second pair ...
Ki(:,2) = D_zero
Kf(:,2) = 2.0d0 * b1 / 3.0d0 + b2 / 3.0d0
!Third pair ...
Ki(:,3) = 2.0d0 * b1 / 3.0d0 + b2 / 3.0d0
Kf(:,3) = 0.5d0 * b1

f = 1.0d0
do p = 1 , pnts

    if( p == pnts ) f = 0.5d0

    ! Direction k ...
    dir  = Kf(:,p) - Ki(:,p)
    a    = dsqrt( dir(1)*dir(1) + dir(2)*dir(2) + dir(3)*dir(3) )
    dir  = dir / a
    step = a / dfloat(ktot)

    do k_ind = 0 , ktot

        K = Ki(:,p) + step * dfloat(k_ind) * dir

        S_rec = C_zero
    
        do j = 1 , N_rec
    
            do i = j , N_rec
    
                do n = 1 , N_uc
    
                    l = ( n - 1 ) * N_rec
    
                    Tn(1) = basis(l+j) % x - basis(l0+i) % x
                    Tn(2) = basis(l+j) % y - basis(l0+i) % y
                    Tn(3) = basis(l+j) % z - basis(l0+i) % z
    
                    ang = dot_product( K , Tn )
    
                    a = dcos( ang )
                    b = dsin( ang )
                    z = a + zi * b
    
                    S_rec( i , j ) = S_rec( i , j ) + z * S( i + l0 , j + l )
    
                end do
    
                ! computing H_rec(i,j) ...
                if( i /= j ) H_rec( i , j ) = X_ij( i , j , basis ) * S_rec( i , j ) 
    
            end do
    
            ! computing H_rec(j,j) ...
            a = X_ij( j + l0 , j + l0 , basis )
            b = X_ij( j + l0 , j , basis )
    
            H_rec( j , j ) = a + b * ( S_rec( j , j ) - D_one )
   
        end do
    
        CALL HEGVD( H_rec , S_rec , erg , 1 , 'N' , 'L' , info )
        if ( info /= 0 ) then
            write(*,*) 'info = ',info,' in HEGVD in EigenSystem '
            stop
        end if
    
        do j = 1 , N_rec
    
            write(1000+j,100) dfloat(p-1) + f * dfloat(k_ind) / dfloat(ktot) , erg(j)
    
        end do

    end do

end do

do j = 1 , N_rec
    close(1000+j)
end do

deallocate( H_rec , S_rec , erg , Ki , Kf )

100 format(3f16.8,ES16.8)
101 format(4a16)

end subroutine  EigenSystem_Reciprocal
!
!
!
!===================================================
subroutine compute_bs( a1 , a2 , a3 , b1 , b2 , b3 )
!===================================================
implicit none
real*8  , intent(in)    :: a1(3) , a2(3) , a3(3)
real*8  , intent(inout) :: b1(3) , b2(3) , b3(3)

! local variables ...
real*8  :: V , f , c(3)

c = cross_product( a2 , a3 )

V = dot_product( a1 , c )

f = two * pi / V

b1 = cross_product( a2 , a3 )
b1 = f * b1

b2 = cross_product( a3 , a1 )
b2 = f * b2

b3 = cross_product( a1 , a2 )
b3 = f * b3

end subroutine compute_bs
!
!
!
!===========================================
 function cross_product( v1 , v2 ) result(v)
!===========================================
implicit none
real*8  , intent(in)    :: v1(3) , v2(3)

! local variables ... 
real*8  :: v(3)

v(1) = v1(2) * v2(3) - v1(3) * v2(2)
v(2) = v1(3) * v2(1) - v1(1) * v2(3)
v(3) = v1(1) * v2(2) - v1(2) * v2(1)

end function cross_product

end module QCModel_Reciprocal_m
