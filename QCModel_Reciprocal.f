module QCModel_Reciprocal_m

    use f95_precision
    use blas95
    use lapack95
    use omp_lib
    use type_m
    use constants_m
    use parameters_m            , only : NUc , Nat_uc , vecK , Rmod
    use Babel_m                 , only : System_Characteristics
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Extended_Cell

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
real*8      , allocatable   :: erg(:,:)
complex*16                  :: z
real*8                      :: rij(3)
real*8                      :: a , b , c , t , kmod , Nuc_inv
integer                     :: i , j , k_ind , N_rec , n , m , i1 , i2 , lm , ll , info

! local parameters ...
integer :: k_tot  = 10000   ! Discretization of k ...

N_rec = 0
do i = 1 , Nat_uc
    N_rec = N_rec + atom( extended_cell%AtNo(i) ) % DOS
end do

allocate( H_rec (     N_rec , N_rec ) )
allocate( S_rec (     N_rec , N_rec ) )
allocate( erg   ( 0 : k_tot , N_rec ) , source = D_zero )

c = two * pi / Rmod
!c = c / dfloat( k_tot )
c = c / 5000.0d0

Nuc_inv = D_one / dfloat( Nuc )

do k_ind = 0 , k_tot

    kmod = c * dfloat( k_ind )

    H_rec = C_zero
    S_rec = C_zero

    do j = 1 , N_rec

        do i = j , N_rec

            do n = 0 , Nuc - 1

                do m = 0 , Nuc - 1

                    i1 = i + m * N_rec
                    i2 = j + n * N_rec

                    ll = merge( i1 , i2 , i1 .le. i2 )
                    lm = merge( i2 , i1 , i1 .lt. i2 )

                    rij(1) = basis( i1 ) % x - basis( i2 ) % x
                    rij(2) = basis( i1 ) % y - basis( i2 ) % y
                    rij(3) = basis( i1 ) % z - basis( i2 ) % z

                    t = vecK(1) * rij(1) + vecK(2) * rij(2) + vecK(3) * rij(3)
                    t = t * kmod

                    a = dcos(t)
                    b = dsin(t)

                    z = a + zi * b

                    H_rec( i , j ) = H_rec( i , j ) + z * H( lm , ll )
                    S_rec( i , j ) = S_rec( i , j ) + z * S( lm , ll )

                end do

            end do

            rij(1) = basis( j ) % x - basis( i ) % x
            rij(2) = basis( j ) % y - basis( i ) % y
            rij(3) = basis( j ) % z - basis( i ) % z

            t = vecK(1) * rij(1) + vecK(2) * rij(2) + vecK(3) * rij(3)
            t = t * kmod

            a = dcos(t) * Nuc_inv
            b = dsin(t) * Nuc_inv

            z = a + zi * b

            H_rec( i , j ) = H_rec( i , j ) * z
            S_rec( i , j ) = S_rec( i , j ) * z

        end do

    end do

    CALL HEGVD( H_rec , S_rec , erg( k_ind , : ) , 1 , 'N' , 'L' , info )
    if ( info /= 0 ) then
        write(*,*) 'info = ',info,' in HEGVD in EigenSystem '
        write(*,*) "k_ind = " , k_ind
        write(*,*) "N = " , N_rec
        stop
    end if

end do

do j = 1 , N_rec
    do k_ind = 0 , k_tot

        kmod = c * dfloat( k_ind )
        write(1000+j,*) kmod , erg(k_ind,j)

    end do
end do

deallocate( H_rec , S_rec , erg )

end subroutine  EigenSystem_Reciprocal

end module QCModel_Reciprocal_m
