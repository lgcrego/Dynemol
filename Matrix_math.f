module Matrix_Math

use constants_m

logical, parameter, public :: return_full = .true.

public Matrix_Symmetrize, &
       Multiply,          &
       syMultiply,        &
       syInvert


interface syInvert
    module procedure Matrix_syInvert_0, Matrix_syInvert_1
end interface syInvert

interface Multiply
    module procedure Matrix_Multiply_0, Matrix_Multiply_1
end interface Multiply

interface syMultiply
    module procedure Matrix_syMultiply_0, Matrix_syMultiply_1, Matrix_syMultiply_2
end interface syMultiply

interface bra_x_Op
    module procedure vec_bra_x_Op_, vec_bra_x_Op_alpha, mat_bra_x_Op_, mat_bra_x_Op_alpha
end interface bra_x_Op

interface Op_x_ket
    module procedure vec_Op_x_ket_, vec_Op_x_ket_alpha, mat_Op_x_ket_, mat_Op_x_ket_alpha
end interface Op_x_ket


contains


!------------------------------------------------------------------
! C = A*B + C
! _s stands for simple
subroutine Matrix_Multiply_0( A, B, C )
    implicit none
    real*8, intent(in)    :: A(:,:), B(:,:)
    real*8, intent(inout) :: C(:,:)

    call Matrix_Multiply_1( A, B, C, 'N', 'N', d_one, d_zero )

end subroutine Matrix_Multiply_0

!------------------------------------------------------------------
! C = alpha*A*B + beta*C
subroutine Matrix_Multiply_1( A, B, C, transA, transB, alpha, beta )
    implicit none
    real*8,    intent(in)    :: A(:,:), B(:,:), alpha, beta
    real*8,    intent(inout) :: C(:,:)
    character, intent(in)    :: transA, transB

    integer :: m, n, k, ldA, ldB, ldC

    ldA = size(A,1)
    ldB = size(B,1)
    ldC = ldA

    if( transA == 'N' .or. transA == 'n') then
        m = ldA
        n = size(B,2)
        k = ldB
    else
        m = size(A,2)
        n = ldB
        k = ldA
    end if

    call xPU_dgemm( transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC )
end subroutine Matrix_Multiply_1


!------------------------------------------------------------------
! C = A*B  - for symmetric matrices
subroutine Matrix_syMultiply_0( A, B, C )
    implicit none
    real*8,    intent(in)  :: A(:,:), B(:,:)
    real*8,    intent(out) :: C(:,:)

    call Matrix_syMultiply_2( 'L', 'U', A, B, C, d_one, d_zero )

end subroutine Matrix_syMultiply_0


!------------------------------------------------------------------
! C = A*B  - for symmetric matrices
subroutine Matrix_syMultiply_1( side, uplo, A, B, C )
    implicit none
    character, intent(in)  :: side, uplo
    real*8,    intent(in)  :: A(:,:), B(:,:)
    real*8,    intent(out) :: C(:,:)

    call Matrix_syMultiply_2( side, uplo, A, B, C, d_one, d_zero )

end subroutine Matrix_syMultiply_1


!------------------------------------------------------------------
! C = alpha*A*B + beta*C   - for symmetric matrices
subroutine Matrix_syMultiply_2( side, uplo, A, B, C, alpha, beta )
    implicit none
    character, intent(in)    :: side, uplo
    real*8,    intent(in)    :: A(:,:), B(:,:), alpha, beta
    real*8,    intent(inout) :: C(:,:)

    integer :: m, n, ldA, ldB, ldC

    ldA = size(A,1)
    ldB = size(B,1)
    ldC = size(C,1)

    m = ldC
    n = size(C,2)

    call xPU_dsymm( side, uplo, m, n, alpha, A, ldA, B, ldB, beta, C, ldC )

end subroutine Matrix_syMultiply_2


!------------------------------------------------------------------
! Symmetrize a matrix
subroutine Matrix_Symmetrize( A, UpLo )
    implicit none
    real*8,    intent(inout) :: A(:,:)
    character, intent(in)    :: UpLo

    integer :: n, i, j, ii, jj
    integer, parameter :: block = 115

    n = size(A,1)

    if( UpLo == 'U' .or. UpLo == 'u' ) then    ! Copy the upper part to the lower
        ! Loop trough blocks
        !$omp parallel do schedule(dynamic) private(jj,ii,j,i) default(shared)
        do jj = 1, n, block
            do ii = 1, jj-1, block
                ! Loop inside block:
                do j = jj, min( n, jj+block )
                do i = ii, min( n, ii+block )
                    A(j,i) = A(i,j)
                end do
                end do
            end do
            ! Diagonal block is the last
            do j = jj+1, min( n, jj+block )
            do i = jj, min( n, j )
                A(j,i) = A(i,j)
            end do
            end do
        end do
        !$omp end parallel do
    else                                       ! Copy the lower part to the upper
        ! Loop trough blocks
        !$omp parallel do schedule(dynamic) private(jj,ii,j,i) default(shared)
        do jj = 1, n, block
            ! Diagonal block is the first
            do j = jj+1, min( n, jj+block )
            do i = jj, min( n, j )
                A(j,i) = A(i,j)
            end do
            end do
            do ii = jj+1, n, block
                ! Loop inside block:
                do j = jj, min( n, jj+block )
                do i = ii, min( n, ii+block )
                    A(j,i) = A(i,j)
                end do
                end do
            end do
        end do
        !$omp end parallel do
    end if
end subroutine Matrix_Symmetrize



!------------------------------------------------------------------
! Invert a symmetric matrix, A = A^-1
subroutine Matrix_syInvert_0( A, full )
    implicit none
    real*8,    intent(inout)        :: A(:,:)
    logical,   intent(in), optional :: full

    integer :: info
    logical :: to_symmetrize
    to_symmetrize = .false.

    if( present(full) ) to_symmetrize = full
    call xPU_syInvert( A, 'U', size(A,1), info )
#ifndef USE_GPU
! For now, only dgetrf/i exist for the GPU, so we only need to symmetrize for cpu
    if( to_symmetrize ) call Matrix_Symmetrize( A, 'U' )
#endif
end subroutine Matrix_syInvert_0

!------------------------------------------------------------------
! Invert a symmetric matrix, A = A^-1
subroutine Matrix_syInvert_1( A, UpLo, full )
    implicit none
    real*8,    intent(inout)        :: A(:,:)
    character, intent(in)           :: UpLo
    logical,   intent(in), optional :: full

    integer :: info
    logical :: to_symmetrize
    to_symmetrize = .false.

    if( present(full) ) to_symmetrize = full
    call xPU_syInvert( A, UpLo, size(A,1), info )
#ifndef USE_GPU
    if( to_symmetrize ) call Matrix_Symmetrize( A, UpLo )
#endif
end subroutine Matrix_syInvert_1


! Don't use gpu cublas for single blas2 operations (memory bounded, slower than cpu due to transfer overhead)
#define xPU_dzgemv  dzgemv

!------------------------------------------------------------------
! Performs <res| = <bra|Op, where bra is a *vector* and Op a matrix
subroutine vec_bra_x_Op_( res, bra, Op )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: bra(:)
    real*8,     intent(in)  :: Op(:,:)

    integer :: n

    n = size(Op, 1)
    call xPU_dzgemv( 'T', n, n, c_one, Op, n, bra, i_one, c_zero, res, i_one )
end subroutine vec_bra_x_Op_

! Performs <res| = alpha*<bra|Op, where bra is a *vector* and Op a matrix
subroutine vec_bra_x_Op_alpha( res, bra, Op, alpha )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: bra(:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: n

    n = size(Op, 1)
    call xPU_dzgemv( 'T', n, n, alpha, Op, n, bra, i_one, c_zero, res, i_one )
end subroutine vec_bra_x_Op_alpha

! Performs <res| = <bra|Op, where bra is a *matrix* and Op a matrix
subroutine mat_bra_x_Op_( res, bra, Op )
    implicit none
    complex*16, intent(out) :: res(:,:)
    complex*16, intent(in)  :: bra(:,:)
    real*8,     intent(in)  :: Op(:,:)

    call mat_bra_x_Op_alpha( res, bra, Op, c_one )
end subroutine mat_bra_x_Op_

! Performs <res| = alpha*<bra|Op, where bra is a *matrix* and Op a matrix
subroutine mat_bra_x_Op_alpha( res, bra, Op, alpha )
    implicit none
    complex*16, intent(out) :: res(:,:)
    complex*16, intent(in)  :: bra(:,:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: m, n

    m = size(bra, 1)  ! nr. of elements (orbitals)
    n = size(bra, 2)  ! nr. of wavefunctions (particles, n_part)
    call xPU_dzgemm( 'T', 'N', m, n, m, alpha, Op, m, bra, m, c_zero, res, m )

end subroutine mat_bra_x_Op_alpha


!------------------------------------------------------------------
! Performs |res> = Op|ket>, where ket is a *vector* and Op a matrix
subroutine vec_Op_x_ket_( res, Op, ket )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: ket(:)
    real*8,     intent(in)  :: Op(:,:)

    integer :: n

    n = size(Op, 1)
    call xPU_dzgemv( 'N', n, n, c_one, Op, n, ket, i_one, c_zero, res, i_one )
end subroutine vec_Op_x_ket_

! Performs |res> = alpha*Op|ket>, where ket is a *vector* and Op a matrix
subroutine vec_Op_x_ket_alpha( res, Op, ket, alpha )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: ket(:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: n

    n = size(Op, 1)
    call xPU_dzgemv( 'N', n, n, alpha, Op, n, ket, i_one, c_zero, res, i_one )
end subroutine vec_Op_x_ket_alpha

! Performs |res> = Op|ket>, where ket is a *matrix* and Op a matrix
subroutine mat_Op_x_ket_( res, Op, ket )
    implicit none
    complex*16, intent(out) :: res(:,:)
    complex*16, intent(in)  :: ket(:,:)
    real*8,     intent(in)  :: Op(:,:)

    call mat_Op_x_ket_alpha( res, Op, ket, c_one )

end subroutine mat_Op_x_ket_

! Performs |res> = alpha*Op|ket>, where ket is a *matrix* and Op a matrix (for ElHl_*)
subroutine mat_Op_x_ket_alpha( res, Op, ket, alpha )
    implicit none
    complex*16, intent(out) :: res(:,:)
    complex*16, intent(in)  :: ket(:,:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: m, n

    m = size(ket, 1)  ! nr. of elements (orbitals)
    n = size(ket, 2)  ! nr. of wavefunctions (particles, n_part)
    call xPU_dzgemm( 'N', 'N', m, n, m, alpha, Op, m, ket, m, c_zero, res, m )

end subroutine mat_Op_x_ket_alpha

end module Matrix_Math