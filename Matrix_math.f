#include "GPU.h"

module Matrix_Math

use constants_m

logical, parameter, public :: return_full = .true.

public Multiply,          &
       syMultiply,        &
       Matrix_Multiply,   &
       Matrix_syMultiply, &
       Matrix_syInvert


interface Multiply
    module procedure Matrix_Multiply_s, Matrix_Multiply
end interface Multiply

interface syMultiply
    module procedure Matrix_syMultiply_s, Matrix_syMultiply
end interface syMultiply

interface bra_x_Op
    module procedure bra_x_Op_, bra_x_Op_alpha
end interface bra_x_Op

interface Op_x_ket
    module procedure Op_x_ket_, Op_x_ket_alpha
end interface Op_x_ket


contains


!------------------------------------------------------------------
! C = alpha*A*B + beta*C
! _s stands for simple
subroutine Matrix_Multiply_s( A, B, C )
    implicit none
    real*8, intent(in)    :: A(:,:), B(:,:)
    real*8, intent(inout) :: C(:,:)
    
    call  Matrix_Multiply( A, B, C, 'N', 'N', d_one, d_zero )
    
end subroutine Matrix_Multiply_s

!------------------------------------------------------------------
! C = alpha*A*B + beta*C
subroutine Matrix_Multiply( A, B, C, transA, transB, alpha, beta )
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
end subroutine Matrix_Multiply


!------------------------------------------------------------------
! C = alpha*A*B + beta*C   - for symmetric matrices
subroutine Matrix_syMultiply_s( side, uplo, A, B, C )
    implicit none
    character, intent(in)  :: side, uplo
    real*8,    intent(in)  :: A(:,:), B(:,:)
    real*8,    intent(out) :: C(:,:)
    
    call  Matrix_syMultiply( side, uplo, A, B, C, d_one, d_zero )
    
end subroutine Matrix_syMultiply_s


!------------------------------------------------------------------
! C = alpha*A*B + beta*C   - for symmetric matrices
subroutine Matrix_syMultiply( side, uplo, A, B, C, alpha, beta )
    implicit none
    character, intent(in)    :: side, uplo
    real*8,    intent(in)    :: A(:,:), B(:,:), alpha, beta
    real*8,    intent(inout) :: C(:,:)

    integer :: m, n, ldA, ldB, ldC

    ldA = size(A,1)
    ldB = size(B,1)
    ldC = size(C,1)
    
    m = ldA
    n = size(C,2)

    call xPU_dsymm( side, uplo, m, n, alpha, A, ldA, B, ldB, beta, C, ldC )
end subroutine Matrix_syMultiply


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
! Invert a symmetric matrix
subroutine Matrix_syInvert( A, UpLo, full )
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

end subroutine Matrix_syInvert

! Don't use gpu cublas for blas2 operations (slower than mkl)
#define xPU_  
!------------------------------------------------------------------
! Performs <res| = <bra|Op, where bra is a vector and Op a matrix
subroutine bra_x_Op_( res, bra, Op )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: bra(:)
    real*8,     intent(in)  :: Op(:,:)

    integer :: n
    
    n = size(Op, 1)
    call xPU_dzgemv( 'T', n, n, c_one, Op, n, bra, i_one, c_zero, res, i_one )
end subroutine bra_x_Op_

function bra_Op_( bra, Op )
    implicit none
    complex*16, intent(in)  :: bra(:)
    real*8,     intent(in)  :: Op(:,:)

    complex*16 :: bra_Op_(size(bra))
    integer    :: n

    n = size(bra)
    call xPU_dzgemv( 'T', n, n, c_one, Op, n, bra, i_one, c_zero, bra_Op_, i_one )
end function bra_Op_

subroutine bra_x_Op_alpha( res, bra, Op, alpha )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: bra(:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: n
    
    n = size(Op, 1)
    call xPU_dzgemv( 'T', n, n, alpha, Op, n, bra, i_one, c_zero, res, i_one )
end subroutine bra_x_Op_alpha

!------------------------------------------------------------------
! Performs |res> = Op|ket>, where ket is a vector and Op a matrix
subroutine Op_x_ket_( res, Op, ket )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: ket(:)
    real*8,     intent(in)  :: Op(:,:)

    integer :: n
    
    n = size(Op, 1)
    call xPU_dzgemv( 'N', n, n, c_one, Op, n, ket, i_one, c_zero, res, i_one )
end subroutine Op_x_ket_

subroutine Op_x_ket_alpha( res, Op, ket, alpha )
    implicit none
    complex*16, intent(out) :: res(:)
    complex*16, intent(in)  :: ket(:), alpha
    real*8,     intent(in)  :: Op(:,:)

    integer :: n
    
    n = size(Op, 1)
    call xPU_dzgemv( 'N', n, n, alpha, Op, n, ket, i_one, c_zero, res, i_one )
end subroutine Op_x_ket_alpha

end module Matrix_Math
