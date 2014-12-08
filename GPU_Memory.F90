!
! Module to manage pinned memory allocation transparently
!
! 
!
#define NR_MEM_PIN_MATRICES 3

module GPU_Memory

    use iso_c_binding, only : C_PTR, c_f_pointer
    use cuda_runtime,  only : cudaMallocHost, cudaFreeHost
    implicit none

    public GPU_Alloc_double, GPU_Dealloc, GPU_Mem_pin_pool_finalize
    private
    
    type matrix_ptr
        real(8), pointer :: matrix(:,:) => null()
    end type matrix_ptr

    logical                :: MemPoolInit = .false.
    integer, parameter     :: N_matrices = NR_MEM_PIN_MATRICES
    integer, save          :: my_matrix = 0
    type(C_PTR), save      :: MemPool_cptr(N_matrices)
    type(matrix_ptr), save :: MemPool(N_matrices)

contains

!-------------------------------------------------------------------
! Initializes the pool, i.e., allocates the memory
subroutine GPU_Mem_pin_pool_init( rows, cols )

    implicit none

    integer :: i, err
    integer, intent(in) :: rows, cols

    do i = 1, N_matrices
        err = cudaMallocHost( MemPool_cptr(i), rows*cols*size_of_double )
        if( err /= 0 ) then
            write(*,*) "Error in cudaMallocHost, GPU_Mem_pin_pool_init", err
            stop
        end if
        call c_f_pointer( MemPool_cptr(i), MemPool(i)%matrix, [rows,cols] )
    end do
    
    MemPoolInit = .true.

end subroutine GPU_Mem_pin_pool_init



!-------------------------------------------------------------------
subroutine GPU_Alloc_double( Matrix_2D , rows, cols )
    implicit none
    real(8), pointer :: Matrix_2D(:,:)
    integer, intent(in) :: rows, cols

    integer, parameter :: size_of_double = sizeof(1.0d0)
    integer :: err

    if(.not. MemPoolInit) call GPU_Mem_pin_pool_init( rows, cols )

    if( my_matrix < N_matrices ) then
        Matrix_2D => MemPool(my_matrix)%matrix
        my_matrix = my_matrix + 1
    else
        allocate( Matrix_2D(rows,cols) )
    end if
end subroutine GPU_Allocate



!-------------------------------------------------------------------
subroutine GPU_Dealloc( A )
    implicit none
    real(8), pointer :: A
    
    integer :: i
    logical :: ass
    ass = .false.

    do i=1,N_matrices
        if( associated( A, target=MemPool(i)%matrix) ) ass = .true.
    end do

    if( ass ) then
        nullify(A)
    else
        deallocate(A)
    end if
end subroutine GPU_Dealloc



!-------------------------------------------------------------------
subroutine GPU_Mem_pin_pool_finalize
    implicit none
    integer :: i, err

    do i=1,N_matrices
        err = cudaFreeHost( MemPool_cptr(i) )
    end do
end subroutine GPU_Mem_pin_pool_finalize

end module GPU_Memory

