!
! Module to interface the CUDA runtime functions
!

module cuda_runtime
!#include "cuda_runtime.fi"

    integer :: cuda_error

    interface

        integer function cudaMallocHost( ptr, size ) bind(C, name = "cudaMallocHost") 
            use, intrinsic :: iso_c_binding
            implicit none 
            type (C_PTR) :: ptr 
            integer (C_SIZE_T), value :: size
        end function cudaHostAlloc

        integer function cudaFreeHost( ptr ) bind(C, name = "cudaFreeHost") 
            use, intrinsic :: iso_c_binding 
            implicit none 
            type (C_PTR), value :: ptr
        end function cudaFreeHost
 
    end interface

end module cuda_runtime
