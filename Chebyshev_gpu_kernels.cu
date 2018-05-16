#ifdef USE_GPU

#include <cuda.h>
#include "cublas_v2.h"
#include "magma_operators.h"

// Array element correspondig to matrix element (i,j) in column major order
#define IDX( i, j, LD ) ((i) + (j)*(LD))


//-------------------
// external variables
extern cublasHandle_t myHandle;
extern cudaStream_t   cublas_default, stream[];
extern const int      nStreams;


//----------------------------------------------
// Prototypes:
__host__ void
acummulate_vec_ax_async(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecA,
    const cuDoubleComplex * __restrict__ vecsX,
    const cuDoubleComplex * __restrict__ vecY,
    const cudaStream_t stream );

// __global__ void
// acummulate_vec_ax_kernel(
//     const int n, const int ld, const int k,
//     const cuDoubleComplex * __restrict__ vecA,
//     const cuDoubleComplex * __restrict__ vecsX,
//     cuDoubleComplex       * __restrict__ vecY );

__global__ void
acummulate_vec_ax_kernel(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecsX,
    cuDoubleComplex       * __restrict__ vecY );

__host__ void
Zvec_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex * y,
    cuDoubleComplex       * z,
    const cudaStream_t stream );

__global__ void
Zvec_sub_yinplace_kernel(
    const int n,
    const cuDoubleComplex * __restrict__ x,
    cuDoubleComplex       * __restrict__ y );

__global__ void
Zvec_sub_kernel(
    const int n,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z );

__host__ void
fused_Zxpby_and_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex   a,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    cuDoubleComplex       * d,
    const cudaStream_t stream );

__global__ void
fused_Zxpby_and_subtract_kernel(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex   a,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    cuDoubleComplex       * d );

__host__ void
Zaxpby_async(
    const int n,
    const cuDoubleComplex & a,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex & b,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    const cudaStream_t stream );

__global__ void
Zaxpby_async_kernel(
    const int n,
    const cuDoubleComplex a,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex b,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z);

__host__ void
hadamard_minus(
    const int n,
    const int m,
    const double * const __restrict__ x,
    const double * const __restrict__ y,
    double       * const __restrict__ z,
    const cudaStream_t stream );

__global__ void
hadamard_minus_kernel(
    const int n,
    const double * const __restrict__ x,
    const double * const __restrict__ y,
    double       * const __restrict__ z);

__host__ void
calculate_A(
    const int n,
    const int ld,
    const cuDoubleComplex * const __restrict__ bra,
    const cuDoubleComplex * const __restrict__ ket,
    double                * const __restrict__ A,
    const cudaStream_t stream );

__global__ void
calculate_rho_kernel(
    const int n,
    const int ld,
    const cuDoubleComplex * const __restrict__ bra,
    const cuDoubleComplex * const __restrict__ ket,
    double                * const __restrict__ rho );


//----------------------------------------------
// Functions / kernels:

//----------------------------------------------
// __constant__ version
// MAX_ORDER must have the same value it has in Fortran code
#define MAX_ORDER 25
__constant__ cuDoubleComplex cA[MAX_ORDER];

__host__ void
acummulate_vec_ax_async(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecA,
    const cuDoubleComplex * __restrict__ vecsX,
    cuDoubleComplex       * __restrict__ vecY,
    const cudaStream_t stream )
{  
    const int Threads = 128;                         // Threads per block ## opt. for SM >= 3.0
    const int Blocks = (n + Threads-1) / Threads;    // We need enough blocks to span all the elements
    
    cudaMemcpyToSymbolAsync( cA, vecA, k, 0, cudaMemcpyDeviceToDevice, stream );

    acummulate_vec_ax_kernel <<< Blocks, Threads, 0, stream >>> (n, ld, k, vecsX, vecY);
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
acummulate_vec_ax_kernel(
    const int n, const int ld, const int k,
//     const cuDoubleComplex * __restrict__ vecA,  // using constant memory now
    const cuDoubleComplex * __restrict__ vecsX,
    cuDoubleComplex       * __restrict__ vecY )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if( i < n )
    {
        cuDoubleComplex res = make_cuDoubleComplex( 0.0, 0.0 );
        for( int j=0; j<k; ++j )
        {
            res += cA[j] * vecsX[i];
            vecsX += ld;
        }

        vecY[i] = res;
    }
}

// //----------------------------------------------
// // __shared__ version
// __host__ void
// acummulate_vec_ax_async(
//     const int n, const int ld, const int k,
//     const cuDoubleComplex * __restrict__ vecA,
//     const cuDoubleComplex * __restrict__ vecsX,
//     cuDoubleComplex       * __restrict__ vecY,
//     const cudaStream_t stream )
// {  
//     const int Threads = 128;                         // Threads per block ## opt. for SM >= 3.0
//     const int Blocks = (n + Threads-1) / Threads;    // We need enough blocks to span all the elements
// 
//     acummulate_vec_ax_kernel <<< Blocks, Threads, 0, stream >>> (n, ld, k, vecA, vecsX, vecY);
// }
// 
// //- - - - - - - - - - - - - - - - - - - - - - - -
// __global__ void
// acummulate_vec_ax_kernel(
//     const int n, const int ld, const int k,
//     const cuDoubleComplex * __restrict__ vecA,
//     const cuDoubleComplex * __restrict__ vecsX,
//     cuDoubleComplex       * __restrict__ vecY )
// {
//     const int i = blockIdx.x*blockDim.x + threadIdx.x;
//     __shared__ cuDoubleComplex A[MAX_ORDER];
//     
//     // load vecA into A
//     if ((threadIdx.x == 0) && (i < k))
//         A[i] = vecA[i];
//     __syncthreads();
// 
//     if( i < n )
//     {
//         cuDoubleComplex res = make_cuDoubleComplex( 0.0, 0.0 );
//         for( int j=0; j<k; ++j )
//         {
//             res += A[j] * vecsX[i];
//             vecsX += ld;
//         }
// 
//         vecY[i] = res;
//     }
// }


//----------------------------------------------
__host__ void
Zvec_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex * y,
    cuDoubleComplex       * z,
    const cudaStream_t stream )
{
    const int Threads = 128;                         // Threads per block ## opt. for SM >= 3.0
    const int Blocks = (n + Threads-1) / Threads;    // We need enough blocks to span all the elements

    if( z == y )
        Zvec_sub_yinplace_kernel <<< Blocks, Threads, 0, stream >>> (n, x, z);
//     else if( z == x )
//         // not used
    else
        Zvec_sub_kernel <<< Blocks, Threads, 0, stream >>> (n, x, y, z);
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
Zvec_sub_yinplace_kernel(
    const int n,
    const cuDoubleComplex * __restrict__ x,
    cuDoubleComplex       * __restrict__ y )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n )
        y[i] = x[i] - y[i];
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
Zvec_sub_kernel(
    const int n,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n )
        z[i] = x[i] - y[i];
}


//----------------------------------------------
// z = x + a*y
// d = z - x
__host__ void
fused_Zxpby_and_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex   a,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    cuDoubleComplex       * d,
    const cudaStream_t stream )
{
    const int Threads = 128;                         // Threads per block ## opt. for SM >= 3.0
    const int Blocks = (n + Threads-1) / Threads;    // We need enough blocks to span all the elements

    fused_Zxpby_and_subtract_kernel <<< Blocks, Threads, 0, stream >>> ( n, x, a, y, z, d );
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
fused_Zxpby_and_subtract_kernel(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex   a,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    cuDoubleComplex       * d )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n )
    {
        const cuDoubleComplex xx = x[i];     // just in case that d == x
        z[i] = cuCfma( a, y[i], xx );        // z[i] = xx + a*y[i];
        d[i] = z[i] - xx;
    }
}


//----------------------------------------------
// z = a*x + b*y
__host__ void
Zaxpby_async(
    const int n,
    const cuDoubleComplex &              a,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex &              b,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    const cudaStream_t stream )
{
    const int Threads = 128;                                // Threads per block ## opt. for SM >= 3.0
    const int Blocks = (n + Threads-1) / Threads;           // We need enough blocks to span all the elements

    Zaxpby_async_kernel <<< Blocks, Threads, 0, stream >>> ( n, a, x, b, y, z );
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
Zaxpby_async_kernel(
    const int n,
    const cuDoubleComplex                a,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex                b,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n )  z[i] = a*x[i] + b*y[i];
}






//==============================================
// Kernels for diabatic-Ehrenfest

//----------------------------------------------
__host__ void
hadamard_minus(
    const int n,
    const int m,
    const double * const __restrict__ x,
    const double * const __restrict__ y,
    double       * const __restrict__ z,
    const cudaStream_t stream )
{
    const int N = n*m;                               // for the sake of hadamard prodruct, pretend matrices are big vectors
    const int threads = 128;                         // Threads per block ## opt. for SM >= 3.0
    const int blocks = (N + threads-1) / threads;    // We need enough blocks to span all the elements

    hadamard_minus_kernel <<< blocks, threads, 0, stream >>> (N, x, y, z);
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
hadamard_minus_kernel(
    const int n,
    const double * const __restrict__ x,
    const double * const __restrict__ y,
    double       * const __restrict__ z )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i < n)
    {
//         const double zz = z[i];
        z[i] = x[i] * y[i] - z[i];
    }
}


//- - - - - - - - - - - - - - - - - - - - - - - -
// ρ = Re{ ket(j,1)*bra(i,1) -  ket(j,2)*bra(i,2) }
__global__ void
calculate_rho_kernel(
    const int n,
    const int ld,
    const cuDoubleComplex * const __restrict__ bra,
    const cuDoubleComplex * const __restrict__ ket,
    double                * const __restrict__ rho )
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y);
    
    __shared__ cuDoubleComplex ket_j[2];
    
    if ((i < n) && (j < n))
    {
        // first thread loads ket elements into shared memory
        // all threads within the block access the same values
        if (threadIdx.x == 0)
        {
            ket_j[0] = ket[ IDX(j, 0, ld) ];
            ket_j[1] = ket[ IDX(j, 1, ld) ];
        }
        __syncthreads();   // wait for ket to be loaded
        
        rho[ IDX(i, j, ld) ] = real( ket_j[0] * bra[ IDX(i, 0, ld) ] )
                             - real( ket_j[1] * bra[ IDX(i, 1, ld) ] );
    }
    /* fortran code:
    do j = 1, N
        ket_j(:) = AO_ket(j,:)
        do i = 1, N
            rho_eh(i,j) = real( ket_j(1)*AO_bra(i,1) ) - real( ket_j(2)*AO_bra(i,2) )
        end do
    end do
    */
}


//----------------------------------------------
__host__ void
calculate_A(
    const int n,
    const int ld,
    const cuDoubleComplex * const __restrict__ bra,
    const cuDoubleComplex * const __restrict__ ket,
    double                * const __restrict__ A,
    const cudaStream_t stream )
{
    // 1) calculate  ρ = Re{ ket(j,1)*bra(i,1) -  ket(j,2)*bra(i,2) }
    // result is stored in A
    // each block works on the same column
    dim3 threads(128, 1);
    dim3 blocks( (n  + threads.x - 1)/threads.x,
                 (ld + threads.y - 1)/threads.y ); 

    calculate_rho_kernel <<< blocks, threads, 0, stream >>> (n, ld, bra, ket, A);
    
    // 2) A = (ρ + ρ^T) / 2
    double alpha = 0.5;
    cublasSetStream( myHandle, stream );
    cublasDgeam(myHandle, CUBLAS_OP_N, CUBLAS_OP_T, n, n, &alpha, A, ld, &alpha, A, ld, A, ld);
}


#endif
