#ifdef USE_GPU

#include <cuda.h>
#include "magma_operators.h"

//----------------------------------------------
// Prototypes:
__host__ void
acummulate_vec_ax_async(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecA,
    const cuDoubleComplex * __restrict__ vecsX,
    const cuDoubleComplex * __restrict__ vecY,
    const cudaStream_t stream );

__global__ void
acummulate_vec_ax_kernel(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecA,
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


//----------------------------------------------
// Functions / kernels:

//----------------------------------------------
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
    
    acummulate_vec_ax_kernel <<< Blocks, Threads, 0, stream >>> (n, ld, k, vecA, vecsX, vecY);
}

//- - - - - - - - - - - - - - - - - - - - - - - -
__global__ void
acummulate_vec_ax_kernel(
    const int n, const int ld, const int k,
    const cuDoubleComplex * __restrict__ vecA,
    const cuDoubleComplex * __restrict__ vecsX,
    cuDoubleComplex       * __restrict__ vecY )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n )
    {
        cuDoubleComplex res = make_cuDoubleComplex( 0.0, 0.0 );
        for( int j=0; j<k; ++j )
        {
            res += vecA[j] * vecsX[i];
            vecsX += ld;
        }

        vecY[i] = res;
    }
}


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
        z[i] = xx + a*y[i];
        d[i] = z[i] - xx;
    }
}

#endif