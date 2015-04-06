/*
    -- KBLAS (version 1.0) --
       Ahmad Abdelfattah, Center of Extreme Computing
	   Hatem Ltaief, Supercomputing Laboratory
	   David Keyes, Center of Extreme Computing
	   King Abdullah University of Science and Technology (KAUST)
       June 2013
	   KBLAS is a subset of BLAS routines highly optimized for NVIDIA GPUs 
*/
/**
	-- Center of Extreme Computing and Supercomputing Laboratory
	-- Division of Applied Mathematics and Computational Science
	-- King Abdullah University of Science and Technology
	-- (C) Copyright 2013

	Redistribution  and  use  in  source and binary forms, with or without
	modification,  are  permitted  provided  that the following conditions
	are met:

	*	Redistributions  of  source  code  must  retain  the above copyright
		notice,  this  list  of  conditions  and  the  following  disclaimer.
	* 	Redistributions  in  binary  form must reproduce the above copyright
		notice,  this list of conditions and the following disclaimer in the
		documentation  and/or other materials provided with the distribution.
	* 	Neither  the  name of the University of Tennessee, Knoxville nor the
		names of its contributors may be used to endorse or promote products
		derived from this software without specific prior written permission.

	THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	''AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
	LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	
**/

/**
	Original sources merged and modified by Alberto Torres to have:
	 * mixed type functions like MKL dzgemv
	 * specialized cases: inc = 1
	 * out-of-place dzgemv
**/

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
//#include "operators.h"
#include "magma_operators.h"
#include <stdio.h>

//-------------------------------------
// From operators.h:

// __device__ static __inline__ float              conj_if(int _if_, float x){return x;}
__device__ static __inline__ double             conj_if(int _if_, double x){return x;}
// __device__ static __inline__ cuFloatComplex     conj_if(int _if_, cuFloatComplex x){if(_if_==0)return x; else return cuConjf(x);}
// __device__ static __inline__ cuDoubleComplex    conj_if(int _if_, cuDoubleComplex x){if(_if_==0)return x; else return cuConj(x);}

/*************************************************************/
/**
 *   Atomic add on double precision, as suggested by the CUDA programming Guide
 **/
__device__ static __inline__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

/**
 *   Atomic add for double complex ('Z' precision)
 **/
__device__ static __inline__ void atomicAdd(cuDoubleComplex* address, cuDoubleComplex val)
{
    atomicAdd( (double*) (&(*address).x) ,val.x);
    atomicAdd( (double*) (&(*address).y) ,val.y);
}

//-------------------------------------
// From gemv2_core.cuh:

// specialized for:
//  * incx == incy := 1
//  * rows == cols := n
template <class T1, class T2, int nb, int tcol, int ept, int width, int ept_>
__global__ void
gemvn(int n, T1 alpha, T1 *A, int lda, T2 *x, T1  beta, T2 *y, int mod_r, int mod_c, int threshold)
{
    const int   tx   = threadIdx.x ;
    const int   ty   = threadIdx.y ;
    const int   blkc = blockIdx.x ;
    const int   by  =   blockIdx.y; 
    
    T2 res_1_   = MAGMA_Z_ZERO; //make_zero<T2>();
    T1 areg[ept];
    T1 breg[ept];
    
    __shared__ T2 la[nb * tcol];
    
    if(blkc == gridDim.x-1)
    {
        if(mod_r > 0) {if(tx >= mod_r) return;}
    }
    
    // number of full blocks to process
    int count = (n/width)/gridDim.y + (by < (n/width)%gridDim.y);
    
    {
        int start = by * ((n/width)/gridDim.y) + min(by, (n/width)%gridDim.y);
        
        // Advance 'A'
        A += nb * blkc;
        A += start * width * lda;
        
        // Advance 'x'
        x += start * width; 
        
        // Advance 'y'
        y += (blkc * nb);
    }
    
    if(by != gridDim.y-1){if(count == 0) return;}
    else {if(count == 0 && mod_c == 0) return;}
    
    const int j = ty * ept * lda + tx;
    
    if(count >= 2)
    {
        // read 1st block
        #pragma unroll
        for(int k = 0; k < ept; k++)
            areg[k] = A[j + k * lda];
        A += width * lda; 
    }
    
    int Vblocks = 0;
    #pragma unroll
    for(Vblocks = 0; Vblocks < (count/2)*2; Vblocks+=2)
    {
        // read 2nd block
        #pragma unroll
        for(int k = 0; k < ept; k++)
            breg[k] = A[j + k * lda];
        A += width * lda;
        
        // compute 1st
        #pragma unroll
        for(int k = 0; k < ept; k++)
            res_1_ += areg[k] * x[(ty * ept + k)]; 
        x += width; 
        
        // prefetch 1st block
        if(Vblocks != ((count/2)*2-2) )
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                areg[k] = A[j + k * lda];
            A += width * lda;
        }
        
        // compute 2nd
        #pragma unroll
        for(int k = 0; k < ept; k++)
            res_1_ += breg[k] * x[(ty * ept + k)]; 
        x += width;
    }
    
    if(count%2 >= 1)
    {
        //if(ty == 0 && tx == 0)printf("hi \n");
        // read the remaining block
        #pragma unroll
        for(int k = 0; k < ept; k++)
            areg[k] = A[j + k * lda];
        A += width * lda;
        
        // process remaining block
        #pragma unroll
        for(int k = 0; k < ept; k++)
            res_1_ += areg[k] * x[(ty * ept + k)];
        x += width;
    }
    
    //if(ty == 0 && tx == 0)printf("(%d, %d): by = %d\n ", tx, ty, by);
    //if(ty == 0 && tx == 0)printf("(%d, %d): mod_c = %d\n", tx, ty, mod_c);
    
    if(by == gridDim.y-1)
    {
        #pragma unroll
        for(int k = 0; k < ept; k++) {breg[k] = 0.0;} //make_zero<T1>();}
        
        //if(ty == 0 && tx == 0)printf("mod_c = %d\n", mod_c);
        if(mod_c != 0)
        {
            if(ty < threshold)
            {
                #pragma unroll
                for(int k = 0; k < ept; k++)
                    breg[k] = A[j + k * lda];
            }
            else if(ty == threshold)
            {
                #pragma unroll
                for(int k = 0; k < ept_; k++)
                    breg[k] = A[j + k * lda];
            }
            
            // compute
            if(ty < threshold)
            {
                #pragma unroll
                for(int k = 0; k < ept; k++)
                    res_1_ += breg[k] * x[(ty * ept + k)];
            }
            else if (ty == threshold)
            {
                #pragma unroll
                for(int k = 0; k < ept_; k++)
                    res_1_ += breg[k] * x[(ty * ept + k)];
            }
            //if(ty == 0 && tx == 0)printf("hi 2\n");
            //int l;
            //#pragma unroll
            //for(l = 0; l < (mod_c/tcol); l++)
            //  breg[l] = A[l*tcol*lda + ty*lda + tx];
            //if(ty < (mod_c%tcol) )
            //  breg[l] = A[l*tcol*lda + ty*lda + tx];
            
            //#pragma unroll
            //for(l = 0; l < (mod_c/tcol); l++)
            //  res_1_ += breg[l] * x[(l*tcol + ty) * incx];
            //if(ty < (mod_c%tcol) )
            //  res_1_ += breg[l] * x[(l*tcol + ty) * incx];
        }
    }
    
    la[ty * nb + tx] = res_1_;
    __syncthreads();
    
    if(ty == 0)
    {   
        res_1_ = MAGMA_Z_ZERO; //make_zero<T2>();
        #pragma unroll
        for(int k = 0; k < tcol; k++)
            res_1_ += la[k * nb + tx];
        // use atomics
        atomicAdd(&y[tx], (alpha*res_1_));
        //y[tx] = alpha * res_1_ ;
    }
}

// specialized for:
//  * incx == incy := 1
//  * rows == cols := n
template <class T1, class T2, int nb, int tcol, int ept, int width, int ept_>
__global__ void
gemvt(int n, T1 alpha, T1 *A, int lda, T2 *x, T1  beta, T2 *y, int mod_r, int mod_c, int threshold, int conj)
{
    const int   tx   = threadIdx.x ;
    const int   ty   = threadIdx.y ;
    const int   blkc = blockIdx.x ;
    const int   by  =   blockIdx.y; 
    
    T2 res[ept] = {0.0}; //make_zero<T2>()};
    T1 areg[ept];
    T1 breg[ept];
    
    __shared__ T2 la[nb * width];
    
    if(blkc == gridDim.x-1 && mod_c != 0){if(ty > threshold) return;}
    
    // number of full blocks to process
    int count = (n/nb)/gridDim.y + (by < (n/nb)%gridDim.y);
    
    {
        int start = by * ((n/nb)/gridDim.y) + min(by, (n/nb)%gridDim.y);
        
        // Advance 'A'
        A += blkc * width * lda;
        A += start * nb;
        
        // Advance 'x'
        x += start * nb; 
        
        // Advance 'y'
        y += (blkc * width);
    }
    
    if(by != gridDim.y-1){if(count == 0) return;}
    
    const int j = ty * ept * lda + tx;
    
    const int irregular = ( (mod_c != 0) && (blkc == gridDim.x-1) && (ty == threshold) );
    
    if(count >= 2)
    {
        //if(blkc == 0 && by == 0 && tx == 0 && ty == 0)printf("hi-1\n");
        // read 1st block
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                areg[k] = A[j + k * lda];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                areg[k] = A[j + k * lda];
        }
        A += nb; 
    }
    
    int Vblocks = 0;
    #pragma unroll
    for(Vblocks = 0; Vblocks < (count/2)*2; Vblocks+=2)
    {
        //if(blkc == 0 && by == 0 && tx == 0 && ty == 0)printf("hi-2\n");
        // read 2nd block
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                breg[k] = A[j + k * lda];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                breg[k] = A[j + k * lda];
        }
        A += nb;
        
        // compute 1st
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                res[k] += conj_if(conj, areg[k]) * x[tx];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                res[k] += conj_if(conj, areg[k]) * x[tx]; 
        }
        x += nb; 
        
        // prefetch 1st block
        if(Vblocks != ((count/2)*2-2) )
        {
            if(irregular)
            {
                #pragma unroll
                for(int k = 0; k < ept_; k++)
                    areg[k] = A[j + k * lda];
            }
            else
            {
                #pragma unroll
                for(int k = 0; k < ept; k++)
                    areg[k] = A[j + k * lda];
            }
            A += nb;
        }
        
        // compute 2nd
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                res[k] += conj_if(conj, breg[k]) * x[tx];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                res[k] += conj_if(conj, breg[k]) * x[tx]; 
        }
        x += nb; 
    }
    
    if(count%2 >= 1)
    {
        //printf("hi from (%d, %d) \n", tx, ty);
        
        // read the remaining block
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                areg[k] = A[j + k * lda];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                areg[k] = A[j + k * lda];
        }
        A += nb;
        
        // process remaining block
        if(irregular)
        {
            #pragma unroll
            for(int k = 0; k < ept_; k++)
                res[k] += conj_if(conj, areg[k]) * x[tx];
        }
        else
        {
            #pragma unroll
            for(int k = 0; k < ept; k++)
                res[k] += conj_if(conj, areg[k]) * x[tx];
        }
        x += nb;
    }
    
    if(by == gridDim.y-1)
    {
        #pragma unroll
        for(int k = 0; k < ept; k++){breg[k] = 0.0;} //make_zero<T1>();}
        
        //if(ty == 0 && tx == 0)printf("mod_c = %d\n", mod_c);
        
        //if(blkc == 0 && by == 0 && tx == 0 && ty == 0)printf("hi-4\n");
        if(tx < mod_r)
        {
            if(irregular)
            {
                #pragma unroll
                for(int k = 0; k < ept_; k++)
                    breg[k] = A[j + k * lda];
            }
            else
            {
                #pragma unroll
                for(int k = 0; k < ept; k++)
                    breg[k] = A[j + k * lda];
            }
            
            // compute
            if(irregular)
            {
                #pragma unroll
                for(int k = 0; k < ept_; k++)
                    res[k] += conj_if(conj, breg[k]) * x[tx];
            }
            else
            {
                #pragma unroll
                for(int k = 0; k < ept; k++)
                    res[k] += conj_if(conj, breg[k]) * x[tx];
            }
        }
        
    }
    
    #pragma unroll
    for(int k = 0; k < ept; k++)
        la[(ty*ept + k)*nb + tx] = res[k];
    __syncthreads();
    
    if(ty == 0 && tx < width)
    {   
        T2 res_1_ = MAGMA_Z_ZERO; //make_zero<T2>();
        #pragma unroll
        for(int k = tx; k < (tx+nb); k++)
            res_1_ += la[tx * nb + k%nb];
        
        // use atomics
        if(mod_c != 0)
        {
            if(blkc == gridDim.x-1) {if(tx < mod_c)atomicAdd(&y[tx], (alpha*res_1_));}
            else atomicAdd(&y[tx], (alpha*res_1_));
        }
        else {atomicAdd(&y[tx], (alpha*res_1_));}
    }
}


//-------------------------------------
// From scal_core.cuh:

template <class T1, class T2>
__global__ void
scal(int n, T1 alpha, T2 *x)
{
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    
    const int gtx = bx * blockDim.x + tx;
    
    if(gtx < n) x[gtx] *= alpha;
}

// specialized for:
//  * incx == incy := 1
template <class T1, class T2>
__global__ void
scal(int n, T1 alpha, T2 * __restrict__ x, T2 * __restrict__ y)
{
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    
    const int gtx = bx * blockDim.x + tx;
    
    if(gtx < n) y[gtx] = alpha * x[gtx];
}

//-------------------------------------
// From zscal.cu:

#define dzscal_nbx      (128)

// int kblas_dzscal_driver(int n, double alpha, cuDoubleComplex *x, int incx, cudaStream_t stream)
// {
//     int gridx = n / dzscal_nbx + (n % dzscal_nbx != 0);
//     
//     dim3 dimBlock(dzscal_nbx, 1);
//     dim3 dimGrid(gridx, 1);
//     
//     scal<double, cuDoubleComplex> <<<dimGrid, dimBlock, 0, stream>>>(n, alpha, x, incx);
//     
//     return 0;
// }

int kblas_dzscal_driver(int n, double alpha, cuDoubleComplex *x, cudaStream_t stream)
{
    int gridx = n / dzscal_nbx + (n % dzscal_nbx != 0);
    
    dim3 dimBlock(dzscal_nbx, 1);
    dim3 dimGrid(gridx, 1);
    
    scal<double, cuDoubleComplex> <<<dimGrid, dimBlock, 0, stream>>>(n, alpha, x);
    
    return 0;
}

int kblas_dzscal_driver(int n, double alpha, cuDoubleComplex *x, cuDoubleComplex *y, cudaStream_t stream)
{
    int gridx = n / dzscal_nbx + (n % dzscal_nbx != 0);
    
    dim3 dimBlock(dzscal_nbx, 1);
    dim3 dimGrid(gridx, 1);
    
    scal<double, cuDoubleComplex> <<<dimGrid, dimBlock, 0, stream>>>(n, alpha, x, y);
    
    return 0;
}

int kblas_dzscal_async(int n, double alpha, cuDoubleComplex *x, cudaStream_t stream)
{
    return kblas_dzscal_driver(n, alpha, x, stream);
}

int kblas_dzscal_async(int n, double alpha, cuDoubleComplex *x, cuDoubleComplex *y, cudaStream_t stream)
{
    return kblas_dzscal_driver(n, alpha, x, y, stream);
}

//-------------------------------------
// From zgemv2.cu:

#if(SM >= 30)

#define dzgemvn_nb       (32)
#define dzgemvn_ntcol    (4)
#define dzgemvn_ept      (2)
#define dzgemvn_width    (dzgemvn_ntcol*dzgemvn_ept)
#define dzgemvn_by       (4)

#define dzgemvt_nb       (32)
#define dzgemvt_ntcol    (4)
#define dzgemvt_ept      (2)
#define dzgemvt_width    (dzgemvt_ntcol*dzgemvt_ept)
#define dzgemvt_by       (4)

#else

#define dzgemvn_nb               (64)
#define dzgemvn_ntcol    		(8)
#define dzgemvn_ept              (2)
#define dzgemvn_width    (dzgemvn_ntcol*dzgemvn_ept)
#define dzgemvn_by               (1)

#define dzgemvt_nb               (64)
#define dzgemvt_ntcol    		(8)
#define dzgemvt_ept              (2)
#define dzgemvt_width    (dzgemvt_ntcol*dzgemvt_ept)
#define dzgemvt_by               (1)
#endif


  
int kblas_dzgemv2_driver(	char trans, int n,
						double alpha, double *dA, int lda, 
						cuDoubleComplex *dX,
						double  beta, cuDoubleComplex *dY,
						cudaStream_t stream)
{	
	if(trans == 'n' || trans == 'N')
	{
		// scaling with beta
		kblas_dzscal_async(n, beta, dY, stream);
		
		int mod_r = n % dzgemvn_nb;
		int mod_c = n % dzgemvn_width;	
		
		int blocks = n/dzgemvn_nb;
		if(mod_r != 0) blocks += 1;
		
		const int thread_x = dzgemvn_nb;
		const int thread_y = dzgemvn_ntcol; 
		const int ept = dzgemvn_ept;
		
		int threshold = mod_c / ept; 
		int ept_ = mod_c % ept;
		dim3 dimBlock(thread_x, thread_y);
		dim3 dimGrid(blocks, dzgemvn_by);
		switch(ept_)
		{
			case 0: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 0><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 1: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 1><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 2: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 2><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 3: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 3><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 4: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 4><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 5: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 5><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 6: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 6><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 7: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 7><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
            case 8: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 8><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold); break;
			default: printf("irregular part %d is not supported, please extend the case statement of dzgemv\n", ept_); exit(1);
		}
	}	// end of non-transpose case
	else if(trans == 't' || trans == 'T' || trans == 'c' || trans == 'C')
	{
		// scaling with beta
		kblas_dzscal_async(n, beta, dY, stream);
		
		int mod_r = n % dzgemvt_nb;
		int mod_c = n % dzgemvt_width;
		
		int blocks = n/dzgemvt_width;
		if(mod_c != 0) blocks += 1;
		
		const int thread_x = dzgemvt_nb;
		const int thread_y = dzgemvt_ntcol;
		const int ept = dzgemvt_ept;
		
		int threshold = mod_c / ept;
		int ept_ = mod_c % ept;
		
		dim3 dimBlock(thread_x, thread_y);
		dim3 dimGrid(blocks, dzgemvt_by);
		
		int conj;
		if(trans == 'c' || trans == 'C')conj = 1;
		else conj = 0;
		//printf("modr = %d, modc = %d, threshold = %d, ept_ = %d \n", mod_r, mod_c, threshold, ept_);
		switch(ept_)
		{
            case 0: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 0><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 1: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 1><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 2: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 2><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 3: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 3><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 4: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 4><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 5: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 5><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 6: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 6><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 7: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 7><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
            case 8: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 8><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dY, mod_r, mod_c, threshold, conj); break;
			default: printf("irregular part %d is not supported, please extend the case statement of dzgemv\n", ept_); exit(1);
		}
	}
	else
	{	
		printf("ZGEMV error: Unrecognized transpose mode %c \n", trans);
		return -1;
	}
	
	return 0;
}


// specialized for:
//  * incx == incy := 1
//  * rows == cols := n
int kblas_dzgemv2_driver(   char trans, int n,
                            double alpha, double *dA, int lda, 
                            cuDoubleComplex *dX,
                            double beta, cuDoubleComplex *dY,
                            cuDoubleComplex *dZ, cudaStream_t stream)
{   
    if(trans == 'n' || trans == 'N')
    {
        // scaling with beta
        kblas_dzscal_async(n, beta, dY, dZ, stream);
        
        int mod_r = n % dzgemvn_nb;
        int mod_c = n % dzgemvn_width;   
        
        int blocks = n/dzgemvn_nb;
        if(mod_r != 0) blocks += 1;
        
        const int thread_x = dzgemvn_nb;
        const int thread_y = dzgemvn_ntcol; 
        const int ept = dzgemvn_ept;
        
        int threshold = mod_c / ept; 
        int ept_ = mod_c % ept;
        dim3 dimBlock(thread_x, thread_y);
        dim3 dimGrid(blocks, dzgemvn_by);
        switch(ept_)
        {
            case 0: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 0><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 1: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 1><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 2: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 2><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 3: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 3><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 4: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 4><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 5: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 5><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 6: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 6><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 7: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 7><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            case 8: gemvn<double, cuDoubleComplex, dzgemvn_nb, dzgemvn_ntcol, ept, dzgemvn_width, 8><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold); break;
            default: printf("irregular part %d is not supported, please extend the case statement of dzgemv\n", ept_); exit(1);
        }
    }   // end of non-transpose case
    else if(trans == 't' || trans == 'T' || trans == 'c' || trans == 'C')
    {
        // scaling with beta
        kblas_dzscal_async(n, beta, dY, dZ, stream);
        
        int mod_r = n % dzgemvt_nb;
        int mod_c = n % dzgemvt_width;
        
        int blocks = n/dzgemvt_width;
        if(mod_c != 0) blocks += 1;
        
        const int thread_x = dzgemvt_nb;
        const int thread_y = dzgemvt_ntcol;
        const int ept = dzgemvt_ept;
        
        int threshold = mod_c / ept;
        int ept_ = mod_c % ept;
        
        dim3 dimBlock(thread_x, thread_y);
        dim3 dimGrid(blocks, dzgemvt_by);
        
        int conj;
        if(trans == 'c' || trans == 'C') conj = 1;
        else conj = 0;
        //printf("modr = %d, modc = %d, threshold = %d, ept_ = %d \n", mod_r, mod_c, threshold, ept_);
        switch(ept_)
        {
            case 0: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 0><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 1: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 1><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 2: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 2><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 3: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 3><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 4: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 4><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 5: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 5><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 6: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 6><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 7: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 7><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            case 8: gemvt<double, cuDoubleComplex, dzgemvt_nb, dzgemvt_ntcol, ept, dzgemvt_width, 8><<<dimGrid, dimBlock, 0, stream>>>(n, alpha, dA, lda, dX, beta, dZ, mod_r, mod_c, threshold, conj); break;
            default: printf("irregular part %d is not supported, please extend the case statement of dzgemv\n", ept_); exit(1);
        }
    }
    else
    {   
        printf("ZGEMV error: Unrecognized transpose mode %c \n", trans);
        return -1;
    }
    
    return 0;
}


extern "C"
int kblas_dzgemv2_async(	char trans, int n,
						double alpha, double *dA, int lda, 
						cuDoubleComplex *dX,
						double  beta, cuDoubleComplex *dY,
						cudaStream_t stream)
{
	return kblas_dzgemv2_driver(	trans, n, alpha, dA, lda, dX, beta, dY, stream);
}

extern "C"
int kblas_dzgemv2_oop_async(char trans, int n,
                            double alpha, double *dA, int lda, 
                            cuDoubleComplex *dX,
                            double  beta, cuDoubleComplex *dY,
                            cuDoubleComplex *dZ, cudaStream_t stream)
{
    return kblas_dzgemv2_driver( trans, n, alpha, dA, lda, dX, beta, dY, dZ, stream);
}
