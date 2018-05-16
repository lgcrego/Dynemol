#ifdef USE_GPU

#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cublas_v2.h"
#include "magma_operators.h"


#ifdef __INTEL_COMPILER
 #define __restrict__ restrict
 #include <mathimf.h>
#endif


#define SAFE(S)  if(cudaSuccess != S) \
                    printf("ERROR(%i): %s\n%s\n\n", __LINE__, cudaGetErrorName(S), cudaGetErrorString(S));

#ifdef GPU_DEBUG
 #define DEBUG(A) A
#else
 #define DEBUG(A) /*Nothing*/
#endif

#ifdef GPU_TIMING
 #include "magma_auxiliary.h"
 #define time_init()  double __my_time = magma_wtime()
 #define time_end(S)  __my_time = magma_wtime()-__my_time; printf("time %s: %.4f\n", S,__my_time); fflush(stdout)
#else
 #define time_init()  /* No timing */
 #define time_end(S)  /* No timing */
#endif

// Array element correspondig to matrix element (i,j) in column major order
#define IDX( i, j, LD ) ((i) + (j)*(LD))

// Short names to copy functions
#define cudaDcopy( DST, SRC, N, DIRECTION ) cudaMemcpy( (void*) DST, (const void*) SRC, (N)*sizeof(double), DIRECTION );
#define cudaDcopy_D2H( DST, SRC, N )        cudaDcopy( DST, SRC, N, cudaMemcpyDeviceToHost );
#define cudaDcopy_H2D( DST, SRC, N )        cudaDcopy( DST, SRC, N, cudaMemcpyHostToDevice );
#define cudaDcopy_D2D( DST, SRC, N )        cudaDcopy( DST, SRC, N, cudaMemcpyDeviceToDevice );

#define cudaDcopyAsync( DST, SRC, N, DIRECTION, STREAM ) cudaMemcpyAsync( (void*) DST, (const void*) SRC, (N)*sizeof(double), DIRECTION, STREAM );
#define cudaDcopyAsync_D2H( DST, SRC, N, STREAM )        cudaDcopyAsync( DST, SRC, N, cudaMemcpyDeviceToHost, STREAM );
#define cudaDcopyAsync_H2D( DST, SRC, N, STREAM )        cudaDcopyAsync( DST, SRC, N, cudaMemcpyHostToDevice, STREAM );
#define cudaDcopyAsync_D2D( DST, SRC, N, STREAM )        cudaDcopyAsync( DST, SRC, N, cudaMemcpyDeviceToDevice, STREAM );

#define cudaZcopy( DST, SRC, N, DIRECTION ) cudaMemcpy( (void*) DST, (const void*) SRC, (N)*sizeof(cuDoubleComplex), DIRECTION );
#define cudaZcopy_D2H( DST, SRC, N )        cudaZcopy( DST, SRC, N, cudaMemcpyDeviceToHost );
#define cudaZcopy_H2D( DST, SRC, N )        cudaZcopy( DST, SRC, N, cudaMemcpyHostToDevice );
#define cudaZcopy_D2D( DST, SRC, N )        cudaZcopy( DST, SRC, N, cudaMemcpyDeviceToDevice );

#define cudaZcopyAsync( DST, SRC, N, DIRECTION, STREAM ) cudaMemcpyAsync( (void*) DST, (const void*) SRC, (N)*sizeof(cuDoubleComplex), DIRECTION, STREAM );
#define cudaZcopyAsync_D2H( DST, SRC, N, STREAM )        cudaZcopyAsync( DST, SRC, N, cudaMemcpyDeviceToHost, STREAM );
#define cudaZcopyAsync_H2D( DST, SRC, N, STREAM )        cudaZcopyAsync( DST, SRC, N, cudaMemcpyHostToDevice, STREAM );
#define cudaZcopyAsync_D2D( DST, SRC, N, STREAM )        cudaZcopyAsync( DST, SRC, N, cudaMemcpyDeviceToDevice, STREAM );  

//-------------------------------------------------------------------
// Useful macros to economize code and improve readability:

// set proper stream
#define setStream( S )           cublasSetStream( myHandle, stream[S ## _stream] )
enum { bra_stream=0, ket_stream=1 };

// copy bras/kets in the appropriate stream
#define braCopy( SRC, DST )      cudaZcopyAsync_D2D( DST, SRC, n, stream[bra_stream] )
#define ketCopy( SRC, DST )      cudaZcopyAsync_D2D( DST, SRC, n, stream[ket_stream] )

// vec3 = alpha*A^op*vec1 + beta*vec2
#define Mat_op_vec5( A, OP, ALPHA, VEC1, BETA, VEC2, VEC3, STREAM ) \
                                 kblas_dzgemv2_oop_async( OP, n, ALPHA, A, ld, VEC1, BETA, VEC2, VEC3, STREAM )
// vec1 = alpha*A^op*vec2 + beta*vec1
#define Mat_op_vec4( A, OP, ALPHA, VEC2, BETA, VEC1, STREAM ) \
                                 kblas_dzgemv2_async( OP, n, ALPHA, A, ld, VEC2, BETA, VEC1, STREAM )

// bra1 = H^T*bra2
// ket1 = H*ket2
#define bra_H( VEC1, VEC2, H )   Mat_op_vec4( H, 't', d_one, VEC2, d_zero, VEC1, stream[bra_stream] )
#define H_ket( VEC1, H, VEC2 )   Mat_op_vec4( H, 'n', d_one, VEC2, d_zero, VEC1, stream[ket_stream] )


// find the (complex) element with maximum real or imaginary parts in VEC
#define findMax( VEC, MAX_PTR ) \
{ \
    cublasIdamax( myHandle, 2*n, (const double *)VEC, 1, i_max); \
    *i_max = (*i_max-1)/2;     /* From Fortran-Real to C-Complex index conversion */ \
    cudaZcopy_D2H( MAX_PTR, &VEC[*i_max], 1 ) \
}

// reset cuda event: destroy and create
#define cudaEventReset( EVENT ) \
    cudaEventDestroy(EVENT);  cudaEventCreateWithFlags( & EVENT, cudaEventDisableTiming );


// structure to hold allocated memory, used in convergence_gpu()
struct bras_kets_ptrs
{
    cuDoubleComplex *bras,     *kets,
                    *old_bra, *old_ket,
                    *new_bra, *new_ket,
                    *diff_bra, *diff_ket;
    cuDoubleComplex *pin_mem_chunk;
};

//-------------------
// external variables
extern cublasHandle_t myHandle;
extern cudaStream_t   cublas_default, stream[];
extern const int      nStreams;

//-------------------
// global constants
const int max_order = 25;
const int i_zero = 0;
//...................
const double h_bar     = 6.58264e-4;    // eV.ps
const double tolerance = 1.0e-8;
const double norm_tolerance = 1.0e-8;
const double d_neg_one = -1.0;
const double d_zero    =  0.0;
const double d_one     =  1.0;
const double d_two     =  2.0;
//...................
const cuDoubleComplex z_zero     = make_cuDoubleComplex( 0.0, 0.0 );
const cuDoubleComplex z_one      = make_cuDoubleComplex( 1.0, 0.0 );
const cuDoubleComplex z_neg_one  = make_cuDoubleComplex(-1.0, 0.0 );
const cuDoubleComplex z_imag     = make_cuDoubleComplex( 0.0, 1.0 );
const cuDoubleComplex z_neg_imag = make_cuDoubleComplex( 0.0,-1.0 );


//-------------------
// external functions
extern void Zvec_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex * y,
    cuDoubleComplex       * z,
    const cudaStream_t stream );

extern void acummulate_vec_ax_async(
    const int n, const int ld, const int k,
    const cuDoubleComplex * const __restrict__ vecA,
    const cuDoubleComplex * const __restrict__ vecsX,
    cuDoubleComplex       * const __restrict__ vecY,
    const cudaStream_t stream );

extern void fused_Zxpby_and_subtract(
    const int n,
    const cuDoubleComplex * x,
    const cuDoubleComplex   a,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    cuDoubleComplex       * d,
    const cudaStream_t stream );

extern void Zaxpby_async(
    const int n,
    const cuDoubleComplex &              a,
    const cuDoubleComplex * __restrict__ x,
    const cuDoubleComplex &              b,
    const cuDoubleComplex * __restrict__ y,
    cuDoubleComplex       * __restrict__ z,
    const cudaStream_t stream );

extern void gpu_dgeInvert(
    double *dA, 
    const int n, 
    const int lddA, 
    cudaStream_t stream);

extern "C" int kblas_dzgemv2_async(
    char trans, int n,
    double alpha, const double *dA, int lda, 
    cuDoubleComplex *dX,
    double  beta, cuDoubleComplex *dY,
    cudaStream_t stream);

extern "C" int kblas_dzgemv2_oop_async(
    char trans, int n,
    double alpha, const double *dA, int lda, 
    cuDoubleComplex *dX,
    double beta, cuDoubleComplex *dY,
    cuDoubleComplex *dZ, cudaStream_t stream);

extern void hadamard_minus(
    const int n,
    const int m,
    const double * const x,
    const double * const y,
    double       * const z,
    const cudaStream_t stream );

extern void calculate_A(
    const int n,
    const int ld,
    const cuDoubleComplex * const __restrict__ bra,
    const cuDoubleComplex * const __restrict__ ket,
    double                * const __restrict__ A,
    const cudaStream_t stream );


//-------------------
// Prototypes

//...................
// exported:
extern "C"
void propagation_gpucaller_(
    const int       * const __restrict__ n,
    double          * const __restrict__ tau,
    double          * const __restrict__ save_tau,
    const double    * const __restrict__ t_init,
    const double    * const __restrict__ t_max,
    cuDoubleComplex * const __restrict__ h_PSI_bra,  // in host
    cuDoubleComplex * const __restrict__ h_PSI_ket,  // in host
    const double    * const __restrict__ h_H);

extern "C"
void propagationelhl_gpucaller_(
    const int       * const __restrict__ N,
    const double    * const __restrict__ h_S,        // in host
    const double    * const __restrict__ h_h,        // in host
    double          * const __restrict__ h_H,        // in host
    cuDoubleComplex * const __restrict__ h_AO_bra,   // in host
    cuDoubleComplex * const __restrict__ h_AO_ket,   // in host
    cuDoubleComplex * const __restrict__ h_PSI_bra,  // in host
    cuDoubleComplex * const __restrict__ h_PSI_ket,  // in host
    const double    * const __restrict__ t_init,
    const double    * const __restrict__ t_max,
    double          * const __restrict__ tau,
    double          * const __restrict__ save_tau );

extern "C"
void ehrenfestkernel_gpu_(
    const int       * const __restrict__ N,
    const double    * const __restrict__ h_H,        // in host
    const double    * const __restrict__ h_A,        // in host
    const double    * const __restrict__ h_X,        // in host
    double          * const __restrict__ h_K );      // in host

extern "C" 
void ehrenfestkernel2_gpu_(
    const int             * const __restrict__ N,
    const cuDoubleComplex * const __restrict__ h_bra,      // in host
    const cuDoubleComplex * const __restrict__ h_ket,      // in host
    const double          * const __restrict__ h_H,        // in host
    const double          * const __restrict__ h_X,        // in host
    double                * const __restrict__ h_K );      // in host

//...................
// internal:
void coefficient( const double tau, const int k_max, cuDoubleComplex * coeff );

void convergence_gpu_alloc( const int ld, struct bras_kets_ptrs * const mem );
void convergence_gpu_dealloc( struct bras_kets_ptrs * const mem );
void convergence_gpu(
    const int n, const int ld,
    const struct bras_kets_ptrs * const __restrict__ mem,
    cuDoubleComplex * const __restrict__ PSI_bra,  // in device
    cuDoubleComplex * const __restrict__ PSI_ket,  // in device
    cuDoubleComplex * const __restrict__ coeff,
    int             * const __restrict__ k_ref,
    const double    tau,
    const double    * const __restrict__ H,        // in device
    const double    norm_ref,
    bool            * const __restrict__ ok );

void chebyshev_gpu(
    const int n, const int ld, double tau,
    double          * const __restrict__ save_tau,
    const double t_max,
    const double t_init,
    cuDoubleComplex * const __restrict__ PSI_bra,  // in device
    cuDoubleComplex * const __restrict__ PSI_ket,  // in device
    const double    * const __restrict__ H);       // in device



//-------------------
// Functions:

//-------------------
// utilities:
template <typename T>
void swap( T & __restrict__ a, T & __restrict__ b )
{
    volatile T tmp = a;
    a = b;
    b = tmp;
};


//-------------------------------------------------------------------
void propagation_gpucaller_(
    const int       * const __restrict__ n,
    double          * const __restrict__ tau,
    double          * const __restrict__ save_tau,
    const double    * const __restrict__ t_init,
    const double    * const __restrict__ t_max,
    cuDoubleComplex * const __restrict__ h_PSI_bra,  // in host
    cuDoubleComplex * const __restrict__ h_PSI_ket,  // in host
    const double    * const __restrict__ h_H)        // in host
{
    time_init();

    const int m = *n;
    const int ld = ((m + 31)/32)*32;  // making leading dimension multiple of 32 for memory coalesced accesses
    double * d_H;
    cuDoubleComplex *d_PSI_bra, *d_PSI_ket;

    cudaMalloc((void **) &d_H, ld*m*sizeof(double));
    cudaMalloc((void **) &d_PSI_bra, ld*sizeof(cuDoubleComplex));
    cudaMalloc((void **) &d_PSI_ket, ld*sizeof(cuDoubleComplex));

    cublasSetMatrix( m, m, sizeof(double), (const void *) h_H, m, (void *) d_H, ld);
    cudaZcopy_H2D( d_PSI_bra, h_PSI_bra, m );
    cudaZcopy_H2D( d_PSI_ket, h_PSI_ket, m );

    chebyshev_gpu( m, ld, *tau, save_tau, *t_max, *t_init, d_PSI_bra, d_PSI_ket, d_H );

    cudaZcopy_D2H( h_PSI_bra, d_PSI_bra, m );
    cudaZcopy_D2H( h_PSI_ket, d_PSI_ket, m );

    cudaFree( d_H );
    cudaFree( d_PSI_bra );
    cudaFree( d_PSI_ket );
    
    time_end("chebyshev_gpucaller");
}

                 
//-------------------------------------------------------------------
void chebyshev_gpu(
    const int n,
    const int ld,
    double tau,
    double          * const __restrict__ save_tau,
    const double t_max,
    const double t_init,
    cuDoubleComplex * const __restrict__ PSI_bra,  // in device
    cuDoubleComplex * const __restrict__ PSI_ket,  // in device
    const double    * const __restrict__ H)        // in device 
{
    time_init();
    DEBUG( printf("chebyshev_gpu: start\n"); fflush(stdout); );

    cudaError_t stat;

    cuDoubleComplex *bras, *kets, *new_bra, *new_ket, *coeff, *d_coeff;
    stat = cudaMalloc( (void **) &bras,   max_order*ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &kets,   max_order*ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &new_bra,          ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &new_ket,          ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &d_coeff,   max_order*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMallocHost( (void **) &coeff, max_order*sizeof(cuDoubleComplex) );  SAFE(stat);
    
    cuDoubleComplex *bra[max_order], *ket[max_order];
    for(int i=0; i<max_order; ++i)
    {
        bra[i] = &bras[i*ld];
        ket[i] = &kets[i*ld];
    }
    
    cudaEvent_t coeffs_copied;
    cudaEventCreateWithFlags( &coeffs_copied, cudaEventDisableTiming );

    // reference norm
    cuDoubleComplex z_norm;
    cublasSetStream( myHandle, stream[2] );
    cublasZdotc( myHandle, n, PSI_bra, 1, PSI_ket, 1, &z_norm );    // implicit synchronization here
    const double norm_ref = cuCabs( z_norm );

    int k_ref = 0;
    bool ok = false;

    struct bras_kets_ptrs mem;
    convergence_gpu_alloc( ld, &mem );
    while( true )
    {
        convergence_gpu( n, ld, &mem, PSI_bra, PSI_ket, coeff, &k_ref, tau, H, norm_ref, &ok );
        if(ok) break;
        tau *= 0.9;
    }
    convergence_gpu_dealloc( &mem );

    *save_tau = tau;
    double t = t_init + tau*h_bar;

    if( t_max - t < tau*h_bar )
    {
        tau = (t_max - t)/h_bar;
        coefficient( tau, max_order, coeff );
    }

    bool must_copy_coeff = true;

    while (t < t_max)
    {
        if( must_copy_coeff )
        {
            cudaZcopyAsync_H2D( d_coeff, coeff, max_order, stream[2] );   // copy coeff to the GPU
            cudaEventRecord( coeffs_copied, stream[2]);
            must_copy_coeff = false;
        }
        
        // Ѱ₀ = c₀|Ѱ⟩ = |Ѱ⟩
        braCopy( PSI_bra, bra[0] );
        ketCopy( PSI_ket, ket[0] );

        for( int k=1; k<k_ref; ++k )
        {
            // Ѱₖ = H Ѱₖ₋₁
            bra_H( bra[k], bra[k-1], H );
            H_ket( ket[k], H, ket[k-1] );
        }
        
        // Ѱⁿᵉʷ = Σ cₖѰₖ 
        setStream( bra );
        cudaStreamWaitEvent( stream[bra_stream], coeffs_copied, 0 );
        cublasZgemv( myHandle, CUBLAS_OP_N, n, k_ref, &z_one, bras, ld, d_coeff, 1, &z_zero, new_bra, 1 );
        setStream( ket );
        cudaStreamWaitEvent( stream[bra_stream], coeffs_copied, 0 );
        cublasZgemv( myHandle, CUBLAS_OP_N, n, k_ref, &z_one, kets, ld, d_coeff, 1, &z_zero, new_ket, 1 );

        // check charge conservation
        // norm = |⟨Ѱⁿᵉʷ|Ѱⁿᵉʷ⟩|²
        setStream( bra );                                               // queue after new_bra
        cudaStreamSynchronize( stream[ket_stream] );                    // wait for new_ket
        cublasZdotc( myHandle, n, new_bra, 1, new_ket, 1, &z_norm );    // implicit synchronization here
        const double norm = cuCabs(z_norm);
        
        if( fabs(norm - norm_ref) < norm_tolerance )
        {
            // Ѱ = Ѱⁿᵉʷ
            braCopy( new_bra, PSI_bra );
            ketCopy( new_ket, PSI_ket );
        }
        else
        {
            convergence_gpu_alloc( ld, &mem );
            do
            {
                tau *= 0.975;
                printf("rescaling tau to %g\n", tau);  fflush(stdout);
                convergence_gpu( n, ld, &mem, PSI_bra, PSI_ket, coeff, &k_ref, tau, H, norm_ref, &ok );
            }
            while( !ok );
            convergence_gpu_dealloc( &mem );
            must_copy_coeff = true;
            cudaEventReset( coeffs_copied );
        }

        t += tau*h_bar;

        if( t_max - t < tau*h_bar )
        {
            tau = (t_max - t)/h_bar;
            coefficient( tau, max_order, coeff );
            must_copy_coeff = true;
            cudaEventReset( coeffs_copied );
        }
    } // while (*t < t_max)

    // clean up
    cublasSetStream( myHandle, cublas_default );
    cudaDeviceSynchronize();   // wait for all calculations to complete (is ths really needed?)

    cudaFree( bras );
    cudaFree( kets );
    cudaFree( new_bra );
    cudaFree( new_ket );
    cudaFree( d_coeff );
    cudaFreeHost( coeff );
    
    cudaEventDestroy( coeffs_copied );

    DEBUG( printf("chebyshev_gpu: exit\n"); fflush(stdout); );
    time_end("chebyshev_gpu");
}


//-------------------------------------------------------------------
void convergence_gpu_alloc( const int ld, struct bras_kets_ptrs * const mem )
{
    cuDoubleComplex *dev_mem_chunk;
    cudaMalloc( (void **) &dev_mem_chunk, (2*max_order + 4)*ld*sizeof(cuDoubleComplex) );

    int offset = 0;
    mem->bras     = &dev_mem_chunk[offset]; offset += max_order*ld;
    mem->kets     = &dev_mem_chunk[offset]; offset += max_order*ld;
    mem->old_bra = &dev_mem_chunk[offset]; offset += ld;
    mem->old_ket = &dev_mem_chunk[offset]; offset += ld;
    mem->new_bra = &dev_mem_chunk[offset]; offset += ld;
    mem->new_ket = &dev_mem_chunk[offset];

    cudaMallocHost((void**) &mem->pin_mem_chunk, sizeof(cuDoubleComplex) + sizeof(int));
}


//-------------------------------------------------------------------
void convergence_gpu_dealloc( struct bras_kets_ptrs * const mem )
{
    cudaFree( mem->bras );
    cudaFreeHost( mem->pin_mem_chunk );
}



//-------------------------------------------------------------------
void convergence_gpu(
    const int n, const int ld,
    const struct bras_kets_ptrs * const __restrict__ mem,
    cuDoubleComplex * const __restrict__ PSI_bra,  // in device
    cuDoubleComplex * const __restrict__ PSI_ket,  // in device
    cuDoubleComplex * const __restrict__ coeff,
    int             * const __restrict__ k_ref,
    const double    tau,
    const double    * const __restrict__ H,        // in device
    const double    norm_ref,
    bool            * const __restrict__ ok )
{
    cuDoubleComplex *bras, *kets, *old_bra, *old_ket, *new_bra, *new_ket, *max;

    bras    = mem->bras;       kets     = mem->kets;
    old_bra = mem->old_bra;    old_ket = mem->old_ket;
    new_bra = mem->new_bra;    new_ket = mem->new_ket;
    max = mem->pin_mem_chunk;
    
#define delta_bra old_bra
#define delta_ket old_ket

    int * const i_max = reinterpret_cast<int*>(mem->pin_mem_chunk + 1);

    cuDoubleComplex *bra[max_order], *ket[max_order];
    for(int i=0; i<max_order; ++i)
    {
        bra[i] = &bras[i*ld];
        ket[i] = &kets[i*ld];
    }

    // Ѱ₀ = c₀|Ѱ⟩ = |Ѱ⟩
    braCopy( PSI_bra, bra[0] );
    ketCopy( PSI_ket, ket[0] );
    
    // Ѱᵒˡᵈ = Ѱ₀
    braCopy( PSI_bra, old_bra );
    ketCopy( PSI_ket, old_ket );

    // get coefficients
    coefficient( tau, max_order, coeff );
    
    // establish number of terms in the series (k_max+1)
    int k_max = max_order;
    for(int k=1; k<max_order; ++k)
    {
        if( cuCabs(coeff[k]) < 1.0e-16 )   // originally: coeff[k]/coeff[0] , but coeff[0] == 1
        {
            k_max = k;
            break;
        }
    }
    *k_ref = k_max;

    *ok = false;
    for( int k=1; k<k_max; ++k )
    {
        // Ѱₖ = H Ѱₖ₋₁
        // r = coeff[k]/coeff[k-1]
        bra_H( bra[k], bra[k-1], H );
        H_ket( ket[k], H, ket[k-1] );

        // Ѱⁿᵉʷ = Ѱᵒˡᵈ + cₖѰₖ      &&      δѰ = Ѱᵒˡᵈ - Ѱⁿᵉʷ   (old_bra <-- δѰ via #define)
        fused_Zxpby_and_subtract( n, old_bra, coeff[k], bra[k], new_bra, delta_bra, stream[bra_stream] );
        fused_Zxpby_and_subtract( n, old_ket, coeff[k], ket[k], new_ket, delta_ket, stream[ket_stream] );

        // find maximum element of 〈δѰ|
        setStream( bra );
        findMax( delta_bra, max );      // implicit synchronization here

        if( cuCabs(*max) < tolerance )  // if bra is converged
        {
            // find maximum element of |δѰ⟩
            setStream( ket );
            findMax( delta_ket, max );      // implicit synchronization here

            if( cuCabs(*max) < tolerance )  // if ket is converged
            {
                // check charge conservation:
                //  norm = |⟨Ѱⁿᵉʷ|Ѱⁿᵉʷ⟩|²
                cuDoubleComplex z_norm;
                cublasZdotc( myHandle, n, new_bra, 1, new_ket, 1, &z_norm );    // implicit synchronization here
                const double norm = cuCabs(z_norm);

                if( fabs(norm - norm_ref) < norm_tolerance ) // if charge is conserved
                {
                    // copy answer back and exit:  Ѱ = Ѱⁿᵉʷ
                    braCopy( new_bra, PSI_bra );
                    ketCopy( new_ket, PSI_ket );
                    *ok = true;
                    break;
                }
            } // ket conv.
        } // bra conv.

        // exchange Ѱⁿᵉʷ <-> Ѱᵒˡᵈ
        swap<cuDoubleComplex*>( old_bra, new_bra );
        swap<cuDoubleComplex*>( old_ket, new_ket );
        
    } // for( k=2; k<max_order; ++k )

    // clean up
    
    cublasSetStream( myHandle, cublas_default );
    cudaDeviceSynchronize();
    
#undef delta_bra
#undef delta_ket
}


//-------------------------------------------------------------------
void coefficient( const double tau, const int k_max, cuDoubleComplex * const coeff )
{
    coeff[0] = z_one;

    for( int k=1; k<k_max; ++k )
        coeff[k] = -z_imag * coeff[k-1] * (tau/k) ;
}



//-------------------------------------------------------------------
void propagationelhl_gpucaller_(
    const int       * const __restrict__ N,
    const double    * const __restrict__ h_S,        // in host
    const double    * const __restrict__ h_h,        // in host
    double          * const __restrict__ h_H,        // in host
    cuDoubleComplex * const __restrict__ h_AO_bra,   // in host
    cuDoubleComplex * const __restrict__ h_AO_ket,   // in host
    cuDoubleComplex * const __restrict__ h_PSI_bra,  // in host
    cuDoubleComplex * const __restrict__ h_PSI_ket,  // in host
    const double    * const __restrict__ t_init,
    const double    * const __restrict__ t_max,
    double          * const __restrict__ tau,
    double          * const __restrict__ save_tau )
{
    time_init();

    static bool first_call = true;
    cudaError_t stat;

    const int n  = *N;
    const int ld = ((n + 31)/32)*32;  // making leading dimension multiple of 32

    cudaEvent_t S_inverted, brakets_copied, H_done;
    cudaEventCreateWithFlags( &S_inverted,     cudaEventDisableTiming );
    cudaEventCreateWithFlags( &brakets_copied, cudaEventDisableTiming );
    cudaEventCreateWithFlags( &H_done,         cudaEventDisableTiming );

    static double *H = nullptr;
    static double *S = nullptr;
    static double *h = nullptr;
    static cuDoubleComplex *bra = nullptr, *ket = nullptr;

    if (first_call)
    {
        stat = cudaMalloc((void **) &bra, ld*sizeof(cuDoubleComplex));  SAFE(stat);
        stat = cudaMalloc((void **) &ket, ld*sizeof(cuDoubleComplex));  SAFE(stat);
        stat = cudaMalloc((void **) &S, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &h, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &H, n*ld*sizeof(double));           SAFE(stat);
    }

    // Copy S and h to the GPU
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_S, n, (void *)S, ld, stream[0]);
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_h, n, (void *)h, ld, stream[1]);

    // Copy bra/ket to the GPU
    // waiting for this copy to finish below
    cudaZcopyAsync_H2D( bra, h_PSI_bra, n, stream[2] );
    cudaZcopyAsync_H2D( ket, h_PSI_ket, n, stream[2] );
    cudaEventRecord( brakets_copied, stream[2]);

//    cudaDeviceSynchronize();

    // Sinv = S^(-1)
    double *Sinv = S;
    gpu_dgeInvert( Sinv, n, ld, stream[0] );  // implicit synchronization here
    cudaEventRecord( S_inverted, stream[0]);  // is this really needed? magma trf/tri are syncronous

    // H = Sinv*h
    cublasSetStream( myHandle, stream[1] );            // queue after h copy
    cudaStreamWaitEvent( stream[1], S_inverted, 0 );   // wait in stream[1] for Sinv to be done
    cublasDsymm( myHandle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, n, n, &d_one, Sinv, ld, h, ld, &d_zero, H, ld);
    cudaEventRecord( H_done, stream[1]);

    // copy H back to the CPU
    cudaStreamWaitEvent( stream[2], H_done, 0 );   // wait in stream[2] for H to be done
    cublasGetMatrixAsync( n, n, sizeof(double), (void *)H, ld, (void *)h_H, n, stream[2]);

    // wait for things to complete before calling chebyshev_gpu()
    cudaEventSynchronize( brakets_copied );  // wait for bra/ket to be copied to the GPU
    cudaEventSynchronize( H_done );          // wait H to be done

    // propagate particle
    chebyshev_gpu( n, ld, *tau, save_tau, *t_max, *t_init, bra, ket, H );    SAFE(cudaGetLastError());
    // in sync. here

    // copy propagated bra/ket back to the CPU
    cudaZcopyAsync_D2H( h_PSI_ket, ket, n, stream[0] );
    cudaZcopyAsync_D2H( h_PSI_bra, bra, n, stream[1] );

    // just a nickname to re-use the ket vector and save one allocation
    cuDoubleComplex *AO_bra = ket;
    
    // AO_bra = bra * S^(-1)
    kblas_dzgemv2_async( 'n', n, d_one, Sinv, ld, bra, d_zero, AO_bra, stream[0] );

    // copy AO_bra to the CPU
    cudaZcopyAsync_D2H( h_AO_bra, AO_bra, n, stream[0] );

    // wait bra/ket to be copied
    cudaStreamSynchronize( stream[1] );
    cudaStreamSynchronize( stream[0] );
 
    // destroy events
    cudaEventDestroy( S_inverted );
    cudaEventDestroy( brakets_copied );
    cudaEventDestroy( H_done );

    cublasSetStream( myHandle, cublas_default );

    first_call = false;
    time_end("Propagation");
}




//===================================================================
//-------------------------------------------------------------------
void ehrenfestkernel_gpu_(
    const int       * const __restrict__ N,
    const double    * const __restrict__ h_H,        // in host
    const double    * const __restrict__ h_A,        // in host
    const double    * const __restrict__ h_X,        // in host
    double          * const __restrict__ h_K )       // in host
{
    time_init();
    
    static bool first_call = true;
    cudaError_t stat;

    const int n  = *N;
    const int ld = ((n + 31)/32)*32;  // making leading dimension multiple of 32
    
    cudaEvent_t done;
    cudaEventCreateWithFlags( &done, cudaEventDisableTiming );
    
    static double *H = nullptr;
    static double *A = nullptr;
    static double *X = nullptr;
    static double *K = nullptr;
    
    if (first_call)
    {
        stat = cudaMalloc((void **) &H, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &A, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &X, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &K, n*ld*sizeof(double));           SAFE(stat);
        first_call = false;
    }
    
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_H, n, (void *)H, ld, stream[0]);
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_A, n, (void *)A, ld, stream[0]);
    
    // B = H' * A   (storing B in K to save memory)
    cublasSetStream( myHandle, stream[0] );
//     cublasDsymm( myHandle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, n, n, &d_one, A, ld, H, ld, &d_zero, K, ld);
    cublasDgemm( myHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &d_one, H, ld, A, ld, &d_zero, K, ld);
    cudaEventRecord( done, stream[0]);
    
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_X, n, (void *)X, ld, stream[1]);
    
    // K = X * A - B
    cudaStreamWaitEvent( stream[1], done, 0 );
    hadamard_minus(n, ld, X, A, K, stream[1]);
    
    cublasGetMatrixAsync( n, n, sizeof(double), (void *)K, ld, (void *)h_K, n, stream[1]);
    
//     cudaStreamSynchronize( stream[0] );
    cudaStreamSynchronize( stream[1] );
    cudaEventDestroy( done );
    
    time_end("EhrenfestKernel");
}


//-------------------------------------------------------------------
void ehrenfestkernel2_gpu_(
    const int             * const __restrict__ N,
    const cuDoubleComplex * const __restrict__ h_bra,      // in host
    const cuDoubleComplex * const __restrict__ h_ket,      // in host
    const double          * const __restrict__ h_H,        // in host
    const double          * const __restrict__ h_X,        // in host
    double                * const __restrict__ h_K )       // in host
{
    time_init();
    
    static bool first_call = true;
    cudaError_t stat;

    const int n  = *N;
    const int ld = ((n + 31)/32)*32;  // making leading dimension multiple of 32
    
    cudaEvent_t A_done, X_done;
    cudaEventCreateWithFlags( &A_done, cudaEventDisableTiming );
    cudaEventCreateWithFlags( &X_done, cudaEventDisableTiming );
    
    cuDoubleComplex *bra = nullptr;
    cuDoubleComplex *ket = nullptr;
    static double *H = nullptr;
    static double *A = nullptr;
    static double *X = nullptr;
    static double *K = nullptr;
    
    if (first_call)
    {
        stat = cudaMalloc((void **) &H, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &A, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &X, n*ld*sizeof(double));           SAFE(stat);
        stat = cudaMalloc((void **) &K, n*ld*sizeof(double));           SAFE(stat);
        first_call = false;
    }
    stat = cudaMalloc((void **) &bra, 2*ld*sizeof(cuDoubleComplex));    SAFE(stat);
    stat = cudaMalloc((void **) &ket, 2*ld*sizeof(cuDoubleComplex));    SAFE(stat);
    
    // copy bra/ket from host (CPU) to device (GPU)
    cublasSetMatrixAsync( n, 2, sizeof(double), (void *)h_bra, n, (void *)bra, ld, stream[0]);
    cublasSetMatrixAsync( n, 2, sizeof(double), (void *)h_ket, n, (void *)ket, ld, stream[0]);
    
    // ρ = Re{ ket(j,1)*bra(i,1) -  ket(j,2)*bra(i,2) }  ∀ i, j
    // A = (ρ + ρ^T) / 2
    calculate_A(n, ld, bra, ket, A, stream[0]);
    cudaEventRecord( A_done, stream[0]);
    
    // copy H from host to device
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_H, n, (void *)H, ld, stream[1]);
    
    // B = H' * A   (storing the result B in matrix K to save memory, that is, K==B)
    cublasSetStream( myHandle, stream[1] );
    cudaStreamWaitEvent( stream[1], A_done, 0 );
//     cublasDsymm( myHandle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, n, n, &d_one, A, ld, H, ld, &d_zero, K, ld);
    cublasDgemm( myHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &d_one, H, ld, A, ld, &d_zero, K, ld);
    
    // copy X from host to device
    cublasSetMatrixAsync( n, n, sizeof(double), (void *)h_X, n, (void *)X, ld, stream[0]);
    cudaEventRecord( X_done, stream[0]);
    
    // K = X * A - B
    cudaStreamWaitEvent( stream[1], X_done, 0 );
    hadamard_minus(n, ld, X, A, K, stream[1]);
    
    // copy K from device to host
    cublasGetMatrixAsync( n, n, sizeof(double), (void *)K, ld, (void *)h_K, n, stream[1]);
    
    cudaStreamSynchronize( stream[1] );
    cudaEventDestroy( A_done );
    cudaEventDestroy( X_done );
    cudaFree( bra );
    cudaFree( ket );
    
    time_end("EhrenfestKernel");
}
#endif
