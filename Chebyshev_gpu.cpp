#ifdef USE_GPU

#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cublas_v2.h"
#include "magma_operators.h"

#endif

#ifdef __INTEL_COMPILER
 #define __restrict__ restrict
 #include <mathimf.h>
#endif

#ifdef USE_GPU

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

// bra1 = alpha*H^T*bra2 - bra3
#define bra_H_m( VEC1, VEC2, H, ALPHA, VEC3 ) \
                                 Mat_op_vec5( H, 't', ALPHA, VEC2, d_neg_one, VEC3, VEC1, stream[bra_stream] )

// ket1 = alpha*H*ket2 - ket3
#define H_ket_m( VEC1, H, VEC2, ALPHA, VEC3 ) \
                                 Mat_op_vec5( H, 'n', ALPHA, VEC2, d_neg_one, VEC3, VEC1, stream[ket_stream] )

// find the (complex) element with maximum real or imaginary parts in VEC
#define findMax( VEC, MAX_PTR ) \
{ \
    cublasIdamax( myHandle, 2*n, (const double *)VEC, 1, i_max); \
    *i_max = (*i_max-1)/2;     /* From Fortran-Real to C-Complex index conversion */ \
    cudaZcopy_D2H( MAX_PTR, &VEC[*i_max], 1 ) \
}


// structure to hold allocated memory, used in convergence_gpu()
struct bras_kets_ptrs
{
    cuDoubleComplex *bras,     *kets,
                    *tmp_bra1, *tmp_ket1,
                    *tmp_bra2, *tmp_ket2,
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
const int max_order = 25;    // note: update zi_pow if changing max_order's value
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
const cuDoubleComplex zi_pow[] = {          // dimension >= max_order+1
    z_imag, z_one, z_neg_imag, z_neg_one,
    z_imag, z_one, z_neg_imag, z_neg_one,
    z_imag, z_one, z_neg_imag, z_neg_one,
    z_imag, z_one, z_neg_imag, z_neg_one,                             
    z_imag, z_one, z_neg_imag, z_neg_one,                             
    z_imag, z_one, z_neg_imag, z_neg_one,
    z_imag, z_one, z_neg_imag, z_neg_one }; // 28 elements here


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
    const int n, const int ld, double tau,
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

    cuDoubleComplex *bras, *kets, *tmp_bra, *tmp_ket, *coeff, *d_coeff;
    stat = cudaMalloc( (void **) &bras,   max_order*ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &kets,   max_order*ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &tmp_bra,          ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &tmp_ket,          ld*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMalloc( (void **) &d_coeff,   max_order*sizeof(cuDoubleComplex) );  SAFE(stat);
    stat = cudaMallocHost( (void **) &coeff, max_order*sizeof(cuDoubleComplex) );  SAFE(stat);
    
    cuDoubleComplex *bra[max_order], *ket[max_order];
    for(int i=0; i<max_order; ++i)
    {
        bra[i] = &bras[i*ld];
        ket[i] = &kets[i*ld];
    }

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
        braCopy( PSI_bra, bra[0] );        // <bra_0| = <PSI_bra|
        ketCopy( PSI_ket, ket[0] );        // |ket_0> = |PSI_ket>

        if( must_copy_coeff )
        {
            cudaZcopyAsync_H2D( d_coeff, coeff, max_order, stream[2] );   // copy coeff to the GPU
            must_copy_coeff = false;
        }

        bra_H( bra[1] ,bra[0], H );        // <bra_1| = <bra_0|H
        H_ket( ket[1], H, ket[0] );        // |ket_1> = H|ket_0>

        for( int k=2; k<k_ref; ++k )
        {
            bra_H_m( bra[k], bra[k-1], H, d_two, bra[k-2] );    // <bra_k| = 2*<bra_(k-1)|H - <bra_(k-2)|
            H_ket_m( ket[k], H, ket[k-1], d_two, ket[k-2] );    // |ket_k> = 2*H|ket_(k-1)> - |ket_(k-2)>
        }

        cudaStreamSynchronize( stream[2] );   // wait for d_coeff to be copied
        // <tmp_bra| = sum(k=0,k_ref) c_k*<bra_k|
        // |tmp_ket> = sum(k=0,k_ref) c_k*|ket_k>
        acummulate_vec_ax_async( n, ld, k_ref, d_coeff, bras, tmp_bra, stream[bra_stream] );
        acummulate_vec_ax_async( n, ld, k_ref, d_coeff, kets, tmp_ket, stream[ket_stream] );

        // check charge conservation
        //  norm = |<tmp_bra|tmp_ket>|^2
        setStream( bra );                              // queue after tmp_bra
        cudaStreamSynchronize( stream[ket_stream] );   // wait for tmp_ket
        cublasZdotc( myHandle, n, tmp_bra, 1, tmp_ket, 1, &z_norm );    // implicit synchronization here
        const double norm = cuCabs(z_norm);

        if( fabs(norm - norm_ref) < norm_tolerance )
        {
            braCopy( tmp_bra, PSI_bra );            // <PSI_bra| = <tmp_bra|
            ketCopy( tmp_ket, PSI_ket );            // |PSI_ket> = |tmp_ket>
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
        }

        t += tau*h_bar;

        if( t_max - t < tau*h_bar )
        {
            tau = (t_max - t)/h_bar;
            coefficient( tau, max_order, coeff );
            must_copy_coeff = true;
        }
    } // while (*t < t_max)

    // clean up
    cublasSetStream( myHandle, cublas_default );
    cudaDeviceSynchronize();   // wait for all calculations to complete

    cudaFree( bras );
    cudaFree( kets );
    cudaFree( tmp_bra );
    cudaFree( tmp_ket );
    cudaFree( d_coeff );
    cudaFreeHost( coeff );

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
    mem->tmp_bra1 = &dev_mem_chunk[offset]; offset += ld;
    mem->tmp_ket1 = &dev_mem_chunk[offset]; offset += ld;
    mem->tmp_bra2 = &dev_mem_chunk[offset]; offset += ld;
    mem->tmp_ket2 = &dev_mem_chunk[offset];

    cudaMallocHost((void**) &mem->pin_mem_chunk, sizeof(cuDoubleComplex) + sizeof(int));
}


//-------------------------------------------------------------------
void convergence_gpu_dealloc( struct bras_kets_ptrs * const mem )
{
    cudaFree( mem->bras );
    cudaFreeHost( mem->pin_mem_chunk );
}
#endif


//-------------------------------------------------------------------
// contribution of 'x' the Bessel's functions without the order's part (it's inverse, actually)
static inline double nakedBessel( const int n, const double x ) { return( (1 << (n - 2)) * (x*x + 4) / pow(x,n) ); };
// for fortran:
extern "C" double nakedbessel_( const int * const n, const double * const x) { return(nakedBessel(*n,*x)); }


#ifdef USE_GPU
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
    cuDoubleComplex *bras, *kets, *tmp_bra1, *tmp_ket1, *tmp_bra2, *tmp_ket2, *diff_bra, *diff_ket, *max;

    bras     = mem->bras;        kets     = mem->kets;
    tmp_bra1 = mem->tmp_bra1;    tmp_ket1 = mem->tmp_ket1;
    tmp_bra2 = mem->tmp_bra2;    tmp_ket2 = mem->tmp_ket2;
    max = mem->pin_mem_chunk;

    int * const i_max = reinterpret_cast<int*>(mem->pin_mem_chunk + 1);

    cuDoubleComplex *bra[max_order], *ket[max_order];
    for(int i=0; i<max_order; ++i)
    {
        bra[i] = &bras[i*ld];
        ket[i] = &kets[i*ld];
    }

    braCopy( PSI_bra, bra[0] );    // <bra_0| = <PSI_bra|
    ketCopy( PSI_ket, ket[0] );    // |ket_0> = |PSI_ket>

    bra_H( bra[1] ,bra[0], H );    // <bra_1| = <bra_0|H
    H_ket( ket[1], H, ket[0] );    // |ket_1> = H|ket_0>

    // get coefficients (in host, potentially overlaping with previous kernel launches)
    coefficient( tau, max_order, coeff );

    // <tmp_bra1| = c_0*<bra_0| + c_1*<bra_1|
    Zaxpby_async( n, coeff[0], bra[0], coeff[1], bra[1], tmp_bra1, stream[bra_stream] );

    // |tmp_ket1> = c_0*|ket_0> + c_1*|ket_1>
    Zaxpby_async( n, coeff[0], ket[0], coeff[1], ket[1], tmp_ket1, stream[ket_stream] );
    
    // new criteria for the number of terms in expansion: coeff[k] lesser than a suitable small value
    int k_max = max_order;
    for(int k=6; k<max_order; ++k)  // minimum of six terms
    {
        if( cuCabs(coeff[k]*nakedBessel(k,tau)) < 1.0e-20 )
        {
            k_max = k;
            break;
        }
    }
    *k_ref = k_max;

    *ok = false;
    for( int k=2; k<k_max; ++k )
    {
        // <bra_k| = 2*<bra_(k-1)|H - <bra_(k-2)|
        // |ket_k> = 2*H|ket_(k-1)> - |ket_(k-2)>
        bra_H_m( bra[k], bra[k-1], H, d_two, bra[k-2] );
        H_ket_m( ket[k], H, ket[k-1], d_two, ket[k-2] );

//      *k_ref = k + 1;   // old criteria for the number of terms in expansion: until full convergence

        // <tmp_bra2| = <tmp_bra1| + c_k*<bra_k|   &&   <diff| = <tmp_bra2| - <tmp_bra1|
        // |tmp_ket2> = |tmp_ket1> + c_k*|ket_k>   &&   |diff> = |tmp_ket2> - |tmp_ket1>
        fused_Zxpby_and_subtract( n, tmp_bra1, coeff[k], bra[k], tmp_bra2, tmp_bra1, stream[bra_stream] );
        fused_Zxpby_and_subtract( n, tmp_ket1, coeff[k], ket[k], tmp_ket2, tmp_ket1, stream[ket_stream] );

        //  find maximum element of |diff>
        setStream( bra );
        findMax( tmp_bra1, max );      // implicit synchronization here

        if( cuCabs(*max) < tolerance ) // if bra is converged
        {
            //  find maximum element of <diff|
            setStream( ket );
            findMax( tmp_ket1, max );      // implicit synchronization here

            if( cuCabs(*max) < tolerance ) // if ket is converged
            {
                // check charge conservation:
                //  norm = |<tmp_bra2|tmp_ket2>|^2
                cuDoubleComplex z_norm;
                cublasZdotc( myHandle, n, tmp_bra2, 1, tmp_ket2, 1, &z_norm );    // implicit synchronization here
                const double norm = cuCabs(z_norm);

                if( fabs(norm - norm_ref) < norm_tolerance ) // if charge is conserved
                {
                    // copy answer back and exit
                    braCopy( tmp_bra2, PSI_bra );    // <PSI_bra| = <tmp_bra2|
                    ketCopy( tmp_ket2, PSI_ket );    // |PSI_ket> = |tmp_ket2>
                    *ok = true;
                    break;
                }
            } // ket conv.
        } // bra conv.

        // No need to copy, just swap pointers!
        swap<cuDoubleComplex*>( tmp_bra1, tmp_bra2 );    // <tmp_bra1| = <tmp_bra2|
        swap<cuDoubleComplex*>( tmp_ket1, tmp_ket2 );    // |tmp_ket1> = |tmp_ket2>
        
    } // for( k=2; k<max_order; ++k )

    // clean up
    cublasSetStream( myHandle, cublas_default );
    cudaDeviceSynchronize();
}


//-------------------------------------------------------------------
void coefficient( const double tau, const int k_max, cuDoubleComplex * const coeff )
{
    coeff[0].x = jn(0, tau);
    coeff[0].y = 0.0;

    for( int k=1; k<k_max; ++k )
        coeff[k] = (d_two * jn(k, tau)) * zi_pow[k+1];
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

#endif
