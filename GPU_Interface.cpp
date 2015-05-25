#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef FORTRAN
  #define FORTRAN
#endif
#ifdef FORTRAN
  #define ADD_
#endif

#ifdef USE_GPU
  #warning "Compiling with GPU support"
  #include <cuda.h>
  #include <cuda_runtime_api.h>
  #include "cublas_v2.h"
  #include "magma.h"
  #include "magma_dbulge.h"
#else
  #warning "Compiling only for CPU"
#endif

// Memory macros
#ifdef GPU_PIN_MEM_WORK
  #define malloc_work( A, DIM, TYPE )   cudaMallocHost( (void **) &A, DIM*sizeof(TYPE) )
  #define free_work( A )                cudaFreeHost( A )
#else
  #define malloc_work( A, DIM, TYPE )   A = (TYPE*) malloc( DIM*sizeof(TYPE) )
  #define free_work( A )                free( A )
#endif

// Timing macros
#ifdef GPU_TIMING
  #include <time.h>
  static inline double get_wtime()
  {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return( t.tv_sec +  t.tv_nsec/1.0e9 );
  };

  #warning "Timing enabled"
  #define time_init()  double __my_time = get_wtime()
  #define time_end()   __my_time = get_wtime()-__my_time; printf("Time %s: %.4f\n", __func__,__my_time); fflush(stdout)
#else
  #warning "Timing disabled"
  #define time_init()  /* No timing */
  #define time_end()   /* No timing */
#endif

// Debug
#ifdef GPU_DEBUG
  #define DEBUG( A ) A
#else
  #define DEBUG( A ) /*Nothing*/
#endif

#define CHECK_INFO(info, F)  if(info!= 0) { printf("error: info = %i in %s: %s\n", info, __func__, F); exit(EXIT_FAILURE); }

#ifdef USE_GPU
  #define SAFE(S) if(S!=cudaSuccess) printf("ERROR(%i): %s\n%s\n\n", __LINE__, cudaGetErrorName(S), cudaGetErrorString(S));
#endif

// Complex type
#ifdef USE_GPU
  #define DoubleComplex cuDoubleComplex
#else
  typedef struct{ double x, y; } DoubleComplex;
#endif

// Matrix (fortran) to array indexing
#define IDX( i, j, LD ) ((i) + (j)*(LD))

// Name mangling for use with fortran:
#ifdef FORTRAN
  #define GPU_Init      gpu_init_
  #define GPU_Finalize  gpu_finalize_
  
  #define GPU_Pin       gpu_pin_
  #define GPU_Unpin     gpu_unpin_

  #define xPU_dzgemv    xpu_dzgemv_

  #define xPU_dgemm     xpu_dgemm_
  #define xPU_dsymm     xpu_dsymm_
  #define xPU_dzgemm    xpu_dzgemm_

  #define xPU_dsygvd    xpu_dsygvd_
  #define xPU_dsygvd2s  xpu_dsygvd2s_
  #define xPU_syInvert  xpu_syinvert_
#endif


//----------------------------------------------------------------
// External functions
extern "C"
{
    int magma_get_parallel_numthreads();
    void dsygvd_(const int* itype, const char* jobZ, const char* UplO, const int* N, 
                 double* A, const int* ldA, double* B, const int* ldB, double* W, 
                 double* work, const int *lwork, int* iwork, const int* liwork, int *INFO);
    void dgemm_(const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                const double * const alpha, double *A, const int *ldA, double *B, const int *ldB, 
                const double * const beta, double *C, const int *ldC);
    void dzgemm_(const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                 const DoubleComplex * const alpha, double *A, const int *ldA, const DoubleComplex *B, const int *ldB, 
                 const DoubleComplex * const beta, DoubleComplex *C, const int *ldC);
    void dsymm_( const char *side, const char *UpLo, const int * const M, const int * const N, 
                 const double * const alpha, double *hA, const int * const LDA,
                 double *hB, const int * const LDB, const double * const beta,
                 double *hC, const int * const LDC );
    void dgemv_( const char * const transA, const int * const M, const int * const N, 
                 const double * const alpha, double * const hA, const int * const LDA,
                 double * const hX, const int * const incX, const double * const beta, 
                 double * const hY, const int * const incY );
    void dzgemv_( const char *transA, const int * const M, const int * const N, 
                  const DoubleComplex *alpha, double *hA, const int * const LDA,
                  DoubleComplex *hX, const int *incX, const DoubleComplex *beta, 
                  DoubleComplex *hY, const int *incY );
    void dsytrf_( const char *UpLo, const int * const M, double *A, const int * const N, int *ipiv, double *work, const int *lwork, int *info );
    void dsytri_( const char *UpLo, const int * const M, double *A, const int * const N, int *ipiv, double *work, int *info);
}


//----------------------------------------------------------------
// Function prototypes:

// 1. Exported:
extern "C"
{
    void GPU_Init( const int *dev, const int *procs_per_dev );
    void GPU_Finalize( void );
    
    void GPU_Pin( void * ptr, int * size );
    void GPU_Unpin( void * ptr );

    void xPU_dzgemv( const char *transA, const int * const M, const int * const N, 
                     const DoubleComplex *alpha, double *hA, const int * const LDA,
                     DoubleComplex *hX, const int *incX, const DoubleComplex *beta, 
                     DoubleComplex *hY, const int *incY );

    void xPU_dgemm( const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                    const double * const alpha, double *A, const int * const LDA,
                    double *B, const int * const LDB, const double * const beta, double *C, const int * const LDC );
    void xPU_dsymm( const char *side, const char *UpLo, const int * const M, const int * const N, 
                    const double * const alpha, double *hA, const int * const LDA,
                    double *hB, const int * const LDB, const double * const beta,
                    double *hC, const int * const LDC );
    void xPU_dzgemm( const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                     const DoubleComplex * const alpha, double *A, const int * const LDA,
                     const DoubleComplex * const B, const int * const LDB, const DoubleComplex * const beta, DoubleComplex *C, const int * const LDC );

    void xPU_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                     double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
    void xPU_dsygvd2s( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                       double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
    void xPU_syInvert( double *A, const char *UpLo, const int * const N, int *info);
}
#ifdef USE_GPU
void gpu_dgeInvert( double *dA, const int n, const int lddA, cudaStream_t stream );
#endif


// 2. Internal:
//  - Helper
static inline int max( const int a, const int b ) { return( a>b ? a:b ); };
static inline int min( const int a, const int b ) { return( a<b ? a:b ); };
#ifdef USE_GPU
static inline magma_vec_t       lapack_to_magma_jobZ( const char jobZ );
static inline magma_uplo_t      lapack_to_magma_uplo( const char uplo );
static inline cublasFillMode_t  lapack_to_cublas_uplo( const char UpLo );
static inline cublasSideMode_t  lapack_to_cublas_side( const char side );
static inline cublasOperation_t lapack_to_cublas_trans( const char trans );
static void cpu_dtranspose_sq_inplace( const int n, double * const A, const int ldA );
static void cpu_dtranspose_naive( const int n, double * const A, const int ldA );
#endif

//  - Compute
static inline void cpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void cpu_syInvert(const char *UpLo, double *A, const int * const N, int *info);
#ifdef USE_GPU
static inline void gpu_dsygvd_2s( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                                  double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_dsygvd_m( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                                 double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_syInvert(const char *UpLo, double *hA, const int * const N, int *info);
#endif



/******************************
 *     Globlal variables
 ******************************/

// Device count and management
int devCount = 0;
int VirtualDevCount = 0;
int IuseGPU = 0; // false
int myGPU = -1;  // none

#ifdef USE_GPU
// CUBLAS handle
cublasHandle_t myHandle;
cudaStream_t cublas_default;

// streams
const int nStreams = 3;
cudaStream_t stream[nStreams];
#endif


/******************************
 *      Helper functions
 ******************************/

//-----------------------------------------
// Initializes the environment
//
// Alberto Torres
void GPU_Init(const int *pid, const int *procs_per_dev)
{
#ifdef USE_GPU
    cudaGetDeviceCount(&devCount);
    VirtualDevCount = devCount * (*procs_per_dev);   // We can pin more than one process to a gpu with procs_per_dev > 1
    if(*pid < VirtualDevCount )
    {
        IuseGPU=1; // true
        magma_init();

        if(*pid==0)
        {
            printf("\nUsing: (%i processes per GPU)\n",*procs_per_dev);
            magma_print_environment();
            printf("\n");
            fflush(stdout);
        }
 
        myGPU = (*pid / *procs_per_dev) % devCount;
        cudaSetDevice( myGPU );                          // Bind to GPUs in a round-robin way according to pid
        cublasCreate( &myHandle );
        cublasSetAtomicsMode( myHandle, CUBLAS_ATOMICS_ALLOWED );   // Enables faster *symm (see cublas manual)
        cublasGetStream( myHandle, &cublas_default );

        for(int i=0; i<nStreams; ++i)
            cudaStreamCreate( &stream[i] );

        printf("Process nr. %i using GPU device nr. %i of %i(%i)\n", *pid, myGPU, devCount, VirtualDevCount);
        fflush(stdout);
    }
#else
    printf("Using only CPU\n");
    fflush(stdout);
#endif
}


//-----------------------------------------
// Finalizes the environment
//
void GPU_Finalize(void)
{
#ifdef USE_GPU
    if(IuseGPU)
    {
        for(int i=0; i<nStreams; ++i)
            cudaStreamDestroy( stream[i] );
        cublasDestroy( myHandle );
        magma_finalize();
    }
#endif
}

//-----------------------------------------
// Memory
//
void GPU_Pin( void * ptr, int * size )
{
#ifdef USE_GPU
    time_init();
    cudaError_t stat = cudaHostRegister( ptr, *size, cudaHostRegisterDefault );     SAFE(stat);
    time_end();
#endif
}

void GPU_Unpin( void * ptr )
{
#ifdef USE_GPU   
    cudaError_t stat = cudaHostUnregister( ptr );     SAFE(stat);
#endif
}


//-----------------------------------------
// Flag conversions
// No parameter checking is done!
#ifdef USE_GPU
// magma:
static inline magma_vec_t lapack_to_magma_jobZ( const char jobZ )
    { return( jobZ == 'V' || jobZ == 'v' ? MagmaVec : MagmaNoVec ); }

static inline magma_uplo_t lapack_to_magma_uplo( const char UpLo )
    { return( UpLo == 'U' || UpLo == 'u' ? MagmaUpper : MagmaLower ); }

// cublas:
static inline cublasFillMode_t lapack_to_cublas_uplo( const char UpLo )
    { return( UpLo == 'U' || UpLo == 'u' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER ); }
    
static inline cublasSideMode_t lapack_to_cublas_side( const char side )
    { return( side == 'L' || side == 'l' ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT ); }

static inline cublasOperation_t lapack_to_cublas_trans( const char trans )
    { return( trans == 'N' || trans == 'n' ? CUBLAS_OP_N : (trans == 'T' || trans == 't' ? CUBLAS_OP_T : CUBLAS_OP_C) ); }
#endif


//-----------------------------------------
// Transpose square matrix
// Half copy blocked algorithm
static void cpu_dtranspose_sq_inplace( const int n, double * const A, const int ldA )
{
    const int bl = 16; // For 32 KB L1 cache

    #pragma omp parallel default(none) shared(bl,n,ldA,A)
    {
        double * const buffer = (double *) alloca( bl*bl*sizeof(double) );

        #pragma omp for schedule(dynamic)
        for(int jj=0; jj<n; jj+=bl)
        {
            const int j_max = min( n, jj+bl );
            
            for(int ii=0; ii<jj; ii+=bl)  // Cycles trhough upper triangle
            {
                const int i_max = min( n, ii+bl );

                // buffer lower block: buffer = A_jj,ii
                for(int i=ii; i<i_max; i++)
                    memcpy( &buffer[ IDX(0,i-ii,bl) ], &A[ IDX(jj,i,ldA) ], (j_max - jj)*sizeof(double) );

                // copy/transpose: block A_jj,ii = A_ii,jj^T
                for(int j=jj; j<j_max; j++)
                    #pragma ivdep
                    for(int i=ii; i<i_max; i++)
                        A[ IDX(j,i,ldA) ] = A[ IDX(i,j,ldA) ];

                    // copy/transpose: block A_ii,jj = buffer^T
                    for(int j=jj; j<j_max; j++)
                        #pragma ivdep
                        for(int i=ii; i<i_max; i++)
                            A[ IDX(i,j,ldA) ] = buffer[ IDX(j-jj,i-ii,bl) ];
            }

            // transpose diagonal block
            for( int j=jj+1; j<min( jj+bl, n ); j++ )
                #pragma ivdep
                for( int i=jj; i<min( j, ldA ); i++ )
                {
                    double volatile tmp = A[ IDX(j,i,ldA) ];
                    A[ IDX(j,i,ldA) ] = A[ IDX(i,j,ldA) ];
                    A[ IDX(i,j,ldA) ] = tmp;
                }
        }

        // free( buffer );
    }
}


static void cpu_dtranspose_naive( const int n, double * const A, const int ldA )
{
    for( int j=0; j<n; j++ )
        for( int i=0; i<j; i++ )
        {
            double volatile tmp = A[ IDX(j,i,ldA) ];
            A[ IDX(j,i,ldA) ] = A[ IDX(i,j,ldA) ];
            A[ IDX(i,j,ldA) ] = tmp;
        }
}



/****************************************
 *     "Compute" functions - BLAS 2
 ****************************************/

//===================================================================
// Matrix-vector multiplication - in CPU or GPU
// mixed precision: A*x = y
// double  : A
// complex : x, y
//
// Alberto Torres
void xPU_dzgemv( const char * const transA, const int * const M, const int * const N, 
                 const DoubleComplex * const alpha, double * const hA, const int * const LDA,
                 DoubleComplex * const hX, const int * const incX, const DoubleComplex * const beta, 
                 DoubleComplex * const hY, const int * const incY )
{
    time_init();
#ifdef USE_GPU
    if(IuseGPU)
    {
        const int m = *M;
        const int n = *N;
        const int ldA = *LDA;
        const int incx = *incX;
        const int incx2 = 2*incx;
        const int incy = *incY;
        const int incy2 = 2*incy;
        const cublasOperation_t opA = lapack_to_cublas_trans( *transA );
        int rows, cols, ldX, ldY;

        if( opA == CUBLAS_OP_N )
        {
            rows = m;
            cols = n;
            ldX = 1+(n-1)*abs(incx);
            ldY = 1+(m-1)*abs(incy);
        }
        else
        {
            rows = n;
            cols = m;
            ldX = 1+(m-1)*abs(incx);
            ldY = 1+(n-1)*abs(incy);
        }

        const int lddA = ( (ldA + 31)/32 )*32;
        const int lddx = ( (ldX + 31)/32 )*32;
        const int lddy = ( (ldY + 31)/32 )*32;

        double *dA, *dx, *dy;

        // Allocate matrices in the GPU
        cudaMalloc( (void**) &dA, lddA*n*sizeof(double) );
        cudaMalloc( (void**) &dx, lddx*sizeof(DoubleComplex) );
        cudaMalloc( (void**) &dy, lddy*sizeof(DoubleComplex) );

        // pointers to the real and imaginary parts (device)
        double * const dx_Re = dx;
        double * const dx_Im = &dx[lddx];
        double * const dy_Re = dy;
        double * const dy_Im = &dy[lddy];
        
        cudaStream_t stream1, stream2;
        cudaStreamCreate( &stream1 );
        cudaStreamCreate( &stream2 );

        // Copy from CPU to GPU
        cublasSetMatrix( m, n, sizeof(double), (void *) hA, ldA, (void *) dA, lddA);
        
        // real part
        cublasSetStream( myHandle, stream1 );
        cublasSetVectorAsync( cols, sizeof(double), (void *) &hX->x, incx2, (void *) dx_Re, incx, stream1 );
        cublasDgemv( myHandle, opA, m, n, (const double *)alpha, dA, lddA, dx_Re, incx, (const double *)beta, dy_Re, incy );
        cublasGetVectorAsync( cols, sizeof(double), (void *) dy_Re, incy, (void *) &hY->x, incx2, stream1 );

        // imag part
        cublasSetStream( myHandle, stream2 );
        cublasSetVectorAsync( cols, sizeof(double), (void *) &hX->y, incx2, (void *) dx_Im, incx, stream2 );
        cublasDgemv( myHandle, opA, m, n, (const double *)alpha, dA, lddA, dx_Im, incx, (const double *)beta, dy_Im, incy );
        cublasGetVectorAsync( cols, sizeof(double), (void *) dy_Im, incy, (void *) &hY->y, incx2, stream2 );

        // Compute
        //cublasDgemm( myHandle, opA, CUBLAS_OP_N, rows, 2, cols, (const double *)alpha, dA, lddA, dx, lddx, (const double *)beta, dy, lddy);

        cudaStreamSynchronize( stream1 );
        cudaStreamSynchronize( stream2 );

        // Deallocate
        cudaFree(dA);
        cudaFree(dx);
        cudaFree(dy);

        cudaStreamDestroy(stream1);
        cudaStreamDestroy(stream2);

        cublasSetStream( myHandle, cublas_default );
    }
    else
#endif
        dzgemv_( transA, M, N, alpha, hA, LDA, hX, incX, beta, hY, incY );

    time_end();
}


 
/****************************************
 *     "Compute" functions - BLAS 3
 ****************************************/

//===================================================================
// Matrix-matrix multiplication - in CPU or GPU
//
// Alberto Torres
void xPU_dgemm( const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                const double * const alpha, double *hA, const int * const LDA,
                double *hB, const int * const LDB, const double * const beta,
                double *hC, const int * const LDC )
{
    time_init();

#ifdef USE_GPU
    if(IuseGPU)
    {
        const int m = *M;
        const int n = *N;
        const int k = *K;

        const int ldA = *LDA;
        const int ldB = *LDB;
        const int ldC = *LDC;

        const int lddA = ( (ldA + 31)/32 )*32;
        const int lddB = ( (ldB + 31)/32 )*32;
        const int lddC = ( (ldC + 31)/32 )*32;

        const cublasOperation_t opA = lapack_to_cublas_trans( *transA );
        const cublasOperation_t opB = lapack_to_cublas_trans( *transB );

        int rA, cA;
        int rB, cB;

        if( opA == CUBLAS_OP_N ) { rA=m; cA=k; }
        else { rA=k; cA=m; }

        if( opB == CUBLAS_OP_N ) { rB=k; cB=n; }
        else { rB=n; cB=k; }

        double *dA, *dB, *dC;

        // Allocate matrices in the GPU
        cudaMalloc( (void**) &dA, lddA*cA*sizeof(double) );
        cudaMalloc( (void**) &dB, lddB*cB*sizeof(double) );
        cudaMalloc( (void**) &dC, lddC*n*sizeof(double) );

        // Copy from CPU to GPU
        cublasSetMatrix( rA, cA, sizeof(double), (void *)hA, ldA, (void *)dA, lddA);
        cublasSetMatrix( rB, cB, sizeof(double), (void *)hB, ldB, (void *)dB, lddB);

        // Compute
        cublasDgemm( myHandle, opA, opB, m, n, k, alpha, dA, lddA, dB, lddB, beta, dC, lddC);

        // Copy from GPU to CPU
        cublasGetMatrix( m, n, sizeof(double), (void *)dC, lddC, (void *)hC, ldC);

        // Deallocate
        cudaFree(dA); cudaFree(dB); cudaFree(dC);
    }
    else
#endif
        dgemm_( transA, transB, M, N, K, alpha, hA, LDA, hB, LDB, beta, hC, LDC);

    time_end();
}


//===================================================================
// Matrix-matrix multiplication - in CPU or GPU
//
// Alberto Torres
void xPU_dsymm( const char * const Side, const char * const UpLo, const int * const M, const int * const N, 
                const double * const alpha, double *hA, const int * const LDA,
                double *hB, const int * const LDB, const double * const beta,
                double *hC, const int * const LDC )
{
    time_init();
#ifdef USE_GPU
    if(IuseGPU)
    {
        const int m = *M;
        const int n = *N;

        const int ldA = *LDA;
        const int ldB = *LDB;
        const int ldC = *LDC;

        const int lddA = ( (ldA + 31)/32 )*32;
        const int lddB = ( (ldB + 31)/32 )*32;
        const int lddC = ( (ldC + 31)/32 )*32;

        const cublasSideMode_t side = lapack_to_cublas_side( *Side );
        const cublasFillMode_t uplo = lapack_to_cublas_uplo( *UpLo );

        int rA, cA;

        if( side == CUBLAS_SIDE_LEFT )  { rA = n; cA = m; }
        else                            { rA = m; cA = n; }

        double *dA, *dB, *dC;

        // Allocate matrices in the GPU
        cudaMalloc( (void**) &dA, lddA*cA*sizeof(double) );
        cudaMalloc( (void**) &dB, lddB*n*sizeof(double) );
        cudaMalloc( (void**) &dC, lddC*n*sizeof(double) );

        // Copy from CPU to GPU
        cublasSetMatrix( rA, cA, sizeof(double), (void *)hA, ldA, (void *)dA, lddA);
        cublasSetMatrix( m, n, sizeof(double), (void *)hB, ldB, (void *)dB, lddB);

        // Compute
        cublasDsymm( myHandle, side, uplo, m, n, alpha, dA, lddA, dB, lddB, beta, dC, lddC);

        // Copy from GPU to CPU
        cublasGetMatrix( m, n, sizeof(double), (void *)dC, lddC, (void *)hC, ldC);

        // Deallocate
        cudaFree(dA); cudaFree(dB); cudaFree(dC);
    }
    else
#endif
        dsymm_( Side, UpLo, M, N, alpha, hA, LDA, hB, LDB, beta, hC, LDC);

    time_end();
}


//===================================================================
// Matrix-matrix multiplication - in CPU or GPU
// mixed type D Z
// Alberto Torres
void xPU_dzgemm( const char *transA, const char *transB, const int * const M, const int * const N, const int * const K, 
                 const DoubleComplex * const alpha, double * const hA, const int * const LDA,
                 const DoubleComplex * const hB, const int * const LDB, const DoubleComplex * const beta,
                 DoubleComplex * const hC, const int * const LDC )
{
    time_init();
// #ifdef USE_GPU
//     if(IuseGPU)
//     {
//         const int m = *M;
//         const int n = *N;
//         const int k = *K;
//         
//         const int ldA = *LDA;
//         const int ldB = *LDB;
//         const int ldC = *LDC;
//         
//         const int lddA = ( (ldA + 31)/32 )*32;
//         const int lddB = ( (ldB + 31)/32 )*32;
//         const int lddC = ( (ldC + 31)/32 )*32;
//         
//         const cublasOperation_t opA = lapack_to_cublas_trans( *transA );
//         const cublasOperation_t opB = lapack_to_cublas_trans( *transB );
//         
//         int rA, cA;
//         int rB, cB;
//         
//         if( opA == CUBLAS_OP_N ) { rA=m; cA=k; }
//         else { rA=k; cA=m; }
//         
//         if( opB == CUBLAS_OP_N ) { rB=k; cB=n; }
//         else { rB=n; cB=k; }
//         
//         double *dA;
//         DoubleComplex *dB, *dC;
//         
//         // Allocate matrices in the GPU
//         cudaMalloc( (void**) &dA, lddA*cA*sizeof(double) );
//         cudaMalloc( (void**) &dB, lddB*cB*sizeof(DoubleComplex) );
//         cudaMalloc( (void**) &dC, lddC*n*sizeof(DoubleComplex) );
//         
//         // Copy from CPU to GPU
//         cublasSetMatrix( rA, cA, sizeof(DoubleComplex), (void *)hA, ldA, (void *)dA, lddA);
//         cublasSetMatrix( rB, cB, sizeof(DoubleComplex), (void *)hB, ldB, (void *)dB, lddB);
//         
//         // Compute
//         cublasDgemm( myHandle, opA, opB, m, n, k, &alpha->x, dA, lddA, dB, lddB, &beta->x, dC, lddC);
//         
//         // Copy from GPU to CPU
//         cublasGetMatrix( m, n, sizeof(DoubleComplex), (void *)dC, lddC, (void *)hC, ldC);
//         
//         // Deallocate
//         cudaFree(dA); cudaFree(dB); cudaFree(dC);
//     }
//     else
// #endif
        dzgemm_( transA, transB, M, N, K, alpha, hA, LDA, hB, LDB, beta, hC, LDC);
    
    time_end();
}



/****************************************
 *     "Compute" functions - LAPACK
 ****************************************/

//===================================================================
// Eigenvalues and eigenvectors - Ax = lambda Bx - dispatcher
//
// Alberto Torres
void xPU_dsygvd(const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                double *A, const int *ldA, double *B, const int *ldB, double *W, int *info)
{
    time_init();

#ifdef USE_GPU
    if(IuseGPU)
    {
#ifdef GPU_SYGVDM_VER
#warning "GPU_dsygvd calls magma_dsygvd_m"
        gpu_dsygvd_m( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, info );
#else
#warning "GPU_dsygvd calls magma_dsygvd"
        gpu_dsygvd( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, info );
#endif
    }
    else
#endif
        cpu_dsygvd( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, info );

    time_end();
}


//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - Dispatchs to CPU or GPU
//
// Alberto Torres
void xPU_dsygvd2s( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                   double *A, const int *ldA, double *B, const int *ldB, double *W, int *info )
{
    time_init();

#ifdef USE_GPU
    if(IuseGPU)
        gpu_dsygvd_2s( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, info );
    else
#endif
        cpu_dsygvd( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, info );

    time_end();
}


//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - on GPU (2 stages), called internally
//
// Alberto Torres
static inline void gpu_dsygvd_2s( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                                  double *A, const int *ldA, double *B, const int *ldB, double *W, int *info )
{
#ifdef USE_GPU
    const int n = *N;
    const int lq2 = magma_dbulge_get_lq2( n, magma_get_parallel_numthreads() );
    const int liwork = 3 + 5*n;
    const int lwork  = 1 + 6*n + 2*n*n + lq2;
    const magma_vec_t  vec  = lapack_to_magma_jobZ( *jobZ );
    const magma_uplo_t uplo = lapack_to_magma_uplo( *UpLo );
    int nZ;
    int *iwork;
    double *work;

    malloc_work( iwork, liwork, int );
    malloc_work( work, lwork, double );

    DEBUG( printf("gpu_dsygvd_2s: n= %i lwork= %i liwork= %i\n", n, lwork, liwork); fflush(stdout); )

    // Only implemented for MagmaLower in magma-1.6.0
    if( uplo ==  MagmaUpper ) cpu_dtranspose_sq_inplace( n, A, *ldA );

    magma_dsygvdx_2stage( *itype, vec, MagmaRangeAll, MagmaLower, n, A, *ldA, B, *ldB, 
                          0.0, 0.0, 1, n, &nZ, W, work, lwork, iwork, liwork, info );

    free_work(work); free_work(iwork);
#endif
}

//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - on GPU (multi), called internally
//
// Alberto Torres
static inline void gpu_dsygvd_m( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                                 double *A, const int *ldA, double *B, const int *ldB, double *W, int *info )
{
#ifdef USE_GPU
    const int n = *N;
    const int liwork = 3 + 5*n;
    const int lwork  = max( 1 + 6*n + 2*n*n, (2 + magma_get_dsytrd_nb(n))*n );
    const magma_vec_t  vec  = lapack_to_magma_jobZ( *jobZ );
    const magma_uplo_t uplo = lapack_to_magma_uplo( *UpLo );
    int *iwork;
    double *work;

    malloc_work( iwork, liwork, int );
    malloc_work( work, lwork, double );

    DEBUG( printf("gpu_dsygvd_m: n= %i lwork= %i liwork= %i\n", n, lwork, liwork); fflush(stdout); )
    magma_dsygvd_m( devCount, *itype, vec, uplo, n, A, *ldA, B, *ldB, W, work, lwork, iwork, liwork, info );

    free_work(work); free_work(iwork);
#endif
}

//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - on GPU, called internally
//
// Alberto Torres
static inline void gpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info )
{
#ifdef USE_GPU
    const int n = *N;
    const int liwork = 3 + 5*n;
    const int lwork  = max( 1 + 6*n + 2*n*n, (2 + magma_get_dsytrd_nb(n))*n );
    const magma_vec_t  vec  = lapack_to_magma_jobZ( *jobZ );
    const magma_uplo_t uplo = lapack_to_magma_uplo( *UpLo );
    int *iwork;
    double *work;

    malloc_work( iwork, liwork, int );
    malloc_work( work, lwork, double );

    DEBUG( printf("gpu_dsygvd: n= %i lwork= %i liwork= %i\n", n, lwork, liwork); fflush(stdout); )
    magma_dsygvd( *itype, vec, uplo, n, A, *ldA, B, *ldB, W, work, lwork, iwork, liwork, info );

    free_work(work); free_work(iwork);
#endif
}

//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - on CPU, called internally
//
// Alberto Torres
static inline void cpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int * const N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info )
{
    const int n = *N;
    const int liwork = 3 + 5*n;
    const int lwork = 1 + 6*n + 2*n*n;

    int *iwork = (int*) malloc( liwork * sizeof(int) );
    double *work = (double*) malloc( lwork * sizeof(double) );

    DEBUG( printf("cpu_dsygvd: n= %i lwork= %i liwork= %i\n", n, lwork, liwork); fflush(stdout); )
    dsygvd_( itype, jobZ, UpLo, N, A, ldA, B, ldB, W, work, &lwork, iwork, &liwork, info );

    free(work); free(iwork);
}



//===================================================================
// Matrix inversion - Dispatchs to CPU or GPU
//
// Alberto Torres
void xPU_syInvert( double *A, const char *UpLo, const int * const N, int *info )
{
    time_init();

#ifdef USE_GPU
    if(IuseGPU)
        gpu_syInvert( UpLo, A, N, info);
    else
#endif
        cpu_syInvert( UpLo, A, N, info );

    time_end();
}


//-------------------------------------------------------------------
// Symmetric Matrix inversion - on GPU
//
// Alberto Torres
static inline void gpu_syInvert( const char *UpLo, double *hA, const int * const N, int *info)
{
#ifdef USE_GPU
    const int n = *N;
    const int lwork = magma_get_dgetri_nb(n)*n;
    const int lddA = ( (n + 31)/32 )*32;  // Make lddA multiple of 32 (faster)
    int *ipiv = (int*) malloc( n*sizeof(int) );
    double *dA, *dwork;

    magmablasSetKernelStream( cublas_default );

    magma_dmalloc( &dA, n*lddA );
    magma_dmalloc( &dwork, lwork );

    magma_dsetmatrix( n, n, hA, n, dA, lddA );

    // dsytrf/i are not impemented in magma-1.6.1
    // Remember to remove the #ifdef in Matrix_Math.F90 -> Matrix_syInvert when changing to sy versions
    DEBUG( printf("gpu_syInvert: n= %i lwork= %i\n", n, lwork); fflush(stdout); )
    magma_dgetrf_gpu( n, n, dA, lddA, ipiv, info );              CHECK_INFO(*info, "magma_dgetrf_gpu");
    magma_dgetri_gpu( n, dA, lddA, ipiv, dwork, lwork, info );   CHECK_INFO(*info, "magma_dgetri_gpu");

    magma_dgetmatrix( n, n, dA, lddA, hA, n );

    magma_free(dA); magma_free(dwork); free(ipiv);
#endif
}


#ifdef USE_GPU
void gpu_dgeInvert( double *dA, const int n, const int lddA, cudaStream_t stream )
{
    time_init();

    const int lwork = magma_get_dgetri_nb(n)*n;
    int info;
    int *ipiv = (int*) malloc( n*sizeof(int) );
    double *dwork;

    magmablasSetKernelStream( cublas_default );  // doesn't work in non-default stream!

    magma_dmalloc( &dwork, lwork );

    magma_dgetrf_gpu( n, n, dA, lddA, ipiv, &info );              CHECK_INFO(info, "magma_dgetrf_gpu");
    magma_dgetri_gpu( n, dA, lddA, ipiv, dwork, lwork, &info );   CHECK_INFO(info, "magma_dgetri_gpu");

    magma_free(dwork); free(ipiv);

    time_end();
}
#endif

//-------------------------------------------------------------------
// Symmetric Matrix inversion - on CPU
//
// Alberto Torres
static inline void cpu_syInvert( const char *UpLo, double *A, const int * const N, int *info )
{
    const int n = *N;
    const int lwork = 64*n;

    int *ipiv = (int*) malloc( n*sizeof(int) );
    double *work = (double*) malloc( lwork*sizeof(double) );

    DEBUG( printf("cpu_syInvert: n= %i lwork= %i\n", n, lwork); fflush(stdout); )
    dsytrf_( UpLo, N, A, N, ipiv, work, &lwork, info );
    dsytri_( UpLo, N, A, N, ipiv, work, info);

    free(work); free(ipiv);
}
