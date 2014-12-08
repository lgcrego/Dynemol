#include <stdlib.h>
#include <stdio.h>

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

  #ifndef USE_GPU
    #include <time.h>
    static inline double get_wtime()
    {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC, &t);
        return( t.tv_sec +  t.tv_nsec/1.0e9 );
    };
  #else
    #define get_wtime magma_wtime
  #endif

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
  #warning "Debug messages enabled"
  #define DEBUG( A ) A
#else
  #warning "Debug messages disabled"
  #define DEBUG( A ) 
#endif

// Name mangling for use with fortran:
#ifdef FORTRAN
  #define GPU_Init      gpu_init_
  #define GPU_Finalize  gpu_finalize_
  #define GPU_dsygvd    gpu_dsygvd_
  #define GPU_dsygvd2s  gpu_dsygvd2s_
  #define GPU_dgemm     gpu_dgemm_
  #define GPU_syInvert  gpu_syinvert_
#endif


//----------------------------------------------------------------
// External functions
extern "C"
{
    int magma_get_parallel_numthreads();
    void dsygvd_(const int* itype, const char* jobZ, const char* UplO, const int* N, 
                 double* A, const int* ldA, double* B, const int* ldB, double* W, 
                 double* work, const int *lwork, int* iwork, const int* liwork, int *INFO);
    void dgemm_(const char *transA, const char *transB, const int *M, const int *N, const int *K, 
                const double *alpha, double *A, const int *ldA, double *B, const int *ldB, 
                const double *beta, double *C, const int *ldC);
    void dsytrf_( const char *UpLo, const int *M, double *A, const int *N, int *ipiv, double *work, const int *lwork, int *info );
    void dsytri_( const char *UpLo, const int *M, double *A, const int *N, int *ipiv, double *work, int *info);
}


//----------------------------------------------------------------
// Function prototypes:

// 1. Exported:
extern "C"
{
    void GPU_Init( const int *dev, const int *procs_per_dev );
    void GPU_Finalize( void );
    void GPU_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                     double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
    void GPU_dsygvd2s( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                       double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
    void GPU_dgemm( const char *transA, const char *transB, const int *M, const int *N, const int *K, 
                    const double *alpha, double *A, const int *LDA,
                    double *B, const int *LDB, const double *beta, double *C, const int *LDC );
    void GPU_syInvert( const char *UpLo, double *A, const int *N, int *info);
}

// 2. Internal:
//  - Helper
static inline int max( const int a, const int b ) { return( a>b ? a:b ); };
//  - Compute
static inline void cpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void cpu_syInvert(const char *UpLo, double *A, const int *N, int *info);
#ifdef USE_GPU
//  - Helper
static inline magma_vec_t lapack_to_magma_jobZ( const char jobZ );
static inline magma_uplo_t lapack_to_magma_uplo( const char uplo );
static inline cublasOperation_t lapack_to_cublas_trans( const char trans );
//  - Compute
static inline void gpu_dsygvd_2s( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                                  double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_dsygvd_m( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                                 double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int *N,
                               double *A, const int *ldA, double *B, const int *ldB, double *W, int *info );
static inline void gpu_syInvert(const char *UpLo, double *hA, const int *N, int *info);
#endif



/******************************
 *     Globlal variables
 ******************************/

// Device count and management
int devCount = 0;
int VirtualDevCount = 0;
int IuseGPU = 0; // false
int myGPU = -1;  // none

// CUBLAS handle
#ifdef USE_GPU
cublasHandle_t myHandle;
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
            printf("\n");
            printf("Using: (%i processes per GPU)\n",*procs_per_dev);
            magma_print_environment();
            printf("\n");
            fflush(stdout);
        }
        
        myGPU = (*pid / *procs_per_dev) % devCount;
        cudaSetDevice( myGPU );   // Bind to GPUs in a round-robin way according to pid
        cublasCreate( &myHandle );

        printf("Process nr. %i using GPU device nr. %i of %i(%i)\n", *pid, myGPU, devCount, VirtualDevCount);
        fflush(stdout);
    }
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
        cublasDestroy( myHandle );
        magma_finalize();
    }
#endif
}


//-----------------------------------------
// Flag conversions
// No parameter checking is performed!
#ifdef USE_GPU
static inline magma_vec_t lapack_to_magma_jobZ( const char jobZ )
    { return( jobZ == 'V' || jobZ == 'v' ? MagmaVec : MagmaNoVec ); }

static inline magma_uplo_t lapack_to_magma_uplo( const char UpLo )
    { return( UpLo == 'U' || UpLo == 'u' ? MagmaUpper : MagmaLower ); }

static inline cublasOperation_t lapack_to_cublas_trans( const char trans )
    { return( trans == 'N' || trans == 'n' ? CUBLAS_OP_N : (trans == 'T' || trans == 't' ? CUBLAS_OP_T : CUBLAS_OP_C) ); }
#endif



/******************************
 *     "Compute" functions
 ******************************/

//===================================================================
// Eigenvalues and eigenvectors - Ax = lambda Bx
//
// Alberto Torres
void GPU_dsygvd(const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
// Eigenvalues and eigenvectors
//
// Alberto Torres
void GPU_dsygvd2s( const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
static inline void gpu_dsygvd_2s( const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
    magma_dsygvdx_2stage( *itype, vec, MagmaRangeAll, uplo, n, A, *ldA, B, *ldB, 
                          0.0, 0.0, 1, n, &nZ, W, work, lwork, iwork, liwork, info );

    free_work(work); free_work(iwork);
#endif
}

//-------------------------------------------------------------------
// Eigenvalues and eigenvectors - on GPU (multi), called internally
//
// Alberto Torres
static inline void gpu_dsygvd_m( const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
static inline void gpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
static inline void cpu_dsygvd( const int *itype, const char *jobZ, const char *UpLo, const int *N,
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
// Matrix multiplication
//
// Alberto Torres
void GPU_dgemm( const char *transA, const char *transB, const int *M, const int *N, const int *K, 
                const double *alpha, double *hA, const int *LDA,
                double *hB, const int *LDB, const double *beta,
                double *hC, const int *LDC )
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
// Matrix inversion
//
// Alberto Torres
void GPU_syInvert( const char *UpLo, double *A, const int *N, int *info )
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
static inline void gpu_syInvert( const char *UpLo, double *hA, const int *N, int *info)
{
#ifdef USE_GPU
    const int n = *N;
    const int lwork = 64*n;
    const int lddA = ( (n + 31)/32 )*32;  // Make lddA multiple of 32 (faster)
    int *ipiv = (int*) malloc( n*sizeof(int) );
    double *dA, *dwork;
    
    magmablasSetKernelStream( NULL );

    magma_dmalloc( &dA, n*lddA );
    magma_dmalloc( &dwork, lwork );

    magma_dsetmatrix( n, n, hA, n, dA, lddA );

    // dsytrf/i are not impemented in magma-1.6.0
    // Remember to remove the #ifdef in Chebishev.f:Invertion_Matrix when changing to sy versions
    DEBUG( printf("gpu_syInvert: n= %i lwork= %i\n", n, lwork); fflush(stdout); )
    magma_dgetrf_gpu( n, n, dA, lddA, ipiv, info );
    magma_dgetri_gpu( n, dA, lddA, ipiv, dwork, lwork, info );
    
    magma_dgetmatrix( n, n, dA, lddA, hA, n );
    
    magma_free(dA); magma_free(dwork);
#endif
}

//-------------------------------------------------------------------
// Symmetric Matrix inversion - on CPU
//
// Alberto Torres
static inline void cpu_syInvert( const char *UpLo, double *A, const int *N, int *info )
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