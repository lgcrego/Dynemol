
#ifdef GPU_SYGVD2S_VER
#  warning "SYGVD calls XPU_dsygvd2s"
#  define SYGVD( A , B , C , D , E , F , G )      XPU_dsygvd2s( D, E, F, size(A,2), A, size(A,1), B, size(B,1), C, G )
#else
#  warning "SYGVD calls XPU_dsygvd"
#  define SYGVD( A , B , C , D , E , F , G )      XPU_dsygvd( D, E, F, size(A,2), A, size(A,1), B, size(B,1), C, G )
#endif

#define DGEMM( A, B, C, D, E, F, G, H, I, J, K, L, M )    XPU_dgemm( A, B, C, D, E, F, G, H, I, J, K, L, M )
