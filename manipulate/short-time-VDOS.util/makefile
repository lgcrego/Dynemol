#

SUFFIX= .f 

FC=ifort
SAFE = -g -traceback
FFLAGS = -xhost -qopenmp -fpp -O3 -free -qmkl $(SAFE)

# MKLROOT  = If MKLROOT is not defined in your environment, edit and uncomment this line
LIB_BLAS   = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP    = -liomp5 -lpthread
INCS_MKL   = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include/fftw


LIB    =  -liomp5 -L$(MKLROOT)/lib/intel64/ 
INCS   = $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

SOURCE = short-time-FFT.o

short-time-FFT: $(SOURCE)  
	rm -f short-time-FFT
	$(FC) $(INCS) $(FFLAGS) -o short-time-FFT $(SOURCE) $(LIB) 
	-rm -f *.log
.f.o:
	$(FC) $(FFLAGS) $(INCS) -c $*$(SUFFIX)
clean: 
	-rm -f *.o *.mod short-time-FFT; touch *.f 
