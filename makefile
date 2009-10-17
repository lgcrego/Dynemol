#

SUFFIX=.f

FC=/usr/opt/intel/Compiler/11.1/046/bin/intel64/ifort 
FREE = -free
#FC=gfortran
#FREE = -ffree-form 

FFLAGS = -O3 $(FREE) 

#LIB    = -L/usr/lib64 -llapack -lblas -L/home/lrego/lib -lmyblas95
#INCS   = -I/home/lrego/lib/blas_95

LIB_BLAS   = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lmkl_blas95 -lmkl_em64t -lguide -lpthread
LIB_LAPACK = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lmkl_lapack95 -lmkl_em64t -lguide -lpthread
INCS_MKL   = -I/opt/intel/mkl/10.0.5.025/lib/em64t/

LIB  = $(LIB_BLAS) $(LIB_LAPACK)
INCS = $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

SOURCE = constants_m.o type_m.o allocation_m.o util.o babel.o EHT_input.o structure.o \
		 pbc.o overlap_D.o STO.o QCModel_Huckel.o projectors.o FMO.o data_output.o film_STO.o  \
		 DOS_m.o multip_core.o oscillator.o QOptics.o dynamics.o rk4.o solvated_M.o main.o

a: $(SOURCE)  
	rm -f a
	$(FC) $(INCS) -o a $(SOURCE) $(LIB) 
	-rm -f *.log
#-rm -f *.o *.mod; touch *.f
.f.o:
	$(FC) $(FFLAGS) $(INCS) -c $*$(SUFFIX)
clean: 
	-rm -f *.o *.mod; touch *.f 
