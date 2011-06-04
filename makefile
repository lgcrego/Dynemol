#
.SUFFIXES: .f .F .for

FC=/opt/intel/Compiler/11.1/064/bin/intel64/ifort 
FREE = -free
#FC=gfortran
#FREE = -ffree-form 

FFLAGS1 = -O3 
FFLAGS2 = -O3 -openmp -parallel $(FREE) 

#LIB    = -L/usr/lib64 -llapack -lblas -L/home/lrego/lib -lmyblas95
#INCS   = -I/home/lrego/lib/blas_95

LIB_BLAS   = -L/opt/intel/mkl/10.2.4.032/lib/em64t -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LIB_LAPACK = -L/opt/intel/mkl/10.2.4.032/lib/em64t -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#LIB_BLAS   = -L/opt/intel/mkl/10.2.4.032/lib/em64t -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lguide -lpthread
#LIB_LAPACK = -L/opt/intel/mkl/10.2.4.032/lib/em64t -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lguide -lpthread
INCS_MKL   = -I/opt/intel/mkl/10.2.4.032/include/em64t/lp64

LIB  = $(LIB_BLAS) $(LIB_LAPACK)
INCS = $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

#INCS1   = comun.inc integcoul.inc m2cdat.inc $(INCS)

SOURCE1 = integ-Coul.o \
	 	  Coul0sim.o \
		  m2caux3-Coul.o \
		  abcpes-Coul.o \
		  ckplm-Coul.o \
		  util-Coul.o

SOURCE2 = constants_m.o  \
		  exec_time.o  \
		  type_m.o  \
		  parameters.o  \
		  allocation_m.o  \
		  util.o  \
		  EHT_input.o  \
		  tuning.o \
		  babel_routines.o  \
		  babel.o  \
		  structure.o  \
		  pbc.o  \
		  overlap_D.o  \
		  STO.o \
		  multip_routines.o  \
		  electron_hole_DP.o	\
		  FMO.o \
		  DP_main.o \
		  td_dp.o \
		  DP_FMO.o \
		  dipole_phi.o \
		  Coulomb.o \
		  QCModel_Huckel.o \
		  data_output.o \
		  film_STO.o \
		  DOS_m.o \
		  oscillator.o \
		  QOptics.o \
		  ga_QCModel.o \
		  ga_routines.o \
		  solvated_M.o \
		  schroedinger.o \
		  rk4.o \
		  diagnostic.o \
		  qdynamics.o \
		  Chebyshev.o \
		  backup.o \
		  AO_adiabatic.o \
		  MO0_adiabatic.o \
		  MOt_adiabatic.o \
		  eigen_driver.o \
		  Chebyshev_driver.o \
		  ga_driver.o \
		  avrg_confgs.o \
 		  main.o


a: $(SOURCE1) $(SOURCE2)
	rm -f a
	$(FC) $(INCS) -o a $(SOURCE1) $(SOURCE2) $(LIB) 
	-rm -f *.log

.F.o:
	$(FC) $(FFLAGS1) $(INCS) -c $*.F

.f.o:
	$(FC) $(FFLAGS2) $(INCS) -c $*.f
 
clean: 
	-rm -f *.o *.mod; touch *.f
