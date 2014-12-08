#
.SUFFIXES: .f .F .for .cpp .F90 

FC=ifort -xHost -ip -fpp
FREE = -free
#FC=gfortran
#FREE = -ffree-form 

FFLAGS1 = -O3
FFLAGS2 = -O1 -openmp -parallel $(FREE) -static 

CXX = icpc -std=c++11
CFLAGS = -O2 -xHost -ip -fno-exceptions  

LIB_BLAS   = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP    = -liomp5 -lmatmul -lpthread
INCS_MKL   = -I$(MKLROOT)/include/intel64/lp64

GPU_DEFS  = -DFORTRAN -DGPU_TIMING #-DUSE_GPU -DGPU_SYGVDM_VER 
#-DGPU_DEBUG -DGPU_SYGVDM_VER -DGPU_SYGVD2S_VER #-DGPU_PIN_MEM -DGPU_PIN_MEM_WORK
CUDADIR   = /usr/local/cuda
LIB_CUDA  = -L$(CUDADIR)/lib64 -lcublas -lcudart
MAGMADIR  = /opt/magma
LIB_MAGMA = $(MAGMADIR)/lib/libmagma.a
#LIB_GPU   = $(LIB_MAGMA) $(LIB_CUDA) -lstdc++
INCS_GPU  = -I$(CUDADIR)/include -I$(MAGMADIR)/include

LIB  = $(LIB_GPU) $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP)
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

SOURCE2 = constants_m.o \
		  exec_time.o \
		  types_EHT.o \
		  types_MM.o \
		  parameters.o \
		  parameters_MM.o \
		  allocation_m.o \
		  util.o \
		  EHT_input.o \
		  tuning.o \
          IdentifyNonBonded.o \
		  gmx2mdflex.o \
		  babel_routines.o \
		  babel.o \
		  structure.o \
		  md_read.o	\
		  md_setup.o \
		  f_intra.o \
		  f_inter.o \
		  md_output.o \
		  pbc.o \
		  overlap_D.o \
		  STO.o \
		  multip_routines.o \
		  electron_hole_DP.o \
		  LCMO_Builder.o  \
		  FMO.o \
		  DP_main.o \
		  td_dp.o \
		  DP_FMO.o \
		  dipole_phi.o \
		  Coulomb.o \
		  polarizability.o \
		  CoulInt_QMMM.o \
		  QCModel_Huckel.o \
		  QCModel_ElHl.o \
		  AlphaPolar.o \
		  data_output.o \
		  verlet.o \
          backup_MM.o \
		  MM_dynamics.o \
		  MM_driver.o \
		  film_STO.o \
		  DOS_m.o \
		  oscillator.o \
		  ga_QCModel.o \
		  cost_tuning.o \
		  CG_class.o \
		  MM_ERG_class.o \
		  nonlinearMM.o \
		  nonlinearCG.o \
		  vibes_driver.o \
		  CG_routines.o \
		  ga_routines.o \
		  solvated_M.o \
		  DOS_tool.o \
		  backup.o \
		  ElHl_schroedinger.o \
		  diagnostic.o \
		  qdynamics.o \
		  ElHl_Chebyshev.o \
		  Chebyshev.o \
		  AO_adiabatic.o \
		  ElHl_adiabatic.o \
		  eigen_driver.o \
		  Chebyshev_driver.o \
		  ga_driver.o \
		  avrg_confgs.o \
 		  main.o

SOURCE_GPU = GPU_Interface.o \
             cuda_runtime.o


a: $(SOURCE1) $(SOURCE2) $(SOURCE_GPU)
	rm -f a
	$(FC) $(INCS) -o a $(SOURCE1) $(SOURCE2) $(SOURCE_GPU) $(LIB) 
	-rm -f *.log

.F.o:
	$(FC) $(FFLAGS1) $(INCS) -c $*.F

.f.o:
	$(FC) $(FFLAGS2) $(INCS) $(GPU_DEFS) -c $*.f

.F90.o:
	$(FC) $(FFLAGS1) $(INCS) $(GPU_DEFS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) $(INCS_GPU) $(GPU_DEFS) -c $<


clean: 
	-rm -f *.o *.mod; touch *.f


