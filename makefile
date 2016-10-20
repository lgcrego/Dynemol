#
.SUFFIXES:
.SUFFIXES: .f .F .for .cpp .F90 .cu .o 

FC=ifort -xHost -ip -fpp
FREE = -free

# use this flag for debugging and coding up
SAFE = #-check all -traceback -fstack-protector -assume protect_parens -implicitnone -warn all 

FFLAGS1 = -O3 -align #array64byte
FFLAGS2 = -O2 -align -openmp -parallel $(FREE) $(SAFE) -static #array64byte 

CXX = icpc -std=c++11
CFLAGS = -O2 -align -xHost -ip -openmp -fno-exceptions -restrict 

# MKLROOT  = If MKLROOT is not defined in your environment, edit and uncomment this line
LIB_BLAS   = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP    = -liomp5 -lpthread
INCS_MKL   = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include/fftw

# Uncomment the lines below when compiling for GPUs
# GPU_DEFS options:
#   -DGPU_TIMING       : Print timings (CPU/GPU)
#   -DUSE_GPU          : Compile with GPU support
#   -DGPU_DEBUG        : Print some debug messages
#   -DGPU_SYGVDM_VER   : Use multi-gpu version of SYGVD (faster than single-gpu even with 1 gpu)
#   -DGPU_SYGVD2S_VER  : Use two stage version of SYGVD (faster, but needs more memory)
#   -DGPU_DONT_PIN_MEM : Don't use pinned memory for faster transfers (in Fortran code)
#   -DGPU_PIN_MEM_WORK : Use pinned memory for work spaces (in C code)
#GPU_DEFS  = -DUSE_GPU
#
#NVCC      = nvcc
#NVCCFLAGS = -O3 -arch=sm_35 -Xcompiler "-fno-strict-aliasing -march=native -fno-exceptions"
#
# CUDA and MAGMA paths:
#CUDADIR   = /usr/local/cuda
#MAGMADIR  = /opt/magma-2.1.0
#
# CUDA and MAGMA libs:
#LIB_CUDA  = -L$(CUDADIR)/lib64 -lcublas -lcusparse -lcudart
#LIB_MAGMA = $(MAGMADIR)/lib/libmagma.a
#
#LIB_GPU   = $(LIB_MAGMA) $(LIB_CUDA) -lstdc++
#INCS_GPU  = -I$(CUDADIR)/include -I$(MAGMADIR)/include

LIB  = $(LIB_GPU) $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP) -lrt
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
		  Matrix_math.o \
		  exec_time.o \
		  types_EHT.o \
		  types_MM.o \
		  OPT_parent.o \
		  parameters.o \
		  parameters_MM.o \
		  allocation_m.o \
		  util.o \
		  EHT_input.o \
		  tuning.o \
          IdentifyNonBonded.o \
		  babel_routines.o \
		  babel.o \
		  gmx2mdflex.o \
		  structure.o \
		  md_read.o	\
		  md_setup.o \
		  f_intra.o \
		  f_inter.o \
		  md_output.o \
		  pbc.o \
		  overlap_D.o \
		  Ehrenfest.o \
		  HuckelForces.o \
		  STO.o \
		  multip_routines.o \
		  electron_hole_DP.o \
		  LCMO_Builder.o \
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
          backup_MM.o \
		  Berendsen.o \
		  NoseHoover.o \
		  NoseHoover_Reversible.o \
          NVE.o \
		  VDOS_m.o \
		  MM_dynamics.o \
		  MM_driver.o \
		  film_STO.o \
		  DOS_m.o \
		  oscillator.o \
		  ga_QCModel.o \
		  cost_tuning_EH.o \
		  cost_tuning_MM.o \
		  nonlinearCG.o \
		  CG_class.o \
		  MM_ERG_class.o \
		  nonlinear-sidekick.o \
		  FF_OPT_class.o \
		  CG_EH_driver.o \
		  ga_routines.o \
		  CG_MM_driver.o \
		  vibes_driver.o \
		  solvated_M.o \
		  DOS_tool.o \
		  backup.o \
		  auto_correlation.o \
		  ElHl_schroedinger.o \
		  diagnostic.o \
		  qdynamics.o \
		  Chebyshev.o \
		  ElHl_Chebyshev.o \
		  AO_adiabatic.o \
		  ElHl_adiabatic.o \
		  Chebyshev_driver.o \
		  eigen_driver.o \
		  ga_driver.o \
		  avrg_confgs.o \
 		  main.o

SOURCE_GPU = GPU_Interface.o \
             Chebyshev_gpu.o

ifneq (,$(findstring USE_GPU,$(GPU_DEFS)))
SOURCE_CUDA= Chebyshev_gpu_kernels.o \
             dzgemv_kernels.o
endif


a: $(SOURCE1) $(SOURCE2) $(SOURCE_GPU) $(SOURCE_CUDA)
	rm -f a
	$(FC) $(INCS) -o a $(SOURCE1) $(SOURCE2) $(SOURCE_GPU) $(SOURCE_CUDA) $(LIB) 
	-rm -f *.log


.F.o:
	$(FC) $(FFLAGS1) $(INCS) -c $<

.f.o:
	$(FC) $(FFLAGS2) $(INCS) $(GPU_DEFS) -c $<

.F90.o:
	$(FC) $(FFLAGS1) $(INCS) $(GPU_DEFS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) $(INCS_GPU) $(GPU_DEFS) -c $<

.cu.o:
	$(NVCC) $(NVCCFLAGS) $(INCS_GPU) $(GPU_DEFS) -c $<


clean: 
	-rm -f a *.o *.mod; touch *.f

depend:
	@echo -en "Searching module dependencies..."
	@chmod +x ./makedepend.bsh
	@./makedepend.bsh > dependencies.txt
	@echo -en " done.\n"


## Dependency list:
-include dependencies.txt

