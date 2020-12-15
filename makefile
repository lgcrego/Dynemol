.SUFFIXES: .f .F .for .cpp .F90 .cu .o

#make a - standard compilation
#make safe - compilation with safe features
#make debug - adds flag -g for debugging
#make serial - remove all parallelization flags
#make gdb - prepare code to GDB (equivalent to debug + serial) analysis
#make vtune - prepare code to intel-Vtune analysis


##########################
# FORTRAN CONFIGURATIONS #
##########################
# Compiler
FC = ifort

# Applied to all fortran files
FC_ALL = -xHost -ip -align

# Parallelization flags
FC_PARALLEL = -qopenmp -parallel

# Others flags for each file type *.f, *.F and *.F90
F_FLAGS =  -O3
F90_FLAGS = $(F_FLAGS)
f_FLAGS = -O2 -static $(FC_PARALLEL)


######################
# CPP CONFIGURATIONS #
######################
# Compiler
CC = icpc -std=c++11

# Applied to all cpp files
CC_ALL = -xHost -ip

# Parallelization flags
CC_PARALLEL = -qopenmp

# Others
CC_FLAGS = -O2 $(CC_PARALLEL)


#######################
# CUDA CONFIGURATIONS #
#######################
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
ifneq (,$(findstring USE_GPU, $(GPU_DEFS)))
# CUDA compiler
NVCC = nvcc
# compute capality (depends on your GPU, check!)
SM = 35
#SAFE_NVCC = -g -lineinfo
NVCCFLAGS = -O3 -gencode arch=compute_${SM},code=sm_${SM} -Xcompiler "-fno-strict-aliasing -march=native -fno-exceptions" $(SAFE_NVCC)
# -fno-strict-aliasing
#
# CUDA and MAGMA paths:
CUDADIR   = /usr/local/cuda
MAGMADIR  = /opt/magma
#
# CUDA and MAGMA libs:
# dynamic linking:
#LIB_CUDA  = -L$(CUDADIR)/lib64 -lcublas -lcusparse -lcudart
# static linking:
LIB_CUDA  = -L$(CUDADIR)/lib64 -lcublas_static -lcusparse_static -lculibos -lcudart_static -ldl
LIB_MAGMA = $(MAGMADIR)/lib/libmagma.a
#
LIB_GPU   = $(LIB_MAGMA) $(LIB_CUDA) -lstdc++
INCS_GPU  = -I$(CUDADIR)/include -I$(MAGMADIR)/include
endif


############################
# LIBARRIES CONFIGURATIONS #
############################
ifndef MKLROOT
$(error 'MKLROOT not set')
endif
# If MKLROOT is not defined, add it's path here
# MKLROOT =
LIB_BLAS = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP = -liomp5 -lpthread

INCLUDES_MKL = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include/fftw


LIBS  = $(LIB_GPU) $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP) -lrt
INCLUDES = $(INCLUDES_MKL)


#########################
# FILES AND DEPENDECIES #
#########################
#INCS1   = comun.inc integcoul.inc m2cdat.inc $(INCLUDES)

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
                  checklist.o \
		  allocation_m.o \
		  util.o \
		  EHT_input.o \
		  tuning.o \
                  IdentifyNonBonded.o \
		  babel_routines.o \
		  babel.o \
		  gmx2mdflex.o \
		  namd2mdflex.o \
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
		  LCMO_Builder.o \
		  Coulomb.o \
		  DP_main.o \
		  td_dp.o \
		  DP_FMO.o \
		  dipole_phi.o \
                  EnvField.o \
		  polarizability.o \
		  hamiltonians.o \
		  QCModel_Huckel.o \
		  HuckelForces.o \
		  Ehrenfest.o \
		  CoulInt_QMMM.o \
		  FMO.o \
		  electron_hole_DP.o \
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
		  backup.o \
		  auto_correlation.o \
		  ElHl_schroedinger.o \
		  diagnostic.o \
		  qdynamics.o \
                  Taylor.o \
		  ElHl_Chebyshev.o \
		  AO_adiabatic.o \
		  Chebyshev_driver.o \
		  eigen_driver.o \
		  ga_driver.o \
		  avrg_confgs.o \
 		  main.o

SOURCE_GPU = GPU_Interface.o \
             Taylor_gpu.o
#             Chebyshev_gpu.o

ifneq (,$(findstring USE_GPU,$(GPU_DEFS)))
SOURCE_CUDA= Chebyshev_gpu_kernels.o \
             dzgemv_kernels.o
endif


#########
# RULES #
#########
a: $(SOURCE1) $(SOURCE2) $(SOURCE_GPU) $(SOURCE_CUDA)
	rm -f a
	$(FC) $(FC_ALL) $(INCLUDES) -o a $(SOURCE1) $(SOURCE2) $(SOURCE_GPU) $(SOURCE_CUDA) $(LIBS)
	-rm -f *.log

# Program runs very slowly with this
safe: FC_ALL += -check all -traceback -fstack-protector -assume protect_parens
safe: CC_ALL += -traceback -fstack-protector
safe: a

# Just adds debug flag to everything
debug: FC_ALL += -g
debug: CC_ALL += -g
debug: a

# Removes parallel flags
serial: FC_PARALLEL =
serial: CC_PARALLEL =
serial: a

# Easiert do debug when there is no threads around
gdb: F_FLAGS = -O0
gdb: f_FLAGS = -O0
gdb: CC_FLAGS = -O0
gdb: debug

# Adds lots of flags and remove static from f_FLAGS
# Flags taken from here:
# https://software.intel.com/en-us/vtune-amplifier-help-compiler-switches-for-performance-analysis-on-linux-targets
vtune: FC_ALL += -debug inline-debug-info -D_DEBUG -qopenmp-link dynamic -parallel-source-info=2
vtune: CC_ALL += -debug inline-debug-info -D_DEBUG -qopenmp-link dynamic -parallel-source-info=2
vtune: f_FLAGS = -O2 $(FC_PARALLEL)
vtune: debug

.F.o:
	$(FC) -fpp $(FC_ALL) $(F_FLAGS) $(INCLUDES) -c $<

.f.o:
	$(FC) -fpp -free $(FC_ALL) $(f_FLAGS) $(INCLUDES) $(GPU_DEFS) -c $<

.F90.o:
	$(FC) -fpp $(FC_ALL) $(F_FLAGS) $(INCLUDES) $(GPU_DEFS) -c $<

.cpp.o:
	$(CC) $(CC_ALL) -align -fno-exceptions -restrict $(CC_FLAGS) $(INCS_GPU) $(GPU_DEFS) -c $<

.cu.o:
	$(NVCC) $(NVCCFLAGS) $(INCS_GPU) $(GPU_DEFS) -c $<


clean:
	-rm -fv a *.o *.mod *__genmod.f90 *.i

depend:
	@echo -en "Searching module dependencies..."
	@chmod +x ./makedepend.bsh
	@./makedepend.bsh > dependencies.txt
	@echo -en " done.\n"


# Dependency list:
-include dependencies.txt

