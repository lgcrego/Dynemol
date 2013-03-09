#
.SUFFIXES: .f .F .for

FC=/opt/intel/bin/ifort 
FREE = -free
#FC=gfortran
#FREE = -ffree-form 

FFLAGS1 = -O3 
FFLAGS2 = -O1 -openmp -parallel $(FREE) -static

LIB_BLAS   = -L/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIB_LAPACK = -L/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIB_OMP    = -L/opt/intel/lib/intel64 -liomp5 -lmatmul
INCS_MKL   = -I/opt/intel/composer_xe_2011_sp1.9.293/mkl/include/intel64/lp64 

LIB  = $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP)
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
		  types_MM.o \
		  parameters.o  \
		  allocation_m.o  \
		  util.o  \
		  EHT_input.o  \
		  tuning.o \
		  babel_routines.o  \
		  babel.o  \
		  structure.o  \
		  md_read.o	\
		  f_intra.o \
		  f_inter.o \
		  md_setup.o \
		  md_output.o \
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
		  polarizability.o \
		  QCModel_Huckel.o \
		  QCModel_ElHl.o \
		  data_output.o \
		  md_dynamics.o \
		  MM_dynamics.o \
		  film_STO.o \
		  DOS_m.o \
		  oscillator.o \
		  ga_QCModel.o \
		  cost_tuning.o \
		  CG_class.o \
		  nonlinearCG.o \
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
