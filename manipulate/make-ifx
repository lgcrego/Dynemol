#

SUFFIX=.f

FC=ifx
SAFE = #-g -check all,nouninit -traceback -fstack-protector -assume protect_parens -implicitnone -warn all,noexternal -fpe-all=0
FFLAGS = -xhost -qopenmp -fpp -O3 -free $(SAFE)

LIB    =  -liomp5
INCS   = 

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

SOURCE = types_m.o \
 		 constants_m.o \
		 diagnosis.o \
		 parms.o \
		 tuning.o \
		 topol.o \
         amber.o \
		 RW.o \
		 gmx.o \
		 EDT.o \
		 crystal.o \
		 solvent.o \
		 elastic.o \
		 correlation.o \
		 stat.o \
         RW_driver.o \
		 trj.o \
		 multiple_trj.o \
		 occupation_bin.o \
         aminoacids.o \
		 edview.o

manipulate: $(SOURCE)  
	rm -f manipulate
	$(FC) $(INCS) -o manipulate $(SOURCE) $(LIB) 
	-rm -f *.log
.f.o:
	$(FC) $(FFLAGS) $(INCS) -c $*$(SUFFIX)
clean: 
	-rm -f *.o *.mod; touch *.f 
