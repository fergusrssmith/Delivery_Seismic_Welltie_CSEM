# makefile called recursively by ../../Makefile

# Your environment needs to have defined CSIRO_OSTYPE, as e.g.
# export CSIRO_OSTYPE=`uname -smp | sed 's/ /-/g;s/-[^-]*$/-unknown/g'`
# Typically a platform-specific OS-identifying string: e.g. $CSIRO_OSTYPE = Linux-i686-unknown

# Platform & compiler-specific things enter via this include....
# i.e. provide one such file for each desired OS
include ../../Makefile.${CSIRO_OSTYPE}

# allow preprocessing (gfortran)
FORFLAGS+=-x f95-cpp-input

# compile as optimized, or debug....
ifeq (debug,$(findstring debug,$(MAKECMDGOALS)))
  FFLAGS=-g  $(FORFLAGS) 
else
  FFLAGS=-w -O3 $(FORFLAGS)
endif

LIBDIR = ../../${CSIRO_OSTYPE}

# New scripps - order of f90 objects probably important here !
#F_OBJS = Occam.o JavaModules.o FilterModules.o Dipole1D_ptrs.o JavaInterface.o CallDipole1D.o  OccamModules.o Occam.o CSEM1D.o  
F_OBJS = JavaModules.o FilterModules.o Dipole1D_ptrs.o JavaInterface.o CallDipole1D.o  OccamModules.o CSEM1D.o  

all: fortran_codes 
.PHONY: all fortran_codes debug fast

fortran_codes: $(addprefix $(LIBDIR)/, $(F_OBJS))
#   OSX prefers libtool, but the below will do...
#	libtool -static -s -o $(LIBDIR)/Scripps.a $(addprefix $(LIBDIR)/, $(F_OBJS))
	ar -r $(LIBDIR)/Scripps.a $(addprefix $(LIBDIR)/, $(F_OBJS))
	$(RANLIB) $(LIBDIR)/Scripps.a 

debug: all
fast: all

# General cleanup

clean: 
	rm *.o *.mod $(LIBDIR)/*.o $(LIBDIR)/Scripps.a

$(LIBDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

