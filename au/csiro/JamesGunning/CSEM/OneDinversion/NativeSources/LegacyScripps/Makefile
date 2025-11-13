# specify either 
# make 
# make debug

# Your environment needs to have defined CSIRO_OSTYPE, as e.g.
# export CSIRO_OSTYPE=`uname -smp | sed 's/ /-/g;s/-[^-]*$/-unknown/g'`
# Typically a platform-specific OS-identifying string: e.g. $CSIRO_OSTYPE = Linux-i686-unknown

# Platform & compiler-specific things enter via this include....
# i.e. provide one such file for each desired OS
include ../../Makefile.${CSIRO_OSTYPE}

# compile as optimized, or debug....
ifeq (debug,$(findstring debug,$(MAKECMDGOALS)))
  FFLAGS=-g  $(FORFLAGS) -fallow-argument-mismatch
else
  FFLAGS=-w -O3 $(FORFLAGS) -fallow-argument-mismatch
endif

LIBDIR = ../../${CSIRO_OSTYPE}

# Legacy scripps native codes
F_OBJS = modules.o hedmods.o hy2fht.o seafloor_JG.o

all: fortran_codes 
.PHONY: all fortran_codes debug fast

fortran_codes: $(addprefix $(LIBDIR)/, $(F_OBJS))
#   ON OSX
#	libtool -static -s -o $(LIBDIR)/LegacyScripps.a $(addprefix $(LIBDIR)/, $(F_OBJS))
	ar -r $(LIBDIR)/LegacyScripps.a $(addprefix $(LIBDIR)/, $(F_OBJS))
	$(RANLIB) $(LIBDIR)/LegacyScripps.a 

debug: all
fast: all

# General cleanup

clean: 
	rm *.o $(LIBDIR)/*.o

$(LIBDIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(LIBDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

