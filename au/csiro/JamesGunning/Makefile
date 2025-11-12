# Master make file for seismic inverter codes and the requisite libraries
# To install, type "make"
#
#
# standard make line
MAKE = make

###############################################################################
#
#	Call 'make' in all sub-directories which have a 'Makefile'.
#
#
###############################################################################

# DIRS		= lib/Numerical lib/DifferentialEvolution lib/Minpack/corejava lib/Minpack/linear_algebra lib/Minpack/optimization lib sulib seismicinverter
DIRS = \
./seismicinverter \
./sulib \
./waveletextraction \
./lib \
./lib/graphutil \
./lib/DifferentialEvolution \
./lib/util \
./lib/stat/clustering \
./lib/stat/segmentation \
./lib/stat/util \
./lib/Minpack/linear_algebra \
./lib/Minpack/optimization \
./lib/Minpack/corejava \
./lib/Numerical \
./lib/Numerical/djbFFT \
./CSEM/OneDinversion \
./CSEM/plotting


#
#	MACROS
#

MAKE_TARGETS	= all clean

.PHONY: $(MAKE_TARGETS)


#
#	TARGETS
#

$(MAKE_TARGETS):
	@if [ -z "$(DIRS)" ]; \
        then \
		dirs=`echo *`; \
	else \
		dirs="$(DIRS)"; \
	fi; \
        for dir in $$dirs; \
	do \
		if [ -f $$dir/Makefile ]; \
		then \
			echo ""; \
			echo "Directory $$dir: make $@"; \
			( cd $$dir; make $@ ); \
		fi; \
	done


make:
	@if [ -z "$(DIRS)" ]; \
        then \
		dirs=`echo *`; \
	else \
		dirs="$(DIRS)"; \
	fi; \
        for dir in $$dirs; \
	do \
		if [ -f $$dir/Makefile ]; \
		then \
			echo ""; \
			echo "Directory $$dir: make "; \
			( cd $$dir; make ); \
		fi;\
	done


world:
	@make make; make all
