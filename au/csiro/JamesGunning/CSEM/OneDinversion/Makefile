# specify either 
# make 
# make debug

# Your environment needs to have defined CSIRO_OSTYPE, as e.g.
# export CSIRO_OSTYPE=`uname -smp | sed 's/ /-/g;s/-[^-]*$/-unknown/g'`
# Typically a platform-specific OS-identifying string: e.g. $CSIRO_OSTYPE = Linux-i686-unknown

# Insert here the location of the top of the delivery install tree: (e.g. /home/username/java)
SEISMIC_INVERTER_HOME = ../../../../..
CLASSPATH=${SEISMIC_INVERTER_HOME}:${SEISMIC_INVERTER_HOME}/jars/colt.jar

# Platform & compiler-specific things enter via this include....
# i.e. provide one such file for each desired OS
include Makefile.${CSIRO_OSTYPE}
#include Makefile.cygwin

# compile as optimized, or debug....

ifeq (debug,$(findstring debug,$(MAKECMDGOALS)))
  FDBGFLAGS= -g
  FFLAGS=-g  $(FORFLAGS)
  COPTIONS = -g 
  JAVAOPTIONS = -g -deprecation
else
  FDBGFLAGS=
  FFLAGS=-w -O3 $(FORFLAGS)
  COPTIONS = -O3 
  JAVAOPTIONS = -O -deprecation
endif

# subdirectories in which Make is invoked
F_LIBRARIES := NativeSources/LegacyScripps NativeSources/Scripps

LIBDIR = ./${CSIRO_OSTYPE}

.SUFFIXES: .java .class

all: interface javaclasses

debug: all
fast: all

.PHONY: all interface $(F_LIBRARIES) javaclasses debug fast

interface: $(F_LIBRARIES) ScrippsLegacyModelInterface.class \
		ScrippsModelInterface.class \
		$(LIBDIR)/ScrippsLegacyModelInterfaceLib.$(JNILIBSUFFIX) \
		$(LIBDIR)/ScrippsModelInterfaceLib.$(JNILIBSUFFIX) 

javaclasses: CSEM1DInverter.class VersionNumber.class DataQCchecker.class \
	CSEM1DBootstrapChecker.class ModelStudyChecker.class DataSetMaker.class ModePairConnector.class

debug: all
fast: all

# Legacy scripps native codes
C_LEGACY_OBJS = ScrippsLegacyModelInterfaceImp.o

# New scripps 
COBJS = ScrippsLegacyModelInterfaceImp.o

# Java parts
CSEM1DInverter.class: CSEM1DInverter.java
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) CSEM1DInverter.java

# VersionNumber is special -script writes the file....
VersionNumber.class: VersionNumber.java
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) VersionNumber.java

# Note full package name!! for header in javac -h  call
ScrippsLegacyModelInterface.class: ScrippsLegacyModelInterface.java
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) ScrippsLegacyModelInterface.java
	javac -h . -classpath ${CLASSPATH} ScrippsLegacyModelInterface.java
# Note here that having the .c files *first* is important for linking
# must have this set for things to work.....
$(LIBDIR)/ScrippsLegacyModelInterfaceLib.$(JNILIBSUFFIX): platform_specifics.h ScrippsLegacyModelInterfaceImp.c $(F_LEGACY_OBJS)
#   Note here gfortran can compile the single ScrippsLegacyModelInterfaceImp.c and link everything to make the shared library
#   => the .c compile below is skipped.
#	cc -c $(CFLAGS) $(COPTIONS) $(CSHAREDLIBFLAG) $(INCLUDEPATH) ScrippsLegacyModelInterfaceImp.c -o ScrippsLegacyModelInterfaceImp.o
	$(FLD) $(FDBGFLAGS) $(FLDFLAGS) $(INCLUDEPATH) ScrippsLegacyModelInterfaceImp.c $(LIBDIR)/LegacyScripps.a \
		-o $(LIBDIR)/libScrippsLegacyModelInterface.$(JNILIBSUFFIX)
#   MacOSX guides suggest making .dylib paths absolute, so this might help if needed...
#		-o $(CURDIR)/$(LIBDIR)/ScrippsLegacyModelInterfaceLib.$(JNILIBSUFFIX)
#	echo "setenv LD_LIBRARY_PATH to include this directory! " 


ScrippsModelInterface.class: ScrippsModelInterface.java
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) ScrippsModelInterface.java
	javac -h . -classpath ${CLASSPATH} ScrippsModelInterface.java


$(LIBDIR)/ScrippsModelInterfaceLib.$(JNILIBSUFFIX): platform_specifics.h ScrippsModelInterfaceImp.c $(addprefix $(LIBDIR)/, $(FOBJS))
#	cc -c $(CFLAGS) $(COPTIONS) $(CSHAREDLIBFLAG) $(INCLUDEPATH) ScrippsModelInterfaceImp.c -o ScrippsModelInterfaceImp.o
	$(FLD) $(FDBGFLAGS) $(FLDFLAGS) $(INCLUDEPATH) ScrippsModelInterfaceImp.c $(LIBDIR)/Scripps.a \
		-o $(LIBDIR)/libScrippsModelInterface.$(JNILIBSUFFIX)


# General cleanup

clean: 
	rm *.o $(LIBDIR)/*Scripps* $(LIBDIR)/*.o *.class

# Implicit rules
$(interface) $(F_LIBRARIES):
	$(MAKE) --directory=$@ $(MAKECMDGOALS)

.java.class:
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?

