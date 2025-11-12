# Usage needs sensible 
# definition of CSIRO_OSTYPE, e.g. I use
# % uname -smp | /usr/bin/tr ' ' '-' |  sed 's/\(.*\)-\([^-]*\)/\1-unknown/g'`
# yielding, say, Linux-x86_64-unknown

# If moving/porting -> note full package names in javac -h  calls!!

# NB if using JAVA_LIBRARY_PATH methods .....
# set LD_LIBRARY_PATH to include ./${CSIRO_OSTYPE}/lib
# New implementation loads library by a full path, so 
# JAVA_LIBRARY_PATH should not be needed

# Platform dependent things in this include....
# i.e. provide one such file for each desired OS
include Makefile.${CSIRO_OSTYPE}

SEISMIC_INVERTER_HOME=../../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar


COMPILER = javac
JAVAOPTIONS = -g -O -deprecation
COPTIONS = -O3 -g ${M32}

CLASSES = FFTW.class \

.SUFFIXES: .java .class 

all: FFTW

FFTW: FFTW.java 
	#javac  -classpath ${CLASSPATH} ${JAVAOPTIONS} FFTW.java
	# Note full package name!! for header
	javac -h . -classpath ${CLASSPATH} FFTW.java
	# Note here that having the .c files *first* is important for linking
	cc ${COPTIONS} $(CFLAGS_UNIVERSAL) ${CSHAREDLIBFLAG} ${INCLUDEPATH} FFTWImp.c -o ./${CSIRO_OSTYPE}/libFFTW.${DYLIBSUFFIX}

# on OSX, may require (leopard) (google away)
# setenv MACOSX_DEPLOYMENT_TARGET 10.5
# COPTIONS += -no_compact_linkedit

# Remove all class files from the project directory. 
clean: 
	rm *.class 

# Implicit rules
.java.class:
	$(COMPILER) -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?

