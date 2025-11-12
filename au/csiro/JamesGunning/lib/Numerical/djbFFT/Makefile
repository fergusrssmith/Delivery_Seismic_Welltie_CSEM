# Usage needs sensible 
# definition of CSIRO_OSTYPE, e.g. I use
# % uname -smp | /usr/bin/tr ' ' '-' |  sed 's/\(.*\)-\([^-]*\)/\1-unknown/g'`
# yielding, say, Linux-x86_64-unknown

# If moving/porting -> note full package names in javac -h  calls!!

# NB if using JAVA_LIBRARY_PATH methods .....
# set LD_LIBRARY_PATH to include ./djbfft-0.76/${CSIRO_OSTYPE}/lib
# New implementation loads library by a full path, so 
# JAVA_LIBRARY_PATH should not be needed

# Platform dependent things in this include....
# i.e. provide one such file for each desired OS
include Makefile.${CSIRO_OSTYPE}

SEISMIC_INVERTER_HOME=../../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar


COMPILER = javac
JAVAOPTIONS = -g -O -deprecation     
COPTIONS = -O2 -g ${M32} 
CC = cc

# CFLAGS_UNIVERSAL now in Makefile.${CSIRO_OSTYPE}
#UNAME=$(shell uname)
# for universal binaries on OSX
#ifeq "$(UNAME)" "Darwin"
#    CFLAGS_UNIVERSAL= -std=gnu99 -arch i386 -arch x86_64
#endif
#else
#    CFLAGS_UNIVERSAL=
#endif

# extra -m32 flag for $COPTIONS: now in Makefile.${CSIRO_OSTYPE}
# M32=$(shell echo ${CSIRO_OSTYPE}  | awk '/i.86/{mflag=" -m32 "}END{print mflag}')


CLASSES = DJBFFT.class \

.SUFFIXES: .java .class 

all: DJBFFT

DJBFFT: NATIVE_LIB DJBFFT.java
	javac  -classpath ${CLASSPATH} ${JAVAOPTIONS} DJBFFT.java
	# Note full package name!! for header
	javac -h . -classpath ${CLASSPATH} DJBFFT.java
	# Note here that having the .c files *first* is important for linking
	${CC} ${COPTIONS} $(CFLAGS_UNIVERSAL) ${CSHAREDLIBFLAG} ${INCLUDEPATH} -I./djbfft-0.76/include DJBFFTImp.c \
	./djbfft-0.76/${CSIRO_OSTYPE}/lib/djbfft.a -o ./djbfft-0.76/${CSIRO_OSTYPE}/lib/libDJBFFT.${DYLIBSUFFIX}

NATIVE_LIB:
	echo "Compiling native djbfft in " ./djbfft-0.76/ 
	cd ./djbfft-0.76/ ; make clean; make; cd ..



# Remove all class files from the project directory. 
clean: 
	cd ./djbfft-0.76/ ; make clean; cd ..
	rm *.class 

# Implicit rules
.java.class:
	$(COMPILER) -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?

