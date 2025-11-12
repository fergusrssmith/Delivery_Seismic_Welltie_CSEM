# Platform dependent things in this include....
# i.e. provide one such file for each desired OS
include Makefile.${CSIRO_OSTYPE}

# CLASSSPATH locally defined
SEISMIC_INVERTER_HOME=../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar

JAVAOPTIONS =   

# CLASSSPATH locally defined
COPTIONS = -O2 

# this magic gets all changed .java files compiled in one javac invocation :-)
# Get listing of *.class targets
JSOURCES = $(wildcard *.java)
CLASSES:=$(patsubst %.java,%.class,$(JSOURCES))
.PHONY: classes
LIST:=
classes: $(CLASSES)
	if [ ! -z "$(LIST)" ] ; then \
		javac ${JAVAOPTIONS} ${LIST} ; \
	fi

$(CLASSES) : %.class : %.java
	$(eval LIST+=$$<)

all: classes MathUtil BlasLapackUtil

MathUtil: MathUtilImp.c
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) MathUtil.java
	# Note full package name!! for header
	javac -h . -classpath ${CLASSPATH} MathUtil.java
	# Note here that having the .c files *first* is important for linking
	cc $(COPTIONS) ${CSHAREDLIBFLAG} ${INCLUDEPATH} MathUtilImp.c -o ./${CSIRO_OSTYPE}/libMathUtil.${DYLIBSUFFIX}

BlasLapackUtil: BlasLapackUtilImp.c
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS) BlasLapackUtil.java
	# Note full package name!! for header  
	javac -h . -classpath ${CLASSPATH}:.. BlasLapackUtil.java
	# Note here that having the .c files *first* is important for linking
	cc $(COPTIONS) ${ATLASOPTS} ${CSHAREDLIBFLAG} ${INCLUDEPATH} BlasLapackUtilImp.c -o ./${CSIRO_OSTYPE}/libBlasLapackUtil.${DYLIBSUFFIX}

.SUFFIXES: .java .class 

# Remove all class files from the project directory. 
clean: 
	rm *.class 

# Implicit rules
.java.class:
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?


