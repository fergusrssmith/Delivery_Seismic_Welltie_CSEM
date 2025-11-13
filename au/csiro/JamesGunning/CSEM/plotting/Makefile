SEISMIC_INVERTER_HOME=../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar:$(SEISMIC_INVERTER_HOME)/jars/FreeHEP.jar

# compile as optimized, or debug....
JAVAOPTIONS = -deprecation   
ifeq (debug,$(findstring debug,$(MAKECMDGOALS)))
  JAVAOPTIONS += -g 
endif


# targets
CLASSES = OneDModelStatsPlotter.class  Plot1DCSEM.class  Plot1DCSEMSection.class  PlotMARE2D.class  PlotTriangle.class

.SUFFIXES: .java .class 

all: $(CLASSES) 
debug: all

# General cleanup

clean: 
	rm *.class

# Implicit rules
.java.class:
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?
