SEISMIC_INVERTER_HOME=../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar

COMPILER = javac
JAVAOPTIONS = -g -O -deprecation    

CLASSES =  TriangleUtil.class PolygonIntersect.class \



.SUFFIXES: .java .class 

all: $(CLASSES) 


# Remove all class files from the project directory. 
clean: 
	rm *.class 


# Implicit rules
.java.class:
	$(COMPILER) -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?

