SEISMIC_INVERTER_HOME=../../../../..
CLASSPATH = $(SEISMIC_INVERTER_HOME):$(SEISMIC_INVERTER_HOME)/jars/colt.jar

COMPILER = javac
JAVAOPTIONS = -g -O -deprecation    

CLASSES = DifferentialEvolution.class \


.SUFFIXES: .java .class 

all: $(CLASSES) 


# Remove all class files from the project directory. 
clean: 
	rm *.class 

# Implicit rule for making Java class files from Java source files. 
# Those who use jikes will prefer leaving off CLASSPATH and having a
# suitable JIKESPATH
#.java.class:
.java.class:
		$(COMPILER) -classpath $(CLASSPATH) $(JAVAOPTIONS) $?


