# CLASSSPATH locally defined
COMPILER = javac
JAVAOPTIONS = -g -O -deprecation    

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

.SUFFIXES: .java .class 

all: classes


# Remove all class files from the project directory. 
clean: 
	rm *.class 

# Implicit rules
.java.class:
	javac -classpath $(CLASSPATH) $(JAVAOPTIONS)  $?

