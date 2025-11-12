
###############################################################################
#
#	Makefile.java
#
#	Targets and macros for making java classes.
#
#	Required variables/macros:
#
#		ROOT_DIR
#		JAVA_SOURCES
#
#	Optional macros:
#
#		JAVA_TITLE
#
###############################################################################



#
#	MACROS
#

JAVA_CLASSES	= $(JAVA_SOURCES:%.java=%.class)
JAVAC		= javac



#
#	ACTION TARGETS
#

.PHONY: all clean


all: $(JAVA_CLASSES)


clean: 
	rm -f *.class



#
#	GENERIC TARGETS
#

%.class: %.java
	$(JAVAC) $<

