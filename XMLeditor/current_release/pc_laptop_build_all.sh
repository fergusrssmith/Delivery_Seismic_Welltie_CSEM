#!/bin/sh

# define elsewhere
DEDITHOME="."
JDKBINPATH="/cygdrive/c/j2sdk1.4.1_01/bin"
#DEDITHOME="/home/ahryat/DEdit"

CLASSPATH="$DEDITHOME/jars/msv.jar;$DEDITHOME/jars/isorelax.jar;$DEDITHOME/jars/relaxngDatatype.jar;$DEDITHOME/jars/xsdlib.jar;$DEDITHOME/jars/servlet-api.jar;$DEDITHOME/jars/jhall.jar"

$JDKBINPATH/javac -classpath $CLASSPATH -d $DEDITHOME/classes/ \
		  $DEDITHOME/src/*/*.java $*

CLASSPATH="$CLASSPATH;$DEDITHOME/classes"

$JDKBINPATH/javac -classpath $CLASSPATH -d $DEDITHOME/plugins/classes/ \
		  $DEDITHOME/plugins/src/*/*.java $*
#needs classpath stuff actually
