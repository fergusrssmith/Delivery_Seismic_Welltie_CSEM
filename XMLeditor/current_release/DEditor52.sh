#!/bin/sh

# define elsewhere
JDKBINPATH="/usr/j2se/bin"
DEDITHOME="/home/ahryat/DEdit"

CLASSPATH=$DEDITHOME/jars/msv.jar:$DEDITHOME/jars/isorelax.jar:$DEDITHOME/jars/relaxngDatatype.jar:$DEDITHOME/jars/xsdlib.jar:$DEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$DEDITHOME/classes

$JDKBINPATH/java -classpath $CLASSPATH model.DEditor $*

