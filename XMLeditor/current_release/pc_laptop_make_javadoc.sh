#!/bin/sh

# define elsewhere
DEDITHOME="."
JDKBINPATH="/cygdrive/c/j2sdk1.4.1_01/bin"

rm -rf $DEDITHOME/javadoc/*

CLASSPATH="$DEDITHOME/jars/msv.jar;$DEDITHOME/jars/isorelax.jar;$DEDITHOME/jars/relaxngDatatype.jar;$DEDITHOME/jars/xsdlib.jar;$DEDITHOME/jars/jhall.jar"

CLASSPATH="$CLASSPATH;$DEDITHOME/classes"

$JDKBINPATH/javadoc -classpath $CLASSPATH -d $DEDITHOME/javadoc src/*/*.java \
    $DEDITHOME/plugins/src/*/*.java
