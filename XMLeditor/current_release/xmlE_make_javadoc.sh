#!/bin/sh

# define elsewhere
#JDKBINPATH="/usr/java/j2sdk1.4.0_01/bin"
#XMLEDITHOME="/home/ahryat/DEdit"

rm -rf $XMLEDITHOME/javadoc/*

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/jhall.jar:$XMLEDITHOME/jars/servlet-api.jar

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes

$JDKBINPATH/javadoc -classpath $CLASSPATH -d $XMLEDITHOME/javadoc \
    src/*/*.java $XMLEDITHOME/plugins/src/*/*.java
