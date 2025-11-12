#!/bin/sh

# define elsewhere
#JDKBINPATH="/usr/java/j2sdk1.4.0/bin"
#DEDITHOME="/home/ahryat/DEdit"

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/jhall.jar:$XMLEDITHOME/jars/servlet-api.jar

$JDKBINPATH/javac -classpath $CLASSPATH -d $XMLEDITHOME/classes/ \
		  $XMLEDITHOME/src/*/*.java $*

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes

$JDKBINPATH/javac -classpath $CLASSPATH -d $XMLEDITHOME/plugins/classes/ \
		  $XMLEDITHOME/plugins/src/*/*.java $*
#needs classpath stuff actually