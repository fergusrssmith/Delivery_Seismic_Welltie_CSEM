#!/bin/sh

# define elsewhere
#XMLJDKBINPATH="/usr/java/j2sdk1.4.0/bin"
#XMLEDITHOME="/home/ahryat/DEdit"

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/classes:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes:$XMLEDITHOME/plugins/classes/:$XMLEDITHOME

$JDKBINPATH/java -classpath $CLASSPATH model.DEditor $*
#java -classpath $CLASSPATH model.DEditor $*

