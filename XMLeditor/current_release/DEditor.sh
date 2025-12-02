#!/bin/sh

# define elsewhere
JDKBINPATH="/usr/java/j2sdk1.4.0_01/bin"
DEDITHOME="/home/ahryat/DEdit"

CLASSPATH=$DEDITHOME/jars/msv.jar:$DEDITHOME/jars/isorelax.jar:$DEDITHOME/jars/relaxngDatatype.jar:$DEDITHOME/jars/xsdlib.jar:$DEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$DEDITHOME/classes:$DEDITHOME/plugins/classes/:$DEDITHOME

$JDKBINPATH/java -classpath $CLASSPATH model.DEditor $*

