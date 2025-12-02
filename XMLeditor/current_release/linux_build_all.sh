#!/bin/sh

# define elsewhere
DEDITHOME="/home/polonius/jsg/java/XMLeditor/current_release"
JDKBINPATH="/usr/java/jdk1.5.0_07/bin/"
#DEDITHOME="/home/ahryat/DEdit"

CLASSPATH=../../jars/javax_servlet.jar:$DEDITHOME/jars/msv.jar:$DEDITHOME/jars/isorelax.jar:$DEDITHOME/jars/relaxngDatatype.jar:$DEDITHOME/jars/xsdlib.jar:$DEDITHOME/jars/jhall.jar

$JDKBINPATH/javac -source 1.4 -nowarn -classpath $CLASSPATH -d $DEDITHOME/classes/ \
		  $DEDITHOME/src/*/*.java $*

CLASSPATH=$CLASSPATH:$DEDITHOME/classes

$JDKBINPATH/javac -source 1.4 -nowarn -classpath $CLASSPATH -d $DEDITHOME/plugins/classes/ \
		  $DEDITHOME/plugins/src/*/*.java $*
#needs classpath stuff actually
