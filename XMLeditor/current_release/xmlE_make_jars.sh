#!/bin/sh

#define elsewhere
#JDKBINPATH="/usr/java/j2sdk1.4.0_01/bin"
#XMLEDITHOME="/home/ahryat/DEdit"

cd $XMLEDITHOME/classes
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/xmleditor.jar *
cd ..
$JDKBINPATH/jar uvf $XMLEDITHOME/jars/xmleditor.jar configFiles/ icons/
cd plugins/classes
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/scriptplugin.jar script/
cd $XMLEDITHOME
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/xsd.jar xsd/
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/xehelp.jar help/
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/scrhelp.jar plugins/help/script/
$JDKBINPATH/jar cvf $XMLEDITHOME/jars/delivhelp.jar plugins/help/delivery/
