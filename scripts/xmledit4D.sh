#!/usr/bin/env bash
#
# This script runs the 4D XML editor code
#

# All CAPITALISED variables imported from this file 
source environment_setup.sh

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/classes:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes:$XMLEDITHOME/plugins/classes/:$XMLEDITHOME

# xsd with source...
XSD=${SEISMIC_INVERTER_HOME}/au/csiro/JamesGunning/seismicinverter4D/xsd/Delivery4D.xsd
if [ ! -f ${XSD} ]
then
    echo "You don't have the 4D distribution in \
${SEISMIC_INVERTER_HOME}/au/csiro/JamesGunning/seismicinverter4D - can't run the xmleditor!"
	exit	
fi

echo "Using xsd  " $XSD
# modified with help for Delivery, and standard delivery schema
java -Dsun.java2d.xrender=false -classpath $CLASSPATH model.DEditor plugin=delivery.DeliveryPlugin xsd=${XSD} $*

# CYGWIN users: modify the above to 
# java -classpath "$(cygpath -wp $CLASSPATH)" model.DEditor plugin=delivery.DeliveryPlugin xsd=${XSD} $*
