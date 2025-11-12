#!/usr/bin/env bash
#
# This script requires JDK >= 1.5, 
# and runs the XML editor code
#

# All CAPITALISED variables imported from this file 
source environment_setup.sh

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/classes:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes:$XMLEDITHOME/plugins/classes/:$XMLEDITHOME

# standard invocation in William's version
# java -classpath $CLASSPATH model.DEditor $*
# older version: XSD=${XMLEDITHOME}/xsd/Delivery.xsd
# xsd with source...
XSD=${SEISMIC_INVERTER_HOME}/au/csiro/JamesGunning/seismicinverter/xsd/Delivery.xsd
echo "Using xsd  " $XSD
# modified with help for Delivery, and standard delivery schema
java -Dsun.java2d.xrender=false -classpath $CLASSPATH model.DEditor plugin=delivery.DeliveryPlugin xsd=${XSD} $*

# CYGWIN users: modify the above to 
# java -classpath "$(cygpath -wp $CLASSPATH)" model.DEditor plugin=delivery.DeliveryPlugin xsd=${XSD} $*
