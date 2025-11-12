#!/usr/bin/env bash
#
# This script requires JDK >= 1.5, 
# and runs the XML editor code for the wavelet extraction
#

# All CAPITALISED variables imported from this file 
source environment_setup.sh

CLASSPATH=$XMLEDITHOME/jars/msv.jar:$XMLEDITHOME/classes:$XMLEDITHOME/jars/isorelax.jar:$XMLEDITHOME/jars/xsdlib.jar:$XMLEDITHOME/jars/relaxngDatatype.jar:$XMLEDITHOME/jars/jhall.jar

CLASSPATH=$CLASSPATH:$XMLEDITHOME/classes:$XMLEDITHOME/plugins/classes/:$XMLEDITHOME

# standard invocation in William's version
# java -classpath $CLASSPATH model.DEditor $*
XSD=$SEISMIC_INVERTER_HOME/au/csiro/JamesGunning/waveletextraction/xsd/WaveletExtraction.xsd
echo "Using xsd  " $XSD
# modified for standard wavelet extractor schema
# java bug workaround -Dsun.java2d.xrender=false
java -Dsun.java2d.xrender=false -classpath $CLASSPATH model.DEditor xsd=${XSD} $*

# CYGWIN users: modify the above to 
# java -classpath "$(cygpath -wp $CLASSPATH)" model.DEditor xsd=${XSD} $*
