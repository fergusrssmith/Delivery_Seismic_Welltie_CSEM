#!/usr/bin/env bash
#
# All CAPITALISED variables imported from this file
source environment_setup.sh

CLASSPATH=$SEISMIC_INVERTER_HOME/jars/colt.jar:$SEISMIC_INVERTER_HOME/jars/FreeHEP.jar:$SEISMIC_INVERTER_HOME/au/csiro/JamesGunning/lib:$SEISMIC_INVERTER_HOME/


java -classpath $CLASSPATH au.csiro.JamesGunning.seismicinverter.InversionAnalyser $*
 
