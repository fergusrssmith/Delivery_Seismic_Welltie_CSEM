#!/usr/bin/env bash
#
# Set up of environmental stuff for all delivery scripts
#
# Sensible default definition of the environmental variables
#
# SEISMIC_INVERTER_HOME and XMLEDITHOME
#
# SEISMIC_INVERTER_HOME is the directory where you install the
# source tree. It behaves rather like CWPROOT in the Seismic Unix (SU) install, for people
# who are familiar with installing SU

# This sets SEISMIC_INVERTER_HOME as the directory at the root of the
# source archive. It contains the directories
# au, jars, scripts. Users can set this variable in their .bashrc or .cshrc
# file if they'd prefer to never have to edit this script

if [ -z "$SEISMIC_INVERTER_HOME" ]
then
	SEISMIC_INVERTER_HOME=$HOME/java
fi

# Libs currently built for only 32 bit java
MODEL=""
#if [[ $OSTYPE == *darwin* ]]; then 
#  MODEL=-d32
#fi 

# An environment variable for the XML editor home area in the 
# Delivery distribution. Users with another XMLeditor tree may override this 
# by setting XMLEDITHOME in their .bashrc or .cshrc
# (again, avoiding having to edit this script)

if [ -z "$XMLEDITHOME" ]
then
  XMLEDITHOME=$SEISMIC_INVERTER_HOME/XMLeditor/current_release
fi



