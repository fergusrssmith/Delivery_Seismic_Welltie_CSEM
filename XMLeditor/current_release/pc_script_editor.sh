#!/bin/sh

#XMLEDITHOME="."
XMLEDITHOME="c:/cygwin/home/williamr/shigoto/rls030712"

xmledit xsd=$XMLEDITHOME/xsd/script.xsd cfg=$XMLEDITHOME/configFiles/scriptconfig.xml -S plugin=script.ScriptPlugin $*
