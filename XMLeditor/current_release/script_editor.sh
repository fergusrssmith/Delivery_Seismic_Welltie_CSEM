#!/bin/sh

xmledit xsd=$XMLEDITHOME/xsd/script.xsd cfg=$XMLEDITHOME/configFiles/scriptconfig.xml -S plugin=script.ScriptPlugin $*
