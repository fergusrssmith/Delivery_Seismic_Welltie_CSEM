#!/bin/sh
XMLJDKBINPATH="/usr/site/jdk-1.4/bin"
XMLEDITHOME="/home/williamr/shigoto/beta"
rice_xmledit.sh xsd=$XMLEDITHOME/xsd/script.xsd cfg=$XMLEDITHOME/configFiles/scriptconfig.xml -S plugin=script.ScriptPlugin $*
