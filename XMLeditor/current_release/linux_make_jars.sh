
JDKBINPATH="/usr/java/j2sdk1.4.0_01/bin"
DEDITHOME="/home/ahglim/prg/bhp_apps/xmledit/beta"

cd $DEDITHOME/classes
$JDKBINPATH/jar cvf $DEDITHOME/jars/xmleditor.jar *
cd ..
$JDKBINPATH/jar uvf $DEDITHOME/jars/xmleditor.jar configFiles/ icons/
cd plugins/classes
$JDKBINPATH/jar cvf $DEDITHOME/jars/scriptplugin.jar script/
$JDKBINPATH/jar cvf $DEDITHOME/jars/deliveryplugin.jar delivery/ 
cd $DEDITHOME
$JDKBINPATH/jar cvf $DEDITHOME/jars/xsd.jar xsd/
$JDKBINPATH/jar cvf $DEDITHOME/jars/xehelp.jar help/
$JDKBINPATH/jar cvf $DEDITHOME/jars/scrhelp.jar plugins/help/script/ 
$JDKBINPATH/jar cvf $DEDITHOME/jars/delivhelp.jar plugins/help/delivery/
