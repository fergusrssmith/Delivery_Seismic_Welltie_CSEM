#!/bin/sh


#
#   Generates a Makefile for producing binaries.
#


. $ROOT_DIR/config/make.functions.sh


echo "Generating Makefile..."

rm -f Makefile
make_bins_makefile > Makefile
touch Makefile.depend
chmod g+rw Makefile Makefile.depend

echo "Generating Makefile...done"


exit 0
