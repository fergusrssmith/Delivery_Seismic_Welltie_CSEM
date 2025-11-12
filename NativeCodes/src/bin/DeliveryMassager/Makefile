
#
#	BASE
#

ROOT_DIR	= ../../..

BIN_SOURCES	= deliveryMassager.c analyseOutputCompacter.c eclipseXsection.c getCholeskySafeTest.c
AUX_SOURCES	= clapack_cblas_wrappers.c errorHandler.c util.c modified_left_Cholesky.c massager.c xmlParser.c variogram.c upperTriangleMatrices.c solver.c decorator.c eclipseHandler.c coarseGridMaker.c su_util.c
BIN_CFLAGS	= -I/usr/include/libxml2 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include   
BIN_CPPFLAGS	= 
BIN_FFLAGS	= 
BIN_LDFLAGS	= -L/usr/lib64/atlas -llapacke -lblas -llapack_atlas -lglib-2.0 -lxml2 -lECLIO -lcommandline -lkdtree -lcgeom -lNR -lm

#include $(ROOT_DIR)/config/Makefile.bins


#
#	MACROS
#

include $(ROOT_DIR)/config/Makefile.ostype.$(CSIRO_OSTYPE)

BIN_FILES	= $(CSIRO_OSTYPE)/deliveryMassager $(CSIRO_OSTYPE)/analyseOutputCompacter $(CSIRO_OSTYPE)/eclipseXsection $(CSIRO_OSTYPE)/getCholeskySafeTest
BIN_OBJS	= $(CSIRO_OSTYPE)/deliveryMassager.o $(CSIRO_OSTYPE)/analyseOutputCompacter.o $(CSIRO_OSTYPE)/eclipseXsection.o $(CSIRO_OSTYPE)/getCholeskySafeTest.o
AUX_OBJS	= $(CSIRO_OSTYPE)/clapack_cblas_wrappers.o $(CSIRO_OSTYPE)/errorHandler.o $(CSIRO_OSTYPE)/util.o $(CSIRO_OSTYPE)/modified_left_Cholesky.o $(CSIRO_OSTYPE)/massager.o $(CSIRO_OSTYPE)/xmlParser.o $(CSIRO_OSTYPE)/variogram.o $(CSIRO_OSTYPE)/upperTriangleMatrices.o $(CSIRO_OSTYPE)/solver.o $(CSIRO_OSTYPE)/decorator.o $(CSIRO_OSTYPE)/eclipseHandler.o $(CSIRO_OSTYPE)/coarseGridMaker.o $(CSIRO_OSTYPE)/su_util.o

INSTALL_DIR	= $(ROOT_DIR)/bin/$(CSIRO_OSTYPE)

CFLAGS		= $(SYS_CFLAGS) $(BIN_CFLAGS) $(AUX_CFLAGS) -I$(ROOT_DIR)/include
CPPFLAGS	= $(SYS_CPPFLAGS) $(BIN_CPPFLAGS) $(AUX_CPPFLAGS) -I$(ROOT_DIR)/include
FFLAGS		= $(SYS_FFLAGS) $(BIN_FFLAGS) $(AUX_FFLAGS)
LDFLAGS		= $(SYS_LDFLAGS) $(BIN_LDFLAGS) $(AUX_LDFLAGS) -L$(ROOT_DIR)/lib/$(CSIRO_OSTYPE)


#
#	ACTION TARGETS
#

.PHONY: all debug optimise fast install clean make world depend hosts


all: $(CSIRO_OSTYPE) $(BIN_FILES)


debug: 
	@make all \
	AUX_CFLAGS="$(SYS_CFLAGS_DEBUG)" \
	AUX_CPPFLAGS="$(SYS_CPPFLAGS_DEBUG)" \
	AUX_FFLAGS="$(SYS_FFLAGS_DEBUG)" \
	AUX_LDFLAGS="$(SYS_LDFLAGS_DEBUG)"


optimise fast: 
	@make all \
	AUX_CFLAGS="$(SYS_CFLAGS_OPTIMISE)" \
	AUX_CPPFLAGS="$(SYS_CPPFLAGS_OPTIMISE)" \
	AUX_FFLAGS="$(SYS_FFLAGS_OPTIMISE)" \
	AUX_LDFLAGS="$(SYS_LDFLAGS_OPTIMISE)"


install: all $(INSTALL_DIR)
	$(CP) $(BIN_FILES) $(INSTALL_DIR)


clean: 
	$(RM) -r $(CSIRO_OSTYPE)


make: 
	@make -f Makefile.base make


world: 
	@make make; make depend; make all


depend: 
	@if ( gcc -v > /dev/null 2> /dev/null );\
	then\
		echo "Making dependencies...";\
		$(RM) Makefile.depend;\
		touch Makefile.depend;\
		gcc -MM $(CFLAGS) $(CPPFLAGS) $(BIN_SOURCES) $(AUX_SOURCES) | sed "s|^.*\.o:|\$$(CSIRO_OSTYPE)\/&|g" >> Makefile.depend;\
	       chmod g+rw Makefile.depend;\
		echo "Making dependencies...done";\
	else\
		echo "Making dependencies...ignored (no gcc)";\
	fi


hosts: 
	@(\
	   if [ -n "$(CSIRO_HOSTS)" ]; \
	   then \
	       targets=`echo $(MAKECMDGOALS) | sed 's|hosts||g'`; \
	       for host in $(CSIRO_HOSTS); \
	       do \
	           echo ""; \
	           echo "Making" $$targets on $$host; \
	           rsh $$host "cd $$PWD; make $$targets"; \
	       done; \
	       echo ""; \
	   fi \
	)


#
#	FILE / DIRECTORY TARGETS
#


$(CSIRO_OSTYPE): 
	mkdir -p $@
	chmod g+rwX $@


$(INSTALL_DIR): 
	mkdir -p $@
	chmod g+rwX $@


$(CSIRO_OSTYPE)/deliveryMassager: $(CSIRO_OSTYPE)/deliveryMassager.o $(AUX_OBJS)
	$(CC) $(CSIRO_OSTYPE)/deliveryMassager.o $(AUX_OBJS) $(LDFLAGS) -o $@
	chmod g+rwX $@


$(CSIRO_OSTYPE)/analyseOutputCompacter: $(CSIRO_OSTYPE)/analyseOutputCompacter.o $(AUX_OBJS)
	$(CC) $(CSIRO_OSTYPE)/analyseOutputCompacter.o $(AUX_OBJS) $(LDFLAGS) -o $@
	chmod g+rwX $@


$(CSIRO_OSTYPE)/eclipseXsection: $(CSIRO_OSTYPE)/eclipseXsection.o $(AUX_OBJS)
	$(CC) $(CSIRO_OSTYPE)/eclipseXsection.o $(AUX_OBJS) $(LDFLAGS) -o $@
	chmod g+rwX $@


$(CSIRO_OSTYPE)/getCholeskySafeTest: $(CSIRO_OSTYPE)/getCholeskySafeTest.o $(AUX_OBJS)
	$(CC) $(CSIRO_OSTYPE)/getCholeskySafeTest.o $(AUX_OBJS) $(LDFLAGS) -o $@
	chmod g+rwX $@


#
#	OBJECT TARGETS
#


$(CSIRO_OSTYPE)/deliveryMassager.o: 
	$(CC) $(CFLAGS) -c deliveryMassager.c -o $(CSIRO_OSTYPE)/deliveryMassager.o


$(CSIRO_OSTYPE)/analyseOutputCompacter.o: 
	$(CC) $(CFLAGS) -c analyseOutputCompacter.c -o $(CSIRO_OSTYPE)/analyseOutputCompacter.o


$(CSIRO_OSTYPE)/eclipseXsection.o: 
	$(CC) $(CFLAGS) -c eclipseXsection.c -o $(CSIRO_OSTYPE)/eclipseXsection.o


$(CSIRO_OSTYPE)/getCholeskySafeTest.o: 
	$(CC) $(CFLAGS) -c getCholeskySafeTest.c -o $(CSIRO_OSTYPE)/getCholeskySafeTest.o


$(CSIRO_OSTYPE)/clapack_cblas_wrappers.o: 
	$(CC) $(CFLAGS) -c clapack_cblas_wrappers.c -o $(CSIRO_OSTYPE)/clapack_cblas_wrappers.o


$(CSIRO_OSTYPE)/errorHandler.o: 
	$(CC) $(CFLAGS) -c errorHandler.c -o $(CSIRO_OSTYPE)/errorHandler.o


$(CSIRO_OSTYPE)/util.o: 
	$(CC) $(CFLAGS) -c util.c -o $(CSIRO_OSTYPE)/util.o


$(CSIRO_OSTYPE)/modified_left_Cholesky.o: 
	$(CC) $(CFLAGS) -c modified_left_Cholesky.c -o $(CSIRO_OSTYPE)/modified_left_Cholesky.o


$(CSIRO_OSTYPE)/massager.o: 
	$(CC) $(CFLAGS) -c massager.c -o $(CSIRO_OSTYPE)/massager.o


$(CSIRO_OSTYPE)/xmlParser.o: 
	$(CC) $(CFLAGS) -c xmlParser.c -o $(CSIRO_OSTYPE)/xmlParser.o


$(CSIRO_OSTYPE)/variogram.o: 
	$(CC) $(CFLAGS) -c variogram.c -o $(CSIRO_OSTYPE)/variogram.o


$(CSIRO_OSTYPE)/upperTriangleMatrices.o: 
	$(CC) $(CFLAGS) -c upperTriangleMatrices.c -o $(CSIRO_OSTYPE)/upperTriangleMatrices.o


$(CSIRO_OSTYPE)/solver.o: 
	$(CC) $(CFLAGS) -c solver.c -o $(CSIRO_OSTYPE)/solver.o


$(CSIRO_OSTYPE)/decorator.o: 
	$(CC) $(CFLAGS) -c decorator.c -o $(CSIRO_OSTYPE)/decorator.o


$(CSIRO_OSTYPE)/eclipseHandler.o: 
	$(CC) $(CFLAGS) -c eclipseHandler.c -o $(CSIRO_OSTYPE)/eclipseHandler.o


$(CSIRO_OSTYPE)/coarseGridMaker.o: 
	$(CC) $(CFLAGS) -c coarseGridMaker.c -o $(CSIRO_OSTYPE)/coarseGridMaker.o


$(CSIRO_OSTYPE)/su_util.o: 
	$(CC) $(CFLAGS) -c su_util.c -o $(CSIRO_OSTYPE)/su_util.o


#
#	DEPENDENCIES
#
include Makefile.depend
