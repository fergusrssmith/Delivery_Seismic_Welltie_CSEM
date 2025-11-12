
#
#	BASE
#

ROOT_DIR	= ../../..

LIB_NAME	= NR
LIB_SOURCES	= nrutil.c
LIB_CFLAGS	= 
LIB_CPPFLAGS	= 
LIB_FFLAGS	= 
LIB_LDFLAGS	= 

#include $(ROOT_DIR)/config/Makefile.lib


#
#	MACROS
#

include $(ROOT_DIR)/config/Makefile.ostype.$(CSIRO_OSTYPE)

LIB_FILE	= $(CSIRO_OSTYPE)/lib$(LIB_NAME).a
LIB_OBJS	= $(CSIRO_OSTYPE)/nrutil.o

LIB_INSTALL_DIR	= $(ROOT_DIR)/lib/$(CSIRO_OSTYPE)
INC_INSTALL_DIR	= $(ROOT_DIR)/include/$(LIB_NAME)

CFLAGS		= $(SYS_CFLAGS) $(LIB_CFLAGS) $(AUX_CFLAGS) -I.. -I$(ROOT_DIR)/include
CPPFLAGS	= $(SYS_CPPFLAGS) $(LIB_CPPFLAGS) $(AUX_CPPFLAGS) -I.. -I$(ROOT_DIR)/include
FFLAGS		= $(SYS_FFLAGS) $(LIB_FFLAGS) $(AUX_FFLAGS)
LDFLAGS		= $(SYS_LDFLAGS) $(LIB_LDFLAGS) $(AUX_LDFLAGS)


#
#	ACTION TARGETS
#

.PHONY: all debug optimise fast install clean make world depend hosts


all: $(CSIRO_OSTYPE) $(LIB_FILE)


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


install: all $(LIB_INSTALL_DIR) $(INC_INSTALL_DIR)
	$(CP) $(LIB_FILE) $(LIB_INSTALL_DIR)
	$(RM) -r $(INC_INSTALL_DIR)
	mkdir -p $(INC_INSTALL_DIR)
	$(CP) *.h $(INC_INSTALL_DIR)
	chmod -R g+rwX $(INC_INSTALL_DIR)


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
		gcc -MM $(CFLAGS) $(CPPFLAGS) $(LIB_SOURCES) | sed "s|^.*\.o:|\$$(CSIRO_OSTYPE)\/&|g" >> Makefile.depend;\
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


$(LIB_INSTALL_DIR): 
	mkdir -p $@
	chmod g+rwX $@


$(INC_INSTALL_DIR): 
	mkdir -p $@
	chmod g+rwX $@


$(LIB_FILE): $(CSIRO_OSTYPE) $(LIB_OBJS)
	$(RM) $@
	$(AR) -q $@ $(LIB_OBJS)
	chmod g+rw $@


#
#	OBJECT TARGETS
#


$(CSIRO_OSTYPE)/nrutil.o: 
	$(CC) $(CFLAGS) -c nrutil.c -o $(CSIRO_OSTYPE)/nrutil.o


#
#	DEPENDENCIES
#
include Makefile.depend
