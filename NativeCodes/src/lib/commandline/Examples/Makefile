
#
#	BASE
#

ROOT_DIR	= ../../../..

BIN_SOURCES	= commandline.c 
AUX_SOURCES	= 
BIN_CFLAGS	= 
BIN_CPPFLAGS	= 
BIN_FFLAGS	= 
BIN_LDFLAGS	= -lcommandline

#include $(ROOT_DIR)/config/Makefile.bins


#
#	MACROS
#

include $(ROOT_DIR)/config/Makefile.ostype.$(CSIRO_OSTYPE)

BIN_FILES	= $(CSIRO_OSTYPE)/commandline
BIN_OBJS	= $(CSIRO_OSTYPE)/commandline.o
AUX_OBJS	=

INSTALL_DIR	= $(ROOT_DIR)/bin/$(CSIRO_OSTYPE)

CFLAGS		= $(BIN_CFLAGS) $(SYS_CFLAGS) $(AUX_CFLAGS) -I$(ROOT_DIR)/include
CPPFLAGS	= $(BIN_CPPFLAGS) $(SYS_CPPFLAGS) $(AUX_CPPFLAGS) -I$(ROOT_DIR)/include
FFLAGS		= $(BIN_FFLAGS) $(SYS_FFLAGS) $(AUX_FFLAGS)
LDFLAGS		= $(BIN_LDFLAGS) $(SYS_LDFLAGS) $(AUX_LDFLAGS) -L$(ROOT_DIR)/lib/$(CSIRO_OSTYPE)


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
	umask 002; $(CP) $(BIN_FILES) $(INSTALL_DIR)


clean: 
	$(RM) -r $(CSIRO_OSTYPE)


make: 
	@umask 002;  make -f Makefile.base make


world: 
	@make make; make depend; make all


depend: 
	@if ( gcc -v > /dev/null 2> /dev/null );\
	then\
		echo "Making dependencies...";\
		$(RM) Makefile.depend;\
		umask 002;  touch Makefile.depend;\
		umask 002;  gcc -MM $(CFLAGS) $(CPPFLAGS) $(BIN_SOURCES) $(AUX_SOURCES) | sed "s|^.*\.o:|\$$(CSIRO_OSTYPE)\/&|g" >> Makefile.depend;\
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
	umask 002; mkdir -p $@


$(INSTALL_DIR): 
	umask 002; mkdir -p $@


$(CSIRO_OSTYPE)/commandline: $(CSIRO_OSTYPE)/commandline.o $(AUX_OBJS)
	umask 002; $(CC) $(CSIRO_OSTYPE)/commandline.o $(AUX_OBJS) $(LDFLAGS) -o $@


#
#	OBJECT TARGETS
#


$(CSIRO_OSTYPE)/commandline.o: 
	umask 002; $(CC) $(CFLAGS) -c commandline.c -o $(CSIRO_OSTYPE)/commandline.o


#
#	DEPENDENCIES
#
include Makefile.depend
