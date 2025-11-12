


#
#   MISC. FUNCTIONS
#


make_target()
{
    target="$1"
    depend="$2"
    shift
    shift
    
    echo ""
    echo ""
    echo "$target: $depend"
    
    for action
    do
        if [ -n "$action" ];
        then
            echo "	$action"
        fi
    done
}




make_directory_target()
{
    make_target "$1" "" "mkdir -p \$@" "chmod g+rwX \$@"

#    make_target "$1" "" "mkdir -p \$@"
}










#
#   SOURCE / OBJECT / BINARY  FUNCTIONS
#


set_compiler()
{
    compiler=""
    flags=""
    
    case "$1" in
        *.c)
            compiler="\$(CC)"
            flags="\$(CFLAGS)"
            ;;
        *.cc)
            compiler="\$(CPP)"
            flags="\$(CPPFLAGS)"
            ;;
        *.cpp)
            compiler="\$(CPP)"
            flags="\$(CPPFLAGS)"
            ;;
        *.f)
            compiler="\$(FC)"
            flags="\$(FFLAGS)"
            ;;
    esac
}




set_compiler_list()
{
    level=0
    compiler="\$(CC)"
    flags="\$(CFLAGS)"
    
    for src
    do
        case "$src" in

            *.f)
                if [ $level -lt 3 ];
                then
                    level=3
                    compiler="\$(FC)"
                    flags="\$(FFLAGS)"
                fi
                ;;

            *.cc)
                if [ $level -lt 2 ];
                then
                    level=2
                    compiler="\$(CPP)"
                    flags="\$(CPPFLAGS)"
                fi
                ;;

            *.cpp)
                if [ $level -lt 2 ];
                then
                    level=2
                    compiler="\$(CPP)"
                    flags="\$(CPPFLAGS)"
                fi
                ;;

            *.c)
                if [ $level -lt 1 ];
                then
                    level=1
                    compiler="\$(CC)"
                    flags="\$(CFLAGS)"
                fi
                ;;
        esac
        
    done
}




set_obj()
{
    obj="$OS_PREFIX"`echo $1 | sed "s|\..*|\.o|g"`
}




set_bin()
{
    bin="$OS_PREFIX"`echo $1 | sed "s|\..*||g"`"$OS_EXE"
}




set_obj_list()
{
    obj_list=""
    
    for src
    do
        set_obj $src
        obj_list="$obj_list $obj"
    done
}




set_bin_list()
{
    bin_list=""
    
    for src
    do
        set_bin $src
        bin_list="$bin_list $bin"
    done
}




make_obj_target()
{
    set_compiler $1
    set_obj $1

    make_target "$obj" "" "$compiler $flags -c $1 -o $obj"
}




make_obj_target_list()
{
    for src
    do
        make_obj_target $src
    done
}




make_bin_target()
{
    set_compiler $1
    set_obj $1
    set_bin $1

    make_target "$bin" "$obj $2" \
        "$compiler $obj $2 \$(LDFLAGS) -o \$@" \
        "chmod g+rwX \$@"

#    make_target "$bin" "$obj $2" \
#        "$compiler $obj $2 \$(LDFLAGS) -o $bin"
}




make_bin_target_list()
{
    aux="$1"
    shift
    
    for src
    do
        make_bin_target $src $aux
    done
}









make_init()
{
    UNAME=`uname`
    OS_PREFIX=""
    OS_SUFFIX=""
    OS_DEPEND=""
    OS_EXT=""
    OS_EXE=""


    if [ "$OSTYPE" = "cygwin" ];
    then
        OS_EXT=".cygwin"
        OS_EXE=".exe"
    fi

    if [ -n "$CSIRO_OSTYPE" ];
    then
        OS_PREFIX="\$(CSIRO_OSTYPE)/"
        OS_SUFFIX="/\$(CSIRO_OSTYPE)"
        OS_DEPEND="\$(CSIRO_OSTYPE) "
        
        if [ -f $ROOT_DIR/config/Makefile.ostype.$CSIRO_OSTYPE ];
        then
            OS_EXT=".\$(CSIRO_OSTYPE)"
        fi
    fi
    
    
    if [ -z "$OS_EXT" ];
    then
        if [ -f $ROOT_DIR/config/Makefile.ostype.$OSTYPE ];
        then
            OS_EXT=".$OSTYPE"
        elif [ -f $ROOT_DIR/config/Makefile.ostype.$UNAME ];
        then
            OS_EXT=".$UNAME"
        else
            OS_EXT=".default"
        fi
    fi
}











#
#   BASE FUNCTIONS
#


make_base_header()
{
    echo ""
    echo "#"
    echo "#	BASE"
    echo "#"
}









#
#   MACROS FUNCTIONS
#


make_macros_header()
{
    echo ""
    echo ""
    echo "#"
    echo "#	MACROS"
    echo "#"
    echo ""
    echo "include \$(ROOT_DIR)/config/Makefile.ostype$OS_EXT"
}









#
#   ACTION TARGET FUNCTIONS
#



make_action_target_header()
{
    echo ""
    echo ""
    echo "#"
    echo "#	ACTION TARGETS"
    echo "#"
    echo ""
    echo ".PHONY: all debug optimise fast install clean make world depend hosts"
}




make_action_target_all()
{
    make_target "all" "$*"
}




make_action_target_debug()
{
    make_target "debug" "" "@make all \\" \
        "AUX_CFLAGS=\"\$(SYS_CFLAGS_DEBUG)\" \\" \
        "AUX_CPPFLAGS=\"\$(SYS_CPPFLAGS_DEBUG)\" \\" \
        "AUX_FFLAGS=\"\$(SYS_FFLAGS_DEBUG)\" \\" \
        "AUX_LDFLAGS=\"\$(SYS_LDFLAGS_DEBUG)\""
}




make_action_target_optimise()
{
    make_target "optimise fast" "" "@make all \\" \
        "AUX_CFLAGS=\"\$(SYS_CFLAGS_OPTIMISE)\" \\" \
        "AUX_CPPFLAGS=\"\$(SYS_CPPFLAGS_OPTIMISE)\" \\" \
        "AUX_FFLAGS=\"\$(SYS_FFLAGS_OPTIMISE)\" \\" \
        "AUX_LDFLAGS=\"\$(SYS_LDFLAGS_OPTIMISE)\""
}




make_action_target_clean()
{
    if [ -z "$CSIRO_OSTYPE" ];
    then
        make_target "clean" "" "\$(RM) $*"
    else
        make_target "clean" "" "\$(RM) -r \$(CSIRO_OSTYPE)"
    fi
}




make_action_target_make()
{
    make_target "make" "" "@make -f Makefile.base make"
}




make_action_target_world()
{
    make_target "world" "" "@make make; make depend; make all"
}




make_action_target_depend()
{
#    if [ -z "$CSIRO_OSTYPE" ];
#    then
#        GCC_COMMAND="gcc -MM \$(CFLAGS) \$(CPPFLAGS) $* >> Makefile.depend 2> /dev/null"
#    else
#        GCC_COMMAND="gcc -MM \$(CFLAGS) \$(CPPFLAGS) $* 2> /dev/null | sed \"s|^.*\\.o:|\\\$\$(CSIRO_OSTYPE)\\/&|g\" >> Makefile.depend"
#    fi

    if [ -z "$CSIRO_OSTYPE" ];
    then
        GCC_COMMAND="gcc -MM \$(CFLAGS) \$(CPPFLAGS) $* >> Makefile.depend"
    else
        GCC_COMMAND="gcc -MM \$(CFLAGS) \$(CPPFLAGS) $* | sed \"s|^.*\\.o:|\\\$\$(CSIRO_OSTYPE)\\/&|g\" >> Makefile.depend"
    fi

    make_target "depend" "" \
        "@if ( gcc -v > /dev/null 2> /dev/null );\\" \
        "then\\" \
        "	echo \"Making dependencies...\";\\" \
        "	\$(RM) Makefile.depend;\\" \
        "	touch Makefile.depend;\\" \
        "	$GCC_COMMAND;\\" \
        "       chmod g+rw Makefile.depend;\\" \
        "	echo \"Making dependencies...done\";\\" \
        "else\\" \
        "	echo \"Making dependencies...ignored (no gcc)\";\\" \
        "fi"
}




make_action_target_hosts()
{
    if [ -n "$CSIRO_HOSTS" ];
    then
        make_target "hosts" "" \
            "@(\\" \
            "   if [ -n \"\$(CSIRO_HOSTS)\" ]; \\" \
			"   then \\" \
            "       targets=\`echo \$(MAKECMDGOALS) | sed 's|hosts||g'\`; \\" \
            "       for host in \$(CSIRO_HOSTS); \\" \
            "       do \\" \
            "           echo \"\"; \\" \
            "           echo \"Making\" \$\$targets on \$\$host; \\" \
            "           rsh \$\$host \"cd \$\$PWD; make \$\$targets\"; \\" \
            "       done; \\" \
            "       echo \"\"; \\" \
            "   fi \\" \
            ")"
    fi
}










#
#   FILE TARGET_FUNCTIONS
#


make_file_target_header()
{
    echo ""
    echo ""
    echo "#"
    echo "#	FILE / DIRECTORY TARGETS"
    echo "#"
}




make_file_target_ostype()
{
    if [ -n "$CSIRO_OSTYPE" ];
    then
        make_directory_target "\$(CSIRO_OSTYPE)"
    fi
}










#
#   OBJECT TARGET_FUNCTIONS
#


make_object_target_header()
{
    echo ""
    echo ""
    echo "#"
    echo "#	OBJECT TARGETS"
    echo "#"
}




make_object_targets()
{
    make_object_target_header
    make_obj_target_list $*
}










#
#   DEPENDECIES FUNCTIONS
#


make_depend_target_header()
{
    echo ""
    echo ""
    echo "#"
    echo "#	DEPENDENCIES"
    echo "#"
    echo "include Makefile.depend"
}




















#
#   LIBRARY MAKEFILE FUNCTIONS
#




make_lib_action_target_install()
{
    make_target "install" "all \$(LIB_INSTALL_DIR) \$(INC_INSTALL_DIR)" \
        "\$(CP) \$(LIB_FILE) \$(LIB_INSTALL_DIR)" \
        "\$(RM) -r \$(INC_INSTALL_DIR)" \
        "mkdir -p \$(INC_INSTALL_DIR)" \
        "\$(CP) *.h \$(INC_INSTALL_DIR)" \
        "chmod -R g+rwX \$(INC_INSTALL_DIR)"
}




make_lib_file_target_directories()
{
    make_directory_target "\$(LIB_INSTALL_DIR)"
    make_directory_target "\$(INC_INSTALL_DIR)"
}




make_lib_file_target()
{
#    make_target "\$(LIB_FILE)" "\$(LIB_OBJS)" \
#        "\$(LD) -r -o \$(LIB_FILE) \$(LIB_OBJS)"

    make_target "\$(LIB_FILE)" "\$(LIB_OBJS)" \
        "\$(LD) -r -o \$@ \$(LIB_OBJS)" \
        "chmod g+rw \$@"
}




make_lib_base()
{
    make_base_header
    
    echo ""
    echo "ROOT_DIR	= $ROOT_DIR"
    echo ""
    echo "LIB_NAME	= $LIB_NAME"
    echo "LIB_SOURCES	= $LIB_SOURCES"
    echo "LIB_CFLAGS	= $LIB_CFLAGS"
    echo "LIB_CPPFLAGS	= $LIB_CPPFLAGS"
    echo "LIB_FFLAGS	= $LIB_FFLAGS"
    echo "LIB_LDFLAGS	= $LIB_LDFLAGS"
    echo ""
    echo "#include \$(ROOT_DIR)/config/Makefile.lib"
}





make_lib_macros()
{
    make_macros_header

    set_obj_list $LIB_SOURCES

    echo ""
    echo "LIB_FILE	= ${OS_PREFIX}lib\$(LIB_NAME).a"
    echo "LIB_OBJS	=$obj_list"
    echo ""
    echo "LIB_INSTALL_DIR	= \$(ROOT_DIR)/lib$OS_SUFFIX"
    echo "INC_INSTALL_DIR	= \$(ROOT_DIR)/include/\$(LIB_NAME)"
    echo ""

    echo "CFLAGS		= \$(SYS_CFLAGS) \$(LIB_CFLAGS) \$(AUX_CFLAGS) -I.. -I\$(ROOT_DIR)/include"
    echo "CPPFLAGS	= \$(SYS_CPPFLAGS) \$(LIB_CPPFLAGS) \$(AUX_CPPFLAGS) -I.. -I\$(ROOT_DIR)/include"
    echo "FFLAGS		= \$(SYS_FFLAGS) \$(LIB_FFLAGS) \$(AUX_FFLAGS)"
    echo "LDFLAGS		= \$(SYS_LDFLAGS) \$(LIB_LDFLAGS) \$(AUX_LDFLAGS)"
}




make_lib_action_targets()
{
    make_action_target_header
    make_action_target_all "$OS_DEPEND\$(LIB_FILE)"
    make_action_target_debug
    make_action_target_optimise
    make_lib_action_target_install
    make_action_target_clean "\$(LIB_FILE) \$(LIB_OBJS)"
    make_action_target_make
    make_action_target_world
    make_action_target_depend "\$(LIB_SOURCES)"
    make_action_target_hosts
}




make_lib_file_targets()
{
    make_file_target_header
    make_file_target_ostype
    make_lib_file_target_directories
    make_lib_file_target
}



make_header_lib()
{
    make_lib_base

    make_macros_header
    echo ""
    echo "INSTALL_DIR	= \$(ROOT_DIR)/include/\$(LIB_NAME)"
    echo ""

    make_action_target_header
    make_target "all debug optimise fast clean depend" ""

#    make_target "install" "\$(INSTALL_DIR)" \
#        "\$(RM) -r \$(INSTALL_DIR)" \
#        "mkdir -p \$(INSTALL_DIR)" \
#        "\$(CP) *.h \$(INSTALL_DIR)"

    make_target "install" "\$(INSTALL_DIR)" \
        "\$(RM) -r \$(INSTALL_DIR)" \
        "mkdir -p \$(INSTALL_DIR)" \
        "\$(CP) *.h \$(INSTALL_DIR)" \
        "chmod -R g+rwX \$(INSTALL_DIR)"

    make_action_target_make
    make_action_target_world
    make_action_target_hosts

    make_file_target_header
    make_directory_target "\$(INSTALL_DIR)"
}




make_lib_makefile()
{
    make_init

    if [ -z "$LIB_SOURCES" ];
    then
        make_header_lib
    else
        make_lib_base
        make_lib_macros
        make_lib_action_targets
        make_lib_file_targets
        make_object_targets $LIB_SOURCES
        make_depend_target_header
    fi
}









#
#   ARCHIVE MAKEFILE FUNCTIONS
#


make_archive_file_target()
{
#    make_target "\$(LIB_FILE)" "$OS_DEPEND\$(LIB_OBJS)" \
#        "\$(RM) \$(LIB_FILE)" \
#        "\$(AR) -q \$(LIB_FILE) \$(LIB_OBJS)"

    make_target "\$(LIB_FILE)" "$OS_DEPEND\$(LIB_OBJS)" \
        "\$(RM) \$@" \
        "\$(AR) -q \$@ \$(LIB_OBJS)" \
        "chmod g+rw \$@"

#        "@for obj in \$(LIB_OBJS);\\" \
#        "do\\" \
#        "	echo \"\$(AR) -q \$(LIB_FILE) \$\$obj\";\\" \
#        "	\$(AR) -q \$(LIB_FILE) \$\$obj;\\" \
#        "done" \
}




make_archive_file_targets()
{
    make_file_target_header
    make_file_target_ostype
    make_lib_file_target_directories
    make_archive_file_target
}




make_archive_makefile()
{
    make_init

    if [ -z "$LIB_SOURCES" ];
    then
        make_header_lib
    else
        make_lib_base
        make_lib_macros
        make_lib_action_targets
        make_archive_file_targets
        make_object_targets $LIB_SOURCES
        make_depend_target_header
    fi
}




















#
#   BINARY MAKEFILE FUNCTIONS
#


make_bin_action_target_install()
{
    make_target "install" "all \$(INSTALL_DIR)" \
        "\$(CP) \$(BIN_FILE) \$(INSTALL_DIR)"
}




make_bin_file_target_directories()
{
    make_directory_target "\$(INSTALL_DIR)"
}




make_bin_file_target()
{
    set_compiler_list $BIN_SOURCES

#    make_target "\$(BIN_FILE)" "$OS_DEPEND\$(BIN_OBJS)" \
#        "$compiler \$(BIN_OBJS) \$(LDFLAGS) -o \$(BIN_FILE)"

    make_target "\$(BIN_FILE)" "$OS_DEPEND\$(BIN_OBJS)" \
        "$compiler \$(BIN_OBJS) \$(LDFLAGS) -o \$@" \
        "chmod g+rwx \$@"
}




make_bin_base()
{
    make_base_header
    
    echo ""
    echo "ROOT_DIR	= $ROOT_DIR"
    echo ""
    echo "BIN_NAME	= $BIN_NAME"
    echo "BIN_SOURCES	= $BIN_SOURCES"
    echo "BIN_CFLAGS	= $BIN_CFLAGS"
    echo "BIN_CPPFLAGS	= $BIN_CPPFLAGS"
    echo "BIN_FFLAGS	= $BIN_FFLAGS"
    echo "BIN_LDFLAGS	= $BIN_LDFLAGS"
    echo ""
    echo "#include \$(ROOT_DIR)/config/Makefile.bin"
}





make_bin_macros()
{
    make_macros_header

    set_obj_list $BIN_SOURCES
    set_bin $BIN_NAME

    echo ""
    echo "BIN_FILE	= $bin"
    echo "BIN_OBJS	=$obj_list"
    echo ""
    echo "INSTALL_DIR	= \$(ROOT_DIR)/bin$OS_SUFFIX"
    echo ""

    echo "CFLAGS		= \$(SYS_CFLAGS) \$(BIN_CFLAGS) \$(AUX_CFLAGS) -I\$(ROOT_DIR)/include"
    echo "CPPFLAGS	= \$(SYS_CPPFLAGS) \$(BIN_CPPFLAGS) \$(AUX_CPPFLAGS) -I\$(ROOT_DIR)/include"
    echo "FFLAGS		= \$(SYS_FFLAGS) \$(BIN_FFLAGS) \$(AUX_FFLAGS)"
    echo "LDFLAGS		= \$(SYS_LDFLAGS) \$(BIN_LDFLAGS) \$(AUX_LDFLAGS) -L\$(ROOT_DIR)/lib$OS_SUFFIX"
}




make_bin_action_targets()
{
    make_action_target_header
    make_action_target_all "$OS_DEPEND\$(BIN_FILE)"
    make_action_target_debug
    make_action_target_optimise
    make_bin_action_target_install
    make_action_target_clean "\$(BIN_FILE) \$(BIN_OBJS)"
    make_action_target_make
    make_action_target_world
    make_action_target_depend "\$(BIN_SOURCES)"
    make_action_target_hosts
}




make_bin_file_targets()
{
    make_file_target_header
    make_file_target_ostype
    make_bin_file_target_directories
    make_bin_file_target
}




make_bin_makefile()
{
    make_init
    make_bin_base
    make_bin_macros
    make_bin_action_targets
    make_bin_file_targets
    make_object_targets $BIN_SOURCES
    make_depend_target_header
}




















#
#   BINARIES MAKEFILE FUNCTIONS
#


make_bins_action_target_install()
{
    make_target "install" "all \$(INSTALL_DIR)" \
        "\$(CP) \$(BIN_FILES) \$(INSTALL_DIR)"
}




make_bins_file_target_directories()
{
    make_directory_target "\$(INSTALL_DIR)"
}




make_bins_file_target()
{
    make_bin_target_list "\$(AUX_OBJS)" $BIN_SOURCES
}




make_bins_base()
{
    make_base_header
    
    echo ""
    echo "ROOT_DIR	= $ROOT_DIR"
    echo ""
    echo "BIN_SOURCES	= $BIN_SOURCES"
    echo "AUX_SOURCES	= $AUX_SOURCES"
    echo "BIN_CFLAGS	= $BIN_CFLAGS"
    echo "BIN_CPPFLAGS	= $BIN_CPPFLAGS"
    echo "BIN_FFLAGS	= $BIN_FFLAGS"
    echo "BIN_LDFLAGS	= $BIN_LDFLAGS"
    echo ""
    echo "#include \$(ROOT_DIR)/config/Makefile.bins"
}





make_bins_macros()
{
    make_macros_header

    echo ""
    
    set_bin_list $BIN_SOURCES
    echo "BIN_FILES	=$bin_list"

    set_obj_list $BIN_SOURCES
    echo "BIN_OBJS	=$obj_list"

    set_obj_list $AUX_SOURCES
    echo "AUX_OBJS	=$obj_list"

    echo ""
    echo "INSTALL_DIR	= \$(ROOT_DIR)/bin$OS_SUFFIX"
    echo ""

    echo "CFLAGS		= \$(SYS_CFLAGS) \$(BIN_CFLAGS) \$(AUX_CFLAGS) -I\$(ROOT_DIR)/include"
    echo "CPPFLAGS	= \$(SYS_CPPFLAGS) \$(BIN_CPPFLAGS) \$(AUX_CPPFLAGS) -I\$(ROOT_DIR)/include"
    echo "FFLAGS		= \$(SYS_FFLAGS) \$(BIN_FFLAGS) \$(AUX_FFLAGS)"
    echo "LDFLAGS		= \$(SYS_LDFLAGS) \$(BIN_LDFLAGS) \$(AUX_LDFLAGS) -L\$(ROOT_DIR)/lib$OS_SUFFIX"
}




make_bins_action_targets()
{
    make_action_target_header
    make_action_target_all "$OS_DEPEND\$(BIN_FILES)"
    make_action_target_debug
    make_action_target_optimise
    make_bins_action_target_install
    make_action_target_clean "\$(BIN_FILES) \$(BIN_OBJS) \$(AUX_OBJS)"
    make_action_target_make
    make_action_target_world
    make_action_target_depend "\$(BIN_SOURCES)" "\$(AUX_SOURCES)"
    make_action_target_hosts
}




make_bins_file_targets()
{
    make_file_target_header
    make_file_target_ostype
    make_bins_file_target_directories
    make_bins_file_target
}




make_bins_makefile()
{
    make_init
    make_bins_base
    make_bins_macros
    make_bins_action_targets
    make_bins_file_targets
    make_object_targets $BIN_SOURCES $AUX_SOURCES
    make_depend_target_header
}
