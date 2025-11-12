/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 

#ifndef GEOSTATS_COMMANDLINE_H
#define GEOSTATS_COMMANDLINE_H


/**
    Organisation:   CSIRO Petroleum Resources
    Library:        geostats
    Author:         James Forbes
    Date:           September 1999
    
    Component:      Command Line
                    Processing of command line arguments.

    Description:
        
        The command line component of the geostats library is a collection of
        objects and functions for parsing command line arguments.
        
        TO BE FINISHED
            - mention short option lists.
    
    
    Example:
        See the program 'commandline.c' in the 'Examples' directory for an
        example of the usage of this component of the library.
*/





/*  The maximum number of arguments an option may have. */

#ifndef CL_MAX_ARGS
#define CL_MAX_ARGS         10
#endif



/*  The value type of an argument of a command line option.  */

#define CL_INT_ARG          1
#define CL_FLOAT_ARG        2
#define CL_STRING_ARG       3



/*  The builtin help option.  This does NOT need to be included in the option
    list.  It is recognised automatically and prints a usage message before
    terminating the program.  It can be overidden by supplying an option in the
    option list.
*/

#define CL_HELP_OPTION      "--help"



/*  The last element in the option list.  Rather than having to specify
    the number of options in the option list, the last option must be this
    one to indicate the end of the list.
*/

#define CL_END_OPTIONS      { '\0', NULL, 0, NULL, NULL }





/*  An argument of a command line option.

    Components:
        type
            the type of the value of the argument as specified above.
        address
            the address of the variable which is assigned the value of the
            argument.  It may be NULL in which case no value is assigned and the
            argument acts only as a placeholder.  The type of the variable
            should correspond to the type of the argument as that is the type of
            value which is assigned.
        name,
        description
            the name and brief description of the argument.  Used for printing
            usage information.
*/

struct CLArgStruct
{
    int     type;
    void*   address;
    char*   name;
    char*   description;
};

typedef struct CLArgStruct  CLArg;






/**
    A command line option.
    
    Components:
        shortName
            the single character short name version of the option.  It may be
            the null character '\0' in which case there is no short version of
            this option, only a long version.
        longName
            the long name version of the option.  It may be NULL in which case
            there is no long version of this option, only a short version.
        argCount
            the number of arguments the option has.
        address
            the address of a boolean (int) variable which is assigned the
            value of 1 if the option is present on the command line.  This
            variable should be initialised to 0 prior to parsing the
            command line.  The address may be NULL in which case no assignment
            is made, however, there would be no way of determining whether the
            option was present or not.  This would only be useful for obselete
            options which are retained for compatability.
        description
            a brief description of the option for printing usage information.
        args
            the list (array) of argument information for this option.
*/

struct CLOptionStruct
{
    char        shortName;
    char*       longName;
    int         argCount;
    int*        address;
    char*       description;
    CLArg       args[CL_MAX_ARGS];
};

typedef struct CLOptionStruct  CLOption;






/**
    Contains all the contextual information regarding the parse of the command
    line options.
    Once the command line has been successfully parsed, the 'argc' and 'argv'
    fields contain the command line arguments and their values as per the
    parameters of the 'main' function. The major difference being these are the
    command line arguments which are left after the options and their arguments
    have been extracted.
*/

struct CLContextStruct
{
    int         argc;
    char**      argv;
    CLOption*   options;
    char*       argsDesc;
    char*       extraHelp;
};

typedef struct CLContextStruct  CLContext;






/**
    INTERFACE:
    
        clInit
                This function initialises a context variable with the list of
                options prior to a parse being performed.
                The 'arguments' parameter is a list of the expected command
                line arguments in addition to the options and their arguments
                as specified in the 'options' list.  This is used for printing
                usage information.
        
        clDone
                This function simply 'cleans up' any memory allocated during
                the parse of the command line and extraction of command line
                arguments.
                It should only be called after the values of the command line
                arguments are nolonger required, but need not be called at
                all.  The memory allocated is minimal and since the command
                line is parsed only once, there is no need to deallocate the
                memory unless the operating system does not do so automatically
                upon exit.
            
        clParse
                This function parses the entire command line as specified by
                the 'options' field of the given 'context' structure and
                returns 1 if successful and 0 if not.
                
                If not successful, then an appropriate error message is
                displayed together with an automatically generated usage
                message.
                
                If the parse is successful, then all the variables associated
                with the options which are present on the command line are
                assigned the appropariate value, and all the remaining command
                line arguments not part of an option, are placed in the 'argv'
                array of the given 'context' structure.  At this point the
                'argc' and 'argv' fields of 'context' contain the reduced
                command line obtained from removing the options and their
                arguments from the given command line.
        
        clPrintUsage
                This function prints an automatically generated usage message to
                stderr.
*/

void    clInit(  CLContext* context, CLOption options[], char* argsDesc, char* extraHelp );
void    clDone(  CLContext* context );
int     clParse( CLContext* context, int argc, char* argv[] );
void    clPrintUsage( CLContext* context );
void    clPrintExtraHelpUsage( CLContext* context );




#endif


