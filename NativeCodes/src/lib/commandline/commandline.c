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
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "commandline/commandline.h"


/*  Useful macros for identifying options.  */

#define AN_OPTION(s)            (((s)[0] == '-') && ((s)[1] != '\0') && !isdigit((s)[1]))
#define A_SHORT_OPTION(s)       (((s)[0] == '-') && ((s)[1] != '-'))
#define A_SHORT_OPTION_LIST(s)  ((s)[2] != '\0')
#define MORE_OPTIONS(c,i)       (((c)->options[(i)].shortName != '\0') ||\
                                 ((c)->options[(i)].longName  != NULL))



static  CLOption   helpOption = 
    { '\0', "help", 0, 0, "print this usage information together with additional help" };



/*  Local Context manipulation functions   */

static  void        contextAllocate( CLContext* context, int argc );
static  void        contextDispose(  CLContext* context );
static  void        contextAppendArg(  CLContext* context, char* arg );
static  CLOption*   contextGetShortOption( CLContext* context, char  name );
static  CLOption*   contextGetLongOption(  CLContext* context, char* name );


/*  Local printing functions  */

static  void        printUsage(   CLContext* context );
static  void        printOptions( CLContext* context );
static  void        printOption(  CLOption* option );
static  void        printExtraHelp( CLContext* context );


/*  Local parse functions   */

static  int         parseArgs( CLContext* context, int argc, char* argv[] );

static  int         parseArg(           CLContext*  context,
                                        int         argc,
                                        char*       argv[],
                                        int*        i       );
                                
static  int         parseOption(        CLContext*  context,
                                        int         argc,
                                        char*       argv[],
                                        int*        i       );
                                    
static  int         parseShortOption(   CLContext*  context,
                                        int         argc,
                                        char*       argv[],
                                        int*        i       );
                                        
static  int         parseLongOption(    CLContext*  context,
                                        int         argc,
                                        char*       argv[],
                                        int*        i       );
                                        
static  int         parseOptionArgs(    CLContext*  context,
                                        CLOption*   option,
                                        int         argc,
                                        char*       argv[],
                                        int*        i       );


/*  Misc. functions */

static  char*       pathBasename( char* path );










/*----------------------------------

        IMPLEMENTATION

  ----------------------------------*/
  

void
clInit( CLContext* context, CLOption options[], char* argsDesc, char* extraHelp )
{
    if ( context != NULL )
    {
        context->argc       = 0;
        context->argv       = NULL;
        context->options    = options;
        context->argsDesc   = argsDesc;
        context->extraHelp  = extraHelp;
    }
}



void
clDone( CLContext* context )
{
    if ( context != NULL )
    {
        contextDispose(context);
    }
}




int
clParse( CLContext* context, int argc, char* argv[] )
{
    if ( (context == NULL) || (argv == NULL) )
    {
        return 0;
    }
    else
    {
        contextDispose(context);
        contextAllocate(context,argc);
        
        
        /*  The first command line arg is the program name so copy into the
            context structure.
        */
            
        contextAppendArg(context,argv[0]);
        
        
        if ( parseArgs(context,argc,argv) )
        {
            return 1;
        }
        else
        {
            clPrintUsage(context);
            
            return 0;
        }
    }
}




void
clPrintUsage( CLContext* context )
{
    if ( context != NULL )
    {
        printUsage(context);
        printOptions(context);
    }
}




void
clPrintExtraHelpUsage( CLContext* context )
{
    if ( context != NULL )
    {
        printUsage(context);
        printOptions(context);
        printExtraHelp(context);
    }
}










/*----------------------------------

        Context Functions

  ----------------------------------*/
  

/*  Allocate memory to hold the command line arguments after the options and
    their arguments have been extracted.
*/

void
contextAllocate( CLContext* context, int argc )
{
    context->argc = 0;
    context->argv = (char**)malloc( sizeof(char*) * argc );
}




/*  Free up the memory allocated for the command line arguments.    */

void
contextDispose( CLContext* context )
{
    if ( context->argv != NULL )
    {
        free(context->argv);
    }
}




/*  As the command line is parsed, if an argument is found which is not an
    option nor an argument of an option, then it is copied to the list of
    command line arguments 'argv' of the 'context'. This functions appends such
    an argument onto the list.
*/

void
contextAppendArg( CLContext* context, char* arg )
{
    context->argv[ context->argc++ ] = arg;
}




/*  Searches the options of a 'context' for the given short option.
    Returns a pointer to the option if found, NULL otherwise.
*/

CLOption*
contextGetShortOption( CLContext* context, char name )
{
    int i = 0;
    
    for ( i = 0; MORE_OPTIONS(context,i); i++ )
    {
        if ( context->options[i].shortName == name )
        {
            return &(context->options[i]);
        }
    }
    
    return NULL;
}




/*  Searches the options of a 'context' for the given long option.
    Returns a pointer to the option if found, NULL otherwise.
*/

CLOption*
contextGetLongOption( CLContext* context, char* name )
{
    int i = 0;
    
    for ( i = 0; MORE_OPTIONS(context,i); i++ )
    {
        if ( context->options[i].longName != NULL )
        {
            if ( strcmp(context->options[i].longName,name) == 0 )
            {
                return &(context->options[i]);
            }
        }
    }
    
    return NULL;
}










/*----------------------------------

        Print Functions

  ----------------------------------*/
  

/*  Prints the first part of the usage information.  This is the name of the
    program and the non-option arguments which are expected.
    eg.
        Usage:
                programname arg1 arg2 ...
*/

void
printUsage( CLContext* context )
{
    fprintf( stderr, "\n" );
    fprintf( stderr, "Usage:\n" );
    
    if ( context->argv != NULL )
    {
        if ( context->argv[0] == NULL )
        {
            fprintf( stderr, "\t(unknown)" );
        }
        else
        {
            fprintf( stderr, "\t%s", pathBasename(context->argv[0]) );
        }
    }

    if ( context->argsDesc != NULL )
    {
        fprintf( stderr, " %s", context->argsDesc );
    }

    fprintf( stderr, "\n" );

}




/*  Prints the usage information for the options.
    eg.
        Options:
                -o,--option arg1 arg2
                        description of option
                        arg1    - description of arg1
                        arg2    - description of arg2
*/

void
printOptions( CLContext* context )
{
    int i = 0;
    
    if ( context->options != NULL )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "Options:\n" );

        for ( i = 0; MORE_OPTIONS(context,i); i++ )
        {
            printOption( &(context->options[i]) );
        }
        
        if ( context->extraHelp != NULL )
        {
            printOption( &helpOption );
        }

        fprintf( stderr, "\n" );
    }
}





/*  Prints the usage information regarding the given option.    */

void
printOption( CLOption* option )
{
    int i = 0;
    
    if ( option->shortName == '\0' )
    {
        fprintf( stderr, "\t--%s", option->longName );
    }
    else if ( option->longName == NULL )
    {
        fprintf( stderr, "\t-%c", option->shortName );
    }
    else
    {
        fprintf( stderr, "\t-%c, --%s", option->shortName, option->longName );
    }
    
    
    for ( i = 0; i < option->argCount; i++ )
    {
        if ( option->args[i].name != NULL )
        {
            fprintf( stderr, " %s", option->args[i].name );
        }
    }
    
    fprintf( stderr, "\n" );

        
    if ( option->description != NULL )
    {
        fprintf( stderr, "\t\t%s\n", option->description );
    }
    
    
    for ( i = 0; i < option->argCount; i++ )
    {
        if ( (option->args[i].name != NULL) &&
             (option->args[i].description != NULL) )
        {
            fprintf( stderr, "\t\t%s\t- %s\n",
                                option->args[i].name,
                                option->args[i].description );
        }
    }
}




void
printExtraHelp( CLContext* context )
{
    if ( context->extraHelp != NULL )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "Additional Help\n" );
        fprintf( stderr, "%s\n", context->extraHelp );
    }
}










/*----------------------------------

        Parse Functions

  ----------------------------------*/
  

/*  Parses the command line arguments and returns 1 if successful, 0 otherwise.
*/

int
parseArgs( CLContext* context, int argc, char* argv[] )
{
    int i = 1;
    
    while ( i < argc )
    {
        if ( !parseArg(context,argc,argv,&i) )
        {
            return 0;
        }
        
        i++;
    }
    
    return 1;
}





/*  Parses the command line option given at index 'i'.
    Returns 1 if successful, 0 otherwise.
*/

int
parseArg( CLContext* context, int argc, char* argv[], int* i )
{
    if ( AN_OPTION(argv[*i]) )
    {
        return parseOption(context,argc,argv,i);
    }
    else
    {
        /*  The argument is not an option so append to the non-option command
            line arguments within the context.
        */
        
        contextAppendArg(context,argv[*i]);
        
        return 1;
    }
}




/*  Parses the option which is at index 'i' on the command line.    */

int
parseOption( CLContext* context, int argc, char* argv[], int* i )
{
    if ( A_SHORT_OPTION(argv[*i]) )
    {
        return parseShortOption(context,argc,argv,i);
    }
    else
    {
        return parseLongOption(context,argc,argv,i);
    }
}





/*  Parses the short option at index 'i' in the command line.
    The argument may be a list of short options such as '-ab' which is
    treated as if there were the two arguments '-a' and '-b'.  In this case
    neither option can have arguments.
*/

int
parseShortOption( CLContext* context, int argc, char* argv[], int* i )
{
    CLOption*   option;
    char*       p;
    

    /*  Step through each of the short (char) options in the list.  */
    
    for ( p = &(argv[*i][1]); *p != '\0'; p++ )
    {
        option = contextGetShortOption(context,*p);
        
        if ( option == NULL )
        {
            /*  No such short option. */
            
            if ( A_SHORT_OPTION_LIST(argv[*i]) )
            {
                fprintf( stderr, "\nError: unknown short option '-%c' within "
                                 "short option list '%s'\n",
                                 *p, argv[*i] );
            }
            else
            {
                fprintf( stderr, "\nError: unknown short option '-%c'\n", *p );
            }
            
            return 0;
        }
        else if ( (option->argCount > 0) && A_SHORT_OPTION_LIST(argv[*i]) )
        {
            /*  If the short option has arguments and is therefore not a flag,
                then it cannot be a member of a list of short options.
            */

            fprintf( stderr, "\nError: short option '-%c' requires %d "
                             "argument(s) and therefore cannot be part of"
                             "\n       the short option list '%s'\n",
                             *p, option->argCount, argv[*i] );
            
            return 0;
        }
        
        /*  Have recognised the option so now parse it arguments, if any.   */
        
        else if ( !parseOptionArgs(context,option,argc,argv,i) )
        {
            return 0;
        }
    }

    return 1;
}





/** Parses the long option at index 'i' in the command line.
    The option may be the builtin help option (CL_HELP_OPTION) in which case
    print the usage information and exit the program.
*/

int
parseLongOption( CLContext* context, int argc, char* argv[], int* i )
{
    CLOption*   option;
    

    option = contextGetLongOption(context,&(argv[*i][2]));

    if ( option == NULL )
    {
        /*  Option not found so check for builtin help option. */
        
        if ( strcmp(argv[*i], CL_HELP_OPTION) == 0 )
        {
            clPrintUsage(context);
            exit(0);
        }
        else
        {
            fprintf( stderr, "\nError: unknown long option '%s'\n", argv[*i] );

            return 0;
        }
    }
    else
    {
        return parseOptionArgs(context,option,argc,argv,i);
    }
}





/*  Parse the arguments of the option at index 'i' in the command line if it
    has any.  If the option does have arguments, then they must follow
    immediately after the option on the command line.
*/

int
parseOptionArgs(    CLContext*  context,
                    CLOption*   option,
                    int         argc,
                    char*       argv[],
                    int*        i       )
{
    int j = 0;
    int k = 0;
    
    
    /*  Flag that the option was encountered on the command line. */
    
    if ( option->address != NULL )
    {
        *(option->address) = 1;
    }

    
    /*  Check that the correct number of arguments have been supplied. */
    
    if ( (*i + option->argCount) >= argc )
    {
        fprintf( stderr, "\nError: insufficient number of arguments for "
                         "option '%s'."
                         "\n       It requires %d argument(s).\n",
                         argv[*i], option->argCount );
        
        return 0;
    }
    
    
    /*  Get the value for each argument. */
    
    for ( j = 0; j < option->argCount; j++ )
    {
        k = *i + j + 1;
        
        if ( AN_OPTION(argv[k]) )
        {
            /*  An option has been found instead of a required argument of the
                given option.
            */

            fprintf( stderr, "\nError: insufficient number of arguments for "
                             "option '%s'."
                             "\n       It requires %d argument(s).\n",
                             argv[*i], option->argCount );
        
            return 0;
        }
        
        if ( option->args[j].address != NULL )
        {
            /*  Assign the variable associated with this argument the value of 
                the command line argument.
            */
            
            switch (option->args[j].type)
            {
                case CL_INT_ARG:
                    *(int*)(option->args[j].address) = atoi(argv[k]);
                    break;

                case CL_FLOAT_ARG:
                    *(float*)(option->args[j].address) = atof(argv[k]);
                    break;

                case CL_STRING_ARG:
                    *(char**)(option->args[j].address) = argv[k];
                    break;
                
                default:
                    break;
            }
        }
    }
    
    *i += option->argCount;
    
    return 1;
}










/*----------------------------------

        Misc. Functions

  ----------------------------------*/
  

/*  Returns the basename of the given path.  That is the filename without any
    preceding directories.
    Simply returns a pointer to where the basename begins within the given
    path.
*/

char*
pathBasename( char* path )
{
    char*   start;
    
    for ( start = path; *path != '\0'; path++ )
    {
        if ( *path == '/' )
        {
            start = path + 1;
        }
    }
    
    return start;
}


