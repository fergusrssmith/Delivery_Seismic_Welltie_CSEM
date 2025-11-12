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
#include "commandline/commandline.h"



/**
    Organisation:   CSIRO Petroleum Resources
    Library:        geostats
    Author:         James Forbes
    Date:           September 1999

    Program:        commandline
                    
                    A sample program showing the use of the command line
                    component of the geostats library.
*/



/*
    'boolean' variables to indicate whether an option is present on the
    command line.
*/

int s = 0;
int f = 0;
int m = 0;



/*
    The variables which are assigned the argument values of the command line
    option 'multi' if it is present.
*/
 
int      mi  = 0;
float    mf  = 0.0f;
char*    ms  = NULL;



/*
    The option array defining the structure of the command line options.
    That is:
        [-s] [--long] [-f|--flag] [-m|--multi int float string]
*/

CLOption options[] =
  {
    { 's',  NULL,   0,  &s,     "flag with no long option" },
    { '\0', "long", 0,  NULL,   "flag with no short option" },
    { 'f', "flag",  0,  &f,     "a flag with both long and short options" },
    { 'm', "multi", 3,  &m,     "an option with multiple (3) arguments",
        {
            { CL_INT_ARG,    &mi, "int",    "an integer argument" },
            { CL_FLOAT_ARG,  &mf, "float",  "a float argument" },
            { CL_STRING_ARG, &ms, "string", "a string argument" }
        }
    },
    CL_END_OPTIONS  /* defined as {0,NULL,..} to mark the end of the list   */
  };





int
main( int argc, char* argv[] )
{
    int         i = 0;
    CLContext   context;


    /*
        Initialise the command line context variable.
        The third 'string' parameter is printed with the usage message
        to indicate what arguments, excluding the options and their arguments,
        are required on the command line.
    */
    
    clInit( &context, options, "argument ...", NULL );
    
    
    
    /*  Parse the command line and if successful, print out its contents.
        If not successful the function prints out the usage message.
    */
    
    if ( clParse( &context, argc, argv ) )
    {
        /*  check which options were present.   */
        
        printf( "\nOptions:\n" );
        
        if ( s )
        {
            printf( "   -s was present.\n" );
        }
        
        if ( f )
        {
            printf( "   -f/--flag was present.\n" );
        }
        
        if ( m )
        {
            printf( "   -m/--multi was present. Arguments are:\n" );
            printf( "      int    arg1 = %d\n", mi );
            printf( "      float  arg2 = %f\n", mf );
            printf( "      string arg3 = '%s'\n", ms );
        }
        
        
        /*
            Print out the remaining command line arguments after the 
            options have been extracted.
            These are contained in the 'argc' and 'argv' fields of the
            CLContext variable.
        */

        printf( "\nArguments:\n" );
        
        for ( i = 0; i < context.argc; i++ )
        {
            printf( "   argv[%d] = %s\n", i, context.argv[i] );
        }
        
        printf( "\n" );
    
    
        /*  Clean up memory allocated by during the parse.  */

        clDone( &context );

        return 0;
    }
    else
    {
        /*  Unsuccessful parse of the command line. */
        
        return 1;
    }
}

