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
#include <unistd.h>
#include <time.h>
#include "commandline/commandline.h"
#include "progress/progress.h"


extern	void		srand48( long int );
extern	long int	lrand48();



/**
    Organisation:   CSIRO Petroleum Resources
    Program:        taskgenerate
    Author:         James Forbes
    Date:           September 1999
    
    Description:    A sample program showing the use of the progress
                    component of the geostats library.  It generates a given
                    number of tasks and randomly updates their progress until
                    completion.  The program 'progress' can be used to monitor
                    the progress of the tasks.
*/



#define OUTPUT(statements)      if ( !silent ) { statements }



int taskCount   = 1;
int silent      = 0;


CLOption options[] = 
    {
        { 't', "tasks", 1, NULL, "number of tasks to generate (default is 1)",
            {
                { CL_INT_ARG, &taskCount, "n", NULL }
            },
        },
        { 's', "silent", 0, &silent, "operate in silent mode (ie. no output)" },
        CL_END_OPTIONS
    };





void
updateTask()
{
    int units = 0;


    /*  Determine a random number between 0 and 9 of work units completed. */

    units = lrand48() % 10;

    OUTPUT(
        printf( "\r   updating by: %d  ", units );
        fflush(stdout);
    )


    /*  Update the progress information.    */

    progressUpdateTaskBy(units);
}





void
generateTask( int id )
{
    char    taskName[100];


    OUTPUT( printf( "Generating task %d\n", id ); )
    sprintf( taskName, "Sample Task %d", id );


    /*  Indicate a new task has begun.
        It requires 100 work units to complete.
    */

    progressNewTask( taskName, 100 );


    /*  Continue updating task until it is completed.   */

    while ( progressTaskStatus() != PROGRESS_COMPLETED )
    {
        /*  Pause briefly for better viewing of progress changes. */

        sleep(1);
        
        updateTask();
    }

    OUTPUT( printf( "\n" ); )
}






void
generateTasks()
{
    int     i = 0;
    char    input[2];


    OUTPUT(
        printf( "\n" );
        printf( "Sample Progress Task Generator\n" );
        printf( "Process id: %d\n", (int)getpid() );
        printf( "\n" );
        printf( "Press return to continue..." );
        fflush(stdout);
        read( 0, input, 1 );
        printf( "\n" );
    )


    /*  Initialisation.    */

    progressInit();
    srand48(time(NULL));


    /*  Sequentially generate the tasks.    */
    
    for ( i = 1; i <= taskCount; i++ )
    {
        generateTask(i);

        /*  Pause between tasks for better viewing of progress changes. */

        sleep(2);
    }

    OUTPUT( printf( "\n" ); )

    progressDone();
}





int
main( int argc, char* argv[] )
{
    CLContext   context;


    clInit( &context, options, NULL, NULL );
    
    if ( clParse( &context, argc, argv ) )
    {
        generateTasks();
        clDone(&context);
        
        return 0;
    }
    else
    {
        clDone(&context);

        return 1;
    }
}


