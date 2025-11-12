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
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>

	
#include "nn.h"
#include "optkd.h"
#include "pqueue.h"
#include "sproullkd.h"

#include "NR/nrutil.h"

/*
  Example driver for kdtree library: creates a bunch of random points
  in 2 dimensions, sets up the kdtree, then does a query to get a list
  of the nearest nNeighborsDesired points. Time to produce the
  points is shown. To graph the output using GNU plotutils ( a handy thing), try

  % example | awk 'NR>1{print $0}'| graph -TX -C -m-1 -S 2

 
  A typical compile might be
  % cc -g main.c -I../../../include -lkdtree -lNR -lm -L../../../lib/Linux-i686-unknown -o Linux-i686-unknown/example

*/

int main(int argc, char *argv[])
{ 
  double **Points,**Points2;
  double** met;
  double str,dt; 
  int k,j;
  optkdNode *OptkdRoot,*OptkdRoot2;
  int dimension;
  int numpoints=10000;           // number of random points in unit square
  int nNeighborsDesired=25;      // number of nearest neighbors desired
  int *neigh, *neigh2;
  long told;
  struct timeval tv;
  struct timezone tz;

  int DO_3D_CASE = 0;              // Set true (1) to try general 3D case
  int DO_SPROULL = 0;              // Set true (1) to try general 3D case with sproull
  int TIME3D = 1;              // Set true (1) to time 3D case
  int DO_TWO_TREES =  0;           // 2D: set to true (1) to do two-tree test

  // seed random points from time of day
  gettimeofday(&tv,&tz);
  srand((int)tv.tv_sec);

  if (!DO_3D_CASE)
	{
	  dimension = 2;
	  printf("2D case\n");
	  if (!DO_TWO_TREES)
		{
		  // make set of random points in unit square
		  Points = (double **)malloc((numpoints)*sizeof(double *));
		  Points[0] = (double *)malloc((numpoints*dimension)*sizeof(double));
		  for (k=0; k < numpoints; k++) 
			{
			  if (k>0) Points[k] =  Points[k-1] + dimension;
			  for (j=0;j<dimension;j++)
				{
				  Points[k][j]=(rand()/(double)RAND_MAX);
				}
			}
		  // set up kd-tree
		  OptkdRoot = BuildOptTree(Points,numpoints,dimension);
		  // make a 2d metric corresp. to a 30 degree rotated ellipse metric
		  setUp2DEuclidenGenMetric(OptkdRoot->Kglobals,30*3.14/180,8);

		  // get nearest neighbor list to some random points Point[j]. Time and average this.
		  double tAv=0;
		  for (k=0; k < MIN(10000,numpoints); k++) 
			{
			  gettimeofday(&tv,&tz);
			  dt=(double)((tv.tv_sec*1e6+tv.tv_usec));
			  neigh = kdOptNNQuery(Points,dimension,
								   Points[(int)(numpoints * rand()/(double)RAND_MAX)],
								   nNeighborsDesired,EUCLIDEAN_GEN,OptkdRoot,1);
			  gettimeofday(&tv,&tz);
			  dt=(double)((tv.tv_sec*1e6+tv.tv_usec)) - dt;
			  tAv += (double)dt;
			  free(neigh);
			}
		  printf("Mean access time for search on %d points for %d nbrs = %f uS\n",
				 numpoints,nNeighborsDesired,(double)tAv/(double)k);

#if 0
		  // For graphing, print points then neighbor list
		  for (k=0; k < numpoints; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points[k][j]);
				}
			  printf("\n");
			}
		  printf("\n");
		  // neighbour list
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points[neigh[k]][j]);
				}
			  printf("\n");
			}
		  printf("\n\n");
		  // centre point
		  for (j=0;j<dimension;j++)
			{
			  printf("%f ",Points[0][j]);
			}
		  printf("\n");

#endif		  

		  /* clean up */
		  KillOptTree(OptkdRoot);
		}
	  else
		{
		  // multiple tree test (thread/user safety)
		  // make 2 sets of random points in unit square
		  Points = (double **)malloc((numpoints)*sizeof(double *));
		  Points[0] = (double *)malloc((numpoints*dimension)*sizeof(double));
		  Points2 = (double **)malloc((numpoints)*sizeof(double *));
		  Points2[0] = (double *)malloc((numpoints*dimension)*sizeof(double));
		  for (k=0; k < numpoints; k++) 
			{
			  if (k>0)
				{
				  Points[k] =  Points[k-1] + dimension;
				  Points2[k] =  Points2[k-1] + dimension;
				}
			  for (j=0;j<dimension;j++)
				{
				  Points[k][j]=(rand()/(double)RAND_MAX);
				  Points2[k][j]=(rand()/(double)RAND_MAX)+1.05;
				}
			}
		  // set up 2 kd-trees
		  OptkdRoot = BuildOptTree(Points,numpoints,dimension);
		  OptkdRoot2 = BuildOptTree(Points2,numpoints,dimension);
		  // make two different metrics
		  // make a 2d metric corresp. to a 30 degree rotated ellipse metric
		  // up
		  setUp2DEuclidenGenMetric(OptkdRoot->Kglobals,30*3.14/180,10);
		  // down
		  setUp2DEuclidenGenMetric(OptkdRoot2->Kglobals,-30*3.14/180,4);

		  // get nearest neighbor list to first point Point[0]. Time this.
		  gettimeofday(&tv,&tz);
		  told=tv.tv_usec;
		  neigh = kdOptNNQuery(Points,dimension,Points[0],nNeighborsDesired,EUCLIDEAN_GEN,OptkdRoot,1);
		  neigh2 = kdOptNNQuery(Points2,dimension,Points2[0],nNeighborsDesired,EUCLIDEAN_GEN,OptkdRoot2,1);
		  gettimeofday(&tv,&tz);
		  told=tv.tv_usec-told;
		  printf("time = %ld uS\n",told);

		  // For graphing, print points then neighbor list
		  for (k=0; k < numpoints; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points[k][j]);
				}
			  printf("\n");
			}
		  printf("\n");
		  for (k=0; k < numpoints; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points2[k][j]);
				}
			  printf("\n");
			}
		  printf("\n");

		  // neighbour list
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points[neigh[k]][j]);
				}
			  printf("\n");
			}
		  printf("\n");
		  for (k=1; k <= nNeighborsDesired; k++) 
			{
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points2[neigh2[k]][j]);
				}
			  printf("\n");
			}
		  printf("\n");

		  // centre point
		  for (j=0;j<dimension;j++)
			{
			  printf("%f ",Points[0][j]);
			}
		  printf("\n\n");
		  for (j=0;j<dimension;j++)
			{
			  printf("%f ",Points2[0][j]);
			}
		  printf("\n");
		}
	}
  else
	{
	  // 3D case. Do single data set only
	  // make set of random points in unit square
	  dimension = 3;
	  int zChop=30;
	  numpoints = 10000;
	  nNeighborsDesired=25;
	  Points = (double **)malloc((numpoints)*sizeof(double *));
	  Points[0] = (double *)malloc((numpoints*dimension)*sizeof(double));
	  for (k=0; k < numpoints; k++) 
		{
		  if (k>0) Points[k] =  Points[k-1] + dimension;
		  for (j=0;j<dimension;j++)
			{
			  // random
			  Points[k][j]=(rand()/((double)RAND_MAX+1));
			}
		  // Points[k][dimension] = ((int)(zChop*Points[k][j]))/(float)zChop;
		}

	  if (!DO_SPROULL)
		{
		  // set up kd-tree
		  OptkdRoot = BuildOptTree(Points,numpoints,dimension);
		  // make a 3d metric corresp. to a 30 degree azimuth and 20 degree dip
		  setUp3DEuclidenGenMetric(OptkdRoot->Kglobals,1.0,1.0,1.0,0,0,0);
		  //	  setUp3DEuclidenGenMetric(OptkdRoot->Kglobals,1.0,1.0,1.0,0,0,0);
		  // get nearest neighbor list to first point Point[0]

		  neigh = kdOptNNQuery(Points,dimension,Points[0],nNeighborsDesired,EUCLIDEAN_GEN,OptkdRoot,1);

		  if (TIME3D)	  
			{ 
			  long tTot=0,nSamp;
			  nSamp= 20;
			  for (k=0; k < nSamp; k++) 
				{
				  gettimeofday(&tv,&tz);
				  told=tv.tv_usec;
				  neigh = kdOptNNQuery(Points,dimension,Points[k],nNeighborsDesired,EUCLIDEAN_GEN,OptkdRoot,1);
				  //neigh = kdOptNNQuery(Points,dimension,Points[k],nNeighborsDesired,EUCLIDEAN,OptkdRoot,1);
				  gettimeofday(&tv,&tz);
				  told=tv.tv_usec-told;
				  tTot += told;
				}
			  printf("Mean access time for search on %d points for %d nbrs = %f uS\n",numpoints,nNeighborsDesired,tTot/(float)nSamp);
			}
		  if (!TIME3D)
			{
			  // For graphing, print points then neighbor list
			  for (k=0; k < numpoints; k++) 
				{
				  for (j=0;j<dimension;j++)
					{
					  printf("%f ",Points[k][j]);
					}
				  printf("\n");
				}
			  printf("\n");
			  // neighbour list
			  for (k=1; k <= nNeighborsDesired; k++) 
				{
				  for (j=0;j<dimension;j++)
					{
					  printf("%f ",Points[neigh[k]][j]);
					}
				  printf("\n");
				}
			  printf("\n\n");
			  // centre point
			  for (j=0;j<dimension;j++)
				{
				  printf("%f ",Points[0][j]);
				}
			  printf("\n");



			  /* Suitable output for trip or chull */
			  /* view it like this...
				 example | awk '/trip/{print $2}' | a2b n1=1 | trip n1=40 n2=40 n3=40
			  */

			  {
				int Ncube=40,Nchull=300;
				int iX,iY,iZ;
				double ***c;
				c = d3tensor(0,Ncube-1,0,Ncube-1,0,Ncube-1);
				// zero c
				memset(&c[0][0][0], 0,Ncube*Ncube*Ncube*sizeof(double));

				// set points in rasterized cube vertices to 1
				for (k=1; k <= nNeighborsDesired; k++) 
				  {
					iX = (int)Ncube*Points[neigh[k]][0];
					iY = (int)Ncube*Points[neigh[k]][1];
					iZ = (int)Ncube*Points[neigh[k]][2];
					c[iX][iY][iZ]=1;
					iX = (int)Nchull*Points[neigh[k]][0];
					iY = (int)Nchull*Points[neigh[k]][1];
					iZ = (int)Nchull*Points[neigh[k]][2];
					printf("chull %d %d %d\n",iX,iY,iZ);
				  }
				for (k=0;k<Ncube*Ncube*Ncube;k++)
				  {
					printf("trip %d\n",(int)c[0][0][k]);
				  }
			  }
			}
		  /* clean up */
		  KillOptTree(OptkdRoot);
		}
	  else
		{
		  // Sproull 3D
		  // set up kd-tree
		  // Tests indicate this is quite slow for small numbers.....

		  sproullkdNode *SprkdRoot;
		  SprkdRoot = BuildSproullkdTree(Points,numpoints,dimension);

		  // get nearest neighbor list to first point Point[0]
		  neigh = kdSproullNNQuery(Points,dimension,Points[0],nNeighborsDesired,EUCLIDEAN,SprkdRoot,1);

		  if (TIME3D)	  
			{ 
			  long tTot=0,nSamp;
			  nSamp= 20;
			  for (k=0; k < nSamp; k++) 
				{
				  gettimeofday(&tv,&tz);
				  told=tv.tv_usec;
				  neigh = kdSproullNNQuery(Points,dimension,Points[k],nNeighborsDesired,EUCLIDEAN,SprkdRoot,1);
				  gettimeofday(&tv,&tz);
				  told=tv.tv_usec-told;
				  tTot += told;
				}
			  printf("Mean access time for sproull search on %d points for %d nbrs = %f uS\n",numpoints,nNeighborsDesired,tTot/(float)nSamp);
			}
		}
	}
  return EXIT_SUCCESS;
}
