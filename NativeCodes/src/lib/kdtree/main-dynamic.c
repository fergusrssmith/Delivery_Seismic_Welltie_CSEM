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
#include "dynamic-tree.h"
#include "NR/nrutil.h"

/*
  Test code for *extremely naive* dynamic kdtree library, where full set 
  of Points and the random path through them is known in advance. Neighbour
  searches can only return nodes that are already "visited".

  The basic idea is to build an initial kdtree based on the first, say, 100
  points in the random path, and rebuild the tree to double this in size when 
  the current tree is fully visited. The process is iterated. This seems 
  to preserve the logarithmic scaling of access times. 

 
  A typical compile might be
  % cc -g main-dynamic.c -I../../../include -lkdtree -lNR -lm -L../../../lib/Linux-i686-unknown -o Linux-i686-unknown/example

*/
// index pair
struct INTPAIR
{
  int indx;
  int val;
};

#define MAX_DOUBLINGS 10

int main(int argc, char *argv[])
{ 
  void formMultigridTraverseOrder(struct INTPAIR* mgIndices, int size);
  
  double **Points;
  int *knownFlag;
  double** met;
  double str,dt=0; 
  int i,k,j,k_d;
  optkdNode *OptkdRoot,*OptkdRoot2;
  int dimension;
  int numpoints=80;          // number of random points in unit square
  int nNeighborsDesired=20;     // number of nearest neighbors desired
  int *neigh, *neigh2;
  double told;
  struct timeval tv;
  struct timezone tz;

  int nKnown,nSeek,nTry,nFound;

  // seed random points from time of day
  gettimeofday(&tv,&tz);
  srand((int)tv.tv_sec);

  dimension = 2;

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
	  
	  if (0)
		{
		  Points[k][0] = (k<=1)?0:(Points[k-2][0] + (((k/2)!=17 && (k/2)!=31)?100:0));
		  Points[k][1] = (k%2==0)?0:100.0;
		  printf("Pt %d %g %g\n",k,Points[k][0],Points[k][1]);
		}
	}

  // Example of how to use library
  int *path = ivector(0,numpoints-1);
  // make traversal path.
  if (0)
	{
	  // direct path
	  for (k=0; k < numpoints; k++) path[k]=k;
	}
  else
	{
	  // multigrid path
	  struct INTPAIR* mgIndices = (struct INTPAIR*)malloc((size_t) ((numpoints * sizeof(struct INTPAIR))));
	  formMultigridTraverseOrder(mgIndices,numpoints);
	  for (k=0; k < numpoints; k++) path[k] = mgIndices[k].val;
	  free(mgIndices);
	}
  // build dynamic tree
  dynamicTree *dTree;

  dTree = buildDynamicTree(Points,numpoints,dimension,path);

  // set up metric: make a 2d metric corresp. to a 30 degree rotated ellipse metric
  //  setUp2DEuclidenGenMetric(dTree->OptkdRoot->Kglobals,30*3.14/180,2);
  setUp2DEuclidenGenMetric(dTree->OptkdRoot->Kglobals,0.0,1.0);
	  
  // loop through traversal path
  dt=0;
  for (k=0; k < numpoints; k++)
	{
	  gettimeofday(&tv,&tz);
	  told=(tv.tv_sec*1e6 + tv.tv_usec);

	  // number of neighbours to seek
	  nSeek = IMIN(k,nNeighborsDesired);   
	  // fetch neighbours
	  neigh = nextDynamicKdOptNNQuery(nSeek, EUCLIDEAN_GEN, dTree);
	  // show list of neighbours if desired		  
	  if (k>0 && 0)
		{
		  printf("neighC %d %g %g\n",k,Points[k][0],Points[k][1]);
		  for(j=1; j<=nSeek; j++) 
			{
			  printf("neigh %d %g %g\n",k,Points[neigh[j]][0],Points[neigh[j]][1]);
			}
		}
	  // time accesses
	  gettimeofday(&tv,&tz);
	  dt += (double)((tv.tv_sec*1e6+tv.tv_usec) - told);
	  if (neigh) free(neigh);
	}
  printf("mean time per query = %g uS\n",dt/(double)numpoints);
  gettimeofday(&tv,&tz);


  /* clean up */
  killDynamicTree(dTree);


  free_ivector(path,0,numpoints-1);
  free(Points[0]);
  free(Points);


  return EXIT_SUCCESS;
}



// reverse level based sort
int intPairCompare(const void *a, const void *b)
{
  struct INTPAIR *A,*B;
  A = (struct INTPAIR*)a;
  B = (struct INTPAIR*)b;
  return ( (*A).indx < (*B).indx ? 1 : ( (*A).indx > (*B).indx ? -1 : 0 ));
}

void formMultigridTraverseOrder(struct INTPAIR* mgIndices, int size)
{
  int j,k,v,sizeP,kstart,level;

  v = size-1;
  sizeP=0;
  while ((v >>= 1) > 0 ) sizeP++;
  // sizeP is now the smallest power of 2 s.t. 2^sizeP >= size
  kstart=(2 << sizeP);   // 2^sizeP
  for (j=0; j < size; j++) 
	{
	  k=kstart;
	  v=j+1;
	  level=sizeP+1;
	  // keep dividing by 2 to find what level this index is in a binary tree
	  while (v%k != 0) {k >>= 1;level--;}
	  mgIndices[j].indx = level;
	  mgIndices[j].val = j;
	  // printf("powers %d %d %d\n",kstart,j,level);
	}
  // sort index list into level-based order
  qsort(&mgIndices[0],(size_t)size,sizeof(struct INTPAIR),intPairCompare);
}

