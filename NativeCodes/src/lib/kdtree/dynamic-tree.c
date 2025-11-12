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
#include <stdlib.h>
#include <string.h>
#include "nn.h"
#include "optkd.h"
#include "pqueue.h"
#include "dynamic-tree.h"
#include "NR/nrutil.h"
 
/*
  *Extremely naive* dynamic kdtree calls, where the full set of Points
  and a complete visiting path through them is known in advance. We
  start with a full set of points, and visit this in an order
  described by the unique indices in path[0....numPoints-1]. At each
  point, neighbour searches can only return nodes that have been already
  "visited" (rather like a sequential simulation routine in geostatistics).

  The basic idea is to build an initial kdtree based on the first, say, 100
  points in the random path, and rebuild the tree to double this in size when 
  the current tree is fully visited. The process is iterated. This seems 
  to preserve the logarithmic scaling of access times.

  Usage involves (1) building the tree using buildDynamicTree()
                 (2) set up some suitable distance metric
			     (3) loop through the nodes using nextDynamicKdOptNNQuery()
			     (4) Either repeat loop using a resetDynamicTree() before (3) again, or
			         clean up with killDynamicTree()
*/

/* 
   Build a new dynamic tree. This is analogous to a BuildOptTree()
   call, but we have to provide a complete path through the Points set
   path[0....numPoints-1] consisting of unique indexes from
   0....numPoints-1.  Make sure path[] is left intact as you operate
   on the tree.
*/
dynamicTree *buildDynamicTree(double** points,int numPoints,int dimension,int *path)
{
  int k,sz;
  dynamicTree *dTree;
  
  dTree = (dynamicTree *)malloc(sizeof(dynamicTree));

  dTree->dimension = dimension;
  dTree->numPoints = numPoints;
  dTree->originalIndexes = path;
  dTree->shuffledPoints = dmatrix(0,numPoints-1,0,dimension-1);
  /* copy points into shuffled array in path's access order */
  sz = dimension*sizeof(double);
  for (k=0; k<numPoints; k++)
	{
	  memcpy(&(dTree->shuffledPoints[k][0]),&(points[path[k]][0]),sz);
	}
  // build initial tree
  dTree->currentSize = MIN(DYNAMIC_TREE_START_SIZE, numPoints);
  dTree->OptkdRoot = BuildOptTree(dTree->shuffledPoints, dTree->currentSize, dTree->dimension);
  dTree->currentIndex=0;
  dTree->knownFlag = calloc(numPoints, sizeof(int));
  dTree->distance_metric = NULL;
  dTree->wksp1 = NULL;
  dTree->wksp2 = NULL;

  return dTree;
}

/* Clean up after finished with tree and neighbour searches */
void killDynamicTree(dynamicTree *dTree)
{
  if (dTree->OptkdRoot != NULL)
	{
	  KillOptTree(dTree->OptkdRoot);
	}
  if (dTree->shuffledPoints != NULL) 
	{
	  free_dmatrix(dTree->shuffledPoints, 0,dTree->numPoints-1,0,dTree->dimension-1);
	}
  if (dTree->knownFlag != NULL) 
	{
	  free(dTree->knownFlag);
	  dTree->knownFlag = NULL;
	}

  if (dTree->distance_metric != NULL)
	{
	  free_dmatrix(dTree->distance_metric,0,dTree->dim-1,0,dTree->dim-1);
	  dTree->distance_metric = NULL;
	}
  if (dTree->wksp1 != NULL)
	{
	  free_dvector(dTree->wksp1,0,dTree->dim-1);
	  dTree->wksp1 = NULL;
	}
  if (dTree->wksp2 != NULL)
	{
	  free_dvector(dTree->wksp2,0,dTree->dim-1);
	  dTree->wksp2 = NULL;
	}
  if (dTree != NULL)
	{
	  free(dTree);
	  dTree = NULL;
	}
}

/* 
   Return a list indx[1....MIN(dTree->currentIndex, nNeighborsDesired)]
   of indexes of the nearest neighbours to the current site from those already visited
   On the first call, there are no neighbors so we return NULL.
*/
int* nextDynamicKdOptNNQuery(int nNeighborsDesired,int Metric,dynamicTree *dTree)
{
  int j,k,nFound,nSeek,nTry;
  int *neigh=NULL;
  INDEX_DIST_PAIR *neighPairs;
  INDXPAIR *desiredSet;

  if (dTree->currentIndex == 0) 
	{
	  dTree->knownFlag[dTree->currentIndex] = 1;
	  dTree->currentIndex++;
	  return NULL;
	}
  // if path complete bail out with null
  if (dTree->currentIndex == dTree->numPoints) return NULL;

  // update current tree if necessary
  if ((dTree->currentIndex == dTree->currentSize) && (dTree->currentSize != dTree->numPoints))
	{
	  // copy metrics before rebuilding kdtree
	  copyMetrics(dTree);
	  KillOptTree(dTree->OptkdRoot);
	  // make tree of double size
	  dTree->currentSize = IMIN(2*dTree->currentSize, dTree->numPoints);
	  dTree->OptkdRoot = BuildOptTree(dTree->shuffledPoints, dTree->currentSize, dTree->dimension);
	  // restore the metrics
	  restoreMetrics(dTree);
	}

  // we *must* be able to find this many known points, if size(tree) > nNeighborsDesired
  nSeek = MIN(dTree->currentIndex, nNeighborsDesired); 
  // start size for search
  // random studies indicate double the nNeighborsDesired is about right.
  nTry = MIN(dTree->currentSize, 2*nNeighborsDesired);
  for(;;)
	{
	  // by definition, query point is dTree->shuffledPoints[dTree->currentIndex]
	  neighPairs = kdOptNNDistQuery(dTree->shuffledPoints, dTree->dimension, dTree->shuffledPoints[dTree->currentIndex], 
									nTry, Metric, dTree->OptkdRoot, 1);
	  for (nFound=0,j=1; j<=nTry; j++) if (dTree->knownFlag[neighPairs->indx[j]]) nFound++;
	  if (nFound >= nSeek) 
		{
		  break;    // range search was big enough
		}
	  else
		{
		  nTry = IMIN(2*nTry,dTree->currentSize); // double size and try again
		  free(neighPairs->indx); neighPairs->indx = NULL;
		  free(neighPairs->dist); neighPairs->dist = NULL;
		  free(neighPairs); neighPairs = NULL;
		}
	}


  // neighPairs now contains at least as many points as we need. 
  // Extract the nSeek nearest ones
  desiredSet = (INDXPAIR*) malloc(nFound*sizeof(INDXPAIR));
  for (j=1,k=0; j<=nTry; j++) 
	{
	  // accumulate subset of known points
	  if (dTree->knownFlag[neighPairs->indx[j]]) 
		{
		  desiredSet[k].indx = neighPairs->indx[j];
		  desiredSet[k++].val = neighPairs->dist[j];
		} 
	}
  // we no longer need neigh and its members
  free(neighPairs->indx); neighPairs->indx = NULL;
  free(neighPairs->dist); neighPairs->dist = NULL;
  free(neighPairs); neighPairs = NULL;


  // sort the desired set in order of increasing distance
  qsort(desiredSet,nFound,sizeof(INDXPAIR),distCompare);
  // and copy back the array of neighbour indexes. The indexes apply to the shuffled set, so
  // map these back to the original ordering
  neigh = (int*) malloc((nSeek+1)*sizeof(int));
  for (j=0; j<nSeek; j++) neigh[j+1] = dTree->originalIndexes[desiredSet[j].indx];
  free(desiredSet);
		  
  dTree->knownFlag[dTree->currentIndex] = 1;

  dTree->currentIndex++;

  // sorted indices neigh[1....MIN(dTree->currentIndex, nNeighborsDesired)] are the nighbour indexes
  return neigh;  
}

/* 
   Reset tree iterator so a loop over nodes in the path can be repeated.
   Every complete loop over the numPoints nodes in the path nextDynamicKdOptNNQuery(...)
   must be followed by this call if you wish to repeat the traversal of the path, fetching
   nearest neighbor sets as you go.
*/
void resetDynamicTree(dynamicTree *dTree)
{
  // copy metrics before rebuilding kdtree
  copyMetrics(dTree);
  KillOptTree(dTree->OptkdRoot);
  // restore tree to initial size
  dTree->currentSize = MIN(DYNAMIC_TREE_START_SIZE, dTree->numPoints);
  // rebuild
  dTree->OptkdRoot = BuildOptTree(dTree->shuffledPoints, dTree->currentSize, dTree->dimension);
  // restore the metrics
  restoreMetrics(dTree);
  // reset the counter & flags
  dTree->currentIndex =0;
  memset(dTree->knownFlag,0,dTree->numPoints*sizeof(int));
}


/*
  Comparator for distance based sort of an INDXPAIR array
*/
int distCompare(const void *a, const void *b)
{
  INDXPAIR *A,*B;
  A = (INDXPAIR*)a;
  B = (INDXPAIR*)b;
  return ( (*A).val > (*B).val ? 1 : ( (*A).val < (*B).val ? -1 : 0 ));
}

void copyMetrics(dynamicTree *dTree)
{
  dTree->dim = dTree->OptkdRoot->Kglobals->dim;
  if (dTree->distance_metric == NULL)
	{
	  dTree->distance_metric = dmatrix(0,dTree->dim-1,0,dTree->dim-1);
	  memcpy(&(dTree->distance_metric[0][0]), &(dTree->OptkdRoot->Kglobals->distance_metric[0][0]), 
			 dTree->dim * dTree->dim * sizeof(double));
	}
  if (dTree->wksp1 == NULL) 
	{
	  dTree->wksp1 = dvector(0,dTree->dim-1);
	  memcpy(dTree->wksp1,dTree->OptkdRoot->Kglobals->wksp1,dTree->dim * sizeof(double));
	}
  if (dTree->wksp2 == NULL && dTree->OptkdRoot->Kglobals->wksp2 != NULL) 
	{
	  dTree->wksp2 = dvector(0,dTree->dim-1);
	  memcpy(dTree->wksp2,dTree->OptkdRoot->Kglobals->wksp2,dTree->dim * sizeof(double));
	}
}

void restoreMetrics(dynamicTree *dTree)
{
  dTree->OptkdRoot->Kglobals->dim = dTree->dim;

  dTree->OptkdRoot->Kglobals->distance_metric = dmatrix(0,dTree->dim-1,0,dTree->dim-1);
  memcpy(&(dTree->OptkdRoot->Kglobals->distance_metric[0][0]), &(dTree->distance_metric[0][0]), 
		 dTree->dim * dTree->dim * sizeof(double));

  dTree->OptkdRoot->Kglobals->wksp1 = dvector(0,dTree->dim-1);
  memcpy(dTree->OptkdRoot->Kglobals->wksp1,dTree->wksp1,dTree->dim * sizeof(double));

  if (dTree->wksp2 != NULL)
	{
	  dTree->OptkdRoot->Kglobals->wksp2 = dvector(0,dTree->dim-1);
	  memcpy(dTree->OptkdRoot->Kglobals->wksp2,dTree->wksp2,dTree->dim * sizeof(double));
	}
}
