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
#ifndef __dynamic_tree_h__
# define __dynamic_tree_h__

#define DYNAMIC_TREE_START_SIZE 50

typedef struct INDXPAIR
{
  int indx;
  double val;
} INDXPAIR;

typedef struct dynamicTree
{
  double **shuffledPoints;
  int *knownFlag;
  int dimension;
  int numPoints;
  int currentSize;
  int currentIndex;
  int *originalIndexes;
  optkdNode *OptkdRoot;
  /*
	Useful stuff for general Euclidean metrics. Copy of equivalent optKglobalsvariables, since tree may be repeatedly
	destroyed and rebuilt, and we want to preserve the metrics
  */
  int dim;                    // dimension
  double** distance_metric;   // symmetric matrix S for distance = dr'.S.dr measures if EUCLIDEAN_GEN is used
  double* wksp1;              // workspace
  double* wksp2;              // workspace
} dynamicTree;

dynamicTree *buildDynamicTree(double** points,int numPoints,int dimension,int *path);
void killDynamicTree(dynamicTree *dTree);
int* nextDynamicKdOptNNQuery(int nNeighborsDesired,int Metric,dynamicTree *dTree);
int distCompare(const void *a, const void *b);
void copyMetrics(dynamicTree *dTree);
void restoreMetrics(dynamicTree *dTree);
void resetDynamicTree(dynamicTree *dTree);

#endif
