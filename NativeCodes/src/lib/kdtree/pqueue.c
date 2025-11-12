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
#include "metric.h"
#include "pqueue.h"
#include <float.h>

/*************************************************************************************/
/*  Code to implement the abstract data type priority queue for use in j nearest     */
/*  neighbor searching.  Actual implementation is done using heaps.                  */
/*                                                                                   */
/*  Adapted from Sedgewick's: Algorithms in C p. 148-160.                            */
/*************************************************************************************/

/* 
   The heap data structure consists of two priority queues.  One for the j-smallest
   distances encountered, one to keep the indexes into the points array of the  
   points corresponding to the j-smallest distances. 
*/



/*************************************************************************************/

void PQupheap_kd(DistArr,FoundArr,k)

	 /*************************************************************************************/

	 double *DistArr;  /* j-smallest distances encountered */

	 int *FoundArr,k;

{
  double v;
  int j;

  v=DistArr[k]; 
  DistArr[0] = DBL_MAX;
  j=FoundArr[k];

  while(DistArr[k/2] <= v) {
	DistArr[k] = DistArr[k/2];
	FoundArr[k] = FoundArr[k/2];
	k=k/2;
  }
  DistArr[k] = v;
  FoundArr[k] = j;
}

/*************************************************************************************/

void PQInsert_kd(double distance,int index,double *DistArr,int *FoundArr)

	 /*************************************************************************************/

{
  FoundArr[0]=FoundArr[0]+1;
  DistArr[FoundArr[0]] = distance;
  FoundArr[FoundArr[0]] = index;
  PQupheap_kd(DistArr,FoundArr,FoundArr[0]);
}



/*************************************************************************************/

void PQdownheap_kd(DistArr,FoundArr,k,index)

	 /*************************************************************************************/

	 double *DistArr;  /* j-smallest distances encountered */

	 int *FoundArr,k,index;

{

  int j,l,N;
  double v;

  v=DistArr[k];

  N = FoundArr[0];  /* tricky patch to maintain the data structure */
  FoundArr[0]=index;

  while (k <= N/2) {
	j=k+k;
	if (j < N && DistArr[j] <DistArr[j+1]) j++;
	if (v>=DistArr[j]) break;
	DistArr[k]=DistArr[j]; 
	FoundArr[k]=FoundArr[j];
	k=j;
  }

  DistArr[k] = v;
  FoundArr[k]= index;
  FoundArr[0]=N;  /* restore data struct */


}

/*************************************************************************************/

void PQreplace_kd(distance,DistArr,FoundArr,index)

/*************************************************************************************/

	 double *DistArr,distance;
	 int *FoundArr,index;

{
  DistArr[0]=distance;
  PQdownheap_kd(DistArr,FoundArr,0,index);
}
