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
#include <math.h>
#include <float.h>
#include "metric.h"



/* module contains implementations of the distance function for the 
   Euclidean, Manhattan, and L_Infinity metric */

/****************************************************************************/

double EuclidDist2(Points,Index,NNQPoint,dimension)

	 /****************************************************************************/

	 /* returns the square of the Euclidean distance of Points[Index] and
		NNQPoint */

	 double **Points,*NNQPoint;
	 int Index,dimension;

{
  double dist,d;
  int j;

  dist=0;
  for (j=0; j<dimension; j++) {
	d=Points[Index][j]-NNQPoint[j];
	dist += d*d;
  }
  return(dist);

}
/****************************************************************************/

double ManhattDist(Points,Index,NNQPoint,dimension)

	 /****************************************************************************/

	 /* returns the Manhattan distance between Points[Index] and
		NNQPoint */

	 double **Points,*NNQPoint;
	 int Index,dimension;

{
  double dist,d;
  int j;

  dist=0;
  for (j=0; j<dimension; j++) {
	d=Points[Index][j]-NNQPoint[j];
	dist = dist + fabs(d);
  }

  return(dist);
}

/****************************************************************************/

double LInfinityDist(Points,Index,NNQPoint,dimension)

	 /****************************************************************************/

	 /* returns the square of the L_inf distance of Points[Index] and
		NNQPoint (farthest distance along all directions) 
	 */

	 double **Points,*NNQPoint;
	 int Index,dimension;

{
  double dist,d;
  int j;

  dist=-DBL_MAX;
  for (j=0; j<dimension; j++) {
	d=fabs(Points[Index][j]-NNQPoint[j]);
	if (dist < d) {
	  dist = d;
	}
  }
  return(dist);
}


/****************************************************************************/

double LGeneralDist(Points,Index,NNQPoint,dimension,MinkP)

	 /****************************************************************************/

	 /* returns the square of the Minkowski distance of Points[Index] and
		NNQPoint */

	 double **Points,*NNQPoint;
	 int Index,dimension,MinkP;

{
  double dist,d;
  int j;

  dist=0;
  for (j=0; j<dimension; j++) {
	dist +=fabs(pow(Points[Index][j]-NNQPoint[j],(double)MinkP));
  }
  return(pow(dist,1.0/(double)MinkP));
}















