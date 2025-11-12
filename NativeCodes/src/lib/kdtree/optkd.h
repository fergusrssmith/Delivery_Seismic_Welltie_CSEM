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
#ifndef __optkd_h__
#define __optkd_h__

#define BUCKETSIZE 50

# ifndef MIN
#  define	MIN(a,b)			(((a) < (b)) ? (a) : (b))
# endif

# ifndef MAX
#  define	MAX(a,b)			(((a) > (b)) ? (a) : (b))
# endif

typedef struct INDEX_DIST_PAIR
{
  int* indx;
  double* dist;
} INDEX_DIST_PAIR;


typedef struct optKglobals
{
  int *perm;        /* permutation array */
  int *optfound;
  int Count;
  double *nndist;
  double (*Distance)();
  /*
	Now, useful stuff for general Euclidean metrics. May be null (unused) in other cases
  */
  int dim;                    /* dimension */
  double** distance_metric;   /* symmetric matrix S for distance = dr'.S.dr measures if EUCLIDEAN_GEN is used */
  double* wksp1;              /* workspace */
  double* wksp2;              /* workspace */
} optKglobals;


typedef struct optkdnode {
  int bucket;  /* true if the node is a bucket node */
  int discrim; /* discriminator of node */
  double cutval; 
  struct optkdnode *loson, *hison;
  int lopt,hipt;
  optKglobals *Kglobals; /* ptrs to global variables */
} optkdNode;

void            Selection                       (double** , int , int , int , int, optkdNode*);
int             findmaxspread                   (int , int , int , double**, optkdNode*);
optkdNode*      BuildkdTree                     (double** , int , int , int, optKglobals*);
optkdNode*      BuildOptTree                    (double** , int , int);
void            rnnEuclidean                    (optkdNode* , double* , double** , int , int, optKglobals*);
void            setUp2DEuclidenGenMetric        (optKglobals* Kglobals,double angle_radians, double r);
void            setUp3DEuclidenGenMetric        (optKglobals* Kglobals,double a_hmax,double a_hmin,double a_vert, 
												 double ang1,double ang2,double ang3);
void            rnnEuclideanGen                 (optkdNode* , double* , double** , int , int, optKglobals*);
void            rnnGeneral                      (optkdNode* , double* , double** , int , int , int, optKglobals*);
int*            kdOptNNQuery                    (double** , int , double* , int , int , optkdNode* , int);
void            KillOptTree                     (optkdNode*);
void            KillMetricSpace                 (optKglobals*);
void            optInRegion                     (optkdNode* , int , double** , double**);
void            optAddRegion                    (optkdNode* , int , double** , double**);
int             optBoundsIntersectRegion        (double* , double** , int);
int             optBoundsContainsRegion         (double* , double** , int);
void            optRangeSearch                  (optkdNode* , double** , int , double** , double*);
int*            kdOptRectQuery                  (optkdNode* , double** , int , double**);
INDEX_DIST_PAIR* kdOptNNDistQuery               (double** , int , double* , int , int , optkdNode* , int);

#endif
