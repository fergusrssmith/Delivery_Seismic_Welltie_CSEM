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
#ifndef __eclio_h__
#define __eclio_h__

#include <stdio.h>

/* Defines */

#define JG_FEATURES 1

/*************************************************************************/
#define D_BIG 1.0E15
#define I_BIG 32000
#define D_LIL 1.0E-15
#define D_TOL 1.E-9

#define N_INC_MAX 10

#define TOKEN_LENGTH 64
#ifndef NAME_LENGTH
#define NAME_LENGTH 128
#endif
#ifndef LINE_LENGTH
#define LINE_LENGTH  256
#endif

#define D_MISSING -999.25
#define I_MISSING -999

#define GENERIC_POINTER void*
#define FREE_ARG char*

#define ECL_WHITESPACE "\n\t\v\b\r\f "
#define ECL_TOKENS "\n\t\v\b\r\f\'\" "
#define ECL_TRINKETS "\n\t\v\b\r\f\'\" /"
#define ONE_LINE_RECORDS " ECHO NOECHO  "
#define MULTI_SLASH_RECORDS " COORDSYS "

/* standard cliches */
/*************************************************************************/
#ifndef  MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef  MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#define APPROX(x,y,tol) (fabs((double)(x)-(double)(y))/(0.5*(fabs((double)(x))+fabs((double)(y)))+fabs((double)(tol)))<= fabs((double)(tol)) ? 1 : 0)
#define AVE(x,y) (0.5*(x +y))

/*Structures */
/*************************************************************************/
typedef struct PINCH /* Eclipse record */
{
  int set;
  double pVal;
  char gap[TOKEN_LENGTH+1];
  double pinchPV;
  char topbot[TOKEN_LENGTH+1];
  char tMeth[TOKEN_LENGTH+1];
} PINCH;

typedef struct SPECGRID /* Eclipse record */
{
  int set;
  int nx;
  int ny;
  int nz;
  int numRes;
  char cylindrical;
} SPECGRID ;

typedef struct COORDSYS /* Eclipse record */
{
  int set;
  int numRes;
  int *k1;
  int *k2;
  char** comp;
  char** join;
} COORDSYS;

typedef struct MODFILES /* Holds queue of files for includes; read and write */
{
  int nInc;
  FILE *fp[N_INC_MAX+1];
  char fn[NAME_LENGTH+1][N_INC_MAX+1];
} MODFILES;

typedef struct GRIDUNIT /* Eclipse record */
{
  int set;
  char units[TOKEN_LENGTH+1];
  char map[TOKEN_LENGTH+1];
} GRIDUNIT;

typedef struct MAPAXES /* Optional Eclipse record */
{
  int set;
  double **axes; /* an array with 
		    axes[i][1] being the x point and [i][2] being the y point
		    1  is on the y axis, 2 is the origin, 3 is on the x-axis
		    NOTE: y is flipped between MAP and GRID systems!! */
  int isReversed; /* flag to indicate that the axis has been flipped */
} MAPAXES;

typedef struct MAPUNITS /* Optional Eclipse record */
{
  int set;
  char units[TOKEN_LENGTH+1];
} MAPUNITS;



typedef struct PROP /* any type property structure */
{
  int set;
  char name[TOKEN_LENGTH+1];
  char name2[TOKEN_LENGTH+1]; 
  /* Used to hold "other" name for covariances OR to hold an "alias" for outputing
     non-Eclipse records in an Eclipse-parsable form */
  int nInd; /* number of indices used for this property; usually 3 BUT
	       may be 4 for ZCORN or other applications */
  int **bnd; /* bnd[1..nInd][1..2] give lower & upper bounds in each direction
		Typically, [][1] = 1 and [][2]=(nx,ny,or nz), but this is not
		the case for COORD bnd = {{0,nx}, {0,ny}, {1,6}} 
		or ZCORN bnd = {{1,nx}, {1,ny}, {1, nz}, {1,8}} */
  char ptype; /* type is 'i', 'f', or 'd' */
  GENERIC_POINTER p; 
  /* "Universal" property array pointer for all types and nInd */   
  GENERIC_POINTER c; /* Holds count data for covariance structure */
  int extremes;
  double pMax;	
  double pMin; 
  int stats;
  double mean;
  double stdDev;
  int bounds;
  int clipped;
  double floor;
  double ceiling;
} PROP;

typedef struct ECL
{	
  char name[NAME_LENGTH+1];
  SPECGRID specGrid;
  MODFILES inFiles;
  MODFILES outFiles;
  GRIDUNIT gridUnit;
  COORDSYS coordSys;
  MAPUNITS mapUnits;
  MAPAXES mapAxes;
  PROP *zCornp;
  PROP *coordp;
  PROP *actnump;
  PINCH pinch;
  double ****zCorn;
  double ***coord;
  int ***actnum;
  int nProp;
  int nPropAlloc;
  PROP *prop;	
} ECL;

#endif
/*************************************************************************/
