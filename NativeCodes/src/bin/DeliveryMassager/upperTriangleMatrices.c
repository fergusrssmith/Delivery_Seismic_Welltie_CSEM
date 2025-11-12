/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "util.h"
#include "platform_specifics.h"
#include "upperTriangleMatrices.h"

/*
	Utilities for Upper-Packed Triangular Storage Mode:  

	Two possible storage formats:

	  (A)
	  When a symmetric matrix is stored in upper-packed storage mode,
	  the upper triangular part of the symmetric matrix is stored,
	  including the diagonal, in a one-dimensional array. The upper
	  triangle is packed by columns. (This is equivalent to packing the
	  lower triangle by rows.) The matrix is packed sequentially column
	  by column in n(n+1)/2 elements of a one-dimensional array. To
	  calculate the location of each element aij of matrix A in an array
	  AP using the upper triangular packed technique, use the following
	  formula:
	  
	  A[i+(j(j-1)/2)] = aij    where j >= i , i,j>=1 (fortranesque)
	  A[i+(j(j+1)/2)] = aij    where j >= i , i,j>=0 (C)
	  
	  (B) 
	  Storage by diagonals for A[0.....n(n+1)/2-1]

	  A[i+(j-i)(2n+1-(j-i))/2] = A_ij   j>=i, i,j>=0 (C-style)

	These routines use format B, so diagonal-only matrices can be stored in n 
	consecutive members without wasted space.

	Usage assumes any multiple instances of getUpperTriMatrix() in concurrent usage have the
	same dimension n. If useDiagOnly=1 called, storage for diagonal only is alloced,
	and subsequent accesses are *not protected*, for speed reasons. Caveat emptor.

	The *_quick methods are for increased effeciency for routines having repeated
	access to the same i,j element. 

*/




static THREAD short needsDiagOnly=0;
static THREAD short twoNplus1 = 1;
static THREAD short mostSigByte = 0;
static THREAD int last_offset = 0;

#if !USE_TWOBYTE_COVARIANCE_MATRIX

/* routines storing entries as full floats */

float* getUpperTriMatrix(int n,int useDiagOnly)
{
  float *A;
  if (useDiagOnly) 
	{
	  needsDiagOnly=1;
	  A = (float*)calloc((size_t)n, sizeof(float));
	}
  else
	{
	  needsDiagOnly = 0;
	  A = (float*)calloc((size_t)((n*(n+1))/2), sizeof(float));
	}
  if (!A) fprintf(stderr,"allocation failure in getUpperTriMatrix()");

  twoNplus1 = 2*n+1;


  return A;
}
 
/* Free an upper triangular matrix */
void freeUpperTriMatrix(float* A, int n)
{
  free(A);
  A=NULL;
}

void copyUpperTriMatrix(float* Adest,float* Asrc,int n)
{
  if (needsDiagOnly) 
	{
	  memcpy(Adest,Asrc,n*sizeof(float));
	}
  else
	{
	  memcpy(Adest,Asrc,((n*(n+1))/2)*sizeof(float));
	}
}

/* get the i,j value. Ensure j>=i */
float get(float* A,int i,int j)
{ 
  // format A:
  // return A[i+(j*(j+1))/2];
  // return A[i+((j*(j+1))>>1)];
  // format B:
  j-=i; // j <- j-i
  return A[i+((j*(twoNplus1-j))>>1)];
}

/* set handy offset index for rapid repeated get_quick(float* A) access */
void set_quick_access(int i,int j)
{
  j-=i; // j <- j-i
  last_offset = i+((j*(twoNplus1-j))>>1);
}

/* get the i,j value based on set_quick_access(int i,int j) */
float get_quick(float* A)
{ 
  return A[last_offset];
}

/*
 no roundup here with full precision
*/
void set_quick_roundup(float *A, float val)
{
  A[last_offset] = val;
}

/* set the i,j value based on set_quick_access(int i,int j) */
void set_quick(float* A,float val)
{
  /*
   A[i+(j*(j+1))/2] = val; 
   A[i+((j*(j+1))>>1)] = val;
  */
  A[last_offset] = val;
}


float getDiag(float* A,int j)
{ 
  // format A:
  //  return A[(j*(j+3))/2];
  // return A[(j*(j+3))>>1];
  // format B:
  return A[j];
}
 
/* set the i,j value. Ensure j>=i */
void set(float* A,float val,int i,int j)
{
  //   A[i+(j*(j+1))/2] = val;
  //  A[i+((j*(j+1))>>1)] = val;
  j-=i; // j <- j-i
  A[i+((j*(twoNplus1-j))>>1)] = val;
}

/* set the i,j value. Ensure j>=i 
   Full precision, so no roundup here 
*/
void set_roundup(float* A,float val,int i,int j)
{
  //   A[i+(j*(j+1))/2] = val;
  //  A[i+((j*(j+1))>>1)] = val;
  j-=i; // j <- j-i
  A[i+((j*(twoNplus1-j))>>1)] = val;
}

/* Accumulate: A[i][j] += val. Ensure j>=i */
void accumulate(float* A,float val,int i,int j)
{ 
  //  int indx = i+(j*(j+1))/2;
  // int indx = i+((j*(j+1))>>1);

  // A[i+((j*(j+1))>>1)] += val;
  j-=i; // j <- j-i
  A[i+((j*(twoNplus1-j))>>1)] += val;
}

/* Divide: A[i][j] /= val. Ensure j>=i */
void divide(float* A,float val,int i,int j)
{ 
  //  int indx = i+(j*(j+1))/2;
  // int indx = i+((j*(j+1))>>1);

  //  A[i+((j*(j+1))>>1)] /= val;
  j-=i; // j <- j-i
  A[i+((j*(twoNplus1-j))>>1)] /= val;
}

#else  // USE_TWOBYTE_COVARIANCE_MATRIX

/**********************************************************************
 
   Memory efficient calls using only two bytes per float 

************************************************************************/

/* Allocate an upper triangular matrix A[i,j], with i,j = 0....n-1, with j>=i */

float* getUpperTriMatrix(int n,int useDiagOnly)
{
  float *A;

  if (useDiagOnly)
	{
	  needsDiagOnly = 1;
	  A = (float *)calloc((size_t)(BYTESPERFLOAT*n),1);
	}
  else
	{
	  needsDiagOnly = 0;
	  A = (float *)calloc((size_t)(BYTESPERFLOAT*((n*(n+1))/2)),1);
	}
  if (!A) fprintf(stderr,"allocation failure in getUpperTriMatrix()");

  twoNplus1 = 2*n+1;

  mostSigByte = is_big_endian()? 0 : 2;

  return A;
}
 
/* Free an upper triangular matrix */
void freeUpperTriMatrix(float* A, int n)
{
  free(A);
  A=NULL;
}

void copyUpperTriMatrix(float* Adest,float* Asrc,int n)
{
  if (needsDiagOnly) 
	{
	  memcpy(Adest,Asrc,n*BYTESPERFLOAT);
	}
  else
	{
	  memcpy(Adest,Asrc,((n*(n+1))/2)*BYTESPERFLOAT);
	}
}

/* get the i,j value. Ensure j>=i */
float get(float* A,int i,int j)
{ 
  //  return A[i+(j*(j+1))/2];
  float val=0;
  int indx;
  char *val_ptr;
  char *A_ptr;

  val_ptr = (char*)&val;
  val_ptr += mostSigByte;

  // format B: A[i+(j-i)(2n+1-(j-i))/2] 
  j-=i; // j <- j-i
  indx = BYTESPERFLOAT*(i+((j*(twoNplus1-j))>>1));
  // little endian code....use [2][3], BE use [0][1]
  A_ptr = (char*)&(((char*)A)[indx]);
  *val_ptr = *A_ptr;
  *(++val_ptr) = *(++A_ptr);
  
  return val;
}

/* set handy offset index for rapid repeated get_quick(float* A) access */
void set_quick_access(int i,int j)
{
  j-=i; // j <- j-i
  last_offset = BYTESPERFLOAT*(i+((j*(twoNplus1-j))>>1));
}

/* get the i,j value. Ensure j>=i */
float get_quick(float* A)
{ 
  float val=0;
  char *val_ptr;
  char *A_ptr;

  val_ptr = &((char*)&val)[mostSigByte];
  A_ptr = (char*)&(((char*)A)[last_offset]);
  *val_ptr = *A_ptr;
  *(++val_ptr) = *(++A_ptr);

  return val;
}

/* set the i,j value based on set_quick_access(int i,int j) */
void set_quick(float* A,float val)
{
  //   A[i+(j*(j+1))/2] = val;
  //  A[i+((j*(j+1))>>1)] = val;
  char *val_ptr;
  char *A_ptr = (char*)A;

  val_ptr = &((char*)&val)[mostSigByte];
  // little endian code....use [2][3], BE use [0][1]
  A_ptr = (char*)&(((char*)A)[last_offset]);
  *A_ptr = *val_ptr;
  *(++A_ptr) = *(++val_ptr);
}

/* set the i,j value based on set_quick_access(int i,int j). 
   Makes sure things are safely rounded up, so the
   recovered A[i][j] >= val on all subsequent recoveries */
void set_quick_roundup(float* A,float val) 
{
  set_quick(A, val*ROUNDUP_MULTIPLIER);
}


float getDiag(float* A,int j)
{ 
  //  return A[(j*(j+3))/2];
  int indx;
  float val=0;
  char *val_ptr;
  char *A_ptr = (char*)A;

  val_ptr = &((char*)&val)[mostSigByte];
  indx = BYTESPERFLOAT*j;
  // little endian code....use [2][3], BE use [0][1]
  *val_ptr = A_ptr[indx];
  *(++val_ptr) = A_ptr[indx+1];

  return val;
  //return get(A,j,j);
}
 
/* set the i,j value. Ensure j>=i */
void set(float* A,float val,int i,int j)
{
  //   A[i+(j*(j+1))/2] = val;
  //  A[i+((j*(j+1))>>1)] = val;
  int indx;
  char *val_ptr;
  char *A_ptr = (char*)A;

  val_ptr = (char*)&val;
  //  return A[i+(j*(j+1))/2];
  j-=i; // j <- j-i
  indx = BYTESPERFLOAT*(i+((j*(twoNplus1-j))>>1));
  // little endian code....use [2][3], BE use [0][1]
  A_ptr[indx] = val_ptr[mostSigByte];
  A_ptr[indx+1] = val_ptr[mostSigByte+1];
}

/* set the i,j value. Ensure j>=i 
   Rounds safely so returned value is always >= val 
*/
void set_roundup(float* A,float val,int i,int j)
{
  set(A, val*ROUNDUP_MULTIPLIER, i, j);
}

/* Accumulate: A[i][j] += val. Ensure j>=i */
void accumulate(float* A,float val,int i,int j)
{ 
  // format A:
  // int indx = i+(j*(j+1))/2;
  // int indx = i+((j*(j+1))>>1);
  // A[i+((j*(j+1))>>1)] += val;
  int indx;
  float oldval=0;
  char *val_ptr;
  char *A_ptr = (char*)A;
 
  val_ptr = (char*)&oldval;

  // format B:
  j-=i; // j <- j-i
  indx = BYTESPERFLOAT*(i+((j*(twoNplus1-j))>>1));
  // little endian code....use [2][3], BE use [0][1]
  val_ptr[mostSigByte] = A_ptr[indx];
  val_ptr[mostSigByte+1] = A_ptr[indx+1];

  // accumulate
  oldval+=val;

  // write
  A_ptr[indx] = val_ptr[mostSigByte];
  A_ptr[indx+1] = val_ptr[mostSigByte+1];
}

/* Divide: A[i][j] /= val. Ensure j>=i */
void divide(float* A,float val,int i,int j)
{ 
  //  int indx = i+(j*(j+1))/2;
  // int indx = i+((j*(j+1))>>1);

  // A[i+((j*(j+1))>>1)] /= val;

  int indx;
  float oldval=0;
  char *val_ptr;
  char *A_ptr = (char*)A;

  val_ptr = (char*)&oldval;

  j-=i; // j <- j-i
  indx = BYTESPERFLOAT*(i+((j*(twoNplus1-j))>>1));
  // little endian code....use [2][3], BE use [0][1]
  val_ptr[mostSigByte] = A_ptr[indx];
  val_ptr[mostSigByte+1] = A_ptr[indx+1];

  // divide
  oldval/=val;

  // write
  A_ptr[indx] = val_ptr[mostSigByte];
  A_ptr[indx+1] = val_ptr[mostSigByte+1];
}

#endif


