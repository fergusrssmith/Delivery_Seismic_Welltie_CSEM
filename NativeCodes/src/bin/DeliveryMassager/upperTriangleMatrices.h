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
#ifndef __upperTriangleMatrices_h__
# define __upperTriangleMatrices_h__

#define USE_TWOBYTE_COVARIANCE_MATRIX 1

/* Although this is a macro, code won't work with values other than 2, if 
   USE_TWOBYTE_COVARIANCE_MATRIX == 1
*/
#define BYTESPERFLOAT 2

#define ROUNDUP_MULTIPLIER 1.0078125
/*
  ROUNDUP_MULTIPLIER = 129/128, closest factor to unity that will reliably roundup with
  2 byte floats

  Obtaining  multiplier in IEEE 754 floating point representation:
  1 for sign, 8 exponent, 23 bits for significand. Killing two bytes gives
  7 for significant, so 2^7=128.

  The hard way...

  float val,valorig,val2;
  char *val_ptr;
  int IMAX=10000000;
  long seed=-1;
  int base=128; // the biggest number that passes the test....

  for (i=1;i<=IMAX;i++)
	{
	  val=(float)i/(float)IMAX;
	  valorig=val;

	  val *= (base+1.0)/((float)base);
	  val_ptr = (char*)&val;	  
	  val_ptr[0]=(char)0;
	  val_ptr[1]=(char)0;

	  if (val < valorig ) {printf("FAILS TEST!\n"); exit(0);}
	}
  printf("PASSES TEST\n");

*/



float* getUpperTriMatrix(int n,int useDiagOnly);
void freeUpperTriMatrix(float* A, int n);
void copyUpperTriMatrix(float* Adest,float* Asrc,int n);
float get(float* A,int i,int j);
float getDiag(float* A,int i);
void set(float* A,float val,int i,int j); 
void accumulate(float* A,float val,int i,int j);
void divide(float* A,float val,int i,int j);

void set_quick_access(int i,int j);
void set_quick(float* A,float val);
void set_quick_roundup(float *A,float val);
float get_quick(float* A);

#endif
