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
#include <stdio.h>
#include "NR/nrutil.h"
#include "bayes/bayes.h"

#include "kdtree/nn.h"
#include "kdtree/optkd.h"

#include "std/std.h"
#include "std/error.h"

#include "massager_types.h"

#ifdef USE_LASPACK
#include "laspack/qmatrix.h"
#include "laspack/vector.h"
#endif

/*
  Various utility functions, chiefly for debugging
*/

void print_dmatrix0(double **m,int nRow,int nCol,char* filename)
{
  int i,j;
  FILE *out;
        
  if ((out = fopen(filename, "w")) == NULL) 
	{
	  printf("Warning: Can't open %s",filename);
	  return;
	}
  printf("Printing matrix to %s\n",filename);
  for (i=0;i<nRow;i++)
	{
	  for (j=0;j<nCol;j++)
		{
		  fprintf(out,"%15.15f\n",m[i][j]);
		}
	}
  fclose(out);
}

void print_dmatrix0_mma(double **m,int nRow,int nCol,char* filename)
{
  int i,j;
  FILE *out;
        
  if ((out = fopen(filename, "w")) == NULL) 
	{
	  printf("Warning: Can't open %s",filename);
	  return;
	}
  printf("Printing matrix to %s\n",filename);
  fprintf(out,"{");
  for (i=0;i<nRow;i++)
	{
	  fprintf(out,"{");
	  for (j=0;j<nCol;j++)
		{
		  fprintf(out,"%15.15f",m[i][j]);
		  if (j<nCol-1)fprintf(out,", ");
		}
	  fprintf(out,"}");
	  if (i<nRow-1) fprintf(out,",");
	}
  fprintf(out,"}");
  fclose(out);
}

void print_dvector0(double *m,int n,char* filename)
{
  int i;
  FILE *out;
  
  if ((out = fopen(filename, "w")) == NULL) 
	{
	  printf("Warning: Can't open %s",filename);
	  return;
	}
  printf("Printing vector to %s\n",filename);
  for (i=0;i<n;i++)
	{
	  fprintf(out,"%g\n",m[i]);
	}
  fclose(out);
}

#ifdef USE_LASPACK

void print_dmatrix0_laspack(Matrix *m,char* filename)
{
  int i,j,nRows,nCols;
  FILE *out;
        
  if ((out = fopen(filename, "w")) == NULL) 
	{
	  printf("Warning: Can't open %s",filename);
	  return;
	}
  printf("Printing matrix to %s\n",filename);
  nRows = M_GetRowDim(m);
  nCols = M_GetClmDim(m);
  
  for (i=1;i<=nRows;i++)
	{
	  for (j=1;j<=nCols;j++)
		{
		  fprintf(out,"%15.15f\n",M_GetEl(m,i,j));
		}
	}
  fclose(out);
}

#endif
