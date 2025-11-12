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
// Uses glib hashing fuctions
#include <glib.h>

#include "laspack/matrix.h"

// sparse library utilities
#include "sparse_util.h"



/* Stores lookup hash table in integer array. Free this when hashing no longer required */

struct LASPACK_HASH buildLaspackHashTable(Matrix *A)
{
  int i,j,k,N,size;
  struct LASPACK_HASH Lhash;
  N = M_GetRowDim(A);
  size=0;
  int key;
  double *result;

  Lhash.RowDim = N;
  for (i=1;i<=N; i++)
	{
	  // size of local system from sparse prior entry count on this row
	  size += M_GetLen(A,i); 
	}
  Lhash.keys = (int*)malloc((size_t) (size* sizeof(int)));
  Lhash.hash_table = g_hash_table_new(g_int_hash,g_int_equal);
  for (i=1,k=0;i<=N; i++)
	{
	  size = M_GetLen(A,i);
	  for (j=0;j<size;j++,k++)
		{
		  // key is unrolled, zero-offset index of equivalent rectangular array
		  Lhash.keys[k]= (i-1)*N + M_GetPos(A, i, j)-1;
		  // data is sparse element
		  g_hash_table_insert(Lhash.hash_table,(gpointer)(&Lhash.keys[k]),
                                             (gpointer)(&(A->El[i][j].Val)));
		}
	}


#if 0
  // test of hash map
  for (i=1;i<=N; i++)
	{
	 size = M_GetLen(A,i);
	  for (j=1;j<=N;j++)
		{ 
		  key = (i-1)*N + j-1;
		  result = g_hash_table_lookup(Lhash.hash_table,&key);
		  double ans = M_GetEl(A,i,j);
		  if (result!=NULL || ans!=0)
			{
			  if ((double)(*result) != ans)
				{
				  printf("error %d %d\n",i,j);
				}
			}
		}
	}
#endif

  return Lhash;
}

void freeLaspackHashTable(struct LASPACK_HASH *Lhash)
{
  free(Lhash->keys);
  g_hash_table_destroy(Lhash->hash_table); 
}

int getHashKey(struct LASPACK_HASH *Lhash,int row,int col)
{
  return (row-1)*Lhash->RowDim + col-1;
}

double* getEntry(struct LASPACK_HASH *Lhash,int row,int col)
{
  int key;
  
  key = getHashKey(Lhash,row,col);
  return g_hash_table_lookup(Lhash->hash_table,&key);
}
