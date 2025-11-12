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
#ifndef __util_h__
# define __util_h__

# ifndef MIN
#  define       MIN(a,b)                        (((a) < (b)) ? (a) : (b))
# endif

# ifndef MAX
#  define       MAX(a,b)                        (((a) > (b)) ? (a) : (b))
# endif

# ifndef SQR
#  define       SQR(x)                          ((x) * (x))
# endif

int is_big_endian();
float inverseG(float p);
void show_malloc_stats();
void f4tensor_fill(float val,float ****f,int Nx,int Ny,int Nz,int Nw);
void cleanKrigingSystemOfDuplicateNeighbours(int N,double** C,double* RHSkrig);
void getCholeskySafe(double **a, int n);
void getLDLTdecompositionSafe(double **a, int n);
float random1(long *seed);
float Gaussian_deviate(long *seed);

int CrossingsTest(double	pgon[][2], int numverts, double	point[2] );
double quick_select(double arr[], int n);

void print_dmatrix0(double **m,int nRow,int nCol,char* filename);
void print_dmatrix0_mma(double **m,int nRow,int nCol,char* filename);
void print_dvector0(double *m,int n,char* filename);
void print_dvector0_mma(double *m,int n,char* filename);

#endif
