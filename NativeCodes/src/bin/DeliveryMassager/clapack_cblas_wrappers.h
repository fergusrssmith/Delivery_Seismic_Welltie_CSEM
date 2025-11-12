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
/*
** massager.h
**
** Automatically generated header file.
** Source file:   deliveryMassager.c
** Last modified: Monday, March 29, 2004 at 15:26:43
*/

#ifndef __clapack_cblas_wrapper_h__
# define __clapack_cblas_wrapper_h__

/* clapack wrappers */
int     clapack_wrapper_dgetrf  (int , int , double* , int , int*);
int     clapack_wrapper_dgetrs  (char , int , int , double* , int , int* , double* , int);
int     clapack_wrapper_spotrf  (char , int , float* , int);
int     clapack_wrapper_dpotrf  (char , int , double* , int);
int     clapack_wrapper_dpotri  (char , int , double* , int);
void    cblas_wrapper_dgemv(char trans, int M, int N,
							double alpha, double *A, int lda,
							double *X, int incX, double beta,
							double *Y, int incY);
void cblas_wrapper_dgemm(char TransA,
                 char TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
#endif
