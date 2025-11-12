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
// This is a translation & adaption of toms 695 by Elizabeth Eskow and 
// Robert B. Schnabel
//
// Author: James Gunning
// Date:   Aug 1/2005

// Includes also updates in "A Revised Modified Cholesky Factorization Algorithm"
// Robert B. Schnabel, Elizabeth Eskow.
// SIAM Journal on Optimization, 9,No 4, Pages 1135 - 1148
//
//*********************************************************************
//
//       subroutine name: modchl
//
//       authors :  Elizabeth Eskow and Robert B. Schnabel
//
//       date    :  August, 1999 (revised from original code of 12/88)
//
//       purpose : perform a modified cholesky factorization
//                 of the form (Ptranspose)AP  + E = L(Ltranspose),
//       where L is stored in the lower triangle of the
//       original matrix A.
//       The factorization has 2 phases:
//        phase 1: Pivot on the maximum diagonal element.
//            Check that the normal cholesky update
//            would result in a positive diagonal
//            at the current iteration, and
//            if so, do the normal cholesky update,
//            otherwise switch to phase 2.
//        phase 2: Pivot on the minimum of the negatives
//            of the lower gerschgorin bound
//            estimates.
//            Compute the amount to add to the
//            pivot element and add this
//            to the pivot element.
//            Do the cholesky update.
//            Update the estimates of the
//            gerschgorin bounds.
//
//       input   : 
//                 n       - dimension of matrix A
//
//                 A       - n*n symmetric matrix (only lower triangular
//            portion of A, including the main diagonal, is used)
//
//                 g       - n*1 work array
//
//                 mcheps - machine precision
//
//                 mu     - tolerance used for determining when to switch
//                          phase 2 via lookahead, or due to the  
//                          relative magnitude of the current 
//                          min. diag. element.
//
//                tau     - tolerance used for determining the maximum
//                          condition number of the final 2X2 submatrix.
//          
//                taubar  - tolerance used for determining the amount
//                          added to the smallest eigenvalue of the 
//                          final 2X2 submatrix.
//
//
//
//       output  : L     - stored in the matrix A (in lower triangular
//                           portion of A, including the main diagonal)
//
//                 P     - a record of how the rows and columns
//                         of the matrix were permuted while
//                         performing the decomposition
//
//                 E     - n*1 array, the ith element is the
//                         amount added to the diagonal of A
//                         at the ith iteration
//
//
//***********************************************************************
// Modifications added by J.G. as per "A revised Modified Cholesky Factorization Algorithm", 1999
// Makes a big difference for approximate covariance matrices
// 
// ***********************************************************************

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "NR/nrutil.h"



#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

// switch on for testing only
#define PERFORM_TESTS 0

// Machine epsilon
#define MCHEPS 1.11e-16
#define TAUBAR 2.3e-11   // MCHEPS^(2/3)
#define TAU 4.8e-6       // MCHEPS^(1/3)
#define MU 0.1 

struct MCLD_GLOBALS
{
  double delta;
  double ming;
  double taugam;
  double taubargam;
  double mugam;
  double* diagScale; 
};

/* private routines not needed outside here */

void modchl(struct MCLD_GLOBALS *mcld_globals,int n,double **A,double* g,double mcheps,double tau,double taubar,double mu,int* P,double* E);
void init(struct MCLD_GLOBALS *mcld_globals, int n,double** A,int* P,double* g,double* E, double taubar,double mu);
void gersch(int , double** , int , double*);
void fin2by2(struct MCLD_GLOBALS*, int , double** , double* , int , double);


/** 
	<p>
	From a symmetric, positive definite matrix <tt>A</tt>,
	constructs a 'safely' positive definite Cholesky Left decomposition
	using the methods described in
	"Algorithm 695 , collected algorithms from ACM.
	Transactions on mathematical software,
	Vol. 17, no. 3, September, 1991, pp. 306-312."

	Includes also updates in "A Revised Modified Cholesky Factorization Algorithm"
    Robert B. Schnabel, Elizabeth Eskow.
    SIAM Journal on Optimization, 9,No 4, Pages 1135 - 1148)

	Guaranteed to rescue
	you if the matrix is ill-conditioned or not quite positive definite.
	</p>
	<p>
	The decomposition is of form P<sup>T</sup>.A.P = L.L<sup>T</sup>
	</p>
	<p>
	Example Usage:<br>
	<pre>

	int n=10;
	// set up some matrix A[0...n-1][0...n-1] ....,then
	int *p = malloc((size_t)n*sizeof(int));
	modifiedCholeskyLeftDecomposition(n,A,p,preScale);

	// L is now lower tri of A
	// p is now inplace in p[]

	</pre>
    </p>	
	
	@param  n   size of square matrix
	@param  A   Square matrix: Positive definite or close thereto.
	@param  P   0...n-1 integer arrays of permutations (call allocated: will be filled).
	@param  preScale. Set to 1 to prescale everything by the diagonals before decomposing)
*/
void modified_Cholesky_left_decomposition(int n, double** A, int* P, int preScale)
{
  struct MCLD_GLOBALS mcld_globals = {0,0,0,0,0,NULL};

  // initialise
  int i,j;
  double *g,*E; 
  double *diagScale,*divScale;
  double s;

#if PERFORM_TESTS
  int k;
  double **B,**Acopy;
  Acopy = dmatrix(0,n-1,0,n-1);
  memcpy(&Acopy[0][0],&A[0][0],n*n*sizeof(double));
#endif

  g = calloc(n,(size_t)sizeof(double));
  E = calloc(n,(size_t)sizeof(double));

  if (preScale)
	{
	  diagScale = malloc((size_t)n*sizeof(double));
	  divScale = malloc((size_t)n*sizeof(double));
	  mcld_globals.diagScale = diagScale;
	  // scaling seems to work best on numbers 
	  // within representable range. Smaller numbers are probably noise, and will
	  // wreck the scaling
	  s=0;
	  for (i=0;i<n;i++) s=MAX(s,fabs(A[i][i]));
	  s *= MCHEPS;
	  for (i=0;i<n;i++) 
		{
		  if (A[i][i]!=0.0 && fabs(A[i][i])>s)
			{
			  diagScale[i] = sqrt(fabs(A[i][i]));
			  divScale[i] = diagScale[i]; 
			}
		  else
			{
			  diagScale[i] = 1.0;
			  divScale[i] = 1.0;
			}
		} 
	  for (i=0;i<n;i++) 
		{
		  for (j=0;j<n;j++) 
			{
			  A[i][j] /= (divScale[i]*divScale[j]); 
			}
		}
	  free(divScale);
	}
  
  modchl(&mcld_globals,n, A, g, MCHEPS, TAU, TAUBAR, MU,P, E);

  // A final fix of Eskow's routine: erase the upper right triangle
  for (i=0;i<n;i++) {
	for (j=i+1;j<n;j++) A[i][j] = 0; 
  }
  
  // and a final unscaling if required
  if (preScale)
	{
	  // L=S.Lhat
	  for (i=0;i<n;i++) 
		{
		  s = diagScale[P[i]];
		  for (j=0;j<=i;j++) 
			{
			  A[i][j] *= s;
			}
		}
	  free(diagScale);
	}
  free(E);
  free(g);



#if PERFORM_TESTS
  double err=0,*sc,maxAbs=0;
  // Check multiplication
  // B = A.A^T
  B = dmatrix(0,n-1,0,n-1);
  sc = dvector(0,n-1);
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=A[i][k]*A[j][k];
	}
	sc[i] = (Acopy[P[i]][P[i]]==0.0)? 0.0 : sqrt(fabs(Acopy[P[i]][P[i]]));
  }
  // permute
  double max=0.0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  err = fabs(B[i][j] - Acopy[P[i]][P[j]]);
	  maxAbs=MAX(maxAbs, err);
	  if (sc[i]>0)
		{
		  if (sc[j]>0) err /= (sc[i] * sc[j]);
		  else err /= (sc[i] * sc[i]);
		}
	  else
		{
		  if (sc[j]>0) err /= (sc[j] * sc[j]);
		  else err= 0.0;
		}
	  max=MAX(max, err);
	}
  }

  printf("\nMax scaled error in modified_Cholesky_left_decomposition = %g, max abs error=%g\n",max,maxAbs);  
  free_dvector(sc,0,n-1);
  free_dmatrix(B,0,n-1,0,n-1); 
  free_dmatrix(Acopy,0,n-1,0,n-1);
#endif

}

//  **********************************************************************
//  
//  subroutine name : modified_Cholesky_solve
//  
//  purpose :  solves (L.L^T)P(x)=P(b), where L is stored in
//  the lower triangle of A, P is the record of the
//  permutations  performed in forming L,b is the rhs vecto
//  and x is the result.
//  L is the result of the modified cholesky factorization
//  which computes Ptranspose(A+E)P=LLtranspose.
//  
//  n     - dimension of matrix
//  
//  A     - n*n array (only lower triangle portion,
//  including the diagonal, is used).
//  
//  P     - n*1 int vector that contains a record
//  of permutations performed in forming L.
//  i.e. each Pi was initialized to i, and if
//  row and columns i & j were interchanged when
//  computing L, the values of Pi and Pj were swapped.
//  
//  b     - n*1 double vector which is the
//  right hand side of the system to be solved.
//  
//  y,z   - n*1 double  work vectors.
//  
//  output:  x     - n*1 double result vector
//  
//  
//  **********************************************************************
//
// Solution checked to be OK; 19 Aug 2003 James Gunning

/*
   After the Cholesky Decomposition is instantiated, we can solve
   the linear problem A.x=b using this call. Solution returned in b.

   @param  n: dimension of system: n x n
   @param  A: Cholesky factor from modifiedCholeskyLeftDecomposition()
   @param  P: permutation from modifiedCholeskyLeftDecomposition()
   @param  b: rhs vector b[0...n-1]. Solution returned in b
*/
void modified_Cholesky_solve(int n, double **A,int *P,double* b)
{
  int i,j;
  double sum;
  double *y,*z;
  
  y = malloc((size_t)n*sizeof(double));
  z = malloc((size_t)n*sizeof(double));
  
  //  
  //  Solve Ly=Pb
  //  
  y[0] = b[P[0]]/A[0][0] ;
  for (i=0 ; i<n ; i++) 
	{
	  sum=0.0;
	  for (j=0 ; j<=i-1 ; j++) 
		{
		  sum += A[i][j] * y[j];
		}
	  y[i] = (A[i][i]!=0) ? (b[P[i]] - sum) / A[i][i] : 0.0;
	}
  //  
  //  Solve Ltranspose z = y
  //  
  z[n-1] = (A[n-1][n-1]!=0) ? y[n-1]/A[n-1][n-1] : 0.0;
  for (i=n-2; i>=0 ; i--) 
	{
	  sum=0.0;
	  for (j=i+1 ; j<n ; j++) 
		{
		  sum += A[j][i]  * z[j];
		}
	  z[i] = (A[i][i]!=0) ? (y[i] - sum) / A[i][i] : 0.0;
	}
  //  
  //  x = Ptranspose z
  //  Solve in-place
  //
  for (i=0 ; i<n ; i++) 
	{
	  b[P[i]] = z[i];
	}
  free(y);
  free(z);
}

/*
   Returns the inverse of the left Cholesky matrix L *in place*. L must
   be a square matrix lower triangular matrix L[0....n-1][0...n-1]
   
*/

void modified_Cholesky_left_inverse(int n, double** L)
{
  int i,j,k;
  double sum,tmp;
  
  // Main loop.
  for (i=0;i<n;i++) 
	{
	  tmp = L[i][i];
	  L[i][i] = (tmp!=0)? 1.0/tmp : 0.0;
	  for (j=i+1; j<n; j++) 
		{
		  sum=0.0;
		  for (k=i; k<j;k++) sum -= L[j][k]*L[k][i];
		  tmp = L[j][j];
		  L[j][i] = (tmp!=0) ? sum/tmp: 0.0;
		}
	}
}

void modified_Cholesky_restore(int n, double **A,int *p)
{
  int i,j,k;
  double** B;

  // B = A.A^T
  B = dmatrix(0,n-1,0,n-1);
  for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=A[i][k]*A[j][k];
	  // symmetry
	  B[j][i] = B[i][j];
	}
  }

  // permute
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  A[p[i]][p[j]] =  B[i][j] ; 
	}
  }
  free_dmatrix(B,0,n-1,0,n-1);
}







/**
   Translation of Eskow's modified Cholesky routine; zero offset.
   Private access: do not call from outside
*/
void modchl(struct MCLD_GLOBALS *mcld_globals,int n,double **A,double* g,
			double mcheps,double tau,double taubar,double mu, 
			int* P,double* E)
{
//*********************************************************************
//
//       subroutine name: modchl
//
//       authors :  Elizabeth Eskow and Robert B. Schnabel
//
//       date    :  August, 1999 (revised from original code of 12/88)
//
//       purpose : perform a modified cholesky factorization
//                 of the form (Ptranspose)AP  + E = L(Ltranspose),
//       where L is stored in the lower triangle of the
//       original matrix A.
//       The factorization has 2 phases:
//        phase 1: Pivot on the maximum diagonal element.
//            Check that the normal cholesky update
//            would result in a positive diagonal
//            at the current iteration, and
//            if so, do the normal cholesky update,
//            otherwise switch to phase 2.
//        phase 2: Pivot on the minimum of the negatives
//            of the lower gerschgorin bound
//            estimates.
//            Compute the amount to add to the
//            pivot element and add this
//            to the pivot element.
//            Do the cholesky update.
//            Update the estimates of the
//            gerschgorin bounds.
//
//       input   : ndim    - largest dimension of matrix that
//                           will be used
//
//                 n       - dimension of matrix A
//
//                 A       - n*n symmetric matrix (only lower triangular
//            portion of A, including the main diagonal, is used)
//
//                 g       - n*1 work array
//
//                 mcheps - machine precision
//
//                 mu     - tolerance used for determining when to switch
//                          phase 2 via lookahead, or due to the  
//                          relative magnitude of the current 
//                          min. diag. element.
//
//                tau     - tolerance used for determining the maximum
//                          condition number of the final 2X2 submatrix.
//          
//                taubar  - tolerance used for determining the amount
//                          added to the smallest eigenvalue of the 
//                          final 2X2 submatrix.
//
//
//
//       output  : L     - stored in the matrix A (in lower triangular
//                           portion of A, including the main diagonal)
//
//                 P     - a record of how the rows and columns
//                         of the matrix were permuted while
//                         performing the decomposition
//
//                 E     - n*1 array, the ith element is the
//                         amount added to the diagonal of A
//                         at the ith iteration
//
//
//***********************************************************************

  char phase1=1;
  int j,iming,i,imaxd,itemp,jp1,k;
  double normj,maxd,mind;
  double delta1,temp,jdmin,tdmin,tempjj;
  
  init(mcld_globals, n, A, P, g, E,taubar,mu);
  //  
  //  check for n=1
  //  
  if (n==1) 
	{
	  mcld_globals->delta = (taubar * fabs(A[0][0])) - A[0][0];
	  if (mcld_globals->delta > 0) E[0]=mcld_globals->delta;
	  if (A[0][0]  == 0) E[0]=taubar;
	  A[0][0] = sqrt(MAX(0,A[0][0] +E[0])); 
	}

  for (j = 0;j < n;j++)
	{
	  //  
	  //  PHASE 1
	  //  
	  if (phase1 ) 
		{
		  //  
		  //  Find index of maximum diagonal element A(i,i) where i>=j
		  //  
		  mind = maxd = A[j][j] ;
		  imaxd = j;
		  for (i = j+1;i<n;i++) 
			{
			  if (maxd < A[i][i] ) 
				{
				  maxd = A[i][i] ;
				  imaxd = i;
				}
			  if (mind > A[i][i]) mind = A[i][i] ;
			}
		  // if (maxd <mcld_globals->taubargam || (mind < -MU*maxd))
		  if (mind < -mu*maxd)
			{
			 phase1 = 0;
			}
		  //  Pivot to the top the row and column with the max diag
		  // 
		  // printf("j=%d, pivot=%d\n",j,imaxd);
		  if (phase1 ) 
			{			  
			  if (imaxd != j) 
				{
				  //  Swap row j with row of max diag
				  for (i=0; i<= j-1;i++) 
					{
					  temp = A[j][i] ;
					  A[j][i] = A[imaxd][i]; 
					  A[imaxd][i] = temp; 
					}
				  //  Swap colj and row maxdiag between j and maxdiag				
				  for (i = j+1;i<=imaxd-1;i++)
					{
					  temp = A[i][j] ;
					  A[i][j] = A[imaxd][i];                 
					  A[imaxd][i] = temp; 
					}
				  //  Swap column j with column of max diag
				  for (i = imaxd+1 ; i < n ; i ++) 
					{
					  temp = A[i][j] ;
					  A[i][j] = A[i][imaxd]; 
					  A[i][imaxd] = temp; 
					}
				  //  Swap diag elements
				  temp = A[j][j] ;
				  A[j][j] = A[imaxd][imaxd]; 
				  A[imaxd][imaxd] = temp; 
				  //  Swap elements of the permutation vector
				  itemp = P[j];
				  P[j]=P[imaxd];
				  P[imaxd]=itemp;
				}
			  //  Check to see whether the normal cholesky update for this;
			  //  iteration would result in a positive diagonal,
			  //  and if not then switch to phase 2.
		  
			  jp1 = j+1;
			  tempjj=A[j][j] ;	

			  if (tempjj > mcld_globals->taubargam) 
				{  		
				  if (j<(n-1))
					{	
					  jdmin=A[jp1][jp1] ;
					  for (i = jp1 ; i < n ; i ++) 
						{
						  temp = A[i][j]  * A[i][j]  / tempjj;
						  tdmin = A[i][i]  - temp;
						  jdmin = MIN(jdmin, tdmin);
						}
					  if (jdmin < -mcld_globals->mugam) phase1 = 0;
					}
				}
			  else
				{
				 phase1 = 0;
				}
			}

		  if (phase1) 
			{
			  //  Do the normal cholesky update if still in phase 1
			  A[j][j] = sqrt(MAX(0.0,A[j][j])); 
			  if (j<(n-1))
				{
				  tempjj = A[j][j] ;
				  for (i = jp1 ; i < n ; i ++) 
					{
					  A[i][j] /= tempjj; 
					}
				  for (i=jp1 ; i<n ; i++) 
					{
					  temp=A[i][j] ;
					  for (k = jp1 ; k <= i ; k++) 
						{
						  A[i][k] -= (temp * A[k][j]); 
						}
					}
				}
			}
		  else 
			{
			  //  Calculate the negatives of the lower gerschgorin bounds
			  gersch(n,A,j,g);
			}
		}

	  //  PHASE 2
	  // printf("j=%d, phase1=%d\n",j,phase1);	// FIXME
	  if (!phase1 && (j == n-1)) 
		{
		  // See algorithm in paper
		  delta1= -A[n-1][n-1] + MAX(-tau*A[n-1][n-1]/(1.0-tau), 
									  mcld_globals->taubargam);
		  mcld_globals->delta = MAX(mcld_globals->delta, delta1);
		  E[n-1] = mcld_globals->delta;
		  A[n-1][n-1] = sqrt(MAX(0,A[n-1][n-1] + E[n-1])); 
		  return;
		}
	  if (!phase1)
		{
		  if (j < (n-2))
			{
			  //  Find the minimum negative gershgorin bound
			  iming=j;
			  mcld_globals->ming = g[j];
			  for (i = j+1 ; i <n ; i ++) 
				{
				  if (mcld_globals->ming > g[i]) 
					{
					  mcld_globals->ming = g[i];
					  iming = i;
					}
				}
			  //  Pivot to the top the row and column with the
			  //  minimum negative gerschgorin bound
			  if (iming != j) 
				{
				  //  Swap row j with row of min gersch bound
				  for (i = 0 ; i<= j-1;i++) 
					{
					  temp = A[j][i] ;
					  A[j][i] = A[iming][i]  ;
					  A[iming][i] = temp ;
					}
				  //  Swap colj with row iming from j to iming
				  for (i = j+1 ; i <=iming-1 ; i ++) 
					{
					  temp = A[i][j] ;
					  A[i][j] = A[iming][i]  ;
					  A[iming][i] = temp ;
					}
				  //  Swap column j with column of min gersch bound
				  for (i = iming+1 ; i < n ; i ++) 
					{
					  temp = A[i][j] ;
					  A[i][j] = A[i][iming]  ;
					  A[i][iming] = temp ;
					}
				  //  Swap diagonal elements
				  temp = A[j][j] ;
				  A[j][j] = A[iming][iming]  ;
				  A[iming][iming] = temp ;
				  //  Swap elements of the permutation vector
				  itemp = P[j];
				  P[j]=P[iming];
				  P[iming]=itemp;
				  //  Swap elements of the negative gerschgorin bounds vecto
				  temp = g[j];
				  g[j] = g[iming];
				  g[iming]=temp;
				}
			  //  Calculate delta and add to the diagonal.
			  //  delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
			  //  where normj=sum of |A(i,j)|,for i=1,n,
			  //  delta_previous is the delta computed at the previous iter
			  //  and taugam is tau1*gamma.
			  normj = 0.0;
			  for (i = j+1 ; i < n ; i ++) 
				{
				  normj += fabs(A[i][j]);
				}
			  temp = MAX(normj,mcld_globals->taubargam);
			  delta1 = temp - A[j][j] ;
			  temp = 0.0;
			  delta1 = MAX(temp, delta1);
			  mcld_globals->delta = MAX(delta1,mcld_globals->delta);
			  E[j]=mcld_globals->delta;
			  A[j][j] += E[j] ;
			  //  Update the gerschgorin bound estimates
			  //  (note: g(i) is the negative of the
			  //  Gerschgorin lower bound.)
			  //  
			  if (A[j][j]  != normj) 
				{
				  temp = (A[j][j] != 0) ? (normj/A[j][j]) - 1.0 : -1.0; // assume norm=0 if A[j][j]==0
				  for (i = j+1 ; i < n ; i ++) 
					{
					  g[i] += fabs(A[i][j])* temp;
					}
				}
			  //  Do the cholesky update
				  
			  A[j][j] = sqrt(MAX(0,A[j][j])) ;
			  tempjj = A[j][j] ;
			  if (tempjj>0)
				{
				  for (i = j+1 ; i < n ; i ++) 
					{
					  A[i][j] /= tempjj ;
					}
				  for (i = j+1 ; i < n ; i ++) 
					{
					  temp = A[i][j] ;
					  for (k = j+1 ; k <= i ; k ++) 
						{
						  A[i][k] -= (temp * A[k][j]) ;
						}
					}
				}
			  else
				{
				  for (i = j+1 ; i < n ; i ++) A[i][j] = 0.0;
				} 
			}
		  else  
			{
			  fin2by2(mcld_globals,n, A, E, j, tau);
			  return;
			}
		}
	}
}

//  ***********************************************************************
//  subroutine name : init
//  
//  purpose : set up for start of cholesky factorization
//  
//  input : n, A, tau1
//  
//  sets: phase1    - short value set to 1 if in phase one,
//  otherwise 0.
//  delta     - amount to add to Ajj at iteration j
//  P,g,E - described above in modchl
//  ming      - the minimum negative gerschgorin bound
//  gamma     - the maximum diagonal element of A
//  taugam  - tau1 * gamma
//  
//  ***********************************************************************
void init(struct MCLD_GLOBALS *mcld_globals,
		  int n,double** A,int* P,double* g,double* E,
		  double taubar,double mu)
{
  int i;
  double gamma;

  mcld_globals->delta = 0.0;
  mcld_globals->ming = 0.0;
  for (i=0 ; i<n ; i++) 
	{
	  P[i]=i;
	  g[i]=0.0;
	  E[i]=0.0;
	}
  //  
  //  Find the maximum magnitude of the diagonal elements.
  //  
  gamma = 0.0;
  for (i=0 ; i<n ; i++)
	{
	  gamma = MAX(gamma,fabs(A[i][i]));
	}
  mcld_globals->taubargam = taubar * gamma;
  mcld_globals->mugam = mu * gamma;
}




//  ***********************************************************************
//  
//  subroutine name : gersch
//  
//  purpose : Calculate the negative of the gerschgorin bounds
//  called once at the start of phase II.
//  
//  input   : n, A, j
//  
//  output  : g - an n vector containing the negatives of the
//  Gerschgorin bounds.
//  
//  ***********************************************************************
void gersch(int n, double** A, int j, double* g)
{
  int i, k;
  double offrow;

  for (i = j ; i < n ; i ++) 
	{
	  offrow = 0.0;
	  for (k = j ; k <= i-1 ; k ++) 
		{
		  offrow += fabs(A[i][k] );
		}
	  for (k = i+1 ; k < n ; k ++) 
		{
		  offrow += fabs(A[k][i] );
		}
	  g[i] = offrow - A[i][i] ;
	}
}

//  ***********************************************************************
//  
//  subroutine name : fin2by2
//  
//  purpose : Handles final 2X2 submatrix in Phase II.
//  Finds eigenvalues of final 2 by 2 submatrix,
//  calculates the amount to add to the diagonal,
//  adds to the final 2 diagonal elements,
//  and does the final update.
//  
//  input : n, A, E, j, tau,
//  delta - amount added to the diagonal in the
//  previous iteration
//  
//  output : A - matrix with complete L factor in the lower triangle,
//  E - n*1 vector containing the amount added to the diagonal
//  at each iteration,
//  delta - amount added to diagonal elements n-1 and n.
//  
//  ***********************************************************************
void fin2by2(struct MCLD_GLOBALS*mcld_globals, int n, double** A, double* E, int j, double tau)
{
  double t1, t2, t3,lmbd1,lmbd2,lmbdhi,lmbdlo;
  double delta1;
  //  
  //  Find eigenvalues of final 2 by 2 submatrix
  //  
  t1 = A[n-2][n-2]  + A[n-1][n-1] ;
  t2 = A[n-2][n-2]  - A[n-1][n-1] ;
  t3 = sqrt(MAX(0,t2*t2 + 4.0*A[n-1][n-2] *A[n-1][n-2]));
  lmbd1 = (t1 - t3)/2.;
  lmbd2 = (t1 + t3)/2.;
  lmbdhi = MAX(lmbd1,lmbd2);
  lmbdlo = MIN(lmbd1,lmbd2);
  //  
  //  Find delta such that:
  //  1.  the l2 condition number of the final
  //  2X2 submatrix + delta*I <= tau
  //  2. delta >= previous delta,
  //  3. lmbdlo + delta >= tau * gamma,
  //  where lmbdlo is the smallest eigenvalue of the final
  //  2X2 submatrix
  //

  // See algorithm in paper
  delta1= MAX(tau*(lmbdhi-lmbdlo)/(1.0-tau), mcld_globals->taubargam);
  delta1 -= lmbdlo;
  delta1 = MAX(delta1, 0.0);
  mcld_globals->delta = MAX(mcld_globals->delta, delta1);

  if (mcld_globals->delta > 0.0) 
	{
	  A[n-2][n-2] += mcld_globals->delta ;
	  A[n-1][n-1] += mcld_globals->delta ;
	  E[n-2] = mcld_globals->delta;
	  E[n-1] = mcld_globals->delta;
	}
  //  
  //  Final update
  //  
  A[n-2][n-2] = sqrt(MAX(0,A[n-2][n-2])) ;
  A[n-1][n-2] = (A[n-2][n-2] !=0 ) ? A[n-1][n-2] /A[n-2][n-2] : 0.0 ;
  A[n-1][n-1] = A[n-1][n-1]  - (A[n-1][n-2] * A[n-1][n-2]);
  A[n-1][n-1] = sqrt(MAX(0,A[n-1][n-1]));
  return;
}


/* A useful debugging func */
void debug_print_matrix(void **m,int n,char type)
{
  int i,j;

  for (i=0;i<n;i++)
	{
	  for (j=0;j<n;j++)
		{
		  printf("%14.8g ",(double)((type=='f') ? ((float**)m)[i][j] : ((double**)m)[i][j]));
		}
	  printf("\n");
	}
 
}

#undef PERFORM_TESTS 
