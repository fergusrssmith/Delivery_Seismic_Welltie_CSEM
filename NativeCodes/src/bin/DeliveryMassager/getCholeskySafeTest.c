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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "NR/nrutil.h"

#include "util.h"
#include "modified_left_Cholesky.h"


/** 
	Test harness for getCholeskySafe()
	Compare with: java au.csiro.JamesGunning.lib.Numerical.ModifiedCholeskyLeftDecomposition [-s]

*/

int main(int argc,char **argv)
{
  void debug_print_matrix(void **m,int n,char type);
  int i,j,k,n,I=0,J=0;
  int I_rel,J_rel;
  double err,max_rel=0,Cr_rel=0,Cc_rel=0;
  double lambda,tmp;
  double **A,**R,**C,**B,*sc,Cr=0,Cc=0;
  int preScale= 0;
  char* filename=NULL;

  if (argc>1)
	{
	  if (!strcmp(argv[1],"-s")) 
		{
		  preScale=1;
		  printf("Using scaling\n");
		}
	  else if (!strcmp(argv[1],"-h")) 
		{
		  printf("Usage:\ngetCholeskySafeTest [-s|-h] [filename].\n");
		  printf("  - no filename supplied will invoke default test matrix\n");
		  printf("  - matrix may be native doubles (nxn), with .bin suffix, otherwise ascii assumed.\n");
		  exit(0);
		}
	  else filename = argv[1];
	}
  if (argc>2) 
	{
	  filename = argv[2];
	}


  if (!filename)
	{
	  n = 9;

	  // make a funky just-positive definite matrix.
	  double scale_eps = 1.0e-12;
	  double corr_range=100.0;
	  A = dmatrix(0,n-1,0,n-1);
	  double* scale = calloc(n,(size_t)sizeof(double));

	  for (i = 0; i < n; i++) 
		{
		  scale[i] = (i<n/3) ? 1.0-(double)i/(10.0*n) : ((i>(0.6666666*n))? 1.0-(double)i/(10.0*n) : scale_eps);
		}
  
	  for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
		  A[i][j] = scale[i]*scale[j]*exp(-fabs((double)i-j)/corr_range); 
		  // A[i][j] = exp(-(i-j)*(i-j)/10000.0); 
		}
	  }
	  
	  printf("Original A matrix:\n");
	  debug_print_matrix((void**)A,n,'d');
	}
  else
	{
	  FILE *in;
	  fpos_t pos;
	  double f,*ptr;
	  int cnt=0,isBin=0;

	  printf("Reading matrix from %s\n",filename);
	  // detect ".bin" suffix
	  for(i=strlen(filename)-1;i>=0;i--) 
		{
		  if (filename[i]=='.') 
			{
			  if (!strncmp(&(filename[i+1]),"bin",3)) 
				{
				  isBin=1;
				  printf("Reading as binary\n");
				}
			  break;
			}
		}
	  if (!isBin)
		{
		  // ascii
		  if (!(in=fopen(filename,"r"))) {printf("Can't open %s\n",filename); exit(0);}
		  printf("Reading as ascii\n");
		  fgetpos(in,&pos);
		  while (fscanf(in,"%lg\n",&f) == 1) cnt++;
		  n = (int)sqrt(cnt+0.1);
		  fsetpos(in,&pos);
		  A = dmatrix(0,n-1,0,n-1);
		  ptr=&(A[0][0]);
		  for (i=0;i<n*n;i++) fscanf(in,"%lg\n",&ptr[i]);
		  fclose(in);
		} 
	  else 
		{
		  // binary
		  if (!(in=fopen(filename,"rb"))) {printf("Can't open %s\n",filename); exit(0);}
		  fgetpos(in,&pos);
		  while (fread(&f,sizeof(double),1,in) == 1) cnt++;
		  n = (int)sqrt(cnt+0.1);
		  fsetpos(in,&pos);
		  A = dmatrix(0,n-1,0,n-1);
		  ptr=&(A[0][0]);
		  for (i=0;i<n*n;i++) fread(&ptr[i],sizeof(double),1,in);
		  fclose(in);
		}
	  if (n<10)
		{
		  printf("Original A matrix:\n");
		  debug_print_matrix((void**)A,n,'d');
		}
	}
  // precleaning?
  for (i=0;i<n;i++) {
	if (A[i][i]<=0)
	  {
		for (j=0;j<n;j++) 
		  {
			A[i][j] = A[j][i] = 0.0;
		  }
	  }
  }



  // take 2 copies
  C = dmatrix(0,n-1,0,n-1);
  double** D = dmatrix(0,n-1,0,n-1);
  memcpy(&C[0][0],&A[0][0],n*n*sizeof(double));
  memcpy(&D[0][0],&A[0][0],n*n*sizeof(double));

  // perform Left Cholesky
  int *p = ivector(0,n-1);
  for (i=0;i<n;i++) p[i]=i;
  //  modified_Cholesky_left_decomposition(n,A,p,preScale);
  getCholeskySafe(A,n);

  if (filename && n<10)
	{
	  printf("L=\n"); 
	  debug_print_matrix((void**)A,n,'d');
	  printf("Permutations:\n");
	  for (i=0;i<n;i++) {
		printf("%3d ",p[i]);
	  }
	  printf("\n");
	}
  
  // Check multiplication
  // B = A.A^T
  B = dmatrix(0,n-1,0,n-1);
  sc = dvector(0,n-1);
  double max=0;
  for (i=0;i<n;i++) {
	max=0;
	for (j=0;j<n;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=A[i][k]*A[j][k];
	  max=DMAX(max, fabs(B[i][j] - C[p[i]][p[j]]));
	}
	sc[i] = sqrt(DMAX(0,C[p[i]][p[i]]));
  }
  // permute
  max=max_rel=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  tmp = B[i][j] ;
	  A[p[i]][p[j]] =  tmp; 
	  err = fabs(tmp - C[p[i]][p[j]]);
	  if (err > max)
		{
		  max=err;
		  Cr = C[p[i]][p[i]];
		  Cc = C[p[j]][p[j]];
		  I=i;
		  J=j;
		}
	  if (sc[i]>0)
		{
		  err /= (sc[i] * ((sc[j]>0)? sc[j] : sc[i]));
		}
	  else err /= ((sc[j]>0)? (sc[j]*sc[j]) : 1e50);
	  if (err > max_rel)
		{
		  max_rel=err;
		  Cr_rel = C[p[i]][p[i]];
		  Cc_rel = C[p[j]][p[j]];
		  I_rel=i;
		  J_rel=j;
		}
	}
  }
  printf("Max abs error after multiplying out = %g, C_ii=%g  C_jj=%g i=%d j=%d\n",max,Cr,Cc,I,J);  
  printf("Max rel error after multiplying out = %g, C_ii=%g  C_jj=%g i=%d j=%d\n",max_rel,Cr_rel,Cc_rel,I_rel,J_rel);  
  

  // Check of inverse L bits

  printf("Check of inverse L pieces...\n");
  // c has original unpolluted matrix.
  getCholeskySafe(C,n);
  // clear upper right tri of L
  for (i=0;i<n;i++) {
	for (j=i+1;j<n;j++) C[i][j] =0.0;
  }

  double** L = dmatrix(0,n-1,0,n-1);
  memcpy(&L[0][0],&C[0][0],n*n*sizeof(double));
  modified_Cholesky_left_inverse(n,C);

  if (n<10)
	{
	  printf("L:\n");
	  debug_print_matrix((void**)L,n,'d');
	  printf("Cholesky L inverse:\n");
	  debug_print_matrix((void**)C,n,'d');
	}
  // form L.L^{-1}
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=L[i][k]*C[k][j];
	}
  }
  for (i=0,max=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  max=DMAX(max,fabs(B[i][j]-((i==j)?1:0)));
	}
  }
  printf("Max. error in L.L^{-1} - I = %lg\n",max);
  
  // form C^{-1} = L^{-1}T.L^{-1}
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=C[k][i]*C[k][j];
	}
  }
  printf("Cholesky C inverse:\n");
  debug_print_matrix((void**)B,n,'d');

  
  // form C^{-1}.C
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  for (C[i][j]=0,k=0;k<n;k++) C[i][j]+=B[i][k]*D[k][j];
	}
  }
  for (i=0,max=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  max=DMAX(max,fabs(C[i][j]-((i==j)?1:0)));
	}
  }
  printf("Max. error in C^{-1}.C - I = %lg\n",max);
  return 0;
}
