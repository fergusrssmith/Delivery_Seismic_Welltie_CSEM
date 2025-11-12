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
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>



// #include <malloc.h> only if malloc_stats below needs to be invoked (tends to be platform specific)


#include "NR/nrutil.h"

#include "platform_specifics.h"
#include "modified_left_Cholesky.h" 
#include "util.h"
#include "errorHandler.h"

static THREAD int bad_matrix_ctr=0;




int is_big_endian()
{
  int x = 1;
  if (*(char *)&x == 1) return 0; // little
  else	return 1;                 // big
}



/* Pade constants in Kennedy and Gentle */
#define p0 -0.322232431088
#define p1 -1.0
#define p2 -0.342242088547
#define p3 -0.0204231210245
#define p4 -0.0000453642210148
#define q0 0.0993484626060
#define q1 0.588581570495
#define q2 0.531103462366
#define q3 0.103537752850
#define q4 0.0038560700634
float inverseG(float p)
/* inverse standard normal cdf, as per Kennedy and Gentle. p is a quantile 0<p<1, 
   returns the normal deviate y with corresponding quantile */
{
  float		y,pp,xp;

  pp=p;
  if (p>0.5) pp=1.0-p;
  y	 = (float)sqrt(log(1. /(pp*pp)));

  xp = (float)(y + ((((y * p4 + p3) * y + p2) * y + p1) * y + p0)
				   / ((((y * q4 + q3) * y + q2) * y + q1) * y + q0));
  if (p==pp) xp = -xp;
  return xp;
} /* inverseG() */
#undef p0
#undef p1
#undef p2
#undef p3
#undef p4
#undef q0
#undef q1
#undef q2
#undef q3
#undef q4


// memory usage info
// uncomment #include <malloc.h> only if malloc_stats below needs to be invoked (tends to be platform specific)
void show_malloc_stats()
{
  // malloc_stats();
}

/*
  method to detect co-located neighbours in a kriging system and adjust
  matrix. Co-located neighbours will have the same conditioning data (or 
  *ought to*, so one of them can be eliminated. This method
  send the duplicated point to infinity (so the off diagonals of C
  attached to that point  ->0, and the same for the RHS entry

  For SK, set N = size of system. For OK with covariance form, set N=size-1
*/
void cleanKrigingSystemOfDuplicateNeighbours(int N,double** C,double* RHSkrig)
{
  int k,kk,deleteRow;
  char cleaned='\0';
  double Cscale;

  // Clean up for co-located neighbours (rank deficiencies)
  for (k=0; k < N; k++) 
	{
	  deleteRow=0;
	  Cscale =C[k][k];
	  for (kk=0; kk < k ; kk++) 
		{
		  if (C[k][kk] >= Cscale)
			{
			  deleteRow=1;
			  cleaned=(char)1;
			  break; 
			}
		}
	  if (deleteRow)
		{
		  // nbr k and kk are redundant: kill k by sending point k->infinity
		  // (put zeros in all off-diagonals attached to kk)
		  for (kk=0; kk < k ; kk++) C[k][kk]=C[kk][k]=0.0; 
		  for (kk=(k+1); kk < N ; kk++) C[k][kk]=C[kk][k]=0.0; 
		  RHSkrig[k]=0.0;
		}
	}
  // if (cleaned) printf("SK cleaned!\n");
  return;
}


/*
  Test routine for cholesky safe methods - not a library call. Used for
  reading matrices from files and examining.
  -1 for binary read
  0 for ascii
  1 for made-up matrix
*/
void TestgetLDLTdecompositionSafe(int n,char *filename) 
{
  void debug_print_matrix(void **m,int n,char type);
  int i,j,k;
  double lambda,tmp;
  double **A,**R,**C,**B;
  int preScale= 0;

  // -1 for binary read
  // 0 for ascii
  // 1 for made-up
  if (n<=0)
	{
	  if (n==0)
		{
		  // read matrix from ascii stream (n*n numbers)
		  double f,*ptr;
		  int cnt=0;
		  FILE *in;
		  fpos_t pos;

		  if (!(in=fopen(filename,"r"))) {printf("Can't open %s\n",filename); return;}
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
		  // read matrix from bin stream (n*n numbers)
		  double f,*ptr;
		  FILE *in;

		  if (!(in=fopen(filename,"rb"))) {printf("Can't open %s\n",filename); return;}
		  fseek(in, 0L, SEEK_END);
		  n = (int)sqrt(ftell(in)/sizeof(double));
		  rewind(in);
		  A = dmatrix(0,n-1,0,n-1);
		  ptr=&(A[0][0]);
		  fread(ptr,sizeof(double),n*n,in);
		  fclose(in);
 		}
	} else {
	  // make a funky just-positive definite matrix.

	  A = dmatrix(0,n-1,0,n-1);
	  double* scale = calloc(n,(size_t)sizeof(double));
	  for (i = 0; i < n; i++) 
		{
		  scale[i] = (i<n/3) ? 1 : ((i>(2*n/3))? 1 : 1e-9);
		  //scale[i]=4;
		}
  
	  for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
		  A[i][j] = scale[i]*scale[j]*exp(-fabs(i-j)/100.0); 
		}
	  }
	}

  if (n<10) 
	{
	  printf("Original A matrix:\n");
	  debug_print_matrix((void**)A,n,'d');
	}

  // take a copy
  C = dmatrix(0,n-1,0,n-1);
  memcpy(&C[0][0],&A[0][0],n*n*sizeof(double));

  // perform Left Cholesky
  getCholeskySafe(A,n);
  for (i=0;i<n;i++) 
	for (j=i+1;j<n;j++) 
	  A[i][j]=0.0;

  if (n<10)
	{
	  printf("L=\n");
	  debug_print_matrix((void**)A,n,'d');
	}
  
  // Check multiplication
  // B = A.A^T
  B = dmatrix(0,n-1,0,n-1);
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=A[i][k]*A[j][k];
	}
  }
  // permute
  double max=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  max=DMAX(max, fabs(B[i][j]  - C[i][j]));
	}
  }
  printf("\nMax abs diff. error after multiplying out = %g\n",max);  
}


#define PERFORM_RESTORE_TESTS 0

void getCholeskySafe(double **a, int n)
{
  int i,j,info;
  int *perm = ivector(0,n-1);

#if PERFORM_RESTORE_TESTS
  double **B,**C; 
  C = dmatrix(0,n-1,0,n-1);
  memcpy(&C[0][0],&a[0][0],n*n*sizeof(double));
#endif
 
  modified_Cholesky_left_decomposition(n,a,perm,1);
  // test if rows were permuted. If not, we are done. 
  for (i=0,info=0; i<n; i++) if (perm[i]!=i) {info=1;break;}
  // If permuted, restore matrix and do conventional bulletproof Cholesky
  if (info) 
	{
	  // Permuted. Reconstruct and use bulletproof method with normal cholesky
	  modified_Cholesky_restore(n,a,perm);
	  getLDLTdecompositionSafe(a,n); 
	} 
  else
	{
	  // zero upper triangle for safety
	  for (i=0;i<n;i++) {
		for (j=i+1;j<n;j++) a[i][j]=0.0;
	  }
	}

  // we are done. Free up and go on 
  free_ivector(perm,0,n-1);
  perm = NULL;

#if PERFORM_RESTORE_TESTS
  int k;
  double err=0,abserr=0,*sc;
  // Check multiplication
  // B = a.a^T
  B = dmatrix(0,n-1,0,n-1);
  sc = dvector(0,n-1);

  // only lower tri of a is usable....
  for (i=0;i<n;i++) 
	for (j=i+1;j<n;j++) 
	  a[i][j]=0.0;
  for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=a[i][k]*a[j][k];
	}
	sc[i] = sqrt(DMAX(0,C[i][i]));
  }
  double max=0;
  for (i=0;i<n;i++) { 
	for (j=0;j<=i;j++) {
	  err = fabs(B[i][j] - C[i][j]);
	  abserr=DMAX(abserr,err);
	  if (sc[i]>0)
		{
		  err /= (sc[i] * ((sc[j]>0)? sc[j] : sc[i]));
		}
	  else err /= ((sc[j]>0)? (sc[j]*sc[j]) : 1e50);
	  max=DMAX(max, err);
	}
  }

  printf("\nMax error in getCholeskySafe (Cholesky after restore): abs = %g , relative = %g\n",abserr,max);  
  free_dvector(sc,0,n-1);
  free_dmatrix(B,0,n-1,0,n-1);
  if (max>0.01 && abserr>0.01)
	{
	  FILE *in;
	  char* name;
	  asprintf_w(&name,"badMatrix.%d.bin",bad_matrix_ctr++);
	  if (!(in=fopen(name,"wb"))) {printf("Can't open %s\n",name); exit(0);}
	  fwrite(&(C[0][0]),sizeof(double),n*n,in);
	  fclose(in); 
	  printf("Wrote %s\n",name);
	  free(name);
	  // exit(0);
	}
  free_dmatrix(C,0,n-1,0,n-1);
#endif

}

#undef PERFORM_RESTORE_TESTS


#define TINY 1.0e-20;
#define PRESCALE 0
#define PERFORM_LDL_TESTS 0

/*
  Zero offset no-pivoting LDL^T decomp routine. Used to get Cholesky factors without
  permutations of matrices that are close to symmetric positive definite. 
  Returns approximate L factor in-place
  Input matrix has only to be approximately positive definite.
  Inputs: a[0...n-1][0...n-1],   symetric pos def. matrix. Only lower triangle required.
  n = dimension of matrix
*/

void getLDLTdecompositionSafe(double **a, int n)
{
  int i,j,k;
  double sum;
  double *d,*d_piv;

  d = dvector(0,n-1); 
  d_piv = dvector(0,n-1);

#if PERFORM_LDL_TESTS
  double **B,**C;
  C = dmatrix(0,n-1,0,n-1);
  memcpy(&C[0][0],&a[0][0],n*n*sizeof(double));
#endif


#if PRESCALE
  double* s;
  s = dvector(0,n-1);

  for (i=0;i<n;i++) 
	{
	  s[i] = (a[i][i]>0)? sqrt(a[i][i]) : 0.0;
	  d[i] = (s[i]!=0)? s[i] : 1.0;
	}
  // must scale all the matrix
  for (i=0;i<n;i++) 
	{
	  for (j=0;j<=i;j++) 
		{
		  a[i][j] /= (d[i]*d[j]); 
		  a[j][i] = a[i][j];
		}
	  a[i][i]=1.0;
	}
#endif

  for (i=0;i<n;i++) 
	{
	  for (j=i;j<n;j++) 
		{
		  for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k]*d[k]; // eqn 4
		  if (i == j) 
			{
			  d[i] = sum;
			  d_piv[i] = (d[i] != 0) ? d[i] : TINY;			  
			}
		  else 
			{
			  a[j][i] = (d[i]!=0) ? sum/d_piv[i] : 0.0;
			}
		}
	}

  // reinstate L <- L.diag(d)^1/2. Set -ve elements of d to 0
  // double max=0,min=0;
  for (i=0;i<n;i++) 
	{
	  // max = DMAX(d[i],max);
	  // min = DMIN(d[i],min);
	  d[i] = (d[i]>0)? sqrt(d[i]) : 0.0;
	  a[i][i] = 1.0; // diag (L) =1.0
	}

  for (i=0;i<n;i++) 
	{
	  for (j=0;j<=i;j++) 
		{
		  a[i][j] *= d[j];
		}
	  // zero upper tri...
	  for (j=i+1;j<n;j++) a[i][j]=0.0;
	}

#if PRESCALE
  for (i=0;i<n;i++) 
	{
	  for (j=0;j<=i;j++) 
		{
		  a[i][j] *= s[i]; 
		}
	}
  free_dvector(s,0,n-1);
  s=NULL;
#endif


  free_dvector(d,0,n-1); d=NULL;
  free_dvector(d_piv,0,n-1); d_piv=NULL;


#if PERFORM_LDL_TESTS
  double err=0,*sc;
  // Check multiplication
  // B = a.a^T
  B = dmatrix(0,n-1,0,n-1);
  sc = dvector(0,n-1);

  // only lower tri of a is usable....
  for (i=0;i<n;i++) 
	for (j=i+1;j<n;j++) 
	  a[i][j]=0.0;
  for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
	  for (B[i][j]=0,k=0;k<n;k++) B[i][j]+=a[i][k]*a[j][k];
	}
	sc[i] = sqrt(DMAX(0,C[i][i]));
  }
  double max=0;
  for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
	  err = fabs(B[i][j] - C[i][j]);
	  if (sc[i]>0)
		{
		  err /= (sc[i] * ((sc[j]>0)? sc[j] : sc[i]));
		}
	  else err /= ((sc[j]>0)? sc[j] *sc[j] : 1e50);
	  max=DMAX(max, err);
	}
  }

  printf("\nMax relative error in getLDLTdecompositionSafe = %g\n",max);  
  free_dvector(sc,0,n-1); sc=NULL;
  free_dmatrix(B,0,n-1,0,n-1); B=NULL;
  free_dmatrix(C,0,n-1,0,n-1); C=NULL;
#endif

}
#undef TINY
#undef PERFORM_LDL_TESTS 

// random nunber generators
/*
Uniform, similar to NR ran1()
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float random1(long *seed)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*seed <= 0 || !iy) {
        if (-(*seed) < 1) *seed=1;
        else *seed = -(*seed);
        for (j=NTAB+7;j>=0;j--) {
            k=(*seed)/IQ;
            *seed=IA*(*seed-k*IQ)-IR*k;
            if (*seed < 0) *seed += IM;
            if (j < NTAB) iv[j] = *seed;
        }
        iy=iv[0];
    }
    k=(*seed)/IQ;
    *seed=IA*(*seed-k*IQ)-IR*k;
    if (*seed < 0) *seed += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *seed;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* 
   Gaussian deviates from N(0,1) 
   Similar to NR gasdev().
*/

float Gaussian_deviate(long *seed)
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*random1(seed)-1.0;
			v2=2.0*random1(seed)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}



#define CONVEX
#define HYBRID 
#define WINDING

#define X	0
#define Y	1

#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#ifndef HUGE
#define HUGE	1.797693134862315e+308
#endif

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* test if a & b are within epsilon.  Favors cases where a < b */
#define Near(a,b,eps)	( ((b)-(eps)<(a)) && ((a)-(eps)<(b)) )


/* ======= Crossings algorithm for point-in-polygon test ====================== */

/* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to
 * compare vertex Y values to the testing point's Y and quickly discard
 * edges which are entirely to one side of the test ray.
 *
 * Input 2D polygon _pgon_ with _numverts_ number of vertices and test point
 * _point_, returns 1 if inside, 0 if outside.	WINDING and CONVEX can be
 * defined for this test.
 */
int CrossingsTest(double	pgon[][2], int numverts, double	point[2] )
{
#ifdef	WINDING
  register int	crossings ;
#endif
  register int	j, yflag0, yflag1, inside_flag, xflag0 ;
  register double ty, tx, *vtx0, *vtx1 ;
#ifdef	CONVEX
  register int	line_flag ;
#endif

  tx = point[X] ;
  ty = point[Y] ;

  vtx0 = pgon[numverts-1] ;
  /* get test bit for above/below X axis */
  yflag0 = ( vtx0[Y] >= ty ) ;
  vtx1 = pgon[0] ;

#ifdef	WINDING
  crossings = 0 ;
#else
  inside_flag = 0 ;
#endif
#ifdef	CONVEX
  line_flag = 0 ;
#endif
  for ( j = numverts+1 ; --j ; ) {

	yflag1 = ( vtx1[Y] >= ty ) ;
	/* check if endpoints straddle (are on opposite sides) of X axis
	 * (i.e. the Y's differ); if so, +X ray could intersect this edge.
	 */
	if ( yflag0 != yflag1 ) {
	  xflag0 = ( vtx0[X] >= tx ) ;
	  /* check if endpoints are on same side of the Y axis (i.e. X's
	   * are the same); if so, it's easy to test if edge hits or misses.
	   */
	  if ( xflag0 == ( vtx1[X] >= tx ) ) {

		/* if edge's X values both right of the point, must hit */
#ifdef	WINDING
		if ( xflag0 ) crossings += ( yflag0 ? -1 : 1 ) ;
#else
		if ( xflag0 ) inside_flag = !inside_flag ;
#endif
	  } else {
		/* compute intersection of pgon segment with +X ray, note
		 * if >= point's X; if so, the ray hits it.
		 */
		if ( (vtx1[X] - (vtx1[Y]-ty)*
			  ( vtx0[X]-vtx1[X])/(vtx0[Y]-vtx1[Y])) >= tx ) {
#ifdef	WINDING
		  crossings += ( yflag0 ? -1 : 1 ) ;
#else
		  inside_flag = !inside_flag ;
#endif
		}
	  }
#ifdef	CONVEX
	  /* if this is second edge hit, then done testing */
	  if ( line_flag ) goto Exit ;

	  /* note that one edge has been hit by the ray's line */
	  line_flag = TRUE ;
#endif
	}

	/* move to next pair of vertices, retaining info as possible */
	yflag0 = yflag1 ;
	vtx0 = vtx1 ;
	vtx1 += 2 ;
  }
#ifdef	CONVEX
 Exit: ;
#endif
#ifdef	WINDING
  /* test if crossings is not zero */
  inside_flag = (crossings != 0) ;
#endif

  return( inside_flag ) ;
}

#undef Near
#undef X
#undef Y
#undef CONVEX
#undef WINDING

void f4tensor_fill(float val,float ****f,int Nx,int Ny,int Nz,int Nw)
{
  int i,N;

  float *fp;

  fp = &(f[0][0][0][0]);

  N = Nx*Ny*Nz*Nw;
  for (i=0; i<N; i++) fp[i] = val;
  
}

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


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */
#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double quick_select(double arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP


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

void print_dvector0_mma(double *m,int n,char* filename)
{
  int i,j;
  FILE *out;
        
  if ((out = fopen(filename, "w")) == NULL) 
	{
	  printf("Warning: Can't open %s",filename);
	  return;
	}
  printf("Printing vector to %s\n",filename);
  fprintf(out,"{");
  for (j=0;j<n;j++)
	{
	  fprintf(out,"%15.15f",m[j]);
	  if (j<n-1)fprintf(out,", ");
	}
  fprintf(out,"}");
  fclose(out);
}

