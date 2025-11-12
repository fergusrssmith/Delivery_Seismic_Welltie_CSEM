/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the same subdirectory or archive
   as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.           

   Modified (enlarged) by James Gunning, CSIRO Petroleum
   Modified (enlarged) by Chris White, LSU Petroleum Engineering


    */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	/* exit(1); */ 
	abort();
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
} 

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a integer 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;

	/* allocate pointers to pointers to rows */
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free an int f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

/* EXTRA CSIRO STUFF */

short *svector(long nl, long nh)
/* allocate a short vector with subscript range v[nl..nh] */
{
	short *v;

	v=(short *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(short)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_svector(short *v, long nl, long nh)
/* free a short vector allocated with svector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

short **smatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a short matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	short **m;

	/* allocate pointers to rows */
	m=(short **) malloc((size_t)((nrow+NR_END)*sizeof(short*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(short *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_smatrix(short **m, long nrl, long nrh, long ncl, long nch)
/* free an short matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

short ***s3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a short 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	short ***t;

	/* allocate pointers to pointers to rows */
	t=(short ***) malloc((size_t)((nrow+NR_END)*sizeof(short**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(short **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(short *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(short)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_s3tensor(short ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh) 
/* free an short 3tensor allocated by s3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


double ***d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh) 
{
  /* allocate a double 3tensor with rang t[nrl..nrh][ncl..nch][ndl..ndh] */
  long l,s, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  double ***t;

  /* allocate pointers to rows */
  t = (double ***) malloc( (size_t)( (nrow+NR_END)*sizeof(double**) ) );
  if( !t ) nrerror( "allocation failure 1 in d3tensor()" );
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if( !t[nrl] ) nrerror( "allocation failure 2 in d3tensor()" );
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double*)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if( !t[nrl][ncl] ) nrerror( "allocation failure 3 in d3tensor()" );
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for( s=ncl+1; s<=nch; s++ ) t[nrl][s] = t[nrl][s-1] + ndep;
  for( l=nrl+1; l<=nrh; l++ ) {
	t[l] = t[l-1] + ncol;
	t[l][ncl] = t[l-1][ncl] + ncol*ndep;
	for( s=ncl+1; s<=nch; s++ ) t[l][s] = t[l][s-1] + ndep;
  }
  /* return pointer to array of pointers to rows */
  return t;
}



#if 0
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    double ***t;

    /* allocate pointers to pointers to rows */
    t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
    if (!t) nrerror("allocation failure 1 in d3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
    if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}
#endif



void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}

/* -EXTRA CSIRO STUFF- */

/* EXTRA LSU STUFF */


char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;

	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof( char *)));
	if (!m) nrerror("allocation failure 1 in cmatrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
		   long npl, long nph)
/* allocate a double 4tensor with range 
   t[nrl..nrh][ncl..nch][ndl..ndh][npl..nph] */

/* Basically the same set-up as NR, but this way of extending seemed
   more natural to me... */
{
    long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,npar=nph-npl+1;
    long n;
    double ****t;
    double *t1;
    double **t2;
    double ***t3;

    /* allocate pointers to pointers to rows (YZP slices) */
    t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
    if (!t) nrerror("allocation failure 1 in d4tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
    if (!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
 
    /* allocate pointers to columns and set pointers to them */
    t[nrl][ncl] = (double **) malloc((size_t)((nrow*ncol*ndep+NR_END)*
					      sizeof(double*)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in d4tensor");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    
    t[nrl][ncl][ndl]=
      (double *) malloc((size_t)((nrow*ncol*ndep*npar+NR_END)*sizeof(double)));
    if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
    t[nrl][ncl][ndl] += NR_END;
    t[nrl][ncl][ndl] -= npl;

    
    /* pointers to rows */
    t3 = &(t[nrl][ncl])-ncl;
    for (i=nrl;i<=nrh;i++) {
      t[i] = t3;
      t3+=ncol;
    }

    /* pointers to columns */
    t2 = &(t[nrl][ncl][ndl])-ndl;
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
	t[i][j]=t2;
	t2+=ndep;
      }
    }

    /* pointers to cells values */
    t1 = t[nrl][ncl][ndl];
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
		for (k=ndl;k<=ndh;k++) {
		  t[i][j][k]=t1;
		  t1+=npar;
		}
      }
    }


    return (t);
}

void free_d4tensor(double ****t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh, long npl, long nph)
/* free a double 4tensor allocated by d4tensor() */
{
  /*    free((FREE_ARG) (t[nrl][ncl][ndl]+npl-NR_END));    */
    free((FREE_ARG) (t[nrl][ncl][ndl]+npl-NR_END));    
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}

float ****f4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
		   long npl, long nph)
/* allocate a double 4tensor with range 
   t[nrl..nrh][ncl..nch][ndl..ndh][npl..nph] */

/* Basically the same set-up as NR, but this way of extending seemed
   more natural to me... */
{
    long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,npar=nph-npl+1;
    long n;
    float ****t;
    float *t1;
    float **t2;
    float ***t3;

    /* allocate pointers to pointers to rows (YZP slices) */
    t=(float ****) malloc((size_t)((nrow+NR_END)*sizeof(float***)));
    if (!t) nrerror("allocation failure 1 in f4tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(float ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float**)));
    if (!t[nrl]) nrerror("allocation failure 2 in f4tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
 
    /* allocate pointers to columns and set pointers to them */
    t[nrl][ncl] = (float **) malloc((size_t)((nrow*ncol*ndep+NR_END)*
					      sizeof(float*)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in f4tensor");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    
    t[nrl][ncl][ndl]=
      (float *) malloc((size_t)((nrow*ncol*ndep*npar+NR_END)*sizeof(float)));
    if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in f4tensor()");
    t[nrl][ncl][ndl] += NR_END;
    t[nrl][ncl][ndl] -= npl;

    
    /* pointers to rows */
    t3 = &(t[nrl][ncl])-ncl;
    for (i=nrl;i<=nrh;i++) {
      t[i] = t3;
      t3+=ncol;
    }

    /* pointers to columns */
    t2 = &(t[nrl][ncl][ndl])-ndl;
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
	t[i][j]=t2;
	t2+=ndep;
      }
    }

    /* pointers to cells values */
    t1 = t[nrl][ncl][ndl];
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
		for (k=ndl;k<=ndh;k++) {
		  t[i][j][k]=t1;
		  t1+=npar;
		}
      }
    }


    return (t);
}

void free_f4tensor(float ****t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh, long npl, long nph)
/* free a float 4tensor allocated by f4tensor() */
{
    free((FREE_ARG) (t[nrl][ncl][ndl]+npl-NR_END));    
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}

int ****i4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
		   long npl, long nph)
/* allocate an int 4tensor with range 
   t[nrl..nrh][ncl..nch][ndl..ndh][npl..nph] */

/* Basically the same set-up as NR, but this way of extending seemed
   more natural to me... */
{
    long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,npar=nph-npl+1;
    long n;
    int ****t;
    int *t1;
    int **t2;
    int ***t3;

    /* allocate pointers to pointers to rows (YZP slices) */
    t=(int ****) malloc((size_t)((nrow+NR_END)*sizeof(int***)));
    if (!t) nrerror("allocation failure 1 in i4tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(int ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int**)));
    if (!t[nrl]) nrerror("allocation failure 2 in i4tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
 
    /* allocate pointers to columns and set pointers to them */
    t[nrl][ncl] = (int **) malloc((size_t)((nrow*ncol*ndep+NR_END)*
					      sizeof(int*)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in i4tensor");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    
    t[nrl][ncl][ndl]=
      (int *) malloc((size_t)((nrow*ncol*ndep*npar+NR_END)*sizeof(int)));
    if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in i4tensor()");
    t[nrl][ncl][ndl] += NR_END;
    t[nrl][ncl][ndl] -= npl;

    
    /* pointers to rows */
    t3 = &(t[nrl][ncl])-ncl;
    for (i=nrl;i<=nrh;i++) {
      t[i] = t3;
      t3+=ncol;
    }

    /* pointers to columns */
    t2 = &(t[nrl][ncl][ndl])-ndl;
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
	t[i][j]=t2;
	t2+=ndep;
      }
    }

    /* pointers to cells values */
    t1 = t[nrl][ncl][ndl];
    for (i=nrl;i<=nrh;i++) {
      for (j=ncl;j<=nch;j++) {
		for (k=ndl;k<=ndh;k++) {
		  t[i][j][k]=t1;
		  t1+=npar;
		}
      }
    }


    return (t);
}

void free_i4tensor(int ****t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh, long npl, long nph)
/* free an int i4tensor() */
{
    free((FREE_ARG) (t[nrl][ncl][ndl]+npl-NR_END));    
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}
/* -EXTRA LSU STUFF- */
