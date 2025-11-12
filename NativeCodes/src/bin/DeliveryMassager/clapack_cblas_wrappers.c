#ifdef __APPLE__
#include   <Accelerate/Accelerate.h>
#else
#include "../../lib/ATLAS/include/cblas.h"
#include "../../lib/ATLAS/include/clapack.h"
#endif

/*
  Set of wrappers to isolate lapack dependencies to this file and enable platform variations.
  e.g. clapack comes in OSX's "Accelerate" framework
  For casting, see also http://developer.apple.com/macosx/64bit.html
  and on OSX,
  /System/Library/Frameworks/vecLib.framework/Headers/clapack.h

  Note also OSX calls are to machine-translated lapack, so these
  calls expect colMajor matrices. All the matrices we use are symmetric, so this
  doesn't matter here, but in general
  
  dgetrf_(M,N,A,lda,ipiv, &info); // OSX clapack
  
  should be equivalent to
  
  clapack_dgetrf(CblasColMajor,M,N,A,lda,ipiv);  // ATLAS
  
*/

/* Row major only dgetrf() */
int clapack_wrapper_dgetrf(int M, int N,
						  double *A, int lda, int *ipiv)
{
#ifdef __APPLE__
  /*
	int dgetrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
	lda, __CLPK_integer *ipiv, __CLPK_integer *info);
  */ 
  __CLPK_integer info;
  return dgetrf_((__CLPK_integer*)&M, (__CLPK_integer*)&N, A, (__CLPK_integer*)&lda,(__CLPK_integer*)ipiv, &info);

#else
  /* ATLAS */
  /*
	int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
	double *A, const int lda, int *ipiv);
  */
 return clapack_dgetrf(CblasRowMajor,M,N,A,lda,ipiv);
#endif
}

/* Row major only dgetrs(): trans should be 'T' or 'N' */
int clapack_wrapper_dgetrs(char trans, int N, int NRHS, double *A, int lda,
						  int *ipiv, double *B, int ldb)
{

#ifdef __APPLE__
  /*
	int dgetrs_(char *trans, __CLPK_integer *n, __CLPK_integer *nrhs, 
	__CLPK_doublereal *a, __CLPK_integer *lda, __CLPK_integer *ipiv, __CLPK_doublereal *b, __CLPK_integer *
	ldb, __CLPK_integer *info);
  */
 __CLPK_integer info;
 return dgetrs_(&trans,(__CLPK_integer*)&N,(__CLPK_integer*)&NRHS, A, (__CLPK_integer*)&lda, 
				(__CLPK_integer*)ipiv, B, (__CLPK_integer*)&ldb, (__CLPK_integer*)&info);

#else
  /* ATLAS */
 /*
   int clapack_dgetrs
   (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const double *A, const int lda,
    const int *ipiv, double *B, const int ldb);
 */
 return clapack_dgetrs(CblasRowMajor, ((trans=='T')||(trans=='t')) ? CblasTrans : CblasNoTrans,
					   N,NRHS,A,lda,ipiv,B,ldb);

#endif  
}

/* Row major spotrf wrapper */
int clapack_wrapper_spotrf(char uplo, int N, float *A, int lda)
{

#ifdef __APPLE__
  /*
	int spotrf_(char *uplo, __CLPK_integer *n, __CLPK_real *a, __CLPK_integer *lda, 
	__CLPK_integer *info);
  */
  __CLPK_integer info;
  return spotrf_(&uplo,(__CLPK_integer*)&N,A,(__CLPK_integer*)&lda,&info);

#else
  /* ATLAS */
  /*
	int clapack_spotrf(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
	const int N, float *A, const int lda);
  */
  return clapack_spotrf(CblasRowMajor, ((uplo=='U')||(uplo=='u')) ? CblasUpper : CblasLower,
						N,A,lda);
#endif
}

/* Row major dpotrf wrapper */ 
int clapack_wrapper_dpotrf(char uplo, int N, double *A, int lda)
{

#ifdef __APPLE__
  /*
	int dpotrf_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
	lda, __CLPK_integer *info);

  */
  __CLPK_integer info;
  return dpotrf_(&uplo,(__CLPK_integer*)&N,A,(__CLPK_integer*)&lda,(__CLPK_integer*)&info);

#else
  /* ATLAS */
  /*
	int clapack_dpotrf(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
	const int N, double *A, const int lda);	
   */
  return clapack_dpotrf(CblasRowMajor, ((uplo=='U')||(uplo=='u')) ? CblasUpper : CblasLower,
						N,A,lda);
#endif
}

/* Row major dpotri wrapper */
int clapack_wrapper_dpotri(char uplo, int N, double *A, int lda)
{

#ifdef __APPLE__
  /*
	int dpotri_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
	lda, __CLPK_integer *info);
  */
  __CLPK_integer info;
  return dpotri_(&uplo, (__CLPK_integer*)&N, A, (__CLPK_integer*)&lda, (__CLPK_integer*)&info);
#else
  /* ATLAS */
  /*
	int clapack_dpotri(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
	const int N, double *A, const int lda);
   */
  return clapack_dpotri(CblasRowMajor,((uplo=='U')||(uplo=='u')) ? CblasUpper : CblasLower,
						N, A, lda);
#endif
}

/* dgemv wrapper, rowMajor only */
void cblas_wrapper_dgemv(char trans, int M, int N,
						double alpha, double *A, int lda,
						double *X, int incX, double beta,
						double *Y, int incY)
{
  /* same on both platforms... */

  cblas_dgemv(CblasRowMajor,((trans=='T')||(trans=='t')) ? CblasTrans : CblasNoTrans,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}
/* dgemm wrapper, rowMajor only */
void cblas_wrapper_dgemm(char transA,
                 char transB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc)
{

  cblas_dgemm(CblasRowMajor, 
			  ((transA=='T')||(transA=='t')) ? CblasTrans : CblasNoTrans,
              ((transB=='T')||(transB=='t')) ? CblasTrans : CblasNoTrans,
			  M,  N,
			  K,  alpha,  A,
			  lda,  B,  ldb,
			  beta, C,  ldc);

}
