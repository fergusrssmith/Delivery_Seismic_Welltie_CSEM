#include <jni.h>
#include <stdio.h>
#include <math.h>

#include "au_csiro_JamesGunning_lib_util_BlasLapackUtil.h"

#ifdef __APPLE__
  #include   <Accelerate/Accelerate.h>
#else
  #include "cblas.h"
#endif


/* set to try to avoid data copies */
#define USE_GETPRIMITIVEARRAYCRITICAL 1

JNIEXPORT jboolean JNICALL Java_au_csiro_JamesGunning_lib_util_BlasLapackUtil_dgemvNative
  (JNIEnv *e, jclass o, jint M, jint N, jboolean transA, jdouble alpha, jdoubleArray A, jdoubleArray X, jdouble beta, jdoubleArray Y)
{
  jboolean isCopy=JNI_FALSE;
  jdouble *x,*y,*a;

#if USE_GETPRIMITIVEARRAYCRITICAL
  // A here is M*N, unrolled
  a = (jdouble*) ((*e)->GetPrimitiveArrayCritical(e,A,&isCopy));
  if (a == NULL) 
	{
	  fprintf(stderr,"FAILED GetPrimitiveArrayCritical() for a in dgemvNative()\n");
	  (*e)->ReleasePrimitiveArrayCritical(e,A,a,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  a = (*e)->GetDoubleArrayElements(e,A, NULL);
  if (a == NULL) return (jboolean)0; /* exception occurred */
#endif



#if USE_GETPRIMITIVEARRAYCRITICAL
  x = (jdouble*) ((*e)->GetPrimitiveArrayCritical(e,X,&isCopy));
  if (x == NULL) 
	{
	  fprintf(stderr,"FAILED GetPrimitiveArrayCritical() for x in dgemvNative()\n");
	  (*e)->ReleasePrimitiveArrayCritical(e,X,x,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  x = (*e)->GetDoubleArrayElements(e,X, NULL);
  if (x == NULL) return (jboolean)0; /* exception occurred */
#endif

#if USE_GETPRIMITIVEARRAYCRITICAL
  y = (jdouble*) ((*e)->GetPrimitiveArrayCritical(e,Y,&isCopy));
  if (y == NULL) 
	{
	  fprintf(stderr,"FAILED GetPrimitiveArrayCritical() for y in dgemvNative()\n");
	  (*e)->ReleasePrimitiveArrayCritical(e,Y,y,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  y = (*e)->GetDoubleArrayElements(e,Y, NULL);
  if (y == NULL) return (jboolean)0; /* exception occurred */
#endif

  
  /*
	void cblas_dgemv(const enum CBLAS_ORDER Order,
	const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
	const double alpha, const double *A, const int lda,
	const double *X, const int incX, const double beta,
	double *Y, const int incY);

	A = m rows x n cols
	lda = "leading" dimension of A = #cols if A is row-major (i.e. usual C-style),
  */

  cblas_dgemv(CblasRowMajor,transA ? CblasTrans : CblasNoTrans, M,N,
			  alpha , a, N, x, 1, beta, y, 1);

#if USE_GETPRIMITIVEARRAYCRITICAL

  (*e)->ReleasePrimitiveArrayCritical(e,A,a,0);
  (*e)->ReleasePrimitiveArrayCritical(e,X,x,0);
  (*e)->ReleasePrimitiveArrayCritical(e,Y,y,0);

#else

  (*e)->ReleaseDoubleArrayElements(e,A,a,0);
  (*e)->ReleaseDoubleArrayElements(e,X,x,0);
  (*e)->ReleaseDoubleArrayElements(e,Y,y,0);

#endif

  return JNI_TRUE;
}

#undef USE_GETPRIMITIVEARRAYCRITICAL 
