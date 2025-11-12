#include <jni.h>
#include "au_csiro_JamesGunning_lib_Numerical_djbFFT_DJBFFT.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>
#include <fftr8.h>

#ifdef __APPLE__
#include   <Accelerate/Accelerate.h>
#endif

// Hack to work around glibc memcpy saga. Sigh
// http://www.win.tue.nl/~aeb/linux/misc/gcc-semibug.html
//__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");

/* set to try to avoid data copies */
#define USE_GETPRIMITIVEARRAYCRITICAL 1

/* set to 1 to test timings */
#define TEST_TIMES 0

#if TEST_TIMES
static double prep_time=0, fft_time = 0;
static int n_calls=0;
#endif

/*
  Convolve method using real convolutions from DJBFFT. Note strange
  interleaving required
*/
JNIEXPORT jboolean JNICALL Java_au_csiro_JamesGunning_lib_Numerical_djbFFT_DJBFFT_convolveNative
  (JNIEnv *env, jclass obj, jdoubleArray A, jint nA, jdoubleArray B, jint nB)
{
#if TEST_TIMES
  struct timeval tv;
  struct timezone tz;
  double dt,tn;
#endif

  extern int errno;
  int i,j;
  real8* b=NULL;
  jboolean isCopy=JNI_FALSE;
  jdouble* a=NULL;
  jdouble* Bptr=NULL;

#if TEST_TIMES
  gettimeofday(&tv,&tz);
  dt = (double)((tv.tv_sec*1e6+tv.tv_usec));
#endif


  if (nA!=nB) return (jboolean)0;
  if ((nA & -nA) != nA) return (jboolean)0; /* i.e., nA is not a power of 2 */
  if ((nB & -nB) != nB) return (jboolean)0; /* i.e., nB is not a power of 2 */
  if (nA > 8192) return (jboolean)0;


  /* test stuff */
 #if 0
	jdouble *ac=NULL;
	ac = (jdouble*) ((*env)->GetPrimitiveArrayCritical(env,A,&isCopy));
	printf("DJBFFT isCopy=%d\n",isCopy);
	(*env)->ReleasePrimitiveArrayCritical(env,A,ac,0);
#endif

  /* We need ptrs to the primitive arrays in the VM */
#if USE_GETPRIMITIVEARRAYCRITICAL
  a = (jdouble*) ((*env)->GetPrimitiveArrayCritical(env,A,&isCopy));
  if (a == NULL) 
	{
	  printf("FAILED GetPrimitiveArrayCritical() for A in convolveNative()\n");
	  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  a = (*env)->GetDoubleArrayElements(env,A, NULL);
  if (a == NULL) return (jboolean)0; /* exception occurred */
#endif
   
#if USE_GETPRIMITIVEARRAYCRITICAL
  Bptr = (jdouble*) ((*env)->GetPrimitiveArrayCritical(env,B,&isCopy));
  if (Bptr == NULL) 
	{
	  printf("FAILED GetPrimitiveArrayCritical() for B in convolveNative()\n");
	  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0); // free a above
	  (*env)->ReleasePrimitiveArrayCritical(env,B,Bptr,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  Bptr = (*env)->GetDoubleArrayElements(env,B, NULL);
  if (Bptr == NULL) return (jboolean)0; /* exception occurred */
#endif


  /* copy of wavelet.... */
  b = (real8 *) malloc(nA*sizeof(real8));
  if (b == NULL)
	{
	  /* memory failure, */
	  free(b);
#if USE_GETPRIMITIVEARRAYCRITICAL
	  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0);
	  (*env)->ReleasePrimitiveArrayCritical(env,B,Bptr,0);
#endif	  
	  printf("FAILED malloc in convolveNative()\n"); 
	  return (jboolean)0; /* malloc failure*/
	}



  /* interleaving. Rearrange A */
  memcpy(b,a,nA*sizeof(real8));
  for (i=j=0;i<nA;i+=2,j++) {a[i]=b[j];}
  for (i=1,j=nA/2;i<nA;i+=2,j++) {a[i]=b[j];}

  /* rearrange B (wavelet) into b. Leave wavelet unmangled */
  for (i=j=0;i<nA;i+=2,j++) {b[i]=Bptr[j];}
  for (i=1,j=nA/2;i<nA;i+=2,j++) {b[i]=Bptr[j];}

#if TEST_TIMES
  gettimeofday(&tv,&tz);
  prep_time += (tn = (double)((tv.tv_sec*1e6+tv.tv_usec))) - dt;
  dt = tn;
#endif

  /* convolution....*/
  switch (nA)
	{
	case 2:
	  {  
		fftr8_2(a);
		fftr8_2(b);
		fftr8_mul2(a,b);
		fftr8_scale2(a);
		fftr8_un2(a);
		break;
	  }
	case 4:
	  {  
		fftr8_4(a);
		fftr8_4(b);
		fftr8_mul4(a,b);
		fftr8_scale4(a);
		fftr8_un4(a);
		break;
	  }
	case 8:
	  {  
		fftr8_8(a);
		fftr8_8(b);
		fftr8_mul8(a,b);
		fftr8_scale8(a);
		fftr8_un8(a);
		break;
	  }
	case 16:
	  {  
		fftr8_16(a);
		fftr8_16(b);
		fftr8_mul16(a,b);
		fftr8_scale16(a);
		fftr8_un16(a);
		break;
	  }
	case 32:
	  {  
		fftr8_32(a);
		fftr8_32(b);
		fftr8_mul32(a,b);
		fftr8_scale32(a);
		fftr8_un32(a);
		break;
	  }
	case 64:
	  {  
		fftr8_64(a);
		fftr8_64(b);
		fftr8_mul64(a,b);
		fftr8_scale64(a);
		fftr8_un64(a);
		break;
	  }
	case 128:
	  {  
		fftr8_128(a);
		fftr8_128(b);
		fftr8_mul128(a,b);
		fftr8_scale128(a);
		fftr8_un128(a);
		break;
	  }
	case 256:
	  {  
		fftr8_256(a);
		fftr8_256(b);
		fftr8_mul256(a,b);
		fftr8_scale256(a);
		fftr8_un256(a);
		break;
	  }
	case 512:
	  {  
		fftr8_512(a);
		fftr8_512(b);
		fftr8_mul512(a,b);
		fftr8_scale512(a);
		fftr8_un512(a);
		break;
	  }
	case 1024:
	  {  
		fftr8_1024(a);
		fftr8_1024(b);
		fftr8_mul1024(a,b);
		fftr8_scale1024(a);
		fftr8_un1024(a);
		break;
	  }
	case 2048:
	  {  
		fftr8_2048(a);
		fftr8_2048(b);
		fftr8_mul2048(a,b);
		fftr8_scale2048(a);
		fftr8_un2048(a);
		break;
	  }
	case 4096:
	  {  
		fftr8_4096(a);
		fftr8_4096(b);
		fftr8_mul4096(a,b);
		fftr8_scale4096(a);
		fftr8_un4096(a);
		break;
	  }
	case 8192:
	  {  
		fftr8_8192(a);
		fftr8_8192(b);
		fftr8_mul8192(a,b);
		fftr8_scale8192(a);
		fftr8_un8192(a);
		break;
	  }
	default:
	  {}
	}
#if TEST_TIMES
  gettimeofday(&tv,&tz);
  fft_time += (tn = (double)((tv.tv_sec*1e6+tv.tv_usec))) - dt;
  dt = tn;
#endif

  /* undo interleaving. Copy result in a back to b */
  memcpy(b,a,nA*sizeof(real8));
  /* then stuff into a the correct ordering */
  for (i=j=0;i<nA;i+=2,j++) {a[j]=b[i];}
  for (i=1,j=nA/2;i<nA;i+=2,j++) {a[j]=b[i];}
  
 MEM_CLEAN:
  /* free mem */
#if USE_GETPRIMITIVEARRAYCRITICAL
  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0);
  (*env)->ReleasePrimitiveArrayCritical(env,B,Bptr,0);
#else
  (*env)->ReleaseDoubleArrayElements(env,A,a,0);
  (*env)->ReleaseDoubleArrayElements(env,B,Bptr,0);
#endif
  free(b);


#if TEST_TIMES
  gettimeofday(&tv,&tz);
  prep_time += (tn = (double)((tv.tv_sec*1e6+tv.tv_usec))) - dt;
  dt = tn;
  n_calls++;
  if (n_calls>10000) 
	{
	  printf("Ratio of fft (n=%d) to prep/cleanup time = %g\n", nA, fft_time/prep_time);
	  n_calls=0;
	}
#endif


  return (jboolean)1;
}

#if 0
/* Complex version..... */
#include <fftc8.h>

JNIEXPORT jboolean JNICALL Java_au_csiro_JamesGunning_lib_Numerical_djbFFT_DJBFFT_convolveNative
  (JNIEnv *env, jclass obj, jdoubleArray A, jint nA, jdoubleArray B, jint nB)
{
  extern int errno;
  int i,na,nb;
  complex8* a;
  complex8* b;
  jdouble* Aptr;
  jdouble* Bptr;

  na=(int)nA;
  nb=(int)nB;
  if (nA!=nB) return (jboolean)0;
  if ((nA & -nA) != nA) return (jboolean)0; /* i.e., nA is not a power of 2 */
  if ((nB & -nB) != nB) return (jboolean)0; /* i.e., nB is not a power of 2 */
  if (nA > 8192) return (jboolean)0;


  /*  We need ptrs to the primitive arrays in the VM */
  Aptr = (*env)->GetDoubleArrayElements(env,A, NULL);
  if (Aptr == NULL) return (jboolean)0; /* exception occurred */
   
  Bptr = (*env)->GetDoubleArrayElements(env,B, NULL);
  if (Bptr == NULL) return (jboolean)0; /* exception occurred */
  

  /*  get real data into complex form */
  a = (complex8 *) calloc(nA,sizeof(complex8));
  if (errno==ENOMEM)
	{
	  free(a);
		/* memory failure, */
	  return (jboolean)0; /* calloc failure*/
	}

  b = (complex8 *) calloc(nB,sizeof(complex8));
  if (errno==ENOMEM)
	{
	  free(b);
	  /* memory failure, */
	  return (jboolean)0; /* calloc failure*/
	}

  /*  copy data to real parts.... */

  for (i=0;i<nA;i++) a[i].re = (real8)Aptr[i];
  for (i=0;i<nB;i++) b[i].re = (real8)Bptr[i];

  /*  convolution.... */

  switch (nA)
	{
	case 2:
	  {  
		fftc8_2(a);
		fftc8_2(b);
		fftc8_mul2(a,b);
		fftc8_scale2(a);
		fftc8_un2(a);
		break;
	  }
	case 4:
	  {  
		fftc8_4(a);
		fftc8_4(b);
		fftc8_mul4(a,b);
		fftc8_scale4(a);
		fftc8_un4(a);
		break;
	  }
	case 8:
	  {  
		fftc8_8(a);
		fftc8_8(b);
		fftc8_mul8(a,b);
		fftc8_scale8(a);
		fftc8_un8(a);
		break;
	  }
	case 16:
	  {  
		fftc8_16(a);
		fftc8_16(b);
		fftc8_mul16(a,b);
		fftc8_scale16(a);
		fftc8_un16(a);
		break;
	  }
	case 32:
	  {  
		fftc8_32(a);
		fftc8_32(b);
		fftc8_mul32(a,b);
		fftc8_scale32(a);
		fftc8_un32(a);
		break;
	  }
	case 64:
	  {  
		fftc8_64(a);
		fftc8_64(b);
		fftc8_mul64(a,b);
		fftc8_scale64(a);
		fftc8_un64(a);
		break;
	  }
	case 128:
	  {  
		fftc8_128(a);
		fftc8_128(b);
		fftc8_mul128(a,b);
		fftc8_scale128(a);
		fftc8_un128(a);
		break;
	  }
	case 256:
	  {  
		fftc8_256(a);
		fftc8_256(b);
		fftc8_mul256(a,b);
		fftc8_scale256(a);
		fftc8_un256(a);
		break;
	  }
	case 512:
	  {  
		fftc8_512(a);
		fftc8_512(b);
		fftc8_mul512(a,b);
		fftc8_scale512(a);
		fftc8_un512(a);
		break;
	  }
	case 1024:
	  {  
		fftc8_1024(a);
		fftc8_1024(b);
		fftc8_mul1024(a,b);
		fftc8_scale1024(a);
		fftc8_un1024(a);
		break;
	  }
	case 2048:
	  {  
		fftc8_2048(a);
		fftc8_2048(b);
		fftc8_mul2048(a,b);
		fftc8_scale2048(a);
		fftc8_un2048(a);
		break;
	  }
	case 4096:
	  {  
		fftc8_4096(a);
		fftc8_4096(b);
		fftc8_mul4096(a,b);
		fftc8_scale4096(a);
		fftc8_un4096(a);
		break;
	  }
	case 8192:
	  {  
		fftc8_8192(a);
		fftc8_8192(b);
		fftc8_mul8192(a,b);
		fftc8_scale8192(a);
		fftc8_un8192(a);
		break;
	  }
	default:
	}

  /*  write result back... */
  for (i=0;i<nA;i++)  Aptr[i] = a[i].re;

  /*  clear mem */
  (*env)->ReleaseDoubleArrayElements(env,A,Aptr,0);
  (*env)->ReleaseDoubleArrayElements(env,B,Bptr,0);
  free(a);
  free(b);

  return (jboolean)1;
}
#endif
