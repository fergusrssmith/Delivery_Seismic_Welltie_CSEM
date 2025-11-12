#include <jni.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>

/* FFTW stuff */
#include <fftw3.h>

// Hack to work around glibc memcpy saga. Sigh
// http://www.win.tue.nl/~aeb/linux/misc/gcc-semibug.html
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");


// statics for fftw plans
int _N,_Nbp,_Nold=-1;
fftw_plan _p=NULL,_q=NULL,_s=NULL;
double *_a_pad, *_b_pad;
fftw_complex *_cc;

void dfftw_convolve_1d(int n, double *a, double *b, double *c);
void dfftw_convolve_1d_tight(int n, double *a, double *b, int nb, int bc, double *c,int nc);
void fftwCleanup();

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* set to try to avoid data copies */
#define USE_GETPRIMITIVEARRAYCRITICAL 1

/* set to 1 to test timings */
#define TEST_TIMES 0

#if TEST_TIMES
static double prep_time=0, fft_time = 0;
static int n_calls=0;
#endif

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_lib_Numerical_FFTW_FFTW_fftwCleanup
(JNIEnv *env, jclass obj)
{
  fftwCleanup();
}

/*
  Convolve method using real convolutions from FFTW. Mainly good if sse happens
  data A[0...nA-1] convolved with B[0...nB-1] to land into C[0....nA-1]
  filter B has bc +ve values in entries B[0...bc-1]
  B[bc....nB-1] are -ve values
*/
JNIEXPORT jboolean JNICALL Java_au_csiro_JamesGunning_lib_Numerical_FFTW_FFTW_convolveNative
(JNIEnv *env, jclass obj, jdoubleArray A, jint nA, jdoubleArray B, jint nB, jint bc, jdoubleArray C, jint nC)
{
#if TEST_TIMES
  struct timeval tv;
struct timezone tz;
  double dt,tn;
#endif

extern int errno;
  int i,j;
  jboolean isCopy=JNI_FALSE;
  jdouble* a=NULL;
  jdouble* Bptr=NULL;
  jdouble* Cptr=NULL;

#if TEST_TIMES
  gettimeofday(&tv,&tz);
  dt = (double)((tv.tv_sec*1e6+tv.tv_usec));
#endif

  if (nC<nA)  return (jboolean)0;

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

#if USE_GETPRIMITIVEARRAYCRITICAL
  Cptr = (jdouble*) ((*env)->GetPrimitiveArrayCritical(env,C,&isCopy));
  if (Cptr == NULL) 
	{
	  printf("FAILED GetPrimitiveArrayCritical() for C in convolveNative()\n");
	  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0); // free a above
	  (*env)->ReleasePrimitiveArrayCritical(env,B,Bptr,0);
	  (*env)->ReleasePrimitiveArrayCritical(env,C,Cptr,0);
	  return (jboolean)0; /* exception occurred */
	}
#else
  Cptr = (*env)->GetDoubleArrayElements(env,C, NULL);
  if (Cptr == NULL) return (jboolean)0; /* exception occurred */
#endif


#if TEST_TIMES
  gettimeofday(&tv,&tz);
  prep_time += (tn = (double)((tv.tv_sec*1e6+tv.tv_usec))) - dt;
  dt = tn;
#endif

  /* FFTW convolution here */
  /* Good
	 dfftw_convolve_1d(nA, a, Bptr,Cptr);
  */
  /* better */
  dfftw_convolve_1d_tight(nA, a, Bptr, nB, bc, Cptr, nC);
	
#if TEST_TIMES
  gettimeofday(&tv,&tz);
  fft_time += (tn = (double)((tv.tv_sec*1e6+tv.tv_usec))) - dt;
  dt = tn;
#endif

 MEM_CLEAN:
  /* free mem */
#if USE_GETPRIMITIVEARRAYCRITICAL
  (*env)->ReleasePrimitiveArrayCritical(env,A,a,0);
  (*env)->ReleasePrimitiveArrayCritical(env,B,Bptr,0);
  (*env)->ReleasePrimitiveArrayCritical(env,C,Cptr,0);
#else
  (*env)->ReleaseDoubleArrayElements(env,A,a,0);
  (*env)->ReleaseDoubleArrayElements(env,B,Bptr,0);
  (*env)->ReleaseDoubleArrayElements(env,C,Cptr,0);
#endif


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

/* 
   See eg http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
   takes 256->256, 255->256
   Works only for 32 bit integer
*/
inline int roundUpToNearestPwrOf2(unsigned int v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

void dfftw_convolve_1d(int n, double *a, double *b, double *c)
{
  int i,Nr,N,Npad;
  double scale, *a_pad, *b_pad, *rc;
  fftw_complex *cc, *ca, *cb;
  fftw_plan p,q,s;
 

  N = 2*n-1;  
  
  // lift N to nearest pwr of 2
  N = roundUpToNearestPwrOf2(N);
  
  Npad = 2*(N/2+1);
  
  /* Allocate arrays and set up FFTW plan. */
  if ((_Nold != N) || (p==NULL))
	{
	  // destroy/cleanup any old stuff

	  a_pad = fftw_malloc(sizeof(double) * Npad);
	  b_pad = fftw_malloc(sizeof(double) * Npad);
	  cc = fftw_malloc(sizeof(fftw_complex) * ((N/2)+1));

	  p = fftw_plan_dft_r2c_1d(N, a_pad, (fftw_complex *)a_pad, FFTW_MEASURE);
	  q = fftw_plan_dft_r2c_1d(N, b_pad, (fftw_complex *)b_pad, FFTW_MEASURE);
	  s = fftw_plan_dft_c2r_1d(N, cc, (double *)cc, FFTW_MEASURE);
	  //fprintf(stderr,"fftw plan: N=%d\n",N);
	  // set cache
	  _Nold =N;
	  _p=p; 
	  _q=q; 
	  _s=s;
	  _a_pad = a_pad;
	  _b_pad = b_pad;
	  _cc = cc;
	  //fprintf(stderr,"fftw plan: N=%d\n",N);
	}
  else
	{
      // reuse plan in static cache
	  p=_p; q=_q; s=_s;
	  a_pad = _a_pad;
	  b_pad = _b_pad;
	  cc = _cc;
	}


  /* Initialize padded input arrays. */
  ca = (fftw_complex *) a_pad;
  memcpy(a_pad, a, n*sizeof(double));
  memset(&a_pad[n], 0, (Npad-n)*sizeof(double));
 
  cb = (fftw_complex *) b_pad;

  // JG: assume coda in [0....n/2-1], precursor in [n/2....n-1]
  // maps to [0....n/2-1............Npad-n/2...Nb-1]
  int Nb_pos = n/2;
  memcpy(b_pad, b, Nb_pos*sizeof(double));  //coda
  memset(&b_pad[Nb_pos], 0, (N-n)*sizeof(double));  //centre
  memcpy(&b_pad[N-(n-Nb_pos)], &b[Nb_pos], (n-Nb_pos)*sizeof(double)); 

  /* Get DFTs of a, b */
  fftw_execute(p); 
  fftw_execute(q); 
 
  scale=1.0/N;
  /* Multiply spectra */
  rc = (double *) cc;
  for (i=0; i<(N/2+1); i++)
	{
	  double ra, ia, rb,ib;
	  ra = ca[i][0];
	  ia = ca[i][1];
	  rb = cb[i][0];
	  ib = cb[i][1];
     
	  cc[i][0] = ra*rb - ia*ib;
	  cc[i][1] = ra*ib + ia*rb;
	}

  /* Compute the convolution by inverse DFT. */
  fftw_execute(s); 

  /* Scale. */
  //for (i=0; i<N; i++)  c[i] = scale*rc[i];
  for (i=0; i<n; i++)  c[i] = scale*rc[i];
 	
  //fftwCleanup();
}

/*
  As above. b a filter with b[0...nb-1]
  values b[0...bc-1] are +ve times. Remainder is wraparound (precursor)
  presumes nb<n
*/
void dfftw_convolve_1d_tight(int n, double *a, double *b, int nb, int bc, double *c,int nc)
{
  int i,Nr,N,Npad;
  double scale, *a_pad, *b_pad, *rc;
  fftw_complex *cc, *ca, *cb;
  fftw_plan p,q,s;
 

  N = n+nb-1;  
  
  // lift N to nearest pwr of 2
  N = roundUpToNearestPwrOf2(N);
  
  Npad = 2*(N/2+1);
  
  /* Allocate arrays and set up FFTW plan. */
  if (_Nold != N)
	{
	  // destroy/cleanup any old stuff
	  //fftwCleanup();

	  a_pad = fftw_malloc(sizeof(double) * Npad);
	  b_pad = fftw_malloc(sizeof(double) * Npad);
	  cc = fftw_malloc(sizeof(fftw_complex) * ((N/2)+1));

	  p = fftw_plan_dft_r2c_1d(N, a_pad, (fftw_complex *)a_pad, FFTW_MEASURE);
	  q = fftw_plan_dft_r2c_1d(N, b_pad, (fftw_complex *)b_pad, FFTW_MEASURE);
	  s = fftw_plan_dft_c2r_1d(N, cc, (double *)cc, FFTW_MEASURE);
	  //fprintf(stderr,"fftw plan: N=%d\n",N);
	  //fprintf(stderr,"fftw nb bc: %d %d\n",nb,bc);
	  // set cache
	  _Nold =N;
	  _p=p; 
	  _q=q; 
	  _s=s;
	  _a_pad = a_pad;
	  _b_pad = b_pad;
	  _cc = cc;
	}
  else
	{
      // reuse plan in static cache
	  p=_p; q=_q; s=_s;
	  a_pad = _a_pad;
	  b_pad = _b_pad;
	  cc = _cc;
	}


  /* Initialize padded input arrays. */
  ca = (fftw_complex *) a_pad;
  memcpy(a_pad, a, n*sizeof(double));
  memset(&a_pad[n], 0, (Npad-n)*sizeof(double));
 


  // coda in [0....bc-1], precursor in [bc....nb-1]
  // maps to [0....bc-1............(N-(nb-bc))...N-1]
  memcpy(b_pad, b, bc*sizeof(double));  //coda
  memset(&b_pad[bc], 0, (N-nb)*sizeof(double));  //centre
  memcpy(&b_pad[N-(nb-bc)], &b[bc], (nb-bc)*sizeof(double)); 

  //for(i=0;i<N;i++)printf("Q %f\n",a_pad[i]);

  /* Get DFTs of a, b */
  fftw_execute(p); 
  fftw_execute(q); 
 
  scale=1.0/N;
  /* Multiply spectra */
  rc = (double *) cc;
  cb = (fftw_complex *) b_pad;
  for (i=0; i<(N/2+1); i++)
	{
	  double ra, ia, rb,ib;
	  ra = ca[i][0];
	  ia = ca[i][1];
	  rb = cb[i][0];
	  ib = cb[i][1];
     
	  cc[i][0] = ra*rb - ia*ib;
	  cc[i][1] = ra*ib + ia*rb;
	}

  /* Compute the convolution by inverse DFT. */
  fftw_execute(s); 

  /* Scale. */
  //for (i=0; i<N; i++)  c[i] = scale*rc[i];

  // must have n+nb-1 components correct, but don't overrun c[] size  
  nc = MIN(n+nb-1,nc);
  for (i=0; i<nc; i++)  c[i] = scale*rc[i];
}

void fftwCleanup()
{
  fftw_destroy_plan(_p);
  fftw_destroy_plan(_q);
  fftw_destroy_plan(_s);
  
  fftw_free(_a_pad); 
  fftw_free(_b_pad); 
  fftw_free(_cc);
  _a_pad = NULL;
  _b_pad = NULL;
  _cc = NULL;
}

