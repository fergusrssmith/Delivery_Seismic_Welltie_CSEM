#include <stdio.h>
#include <fftc8.h>
#include <math.h>
/*
  Test of complex convolution interface for djbFFT library. To compile, e.g.
  cc -g -I./djbfft-0.76/include cmplxConvolutionsTest.c ./djbfft-0.76/Linux-i686-unknown/lib/djbfft.a -lm -o cmplxConvolutionsTest

  Pipe output to graphing program, e.g.
  cmplxConvolutionsTest|gr

  Performs a real convolution by using cmplx no's with zero im. parts...
*/
int main ()
{
  int i,j,m,spike,spike2,bspike=0;
  int nA=128;
  int nB=nA;
  double t;
  complex8 a[nA];
  complex8 b[nA];
  complex8 z4[nA];
  for(i=0;i<nA;i++) b[i].re=b[i].im=a[i].re=a[i].im=z4[i].re=z4[i].im=0.0;
  for(i=-nA/2;i<nA/2;i++)
	{
	  t=(i<=0)? exp(-i*i/100.0) : exp(-i*i/20.0);
	  b[((i+bspike)+2*nA)%nA].re=t;
	}
  for(i=0;i<nA;i++) printf("%g\n",b[i]); printf("\n");
  spike=15;
  spike2=56;
  a[spike].re=1; 
 a[spike2].re=1.6; 
  //a[57]=1; 
  // a[75]=1.5;
  for(i=0;i<nA;i++) printf("%g\n",b[((i-spike)>=0)?(i-spike):(nA+i-spike)].re+1.6*b[((i-spike2)>=0)?(i-spike2):(nA+i-spike2)].re); 
  printf("\n");
  // convolution....
  switch (nA)

	{
	case 128:
	  {  
		fftc8_128(a);
		fftc8_128(b);
		fftc8_mul128(a,b);
		fftc8_scale128(a);
		fftc8_un128(a);
	  }
	default:
	}
  for(i=0;i<nA;i++) printf("%g\n",a[i].re+0.001); printf("\n");
  
  return 0;
}
