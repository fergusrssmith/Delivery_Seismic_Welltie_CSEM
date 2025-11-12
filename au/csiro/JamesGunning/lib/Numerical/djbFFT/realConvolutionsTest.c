#include <stdio.h>
#include <fftr4.h>
#include <math.h>
/*
  Test of real convolution interface for djbFFT library. To compile...
  cc -g -I./djbfft-0.76/include realConvolutionsTest.c ./djbfft-0.76/Linux-i686-unknown/lib/djbfft.a -lm -o realConvolutionsTest

  Pipe output to graphing program, e.g.
  realConvolutionsTest|gr

  Performs a real convolution using unusual ordering needed by interface
*/
#ifdef __APPLE__
#include   <Accelerate/Accelerate.h>
#endif

int main ()
{
  int i,j,m,spike,bspike=0;
  int inter=1;
  int nA=256;
  int nB=nA;
  double t;
  real4 a[nA];  
  real4 aa[nA];
  real4 b[nA];  
  real4 bb[nA];
  real4 z4[nA];
  for(i=0;i<nA;i++) b[i]=a[i]=z4[i]=0.0;
  /* b is Gaussian response */
  for(i=-nA/2;i<nA/2;i++)
	{
	  t=((i>0)?1.0/i:1)*exp(-i*i/100.0);
	  b[((i+bspike)+2*nA)%nA]=t;
	}
  for(i=0;i<nA;i++) printf("%g\n",b[i]); printf("\n");
  spike=15;
  /* a[] is a set of spikes */
  a[spike]=1; 
  a[65]=1.5;  
  a[125]=-0.9;
  //a[75]=1.5;
  for(i=0;i<nA;i++) printf("%g\n",a[i]); printf("\n");
  // convolution....

#ifdef __APPLE__
  real4 a_sm[2*nA-1];
  real4 b_sm[nA];

  for(i=0;i<nA/2;i++) a_sm[i]=0;              // zero-pad filter precursor length
  for(i=0,j=nA/2-1;i<nA;i++,j++) a_sm[j]=a[i];  // append signal
  for(i=j;i<(2*nA-1);i++) a_sm[i]=0;          // zero-pad filter coda length

  for(i=0;i<nA/2; i++) b_sm[i]=b[i+nA/2];
  for(j=0;j<nA/2;i++,j++) b_sm[i]=b[j];

  vDSP_conv(a_sm,1, // signal
			&b_sm[nA-1],-1, // filter
			aa, 1, // answer
			nA, // len(aa)
			nA); // filter width
				 
  for(i=0;i<nA;i++) printf("%g\n",aa[i]); printf("\n");

  //exit(0);

#endif


  if (inter)
	{
	  //input reordering...
	  // this works..

	  for (i=0;i<nA;i++) {aa[i]=a[i];bb[i]=b[i];} // copy
		for (i=j=0;i<nA;i+=2,j++) {a[i]=aa[j];b[i]=bb[j];}
		for (i=1,j=nA/2;i<nA;i+=2,j++) {a[i]=aa[j];b[i]=bb[j];}

	}
  if (0)
	{
	  // explicit formulas from djbfft accuracy code 
	  m=nA/2;
	  for (i = 0;i < m; i++) {
		for (j = 0;j < m; j++) {
		 if (i + j < m) {
		   z4[2*i + 2*j] += a[2*i] * b[2*j] + a[2*i+1] * b[2*j+1];
		   z4[2*i + 2*j+1] += a[2*i] * b[2*j+1] +  a[2*i+1] * b[2*j];
		 }
		 else {
		   z4[2*i + 2*j - 2*m+1] += a[2*i] * b[2*j] + a[2*i+1] * b[2*j+1];
		   z4[2*i + 2*j - 2*m] += a[2*i] * b[2*j+1] + a[2*i+1] * b[2*j];
		 }
		}
	  }
	  if (inter)
		{// undo interleaving ??
		  for (i=j=0;i<nA;i+=2,j++) {aa[j]=z4[i];}
		  for (i=1,j=nA/2;i<nA;i+=2,j++) {aa[j]=z4[i];}
		} else {  for(i=0;i<nA;i++) aa[i]=z4[i];}
	   for(i=0;i<nA;i++) printf("%g\n",aa[i]+0.1); printf("\n");
	}
  else
	{
	  
	  switch (nA)
		{
		case 256:
		  {
			fftr4_256(b);
			fftr4_256(a); 
			fftr4_mul256(a,b);
			fftr4_scale256(a);
			fftr4_un256(a);
		  }
		default:
		  {}
		}
	  if (inter)
		{// undo interleaving ??
		  for (i=j=0;i<nA;i+=2,j++) {aa[j]=a[i];}
		  for (i=1,j=nA/2;i<nA;i+=2,j++) {aa[j]=a[i];}
		} else {  for(i=0;i<nA;i++) aa[i]=a[i];}
	  for(i=0;i<nA;i++) printf("%g\n",aa[i]); printf("\n");
	


	}
  
  return 0;
}
