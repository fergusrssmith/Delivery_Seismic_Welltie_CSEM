/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 
// USE_NEW_COPYRIGHT
package au.csiro.JamesGunning.lib.Numerical.djbFFT;

import java.io.File;
import java.util.Arrays;

import au.csiro.JamesGunning.lib.util.*; 
import au.csiro.JamesGunning.lib.Numerical.FFT;

/**
   Class to provide JNI interface to the c-library djbfft fast Fourier routines 
   (http://cr.yp.to/djbfft.html). These
   seem to be the fastest FFTs around for small (2-8192) 1D systems, and are
   claimed to be better than FFTW routines. (Both are better than NR).
   <p>
   The convolve method will also default to the native FFT class convolution
   routines if it can't resolve the native dlbfft library.
*/


public final class DJBFFT 
{
  /**
	 Location of native builds of djbfft wrt root of java packages.
	 This location contains builds in the directories<p>

	 DJBRelativeDir/Linux-i686-unknown/lib/libDJBFFT.so,<p>
	 DJBRelativeDir/Linux-x86_64-unknown/lib/libDJBFFT.so,<p>
	 DJBRelativeDir/Darwin-Power-Macintosh-unknown/lib/libDJBFFT.jnilib<p>

	 etc etc
  */
  public static final String DJBRelativeDir = "au/csiro/JamesGunning/lib/Numerical/djbFFT/djbfft-0.76";

  private static boolean haveNativeLibrary = true;

  public static native boolean convolveNative(double[] a, int na, double[] b, int nb);

  static 
  {
	String path;
	
	path = getFullDJBLibraryPath();
	try
	  {
		// System.loadLibrary("DJBFFT"); // old style
		System.load(path);
	  }
	catch (SecurityException se) 
	  {
		haveNativeLibrary=false;
	  }
 	catch (UnsatisfiedLinkError le) 
	  {
		haveNativeLibrary=false;
		System.err.println("le: "+le);
	  }
  }

  /*
	Method to return full native library path of the correct
	platform-specific libDJBFFT.* library.  On 64bit Linux platforms
	("Linux-x86_64-unknown") running 32bit jvms, the path will be the
	32bit build in "Linux-i686-unknown"

	Typical full library names are:
	ROOT/au/csiro/JamesGunning/lib/Numerical/djbFFT/djbfft-0.76/Linux-x86_64-unknown/lib/libDJBFFT.so
  */
  public static String getFullDJBLibraryPath()
  {
	String djbPath=null;

	// Resolve absolute name to this class

	try
	  {
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		if (true)
		  {
			djbPath = loader.getResource(DJBRelativeDir).toURI().getPath();
			djbPath = djbPath + File.separatorChar + JNIPathUtil.getJNILibraryPath("libDJBFFT","lib");
		  }
		else
		  {
			// Experimental ideas here for loading from jar files. Suspect it can't
			// be done for jni libraries: probably has to be unpacked in some way 
			// like http://java.itags.org/java-core-apis/32891
			Class c= new Object() {}.getClass().getEnclosingClass();
			djbPath = "/"+DJBRelativeDir + File.separatorChar +
			  JNIPathUtil.getJNILibraryPath("libDJBFFT","lib");
			java.net.URL url = c.getResource(djbPath);
		  }
	  }
	catch (Exception e) 
	  {
		System.err.println("ClassLoader-error: "+e);
		djbPath="";
	  }

	return djbPath;
  }

  /**
	 Convolution of data[0...(n-1)] with respns[0...(n-1)]. n *must* be a power
	 of 2. Answer is placed in first 0...(n-1) elements of ans, but ans
	 must be of length 2n
  */
  public static final void convolve(double[] data, int n, double[] respns,double[] ans)
  {
	if (haveNativeLibrary && (n>=2) && (n<=8192))
	  {
		//System.err.println("Native");
		if (!convolveNative(data,n,respns,n))
		  {
			throw new RuntimeException("Native FFT routine failed");
		  }
		System.arraycopy(data,0,ans,0,n);
	  }
	else
	  {
		// Java fallbacks
		// 1) FFT class, zero offset convolution
		// FFT.convlvZO_clean(data,n,respns,1,ans);
		// 2) CWP convolve fallback (faster than FFT.convlvZO_clean): tests OK wrt FFT.convlvZO_clean
		FFT.convolve_CWP(data,n,respns,ans);
	  }
  }

  /**
	 Convolution of data[0...(n-1)] with respns[0...(n-1)]. Arbitrary n
	 Answer is placed in first 0...(n-1) elements of ans, but ans
	 must be of length 2n
  */
  public static final void convolve_gen(double[] data, int n, double[] respns,double[] ans)
  {
	if (haveNativeLibrary && (n>=2) && (n<=8192))
	  {
		double[] _data=data;
		double[] _respns=respns;
		int N;
		// round Up To Nearest Pwr Of 2...   
		N = FFT.getNearestPowerOfTwo(n);
		if (N!=n)
		  {
			_data = new double[N];
			System.arraycopy(data,0,_data,0,n);
			_respns = new double[N];
			int nc = (n+1)/2;
			System.arraycopy(respns,0,_respns,0,nc);
			System.arraycopy(respns,nc,_respns,N-nc,n-nc);
		  }
		if (!convolveNative(_data,N,_respns,N))
		  {
			throw new RuntimeException("Native FFT routine failed");
		  }
		System.arraycopy(_data,0,ans,0,n);
	  }
	else
	  {
		// Java fallbacks
		// 1) FFT class, zero offset convolution
		// FFT.convlvZO_clean(data,n,respns,1,ans);
		// 2) CWP convolve fallback (faster than FFT.convlvZO_clean): tests OK wrt FFT.convlvZO_clean
		FFT.convolve_CWP(data,n,respns,ans);
	  }
  }


  /** 
	means to interrogate library type
  */
  public static boolean haveNativeFFT()
  {
	return haveNativeLibrary;
  }

  /** 
	Switch off access to native FFT
  */
  public static void disableNativeFFT()
  {
	haveNativeLibrary = false;
  }

  /** 
	  Dummy testbed for method. Shows convolution of two impulses.
	  furnish also N lambda as commandline args (N=2^p) for variety, e.g.
	  java DJBFFT 64 20 | graphing-program
  */
  public static void main(String[] args) 
  {
	  double[] a;
	  double[] b;
	  double[] ans;

	  int i,j,N=128;
	  double lambda=20.0;
	  boolean tryCWP=false;
	  boolean testCorrectness=false, testTiming=false, help=false;

	  if (args.length>0) 
		{
		  for(i=0; i<args.length; i++)
			{
			  if (args[i].equals("-CWP")) {tryCWP=true; continue;}
			  if (args[i].equals("-t")) {testCorrectness=true; continue;}
			  if (args[i].equals("-f")) {testTiming=true; continue;}
			  if (args[i].equals("-h")) {help=true; continue;}
			  N=Integer.parseInt(args[i]);
			  if (args.length>2) lambda=Double.parseDouble(args[i+1]);
			}
		}
	  if (help || N<16) 
		{
		  System.err.println("usage: java DJBFFT [-CWP| -t | -f] [N lambda],    n >= 16 ");
		  System.err.println("-t tests accuracy cf CWP, FFTW methods");
		  System.err.println("-f tests speed cf CWP, FFTW methods");
		  System.exit(0);
		}
	  System.err.println("Attempting to load: "+getFullDJBLibraryPath());
	  System.err.println("Have native library = "+haveNativeFFT());
	  System.err.println("Using CWP version = "+tryCWP);

	  if(testCorrectness)
		{
		  DJBFFT.correctness_test();
		  System.exit(0);
		}
	  else if (testTiming)
		{
		  DJBFFT.timing_test(N);
		  System.exit(0);
		}
	  if ((N & -N) != N) 
		{
		  System.err.println("usage: java DJBFFT [-CWP | -t | -f] N lambda, N a power of 2");
		  System.exit(0);
		}

	  a = new double[2*N];
	  b = new double[2*N];
	  ans = new double[4*N];
	  
	  // example wavelet
	  double spacing = lambda/(2.0*Math.PI);
	  double mx=1,eps=1e-3, ls = spacing*spacing;
	  for (i=0;i<N;i++)
		{
		  double x=i*i;
		  b[i]=(1.0-x/ls)*Math.exp(-0.5*x/ls);
		  if(i==0) mx = Math.abs(b[i])*eps;
		  else if (Math.abs(b[i])<mx) b[i]=0.0;
		}
	  // wraparound bit (weaker precursor)
	  for (j=2*N-1,i=1;i<N;i++,j--)
		{
		  b[j]=b[i]*(6.0/(5.0+i));
		}
	  // impulses

	  a[N/4]=1.0;
	  a[N/2]=1.5;
	  a[(int)(0.05*N)]=-0.5;
	  a[(int)(0.95*N)]=0.8;

	  for (i=0;i<N;i++) System.out.println(a[i]); System.out.println("");
	  for (i=0;i<N;i++) System.out.println(b[i]); System.out.println("");
	  if(!tryCWP)
		DJBFFT.convolve(a,a.length,b,ans);
	  else
		FFT.convolve_CWP(a,a.length,b,ans);
	  for (i=0;i<N;i++) System.out.println(ans[i]); System.out.println("");
	  // thrashing test...
	  int Ntry=20000;
	  long t=System.currentTimeMillis();
	  for (i=0;i<Ntry;i++)
		{
		  Arrays.fill(a,0);  Arrays.fill(b,0);  Arrays.fill(ans,0);
		  // impulses
		  a[N/4]=1.0;
		  a[N/2]=1.5; 
		  if (!tryCWP) DJBFFT.convolve(a,a.length,b,ans);
		  else FFT.convolve_CWP(a,a.length,b,ans);
		}
	  System.err.println("Mean time per conv (us) = "+1e3*(System.currentTimeMillis()-t)/(double)Ntry);
  }

  static void correctness_test()
  {
	int i,j,k,N,nc,Ntest;
	double a[],b[],bl[],al[],ans[],ansd[],ansfd[];
	double mx,u,mw=0,mwf=0;
	
	Ntest=1000;
	for(k=0;k<Ntest;k++)
	  {
		N = (int)(4+1000*Math.random());
		// lift to pwr of 2
		//for(i=1;i<N;i<<=1){}
		//N=i;
		
		a = new double[N];
		b = new double[N];
		al = new double[2*N];  // first half has data
		bl = new double[2*N];  // must we wraparound form in size 2N
		ans = new double[4*N];
		ansd = new double[4*N];
		ansfd = new double[4*N];
		
		// "wavelet"
		for (i=0;i<N;i++) b[i]=Math.random();
		nc=(N+1)/2;
		System.arraycopy(b,0,bl,0,nc);
		System.arraycopy(b,nc,bl,N+nc,N-nc); // bl is length 2N version, wraparound form

		// data
		for (i=0;i<N;i++) a[i] = al[i] = Math.random(); // al is length 2N version, second half zeroed

		FFT.convolve_CWP(a,N,b,ans);
		for(i=0,mx=0;i<N;i++) mx = Math.max(mx,Math.abs(ans[i]));

		// FFTW
		au.csiro.JamesGunning.lib.Numerical.FFTW.FFTW.convolve(a,N,b,ansfd);
		for(i=0,u=0;i<N;i++) u = Math.max(u,Math.abs(ans[i]-ansfd[i]));
		mwf = Math.max(mwf,u/mx);
		
		// DJBFFT
		DJBFFT.convolve_gen(al,2*N,bl,ansd);
		for(i=0,u=0;i<N;i++) u = Math.max(u,Math.abs(ans[i]-ansd[i]));
		mw = Math.max(mw,u/mx);

		//for (i=0;i<N;i++) System.out.printf("%f %f%n",ans[i],ansd[i]);		
	  }
	System.out.printf("MAX fractional error DJBFFT: %18g,  FFTW: %18g%n",mw,mwf);		
  }

  static void timing_test(int N)
  {
	int i,j,k,Ntest,nc;
	double a[],b[],bl[],al[],ans[],ansd[];
	double mx,u,mw=0;

	a = new double[N];
	b = new double[N];
	al = new double[2*N];  // first half has data
	bl = new double[2*N];  // must we wraparound form in size 2N
	ans = new double[4*N];
	ansd = new double[4*N];
		
	// "wavelet"
	for (i=0;i<N;i++) b[i]=Math.random();
	// data
	for (i=0;i<N;i++) a[i] = Math.random(); // al is length 2N version, second half zeroed

	long t0=0,t1=0,t2=0,tp=0,tStart;

	Ntest=100000;
	au.csiro.JamesGunning.lib.Numerical.FFTW.FFTW.convolve(a,N,b,ansd); // warmup FFTW
	for(k=0;k<Ntest;k++)
	  {
		// data
		Arrays.fill(al,0);
		System.arraycopy(a,0,al,0,N);
		Arrays.fill(bl,0);
		nc=(N+1)/2;
		System.arraycopy(b,0,bl,0,nc);
		System.arraycopy(b,nc,bl,N+nc,N-nc); // N+nc = (2N-(N-nc))
		
		tStart = System.nanoTime();
		FFT.convolve_CWP(a,N,b,ans);
		t0 += ((tp = System.nanoTime()) - tStart); tStart = tp;

		// FFTW
		au.csiro.JamesGunning.lib.Numerical.FFTW.FFTW.convolve(a,N,b,ansd);
		t2 += ((tp = System.nanoTime()) - tStart);	tStart = tp;

		// DJBFFT 
		DJBFFT.convolve_gen(al,2*N,bl,ansd);
		t1 += System.nanoTime() - tStart;
	  }
	System.out.printf("Relative timings:%nCWP %d (1.0)  DJBFFT %d (%f)  FFTW %d (%f): %n",t0,t1,t1/(double)t0,t2,t2/(double)t0);
  }
}

