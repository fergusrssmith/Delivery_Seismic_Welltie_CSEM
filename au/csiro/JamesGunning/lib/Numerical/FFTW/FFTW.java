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
package au.csiro.JamesGunning.lib.Numerical.FFTW;

import java.io.File;
import java.util.Arrays;
import static java.lang.Math.*;

import au.csiro.JamesGunning.lib.util.*; 
import au.csiro.JamesGunning.lib.Numerical.FFT;


/**
   Class to provide JNI interface to the c-library FFTW3 
   <p>
   The convolve method will also default to the native FFT class convolution
   routines if it can't resolve the native fftw library.
*/

public final class FFTW 
{
  /**
	 Location of native builds of fftw wrt root of java packages.
	 This location contains builds in the directories<p>

	 FFTWRelativeDir/Linux-i686-unknown/lib/libFFTW.so,<p>
	 FFTWRelativeDir/Linux-x86_64-unknown/lib/libFFTW.so,<p>
	 FFTWRelativeDir/Darwin-Power-Macintosh-unknown/lib/libFFTW.jnilib<p>

	 etc etc
  */
  public static final String FFTWRelativeDir = "au/csiro/JamesGunning/lib/Numerical/FFTW";
 
  private static boolean haveNativeLibrary = true;

  public static native boolean convolveNative(double[] a, int na, double[] b, int nB, int bc,double[] c, int nC);
  public static native void fftwCleanup();

  static 
  {
	String path;
	
	path = getFullFFTWLibraryPath();
	try
	  {
		// System.loadLibrary("FFTW"); // old style
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
	platform-specific libFFTW.* library.  On 64bit Linux platforms
	("Linux-x86_64-unknown") running 32bit jvms, the path will be the
	32bit build in "Linux-i686-unknown"

	Typical full library names are:
	ROOT/au/csiro/JamesGunning/lib/Numerical/FFTW/Linux-x86_64-unknown/libFFTW.so
  */
  public static String getFullFFTWLibraryPath()
  {
	String fftwPath=null;

	// Resolve absolute name to this class

	try
	  {
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		
		fftwPath = loader.getResource(FFTWRelativeDir).toURI().getPath();
		fftwPath = fftwPath + File.separatorChar + JNIPathUtil.getJNILibraryPath("libFFTW",null);		 
	  }
	catch (Exception e) 
	  {
		System.err.println("ClassLoader-error: "+e);
		fftwPath="";
	  }
	return fftwPath;
  }

  /**
	 Convolution of data[0...(n-1)] with respns[0...(n-1)]. 
	 Answer is placed in first 0...(n-1) elements of ans, but ans
	 must be of length 2n. +ve part of respns presumed to be elements [0...(n+1)/2-1]
	 (Note implied integer round-down). Works for odd n.
  */
  public static final void convolve(double[] data, int n, double[] respns,double[] ans)
  {
	if (haveNativeLibrary)
	  {
		//System.err.println("Native");
		if (!convolveNative(data,n,respns,n,(n+1)/2,ans,2*n))
		  {
			throw new RuntimeException("Native FFT routine failed");
		  }
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
	 Convolution of data[0...(n-1)] with respns[0...(nr-1)]. 
	 Respns has +ve (coda) elements respns[0...nrc] (i.e. nrc+1 of these)
	 Answer is placed in first 0...(n-1) elements of ans.
  */
  public static final void convolve(double[] data, int nd, double[] respns,int nr,int nrc,double[] ans, int na)
  {
	if (haveNativeLibrary)
	  {
		//System.err.println("Native");
		if (!convolveNative(data,nd,respns,nr,nrc+1,ans,na))
		  {
			throw new RuntimeException("Native FFT routine failed");
		  }
	  }
	else
	  {
		// Java fallbacks
		/**
		   -kx = sample index of start of precursor in filt[] in indexing: nr = kx+lx+1 = kx + nrc+1
		   lx = sample index of last nonzero sample in filt[] coda
		   filter is in wrap-around form. nr=kx+nrc;
		*/		
		FFT.convolve_CWP(data,nd,respns,nr-nrc-1,nrc,ans,na);
	  }
  }

  /**
   Cleanup of fftw plans etc
  */
  public static final void cleanupFFTW()
  {
	if (haveNativeLibrary) fftwCleanup();
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
	  Dummy testbed for method. Shows convolution of some impulses.
	  furnish also N lambda as commandline args (N=2^p) for variety, e.g.
	  java FFTW 64 20 | graphing-program
  */
  public static void main(String[] args) 
  {
	  double[] a;
	  double[] b;
	  double[] ans;

	  int i,j,Nend,bc=0,N=128;
	  double x,y,ls,spacing,mx,lambda=20.0;
	  boolean tryCWP=false;
	  boolean compactFilter = false;

	  if (args.length>0) 
		{
		  for(i=0; i<args.length; i++)
			{
			  if (args[i].equals("-CWP")) {tryCWP=true; continue;}
			  if (args[i].equals("-c")) {compactFilter=true; continue;}
			  if (args[i].equals("-h")) 
				{
				  System.err.println("usage: java FFTW [-CWP] [-c] [N lambda], N a power of 2");
				  System.exit(0);
				}
			  N=Integer.parseInt(args[i]);
			  lambda=Double.parseDouble(args[i+1]);
			  i++;
			  if ((N & -N) != N) 
				{
				  System.err.println("usage: java FFTW [-CWP] [-c] [N lambda], N a power of 2");
				  System.exit(0);
				}
			}
		}
	  if (N<16) 
		{
		  System.err.println("usage: java FFTW [-CWP] [-c] [N lambda],    n >= 16 ");
		  System.exit(0);
		}
	  System.err.println("Attempting to load: "+getFullFFTWLibraryPath());
	  System.err.println("Have native library = "+haveNativeFFT());
	  System.err.println("Using CWP version = "+tryCWP);
	  System.err.println("Compact filter= "+compactFilter);

	  a = new double[N];
	  ans = new double[N];

	  // get bc
	  // example wavelet
	  double eps=1e-3;	
	  spacing = lambda/(2.0*Math.PI);
	  ls = spacing*spacing;
	  for (i=N/2-1;i>=0;i--)
		{
		  x=i*i;
		  y=(1.0-x/ls)*exp(-0.5*x/ls);
		  if (abs(y)>eps) {bc=i; break;}
		}
	  int Nb = 2*bc+1;
	  if (compactFilter)
		{
		  b = new double[Nb];
		  Nend = Nb;
		}
	  else
		{
		  b = new double[N];
	  	  Nend=N;
		}
	  for (i=0,mx=1; i<=bc; i++)
		{
		  x=i*i;
		  b[i]=(1.0-x/ls)*Math.exp(-0.5*x/ls);
		  if(i==0) mx = Math.abs(b[i])*eps;
		  else if (Math.abs(b[i])<mx) b[i]=0.0;
		}
	  // wraparound bit (weaker precursor)
	  
	  // asymmetry...
	  for (j=Nend-1,i=1;i<=bc;i++,j--)
		{
		  b[j]=b[i]*(6.0/(5.0+i));
		}

	  // impulses

	  a[N/4]=1.0;
	  a[N/2]=1.5;
	  a[(int)(0.05*N)]=-0.5;
	  a[(int)(0.95*N)]=0.8;

	  // show data
	  for (i=0;i<N;i++) System.out.println(a[i]); System.out.println("");
	  // show wavelet
	  if (compactFilter)
		{
		  for (i=0;i<=bc;i++) System.out.printf("%f%n",b[i]);
		  for (j=0;j<(N-Nb);j++) System.out.printf("0%n");
		  for (;i<Nb;i++) System.out.printf("%f%n",b[i]);
		}
	  else
		{
		  for (i=0;i<N;i++) System.out.println(b[i]); 
		}
	  System.out.printf("%n");

	  if(!tryCWP)
		{
		  //FFTW.convolve(a,N,b,ans);
		  FFTW.convolve(a,N,b,Nend,bc,ans,N);		
		}
	  else
		{
		  if (compactFilter)
			{
			  // nr-nrc,nrc-1
			  FFT.convolve_CWP(a,N,b, Nb-bc,bc-1, ans,N);			  
			}
		  else
			{
			  FFT.convolve_CWP(a,N,b,ans);
			}
		}
	  for (i=0;i<N;i++) System.out.println(ans[i]); 
	  System.out.println("");

	  
	  // thrashing test...
	  int Ntry=200000;
	  long t=System.currentTimeMillis();
	  for (i=0;i<Ntry;i++)
		{
		  Arrays.fill(a,0);  Arrays.fill(b,0);  Arrays.fill(ans,0);
		  // impulses
		  a[N/4]=1.0;
		  a[N/2]=1.5; 
		  if (!tryCWP) 
			{
			  FFTW.convolve(a,N,b,Nend,bc,ans,N);		
			}
		  else
			{
			  if (compactFilter)
				{
				  // nr-nrc,nrc-1
				  FFT.convolve_CWP(a,N,b,Nb-bc,bc-1,ans,N);
				}
			  else
				{
				  FFT.convolve_CWP(a,N,b,ans);
				}
			} 
		}
	  System.err.println("Mean time per conv (us) = "+1e3*(System.currentTimeMillis()-t)/(double)Ntry);

	  // final fftw junk cleanup
	  FFTW.cleanupFFTW();
  }
}

