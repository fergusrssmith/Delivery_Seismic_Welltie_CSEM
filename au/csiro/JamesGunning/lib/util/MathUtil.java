package au.csiro.JamesGunning.lib.util;

/**
Basic class of static public functions for elementary maths, plus
useful constants.

Class made final to assist compiler inlining
JNI methods added 2008/2009.

@author james.gunning@dpr.csiro.au
@version 2.0, 1/7/2008 
*/

import java.util.*;
import java.io.File;

public final class MathUtil
{
  /**
	 Location of native builds of this library wrt root of java packages.
  */
  public static final String relativeDir = "au/csiro/JamesGunning/lib/util";
  private static boolean haveNativeLibrary = false;

  public static native double expNative(double x);
  public static native double logNative(double x);
  public static native double sinNative(double x);
  public static native double cosNative(double x);
  public static native double tanNative(double x);

  


  /*
  static 
  { 
	String path;	
	path = getFullLibraryPath();

	try
	  {
		//System.err.println("PATH: "+path);
		System.load(path);
		haveNativeLibrary=true;
	  }
	catch (Exception e) 
	  {
		haveNativeLibrary=false;
		System.err.println("e: "+e);
	  }
	 System.err.println("Using native MathUtil = "+haveNativeLibrary);
  }
  
  */  

  /*
	Method to return full native library path of the correct
	platform-specific libMathUtil.* library.  On 64bit Linux platforms
	("Linux-x86_64-unknown") running 32bit jvms, the path will be the
	32bit build in "Linux-i686-unknown"
  */
  public static String getFullLibraryPath()
  {
	String path;

	// Resolve absolute name to this class
	try
	  {
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		path = loader.getResource(relativeDir).toURI().getPath();
		path = path + File.separatorChar + JNIPathUtil.getJNILibraryPath("libMathUtil",null);
	  }
	catch (Exception e) 
	  {
		System.err.println("ClassLoader-error: "+e);
		path="";
	  }
	return path;
  }

  /** 
	Sample tester for fast ceil and floor methods
  */
  public static void main(String argv[])
  {
	int i,j,k,M=0,N=0,tt;
	double x,y,z,MAX=100.0;
	boolean doInt;
	long tStart,tE,tJ;
	double tmin,dt,dt0;

	System.out.println("Have native = "+haveNativeLibrary);

	N=10000000;
	x=2.05;
	System.out.println("Correctness tests:");
	System.out.printf("Native log(2.0) test: %g, native %g%n",Math.log(x),log(x));
	System.out.printf("Native exp(2.0) test: %g, native %g%n",Math.exp(x),exp(x));

	System.out.println("Speed tests:");

	for (k=0;k<1;k++)
	  {
		y=10.0/N;
		tStart = System.currentTimeMillis();
		for (i=0,x=0.2;i<N;i++)
		  {
			x+=y;
			z = Math.sin(x);
		  }
		dt0 = (System.currentTimeMillis()-tStart)/(double)N;
		System.out.println("Mean call-time java.lang.Math.sin = "+dt0);
		tStart = System.currentTimeMillis();
		for (i=0,x=0.2;i<N;i++)
		  {
			x+=y;
			z = sin(x);
		  }
		dt = (System.currentTimeMillis()-tStart)/(double)N;
		System.out.println("Mean call-time MathUtil.sin = "+dt+"  java/JNI ratio = "+dt0/dt);

		y = 100.0/N;
		tStart = System.currentTimeMillis();
		for (i=0,x=2.0;i<N;i++)
		  {
			//x+=y;
			z = Math.exp(x);
		  }
		dt0 = (System.currentTimeMillis()-tStart)/(double)N;
		System.out.println("Mean call-time java.lang.Math.exp = "+dt0);
		for (i=0,x=2.0;i<N;i++)
		  {
			//x+=y;
			z = exp(x);
		  }
		dt = (System.currentTimeMillis()-tStart)/(double)N;
		System.out.println("Mean call-time MathUtil.exp = "+dt+"  java/JNI ratio = "+dt0/dt);
	  }
	System.out.println("Ceil, floor tests");
	for (k=0;k<2;k++)
	  {
		for (i=0;i<N;i++)
		  {
			x = 1e7*(Math.random()-0.5);
			doInt = (Math.random() < 0.2);
			if (doInt) x = (double)(int)x;
			if (fFloor(x) != Math.floor(x)) System.out.println("FAILED on "+x);
			if (fCeil(x) != Math.ceil(x)) System.out.println("FAILED on "+x);
		  }

		tStart = System.currentTimeMillis();
		x=y=0;
		for (i=0;i<N;i++)
		  {
			x+=0.777;
			y = Math.floor(x);
		  }
		System.out.println("Time in strict = "+(System.currentTimeMillis()-tStart));
		tStart = System.currentTimeMillis();
		x=y=0;
		for (i=0;i<N;i++)
		  {
			x+=0.777;
			y = fFloor(x);
		  }
		System.out.println("Time in fFloor = "+(System.currentTimeMillis()-tStart));
	  }
  }


  /**
	 0.5 Log(2 Pi)
  */
  public static final double HALF_LOG_TWOPI = 0.9189385332046727418;
  /**
	 Log(10)
  */
  public static final double LOG10  = 2.3025850929940456840;


  public final static double exp(double x)
  {
	return haveNativeLibrary ? expNative(x) : Math.exp(x);
  }

  public final static double log(double x)
  {
	return haveNativeLibrary ? logNative(x) : Math.log(x);
  }

  public final static double sin(double x)
  {
	return haveNativeLibrary ? sinNative(x) : Math.sin(x);
  }

  public final static double cos(double x)
  {
	return haveNativeLibrary ? cosNative(x) : Math.cos(x);
  }

  public final static double tan(double x)
  {
	return haveNativeLibrary ? tanNative(x) : Math.tan(x);
  }


  // return base^p. Must be called with p>=0
  public static final int ipow(int base,int p)
  {
	int r = 1;
	while (p!=0)
	  {
        if ((p & 1) != 0) r *= base;
        p >>= 1;
        base *= base;
	  }
	
    return r;
  }

  public static final int sqr(int x)
  {
	return x*x;
  }

  public static final float sqr(float x)
  {
	return x*x;
  }

  public static final double sqr(double x)
  {
	return x*x;
  }

  public static final int cube(int x)
  {
	return x*x*x;
  }

  public static final float cube(float x)
  {
	return x*x*x;
  }

  public static final double cube(double x)
  {
	return x*x*x;
  }

  public static final float pow4(float x)
  {
	float y=x*x;
	return y*y;
  }

  public static final double pow4(double x)
  {
	double y=x*x;
	return y*y;
  }

  public static final float fpow(float a,float x)
  {
	return (float) Math.pow(a,x);
  }

  public static final float fpow10(float x)
  {
	return (float) Math.pow(10.0,x);
  }

  public static final float fpow10(double x)
  {
	return (float) Math.pow(10.0,x);
  }

  public static final float flog10(float x)
  {
	return (float) Math.log10(x);
  }

  public static final float flog10(double x)
  {
	return (float) Math.log10(x);
  }

 public static final int factorial(int i)
  {
	int f=1;
	if (i<0) return 0;

	for (int j=i;j>=2;j--) f*=j;
	return f;
  }




  /** 
	  A faster (sloppy) version of StrictMath.floor(). valid over the rang of "longs"
  */
  public static final double fFloor(double y)
  {
	if (y>=0) return (double)(long)y;
	else 
	  {
		if (((double)(int)y) == y) return y;   // exact -ve integers
		else return (double)(-(long)(-y) - 1);
	  }
  }

 /** 
	  A faster (sloppy) version of StrictMath.floor(), valid over the
	  int range
  */
  public static final int fFloor_int(double y)
  {
	if (y>=0) return (int)y;
	else 
	  {
		if (((double)(int)y) == y) return (int)y;   // exact -ve integers
		else return (-(int)(-y) - 1);
	  }
  }

  /** 
	  A faster (sloppy) version of StrictMath.ceil()
	  returning int.
  */
  public static final int fCeil_int(double y)
  {
	return -fFloor_int(-y);
  }

  /** 
	  A faster (sloppy) version of StrictMath.ceil()
  */
  public static final double fCeil(double y)
  {
	return -fFloor(-y);
  }

  /**
	 This subroutine finds the combination set of n things<br>
	 taken p at a time for a given lexicographical index.<br>
	 n - number of things in the set<br>
	 p - number of things in each combination<br>
	 l - lexicographical index of combination sought<br>
	 output array containing the combination set<br>
	 the following relationships must exist among the input<br>
	 variables.  l must be greater than or equal to 1 and less<br>
	 than or equal to the maximum lexicographical index.<br>
	 p must be less than or equal to n and greater than zero.<br>
	 @param  n number of things in the set<br>
	 @param  p number of things in each combination<br>
	 @param  l lexicographical index of combination sought<br>
	 @returns subset array (unit offset)
  */
  public static final int[] getSubsets(int n, int p, int l)
  {
	
	int i,k, r, p1;
	int[] c;
	c = new int[p];
	//initialize lower bound index at zero
	k = 0;
	//loop to select elements in ascending order
	p1 = p - 1;
	//for (i=1;i<=p1;i++)
	for (i=0;i<p1;i++)
	  {
		//set lower bound element number for next element value
		c[i] = 0;
		if (i!=0) c[i] = c[i-1];
		//loop to check validity of each element value
		do
		  {
			c[i]++;
			r = binomial(n-c[i],p-(i+1));
			k += r;
		  }
		while (k<l);
		k = k - r;
	  }
	c[p-1] = c[p1-1] + l - k;
	return c;
  }
  
  /**
   Binomial function. ACM algorithm 160 translated to java. Calculates the
   number of combinations of m things taken n at a time.

   @param  m integer
   @param  n integer
   @returns  number of combinations
  */

  public static final int binomial(int m, int n)
  {
	int p, i, n1, r;
	n1 = n;
	p = m - n1;
	if (n1<p) 
	  {
		p = n1;
		n1 = m - p;
	  }
	r = n1 + 1;
	if (p==0) r = 1;
	if (p>=2)
	  {
		for (i=2;i<=p;i++)
		  {
			r = (r*(n1+i))/i;
		  }
	  }
	return r;
  }
}








