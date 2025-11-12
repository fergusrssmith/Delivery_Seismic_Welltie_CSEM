package au.csiro.JamesGunning.lib.util;

/**
Basic class of static public functions for blas/lapack wrappers

@author James.Gunning@csiro.au
@version 1.0, 1/7/2009 
*/

import java.util.*;
import java.io.File;

public final class BlasLapackUtil
{
  /**
	 Location of native builds of this library wrt root of java packages.
  */
  public static final String relativeDir = "au/csiro/JamesGunning/lib/util";
  private static boolean haveNativeLibrary = false;


  // dgemv: y = alpha*op(A)*x + beta*y
  // op = ^T if transA. 
  // A = double[M][N], unrolled as a N*M 1D array, row-major, so columns fastest
  // x = double[columns of op(A)]
  // y = double[rows of op(A)]
  public static native boolean dgemvNative(int M,int N,boolean transA,double alpha,
										   double[] A,double[] x,double beta,double[] y);
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
	catch (UnsatisfiedLinkError e)
	  {
		System.err.println("Using use native BlasLapackUtil library. Defaulting to java methods");
		haveNativeLibrary = false;
	  }
	catch (Exception e) 
	  {
		haveNativeLibrary=false;
		System.err.println("e: "+e);
	  }
	 System.err.println("Using native BlasLapackUtil = "+haveNativeLibrary);
  }
  

  /*
	Method to return full native library path of the correct
	platform-specific libBlasLapackUtil.* library.  On 64bit Linux platforms
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
		path = path + File.separatorChar + JNIPathUtil.getJNILibraryPath("libBlasLapackUtil",null);
	  }
	catch (Exception e) 
	  {
		System.err.println("ClassLoader-error: "+e);
		path="";
	  }
	return path;
  }

  /** 
	Sample tester for blas functions
  */
  public static void main(String argv[])
  {
	int i,j,k,M=0,N=0,tt;
	double x,y,z,MAX=100.0;
	boolean doInt;
	long tStart,tE=0,tJ;
	double tmin;

	System.out.println("Have native = "+haveNativeLibrary);

	System.out.println("DGEMV tests: (1) coding accuracy");

	// DGEMV array test
	int m=4,n=7;
	double alpha=.666,beta=0.1;

	double[] A = new double[m*n];
	for(i=0,k=0;i<m;i++)
	  {
		for(j=0;j<n;j++,k++) A[k]=k;
	  }
	for(tt=0;tt<=1;tt++)
	  {
		if (tt==0)
		  {
			N=n;
			M=m;
		  }
		else
		  {
			N=m;
			M=n;
		  }
		double[] X = new double[N];
		double[] Y = new double[M];
		double[] Y0 = new double[M];
		for(j=0;j<N;j++) X[j]=1-j; 
		for(j=0;j<M;j++) Y0[j]=j*j;

		System.out.println("JAVA:"+((tt==1)?" (A trans)":""));	
		for(j=0;j<M;j++) Y[j]=Y0[j];
		for(i=0;i<M;i++) 
		  {
			if (tt==0) for(j=0,x=0;j<N;j++) x += A[i*n+j]*X[j];
			else for(j=0,x=0;j<N;j++) x += A[j*n+i]*X[j];
			Y[i] = alpha*x + beta*Y[i];
			System.out.printf("%6g ",Y[i]);
		  }
		System.out.println("");

		System.out.println("LIBRARY"+((tt==1)?" (A trans)":""));
		for(j=0;j<M;j++) Y[j]=Y0[j];
		dgemv_kernel(m,n,(tt==0)?false:true,alpha,A,X,beta,Y);
		for(i=0;i<M;i++) System.out.printf("%6g ",Y[i]);
		System.out.println("");

		if(haveNativeLibrary) 
		  {
			for(j=0;j<M;j++) Y[j]=Y0[j];
			// y = alpha*op(A)*x + beta*y
			dgemvNative(m,n,(tt==0)?false:true,alpha,A,X,beta,Y);
			System.out.println("Native:"+((tt==1)?" (A trans)":""));	
			for(i=0;i<M;i++) System.out.printf("%6g ",Y[i]);
			System.out.println("");
		  }
	  }

	// Benchmarks
	System.out.printf("DGEMV benchmarks:%n");
	tmin=500; // min time per size, ms;
	int k2;
	for(m=100;m<2000;m+=200)
	  {
		n=m;
		A = new double[m*n];
		for(i=0,k=0;i<m;i++)
		  {
			for(j=0;j<n;j++,k++) A[k]=0.001*k;
		  }
		double[] X = new double[n];
		double[] Y = new double[m];
		for(j=0;j<n;j++) X[j]=1-j; 
		for(j=0;j<m;j++) Y[j]=0.2*j;
		x = 1e-3*((double)(n*m));
		if (haveNativeLibrary)
		  {
			tStart = System.currentTimeMillis();
			for(k=1;;k++)
			  {
				dgemvNative(m,n,false,alpha,A,X,beta,Y);
				tE = System.currentTimeMillis() -tStart;
				if (tE > tmin) break;
			  }
		  }
		
		for(j=0;j<n;j++) X[j]=1-j; 
		for(j=0;j<m;j++) Y[j]=0.2*j;
		tStart=System.currentTimeMillis();
		for(k2=1;;k2++)
		  {
			dgemv_kernel(m,n,false,alpha,A,X,beta,Y);
			tJ = System.currentTimeMillis() -tStart;
			if (tJ > tmin) break;
		  }
		String s = haveNativeLibrary ? (""+2*k*x/tE) : "not-available";
		System.out.printf("Mflop rate square dgemv, m=%d, java=%4g native=%s%n",m,2*k2*x/tJ,s);
	  }
  }


  /** 
	  DGEMV wrapper to BLAS
	  Matrix A = m-rows x n-columns is is 
	  passed in zero-offset, row-major (i.e. columns fastest)
	  (i.e. in C-style, "A[i][j]" = A[i*n + j])
	  x and y are presumed zero-offset, and big enough for the operation, i.e.
	  if (!transA)
	     y = double[m]
         x = double[n]
	  else
	     y = double[n]
         x = double[m]

	  Returned result is

	  @returns  (in situ) y = alpha * op(A) * x + beta * y
  */
  public static final void dgemv(int m,int n,boolean transA, double alpha,
					 double[] A,double[] x,double beta,double[] y)
  {
	if (haveNativeLibrary)
	  {
		dgemvNative(m,n,transA,alpha,A,x,beta,y);
	  }
	else
	  {
		dgemv_kernel(m,n,transA,alpha,A,x,beta,y);
	  }
  }

  private static final void dgemv_kernel(int m,int n,boolean transA, double alpha,
										 double[] A,double[] x,double beta,double[] y)
  {
	// dopy literal version
	int i,ii,j,k;
	double sum;
	if (!transA)
	  {		
		for(i=0;i<m;i++) 
		  {
			k= i*n;
			for(j=0,sum=0;j<n;j++) sum += A[k+j] * x[j];
			y[i] = alpha*sum+ beta*y[i];
		  }
		
		/* the JVM optimises above better than unrolling manually.. as here
			// unrolled
			for(i=0;i<m;i++) 
			  {
				ii=i*n;
				for(k=n%4,sum=0,j=0;--k>=0;) 
				  {
					sum += A[ii++] * x[j++]; 
				  }
				for(k=n/4;--k>=0;) 
				  {
					sum += A[ii++] * x[j++];
					sum += A[ii++] * x[j++];
					sum += A[ii++] * x[j++];
					sum += A[ii++] * x[j++];
				  }
				y[i] = alpha*sum+ beta*y[i];
			  }
		*/
	  }
	else
	  {
		for(i=0;i<n;i++) 
		  {
			for(j=0,sum=0;j<m;j++) sum += A[j*n+i] * x[j];
			y[i] = alpha*sum + beta*y[i];
		  }
	  }
  }
}








