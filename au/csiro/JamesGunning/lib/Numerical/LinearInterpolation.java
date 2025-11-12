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


package au.csiro.JamesGunning.lib.Numerical;

import au.csiro.JamesGunning.lib.*;


/**
   Class with useful methods for linear data interpolation
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public final class LinearInterpolation
{
  private double xa[];
  private double ya[];
  private boolean flatExtrap = false;

  // local variables used in fast table location stuff
  private int last_indx=0;
  private int n;

  ErrorHandler EH;

  /**
	 LinearInterpolation constructor. Set up an interpolatorfor data x[0...n-1]
	 y[0...n-1]. Ensure data has appropriate range, and x[] has no
	 repetitions. No checks in place.  Subsequent calls to interpolate use
	 the dieWithStackTrace call of the ErrorHandler, so calls outside the data
	 range should be avoided.
  */
  public LinearInterpolation(double x[], double y[], int n, ErrorHandler EH)
  {
	xa = x;
	ya = y;
	this.n=n;
	this.EH = EH;
  }

  /**
	 LinearInterpolation constructor with controllable
	 extrapolation behaviour. Set up an interpolatorfor data x[0...n-1]
	 y[0...n-1]. Ensure data has appropriate range, and x[] has no
	 repetitions. No checks in place. 
	 
	 If flatExtrap == true, extrapolation returns tolerant endpoint, not fatal error
  */
  public LinearInterpolation(double x[], double y[], int n, ErrorHandler EH,boolean flatExtrap)
  {
	xa = x;
	ya = y;
	this.n=n;
	this.EH = EH;
	this.flatExtrap = flatExtrap;
  }
  

  /**
	 Test harness
  */
  public static void main(String argv[])
  {
	int i,j,N=3,n=10;
	double X,dX;
	double[] x = new double[N];
	double[] y = new double[N];

	// generate & print true data
	for(i=0;i<N;i++)
	  {
		y[i] = Math.random();
		x[i] = (i>0) ? x[i-1]+1+Math.random(): 0;
		System.out.println(x[i]+ " " + y[i]);
	  }
	LinearInterpolation sp = new LinearInterpolation(x,y,N,new ErrorHandler());

	// now compute & print some intermediate points
	System.out.println("");
	dX=0.99999*(x[N-1]-x[0])/(n-1);

	for(i=0,X=x[0];i<n;i++,X+=dX)
	  {
		System.out.println(X+" " + sp.interpolate(X,true));
	  }

	// timing runs
	n=100000;
	dX=0.99999*(x[N-1]-x[0])/(n-1);
	long t = System.currentTimeMillis();

	for(i=0,X=x[0];i<n;i++,X+=dX)
	  {
		sp.interpolate(X,false);
	  }
	System.err.println("Slow Bisection time: "+(System.currentTimeMillis()-t));
	t = System.currentTimeMillis();
	for(i=0,X=x[0];i<n;i++,X+=dX)
	  {
		sp.interpolate(X,true);
	  }
	System.err.println("Fast bisection time: "+(System.currentTimeMillis()-t));
 }

 /**
	 return the spline--interpolated value y associated 
	 with the data xa[], ya[], whose
	 curvature is y2a[], at coordinate x.
	 @param x Value at which interpolated y is sought.
	 @param useLastValueInSearch Set to true if x is likely to be close to the last x used when this
	 method was called (performs more optimal bracketing search)
  */
  public double interpolate(double x, boolean useLastValueInSearch)
  {
	int klo,khi,k;
	double h,y,wt;

	if (x<xa[0] || x>xa[n-1]) 
	  {
		if (flatExtrap) return (x<xa[0]) ? ya[0] : ya[n-1];		  
		showState();
		throw new RuntimeException("Bad x input (x="+x+") to routine interpolate");
	  }

	if (useLastValueInSearch)
	  {
		// fast table searching method starting from last call
		klo = hunt(x);
		khi = (klo<n-1) ? klo+1 : klo;
	  }
	else
	  {
		// bisection to find last value
		klo = 0;
		khi = n-1;
		while (khi-klo > 1) {
		  k=(khi+klo) >> 1;
		  if (xa[k] > x) khi=k;
		  else klo=k;
		}
	  }
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
	  {
		throw new RuntimeException("Bad x input (x="+x+") to routine interpolate");
	  }

	// Ok now interpolate
	wt = (x-xa[klo])/h;
	y = wt*ya[khi]+(1-wt)*ya[klo];

	last_indx = klo;

	return y;
  }

  /**
	 Debugging
  */
  private void showState()
  {
	EH.println("Current data range");
	EH.println("X = "+xa[0]+"..."+xa[n-1]+", Y = "+ya[0]+"..."+ya[n-1]);
  }

  /** 
	  fast hunting method for searching table
  */
  private int hunt(double x)
  {
	int jlo=0,jm,jhi,inc;
	boolean ascnd;

	ascnd = (xa[n-1] >= xa[0]);
	if (jlo < 0 || jlo > (n-1)) 
	  {
		jlo=-1;
		jhi=n;
	  } 
	else 
	  {
		inc=1;
		if ((x >= xa[jlo]) == ascnd) 
		  {
			if (jlo == n-1) return jlo;
			jhi=jlo+1;
			while ((x >= xa[jhi]) == ascnd) 
			  {
				jlo=jhi;
				inc += inc;
				jhi=jlo+inc;
				if (jhi > (n-1)) 
				  {
					jhi=n;
					break;
				  }
			  }
		  } 
		else 
		  {
			if (jlo == 0) 
			  {
				jlo=-1;
				return jlo;
			  }
			jhi=jlo--;
			while ((x < xa[jlo]) == ascnd) 
			  {
				jhi=jlo;
				inc <<= 1;
				if (inc >= jhi) 
				  {
					jlo=-1;
					break;
				  }
				else jlo=jhi-inc;
			  }
		  }
	  }
	while (jhi-jlo != 1) 
	  {
		jm=(jhi+jlo) >> 1;
		if ((x >= xa[jm]) == ascnd)
		  jlo=jm;
		else
		  jhi=jm;
	  }
	if (x == xa[n-1]) jlo=n-2;
	if (x == xa[0]) jlo=0;
	return jlo;
  }
}

