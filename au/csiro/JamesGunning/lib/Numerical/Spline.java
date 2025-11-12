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

import java.io.*;
import java.util.*;
import au.csiro.JamesGunning.lib.*;


/**
   Class with useful methods for data splining, regular
   cubic splines.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public final class Spline
{
  public double xa[];
  public double ya[];
  public double y2a[];
  // local variable used in fast table location stuff
  private int last_indx=0;
  public int n;

  ErrorHandler EH;
  /**
	 Set endpoint conditions to zero curvature (natural spline)
	 if yp1 or ypn exceed this threshold
	 in the constructor
  */
  public static double ENDPOINT_THRESHOLD_FLAG = 0.99e30;
  /**
	 Spline constructor. Set up a spline for data x[0...n-1]
	 y[0...n-1]. Ensure data has appropriate range, and x[] has no
	 repetitions. No checks in place.  Subsequent calls to splint use
	 the dieWithStackTrace call of the ErrorHandler, so calls outside the data
	 range should be avoided. If either<br>
	 yp1>ENDPOINT_THRESHOLD_FLAG<br>
	 or<br>
	 ypn>ENDPOINT_THRESHOLD_FLAG<br>
	 natural splines (zero curvature) are
	 used at that spline endpoint

	 @param x array of x values (increasing!)
	 @param x array of y values (matching length of x)
	 @param n length of data arrays x[0....n-1], y[0....n-1]
	 @param yp1 derivative at leftmost endpoint
	 @param ypn derivative at rightmost endpoint
  */
  public Spline(double x[], double y[], int n, double yp1, double ypn, ErrorHandler EH)
  {
	xa = x;
	ya = y;
	this.n=n;
	this.EH = EH;
	initialiseSpline(yp1,ypn);
  }

  /**
	 Test harness
  */
  public static void main(String argv[])
  {
	int i,j,N=3,n=10,na=argv.length;
	boolean time=false;
	String fin=null,fout=null;
	double X,dX,u,v;
	ErrorHandler EH = new ErrorHandler();
  
	for(i=0;i<na;i++)
	  {
		if(argv[i].equals("-t")) time=true;
		if(argv[i].equals("-i") && (i+1)<na) fin= argv[++i];
		if(argv[i].equals("-o") && (i+1)<na) fout= argv[++i];
		if(argv[i].equals("-h"))
		  {
			EH.printf("Usage Spline [-i xy.txt] [-o xnew.txt]%n");
			EH.printf("  where xy is  ascii 2 column, xnew is 1col; new values to be interpolated%n");
			EH.printf("  -t runs a demo test of self generated data with timing%n");
			return;
		  }
	  }
	if (fin!=null && fout!=null)
	  {
		BufferedReader in=null;
		String s;
		try 
		  {
			in = new BufferedReader(new FileReader(fin));
		  }
		catch (Exception e)
		  {
			EH.warn("Can't open file "+fin);
		  }
		ArrayList<Float> xi = new ArrayList<Float>();
		ArrayList<Float> yi = new ArrayList<Float>();
		try
		  {
			while ((s = in.readLine()) != null) 
			  { 
				if (s.startsWith("#")) continue;
				String[] st = s.trim().split("\\s+");
				xi.add(Float.parseFloat(st[0]));
				yi.add(Float.parseFloat(st[1]));
			  }
		  }
		catch (Exception e)
		  {
			EH.warn("Can't parse num num"); 
		  }	
		try
		  {
			in.close();
		  }
		catch (Exception e)
		  {
			EH.warn("Can't close file "+fin);
		  }
		double[] _x = new double[N=xi.size()];
		double[] _y = new double[N];
		for(i=0;i<N;i++) {_x[i]=xi.get(i); _y[i]=yi.get(i);}
		Spline sp = new Spline(_x,_y,N,Double.MAX_VALUE,Double.MAX_VALUE, new ErrorHandler());
		try 
		  {
			in = new BufferedReader(new FileReader(fout));
		  }
		catch (Exception e)
		  {
			EH.warn("Can't open file "+fout);
		  }
		try
		  {
			while ((s = in.readLine()) != null) 
			  { 
				if (s.startsWith("#")) continue;				
				String[] st = s.trim().split("\\s+");
				u = Double.parseDouble(st[0]);
				if      (u<=_x[0])   v=_y[0];
				else if (u>=_x[N-1]) v=_y[N-1];
				else v = sp.splint(u,false);
				EH.printf("%f%n",v);			  
			  }
		  }
		catch (Exception e)
		  {
			EH.warn("Can't parse num in "+fout+" "+e); 
		  }
			
		try
		  {
			in.close();
		  }
		catch (Exception e)
		  {
			EH.warn("Can't close file "+fin);
		  }
	  }
	if(time)
	  {		
		double[] x = new double[N];
		double[] y = new double[N];
		// generate & print true data
		for(i=0;i<N;i++)
		  {
			y[i] = Math.random();
			x[i] = (i>0) ? x[i-1]+1+Math.random(): 0;
			EH.println(x[i]+ " " + y[i]);
		  }
		Spline sp = new Spline(x,y,N,Double.MAX_VALUE,Double.MAX_VALUE, new ErrorHandler());

		// now compute & print some intermediate points
		EH.println("");
		dX=0.99999*(x[N-1]-x[0])/(n-1);

		for(i=0,X=x[0];i<n;i++,X+=dX)
		  {
			EH.println(X+" " + sp.splint(X,true));
		  }

		// timing runs

		n=100000;
		dX=0.99999*(x[N-1]-x[0])/(n-1);
		long t = System.currentTimeMillis();
		
		for(i=0,X=x[0];i<n;i++,X+=dX)
		  {
			sp.splint(X,false);
		  }
		EH.warn("Slow Bisection time: "+(System.currentTimeMillis()-t));
		t = System.currentTimeMillis();
		for(i=0,X=x[0];i<n;i++,X+=dX)
		  {
			sp.splint(X,true);
		  }
		EH.warn("Fast bisection time: "+(System.currentTimeMillis()-t));
	  }
  }



  /**
	 Sets up spline curvature terms.....
  */
  private void initialiseSpline(double yp1, double ypn)
  {
	int i,k;
	double p,qn,sig,un;
	double[] u;

	u = new double[n];
	y2a = new double[n];

	if (yp1 > ENDPOINT_THRESHOLD_FLAG)
		y2a[0] = u[0] = 0.0;
	else 
	  {
		y2a[0] = -0.5;
		u[0] = (3.0/(xa[1]-xa[0]))*((ya[1]-ya[0])/(xa[1]-xa[0])-yp1);
	  }
	for (i=1;i<n-1;i++) 
	  {
		sig = (xa[i]-xa[i-1])/(xa[i+1]-xa[i-1]);
		p = sig*y2a[i-1]+2.0;
		y2a[i] = (sig-1.0)/p;
		u[i] = (ya[i+1]-ya[i])/(xa[i+1]-xa[i]) - (ya[i]-ya[i-1])/(xa[i]-xa[i-1]);
		u[i] = (6.0*u[i]/(xa[i+1]-xa[i-1])-sig*u[i-1])/p;
	  }
	if (ypn > ENDPOINT_THRESHOLD_FLAG)
		qn = un = 0.0;
	else 
	  {
		qn = 0.5;
		un = (3.0/(xa[n-1]-xa[n-2]))*(ypn-(ya[n-1]-ya[n-2])/(xa[n-1]-xa[n-2]));
	  }
	y2a[n-1] = (un-qn*u[n-2])/(qn*y2a[n-2]+1.0);
	for (k=n-2;k>=0;k--) y2a[k] = y2a[k]*y2a[k+1]+u[k];
	u = null;
  }




  /**
	 return the spline--interpolated value y associated 
	 with the data xa[], ya[], whose
	 curvature is y2a[], at coordinate x.
	 @param x Value at which interpolated y is sought.
	 @param useLastValueInSearch Set to true if x is likely to be close to the last x used when this
	 method was called (performs more optimal bracketing search)
  */
  public double splint(double x, boolean useLastValueInSearch)
  {
	int klo,khi,k;
	double h,b,a,y;

	if (x<xa[0] || x>xa[n-1]) 
	  {
		throw new RuntimeException("Bad x input to routine splint");
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
		throw new RuntimeException("Bad x input to routine splint");
	  }
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

	last_indx = klo;

	return y;
  }

  /**
	 Return the value of x at the absolute maximum or minimum of the function.
	 O(n) runtime.
  */
  public double getXatMinOrMax(boolean restrictToSpecifiedSign, boolean seekMax)
  {
	int j,klo,khi,jMax=0;
	boolean maximise;
	double h,h2,h4,b,a,y,x,det,tmp,A;
	double max=0, min=0, xMin, xMax;

	if (!restrictToSpecifiedSign)
	  {
		// find largest |sample|
		for (j=0; j<n; j++) 
		  {
			if (Math.abs(ya[j]) > max)
			  {
				jMax=j;
				max = Math.abs(ya[j]);
			  }
		  }
	  }
	else
	  {
		if (seekMax)
		  {
			max = ya[0];
			for (j=0; j<n; j++) 
			  {
				if (ya[j] > max)
				  {
					jMax=j;
					max = ya[j];
				  }
			  }
		  }
		else
		  {
			min = ya[0];
			for (j=0; j<n; j++) 
			  {
				if (ya[j] < min)
				  {
					jMax=j;
					min = ya[j];
				  }
			  }
		  }
	  }

	// if it's internal, get max or min from the adjacent two intervals
	if (jMax<=(n-2) && jMax>=1)
	  {
		// Peak appears to be internal. Find its time using 
		// analytical methods based on spline

		maximise = ((ya[jMax] >= ya[jMax+1]) && (ya[jMax] >= ya[jMax-1]));

		// loop over two intervals on either side of max
		max = -Double.MAX_VALUE;
		min = Double.MAX_VALUE;
		xMin = xMax = xa[jMax]; // overwritten almost certainly

		for (klo=jMax-1;klo<=jMax;klo++)
		  {
			khi = klo+1;	
			if (klo<0 || khi>(n-1)) continue;

			h = xa[khi]-xa[klo];
			if (h == 0.0) 
			  {
				throw new RuntimeException("Bad x in getXatMinOrMax()");
			  }
			h2 = h*h;
			h4 = h2*h2;

			// in this interval, the cubic spline has a parabolic derivative. Find its zeros
			tmp = y2a[klo] - y2a[khi];
			det = 36.0 * h4 * y2a[khi] * y2a[khi] + 
			  12.0 * h2 * tmp *(-6.0*ya[klo] + 6.0*ya[khi] + h2*(y2a[klo] + 2.0*y2a[khi]));
			if (det>=0 && tmp != 0)
			  {
				det = Math.sqrt(det);
				// loop over two roots of quadratic: j is the +/- sign.
				for (j=-1;j<=1;j+=2)
				  {
					A = j*((-j)*6.0*h2*y2a[khi] + det)/(6*h2*tmp);
					if (A>=0 && A<=1) // root is within interval
					  {
						x = xa[khi] - h*A;
						a = (xa[khi]-x)/h;
						b = (x-xa[klo])/h;
						y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
						if (y>max)
						  {
							max = y; xMax=x;
						  }
						if (y<min)
						  {
							min = y; xMin = x;
						  }
					  }
				  }
			  }
		  }
		return maximise ? xMax : xMin;
	  }
	else
	  {
		// endpoint
		return xa[jMax];
	  }
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

