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
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.abs;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import java.lang.Math.*;
import au.csiro.JamesGunning.lib.util.*;


/**
   Basic class of static public function for performing FFT's using
   <i>exactly</i> the notation of <b>Numerical Recipes in C, Chs. 12 & 13
   (Press et al)</b>, henceforth denoted NR. The calls for four1(),
   realft(), twofft() and convlv() are reproduced here. They are
   overloaded so as to work with common double[] arrays, but also the
   DoubleMatrix1D arrays from the cern.colt.matrix class libraries.

   Tests so far indicate the FFT's run significantly faster on 
   native java primitive-type arrays, like double[], compared to the
   CERN DoubleMatrix1D arrays. This may not hold for the 
   case  of 2D FFTs, but these are not yet in the library.   If you
   don't want the CERN stuff, comment it out and remove the
   import cern.colt.matrix.* line.
   
   
   Note that NR's unit--offset convention is used throughout, so most data vectors use
   elements data[1..n], where n=data.size()-1 if data is of type
   DoubleMatrix1D. The element data[0] is generally unused.

   Zero-offset equivalents for data in the range data[0..(n-1)] are embodied
   in methods with the suffix ZO (zero offset).
   
   
   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/8/2001 
*/

public final class FFT
{
  /**
	 Simple test functions to convolve ascii stdin data with smoother via FFT..
	 You don't ever need  main(), its just a test...
	 A variety of different tests are encoded here, just replace
	 the commented sections accordingly.....
  */
  
  public static void main(String[] args)
  {
	
	/* from the stdin ascii stream of numbers, */
	/* do some FFT stuff with it... */
	
	BufferedReader stdin;
	String line=null,tok;		
	int j,jmax,i;
	int ndata,nn,ng;
	double y;
	

	if (args.length!=1 && args.length!=2)
	{
	  System.out.println("usage: stdin | java FFT [file] ndata(pwr of 2)");
	  return;
	}
	
	if (args.length==1)
	{
	  nn = ndata = Integer.parseInt(args[0]);
	  stdin  = new BufferedReader(new InputStreamReader(System.in));
	}
	else
	{
	  try
	  {
		stdin = new BufferedReader(new FileReader(args[0]));
	  }
	  catch  (IOException e) {System.out.println("couldn't open " + args[0]);return;}
	  nn = ndata = Integer.parseInt(args[1]);
	}
	ng=0;
	while ((nn >>= 1) > 0 ) ng++;
	if (ndata != (1<<ng))
	  {
		System.out.println("ndata must be a power of 2");
		return;
	  }


	DoubleMatrix1D data = DoubleFactory1D.dense.make(ndata+1);

	String[] firstcolumns = new String[ndata+1];
	i=1;
	try
	  {
		while (((line=stdin.readLine()) != null) && (i<=ndata))
		  {
			firstcolumns[i]=""; 
			StringTokenizer st = new StringTokenizer(line);
			jmax=st.countTokens()-1;
			for (j=1;j<=jmax;j++)
			  {
				tok=st.nextToken();
				firstcolumns[i] += " " + tok;
			  } 
			firstcolumns[i] += " ";
			tok=st.nextToken();
			data.setQuick(i,java.lang.Double.parseDouble(tok));
			i++;
		  }
	  }
	// file i/o exception
	catch (IOException e) {System.out.println("couldn't read stdin");}
	

	
	/* start 2ffttest */
//  DoubleMatrix1D data3 = DoubleFactory1D.dense.make(ndata+1);
//  	DoubleMatrix1D fft= DoubleFactory1D.dense.make(2*ndata+1);
//  	DoubleMatrix1D fft3 = DoubleFactory1D.dense.make(2*ndata+1);
//  	for (i=1;i<=(data.size()-1);i++)
//  	  {
//  		data3.setQuick(i,java.lang.Math.exp(-(i-1)*(i-1)/100.0));
//  	  }
//  	twofft(data,data3,fft,fft3,ndata);

//  	for (i=1;i<=data.size()-1;i++)
//  	  {
//  		System.out.println(data.getQuick(i) + " " + data3.getQuick(i) + " " + fft.getQuick(i) + " "
//  						   + fft3.getQuick(i));
//  	  }

/* start convlv test*/
	/*
	  DoubleMatrix1D wave = DoubleFactory1D.dense.make(ndata+1);
	  DoubleMatrix1D smoothed = DoubleFactory1D.dense.make(2*ndata+1);
	  double[] wave2 = new double[ndata+1];
	  double[] wave2ZO = new double[ndata];
	  double[] smoothed2 = new double[2*ndata+1];
	  double[] smoothed2ZO = new double[2*ndata];
	  double[] data2 = data.toArray();
	  double[] data2ZO = new double[ndata]; 
	  for (i=1;i<=ndata;i++)
		{
		  data2ZO[i-1]=data2[i];
		}
	  Util.print(data2,4);	
	  // make Gaussian filter
	  for (i=1;i<=(data.size()-1)/2;i++)
		{
		  wave.setQuick(i,java.lang.Math.exp(-(i-1)*(i-1)/5.0));
		  wave2[i]=java.lang.Math.exp(-(i-1)*(i-1)/5.0);
		  wave2ZO[i-1]=wave2[i];
		}
	  for (j=(data.size()-1),i=2;i<=(data.size()-1)/2;i++,j--)
		{
		  wave.setQuick(j,java.lang.Math.exp(-(i-1)*(i-1)/5.0));
		  wave2[j]=java.lang.Math.exp(-(i-1)*(i-1)/5.0);
		  wave2ZO[j-1]=wave2[j];
		}
	  //	convlv(data,data.size()-1,wave,data.size()-2,1,smoothed);
	  convlv(data2,ndata,wave2,ndata-1,1,smoothed2);
	  convlvZO(data2ZO,ndata,wave2ZO,ndata-1,1,smoothed2ZO);
	  
	  for (i=1;i<=data.size()-1;i++)
		{
		  //		System.out.println(data.getQuick(i) + " " + wave.getQuick(i) + " " + smoothed.getQuick(i));
		  System.out.println(data2[i] + " " + wave2[i] + " " + smoothed2[i] + " " +  smoothed2ZO[i-1] + " diff = " + (smoothed2[i] - smoothed2ZO[i-1]));
		}
	
	*/
	  //	Start of zero offset test....
	  /*
		ndata=32;
		data = DoubleFactory1D.dense.make(ndata+1);
		double[] data2 = new double[ndata+1]; 
		double[] data2ZO = new double[ndata]; 
		for (i=1;i<=ndata;i++)
		{
		data.setQuick(i,(double)i/12.5);
		data2[i]=data.getQuick(i);
		data2ZO[i-1]=data2[i];
		}
		realft(data2,ndata,1);
		realftZO(data2ZO,ndata,1);
		for (i=1;i<=ndata;i++)
		{
		System.out.println(data2[i] + " " + data2ZO[i-1] + " diff = " + (data2[i]-data2ZO[i-1]));
		}
		System.out.println("After");
	
		realft(data2,ndata,-1);
		realftZO(data2ZO,ndata,-1);

		for (i=1;i<=ndata;i++)
		{
		System.out.println(data2[i] + " " + data2ZO[i-1] + " diff = " + (data2[i]-data2ZO[i-1]));
		}

	  */

	
	  // speed test
	  ndata=2097152;
	  data = DoubleFactory1D.dense.make(ndata+1);
	  double[] data2 = new double[ndata+1]; 
	  double[] data2ZO = new double[ndata]; 
	  for (i=1;i<=ndata;i++)
		{
		  data.setQuick(i,(double)i/12.5);
		  data2[i]=data.getQuick(i);
		  data2ZO[i-1]=data2[i];
		}
	  System.out.println("start CERN");
	  long st,dt;
	  st = System.currentTimeMillis();
	  four1(data,ndata/2,-1);
	  dt = System.currentTimeMillis() - st;
	  System.out.println("new took " + (double)dt/1000.0);
	  System.out.println("start ZO");

	  st = System.currentTimeMillis();
	  four1ZO(data2ZO,ndata/2,-1);
	  dt = System.currentTimeMillis() - st;
	  System.out.println("ZO took " + (double)dt/1000.0);

	  dt = System.currentTimeMillis() - st;
	  st = System.currentTimeMillis();	
	  four1(data2,(data.size()-1)/2,-1);
	  dt = System.currentTimeMillis() - st;
	  System.out.println("old took " + (double)dt/1000.0);
	
  }

 
  /** 
		Replaces data[1..2*nn] by its DFT, if isign=1, otherwise
		replaces data[1..2*nn] by nn times its inverse DFT, if
		isign=-1. <br>
        data[1..2*nn] is a complex array of length nn.  nn <b>MUST</b> be a
        power of 2.<br>
		N.B. data.size() = 2*nn+1, since by sticking to the NR 
		unit-offset convention, we don't use element data[0] 
  */
  public static void four1(DoubleMatrix1D data,int nn,int isign)
  {
	
	int n,mmax,m,istep;
	int i,j;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;
	double swap;

	n=  nn << 1;
	j=1;
	for (i=1;i<n;i+=2) 
	  {
		if (j > i) 
		  {
			swap=data.getQuick(j);
			data.setQuick(j,data.getQuick(i));
			data.setQuick(i,swap);
			swap=data.getQuick(j+1);
			data.setQuick(j+1,data.getQuick(i+1));
			data.setQuick(i+1,swap);
		  }
		m=n >> 1;
		while (m >= 2 && j > m) 
		  {
			j -= m;
			m >>= 1;
		  }
		j += m;
	  }
	mmax=2;
	while (n > mmax) 
	  {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=Math.sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=java.lang.Math.sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) 
		  {
			for (i=m;i<=n;i+=istep) 
			  {
				j=i+mmax;
				tempr=wr*data.getQuick(j)-wi*data.getQuick(j+1);
				tempi=wr*data.getQuick(j+1)+wi*data.getQuick(j);
				data.setQuick(j,data.getQuick(i)-tempr);
				data.setQuick(j+1,data.getQuick(i+1) - tempi);
				data.setQuick(i , data.getQuick(i) + tempr);
				data.setQuick(i+1,data.getQuick(i+1) + tempi);
			  }
			wtemp=wr;
			wr=wtemp*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		  }
		mmax=istep;
	  }
  }
   
  /** 
		Given data1[1..nn] and data2[1..n], this calls four1 and
		returns two complex output arrays fft1[1..2n] and fft2[1..2n],
		which contain the DFT of the respective data arrays.<br> 
		n <b>MUST</b> be a power of 2.<br> 
		N.B. data1.size() =
		data2.size() = n+1 and fft1.size() = fft2.size() = 2n+1, since
		by sticking to the NR unit-offset convention, we don't use
		element data[0] 
  */
  public static void twofft(DoubleMatrix1D data1, DoubleMatrix1D data2, DoubleMatrix1D fft1, 
							DoubleMatrix1D fft2, int n)
  {
    int  nn3,nn2,jj,j;
	double rep,rem,aip,aim;
	
	if (data2.size()!=data1.size())
	  {
		throw new IllegalArgumentException("data1 must be of same size as dat2 in twofft()");
	  }
	if (fft2.size()!=fft1.size())
	  {
		throw new IllegalArgumentException("fft1 must be of same size as fft2 in twofft()");
	  }
	if (2*(data1.size()-1) != (fft1.size()-1))
	  {
		throw new IllegalArgumentException("data1 must be of half size of fft1 in twofft()");
	  }

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
	  fft1.setQuick(jj-1,data1.getQuick(j));
	  fft1.setQuick(jj,data2.getQuick(j));
	}
	four1(fft1,n,1);
	fft2.setQuick(1,fft1.getQuick(2));
	fft1.setQuick(2,0.0);
	fft2.setQuick(2,0.0);
	for (j=3;j<=n+1;j+=2) {
	  rep=0.5*(fft1.getQuick(j)+fft1.getQuick(nn2-j));
	  rem=0.5*(fft1.getQuick(j)-fft1.getQuick(nn2-j));
	  aip=0.5*(fft1.getQuick(j+1)+fft1.getQuick(nn3-j));
	  aim=0.5*(fft1.getQuick(j+1)-fft1.getQuick(nn3-j));
	  fft1.setQuick(j,rep);
	  fft1.setQuick(j+1,aim);
	  fft1.setQuick(nn2-j,rep);
	  fft1.setQuick(nn3-j, -aim);
	  fft2.setQuick(j,aip);
	  fft2.setQuick(j+1, -rem);
	  fft2.setQuick(nn2-j,aip);
	  fft2.setQuick(nn3-j,rem);
	}
  }

  /** 
		Computes the DFT of a set of real data in data1[1..n], going forwards or
		backwards according to isign=1 or -1 respectively. the inverse transform must be
		multiplied by 2/n.<br>
		n <b>MUST</b> be a power of 2.<br>
		N.B. data.size() = n+1  and we don't use data[0].<br>
		See NR section 12.3 for all the gory detail.
  */

 public static void realft(DoubleMatrix1D data, int n, int isign)
  {
	int i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) 
	  {
		c2 = -0.5;
		four1(data,n>>1,1);
	  } 
	else 
	  {
		c2=0.5;
		theta = -theta;
	  }
	wtemp=Math.sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=Math.sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) 
	  {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data.getQuick(i1)+data.getQuick(i3));
		h1i=c1*(data.getQuick(i2)-data.getQuick(i4));
		h2r = -c2*(data.getQuick(i2)+data.getQuick(i4));
		h2i=c2*(data.getQuick(i1)-data.getQuick(i3));
		data.setQuick(i1, h1r+wr*h2r-wi*h2i);
		data.setQuick(i2, h1i+wr*h2i+wi*h2r);
		data.setQuick(i3, h1r-wr*h2r+wi*h2i);
		data.setQuick(i4, -h1i+wr*h2i+wi*h2r);
		wtemp=wr;
		wr=wtemp*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	  }
	if (isign == 1) 
	  {
		h1r=data.getQuick(1);
		data.setQuick(1,h1r+data.getQuick(2));
		data.setQuick(2,h1r-data.getQuick(2));
	  } 
	else 
	  {
		h1r=data.getQuick(1);
		data.setQuick(1,c1*(h1r+data.getQuick(2)));
		data.setQuick(2,c1*(h1r-data.getQuick(2)));
		four1(data,n>>1,-1);
	  }
  }

  /** 
		Convolves data[1..n] with the response function respns[1..m], with m odd and m<=n<br>
		n <b>MUST</b> be a power of 2.<br>
		N.B. Pollutes the response function<br>
		N.B. data.size() = n+1  and we don't use data[0].<br>
		See NR section 13.1 for all the gory detail.
  */

  public static void convlv(DoubleMatrix1D data, int n,DoubleMatrix1D respns, 
							int m, int isign, DoubleMatrix1D ans)
  {
	int i,no2;
	double dum,mag2;

	if (n!=data.size()-1)
	  {
		throw new IllegalArgumentException(" n != data.size() in convlv()");
	  }
	if (n!=respns.size()-1)
	  {
		throw new IllegalArgumentException(" n != respns.size() in convlv()");
	  }
	if (m>n || m<1)
	  {
		throw new IllegalArgumentException(" m>n || m<1 forbidden  in convlv()");
	  }
	DoubleMatrix1D fft = DoubleFactory1D.dense.make((n<<1)+1);

	for (i=1;i<=(m-1)/2;i++)
	  {
		respns.setQuick(n+1-i,respns.getQuick(m+1-i));
	  }
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
	  {
		respns.setQuick(i,0.0);
	  }
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) 
	  {
		if (isign == 1) 
		  {
			dum=ans.getQuick(i-1);
			ans.setQuick(i-1 ,(fft.getQuick(i-1)*dum-fft.getQuick(i)*ans.getQuick(i))/no2);
			ans.setQuick(i, (fft.getQuick(i)*dum+fft.getQuick(i-1)*ans.getQuick(i))/no2);
		  } 
		else if (isign == -1) 
		  {
			if ((mag2=(ans.getQuick(i-1)*ans.getQuick(i-1)+ans.getQuick(i)*ans.getQuick(i))) == 0.0)
			  {
				throw new IllegalArgumentException("Deconvolving at response zero in convlv");
			  }
			dum=ans.getQuick(i-1);
			ans.setQuick(i-1, (fft.getQuick(i-1)*dum+fft.getQuick(i)*ans.getQuick(i))/mag2/no2);
			ans.setQuick(i, (fft.getQuick(i)*dum-fft.getQuick(i-1)*ans.getQuick(i))/mag2/no2);
		  } 
		else 
		  {
			throw new IllegalArgumentException("No meaning for isign in convlv");
		  }
	  }
	ans.setQuick(2, ans.getQuick(n+1));
	realft(ans,n,-1);
	fft=null;   /* force  clean up */
  }



  /** 
		Replaces data[1..2*nn] by its DFT, if isign=1, otherwise
		replaces data[1..2*nn] by nn times its inverse DFT, if
		isign=-1. data[1..2*nn] is a complex array of length nn.  <br>
		nn <b>MUST</b> be a
        power of 2.<br>
		N.B. we don't use data[0].
  */
  public static void four1(double data[], int nn, int isign)
  {
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi,swap;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) 
	  {
		if (j > i) 
		  {
			swap=data[j];data[j]=data[i];data[i]=swap;
			swap=data[j+1];data[j+1]=data[i+1];data[i+1]=swap;
		  }
		m=n >> 1;
		while (m >= 2 && j > m) 
		  {
			j -= m;
			m >>= 1;
		  }
		j += m;
	  }
	mmax=2;
	while (n > mmax) 
	  {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=Math.sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=Math.sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) 
		  {
			for (i=m;i<=n;i+=istep) 
			  {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
			 	data[i+1] += tempi;
			  }
			wtemp=wr;
			wr=wtemp*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		  }
		mmax=istep;
	  }
  }

  /** 
		Zero-offset version<br>.
		Replaces data[0..2*nn-1] by its DFT, if isign=1, otherwise
		replaces data[0..2*nn-1] by nn times its inverse DFT, if
		isign=-1. data[0..2*nn-1] is a complex array of length nn.  <br>
		nn <b>MUST</b> be a
        power of 2.<br>
  */
  public static void four1ZO(double data[], int nn, int isign)
  {
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi,swap;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) 
	  {
		if (j > i) 
		  {
			swap=data[j-1];data[j-1]=data[i-1];data[i-1]=swap;
			swap=data[j];data[j]=data[i];data[i]=swap;
		  }
		m=n >> 1;
		while (m >= 2 && j > m) 
		  {
			j -= m;
			m >>= 1;
		  }
		j += m;
	  }
	mmax=2;
	while (n > mmax) 
	  {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=Math.sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=Math.sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) 
		  {
			for (i=m;i<=n;i+=istep) 
			  {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
			 	data[i] += tempi;
			  }
			wtemp=wr;
			wr=wtemp*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		  }
		mmax=istep;
	  }
  }










  /** 
		Computes the DFT of a set of real data in data1[1..n], going forwards or
		backwards according to isign=1 or -1 respectively. the inverse transform must be
		multiplied by 2/n.<br>
		n <b>MUST</b> be a power of 2.<br>
		N.B. we don't use data[0].<br>
		See NR section 12.3 for all the gory detail.
  */
  public static void realft(double data[], int n, int isign)
	  {
		int i,i1,i2,i3,i4,np3;
		double c1=0.5,c2,h1r,h1i,h2r,h2i;
		double wr,wi,wpr,wpi,wtemp,theta;

		theta=3.141592653589793/(double) (n>>1);
		if (isign == 1) 
		  {
			c2 = -0.5;
			four1(data,n>>1,1);
		  } 
		else 
		  {
			c2=0.5;
			theta = -theta;
		  }
		wtemp=Math.sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=Math.sin(theta);
		wr=1.0+wpr;
		wi=wpi;
		np3=n+3;
		for (i=2;i<=(n>>2);i++) 
		  {
			i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
			h1r=c1*(data[i1]+data[i3]);
			h1i=c1*(data[i2]-data[i4]);
			h2r = -c2*(data[i2]+data[i4]);
			h2i=c2*(data[i1]-data[i3]);
			data[i1]=h1r+wr*h2r-wi*h2i;
			data[i2]=h1i+wr*h2i+wi*h2r;
			data[i3]=h1r-wr*h2r+wi*h2i;
			data[i4] = -h1i+wr*h2i+wi*h2r;
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		  }
		if (isign == 1) 
		  {
			data[1] = (h1r=data[1])+data[2];
			data[2] = h1r-data[2];
		  } 
		else 
		  {
			data[1]=c1*((h1r=data[1])+data[2]);
			data[2]=c1*(h1r-data[2]);
			four1(data,n>>1,-1);
		  }
	  }
  



  /** 
	  Zero offset version.<br>
	  Computes the DFT of a set of real data in data1[0..n-1], going forwards or
	  backwards according to isign=1 or -1 respectively. the inverse transform must be
	  multiplied by 2/n.<br>
	  n <b>MUST</b> be a power of 2.<br>
	  See NR section 12.3 for all the gory detail.
  */
  public static void realftZO(double data[], int n, int isign)
  {
	int i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) 
	  {
		c2 = -0.5;
		four1ZO(data,n>>1,1);
	  } 
	else 
	  {
		c2=0.5;
		theta = -theta;
	  }
	wtemp=Math.sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=Math.sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) 
	  {
		// unit-offset line... i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		i1=i+i-2;
		i2=i1+1;		
		i3=np3-i2-2;		
		i4=i3+1;
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	  }
	if (isign == 1) 
	  {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	  } 
	else 
	  {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1ZO(data,n>>1,-1);
	  }
  }






 /** 
		Given data1[1..nn] and data2[1..n], this calls four1 and
		returns two complex output arrays fft1[1..2n] and fft2[1..2n], which contain the DFT of the
		respective data arrays.<br>
		n <b>MUST</b> be a power of 2.<br>
		N.B. we don't use element data1[0], data2[0], fft1[0], fft2[0]
  */
  public static void twofft(double data1[], double data2[], double fft1[], double fft2[],int n)
  {
	int nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) 
	  {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	  }
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) 
	  {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	  }
  }


  /** 
		Zero-offset version of above<br>
		Given data1[0..nn-1] and data2[0..n-1], this calls four1 and
		returns two complex output arrays fft1[0..2n-1] and fft2[0..2n-1], which contain the DFT of the
		respective data arrays.<br>
		n <b>MUST</b> be a power of 2.<br>
  */
  public static void twofftZO(double data1[], double data2[], double fft1[], double fft2[],int n)
  {
	int nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	//nn3=1+(nn2=2+n+n);
	nn2=1+n+n;
	nn3=1+nn2;
	for (j=1,jj=2;j<=n;j++,jj+=2) 
	  {
		fft1[jj-2]=data1[j-1];
		fft1[jj-1]=data2[j-1];
	  }
	four1ZO(fft1,n,1);
	fft2[0]=fft1[1];
	fft1[1]=fft2[1]=0.0;
	for (j=3;j<=n+1;j+=2) 
	  {
		rep=0.5*(fft1[j-1]+fft1[nn2-j]);
		rem=0.5*(fft1[j-1]-fft1[nn2-j]);
		aip=0.5*(fft1[j]+fft1[nn3-j]);
		aim=0.5*(fft1[j]-fft1[nn3-j]);
		fft1[j-1]=rep;
		fft1[j]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j-1]=aip;
		fft2[j] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	  }
  }



  /** 
		Convolves data[1..n] with the response function respns[1..m], with m odd and m<=n<br>
		n <b>MUST</b> be a power of 2.<br>
		N.B. Pollutes the response function<br>
		N.B. We don't use data[0].<br>
		See NR section 13.1 for all the gory detail.
  */
  public static void convlv(double data[], int n, double respns[], int m, int isign, double ans[])
  {
	int i,no2;
	double dum,mag2;

	double[] fft=new double[1+n<<1];
	for (i=1;i<=(m-1)/2;i++)
	  respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
	  respns[i]=0.0;
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) 
	  {
		if (isign == 1) 
		  {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		  } 
		else if (isign == -1) 
		  {
			if ((mag2=(ans[i-1]*ans[i-1]+ans[i]*ans[i])) == 0.0)
			  {
				throw new IllegalArgumentException("Deconvolving at response zero in convlv");
			  }
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		  } 
		else 
		  {
			throw new IllegalArgumentException("No meaning for isign in convlv");
		  }
	  }
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	fft=null;
  }
  
  /** 
	  Zero-offset version of above<br>
	  Convolves data[0..n-1] with the response function respns[0..m-1], with m odd and m<=n<br>
	  n <b>MUST</b> be a power of 2.<br>
	  N.B. Pollutes the response function<br>
	  See NR section 13.1 for all the gory detail.
  */
  public static void convlvZO(double data[], int n, double respns[], int m, int isign, double ans[])
  {
	int i,no2;
	double dum,mag2; 

	double[] fft=new double[1+n<<1];
	for (i=1;i<=(m-1)/2;i++)
	  respns[n-i]=respns[m-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
	  respns[i-1]=0.0;
	twofftZO(data,respns,fft,ans,n);
	no2=n>>1; 
	for (i=1;i<=n+1;i+=2) 
	  {
		if (isign == 1) 
		  {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		  } 
		else if (isign == -1) 
		  {
			if ((mag2=(ans[i-1]*ans[i-1]+ans[i]*ans[i])) == 0.0)
			  {
				throw new IllegalArgumentException("Deconvolving at response zero in convlv");
			  }
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		  } 
		else 
		  {
			throw new IllegalArgumentException("No meaning for isign in convlv");
		  }
	  }
	ans[1]=ans[n];
	realftZO(ans,n,-1);
	fft=null;
  }




/** 
		Convolves data[1..n] with the response function respns[1..n]<br>
		n <b>MUST</b> be a power of 2.<br>
		This is the same as convlv(), but leaves the response function <i>intact</i>. <br>
		The response must be stored in respns[1..n] in wrap-around order, with
		user-supplied zero padding in the middle.
		N.B. We don't use data[0].<br>
		See NR section 13.1 for all the gory detail.
  */
  public static void convlv_clean(double data[], int n, double respns[], int isign, double ans[])
  {
	int i,no2;
	double dum,mag2;

	double[] fft=new double[1+n<<1];
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) 
	  {
		if (isign == 1) 
		  {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		  } 
		else if (isign == -1) 
		  {
			if ((mag2=(ans[i-1]*ans[i-1]+ans[i]*ans[i])) == 0.0)
			  {
				throw new IllegalArgumentException("Deconvolving at response zero in convlv");
			  }
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		  } 
		else 
		  {
			throw new IllegalArgumentException("No meaning for isign in convlv");
		  }
	  }
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	fft=null;
  }





  /** 
		Zero offset version of above<br>
		Convolves data[0..n-1] with the response function respns[0..n-1]<br>
		n <b>MUST</b> be a power of 2.<br>
		This is the same as convlv(), but leaves the response function <i>intact</i>. <br>
		The response must be stored in respns[0..n-1] in wrap-around order, with
		user-supplied zero padding in the middle.

		See NR section 13.1 for all the gory detail.
  */
  public static void convlvZO_clean(double data[], int n, double respns[], int isign, double ans[])
  {
	int i,no2;
	double dum,mag2;

	double[] fft=new double[1+n<<1];
	twofftZO(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=1;i<=n+1;i+=2) 
	  {
		if (isign == 1) 
		  {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		  } 
		else if (isign == -1) 
		  {
			if ((mag2=(ans[i-1]*ans[i-1]+ans[i]*ans[i])) == 0.0)
			  {
				throw new IllegalArgumentException("Deconvolving at response zero in convlv");
			  }
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		  } 
		else 
		  {
			throw new IllegalArgumentException("No meaning for isign in convlv");
		  }
	  }
	ans[1]=ans[n];
	realftZO(ans,n,-1);
	fft=null;
  }

  /**
	 <i>Zero-offset</i> version of fourn
  */
  public static void fourn(double data[],  int nn[], int ndim, int isign)
  {
	int idim;
	int i1,i2,i3,i2rev,i3rev;
	int ip1,ip2,ip3,ifp1,ifp2;
	int ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot = 1,idim = 1;idim<= ndim;idim++) ntot *= nn[idim-1];
	nprev = 1;
	for (idim = ndim;idim>=1;idim--) {
	  n = nn[idim-1];
	  nrem = ntot/(n*nprev);
	  ip1 = nprev << 1;
	  ip2 = ip1*n;
	  ip3 = ip2*nrem;
	  i2rev = 1;
	  for (i2 = 1;i2<=ip2;i2+=ip1) {
		if (i2 < i2rev) {
		  for (i1 = i2;i1<=i2+ip1-2;i1+=2) {
			for (i3 = i1;i3<=ip3;i3+=ip2) {
			  i3rev = i2rev+i3-i2;
			  tempr = data[i3-1]; data[i3-1] = data[i3rev-1]; data[i3rev-1] = tempr;
			  tempr = data[i3]; data[i3] = data[i3rev]; data[i3rev] = tempr;
			}
		  }
		}
		ibit = ip2 >> 1;
		while (ibit >= ip1 && i2rev > ibit) {
		  i2rev -= ibit;
		  ibit >>= 1;
		}
		i2rev += ibit;
	  }
	  ifp1 = ip1;
	  while (ifp1 < ip2) {
		ifp2 = ifp1 << 1;
		theta = isign*2.0*Math.PI/(ifp2/ip1);
		wtemp = Math.sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = Math.sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (i3 = 1;i3<=ifp1;i3+=ip1) {
		  for (i1 = i3;i1<=i3+ip1-2;i1+=2) {
			for (i2 = i1;i2<=ip3;i2+=ifp2) {
			  k1 = i2;
			  k2 = k1+ifp1;
			  tempr = (double)wr*data[k2-1]-(double)wi*data[k2];
			  tempi = (double)wr*data[k2]+(double)wi*data[k2-1];
			  data[k2-1] = data[k1-1]-tempr;
			  data[k2] = data[k1]-tempi;
			  data[k1-1] += tempr;
			  data[k1] += tempi;
			}
		  }
		  wr = (wtemp = wr)*wpr-wi*wpi+wr;
		  wi = wi*wpr+wtemp*wpi+wi;
		}
		ifp1 = ifp2;
	  }
	  nprev *= n;
	}
  }

  /**
	 Alternative convolution version, based on Dave Hale's SU convolve
	 routines. Autodetect of filter length. (Better to provide
	 explicitly externally)
  */
  public static final void convolve_CWP(double[] data, int n, double[] respns,double[] ans)
  {
	int i,j,kx,lx,js;
	double mx=0;
	boolean trim = false;

	if (trim)
	  {
		// coda
		for(lx=0; respns[lx]!=0 && lx<n; lx++) {mx = max(mx,abs(respns[lx]));}
		// precursor
		for(i=n-1,kx=0; respns[i]!=0 && i>=0; i--,kx++) {mx = max(mx,abs(respns[i]));}
	  }
	// coda. Assume half of filter (4 samples of 8, or 4 samples of 7
	js = max((n-1)/2,0);
	for(lx=js; lx>0; lx--) {if (respns[lx]!=0) break;} 	

	// precursor
	// for(i=n-1,kx=0; respns[i]!=0 && i>=0; i--,kx++) {mx = max(mx,abs(respns[i]));}
	for(i=js+1; i<n; i++) {if (respns[i]!=0) break;}	
	kx = n-i;

	//System.err.printf("CC %d %d lx=%d kx=%d\n",i,n,lx,kx);
	/* shrinking of filter length for weak tails... */
	if (trim)
	  { 
		mx*=1e-3;
		for(; abs(respns[lx])<mx && lx>=0; lx--){};
		for(i=n-1-kx; abs(respns[i])<mx && i<n; i++,kx--){};
	  }

	double[] b = new double[lx+kx+1];
	for (i=-kx,j=n-kx; i<0;i++,j++)
	  {
		b[i+kx] = respns[j];
	  }
	for(i=0,j=0; i<=lx; i++,j++) b[i+kx] = respns[j];
	convolve_kernel(lx+kx+1,-kx,b,n,0,data,n,0,ans);
  }

 /**
	 Alternative convolution version, based on Dave Hale's SU convolve
	 routines. Filter lengths supplied:
	 -kx = sample index of start of precursor in filt[] in indexing
	 lx = sample index of last nonzero sample in filt[] coda
	 filter is in wrap-around form
  */
  public static final void convolve_CWP(double[] data, int nd, double[] filt,int kx,int lx, double[] ans, int na)
  {
	int i,j,Nf=1+kx+lx;
	double[] b = new double[Nf];
	for (i=-kx,j=Nf-kx; i<0;i++,j++)
	  {
		b[i+kx] = filt[j];
	  }
	for(i=0,j=0; i<=lx; i++,j++) b[i+kx] = filt[j];
	convolve_kernel(lx+kx+1,-kx,b, nd,0,data, na,0,ans);
  }

 /**
	 Aggressive version, based on convolve_kernel(). Filter lengths supplied:
	 -kx = sample index of start of precursor in filt[] in indexing
	 lx = sample index of last nonzero sample in filt[] coda
	 filter is in wrap-around form
	 ans[] has desired length na
	 data[] has valid samples 0...nd-1 (but maybe bigger storage).
	 trim = true cuts any convolution elements |filt[i]|<1e-3 max{|filt|} out to compress the filter.
	 Useful also if nd < n
  */
  public static final void convolve_CWP(double[] data, int nd, double[] filt,int kx,int lx, double[] ans,int na, boolean trim)
  {
	int i,j,Nf,Ne=0,Nfs;
	Nf=1+kx+lx;

	if(trim)
	  {
		double mx=0;
		for(i=0; i<=lx; i++) {mx = max(mx,abs(filt[i]));}
		mx*=1e-3;
		for(j=lx; j>0; j--) {if (abs(filt[j])>mx) {lx=j;break;}}
		for(j=lx+1; j<Nf; j++) {if (abs(filt[j])>mx) {Ne=j;break;}}
		kx=Nf-Ne+1;
		Nfs = lx+kx+1; 
	  }
	else Nfs=Nf;

	double[] b = new double[Nfs];
	if(!trim)
	  {
		for (i=-kx,j=Nf-kx; i<0;i++,j++) b[i+kx] = filt[j];
		for(i=0,j=0; i<=lx; i++,j++) b[i+kx] = filt[j];
	  }
	else
	  {
		for (i=-kx,j=Nf-kx; i<0;i++,j++) b[i+kx] = filt[j];
		for(i=0,j=0; i<=lx; i++,j++) b[i+kx] = filt[j];
	  }
	convolve_kernel(lx+kx+1,-kx,b, nd,0,data, na,0,ans);
  }

  // From Dave Hale's BSD-licenced SU version CWP
  /**
   * Computes the convolution (or cross-correlation) of two sequences.
   * <p>
   * Convolution of one-dimensional sequences x and y is defined generically 
   * by the following sum:
   * <pre><code>
   *   z[i] =  sum x[j]*y[i-j]
   *            j
   * </code></pre>
   * In practice, the sequences x, y, and z are non-zero for only finite 
   * ranges of sample indices i and j, and these ranges determine limits 
   * on the summation index j.
   * <p>
   * Specifically, the sequences x, y, and z are stored in arrays with 
   * zero-based indexing; e.g., x[0], x[1], x[2], ..., x[lx-1], where lx
   * denotes the length of the array x. Sequences are assumed to be zero
   * for indices outside the bounds of these arrays.
   * <p>
   * Note that an array index need not equal its corresponding sample index. 
   * For each sequence, we must specify the sample index of the first sample 
   * in the array of sample values; e.g., kx denotes the sample index of x[0]. 
   * With this distinction between sample and array indices in mind, in terms 
   * of arrays x, y, and z, the convolution sum may be rewritten as
   * <pre><code>
   *             jhi
   *   z[i-k] =  sum  x[j]*y[i-j] ; i = k, k+1, ..., k+lz-1
   *             j=jlo
   * </code></pre>
   * where k = kz-kx-ky, jlo = max(0,i-ly+1), and jhi = min(lx-1,i). The 
   * summation limits jlo and jhi ensure that array indices are always in 
   * bounds. The effect of the three first-sample indices is encoded in the 
   * single shift k.
   * <p>
   * For example, if sequence z is to be a weighted average of the nearest 
   * five samples of sequence y, one might use 
   * <pre><code>
   *   ...
   *   x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
   *   conv(5,-2,x,ly,0,y,ly,0,z);
   *   ...
   * </code></pre>
   * In this example, the sequence x is symmetric about the origin, with 
   * first-sample index kx = -2.
   * <p>
   * Cross-correlation is similar to convolution. (Indeed, cross-correlation
   * of x and y equals the convolution of x-reversed and y. The generic 
   * definition of cross-correlation is
   * <pre><code>
   *   z[i] =  sum x[j]*y[i+j]
   *            j
   * </code></pre>
   * Unlike convolution, cross-correlation is not commutative. In other words,
   * the cross-correlation of x and y does not equal the cross-correlation of 
   * y and x.

   */
  public static final void convolve_kernel(
    int lx, int kx, double[] x,
    int ly, int ky, double[] y,
    int lz, int kz, double[] z)
  {
    // If necessary, swap x and y so that x is the shorter sequence.
    // This simplifies the logic below.
    if (lx>ly) {
      int lt = lx;  lx = ly;  ly = lt;
      int kt = kx;  kx = ky;  ky = kt;
      double[] t = x;  x = y;  y = t;
    }

    // Bounds for index i.
    int imin = kz-kx-ky;
    int imax = imin+lz-1;

    // Variables that we expect to reside in registers.
    int i,ilo,ihi,j,jlo,jhi,iz;
    double sa,sb,xa,xb,ya,yb;

    // Off left: imin <= i <= -1
    ilo = imin;
    ihi = min(-1,imax);
    for (i=ilo,iz=i-imin; i<=ihi; ++i,++iz)
      z[iz] = 0.0;

    // Rolling on: 0 <= i <= lx-2 and 0 <= j <= i
    ilo = max(0,imin);
    ihi = min(lx-2,imax);
    jlo = 0;
    jhi = ilo;
    for (i=ilo,iz=i-imin; i<ihi; i+=2,iz+=2,jhi+=2) {
      sa = 0.0;
      sb = 0.0;
      yb = y[i-jlo+1];
      for (j=jlo; j<jhi; j+=2) {
        xa = x[j];
        sb += xa*yb;
        ya = y[i-j];
        sa += xa*ya;
        xb = x[j+1];
        sb += xb*ya;
        yb = y[i-j-1];
        sa += xb*yb;
      }
      xa = x[j];
      sb += xa*yb;
      if (j==jhi) {
        ya = y[i-j];
        sa += xa*ya;
        xb = x[j+1];
        sb += xb*ya;
      }
      z[iz  ] = sa;
      z[iz+1] = sb;
    }
    if (i==ihi) {
      jlo = 0;
      jhi = i;
      sa = 0.0;
      for (j=jlo; j<=jhi; ++j)
        sa += x[j]*y[i-j];
      z[iz] = sa;
    }

    // Middle: lx-1 <= i <= ly-1 and 0 <= j <= lx-1
    ilo = max(lx-1,imin);
    ihi = min(ly-1,imax);
    jlo = 0;
    jhi = lx-1;
    for (i=ilo,iz=i-imin; i<ihi; i+=2,iz+=2) {
      sa = 0.0;
      sb = 0.0;
      yb = y[i-jlo+1];
      for (j=jlo; j<jhi; j+=2) {
        xa = x[j];
        sb += xa*yb;
        ya = y[i-j];
        sa += xa*ya;
        xb = x[j+1];
        sb += xb*ya;
        yb = y[i-j-1];
        sa += xb*yb;
      }
      if (j==jhi) {
        xa = x[j];
        sb += xa*yb;
        ya = y[i-j];
        sa += xa*ya;
      }
      z[iz  ] = sa;
      z[iz+1] = sb;
    }
    if (i==ihi) {
      sa = 0.0;
      for (j=jlo; j<=jhi; ++j)
        sa += x[j]*y[i-j];
      z[iz] = sa;
    }

    // Rolling off: ly <= i <= lx+ly-2 and i-ly+1 <= j <= lx-1
    ilo = max(ly,imin);
    ihi = min(lx+ly-2,imax);
    jlo = ihi-ly+1;
    jhi = lx-1;
    for (i=ihi,iz=i-imin; i>ilo; i-=2,iz-=2,jlo-=2) {
      sa = 0.0;
      sb = 0.0;
      yb = y[i-jhi-1];
      for (j=jhi; j>jlo; j-=2) {
        xa = x[j];
        sb += xa*yb;
        ya = y[i-j];
        sa += xa*ya;
        xb = x[j-1];
        sb += xb*ya;
        yb = y[i-j+1];
        sa += xb*yb;
      }
      xa = x[j];
      sb += xa*yb;
      if (j==jlo) {
        ya = y[i-j];
        sa += xa*ya;
        xb = x[j-1];
        sb += xb*ya;
      }
      z[iz  ] = sa;
      z[iz-1] = sb;
    }
    if (i==ilo) {
      jlo = i-ly+1;
      jhi = lx-1;
      sa = 0.0;
      for (j=jhi; j>=jlo; --j)
    	sa += x[j]*y[i-j];
      z[iz] = sa;
    }
	
    // Off right: lx+ly-1 <= i <= imax
    ilo = max(lx+ly-1,imin);
    ihi = imax;
    for (i=ilo,iz=i-imin; i<=ihi; ++i,++iz)
      z[iz] = 0.0;
  }

  /** 
	 <i>Zero-offset</i> version of rlft3. Everything as per the book,
	 but data and speq are zero-offset. data must have powers-of 2
	 dimensions m ~ slices x rows x columns. Method must operate on
	 au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D types as
	 exposure to internal elements is necessary.
  */
  public static void rlft3(au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D data, double speq[][], int isign)
  {
	int i1,i2,i3,j1,j2,j3,ii3;
	int nn1,nn2,nn3;
	int[] nn = new int[3];
	double theta,wi,wpi,wpr,wr,wtemp;
	double c1,c2,h1r,h1i,h2r,h2i;
	double[] dataElements;

	nn1 = data.slices();
	nn2 = data.rows();
	nn3 = data.columns();
	if (getNearestPowerOfTwo(nn1) != nn1) throw new IllegalArgumentException("rlft3: dim1 not a power of 2");
	if (getNearestPowerOfTwo(nn2) != nn2) throw new IllegalArgumentException("rlft3: dim2 not a power of 2");
	if (getNearestPowerOfTwo(nn3) != nn3) throw new IllegalArgumentException("rlft3: dim3 not a power of 2");


	c1 = 0.5;
	c2 = -0.5*isign;
	theta = isign*(2.0*Math.PI/nn3);
	wtemp = Math.sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = Math.sin(theta);
	nn[0] = nn1;
	nn[1] = nn2;
	nn[2] = nn3 >> 1;
	if (isign == 1) {
	  // get raw elements from data cube and point fourn at it
	  dataElements = data.viewElements();
	  fourn(dataElements,nn,3,isign);
	  for (i1 = 1;i1<=nn1;i1++)
		for (i2 = 1,j2=0;i2<=nn2;i2++) {
		  speq[i1-1][(++j2)-1] = data.getQuick(i1-1,i2-1,0);
		  speq[i1-1][(++j2)-1] = data.getQuick(i1-1,i2-1,1);
		}
	}
	for (i1 = 1;i1<=nn1;i1++) {
	  j1 = (i1 != 1 ? nn1-i1+2 : 1);
	  wr = 1.0;
	  wi = 0.0;
	  for (ii3 = 1,i3 = 1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
		for (i2 = 1;i2<=nn2;i2++) {
		  if (i3 == 1) {
			j2 = (i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
			h1r = c1*(data.getQuick(i1-1,i2-1,0)+speq[j1-1][(j2)-1]);
			h1i = c1*(data.getQuick(i1-1,i2-1,1)-speq[j1-1][j2]);
			h2i = c2*(data.getQuick(i1-1,i2-1,0)-speq[j1-1][(j2)-1]);
			h2r =  -c2*(data.getQuick(i1-1,i2-1,1)+speq[j1-1][j2]);
			data.setQuick(i1-1,i2-1,0, h1r+h2r);
			data.setQuick(i1-1,i2-1,1, h1i+h2i);
			speq[j1-1][(j2)-1] = h1r-h2r;
			speq[j1-1][j2] = h2i-h1i;
		  } else {
			j2 = (i2 != 1 ? nn2-i2+2 : 1);
			j3 = nn3+3-(i3<<1);
			h1r = c1*(data.getQuick(i1-1,i2-1,ii3-1)+data.getQuick(j1-1,j2-1,j3-1));
			h1i = c1*(data.getQuick(i1-1,i2-1,ii3)-data.getQuick(j1-1,j2-1,j3));
			h2i = c2*(data.getQuick(i1-1,i2-1,ii3-1)-data.getQuick(j1-1,j2-1,j3-1));
			h2r = -c2*(data.getQuick(i1-1,i2-1,ii3)+data.getQuick(j1-1,j2-1,j3));
			data.setQuick(i1-1,i2-1,ii3-1,  h1r+wr*h2r-wi*h2i);
			data.setQuick(i1-1,i2-1,ii3, h1i+wr*h2i+wi*h2r);
			data.setQuick(j1-1,j2-1,j3-1, h1r-wr*h2r+wi*h2i);
			data.setQuick(j1-1,j2-1,j3,  -h1i+wr*h2i+wi*h2r);
		  }
		}
		wr = (wtemp = wr)*wpr-wi*wpi+wr;
		wi = wi*wpr+wtemp*wpi+wi;
	  }
	}
	if (isign == -1)
	  {
		dataElements = data.viewElements();
		fourn(dataElements,nn,3,isign);
	  }
  }







  /**
	 Useful public static method for getting nearest power of 2 equal
	 to or bigger than N given, suitable for most zero
	 padding applications. for 32 bit N only, also N>0
  */
  public final static int getNearestPowerOfTwo(int N)
  {
	N--;
	N |= N >> 1;
	N |= N >> 2;
	N |= N >> 4;
	N |= N >> 8;
	N |= N >> 16;
	N++;
	return N;
  }

  //old
  /* FFT power of 2 stuff
	 int n = N;
	 int ng = 0;
	 
	 while ((n >>= 1) > 0 ) ng++;
	 // In general N not a power of 2, so find the right n = 2^ng  bigger than data.length
	 if (N != (1<<ng))
	 {
	 // values.length not a power of 2
	 n=(1<<(ng+1));
	 }
	else
	{
	n = N; // it was already a power of 2
	}
	return n;
  }
  */

  /**
	 A public domain FFT, obtained from http://www.nauticom.net/www/jdtaft/JavaFFT.htm
	 which we include here for diversion.......
  */

  public static double[][] fft_1d(  double[][] array )
  {
	double  u_r,u_i, w_r,w_i, t_r,t_i;
	int     ln, nv2, k, l, le, le1, j, ip, i, n;

	n = array.length;
	ln = (int)( Math.log( (double)n )/Math.log(2) + 0.5 );
	nv2 = n / 2;
	j = 1;
	for (i = 1; i < n; i++ )
	  {
		if (i < j)
		  {
			t_r = array[i - 1][0];
			t_i = array[i - 1][1];
			array[i - 1][0] = array[j - 1][0];
			array[i - 1][1] = array[j - 1][1];
			array[j - 1][0] = t_r;
			array[j - 1][1] = t_i;
		  }
		k = nv2;
		while (k < j)
		  {
			j = j - k;
			k = k / 2;
		  }
		j = j + k;
	  }

	for (l = 1; l <= ln; l++) /* loops thru stages */
	  {       
		le = (int)(Math.exp( (double)l * Math.log(2) ) + 0.5 );
		le1 = le / 2;
		u_r = 1.0;
		u_i = 0.0;
		w_r =  Math.cos( Math.PI / (double)le1 );
		w_i = -Math.sin( Math.PI / (double)le1 );
		for (j = 1; j <= le1; j++) /* loops thru 1/2 twiddle values per stage */
		  {
			for (i = j; i <= n; i += le) /* loops thru points per 1/2 twiddle */
			  {
				ip = i + le1;
				t_r = array[ip - 1][0] * u_r - u_i * array[ip - 1][1];
				t_i = array[ip - 1][1] * u_r + u_i * array[ip - 1][0];

				array[ip - 1][0] = array[i - 1][0] - t_r;
				array[ip - 1][1] = array[i - 1][1] - t_i; 

				array[i - 1][0] =  array[i - 1][0] + t_r;
				array[i - 1][1] =  array[i - 1][1] + t_i;  
			  } 
			t_r = u_r * w_r - w_i * u_i;
			u_i = w_r * u_i + w_i * u_r;
			u_r = t_r;
		  } 
	  }  
	return array;
  } /* end of FFT_1d */
}
















