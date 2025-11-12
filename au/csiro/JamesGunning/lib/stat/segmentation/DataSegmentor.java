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


package au.csiro.JamesGunning.lib.stat.segmentation;

import java.io.*;
import java.util.*; 

import cern.colt.list.*;
import cern.jet.stat.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.lib.util.*;
 
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;

/**
   Class with useful static methods for data segmentation, as per 
   D. M. Hawkins and J. A. ten Krooden, "A Review of Several 
   Methods of Segmentation", Geomathematical and Petrophysical 
   Studies in Sedimentology.

   Also, dynamic programming methods for Maximum likelihood methods
   added Sep 2004, inspired from D. M. Hawkins, "Fitting multiple
   change-point models to data", Computational and Statistical Data
   Analysis, 2001.

   All methods that return IntArrayLists are lists of zero-offset
   changepoints, with the convention that the subsets 
   data[b_i + 1,b_i + 2,......b_{i+1}] inclusive are "segments", and the IntArrayLists
   contain the set of changepoints {b_i}. The first should always be 0, the last data.length-1
   
   @author James.Gunning@csiro.au
   @version 2.0, 12/9/2004
*/

public final class DataSegmentor
{
  // set this to generate full range of optimal segmentations up to the BIC optimum
  private boolean generateAllDPSegmentations = false;
  // set this to print results
  private boolean printAllDPSegmentations = false;

  private IntArrayList[] changepts_array;
  private IntArrayList changepts;   // will point at useful member of changepts_array as required
  private double[] data;
  private ArrayList <SortableNumbers> decimationLookups=null;
  private int minSegSize = 0;
  private int prescribedNSegments = -1; // set this >0 to show segmentation with prescribed number
  private float minStdDevFrac =0.0f;
  private boolean usingDecimation = false;
  boolean doPieceWiseLinear = false;
  private boolean verb;

  public static	int MAX_DATA_SIZE_FOR_ML_METHOD = 1400; // MAX IS REALLY 1400;

  private static float BAYES_PRIOR_SEG_PENALTY = 1e20f;
  
  /**
	 Make a new default segmentor based on the given data array
  */
  public DataSegmentor(double[] data)
  {
	this.data=data;
  }

  /**
	 Make a new segmentor based on the given data array, with specified minimum segment size, and
	 fractional acceptable std_deviation. 
 	 @param minStdDevFrac least acceptable segment standard deviation as fraction of data range (default 0)
	 @param prescribedNSegments 
	 Set prescribedNSegments>1 to fix number of segments.
	 Set prescribedNSegments=-1 for autmatic selection of number of segments.
	 @param minSegSize set to restrict segments to at least this size (measured in number of samples)
	 @param verb set true for verbose (stderr) behaviour
  */
  public DataSegmentor(double[] data,float minStdDevFrac,int minSegSize,int prescribedNSegments,boolean verb)
  {
	this.data=data;
	this.minStdDevFrac = minStdDevFrac;
	this.minSegSize = minSegSize;
	this.prescribedNSegments = prescribedNSegments;
	this.verb=verb;
  }

  /**
	 Make a new segmentor based on the given data array, with specified minimum segment size, and
	 fractional acceptable std_deviation. 
 	 @param minStdDevFrac least acceptable segment standard deviation as fraction of data range (default 0)
	 @param prescribedNSegments 
	 Set prescribedNSegments>1 to fix number of segments.
	 Set prescribedNSegments=-1 for autmatic selection of number of segments.
	 @param minSegSize set to restrict segments to at least this size (measured in number of samples)
	 @param verb set true for verbose (stderr) behaviour
	 @param doPieceWiseLinear set true for piecewise-linear fitting
  */
  public DataSegmentor(double[] data,float minStdDevFrac,int minSegSize,int prescribedNSegments,boolean verb, boolean doPieceWiseLinear)
  {
	this.data=data;
	this.minStdDevFrac = minStdDevFrac;
	this.minSegSize = minSegSize;
	this.prescribedNSegments = prescribedNSegments;
	this.verb=verb;
	this.doPieceWiseLinear = doPieceWiseLinear;
  }

  /**
	 Test harness. Produces stream of numbers x,y of a
	 quasi segmented process, followed by a stream of change-points
	 and associated cluster-means for the segment. 
	 Suggest: pipe to a graphing program

	 stdin | java au.csiro.JamesGunning.lib.stat.segmentation.DataSegmentor [-d]| graph -TX

	 If stdin is supplied, an input stream of ascii numbers (either
	 one per line or white space separated x,y pairs per line) is read
	 and segmented. If stdin not supplied, an internal fake sequence
	 is generated and used.  Supplying -d forces use of maximum
	 likelihood methods (dynamic programming).
  */
  public static void main(String argv[])
  {
	int i,j,minSeg=0,nSeg=-1,N=200;
	double[] dat,xdat=null;
	BufferedReader inReader;
	boolean useFakeData=false;
	boolean dyn=false;
	boolean doPieceWiseLinear = false;
	double cutoffscale=1.0;
	float minStd=0f;
	String filename=null;
	boolean verb=false;
	boolean cascade=false;
	ErrorHandler EH = new ErrorHandler();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("DataSegmentor",EH);

   	options.addMandatoryArg(1,"data_file","data file (single column 'y' or pairs 'x y', per line, ascii). Use single '-' for "
							+"stdin. If -  (stdin) specified and no input found, fake internal data is used. "
							+"Original data and segmented version are written to stdout");	
	options.addOption("-d",null,0,null,null,"Use dynamic programming (DP), i.e. strict max. likelihood method");
	options.addOption("-c",null,0,null,null,"Show (plotutils) graph-friendly cascade of segmentations (up to BIC optimum) when DP used");
	options.addOption("-m",null,1,"int","size","Minimum segment size (-d method only) (in samples). Default 0.");
	options.addOption("-n",null,1,"int","size","Use specified number of segments (n>0) (DP only)");
	options.addOption("-s",null,1,"float","fraction","Least acceptable segment standard deviation "
					  +"(DP method only) (as fraction of data range). Default 0.");
	options.addOption("-l",null,0,null,null,"Use piecewise-linear segments");
	options.addOption("--blocking-cutoff",null,1,"float","cutoff","specify aggregative method variance cutoff (recommended value = 1)");
	options.addOption("-v",null,0,null,null,"Be verbose to stderr");

	options.parseOptions(argv,true);

	// data file
	filename = options.getMandatoryArg(1);

	if (options.isOptionSet("-d")) dyn=true;
	if (options.isOptionSet("-c") && options.isOptionSet("-d")) cascade=true;
	if (options.isOptionSet("-v")) verb=true;

	if (options.isOptionSet("-m")) 
	  {
		minSeg = ((Integer)options.getOptionValues("-m")[0]).intValue();
		if (verb) System.err.println("Using minimum seg size "+minSeg);
	  }
	if (options.isOptionSet("-n")) 
	  {
		nSeg = ((Integer)options.getOptionValues("-n")[0]).intValue();
		if (nSeg<1)
		  {
			System.err.println("Number of segments must be >0");
			System.exit(0);
		  }
		if (verb) System.err.println("Using fixed number of segments "+nSeg);
	  }
	if (options.isOptionSet("-s")) 
	  {
		minStd = ((Float)options.getOptionValues("-s")[0]).floatValue();
	  }
	if (options.isOptionSet("-l")) 
	  {
		doPieceWiseLinear = true;
	  }
	if (options.isOptionSet("--blocking-cutoff")) 
	  {
		cutoffscale = ((Float)options.getOptionValues("--blocking-cutoff")[0]).floatValue();
	  }

	if (!filename.equals("-"))
	  {
		// from file
		try
		  {
			inReader = new BufferedReader(new FileReader(filename));
		  }
		catch  (IOException e) 
		  {
			System.err.println("couldn't open " + filename);return;
		  }
	  }
	else
	  {
		// try stdin?
		InputStreamReader iReader = new InputStreamReader(System.in);
		inReader  = new BufferedReader(iReader);
		try
		  {
			// Try to detect stdin in 1000ms
			boolean haveStdIn=false;
			long tEnd = System.currentTimeMillis()+1000;
			for (;;)
			  {
				if (iReader.ready())
				  {
					haveStdIn=true;
					break;
				  }
				if (System.currentTimeMillis()>tEnd) break;
			  }
			if (!haveStdIn) useFakeData = true;
		  }
		catch (IOException e) {throw new RuntimeException("couldn't check stdin bytes available");}
	  }

	if (useFakeData)
	  {
		if (verb) System.err.println("Using fake data");
		int i_s,ns=0;
		double scale=1,ol=0,offset;
		dat = new double[N];
		Random gausGen = new Random();
		for(i=0;i<N;)
		  {
			int block_size = 30+(int)(20*Math.random());
			do
			  {
				offset= 4*Math.random();
			  } while (Math.abs(offset-ol)<1);
			ol=offset;
			ns++;
			if (!doPieceWiseLinear)
			  {
				i_s=i;
				for(j=0;j<block_size && i<N;j++,i++)
				  {
					dat[i]=scale*(offset+gausGen.nextGaussian());
				  }
				System.err.println("True segment: "+i_s+" "+(i-1));
			  }
			else
			  {
				i_s=i;
				for(j=0;j<block_size && i<N;j++,i++)
				  {
					dat[i] = scale*(offset+(offset*j/(double)N)+gausGen.nextGaussian());
				  }
				System.err.println("True segment: "+i_s+" "+(i-1));
			  }
		  }

		if (verb) System.err.println("Fake data has "+ns+" segments");
	  }
	else
	  {
		// read ascii data inReader
		if (verb) System.err.println("Using user's data");
		String line=null;
		String ws="[ \t]";
		String[] words;
		boolean twoCol=false;
		DoubleArrayList dList = new DoubleArrayList();
		DoubleArrayList xList = new DoubleArrayList();
		i=0;
		try
		  {
			while (((line=inReader.readLine()) != null))
			  {
				words = line.split(ws);
				if (i==0)
				  {
					if (words.length>1) twoCol=true;
				  }
				if (twoCol)
				  {
					xList.add(java.lang.Double.parseDouble(words[0]));
					dList.add(java.lang.Double.parseDouble(words[1]));
				  }
				else dList.add(java.lang.Double.parseDouble(words[0]));
			  }
		  }
		// file i/o exception
		catch (IOException e) {System.out.println("couldn't read data");}
		dList.trimToSize();
		dat = dList.elements();
		if (twoCol) xdat = xList.elements();
	  }

	// segment it

	DataSegmentor segmentor = new DataSegmentor(dat,minStd,minSeg,nSeg,verb,doPieceWiseLinear);
	if (!dyn)
	  {
		//segmentor.aggregativeSegment(false,1.0);
		segmentor.aggregativeSegment(false,cutoffscale);
		segmentor.debugPrintSegmentationResults(xdat);
	  }
	else 
	  {
		// segmentor.dynamicProgrammingMLSegment();
		segmentor.showCascade(cascade);
		segmentor.dynamicProgrammingMLSegmentWithDecimationIfRequired(cascade);
		segmentor.debugPrintSegmentationResults(xdat);
	  }
  }

  /**
	 Split-moving window method. Webster 1973. Suggested window=20.
	 Returns list of suggested change points.
	 @param window. Half window height as per Hawkins' discussion.
  */
  public IntArrayList splitWindowSegment(int window)
  {
	int i;
	int h;
	int N=data.length;
	double deltaQ;
	double cutoff=0.1;

	h=(int)Math.abs(window); // window half-size

	DoubleArrayList list = new DoubleArrayList(data);
	
	// suitably scaled cutoff
	cutoff=0.5*h*Descriptive.sampleVariance(list,Descriptive.mean(list));

 	changepts = new IntArrayList();
	changepts.add(0);
	for (i=h;i<N-h;i++)
	  {
		deltaQ = segmentVariance(i-h,i+h,data)
		  -segmentVariance(i-h,i,data)
		  -segmentVariance(i,i+h,data);
		if (deltaQ > cutoff) changepts.add(i);
	  }
	changepts.add(N);
	changepts.trimToSize();

	return changepts;
  }
  /**********************************************/

  /**
	 Aggregational upscale, cf p 120 Hawkins
  */
  public IntArrayList aggregativeSegment(boolean useDefaultCutoff, double cutoffscale)
  {
	int i,j,m,min_indx,max_indx,max_m_indx=0;
	int j_delete=0;
	int N=data.length;
	double deltaQ=0;
	double cutoff = 0.1;
	double split_cutoff;
	double min,max,max_m;
	boolean splitPerformed,mergePerformed;
	DoubleArrayList dQ;
	DoubleArrayList list = new DoubleArrayList(data);
	
	// suitably scaled cutoff

	if (useDefaultCutoff)
	  {
		cutoff = Descriptive.sampleVariance(list,Descriptive.mean(list));
	  }
	else
	  {
		// JG - (12/01/2009) - keep more detail (e.g. extra 0.5)
		cutoff = cutoffscale*Descriptive.sampleVariance(list,Descriptive.mean(list));
	  }
	split_cutoff=10*cutoff;

 	changepts = new IntArrayList();

	// make full suite of change pts.
	//for (i=0;i<N;i++) changepts.add(i);
	for (i=-1;i<N;i++) changepts.add(i); // left edge must be -1 (JG Mar 09)??

	changepts.trimToSize();
	
	dQ = new DoubleArrayList(changepts.size());
	dQ.add(0);      // dummy placeholder for element 0
	for (j=1;j<changepts.size()-1;j++)
	  {
		dQ.add(segmentVariance(changepts.get(j-1),changepts.get(j+1),data)
			   -segmentVariance(changepts.get(j-1),changepts.get(j),data)
			   -segmentVariance(changepts.get(j),changepts.get(j+1),data));
	  }

	// now start deleting them
	for (int nseg=1;;nseg++)
	  {
		// MERGING
		mergePerformed=false;
		min=Double.MAX_VALUE;
		min_indx=0;
		for (j=1;j<changepts.size()-1;j++)
		  {
			if (dQ.get(j) < min) 
			  {
				min = dQ.get(j);
				min_indx = changepts.get(j);
				j_delete = j;
			  }
		  }
		if (min < cutoff)
		  {
			changepts.remove(j_delete);
			// update dQ array and recompute elemnts of list affected by this changepoint removal
			if (j_delete-1 >=0 ) dQ.remove(j_delete-1);
			for (j=Math.max(1,j_delete-1); j<=Math.min(changepts.size()-2,j_delete); j++)
			  {
				deltaQ =  segmentVariance(changepts.get(j-1),changepts.get(j+1),data)
				  -segmentVariance(changepts.get(j-1),changepts.get(j),data)
				  -segmentVariance(changepts.get(j),changepts.get(j+1),data);
				dQ.set(j,deltaQ);
			  }
			mergePerformed=true;
		  }


		// SPLITTING
		splitPerformed=false;
		// not working yet.....FIXME sometime!!!!!!
		if (false)
		  {
			max=-Double.MAX_VALUE;
			max_indx=0;
			int insert_index=0;
			for (j=1;j<changepts.size()-1;j++)
			  {
				max_m=-Double.MAX_VALUE;
				max_m_indx=0;
				// possible splits of each segment
				for (m=changepts.get(j)+1;m<changepts.get(j+1);m++)
				  {
					deltaQ = Math.abs(
									  segmentVariance(changepts.get(j-1),changepts.get(j),data)
									  -segmentVariance(changepts.get(j-1),m,data)
									  -segmentVariance(m,changepts.get(j),data));
					if (deltaQ > max_m) 
					  {
						max_m=deltaQ;
						max_m_indx=m;
					  }
				  }
				// keep track of best possible split of all possible splits of
				// all segments
				if (max_m > max) 
				  {
					max = max_m;
					max_indx = max_m_indx;
					insert_index=j;
				  }
			  }
			//if (verb) System.err.println("best split "+max+" > "+cutoff+" at "+max_indx);
			if (max>split_cutoff)
			  {
				if (verb) System.err.println(max+" "+cutoff+" "+changepts.size()+" added: "+max_m_indx);
				changepts.beforeInsert(insert_index,max_indx);
				splitPerformed=true;
			  }
		  }
		if (!mergePerformed && !splitPerformed) break;
	  }
	// endpoint
	//for(i=0;i<Math.min(10,changepts.size());i++) System.err.println("CP: "+changepts.get(i));
	if ((changepts.get(0)==-1) && (changepts.get(1)==0)) changepts.remove(0);
	if (changepts.get(0) != 0) changepts.set(0,0);
	if (!changepts.contains(N-1)) changepts.add(N-1);

	changepts.trimToSize();
 
	return changepts;
  }

  private double xValue(int j, double[] xarray)
  {
	return (xarray==null)? (double)j : xarray[j];
  }

  /**
	 Print to stdout successive streams of the raw data and the
	 upscaled (segmented) data. Suitable for graphing, debugging etc
  */
  public void debugPrintSegmentationResults(double[] xarray)
  {
	int i,j,js,je;
	double xs,xe,min,max,s=0,v=0;

	if (printAllDPSegmentations)
	  {
		min=Double.MAX_VALUE;
		max=-Double.MAX_VALUE;
		for (i=0;i<data.length; i++) 
		  {
			min=Math.min(min,data[i]);
			max=Math.max(max,data[i]);
		  }
		v= (max-min)/30.0;
		s = min;
	  }
	// Diagnostic output
	// print orig data
	for (j=0; j<data.length; j++) System.out.println(xValue(j,xarray)+" "+data[j]);
	System.out.println("");

	System.out.println("#segmentation:");

	// now Arithmetically upscaled segments, or linear trends
	if (doPieceWiseLinear)
	  {
		int Ns;
		double Sy,Syy,Sxy,Sxx,Sx,m,b;
		for (j=0; j<changepts.size()-1; j++)
		  {
			js = changepts.get(j)+1; 
			if(js<0) continue;
			je = changepts.get(j+1);
			System.err.println("Guess segment: "+js+" "+je);
			for (i=js,Sy=0,Syy=0,Sxy=0; i<=je; i++)
			  {
				Sy += data[i];
				Syy += data[i]*data[i];
				Sxy += i*data[i];
			  }
			Ns = je-js+1;
			Sx = Ns *(je + js)/2.0;			
			Sxx = Ns *(je*je - 2*je*(js-1) + (js-2)*js)/12.0;
			Syy -= Sy*Sy/Ns;
			Sxy -= Sx*Sy/Ns;
			m = Sxy/Sxx;
			b = (Sy/Ns) - ((je+js)/2.0)*m;
			//xs = (js>0) ? 0.5*(xValue(js,xarray)+xValue(js-1,xarray)) : xValue(js,xarray);
			//xe = (je>0) ? 0.5*(xValue(je,xarray)+xValue(je-1,xarray)) : xValue(je,xarray);
			xs = xValue(js,xarray);
			xe = xValue(je,xarray);
			System.out.println(xs+" "+(b+js*m));
			System.out.println(xe+" "+(b+je*m));
			System.out.println("");
		  }
	  }
	else
	  {
		int nd = data.length;
        for (j=0; j<changepts.size()-1; j++)
		  {
			double sum=0;
			js = (j>0)? changepts.get(j): -1;
			je = changepts.get(j+1);
			for (i=js+1; i<=je; i++) sum+=data[i];
			sum/=(double)(je-js);
			System.err.println("Guess constant-seg: "+(js+1)+" "+je);
            xs = (js==-1) ? xValue(0,xarray) : ((js>0) ? 0.5*(xValue(js,xarray)+xValue(js+1,xarray)) : xValue(js,xarray));
            xe = (je==nd-1) ? xValue(je,xarray) : ((je>0) ? 0.5*(xValue(je,xarray)+xValue(je+1,xarray)) : xValue(je,xarray));
            System.out.printf("%f %f%n%f %f%n",xs,sum,xe,sum);
		  }
	  }
	System.out.println("");

	if (printAllDPSegmentations)
	  {
		// insert gap
		s-=v;
		// loop over optimal segmentations up to BIC
		for (j=0;j<changepts_array.length;j++)
		  {
			// useful reference
			changepts = changepts_array[j];
			// Graph some details of the DP solution: this shows a graph-friendly cascade of the 
			// optimal segmentation up to the BIC-optimal choice
			for(i=0;i<changepts.size();i++)
			  {
				// segment
				js = changepts.get(i); // -0.5;
				xs = (js>0) ? 0.5*(xValue(js,xarray)+xValue(js-1,xarray)) : xValue(js,xarray);
				System.out.println(xs+" "+s);
				System.out.println(xs+" "+(s-v));
				System.out.println("");
			  }
			s -=1.5*v;
		  }
	  }
  }

  // Method to get segment variance. (Hawkins Q_{i_start,i_end}).
  private double segmentVariance(int i_start,int i_end,double[] data)
  {
	int i;
	double S=0,xbar=0,tmp;

	for (i=i_start+1; i<=i_end; i++) xbar+=data[i];
	xbar/=(double)(i_end-i_start);
	for (i=i_start+1; i<=i_end; i++) S+=(tmp=(data[i]-xbar))*tmp;
	return S; 
  }


  /** 
	  Dynamic programming methods of Hawkins, 
	  "Fitting multiple change-point models to data", Computational and Statistial Data Analysis, 2001.

	  Wrapping call to segment data using  recursive calls until BIC is satisfied. Start with guess of
	  Min(10,data.size/50) segments.

	  Returns null if there are memory exceptions (these methods are O(n^2) on memory, n=data size)
  */

  public IntArrayList dynamicProgrammingMLSegment()
  {
	int currentMaxSegments,suggestedSize,newSize;


	suggestedSize = Math.max((int)(data.length/50.0),10);
	if (prescribedNSegments > suggestedSize)
	  {
		// specified segmentation
		suggestedSize = prescribedNSegments;
		if (prescribedNSegments >= data.length-1)
		  {
			System.err.println("Suggested number of segments "+suggestedSize+" too many for data size");
			return null;
		  }
	  }
	try
	  {
		for(;;)
		  {
			newSize = dynamicProgrammingMLSegmentUsingBIC(suggestedSize);
			if (newSize < (suggestedSize-1))
			  {
				// must have found BIC solution. Break
				break;
			  }
			else
			  {
				// increase and call again
				suggestedSize = Math.min(data.length,suggestedSize*2);
				if (verb) System.err.println("Trying up to "+suggestedSize+" segments on "+data.length+" data");
			  }
		  }
	  }
	catch (OutOfMemoryError e)
	  {
		// Big systems will require O(n^2) memory. Catch this
		return null;
	  }
	return changepts;
  }
  
  /** 
	  Segment data using ML methods based on at most nSegments segments. Returns the 
	  BIC-calculated suggested number of segments
  */
  private int dynamicProgrammingMLSegmentUsingBIC(int nSegments)
  {
	int i,j=0,h,m,n,r,segSize=0;
	double s=0,v,min,max,var,min_variance;
	DoubleMatrix2D W,F;
	int[][] H;
	int[] fullChPoints;
	int rmin,count;

	n = data.length;
	rmin = nSegments;

	if (doPieceWiseLinear)
	  {
		// piecewise-linear algorithm from Kehagias
		
		// (1) segment variances
		int k,kL,kR,kmax,kmin;
		int[][] z;
		double[][] d,c,Qij,Qijp1;
		double[] Rij,Rijp1,Aij,Aijp1,dA,z1,d_ptr,QdA;
		double xi,e_ji;
		d = new double[n+1][n+1];	
		Rij = new double[2];
		Rijp1 = new double[2];
		Aij = new double[2];
		Aijp1 = new double[2];
		dA = new double[2];
		Aij = new double[2];
		Qij = new double[2][2];
		Qijp1 = new double[2][2];
		QdA = new double[2];
		double[][][] A=null;

		if (verb) A = new double[n+1][n+1][2]; // stored regression coeffs (not in Kehagias)

		for(i=1;i<=n;i++)
		  {
			xi = data[i-1]; // data shifted to unit-offset ****
			Qij[0][0] = 1; Qij[0][1] = i; Qij[1][0] = i; Qij[1][1] = i*i;
			Rij[0] = xi;
			Rij[1] = i * xi;
			// should be: linSolve2x2(Qij,Rij,Aij);
			// but system is degenerate with 2 unknowns & 1 data point, so impose simple solution...
			Aij[0] = xi; /* =Rij[0]/Qij[0][0] */; Aij[1]=0.0;

			d_ptr = d[i]; 
			d_ptr[i]=0;
			for (j=i; j<n; j++)
			  {
				Qijp1[0][0] = Qij[0][0] + 1;
				Qijp1[0][1] = Qij[0][1] + (j+1);
				Qijp1[1][0] = Qijp1[0][1];
				Qijp1[1][1] = Qij[1][1] + (j+1)*(j+1);

				Rijp1[0] = Rij[0] + data[j];
				Rijp1[1] = Rij[1] + (j+1)*data[j];

				linSolve2x2(Qijp1,Rijp1,Aijp1);

				dA[0] = Aijp1[0] - Aij[0];
				dA[1] = Aijp1[1] - Aij[1];

				QdA[0] = Qij[0][0] * dA[0] + Qij[0][1] * dA[1];
				QdA[1] = Qij[1][0] * dA[0] + Qij[1][1] * dA[1];

				d_ptr[j+1] = d_ptr[j] + (dA[0]*QdA[0] + dA[1]*QdA[1])
				  + MathUtil.sqr(data[j]-(Aijp1[0] + (j+1)*Aijp1[1]));
				// store regression coeffs
				if (verb) 
				  {
					A[i][j][0] = Aij[0];
					A[i][j][1] = Aij[1];
				  }
				// save for next pass in loop..
				ArrayUtils.copy(Qijp1,Qij);

				Aij[0] = Aijp1[0]; Aij[1] = Aijp1[1];
				Rij[0] = Rijp1[0]; Rij[1] = Rijp1[1];
			  }
			// store regression coeffs
			if (verb) 
			  {
				A[i][n][0] = Aij[0];
				A[i][n][1] = Aij[1];
			  }
		  }

		// (2) Dynamic programming bit
		// initialization
		kmax = nSegments+1;
		z = new int[kmax+1][n+1];
		c = new double[kmax+1][n+1];  // initializes as zero
		double[] c_ptr,ckm1_ptr;
		c_ptr = c[1];
		for (i=1;i<=n;i++) c_ptr[i] = d[1][i];

		// minimization
		int jmin=1;
		int[] z_ptr;
		for (k=2; k<=kmax; k++)
		  {
			c_ptr = c[k];
			ckm1_ptr = c[k-1];
			z_ptr = z[k];
			c_ptr[0]=0.0;
			for (i=1; i<=n;i++)
			  {
				for (j=1, min=Double.MAX_VALUE; j<=i-1; j++) 
				  {
					e_ji = ckm1_ptr[j] + d[j+1][i]; 
					if (e_ji < min)
					  {
						min = e_ji;
						jmin = j;
					  }
				  }
				c_ptr[i] = min;
				z_ptr[i] = jmin;
			  }
		  }
		// backtracking
		int[][] t_hat = new int[kmax+1][kmax+1];
		for (k=1; k<=kmax; k++)
		  {
			t_hat[k][k] = n;
			for (j=k; j>=1; j--) t_hat[k][j-1] = z[j][t_hat[k][j]];
		  } 
		// diagnostics
		int k_scan_min,k_scan_max;

		if (prescribedNSegments > 0)
		  {
			// scanning minima in loop below not required. Compute for specified number directly
			k_scan_min = k_scan_max = prescribedNSegments;
		  }
		else
		  {
			k_scan_min = 1;
			k_scan_max = kmax;
		  }
		v = Math.log(n);
		kmin=0;
		min = Double.MAX_VALUE;
		for (k = k_scan_min; k <= k_scan_max; k++)
		  {
			// Log lik when we write term from variance estimate
			var = c[k][n] / ((n>1) ? (double)(n-1) : 1.0);
			// 3k parameters => 1.5 = 3/2. See p 4 Kehagias (note I work with extra factor of 0.5)
			s = ((var>0) ? 0.5*n*Math.log(var) : 0) + 1.5*k*v;
			if (s < min)
			  {
				kmin = k;
				min = s;
			  }
			if (verb) System.err.println("BIC: segs= "+k+" bic= "+s+
										 " loglik= "+((var>0) ? 0.5*n*Math.log(var) : 0)
										 +" var= "+var);
		  }
		if (verb) System.err.println("Suggest "+(kmin)+" segments");
		rmin = kmin; // for return below

		// Store segmentation results suitably
		int kMaxModelRange,kMinModelRange,model_indx;
		if (generateAllDPSegmentations)
		  {
			kMaxModelRange = kmin;
			kMinModelRange = 1;
		  }
		else
		  {
			kMaxModelRange = kmin;
			kMinModelRange = kmin;
		  }
		changepts_array = new IntArrayList[kMaxModelRange-kMinModelRange+1];


		for (i = kMinModelRange,model_indx=0; i<=kMaxModelRange; i++,model_indx++)
		  {
			changepts_array[model_indx] = new IntArrayList();
			// useful reference
			changepts = changepts_array[model_indx];
			// flush list if required
			if (i >kMinModelRange) changepts.removeAll(changepts);
			// Form the chpt list. Add start and end of array
			for (k=0; k<i; k++)
			  {
				// start of segment x[kL+1,...kR] = data[kL....kR-1]
				// kR = t_hat[i][k+1];
				kL = t_hat[i][k];
				if(kL>=0) changepts.add(kL-1);
			  }
			// .. and add and end of array if needed
			if (!changepts.contains(n-1))	changepts.add(n-1);
			// sort
			changepts.sort(); // may not be needed
			changepts.trimToSize();
		  }
		if (false)
		  {
			// debugging....
			min=Double.MAX_VALUE;
			max=-Double.MAX_VALUE;
			for (i=0;i<n; i++) 
			  {
				min=Math.min(min,data[i]);
				max=Math.max(max,data[i]);
			  }
			// print data
			for (i=0;i<n; i++) System.out.println(i+" "+data[i]);
			System.out.println("");
			// then piecewise regressions..
			System.out.println("#m=2");
			for (k=0; k<kmin; k++)
			  {
				kL = t_hat[kmin][k]; // start of segment- x[kL+1....kR] = data[kL....(kR-1)]
				kR = t_hat[kmin][k+1]; // end
				double a = A[kL][kR][0]; // (I think)
				double b = A[kL][kR][1];
				System.out.printf("%d %g%n%d %g%n",kL,a+b*(kL),kR-1,a+b*(kR-1)); // to align with zero-offset data
			  }
			// System.exit(0);
		  }
	  }
	else
	  {
		// Piecewise constant models. Originally from  Hawkins 
		// Precompute and save W array
		double x;
		W = new DenseDoubleMatrix2D(n,n+1);
		min = Double.MAX_VALUE;
		max = -Double.MAX_VALUE;
		for (h=0;h<n;h++)
		  {
			s=v=0.0;
			// ORIG for (m=h+1;m<n;m++)
			for (m=h+1;m<=n;m++)
			  {
				x = data[m-1]; // data offset so changepoints come out unit-offset ****
				s += x;
				v += x * x;
				W.setQuick(h,m, v - s*s/(double)(m-h));
			  }
			// System.out.println("WVAR: "+W.get(h,n-1)/((n-1)-(h+1)+1));
			min = Math.min(data[h],min);
			max = Math.max(data[h],max);
		  }
		min_variance=(max-min)*minStdDevFrac; 
		min_variance *= min_variance; // convert to variance

		// Step 1
		F = new DenseDoubleMatrix2D(n,n);
		for (m=0;m<n;m++) F.setQuick(0,m,W.getQuick(0,m));
	
		// Step 2
		H = new int[nSegments][n];
		for (r=1;r<nSegments;r++)
		  {
			for (m=0;m<n;m++)
			  {
				min=Double.MAX_VALUE;
				// HACK
				boolean hasTooLongSegment = false;
				// ...HACK
				for (h=0;h<m;h++)
				  {
					// model with Q = W
					s = F.getQuick(r-1,h) + W.getQuick(h,m);

					// Bayesian hack for short segments
					if (minSegSize>0)
					  {
						if (!usingDecimation) 
						  {
							segSize = m-h;
						  }
						else
						  {
							segSize = decimationLookups.get(m).index - decimationLookups.get(h).index;
						  }
						if (segSize <= minSegSize) 
						  {
							// System.err.println("Using SEG penalty "+m +" "+h+" "+(m-h));
							s += BAYES_PRIOR_SEG_PENALTY;
						  }
					  }
					// ...HACK
					if (s < min)
					  {
						min = s;
						H[r][m] = h;
					  }
				  }
				F.setQuick(r,m,min);
			  }
		  }

		// First, get optimum number of changepoints from a BIC like argument
		// See Yao 1988 or Kehagias etc...
		min = Double.MAX_VALUE;
		v = Math.log(n);
		int r_scan_min,r_scan_max;
		// BIC = 0.5 * Np * log(n), where 
		// Np = N(means) + N(ch.points) +1 x variance = ((r+1) + r + 1) = 2r+2
		if (prescribedNSegments > 0)
		  {
			// scanning minima in loop below not required. Compute for specified number directly
			r_scan_min = r_scan_max = prescribedNSegments-1;
		  }
		else
		  {
			r_scan_min = 1;
			r_scan_max = nSegments-1;
		  }

		for (r=r_scan_min; r<=r_scan_max; r++)
		  {
			// Hawkins - known variance.
			// s = 0.5*F.getQuick(r,n-1)+(r+1)*v;
			// Log lik when we write term from variance estimate
			var = F.getQuick(r,n-1)/((n>1) ? (double)(n-1) : 1.0);
			// Bayesian hack
			if (var < min_variance) 
			  {
				var = min_variance;
			  }
			// C.F. Yao 1988: r= # changepoints, so
			// s = 0.5 n log(var) + r log(n)
			s = ((var>0) ? 0.5*n*Math.log(var) : 0) + r*v;
			if (s<min)
			  {
				rmin=r;
				min=s;
			  }
			if (verb) System.err.println("BIC: segs= "+(r+1)+" bic= "+s+
										 " loglik= "+((var>0) ? 0.5*n*Math.log(var) : 0)
										 +" var= "+var);
		  }
		if (verb) System.err.println("Suggest "+(rmin+1)+" segments");


		// Step 3. DP backtrack for changepoint boundaries. Use m for Hawkins' tau
		// Backtrack here depends on rmin
		int rMaxModelRange,rMinModelRange,model_indx;

		if (generateAllDPSegmentations)
		  {
			rMaxModelRange = rmin;
			rMinModelRange = 1;
		  }
		else
		  {
			rMaxModelRange = rmin;
			rMinModelRange = rmin;
		  }
		changepts_array = new IntArrayList[rMaxModelRange-rMinModelRange+1];


		for (rmin = rMinModelRange,model_indx=0; rmin<=rMaxModelRange; rmin++,model_indx++)
		  {
			changepts_array[model_indx] = new IntArrayList();
			// useful reference
			changepts = changepts_array[model_indx];
			// flush list if required
			if (rmin >rMinModelRange) changepts.removeAll(changepts);
			m = n-1;
			// DP backtrack
			// Form the chpt list. Add start and end of array
			changepts.add(0);
			for (i=1,r=rmin-1;r>=0;r--,i++)
			  {
				m = H[r+1][m]; // as computed on unit-offset array data[m-1] in **** above
				changepts.add(m-1);  // now saved as zero-offset
			  }
			// .. and add and end of array if needed
			if (!changepts.contains(n-1))	changepts.add(n-1);
			// sort
			changepts.sort();
			changepts.trimToSize();
		  }
	  }
	// clean up
	W = null;
	F = null;
	H = null;

	return rmin;
  }

  // Solves A.x = b
  void linSolve2x2(double[][] A,double[] b,double[] x)
  {
	double det;

	det = (-(A[0][1]*A[1][0]) + A[0][0]*A[1][1]);

	if (Math.abs(det) != 0.0)
	  {
		x[0] = (-(b[1]*A[0][1]) + b[0]*A[1][1])/det;
		x[1] = (b[1]*A[0][0] - b[0]*A[1][0])/det;
	  }
	else
	  {
		if (A[0][0] != 0.0)
		  {
			x[0] = b[0]/A[0][0];
			x[1] = 0.0;
		  }
	  }
  }
  
  /**
	 Sets variable to print to stdout the cascade of model segmentations.
	 Useful in testing mode
  */
  void showCascade(boolean cascade)
  {
	printAllDPSegmentations = cascade;
  }
  
  /** 
	Exact Dynamic Programming method, with
	very approximate decimation method to reduce data sequences to
	manageable size based on filtering small increments. Pairs
	(y_i,y_i+1) are replaced with (y_i+y_{i+1})/2 in a sequence
	computed from the sorted array { |y_{i+1}-y_i| } (i.e. "merging"
	small increments.) This is continued until the array is small
	enough for the dynamicProgrammingMLSegment() methods, which are
	O(n^2) in memory.

	If makeFullHierarchy==true, an array of segmentations list = IntArrayListt[] 
	is produced, each entry being the
	maximum likelihood segmentation for i = 1,2,..list.length segments. 
	The last in the list is the most likely model selected using the BIC criterion.
	If makeFullHierarchy==false, the array has length 1, with list[0] corresponding
	to the BIC-selected ML model.
*/

  public IntArrayList[] dynamicProgrammingMLSegmentWithDecimationIfRequired(boolean makeFullHierarchy)
  {
	int i,j,k,N;
	SortableNumbers[] sN;
	SortableNumbers num;
	double[] olddata;
	double val;

	N = data.length;

	generateAllDPSegmentations = makeFullHierarchy;

	if (prescribedNSegments > N)
	  {
		System.err.println("Too many segments specified for data size "+N);
		return null;
	  }

	if (N>MAX_DATA_SIZE_FOR_ML_METHOD)
	  {
		if (verb) System.err.println("Decimating to size "+MAX_DATA_SIZE_FOR_ML_METHOD);
		usingDecimation = true;

		sN = new SortableNumbers[N-1];
		decimationLookups = new ArrayList<SortableNumbers>(N);
		for (i=0;i<N-1;i++)
		  {
			sN[i] = new SortableNumbers(Math.abs(data[i+1]-data[i]),i);
			decimationLookups.add(i,new SortableNumbers(data[i],i));
			// System.out.println(i+" "+data[i]);
		  }
		// last entry
		decimationLookups.add(N-1,new SortableNumbers(data[N-1],N-1));
		// System.out.println(i+" "+data[i]+"\n");

		// sort by difference size. Leading part of array is then list of best candidates to drop
		// (smallest differences)
		Arrays.sort(sN);
		// sort this leading part of the candidate list then by index
		Arrays.sort(sN, 0, N-MAX_DATA_SIZE_FOR_ML_METHOD, new SortableNumbers.IndexSorter());
	
		// remove smallest entries running backwards in array
		// and replace first element by average of pair
		for (i=N-MAX_DATA_SIZE_FOR_ML_METHOD-1;i>=0;i--)
		  {
			// set value to mean
			val = 0.5*((decimationLookups.get(sN[i].index)).value + (decimationLookups.get(sN[i].index+1)).value);
			(decimationLookups.get(sN[i].index)).value = val;
			// remove second of pair
			decimationLookups.remove(sN[i].index+1);
		  }

		// redefine data
		olddata = data;

		data = new double[MAX_DATA_SIZE_FOR_ML_METHOD];
		int lastIndx=0;
		double lastVal=data[0];
		for (i=0;i<MAX_DATA_SIZE_FOR_ML_METHOD;i++)
		  {
			num = decimationLookups.get(i);
			data[i] = num.value;
			if (false)
			  {
				for(j=lastIndx+1; j<num.index;j++) System.out.println(j+" "+lastVal);
				System.out.println(num.index+" "+data[i]);
				lastIndx = num.index;
				lastVal = num.value;
			  }
		  }
		// System.out.println("");

		// and redo segmentation on reduced data set
		changepts = dynamicProgrammingMLSegment();
		if (changepts==null) return null;
		// then realign segments, for each set of changepoints in the cascade
		for (i=0;i<changepts_array.length;i++)
		  {			
			changepts = changepts_array[i];
			for (j=0;j<changepts.size();j++)
			  {
				k = changepts.get(j);
				if(k>=0)
				  {
					num = decimationLookups.get(k);
					changepts.set(j,num.index);
				  }
				// else case is k=-1. Leave as-is
			  }
		  }
		data = olddata;
	  }
	else
	  {
		usingDecimation=false;
		dynamicProgrammingMLSegment();
	  }
	return changepts_array;
  }
}
