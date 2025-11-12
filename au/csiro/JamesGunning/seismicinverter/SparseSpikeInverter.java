/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT ---------------------------------  

/*
  This is a simple "demo type" threaded linear 
  sparse spike inversion code based on the LARS approach to
  least-squares with L1 regularization. It's roughly
  a hybrid LARS-OLS regression method, at present.

  The reflectivity model is linearized, and traces are independent.
  
  @author James.Gunning@csiro.au
  @version 1.1, 1/7/2012 
*/
package au.csiro.JamesGunning.seismicinverter;


import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.Numerical.FFT;
import au.csiro.JamesGunning.lib.util.*;
import static au.csiro.JamesGunning.lib.util.MathUtil.sqr;

import au.csiro.JamesGunning.lib.graphutil.*;

import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import cern.colt.map.OpenLongObjectHashMap;

import static java.lang.Math.*;
import java.awt.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;


public class SparseSpikeInverter 
{
  static final boolean USE_CHOL_UPDATE = true;
  static final boolean USE_LOOKUP = true;

  // options
  int verbosity=0;

  // IRLS working variables
  boolean useIRLS=false;
  // Use steepest descent
  boolean useSD = false; 
  // ITS method
  boolean useITS=false;

  boolean addSmoothingTerm = true;
  String inline_word=null,xline_word=null;
  ArrayList<float[]>Ytr;
  int[][] XY; // geometry here: {x,y,east-nbr-indx,nth-nbr-indx}
  double sqrt_lambda=1.0, sqrt_mu_smooth=1.0;
  double[] y,e,resid,y_orig;
  int Np,Nd,KDIP;
  long DTT,DTTS;

  ErrorHandler EH;
  float tmin,tmax;

  boolean fromLE= false;
  float est_noise = 0 ;
  float est_RMS_noise = 0;
  boolean show_synth=false;
  int verb=0;
  int count=-1;
  boolean append =false; 
  boolean useBIC =false;
  int N,Nstacks;
  StackAndWavInfo[] sI = null;
  double[][] Xt;
  DenseDoubleMatrix2D sG_A_lookup;

  byte[] free_ind=null;
  boolean do_PN=false;
  double eps_OLS=1.0;

  // Experimental MRF stuff
  boolean doGraphics=false;
  boolean useMRF =false;
  ArrayList<float[]> Slist=null;

  String oFile=null;
  SuWriter suwriter= new SuWriter();
  SuSegyFormat tr_out = new SuSegyFormat();
  FileOutputStream out_strm=null;

  // Thread management
  Vector<SparseSpikeInverter.Job> freeJobList = new Vector();
  int Nthreads=1;
  volatile int i_consumed,i_delivered; 
  volatile int m; // tracks work
  PriorityBlockingQueue<SparseSpikeInverter.Job> PQ;


  SparseSpikeInverter()
  {
	EH = new ErrorHandler();
  }

  public static void main(String[] args) 
  {
	SparseSpikeInverter SSI = new SparseSpikeInverter();

	SSI.parse(args);
	if (!SSI.useIRLS && !SSI.useSD && !SSI.useITS) 
	  SSI.runLARS();
	else
	  SSI.runIRLS_or_SD_or_ITS();
  }


  void parse(String argv[])
  {
	int i;
	String s1;
	String[] opts;

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("SparseSpikeInverter",EH);

	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber); 

	// some decoration to make help more intelligible
	options.addMainCommandlineText(" [ < stdin ]");
	options.addSpielText("Simple sparse spike inversion code based on a hybrid LARS-OLS algorithm and complete "
						 +"linearization of the reflectivity. Multithreading should help performance "
						 +"on SMP/multi-core machines. Mainly a demo code on java multithreading, but "
						 +"may be useful for other purposes. Example (on big-endian machine): %n"
						 +"deliverySS near-stack.su wavelet.su 2,3 --BIC --synth| suxwigb");



   	options.addMandatoryArg(1,"stack1.su,stack2.su,...","stacks for inversion: use \"stdin\" to indicate std input");
   	options.addMandatoryArg(2,"wavelet1.su,wavelet2.su,...","wavelets for inversion");
   	options.addMandatoryArg(3,"tmin,tmax","time range (seconds) of inversion");

	options.addOption("-h","--help",0,null,null,"Show this help");
	options.addOption("-v",null,1,"int","level","Execute with stated level of verbosity");
	options.addOption("-o",null,1,"String","filename.su","Output file other than stdout");
	options.addOption("-n",null,1,"float","amplitude","Supply absolute noise level");
	options.addOption("-nRMS",null,1,"float","amplitude","Supply noise-to-signal ratio (relative to RMS, per trace)");
	options.addOption("--synth",null,0,null,null,"Dump out synthetic model, not reflectivity");
	options.addOption("--count",null,1,"int","number","Do only \"count\" traces");
	options.addOption("--append",null,0,null,null,"Append model to original amplitude data in output "
					  +"(makes e.g. , using --synth, {true,synthetic} SU pairs)");
	options.addOption("--BIC",null,0,null,null,"Use BIC to determine number of spikes, not default discrepancy principle");
	options.addOption("--MRF",null,0,null,null,"Try experimental MRF spike-cleaner (incomplete/untested)");
	options.addOption(null,"--threads",1,"int","level","Use nthreads suggested (default from hardware)");
	options.addOption("--LE",null,0,null,null,"Flip seismic & wavelet data from little endian when reading");
	options.addOption("--IRLS",null,1,"String","inline_word,xline_word,lambda,mu",
					  "Try global iteratively reweighted least squares (IRLS) L1 method, not LARS (e.g. ep,cdp,10,1).");
	options.addOption("--SD",null,1,"String","inline_word,xline_word,lambda,mu",
					  "Try projected steepest-descent L1 method, not LARS (e.g. ep,cdp,10,1).");
	options.addOption("--ITS",null,1,"String","inline_word,xline_word,lambda,mu",
					  "Try global iterative soft thresholding (ITS) L1 method, not LARS (e.g. ep,cdp,10,1).");
	options.addSpielText("The IRLS, SD and ITS options computes the global reflectivity solution %nR = argmin |W*R-S|_2 + 100 mu |Dt.R|_2 + lambda |R|_1 %nwith S,W rescaled by global seismic RMS. Dt is transverse difference operator acting on reflectivities. The smoothness term mu has an extra 100 factor since typical reflectivities are 0.1");

	// leave this to last!
	//	options.addExtraHelpOption("-X","Show help with all options (there are lots more!)");

	options.parseOptions(argv,true);

	if (options.isOptionSet("-h"))
	  {
		options.Usage();
	  }

  	// stacks
	s1 = options.getMandatoryArg(1);
	opts = s1.split(",");
	Nstacks = opts.length;
	sI = new StackAndWavInfo[Nstacks];
	for(i=0;i<Nstacks;i++) 
	  {
		sI[i] = new StackAndWavInfo(this);
		sI[i].sFile = opts[i];
	  }
	// wavelets
	s1 = options.getMandatoryArg(2);
	opts = s1.split(",");
	if (Nstacks != opts.length)
	  {
		EH.die_err("Must have matching number of stacks and wavelets: "+
				   options.getMandatoryArg(1)+ " : "+options.getMandatoryArg(2));
	  }
	for(i=0;i<Nstacks;i++) sI[i].wFile = opts[i];

	if (Nstacks!=1) EH.die_err("Not coded for Nstacks>1 yet");

	// Time gate
	s1 = options.getMandatoryArg(3);

	opts = s1.split(",");
	if (opts.length!=2) EH.die_err("Must supply tmin,tmax");
	
	try {tmin = Float.parseFloat(opts[0]);}
	catch (NumberFormatException e) {EH.die_err("tmin not a number");}
	tmin = Math.max(tmin, 0);
	
	try {tmax = Float.parseFloat(opts[1]);}
	catch (NumberFormatException e) {EH.die_err("tmax not a number");}
	if (tmax<=tmin) EH.die_err("tmax <= tmin ??");
  
	
	if (options.isOptionSet("-v")) verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();
	if (options.isOptionSet("-o")) oFile = ((String)options.getOptionValues("-o")[0]);

	if (options.isOptionSet("-n")) est_noise = ((Float)options.getOptionValues("-n")[0]).floatValue();
	if (options.isOptionSet("-nRMS")) 
	  {
		est_RMS_noise = ((Float)options.getOptionValues("-nRMS")[0]).floatValue();
		if (est_RMS_noise <0.01 || est_RMS_noise>10)
		  {
			EH.die_err("-nRMS "+est_RMS_noise+" should be in range : 0.01 < est_rms_noise < 10");	
		  }
	  }
	if (options.isOptionSet("--synth")) show_synth = true;
	if (options.isOptionSet("--LE")) fromLE = true;
	if (options.isOptionSet("--BIC")) useBIC = true;
	if (options.isOptionSet("--MRF"))
	  {
		useMRF = true;
	  }
	if (options.isOptionSet("--append")) append = true;

	if (options.isOptionSet("--threads"))
	  {
		Nthreads = ((Integer)options.getOptionValues("--threads")[0]).intValue();
		if (Nthreads <1) EH.die_err("Can't use nthreads < 1");
	  }
	if (options.isOptionSet("--count"))
	  {
		count = ((Integer)options.getOptionValues("--count")[0]).intValue();
		if (count < 1) EH.die_err("Can't use count < 1");
	  }
	s1=null;
	if (options.isOptionSet("--IRLS")) 
	  {
		s1 = "--IRLS";
		useIRLS=true;
	  }
	else if (options.isOptionSet("--SD")) 
	  {
		s1 = "--SD";
		useSD=true;
	  }
	else if (options.isOptionSet("--ITS"))
	  {
		s1 = "--ITS";
		useITS=true;
	  }
	if (s1 != null)
	  {
		// parse --IRLS, --ITS options
		opts = options.getOptionString(s1,0).split(",");
		if (opts.length!=4) EH.die_err("Check --IRLS, --SD, --ITS usage!");
		inline_word = opts[0];
		xline_word = opts[1];
		try
		  {
			sqrt_lambda = Double.parseDouble(opts[2]);
			sqrt_mu_smooth = Double.parseDouble(opts[3]);
			if (sqrt_lambda<=0 ||sqrt_mu_smooth<0) throw new NumberFormatException("");
			sqrt_lambda = sqrt(sqrt_lambda);
			sqrt_mu_smooth = 10.0*sqrt(sqrt_mu_smooth);
		  }
		catch (NumberFormatException e)
		  {
			EH.die_err("--IRLS, --ITS usage needs numbers lamba,mu>0");
		  }
	  }
  }

  /** 
	  Global L1 inversion scheme, minimising
	  objective

	  chisq = ||S-w*R||^2 + lambda |R||^1 + 100 mu ||D.R||^2 

	  using either
	  a) IRLS with LSQR as the main CG engine for each OLS reweighted problem. D=transverse difference operator
	  b) a projected steepest descent method
	  c) iterative thresholding
  */
  

  void runIRLS_or_SD_or_ITS()
  {
	int i,i_outer,j,k,kk,i_s=0,Ntr;
	double var,chisq,chisqO,u=0,tmp;
	double eps,sqrteps=0;
	double[] R=null,S=null;
	byte[] hdr_bytes=null;

	if (Nstacks!=1) EH.die_err("IRLS/SD/ITS not yet coded properly for multistack: use 1 stack");
	// read wavelets, open seismic
	for(i=0;i<Nstacks;i++)
	  {
		sI[i].readWaveletSetUpX();
		sI[i].openSeismic();
	  }
	// 
	try 
	  {
		if (oFile==null)
		  {
			out_strm = (FileOutputStream) new FileOutputStream(FileDescriptor.out);
		  }
		else
		  {
			out_strm = (FileOutputStream) new FileOutputStream(oFile);
		  }
	  }
	catch (Exception e) 
	  {
		throw new RuntimeException("couldn't open output stream");
	  }

	// List of traces
	ArrayList<float[]>Ytr = new ArrayList<float[]>();
	ArrayList<byte[]>headers = new ArrayList<byte[]>();
	// XY lookups
	ArrayList<int[]>XYlist = new ArrayList<int[]>();
	OpenLongObjectHashMap XYmap = new OpenLongObjectHashMap();


	// read & store full data
	int i_tr,X,Y,Xlast,Ylast,DX,DY;
	i_tr=Xlast=Ylast=0;
	DX=DY=Integer.MAX_VALUE;
	FULLPASS:
	for(i=0,var=0.0,Np=0;;i++)
	  {		
		for (j=0;j<Nstacks;j++) 
		  {						
			sI[j].getSeismic();			
			if (sI[j].stackExhausted) break FULLPASS;
			n = (int)(1+(tmax-tmin)/(sI[j].tr.get_short("dt")*1e-6));
			if (j==0 && i==0) 
			  {
				N=n;
				if (verbosity>0) EH.warn("Vertical system size (N samples in window)="+N);
				i_s = sI[j].initWindowA();
				tr_out.set("ns",sI[j].tr.get_short("ns"));
				tr_out.set("dt",sI[j].tr.get_short("dt"));				
			  }
			else if (n!=N)
			  {
				EH.die_err("dt appears to vary between stacks??? resample to same!");
			  }
			float[] data = new float[N];
			for (k=0,kk=sI[j].i_s; k<N; k++,kk++) 
			  {
				data[k] = sI[j].tr.get_data_value(kk);
				var += sqr(data[k]);
			  }
			Ytr.add(data);
			// take copy of hdrs, set neighbourhood etc, from stack 0
			if (j==0) 
			  {
				hdr_bytes= new byte[SuSeismicConstants.HDRBYTES];
				sI[0].tr.copyHeader(hdr_bytes);
				headers.add(hdr_bytes);
				X = sI[j].tr.get_int(inline_word);
				Y = sI[j].tr.get_int(xline_word);
				if (i_tr != 0) 
				  {
					if (Xlast!=X) DX = min(abs(X - Xlast),DX);
					if (Ylast!=Y) DY = min(abs(Y - Ylast),DY);
				  }
				XYlist.add(new int[] {X,Y,-1,-1});   // will store x,y,e-indx,n-indx
				long key = (((long)X << 32) | (Y & 0xFFFFFFFFL)); // pack x,y ints into a long
				XYmap.put(key,new Integer(i_tr));
				Xlast=X; 
				Ylast=Y;
				i_tr++;
			  }  
			Nd+=n;
		  }
	  }
	if (DX==Integer.MAX_VALUE) DX=-1;  // plane
	if (DY==Integer.MAX_VALUE) DY=-1;  //plane
	EH.warn("Geometry: inferred DX=%s DY=%s\n",((DX>0)?DX:"none"),(DY>0)?DY:"none");
	var = sqrt(var/Nd);  // we'll work with scaled data. Use sample stddev as scale..
	Ytr.trimToSize();
	XYlist.trimToSize();
	// sort out adjacent pairs in spatial pattern
	Ntr = Ytr.size();
	for(i=0,k=0;i<Ntr;i++)
	  {
		int[] xy = XYlist.get(i);
		// East: nbr of target pt (x+dx,y)
		long key = (((long)(xy[0]+DX) << 32) | (xy[1] & 0xFFFFFFFFL));
		Integer I = (Integer)XYmap.get(key);
		// save trace number
		if (I!=null) xy[2]=I.intValue(); // trace number of east neighbour
		// North: nbr of target pt (x,y+dy)
		key = (((long)(xy[0]) << 32) | ((xy[1]+DY) & 0xFFFFFFFFL));
		I = (Integer)XYmap.get(key);
		if (I!=null) xy[3]=I.intValue(); // trace number of north neighbour
	  }
	XYmap=null; // CG

	// global array useful for smoothing operator
	XY = (int[][]) XYlist.toArray(new int[XYlist.size()][4]);


	Np = Nd;  // for now. single stack


	y = new double[Nd+Np + (addSmoothingTerm ?Np : 0)];
	for(i=0,k=0;i<Ntr;i++)
	  {
		float[] _y = Ytr.get(i);
		for(j=0;j<n;j++,k++) y[k] =(double)_y[j]/var;
	  }

	// Note: last Np-sized block of RHS y is 0  !

	// pre scale wavelets here too
	for(i=0;i<Nstacks;i++)
	  {
		for(k=0;k<sI[i].Xw.length;k++) sI[i].Xw[k]/=var;
		sI[i].setUpAdjoint();
	  }
	// estimate oversampling discount rate based on first wavelet
	double delta_t_indep_samples = getDiscountedSamplingRate(sI[0].Xw);

	double[] x = new double[Np]; // unknown reflectivities
	e = new double[Np];          // scaling errors

	double EPS_EM=0.0001; // outer loop of IRLS reweighting is an EM algorithm. 
	KDIP=1;  // allowable dip-sample-shift in smoothness constraint

	StackAndWavInfo _SI = sI[0];
	int Nc = _SI.n_coda - _SI.n_precursor+1;
	if (useIRLS || useSD)
	  {
		double chisqReg,scaling_noise=1.0,lam = sqr(sqrt_lambda);
		if (est_RMS_noise > 0) 
		  {
			scaling_noise = est_RMS_noise;
		  }
		// IRLS method. Sequence of OLS solves
		S = new double[n];
		R = new double[n];

		// parameters

		eps=0.0001;  // Huber-like truncation of residuals of reflectivity: resid_k = max(|R_k|,eps)
		sqrteps=sqrt(eps);
		Arrays.fill(e,1.0);

		if (useIRLS)
		  {
			// IRLS method
			// This chunk correct and debugged.
			// 
			// Note lambda sum_i |m_i| term goes to (lambda/2) sum_i (m_{i,old}^2 + m_i^2/|m_{i,old}|) under IRLS
			// Bailout when successive chisq values match at this relative level
			EH.warn("IRLS method: CG system size = "+Np);
			// use truncated Newton ? 
			useTruncatedNewtonMethod(true,0.3);

			int N_CG_MAX=1000;   // absolute ceiling on #CG its
			int N_CG_MIN = 10;   // approximate limit for early on
			int N_CG = 5;        // 10? N_CG_MIN;
		
			DTT=DTTS=0;
			for(i_outer=0,chisqO=1.0; i_outer<100; i_outer++)
			  {
				// OLS solve
				CG_solve(i_outer!=0,y,x,N_CG);
				// get residuals
				long t0=System.nanoTime();
				for(k=0,chisq=0; k<Nd; k+=n)
				  {
					// loop over traces
					System.arraycopy(x,k,R,0,n);  // copy trace from x into wkspace R
					FFT.convolve_kernel(Nc, _SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
					for(j=0,kk=k;j<n;j++,kk++) chisq += sqr(S[j]-y[kk]);
				  }
				DTT += System.nanoTime()-t0;
				u = chisq; // L2 data misfit bit for diagnostics below
				// L1 bit
			    // for (j=0,chisqReg=0;j<Np;j++) chisqReg += lam * abs(x[j]);
				// IRLS equivalent - see sec 2.2 Schmidt etal, "Fast Optimization Methods for L1..."
				for (j=0,chisqReg=0;j<Np;j++) chisqReg += sqr(x[j]/e[j]) + sqr(e[j]);
				chisqReg *= lam/2;

				// Add 3rd term smoothing energy
				tmp = getSmoothingEnergy(x);
				chisq = u + chisqReg + tmp;
				EH.warn("IRLS: chisq_MS=%7g chisq_data=%7g chisq_L1=%7g chisq_sm=%7g chisqTot=%7g%n",
						chisq/(scaling_noise*Nd),u,chisqReg,tmp,chisq);

				if ((i_outer>0) && (abs(chisq-chisqO)/chisqO < EPS_EM)) 
				  {
					if (N_CG >= N_CG_MAX) break; // infallible success ;-)
					else
					  {
						N_CG *= 2;
						EH.warn("Nmax Rescale: "+N_CG);
					  }
				  }
				// estimate coefficient residuals: sort of "Huberise"
				for (j=0;j<Np;j++) 
				  {
					if ((tmp=abs(x[j])) < eps) e[j] = sqrteps;
					else e[j] = sqrt(tmp);
				  }  
				chisqO=chisq;			
			  }
		  } // IRLS
		else
		  {
			// projected Steepest descent, or projected Newton (Newton not finished)
			double EPSG=0.0005;      // close-to-zero epsilon to detect non-free variables
			double ALPHA_MIN=1e-5;   // smallest backtrack epsilon
			double EPS_SD = 1e-7;    // outer termination fractional change in chisq, SD method
			double MAX_STEP=0.05;    // uppermost step-length limit in SD method (for reflectivities)
			double EPS_SLOPE=1e-4;   // slope condition in SD line-search
			double alpha=0,s;
			double[] grad=null;
			// for Brent line-search
			SDinfo sdinfo=null;

			// workspace

			// free-variable 0,1 indicators 
			free_ind = new byte[Np]; // unknown reflectivities
			Arrays.fill(free_ind,(byte)1);
			double[] xold = (double[])x.clone();
			
			if (useSD) 
			  {
				// wkspace for proj-SD
				grad = (double[])x.clone();
				sdinfo = new SDinfo(Nc,_SI.n_precursor,R,x,y,xold,grad,S,_SI.Xw,lam);
			  }
			else
			  {
				do_PN=true;  // Switches on correct behavior of AV_product
				// wkspace for PN
				y_orig=(double[])y.clone(); // copy of true scaled amplitudes
			  }
			OUTER:
			for(i_outer=0,chisqO=Double.MAX_VALUE; i_outer<10000; i_outer++)
			  {
				if (useSD)
				  {
					// Projected Steepest Descent algorithm
					// get gradient at xold, set free variables, also max|grad| over free vars

					// Compute gradient part from 3rd term: grad = 2*Ds^T.Ds.xold
					Arrays.fill(grad,0);
					doSmoothingMult(false,0,xold,x); // x as workspace
					doSmoothingMult(true,0,x,grad); // this is additive, hence grad=0 above

					// now grad from 1st and 2nd term
					double xGM=0,free_grad_max=0;
					for(k=0,i=0; k<Nd; k+=n)
					  {
						// loop over traces
						System.arraycopy(xold,k,R,0,n);  // copy trace from xold into wkspace R
						FFT.convolve_kernel(Nc, _SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
						for(j=0,kk=k;j<n;j++,kk++) R[j] = y[kk] - S[j]; // R[] has residuals y-Xm chunk
						// adjoint mult, so S[] has chunks of X^T.(y-Xm)			
						FFT.convolve_kernel(Nc,-_SI.n_coda,_SI.Xwt,n,0,R,n,0,S);
						
						for(j=0,kk=k;j<n;j++,kk++) 
						  {
							u = -2*S[j] +2*grad[kk]; // smooth gradient from first term + 3rd term here
							
							// set free vars
							free_ind[kk]=(byte)1;
							if (abs(xold[kk]) < EPSG)
							  {
								// test for "V" profile indicating trapping
								if (lam > abs(u))
								  {
									// mark as not-free
									free_ind[kk]=0; // free_ind constraint 							  
									i++;
								  }
							  }
							if (xold[kk] == 0.0)
							  {
								// Won't matter for non-free components
								// for active guys at exactly zero, set to larger gradient
								// if(abs(u+lam) > abs(u-lam)) grad[kk] = u+lam;
								if (u>0) grad[kk] = u-lam;  // we're heading -ve. Use left-side grad
								else grad[kk]=u+lam;        // conversely
							  }
							else
							  {
								grad[kk] = u + ((xold[kk]>=0) ? lam : -lam);   // piece to be added
							  }
							if (free_ind[kk]!=0) 
							  {
								if (free_grad_max < abs(grad[kk]))
								  {
									free_grad_max = abs(grad[kk]);
									xGM=xold[kk];
									i_s=kk;
								  }
							  }
						  }
					  }

					if (i_outer%10==0) EH.warn("It: "+i_outer+" Nfree= "+(Np-i)+" Max free grad="+free_grad_max);

				  
					// scale gradient step x := (xold - grad) to maximum of MAX_STEP
					// zero grad components for non-free components
					tmp = MAX_STEP/free_grad_max;
					for(k=0; k<Nd; k++) 
					  {
						if (free_ind[k]== 1) grad[k] *= tmp;
						else grad[k] = 0;
					  }
					
					// linesearch with projection: components crossing 0 get parked at 0 en-route  
					/* 
					   double alpha_ratio=0.9;  // backtrack ratio decrease in SD methods
					   for(alpha=1; alpha>ALPHA_MIN; alpha*=alpha_ratio)
					  {
						chisq = chisqSD(alpha,sdinfo);

						if (chisq < chisqO) // FIX
						  {
							// take copy of new point & chisq
							System.arraycopy(x,0,xold,0,Np);
							EH.warn("SD: alpha= %5f, chisq_data=%7g chisq_L1=%7g chisq_sm=%7g chisqTot=%7g%n",
									alpha,sdinfo.E_data,sdinfo.E_L1,sdinfo.E_sm,chisq);
							if ((i_outer>0) && (abs(chisq-chisqO)/chisqO < EPS_SD)) 
							  {
								chisqO=chisq;							
								//EH.warn("Termination: Max free ref.grad=%6g at x=%6g, g=%6g%n,",free_grad_max,xGM,grad[i_s]*free_grad_max/MAX_STEP);	
								//for(k=0,tmp=0;k<Np;k++) if (free_ind[k]!=0)EH.printf("GG %f\n",abs(grad[k])*free_grad_max/MAX_STEP);
								//System.exit(0);
								break OUTER; // infallible success ;-)
							  }
							chisqO=chisq;
							break;    // backtrack loop
						  }
					  }
*/
					// Minpack-like backtrack (roughly, Schnabel's lnsrch)
					double slp,pfpls=0.0,plmbda = 0.0,alpha_new;
					int terminate=0;
					// at this point, grad <- tmp*grad, and the proposed full step is Dx = -grad (non-zero components)
					for(k=0,slp=0; k<Nd; k++) slp += sqr(grad[k]);					 
					slp=-slp/tmp;
					// here, slp = Dx.{grad(f)}, where Dx is the full step for alpha=1, which is a scaled -grad step
					for(alpha=1.0;;)
					  {
						chisq = chisqSD(alpha,sdinfo);
						if (chisq <= chisqO + slp*EPS_SLOPE*alpha) terminate = 2;
						else
						  {
							if (alpha==1.0)
							  {
								// first backtrack: quadratic fit
								alpha_new = 2.0*(chisq - chisqO - slp);
								if (alpha_new!=0.0) 
								  {
									alpha_new = -slp/alpha_new; 
								  }
								else
								  {
									alpha_new=0;
									terminate=4;
								  }
							  } 
							else 
							  {								  
								// all subsequent backtracks: cubic fit
								double t1,t2,t3,a3,b,disc;
								t1 = chisq - chisqO - alpha*slp;
								t2 = pfpls - chisqO - plmbda*slp;
								t3 = 1.0/(alpha - plmbda);
								a3 = 3.0*t3*(t1/(alpha*alpha) - t2/(plmbda*plmbda));
								b = t3*(t2*alpha/(plmbda*plmbda) - t1*plmbda/(alpha*alpha));
								disc = b*b - a3*slp;									
								if (disc > b*b) 
								  {										
									// only one positive critical point, must be minimum			
									alpha_new = (-b + signum(a3) * sqrt(disc))/a3;
								  } 
								else 
								  {										
									// both critical points positive, first is minimum
									alpha_new = (-b - signum(a3) * sqrt(disc))/a3;	
								  }
								if (alpha_new > 0.5*alpha) alpha_new = 0.5*alpha;                  
							  }
							plmbda = alpha;
							pfpls = chisq;
							// additional JG safeguards
							if (alpha_new<=0 || alpha_new>1.0)
							  {
								EH.warn("Termination: wierd line-search result, alpha=%6g%n,",alpha_new);	
								terminate=3;
							  }
							// below, as per uncmin/minpack
							else if (alpha_new < 0.1*alpha) 
							  {
								alpha *= 0.1;
							  } 
							else 
							  {
								alpha = alpha_new;
							  }
						  }		
						if (terminate != 0 )
						  {
							// take copy of new point & chisq
							System.arraycopy(x,0,xold,0,Np);
							EH.warn("SD: alpha= %5f, chisq_data=%7g chisq_L1=%7g chisq_sm=%7g chisqTot=%7g%n",
									alpha,sdinfo.E_data,sdinfo.E_L1,sdinfo.E_sm,chisq);
							if ((i_outer>0) && (abs(chisq-chisqO)/chisqO < EPS_SD)) 
							  {
								chisqO=chisq;							
								//EH.warn("Termination: Max free ref.grad=%6g at x=%6g, g=%6g%n,",free_grad_max,xGM,grad[i_s]*free_grad_max/MAX_STEP);	
								//for(k=0,tmp=0;k<Np;k++) if (free_ind[k]!=0)EH.printf("GG %f\n",abs(grad[k])*free_grad_max/MAX_STEP);
								break OUTER; // infallible success ;-)
							  }
							chisqO=chisq;
							break;    // backtrack loop
						  }				  
						if (alpha<ALPHA_MIN)
						  {
							EH.warn("SD: line-search degenerate break!");
							break OUTER; 
						  }
					  }
				  }		
				else
				  {
					// Start of Projected Newton method.
					// FIXME This guff is completely wrong at present (=> disabled)....
					eps_OLS=.01;
					// set RHS with right gradient
					// block 1
					for(k=0; k<Nd; k+=n)
					  {
						System.arraycopy(xold,k,R,0,n);  // copy trace from xold into wkspace R
						for(j=0,kk=k;j<n;j++,kk++) if (free_ind[kk]==(byte)0) R[j]=0;
						FFT.convolve_kernel(Nc, _SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
						for(j=0,kk=k;j<n;j++,kk++) y[kk] = y_orig[kk] - S[j]; // residuals y-Xm chunk, this trace
					  }
					// block 2
					tmp = -lam/(2*eps_OLS);
					for(i=Nd,k=0;i<(Nd+Np);i++,k++)
					  {
						if (abs(u=xold[k]) > EPSG) y[i] = tmp * free_ind[k]* ((u>0) ?  1 : -1);
						else y[i]=0;
					  }
					// block 3
					if (addSmoothingTerm)
					  {					
						// FIXME HERE!! doSmoothingMult(true,0,res,x);
						for(i=Nd+Np; i<(Nd+2*Np); i++) y[i]=0;
					  }
					// OLS solve for Newton updates, current free variables
					int N_CG_PN=500;
					CG_solve(false,y,x,N_CG_PN);
					// REST OF PN algorithm wrong/corrupt! FIXME
				  }
			  } // PN or SD
			// Profile check....
			// for(alpha=-.1;alpha<=.1;alpha+=.005) EH.printf("XX %12f %12f%n",alpha,chisqSD(alpha,sdinfo));
		  } // SD/PN
	  } // IRLS or SD/PN
	else
	  {
		// ITS = iterative soft thresholding/Landweber method
		// see, e.g. "On the performance of algorithms for the 
		// minimization of l_1 -penalized functionals", Loris, 2009, Eqn (8)
		double[] res = new double[Nd+Np];          // residuals
		double[] dx = new double[Np];               // IST model changes
		double lambda,eMax=0;

		if (false)
		  {
			// junk: test of correctness of sparse A.x, A^T.x operations
			n=6;
			Nd=4*n;
			Np=Nd;
			res = new double[Nd+Np];          // residuals
			x = new double[Np]; 
			DenseDoubleMatrix2D A = new DenseDoubleMatrix2D(Nd,Np);

			for(i=0;i<Np;i++)
			  {
				Arrays.fill(x,0);
				x[i]=1;
				doSmoothingMult(false,0,x,res);
				for(j=0;j<Nd;j++) A.s(j,i,res[j]);
			  }
			EH.println("A="+A);
			A.assign(0);
			for(i=0;i<Nd;i++)
			  {
				Arrays.fill(res,0);
				res[i]=1;
				Arrays.fill(x,0);
				doSmoothingMult(true,0,res,x);
				for(j=0;j<Np;j++) A.s(i,j,x[j]);
			  }
			EH.println("A="+A.viewDice());
			System.exit(0);
		  }
		// (1) Estimate largest eigenvalue of A^T.A by Rayleigh quotient
		for(j=0,u=0;j<Np;j++) {x[j]=Math.random(); u+=x[j];}
		u = sqrt(u);
		for(j=0;j<Np;j++) x[j]/=u;
		for(i_outer=0; i_outer<10; i_outer++)
		  {
			// form x' = A^T.A.x
			formConvolution(false,x,res);
			doSmoothingMult(false,Nd,x,res);

			formConvolution(true,res,dx);
			doSmoothingMult(true,Nd,res,dx);
			for(j=0,u=0;j<Np;j++) u += sqr(dx[j]);
			u=sqrt(u);
			for(j=0,eMax=0;j<Np;j++) 
			  {	   
				eMax += dx[j]*x[j];
				dx[j]/=u;
			  }
			if (verbosity>0) EH.printf("eMax=%g%n",eMax);
			for(j=0;j<Np;j++) x[j]=dx[j];
		  }
		eMax *=2.0;   // factor of safety
		if (verbosity>0) EH.printf("Taking (approx 2 x max{eigenvalue}) eMax=%g%n",eMax);
		// below sqrt_lambda^2 = 2 lambda in Loris
		lambda = 0.5*sqr(sqrt_lambda)/eMax;


		// (2) IST method
		Arrays.fill(x,0);
		for(i_outer=0,chisqO=1.0; i_outer<2000; i_outer++)
		  {
			// res = y - A.x
			Arrays.fill(res,0);
			formConvolution(false,x,res);
			doSmoothingMult(false,Nd,x,res);
			// data residuals
			for(j=0,chisq=0;j<Nd;j++) 
			  {
				u = y[j]-res[j];
				res[j] = u/eMax;
				chisq += u*u;
			  }
			// smoothness residuals (y_data=0)
			for(j=Nd;j<Nd+Np;j++) 
			  {
				u = -res[j];
				res[j] = u/eMax;
			  }
			if (verbosity>0) EH.printf("Chisq= %g%n",chisq);
			Arrays.fill(dx,0);
			formConvolution(true,res,dx);
			doSmoothingMult(true,Nd,res,dx);
			for(j=0;j<Np;j++)
			  {
				u = x[j]+dx[j];
				if (u >= lambda) x[j] = u-lambda;
				else if (u <= -lambda) x[j] = u+lambda;
				else x[j]=0;
			  }
			if(i_outer>0 && ((abs(chisq-chisqO)/chisqO)<EPS_EM)) break;
			chisqO=chisq;
		  } // i_outer
	  } // ITS

	// IO

	// write result
	if (show_synth)
	  {
		R = new double[n];
		S = new double[n];
	  }

	// IO
	int tracl=0;
	for(i=0;i<Nd;i+=n,tracl++)
	  {
		// set correct hdr geom keywords etc
		hdr_bytes = headers.get(tracl);
		tr_out.setHeader(hdr_bytes);
		// but overwrite these
		tr_out.set("ns",(short)n);
		tr_out.set("f1",tmin);
		tr_out.set("delrt",(short)rint(tmin*1000));

		// tr_out.set("tracl",tracl);
		if (append)
		  {
			// rescale to real units
			kk=i;
			for(k=0;k<n;k++,kk++) tr_out.set_data_value((float)(var*y[kk]),k);
			// write out
			suwriter.fputtr(out_strm,tr_out);
		  }
		if (!show_synth)
		  {
			kk=i;
			for(k=0;k<n;k++,kk++) tr_out.set_data_value((float)x[kk],k);
			// write out
			suwriter.fputtr(out_strm,tr_out);
		  }
		else
		  {
			System.arraycopy(x,i,R,0,n);  // copy trace from x into wkspace R
			FFT.convolve_kernel(Nc,_SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
			// rescale to real units
			for(j=0;j<n;j++) tr_out.set_data_value((float)(var*S[j]),j);
			// write out
			suwriter.fputtr(out_strm,tr_out);
		  }
	  }
	  

	// finish up
	for(i=0;i<Nstacks;i++) sI[i].closeSeismic();
	suwriter.flush(out_strm);
	try {out_strm.close();}
	catch (IOException err) 
	  {
		throw new RuntimeException("couldn't close SU out stream");
	  }
  }

  private final class SDinfo
  {
	// in
	int Nc,n_precursor;
	double[] R,x,y,xold,grad,S,Xw;
	double lam;
	// out
	double E_data,E_L1,E_sm;

	public SDinfo(int _Nc,int _n_precursor,double[] _R,double[] _x,double[] _y,double[] _xold,double[] _grad,
						double[] _S,double[] _Xw,double _lam)
	{
	  Nc = _Nc;
	  n_precursor = _n_precursor;
	  R=_R;
	  x=_x;
	  y=_y;
	  xold=_xold;
	  grad=_grad;
	  S=_S;
	  Xw=_Xw;
	  lam=_lam;
	}
  }

  private final double chisqSD(double alpha,SDinfo info)
  {
	int j,k,kk;
	double s,u,tmp,chisq;
	int Nc = info.Nc;
	int n_precursor = info.n_precursor;
	double[] R = info.R;
	double[] x = info.x;
	double[] y = info.y;
	double[] S = info.S;
	double[] xold = info.xold;
	double[] grad = info.grad;
	double[] Xw = info.Xw;

	// loop over traces
	for(k=0,tmp=0,s=0; k<Nd; k+=n) // s = |m|_1 norm, tmp = |y-Xm|^2 norm
	  {
		for(j=0,kk=k;j<n;j++,kk++) 
		  {
			// backtrack
			// new point
			u  = xold[kk] - alpha*grad[kk];
			// project
			if (xold[kk]>0 && u<0) u=0;      // parking
			else if (xold[kk]<0 && u>0) u=0; // parking
			x[kk] = R[j]=u; // projected, backtracked new point
		  }
		FFT.convolve_kernel(Nc,n_precursor,Xw,n,0,R,n,0,S);
		for(j=0,kk=k;j<n;j++,kk++) 
		  {
			tmp += sqr(y[kk] - S[j]);
			s += abs(x[kk]);
		  }
	  }
	info.E_data = tmp; 	
	info.E_L1 = s * info.lam;   // L1 penalty
	// 3rd term smoothing energy
	info.E_sm = getSmoothingEnergy(x);

	chisq = info.E_data + info.E_L1 + info.E_sm;  // 3 terms
	return chisq;
  }


  double getDiscountedSamplingRate(double[] signal)
  {
	// information content guess for wavelet 0
	int ng,Nw,i_dft,i_max=0,i;
	double MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER=1.124;
	double SAMPLE_TIME_MAGIC_CONSTANT = 0.253;
	double delta_t_indep_samples;

	Nw = signal.length;
	ng=0;
	while ((Nw >>= 1) > 0 ) ng++;
	// In general Nw not a power of 2, so find the right n = 2^ng  bigger than Nw
	if (signal.length != (1<<ng))
	  {
		// Nw not a power of 2
		Nw =(1<<(ng+1));
	  }
	else Nw= signal.length; // it was already a power of 2
	// now pad out W
	double[] padded_wavelet = new double[Nw];
	// copy wavelet
	System.arraycopy(signal,0,padded_wavelet,0,signal.length);
	FFT.realftZO(padded_wavelet,Nw,1);
	// set the peak frequency from the fft-weighted index (first moment)
	// this is the "mean" frequency. If the spectrum is hump shaped, this will be a good guess	
	double delta_t=1,i_peak,f_peak,spec,max=0,num=0,den=0, freq_den=((double)n * delta_t);

	// OLD for (i=2,i_dft=1;i<=Nw/2;i+=2,i_dft++)
	for (i=2,i_dft=1;i<=Nw-2;i+=2,i_dft++)
	  {
		spec = sqr(padded_wavelet[i])+sqr(padded_wavelet[i+1]);
		if (spec>=max)
		  {
			max = spec;
			i_max=i_dft;
		  }
		num+=spec*i_dft;  // first moment
		den+=spec;    // normalisation
	  }
	i_peak = (num/den)/MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER;
	//EH.warn("Wavelet mean-peak detected at effective index %4g",i_peak/2.0);
		
	f_peak = i_peak/((double)Nw * delta_t);
	//EH.warn("Wavelet peak frequency=%3g Hz",f_peak);

	// Set error-sampling time of all synthetics
	delta_t_indep_samples = SAMPLE_TIME_MAGIC_CONSTANT/f_peak;
	return delta_t_indep_samples;
  }

  /*
	Main IRLS conjugate gradient product routine
	returns a = A.v or a=A^T.v 
	System matrix is blocks like
A = A1   <- convolution
    A2   <- IRLS matrix for current l2 norm of reflectivities R_k sqrt_lambda diag (e_kk)
    A3   <- Smoothing matrix for difference terms R_k - (sum_j R_{k,j})/norm on adjacent trace j~k.

  A^T = (A1^T| A2^T| A3^T)

  A3^T (case KDIP=1):
0
   1 
      1

0 -.33
0 -.33 -.33
0 -.33 -.33 -.33
0  0   -.33 -.33 ..
*/
  private final void Av_product(boolean transp,double[] v,double[] a)
  {
	int i,k,j,kk;
	double u,norm,tmp;

	if (!do_PN)
	  {
		tmp = sqrt_lambda/sqrt(2);
		if (!transp)
		  {
			// top-block: convolution for all traces		
			formConvolution(false,v,a);

			// next block: reweighted prior
			for (k=0,kk=Nd; k<Np; k++,kk++)
			  {
				a[kk] = tmp * v[k]/e[k];
			  }
			// smoothing block (2D x-section)
			doSmoothingMult(false,2*Nd,v,a);
		  }
		else
		  {
			// left-block: reverse convolution
			formConvolution(true,v,a);
			// block (1,2): add reweighted prior
			for(k=0,kk=Nd; k<Np; k++,kk++)
			  {
				a[k] += tmp * v[kk]/e[k];
			  }
			// smoothing block (2D x-section)
			doSmoothingMult(true,2*Nd,v,a);
		  }
	  }
	else
	  {
		u = eps_OLS;
		if (!transp)
		  {
			// top-block: convolution for all traces		
			formConvolution(false,v,a);
			// next block: limiting diagonal
			for (k=0,kk=Nd; k<Np; k++,kk++)
			  {
				a[kk] = u * v[k];
			  }
			// smoothing block (2D x-section)
			doSmoothingMult(false,2*Nd,v,a);
		  }
		else
		  {
			// left-block: reverse convolution
			formConvolution(true,v,a);
			// block (1,2): add limiting diagonal
			for(k=0,kk=Nd; k<Np; k++,kk++)
			  {
				a[k] += u * v[kk];
			  }
			// smoothing block (2D x-section)
			doSmoothingMult(true,2*Nd,v,a);
		  }
	  }
  }


  final void formConvolutionChunk(boolean transp, int kmin,int kmax, double[] v, double[] a)
  {
	int k,Nc;
	double[] S = new double[n];
	double[] R = new double[n];

	StackAndWavInfo _SI = sI[0];
	Nc = _SI.n_coda - _SI.n_precursor+1;
	if (!do_PN)
	  {
		if (!transp)
		  {
			for(k=kmin; k<kmax; k+=n)
			  {
				System.arraycopy(v,k,R,0,n);  // copy trace from v into wkspace R
				FFT.convolve_kernel(Nc,_SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
				System.arraycopy(S,0,a,k,n);			
			  }
		  }
		else
		  {
			for(k=kmin; k<kmax; k+=n)
			  {
				System.arraycopy(v,k,R,0,n);
				FFT.convolve_kernel(Nc,-_SI.n_coda,_SI.Xwt,n,0,R,n,0,S);
				System.arraycopy(S,0,a,k,n);
			  }
		  }	
	  }
	else
	  {
		int i,j,ke;
		// same as above, but non-free components actively zeroed
		if (!transp)
		  {
			for(k=kmin; k<kmax; k+=n)
			  {
				System.arraycopy(v,k,R,0,n);  // copy trace from v into wkspace R
				ke=k+n;
				for(i=k,j=0; i<ke; i++,j++) if(free_ind[i]==(byte)0) R[j]=0;
				FFT.convolve_kernel(Nc,_SI.n_precursor,_SI.Xw,n,0,R,n,0,S);
				System.arraycopy(S,0,a,k,n);			
			  }
		  }
		else
		  {
			for(k=kmin; k<kmax; k+=n)
			  {
				System.arraycopy(v,k,R,0,n);
				ke=k+n;
				for(i=k,j=0; i<ke; i++,j++) if(free_ind[i]==(byte)0) R[j]=0;
				FFT.convolve_kernel(Nc,-_SI.n_coda,_SI.Xwt,n,0,R,n,0,S);
				System.arraycopy(S,0,a,k,n);
			  }
		  }	
	  }
  }

  // usually about .1-1ms/trace
  final void formConvolution(boolean transp, double[] v, double[] a)
  {
	int k,Nc;

	long t0 = System.nanoTime();
	if (Nthreads<=1)
	  {
		if (!transp)
		  {
			formConvolutionChunk(transp,0,Nd,v,a);
		  }
		else
		  {
			formConvolutionChunk(transp,0,Np,v,a);
		  }
	  }
	else
	  {
		int off,i_th,kchunk,dk;
		Thread[] threads = new Thread[Nthreads];
		final boolean _transp=transp;
		final double[] _v=v, _a=a;
		if (!transp)
		  {			
			kchunk = (int)((Nd/n)/Nthreads) + 1;
			dk = n*kchunk;
			for (i_th=0, off=0; i_th<Nthreads; i_th++, off+=dk)
			  {
				final int kmin = off;
				final int kmax = min(dk + off,Nd);				
				threads[i_th] = new Thread(new Runnable() 
				  {
					public void run() {
					  formConvolutionChunk(_transp,kmin,kmax,_v,_a);}
				  });
			  }
		  }
		else
		  {
			kchunk = (int)((Np/n)/Nthreads) + 1;
			dk = n*kchunk;
			for (i_th=0, off=0; i_th<Nthreads; i_th++, off+=dk)
			  {
				final int kmin = off;
				final int kmax = min(dk + off,Np);
				threads[i_th] = new Thread(new Runnable() 
				  {
					public void run() {formConvolutionChunk(_transp,kmin,kmax,_v,_a);}
				  });
			  }
		  }
	  for (i_th=0; i_th<Nthreads; i_th++) threads[i_th].start();
	  try 
		{
		  for (i_th=0; i_th<Nthreads; i_th++) threads[i_th].join();
		}
	  catch (InterruptedException ie) 
		{
		  throw new RuntimeException(ie);
		}
	  } // threaded

	DTT += System.nanoTime()-t0;
	//EH.warn("Conv time=%5g, per-trace=%5g (uS)%n",1e-3*DT,(1e-3*DT)/n);
  }
  
  // Not used (obsolete)....
  final void doSmoothingMultOLD(boolean transp, int kstart,double[] v,double[] a)
  {
	if (addSmoothingTerm && (sqrt_mu_smooth>0))
	  {
		int j,je,js,k,kk,ks,ke=n-1-KDIP;
		double u,tmp,norm = 1.0/(2*KDIP+1);
		long t0 = System.nanoTime();
		if (!transp)
		  {
			je=n+KDIP;   // last "adjacent" sample, next trace
			js = n-KDIP;
			for (k=0,kk=kstart; k<(Np-n); k++,kk++,je++,js++)
			  {
				ks = k%n;
				if ((ks<KDIP) || (ks>ke)) {a[kk]=0.0; continue;}
				//je=k+n+KDIP;  // last "adjacent" sample, next trace
				//for(j=k+n-KDIP,u=0.0; j<=je; j++) u+=v[j];
				for(j=js,u=0.0; j<=je; j++) u+=v[j];
				u*=norm;
				a[kk] = sqrt_mu_smooth*(v[k] - u);
			  }
			for(;k<Np;k++,kk++) a[kk]=0.0;			
		  }
		else
		  {
			// smoothing block (2D x-section)
			// columnwise-driven
			tmp = sqrt_mu_smooth * norm;
			je = n+KDIP;
			js = n-KDIP;
			int kV=kstart;
			for(k=0;k<Np-n;k++,js++,je++,kV++)
			  {
				ks = k%n;
				if ((ks<KDIP)  || ks>ke) continue;	
				// nonzeros
				u = v[kV];
				a[k] += sqrt_mu_smooth*u;
				u *= tmp;						
				for(j=js; j<=je;j++) a[j] -= u;
			  }
		  }
		DTTS += System.nanoTime()-t0;
	  }
  }
  
  
  final void doSmoothingMult(boolean transp, int kstart,double[] v,double[] a)
  {
	if (addSmoothingTerm && (sqrt_mu_smooth>0))
	  {
		int j,je,js,k,k_nbr,kk,ks,ke=n-1-KDIP,nnbr;
		int i_tr,i_EN,Ntr=XY.length;
		int[] js_EN={0,0};
		double u=0,s;

		if (!transp)
		  {
			double norm = 1.0/(2*KDIP+1);
			double norm2 = norm/2.0;
			Arrays.fill(a,kstart,kstart+Np,0);			

			for(i_tr=0,k=0; i_tr<Ntr; i_tr++,k+=n)
			  {
				int[] _XY = XY[i_tr];
				// lookups
				js_EN[0] = (_XY[2]>0) ? (n*_XY[2]) : -1;
				js_EN[1] = (_XY[3]>0) ? (n*_XY[3]) : -1;
				for(i_EN=0,nnbr=0; i_EN<2; i_EN++) if (js_EN[i_EN]>=0) nnbr++;
				if (nnbr==0) continue;
				if (nnbr==1) {s=norm;} else {s=norm2;}
				// ks=relative index, this trace. k_nbr = relative indx, nbr trace				
				for(kk=k+KDIP,ks=KDIP,k_nbr=0; ks<=ke; kk++,ks++,k_nbr++)
				  {
					// loop over E then N nbrs					
					for(i_EN=0,u=0.0; i_EN<2; i_EN++)
					  {
						j = js_EN[i_EN];  // trace indx of E or N nbr
						if (j<0) continue; // no nbr 
						js = j+k_nbr;
						je = js+2*KDIP;
						for(j=js; j<=je; j++) u += v[j];
					  }
					u*=s;
					// only here for nnbr=1,2
					a[kk+kstart] = sqrt_mu_smooth*(v[kk] - u);
				  }
			  }
		  }
		else 
		  {
			// A^T version
			double norm = sqrt_mu_smooth/(2*KDIP+1);
			double norm2 = norm/2.0;
			for(i_tr=0,k=0; i_tr<Ntr; i_tr++,k+=n)
			  {
				int[] _XY = XY[i_tr];
				// lookups
				js_EN[0] = (_XY[2]>0) ? (n*_XY[2]) : -1;
				js_EN[1] = (_XY[3]>0) ? (n*_XY[3]) : -1;
				for(i_EN=0,nnbr=0; i_EN<2; i_EN++) if (js_EN[i_EN]>=0) nnbr++;
				if (nnbr==0) continue;
				if (nnbr==1) {u = -norm;} else {u = -norm2;}
				// ks=relative index, this trace. k_nbr = indx, nbr trace				
				for(kk=k+KDIP,ks=KDIP,k_nbr=0; ks<=ke; kk++,ks++,k_nbr++)
				  {
					s = v[kstart+ kk];
					// loop over E then N nbrs					
					for(i_EN=0; i_EN<2; i_EN++)
					  {
						j = js_EN[i_EN];  // trace indx of E or N nbr
						if (j<0) continue; // no nbr 
						js = j+k_nbr;
						je = js+2*KDIP;
						for(j=js; j<=je; j++) 
						  {
							// and here: wt of nbrs
							// For matrix in standard view: A[kk+kstart][j] = u;
							// a[j] += u * v[kstart+ kk];
							a[j] += u * s;
						  }
					  }
					// here
					// A[kk+kstart][kk] = sqrt_mu_smooth; 
					// a[kk] += sqrt_mu_smooth * v[kk+kstart];
					a[kk] += sqrt_mu_smooth * s;
				  }
			  }
		  }
	  }
  }


  final double getSmoothingEnergy(double[] v)
  {
	double chisq=0;
	if (addSmoothingTerm && (sqrt_mu_smooth>0))
	  {
		int j,je,js,k,kk,ks,ke=n-1-KDIP,k_nbr,nnbr;
		int i_tr,i_EN,Ntr=XY.length;
		int[] js_EN={-1,-1};
		double u,s;
		double norm = 1.0/(2*KDIP+1);
		double norm2 = norm/2.0;

		chisq=0.0;
		for(i_tr=0,k=0; i_tr<Ntr; i_tr++,k+=n)
		  {
			int[] _XY = XY[i_tr];
			// neighbor trace lookups
			js_EN[0] = (_XY[2]>0) ? (n*_XY[2]) : -1;
			js_EN[1] = (_XY[3]>0) ? (n*_XY[3]) : -1;
			for(i_EN=0,nnbr=0; i_EN<2; i_EN++) if (js_EN[i_EN]>=0) nnbr++;
			if (nnbr==0) continue;
			if (nnbr==1) {s=norm;} else {s=norm2;}
			// ks=relative index, this trace. k_nbr = indx, nbr trace				
			for(kk=k+KDIP,ks=KDIP,k_nbr=0; ks<=ke; kk++,ks++,k_nbr++)
			  {
				// loop over E then N nbrs					
				for(i_EN=0,u=0.0; i_EN<2; i_EN++)
				  {
					j = js_EN[i_EN];  // trace indx of E or N nbr
					if (j<0) continue; // no nbr 
					js = j+k_nbr;
					je = js+2*KDIP;
					for(j=js; j<=je; j++) u += v[j];					
				  }				
				chisq += sqr(v[kk] - u*s);
			  }
		  }
		chisq *= sqr(sqrt_mu_smooth);
	  }
	return chisq;
  }

  void runLARS()
  {
	int i;
	long t0,tt;

  
	EH.warn("Main thread..");


	// read wavelets, open seismic
	for(i=0;i<Nstacks;i++)
	  {
		sI[i].readWaveletSetUpX();
		sI[i].openSeismic();
	  }
	// 
	try 
	  {
		if (oFile==null)
		  {
			out_strm = (FileOutputStream) new FileOutputStream(FileDescriptor.out);
		  }
		else
		  {
			out_strm = (FileOutputStream) new FileOutputStream(oFile);
		  }
	  }
	catch (Exception e) 
	  {
		throw new RuntimeException("couldn't open output stream");
	  }

	if (useMRF) initialiseMRF();


	// A list of stuff to do. termination criteria later on examine
	// SSI.m - this could be altered suitably with work streaming in....
	m = 100000;

	PQ = new PriorityBlockingQueue();

	// Number of threads. Set default to N processors
	if (Nthreads < 1) Nthreads = Runtime.getRuntime().availableProcessors();
	EH.warn("Using "+Nthreads+" thread(s)");

	Task t[] = new Task[Nthreads];
	
	i_consumed=0;
	i_delivered=0;
	t0 = System.nanoTime();
	for (i=0; i<Nthreads; i++) 
	  {			
		t[i] = new Task(i,this);
		t[i].initialise();
	  }

	for (i=0; i<Nthreads; i++) 
	  {	
		t[i].start();
	  }
	for (i=0; i<Nthreads; i++) 
	  {
		try 
		  {
			t[i].join();
		  } 
		catch (java.lang.InterruptedException e) {};
	  }
	// for (i=0; i<Nthreads; i++)  EH.println ("Thread "+i+" ans = "+delay[i]);
	EH.warn("Gather in main thread");
	
	for(i=0;i<Nstacks;i++)
	  {
		sI[i].closeSeismic();
	  }
	suwriter.flush(out_strm);
	try
	  {
		out_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close SU out stream");
	  }

	tt = (System.nanoTime()-t0); 
	EH.warn("Processed "+i_consumed+" trace-groups");
	EH.warn("Total-time = "+(1e-9*tt)+", Time per op = "+1e-6*(tt/(double)i_consumed)+" ms");
	
	if (useMRF)
	  {
		performMRFcleanup();
	  }
 }


  // matrices etc shared across threads
  void initGlobals(int N)
  {
	this.N = N;
	Xt = new double[N][N];
  }  

  // place holder to return true if the
  // supplied index i indicates all the work will have been done
  boolean haveFinished(int i)
  {
	int j;

	if (count>0 && (i_consumed >= count))
	  {
		return true;
	  }
	for(j=0;j<Nstacks;j++)
	  {
		if (sI[j].stackExhausted) return true;
	  }
	return false;
	// return (i>=m);
  }

  synchronized SparseSpikeInverter.Job getTask(int i_th, Task T)
  {
	int i,k;
	SparseSpikeInverter.Job J;

	if (haveFinished(i_consumed)) // some suitable "run-out-of-work" criterion
	  {
		return null;
	  }

	if (freeJobList.isEmpty())
	  {
		J = new SparseSpikeInverter.Job(i_consumed,Nstacks,N,append);
	  }
	else 
	  {
		// reuse an old object, just give it a fresh number
		J = freeJobList.remove(freeJobList.size()-1);
		J.i = i_consumed;
	  }
	float[] Jyi;
	SuSegyFormat _tr;
	int i_s;

	for(i=0;i<Nstacks;i++)
	  {
		// seismic already in place from initialisation
		if (i_consumed > 0) sI[i].getSeismic();
		if (sI[i].stackExhausted) return null; // suitable "run-out-of-work" criterion
		_tr = sI[i].tr;
		Jyi = J.y[i];
		i_s = sI[i].i_s;
		// stuff data into Job
		for(k=0;k<N;k++) Jyi[k] = _tr.get_data_value(i_s+k);
		// keep a copy if we're appending
		if (append)
		  {
			System.arraycopy(Jyi,0,J.y_orig[i],0,N);
		  }

		// and copy of header, first stack
		if (i==0) _tr.copyHeader(J.hdr_bytes);
		
	  }

	i_consumed++;
	return J;
  }
 

  synchronized void doIO(boolean all,int i_th)
  {
	for(;;)
	  {
		tryNextIO(i_th);
		if (!all) break;
		else
		  {
			// Again
			if (haveFinished(i_delivered)) break;
		  }
	  }
  }


  synchronized void tryNextIO(int i_th)
  {
	SparseSpikeInverter.Job J;

	// dump out all the next chunk that's in order
	for(;;)
	  {
		if ((J = PQ.peek()) !=null)
		  {
			if (J.i == i_delivered)
			  {
				int i,k,i_s;
				float[] y;


				// pull job off queue
				PQ.poll();

				// output
				// copy header across
				tr_out.setHeader(J.hdr_bytes);
				if (append)
				  {
					// prepend true seismic
					i_s = sI[0].i_s;
					y = J.y_orig[0];
					// offset it back:
					for(k=0;k<N;k++) tr_out.set_data_value(y[k], i_s + k);
					// write out
					suwriter.fputtr(out_strm,tr_out);
				  }
				// assume desired seismic model now stored in y[0]
				y = J.y[0];
				// offset it back:
				i_s = sI[0].i_s;
				for(k=0;k<N;k++) tr_out.set_data_value(y[k], i_s + k);
				// write out
				suwriter.fputtr(out_strm,tr_out);
				if (useMRF)
				  {
					// keep copy of spikes for later MRF operations
					Slist.add((float[])y.clone());
				  }
				i_delivered++;
				// free up space of old job for a new one
				freeJobList.add(J);
			  }
			else break;
		  }
		else break;
	  }
  }

  /* 
	 Stuff to try MRF cleanup on slice of picks
  */
  void initialiseMRF()
  {
	Slist = new ArrayList<float[]>();
  }

  void performMRFcleanup()
  {
	int i,j,k,n,iq,q,it,Nt,Nsw,NMAX;
	PrintWriter out = null;
	float alpha,beta,scale,gamma;
	float[][] wt;
	float[][] S = 	Slist.toArray(new float[0][N]);
	float[] U = new float[2];
	int[] p = new int[2];
	int[][] P;
	float temp0,tempF,temp=0,tempR=0;
	boolean asciiIO=false;

	// Graphics
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDShapePlot plot = null;
	CanvasPlotter CP=null;
	ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups=null;
	double[] x,y,x1,y1,col=null;
	int[] gcol;

	Nt = S.length;
	Nsw = Nt*N;
	doGraphics=true;
	if (doGraphics)
	  {
		CP = new CanvasPlotter(true);
		Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		plot = new CanvasPlotter.TwoDShapePlot();
		x = new double[Nsw];
		y = new double[Nsw];
		x1 = new double[Nsw];
		y1 = new double[Nsw];
		col = new double[Nsw];
		gcol = new int[Nsw];
		for(i=0,j=0; i<Nt; i++)
		  {
			for(k=0;k<N;k++,j++) 
			  {
				x[j] = i+1;
				x1[j] = i+2;
				y[j] = k;
				y1[j] = k+1;
				gcol[j]=-1;
			  }
		  }
		plot.x = x;
		plot.x1 = x1;
		plot.y = y;
		plot.y1 = y1;
		plot.col = col;
		plot.gcol=gcol;
	  }
	if (asciiIO)
	  {
		try
		  {
			out = new PrintWriter(new FileWriter("Picks.txt"));
		  }
		catch (IOException e) {};
	  }

	// (1) set weights
	for(i=0,scale=0.0f; i<Nt; i++)
	  {
		for(j=0;j<N;j++) scale = max(scale,abs(S[i][j]));
	  }
	wt = new float[Nt][N];


	// Coupling params

	alpha=.3f;   // prize for vertical continuity across missing spike
	gamma=1.4f;   // prize for state flip vertically across spike
	beta=0.5f;    // prize for transverse continuity

	for (i=0;i<Nt;i++)
	  {
		for (j=0;j<N;j++)
		  {
			// vertical couplings
			if (S[i][j]==0.0) 
			  {
				wt[i][j] = -alpha;
			  }
			else 
			  {
				wt[i][j] = gamma * abs(S[i][j])/scale;
			  }
		  }
	  }

	// initialise Ising/Potts model
	P = new int[Nt][N];

	boolean doICM = false;
	boolean doFFMAP = true;

	if (doFFMAP)
	  {
		getMAPstateFF(Nsw,N,Nt,gamma, wt,beta,P);
	  }
	else
	  {
		for(i=0; i<Nt; i++)
		  {
			for(k=0,j=-1;k<N;k++) 
			  {
				// random start initialise as white noise
				// P[i][k] = (random()<0.5) ? -1 : 1;
			
				// set by counting spikes: start at -1
				if (S[i][k]!=0.0) j=-j;  // flip from -1 to 1 at reflection
				P[i][k] = j;
			  }
		  }

		if (doICM)
		  {
			NMAX=20;
		  }
		else
		  {
			// Annealing
			NMAX=10000;
			temp0=10f;
			tempF=0.1f;
			temp=temp0;
			tempR = (float)pow(tempF/temp0,1.0/NMAX);
		  }
	  

		for(it=0; it<NMAX;it++)
		  {
			if (!doICM && ((it%(NMAX/10))==0)) EH.warn("it="+it);
			for(n=0;n<Nsw;n++)
			  {
				// random pixel
				i = min((int)(random()*Nt),Nt-1);
				k = min((int)(random()*N),N-1);
				for(iq=0;iq<2;iq++)                    // states +/-1 =of ij pixel
				  {
					q=2*iq-1;  // +/- 1
					if (k>0)    U[iq] =  P[i][k-1] * q * wt[i][k];    // neighbour above
					if (k<N-1)  U[iq] += P[i][k+1] * q * wt[i][k+1];  // below
					if (i<Nt-1) U[iq] -= P[i+1][k] * q * beta;  // right
					if (i>0)    U[iq] -= P[i-1][k] * q * beta;  // left
				  }
				if (doICM)
				  {
					// greedy flip to lowest energy state
					P[i][k] = (U[0] < U[1]) ? -1 : 1;
				  }
				else
				  {
					// Annealing
					P[i][k] = (random() < 1.0/(1.0 + exp(-(U[1] - U[0])/temp))) ? -1 : 1;
				  }
			  }
			if (!doICM) 
			  {
				// annealing cooling
				temp *= tempR;
			  }
		  }
	  }

	if (asciiIO)
	  {
		for(i=0;i<Nt;i++)
		  {
			//		for(j=0;j<N;j++)
			for(j=N-1;j>=0;j--)
			  {
				out.printf("%d %d %d%n",i,j,P[i][j]);
			  }
		  }
		out.close();
	  }
	if (doGraphics)
	  {
		for(i=0,j=0;i<Nt;i++)
		  {
			for(k=0;k<N;k++,j++) col[j] = P[i][k];
		  }
		pg.GG.reverseY = true;
		pg.GG.computeRanges(plot.x,plot.y);
		pg.GG.computeRanges(plot.x1,plot.y1);
		pg.GG.setTickNumbersEngStyle(!true,!true);
		plot.normalize(); // rescales colors to [0-1]
		plot.setColorScheme(Color.white , Color.black);

		pg.plots.add(plot);
		pg.GG.computeRanges(plot.x,plot.y);
		pg.GG.Xword = "traces";
		pg.GG.Yword = "sample";

		CP.windowTitle = "MRF pick classification";
		CP.run(Pgroups); 
	  }
  }


  // Ford-Fulkerson MAP detection
  // Basic Definitions
  public static final int WHITE=0;
  public static final int GRAY=1;
  public static final int BLACK=2;

  int n;  // number of nodes
  double[][] capacity;  // capacity matrix
  double[][] flow;      // flow matrix
  int[] color;          // needed for breadth-first search               
  int[] pred;           // array to store augmenting path


  // A Queue for Breadth-First Search

  int head,tail;
  int[] q; // array size[n+2];

  void enqueue (int x) 
  {
	q[tail] = x;
	tail++;
	color[x] = GRAY;
  }

  int dequeue () 
  {
	int x = q[head];
	head++;
	color[x] = BLACK;
	return x;
  }

  
  /*
	N=total pixels
	Nz= # vertical pixels
	Nr= # horix pixels
	wt[i][k] = edge coupling, pixel ik
	beta = edge coupling horiontally
  */


  void getMAPstateFF(int N,int Nz,int Nt,float beta_max, float[][] wt, float beta,int[][] P)
  {
	int a,b,c,i,j,k,ii,jj,e;
	double s;

	n = N+2;              // # vertices in network [0...N-1]=pixels. [N]=source [N+1]=sink
	e = N + 2*(N-1);      // # edges (not needed)

	// workspace

	capacity = new double[n][n]; // initializes capacity matrix to 0
	flow = new double[n][n];
	color = new int[n];
	pred = new int[n];
	q = new int[n+2];

	// set edge capacities. Indices are: [0...N-1]=pixels. [N]=source [N+1]=sink
	// source, sink to nodes
	// in this problem, external field=0: lambda[i]=0

	if (!false)
	  {
		for (i=0;i<N;i++)
		  {
			s = 0.00001*(i%2-0.5);   // lambda[i]
			if (s>0) capacity[N][i] = s;
			else capacity[i][N+1] = -s;
		  }
	  }
	// internal edges (MRF couplings). beta_max added to make all capacities >=0
	for(i=0,ii=0;i<Nt;i++)
	  {
		for(k=0;k<Nz;k++,ii++)
		  {
			if (k>0) 
			  {
				capacity[ii][ii-1] = capacity[ii-1][ii] = (beta_max - wt[i][k]);    // neighbour above
			  }
			if (i>0)   
			  {
				jj = ii-Nz;
				capacity[ii][jj] = capacity[jj][ii] = beta_max + beta;                   // neighbour left
			  }
		  }
	  }


	// get max flow from source to sink
	max_flow(N,N+1);

	// set pixels
	for(i=0,ii=0;i<Nt;i++)
	  {
		for(k=0;k<Nz;k++,ii++)
		  {
			P[i][k] = (color[ii]==BLACK) ? -1 : 1;
		  }
	  }
	// clean up
	capacity = null;
	flow=null;
	color=null;
	pred=null;
	q=null;
  }
   
  // Breadth-First Search for an augmenting path
  boolean bfs (int start, int target) 
  {
	int u,v;
	for (u=0; u<n; u++) color[u] = WHITE;
	head = tail = 0;
	enqueue(start);
	pred[start] = -1;
	while (head!=tail) 
	  {
		u = dequeue();
		// Search all adjacent white nodes v. If the capacity
		// from u to v in the residual network is positive,
		// enqueue v.
		for (v=0; v<n; v++)
		  {
			if ((color[v]==WHITE) && ((capacity[u][v]-flow[u][v])>0))
			  {
				enqueue(v);
				pred[v] = u;
			  }
		  }
	  }
    // If the color of the target node is black now,
    // it means that we reached it.
	return (color[target]==BLACK) ? true : false;
  }

  // Ford-Fulkerson Algorithm
  double max_flow(int source, int sink) 
  {
    int i,j,u;
	double increment;

    // Initialize empty flow.
    double max_flow = 0.0;
    for (i=0; i<n; i++) 
	  {
		for (j=0; j<n; j++) flow[i][j] = 0;
	  }
    
    // While there exists an augmenting path,
    // increment the flow along this path.
    while (bfs(source,sink)) 
	  {
        // Determine the amount by which we can increment the flow.
		increment = Double.MAX_VALUE;
		for (u=n-1; pred[u]>=0; u=pred[u]) 
		  {
			increment = min(increment,capacity[pred[u]][u]-flow[pred[u]][u]);
		  }
        // Now increment the flow.
		for (u=n-1; pred[u]>=0; u=pred[u]) 
		  {
			flow[pred[u]][u] += increment;
			flow[u][pred[u]] -= increment;
		  }
		max_flow += increment;
	  }
    // No augmenting path anymore. We are done.
    return max_flow;
  }

  // variables for CG method
  boolean doTruncatedNewton = false;
  private double ATr_norm_lastuse = Double.MAX_VALUE;
  private double eta_TN = 0.5;
 /**
	 Set the method to use truncated Newton. Will terminate at criterion
	 determined by previous successful exit (see eg Dembo-Steihaug (1983)
	 Essentially imposes geometric decrease of termination residual
	 r = ||A^T.(b-Ax)||
	 on successive calls. Useful when LSQR is the inner problem in an 
	 outer optimisation loop.
	 Invoke this after constructor if desired.

	 @param use sets Truncated Newton on or off
	 @param eta sets gain ratio 0<eta<1 required for TN method to get linear convergence of rate eta
  */
  void useTruncatedNewtonMethod(boolean use,double eta)
  {
	this.doTruncatedNewton = use;
	this.eta_TN = eta;
  }
  /**
	 Solve the ordinary least squares problem<p>
	 x = argmin ||Ax - b||^2<p>

	 The matrix A is implicitly represented by Av_product.
	 @param x unknown vector of least squares system
	 @param b data or RHS vector of least squares system
	 @return solution is placed into x[]
 */
  private void CG_solve(boolean doWarmStart,double[] b,double[] x,int Nmax)
  {
	// CG method of Paige & Saunders (LSQR)
	// relies on calling class's Av_product() 
	// Unknown is of size N, rhs of size Nd
	int i,j,k,N,Nd;
	double[] _u,_v,_u_next,_v_next;
	double[] u,u_next,v,v_next,w,tmp;
	double[] beta={0},alpha={1};
	double[] beta_next={0},alpha_next={0};
	double s,phi,phibar,rho,rhobar,c,theta_next,rhobar_next,phibar_next,bbnorm;
	double ATy_norm=0.0;

	double eps_CG = 1e-6;
	double resid_tol_match=1e-14; 

	boolean verbose = verbosity>2;

	// parallel stuff
	ExecutorService es=null;
	CompletionService<Void> cs=null;
	if (Nthreads>1)
	  {
		es = Executors.newFixedThreadPool(Nthreads);
		cs = new ExecutorCompletionService<Void>(es);
	  }

	N = x.length;
	Nd = b.length;

	// actual memory
	_u = new double[Nd];
	_u_next = new double[Nd];
	_v = new double[N];
	_v_next = new double[N];
	w = new double[N];
	//	x = double[N], passed in
	// must initialise x
	// swapping pointers
	u=_u; 
	v=_v;
	u_next = _u_next;
	v_next = _v_next;

	// (1) Initialize
	if (doWarmStart)
	  {
		Av_product(false,x,u);   
		for(i=0; i<Nd; i++) u[i] = b[i] - u[i]; // u = b-A.x
	  }
	else
	  {
		System.arraycopy(b,0,u,0,Nd);		
		Arrays.fill(x,0.0);
	  }


	normalize(cs,u,beta);

	Av_product(true,u,v);
	normalize(cs,v,alpha);
	if (doTruncatedNewton) ATy_norm = alpha[0]*beta[0];
	System.arraycopy(v,0,w,0,N);

	phibar=beta[0];
	rhobar=alpha[0];
	// matrix norm
	bbnorm = 0.0;

	// iterates	
	for(j=0;j<min(100*N,Nmax);j++)
	  {
		// P & S step 3: bi-diagonalization
		Av_product(false,v,u_next);
		s=alpha[0];
		if (Nthreads==1) for(i=Nd-1; i>=0; i--) u_next[i] -= s * u[i];
		else psum(cs,true,Nd,u_next,null,u,-s);

		normalize(cs,u_next,beta_next);

		Av_product(true,u_next,v_next);
		s = beta_next[0];
		if (Nthreads==1) for(i=N-1; i>=0; i--) v_next[i] -= s * v[i];
		else psum(cs,true,N,v_next,null,v,-s);

		normalize(cs,v_next,alpha_next);

		// accumulate this quantity to estimate Frobenius norm of matrix A
		bbnorm += sqr(alpha_next[0])+ sqr(beta_next[0]);
		// next orth. transform.
		rho = hypot(rhobar,beta_next[0]);
		//rho = sqrt(sqr(rhobar)+sqr(beta_next[0]));
		c =  rhobar/rho;
		s = beta_next[0]/rho;
		theta_next = s*alpha_next[0];
		rhobar_next = -c*alpha_next[0];
		phi = c*phibar;
		phibar_next  = s* phibar;
	  
		// update x,w
		s = phi/rho;
		if (Nthreads==1) for(i=N-1; i>=0; i--) x[i] += s*w[i];
		else psum(cs,true,N,x,null,w,s);
		s = theta_next/rho;
		if (Nthreads==1) for(i=N-1; i>=0; i--) w[i] = v_next[i] - s*w[i];
		else psum(cs,false,N,w,v_next,w,-s);

		// test for convergence...
		// ATr_norm = ||A^T r_i|| = eqn 5.4 P&S
		// here ||r|| = phibar_next (5.2), P & S
		double ATr_norm = phibar_next*alpha_next[0]*abs(c);
		// Nash
		if (doTruncatedNewton)
		  {
			double QR;
			//QR = (j+1)*ATr_norm/ATy_norm;   // Based on Dembo-Steihaug 1983 ??
			if (ATr_norm_lastuse == Double.MAX_VALUE)
			  {
				// init cycle
				QR = ATr_norm/ATy_norm;
			  }
			else
			  {
				QR = ATr_norm/ATr_norm_lastuse;
			  }
			if (QR < eta_TN) 
			  {
			    if (verbose) System.out.printf("Trunc. Newt. break ||r||=%8g, ATr_norm=%8g%n",phibar_next,ATr_norm);
				ATr_norm_lastuse= ATr_norm;
				break;
			  }
		  }

		if (verbose && (j%100==0)) System.err.printf("j=%d, ||r||=%8g%n",j,phibar_next);

		// Stopping criteria S1
		// S1 = phibar_next/bbnorm;
		s = (bbnorm!=0) ? (phibar_next/bbnorm) : 0.0;
		if (s< resid_tol_match)
		  {
			if (verbose) System.out.printf("Broke at j=%d, S1=%8g%n",j,s);
			break;
		  }



		// stopping criterion S2 ||A^T.r_k||/||A||||r_k|| < eps
		double S2 = ATr_norm/(sqrt(bbnorm)*phibar_next);

		if (S2 < eps_CG) 
		  {
			if (verbose) EH.warn("Broke at j=%d, ||S2||=%8g%n",j,S2);
			break;
		  }
		// reset variables to point at correct place 
		alpha[0] = alpha_next[0];
		rhobar = rhobar_next;
		phibar =  phibar_next;
		tmp = u_next; u_next = u; u = tmp;
		tmp = v_next; v_next = v; v = tmp;
	  }
	if (Nthreads>1) es.shutdown();
  }
  /*
	Parallel running sum update, over vector a
	If addStimesConly==true
	  a += s*c
	else
	  a = b+s*c
  */
  final void psum(CompletionService<Void> cs,boolean addStimesConly,int Nsum,double [] a,double[] b,double[] c,double s)
  {
	int k,Nch;
	Nch = Nsum/Nthreads+1;
	for (k=0; k<Nthreads; k++) 
	  {
		final int k1=k*Nch;
		final int k2=min(Nsum,k1+Nch);
		final double[] _a=a;
		final double[] _c=c;
		final double _s=s;
		if (addStimesConly) 
		  {
			cs.submit(new Runnable() {public void run() {sum1(k1,k2,_a,_c,_s);}},null);
		  }
		else
		  { 
			final double[] _b=b;
			cs.submit(new Runnable() {public void run() {sum2(k1,k2,_a,_b,_c,_s);}},null);
		  }
	  }
	try 
	  {
		for (k=0; k<Nthreads; k++) cs.take();
	  } 
	catch (InterruptedException ie) 
	  {
		throw new RuntimeException(ie);
	  }
  }

  private final void sum1(int k1, int k2,double[] x,double[] y,double s)
  {
	for(int k=k1;k<k2;k++) x[k] += s * y[k];
  }

  private final void sum2(int k1, int k2,double[] x,double[] y,double[] z,double s)
  {
	for(int k=k1;k<k2;k++) x[k] = y[k] + s * z[k];
  }

  private final void normalize(CompletionService<Void> cs,double[] u,double[] norm)
  {
	if (Nthreads==1)
	  {
		// same-thread version
		normalizeST(u,norm);
	  }
	else
	  {
		// parallelized
		pnormalize(cs,u,norm);
	  }
  }

  private final void normalizeST(double[] u,double[] norm)
  {
	int i,n=u.length;
	double s=0.0;
	for (i=n-1; i>=0; i--) s += u[i]*u[i];
	norm[0] = s = sqrt(s);
	if (s>0) for(i=n-1; i>=0; i--) u[i]/=s;
  }

  // variables for summations
  private double __SS;
  private synchronized void accrueS(double s) {__SS+=s;}

  final void pnormalize(CompletionService<Void> cs,double[] a,double[] norm)
  {
	int k,Nch;
	int Nsum = a.length;
	Nch = Nsum/Nthreads+1;

	__SS=0;
	final double[] _a=a;
	for (k=0; k<Nthreads; k++) 
	  {
		final int k1=k*Nch;
		final int k2=min(Nsum,k1+Nch);
		cs.submit(new Runnable() {public void run() 
			{
			  double s=0.0;
			  for (int i=k1;i<k2;i++) s += _a[i]*_a[i];
			  accrueS(s);
			}
		  },null);
	  }
	
	try {for (k=0; k<Nthreads; k++) cs.take();}
	catch (InterruptedException ie) {throw new RuntimeException(ie);}

	norm[0] = sqrt(__SS);
	
	final double s = norm[0];
	if (s>0) 
	  {
		for (k=0; k<Nthreads; k++) 
		  {
			final int k1=k*Nch;
			final int k2=min(Nsum,k1+Nch);
			cs.submit(new Runnable() {public void run() {for (int i=k1;i<k2;i++) _a[i]/=s;}},null); 
		  }
		try {for (k=0; k<Nthreads; k++) cs.take();}
		catch (InterruptedException ie) {throw new RuntimeException(ie);}
	  }
  }

  
  public static class Job implements Comparable
  {
	int i=0;
	// header bytes for use
	byte[] hdr_bytes = new byte[SuSeismicConstants.HDRBYTES];
	// data copy over stacks, also working variables
	float[][] y;
	// data copy over stacks
	float[][] y_orig;

	public Job(int j,int Nstacks,int N,boolean takeCopy)
	{
	  i = j;
	  y = new float[Nstacks][N];
	  if (takeCopy)
		{
		  y_orig = new float[Nstacks][N];
		}
	}

	public int compareTo(Object o) 
	{ 
	  int i,l;
	
	  Job m1 = this;
	  Job m2 = (Job)o;

	  return (m1.i < m2.i) ? -1 : ((m1.i > m2.i) ? 1 : 0);
	}
  }

  public static class StackAndWavInfo
  {
	SparseSpikeInverter _SSI;

	// wavelet stuff
	String wFile;
	int ns,idel1;
	float f1;
	int n_precursor,n_coda;
	double[] Xw;

	// optional entry if wavelet adjoint used (IRLS)
	double[] Xwt;

	float noise;
	SuReader sureader = new SuReader();
	SuSegyFormat tr = new SuSegyFormat();
	FileInputStream su_in_strm=null;

	// stack info
	String sFile;
	int i_s;
	boolean stackExhausted=false;

	public StackAndWavInfo(SparseSpikeInverter SSI)
	{
	  this._SSI = SSI;
	}

	void readWaveletSetUpX()
	{
	  int i,off;
	  double tmp;

	  // open wavelets
	  try 
		{
		  su_in_strm = (FileInputStream) new FileInputStream(wFile);			
		} 
	  catch (IOException e) 
		{
		  _SSI.EH.die_err("Couldn't open " +wFile);
		}
	  
	  i = sureader.fgettr(su_in_strm,tr,_SSI.fromLE);

	  if (i == 0)
		{
		  _SSI.EH.die_err("Couldn't read "+wFile);
		}
	  try 
		{
		  su_in_strm.close();
		} 
	  catch (IOException e) 
		{
		  throw new RuntimeException("couldn't close " + wFile);
		}

	  ns = tr.get_short("ns");  
	  if (_SSI.verbosity > 0) _SSI.EH.warn("wavelet: ns="+ns);
	  f1 = tr.get_float("f1");
	  if (f1==0) f1 = tr.get_short("delrt") * 1e-3f;
	  // idel1 = -ve index denoting first precursor sample
	  idel1 = (int)Math.floor(f1/(tr.get_short("dt")*1e-6)+0.5);
	  Xw = new double[ns];
	  // copy
	  for (i=0;i<ns;i++) Xw[i] = (double)(tr.get_data_value(i));
	  // wavelet zero-time sample is at Xw[-idel1]

	  // make as compact as possible; threshold small values relative to peak
	  double WAV_EPS=0.001f;
	  for (i=0,tmp=0;i<ns;i++) tmp = max(tmp,abs(Xw[i]));
	  tmp *= WAV_EPS;
	  // precursor
	  n_precursor = idel1;
	  for (i=0;i<min(-idel1,ns);i++)
		{
		  if (abs(Xw[i]) > tmp) break;
		  n_precursor++;
		}
	  // precursor useful samples are now i=[wM.n_precursor ...-1]

	  // coda
	  n_coda = ns+idel1-1;
	  for (i=ns-1; i>=0; i--)
		{
		  if (abs(Xw[i]) > tmp) break;
		  n_coda--;
		}
	  // coda useful samples are now i=[0...wM.n_coda]

	  // make clean copy corresponding to this coda etc
	  // reset ns, idel1
	  ns = n_coda -n_precursor+1;
	  double[] Xw_new = new double[ns];
	  off = n_precursor-idel1;
	  for(i=0;i<ns;i++) Xw_new[i] = Xw[i+off];
	  idel1 = n_precursor;
	  Xw = Xw_new;
	}

	int initWindowA()
	{
	  double f1s,dt,eps=1e-6;  

	  if ((f1s = tr.get_float("f1"))==0)	f1s = (tr.get_short("delrt") * 1e-3);
	  if (_SSI.tmin<(f1s-eps))
		{
		  _SSI.EH.die_err("first trace has no samples as early as tmin="+_SSI.tmin);
		}
	  dt = tr.get_short("dt")*1e-6;
	  int ns = tr.get_short("ns");
	  i_s = (int)rint((_SSI.tmin - f1s)/dt);
	  if (_SSI.tmax > (f1s+(ns-1)*dt + eps))
		{
		  _SSI.EH.die_err("first trace has no samples as late as tmax="+_SSI.tmax);
		}
	  return i_s;
	}

	void initWindow()
	{
	  int i,j,k,_N=0;
	  double tmp;

	  initWindowA();
		
	  // initialise output
				
	  _SSI.tr_out.set("ns",(short)ns);
	  _SSI.tr_out.set("dt",tr.get_short("dt"));


	  // fill Xt
	  // FIXME, for multi stacks this will of course be a block matrix...
	  // DenseDoubleMatrix2D Xt = _SSI.Xt;
	  double[][] Xt = _SSI.Xt;
	  _N = _SSI.N;
	  for (j=0;j<_N;j++)
		{
		  for (i=0;i<_N;i++)
			{
			  k = i-j;
			  if (k>=n_precursor && k<=n_coda) 
				{
				  // Xt.setQuick(j,i,Xw[k-idel1]);
				  Xt[j][i] = Xw[k-idel1];
				}
			}
		}

	  if (_SSI.USE_LOOKUP)
		{
		  // complete sG_A_lookup = X^T.X multiply, so we can lookup active components of these as we please
		  _SSI.sG_A_lookup = new DenseDoubleMatrix2D(_N,_N);
		  // SeqBlas.seqBlas.dgemm(false,true,1,Xt,Xt,0,_SSI.sG_A_lookup); 
		  for (i=0;i<_N;i++)
			{
			  for (j=i;j<_N;j++)
				{
				  for(k=0,tmp=0;k<_N;k++) tmp += Xt[i][k]*Xt[j][k];
				  _SSI.sG_A_lookup.setQuick(i,j,tmp);
				  if(j!=i) _SSI.sG_A_lookup.setQuick(j,i,tmp);
				}
			}
		}
	}

	void setUpAdjoint()
	{
	  int n,i,j;
	  Xwt = Xw.clone();
	  n = Xw.length;

	  for(i=0,j=n-1;i<n;i++,j--) Xwt[j] = Xw[i]; 
	}


	void getSeismic()
	{
	  int n;
	  n = sureader.fgettr(su_in_strm,tr,_SSI.fromLE);
	  if (n == 0)
		{
		  _SSI.EH.warn("Exhausted seismic "+sFile);
		  stackExhausted = true;
		}
	}

	void openSeismic()
	{
	  // open seismic
	  try 
		{
		  if (sFile.equals("stdin"))
			{
			  su_in_strm  = (FileInputStream) new FileInputStream(FileDescriptor.in);
			}
		  else
			{
			  su_in_strm = (FileInputStream) new FileInputStream(sFile);			
			} 
		}
	  catch (IOException e) 
		{
		  _SSI.EH.die_err("Couldn't open seismic " +sFile);
		}
	}

	void closeSeismic()
	{
	  try 
		{
		  su_in_strm.close();
		} 
	  catch (IOException e) 
		{
		  throw new RuntimeException("couldn't close " + sFile);
		}
	}
  }
}


class Task extends Thread 
{
  int i_th;
  SparseSpikeInverter _SSI;

  // work variables for LARS
  double[] mu,beta,c_hat,a,u_A,wk,y_opt,beta_full,beta_opt;
  int[] s,act;
  boolean[] A,Alast;
  DenseDoubleMatrix2D sG_A_L; // _sG_A_lookup;

  // small stuff for profiling
  boolean profile=false; 
  long[] tT0,tT;


  public Task(int i,SparseSpikeInverter SSI) 
  {
	_SSI = SSI;
	i_th = i;
	if (profile)
	  {
		tT0 = new long[6];
		tT = new long[6];
	  }
  }
  
  // called once, for each thread, reading first seismic trace-group on thread 0,
  // but doing thread-local allocs for all threads
  void initialise()
  {
	int i,j,k,N=0;
	float f1;

	// do some initialisations
	// initialise on first trace
	SuSegyFormat _tr = _SSI.sI[0].tr;
	if (i_th==0)
	  {
		// read seismic only once, on "first" thread
		// set up some stuff from it
		for (i=0;i<_SSI.Nstacks;i++) 
		  {			
			_SSI.sI[i].getSeismic();			
			N = (int)( (_SSI.tmax-_SSI.tmin)/(_tr.get_short("dt")*1e-6));
			if (i>0 && N!=_SSI.N)
			  {
				_SSI.EH.die_err("dt appears to vary between stacks??? resample to same!");
			  }
			if (i==0) _SSI.initGlobals(N); // sets _SSI.N
			if (_SSI.verbosity>0 && (i==0)) _SSI.EH.warn("Vertical system size (N samples in window)="+N);
			_SSI.sI[i].initWindow();
		  }
	  }
	else
	  {
		N = _SSI.N;
	  }

	// thread-locals
	mu = new double[N];
	if (!_SSI.show_synth) 
	  {
		beta = new double[N];
		beta_full = new double[N];
	  }
	c_hat = new double[N];
	s = new int[N]; // zero s
	act = new int[N];
	wk = new double[N];
	a = new double[N];
	u_A = new double[N];
	A = new boolean[N];
	Alast = new boolean[N];
	// N_fbytes = N*FSIZE;

	if (_SSI.USE_CHOL_UPDATE) sG_A_L = new DenseDoubleMatrix2D(N,N);
	if (_SSI.useBIC)
	  {
		y_opt = new double[N];
		beta_opt = new double[N];
	  }		

	//	_sG_A_lookup = (DenseDoubleMatrix2D)_SSI.sG_A_lookup.like().assign(_SSI.sG_A_lookup);
  }



  public void run() 
  {
	double delay;
	long tt,t0;
	SparseSpikeInverter.Job J;
	boolean done=false;
	int i;

	_SSI.EH.warn("Running thread " + i_th);

	for(;;)
	  {
		// get a job
		if ((J = _SSI.getTask(i_th,this)) == null) 
		  {
			done=true;
			break;
		  }

		doWork(J);

		// dump to queue
		_SSI.PQ.add(J); 

		// _SSI.EH.println("QSIZE: "+_SSI.PQ.size());
		
		// handle some IO as well
		_SSI.doIO(false,i_th);	
	  }
	// at end of all the work, finish off IO
	if (done)
	  {
		_SSI.doIO(true,i_th);
	  }

	if (profile)
	  {
		for(i=0;i<tT.length;i++)
		  {
			_SSI.EH.warn("Profile: chunk %d -time(ms)=%4g%n",i,(tT[i]*1e-6)); // profiling
		  }
	  }
  }

  /**
	 Work done by thread on each "job". Here, "sparse-spike" inversion
	 via the LARS algorithm
  */
  void doWork(SparseSpikeInverter.Job J)
  {
	int i,j,k,it,A_info,info=0,N,Na,Na_opt;
	double t=0;
	double tmp,chiBIC_min,chiBIC_last,min,A_A,C_hat,epsC_hat,gamma_hat,chsq,tmp1,x;

	boolean do_break;
	double[] rhs;
	double[][] Xat=null,sG_A;
	double[] Xat_ur=null;
	float[] y;

	int ITMAX=500;
	double EPS=1e-5;
	double est_noise;
	boolean unroll_Xat = true;


	N = _SSI.N;


	if (_SSI.Nstacks!=1) _SSI.EH.die_err("Not coded for Nstacks>1 yet");

	// LARS algorithm
	// data window: do first stack for now
	y = J.y[0];

	if (_SSI.est_RMS_noise <= 0) 
	  {
		// use absolute level
		if (_SSI.est_noise <= 0)
		  {
			// use default noise level 0.15 max{|y|}
			est_noise=0;
			for (i=0;i<N;i++) est_noise = max(abs(y[i]),est_noise);
			est_noise *= 0.15;
		  }	
		else
		  {
			est_noise = _SSI.est_noise;
		  }
	  }
	else 
	  {
		// guess suitable trace-local noise level from supplied RMS fraction, for this trace
		for (i=0,tmp=0;i<N;i++) tmp += y[i]*y[i];
		est_noise = _SSI.est_RMS_noise * (float)sqrt(tmp/N);
	  }
		
	// LARS ALGORITHM		
	// zero mu
	Arrays.fill(mu,0);
	if (!_SSI.show_synth) 
	  {
		Arrays.fill(beta,0); // zero beta
		Arrays.fill(beta_full,0); // zero beta
	  }
	Arrays.fill(Alast,false);
	Na=0;
	A_info=0;
	chiBIC_min = Float.MAX_VALUE;
	Arrays.fill(A,false);
	for (it=0;it<ITMAX;it++)
	  {
		// Reasonably literal version from Efron et al 
		// get c_hat, C_hat, wk = y-mu
		C_hat = getC_hat(N,y,mu,c_hat,wk,s);
		// fprintf(stderr,"First C_hat=%g\n",C_hat); 
		// grow the active set
		epsC_hat = EPS*C_hat;
		if (true)
		  {
			// Forcibly search for nearest new member to C_hat
			for(i=0,k=0,min=Float.MAX_VALUE; i<N; i++) 
			  {
				if (!A[i])
				  {
					if ((tmp = abs(abs(c_hat[i]) - C_hat)) < min)
					  {
						k = i;
						min=tmp;
					  }
				  }
			  }
			tmp = min/epsC_hat;
			if (tmp>1)
			  {
				_SSI.EH.warn("Active set has numerical roundoff issues: added value which missed threshold by "+tmp);
				// I think, once roundoff error really sets in, it's best to give up....
				break;
			  }
			A[k]=true;
			j = Na+1;
		  }
		else
		  {
			// Closer to Efron's "c_hat_j=C_max" numerical criterion. B/c of roundoff, it is 
			// hard to get always increments of
			// exactly one new member in the active set
			Arrays.fill(A,false); // zero A
			// Naive simple search for members in window:   |c_hat[i]) - C_hat| < epsC_hat
			for(i=0,j=0;i<N;i++) 
			  {
				if ((abs(abs(c_hat[i]) - C_hat) < epsC_hat))
				  {
					A[i]=true; 
					j++;
				  } 
			  }
		  }

		if (j==(Na+1))
		  {
			// normal increment
			Na++;
		  }
		else if (j==Na)
		  {
			// One new guy didn't quite make it in. Search inactives for nearest miss j
			for(i=0,min=Float.MAX_VALUE; i<N; i++) 
			  {
				if (!A[i])
				  {
					if ((tmp = abs(abs(c_hat[i]) - C_hat)) < min)
					  {
						k = i;
						min=tmp;
					  }
				  }
			  }
			_SSI.EH.warn("Active set not incrementing by one. Added value which missed threshold by "+min/epsC_hat);

			A[k]=true;
		
			Na++;
		  }	
		else
		  {
			_SSI.EH.warn("Added "+(j-Na)+" actives?? BUG? Numerics?");
			Na = j;
			A_info=1; // flag disaster. Enshrined for rest of algorithm	
		  }

		if (!true)
		  {
			System.err.printf("Actives: ");
			for(j=0;j<N;j++) System.err.printf("%d ",A[j]?1:0);
			System.err.println("");
		  }
		if (A_info==0)
		  {
			// Normal behaviour; increment active components: tack onto end the new one
			// _SSI.EH.warn("Normal\n");
			for(i=0;i<N;i++)    
			  {
				if (A[i] && !(Alast[i]))
				  {
					act[Na-1] = i;
					Alast[i] = true;
					break;
				  } 
			  }
		  }
		else
		  {
			//_SSI.EH.warn("Reform\n");
			// re-form list afresh
			for(i=0,j=0;i<N;i++) 
			  {
				if (A[i])
				  {
					act[j]=i;
					j++;
				  }
			  }
			System.arraycopy(A,0,Alast,0,N);
		  }
		
		if ((Na>0) && (Na<=N))
		  {
			if (_SSI.verb>1) _SSI.EH.warn("Na="+Na);
		  }
		else 
		  {
			// terminate
			break;
		  }
		// LS on active subset
		if (!unroll_Xat)
		  {
			Xat = new double[Na][N]; // Xat = matrix(0,Na-1,0,N-1);
		  }
		else
		  {
			Xat_ur = new double[N*Na]; // Xat^T, unrolled
		  }
		rhs = new double[Na];    // alloc1float(Na);


		// fill Xa^T:
		if (profile) tT0[0] = System.nanoTime();
		double[][] _Xt = _SSI.Xt;
		for(j=0;j<Na;j++)
		  {
			i = act[j];
			// _SSI.EH.warn("j=%d i=%d\n",j,i);
			if (unroll_Xat) 
			  {
				double[] p1 = _Xt[i];
				System.arraycopy(_Xt[i],0,Xat_ur,j*N,N);
				// for(k=0;k<N;k++) Xa_ur[k*Na+j] = p1[k];   // Xa[k][j]
			  }
			else System.arraycopy(_Xt[i],0,Xat[j],0,N);
			// form S.1_A while we're here
			rhs[j] = s[i];
		  } 
		if (j!=Na) _SSI.EH.warn("BUG j!=Na="+Na);



		// form sG_A
		sG_A = new double[Na][Na]; // sG_A = Xa^T.Xa

		if (profile) tT[0] += (System.nanoTime()-tT0[0]);

		if (_SSI.USE_LOOKUP)
		  {
			// use lookups
			int ii,jj;
			DenseDoubleMatrix2D _sG_A_lookup = _SSI.sG_A_lookup;
			double[] p1;
			for (i=0;i<Na;i++)
			  { 
				ii = act[i]; 
				p1 = sG_A[i];
				for(j=0;j<=i;j++)
				  {
					jj = act[j];
					p1[j] = _sG_A_lookup.getQuick(ii,jj); // [ii][jj];							
					//sG_A[j][i] = p1[j]; // omit upper triangle for now
				  }							
			  }
		  }
		else
		  {
			_SSI.EH.die_err("FIX Implementation 1");
			// Not implemented yet
			// cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,Na,Na,N, 1.0, &(Xat[0][0]),N, &(Xat[0][0]),N, 0.0,&(sG_A[0][0]),Na);
		  }


		// solve sG_A^-1  . rhs
		/*	int clapack_sposv(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
			const int N, const int NRHS, float *A, const int lda,
			float *B, const int ldb)
		*/

		// info=0 here means conventional incremental update
		if (profile) tT0[1] = System.nanoTime();
		if (_SSI.USE_CHOL_UPDATE && (A_info==0)) 
		  {
			// assumes here that Na grows incrementally from 1...Na
			// Cholesky update			
			info = do_cholesky_update(Na, sG_A_L, sG_A, rhs, wk);
		  }
		else 
		  {
			_SSI.EH.die_err("Full Cholesky not coded yet...");
			// FINISH:
			// info = clapack_sposv(CblasRowMajor,CblasUpper,Na,1,sG_A[0],Na,rhs,Na);				
		  }  
		if (profile) tT[1] += (System.nanoTime()-tT0[1]);

		if (info!=0)
		  {
			_SSI.EH.die_err("Cholesky failed! Exiting");
		  }

		// get A_A
		for (j=0,A_A=0;j<Na;j++) A_A += s[act[j]] * rhs[j];
		A_A = 1.0/sqrt(A_A);
		// form u_A = A_A X_A.rhs
		
		if (profile) tT0[2] =System.nanoTime(); // profiling
		if (false)
		  {
			// FIXME finish
			// cblas_sgemv(CblasRowMajor, CblasTrans,Na,N,A_A, &(Xat[0][0]),N, &(rhs[0]),1, 0.0,&(u_A[0]),1);
		  }
		else
		  {
			// hand coded version
			// Slow here (25% time?). Parallelizes badly

			if (!unroll_Xat)
			  {
				for(i=0;i<N;i++)
				  {
					for(j=0,tmp=0;j<Na;j++) tmp += Xat[j][i]*rhs[j];
					u_A[i] = A_A * tmp;
				  }
			  }
			else 
			  {
				if (true)
				  {
					// possible blas use here Xat is m x n = Na x N, unrolled
					BlasLapackUtil.dgemv(Na,N,true,A_A,Xat_ur,rhs,0.0,u_A);
				  }
				else
				  {
					for(i=0;i<N;i++)
					  {
						for(j=0,tmp=0;j<Na;j++) tmp += Xat_ur[j*N+i]*rhs[j]; // tmp += Xat[j][i]*rhs[j];
						u_A[i] = A_A * tmp;
					  }
				  }
			  }
		
		  }
		if (profile)
		  {
			tT[2] += (System.nanoTime()-tT0[2]); // profiling		  
		  }

		// form a
		if (profile) tT0[3] =System.nanoTime(); // profiling
		get_Xt_vec_prod(N,a,u_A);
		if (profile) tT[3] += (System.nanoTime()-tT0[3]); // profiling


		// get gamma_hat
		if (profile) tT0[4] =System.nanoTime(); // profiling
		gamma_hat = Float.MAX_VALUE;
		if (Na>0)
		  {
			if (Na<N)
			  {
				int jA=-1; 
				for(j=0;j<N;j++)
				  {
					// search the inactive space
					if (!A[j])
					  {
						tmp = (C_hat-c_hat[j])/(A_A - a[j]);
						if ((tmp>0) && (tmp<gamma_hat)) {gamma_hat = tmp; jA=j;}
						tmp = (C_hat+c_hat[j])/(A_A + a[j]);
						if ((tmp>0) && (tmp<gamma_hat)) {gamma_hat = tmp; jA=j;}
					  }
				  }
			  }
			else 
			  {
				// FIXME need also case here where problem uses all components: gamma_hat = C_hat/A_A in Efron
				gamma_hat = C_hat/A_A;
			  }
			// _SSI.EH.warn("New active j=%d, new C_hat=%g\n",jA,C_hat - gamma_hat*A_A);
		  }
		else gamma_hat=1.0;

		// "finalised" model for use in BIC if required
		// update mu as if last active entry was the final in the basiss set
		tmp1 = C_hat/A_A; // " = gamma_hat_finish"
		for(j=0; j<N; j++) wk[j] = mu[j] + tmp1 * u_A[j];
		// wk now contains mu_full

		// update mu
		for(j=0; j<N; j++) mu[j] += gamma_hat * u_A[j];


		if (!_SSI.show_synth)
		  {
			// keep track of spike values
			tmp = A_A*gamma_hat;
			// full & partial updates
			for(j=0;j<Na;j++)
			  {
				i = act[j];
				x = s[i] * rhs[j];
				beta_full[i] = beta[i] + tmp1 * x;
				beta[i] += tmp * x;
			  }
		  }
		// check fit, exit at suitable noise level match

		do_break=false;
		chsq=0;
		if (true)
		  {
			// use "finished" model to make chisq comparisons
			for(i=0;i<N;i++) chsq += ((tmp=(y[i]-wk[i])) * tmp);
		  }
		else
		  {
			for(i=0;i<N;i++) chsq += ((tmp=(y[i]-mu[i])) * tmp);
		  }
		if (_SSI.verbosity > 1) _SSI.EH.warn("FIT: DOF=%d, (target-noise=%g, est-noise=%g), chisq=%g\n",N-Na,
											 est_noise,sqrt(chsq/(N-Na)),chsq/est_noise);
		chsq /= (est_noise*est_noise);
		if (!_SSI.useBIC)
		  {
			// Break criterion based on DOF-adjusted discrepancy principle
			if (chsq <= (N-Na)) 
			  {
				do_break=true;
				// "Finished solution"....
				if (_SSI.show_synth)
				  {
					System.arraycopy(wk,0,mu,0,N); // memcpy(mu,wk,N_fbytes);
				  }   
				else
				  {
					System.arraycopy(beta_full,0,mu,0,N);
				  }
			  }
		  }
		else
		  {
			// BIC method: add model penalty.
			// curve is a bit zigzaggy, so allow a little tolerance, i.e. overrun the minimum a little
			tmp = Na*log(N);
			if (_SSI.verbosity>1) _SSI.EH.warn("Na= %d Chisq= %g BIC= %g tot= %g\n",Na, chsq,tmp,chsq+tmp);
			chsq += tmp;
			if (chsq < chiBIC_min)
			  {
				// Improvement: save state
				Na_opt = Na;
				chiBIC_min = chsq;
				// save minima 
				if (_SSI.show_synth)
				  {
					System.arraycopy(wk,0,y_opt,0,N);
				  }
				else
				  {
					System.arraycopy(beta_full,0,beta_opt,0,N);
				  }
			  }
			chiBIC_last = chsq;
			// after climbing uphill for a bit, give up. Allow a little tolerance
			// chsq should scale as N, so make this, e.g. 0.1*N
			if (chsq > chiBIC_min + 0.1*N) 
			  {
				// Termination. Copy back in best solution to mu
				if (_SSI.show_synth) 
				  {
					System.arraycopy(y_opt,0,mu,0,N);
				  }
				else
				  {
					System.arraycopy(beta_opt,0,mu,0,N);
				  }
				do_break = true;
			  }
		  }
		if (profile) tT[4] += (System.nanoTime()-tT0[4]); // profiling

		// CLEANUP
		Xat=null;
		rhs=null;

		// break out, with solution in mu, either synthetics, or refl coeffs
		if (do_break) break;
	
	  }
	// END LARS
	if (!_SSI.show_synth) 
	  {
		// make true reflection coefficients (i.e. undo sign flip). mu[] contains the OLS beta[] at last model
		for(j=0;j<N;j++) mu[j] *= s[j];
	  }
	// Copy solution mu into Job variables for queue: FIXME stack1
	y = J.y[0]; // a reminder
	for(i=0;i<N;i++) y[i] = (float)mu[i];
	
	return;
  }

  double getC_hat(int N,float[] y,double[] mu,double[] c_hat,double[] wk,int[] s)
  {
	int i,j,imax=-1;
	double tmp,C_hat;

	boolean verb=false;

	for(i=0;i<N;i++) wk[i] = y[i] - mu[i];

	get_Xt_vec_prod(N,c_hat,wk);

	C_hat=0;
	for(i=0;i<N;i++)
	  {
		//		for(j=0,tmp=0;j<N;j++) tmp += Xt[i][j] * wk[j];
		// c_hat[i] = tmp;
		if ((tmp = abs(c_hat[i])) > C_hat)
		  {
			C_hat = tmp;
			imax=i;
		  }
		s[i] = (c_hat[i]>=0)? 1 : -1;
	  }
	if (verb) _SSI.EH.warn("Max at i="+imax);
	return C_hat;
  }

  int do_cholesky_update(int Na,DenseDoubleMatrix2D sG_A_L, double[][] sG_A, double[] rhs, double[] wk)
  {
	int i,j;
	double tmp, lsq;
	// double[] p1,p2;
	double[] p3;

	lsq=0;
	p3 = sG_A[Na-1];
	// p2 = sG_A_L[Na-1];
	for (i=0;i<(Na-1);i++)
	  {
		// p1 = sG_A_L[i];
		for (j=0,tmp=0; j<i ;j++) tmp += sG_A_L.getQuick(i,j) * sG_A_L.getQuick(Na-1,j); // tmp += p1[j] * p2[j];
		// sG_A_L[Na-1][i] = (sG_A[Na-1][i] - tmp)/p1[i];
		tmp = (p3[i] - tmp)/sG_A_L.getQuick(i,i);
		sG_A_L.setQuick(Na-1,i,tmp);
		lsq += tmp*tmp; // p2[i]*p2[i];
	  }
	lsq = sG_A[Na-1][Na-1] - lsq;

	if (lsq<0) 
	  {
		_SSI.EH.warn("Cholesky-update failed! Exiting");
		return 1;
	  }
	else 
	  {
		//sG_A_L[Na-1][Na-1] = sqrt(lsq);
		sG_A_L.setQuick(Na-1, Na-1, sqrt(lsq));
	  }
	// backsolve
	for (i=0;i<Na;i++) 
	  {
		// p1 = sG_A_L[i];
		for (tmp = rhs[i],j=i-1; j>=0; j--) tmp -= sG_A_L.getQuick(i,j)*wk[j];  // tmp -= p1[j]*wk[j];
		wk[i] = tmp/sG_A_L.getQuick(i,i); // tmp/p1[i];
	  }
	for (i=Na-1;i>=0;i--) 
	  {
		for (tmp = wk[i],j=i+1; j<Na; j++) tmp -= sG_A_L.getQuick(j,i)*wk[j];  //tmp -= sG_A_L[j][i]*wk[j];
		wk[i] = tmp/sG_A_L.getQuick(i,i); // tmp/sG_A_L[i][i];
	  }

	// overwrite into rhs
	for(i=0;i<Na;i++) rhs[i] = wk[i];
 	
	return 0;
  }

  void get_Xt_vec_prod(int N,double[] ans,double[] vec)
  {
	int i,j,jl,ju,p,c;
	double tmp;
	double[] X;

	double[][] Xt = _SSI.Xt;
	if (false)
	  {
		// full dense multiply
		for(i=0;i<N;i++)
		  {
			for(j=0,tmp=0;j<N;j++) tmp += Xt[i][j] * vec[j];
			ans[i] = tmp;
		  }
	  }
	else
	  {
		// band multiply
		// case of stack 1 only (FIXME)
		
		p = _SSI.sI[0].n_precursor; // wM->n_precursor; // < 0
		c = _SSI.sI[0].n_coda;      // wM->n_coda;      // > 0
		for (i=0;i<N;i++)
		  {
			jl = max(0,i-c);
			ju = min(N-1,i-p);
			X = Xt[i];
			for(j=jl,tmp=0; j<=ju; j++) tmp += X[j] * vec[j];
			ans[i] = tmp;
		  }
	  }
  }
}
