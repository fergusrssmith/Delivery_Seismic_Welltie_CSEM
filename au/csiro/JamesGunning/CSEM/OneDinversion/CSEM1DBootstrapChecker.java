/**********************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion
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
package au.csiro.JamesGunning.CSEM.OneDinversion;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;
import java.awt.*;

import cern.jet.stat.*;
import cern.colt.list.*;
import cern.jet.random.engine.*;
import hep.aida.*;
import hep.aida.ref.*;

import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.graphutil.*;

import au.csiro.JamesGunning.lib.stat.util.MCMCConvergenceMonitor;


/**
   Bootstrap checker. Code for plotting 
   -2log(posterior) sample fluctuation histograms from either
   Bootstrap or MCMC vs likely sampling distributions chisq_p
   known from linear theory.


   @author James.Gunning@csiro.au
   @version 1.0, 1/5/2009

*/

public class CSEM1DBootstrapChecker
{
  ErrorHandler EH = new ErrorHandler();
  int verbosity;
  boolean showHist=false;
  boolean checkCorr=false;
  boolean stdoutGraphics=false;
  boolean normalVideo=false;

  String filename = "PB_RMS.txt";
  ArrayList<PBSample> samples=null;
  String depths_header=null;
		  
  int nd=-1,np=-1,Nptot=-1,np_eff;
  float chisq_RMS,chsqNLPmin,chsqNLPminAct=-1f;
  double fac;

  // data from BS samples
  int Nsamp=0;
  double[] chsqTot;

  /** 
	  main driver from command line
  */
  public static void main(String argv[]) 
  {
	CSEM1DBootstrapChecker BC = new CSEM1DBootstrapChecker();
	ErrorHandler _EH = BC.EH;


	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("CSEM1DBootstrapChecker",BC.EH);
	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber);

   	options.addMandatoryArg(1,"PBsamples.txt","bootstrap or MCMC samples file, or - for stdin");

	options.addOption("-v",null,1,"int"," verbosity","level of verboseness");	
	options.addOption("-H",null,0,null,null,"Plot histogram of total misfit for samples, "
					  +"plus (linear-theory) offset-chisquare_np distribution");	
	options.addOption("--GO",null,0,null,null,"dump graphics numbers to 2-column ascii stdout (i.e. suppress java graphics)");	
	options.addOption("-N",null,0,null,null,"Normal-video graphics (i.e. black-on-white)");
	options.addOption("--Raftery",null,0,null,null,"show Raftery-style correlation statistics");

	// parse options
	options.parseOptions(argv,true);

	BC.filename = options.getMandatoryArg(1);

	if (options.isOptionSet("-v")) 
	  {
		BC.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();	
	  }
	if (options.isOptionSet("--GO")) BC.stdoutGraphics = true;
	if (options.isOptionSet("-N")) BC.normalVideo = true;

	if (options.isOptionSet("--Raftery")) BC.checkCorr=true;  

	if (options.isOptionSet("-H")) 
	  {
		BC.showHist=true;
	  }	

	BC.runBC();
  }
  
  /**
	 Main routine to read existing bootstrap samples and do stuff.
  */
  void runBC()
  {
	int i,j,k,io,i_theory,Ntheory;
	double x,xo=0,xOff;
	double mx=0,mn=0,peak=0;
	double[] m;
	double ChsqOffset;

   	// read samples...
	getData();


	if (showHist)
	  {
		// check fit against histogram
		DoubleArrayList Xl=null,Yl=null;
		CanvasPlotter.TwoDPlot plot = null;
		ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups=null;
		CanvasPlotter CP = null;
		CanvasPlotter.TwoDPlotGroup pg = null;

		if (!stdoutGraphics)
		  {
			Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
			CP = new CanvasPlotter(true);
			pg = new CanvasPlotter.TwoDPlotGroup(CP);
			Pgroups.add(pg);

			CP.windowTitle = "Sample misfit histogram";
			pg.GG.Xword = "-2 log(posterior) (prior & data)";
			pg.GG.Yword = "count";
			CP.reverseVideo = !normalVideo;
		  }

		// form histogram with some % extra range at either end (say 10%)
		// Ensure range is all positive too, *unless* we're examining reflection coeffs

		System.err.printf("Equivalent linearized target distribution: (offset=%-6g) chisq_(nu=%d)%n",chsqNLPmin,Nptot);

		double y,expMax,h,hMax=0;
		if (Nsamp>1)
		  {
			mx = mn = chsqTot[0];
			for (i=1;i<Nsamp;i++)
			  {
				mx = max( chsqTot[i],mx);
				mn = min( chsqTot[i],mn);
			  }
		


			IHistogram1D h1;
			int Nbins = Math.min(50,Math.max(4,(int)Math.rint(2*Math.sqrt(Nsamp))));

			h1 = new Histogram1D(null,Nbins,mn-0.1*(mx-mn),mx+0.1*(mx-mn));

			// populate histogram
			for (i=0; i < Nsamp; i++) h1.fill(chsqTot[i]);
			// Generate 2 column output: bin_centre, count_in_bin.
			IAxis ax = h1.xAxis();
	
			int nB = ax.bins();
			double[] X = new double[2*nB+2];
			double[] Y = (double[])X.clone();

			X[0] = ax.binLowerEdge(0);
			Y[0] = 0;
			h = hMax = 0;

			for (i=0,j=1; i < nB; i++,j++) 
			  {
				h = h1.binHeight(i);
				X[j] = ax.binLowerEdge(i);
				Y[j] = h;
				j++;
				X[j] = ax.binUpperEdge(i);
				Y[j] = h;
				hMax = max(hMax,h);
				if (stdoutGraphics) EH.printf("%g %g%n%g %g%n",X[j-1],Y[j-1],X[j],Y[j]);
			  }
			X[j] = ax.binUpperEdge(nB-1);     
			Y[j] = 0;
			EH.println("");

			if (!stdoutGraphics)
			  {
				plot = new CanvasPlotter.TwoDPlot();
				plot.x = X;
				plot.y = Y;
				pg.plots.add(plot);
				pg.GG.computeRanges(plot.x,plot.y);
			  }
		  }
		else
		  {
			// desperate defaults, just to draw sane distro from the header
			mx = chsqNLPmin + 3*Nptot;
			mn = chsqNLPmin;
			hMax=1;
		  }

		// and plot target distribution offset chisq(np) 
		// (1) theory, correct noise from DOF, (2) theory, actual (probably) wrong noise

		// get peak, for sane scaling/overflow prevention etc
		// mode of chisq distn is at x-2, for p>=2
		xOff = (Nptot>=2) ? Nptot-2 : 0;
		expMax = (Nptot/2.0 -1)*log(xOff)-xOff/2.0;

		Ntheory = 100;
		for (i_theory =0; i_theory<=1; i_theory++)
		  {
			if (i_theory==1 && chsqNLPminAct<0) break; // can't plot second theory curve		  {
			if (!stdoutGraphics)
			  {
				Xl = new DoubleArrayList(Ntheory);
				Yl = new DoubleArrayList(Ntheory);
			  }

			// now plot
			ChsqOffset = (i_theory==0) ? chsqNLPmin : chsqNLPminAct;
			for(i=0,x=mn-0.1*(mx-mn); x<=mx+0.1*(mx-mn) && i<Ntheory; x+=(mx-mn)/50.0, i++)
			  {
				xOff = x - ChsqOffset;
				h = ((Nptot/2.0 -1)*log(xOff)-xOff/2.0) - expMax;
				y = (x <= ChsqOffset)? 0 :  (hMax * exp(h));
				if (stdoutGraphics)	EH.printf("%g %g%n",x,y);
				if (!stdoutGraphics)
				  {
					Xl.add(x);
					Yl.add(y);
				  }
			  }
			if (stdoutGraphics)	EH.printf("%n");		  
			if (!stdoutGraphics)
			  {	
				Xl.trimToSize();
				Yl.trimToSize();
				plot = new CanvasPlotter.TwoDPlot();
				plot.x = Xl.elements();
				plot.y = Yl.elements();
				plot.c = (i_theory==0)? Color.red : Color.blue;
				pg.plots.add(plot);
				pg.GG.computeRanges(plot.x,plot.y);
			  }
		  }
		if (!stdoutGraphics)
		  {	
			CP.run(Pgroups); 
		  }
	  }


	if (checkCorr)
	  {
		// Check of correlation lengths in samples files (for all model columns)
		// using Raftery-style diagnostics: see notes in MCMCConvergenceMonitor class
		MCMCConvergenceMonitor CM;

		CM = new MCMCConvergenceMonitor(EH,2,false,true);

		CM.initialiseRafteryStyleConvergenceDiagnostics(Nptot,Nsamp);
		for(i=0;i<Nsamp;i++)
		  {
			m = samples.get(i).m;
			for(j=0;j<Nptot;j++) CM.setEstimand(j,m[j]); // value of model in (i,j)th block
			CM.completeScalarEstimandsUpdate(i);
		  }
		CM.finishRafteryStyleConvergenceDiagnostics();
		EH.printf("Suggested thinning rates:%nparam  rate%n");
		for(j=0;j<Nptot;j++) EH.printf("%3d    %d%n",(j+1),CM.Raftery_kthin[j]);
	  }
  }


  /**
	 Parse file for data format

	 #DEPTHS     915.60     965.60     1023.1     1089.2     1165.3     1252.7 ...
	 #nd= 138
	 #np= 18
	 #Nptot= 22
	 #error-rescaling-factor= xxx
	 # m_1 m_2 ...m_Nptot rms_misfit  chisq_data chisq_prior chisq_total rms_misfit_orig_noise rms_misfit_orig_noise_pre_dilution
	 ...ascii columns as per last line
  */

  void getData()
  {
	int i,j;
	PBSample samp;
	boolean checkedHdr=false;
	samples = new ArrayList<PBSample>();


	// real data
	BufferedReader in;
	String s;
	String[] st;

	samples = new ArrayList<PBSample>();

	try
	  {
		// use stdin
		if (filename.equals("-"))
		  {
			in = new BufferedReader(new FileReader(FileDescriptor.in));
		  }
		else
		  {
			in = new BufferedReader(new FileReader(filename));
		  }
		while ((s = in.readLine()) != null)
		  {
			s = s.trim();
			if (s.startsWith("#DEPTHS")) depths_header = s; 
			if (s.startsWith("#np=")) 
			  {
				try {np = Integer.parseInt(s.replace("#np=","").trim());} 
				catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for np");}
			  }
			else if (s.startsWith("#Nptot=")) 
			  {
				try {Nptot = Integer.parseInt(s.replace("#Nptot=","").trim());} 
				catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for np");}
			  }
			else if (s.startsWith("#nd=")) 
			  {
				try {nd = Integer.parseInt(s.replace("#nd=","").trim());} 
				catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for nd: ");}
			  }
			else if (s.startsWith("#expected-min-chisq_total="))
			  {
				try {chsqNLPmin = Float.parseFloat(s.replace("#expected-min-chisq_total=","").trim());} 
				catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for expected-min-chisq_total ");}
			  }
			else if (s.startsWith("#actual-min-chisq_total="))
			  {
				try {chsqNLPminAct = Float.parseFloat(s.replace("#actual-min-chisq_total=","").trim());} 
				catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for actual-min-chisq_total ");}
			  }
			else if (s.startsWith("#")) 
			  {
				continue; // skip comments
			  }

			st = s.split("\\s+");

			if (Nptot==-1) Nptot=np; // legacy formats

			// Expect format
			// # m_1 m_2 ... m_np rms_misfit  chisq_data chisq_prior chisq_total rms_misfit_orig_noise rms_misfit_orig_noise_pre_dilution

			if (st.length >= (Nptot+4)) 
			  {
				if (!checkedHdr)
				  {
					// must have these
					if (nd<0 || Nptot <0 || chsqNLPmin<0)
					  {
						EH.die_err("File "+filename+" missing header info of np, "
								   +"nd and expected-min-chisq_total. Can't proceed");
					  }
					else checkedHdr=true;
				  }
				samp = new PBSample(Nptot);
				try
				  {
					// next np entries are: model
					for (j=0;j<Nptot;j++) samp.m[j] = Double.parseDouble(st[j]);
					// then .....
					samp.chisq_RMS = Double.parseDouble(st[j++]);
					samp.chisq_data = Double.parseDouble(st[j++]);
					samp.chisq_prior = Double.parseDouble(st[j++]);
					samp.chisq_tot = Double.parseDouble(st[j]);
					// don't need rms_misfit_orig_noise
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Line containing "+s+" does not have numbers??");
				  }
				samples.add(samp);
			  }
		  }
		in.close();
		samples.trimToSize();
		Nsamp = samples.size();
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't find or use file "+filename);
	  }
	  
	chsqTot = new double[Nsamp];
	for(i=0;i<Nsamp;i++) 
	  {
		chsqTot[i] = samples.get(i).chisq_tot; 
		//EH.println("%g", chsqTot[i]);
	  }
  }
  
  static class PBSample
  {
	double[] m;
	double chisq_RMS, chisq_data, chisq_prior, chisq_tot;

	public PBSample(int n)
	{
	  m = new double[n];
	}
  }
}