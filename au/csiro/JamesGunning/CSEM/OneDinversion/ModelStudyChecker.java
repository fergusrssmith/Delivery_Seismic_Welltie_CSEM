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
import au.csiro.JamesGunning.lib.Numerical.*;

import au.csiro.JamesGunning.lib.stat.util.MCMCConvergenceMonitor;

/**
   Code for assembling marginal distributions from 
   ModelStudy files.

   @author James.Gunning@csiro.au
   @version 1.0, 1/5/2009
*/

public class ModelStudyChecker implements Comparator<double[]>
{
  ErrorHandler EH = new ErrorHandler();
  int verbosity;
  boolean showMarg=false;
  boolean showQuantiles=false;
  float quantile=0;
  boolean checkCorr=false;

  String filename = "ModelStudy.txt";
  ArrayList<double[]>data = null;
  ArrayList<String>colnames = null;
  int Nstud=0;
  int rmsCol=0;
  int MMLcol=0;
  int rhoStartCol;

  int dataSortIndx=0;

  /** 
	  main driver from command line
  */
  public static void main(String argv[]) 
  {
	ModelStudyChecker MSC = new ModelStudyChecker();
	ErrorHandler _EH = MSC.EH;


	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("ModelStudyChecker",MSC.EH);
	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber);

   	options.addMandatoryArg(1,"ModelStudy.txt","model study file, or - for stdin");

	options.addOption("-v",null,1,"int"," verbosity","level of verboseness (adds extra data to marginal files too)");	
	options.addOption("-M",null,0,null,null,"Plot marginals");
	options.addOption("-Q",null,1,"float","quantile(%%)","Form quantiles from marginals");
	options.addOption("--GO",null,0,null,null,"dump graphics numbers to 2-column ascii stdout (i.e. suppress java graphics)");	
	options.addOption("-N",null,0,null,null,"Normal-video graphics (i.e. black-on-white)");

	// parse options
	options.parseOptions(argv,true);

	MSC.filename = options.getMandatoryArg(1);

	if (options.isOptionSet("-v")) 
	  {
		MSC.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();	
	  }
	if (options.isOptionSet("-M")) 
	  {
		MSC.showMarg=true;
	  }	
	if (options.isOptionSet("-Q")) 
	  {
		MSC.showQuantiles =true;
		MSC.quantile = ((Float)options.getOptionValues("-Q")[0]).floatValue();
		if (MSC.quantile<=0 || MSC.quantile>=100)
		  {
			_EH.die_err("Quantile "+ MSC.quantile + " must be in (0,100)");
		  }
	  }

	MSC.runMSC();
  }
  
  /**
	 Parse file for data format like....
	 # m_1 m_2 rho_1 rho_2 logMMLapprox rmsMisfit d_0 d_1 d_2 d_3 d_4 log10rho_0 log10rho_1 log10rho_2 log10rho_3 
	 columns.....

  */

  void getData()
  {
	int i,j,k,Ncol;
	boolean hdrSet=false;

	// real data
	BufferedReader in;
	String s;
	String[] st;
	double[] d;
	Nstud=0;
	Ncol=0;
	rmsCol=0;
	MMLcol=0;
	rhoStartCol=-1;

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
			if (s.startsWith("# m")) 
			  {
				// 2 sets of Nstud columns + MML + rms
				// e.g. # m_1 m_2 rho_1 rho_2 logMMLapprox rmsMisfit d_0 d_1 d_2 d_3 d_4 log10rho_0 log10rho_1 log10rho_2 log10rho_3 
				colnames = new ArrayList();
				st = s.split("\\s+");
				for (j=1,k=0;j<st.length;j++,k++)
				  {
					if (st[j].startsWith("m_")) 
					  {
						Nstud++;
						colnames.add(st[j]);
					  }
					else if (st[j].startsWith("rms")) rmsCol=k;
					else if (st[j].startsWith("logMML")) MMLcol=k;
					else if (st[j].startsWith("rho") && (rhoStartCol==-1)) rhoStartCol=k;
				  }
				Ncol = st.length-1; 
				data = new ArrayList<double[]>();

				if (Nstud<1) EH.die_err("File "+filename+" has no model m_ columns???");
				hdrSet = true;			   
			  }
			else if (s.startsWith("#")) continue; // skip comments
			else
			  {
				if (!hdrSet) EH.die_err("File "+filename+" has no header # m_....");

				st = s.split("\\s+");
				
				// format
				// # m_1 m_2 ..rho_1 rho_2 ...logMMLapprox rmsMisfit d_0 d_1 ... log10rho_0 log10rho_1 ...  
				
				if (st.length != Ncol) EH.die_err("File "+filename+" has incomplete line: "+s); 
				d = new double[Ncol];
				try
				  {
					for (j=0;j<Ncol;j++) d[j] = Double.parseDouble(st[j]);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Line containing "+s+" does not have numbers??");
				  }
				data.add(d);
			  }
		  }
		in.close();
		data.trimToSize();
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't find or use file "+filename);
	  }
  }

  /**
	 Comparator to sort by field data[dataSortIndx]
  */
  public int compare(double[] data1, double[] data2) 
  { 
	double x = data1[dataSortIndx];
	double y = data2[dataSortIndx];
	return (x<y)? -1: ((x>y)?1:0);
  }

  /**
	 Main routine to form marginals etc
  */
  void runMSC()
  {
	int i,j,k,pass,type,N,N_table;
	double min,sum=0,norm,xMin,xMax,cdf=0,eps;
	double x,x0=0,y,p,xp,q;
	double[] X=null,Y=null;
	PrintWriter out=null;
	String filenameO=null;
	LinearInterpolation LI;

 	// read samples...
	getData();


	if (showMarg || showQuantiles)
	  {
		double[][] arrs = (double[][]) data.toArray(new double[0][0]);
		N = arrs.length;
		// get min
		min = Double.MAX_VALUE;
		for (j=0;j<N;j++) min = Math.min(min,arrs[j][MMLcol]);
		// overwrite last col with scaled probs
		for (j=0;j<N;j++) 
		  {
			p = Math.exp(-(arrs[j][MMLcol]-min));
			arrs[j][rmsCol] = p;
		  } 
		// form all marginals
	
		double[][] castTarget = new double[0][0];
		for(i=0;i<Nstud;i++)
		  {
			dataSortIndx = i;
			// sort data so variable of interest  goes slowest
			Collections.sort(data,this);
			arrs = (double[][]) data.toArray(castTarget);
			if (false)
			  { 
				// show ordered data
				EH.printf("Resorted data:%n");
				for (j=0;j<N;j++) 
				  {
					for(k=0;k<Nstud;k++) EH.printf("%10g ", arrs[j][k]); 
					EH.printf("%n");
				  }
			  }
			if (showMarg)
			  {
				filenameO = (filename.equals("-") ? ("Marginals."+colnames.get(i)+".txt") : (filename+"_marginals."+colnames.get(i)+".txt"));
				try 
				  {
					out = new PrintWriter(new BufferedWriter(new FileWriter(filenameO)));
				  } 
				catch (IOException ioe) 
				  {
					throw new RuntimeException("couldn't open " + filenameO);
				  }
	
				out.printf("#Marginal %d%n",i);
			  }

			for (pass=0,norm=0,N_table=1; pass<=1; pass++)
			  {
				x0 = arrs[0][dataSortIndx];				
				if (showQuantiles && pass==1)
				  {
					X = new double[N_table];
					Y = new double[N_table];
					LI = new LinearInterpolation(X,Y,N_table,EH);
					X[0]=0;
					Y[0]=0;
				  }
				cdf=0;
				sum=0;
				for (j=0,k=1;j<N;j++) 
				  {
					x = arrs[j][dataSortIndx]; 
					if (x==x0)
					  {
						sum += arrs[j][rmsCol]; 
					  }
					else
					  {
						if (pass==0)
						  {
							norm +=sum;
							N_table++;
						  }
						else 
						  {
							cdf += sum/norm;
							if (showMarg) 
							  {
								out.printf("%g %g%n",x0,sum/norm); // PDF
								//out.printf("%g %g%n",x0,cdf); // CDF
							  }
							if (showQuantiles)
							  {
								X[k] = x0;
								Y[k]= cdf;
								k++;
							  }
						  }
						x0=x;
						sum=0;
					  }
				  }
			  }
			if (showMarg)
			  {
				// final point 
				out.printf("%g %g%n",x0,sum/norm); // CDF
				out.close();
				EH.println("Wrote: "+filenameO);
			  }
			if (showQuantiles)
			  {
				// CDF: for(k=0; k<N_table; k++) EH.println("%g %g",X[k],Y[k]);

				q = quantile/100.0;
				EH.printf("%g-Quantiles, %s: ",quantile, colnames.get(i));
				for(pass=0;pass<=1;pass++)
				  { // q, then 1-q
					for(k=0; k<N_table-1; k++)
					  {
						// dumb linear search: q
						if ((Y[k]<=q) && (Y[k+1]>q))
						  {
							x = X[k] + ((q-Y[k])/(Y[k+1]-Y[k]))*(X[k+1]-X[k]);
							EH.printf("%g ",x);
							break;
						  }					
					  }
					q = 1.0-q;
				  }
				EH.println("");
			  }
		  }

		// RTP and R histograms, assuming equal weight on all models
		for (type=0; type<=1; type++)
		  {
			String prefix = (type==0)? "R." : "RTP.";
			for(i=0;i<Nstud;i++)
			  {
				dataSortIndx = i+Nstud;      // d[dataSortIndx] = d[i+Nstud]  = resistivity
				double[] d;
				if (type==1)
				  {
					for (j=0;j<N;j++) 
					  {
						d = data.get(j);
						d[dataSortIndx] *= d[i]; // form (insitu) RTP = d[i] (thick) * d[i+Nstud] (res)
					  }
				  }
				// sort data so variable of interest  goes slowest
				Collections.sort(data,this);
				arrs = (double[][]) data.toArray(castTarget);
				if (showMarg)
				  {
					filenameO = (filename.equals("-") ? ("Marginal_"+prefix+colnames.get(i)+".txt") : 
								 (filename+"_marginal_"+prefix+colnames.get(i)+".txt"));
					try 
					  {
						out = new PrintWriter(new BufferedWriter(new FileWriter(filenameO)));
					  } 
					catch (IOException ioe) 
					  {
						throw new RuntimeException("couldn't open " + filenameO);
					  }


					if (verbosity>0)
					  {
						out.printf("#%s Marginal %d%n",prefix,i);
						for (j=0;j<N;j++) 
						  {
							// Show RTP, prob
							out.printf("%g %g%n",arrs[j][dataSortIndx],arrs[j][rmsCol]);
						  }
					  }
				  }


				// Assemble raw CDF
				if (showMarg && verbosity>0) out.printf("#Raw %scdf%n",prefix);
			
				xMin = arrs[0][dataSortIndx];
				xMax = arrs[N-1][dataSortIndx];
				norm = 0;
				N_table = 1;
				for(pass=0; pass<=1; pass++)
				  {
					if (pass==1)
					  {
						X = new double[N_table];
						Y = new double[N_table];
						LI = new LinearInterpolation(X,Y,N_table,EH);
						X[0]=xMin;
						Y[0]=0;
					  }
					cdf=0;
					x0 =  xMin;
					for (j=0,k=1,sum=0;j<N;j++)
					  {
						x = arrs[j][dataSortIndx]; // R, or RTP
						y = arrs[j][rmsCol];       // this is the scaled probability by this point ( ~ exp(-logMMML))
						if (x <= x0) sum += y;
						else
						  {
							// finish bin: print midpoint, fraction, setup new bin
							if (pass==0)
							  {
								norm+=sum;
								N_table++;
							  }
							else
							  {
								cdf += sum;
								X[k]=x;
								Y[k]=cdf/norm;
								if (showMarg && verbosity>0)  out.printf("%g %g%n",x,Y[k]);
								k++;
							  }
							sum = y;
						  }
					  }
				  }
				// X[], Y[] are now the "jittery" posterior marginal CDF of the RTP

				// Make smoothed histogram by splining CDF into Nbins equal chunks
				int Nbins = (int)(sqrt(N));
				double dX = (xMax-xMin)/Nbins;
				double[] xS = new double[Nbins];
				double[] yS = new double[Nbins];

				for(x=xMin,j=0,k=0; j<Nbins; j++,x+=dX)
				  {
					// dumb linear search to bracket x in CDF: 
					for(;k<N-1;k++)
					  {
						if ((X[k] <= x) && (X[k+1]>x))
						  {
							xS[j] = x;
							q = (X[k+1]>X[k]) ? (x - X[k])/(X[k+1]-X[k]) : 0;
							yS[j] = (1-q)*Y[k] + q*Y[k+1];
							// out.printf("SCDF %g %g%n",xS[j],yS[j]);
							break;
						  }
					  }
				  }
				// just to be sure:
				xS[Nbins-1]=xMax;
				yS[Nbins-1]=1.0;
			
				// Now plot PDF by numerical derivatives of splined CDF
				Spline sp = new Spline(xS,yS,Nbins,0,0,EH);
				// 100 points in curve
				dX = (xMax-xMin)/100.0;

				if (showMarg)	out.println("#smoothed histogram");
				for(x=xMin;x<xMax;x+=dX)
				  {
					eps = 0.0001*x;
					// numerical derivative
					xp = x+eps;
					if (xp>xMax) xp = x-eps;
					//y = sp.splint(xp,j!=0);  // CDF
					y = (sp.splint(xp,j!=0) - sp.splint(x, j!=0))/(xp-x); // PDF
					if (showMarg) out.printf("%g %g%n",x,max(y,0));
				  }

				if (showMarg)
				  {
					out.close();
					EH.println("Wrote: "+filenameO);
				  }


				if (showQuantiles)
				  {
					q = quantile/100.0;
					EH.printf("%g-Quantiles, %s%s: ",quantile,prefix,colnames.get(i));
					for(pass=0;pass<=1;pass++)
					  { // q, then 1-q
						for(k=0; k<N; k++)
						  {
							// dumb linear search: q
							if ((Y[k]<=q) && (Y[k+1]>q))
							  {
								x = X[k] + ((q-Y[k])/(Y[k+1]-Y[k]))*(X[k+1]-X[k]);
								EH.printf("%g ",x);
								break;
							  }					
						  }
						q = 1.0-q;
					  }
					EH.println("");
				  }
			  }
		  } // type
	  }
  }
}
