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
package au.csiro.JamesGunning.CSEM.plotting;

import java.io.*;
import java.nio.*;
import java.util.*;
import java.awt.*;

import cern.colt.map.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;  
import cern.colt.list.*;
import cern.jet.stat.*;

 

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.Numerical.DenseFloatMatrix2D;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.util.Util.*;
import au.csiro.JamesGunning.lib.graphutil.*;



/**
   Analysing code for MCMC or bootstrap outputs of
   Java Bayesian 1D CSEM code.

   @author James.Gunning@csiro.au
   @version 1.0, 13/3/2009
*/ 

public class OneDModelStatsPlotter
{
  ErrorHandler EH;
  int np=-1,ns,Nptot;
  double[] depths;
  double d_base=-1.0;
  ArrayList<double[]>data;
  String filename=null;

  boolean plotQuantiles=false;
  float quantile=50f;
  boolean plotModels = false;
  boolean normalVideo=false;
  boolean graphics = true;
  boolean asciiOut=false;
  CanvasPlotter CP;

  OneDModelStatsPlotter()
  {
	EH = new ErrorHandler();
  }

  public static void main(String argv[]) 
  {
	OneDModelStatsPlotter PL = new OneDModelStatsPlotter();
	ErrorHandler _EH = PL.EH;
	
	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("OneDModelStatsPlotter",_EH);
	options.addLicense(DeliveryLicense.license_text);
	options.addSpielText("Plots useful stuff from MCMCsamples or PBsamples files from CSEM1DInverter");
	
	options.addMandatoryArg(1,"Samples.txt","MCMC or Bootstrap samples file");	

	options.addOption("-Q",null,1,"float","quantile","plot matched q (0<q<100) quantiles from MCMC or Bootstrap ensemble");
	options.addOption("-m",null,0,null,null,"plot models from MCMC/Bootstrap ensemble");
	options.addOption("-N",null,0,null,null,"normal video (white background)");
	options.addOption("-O",null,0,null,null,"dump ascii versions of graphs to stdout");

	// parse options
	options.parseOptions(argv,true);

	PL.filename = options.getMandatoryArg(1);
	if (options.isOptionSet("-Q")) 
	  {
		PL.plotQuantiles = true;
		PL.quantile = ((Float)options.getOptionValues("-Q")[0]).floatValue();
		if (PL.quantile<=0 || PL.quantile>=100) _EH.die_err("Quantile must be in 0<q<100");
	  }
	if (options.isOptionSet("-m")) 
	  {
		PL.plotModels = true;
	  }
	if (options.isOptionSet("-N")) 
	  {
		PL.normalVideo = true;
	  }
	if (options.isOptionSet("-O")) 
	  {
		PL.asciiOut = true;
	  }

	PL.run();
  }


  void run()
  {
	int i,j,k,N;
	String s;
	BufferedReader in=null;
	String[] st;
	double[] v;
	double val;

	data = new ArrayList();

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
			if (s.startsWith("#"))
			  {
				if (s.startsWith("#DEPTHS")) 
				  {
					st = s.split("\\s+");
					depths = new double[st.length-1];
					for (i=1; i<st.length; i++)
					  {
						try { depths[i-1] = Double.parseDouble(st[i]);}
						catch (NumberFormatException e)
						  {
							EH.die_err("Depth line containing "+s+" contains non-number. ");					
						  }
					  } 
					continue;
				  }
				else if (s.startsWith("#BASE_DEPTH")) 
				  {
					st = s.split("\\s+");
					if (st.length>1)
					  {
						try {d_base = Double.parseDouble(st[1]);}
						catch (NumberFormatException e)
						  {
							EH.die_err("Depth line containing "+s+" contains non-number. ");					
						  }
					  }
					continue;
				  }
				else if (s.startsWith("#np= ")) 
				  {
					// We have np (number geol params)
					try {np = Integer.parseInt(s.replace("#np=","").trim());} 
					catch (NumberFormatException e) {EH.die_err("Line containing "+s+": can't parse for np");}
					continue; // skip other comments
				  }
				else if (s.startsWith("# m")) 
				  {
					// count # params
					st = s.split("\\s+");
					for (i=0,Nptot=0; i<st.length; i++) if (st[i].startsWith("m_")) Nptot++; 
					continue; // skip other comments
				  }
				else 
				  {
					continue; // skip other comments
				  }
			  }
			else
			  {
				st = s.split("\\s+");
				if (Nptot==0)
				  {
					Nptot = st.length;
					EH.warn("Warning: no header - guessing Nptot="+Nptot);
				  }
				double[] m = new double[Nptot];
				for (i=0; i<Nptot; i++)
				  {
					try
					  {
						m[i] = Double.parseDouble(st[i]);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Line containing "+s+" contains non-number.");
					  }
				  }
				data.add(m);
			  }
		  }
		in.close();
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open or close " + filename);
	  }
	ns = data.size();
	// set np
	if (np==-1) 
	  {
		// when #np = missing from header
		np = Nptot;
	  }
	// dummy depths if needed
	if (depths==null)
	  {
		depths = new double[np+1];
		for(i=0;i<np+1;i++) depths[i]=i;
	  }
	// Set up graphics
	ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter CP=null;
	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
	CP = new CanvasPlotter(true); // single kill needed for entire App.
	if (normalVideo) CP.reverseVideo=false;

	pg = new CanvasPlotter.TwoDPlotGroup(CP);
	Pgroups.add(pg);

	// placeholder for depth if not present
	if (d_base == -1.0) d_base = (np > 1) ? (2*depths[np-1] - depths[np-2]) : ( 1.1* depths[np-1]); 

	// sundry tasks
	if (plotQuantiles)
	  {
		// stuff data into DoubleArrayList s for each param
		DoubleArrayList[] trm = new DoubleArrayList[np];

		for(i=0;i<np;i++)
		  {
			trm[i] = new DoubleArrayList();
			DoubleArrayList _t =	trm[i];
			for(j=0;j<ns;j++)
			  {
				_t.add(((double[])(data.get(j)))[i]);
			  }
		  }

		double q;

		CanvasPlotter.TwoDPlot plot;
		pg.GG.Yword ="value";
		pg.GG.Xword = "depth";

		// line graph of "quantile"...
		plot = new CanvasPlotter.TwoDPlot();
		plot.x = new double[2*np];
		plot.y = new double[2*np];
		for (i=0,k=0;i<np;i++)
		  {
			DoubleArrayList dl = trm[i];
			dl.sort();
			q =Descriptive.quantile(dl,0.01*quantile);
			// EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",depths[i],q,(i!=(np-1)) ? depths[i+1] : d_base,q);
			plot.x[k] = depths[i];
			plot.x[k+1] = (i!=(np-1)) ? depths[i+1] : d_base;
			plot.y[k] =  plot.y[k+1] = q;
			if (asciiOut) EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",plot.x[k],q,plot.x[k+1],q);
			k+=2;
		  }
		pg.plots.add(plot);
		pg.GG.computeRanges(plot.x,plot.y);

		// Median (p50 quantile)
		plot = new CanvasPlotter.TwoDPlot();
		plot.c = Color.red;
		plot.x = new double[2*np];
		plot.y = new double[2*np];
		if (asciiOut) EH.printf("%n");
		for (i=0,k=0;i<np;i++)
		  {
			DoubleArrayList dl = trm[i];
			q=Descriptive.quantile(dl,0.5);
			// EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",depths[i],q,(i!=(np-1)) ? depths[i+1] : d_base,q);
			plot.x[k] = depths[i];
			plot.x[k+1] = (i!=(np-1)) ? depths[i+1] : d_base;
			plot.y[k] =  plot.y[k+1] = q;
			if (asciiOut) EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",plot.x[k],q,plot.x[k+1],q);
			k+=2;
		  }
		pg.plots.add(plot);
		pg.GG.computeRanges(plot.x,plot.y);

		// Opposite quantile i.e. (100%-quantile)
		plot = new CanvasPlotter.TwoDPlot();
		plot.x = new double[2*np];
		plot.y = new double[2*np];
		if (asciiOut) EH.printf("%n");
		for (i=0,k=0;i<np;i++) 
		  {
			DoubleArrayList dl = trm[i];
			q =Descriptive.quantile(dl,1.0-0.01*quantile);
			// EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",depths[i],q,(i!=(np-1)) ? depths[i+1] : d_base,q);
			plot.x[k] = depths[i];
			plot.x[k+1] = (i!=(np-1)) ? depths[i+1] : d_base;
			plot.y[k] =  plot.y[k+1] = q;
			if (asciiOut) EH.printf("%15.7g %15.7g%n%15.7g %15.7g%n",plot.x[k],q,plot.x[k+1],q);
			k+=2;
		  }
		pg.plots.add(plot);
		pg.GG.computeRanges(plot.x,plot.y);
		CP.windowTitle = "+/-"+quantile+"%-Quantiles: "+filename;
		if (asciiOut) 
		  {
			EH.printf("%n");
		  }
		else
		  {
			CP.run(Pgroups); 
		  }
	  }

	if (plotModels)
	  {
		// "Pixel" plots of models, vertically by depth, horizontally by "sample"
		// block colour informed by model value

		pg.GG.reverseY = true;
		pg.GG.Yword = "depth";
		pg.GG.Xword = "sample";
		
		pg.GG.setTickNumbersEngStyle(false,false);

		CanvasPlotter.TwoDShapePlot plotS = new CanvasPlotter.TwoDShapePlot();
	  
		double[] x,y,x1,y1,col;
		int[] gcol;

		N=np*ns;
		x = new double[N];
		y = new double[N];
		x1 = new double[N];
		y1 = new double[N];
		col = new double[N];
		gcol = new int[N];
		 
		// search for a safe colour value
		double goodCol=0;
		int blankGcol = CP.reverseVideo ? 0:1;
		HERE:
		for(j=0;j<ns;j++)
		  {
			v = data.get(j); // value of model in (i,j)th block
			for(i=0;i<np;i++)
			  {				
				if (!Double.isNaN(v[i])) 
				  {
					goodCol = v[i]; break HERE;
				  }
			  }
		  }
		// now stuff
		for(j=0,k=0;j<ns;j++)
		  {
			v = data.get(j); // value of model in (i,j)th block
			for(i=0;i<np;i++)
			  {
				x[k]=j;
				x1[k]=j+1;
				y[k]=depths[i];
				y1[k]= (i<np-1) ? depths[i+1] : d_base;
				col[k] = v[i];
				if (!Double.isNaN(v[i])) 
				  {
					gcol[k] = -1; // no "palette" override 
				  }
				else
				  {
					gcol[k] = blankGcol;
					col[k] = goodCol;
				  }
				k++;
			  }
		  }
		
		plotS.x = x;
		plotS.x1 = x1;
		plotS.y = y;
		plotS.y1 = y1;
		plotS.col = col;
		plotS.gcol=gcol;
		pg.GG.computeRanges(plotS.x,plotS.y);
		pg.GG.computeRanges(plotS.x1,plotS.y1);
		plotS.normalize(); // rescales colors to [0-1]
		plotS.setColorScheme(Color.black,Color.white);
		pg.plots.add(plotS);
		CP.windowTitle = "OneDModelStatsPlotter: "+filename;
		// render plot
		CP.run(Pgroups); 
	  }
  
	if (asciiOut) System.exit(0);
  }
}
