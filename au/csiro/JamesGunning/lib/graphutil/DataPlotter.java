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

package au.csiro.JamesGunning.lib.graphutil;

import java.io.*;
import java.util.*;
import java.awt.*;

import cern.colt.map.*;
import cern.colt.list.*;
import hep.aida.*;
import hep.aida.ref.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 

/**
   A very simple class to show sundry scatterplots & histograms from
   multi-column ascii data using, optionally, plotutils' graph routine
   to do the plotting via system calls.  Some native java graphics
   also supplied.

   @author James.Gunning@csiro.au
   @version 1.0, 30/2/2008
*/

public class DataPlotter implements Comparator<float[]>
{
  public boolean USE_JAVA_PLOTTING = true;
  CanvasPlotter CP;
  ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;

  public static String grcmdDefault = 
  "graph -TX --blankout 1.05 -m 0 -F HersheySans --frame-color wheat "
	+"--bg-color black -w 0.7 -h 0.8 -r 0.18 -u 0.18 -C";
  String grcmdPlotSizeCmd = "--bitmap-size "; // append something like...900x700  
  String grcmdExtraOpts="";


  String grcmd = DataPlotter.grcmdDefault;
  String grout=null;
  ErrorHandler EH;
  String filenames;
  boolean colsReqd=false;
  boolean superpose=false;
  boolean doLabels=false;
  boolean normalVideo=false;
  boolean sameScale = false;
  boolean doHistograms = false;
  boolean doBlockImage = false;
  double bClip=Double.MAX_VALUE, wClip=-Double.MAX_VALUE;
  boolean bClipRel=false,wClipRel=false;
  String[] quantileFiles = null;
  ArrayList quantileList=null;
  int ncmax;
  int[] fileLengths;
  int[] columns;
  int[] linecols=null;
  String[] column_names=null;
  int propNamesLine=0;
  ArrayList d = new ArrayList();
  ArrayList<String> dLabels = new ArrayList<String>();
  String[] cnames=null;
  float data[][];
  String labels[];
  boolean noModel=false,noData=false;
  cern.colt.map.OpenDoubleIntHashMap fmap;
  double symSize=0.005;
  String fontAdj=null;
  float fontAdjS1=0f,fontAdjS2=0f;
  String grTitle=null;
  String Xlabel=null,Ylabel=null;
  boolean usePlotUtils=false;

  boolean engStyleNumbers=true;

  public DataPlotter()
	{
	  // handle errors: can change how these are handled if desired
	  EH = new ErrorHandler(); 
	}
  
  /** 
	  Example  test driver to show how to view extraction info from typical
	  MLgraphing-file.*.txt files produced by the WabeletExtractor.
	  The display could also be driven from some other class that has a
	  fully populated MLdata object, and then calling the display() method.
  */
  public static void main(String argv[]) 
  {
	int i,j;
	String arg;
	double x;
	String[] st;


	DataPlotter DP = new DataPlotter();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("DataPlotter",DP.EH);

	options.addSpielText("General plotter for multicolumn data in ascii format. "
						 +"Be sure to select columns (-c). Default mode is scatterplotting. "
						 +"Can also show histograms (-H) of columns, or block-raster plots (-B) of (x,y,value) triples from named columns. "
						 +"Block plots assume (x,y) values live on some regular, x,y aligned sort of grid. "
						 +"%nOptionally uses the GNU plotutils for graphics (default java native). The \"graph\" "
						 +"command from this package must be in your path to use this");

   	options.addMandatoryArg(1,"file(s).txt","Ascii columnated data file (# comments allowed). Use - for stdin");	
	options.addOption("-c",null,1,"String","c1,c2...","restrict to numbered columns c1,c2.. (e.g. 2,3,4)");
	options.addOption(null,"--names",1,"int","line-no.","Use property names on axes, parsed from supplied line in (first) file");
	options.addOption("-C",null,0,null,null,"Use first-line #name comment to label point-group from each file supplied (usually with -s, or just 2 columns)");
	options.addOption("-m",null,1,"String","col1,col2...","use named colors instead of auto-bumping (e.g. 1,1,3)");
	options.addOption("-s",null,0,null,null,"Superpose scatterplots, change symbols & colours for each dataset");
	options.addOption("-H",null,0,null,null,"Show histograms of named columns");
	options.addOption("-B",null,0,null,null,"Show block-model image of three supplied columns (taken as x,y,value). X,y must be on some sort of grid");
	options.addOption("--bclip",null,1,"String","value","in block-model (-B) images values outside bclip are clipped. Use r:val for relative clip");
	options.addOption("--wclip",null,1,"String","value","in block-model (-B) images values outside wclip are clipped. Use r:val for relative clip");

	options.addOption(null,"--HQ",1,"String","quantiles-filename","Superpose on histograms quantiles from file");
	options.addOption("-S",null,1,"float","size","Symbol size (default 0.005)");
	options.addOption(null,"--SS",0,null,null,"Render all scatterplots on same scales");
	options.addOption(null,"--title",1,"String","title","title of scatterplot graphs");
	options.addOption(null,"--fs",1,"String","s1,s2","font adjustments: scale axes fonts(s1), and title fonts(s2). Defaults=1.");
	options.addOption("--NE",null,0,null,null,"Not-Engineering style numbers (xxx.yy)");
	options.addOption("-N",null,0,null,null,"Normal-video (i.e. black-on-white)");
	options.addOption(null,"-X",1,"String","Xlabel","X-axis label");
	options.addOption(null,"-Y",1,"String","Ylabel","Y-axis label");
	options.addOption("--PU",null,0,null,null,"use plotutils, not java graphics");
	options.addOption("-g",null,1,"String","graph command","If using plotutils, user supplied graph "+
					  "command alternative to: "+DataPlotter.grcmdDefault);

	options.addOption("-o",null,1,"String","filename","stdout file to capture graph stdout if required (e.g. plot.ps), used with e.g. -e '-T ps'");
	options.addOption("-e",null,1,"String","options","extra options to graph (quote protect: e.g. -e '-T ps --frame-color black --bg-color white'");

	// parse options
	options.parseOptions(argv,true);

	DP.filenames = options.getMandatoryArg(1);

	if (options.isOptionSet("-s")) DP.superpose = true;
	if (options.isOptionSet("-C")) DP.doLabels = true;
	if (options.isOptionSet("-m")) 
	  {
		arg = options.getOptionString("-m",0);
		st = arg.split(",");
		DP.linecols = new int[st.length];
		for(i=0;i<st.length;i++)
		  {
			try
			  {
				DP.linecols[i] = Integer.parseInt(st[i]); 
			  }
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("-m "+arg +" contains a non-integer");
			  }
			if (DP.linecols[i]  <= 0)
			  {
				throw new RuntimeException("-m "+arg+" contains a non-positive color");
			  }
		  }												 
	  }
  

	if (options.isOptionSet("--PU")) DP.usePlotUtils = true;

	if (options.isOptionSet("-N")) DP.normalVideo = true;

	if (options.isOptionSet("-S")) 
	  {
		DP.symSize = (double)((Float)options.getOptionValues("-S")[0]).floatValue();
	  }
	if (options.isOptionSet("--SS")) DP.sameScale = true;

	if (options.isOptionSet("--title")) 
	  {
		DP.grTitle = (String)options.getOptionValues("--title")[0];
	  }
	if (options.isOptionSet("-X")) DP.Xlabel= (String)options.getOptionValues("-X")[0]; 
	if (options.isOptionSet("-Y")) DP.Ylabel= (String)options.getOptionValues("-Y")[0]; 

	if (options.isOptionSet("--fs")) 
	  {
		DP.fontAdj  = (String)options.getOptionValues("--fs")[0];
		st = DP.fontAdj.split("[,;:]");
		for (i=0; i< st.length; i++)
		  {
			try
			  {
				if (i==0) DP.fontAdjS1 = Float.parseFloat(st[i]);
				else DP.fontAdjS2 = Float.parseFloat(st[i]);
			  }
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("--fs "+DP.fontAdj +" contains a non-number");
			  }
		  }
	  }
	if (options.isOptionSet("--NE")) DP.engStyleNumbers = false;

	if (options.isOptionSet("-c")) 
	  {
		DP.colsReqd = true;
		arg = (String)options.getOptionValues("-c")[0];
		st = arg.split("[,;:]");
		DP.columns = new int[st.length];
		for (i=0; i< st.length; i++)
		  {
			try
			  {
				j = Integer.parseInt(st[i]);
			  }
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("-t "+arg+" contains a non-number");
			  }
			if (j <= 0)
			  {
				throw new RuntimeException("-t "+DP.columns+" contains no-positive type");
			  }
			DP.columns[i] = j;
		  }
	  }
	if (options.isOptionSet("--names")) 
	  {
		DP.propNamesLine = options.getOptionInt("--names",0);
		if (DP.propNamesLine <= 0)
		  {
			throw new RuntimeException("--names "+DP.propNamesLine+" must be +ve");
		  }
	  }
	if (options.isOptionSet("-g")) 
	  {
		DP.grcmd = (String)options.getOptionValues("-g")[0];
	  }
	if (options.isOptionSet("-o")) 
	  {
		DP.grout = (String)options.getOptionValues("-o")[0];
	  }
	if (options.isOptionSet("-e")) 
	  {
		DP.grcmdExtraOpts = (String)options.getOptionValues("-e")[0];
	  }
	if (options.isOptionSet("-H")) DP.doHistograms = true;
	if (options.isOptionSet("-B")) 
	  {
		DP.doBlockImage = true;

		if (options.isOptionSet("--bclip")) 
		  {
			arg = ((String)options.getOptionValues("--bclip")[0]);
			if (arg.startsWith("r:"))
			  {
				DP.bClipRel=true;
				arg = arg.replaceAll("r:","") ;
			  }
			try {DP.bClip = Float.parseFloat(arg);}
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("--bClip r:"+arg+" contains a non-number");
			  }			
		  }
		if (options.isOptionSet("--wclip")) 
		  {
			arg = ((String)options.getOptionValues("--wclip")[0]);
			if (arg.startsWith("r:"))
			  {
				DP.wClipRel=true;
				arg = arg.replaceAll("r:","") ;
			  }
			try {DP.wClip = Float.parseFloat(arg);}
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("--wClip r:"+arg+" contains a non-number");
			  }
		  }
	  }

	if (options.isOptionSet("--HQ")) 
	  {
		DP.quantileFiles = ((String)options.getOptionValues("--HQ")[0]).split("[,;:]");
	  }
	DP.run();
  }

  void run()
  {
	readData();
	determineColumns();
	if (doHistograms)
	  {
		plotHistograms();
	  }
	else if (doBlockImage)
	  {
		plotBlockXimageJava();
	  }
	else
	  {
		plotScatterPlots();
	  }
  }

  void determineColumns()
  {
	int ic;

	if (!colsReqd)
	  {
		// default: do all
		columns = new int[ncmax];
		for (ic=0; ic <columns.length; ic++) columns[ic] = ic+1; 
	  }
	else
	  {
		// check supplied columns
		for (ic=0; ic <columns.length;ic++) 
		  { 
			if (columns[ic] > ncmax) EH.die_err("Specified column larger than data set contains");
		  }
	  }

	if (!doHistograms && columns.length <=1 ) EH.die_err("Need more than one specified column");
	if (doBlockImage && (columns.length!=3))
	  {
		EH.die_err("Block plot: Must specify 3 columns only");
	  }

	if (cnames!=null)
	  {
		for (ic=0; ic <columns.length;ic++) 
		  { 
			if (columns[ic] > cnames.length) 
			  EH.die_err("--names points to file-line with insufficient names for requested columns!"); 
		  }
	  }
	  
  }

  void plotHistograms()
  {
	if (USE_JAVA_PLOTTING && !usePlotUtils) 
	  {
		plotHistogramsJava();
	  }
	else 
	  {
		plotHistogramsPlotutils();
	  }
  }


  void plotHistogramsJava()
  {
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDPlot plot = null,plotq=null;
	CanvasPlotter CP;
	ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups;


	int i,j,ic,jc,iCol,iF,nHistograms,nRows,nCols,nB;
	double h,hMax,hMaxTot=0;

	// read quantiles 
	readQuantiles();

	// number of histograms
	nHistograms = columns.length;

	nRows = (int)Math.floor(Math.sqrt(columns.length));
	nCols = (int)Math.ceil(columns.length/(double)nRows);

	CP = new CanvasPlotter(true);
	CP.reverseVideo = !normalVideo;

	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(nHistograms);

	int Nd = d.size();
	float[] Z = new float[Nd];
	for (iF=0; iF <columns.length;iF++) 
	  {
		ic = iF%nCols;
		jc = iF/nCols;

		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = nCols; 
		pg.GG.plotGroupNY = nRows; 
		pg.GG.plotGroupXindx = ic; 
		pg.GG.plotGroupYindx = jc;

		if (fontAdjS1 > 0) pg.GG.fontScale = fontAdjS1;
		if (fontAdjS2 > 0) pg.GG.fontTitleScale = fontAdjS2;
		  

		// copy data
		iCol = columns[iF]-1;
		for (i=0;i<Nd;i++) Z[i] = data[i][iCol];
		// sort
		java.util.Arrays.sort(Z); 
		// Histogram output
		// Let's take sqrt(2N) bins (same default rule as in |stat), 
		// but no less than 4, and no more than 50.
		int Nbins = Math.min(50,Math.max(4,(int)Math.rint(2*Math.sqrt(Nd))));
		double min = Z[0];
		double max = Z[Nd-1];
		// In wierd cases, if data is a spike (all constant), 
		// we need to ensure we have a non-zero range
		if (max==min) 
		  {
			if (max!=0.0) 
			  {
				max*=1.05;
				min/=1.05;
			  } 
			else
			  {
				max=1.0;
				min=-1.0;
			  }
		  }
		// form histogram with some % extra range at either end (say 10%)
		// Ensure range is all positive too, *unless* we're examining reflection coeffs
		IHistogram1D h1;
		double dx = (max-min)/Nbins;
		int nExtra = (int)Math.max(1,Math.rint(0.1*Nbins));
		Nbins += 2*nExtra;
		h1 = new Histogram1D(null,Nbins,min-nExtra*dx,max+nExtra*dx);
		// populate histogram
		for (i=0; i < Nd; i++) h1.fill(Z[i]);
		// Generate 2 column output: bin_centre, count_in_bin.
		IAxis ax = h1.xAxis();
		
		nB = ax.bins();
		double[] X = new double[2*nB+2];
		double[] Y = (double[])X.clone();

		X[0] = ax.binLowerEdge(0);
		Y[0] = 0;
		hMax = 0;
		for (i=0,j=1; i < nB; i++,j++) 
		  {
			h = h1.binHeight(i);
			X[j] = ax.binLowerEdge(i);
			Y[j] = h;
			j++;
			X[j] = ax.binUpperEdge(i);
			Y[j] = h;
			hMax = Math.max(hMax,h);
			//EH.printf("%g %g%n%g %g%n",X[j-1],Y[j-1],X[j],Y[j]);
		  }
		hMaxTot = Math.max(hMax,hMaxTot); // overall max range
 
		X[j] = ax.binUpperEdge(nB-1);     
		Y[j] = 0;

		// plot quantiles over the top
		if (quantileFiles != null)
		  {
			for (j=0;j<quantileList.size();j++)
			  {
				float[] v = (float[]) (quantileList.get(j));
				plotq = new CanvasPlotter.TwoDPlot();
				plotq.x = new double[2];
				plotq.y = new double[2];
				plotq.x[0] = plotq.x[1] = v[iF];
				plotq.y[0] = 0; plotq.y[1] = hMax/2.0;
				plotq.c = Color.green;
				pg.plots.add(plotq);
				//out.printf("#m=2,S=0%n %10.8g 0 %10.8g %10.8g%n",v[iF],v[iF],hMax/2.0);
				//out.println("");
			  }			
		  }
		//out.close();

		plot = new CanvasPlotter.TwoDPlot();
		plot.x = X;
		plot.y = Y;
		// color
		plot.c = Color.red;
		pg.plots.add(plot);
		pg.GG.computeRanges(plot.x,plot.y);
		if (cnames==null)
		  {
			pg.GG.Xword = Xlabel;
			pg.GG.title = ((grTitle!=null) ? grTitle: "col")+"_"+columns[iF];
		  }
		else
		  {
			pg.GG.Xword = cnames[iCol];
			pg.GG.title = (grTitle!=null) ? grTitle: pg.GG.Xword;
		  }
		pg.GG.Yword = Ylabel;		 
	  	pg.GG.setTickNumbersEngStyle(engStyleNumbers,engStyleNumbers);
	  }
	if (sameScale)
	  {
		for (i=0;i<Pgroups.size();i++)
		  {
			(Pgroups.get(i)).GG.setYRange(hMaxTot,0.0);
		  }
	  }
	CP.windowTitle = "Histograms";
	CP.run(Pgroups); 
  }

  /**
	 Plots a blocky raster plot of the values in c3=column[2], given
	 block-centres (x,y) read from c1=columns[0] and c2=columns[1]. Quite
	 like SU's ximage blockinterp=1, but using ascii data.
	 Assumes columns c1,c2 are on a strict (maybe stretched) grid, i.e.
	 does not plot general spatial samples (a TODO)
  */
  void plotBlockXimageJava()
  {
	int i,j,k;
	double[] x,y,x1,y1,col;
	int[] gcol;

	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDShapePlot plot = null;
	CanvasPlotter CP;
	ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups;


	if (columns.length!=3)
	  {
		EH.die_err("Block plot: Must specify 3 columns only");
	  }

	CP = new CanvasPlotter(true);
	CP.reverseVideo = !normalVideo;

	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
	pg = new CanvasPlotter.TwoDPlotGroup(CP);
	Pgroups.add(pg);

	int Nd = d.size();
	
	// sort data into ascending x, y (columns 1,2)

	java.util.Arrays.sort(data,this);

	plot = new CanvasPlotter.TwoDShapePlot();
	x = new double[Nd];
	y = new double[Nd];
	x1 = new double[Nd];
	y1 = new double[Nd];
	col = new double[Nd];
	gcol = new int[Nd];

	int blankGcol = CP.reverseVideo ? 0:1;
	int iX = columns[0]-1;
	int iY = columns[1]-1;
	int iZ = columns[2]-1; // color from "col 3"
	double wX,wY;

	// form x & y ordered arrays of block centres
	DoubleArrayList xList = new DoubleArrayList();
	DoubleArrayList yList =  new DoubleArrayList();
	double xL,yL;
	xL = yL = -Double.MAX_VALUE;
	for(i=0;i<Nd;i++)
	  {
		if ((wX = data[i][iX]) > xL)
		  {
			xList.add(wX);
			xL = wX;
		  }
		if ((wY = data[i][iY]) > yL)
		  {
			yList.add(wY);
			yL = wY;
		  }		
	  }
	xList.trimToSize();
	yList.trimToSize();
	double[] xLe = xList.elements();
	double[] yLe = yList.elements();
	int Nxlist = xLe.length-1;
	int Nylist = yLe.length-1;

	double defaultXw=1, defaultYw=1,defaultCol=0;
	if (xLe.length==1 && (yLe.length>1))
	  {
		defaultXw = (yLe[Nylist] - yLe[0])/yLe.length;
	  }
	else if (yLe.length==1 && (xLe.length>1))
	  {
		defaultYw = (xLe[Nxlist] - xLe[0])/xLe.length;
	  }


	float cMax=1,cMin=0;
	if (wClipRel || bClipRel)
	  {
		cMax = cMin = data[0][iZ];
		for(i=1; i<Nd; i++)
		  {
			cMax = Math.max(cMax,data[i][iZ]);
			cMin = Math.min(cMin,data[i][iZ]);
		  }

		if (bClipRel) 
		  {
			bClip = cMin + Math.abs(bClip);
		  }
		else  if (wClipRel) 
		  {
			wClip = cMax - Math.abs(wClip);
		  }
	  }
	if (bClip != Double.MAX_VALUE) defaultCol = bClip;
	else if (wClip != -Double.MAX_VALUE) defaultCol = wClip;

	System.err.printf("bClip=%6g wClip=%6g%n",bClip,wClip);

	for(i=0,j=0,k=0;i<Nd;i++)
	  {
		xL = data[i][iX];
		// search for edge of box to right. xL should appear in sorted increasing order
		for(;;j++) {if (xL == xLe[j]) break;}
		if (Nxlist>1)
		  {
			wX = (j != Nxlist) ?  (xLe[j+1] - xLe[j]) :  (xLe[j] - xLe[j-1]);
		  }
		else wX = defaultXw;

		x[i] = xL-wX/2.0;
		x1[i] = xL+wX/2.0;

		yL = data[i][iY];

		// search. start from 0 if overrun..
		if (yL==yLe[0]) k=0;
		for(; k<=Nylist; k++) {if (yL == yLe[k]) break;}
		if (Nylist>1)
		  {
			wY = (k != Nylist) ?  (yLe[k+1] - yLe[k]) :  (yLe[k] - yLe[k-1]);
		  }
		else wY = defaultYw;
		y[i] = yL+wY/2.0;
		y1[i] = yL-wY/2.0;

		// color..
		col[i] = data[i][iZ];

		if ((col[i] > wClip) && (col[i] < bClip))
		  {
			gcol[i] = -1;      // no "override": use true colour
		  }
		else
		  {
			col[i] = defaultCol; // dummy value, but in range
			gcol[k] = blankGcol;
		  }
	  }
	plot.x = x;
	plot.x1 = x1;
	plot.y = y;
	plot.y1 = y1;
	plot.col = col;
	plot.gcol=gcol;

	if (fontAdjS1 > 0) pg.GG.fontScale = fontAdjS1;
	if (fontAdjS2 > 0) pg.GG.fontTitleScale = fontAdjS2;
		  
	pg.GG.computeRanges(plot.x,plot.y);
	pg.GG.computeRanges(plot.x1,plot.y1);
	plot.normalize(); // rescales colors to [0-1]
	plot.setColorScheme(Color.white , Color.black);

	pg.plots.add(plot);
	pg.GG.computeRanges(plot.x,plot.y);
	pg.GG.title = (grTitle!=null) ? grTitle : ( "density,col_"+columns[2]+" vs cols "+columns[0]+","+columns[1]); 
	pg.GG.setTickNumbersEngStyle(engStyleNumbers,engStyleNumbers);
	pg.GG.Xword = Xlabel;
	pg.GG.Yword = Ylabel;

	CP.windowTitle = filenames;
	CP.run(Pgroups); 
  }

  void plotHistogramsPlotutils()
  {
	int i,j,nRows,nCols,nHistograms,iCol,ic,jc,iF,nB;
	float[] v;
	double h,hMax;
	double DXG=1.0,DYG=1.0,XSIZE,fSize;
	String name="",file_commands="";
	File[] temp;
	PrintWriter out=null; 

	// read quantiles 
	readQuantiles();

	// number of histograms
	nHistograms = columns.length;
	temp = new File[nHistograms];


	nRows = (int)Math.floor(Math.sqrt(columns.length));
	nCols = (int)Math.ceil(columns.length/(double)nRows);
	EH.println("CR="+nRows+" "+nCols);
	DXG = 1.0/(nCols);
	DYG = 1.0/(nRows);
	XSIZE = Math.min(DXG,DYG);

	// aim for constant font size...
	fSize = (columns.length==1)? 0.03 :  Math.min(1.0,0.02*Math.max(nCols,nRows));
	String fontSize = " -f "+fSize;

	int Nd = d.size();
	float[] Z = new float[Nd];
	for (iF=0; iF <columns.length;iF++) 
	  {
		ic = iF%nCols;
		jc = iF/nCols;
		double xGshift = ic*XSIZE;
		double yGshift = jc*XSIZE;
		// Open tmp files
		try 
		  {
			name = new String("junk"+iF);
			temp[iF] = File.createTempFile(name,".txt");				
			// repositioning
			file_commands = file_commands.concat(" --reposition "+xGshift+" "+yGshift+" "+XSIZE+" ");
			// labels
			file_commands = file_commands.concat(" -X col_"+ columns[iF]+ " ");
			// data
			file_commands = file_commands.concat(temp[iF].getAbsolutePath()+" ");
			
			temp[iF].deleteOnExit();
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + name);
		  }
		try 
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(temp[iF])));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + name);
		  }
		iCol = columns[iF]-1;
		// EH.println("Doing cols="+iFol+" "+jCol+": "+file_commands);
		// copy data
		for (i=0;i<Nd;i++) Z[i] = data[i][iCol];
		// sort
		java.util.Arrays.sort(Z); 
		// Histogram output
		// Let's take sqrt(2N) bins (same default rule as in |stat), 
		// but no less than 4, and no more than 50.
		int Nbins = Math.min(50,Math.max(4,(int)Math.rint(2*Math.sqrt(Nd))));
		double min = Z[0];
		double max = Z[Nd-1];
		// In wierd cases, if data is a spike (all constant), 
		// we need to ensure we have a non-zero range
		if (max==min) 
		  {
			if (max!=0.0) 
			  {
				max*=1.05;
				min/=1.05;
			  } 
			else
			  {
				max=1.0;
				min=-1.0;
			  }
		  }
		// form histogram with some % extra range at either end (say 10%)
		// Ensure range is all positive too, *unless* we're examining reflection coeffs
		IHistogram1D h1;
		h1 = new Histogram1D(null,Nbins,min-0.1*(max-min),max+0.1*(max-min));

		// populate histogram
		for (i=0; i < Nd; i++) h1.fill(Z[i]);
		// Generate 2 column output: bin_centre, count_in_bin.
		IAxis ax = h1.xAxis();
		nB = ax.bins();

		out.printf("%10.8g 0%n",ax.binLowerEdge(0));
		hMax = 0;
		for (i=0; i < nB; i++) 
		  {
			h = h1.binHeight(i);
			out.printf("%10.8g %10.8g%n%10.8g %10.8g%n",ax.binLowerEdge(i), h ,ax.binUpperEdge(i), h);
			hMax = Math.max(hMax,h);
		  }
		out.printf("%10.8g 0%n",ax.binUpperEdge(0));
		out.println("");
		// plot quantiles over the top
		if (quantileFiles != null)
		  {
			for (j=0;j<quantileList.size();j++)
			  {
				v = (float[]) (quantileList.get(j));
				out.printf("#m=2,S=0%n %10.8g 0 %10.8g %10.8g%n",v[iF],v[iF],hMax/2.0);
				out.println("");
			  }			
		  }
		out.close();
	  }
	int size = (superpose || (columns.length <= 2)) ? 700 : 900;
	i = Math.max(nRows,nCols);
	name = String.format("%s %s %dx%d %s %s -m 1 -S 0 %s",grcmd,grcmdPlotSizeCmd,size,size,
						 // size*nCols/i,size*nRows/i, 
						 fontSize, grcmdExtraOpts,file_commands);
 	callExternalPlot(name);
  }


  void plotScatterPlots()
  {
	if (USE_JAVA_PLOTTING && !usePlotUtils) 
	  {
		plotScatterPlotsJava();
	  }
	else 
	  {
		plotScatterPlotsPlotutils();
	  }
  }

  void plotScatterPlotsJava()
  {
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDPlot plot = null,plotq=null;
	CanvasPlotter CP;
	ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups;

	int i,j,k,ic,jc,nr,iCol,jCol,iF,jF,nFiles,IC;
	float[] v=null;

	// scatter plot matrix will have nRows = (columns.length-1) rows
	int nRows = columns.length-1;
	nFiles = ((nRows) * columns.length)/2;

	CP = new CanvasPlotter(true);
	CP.reverseVideo = !normalVideo;
	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(superpose ? 1 : nFiles);


	double xmin = Double.MAX_VALUE;
	double ymin = Double.MAX_VALUE;
	double xmax = -Double.MAX_VALUE;
	double ymax = -Double.MAX_VALUE;

	iF=0;
	for (ic=1; ic <columns.length;ic++) 
	  {
		for (jc=0; jc <ic; jc++) 
		  {
			iCol = columns[ic]-1;
			jCol = columns[jc]-1;

			if ((iF==0) || !superpose) pg = new CanvasPlotter.TwoDPlotGroup(CP);
			Pgroups.add(pg);
			pg.GG.plotGroupNX = superpose ? 1 : nRows; 
			pg.GG.plotGroupNY = superpose ? 1 : nRows; 
			pg.GG.plotGroupXindx = superpose ? 0 : jc;
			pg.GG.plotGroupYindx = superpose ? 0 : (ic-1);
			if (fontAdjS1 > 0) pg.GG.fontScale = fontAdjS1;
			if (fontAdjS2 > 0) pg.GG.fontTitleScale = fontAdjS2;

			if (grTitle!=null) pg.GG.title = grTitle+"_"+(iCol+1)+"_vs_"+grTitle+"_"+(jCol+1);
			if (engStyleNumbers==false) pg.GG.setTickNumbersEngStyle(false,false);
			
			pg.GG.Xword = Xlabel;
			pg.GG.Yword = Ylabel;


			nr = d.size();
			int i_file = 0;

			// new data set for each input file (as well as each crossplot) also
			for (jF=0,k=0; jF < fileLengths.length; jF++,iF++)
			  {
				int Nd = fileLengths[jF];
				double[] X = new double[Nd];
				double[] Y = new double[Nd];
				for (i=0; i< Nd; i++,k++)
				  {
					v = data[k];
					X[i] = v[jCol];
					Y[i] = v[iCol];
				  }
				plot = new CanvasPlotter.TwoDPlot();
				plot.x = X;
				plot.y = Y;
				plot.drawline = false;
				plot.symbol='r';
				if (doLabels) 
				  {
					plot.label=labels[jF];
				  }
				if (superpose)
				  {
					if (linecols!=null)
					  {
						IC = (iF >= linecols.length)? (linecols[linecols.length-1]-1): (linecols[iF]-1);
						plot.c = CanvasPlotter.defaultColors[IC%(CanvasPlotter.defaultColors.length)];
					  }
					else IC=iF;
				  }
				else
				  {
					// change cols each file only
					if (linecols!=null)
					  {
						IC = (jF >= linecols.length)? (linecols[linecols.length-1]-1): (linecols[jF]-1);
					  }
					else IC=jF;
				  }
				plot.c = CanvasPlotter.defaultColors[IC%(CanvasPlotter.defaultColors.length)];
				plot.symSize = symSize;
				pg.plots.add(plot);
				pg.GG.computeRanges(plot.x,plot.y);
				if (sameScale)
				  {
					xmin = Math.min(pg.GG.xMin,xmin);
					ymin = Math.min(pg.GG.yMin,ymin);
					xmax = Math.max(pg.GG.xMax,xmax);
					ymax = Math.max(pg.GG.yMax,ymax);
				  }
				if (cnames!=null)
				  {
					pg.GG.Xword = cnames[jCol];
					pg.GG.Yword = cnames[iCol];
				  }
			  }
		  }
	  }
	if (sameScale)
	  {
		for (i=0;i<Pgroups.size();i++)
		  {
			pg = Pgroups.get(i);
			for (j=0;j<pg.plots.size();j++)
			  {
				plot = (CanvasPlotter.TwoDPlot)pg.plots.get(j);
				pg.GG.xMin = xmin;
				pg.GG.yMin = ymin;
				pg.GG.xMax = xmax;
				pg.GG.yMax = ymax;
				pg.GG.updateTicksBracket();
			  }
		  }
	  }

	CP.windowTitle = "Scatterplots";
	CP.run(Pgroups); 	
  }


  void plotScatterPlotsPlotutils()
  {
	//	plot
	int i,j,ic,jc,nr,iCol,jCol,iF,nFiles,symIndx;
	boolean allow=true;
	float[] v;
	float x,y;
	double DXG=0.0;
	String file_commands="",fontSize="",name=null;
	File[] temp;
	PrintWriter out; 

	// scatter plot matrix will have nRows = (columns.length-1) scatterplots
	int nRows = columns.length-1;
	nFiles = ((nRows) * columns.length)/2;
	temp = new File[nFiles];

	// aim for constant font size...
	fontSize = " -f "+Math.min(1.0,(0.02*(nRows)));


	symIndx=1;
	if (superpose) 
	  {
		DXG = 0.0;
		fontSize = " -f 0.03 ";
	  } 
	else
	  {
		DXG = 1.0/(nRows);
		// aim for constant font size...
		fontSize = " -f "+Math.min(1.0,(0.02*(nRows)));
	  }

	iF=0;
	for (ic=1; ic <columns.length;ic++) 
	  {
		for (jc=0; jc <ic; jc++,iF++) 
		  {
			if (!superpose) symIndx = 1;

			double xGshift = jc*DXG;
			double yGshift = (1.0-DXG) -(ic-1)*DXG;
			// Open tmp files
			try 
			  {
				name = new String("junk"+ic+"."+jc);
				temp[iF] = File.createTempFile(name,".txt");				
				if (!superpose) file_commands = file_commands.concat(" --reposition "+xGshift+" "+yGshift+" "+DXG+
																	 " -X col_"+columns[jc]+" -Y col_"+columns[ic]+" ");
				file_commands = file_commands.concat(temp[iF].getAbsolutePath()+" ");

				temp[iF].deleteOnExit();
			  } 
			catch (IOException ioe) 
			  {
				throw new RuntimeException("couldn't open " + name);
			  }
			try 
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(temp[iF])));
			  } 
			catch (IOException ioe) 
			  {
				throw new RuntimeException("couldn't open " + name);
			  }
			iCol = columns[ic]-1;
			jCol = columns[jc]-1;
			//EH.println("Doing cols="+iCol+" "+jCol+": "+file_commands);
			
			if (superpose)
			  {
				out.printf("#m=%d,S=%d%n",-symIndx,symIndx);
				symIndx++;
			  }
			nr = d.size();
			int i_file = 0;
			for (i=0; i< nr; i++)
			  {
				// change symbol for new files...
				if ((fileLengths[i_file] == i) && !superpose)
				  {
					symIndx++;
					out.printf("#m=%d,S=%d%n",-symIndx,symIndx);
					i_file++;
				  }
				v = data[i];
				x = v[jCol];
				y = v[iCol];
				out.printf("%f %f%n",x,y);
			  }
			out.close();
			
			
		  }
	  }

  }

  void callExternalPlot(String name)
  {
	PrintWriter out;
	try
	  {
		EH.println("Attempting: "+name);
		Runtime runt = Runtime.getRuntime();
		Process p = runt.exec(name); 
		// Process p = runt.exec("uname"); 
		if (grout != null)
		  {
			EH.println("Attempting to capture in: "+grout);
			try 
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(grout)));
			  } 
			catch (IOException ioe)
			  {
				throw new RuntimeException("couldn't open " + grout);
			  }
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			while ((line = br.readLine()) != null) 
			  {
				out.println(line);
			  }
			out.close();
		  }
		// clean up
		try 
		  {
			p.waitFor();
		  }
		catch (InterruptedException ie) 
		  {
			throw new RuntimeException("Got interrupted!");
		  }
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't exec command "+name);
	  }
  }


  void readQuantiles()
  {
	int i,j,ic,nc,lineNo,dataLineNo,maxNumberedColumn=0;
	float[] v;
	String filename,s,word;
	BufferedReader in=null;

	if (quantileFiles == null) return;

	for (ic=0; ic <columns.length; ic++) 
	  {
		maxNumberedColumn = Math.max(maxNumberedColumn, columns[ic]);
	  }
	quantileList = new ArrayList();
	for (i=0; i<quantileFiles.length; i++)
	  {
		filename = quantileFiles[i];
		try
		  {
			in = new BufferedReader(new FileReader(filename));
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't open file "+filename);
		  }
		lineNo=0;
		dataLineNo=0;
		try
		  {
			while ((s=in.readLine()) != null)
			  {
				lineNo++;
				j=0;
				boolean isData = true;
				StringTokenizer st = new StringTokenizer(s);
				nc = st.countTokens();
				if (nc < maxNumberedColumn)
				  {
					EH.die_err("Quantile file has no column large enough to match requested column");
				  }
				v = new float[nc];
				while (st.hasMoreTokens()) 
				  {
					word = st.nextToken();
					if ((j==0) && word.startsWith("#")) isData = false; // skip comments
					if (isData)
					  {
						try
						  {
							v[j] = Float.parseFloat(word);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
						  }			
					  }
					j++;	
				  }
				if (isData)
				  {
					float[] quantiles = new float[columns.length];
					for (ic=0; ic <columns.length; ic++) 
					  {
						quantiles[ic] = v[columns[ic]-1];
					  }
					quantileList.add(quantiles);
				  }
			  }
			in.close();
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't close "+filename);
		  }
	  }
  }


  void readData()
  {
	int i,j,col,nc,lineNo,dataLineNo;
	String s,word,filename;
	String[] fn;
	StringTokenizer st;
	float[] v;
	ArrayList fd;

	BufferedReader in=null;
	fmap = new OpenDoubleIntHashMap();
	col=1;


	ncmax=0;
	// loop over files
	fn = filenames.split(",");
	fileLengths = new int[fn.length];
	for (i=0; i< fn.length; i++)
	  {
		filename = fn[i];
		try
		  {
			if (filename.equals("-"))
			  {
				// use stdin
				in = new BufferedReader(new FileReader(FileDescriptor.in));
			  }
			else
			  {
				in = new BufferedReader(new FileReader(filename));
			  }
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't open file "+filename);
		  }
		fd = new ArrayList();
		lineNo=0;
		dataLineNo=0;
		try
		  {
			while ((s=in.readLine()) != null)
			  {
				lineNo++;
				j=0;
				boolean isData = true;
				st = new StringTokenizer(s);
				nc = st.countTokens();
				if (lineNo==1)
				  {
					if (s.startsWith("#"))
					  dLabels.add(s.substring(1));
					else 
					  dLabels.add(null);
				  }
				if (i==0 && lineNo==propNamesLine)
				  {
					// pull property names off first file
					if (s.startsWith("#")) word = s.substring(1);
					else word=s;
					cnames = word.trim().split("\\s+");
					// take to be comments also
					isData = false; // 
				  }
				if (nc==0 || s.startsWith("#")) isData = false; // skip blank lines
				if (isData) ncmax = Math.max(nc,ncmax);           // track max no columns
				v = new float[nc];
				while (st.hasMoreTokens()) 
				  {
					word = st.nextToken();
					if (isData)
					  {
						try
						  {
							v[j] = Float.parseFloat(word);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
						  }			
					  }
					j++;	
				  }
				if (isData) 
				  {
					dataLineNo++;
					fd.add(v);
				  }
			  }
			in.close();
			fileLengths[i] = dataLineNo;
			// append to master collection
			d.addAll(fd);
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't close "+filename);
		  }
	  } // loop over files

	// dump into matrix
	data = new float[d.size()][];
	for (i=0;i<d.size();i++)
	  {
		v = (float[]) d.get(i);
		data[i] = v;
	  }
	labels = (String[])dLabels.toArray(new String[0]);
  }


  // comparator for data: sort by columns 1 && 2. 
  // called only when these are guaranteed to exist
  public int compare(float a[], float b[])
  {
	int iX = columns[0]-1;
	int iY = columns[1]-1;

	if (a[iX]==b[iX])
	  {
		if (a[iY]==b[iY])
		  {
			return 1; // don't care
		  }
		else return (a[iY] < b[iY]) ? -1 : 1;
	  }
	else return (a[iX] < b[iX]) ? -1 : 1;
  }
}