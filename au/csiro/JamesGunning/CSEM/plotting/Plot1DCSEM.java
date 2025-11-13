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
import java.util.*;
import java.awt.*;
import static java.lang.Math.*;

import cern.colt.map.*;
import cern.colt.list.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.util.*; 
import au.csiro.JamesGunning.lib.graphutil.*;

/**
   A very simple class to show occam response files. 
   Parses occam output files, does some sorting etc,

   @author James.Gunning@csiro.au
   @version 1.0, 30/1/2008
*/

public class Plot1DCSEM 
// implements MouseListener, MouseMotionListener, KeyListener,	ComponentListener, ActionListener, Printable
{
  public static String grcmdDefault = 
  "graph -TX -S 1 .01  -F HersheySans --frame-color wheat "
	+"--bg-color black -f 0.03 -w 0.82 -h 0.85 -r 0.13 -u 0.1 -C --bitmap-size 700x700 ";
  String grcmdExtraOpts="";
  String grTitle=null;

  String grcmd = Plot1DCSEM.grcmdDefault;
  String grout=null;
  ErrorHandler EH;
  String filenames[];
  boolean restrictFreqs=false;
  double[] restrictedFreqs;
  boolean restrictTypes=false;
  int[] restrictedTypes;
  boolean setSymSize = false;
  float symSize;
  boolean joinData = false;

  ArrayList d = new ArrayList();
  boolean noModel=false,noData=false;
  cern.colt.map.OpenDoubleIntHashMap fmap;
  boolean plotEAmps,plotEPhases,plotBAmps,plotBPhases;
  boolean combineGraphs;
  boolean normalVideo=false;

  // I/O
  BufferedReader in=null;

  // Under, here, all the entities for java graphics
  public boolean USE_JAVA_PLOTTING = true;
  int N_Pgroups=1;
  ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;
  CanvasPlotter CP;

  public Plot1DCSEM()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
  }
  
  /** 
	  Show CSEM graphs for named response files, as per Usage()
  */
  public static void main(String argv[]) 
  {
	int i,j;
	String arg;
	double x;
	String[] st;


	Plot1DCSEM PO = new Plot1DCSEM();


	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("csem1Dplotter",PO.EH);

   	options.addMandatoryArg(1,"resp_file(s).txt","response file(s) (space separated), or data files. Format is resp* or *.dat files, "
							+ "both new [EMData_1.1] and legacy Scripps formats");	
	options.addOption("-f",null,1,"String","f1,f2... ","restrict to frequencies named");
	options.addOption("-t",null,1,"String","types t1,t2... ","restrict to data-types (integers) named");
	options.addOption("-g",null,1,"String","graph command","user supplied graph command alternative to: %n "+Plot1DCSEM.grcmdDefault);
	options.addOption("-o",null,1,"String","filename","stdout file to capture graph stdout if required (e.g. plot.ps)");
	options.addOption("-e",null,1,"String","options","extra options to graph (quote protect)");
	options.addOption("-T",null,1,"String","title","graph title");
	options.addOption(null,"--Nm",0,null,null,"Don't plot model");
	options.addOption(null,"--Nd",0,null,null,"Don't plot data");
	options.addOption("-S",null,1,"float","size","datum symbol size (default 0.01)");
	options.addOption("-J",null,0,null,null,"plot data points as if modelled (no error bars & joined)");
	options.addOption(null,"--comb",0,null,null,"Combine phase and amplitude graphs to one window (with -P only)");
	options.addOption("-N",null,0,null,null,"normal video (white background)");
	options.addOption("-P",null,0,null,null,"use plotutils graphics. Supports only legacy Scripps formats");

	options.addSpielText("Optionally uses the GNU plotutils for graphics. The \"graph\" "
					      +"command from this package must be in your path to use this");



	// parse options. Can be many files
	options.parseOptions(argv,false);

	PO.filenames = options.getMandatoryArg(1).split(",");

	if (options.isOptionSet("-N")) PO.normalVideo=true; 

	if (options.isOptionSet("-P")) 
	  {
		PO.USE_JAVA_PLOTTING = false;
	  }
	if (options.isOptionSet("-f")) 
	  {
		PO.restrictFreqs = true;
		arg = (String)options.getOptionValues("-f")[0];
		st = arg.split("[,;:]");
		PO.restrictedFreqs = new double[st.length];
		for (i=0; i< st.length; i++)
		  {
			try
			  {
				x = Double.parseDouble(st[i]);
			  }
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("-f "+arg+" contains a non-number");
			  }
			if (x<=0)
			  {
				throw new RuntimeException("-f "+PO.restrictedFreqs+" contains -ve frequency");
			  }
			PO.restrictedFreqs[i]=x;
		  }
	  }
	if (options.isOptionSet("-t")) 
	  {
		PO.restrictTypes = true;
		arg = (String)options.getOptionValues("-t")[0];
		st = arg.split("[,;:]");
		PO.restrictedTypes = new int[st.length];
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
				throw new RuntimeException("-t "+PO.restrictedTypes+" contains no-positive type");
			  }
			PO.restrictedTypes[i] = j;
		  }
	  }
	if (options.isOptionSet("-g")) 
	  {
		PO.grcmd = (String)options.getOptionValues("-g")[0];
	  }
	if (options.isOptionSet("-o")) 
	  {
		PO.grout = (String)options.getOptionValues("-o")[0];
	  }
	if (options.isOptionSet("-e")) 
	  {
		PO.grcmdExtraOpts = (String)options.getOptionValues("-e")[0];
	  }
	if (options.isOptionSet("-T")) 
	  {
		PO.grTitle = (String)options.getOptionValues("-T")[0];
	  }
	if (options.isOptionSet("--Nm")) PO.noModel = true; 
	if (options.isOptionSet("--Nd")) PO.noData = true; 
	if (options.isOptionSet("-S")) 
	  {
		PO.setSymSize = true;
		PO.symSize = Math.abs(((Float)options.getOptionValues("-S")[0]).floatValue());
	  }
	if (options.isOptionSet("-J")) 
	  {
		PO.joinData = true;
	  }
	if (options.isOptionSet("--comb")) PO.combineGraphs = true; 


	if (PO.USE_JAVA_PLOTTING)
	  {
		PO.readAndJavaPlot();
	  }
	else
	  {
		PO.readAndGraphPlot();
	  }
	
  }

  void readAndGraphPlot()
  {
	readData();

	//	plot
	int i,j,i_set,pass,type;
	boolean allow=true;
	float[] v;
	float r=0,rl=0,yd,ym,e,re;
	String fmod,fdat,name=null;
	int nFiles;
	File[] temp;
	String[] fnames;
	PrintWriter out; 
	boolean firstDataSet,dataOnly;

	setPlotType();

	nFiles = (plotEPhases && plotEAmps) ? 4 : 2;
	temp = new File[nFiles];
	fnames = new String[nFiles];
	// may be up to 4:
	// fnames[0] = model line (amp), plotted on log scale
	// fnames[1] = date+err (amp)
	// fnames[2] = model line (phase), linear scale
	// fnames[3] = data+err (phase)
	// Open tmp files
	try 
	  {
		for(i=0;i<nFiles;i++)
		  {
			temp[i] = File.createTempFile("occam_junk"+i,".txt");
			//if (i==0) fmod = temp[i].getAbsolutePath();
			//if (i==1) fdat = temp[i].getAbsolutePath();
			fnames[i] = temp[i].getAbsolutePath();
			temp[i].deleteOnExit();
		  }
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + name);
	  }
	firstDataSet=true;
	for (pass=0;pass<nFiles;pass++)
	  {
		// name = (pass==0)?fmod:fdat;
		name = fnames[pass];
		try 
		  {
			// out = new PrintWriter(new BufferedWriter(new FileWriter((pass==0)?temp[1]:temp2)));
			out = new PrintWriter(new BufferedWriter(new FileWriter(temp[pass])));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + name);
		  }
		// set colour
		i_set = 1;
		if (pass%2 == 0)
		  {
			if (noModel) continue;
		  }
		else
		  {// pass=1,3
			if (noData) continue;
		  }
		for (i=0; i< d.size(); i++)
		  {
			allow=true;
			v = (float[])d.get(i);
			if (restrictFreqs)
			  {
				for (j=0,allow=false; j<restrictedFreqs.length; j++) 
				  {
					if (Math.abs(v[2] - restrictedFreqs[j]) < 1e-3) 
					  {
						allow=true;break;
					  }
				  }
			  }
			if (!allow) continue;

			type = (int)v[0];

			if (restrictTypes)
			  {
				for (j=0,allow=false; j<restrictedTypes.length;j++) 
				  if (type == restrictedTypes[j]) 
					{
					  allow=true;break;
					}
			  }
			if (!allow) continue;
			// save phase info for pass 2,3
			if (plotEAmps && ((type==2) || (type==4)) && (pass<2)) continue;
			// skip amps for phase plot
			if (plotEPhases && !((type==2) || (type==4)) && (pass>1)) continue;

			r = v[1];
			if (!((type==2) || (type==4)))
			  {
				ym = Math.abs(v[7]);
				yd = Math.abs(v[6]);
			  }
			else
			  {
				ym = v[7];
				yd = v[6];
			  }
			re = v[8];

			// detect data-only entries;
			dataOnly = (v[9] > 0);
			e = dataOnly ? re : ((re!=0.0) ? Math.abs((ym-yd))/re : 0);

			if ((r < rl) || firstDataSet) 
			  {
				// new data set
				if (!firstDataSet) 
				  {
					out.printf("%n");
				  }
				i_set = fmap.get((double)v[2]);
				if (pass%2 == 0)
				  {
					if (type==2 || type==15) 
					  {
						// imag or phase
						out.printf("#m=%d,S=0%n",i_set+20);  // longdashed
					  } 
					else 
					  {
						out.printf("#m=%d,S=0%n",i_set);  
					  }
				  }
				else
				  {
					// pass=1,3
					if (type==2 || type==15) 
					  {
						// imag or phase = circle
						out.printf("#m=%d,S=4%n",-i_set); 
					  } 
					else 
					  {
						// filled square
						out.printf("#m=%d,S=17%n",-i_set);  
					  } 
				  }
				firstDataSet=false;
			  }
			if (pass%2 == 0) 
			  {
				// fitted model
				out.printf("%10.5g %10.5g%n",r,ym);				
			  }
			else
			  {
				// data
				out.printf("%10.5g %10.5g %10.5g%n",r,yd,e);				
			  }
			rl=r;
		  }
		out.close();	
	  }

	for (i=0; i < ((!combineGraphs && plotEAmps && plotEPhases) ? 2 : 1) ; i++)
	  {
		try
		  {
			if (combineGraphs && plotEAmps && plotEPhases)
			  {
				name = grcmd+" --bitmap-size 1300x400 --reposition 0 0 0.5 -f 0.06 -w 0.7 -h 1.6 -u 0.2 " + ((i==0)? " -ly ": "")
				  +grcmdExtraOpts+
				  ((grTitle!=null) ? (" -h 0.82 -L "+ grTitle+" ") : " ") 
				  +fnames[2*i]+" -Ie "+fnames[2*i+1];

				name = name + " --blankout 1 --reposition 0.5 0 0.5 -Ia -ly " +fnames[2]+" -Ie "+fnames[3];
			  }
			else
			  {
				name = grcmd+ ((plotEAmps && i==0)? " -ly ": "")
				  +grcmdExtraOpts+
				  ((grTitle!=null) ? (" -h 0.82 -L "+ grTitle+" ") : " ") 
				  +fnames[2*i]+" -Ie "+fnames[2*i+1];
			  }
			EH.println("Attempting: "+name);
			Runtime runt = Runtime.getRuntime();
			Process p = runt.exec(name); 
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
  }


  void readAndJavaPlot()
  {
	int i,j,code=0,code_last=0,i_set,pass,type=1,N,i_data;
	int NplotsY, NplotsX;
	boolean allow=true,legacyFormat=true,logPlot=true;
	float[] v=null;
	float r=0,rl=0,yd=1,ym=1,e=1,re;
	
	boolean finalize = false,dataOnly = false;
	CanvasPlotter.TwoDPlot pMod = null,pDat = null;
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDPlotGroup _pg=null, _pgE = null,_pgEph=null, _pgB=null, _pgBph=null;
	DoubleArrayList X=null,Y=null,Ym=null,Ye=null;

	readData();
	setPlotType();

	N_Pgroups = 0;
	if (plotEAmps) N_Pgroups++;
	if (plotEPhases) N_Pgroups++;
	if (plotBAmps) N_Pgroups++;
	if (plotBPhases) N_Pgroups++;

	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(N_Pgroups);
	CP = new CanvasPlotter(true);
	if (normalVideo) CP.reverseVideo=false;
	NplotsY = 1;
	NplotsX = 1;
	// top row of plots is E  Ephase
	if ((plotEAmps && plotEPhases) || (plotBAmps && plotBPhases)) NplotsX = 2;
	// next row is B, Bphase (if there is any E)
	if ((plotEAmps || plotEPhases) && (plotBAmps || plotBPhases)) NplotsY = 2;

	if (plotEAmps)
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = NplotsX; 
		pg.GG.plotGroupNY = NplotsY; 
		pg.GG.setTickNumbersEngStyle(false,true);	
		pg.GG.title = grTitle;
		pg.GG.Yword = "|E|,Re,Im(E_i)";
		pg.GG.Xword = "offset";
		_pgE = pg;
	  }
	if (plotEPhases) 
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupXindx = plotEAmps ? 1 : 0;
		pg.GG.plotGroupNX = NplotsX; // amps and phases side-by-side	  
		pg.GG.plotGroupNY = NplotsY; 
		pg.GG.setTickNumbersEngStyle(false,false);
		if (grTitle!=null) pg.GG.title = grTitle+" E phase";
		pg.GG.Yword = "|E|, E_i phase";
		pg.GG.Xword = "offset";
		_pgEph = pg;
	  }
	if (plotBAmps)
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = NplotsX; 	  
		pg.GG.plotGroupNY = NplotsY; 
		pg.GG.plotGroupYindx = NplotsY-1;
		pg.GG.setTickNumbersEngStyle(false,true);	
		pg.GG.title = grTitle;
		pg.GG.Yword = "|B|,Re,Im(B_i)";
		pg.GG.Xword = "offset";
		_pgB = pg;
	  }
	if (plotBPhases) 
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupXindx = plotBAmps ? 1 : 0;
	    pg.GG.plotGroupNX = NplotsX; 	  
		pg.GG.plotGroupNY = NplotsY; 
		pg.GG.plotGroupYindx = NplotsY-1;
		pg.GG.setTickNumbersEngStyle(false,false);
		if (grTitle!=null) pg.GG.title = grTitle+" B phase";
		pg.GG.Yword = "|B|,B_i phase";
		pg.GG.Xword = "offset";
		_pgBph =  pg;
	  }

	
	int nFiles = filenames.length;


	pass=0;
	
	// set colour
	i_set = 1;
	i_data = -1;
	rl = Float.MAX_VALUE;

	N = d.size();
	for (i=0; i<= N ; i++)
	  {
		if (i<N)
		  {
			v = (float[])d.get(i);
			if (!useThisDatum(v)) continue;

			type = (int)v[0];
			code = (int)v[10];
			legacyFormat = (code <= 1); 

			// code =v[10] format: 0,1 = Legacy {E, Ephase}. 2,3,4,5 = NewScripps = {E,Ephase, B, Bphase}
			r = v[1];
			// detect data-only entries;
			dataOnly = (v[9] > 0);

			if (legacyFormat) 
			  {		
				// error
				re = v[8];
		
				if (code==1)
				  { // phase
					ym = v[7];
					yd = v[6];
				  }
				else
				  { // amplitude
					if (!(type==5 || type==6))
					  {
						ym = Math.abs(v[7]);
						yd = Math.abs(v[6]);
					  }
					else
					  {
						// type = 5, 6 is log10(Erho, Ephi). Convert to absolute amplitude
						ym = (float)pow(10,v[7]);
						yd = (float)pow(10,v[6]);
						// rescale error
						if (dataOnly) 
						  {
							re = yd * (float)(v[8] * MathUtil.LOG10);						
						  }
						else
						  {
							re = (float)(v[8]);				
						  }
					  }
				  }
			  }
			else
			  {
				if ((code == 3) || (code == 5))  // E or B phase
				  {
					ym = v[7];
					yd = v[6];
				  }
				else
				  {
					ym = Math.abs(v[7]);
					yd = Math.abs(v[6]);
				  }
				re = v[8];
			  }
			
			e = dataOnly ? re : ((re!=0.0) ? Math.abs((ym-yd))/re : 0);
		  }
		else 
		  {
			finalize=true; // last dummy run to populate graph entities
		  }
		if ((r < rl) || finalize) 
		  {
			if (!finalize) i_data++;
			// use code as per last data set: code_last
			// code format: 0,1 = Legacy {E, Ephase}. 2,3,4,5 = NewScripps = {E,Ephase, B, Bphase}
			if (legacyFormat)
			  {
				if (code_last == 0) 
				  {
					_pg = _pgE ;
					logPlot=true;
				  }
				else
				  {
					_pg = _pgEph;
					logPlot=false;
				  }
			  }
			else
			  {
				switch (code_last)
				  {
				  case 0: // only occurs at initialisation (means nothing)
				  case 2:
					_pg = _pgE; logPlot=true; break;
				  case 3:
					_pg = _pgEph; logPlot=false; break;
				  case 4:
					_pg = _pgB; logPlot=true; break;
				  case 5:
					_pg = _pgBph; logPlot=false; break;
				  default:
					EH.die_err("Illegal code! Bug! Fixme!");
				  }
			  }
			if (pMod != null && _pg != null)
			  {
				// finish up last set 
				X.trimToSize();
				Y.trimToSize();
				Ym.trimToSize();
				Ye.trimToSize();				
				pMod.x = X.elements(); 
				pMod.y = Ym.elements();
				pMod.e = null;
				pMod.drawline = true;
				// pMod.c = Color.black;
				
				pDat.x = pMod.x;
				pDat.y = Y.elements();
				pDat.e = Ye.elements();
				// pDat.c = Color.blue;
				if (!joinData)
				  {
					pDat.symbol = 'r';
					pDat.drawline = false;
				  }
				else
				  {
					pDat.drawline = true;	
					pDat.e = null;
				  }
				if (setSymSize) pDat.symSize = symSize;

				_pg.GG.logPlotY = logPlot;
				if (!noData) 
				  {
					_pg.plots.add(pDat);
					_pg.GG.computeRanges(pDat.x,pDat.y);
				  }
				if (!noModel && !dataOnly)
				  {
					_pg.plots.add(pMod);
					_pg.GG.computeRanges(pMod.x,pMod.y);
				  }
			  }
			if (!finalize)
			  {
				pMod = new CanvasPlotter.TwoDPlot();
				pDat = new CanvasPlotter.TwoDPlot();
				X = new DoubleArrayList();
				Ym = new DoubleArrayList();
				Y = new DoubleArrayList();
				Ye = new DoubleArrayList();
			  }
		  }

		if (!finalize)
		  {
			i_set = fmap.get((double)v[2]) - 1;  // hashmapped colors are 1,2,3,.....Make zero-offset

			Color c = CanvasPlotter.defaultColors[i_set % CanvasPlotter.defaultColors.length];
			pMod.c = c;
			pDat.c = c;
			X.add(r); 
			Ym.add(ym);  
			Y.add(yd); 
			Ye.add(e);
			rl = r;
			code_last = code;
		  }
	  }


	CP.windowTitle = "csem1DPlotter";
	// CP.run(Pgroups);
	if (grout == null) 
	  { 
		CP.run(Pgroups);
	  }
	else
	  {
		CP.runDisplayAsEPS(Pgroups,grout);
		System.exit(0);
	  }
  }

  private boolean useThisDatum(float v[])
  {
	int j,type;
	boolean	allow=true;

	if (restrictFreqs)
	  {
		for (j=0,allow=false; j<restrictedFreqs.length; j++) 
		  {
			if (Math.abs(v[2] - restrictedFreqs[j]) < 1e-3) 
			  {
				allow=true;break;
			  }
		  }
	  }
	if (!allow) return false;
	
	type = (int)v[0];
	
	if (restrictTypes)
	  {
		for (j=0,allow=false; j<restrictedTypes.length;j++) 
		  if (type == restrictedTypes[j]) 
			{
			  allow=true;break;
			}
	  }
	if (!allow) return false;

	return true;
  }


  void setPlotType()
  {
	int i,j,N,type;
	boolean isEPhase,isBPhase,isEAmp,isBAmp,allow;
	boolean legacyFormat;
	float[] v;

	plotEPhases = false;
	plotEAmps = false;
	plotBPhases = false;
	plotBAmps = false;

	isEPhase = isBPhase = isEAmp = isBAmp = false;
	// equiv to pass==1 call
	N = d.size();
	for (i=0; i< N; i++)
	  {
		allow = true;
		v = (float[])d.get(i);
		if (restrictTypes)
		  {
			for (j=0,allow=false; j<restrictedTypes.length;j++) 
			  {
				if (v[0] == restrictedTypes[j]) 
				  {
					allow=true;break;
				  }
			  }
		  }
		if (!allow) continue;
		type = (int)v[0];
		// Coming in to this routine: Format: v[10] = : 0 = Legacy, 1 = NewScripps
		legacyFormat = (v[10] == 0); 
		if (legacyFormat) 
		  {
			// legacy scripps
			isEPhase = ((type==2) || (type==4));
			isEAmp = ((type==1) || (type==3)           // Erho, Ephi
					  || ((type>=14) && (type<=17)));  // Re,Im{Erho,Ephi}
			isEAmp |= (type==5||type==6);   // log10(Erho), log10(Ephi)
			v[10] = isEPhase ? 1 : 0;
		  }
		else
		  {
			// new Scripps
			isEAmp = (((type==21) || (type==23) || (type==25))           // |Ex,Ey,Ez|
					  || ((type>=1) && (type<=6)));                      // Re,Im{Ex,Ey,Ez}
			isEPhase = ((type==22) || (type==24) || (type==26));
			isBAmp= ((type>=11 && type<=16)                              // Re,Im { Bx,By,Bz}
					 || ((type==31) || (type==33) || (type==35)));       // |{ Bx,By,Bz}|
			isBPhase = ((type==32) || (type==34) || (type==36));

			if (isEAmp)        v[10] = 2;
			else if (isEPhase) v[10] = 3;
			else if (isBAmp  ) v[10] = 4;
			else if (isBPhase) v[10] = 5;
		  }

		// After we exit:
		// v[10] = 0    legacy, amp or component 
		// v[10] = 1    legacy, phase
		// v[10] = 2    new, E amp or component 
		// v[10] = 3    new, E phase
		// v[10] = 4    new, B amp or component 
		// v[10] = 5    new, B phase


		if (isEPhase) plotEPhases=true;
		if (isBPhase) plotBPhases=true;
		if (isEAmp) plotEAmps=true;
		if (isBAmp) plotBAmps=true;
		if (!plotEPhases && !plotBPhases && !plotEAmps && !plotBAmps)
		  {
			EH.die_err("No data to show. Check type code in data file, or -t, -f exclusion flags...");
		  }
	  }
  }

  void readData()
  {
	int i,j,k,col,lineNo;
	String s,word,ss,filename;
	String sCode,sValue;
	String[] sa=null;
	String[] fn;
	StringTokenizer st;
	float[] v=null,vtmp=null;
	ArrayList fd;
	DataComparator cmp = new DataComparator();
	boolean isDataFile,isNewScripps=false;

	// New scripps..
	int iTx=0,iRx=0,iFreq=0;
	LineSplitter LS = null, LS1=null;
	EMDataNewScrippsInfo info =null;

	fmap = new OpenDoubleIntHashMap();
	col=1;

	// loop over files
	fn = filenames;
	for (i=0; i< fn.length; i++)
	  {
		filename = fn[i];
		try
		  {
			in = new BufferedReader(new FileReader(filename));
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't open file "+filename);
		  }
		fd = new ArrayList();
		lineNo=0;
		isDataFile = false;
		try
		  {
			while ((s=in.readLine()) != null)
			  {
				// try to cope also with data-files
				if (lineNo==0)
				  {
					ss = s.trim().toLowerCase();
					if (ss.startsWith("format:")) 
					  {
						isDataFile = true; 
						if (ss.contains("emdata_1")) {isNewScripps = true; isDataFile = true;}
						else if (ss.contains("emresp_1")) {isNewScripps = true; isDataFile = false;}
						if (isNewScripps)
						  {
							EH.warn("Expecting new Scripps: "+s);
							LS = new LineSplitter();
							LS1 = new LineSplitter();
						  }
						else
						  {
							EH.warn("Expecting legacy Scripps: "+s);
							vtmp = new float[8];
						  }
					  }
				  }
				if (!isNewScripps)
				  {					
					if (isDataFile)
					  {
						// after 4 line header, data files have 8 cols
						lineNo++;
						if (lineNo > 4)
						  {	
							v = new float[11];
							v[9] = 1; // flags data
							// v[10] = 0; // flags legacy scripps
							j=0;
							st = new StringTokenizer(s);
							while (st.hasMoreTokens()) 
							  {
								if (j>=8) throw new RuntimeException("Too many numbers, line "+lineNo+" of "+filename);
								word = st.nextToken();
								try
								  {
									vtmp[j] = Float.parseFloat(word);
									j++;
								  }
								catch (NumberFormatException e)
								  {
									EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers: "+word);
								  }				
							  }
							mapDataToResp(vtmp,v);
							// collect frequency hashmap while we're here
							if (!fmap.containsKey((double)v[2]))
							  { 
								fmap.put((double)v[2],col);
								col++;
							  }
							fd.add(v);
						  }
					  }
					else
					  {
						// Legacy response file
						v = new float[11]; 
						// v[9]=0; // flags response
						// v[10] = 0; // flags legacy scripps
						j=0;
						st = new StringTokenizer(s);
						lineNo++;
						while (st.hasMoreTokens()) 
						  {
							if (j>=9) throw new RuntimeException("Too many numbers, line "+lineNo+" of "+filename);
							word = st.nextToken();
							try
							  {
								v[j] = Float.parseFloat(word);
								j++;
							  }
							catch (NumberFormatException e)
							  {
								EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
							  }				
						  }
						// collect frequency hashmap while we're here
						if (!fmap.containsKey((double)v[2]))
						  { 
							fmap.put((double)v[2],col);
							col++;
						  }
						fd.add(v);
					  }
				  } // old scripps
				else
				  {
					// New scripps
					LS.parse(s);
					if (LS.bComment) continue;
					sCode = LS.sCode;
					sValue = LS.sValue;
					if (sCode.equals("# transmitters"))
					  {
						info = new EMDataNewScrippsInfo();
						try
						  {
							info.nTxIn = Integer.parseInt(sValue);
							// current restriction:
							if (info.nTxIn != 1)
							  {
								//EH.die_err("Can't plot data with >1 transmitter currently. Sorry.");
							  }
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Number of transmitters not a positive integer");
						  }
						info.TxInfo = new double[info.nTxIn][5];
						//
						// Now read in block of transmitters, skipping any comment lines:
						//
						for (j=0; j < info.nTxIn;) 
						  { 
							// Get the next code/value pair
							try
							  {
								if ((s = in.readLine()) == null) break;
							  }
							catch (IOException e)
							  {
								EH.die_err("Insufficient transmitters");
							  }
        
							// Check if line is a comment, if not strip off comments at end of line:
							LS.parse(s);
							if (LS.bComment ) continue;
							// Else it is not a comment, so read in the Tx parameters from s
							sa = s.trim().split("\\s+");            
							if (sa.length >= 5)
							  {
								try
								  {
									// k= 0,1,2 = x,y,z; 3,4=Azimuth Dip
									for (k=0;k<5;k++) info.TxInfo[j][k] = Util.parseDouble(sa[k]);
								  }
								catch (NumberFormatException e)
								  {
									EH.die_err("X Y Z Azimuth Dip contains a non-number: "+s);
								  }
							  }
							else
							  {
								EH.die_err("Not 5 number entries on transmitter line "+s);
							  }
							j++;   
						  }            
					  }
					else if (Util.matches(sCode,"# frequencies","nfreq","#freq","#freqs", "# freq"))
					  {
						// Get the next code/value pair
						try
						  {
							info.nFreqIn = Integer.parseInt(sValue);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Indecipherable # frequencies - not an integer");
						  }
						info.fTxIn = new double[info.nFreqIn];

						//
						// Now read in block of frequencies, skipping any comment lines:
						//
						for (j = 0; j < info.nFreqIn;) 
						  {
							// Get the next code/value pair
							try
							  {
								if ((s = in.readLine()) == null) break;
							  }
							catch (IOException e)
							  {
								EH.die_err("Insufficient frequencies");
							  }
							// Check if line is a comment, if not strip off comments at end of line:
							LS1.parse(s);
							if (LS1.bComment ) continue;

							// Else it is not a comment, so read in the Tx parameters from s
							try  
							  {
								info.fTxIn[j] = Util.parseDouble(s);
							  }
							catch (NumberFormatException e)
							  {
								EH.die_err("Frequency not a number");
							  }
							j++;
						  }        
						//
						// Check and make sure all the frequencies were read in:
						//
						if (j != info.nFreqIn) 
						  {
							EH.die_err(" Error reading frequencies: iFreqRead!=nFreqIn: "+j+" "+info.nFreqIn);
						  }
					  }
					else if (Util.matches(sCode,"# receivers","nrx","#receivers"))
					  {
						try  
						  {
							info.nRxIn = Integer.parseInt(sValue);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("# receivers not an integer");
						  }

						// Allocate the field arrays: 
						//
						info.RxInfo = new double[info.nRxIn][3];
						// Now read in block of receiver locations, skipping any comment lines:
						//
						for (j = 0;j < info.nRxIn;) 
						  {
							// Get the next code/value pair
							try
							  {
								if ((s = in.readLine()) == null) break;
							  }
							catch (IOException e)
							  {
								EH.die_err("Insufficient receivers");
							  }
        
							// Check if line is a comment, if not strip off comments at end of line:
							LS1.parse(s);
							if (LS1.bComment ) continue;
							// Else it is not a comment, so read in the receivers from s
							sa = s.trim().split("\\s+");            
							if (sa.length >= 3)
							  {
								try
								  {
									for (k=0;k<3;k++) info.RxInfo[j][k] = Util.parseDouble(sa[k]);
								  }
								catch (NumberFormatException e)
								  {
									EH.die_err(s+": Not numbers for receivers");
								  }
							  }
							else
							  {
								EH.die_err("Not 3 columns of numbers for receivers");
							  }
							j++;   
						  }
						//
						// Check and make sure all the layers were read in:
						//
						if (j!=info.nRxIn)
						  {
							EH.die_err("Error reading receivers: j != nRxIn: "+j+" "+info.nRxIn);
						  }
					  }
					else if (Util.matches(sCode,"# data","#data","ndata"))
					  {
						//
						// Lastly we need to read in the DATA block.
						//
						try  
						  {
							info.nd = Integer.parseInt(sValue);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Number of data not an integer");	
						  }
						// loop over data and add to collections
						for (j = 0; j < info.nd;) 
						  {
							// Get the next code/value pair
							try
							  {
								if ((s = in.readLine()) == null) break;
							  }
							catch (IOException e)
							  {
								EH.die_err("Insufficient data");
							  }        
							// Check if line is a comment, if not strip off comments at end of line:
							LS1.parse(s);
							if (LS1.bComment ) continue;
							// Else it is not a comment, so read in the Tx parameters from s
							sa = s.trim().split("\\s+");            
							if (sa.length >= 6)
							  {								
								v = new float[11];
								if (isDataFile) v[9] = 1; // flags data
								v[10] = 1; // flags new Scripps
								try
								  {
									v[0] = Integer.parseInt(sa[0]); // type
									iFreq = Integer.parseInt(sa[1]);
									iTx = Integer.parseInt(sa[2]);  // transmitter #
									iRx = Integer.parseInt(sa[3]);  // receiver #
									v[6] = Util.parseFloat(sa[4]);  // data
									v[7] = v[6];                    // response (default)
									v[8] = Util.parseFloat(sa[5]);  // error
								  
									if (!isDataFile && sa.length >= 8)
									  {
										// Resp file: get extra fields
										v[9] = 0; // flags response
										v[7] = Util.parseFloat(sa[6]);  // data
										v[8] = Util.parseFloat(sa[7]);  // err
									  }
									  
									// restore geometries etc
									// range -> use euclidean TX-Rx distance (got a better idea?)
									if (iRx<1 || iRx>info.nRxIn)
									  {
										EH.die_err("Illegal receiver number in data section, data-line: "+s);
									  }
									if (iTx<1 || iTx>info.nTxIn)
									  {
										EH.die_err("Illegal transmitter number in data section, data-line: "+s);
									  }
									v[1] = (float)Math.hypot( info.RxInfo[iRx-1][0] - info.TxInfo[iTx-1][0], 
															  info.RxInfo[iRx-1][1] - info.TxInfo[iTx-1][1]);
									// freq
									v[2] = (float)info.fTxIn[iFreq-1]; 
								  }
								catch (NumberFormatException e)
								  {
									EH.die_err("Not numbers in data section, data-line "+s);
								  }
							  }
							else
							  {
								EH.die_err("Not 6 columns of numbers for data");
							  }
							// add to collection/array
							if (!fmap.containsKey((double)v[2]))
							  { 
								fmap.put((double)v[2],col);
								col++;
							  }
							fd.add(v);
							j++;
						  }
			
						// Check and make sure all the data was read in:
						if (j != info.nd)
						  {
							EH.die_err("Error reading receivers: j!=nd: "+j+ " "+info.nd);
						  }
					  }				
				  }
			  }
			in.close();
			// sort this file's data by ascending range...
			Collections.sort(fd,cmp);
			// append to master collection
			d.addAll(fd);
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't close "+filename);
		  }

	  }
  }

  // For handling legacy formats: stuff cols from data file to equiv place in resp file
  void mapDataToResp(float[] vtmp,float[] v)
  {
	int i,j;

	v[0] = vtmp[0];                   // type
	for(i=1,j=3;i<=5;i++,j++) v[i] = vtmp[j]; // range
	v[6] = vtmp[1]; // modelled
	v[7] = vtmp[1]; // data
	v[8] = vtmp[2]; // simply save actual here
  }

  // Utility class for sorting collected data into chunks by
  // frequency, type, and offset
  class DataComparator implements Comparator<float[]> 
  {
	// sort by columns 3, 1, 2 = { freq, type, range}
	public int compare(float[] c1, float[] c2) 
	{ 
	  if (c1[2] != c2[2])
		{
		  return (c1[2] < c2[2]) ? -1 : 1;
		}
	  else if (c1[0] != c2[0])
		{
		  return (c1[0] < c2[0]) ? -1 : 1;
		}
	  else if (c1[1] != c2[1])
		{
		  return (c1[1] < c2[1]) ? -1 : 1;
		}
	  else
		{
		  // don't care
		  return 1;
		}
	}
  }

  // Small utility class for handling new Scripps formats
  class EMDataNewScrippsInfo
  {
	int nTxIn;
	double [][] TxInfo;
	int nFreqIn;
	double[] fTxIn;
	int nRxIn;
	double [][] RxInfo;
	int nd; // n data
  }
}