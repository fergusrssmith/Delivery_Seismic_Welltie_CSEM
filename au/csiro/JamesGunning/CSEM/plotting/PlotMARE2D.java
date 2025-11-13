/**********************************************************************
Copyright (C) 2010, CSIRO,  all rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of CSIRO nor the names of its
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

package au.csiro.JamesGunning.CSEM.plotting;

import java.io.*;
import java.util.*;
import java.awt.*;

import java.awt.event.*; 

import static java.lang.Math.*;

import cern.colt.map.OpenDoubleIntHashMap;


import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.util.*; 
import au.csiro.JamesGunning.CSEM.TwoDinversion.*;
import au.csiro.JamesGunning.lib.graphutil.*;

/**
   Plot MARE2D fields in panel
   @author James.Gunning@csiro.au
   @version 1.0, 30/4/2010
*/

public class PlotMARE2D 
// implements MouseListener, MouseMotionListener, KeyListener,	ComponentListener, ActionListener, Printable
{
  ErrorHandler EH;
  String filename;
  float minClip=-Float.MAX_VALUE, maxClip=Float.MAX_VALUE;
  boolean plotPhase=false;
  boolean restrictTrans = false;
  int transNumber = 0;
  boolean addA2Dfile=false;
  String A2Dfilename=null;
  String linecolors=null;
  boolean plot2DfileOnly=false;
  int freqNumber=0;
  boolean restrictFreq=false;
  boolean normalVideo=false;

  // Under, here, all the entities for java graphics
  int N_Pgroups=1;
  ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;
  CanvasPlotter CP;
  CanvasPlotter.TwoDPlotGroup pg = null;
  CanvasPlotter.TwoDPlot plot = null;
  

  MARE2DFields MAREfields = null;
  ArrayList<A2D> A2Dlist=null;


  public PlotMARE2D() 
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
  } 
  
  /** 
	  Show triangle as canvas
  */
  public static void main(String argv[]) 
  {
	PlotMARE2D PM = new PlotMARE2D();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("PlotMARE2D",PM.EH);

   	options.addMandatoryArg(1,"MARE2D.csem", "MARE2D response file");
	
	options.addOption(null,"--clip",2,"float","lower upper","clipping range");
	options.addOption(null,"-l",0,null,null,"use log10 scale for color");
	options.addOption(null,"-p",0,null,null,"plot phase");
	options.addOption(null,"-t",1,"int","number","plot numbered transmitter");
	options.addOption(null,"-F",1,"int","number","plot numbered frequency");
	options.addOption(null,"--A2D",1,"String",".A2D_file","plot approximate modelled fields");
	options.addOption(null,"--A2Donly",0,null,null,"plot .A2D file only");
    options.addOption(null,"--cols",1,"String","col1,col2,..","override default linecolors. Supply eg red,blue,black,..(named in java.awt.Color)");
	options.addOption("-N",null,0,null,null,"normal video (white background)");

	// parse options. Can be many files
	if (!options.parseOptions(argv,false))
	  {
		options.Usage();
	  }

	PM.filename = options.getMandatoryArg(1);
	if (options.isOptionSet("-N")) PM.normalVideo=true; 
	if (options.isOptionSet("-p")) PM.plotPhase=true;
	if (options.isOptionSet("-t")) 
	  {
		PM.restrictTrans = true;
		PM.transNumber = options.getOptionInt("-t",0);
	  }
	if (options.isOptionSet("-F")) 
	  {
		PM.restrictFreq = true;
		PM.freqNumber = options.getOptionInt("-F",0);
	  }
	if (options.isOptionSet("--A2D")) 
	  {
		PM.addA2Dfile = true;
		PM.A2Dfilename  = options.getOptionString("--A2D",0);
		if (options.isOptionSet("--A2Donly")) PM.plot2DfileOnly=true;
	  }
    if (options.isOptionSet("--cols")) 
	  {
		PM.linecolors  = options.getOptionString("--cols",0);
	  }

	PM.run();	
  }


  void run()
  {
	MAREfields = new MARE2DFields();
	if (MAREfields.readMARE2Dfile(filename,1))
	  {
		plot();
	  }
  }

  void plot()
  {
	int i,j,j0,jj,c=0,k,iF,N,NperF,Nsub,isAir;
	float freq;
	double[] f,y;
	OpenDoubleIntHashMap map = new OpenDoubleIntHashMap();

	// start & fire off plotter

	CP = new CanvasPlotter(true);
	if (normalVideo) CP.reverseVideo=false;
	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
	pg = new CanvasPlotter.TwoDPlotGroup(CP);
	Pgroups.add(pg);
	pg.GG.logPlotY = plotPhase? false : true;
    ArrayList<Color> ca=null;
    if(linecolors!=null)
      {
        ca = new ArrayList<Color>();
        String[] st = linecolors.split(",");
        for (i=0; i<st.length; i++)
          {
            Color color;
            try
              {
                java.lang.reflect.Field field = Class.forName("java.awt.Color").getField(st[i]);
                color = (Color)field.get(null);
              }
            catch (Exception e)
              {
                EH.println("No such color "+st[i]+": using defaults");
                color = CanvasPlotter.defaultColors[i%CanvasPlotter.defaultColors.length];
              }
            ca.add(color);
          }
      }
	if (!plot2DfileOnly)
	  {
		for(i=0;i<MAREfields.Nf;i++)
		  {
			map.put(MAREfields.frequencies[i],i+1);
			for (j=0;j<MAREfields.Nt;j++)
			  {
				if (restrictTrans && (j+1)!=transNumber) continue;
				if (restrictFreq && (i+1)!=freqNumber) continue;

				// make a line plot
				plot = new CanvasPlotter.TwoDPlot();
				f = new double[MAREfields.Nr];
				for (k=0;k<MAREfields.Nr;k++) 
				  {
					if (plotPhase)
					  {
						f[k] = toDegrees(atan2(MAREfields.Ey_im[i][j][k],MAREfields.Ey_re[i][j][k]));					
					  }
					else
					  {
						f[k] = hypot(MAREfields.Ey_re[i][j][k],MAREfields.Ey_im[i][j][k]);
					  }
				  }
				plot.x = MAREfields.Yr;
				plot.y = f;
				// colour by frequency 
				if(ca!=null && i<ca.size())
                  plot.c = ca.get(i);
                else
                  plot.c = CanvasPlotter.defaultColors[(i+1) % CanvasPlotter.defaultColors.length];                
				pg.plots.add(plot);  
				// get everything in range  
				pg.GG.computeRanges(plot.x,plot.y);
			  }
		  }
	  }

	pg.GG.Xword = "y";
	pg.GG.Yword = plotPhase? "phase": "|Ey|";
	CP.windowTitle = filename;

	if (addA2Dfile)
	  {
		A2D a2d;
		readA2Dfile(A2Dfilename);
		N = A2Dlist.size();
		int nFreq = MAREfields.Nf;
		for (iF=0,k=0;k<N;k++) 
		  {
			a2d = A2Dlist.get(k);
			freq = a2d.f;
			for (j=k+1; j<N;j++) if (A2Dlist.get(j).f != freq) break;
			j--;
			j0=k;
			// frequency group now [j0....j]
			NperF = (j-j0+1);


			// make line plots: clean data, then air data...
			for (isAir=0; isAir<=1; isAir++)
			  {
				if (isAir==0) 
				  {
					for(jj=j0,c=0; jj<=j; jj++) {if (A2Dlist.get(jj).isAir == 0) c++;}
				  }
				// c = # "clean" data, this freq
				Nsub = (isAir==0) ? c : (NperF-c); 
				if (Nsub>0)
				  {
					plot = new CanvasPlotter.TwoDPlot();
					f = new double[Nsub];
					y = new double[Nsub];
					for (i=0,k=j0; k<=j; k++)
					  {
					   	a2d = A2Dlist.get(k);
						if (a2d.isAir != isAir) continue;
						y[i] = a2d.Yr;
						if (plotPhase)
						  {
							f[i] = a2d.phase;
						  }
						else
						  {
							f[i] = a2d.Ey;
						  }
						i++;
					  }
					plot.x = y;
					plot.y = f;
					// colour by frequency
					iF = map.get(freq);
					if (iF==0) map.put(freq,++nFreq);
					plot.c = CanvasPlotter.defaultColors[iF% CanvasPlotter.defaultColors.length];

					plot.symbol = (isAir==0)? 'r' : 't';
					plot.drawline = (isAir==0) ? true : false;
					if (!restrictFreq || (iF==freqNumber))
					  {
						pg.plots.add(plot);
						// get everything in range
						pg.GG.computeRanges(plot.x,plot.y);
					  }
				  }
			  }
 			k--;
		  }
	  }
	// pg.GG.setXRange(pg.GG.xMax,0);
	CP.run(Pgroups); 
	
  }

  /**
	 Read & superpose output file "Approx2Dresponse.A2D" from Approx2DModeller class
  */
  boolean readA2Dfile(String filename)
  {
	int j;
	String s;
	BufferedReader in=null;
	String[] sa;

	try
	  {
		filename = filename;
		in = new BufferedReader(new FileReader(filename));
	  }
	catch (Exception e) 
	  {
		EH.warn(("Can't open file "+filename)); return false;
	  }
	try
	  {
		j=0; // line num
		A2Dlist = new ArrayList<A2D>();
		while ((s=in.readLine()) != null)
		  {
			j++;
			s=s.trim();

			if (s.startsWith("#")) continue;
			sa = s.split("\\s+");     // split just at white space
			A2D a2d = new A2D();
			try
			  {
				a2d.Yr = Float.parseFloat(sa[0]);
				a2d.Yt = Float.parseFloat(sa[1]);
				a2d.f = Float.parseFloat(sa[2]);
				a2d.isAir = Integer.parseInt(sa[3]);
				a2d.Ey= Float.parseFloat(sa[4]);
				a2d.phase = Float.parseFloat(sa[5]);
			  }
			catch (NumberFormatException e)
			  {
				EH.warn(("Line "+j+" of file " + filename+" has non-numbers")); return false;
			  }	
			A2Dlist.add(a2d);
		  }
		in.close();
		EH.warn("read "+filename);
	  }
	catch (Exception e) 
	  {
		EH.warn(("Can't open, decipher or close file "+filename)); return false;
	  }
	A2Dlist.trimToSize();
	return true;
  }

  public static class A2D
  {
	int isAir;
	float Yt,Yr,f,Ey,phase;
  }
} 
