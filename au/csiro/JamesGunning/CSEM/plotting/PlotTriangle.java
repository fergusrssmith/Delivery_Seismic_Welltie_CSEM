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


import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.cgeom.*;
import au.csiro.JamesGunning.lib.util.*; 
import au.csiro.JamesGunning.lib.graphutil.*;

/**
   Plot triangulation grids from triangle output, as for MARE2D
   @author James.Gunning@csiro.au
   @version 1.0, 30/4/2010
*/

public class PlotTriangle implements KeyListener
// implements MouseListener, MouseMotionListener, KeyListener,	ComponentListener, ActionListener, Printable
{
  ErrorHandler EH;
  int verbosity=1;
  public String filename;
  public float minClip=-Float.MAX_VALUE, maxClip=Float.MAX_VALUE;
  public boolean logScale=false;
  public boolean normalVideo=false;
  public boolean engStyleNumbers=true;

  // I/O
  BufferedReader in=null;

  boolean doXrange=false, doYrange=false;
  float Xmin,Xmax,Ymin,Ymax;

 // Under, here, all the entities for java graphics
  int N_Pgroups=1;
  ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;
  CanvasPlotter CP;
  CanvasPlotter.TwoDPlotGroup pg = null;
  CanvasPlotter.TwoDShapePlot plot = null;
  


  TriangleUtil TU=null;
  /*
  int Nt,Na;
  int[][] nodes=null;
  double[][] xys=null;
  double values[]=null;
  */

  public PlotTriangle()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
  }
  
  /** 
	  Show triangle as canvas
  */
  public static void main(String argv[]) 
  {
	PlotTriangle PT = new PlotTriangle();


	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("PlotTriangle",PT.EH);

	options.addSpielText("Toggle triangles wih 't' key");
   	options.addMandatoryArg(1,"poly_file.[poly]","triangle fill to plot with matching .ele and .node files");	
	
	options.addOption(null,"--clip",2,"float","lower upper","clipping range");
	options.addOption(null,"-x",2,"float","Xmin Xmax","X plotting range");
	options.addOption(null,"-y",2,"float","Ymin Ymax","Y plotting range");
	options.addOption(null,"-l",0,null,null,"use log10 scale for color");
	options.addOption("-N",null,0,null,null,"normal video (white background)");
	options.addOption("--NE",null,0,null,null,"Not-Engineering style numbers (xxx.yy)");
	// parse options. Can be many files
	if (!options.parseOptions(argv,false))
	  {
		options.Usage();
	  }

	PT.filename = options.getMandatoryArg(1);
	if (options.isOptionSet("--clip"))
	  {
		PT.minClip = ((Float)options.getOptionValues("--clip")[0]).floatValue();
		PT.maxClip  = ((Float)options.getOptionValues("--clip")[1]).floatValue();
	  }
	if (options.isOptionSet("-x"))
	  {
		PT.doXrange = true;
		PT.Xmin = options.getOptionFloat("-x",0);
		PT.Xmax = options.getOptionFloat("-x",1);
	  }
	if (options.isOptionSet("-y"))
	  {
		PT.doYrange = true;
		PT.Ymin = options.getOptionFloat("-y",0);
		PT.Ymax = options.getOptionFloat("-y",1);
	  }
	if (options.isOptionSet("-l")) PT.logScale=true;
	if (options.isOptionSet("-N")) PT.normalVideo = true;
	if (options.isOptionSet("--NE")) PT.engStyleNumbers = false;

	PT.run();	
  }


  public void run()
  {
	TU = new TriangleUtil();

	if (TU.readFromEleNodeFilenamePair(filename,verbosity,EH))
	  {
		plot(filename);
	  }
  }

  public void setTU(TriangleUtil TU)
  {
	this.TU=TU;
  }

  public void plot(String windowTitle)
  {
	int i;
	int[] npi=null;
	double[][] _xys;
	double[] _values;

	// Set up arrays needed by triangle plotting library
	double[] X = new double[TU.Nt];
	double[] Y = new double[TU.Nt];
	double[] X1 = new double[TU.Nt];
	double[] Y1 = new double[TU.Nt];
	double[] X2 = new double[TU.Nt];
	double[] Y2 = new double[TU.Nt];
	int[] gcol = new int[TU.Nt];

	_xys = TU.xys;
	for(i=0;i<TU.Nt;i++)
	  {
		npi = TU.nodes[i];
		// triangle corners are {X[i],Y[i]}, {X1[i],Y1[i]}  {X2[i],Y2[i]} and colour is TU.values[i]
		X[i] = _xys[npi[0]-1][0];
		Y[i] = _xys[npi[0]-1][1];
		X1[i] = _xys[npi[1]-1][0];
		Y1[i] = _xys[npi[1]-1][1];
		X2[i] = _xys[npi[2]-1][0];
		Y2[i] = _xys[npi[2]-1][1];
	  }

	Arrays.fill(gcol,-1);


	_values = TU.values;
	// threshold colors
	if (logScale)
	  {
		double vmin = Double.MAX_VALUE,vmax = -Double.MAX_VALUE,vU,vL,eps; 
		for (i=0;i<TU.Nt;i++) 
		  {
			vmin = min(vmin,_values[i]);
			vmax = max(vmax,_values[i]);
		  }
		EH.warn("Data value-range: %g %g%n",vmin,vmax);
		if (vmin>0)
		  {
			for (i=0;i<TU.Nt;i++) _values[i] = max(min(log10(_values[i]),maxClip),minClip);
			vmin = max(log10(vmin),minClip);
			vmax = min(log10(vmax),maxClip);
			eps = 0.1;
			vL = vmin-eps*(vmax-vmin);
			vU = vmax+eps*(vmax-vmin);
			for (i=0;i<TU.Nt;i++) _values[i] = (_values[i]-vL)/(vU-vL);
			EH.warn("Log-clipping: %g %g%n",vL,vU);
		  }
	  }
	else
	  {
		for (i=0;i<TU.Nt;i++) _values[i] = max(min(_values[i],maxClip),minClip);
	  }

	// start & fire off plotter
	CP = new CanvasPlotter(true);
	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
	pg = new CanvasPlotter.TwoDPlotGroup(CP);
	Pgroups.add(pg);

	// make a triangle shape plot with lines
	plot = new CanvasPlotter.TwoDShapePlot(true,true);
	// set its elements
	plot.x = X;
	plot.y = Y;
	plot.x1 = X1;
	plot.y1 = Y1;
	plot.x2 = X2;
	plot.y2 = Y2;
	plot.col = _values;
	plot.gcol = gcol;
	pg.GG.reverseY=true; // flip vertical axis

	// get everything in range
	pg.GG.computeRanges(plot.x,plot.y);
	pg.GG.computeRanges(plot.x1,plot.y1);
	pg.GG.computeRanges(plot.x2,plot.y2);

	// These should act via clip methods in the renderer
	if (doXrange) pg.GG.setXRange(Xmax,Xmin);
	if (doYrange) pg.GG.setYRange(Ymax,Ymin);

	plot.normalize(); // rescales colors to [0-1]
	plot.setColorScheme(Color.black,Color.white);  // large numbers go white...

	pg.plots.add(plot);
	pg.GG.Xword = "Y";
	pg.GG.Yword = "Z";
	pg.GG.setTickNumbersEngStyle(engStyleNumbers,engStyleNumbers);
	CP.windowTitle = windowTitle;// filename;

	if (normalVideo) CP.reverseVideo=false;
	CP.addKeyListener(this);  // for triangle toggle 't' key below
	CP.run(Pgroups); 
  }

  public void keyTyped(KeyEvent e) 
  {
	// inherit methods from CanvasPlotter
	CP.keyTyped(e);


	char key = e.getKeyChar();
	if (key == 't')
	  {   
		plot.showLinesOnTrianglePlot = !plot.showLinesOnTrianglePlot;
		CP.repaint();
	  }
  }
  public void keyReleased(KeyEvent e) {}
  public void keyPressed(KeyEvent e) {}
} 