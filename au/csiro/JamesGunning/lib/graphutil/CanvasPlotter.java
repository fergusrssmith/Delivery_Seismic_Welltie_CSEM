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
import java.text.*;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.Math.min;
import static java.lang.Math.max;

import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*; 
import javax.swing.JFrame;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.graphutil.*;
 
// Classes to get eps output...
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ps.PSGraphics2D;  

 
/** 
	A very simple class to plot 
	groups of common data as line or point data
   
	@author James.Gunning@csiro.au  
	@version 1.0, 25/5/2008
*/ 

public class CanvasPlotter extends ZoomableCanvas implements KeyListener
															 // implements MouseListener, MouseMotionListener, KeyListener,	ComponentListener, ActionListener, Printable
{
  /**
	 Default startup size of a single graph
  */
  public static final int DEFAULT_SINGLE_GRAPH_SIZE=700;
  /**
	 Default maximum dimension of a group-of-graphs (i.e. more than 1)
  */
  public static final int DEFAULT_MAX_SIZE=975;

  public static final Color[] defaultColors = {Color.red, Color.green, Color.blue, Color.magenta, Color.cyan};

  public static final Color[] defaultPalette = {Color.black,            // 0
												Color.white,            // 1
												Color.red,              // 2
												Color.green,            // 3
												Color.blue,             // 4
												Color.magenta,          // 5
												Color.cyan,             // 6
												new Color(101, 67, 33), // 7 (Brown)
												Color.yellow,           // 8 
												Color.orange,           // 9
												Color.pink};            // 10

  public String windowTitle = "";
  public boolean reverseVideo = true; // set white on black


  public boolean windowQuitForcesSystemQuit = false;

  ArrayList<TwoDPlotGroup>_plotgroups;

  ErrorHandler EH;
  public CanvasPlotter()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
	addKeyListener(this); 
  } 

  public CanvasPlotter(boolean enableSystemQuit)
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
	addKeyListener(this); 
	windowQuitForcesSystemQuit = enableSystemQuit;
  }


  public void enableSystemQuit()
  {
	windowQuitForcesSystemQuit = true;
  }

  
  /** 
	  Show example graph array
  */
  public static void main(String argv[]) 
  {
	int i,j,N=30;

	CanvasPlotter CP = new CanvasPlotter();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("CanvasPlotter",CP.EH);

	// parse options. Can be many files
	options.parseOptions(argv,false);

	ArrayList<TwoDPlotGroup>Pgroups = new ArrayList<TwoDPlotGroup>(1);

	CanvasPlotter.TwoDPlotGroup pg = new CanvasPlotter.TwoDPlotGroup(CP);
	Pgroups.add(pg);
	pg.GG.plotGroupNX=1;
	pg.GG.plotGroupXindx=0;
	for (int idata = 0; idata<2; idata++)
	  {
		CanvasPlotter.TwoDPlot p = new CanvasPlotter.TwoDPlot();
		// make fake data
		N = (int)(50*(1.0 + Math.random()));
		double[] x = new double[N];
		double[] y = new double[N];
		double[] e=null;
		if (idata==0) e = new double[N];
		for (i=0;i<N;i++)
		  {
			x[i] = i+1.0;
			y[i] = 5*Math.random();
			if (idata==0) e[i] = Math.abs(y[i]/10.0);
		  }
		p.x = x;
		p.y = y;
		p.e = e;
		p.c = (idata==0) ? Color.green  : Color.red;
		p.symbol = (idata==0) ? 'r': 't';
		p.fillSymbols = true;
		pg.plots.add(p);
		pg.GG.logPlotY = true;
		pg.GG.logPlotX = true;
		pg.GG.computeRanges(x,y);
	  }
	CP.windowTitle = "CanvasPlotterDemo";
	CP.run(Pgroups);
  }

  public void run(ArrayList<TwoDPlotGroup>plotgroups)
  {
	this._plotgroups = plotgroups; 

	// set initial canvas size guess by matrix of plots
	computeCanvasSizes();

	display(windowTitle);
	
	//frame.addWindowListener(this);
  }

  /*
	run, but displaying to supplied file, as eps
  */
  public void runDisplayAsEPS(ArrayList<TwoDPlotGroup>plotgroups, String epsFilename)
  {
	this._plotgroups = plotgroups; 
	
	// set initial canvas size guess by matrix of plots
	computeCanvasSizes();

	display(windowTitle,epsFilename);
	dumpEps(epsFilename);
  }
  
  void computeCanvasSizes()
  {
	// set initial canvas size guess by matrix of plots
	int Nxmax=1,Nymax=1;
	double max;
	TwoDPlotGroup pg;
	for (int i=0;i < _plotgroups.size(); i++)
	  {
		pg = _plotgroups.get(i);
		Nxmax = Math.max(Nxmax, pg.GG.plotGroupNX);
		Nymax = Math.max(Nymax, pg.GG.plotGroupNY);
	  }
	if (Nxmax==1 && (Nxmax==Nymax))
	  {
		canvasWidth = canvasHeight = DEFAULT_SINGLE_GRAPH_SIZE;
	  }
	else if (Nxmax >= Nymax)
	  {
		canvasWidth = DEFAULT_MAX_SIZE;
		canvasHeight = (int) ((canvasWidth*Nymax)/Nxmax);
	  }
	else
	  {
		canvasHeight = DEFAULT_MAX_SIZE;
		canvasWidth = (int) ((canvasHeight*Nxmax)/Nymax);
	  }
	// Try and not go over smaller screens....
	if (true)
	  {
		Toolkit tk = Toolkit.getDefaultToolkit();
		max = tk.getScreenSize().getHeight();

		JFrame dummyFrame = new JFrame();
		dummyFrame.addNotify();
		Insets	ins = dummyFrame.getInsets();
		max -= (ins.bottom + ins.top);
		if (max < canvasHeight) canvasHeight = (int)max-30; // 50 seems needed on OSX/Fink???
		dummyFrame.dispose();
	  }
  }

 public void paint(Graphics g) 
  {
	Graphics2D g2 = (Graphics2D) g;
 
	setBackground(reverseVideo ? Color.black : Color.white);
	//setBackground(Color.gray);

	// save context Xform
	AffineTransform saveXform = g2.getTransform();
    // apply new one
	g2.setTransform(at);
	// then draw away....

	
	// here , the main graphs we want to see
	paintContent(g2);

	// paint Zoom & mouse 3 in old coords
	g2.setTransform(saveXform);

	// zoom box 
	super.paintZoomBox(g2);

  } 

  private void paintContent(Graphics2D g2)
  {
	TwoDPlotGroup pg;
	Color fgColor = reverseVideo ? Color.white : Color.black;

	int Nmax=1;
	for (int i=0;i < _plotgroups.size(); i++)
	  {
		pg = _plotgroups.get(i);
		Nmax = Math.max(1,Math.max(pg.GG.plotGroupNX, pg.GG.plotGroupNY));
	  }

	if (Nmax>1)
	  {
		double scale =  Math.max(0.4,1.9/Nmax);
		double sX,sY,s;
		sX = canvasHeight/(double)DEFAULT_MAX_SIZE;
		sY = canvasWidth/(double)DEFAULT_MAX_SIZE ;
		s = Math.max(sY,sX);
		if (s > 1) scale*=s;

		if (!true)
		  {
			AffineTransform at = new AffineTransform();
			at.scale(scale,scale);		
			g2.setFont(g2.getFont().deriveFont(at));
		  }
		else
		  {
			Font f = g2.getFont();
			g2.setFont(new Font(f.getName(),f.getStyle(),(int)Math.ceil(f.getSize()*(float)scale)));
		  }
	  }

	for (int i=0;i < _plotgroups.size(); i++)
	  {
		pg = _plotgroups.get(i);
		pg.paintContent(g2,fgColor);
	  }
  } 

 
 public void keyTyped(KeyEvent e) 
  {
	char key = e.getKeyChar();
	if (key == 'q')
	  {   
		frame.dispose(); 
		frame=null;
		if (windowQuitForcesSystemQuit) System.exit(0); // ought not to be needed if frames.dispose() used properly
	  }
	else if (key == 'e')
	  {  
		dumpEps();
	  }
  }

  /** 
	  Dump eps version of screen using HEP libraries, with default filename
  */
  public void dumpEps()
  {
	dumpEps(null);
  }
  /** 
	  Dump eps version of screen using HEP libraries, with supplied filename
  */
  public void dumpEps(String filename)
  {
	if (filename==null) filename = new String("EPSdump.ps");
	Properties p = new Properties();
	// p.setProperty("PageSize","A4");
	VectorGraphics g = null;
	try 
	  {
		g = new PSGraphics2D(new File(filename),this); 
	  } 
	catch (IOException ioe) 
	  {
		EH.warn("Couldn't open " + filename + " for eps output");
	  }
	// make sure text comes out as fonts
	p.setProperty(org.freehep.graphicsio.AbstractVectorGraphicsIO.TEXT_AS_SHAPES,"false");
	g.setProperties(p);

	g.startExport(); 
	this.print(g); 
	g.endExport();    
	EH.warn("Wrote EPS screen dump to "+filename);
  }



  // KeyListener
  public void keyReleased(KeyEvent e) {}
  public void keyPressed(KeyEvent e) {}
  
  public static class TwoDPlotGroup
  {
	public ArrayList<TwoDPlotType>plots = new ArrayList<TwoDPlotType>();
	public Rectangle R;
	public GraphGroup GG;
	Canvas _C;

	public TwoDPlotGroup(Canvas C)
	{
	  _C = C;
	  GG = new GraphGroup();
	}

	void paintContent(Graphics2D g2,Color fgColor)
	{
	  double xR,yR,wR,hR;
	  TwoDPlotType plot;

	  // Compute graph pos wrt canvas rectangle
	  Rectangle Rc = _C.getBounds();
	  wR = Rc.width/GG.plotGroupNX;
	  xR = Rc.x + GG.plotGroupXindx * wR;
	  hR = Rc.height/GG.plotGroupNY;
	  yR = Rc.y + GG.plotGroupYindx * hR;
	  
	  R = new Rectangle((int)xR,(int)yR,(int)wR,(int)hR);


	  // common set of axes, labels etc..
	  GG.paintContent(g2,R,fgColor);

	  // and set of data, lines, objects etc in this group...
	  for (int i=0; i < plots.size(); i++)
		{
		  plot = plots.get(i);
		  plot.paintContent(g2,R,GG);
		}
	}
  }

  public static class GraphGroup
  {
	public static final double DEFAULT_FRACTIONAL_MARGIN = 0.1;
	public static int DEFAULT_N_TICKS = 10;
	public static double TICK_SIZE_FRAC = 0.1;

	public double xMin=0,xMax=1,yMin=0,yMax=1;  // properties of data
	public double xMinT=0,xMaxT=1,yMinT=0,yMaxT=1;  // properties of ticks
	public double xMargLeft;
	public double xMargRight;
	public double yMargTop;
	public double yMargBot;

	public String title=null;
	public boolean logPlotX = false;
	public boolean logPlotY = false;
	public boolean reverseY = false;
	public int plotGroupNX=1;
	public int plotGroupNY=1;
	public int plotGroupXindx=0;
	public int plotGroupYindx=0;
	// flags to suppress ticks+words if desired
	public boolean doXticks=true;
	public boolean doYticks=true;

	public String Xword=null, Yword=null;

	boolean haveRanges = false;

	// formatting for numbers
	java.text.DecimalFormat dFormatX,dFormatY;

	// possible needed font scale fudges
	public float fontScale = 0f;
	public float fontTitleScale = 0f;

	TickMarks ticks;
	ErrorHandler EH;

	public GraphGroup()
	{
	  xMargLeft = yMargBot = DEFAULT_FRACTIONAL_MARGIN;
	  xMargRight = yMargTop = DEFAULT_FRACTIONAL_MARGIN/3.0;

	  EH = new ErrorHandler();
	  ticks = new TickMarks(EH);
	  // Format for text numbers for text
	  setTickNumbersEngStyle(true,true);
	}

	public void setTickNumbersEngStyle(boolean x,boolean y)
	{	  
	  dFormatX = new DecimalFormat(x ? "0.#E0" : "0.######");
	  dFormatY = new DecimalFormat(y ? "0.#E0" : "0.######");
	  dFormatX.setMaximumFractionDigits(3);
	  dFormatY.setMaximumFractionDigits(3);
	}

	void updateTicksBracket()
	{
	  double xmin,xmax,ymin,ymax;
	  double tick_spacing;

	  // Xs
	  if (logPlotX)
		{
		  // rount up and down to nearest integer
		  xmin = Math.floor(log10(xMin));
		  xmax = Math.ceil(log10(xMax));
		}
	  else
		{
		  xmin = xMin;
		  xmax = xMax;
		}
	  tick_spacing = ticks.getTickSpacing(xmin,xmax,DEFAULT_N_TICKS);
	  xMinT = ticks.getMinTickCount() * tick_spacing;
	  xMaxT = ticks.getMaxTickCount() * tick_spacing;

	  if (logPlotX)
		{
		  xMinT = pow(10,xMinT);
		  xMaxT = pow(10,xMaxT);
		}
	  

	  // Ys
	  if (logPlotY)
		{
		  ymin = Math.floor(log10(yMin));
		  ymax = Math.ceil(log10(yMax));
		}
	  else
		{
		  ymin = yMin;
		  ymax = yMax;
		}
	  tick_spacing = ticks.getTickSpacing(ymin,ymax,DEFAULT_N_TICKS);
	  yMinT = ticks.getMinTickCount() * tick_spacing;
	  yMaxT = ticks.getMaxTickCount() * tick_spacing;
	  if (logPlotY)
		{
		  yMinT = pow(10,yMinT);
		  yMaxT = pow(10,yMaxT);
		}
	}

	public void clearRanges()
	{
	  haveRanges = false;
	}

	public void setXRange(double xmax,double xmin)
	{
	  xMin = xmin;
	  xMax = xmax;
	  updateTicksBracket();	 
	}

	public void setYRange(double ymax,double ymin)
	{
	  yMin = ymin;
	  yMax = ymax;
	  updateTicksBracket();	 
	}

	public void computeRanges(double[] x, double[] y)
	{
	  int i,n;
	  n = x.length;
	  if (n>0)
		{
		  if (!haveRanges) xMin = xMax = x[0];
		  for (i=0; i < n; i++)
			{ 
			  xMin = Math.min(x[i],xMin);
			  xMax = Math.max(x[i],xMax);
			}
		}
	  n = y.length;
	  if (n>0)
		{
		  if (!haveRanges) yMin = yMax = y[0];
		  for (i=0; i < n; i++)
			{ 
			  yMin = Math.min(y[i],yMin);
			  yMax = Math.max(y[i],yMax);
			}
		}
	  haveRanges = true;
	  updateTicksBracket();	  
	}


	/** 
		Paints axes plus labels, ticks
	*/
	void paintContent(Graphics2D g2, Rectangle R, Color fgColor)
	{
	  int i,n;
	  double X,w,Y,Xl,Yl,x,y,strDim,strDimCmp=0,xN,yN,yBase,zmin,zmax;
	  double aX,bX,aY,bY;
	  double xScale,yScale;
	  double tick_spacing, range;
	  double lastXtickLabel;
	  int min_tick_count,max_tick_count;
	  boolean drawTick;
	  String word;


	  xN = R.width*xMargLeft + R.x;
	  xScale = R.width*(1.0 - xMargLeft - xMargRight);

	  yBase = R.height*(1.0-yMargBot) + R.y;
	  if (!reverseY)
		{
		  yN = R.height*(1.0-yMargBot) + R.y;
		  yScale = R.height*(1.0 - yMargTop - yMargBot);
		}
	  else
		{
		  yN = R.height*yMargTop + R.y;
		  yScale = -R.height*(1.0 - yMargTop - yMargBot);
		}

	  if (fontScale > 0)
		{
		  Font f = g2.getFont();
		  g2.setFont(f.deriveFont(AffineTransform.getScaleInstance(fontScale,fontScale)));
		}

	  g2.setPaint(fgColor);
	  // Axes
	  g2.draw(new Line2D.Double(xN + xScale, yBase,  xN, yBase));
	  g2.draw(new Line2D.Double(xN, yN, xN, yN - yScale));

	  // zero-range protection
	  if (yMaxT <= yMinT) yMaxT = yMinT+1.0;
	  if (xMaxT <= xMinT) xMaxT = xMinT+1.0;

	  // scalars for conversion
	  if (logPlotX)
		{
		  aX = 1.0/(log10(xMaxT) - log10(xMinT));
		  bX = -log10(xMinT) * aX;
		}
	  else
		{
		  aX = 1.0/(xMaxT - xMinT);
		  bX = -xMinT * aX;
		}
	  if (logPlotY)
		{
		  aY = 1.0/(log10(yMaxT) - log10(yMinT));
		  bY = -log10(yMinT) * aY;
		}
	  else
		{
		  aY = 1.0/(yMaxT - yMinT);
		  bY = -yMinT * aY;
		}

	  // X Ticks
	  if (logPlotX)
		{
		  zmin = log10(xMin); 
		  zmax = log10(xMax);
		}
	  else
		{
		  zmin = xMin; 
		  zmax = xMax;
		}
	  // zero-range protection
	  if (zmax > zmin) 
		{
		  tick_spacing = ticks.getTickSpacing(zmin,zmax,DEFAULT_N_TICKS);
		}
	  else
		{
		  // old zmax=zmin+1.0;
		  zmax=zmin;
		  tick_spacing = ticks.getTickSpacing(logPlotX ? log10(xMinT) : xMinT,
											  logPlotX ? log10(xMaxT) : xMaxT,
											  DEFAULT_N_TICKS);
		}

	  min_tick_count = ticks.getMinTickCount();
	  max_tick_count = ticks.getMaxTickCount();
	  
	  // X axis
	  FontMetrics metrics = g2.getFontMetrics();
	  double fHeight = metrics.getAscent();
	  double txtMargin = 0.2 * metrics.getAscent();
	  Y = yBase;
	  Yl = yBase - TICK_SIZE_FRAC * R.height*yMargBot;
	  lastXtickLabel = -Double.MAX_VALUE;

	  if (doXticks)
		{ 
		  for (i=min_tick_count; i<=max_tick_count; i++) 
			{
			  w = tick_spacing * i;
			  //		  X = (X - xMin)/(xMax - xMin);
			  X = w * aX + bX;
			  X = X*xScale + xN;
			  drawTick = (!logPlotX || (logPlotX && (Math.abs(w - Math.rint(w))<0.001)));
			  if (drawTick)
				{
				  g2.draw(new Line2D.Double(X,Y,X,Yl));
				  if (logPlotX) w = pow(10,w);
				  word = dFormatX.format(w);
				  strDim= metrics.stringWidth(word)/2.0;

				  // only draw word if if doesn't collide with last
				  /*
					x = X-strDim;
					if (x > lastXtickLabel+(strDim/2.0)) // force 1/4 word spacing between labels
					{
					g2.drawString(word,(int)x,(int)(Y+2*fHeight));
					lastXtickLabel = X+strDim;
					}
				  */
				  if (i==min_tick_count) strDimCmp = strDim;
				  if ((X-strDimCmp) > lastXtickLabel+(strDimCmp/2.0)) // force 1/4 word spacing between labels
					{
					  g2.drawString(word,(int)(X-strDim),(int)(Y+2*fHeight));
					  lastXtickLabel = X+strDimCmp;
					}
				}
			}
		  // special non-linear ticks for log scales. These always span entire decades.
		  if (logPlotX)
			{
			  Yl = yN - (reverseY ? -1:1) * 0.5*TICK_SIZE_FRAC * R.height*yMargBot;
			  for (i=min_tick_count; i<=max_tick_count; i++) 
				{
				  w = tick_spacing * i;
				  if ((Math.abs(w - Math.rint(w))<0.001) && (i<max_tick_count))
					{
					  // draw intermediate ticks
					  for (int j=0; j<TickMarks.LOG_SMALL_TICKS.length; j++)
						{
						  X = w + TickMarks.LOG_SMALL_TICKS[j];
						  X = X * aX + bX;
						  X = X*xScale + xN;
						  g2.draw(new Line2D.Double(X,Y,X,Yl));
						}
					}
				}
			}
		}
	  // Xword
	  if (Xword!=null)
		{
		  // Axes
		  // yBase = R.height*(1.0-yMargBot) + R.y;
		  //g2.draw(new Line2D.Double(xN + xScale, yBase,  xN, yBase));
		  //g2.draw(new Line2D.Double(xN, yN, xN, yN - yScale));
		  
		  strDim= metrics.stringWidth(Xword)/2.0;
		  X = xN+xScale/2.0 - strDim;
		  Y = R.height + R.y - 0.5*fHeight;  // 1/2 word thickness off bottom edge
		  g2.drawString(Xword,(int)X,(int)Y);
		  //g2.draw(new Line2D.Double(X,Y,X+20,Y));
		  //g2.draw(new Line2D.Double(X,Y,X,Y-20));
		}

	  // Y Ticks
	  if (logPlotY)
		{
		  zmin = log10(yMin); 
		  zmax = log10(yMax);
		}
	  else
		{
		  zmin = yMin; 
		  zmax = yMax;
		}
	  // zero-range protection
	  if (zmax >zmin) 
		{
		  tick_spacing = ticks.getTickSpacing(zmin,zmax,DEFAULT_N_TICKS);
		}
	  else
		{
		  // old zmax=zmin+1.0;
		  zmax=zmin;
		  tick_spacing = ticks.getTickSpacing(logPlotY ? log10(yMinT) : yMinT,
											  logPlotY ? log10(yMaxT) : yMaxT,
											  DEFAULT_N_TICKS);
		}
	  min_tick_count = ticks.getMinTickCount();
	  max_tick_count = ticks.getMaxTickCount();
	  X = xN;
	  Xl = xN + TICK_SIZE_FRAC * R.width * xMargLeft;
	  if (doYticks)
		{
		  for (i=min_tick_count; i<=max_tick_count; i++) 
			{
			  Y = w = tick_spacing * i;
			  // Y = (Y - yMin)/(yMax - yMin):
			  Y = Y * aY + bY;
			  // now 0<Y<1
			  // map into internal rectangle
			  // 
			  Y = -Y*yScale + yN;
			  drawTick = (!logPlotY || (logPlotY && (Math.abs(w - Math.rint(w))<0.001)));
			  if (drawTick)
				{
				  g2.draw(new Line2D.Double(Xl,Y,X,Y));
				  if (logPlotY) w = pow(10,w);
				  word = dFormatY.format(w);
				  g2.drawString(word,(int)(X-metrics.stringWidth(word)-txtMargin),(int)(Y + 0.5*fHeight));
				}
			}
		  // special non-linear ticks for log scales. These always span entire decades.
		  if (logPlotY)
			{
			  Xl = xN + 0.5*TICK_SIZE_FRAC * R.width * xMargLeft;
			  for (i=min_tick_count; i<=max_tick_count; i++) 
				{
				  w = tick_spacing * i;
				  if ((Math.abs(w - Math.rint(w))<0.001) && (i<max_tick_count))
					{
					  // draw intermediate ticks
					  for (int j=0; j<TickMarks.LOG_SMALL_TICKS.length; j++)
						{
						  Y = w + TickMarks.LOG_SMALL_TICKS[j];
						  Y = Y * aY + bY;
						  Y = -Y*yScale + yN;
						  g2.draw(new Line2D.Double(Xl,Y,X,Y));
						}
					}
				}
			}
		}
	  // Yword
	  if (Yword!=null)
		{
		  strDim= metrics.stringWidth(Yword)/2.0;
		  AffineTransform aT,aTrot;
		  aT = g2.getTransform();
		  aTrot = g2.getTransform(); // copy of aT
		  g2.setTransform(aTrot);
		  // rotate about point
		  X = R.x + 1.1*fHeight;
		  Y = yN-yScale/2.0 + strDim;
		  //g2.draw(new Line2D.Double(X,Y,X+20,Y));
		  //g2.draw(new Line2D.Double(X,Y,X,Y-20));

		  g2.rotate(-Math.PI/2.0,(int)X,(int)Y);
		  g2.drawString(Yword,(int)X,(int)Y);
		  // set back to normal
		  g2.setTransform(aT);
		}

	  // Title
	  if (title!=null)
		{
		  FontMetrics metrics2 = metrics;
		  Font f = null;
		  double fHeight2 = fHeight;
		  if (fontTitleScale > 0)
			{
			  f = g2.getFont();
			  g2.setFont(f.deriveFont(AffineTransform.getScaleInstance(fontTitleScale,fontTitleScale)));
			  metrics2 = g2.getFontMetrics();
			  fHeight2 = metrics.getAscent();
			}
		  i = (int)(R.x+R.width/2.0-metrics2.stringWidth(title)/2.0);
		  int j = (int)(R.y + fHeight + 0.5*R.height*yMargTop);
		  g2.drawString(title,i,j);

		  if (fontTitleScale > 0) g2.setFont(f);
		}

	}
  }

  public interface TwoDPlotType
  {
	public void paintContent(Graphics2D g2, Rectangle R,GraphGroup GG);	
  }

  /**
	 Class for a line plot, with possible error bars, colors, symbols etc
  */
  public static class TwoDPlot implements TwoDPlotType
  {
	public static double SYMBOL_SIZE = 0.01;
	public double x[];
	public double y[];
	public double e[]=null;

	public Color c;
	public boolean drawline = true;
	public float lineWidth = 1.0f;
	public boolean fillSymbols= true;
	public double symSize = SYMBOL_SIZE;
	public char symbol='\0';
	public String label=null;


	double xScale,yScale,tmp;
	Rectangle _R;
	GraphGroup _GG;
	Graphics2D _g2;

	Stroke defaultStroke  = new BasicStroke(1,BasicStroke.CAP_ROUND,BasicStroke.JOIN_MITER);

	public void paintContent(Graphics2D g2, Rectangle R,GraphGroup GG)
	{
	  byte eP,eM;
	  int i,n;
	  double X,Y,YeP,YeM,X0=0,Y0=0,xN,yN,xx,yy;
	  double aX,bX,aY,bY;

	  _R = R;
	  _GG = GG;
	  _g2 = g2;

	  if (x.length != y.length)
		{
		  throw new RuntimeException("x[] and y[] different lengths in TwoDPlot.paintContent()");
		}
	  if ((e!=null) && (e.length != x.length))
		{
		  throw new RuntimeException("x[] and e[] different lengths in TwoDPlot.paintContent()");
		}
	  xN = R.width*_GG.xMargLeft + R.x;
	  xScale = R.width*(1.0 - _GG.xMargLeft - _GG.xMargRight);
	  if (!_GG.reverseY)
		{
		  yN = R.height*(1.0-_GG.yMargBot) + R.y;
		  yScale = R.height*(1.0 - _GG.yMargTop - _GG.yMargBot);
		}
	  else
		{
		  yN = R.height*_GG.yMargTop + R.y;
		  yScale = -R.height*(1.0 - _GG.yMargTop - _GG.yMargBot);
		}


	  n = Math.min(x.length,y.length);

	  Paint oldPaint = g2.getPaint();
	  // Stroke oldStroke = g2.getStroke();
	  Stroke oldStroke  = defaultStroke;
	  g2.setPaint(c);
	  g2.setStroke((lineWidth==1.0f) ? oldStroke : new BasicStroke(lineWidth,BasicStroke.CAP_ROUND,BasicStroke.JOIN_MITER));


	  if (GG.logPlotX)
		{
		  aX = 1.0/(log10(GG.xMaxT) - log10(GG.xMinT));
		  bX = -log10(GG.xMinT) * aX;
		}
	  else
		{
		  aX = 1.0/(GG.xMaxT - GG.xMinT);
		  bX = -GG.xMinT * aX;
		}
	  if (GG.logPlotY)
		{
		  aY = 1.0/(log10(GG.yMaxT) - log10(GG.yMinT));
		  bY = -log10(GG.yMinT) * aY;
		}
	  else
		{
		  aY = 1.0/(GG.yMaxT - GG.yMinT);
		  bY = -GG.yMinT * aY;
		}

	  //  xx = X*xScale + xN;   for 0<X<1 is the span of screen xx coordinates
	  //  yy = yN - Y*yScale;   for 0<Y<1 is the span of screen yy coordinates
	  // Don't draw outside the axes...
	  Rectangle2D.Double clip = new Rectangle2D.Double(xN,GG.reverseY ? yN : yN-yScale,Math.abs(xScale),Math.abs(yScale));
	  g2.setClip(clip);


	  if (label!=null)
		{
		  // label at "midpoint"
		  i=x.length/2;
		  // X = (x[i] - GG.xMin)/(GG.xMax - GG.xMin);         //  0 <= X <= 1
		  X = (_GG.logPlotX ? log10(x[i]) : x[i]) * aX + bX;
		  // Y = (y[i] - GG.yMin)/(GG.yMax - GG.yMin);
		  Y = (_GG.logPlotY ? log10(y[i]) : y[i]) * aY + bY;
		  xx = X*xScale + xN;
		  yy = yN - Y*yScale;
		  g2.drawString(label,(int)xx,(int)yy);
		}
	  for (i=0;i<x.length;i++)
		{
		  // X = (x[i] - GG.xMin)/(GG.xMax - GG.xMin);         //  0 <= X <= 1
		  X = (_GG.logPlotX ? log10(x[i]) : x[i]) * aX + bX;
		  // Y = (y[i] - GG.yMin)/(GG.yMax - GG.yMin);
		  Y = (_GG.logPlotY ? log10(y[i]) : y[i]) * aY + bY;
		  xx = X*xScale + xN;
		  yy = yN - Y*yScale;
		  if (drawline && (i>0))
			{			  
			  g2.draw(new Line2D.Double(X0*xScale + xN, yN - Y0*yScale,xx, yy));
			}	
		  if (symbol != '\0')
			{
			  paintSymbol(xx,yy);
			}
		  if (e!=null)
			{
			  eP=eM=0;
			  if (_GG.logPlotY)
				{
				  if ((tmp = y[i]+e[i]) > 0)
					{
					  YeP = log10(tmp)  * aY + bY;
					}
				  else
					{
					  YeP=0;
					  eP=1;
					}
				  if ((tmp = y[i]-e[i])>0)
					{
					  YeM = log10(tmp)  * aY + bY;
					}
				  else
					{
					  YeM=0;
					  eM=1;
					}
				}
			  else
				{
				  YeP = (y[i]+e[i]) * aY + bY;
				  YeM = (y[i]-e[i]) * aY + bY;
				}
			  if (YeP > 1.0) {YeP=1;eP=1;}
			  if (YeP < 0) {YeP=0;eP=1;}
			  if (YeM > 1.0) {YeM=1;eM=1;}
			  if (YeM < 0) {YeM=0;eM=1;}
	  
			  // omit off-scale error-bars ??
			  // if (eP!=0 || eM!=0)
				{	

				  YeP = yN - YeP*yScale;
				  YeM = yN - YeM*yScale;

				  paintErrorSymbol(xx,yy,YeP,eP,YeM,eM);
				}
			}
		  X0=X;
		  Y0=Y;
		}
	  // restore Stroke & Paint
	  g2.setStroke(oldStroke);
	  g2.setPaint(oldPaint);

	  // clear clip
	  g2.setClip(null);
	}

	/**
	   Symbol painter. Currently accepts symbol=
	   'r' rectangle
	   'c' circle
	   't' triangle
	*/
	void paintSymbol(double x,double y)
	{
	  double X,Y,wx,wy;
	  Shape s=null;

	  wx=xScale*symSize;
	  wy=yScale*symSize;
	  switch (symbol)
		{
		case 'r':
		  {
			s = new Rectangle2D.Double(x-wx/2.0, y-wy/2.0, wx,wy);
			break;
		  }
		case 'c':
		  {
			s = new Ellipse2D.Double(x-wx/2.0, y-wy/2.0, wx,wy);
			break;
		  }
		case 't':
		  {
			// Fill polygons (Triangles using GeneralPath)
			GeneralPath p = new GeneralPath(GeneralPath.WIND_EVEN_ODD);
			p.moveTo(x-0.5*wx, y+wy*0.288);
			p.lineTo(x,        y-wy*0.711);
			p.lineTo(x+0.5*wx, y+wy*0.288);
			s=p;
			break;
		  }
		default:
		  break;
		}
	  if (s!=null)
		{
		  if (fillSymbols)
			{
			  _g2.fill(s);
			}
		  else
			{
			  _g2.draw(s);
			}
		}
	}
	
	void paintErrorSymbol(double x,double y,double yeP,byte eP,double yeM,byte eM)
	{
	  double X,Y,wx,wy;
	  Shape s=null;
	  
	  wx=xScale*symSize/2.0;
	  
	  _g2.draw(new Line2D.Double(x,yeP,x,yeM));
	  if (eP==0) _g2.draw(new Line2D.Double(x-wx,yeP,x+wx,yeP));
	  if (eM==0) _g2.draw(new Line2D.Double(x-wx,yeM,x+wx,yeM));
	}
  }

  /**
	 Class for drawing coloured rectangles or triangles etc in place
  */
  public static class TwoDShapePlot implements TwoDPlotType
  {
	public boolean showLinesOnTrianglePlot = true;

	boolean isTrianglePlot = false;


	/** 
		if (!isTrianglePlot)
		    x,x1,y,y1 = left & right edges of box, top & bottom of boxes. x2  y2 unused
		else 
		    x,y  x1,y1   x2,y2 = vertices of triangles
	*/
	public double x[],x1[],x2[];
	public double y[],y1[],y2[];
	/* box color */
	public double col[];
	/* standard "library" color override if gcol[i]>=0. From list in defaultPalette above */
	public int gcol[];

	double minCol=0,maxCol=1;

	private float[] minColComps={0,0,0,255};
	private float[] maxColComps={255,255,255,255};

	private Color maxColor=null;

	public TwoDShapePlot() {} // default

	public TwoDShapePlot(boolean isTrianglePlot,boolean showLinesOnTrianglePlot)
	{
	  this.isTrianglePlot = isTrianglePlot;
	  this.showLinesOnTrianglePlot = showLinesOnTrianglePlot;
	}

	public void paintContent(Graphics2D g2, Rectangle R,GraphGroup GG)
	{
	  int i,ci,n,pass;
	  double X,Y,X1,Y1,xN,yN,c=0,xx,yy,xx1,yy1;
	  double aX,bX,aY,bY;
	  double xScale,yScale;

	  double X2,Y2,xx2=0,yy2=0; // triangles only
	  Stroke oldStroke,lStroke=null;

	  xN = R.width*GG.xMargLeft + R.x;
	  xScale = R.width*(1.0 - GG.xMargLeft - GG.xMargRight);



	  if (!GG.reverseY)
		{
		  yN = R.height*(1.0-GG.yMargBot) + R.y;
		  yScale = R.height*(1.0 - GG.yMargTop - GG.yMargBot);
		}
	  else
		{
		  yN = R.height*GG.yMargTop + R.y;
		  yScale = -R.height*(1.0 - GG.yMargTop - GG.yMargBot);
		}

	  //  xx = X*xScale + xN;   for 0<X<1 is the span of screen xx coordinates
	  //  yy = yN - Y*yScale;   for 0<Y<1 is the span of screen yy coordinates
	  // Don't draw outside the axes...
	  Rectangle2D.Double clip = new Rectangle2D.Double(xN,GG.reverseY ? yN : yN-yScale,Math.abs(xScale),Math.abs(yScale));
	  g2.setClip(clip);

	  n = Math.min(x.length,y.length);

	  Paint oldPaint = g2.getPaint();
	  oldStroke = g2.getStroke();


	  if (GG.logPlotX)
		{
		  aX = 1.0/(log10(GG.xMaxT) - log10(GG.xMinT));
		  bX = -log10(GG.xMinT) * aX;
		}
	  else
		{
		  aX = 1.0/(GG.xMaxT - GG.xMinT);
		  bX = -GG.xMinT * aX;
		}
	  if (GG.logPlotY)
		{
		  aY = 1.0/(log10(GG.yMaxT) - log10(GG.yMinT));
		  bY = -log10(GG.yMinT) * aY;
		}
	  else
		{
		  aY = 1.0/(GG.yMaxT - GG.yMinT);
		  bY = -GG.yMinT * aY;
		}

	  if (isTrianglePlot && showLinesOnTrianglePlot)
		{
		  oldStroke  = g2.getStroke();
		  lStroke = new BasicStroke(.1f);		  
		}

	  for (pass=0; pass< ((isTrianglePlot && showLinesOnTrianglePlot) ? 2 : 1); pass++)
		{
		 if (pass==1)
		   {
			 // isTrianglePlot && showLinesOnTrianglePlot
			 g2.setPaint(Color.red);
			 g2.setStroke(lStroke);
		   }
		 for (i=0;i<x.length;i++)
			{
			  // X = (x[i] - GG.xMin)/(GG.xMax - GG.xMin);         //  0 <= X <= 1
			  X = (GG.logPlotX ? log10(x[i]) : x[i]) * aX + bX;
			  // Y = (y[i] - GG.yMin)/(GG.yMax - GG.yMin);
			  Y = (GG.logPlotY ? log10(y[i]) : y[i]) * aY + bY;
			  xx = X*xScale + xN;
			  yy = yN - Y*yScale;

			  X1 = (GG.logPlotX ? log10(x1[i]) : x1[i]) * aX + bX;
			  // Y = (y[i] - GG.yMin)/(GG.yMax - GG.yMin);
			  Y1 = (GG.logPlotY ? log10(y1[i]) : y1[i]) * aY + bY;
			  xx1 = X1*xScale + xN;
			  yy1 = yN - Y1*yScale;

			  if (isTrianglePlot)
				{	
				  X2 = (GG.logPlotX ? log10(x2[i]) : x2[i]) * aX + bX;
				  // Y = (y[i] - GG.yMin)/(GG.yMax - GG.yMin);
				  Y2 = (GG.logPlotY ? log10(y2[i]) : y2[i]) * aY + bY;
				  xx2 = X2*xScale + xN;
				  yy2 = yN - Y2*yScale;
				}

			  // DRAW SHAPES HERE
			  if (pass==0)
				{
				  if (gcol[i] < 0)
					{
					  c = (col[i]-minCol)/(maxCol-minCol);
					  g2.setPaint(blendColor((float)max(min(c,1.0),0.0)));
					}
				  else
					{
					  // override with prescribed color from palette
					  g2.setPaint(defaultPalette[gcol[i]]); 
					}
				}
			  if (!isTrianglePlot)
				{
				  g2.fill(new Rectangle2D.Double(xx,yy,xx1-xx,Math.abs(yy1-yy)));
				}
			  else 
				{
				  if (pass==0)
					{
					  // Fill polygons (Triangles using GeneralPath)
					  GeneralPath p = new GeneralPath(GeneralPath.WIND_EVEN_ODD);
					  p.moveTo(xx,yy);
					  p.lineTo(xx1,yy1);
					  p.lineTo(xx2,yy2);
					  g2.fill(p);
					}
				  else
					{
					  g2.draw(new Line2D.Double(xx,yy,xx1,yy1));
					  g2.draw(new Line2D.Double(xx1,yy1,xx2,yy2));
					  g2.draw(new Line2D.Double(xx2,yy2,xx,yy));
					}
				}
			}
		}
	  // restore Stroke & Paint
	  g2.setStroke(oldStroke);
	  g2.setPaint(oldPaint);
	  // clear clip
	  g2.setClip(null);
	}

	/**
	   return color that is a blend-mixture of the minColor and
	   maxColor endpoints set by setColorScheme(), where
	   0 <= blend <= 1
	*/
	private Color blendColor(float blend)
	{
	  float x;

	  x = (1f-blend);

	  // darken slightly for white background
	  if (maxColor == Color.white) blend *= 0.9;

	  return new Color(blend*maxColComps[0] + x*minColComps[0],
					   blend*maxColComps[1] + x*minColComps[1],
					   blend*maxColComps[2] + x*minColComps[2]);
	}

	public void setColorScheme(Color mincolor,Color maxcolor)
	{ 
	  maxColor = maxcolor;
	  maxcolor.getComponents(maxColComps);
	  mincolor.getComponents(minColComps);
	}

	/* normalize color range */
	public void normalize()
	{
	  int i;

	  if (col.length < 1) return;
	  maxCol = minCol = col[0];
	  for (i=1;i<col.length;i++)
		{
		  maxCol = Math.max(maxCol,col[i]);
		  minCol = Math.min(minCol,col[i]);
		}
	  if (minCol==maxCol) maxCol = (1+1e-14)*minCol+1e-14;
	}

	public void setColorRange(double minColor,double maxColor)
	{
	  minCol = minColor;
	  maxCol = maxColor;
	}
  }
}
