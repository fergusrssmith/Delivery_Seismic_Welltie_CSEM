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
import java.awt.geom.*;
import java.awt.event.*;

import cern.colt.map.*;
import cern.colt.list.*;
import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.graphutil.*;


// Classes to get eps output...
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ps.PSGraphics2D;  


/**
   A very simple class to show a 2D grayscale rendering of 
   CSEM1DInverter profile files 

   Supplied commandline file CMP.txt loosk like, e.g.

CMP01       636350.29      4104460.55         3129.34
CMP02       635389.21      4106856.89         3233.40
CMP03       634573.47      4108808.22         3511.20
CMP04       633753.50      4110806.57         4575.51
CMP05       633114.29      4112423.09         4699.05
CMP06       632337.02      4114440.49         3559.32
CMP07       631657.72      4116250.92         3555.41
CMP08       630925.01      4118089.26         3832.56

Code extracts profile_final.txt from directories
CMP01, CMP02 etc and assembles 2D graph

   @author James.Gunning@csiro.au
   @version 1.0, 25/5/2008
*/

public class Plot1DCSEMSection extends ZoomableCanvas implements KeyListener
// implements MouseListener, MouseMotionListener, KeyListener,	ComponentListener, ActionListener, Printable
{
  public static int DEFAULT_N_Z_TICKS = 10;
  java.text.DecimalFormat dFormat,cbFormat;

  TickMarks ticks,ticksCB;

  String geom_file;
  String pFilename = "profile_final.txt";
  boolean noRMS=false;
  ErrorHandler EH;
  ArrayList<Geom>geoms;
  float zMax,zMin,vMax,vMin;
  float line_max;

  public Plot1DCSEMSection()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler(); 
	addKeyListener(this);

	// formatting for numbers
	dFormat = new java.text.DecimalFormat("###.##");
	dFormat.setMaximumFractionDigits(3);
	// colorbar
	cbFormat = new java.text.DecimalFormat("###.##");
	cbFormat.setMaximumFractionDigits(1);
	// tickmarks
	ticksCB = new TickMarks(EH);
	ticks = new TickMarks(EH);
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


	Plot1DCSEMSection PO = new Plot1DCSEMSection();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("Plot1DCSEMSection",PO.EH);
	options.addSpielText("Plots 2D section of inversion using 4-column master file CMP.txt to collect inversion info. "+
						 "CMP.txt is a simple ascii file comprising lines of"+
						 "%nDIR_1 X1 Y1 Z1 %nDIR_2 X2 Y2 Z2 %n..."+
						 "%nwhere Xi Yi are the (x,y) locations of the CMP inversion performed in directory DIR_i. "+
						 "Models and responses are collected from the files DIR_i/profile_final.txt and DIR_i/iter_final");


   	options.addMandatoryArg(1,"CMP.txt","CMP geom file");	
	options.addOption("-p",null,1,"String","profile_filename","alternative profile file to seek in each CMP directory");
	options.addOption("--no-RMS",null,0,null,null,"don't show RMS values");

	// parse options. Can be many files
	if (!options.parseOptions(argv,true)) options.Usage();

	PO.geom_file = (options.getAllMandatoryArgs())[0];

	if (options.isOptionSet("-p")) PO.pFilename = ((String)options.getOptionValues("-p")[0]);
	if (options.isOptionSet("--no-RMS")) PO.noRMS=true;
	  

	PO.run();
  }

  void run()
  {
	readGeom();
	readData();
	display("Plot1DCSEMSection");
  }

 public void paint(Graphics g) 
  {
	Graphics2D g2 = (Graphics2D) g;
 

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
	int i,j,Nx,Ny,c;
	float x,x1,dl,dr,y=0,y1,ws,hs,w,h,ymin=0f,ymax=0f;
	double z,y_z;
	Geom g,gR,gL;
	String word;

	// Here, scale everything to canvasHeight, canvasWidth
	// top LH corner
	hs = 0.15f*canvasHeight; // down a bit to leave room for words...
	ws = 0.1f*canvasWidth;
	// image size
	h = 0.8f*canvasHeight;
	w = 0.8f*canvasWidth;

	// scale Fonts to size
	double sX,sY,scale;
	sX = canvasHeight/(double)DEFAULT_SIZE;
	sY = canvasWidth/(double)DEFAULT_SIZE ;
	scale = Math.max(sY,sX);
	Font f = g2.getFont();
	g2.setFont(new Font(f.getName(),f.getStyle(),(int)Math.ceil(f.getSize()*(float)scale)));

	setBackground(Color.white);

	// outside box
	g2.drawRect(rint(ws),rint(hs),rint(w),rint(h));

	Nx = geoms.size();
	for (i=0;i<Nx;i++)
	  {
		g = geoms.get(i);
		Ny = g.z.length;

		gR = (i < (Nx-1)) ? (geoms.get(i+1)) : g;
		gL = (i>0) ? (geoms.get(i-1)) : g;
		dr = (float)(0.5 * (gR.line_dist - g.line_dist));
		dl = (float)(0.5 * (g.line_dist - gL.line_dist));
		x = ws + w * (g.line_dist-dl)/line_max;
		x1 = ws + w * (g.line_dist+dr)/line_max;			
		
		for(j=0;j<(Ny-1);j++)
		  {
			y = hs + h*(g.z[j] - zMin)/(zMax-zMin);
			y1 = hs + h*(g.z[j+1] - zMin)/(zMax-zMin);
			c =(int)( 254.9999f*((g.v[j]-vMin)/(vMax-vMin)));
			g2.setPaint(new Color(c,c,c));
			g2.fill(new Rectangle2D.Float(x,y,x1-x,y1-y));
			g2.setPaint(Color.black);
   			g2.draw(new Line2D.Float(x,y,x,y1));
			g2.draw(new Line2D.Float(x,y1,x1,y1));
			g2.draw(new Line2D.Float(x1,y1,x1,y));
			g2.draw(new Line2D.Float(x1,y,x,y));
		  }
	  }


	double tick_spacing = ticks.getTickSpacing(zMin,zMax,DEFAULT_N_Z_TICKS);
	double range = ticks.getRange();
	int min_tick_count = ticks.getMinTickCount();
	int max_tick_count = ticks.getMaxTickCount();


	// depth scale
	FontMetrics metrics = g2.getFontMetrics();
	double fHeight = metrics.getAscent()/2.0;
	double txtMargin = metrics.getAscent()/5.0;
	for (i=min_tick_count; i<=max_tick_count; i++) 
	  {
		z = tick_spacing * i;
		//y_z =  zMin + h*(z-zMin)/(zMax-zMin); // OLD BUG :-(
		y_z =  hs + h*(z-zMin)/(zMax-zMin);
		if (y_z < hs|| y_z > (h+hs)) continue;
		// coloured time dashes
		g2.setPaint(Color.blue);
		g2.draw(new Line2D.Double(ws,y_z,ws+w,y_z));
		word = dFormat.format(z);
		// black words
		g2.setPaint(Color.black);
		g2.drawString(word,(int)(ws-metrics.stringWidth(word)-txtMargin),(int)(y_z+fHeight));
	  }

	// CMP words
	g2.setPaint(Color.black);
	AffineTransform aT = g2.getTransform();
	AffineTransform aTrot;
	for (i=0;i<Nx;i++)
	  {
		g = geoms.get(i);
		gR = (i < (Nx-1)) ? (geoms.get(i+1)) : g;
		gL = (i>0) ? (geoms.get(i-1)) : g;
		dr = (float)(0.5 * (gR.line_dist - g.line_dist));
		dl = (float)(0.5 * (g.line_dist - gL.line_dist));
		x = ws + w * (g.line_dist-dl)/line_max;
		x1 = ws + w * (g.line_dist+dr)/line_max;		
		word = g.name;
		if (!noRMS) word = word + ", "+dFormat.format(g.misfitRMS);
		// rotate about tick point
		aTrot = g2.getTransform(); // copy of aT
		g2.setTransform(aTrot);
		// x = (float)((x+x1)/2.0 - metrics.stringWidth(word)-txtMargin);
		// x = (float)((x+x1)/2.0);
		g2.rotate(-Math.PI/2.0,(int)(x+x1)/2.0,(int)(hs));
		g2.drawString(word,(int)((x+x1)/2.0 + 5),(int)(hs));			
		//g2.drawString("A"+i,(int)x,(int)hs);			
		// set back to normal
		g2.setTransform(aT);
	  }

	// colorbar legend
	float cb_start,cb_height,cb_block_ht,cb_marg,xR,wordH;

	tick_spacing = ticksCB.getTickSpacing(vMin,vMax,15);
	min_tick_count = ticksCB.getMinTickCount();
	max_tick_count = ticksCB.getMaxTickCount();  

	xR = ws+w;
	cb_marg = 0.15f*(canvasWidth - xR);
	x = xR + cb_marg;
	x1 = canvasWidth - cb_marg;

	Ny = max_tick_count - min_tick_count+1;
	cb_start = canvasHeight/4.0f;
	cb_height = canvasHeight/2.0f;
	cb_block_ht = cb_height/Ny;
	wordH = metrics.getHeight();
	for (i=min_tick_count; i<=max_tick_count; i++) 
	  {
		z =  tick_spacing * i;
		y_z = (z-vMin)/(vMax-vMin);
		c =(int)(254.9999f * y_z);
		if (c<0 || c>255) continue;
		g2.setPaint(new Color(c,c,c));
		y = cb_height*(i-min_tick_count)/(float)Ny + cb_start;
		g2.fill(new Rectangle2D.Float(x,y,x1-x,cb_block_ht));
		//		word = cbFormat.format(Math.pow(10.0,z));
		word = cbFormat.format(z);
		g2.setPaint(Color.black);
		g2.drawString(word,(int)((x+x1 - metrics.stringWidth(word))/2.0),
								 (int)(y+(wordH+cb_block_ht)/2.0));
		if (ymin<=0) ymin = y;
	  }
	ymax = y+cb_block_ht;

	// legend box
	g2.draw(new Line2D.Double(x,ymin,x,ymax));
	g2.draw(new Line2D.Double(x,ymax,x1,ymax));
	g2.draw(new Line2D.Double(x1,ymax,x1,ymin));
	g2.draw(new Line2D.Double(x1,ymin,x,ymin));
  }

 void readGeom()
  {	
	int j,lineNo;
	String s,word,filename = geom_file;
	StringTokenizer st;
	float a=0;
	BufferedReader in = null;

	try
	  {
		in = new BufferedReader(new FileReader(filename));
	  }
	catch (Exception e) 
	  {
		EH.die_err("Can't open file "+filename);
	  }
	geoms = new ArrayList();
	lineNo=0;
	try
	  {
		while ((s=in.readLine()) != null)
		  {
			Geom g = new Geom();
			j=0;
			st = new StringTokenizer(s);
			while (st.hasMoreTokens()) 
			  {
				lineNo++;
				word = st.nextToken();
				if (j>0) 
				  {
					if (j>=4) throw new RuntimeException("Too many numbers, line "+lineNo+" of "+filename);
					try
					  {
						a = Float.parseFloat(word);						
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
					  }
				  }
				switch (j)
				  {
				  case 0:
					g.name = word; break;
				  case 1:
					g.x = a; break;
				  case 2:
					g.y = a; break;
				  case 3:
					g.u = a; break;
				  }
				j++;				
			  }
			geoms.add(g);
		  }
		in.close();
	  }
	catch (Exception e) 
	  {
		EH.die_err("Can't close "+filename);
	  }
	// sort into line order by x,y
	Collections.sort(geoms);
  }

  void readData()
  {
	int i,j,col,lineNo,nGood;
	String s,word,filename;
	StringTokenizer st;
	Geom g;
	BufferedReader in=null;
	col=1;
	float X0,Y0,a=0f;
	FloatArrayList z,v;

	// loop over files
	for (i = 0;i< geoms.size();i++)
	  {
		g = geoms.get(i);
		filename = g.name + File.separator + pFilename;
		// EH.println("Reading "+filename);
		try
		  {
			in = new BufferedReader(new FileReader(filename));
		
			z = new FloatArrayList();
			v = new FloatArrayList();
			lineNo=0;
			while ((s=in.readLine()) != null)
			  {
				j=0;
				st = new StringTokenizer(s);
				while (st.hasMoreTokens()) 
				  {
					lineNo++;
					if (j>=2) throw new RuntimeException("Too many numbers, line "+lineNo+" of "+filename);
					word = st.nextToken();
					try
					  {
						a = Float.parseFloat(word);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
					  }	
					switch (j)
					  {
					  case 0:
						z.add(a); break;
					  case 1:
						v.add(a); break;
					  }
					j++;	
				  }
			  }
			in.close();
			z.trimToSize();
			v.trimToSize();
			g.z = z.elements();
			g.v = v.elements();
		  }
		catch (FileNotFoundException  e) 
		  {
			EH.warn("Can't open file "+filename+" inventing dummy data");
			g.z = new float[2];
			g.v = new float[2];
			g.isGood = false;
		  }
		catch (Exception e) 
		  {
			EH.die_err("Exception e="+e);
		  }
		
		// Now read iter_final files for tlast

		filename = g.name+"/iter_final";
		// EH.println("Reading "+filename);  
		try
		  {
			in = new BufferedReader(new FileReader(filename));

			lineNo=0;
			OUT:
			while ((s=in.readLine()) != null)
			  {
				if (s.startsWith("tlast:"))
				  {
					st = new StringTokenizer(s);
					word = st.nextToken();
					if (word.equals("tlast:") && st.hasMoreTokens())
					  {
						word = st.nextToken();
						try
						  {
							g.misfitRMS = Float.parseFloat(word);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Line "+lineNo+" of file " + filename+" not tlast: number");
						  }	
						break OUT;
					  }
					else
					  {
						EH.die_err("File " + filename+" doesn't have \"tlast: number\" line");
					  }
				  }
			  }
			in.close();
		  }
		catch (FileNotFoundException  e) 
		  {
			EH.warn("Can't open file "+filename+" inventing dummy RMS");
			g.misfitRMS=0.0f;
		  }
		catch (Exception e) 
		  {
			EH.die_err("Exception e="+e);
		  }
	  }
	zMax = vMax = -Float.MAX_VALUE;
	zMin = vMin = Float.MAX_VALUE;
	if (geoms.size()<1)
	  {
		EH.die_err("Can't read CMP directories to get data");
	  }
	
	g = geoms.get(0);
	X0 = g.x;
	Y0 = g.y;
	g = geoms.get(geoms.size()-1);
	line_max = (float)Math.hypot(g.y - Y0,g.x-X0);
	nGood=0;
	for (i=0;i<geoms.size();i++)
	  {
		g = geoms.get(i);
		if (g.isGood)
		  {
			for (j=0;j<g.z.length;j++)
			  {
				zMax = Math.max(g.z[j],zMax);
				zMin = Math.min(g.z[j],zMin);
				vMax = Math.max(g.v[j],vMax);
				vMin = Math.min(g.v[j],vMin);
			  }
			nGood++;
		  }
		g.line_dist = (float)Math.hypot(g.y - Y0, g.x-X0);
		line_max = (float)Math.max(g.line_dist, line_max);		
	  }
	if (nGood<1)
	  {
		EH.die_err("No good inversion to show");
	  }
	else
	  {
		for (i=0;i<geoms.size();i++)
		  {
			g = geoms.get(i);
			EH.println("CMP: "+g.name);
			if (!(g.isGood))
			  {	
				g.z[0]=zMin;
				g.z[1]=zMax;
			  }
		  }
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

  public void keyTyped(KeyEvent e) 
  {
	char key = e.getKeyChar();
	if (key == 'q')
	  {
		System.exit(0);
	  }
	else if (key == 'e')
	  {
		dumpEps();
	  }
  }

  public void keyReleased(KeyEvent e) {}
  public void keyPressed(KeyEvent e) {}

 
  static final int rint(float x)
  {
	return (int)Math.rint((double)x);
  }

  public static class Geom implements Comparable
  {
	boolean isGood=true;
	String name;
	float x,y,u;
	float misfitRMS=-1.0f;

	// distance from leftmost CMP x,y
	float line_dist;

	// data
	float[] z;
	float[] v;

	/** sort by x, then y */
	public int compareTo(Object o) 
	{ 
	  int i,l;
	
	  Geom g1 = this;
	  Geom g2 = (Geom)o;
	
	  if (g1.x != g2.x)
		{
		  return (g1.x < g2.x) ? -1 : 1;
		}
	  else
		{
		  return (g1.y <= g2.y) ? -1 : 1;
		}
	}
  }
}