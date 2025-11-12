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


package au.csiro.JamesGunning.waveletextraction;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.*;
import java.awt.print.PrinterJob;
import java.awt.geom.*;
import javax.swing.*;
import java.awt.print.*;

import java.awt.image.*;
import javax.imageio.ImageIO;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.sulib.*;  
import au.csiro.JamesGunning.lib.graphutil.*; 
import au.csiro.JamesGunning.lib.Numerical.*;
 
import cern.colt.list.*; 

// Classes to get eps output...
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ps.PSGraphics2D;
import org.freehep.graphicsio.gif.GIFGraphics2D; 

/**
   A very simple viewing class to show extracted wavelets
   with log data and seismic on similar scales. Features Zoomable
   enlargements. designed to be used standalone with MLgraphing 
   files produced by the WaveletExtractor. 
   <p>Typical usage:<br>
   <i>(1) Extract wavelet</i>
   <p>
   java au.csiro.JamesGunning.waveletextraction.WaveletExtractor extractionfile.xml ... --dump-ML-parameters
   <p>
   <i>(2) View the maximum likelihood tie</i>
   <p>
   java au.csiro.JamesGunning.waveletextraction.ExtractionViewer [-W] MLgraphing-file.MOST_LIKELY.well-file-name.txt
   <p>

   @author James.Gunning@csiro.au
   @version 1.0, 30/1/2004
*/

public class ExtractionViewer extends Canvas implements MouseListener, MouseMotionListener, KeyListener, ComponentListener, ActionListener, Printable 
{
  // graphics stuff
  public static int DEFAULT_SIZE = 900;	
  public static int GRAPHING_XMARGIN = 15; 
  public static float XPLOT_TO_XRECT_RATIO_MIN = 0.7f; 
  public static float WAVELET_GRAPH_WIDTH_FRACTION = 0.66f;
  public static int WAV_MOV_INC=5;
  private JFrame _jframe;
  // Screen size;
  Dimension Dscreen;
  int sz;                                            // component height
  int sz_width;                                      // component width
  public static int DEFAULT_NUM_GRAPHS=5;            // rho vp vint Strue Sseis
  int yMarginTop=(int)(DEFAULT_SIZE/10.0);           // margin on top of rectangles for logs etc for text
  int yMarginBot=(int)(DEFAULT_SIZE/50.0);           // margin on bottom  of rectangles for logs etc for stuff
  int simpleAxisWidth = (int)(DEFAULT_SIZE/20.0);
  int simpleAxisMargin = (int)(simpleAxisWidth/10.0);
  private int yMin,yMax,graphsYorigin;
  private int xMin;
  private boolean plotVint = false;
  private boolean showGammaAsVsh=false;
  private double vIntFactor = 1.0;
  private boolean plotAI=true;
  private boolean plotDT_blocked=false;
  private boolean plotDTS=true;
  private boolean plotGamma=true;
  private boolean haveFacies=true;
  private boolean expandGraphRangesFully=true;
	
  private boolean blockInterpret = false;
  DoubleArrayList blockInterpList = null;

  // Rectangles for the plots
  GraphInfo gammaRect;
  GraphInfo densityRect;
  GraphInfo sonicRect;
  GraphInfo DTSRect; 
  GraphInfo intVelocityRect;

  Rectangle tvdRect;
  Rectangle mdRect;
  Rectangle timeGridRect;
  GraphInfo[] synthRect;
  GraphInfo[] seisRect;

  // Wavelet stuff
  private boolean showWavelet=false;
  private boolean autoShowWavelet=false;
  String[] wavelet_filenames = null;
  private WaveletData[] waveletData=null;
  private int[] waveletYOffsets;
  private boolean doWaveletDrag;
  private int[] waveletDrag;
  GraphInfo[] waveletRect;

  
  Font axesfont = new Font("Arial", Font.PLAIN, 9);
  Font titlefont = new Font("Arial", Font.PLAIN, 10);
  public static int DEFAULT_N_DEPTH_TICKS = 40;
  public static int DEFAULT_N_TIME_TICKS = 40;


  private double tMax,tMin;
  volatile boolean dumpScreenShot = false;
  private boolean dumpEPSScreenShotAndDie =false;
  private boolean dumpScreenShotAndDie = false;

  // Data for graphing
  int nStacks=1;
  MLgraphingData MLdata;
  private String graphTitle = "Wavelet Extraction";
  private int indxDispMin,indxDispMax;
  private boolean DTnotMicroseconds = false;
  private double DT_scaling=1.0;
  private boolean useZoom=false;
  Petrophysics petro;

  // I/O & error stuff
  ErrorHandler EH;
  String MLgraph_filename=null;
  
  private boolean randomiseDepthsandTimes = false;

  // Stuff for Zoom managing
  private boolean mouseDown=false;
  private Rectangle currentRect;
  private Rectangle rectToDraw = null;
  private Rectangle previousRectDrawn = new Rectangle();
  private AffineTransform at = new AffineTransform();  // current snapshot
  private AffineTransform xAT = new AffineTransform(); // for preconcats...
  private BasicStroke stroke;
  private double currentScale=1.0;

  // mouse 3
  private Rectangle mouseTextRect=null;
  private String mouseText = null;
  private TextLayout mouse3layout=null;
  public static int m3boxMargin=3;

  // Menus
  JMenu menu;

  // scale interpolations..
  LinearInterpolation TinterpTVD;
  LinearInterpolation TinterpMD;
  LinearInterpolation TVDinterpT;
  LinearInterpolation MDinterpT;

  // text formatting
  String dFormatStr = "####.##";    // fixed
  String dFormatStrE = "##0.###E0"; // engineering


  java.text.DecimalFormat dFormat, dFormatE, dSeisTimeFormat;


  public ExtractionViewer()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();

	// Default sizes
	// Get Screen size;
	Dscreen = Toolkit.getDefaultToolkit().getScreenSize();
	sz = min(DEFAULT_SIZE,Dscreen.height-50);         // component height
	sz_width = min(DEFAULT_SIZE,Dscreen.width);    // component width
	// Add mouse handlers
	addMouseMotionListener(this);
	addMouseListener(this);
	// etc
	addComponentListener(this);
	addKeyListener(this);
	// setup Zoom transforms
	at.setToIdentity();
  }


  /** 
	  Example  test driver to show how to view extraction info from typical
	  MLgraphing-file.*.txt files produced by the WabeletExtractor.
	  The display could also be driven from some other class that has a
	  fully populated MLdata object, and then calling the display() method.
  */
  public static void main(String argv[]) 
  {
	ExtractionViewer EV = new ExtractionViewer();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("ExtractionViewer",EV.EH);

   	options.addMandatoryArg(1,"MLgraphing-file.txt","MLgraphing file generated by WaveletExtractor");	
	options.addOption("-c",null,0,null,null,"Add random offset to times and depths for concealment purposes");
	options.addOption("-z",null,0,null,null,"Start up with zoom applied");
	options.addOption("-s",null,0,null,null,"Dump screen shot and and exit");
	options.addOption("-ps",null,0,null,null,"Dump eps screen shot and and exit");
	options.addOption("-w",null,1,"String","wavelet_filename.su","Show supplied wavelet between stacks");
	options.addOption("-W",null,0,null,null,"Autoshow wavelet(s) referenced in MLgraphing file");
	options.addOption("-Vsh",null,0,null,null,"Interpret gamma-log as Vshale");
	options.addSpielText("Use mouse-1 drag to zoom, mouse-2/3 to show values, shift-mouse-1 to "
						 +"manually add interpreted layer-tops, cntrl-mouse-1 to remove interpreted layer-tops");


	// parse options
	options.parseOptions(argv,true);
	if (options.isOptionSet("-c")) 
	  {
		EV.randomiseDepthsandTimes = true;
	  }
	if (options.isOptionSet("-s")) 
	  {
		EV.dumpScreenShotAndDie = true;
	  }
	if (options.isOptionSet("-ps")) 
	  {
		EV.dumpEPSScreenShotAndDie = true;
	  }
	if (options.isOptionSet("-z")) 
	  {
		EV.useZoom = true;
	  }
	if (options.isOptionSet("-w")) 
	  {
		EV.showWavelet=true;
		EV.wavelet_filenames = new String[1];
		EV.wavelet_filenames[0] = (String)options.getOptionValues("-w")[0];
	  }
	if (options.isOptionSet("-W")) 
	  {
		EV.autoShowWavelet=true;
	  }
	if (options.isOptionSet("-Vsh")) 
	  {
		EV.showGammaAsVsh=true;
	  }

	EV.MLgraph_filename = options.getMandatoryArg(1);

	// parse data
	EV.loadASCIIdata();

	EV.loadWavelets();

	// set display range
	EV.setDisplayWindowRange();

	EV.setGraphTitle(EV.MLgraph_filename);

	if (EV.dumpScreenShotAndDie || EV.dumpEPSScreenShotAndDie) 	EV.dumpScreenShot=true;

	// view everything. This can be called as soon as MLdata is fully populated
	EV.display();

	
	if (EV.dumpScreenShotAndDie || EV.dumpEPSScreenShotAndDie) // once only dump and exit
	  {
		EV.dumpScreenShot=true;
		//EV.update((Graphics2D)EV.getGraphics()); 
		EV.grabAndDump(true);
	  }
  }



  // Ascii loading stuff to load from specified file MLgraph_filename
  private void loadASCIIdata()
  {
	int i;
	DoubleArrayList[] SsynthList,StrueList;

	// get a general ASCII reader
	ASCIIfileReader ASCReader = new ASCIIfileReader(EH,false,false);

	// get nStacks & wavelet filenames if available
	nStacks = extractNStacksPlusExtraInfoFromMLGraphingFile(MLgraph_filename);

	// EH.println("Expecting "+nStacks+" stacks");
	if (nStacks<1) EH.die_err("Incorrect Nstacks from "+ 
							  MLgraph_filename +": "+nStacks);

	// form list of required things
	// data buckets
	DoubleArrayList tList = new DoubleArrayList();
	DoubleArrayList tMinicubeList = new DoubleArrayList();
	DoubleArrayList tSeismicList = new DoubleArrayList();
	DoubleArrayList TVDList = new DoubleArrayList();
	DoubleArrayList MDList = new DoubleArrayList();
	DoubleArrayList rhoList = new DoubleArrayList();
	DoubleArrayList DTList = new DoubleArrayList();
	DoubleArrayList DTSList = new DoubleArrayList();
	DoubleArrayList VintList = new DoubleArrayList();
	DoubleArrayList DTblockedList = new DoubleArrayList();
	DoubleArrayList gammaList = new DoubleArrayList();
	DoubleArrayList faciesList = new DoubleArrayList();
	SsynthList = new DoubleArrayList[nStacks];
	StrueList = new DoubleArrayList[nStacks];
	for (i=0; i<nStacks; i++)
	  {
		SsynthList[i] = new DoubleArrayList();
		StrueList[i] = new DoubleArrayList();
	  }
	ASCReader.addToContainerList("t",tList);
	ASCReader.addToContainerList("t-minicube",tMinicubeList,false); // optional: absent on old formats
	ASCReader.addToContainerList("t_seismic",tSeismicList,false);   // optional different name to above, absent on old formats
	ASCReader.addToContainerList("TVD",TVDList);
	ASCReader.addToContainerList("MD",MDList);
	ASCReader.addToContainerList("rho",rhoList);
	ASCReader.addToContainerList("DT",DTList);
	ASCReader.addToContainerList("V-interval",VintList);
	ASCReader.addToContainerList("DT_blocked",DTblockedList,false); // optional: absent on old formats
	ASCReader.addToContainerList("DTS",DTSList,false);              // shear optional: absent on old formats
	ASCReader.addToContainerList("GR",gammaList,false); // optional: may not be there
	ASCReader.addToContainerList("facies_indicator",faciesList,false); // optional: may not be there
	ASCReader.addToContainerList("S_synth",SsynthList[0]);
	ASCReader.addToContainerList("S_true",StrueList[0]);
	if (nStacks>1)
	  {
		for (i=2; i<=nStacks; i++)
		  {
			ASCReader.addToContainerList("S_synth_"+i,SsynthList[i-1]);
			ASCReader.addToContainerList("S_true_"+i,StrueList[i-1]);
		  }
	  }
	// now get the goodies
	ASCReader.extractData(MLgraph_filename);

	// Populate relevent fields of MLdata class
	MLdata = new MLgraphingData();

	MLdata.t = tList.elements();

	// seismic time. Might be missing on old files. 
	// Might also have names "t_seismic" or "t-minicube". Look after all
	// these possibilities
	MLdata.t_minicube = tMinicubeList.elements();
	if (tMinicubeList.size()==0) MLdata.t_minicube = null; // old data formats
	if (MLdata.t_minicube == null)
	  {
		// might be labelled t_seismic - check this
		MLdata.t_minicube = tSeismicList.elements();
		if (tSeismicList.size()==0) MLdata.t_minicube = null; // old data formats
	  }


	MLdata.TVD = TVDList.elements();
	MLdata.MD = MDList.elements();
	MLdata.rho = rhoList.elements();
	MLdata.DT = DTList.elements();
	MLdata.DTS = DTSList.elements();
	if (DTSList.size()==0) 
	  {
		MLdata.DTS = null; // not present
		plotDTS = false;
	  }
	MLdata.Vint = VintList.elements();
	MLdata.DT_blocked = DTblockedList.elements();
	if (DTblockedList.size()==0) MLdata.DT_blocked=null; // old data formats
	MLdata.gamma = gammaList.elements();
	if (gammaList.size()==0) 
	  {
		MLdata.gamma = null; // not present
		plotGamma = false;
	  }
	MLdata.facies_indicator = faciesList.elements();
	if (faciesList.size()==0) 
	  {
		MLdata.facies_indicator = null; // not present
		haveFacies = false;
	  } else haveFacies=true;

	MLdata.Ssynth = new double[nStacks][];
	MLdata.Strue = new double[nStacks][];
	for (i=0; i<nStacks; i++)
	  {
		MLdata.Ssynth[i] = SsynthList[i].elements();
		MLdata.Strue[i] = StrueList[i].elements();
	  }
	// randomisation
	if (randomiseDepthsandTimes)
	  {
		double tOffset = 1+Math.random();
		double zOffset = 1000*(1+Math.random());
		for (i=0;i<MLdata.t.length;i++)
		  {
			MLdata.t[i]+=tOffset;
			if (MLdata.t_minicube !=null) MLdata.t_minicube[i] += tOffset;
			MLdata.TVD[i]+=zOffset;
			MLdata.MD[i]+=zOffset;
		  }
	  }
	// determination of useful data range (from density)
	double[] rho = MLdata.rho;
	MLdata.usefulDataStartIndx=0;
	MLdata.usefulDataEndIndx = rho.length-1;
	for (i=0;i<rho.length;i++)
	  {
		if (rho[i]!=0.0)
		  {
			// add small margin for visual pleasure :-)
			MLdata.usefulDataStartIndx = Math.max(0,i-5); 
			break;
		  }
	  }
	for (i=rho.length-1; i>MLdata.usefulDataStartIndx; i--)
	  {
		if (rho[i]!=0.0)
		  {
			// add small margin for visual pleasure :-)
			MLdata.usefulDataEndIndx = Math.min(rho.length-1,i+5); 
			break;
		  }
	  }
  }

  private void loadWavelets()
  {
	int i,j,k;
	FileInputStream in;
	SuReader sureader = new SuReader();
	SuSegyFormat trW = new SuSegyFormat();


	if (!(showWavelet || autoShowWavelet)) return;

	waveletData = new WaveletData[wavelet_filenames.length];


	for (k=0;k<wavelet_filenames.length;k++)
	  {
		waveletData[k] = new WaveletData();
		// open the wavelet file.
		try
		  {
			in = new FileInputStream(wavelet_filenames[k]);
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + wavelet_filenames[k]);
		  }
		if (sureader.fgettr(in,trW)<1) 
		  throw new RuntimeException("No wavelet data in " + wavelet_filenames[k]);
		// simple tests on critical fields
		if ((waveletData[k].ns = (int) trW.get_short("ns"))<2)
		  throw new RuntimeException("Wavelet is unacceptably short " + wavelet_filenames[k]);
		if (trW.get_short("dt")==0)
		  throw new RuntimeException("Wavelet has weird dt=0 in " + wavelet_filenames[k]);

		waveletData[k].delta_t = 1e-6*trW.get_short("dt"); // in secs

		if (trW.get_short("delrt") != 0)
		  {
			// looks like delrt field defines window start (ms)
			waveletData[k].f1 = (double)trW.get_short("delrt")*1e-3; // to secs
		  }
		else
		  {
			waveletData[k].f1 = (double) trW.get_float("f1");
		  }
		waveletData[k].values = new double[waveletData[k].ns];
	
		// Get wavelet samples.
		// In trace, wavelet is in ordinary form.
		double _values[] = waveletData[k].values;
		for (i=0; i<waveletData[k].ns; i++)
		  {
			_values[i]  = (double) trW.get_data_value(i);
		  }
		// Thin down to eliminate zero paddings
		int ns_new = waveletData[k].ns;
		// coda
		for (i=waveletData[k].ns-1; i>=0 ; i--)
		  {
			if (_values[i] == 0.0)
			  {
				ns_new--;
			  }
			else break;
		  }
		// precursor
		for (i=0; i<waveletData[k].ns; i++)
		  {
			if (_values[i] == 0.0)
			  {
				waveletData[k].f1 += waveletData[k].delta_t;
				ns_new--;
			  }
			else break;
		  }
		if ((i!=0) && (i!= (waveletData[k].ns-1)))
		  {
			double[] new_values = new double[ns_new];
			System.arraycopy(_values,i,new_values,0,ns_new); 
			waveletData[k].values = new_values;
			waveletData[k].ns = ns_new;
		  }
		// .....end of padding elimination

		waveletData[k].tMax = waveletData[k].f1 +(waveletData[k].ns-1)*waveletData[k].delta_t;
	
		waveletYOffsets = new int[nStacks];
		waveletDrag = new int[nStacks];
		Arrays.fill(waveletDrag,-1);


		// close file
		try
		  {
			in.close();
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close " + wavelet_filenames[k]);
		  }
	  }
	// GC
	trW=null;
	sureader=null;
  }

  private void setDisplayWindowRange()
  {
	if (useZoom)
	  {
		indxDispMin = MLdata.usefulDataStartIndx;
		indxDispMax = MLdata.usefulDataEndIndx;
	  }
	else
	  {
		indxDispMin = 0;
		indxDispMax = MLdata.t.length-1;
	  }
	// and set graph mins/maxes
	tMin = MLdata.t[indxDispMin];
	tMax = MLdata.t[indxDispMax];
  }

  /**
	 Extract number of stacks by careful examination of file.
	 (Counts no occurrences of "S_synth*" in 2nd line of filename)
	 Also constructs petrophysics info, plus wavelet filenames, plus establishes DT scaling
  */
  private int extractNStacksPlusExtraInfoFromMLGraphingFile(String filename)
  {
	String[] wC;
	String s;
	int i,cols=0;
	int nStacks=0; 
	try 
	  {
        // f = new RandomAccessFile(filename, "r");
		BufferedReader in = new BufferedReader(new FileReader(filename));		
		int lineNo=0;
		while ((s=in.readLine()) != null)
		  {
			// strip leading & trailing white space
			s = s.trim();
			lineNo++;

			// lineNo==1 is title: look for title keyword
			if (lineNo==1)
			  {
				// (A) look for petrophysics string - strip leading "....Petrophysics:[" then trailing "[...."
				// Petrophysics:[version=1:true,true,70.0,0.0,true,2.2,0.0]
				if (s.contains("Petrophysics:")) 
				  {
					String petrophysics_code = s.replaceFirst(".*Petrophysics:\\[", "").replaceFirst("\\].*","");
					petro = new Petrophysics(petrophysics_code,EH);
				  }

				// (B) look for DT unit scale string - strip leading "....Units:[" then trailing "[...."
				// Units:[rho,dt,dts]
				if (s.contains("Units:")) 
				  {
					String units_code = s.replaceFirst(".*Units:\\[", "").replaceFirst("\\].*","");
					String[] units = units_code.split("\\s+"); // unitsDT unitsDTS unitsRho
					if (units.length>0)
					  {
						if (units[0].toLowerCase().equals("s/m")) 
						  {
							// EH.println("DT in s/m!");
							DTnotMicroseconds = true;
						  }
					  }
				  }
				DT_scaling = DTnotMicroseconds ? 1.0 : FundamentalConstants.SONIC_LOG_TIME_SCALING;

				// (C) look for Wavelets string - strip leading "....Wavelets:[" then trailing "[...."
				// Wavelets:[filename1,filename2..]
				if (!showWavelet)
				  {
					// wavelet filename not user-supplied					
					if (s.contains("Wavelets:")) 
					  {
						String wavelets = s.replaceFirst(".*Wavelets:\\[", "").replaceFirst("\\].*","");
						wavelet_filenames = wavelets.split(",");
					  }
					else
					  {
						// Not available in MLgraphing file. No idea what to do.
						if (autoShowWavelet)
						  {
							EH.println("No wavelets named in "+filename+". No wavelets will be shown. Try -w flag?");			 
							wavelet_filenames = null;
							autoShowWavelet = false;
						  }
					  }
				  }
			  }


			if (lineNo==2)
			  {
				wC = s.split("\\s+");
				cols = wC.length;
				for (i=0;i<cols;i++)
				  {
					if (wC[i].startsWith("S_synth")) nStacks++;
				  }
				wC=null;
				// Should be either 1 wavelet or nStacks wavelets
				if (autoShowWavelet && !(wavelet_filenames.length==1 || wavelet_filenames.length==nStacks))
				  {
					EH.die_err("Number of wavelets (see line 1) not either 1 or the number of stacks, in "+filename);
				  }
			  }
			if (lineNo>2)
			  {
				break;
			  }
		  }	
		in.close();
	  }
	catch (Exception e) 
	  {
		EH.die_err("Can't open file "+filename);
	  }

	return nStacks;
  }


  void display()
  {
	JFrame frame;
	JMenuBar menuBar; // FM
	JMenuItem menuItem;

	precomputeUsefulQuantities();

	frame = new JFrame(MLgraph_filename);
	_jframe = frame;
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	Dimension D = new Dimension (sz_width,sz);
	setSize(D);
	setBackground(Color.lightGray);
	setForeground(Color.black);
	frame.getContentPane().add(this);

	//Create the menu bar.
	menuBar = new JMenuBar();

	// Build the first menu.
	// menus are lightweight and canvas is heavyweight, so we need to toggle this...
	// (Not a good mix in general, but OK for this)
	JPopupMenu.setDefaultLightWeightPopupEnabled(false);
	menu = new JMenu("Actions");

	// Add action Menu 
	menuItem = new JMenuItem("toggle sonic Vp/DT (also Vs/DTS)",KeyEvent.VK_V);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('v'));
	menuItem.addActionListener(this);
	menu.add(menuItem);
	menuItem = new JMenuItem("toggle density / AI",KeyEvent.VK_A);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('a'));
	menuItem.addActionListener(this);
	menu.add(menuItem);
	menuItem = new JMenuItem("toggle Vint / blocked sonic DT",KeyEvent.VK_B);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('b'));
	menuItem.addActionListener(this);
	menu.add(menuItem);
	menuItem = new JMenuItem("toggle detail scale of logs",KeyEvent.VK_L);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('l'));
	menuItem.addActionListener(this);
	menu.add(menuItem);
	menuItem = new JMenuItem("zoom to window over non-zero seismic",KeyEvent.VK_Z);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('z'));
	menuItem.addActionListener(this);
	menu.add(menuItem);

	menuItem = new JMenuItem("dump screen shot to raster file",KeyEvent.VK_S);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('s'));
	menuItem.addActionListener(this);
	menu.add(menuItem);

	menuItem = new JMenuItem("dump EPS file",KeyEvent.VK_E);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('e'));
	menuItem.addActionListener(this);
	menu.add(menuItem);

	menuItem = new JMenuItem("Clear all manual interpretation lines",KeyEvent.VK_C);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('c'));
	menuItem.addActionListener(this);
	menu.add(menuItem);
	menuItem = new JMenuItem("Print",KeyEvent.VK_P);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('p'));
	menuItem.addActionListener(this);
	menu.add(menuItem);	
	
	menu.add(menuItem);
	menuItem = new JMenuItem("quit",KeyEvent.VK_Q);
	menuItem.setAccelerator(KeyStroke.getKeyStroke('q'));
	menuItem.addActionListener(this);
	menu.add(menuItem);

	menuBar.add(menu);

	frame.setJMenuBar(menuBar);

	// frame.setLocation(0,0) ;
	frame.pack();
	frame.setVisible(true);
  } 

  private void setGraphTitle(String title)
  {
	graphTitle = title;
  }

  private void precomputeUsefulQuantities()
  {
	int i;
	// some scales we need once & once only..
	double[] t = MLdata.t;

	DoubleArrayList list = new DoubleArrayList(t);

	// interpolation objects
	TinterpTVD = new LinearInterpolation(MLdata.TVD,MLdata.t,MLdata.TVD.length,EH);
	TinterpMD = new LinearInterpolation(MLdata.MD,MLdata.t,MLdata.MD.length,EH);

	TVDinterpT = new LinearInterpolation(MLdata.t,MLdata.TVD,MLdata.t.length,EH);
	MDinterpT = new LinearInterpolation(MLdata.t,MLdata.MD,MLdata.t.length,EH);
	  
	// velocity scaling
	int cnt=0,cnt2=0;
	double VintBar=0,VlogBar=0;

	
	for (i=0;i<MLdata.Vint.length;i++)
	  {
		if (MLdata.Vint[i]>0)
		  {
			VintBar += MLdata.Vint[i];
			cnt++;
		  }
		if (MLdata.DT[i]>0)
		  {			
			VlogBar += DT_scaling/MLdata.DT[i];
			cnt2++;
		  }
	  }
	VintBar /= Math.max(cnt,1); // will be either ft/s or m/s
	VlogBar /= Math.max(cnt2,1); // will be in ft/s

	if (VlogBar > 1e3*VintBar)
	  {
		EH.warn("Suspect sonic log data not in uS units & not detected in extraction! Converting!");
		DTnotMicroseconds = true;
		DT_scaling = 1.0;
		VlogBar /=  FundamentalConstants.SONIC_LOG_TIME_SCALING;
	  }
	if (VlogBar > 2.0*VintBar)
	  {
		// VintBar in m/s
		vIntFactor = FundamentalConstants.METRES_PER_FOOT;
	  }
	else
	  {
		vIntFactor = 1.0;
	  }

	// formatting for numbers
	dFormat = new java.text.DecimalFormat(dFormatStr); 
	dFormat.setMaximumFractionDigits(3);

	dFormatE = new java.text.DecimalFormat(dFormatStrE); 

	dSeisTimeFormat = new java.text.DecimalFormat("###.##");
	dSeisTimeFormat.setMaximumFractionDigits(0);

  }

  public void paint(Graphics g) 
  //public void paintComponent(Graphics g) 
  {
	Graphics2D g2 = (Graphics2D) g;
	/*
	g.setColor(Color.lightGray);           // Set background color
	g.fillRect(0, 0, getWidth(), getHeight());  // Fill area with background.
	g.setColor(Color.BLACK);             // Restore color for drawing.
	*/
	// set basic sizes
	setBasicSizes();

	// save context Xform
	AffineTransform saveXform = g2.getTransform();
    // apply new one
	g2.setTransform(at);
	// then draw away....

	
	// here , the main graphs we want to see
	paintGraphs(g2);


	// paint Zoom & mouse 3 in old coords
	g2.setTransform(saveXform);

	// mouse 3 stuff
	paintMouse3(g2);

	// zoom box
	paintZoomBox(g2);

	if (dumpScreenShot) 
	  {
		pushGrab();
	  }
  }

  /**
	 Queue-like flagging of screen dump, so image is only drawn when repaint is completed
  */  
  public synchronized void pushGrab()
  {
	notifyAll();
	dumpScreenShot = false;
  }

  /**
	 Call to wait till screen paint is finished, get screen dump, then exit if exit==true.
  */
  public synchronized void grabAndDump(boolean exit)
  {
	while (dumpScreenShot)
	  {
		try { wait();} catch (InterruptedException e) {}
	  }	
	grab();
	if (exit) System.exit(0);
  }


  private void setBasicSizes()
  {
	yMarginTop=(int)(sz/10.0); // margin on top of rectangles for logs etc for text
	yMarginBot=(int)(sz/50.0);  // margin on bottom  of rectangles for logs etc for stuff
	simpleAxisWidth = (int)(sz_width/20.0);
	simpleAxisMargin = (int)(simpleAxisWidth/10.0);
  }

  private void paintMouse3(Graphics2D g2)
  {
	if (mouseTextRect!=null)
	  {	
		g2.setPaint(Color.white);
		g2.fill(mouseTextRect);
		g2.setPaint(Color.red);
		g2.draw(mouseTextRect);
		// txt on top
		g2.setPaint(Color.black);
		float y = (float)(mouseTextRect.getY()+mouseTextRect.getHeight() - m3boxMargin/2.0);
		mouse3layout.draw(g2,(float)mouseTextRect.getX(),y);
	  }
  }

  private void paintGraphs(Graphics2D g2)
  {
	double suggestedScale;
	int i,graphWidth;
	int numGraphs;

	numGraphs = DEFAULT_NUM_GRAPHS;
	numGraphs += 2*(nStacks-1); // another 2 for each extra stack
	if (plotGamma) numGraphs++;
	if (plotDTS) numGraphs++;

	// pane is 3 scales (t,MD,TVD) plus numGraphs rectangles
	graphWidth = (sz_width - 3*simpleAxisWidth)/numGraphs;

	// global registration  values
	graphsYorigin = 0;
	// bounds 
	yMin = graphsYorigin+yMarginTop;
	yMax = graphsYorigin+sz-yMarginBot;
	xMin = 0;

	// optional gamma
	if (plotGamma)
	  {
		gammaRect = new GraphInfo(simpleAxisWidth,0,graphWidth,sz);
		paintLithology(g2,MLdata.t, gammaRect);
		String s = showGammaAsVsh ? "VSH" : "GR";
		paintVerticalGraph(g2,MLdata.gamma,MLdata.t,gammaRect,true,0,false,false,s);
	  }
	else
	  {
		gammaRect = new GraphInfo(simpleAxisWidth,0,0,sz);
	  }

	  
	// density
	// densityRect = new GraphInfo(simpleAxisWidth,0,graphWidth,sz);
	// paintVerticalGraph(g2,MLdata.rho,MLdata.t,densityRect,true,0,false,false,"rho");
	densityRect = new GraphInfo(gammaRect.boundingBox.x+gammaRect.boundingBox.width,0,graphWidth,sz);
	if (!plotGamma) 
	  {
		// stick lithology over density log if gamma not available
		paintLithology(g2,MLdata.t, densityRect);
	  }
	paintVerticalGraph(g2,MLdata.rho,MLdata.t,densityRect,true,0,false,false,"rho");

	// sonic log
	sonicRect = new GraphInfo(densityRect.boundingBox.x+densityRect.boundingBox.width,0,graphWidth,sz);
	paintVerticalGraph(g2,MLdata.DT,MLdata.t,sonicRect,true,0,false,false,"DT");

	// interval velocity
	intVelocityRect = new GraphInfo(sonicRect.boundingBox.x+sonicRect.boundingBox.width,0,graphWidth,sz);
	paintVerticalGraph(g2,MLdata.Vint,MLdata.t,intVelocityRect,true,0,false,true,"V-interval");

	// optional shear sonic log: DTS
	if (plotDTS)
	  {
		DTSRect = new GraphInfo(intVelocityRect.boundingBox.x+intVelocityRect.boundingBox.width,0, graphWidth,sz);
		paintVerticalGraph(g2,MLdata.DTS,MLdata.t,DTSRect,true,0,false,false,"DTS");
	  }
	else
	  {
		DTSRect = new GraphInfo(intVelocityRect.boundingBox.x+intVelocityRect.boundingBox.width,0, 0,sz);
	  }

	// TVD depth axes
	tvdRect = new Rectangle(DTSRect.boundingBox.x+DTSRect.boundingBox.width,0,simpleAxisWidth,sz);
	paintDepthGrid(g2,tvdRect,MLdata.TVD);

	// MD depth axes
	mdRect = new Rectangle(tvdRect.x+tvdRect.width,0,simpleAxisWidth,sz);
	paintDepthGrid(g2,mdRect,MLdata.MD);

	// Synthetics and true seismic. Paint on same scales
	synthRect = new GraphInfo[nStacks];
	seisRect = new GraphInfo[nStacks];

	synthRect[0] = new GraphInfo(mdRect.x+mdRect.width,0,graphWidth,sz);
	suggestedScale = Math.max(getMaxAbs(MLdata.Ssynth[0]),getMaxAbs(MLdata.Strue[0])); 
	paintVerticalGraph(g2,MLdata.Ssynth[0],MLdata.t,synthRect[0],false,suggestedScale,true,true,"synthetic");

	seisRect[0] = new GraphInfo(synthRect[0].boundingBox.x+synthRect[0].boundingBox.width,0,graphWidth,sz);
	paintVerticalGraph(g2,MLdata.Strue[0],MLdata.t,seisRect[0],false,suggestedScale,true,true,"seismic");


	if (nStacks>1)
	  {
		for (i=1; i<nStacks; i++)
		  {
			synthRect[i] = new GraphInfo(seisRect[i-1].boundingBox.x+seisRect[i-1].boundingBox.width,0,graphWidth,sz);
			suggestedScale = Math.max(getMaxAbs(MLdata.Ssynth[i]),getMaxAbs(MLdata.Strue[i])); 
			paintVerticalGraph(g2,MLdata.Ssynth[i],MLdata.t,synthRect[i],false,suggestedScale,true,true,"synth stack "+(i+1));
			
			seisRect[i] = new GraphInfo(synthRect[i].boundingBox.x+synthRect[i].boundingBox.width,0,graphWidth,sz);
			paintVerticalGraph(g2,MLdata.Strue[i],MLdata.t,seisRect[i],false,suggestedScale,true,true,"seis stack "+(i+1));
		  }
	  }

	if (showWavelet||autoShowWavelet)
	  {
		waveletRect = new GraphInfo[nStacks];
		for (i=0; i<nStacks; i++)
		  {
			int widthW =  (int)(WAVELET_GRAPH_WIDTH_FRACTION*synthRect[i].boundingBox.width);
			int xCentre = synthRect[i].boundingBox.x+synthRect[i].boundingBox.width;

			// Set time scales correctly from wavelet span compared to full time scale
			WaveletData wData = (wavelet_filenames.length>1)? waveletData[i] :  waveletData[0];
			int wBoxHeight = (int)((yMax-yMin) * ((wData.tMax - wData.f1)/(tMax-tMin)));

			waveletRect[i] = new GraphInfo(xCentre - widthW/2, 
										   waveletYOffsets[i] + (yMin + yMax- wBoxHeight)/2,
										   widthW,wBoxHeight);
			paintWaveletGraph(g2,waveletRect[i],wData);
		  }
	  }

	// global time grid
	timeGridRect = new Rectangle(0,0,sz_width,sz);
	// any interpreted events
	paintInterpEvents(g2);
	// time grid
	paintTimeGrid(g2,timeGridRect);
	// title
	paintTitle(g2);
  }

  private void paintTitle(Graphics2D g2)
  {
	FontMetrics metrics = g2.getFontMetrics();
	String word = getReducedFilename();
	
	int x = (int)(sz_width - metrics.stringWidth(word))/2;
	int y = (int)(1.5*metrics.getAscent());
	
	g2.setFont(titlefont);

	g2.drawString(word,x,y);
  }

  private void paintInterpEvents(Graphics2D g2)
  {
	if (blockInterpList != null)
	  {
		double t,Y;
		int i,xAxis,rX;
		g2.setPaint(Color.white);
		xAxis = timeGridRect.x + simpleAxisWidth;
		if (MLdata.t_minicube != null)
		  {
			// move a little left to allow for red time words
			xAxis -= (int)(0.35*simpleAxisWidth);
		  }
		rX = timeGridRect.x + timeGridRect.width;
		for (i=0;i<blockInterpList.size(); i++)
		  {
			t = blockInterpList.get(i);
			Y = yMin + (t - tMin)*(yMax-yMin)/(tMax-tMin);
			if (Y>=yMin && Y <=yMax)
			  {
				g2.draw(new Line2D.Double(xAxis,Y,rX,Y));
			  }
		  }
	  }
  }

  // here, depth must be one of MLdata.MD or TVD
  private void paintDepthGrid(Graphics2D g2,Rectangle boundingBox,double depth[])
  {
	int i,oX,graphsYorigin;
	double t,Y;
	double z,zMin,zMax;

	oX = boundingBox.x;
	int rX = oX + boundingBox.width;

	// suitable depth range
	zMin = depth[indxDispMin];
	zMax=  depth[indxDispMax];

	TickMarks ticks = new TickMarks(EH);
	int nDepthTicks = DEFAULT_N_DEPTH_TICKS;
	double tick_spacing = ticks.getTickSpacing(zMin,zMax,nDepthTicks);
	double range = ticks.getRange();
	int min_tick_count = ticks.getMinTickCount();
	int max_tick_count = ticks.getMaxTickCount();
 
	// vertical line
	int xAxis = oX+simpleAxisWidth - simpleAxisMargin;
	setNewStrokeSize(g2,1.0);
	g2.draw(new Line2D.Double(xAxis,yMin,xAxis,yMax));

	// ticks and time labels
	// label
	g2.setFont(titlefont);
	FontMetrics metrics = g2.getFontMetrics();
	String word="TVD";
	if (depth == MLdata.MD) word="MD";
	g2.drawString(word,(int)((oX+rX-metrics.stringWidth(word))/2.0),yMarginTop/2);

	//	Font font = new Font("Arial", Font.PLAIN, 9);
	double fHeight = metrics.getAscent()/2.0;
	double txtMargin = metrics.getAscent()/5.0;
	g2.setFont(axesfont);

	int tickWidth=5;

	// Tinterp = new LinearInterpolation(depth,MLdata.t,depth.length,EH);
	LinearInterpolation TinterpZ = (depth==MLdata.TVD)? TinterpTVD : TinterpMD;
	for (i=min_tick_count; i<=max_tick_count; i++) 
	  {
		// depth
		z = tick_spacing * i;
		if (z<zMin || z>zMax) continue;
		// find corresponding time - draw tick
		t = TinterpZ.interpolate(z,i>min_tick_count);

		Y = yMin + (t - tMin)*(yMax-yMin)/(tMax-tMin);
		g2.draw(new Line2D.Double(xAxis,Y,xAxis+tickWidth,Y));
		String tickWord = dFormat.format(z); // display depth
		//		metrics.stringWidth(tickWord);
		g2.drawString(tickWord,(int)(xAxis-metrics.stringWidth(tickWord)-txtMargin),(int)(Y+fHeight));
	  }

  }

  private void paintTimeGrid(Graphics2D g2,Rectangle boundingBox)
  {
	int i,oX,graphsYorigin,I;
	double tStart,tEnd,dt,t,Y;
	double[] times = MLdata.t;

	oX = boundingBox.x;
	int rX = oX + boundingBox.width;

	TickMarks ticks = new TickMarks(EH);
	int nTimeTicks = DEFAULT_N_TIME_TICKS;
	double tick_spacing = ticks.getTickSpacing(tMin,tMax,nTimeTicks);
	double range = ticks.getRange();
	int min_tick_count = ticks.getMinTickCount();
	int max_tick_count = ticks.getMaxTickCount();

	// vertical line
	int xAxis = oX+simpleAxisWidth;
	FontMetrics metrics = g2.getFontMetrics();
	if (MLdata.t_minicube != null)
	  {
		// move a little left to allow for red time words
		xAxis -= (int)(0.35*simpleAxisWidth);
	  }
	setNewStrokeSize(g2,1.0);
	g2.setPaint(Color.black);
	g2.draw(new Line2D.Double(xAxis,yMin,xAxis,yMax));

	// ticks and time labels
	g2.setFont(titlefont);
	int fHeight = (int)(metrics.getAscent()/2.0);
	double txtMargin = metrics.getAscent()/5.0;
	String word = "chksht time";
	I = oX+(int)Math.max(0,simpleAxisWidth/2.0-metrics.stringWidth(word)/2.0);
	g2.drawString(word,I, yMarginTop/2);
	word = "(ms)";
	I = oX+(int)Math.max(0,simpleAxisWidth/2-metrics.stringWidth(word)/2.0);
	g2.drawString(word, I,(int)(yMarginTop/2+2.5*fHeight));

	if (MLdata.t_minicube != null)
	  {
		g2.setPaint(Color.red);
		word = "seismic time";
		g2.drawString(word,(oX+Math.max(0,simpleAxisWidth/2-metrics.stringWidth(word)/2)),
					  (int)(yMarginTop/2 + 5*fHeight));
		g2.setPaint(Color.black);
	  }

	g2.setFont(axesfont);
	setNewDashedStrokeSize(g2,1,1);


	for (i=min_tick_count; i<=max_tick_count; i++) 
	  {
		t = tick_spacing * i;
		Y = yMin + (t - tMin)*(yMax-yMin)/(tMax-tMin);
		if (Y<yMin || Y>yMax) continue;
		// coloured time dashes
		g2.setPaint(Color.blue);
		g2.draw(new Line2D.Double(xAxis,Y,rX,Y));
		word = dFormat.format(t*1000); // display in ms
		// black words
		g2.setPaint(Color.black);
		g2.drawString(word,(int)(xAxis-metrics.stringWidth(word)-txtMargin),(int)(Y+fHeight));
	  }


	if (MLdata.t_minicube != null)
	  {
		// seismic time scale in red
		g2.setPaint(Color.red);
		double tOffset = (MLdata.t_minicube[0] - MLdata.t[0]);
		for (i=min_tick_count; i<=max_tick_count; i++) 
		  {
			t = tick_spacing * i;
			Y = yMin + (t - tMin)*(yMax-yMin)/(tMax-tMin);
			if (Y<yMin || Y>yMax) continue;
			word = dSeisTimeFormat.format((t+tOffset)*1000); // display in ms
			// red minicube time 
			g2.drawString(word,(int)(xAxis+txtMargin),(int)(Y+fHeight));
		  }
	  }

	// restore stroke & black color
	g2.setPaint(Color.black);
	setNewStrokeSize(g2,1.0);
  }


  /**
	 Draw Lithology colouring behind GI.
  */
  private void paintLithology(Graphics2D g2, double[] t,GraphInfo GI)
  {
	int i;
	double iY;


	// skip all this if no petrophysics available, 
	// or anisotropy terms not used in model
	if ((petro==null || !petro.hasAnisotropy) && !haveFacies) return;

	double[] GR;
	double[] rho;
	double[] zTVD;

	GR = (plotGamma) ? MLdata.gamma: MLdata.rho; 
	// rho always exists as fallback: if (!plotGamma), GR not used in isShale() below

	rho = MLdata.rho;
	zTVD = MLdata.TVD;

	double widthL = (GI.xMax-GI.xMin);
	double yFact = (yMax-yMin)/(tMax-tMin);

	boolean shale=true,shale_last=true;
	Color brown = new Color(150, 75, 0);

	// first pt
	double iYo = (yMin + (t[indxDispMin]-tMin)*yFact);
	shale_last = (petro!=null) ? petro.isShale(plotGamma, zTVD[indxDispMin], GR[indxDispMin], rho[indxDispMin]) : false;
	for (i=indxDispMin+1; i<= indxDispMax ; i++) 
	  {
		if (true)
		  {
			if (petro!=null)
			  {
				if  (petro.version==1) 
				  {
					shale = (MLdata.facies_indicator[i] == 1.0);
				  }
			  }
			else if (haveFacies)
			  {
				// using facies_indicator for something else perhaps
				shale = (MLdata.facies_indicator[i] == 1.0);
			  }
		  }
		else
		  {
			// version from logs
			shale = petro.isShale(plotGamma, zTVD[i], GR[i], rho[i]);
		  }
		iY = (yMin + (t[i]-tMin)*yFact);
		if (shale != shale_last || (i== indxDispMax)) // at least 1 rectangle...
		  {
			g2.setPaint(shale_last ? brown  : Color.yellow);
			g2.fill(new Rectangle2D.Double(GI.xMin,iYo,widthL, iY - iYo));
			// save old pts
			iYo = iY;
			shale_last = shale;
		  }
	  }

	g2.setPaint(Color.black);
  }



  /**
	 Draw graph of data. Use positive scale only if autoScalePositive==true.
	 Use suggested scale if suggestedScale>0.
	 Use fewer ticks if sparseTicks==true.
  */
  private void paintVerticalGraph(Graphics2D g2, double[] data, double[] t,GraphInfo GI, 
								  boolean autoScalePositive, double suggestedScale,
								  boolean drawZeroFAxis,
								  boolean sparseTicks,
								  String title)
  {
	int i,graphsYorigin;
	double iX,iY,u;
	double[] f,fNew;
	boolean outOfRange = false;

	// bounding rectangle
	setNewStrokeSize(g2,1.0);
	//g2.drawRect(boundingBox.x, graphsYorigin,boundingBox.width,boundingBox.height);
	
	// bounds
	// int xMin = boundingBox.x+GRAPHING_XMARGIN;
	// int xMax = boundingBox.x+boundingBox.width-GRAPHING_XMARGIN;

	// possible data views
	f = data;

	if (data == MLdata.rho)
	  {
		if (plotAI)
		  {
			// show pseudo AI log instead of rho
			double s = vIntFactor*DT_scaling;
			fNew = f.clone();
			for (i=0;i<fNew.length;i++) 
			  {
				if (f[i]>0 && MLdata.rho[i]>0)
				  {
					// Phi = rho/DT
					fNew[i] = MLdata.rho[i]*(s/MLdata.DT[i]);
				  }
				else	fNew[i]=0;
			  }
			f=fNew;
			title="AI";
			autoScalePositive=true;
		  }
		else f=data;
	  }
	else if (data == MLdata.DT)
	  {
		if (!plotVint)
		  {
			// show pseudo DT log
			fNew = f.clone();
			for (i=0;i<fNew.length;i++) 
			  {
				if (f[i]>0)
				  {
					// uS/ft scale for DT
					// FIXME! for s/m
					fNew[i] = vIntFactor*DT_scaling/f[i];
				  }
				else fNew[i]=0;
			  }
			f=fNew;
			title="sonic Vp";
			autoScalePositive=true;
		  }
		else f=data;
	  }
	else if (data == MLdata.DTS)
	  {
		if (!plotVint)
		  {
			// show pseudo DT log
			fNew = f.clone();
			for (i=0;i<fNew.length;i++) 
			  {
				if (f[i]>0)
				  {
					// uS/ft scale for DT
					fNew[i] = vIntFactor*
					  DT_scaling/f[i];
				  }
				else fNew[i]=0;
			  }
			f=fNew;
			title="sonic Vs";
			autoScalePositive=true;
		  }
		else f=data;
	  }
	else if (data == MLdata.Vint)
	  {
		if (plotDT_blocked && MLdata.DT_blocked!=null) // DT_blocked poss. absent on old sets
		  {
			// show blocked sonic log
			fNew = f.clone();
			for (i=0;i<fNew.length;i++) 
			  {
				if (f[i]>0)
				  {
					// uS/ft scale for DT
					fNew[i] = MLdata.DT_blocked[i];
				  }
				else fNew[i]=0;
			  }
			f=fNew;
			title="DT-blocked";
			autoScalePositive=true;
		  }
		else f=data;
	  }
	else
	  {
		f=data;
	  }
	// ptr for object
	GI.f = f;

	// scales
	if (suggestedScale > 0)
	  {
		GI.fMax = suggestedScale;
		GI.fMin=0;
		if (autoScalePositive) GI.fMin=0;
		else GI.fMin = -suggestedScale;
	  }
	else
	  {
		DoubleArrayList list = new DoubleArrayList(f);
		GI.fMax = cern.jet.stat.Descriptive.max(list);


		GI.fMin=0;
		if (!autoScalePositive)
		  {
			GI.fMin = cern.jet.stat.Descriptive.min(list);
		  }
		if (expandGraphRangesFully)
		  {
			GI.fMin = getPositiveMin(f);
		  }
	  }
	// convenience factors
	double xFact, yFact;
	if (GI.fMax>GI.fMin)
	  {
		xFact = (GI.xMax-GI.xMin)/(GI.fMax-GI.fMin);
	  }
	else
	  {
		xFact=0.0;
	  }
	if (tMax>tMin)
	  {
		yFact = (yMax-yMin)/(tMax-tMin);
	  }
	else
	  {
		yFact=0.0;
	  }

	// now the curve
	// adjust linewidths
	setNewStrokeSize(g2,1.0);
	// first pt
	double iXo,iYo;
	iXo = GI.xMin + Math.max(f[indxDispMin]-GI.fMin, 0)*xFact;
	iYo = yMin + (t[indxDispMin]-tMin)*yFact;
	for (i=indxDispMin+1; i<= indxDispMax ; i++) 
	  {
		iX = (f[i]-GI.fMin);
		if (expandGraphRangesFully)
		  {
			if (iX<0)
			  {
				iX=0; // contain -ve excursions if needed
				g2.setPaint(Color.red);
				outOfRange=true;
			  }
			else if (outOfRange)
			  {
				outOfRange = false;
				g2.setPaint(Color.black);
			  }
		  }
		iX = (GI.xMin + iX*xFact);
		iY = (yMin + (t[i]-tMin)*yFact);
		// g2.drawLine(iXo,iYo,iX,iY);
		g2.draw(new Line2D.Double(iXo,iYo,iX,iY));
		// save old pts
		iXo = iX;
		iYo = iY;
	  }
	// restore black to be sure
	if (expandGraphRangesFully) g2.setPaint(Color.black);

	// f-axis
	setNewStrokeSize(g2,2.0);
	//BasicStroke curStroke = ((BasicStroke)g2.getStroke());
	// System.out.println("LW: "+curStroke.getLineWidth()*currentScale);
	

	g2.draw(new Line2D.Double(GI.xMin,yMin,GI.xMax,yMin));



	// f-axis tick marks
	setNewStrokeSize(g2,1.0);
	int tick_size = (int)Math.max(10,GI.boundingBox.height/50.0);

	TickMarks ticks = new TickMarks(EH);
	double tick_spacing = ticks.getTickSpacing(GI.fMin,GI.fMax,sparseTicks ? 3 : 3);
	double range = ticks.getRange();
	int min_tick_count = ticks.getMinTickCount();
	int max_tick_count = ticks.getMaxTickCount(); 

	g2.setFont(axesfont);
	FontMetrics metrics = g2.getFontMetrics();
	if (sparseTicks) dFormat.applyPattern("0.###E0##");

	AffineTransform aT = g2.getTransform();
	AffineTransform aTrot;
	int wordHeight = metrics.getHeight();
	for (i=min_tick_count; i<=max_tick_count; i++) 
	  {
		// f = tick_spacing * i;
		iX = (int)(GI.xMin + (tick_spacing*i - GI.fMin)*xFact);
		if (iX<GI.xMin ||iX>GI.xMax) continue;
		iY = yMin;

		g2.draw(new Line2D.Double(iX,yMin,iX, yMin + tick_size));
		u = tick_spacing*i;
		String tick = sparseTicks ? dFormat.format(u) : dFormatE.format(u);

		// int wordWidth = metrics.stringWidth(tick);
		aTrot = g2.getTransform(); // copy of aT
		g2.setTransform(aTrot);
		// rotate about tick point
		g2.rotate(-Math.PI/2.0,(int)iX,yMarginTop);
		g2.drawString(tick,(int)iX+wordHeight/3,yMarginTop+wordHeight/2);
		// set back to normal
		g2.setTransform(aT);
	  }


	if (sparseTicks) dFormat.applyPattern(dFormatStr);

	// t-axis
	if (drawZeroFAxis)
	  {
		setNewStrokeSize(g2,1.0);
		iX = (GI.xMin -GI.fMin*xFact);
		g2.draw(new Line2D.Double(iX,yMin,iX,yMax));
	  }

	// Write a title
	// some words
	if (title!=null)
	  {
 		g2.setFont(titlefont);
		metrics = g2.getFontMetrics();
		g2.drawString(title,(GI.xMin+(GI.xMax-GI.xMin)/2-metrics.stringWidth(title)/2),yMarginTop/2);
	  }
  }


  /**
	 Paint on wavelet
  */

  void  paintWaveletGraph(Graphics2D g2,GraphInfo GI, WaveletData wData)
  {
	int i;
	double iX,iY;
	double[] f;
		
	// bounds
	// int xMin = boundingBox.x+GRAPHING_XMARGIN;
	// int xMax = boundingBox.x+boundingBox.width-GRAPHING_XMARGIN;

	// possible data views
	f = wData.values;

	// ptr for object
	GI.f = f;

	DoubleArrayList list = new DoubleArrayList(f);
	GI.fMax = cern.jet.stat.Descriptive.max(list);
	GI.fMin = cern.jet.stat.Descriptive.min(list);
	// Symmetric scale for wavelet
	GI.fMax = Math.max(GI.fMax,-GI.fMin);
	GI.fMin = -GI.fMax;

	// Set time scales correctly from wavelet span compared to full time scale
	double yWavMin = GI.boundingBox.y;
	double yWavMax = yWavMin + (yMax-yMin) * ((wData.tMax - wData.f1)/(tMax-tMin));
	double dy = (wData.ns>1)? (yWavMax - yWavMin)/(wData.ns-1.0): 0;


	// convenience factors
	double xFact=0;
	if (GI.fMax>GI.fMin)
	  {
		xFact = (GI.xMax-GI.xMin)/(GI.fMax-GI.fMin);
	  }
	// now the curve
	
	// red
    g2.setPaint(Color.red);
	// adjust linewidths
	setNewStrokeSize(g2,1.0);
	// first pt
	double iXo = (GI.xMin + (f[0]-GI.fMin)*xFact);
	double iYo = yWavMin;
	for (i=0; i<wData.ns; i++) 
	  {
		iX = (GI.xMin + (f[i]-GI.fMin)*xFact);
		iY = yWavMin + i*dy;
		if (iY>=yMin && iY<=yMax)
		  {
			g2.draw(new Line2D.Double(iXo,iYo,iX,iY));
		  }
		// save old pts
		iXo = iX;
		iYo = iY;		
	  }

	// black axes
	g2.setPaint(Color.black);
	// t-axis
	iX = (GI.xMin -GI.fMin*xFact);
	g2.draw(new Line2D.Double(iX,Math.max(yMin, yWavMin),
							  iX,Math.min(yMax,yWavMax)));
	// 0-time axis
	iY = yWavMin - (wData.f1/wData.delta_t)*dy;
	g2.draw(new Line2D.Double(GI.xMin,iY, GI.xMax,iY));

  }

  /**
	 Joint movement of all wavelets called by keybinding UP/DOWN arrows 
  */
  void moveAllWavelets(boolean up)
  {
	if (!(showWavelet || autoShowWavelet)) return;

	int i;
	int dy = up ? -WAV_MOV_INC : WAV_MOV_INC;
	int yNew;

	for (i=0; i<nStacks; i++)
	  {
		yNew = waveletYOffsets[i] + (yMin + yMax - waveletRect[i].boundingBox.height)/2 + dy;
		if (yNew>=yMin && (yNew + waveletRect[i].boundingBox.height <= yMax))
		  {
			waveletYOffsets[i] += dy;
		  }
		if (currentScale==1.0)
		  {
			Rectangle Rold = new Rectangle(waveletRect[i].boundingBox);
			Rectangle Rnew = new Rectangle(waveletRect[i].boundingBox);
			Rnew.translate(0,dy) ;
			Rold = Rnew.union(Rold);
			repaint(Rold.x, Rold.y,	Rold.width+1, Rold.height+1);
		  }
		else repaint();
	  }
  }

  /**
	 Individual wavelet motion
  */
  void moveWavelets(int mY)
  {
	int i;
	int dy;
	int yNew;

	if (!(showWavelet||autoShowWavelet)) return;

	for (i=0; i<nStacks; i++)
	  {
		if (waveletDrag[i] >= 0)
		  {
			double xOrig, yOrig;
			
			// get X,Y in original coords..(robust to current Xform)
			// xOrig = (mX - at.getTranslateX())/at.getScaleX();
			yOrig = (mY - at.getTranslateY())/at.getScaleY();

			dy = (int)Math.rint(yOrig - waveletDrag[i]);
			yNew = waveletYOffsets[i] + (yMin + yMax - waveletRect[i].boundingBox.height)/2 + dy;
			if (yNew>=yMin && (yNew + waveletRect[i].boundingBox.height <= yMax))
			  {
				waveletYOffsets[i] += dy;
			  }

			if (currentScale==1.0)
			  {
				Rectangle Rold = new Rectangle(waveletRect[i].boundingBox);
				Rectangle Rnew = new Rectangle(waveletRect[i].boundingBox);
				Rnew.translate(0,dy) ;
				Rold = Rnew.union(Rold);
				repaint(Rold.x, Rold.y,	Rold.width+1, Rold.height+1);
			  }
			else
			  {
				repaint();
			  }
		  }
	  }
  }

  private void checkForWaveletMove(int mX,int mY)
  {

	if (!(showWavelet||autoShowWavelet)) return;

	int i;
	double xOrig, yOrig;

	// get X,Y in original coords..(robust to current Xform)
	xOrig = (mX - at.getTranslateX())/at.getScaleX();
	yOrig = (mY - at.getTranslateY())/at.getScaleY();

	if (!doWaveletDrag)
	  {
		// start of new drag...
		for (i=0; i<nStacks; i++)
		  {
			waveletDrag[i] = -1;
			if (waveletRect[i].boundingBox.contains(xOrig,yOrig))
			  {
				waveletDrag[i] = (int)Math.rint(yOrig);
				doWaveletDrag = true;
			  }
		  }
	  }
	else
	  {
		// already dragging: keep track
		for (i=0; i<nStacks; i++)
		  {
			if (waveletDrag[i] >=0) waveletDrag[i] = (int)Math.rint(yOrig);
		  }
	  }
  }

  // Shows log values in region pointed to
  private void showLogValues(int mX,int mY)
  {
	int i;
	double t,tOffset, xOrig, yOrig,zTVD,zMD,val=0,x;
	boolean showTseis;

	GraphInfo curGI=null;
	// get X,Y in original coords..(robust to current Xform)
	xOrig = (mX - at.getTranslateX())/at.getScaleX();
	x = (xOrig-xMin);
	yOrig = (mY - at.getTranslateY())/at.getScaleY();
	t = (yOrig-yMin)*(tMax-tMin)/(yMax-yMin) + tMin;

	if (t>=tMin && t<=tMax)
	  {
		zTVD = TVDinterpT.interpolate(t,false);
		zMD = MDinterpT.interpolate(t,false);
		// test which we are in and print out nearest value
		int indx = (int)Math.rint(MLdata.t.length*(t-tMin)/(tMax-tMin));

		if (gammaRect.boundingBox.contains(xOrig,yOrig))
		  {
			curGI = gammaRect;
		  }
		else if (densityRect.boundingBox.contains(xOrig,yOrig))
		  {
			curGI = densityRect;
		  }
		else if (sonicRect.boundingBox.contains(xOrig,yOrig))
		  {
			curGI = sonicRect;
		  }
		else if (intVelocityRect.boundingBox.contains(xOrig,yOrig))
		  {
			curGI = intVelocityRect;
		  }
		else if (DTSRect.boundingBox.contains(xOrig,yOrig))
		  {
			curGI = DTSRect;
		  }
		else 
		  {
			for (i=0; i< nStacks; i++)
			  {
				if (synthRect[i].boundingBox.contains(xOrig,yOrig))
				  {
					curGI = synthRect[i];
				  }
				else if (seisRect[i].boundingBox.contains(xOrig,yOrig))
				  {
					curGI = seisRect[i];
				  }
			  }
		  }

		val = (curGI!=null)? (curGI.fMax-curGI.fMin)*(xOrig-curGI.xMin)/(curGI.xMax-curGI.xMin)+curGI.fMin : Double.MAX_VALUE;

		// and on panel
		Graphics2D g2 = (Graphics2D)getGraphics();
		showTseis = (MLdata.t_minicube != null);
		tOffset = showTseis ? (MLdata.t_minicube[0] - MLdata.t[0]) : 0.0;
		String seisTtext = (showTseis) ? String.format(" t_seis=%-14.6g",(t+tOffset)*1000.0) : "";
		if (val != Double.MAX_VALUE)
		  {
			mouseText = String.format(" t=%-14.6g%s TVD=%-14.6g MD=%-14.6g val=%-14.6g",
									  t*1000.0,seisTtext, zTVD, zMD, val);
		  }
		else
		  {
			mouseText = String.format(" t=%-14.6g%s TVD=%-14.6g MD=%-14.6g", t*1000.0, seisTtext, zTVD, zMD);
		  }
		FontRenderContext frc = g2.getFontRenderContext();
		mouse3layout = new TextLayout(mouseText, titlefont, frc);
		// white background
		Rectangle2D bgBox = mouse3layout.getBounds();
		int height = (int)Math.ceil(bgBox.getHeight()+m3boxMargin);
		int width = (int)Math.ceil(bgBox.getWidth())+m3boxMargin;
		if (mX+width < sz_width)
		  {
			// fits ok
			mouseTextRect = new Rectangle(mX,mY-height,width,height);
		  }
		else
		  {
			// fit left
			mouseTextRect = new Rectangle(mX-width,mY-height,width,height);
		  }
		repaint(mouseTextRect.x, mouseTextRect.y,mouseTextRect.width+1, mouseTextRect.height+1);
	  }
  }



  // Print out t, TVD, MD corresponding to current mouse click
  // accepts current mouse coords
  private void showTimesAndDepths(int mX,int mY)
  {
	double t, yOrig,zTVD,zMD;

	// get Y in original coords..(robust to current Xform)
	yOrig = (mY - at.getTranslateY())/at.getScaleY();
	t = (yOrig-yMin)*(tMax-tMin)/(yMax-yMin) + tMin;

	if (t>=tMin && t<=tMax)
	  {
		zTVD = TVDinterpT.interpolate(t,false);
		zMD = MDinterpT.interpolate(t,false);

		// show coords at stdout
		EH.println("%-14.6g %-14.6g %-14.6g",t*1000.0, zTVD, zMD);

		// and on panel
		Graphics2D g2 = (Graphics2D)getGraphics();
		mouseText = String.format(" t=%-14.6g TVD=%-14.6g MD=%-14.6g ", t*1000.0, zTVD, zMD);
		FontRenderContext frc = g2.getFontRenderContext();
		mouse3layout = new TextLayout(mouseText, titlefont, frc);
		// white background
		Rectangle2D bgBox = mouse3layout.getBounds();
		int height = (int)Math.ceil(bgBox.getHeight()+m3boxMargin);
		int width = (int)Math.ceil(bgBox.getWidth())+m3boxMargin;

		mouseTextRect = new Rectangle(mX,mY-height,width,height);
		if (mX+width < (sz_width))
		  {
			// fits ok
			mouseTextRect = new Rectangle(mX,mY-height,width,height);
		  }
		else
		  {
			// fit left
			mouseTextRect = new Rectangle(mX-width,mY-height,width,height);
		  }

		//		System.out.println("BOX "+mouseTextRect.getWidth()+" "+mouseTextRect.getHeight());
		repaint(mouseTextRect.x, mouseTextRect.y,mouseTextRect.width+1, mouseTextRect.height+1);
	  }
  }

  private void addBlockInterpEvent(int mX,int mY)
  {
	double t, yOrig,zTVD,zMD;
		
	// get Y in original coords..(robust to current Xform)
	yOrig = (mY - at.getTranslateY())/at.getScaleY();
	t = (yOrig-yMin)*(tMax-tMin)/(yMax-yMin) + tMin;

	if (t>=tMin && t<=tMax)
	  {
		if (blockInterpList==null) blockInterpList = new DoubleArrayList();

		blockInterpList.add(t);
		blockInterpList.sort(); // maintain order
		if (blockInterpList.size() > 1)
		  {
			// show averages from last section
			zTVD = TVDinterpT.interpolate(t,false);
			zMD = MDinterpT.interpolate(t,false);
			// print interp-blocked current model
			showInterpBlockedModel();
		  }
	  }

	//		System.out.println("BOX "+mouseTextRect.getWidth()+" "+mouseTextRect.getHeight());
	repaint(); // mouseTextRect.x, mouseTextRect.y,mouseTextRect.width+1, mouseTextRect.height+1);
  }

  private void removeBlockInterpEvent(int mX,int mY)
  {
	int i;
	double t, yOrig,zTVD,zMD;
	double DELETE_T_EPS = 2.0e-3;
	double tC;
	boolean renew=false;

	if (blockInterpList==null) return;

	// get Y in original coords..(robust to current Xform)
	yOrig = (mY - at.getTranslateY())/at.getScaleY();
	t = (yOrig-yMin)*(tMax-tMin)/(yMax-yMin) + tMin;

	if (t>=tMin && t<=tMax)
	  {
		for (i=0;i<blockInterpList.size();i++)
		  {
			tC = blockInterpList.get(i);
			if (Math.abs(t - tC) < DELETE_T_EPS) // in MS
			  {
				blockInterpList.delete(tC);
				blockInterpList.sort();
				renew=true;
				break;
			  }
		  }
	  }

	if (renew) repaint(); // mouseTextRect.x, mouseTextRect.y,mouseTextRect.width+1, mouseTextRect.height+1);
  }

  void showInterpBlockedModel()
  {
	int i,j,jS,jE,count;
	double t,tE,zMD,zTVD,tOffset;
	double vpbar,vsbar,rhobar;

	if (blockInterpList.size() > 1)
	  {
		tOffset = (MLdata.t_minicube != null) ? (MLdata.t_minicube[0] - MLdata.t[0]) : 0.0;

		EH.printf("Blocking means over marked intervals: Tchksht(ms)_top Tseismic(ms)_top TVD_top MD_top  <Vp> %s <rho>%n",
				   (MLdata.DTS != null) ? "<Vs>" : "");
		for (i=0;i<blockInterpList.size()-1; i++)
		  {
			t = blockInterpList.get(i);
			zTVD = TVDinterpT.interpolate(t,false);
			zMD = MDinterpT.interpolate(t,false);
			// gets indices over interval
			jS = (int)Math.rint(MLdata.t.length*(t-tMin)/(tMax-tMin));
			tE = blockInterpList.get(i+1);
			jE = (int)Math.rint(MLdata.t.length*(tE-tMin)/(tMax-tMin)) - 1;
			vpbar = 0.0;
			vsbar = 0.0;
			rhobar = 0.0;
			for (j=jS;j<=jE;j++)
			  {
				vpbar += 1.0/MLdata.DT[j];
				rhobar += MLdata.rho[j];
			  }
			if (MLdata.DTS != null)
			  {
				for (j=jS;j<=jE;j++) vsbar += 1.0/MLdata.DTS[j];
			  }
			// fix sonic units
			vpbar *= vIntFactor*DT_scaling;
			vsbar *= vIntFactor*DT_scaling;
			// average
			count = jE-jS+1;
			if (count > 0)
			  {
				vpbar /= count;
				vsbar /= count;
				rhobar /= count;
			  }
			EH.printf("%-14.6g %-14.6g %-14.6g %-14.6g %-14.6g ",t*1000.0, (t+tOffset)*1000.0, zTVD, zMD,vpbar);
			if (MLdata.DTS != null) EH.printf("%-14.6g ",vsbar);  // shear if present
			EH.printf("%-14.6g%n",rhobar);
		  }
	  }
  }

  // get suitable scale for a data set
  private double getMaxAbs(double[] data)
  {
	DoubleArrayList list = new DoubleArrayList(data);
	double min_entry = cern.jet.stat.Descriptive.min(list); 

	return Math.max(cern.jet.stat.Descriptive.max(list),
					-cern.jet.stat.Descriptive.min(list));
  }
  
  // get suitable minima for an expected-positive data set
  private double getPositiveMin(double[] data)
  {
	int i;
	double min=0.0;

	for (i=0;i<data.length;i++)
	  {
		if (data[i]>0) 
		  {
			min = data[i];
			break;
		  }
	  } 
	for (;i<data.length;i++)
	  {
		if (data[i]>0) min = Math.min(min,data[i]);
	  }
	return min;
  }

  // sets appropriate stroke size for current scale
  private void setNewStrokeSize(Graphics2D g2,double size)
  {
	g2.setStroke(new BasicStroke((float)(size/currentScale)));
  }

  /**
	 sets appropriate stroke size for current scale
	 type = 0; equispaced dashes
	 type = 1: light 1:4 ratio dashes
  */
  private void setNewDashedStrokeSize(Graphics2D g2,double width, int type)
  {
	BasicStroke st=null;


	switch (type)
	  {
	  case 0:
		float[] dash = {10.0f};
		st = new BasicStroke((float)(width/currentScale),
							 BasicStroke.CAP_SQUARE, 
							 BasicStroke.JOIN_MITER,
							 10.0f, dash, 0);
		break;
	  case 1:
		float[] sdash = {1.0f,4.0f};
		st = new BasicStroke((float)(width/currentScale),
							 BasicStroke.CAP_SQUARE, 
							 BasicStroke.JOIN_MITER,
							 10.0f, sdash, 0);
		break;
	  default:
		st = new BasicStroke((float)(width/currentScale));
	  }
	g2.setStroke(st);
  }


  
  public void mousePressed(MouseEvent e) 
  {
	int modifiers = e.getModifiers();
	if ((modifiers & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK)
	  {
		if (e.isShiftDown())
		  {
			addBlockInterpEvent(e.getX(), e.getY());
		  }
		else if (e.isControlDown())
		  {
			removeBlockInterpEvent(e.getX(), e.getY());
		  }
		else
		  {
			doWaveletDrag=false;
			checkForWaveletMove(e.getX(), e.getY());
			if (!doWaveletDrag)
			  {
				// mouse 1 zoom
				currentRect = new Rectangle(e.getX(), e.getY(), 0, 0);
				updateDrawableRect(getWidth(), getHeight());
			  }
		  }
	  }
	if ((modifiers & InputEvent.BUTTON2_MASK) == InputEvent.BUTTON2_MASK)
	  {
		// mouse 2: show log values
		showLogValues(e.getX(), e.getY());
	  }
	if ((modifiers & InputEvent.BUTTON3_MASK) == InputEvent.BUTTON3_MASK)
	  {
		// mouse 3: show current t, TVD, MD
		showTimesAndDepths(e.getX(), e.getY());
	  }
  }
  
  public void mouseDragged(MouseEvent e) 
  {
	if (((e.getModifiers() & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK)
		&& !e.isShiftDown() && !e.isControlDown())    // not shift or cntrl
	  { 
		if (doWaveletDrag)
		  {
			moveWavelets(e.getY());
			checkForWaveletMove(e.getX(), e.getY());
		  }
		else if (!blockInterpret)
		  {	// mouse 1 zoom
			// update canvas scale
			updateSize(e);
			mouseDown=true;
		  }
	  }
  }
    
  public void mouseReleased(MouseEvent e) 
  {
	int mods;
	mods = e.getModifiers();
	if (((mods & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK) // button 1
		&& !e.isShiftDown()       // not shift
		&& !e.isControlDown())    // or cntrl
	  {
		if (doWaveletDrag)
		  {
			moveWavelets(e.getY());
			doWaveletDrag=false;
		  }
		else 
		  {
			// mouse 1 zoom
			updateSize(e);
			if (mouseDown)
			  {
				// end of drag. Do A ZoomBox
				// translation defined by left edge of box
				if (rectToDraw!=null)
				  {
					xAT.setToTranslation(-rectToDraw.x,-rectToDraw.y);
					at.preConcatenate(xAT);
		
					// and scale by box ratios
					double Sx = (double)getWidth()/Math.max((double)rectToDraw.width,1);
					double Sy = (double)getHeight()/Math.max((double)rectToDraw.height,1);
		
					xAT.setToScale(Sx,Sy);
					at.preConcatenate(xAT);
					currentScale *= Math.sqrt(Sx*Sx+Sy*Sy);
				  }
			  }
			else
			  {
				// restore original
				at.setToIdentity();
				currentScale=1.0;
			  }
			mouseDown=false;

			// to clear box...
			currentRect=null;
			repaint();
		  }
	  }
	if ((mods & InputEvent.BUTTON2_MASK) == InputEvent.BUTTON2_MASK)
	  {
		// mouse 2 depth box
		if (mouseTextRect!=null) repaint(mouseTextRect.x, mouseTextRect.y,	mouseTextRect.width+1, mouseTextRect.height+1);
		mouseTextRect = null;
	  }
	else if ((mods & InputEvent.BUTTON3_MASK) == InputEvent.BUTTON3_MASK)
	  {
		// mouse 3 depth box
		if (mouseTextRect!=null) repaint(mouseTextRect.x, mouseTextRect.y,	mouseTextRect.width+1, mouseTextRect.height+1);
		mouseTextRect = null;
	  }
  }
  
  /* 
   * Update the size of the current rectangle
   * and call repaint.  Because currentRect
   * always has the same origin, translate it
   * if the width or height is negative.
   * 
   * For efficiency (though
   * that isn't an issue for this program),
   * specify the painting region using arguments
   * to the repaint() call.
   * 
   */
  private void updateSize(MouseEvent e) 
  {
	if (currentRect!=null)
	  {
		currentRect.setSize(e.getX() - currentRect.x, e.getY() - currentRect.y);
		updateDrawableRect(getWidth(), getHeight());
		Rectangle totalRepaint = rectToDraw.union(previousRectDrawn);
		repaint(totalRepaint.x, totalRepaint.y,	totalRepaint.width, totalRepaint.height);
	  }
  }

  private void paintZoomBox(Graphics g2) 
  {
	//If currentRect exists, paint a box on top.
	if (currentRect != null) 
	  {
		//Draw a rectangle on top of the image.
		// g.setXORMode(Color.white); //Color of line varies
		//depending on image colors
		g2.drawRect(rectToDraw.x, rectToDraw.y, 
				   rectToDraw.width - 1, rectToDraw.height - 1);
		
		updateLabel(rectToDraw);
	  }
  }

 
  private void updateLabel(Rectangle rect) {
	int width = rect.width;
	int height = rect.height;
	
	//Make the coordinates look OK if a dimension is 0.
	if (width == 0) {
	  width = 1;
	}
	if (height == 0) {
	  height = 1;
	}
	
	/*
	  label.setText("Rectangle goes from ("
	  + rect.x + ", " + rect.y + ") to ("
	  + (rect.x + width - 1) + ", "
	  + (rect.y + height - 1) + ").");
	*/
  }


  private void updateDrawableRect(int compWidth, int compHeight) 
  {
	int x = currentRect.x;
	int y = currentRect.y;
	int width = currentRect.width;
	int height = currentRect.height;
    
	//Make the width and height positive, if necessary.
	if (width < 0) {
	  width = 0 - width;
	  x = x - width + 1; 
	  if (x < 0) {
		width += x; 
		x = 0;
	  }
	}
	if (height < 0) {
	  height = 0 - height;
	  y = y - height + 1; 
	  if (y < 0) {
		height += y; 
		y = 0;
	  }
	}
    
	//The rectangle shouldn't extend past the drawing area.
	if ((x + width) > compWidth) {
	  width = compWidth - x;
	}
	if ((y + height) > compHeight) {
	  height = compHeight - y;
	}
          
	//Update rectToDraw after saving old value.
	if (rectToDraw != null) {
	  previousRectDrawn.setBounds(
								  rectToDraw.x, rectToDraw.y, 
								  rectToDraw.width, rectToDraw.height);
	  rectToDraw.setBounds(x, y, width, height);
	} else {
	  rectToDraw = new Rectangle(x, y, width, height);
	}
  }




  // This method required by MouseListener.
  public void mouseMoved(MouseEvent e){}
  
  // These methods are required by MouseMotionListener.
  public void mouseClicked(MouseEvent e){}
  public void mouseExited(MouseEvent e){}
  public void mouseEntered(MouseEvent e){}
  


  // Key Listeners
  /*
  public void keyTyped(KeyEvent e) 
  {
  char key = e.getKeyChar();
  if (key == 'q')
  {
  System.exit(0);
  }
  else if (key == 'v')
  {
  plotVint = !plotVint;
  }
  else if (key == 'a')
  {
  plotAI = !plotAI;
  }
  else if (key == 'b')
  {
  plotDT_blocked = !plotDT_blocked;
  }
  else if (key == 'p')
  {
  // print
  PrinterJob printJob = PrinterJob.getPrinterJob();
  // FIXME imaging area problems remain here.....
  printJob.setPrintable(this);
  if (printJob.printDialog()) 
  {
  try 
  {
  printJob.print();
  paint((Graphics2D)getGraphics());
  } 
  catch (Exception ex) 
  {
  ex.printStackTrace();
  }
  }
  }
  else if (key == 's')
  {
  grab(false);
  }
  }
  public void keyPressed(KeyEvent e) {}
  public void keyReleased(KeyEvent e)
  {
	if (e.getKeyChar() == 'v' || e.getKeyChar() == 'a' || e.getKeyChar() == 'b')
	  {
		repaint();
	  }
  }
  */

  // Naive printing interface
  public int print(Graphics g, PageFormat pf, int pi) throws PrinterException 
  {
	int sz_old;
    if (pi >= 1) 
	  {
		return Printable.NO_SUCH_PAGE;
	  }
	// this doesn't work properly yet
	Graphics2D g2 = (Graphics2D) g;
	

	if (!true)
	  {
		paint(g);
	  } 
	else 
	  {
		// set size according to paper width
		sz_old = sz;
		//sz = 1500; //(int) pf.getImageableWidth();
		//System.out.println("width = "+	pf.getImageableWidth()+" offset="+pf.getImageableX());
		//g2.translate((int)(pf.getImageableX()), pf.getImageableY());


		// set basic sizes
		setBasicSizes();


		// save context Xform
		AffineTransform saveXform = g2.getTransform();
		// apply new one
		g2.setTransform(at);


		// printing requires adjustment to the page size
		double xScale = at.getScaleX();
		double yScale = at.getScaleY();

		boolean scaleToX = (pf.getImageableWidth() < pf.getImageableHeight());
		if (scaleToX)
		  {
			g2.translate((int)(pf.getImageableWidth()/at.getScaleX()),
						 (int)(pf.getImageableWidth()/at.getScaleY()));
		  }
		else
		  {
			g2.translate((int)(pf.getImageableHeight()/at.getScaleX()),
						 (int)(pf.getImageableHeight()/at.getScaleY()));
		  }
		// then draw away....

		// here , the main graphs we want to see
		paintGraphs(g2);


		// paint Zoom & mouse 3 in old coords
		//g2.setTransform(saveXform);

		// mouse 3 stuff
		//paintMouse3(g2);

		// zoom box
		//paintZoomBox(g2);

		//if (dumpScreenShot) pushGrab();
	  }

	sz = sz_old;
    return Printable.PAGE_EXISTS;
  }


  // menu listener
  public void actionPerformed(ActionEvent e) 
  {
	//...Get information from the action event...
	// match it on the basis of the key accelerator
	KeyStroke str = ((JMenuItem)(e.getSource())).getAccelerator();
	if (str.equals(KeyStroke.getKeyStroke('z'))) 
	  {
		useZoom = !useZoom;
		setDisplayWindowRange();
		repaint();
	  }
	else if (str.equals(KeyStroke.getKeyStroke('b'))) 
	  {
		plotDT_blocked = !plotDT_blocked;
		repaintGraphInfo(intVelocityRect);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('q'))) 
	  {
		System.exit(0);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('v'))) 
	  {
		plotVint = !plotVint;
		repaintGraphInfo(sonicRect);
		repaintGraphInfo(DTSRect);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('a'))) 
	  {
		plotAI = !plotAI;
		repaintGraphInfo(densityRect);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('e'))) 
	  {
		// experimental eps driver
		printToEPS(); 
	  }
	else if (str.equals(KeyStroke.getKeyStroke('p'))) 
	  {
		// print
		PrinterJob printJob = PrinterJob.getPrinterJob();
		// FIXME imaging area problems remain here.....
		printJob.setPrintable(this);
		if (printJob.printDialog()) 
		  {
			try 
			  {
				printJob.print();
				//paint((Graphics2D)getGraphics());
			  } 
			catch (Exception ex) 
			  {
				ex.printStackTrace();
			  }
		  }
	  }
	else if (str.equals(KeyStroke.getKeyStroke('l'))) 
	  {
		expandGraphRangesFully = !expandGraphRangesFully;
		if (plotGamma) repaintGraphInfo(gammaRect);
		repaintGraphInfo(intVelocityRect);
		repaintGraphInfo(sonicRect);
		repaintGraphInfo(DTSRect);
		repaintGraphInfo(densityRect);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('s'))) 
	  {
		dumpScreenShot=true;
		// fresh paint here: a repaint doesn't erase the menu...
		paint((Graphics2D)getGraphics());
		grabAndDump(false);
	  }
	else if (str.equals(KeyStroke.getKeyStroke('c'))) 
	  {
		blockInterpret = false;
		blockInterpList = null;
		// fresh paint here: a repaint doesn't erase the menu...
		repaint();
	  }
  }

  // KeyListener methods
  public void keyTyped(KeyEvent e) {}
  public void keyPressed(KeyEvent e) {}
  public void keyReleased(KeyEvent e)
  {
	if (e.getKeyCode() == KeyEvent.VK_UP)
	  {
		moveAllWavelets(true);
	  }
	else if (e.getKeyCode() == KeyEvent.VK_DOWN)
	  {
		moveAllWavelets(false);
	  }
  }

  private void repaintGraphInfo(GraphInfo gi)
  {
	if (currentScale==1.0)
	  {
		repaint(gi.boundingBox.x, gi.boundingBox.y,
				gi.boundingBox.width+1, gi.boundingBox.height+1);
	  }
	else repaint();
  }

  // Component Listener
  public void componentResized(ComponentEvent e) 
  {
	Rectangle R = e.getComponent().getBounds();
	// resize to smallest dim.
	// sz = (int)Math.min(R.getWidth(),R.getHeight());
	sz = (int)R.getHeight();
	sz_width = (int)R.getWidth();

	// scale fonts accordingly
	axesfont = new Font("Arial", Font.PLAIN, (int)Math.rint(Math.max(1,9.0*sz/DEFAULT_SIZE)));
	titlefont = new Font("Arial", Font.PLAIN, (int)Math.rint(Math.max(1,10.0*sz/DEFAULT_SIZE)));
	repaint();
  }


  public void componentHidden(ComponentEvent e) {}
  
  public void componentMoved(ComponentEvent e) {}

  public void componentShown(ComponentEvent e) {}

  private String getReducedFilename()
  {
	String fName = new String(MLgraph_filename);
	fName = fName.replaceAll("MLgraphing-file.","");
	fName = fName.replaceAll("MOST_LIKELY.","");
	fName = fName.replaceAll(".las","");
	fName = fName.replaceAll(".txt","");
	return fName;
  }

  public BufferedImage captureComponentUsingImageBuffer(Component c) 
  {
	Dimension preferredSize = c.getPreferredSize();
	BufferedImage img = new BufferedImage(preferredSize.width, preferredSize.height, BufferedImage.TYPE_INT_ARGB);
	c.paint(img.getGraphics());
	return img;
  }

  /**
	 Screen dumping stuff
  */
  public void grab()
  {
	if (dumpEPSScreenShotAndDie)
	  {
		printToEPS();
	  }
	else 
	  {
		// OSX & XQuartz has XgetImage bug: detect & workaround
		boolean bug = detectOSXscreencapturebug();
		if (bug)
		  {
			EH.println("OSX screendump bug, trying ghostscript...");
			if (printToPNG_systemhack()) return; // system hack
		  }
		// default
		printToPNG();
	  }
  }


  /** 
	  Dump eps version of screen using HEP libraries
  */
  void printToEPS()
  {
	String fileName = getReducedFilename();
	fileName = fileName+".eps";

	VectorGraphics g = null;
	try 
	  {
		g = new PSGraphics2D(new File(fileName),this); 
	  }
	catch (IOException ioe) 
	  {
		EH.warn("Couldn't open " + fileName + " for eps output");
	  }
	// make sure text comes out as fonts
	Properties p = new Properties();
	p.setProperty(org.freehep.graphicsio.AbstractVectorGraphicsIO.TEXT_AS_SHAPES,"false");
	g.setProperties(p);

	g.startExport(); 
	this.print(g); 
	g.endExport();   
 
	EH.println("Wrote EPS screen dump to "+fileName);
  }

  // OSX has dud createScreenCapture calls via Xquartz. Detect such a case
  // with this hack
  boolean detectOSXscreencapturebug()
  {
	int i,j,k,N=10;
	int Kbug = -16777216; // bad value returned by failed XGetimage
	try
	  {
		Robot rb = new Robot();
		BufferedImage img = rb.createScreenCapture(new Rectangle(0,0,N,N));
		for(i = 0; i < N ;i++ )
		  for(j = 0; j < N; j++ )
			{
			  k = img.getRGB(i,j);
			  if (k != Kbug) return false;
			}
	  }
	catch (Exception e){}
	return true;
  }
  
  /** 
	  Dump gif version of screen using Robot
  */

  void printToPNG()
  {	
	String fName = getReducedFilename();
	try 
	  {
	   
		Robot rb = new Robot();
		Rectangle R = getBounds();
		Point P = getLocationOnScreen();
		R.setLocation(P);		

		BufferedImage img = rb.createScreenCapture(R);
		// Save as jpg, gif etc
		String type="png";
		fName = fName+"."+type;
		EH.println("Trying screen dump to "+fName);
		File file = new File(fName);
		if (file!=null)
		  {
			try
			  {
				ImageIO.write(img,type, file);
			  }
			catch (Exception e)
			  {
				EH.println("Can't write to "+fName);
			  }
		  }
		else
		  {
			EH.println("Can't write to "+fName);
		  }			  
	  }
	catch (Exception e)
	  {
		e.printStackTrace();
	  }
	EH.println("Wrote "+fName);
  }
  
  // use EPS->png conversion via ghostscript if needed
  boolean printToPNG_systemhack()
  {
	File f;
	int exitVal;
	// test if gs is on system
	try
	  {
		Process proc = Runtime.getRuntime().exec("gs -dBATCH");
		proc.waitFor();
		exitVal = proc.exitValue();
	  }
	catch (Exception e)
	  {
		EH.warn("Invoking gs problem: "+e);
		return false;
	  }
	if (exitVal<0) return false;
	try
	  {
		f = File.createTempFile("_WE",null);
		//f.deleteOnExit();
	  }
	catch (Exception e)
	  {
		EH.warn("Invoking gs problem: "+e);
		return false;
	  }
	String fileName = f.getAbsolutePath();
	VectorGraphics g = null;
	try 
	  {
		g = new PSGraphics2D(new File(fileName),this); 
	  }
	catch (IOException ioe) 
	  {
		EH.warn("Invoking gs problem: "+ioe);
		EH.warn("Couldn't open " + fileName + " for eps output");
	  }
	
	// Dump eps version to tmp file
	// make sure text comes out as fonts
	Properties p = new Properties();
	p.setProperty(org.freehep.graphicsio.AbstractVectorGraphicsIO.TEXT_AS_SHAPES,"false");
	g.setProperties(p);

	g.startExport(); 
	this.print(g); 
	g.endExport();

	String fileNameG = getReducedFilename()+".png";
	// ghostscript conversion from eps to png
	// resolution r300 = 300dpi
	String cmd = "gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 -sOutputFile="+fileNameG+" "+fileName;
	try
	  {
		EH.println("Trying: "+cmd);
		Process proc = Runtime.getRuntime().exec(cmd);
		BufferedReader sout=new BufferedReader(new InputStreamReader(proc.getInputStream()));
        String s;
        while((s=sout.readLine())!=null){EH.println(s);}
		proc.waitFor();
		exitVal = proc.exitValue();
		EH.println(cmd+"%nWrote "+fileNameG+" from "+fileName);
		EH.println("gs process exitValue: " + exitVal);
	  }
	catch (Exception e)
	  {
		e.printStackTrace();
	  }
	if(exitVal>=0)
	  {
		EH.println("Wrote raster screen dump to "+fileNameG);
		return true;
	  }
	return false;
  }

	
  /**
	 Class containing most of the useful info for each graph
  */
  private static class GraphInfo
  {
	Rectangle boundingBox;
	// graphing x bounds
	int xMin;
	int xMax;
	// bounds of current plotted function f
	double fMin,fMax;
	// ptr to current data
	double[] f;

	private GraphInfo(int x,int y,int width,int height) 
	{
	  boundingBox = new Rectangle(x,y,width,height);
	  // bounds
	  // Desire ratio of graph excursions to rectangle width to be not less than
	  // XPLOT_TO_XRECT_RATIO_MIN, but also margin not more than 
	  // ExtractionViewer.GRAPHING_XMARGIN pixels.

	  int xMargin = (int) Math.min((1.0-ExtractionViewer.XPLOT_TO_XRECT_RATIO_MIN)*width/2.0, ExtractionViewer.GRAPHING_XMARGIN);

	  xMin = boundingBox.x+xMargin;
	  xMax = boundingBox.x+boundingBox.width-xMargin;
	}
  }

  /**
	 Simple data class for wavelet data
  */
  private static class WaveletData
  {
	int ns;
	double delta_t,f1,tMax;
	double[] values;
  }
}
