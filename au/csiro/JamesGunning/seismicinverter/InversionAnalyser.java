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


package au.csiro.JamesGunning.seismicinverter;

import java.util.*;
import java.io.*;
import java.awt.*;
import static java.lang.Math.*;

import cern.colt.list.*; 
import cern.jet.stat.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import cern.colt.map.*;
import hep.aida.*;
import hep.aida.ref.*;

import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.stat.util.MCMCConvergenceMonitor;
import au.csiro.JamesGunning.lib.graphutil.*;


/**
   Basic class for reading output traces of
   seismic inverter and extracting relevant properties
   
   @author james.gunning@csiro.au
   @version 1.0, 15/2/2002
*/

public class InversionAnalyser
{
  ErrorHandler EH;

  // runtime parameters
  boolean graphics = false;
  private boolean layerCake = false;
  private boolean layerCakeOfFluids = false;
  private boolean layerCakeOfNG = false;
  private boolean layerCakeOfZ_eff = false;
  private boolean layerCakeOfProp = false;
  private String  layerCakeProp = null;
  private boolean layerCakesInDepth = false;
  private boolean layerRangeSpecified = false;
  private String  layerColors = null;
  private boolean layersMakeFirstTraceInRealisationsBlack = false;
  private ArrayList layerColorArray = null;
  private boolean fluidContentMap = false;
  private boolean showHorizonsAndSeismic=false;
  private boolean showOverLayerCake=false;
  private boolean propertySUMap = false;
  private boolean propertyMap = false;
  private boolean propertyMap2 = false;
  private boolean propertyMap2gen = false;
  private boolean makePropertyAsciiSummary = false;
  private boolean makeSyntheticSeismic = false;
  private boolean performTimeToDepthRemapping = false;
  private boolean computeMean = false;
  private boolean computeStddev = false;
  private boolean computeFluidProb = false;
  private boolean computeQuantile = false;
  private double  quantile=0.50; // default quantile
  private boolean makeFullSummaryTraces = false;
  private boolean makeFullMeanSummaryTraces = false;
  private boolean makeFullStdDevSummaryTraces = false;
  private boolean makeFullQuantileSummaryTraces = false;
  private boolean makeFullAsciiSummary = false;
  private boolean makeTimesAsciiSummary = false;
  private boolean applyFilter=false;
  private boolean filterPad=false;
  private TraceFilter traceFilter=null;
  private boolean doBHPcommand=false;
  private boolean makePriors=false;
  private boolean massageAnalyse=false;
  private boolean massageAnalyseAscii=false;
  private boolean computeHistogram=false;
  private boolean computeRafteryThinning=false;
  private boolean waitForStdin=false;
  boolean useJNI=true;

  int verbosity=0;

  String inversions_su_filenames = null;
  String realisation_filename= null;
  String wavelet_filename = null;
  String seismic_filename = null;

  public static final String NET_SAND_KEYWORD = "net-sand";
  public static final String NET_HYDROCARBON_KEYWORD = "net-hydrocarbon";
  public static final String NET_OIL_KEYWORD = "net-oil";
  public static final String NET_LSG_KEYWORD = "net-lsg";
  public static final String NET_GAS_KEYWORD = "net-gas";
  public static final String THICKNESS_KEYWORD = "thickness";
  public static final String THINNING_KEYWORD = "sample-correlation";
  // Keywords to detect fluid probabilities
  public static final String BRINE_PROB_KEYWORD = "fluid_type_0";
  public static final String OIL_PROB_KEYWORD   = "fluid_type_1";
  public static final String LSG_PROB_KEYWORD   = "fluid_type_2";
  public static final String GAS_PROB_KEYWORD   = "fluid_type_3";

  public static final String NET_HARD_KEYWORD = "net-hard";
  public static final String PHIE_KEYWORD = "phi_e";
  
  /**
	 Useful keys for fast mapping
  */
  public static final int NET_HYDROCARBON_KEY=0;
  public static final int NET_SAND_KEY=1;
  public static final int THICKNESS_KEY=2;
  public static final int THINNING_KEY=3;
  public static final int BRINE_PROB_KEY = 4;
  public static final int OIL_PROB_KEY   = 5;
  public static final int LSG_PROB_KEY   = 6;
  public static final int GAS_PROB_KEY   = 7;
  public static final int NET_OIL_KEY =8;
  public static final int NET_LSG_KEY =9;
  public static final int NET_GAS_KEY =10;
  public static final int NET_HARD_KEY=11;
  public static final int PHIE_KEY=12;
  
  public static final float LAYER_CAKE_DEFAULT_DT = 0.001f;


  // Member variables
  private ArrayList keyWords;
  private Map keyWordsMap;
  private ArrayList completeWords;

  // massaging entities
  /**
	 Flag for writing efficient binary outputs
  */
  public static final boolean WRITE_SCALE_COMPRESSED_BINARY_MOMENT_ANALYSIS = true;
  private int compressed_bytes_per_float=2;

  private ArrayList massagingKeywords;
  private String[] massagingHeaderWords;
  private String[] massagingHeaderWordsValues;
  private String massagingOutputFileName;
  private SuReader sureader;
  private SuWriter suwriter;
  DoubleMatrix2D sampleCovNormed;
  DoubleMatrix2D sampleCovNormedDiff=null;
  private byte[] work_bytes = new byte[4];


  // Prior traces
  private String[] priorKeywords;
  private ArrayList priorKeywordsArray;



  /**
	 model traces
  */
  private SuSegyFormat tr;
  private SuSegyFormat tr_out;

  /** 
	  input stream for inverted traces
  */
  private FileInputStream inversion_in_strm;
  private FileOutputStream out_strm;
  private PrintWriter aux_out_writer;
  private FileOutputStream aux_out_binary_writer;
  private BufferedOutputStream aux_out_binary_writer_buf;

  // useful variables for layer cake viewing
  private int ns, new_ns, n_layers=1;
  private float t_start, d_start, t_end, d_end;
  private float delta_t = LAYER_CAKE_DEFAULT_DT;
  private float[] layer_color_step_values=null;
  private float largest_layer_color=-Float.MAX_VALUE;


  // property fetchers 
  TracePropertyFetcher propFetcher0;   // for property 0
  TracePropertyFetcher propFetcher1;   // for property 1
  TracePropertyFetcher propFetcher2;   // for property 2
  TracePropertyFetcher propFetcher3;   // for property 3
  TracePropertyFetcher propFetcher4;   // for property 4
  TracePropertyFetcher filterPropertyFetcher0;   // for property 2 (in filtering)

  FetcherSet synthSeisFetchers=null;

  // useful for most modes

  private int layerNumber;
  private int layerNumber2; // only used for propertyMap2gen
  private int refIndx0;     // generally useful indexes for positions in model trace
  private int refIndx1;  
  private int refIndx_t;

  // useful for property viewing
  private String property_name1;
  private String property_name2;
  private ArrayList columnTitles=null;

  // for filtering
  private int     filterLayerNumber;
  private char    filter_operation='=';
  private double  filter_value=0;
  private String  filter_property_name;
  private int     filterIndx;

  // for fluid probabilities
  private int desiredFluidType=0;

  // for histograms: force common binning (use widest stats first)
  boolean isFirstHistogram=true;
  double histogram_min,histogram_max;

  // BHP commands to be forked as processes, formed partially from this string
  private String BHPcommandstring="";
  private String BHPcommand_filename;
  
  // generally useful things
  private int trace_number=0;
  private int accepted_trace_number=0;
  private int accepted_tracegroup_number=0;
  private int Ns=0;
  // header entry of tr associated with the realisation_number, and MF.realisationNumberHeaderField
  private SuSegyFormat.SuHeaderEntry realisationNumberSuHeaderMapEntry;


  protected ModelFormat MF;
  private StringBuffer propertyWords;

  // Useful variables for synthetic seismics
  private double tSyntheticStart;    // start of synthetic seismics
  private double tSyntheticEnd;      // end of synthetic seismics
  private OneDSeismicConvolutionalModel SM;
  private OneDSeismicConvolutionalModel.Wavelet wavelet;

  // Useful variables for acquiring stats across multiple traces
  boolean realisationGroupIsNull=false;
  StringBuilder hwb = new StringBuilder();
  private DoubleArrayList dataList;
  private int lastRealisationNumberHeaderField=-1;
  private int groupOffsetNumber = 0;
  // and similarly for simultaneous averages
  private DoubleArrayList[] dataLists;

  // and stored reflection/time data 
  private HashMap SyntheticSeismicDataMap=null;

  private boolean gotFirstTrace = false;

  // Material for graphics
  SeismicGraphicInfo SGI=null;
  OpenIntObjectHashMap GraphingDataMap = null;

  protected InversionAnalyser()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();

	// Default ***BASIC*** Model formats. Overwritten by setUpModelFormats() based on
	// first trace header info later on if needed
	MF = new ModelFormat(true,0);	
	assembleWordsAndMaps();

	// Always need a trace filter object
	traceFilter = new TraceFilter(EH,0); // update verbosity later
  }

  private void setUpModelFormats(int format_version)
  {
	MF = new ModelFormat(false, format_version);	
	assembleWordsAndMaps();
  }

  private void assembleWordsAndMaps()
  {
	int i;

	// Assemble help text for property words in Model Format
	propertyWords = new StringBuffer();
	for (i=0;i<MF.reducedPropertyList.size();i++) 
	  {
		propertyWords.append((String)MF.reducedPropertyList.get(i)+" ");
	  }
	
	// Extra things we'd like to see
	keyWords=new ArrayList();
	keyWordsMap = new HashMap();

	// Raftery thinning
	keyWords.add(THINNING_KEYWORD);
	keyWordsMap.put(THINNING_KEYWORD,new Integer(THINNING_KEY));

	// layer thickness
	if (MF.reducedPropertyList.contains("vp_eff"))
	  {
		keyWords.add(THICKNESS_KEYWORD);
		keyWordsMap.put(THICKNESS_KEYWORD,new Integer(THICKNESS_KEY));
	  }
	// net sand per layer
	if (MF.reducedPropertyList.contains("vp_eff") && MF.reducedPropertyList.contains("NG"))
	  {
		keyWords.add(NET_SAND_KEYWORD);
		keyWordsMap.put(NET_SAND_KEYWORD,new Integer(NET_SAND_KEY));
	  }
	// net hard-rock per layer
	if (MF.reducedPropertyList.contains("vp_eff") && MF.reducedPropertyList.contains("NG_hard"))
	  {
		keyWords.add(NET_HARD_KEYWORD);
		keyWordsMap.put(NET_HARD_KEYWORD,new Integer(NET_HARD_KEY));
	  }
	// phi_e
	if (MF.reducedPropertyList.contains("NG") && MF.reducedPropertyList.contains("phi_s"))
	  {
		// phi_e
		keyWords.add(PHIE_KEYWORD);
		keyWordsMap.put(PHIE_KEYWORD,new Integer(PHIE_KEY));
	  }
	// net hydrocarbon per layer
	if (MF.reducedPropertyList.contains("vp_eff") && MF.reducedPropertyList.contains("NG")
		&& MF.reducedPropertyList.contains("phi_s"))
	  {
		if (MF.reducedPropertyList.contains("S_h"))
		  {
			keyWords.add(NET_HYDROCARBON_KEYWORD);
			keyWordsMap.put(NET_HYDROCARBON_KEYWORD,new Integer(NET_HYDROCARBON_KEY));
		  }
		// Subsequent stuff for contacts model....
		if (MF.reducedPropertyList.contains("S_oil"))
		  {
			keyWords.add(NET_OIL_KEYWORD);
			keyWordsMap.put(NET_OIL_KEYWORD,new Integer(NET_OIL_KEY));
		  }
		if (MF.reducedPropertyList.contains("S_lsg"))
		  {
			keyWords.add(NET_LSG_KEYWORD);
			keyWordsMap.put(NET_LSG_KEYWORD,new Integer(NET_LSG_KEY));
		  }
		if (MF.reducedPropertyList.contains("S_gas"))
		  {
			keyWords.add(NET_GAS_KEYWORD);
			keyWordsMap.put(NET_GAS_KEYWORD,new Integer(NET_GAS_KEY));
		  }
	  }

	if (MF.reducedPropertyList.contains("fluid_type"))
	  {
		// Brine probability
		keyWords.add(BRINE_PROB_KEYWORD);
		keyWordsMap.put(BRINE_PROB_KEYWORD,new Integer(BRINE_PROB_KEY));
		// Oil probability
		keyWords.add(OIL_PROB_KEYWORD);
		keyWordsMap.put(OIL_PROB_KEYWORD,new Integer(OIL_PROB_KEY));
		// LSG probability
		keyWords.add(LSG_PROB_KEYWORD);
		keyWordsMap.put(LSG_PROB_KEYWORD,new Integer(LSG_PROB_KEY));
		// Gas probability
		keyWords.add(GAS_PROB_KEYWORD);
		keyWordsMap.put(GAS_PROB_KEYWORD,new Integer(GAS_PROB_KEY));
	  }

	// tack these on to the list of properties we can use
	for (i=0;i<keyWords.size();i++) 
	  {
		propertyWords.append((String)keyWords.get(i)+" ");
	  }
	// full list of words in Model format plus special keywords
	completeWords = new ArrayList(keyWords);
	completeWords.addAll(MF.reducedPropertyList);
  }
  
  /** 
	  Sample driver from input streams etc
  */
  public static void main(String argv[])
  {

	InversionAnalyser IA = new InversionAnalyser();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("InversionAnalyser",IA.EH);

	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber);

	// some decoration to make help more intelligible
	options.addMainCommandlineText(" < stdin");
	options.addSpielText("Conceivable properties for PROP are:%n"+IA.propertyWords+ "%n "+
						 "- sample-correlation is an indication of the correlation length of MCMC realisations "+
						 "in multi-realisation runs (should be small for good statistics).%n "+
						 " - fluid_type_0,1,2,3 are 0,1 indicators for brine, oil, lsg and gas. "+
						 "Their means are the probabilities of these fluids.%n "+
						 " - fluid_type > 3 values are encodings for layers with internal fluid contacts" +
						 " - Note 1: seismic_chisquare from the *last* layer is the full -ve log-posterior (i.e. not scaled)"+
						 " - Note 2: R_[near|far][_fc] means the associated R_ps if the associated stack has isPSdata = true%n");


	// program options
	options.addOption("-h","--help",0,null,null,"Show this help");
	options.addOption("-G",null,0,null,null,"Use inbuilt graphics to show results. Available on commands marked [G]. Suppresses usual output");
	options.addOption("-v",null,1,"int","level","Execute with stated level of verbosity (to stderr)");
	options.addOption("-i",null,1,"String","realisations.su","Use named file as SU realisations file, not stdin");
	options.addOption("-l",null,0,null,null,"Produce layer-cake image from mean sample times (SU format to stdout) [G]");
	options.addOption("-lf",null,0,null,null,"Produce layer-cake image from mean sample times of fluid content (SU format to stdout) [G]");
	options.addOption("-lNG",null,0,null,null,"Produce layer-cake image from mean sample times of NG content (SU format to stdout) [G]");
	options.addOption("-lZ_eff",null,0,null,null,"Produce layer-cake image from mean sample times of normal impedance content (SU format to stdout) [G]");
	options.addOption("-lprop",null,1,"String","prop","Produce layer-cake image of property prop (SU format to stdout) [G]");
	options.addOption("-Sl",null,1,"String","seismic_filename.su",
					  "Display horizons on traces in (big-endian) seismic_filename.su. Java graphics (-G). "
					  +"Use --layer to adjust range if required. Shows model picks as little bars");
	options.addOption("-Slo",null,1,"String","seismic_filename.su",
					  "Same as -Sl, but overlays seismic on model");
	options.addOption(null,"--layer",2,"float","t_start(sec) t_end(sec)","Produce layer-cake image from sample "+
					  "times over specified range (SU format to stdout) [G]");
	options.addOption("-ldepth",null,0,null,null,"Write layer-cake type images (e.g. -l, -lNG) in depth, not time [G]");
	options.addOption("--layers-add-ensemble-boundaries",null,0,null,null,"Add \"black\" trace for first "
					  +"realisation when doing --layers");
	options.addOption("--layer-colors",null,1,"String","colors","with -l option, specifies successive colors "+
					  "in grayscale (SU output stream only), e.g. 1,2,1,3");
	options.addOption("-f",null,0,null,null,"Produce layer-cake fluid-content map (SU format to stdout)");
	options.addOption("--time-to-depth",null,0,null,null,"remap time model traces to depth-model traces (SU format to stdout). "+
					  "Must be used stand-alone as a filter (output can be piped into another analysis call)");
	options.addOption("--full",null,0,null,null,"Produce full posterior summary model traces of everything of "+
					  "interest (SU format to stdout)");	
	options.addOption("--full-means",null,0,null,null,"Produce full mean summary model traces of everything of "+
					  "interest (SU format to stdout)");	
	options.addOption("--full-stddevs",null,0,null,null,"Produce full stddev summary model traces of everything "+
					  "of interest (SU format to stdout)");
	options.addOption("--full-quantile",null,1,"String","Q","Produce full Qth [0-100] quantile summary model "+
					  "traces of everything of interest (SU format to stdout)");
	options.addOption("-p",null,2,"String","PROP N","Extract property PROP in layer N (SU format to stdout, "+
					  "1 sample per trace)"+
					  ". See available properties below.");
	options.addOption("-p-ascii",null,2,"String","PROP N","Print property PROP in layer N (Ascii to stdout)"+
					  ". See available properties below. [G]");
	options.addOption("-p-ascii2",null,3,"String","PROP1 PROP2 N","Produce scatterplot of property PROP1 vs PROP2 "+
					  "in layer N (2 column ascii to stdout). "+
					  "2 Column ascii output stream. See available properties below. [G]");
	options.addOption("-p-ascii2-gen",null,4,"String","PROP1 PROP2 N1 N2","Produce scatterplot of property "+
					  "PROP1 (layer N1) vs PROP2 (layer N2) (2 column ascii to stdout)."+
					  ". See available properties below. [G]");
	options.addOption("-p-ascii-t",null,0,null,null,"Print all layer times out in ascii form, including t_base (ascii to stdout)");
	options.addOption("-p-ascii-prop",null,1,"String","PROP","Print table of "+
					  "property PROP, one column per layer (ascii to stdout)");	
	options.addOption("-p-ascii-full",null,0,null,null,"Produce full posterior summary ascii table of "+
					  "everything of interest (ascii to stdout)");
	options.addOption(null,"--mean",2,"String","PROP N","Compute average of property PROP in layer N (Ascii to stdout)"+
					  ". See available properties below.");
	options.addOption(null,"--stddev",2,"String","PROP N","Compute std. deviation of property PROP in layer N (Ascii to stdout)"+
					  ". See available properties below.");
	options.addOption(null,"--quantile",3,"String","PROP N Q","Compute Qth [0-100] quantile of property "+
					  "PROP in layer N (Ascii to stdout)"+
					  ". See available properties below.");
	options.addOption(null,"--histogram",2,"String","PROP N","Produce histogram of property PROP in layer N "+
					  "(2 column Ascii to stdout)"+
					  ". See available properties below. [G]");
	options.addOption(null,"--fluid-prob",2,"String","F N","Compute probability of fluid_type F in layer N (Ascii to stdout). "+
					  "F must be 0,1,2,3 for (brine, oil, low-sat-gas, gas) respectively.");
	options.addOption(null,"--Raftery",2,"String","PROP N","Compute Raftery thinning-parameter of property "+
					  "PROP in layer N (Ascii to stdout)"+
					  ". See available properties below.");

	options.addOption("-s",null,4,"String","t_0(sec) t_last(sec) R_[near|far][_fc] wavelet_filename.su",
					  "Produce synthetic seismic images from  refl. coefficients (including fluid contacts). SU to stdout. "+
					  "With _fc, show fluid-contact synthetics only [G]. %nExample: -s 1.2 1.8 R_near near_wavelet.su");
	options.addOption("-S",null,5,"String","t_0(sec) t_last(sec) R_[near|far][_fc] wavelet_filename.su seismic_filename.su",
					  "Same as -s, but appends traces in (big-endian) seismic_filename.su as well, gated to same time window [G]");
	options.addOption(null,"--realisation-header-word",1,"String","name","Use this header field for the realisation number");
	options.addOption(null,"--filter",4,"String","PROP [=|>|<] value N","Restrict analysis to traces "+
					  "satisfying stated [in]equality on property PROP in layer N. "+
					  "(Shell protect <>= if required). E.g. --filter 'phi_s > 0.2 2'. See available properties below");
	options.addOption(null,"--filter-pad",4,"String","PROP [=|>|<] value N","Same as --filter, except that a fully-exclusive " +
					  "filter (i.e. one that excludes all traces of a realisation-group) will generate a -999 filled single \"null\" trace. " +
					  "If used, this filter operation must be the only deliveryAnalyser operation (i.e. chain your filters).");
	options.addOption("--trace-filter",null,1,"String","Expression",
					  "operate only on traces satisfying shell-protected, comma-separated header word filter "+
					  "expressions involving comparisons <,>,=,!=,>=,<=,%% e.g. \"fldr>100,ep=234\"");
	options.addOption("--BHPcommand",null,1,"String","filename",
					  "Write to filename the string \"properties=x,x..\", which lists model-properties associated with the outputs from --full, --full-means, --full-stddevs, --full-quantile");
	options.addOption(null,"--wait",0,null,null,"Wait indefinitely for stdin traces");
	/*
	  options.addOption(null,"--linear-comb",6,"String","A B C PROP1 PROP2 N","Form linear combination A + B PROP1 + "
	  +"C PROP2 in layer n. See available properties below.");
	*/

	options.addOption("--make-prior-traces",null,1,"String",
					  " prop1,prop2,...",
					  "Make an SU model prior file from p50 quantile of named properties "+
					  "Example: --make-prior-traces NG,t");

	options.addOption("--massage-analyse",null,4,"String","header_word1,header_word2.. layer-number1,layer-number2... "
					  +" prop1,prop2,... output_file_name",
					  "Form summary statistics for massaging. Write ascii-header+binary-data file. Example: --massage-analyse cdp,ep 2 net-sand,NG massage-analysis.mab");
	options.addOption("--massage-analyse-ascii",null,4,"String","header_word1,header_word2.. layer-number1,layer-number2... "
					  +" prop1,prop2,... output_file_name",
					  "Same as --massage-analyse, but with bloated ascii header+data output");
	options.addOption(null,"--wait",0,null,null,"Wait indefinitely for stdin traces",1);
	options.addOption(null,"--no-jni",0,null,null,"Disable native libraries (use java fallbacks)");

	// Parse allowing for arbitrary numbers of commandline files
	if (!options.parseOptions(argv,true))
	  {
		 IA.EH.die_err("");	
	  }

	if (options.isOptionSet("--no-jni")) IA.useJNI = false;  

	// (1) try and get a trace to determine format variations if any from 1st trace
	if (options.isOptionSet("-i"))
	  {
		IA.realisation_filename = ((String)options.getOptionValues("-i")[0]);
	  }
	if (options.isOptionSet("-G"))
	  {
		IA.graphics = true;
	  }
	if (options.isOptionSet("--wait"))
	  {
		IA.waitForStdin=true;
	  }
	// (2) Make sure stdin is there
	// test we have something, or show Usage. Wait 1 sec to get something before
	// bailing out.
	if (IA.realisation_filename==null)
	  {
		// use stdin
		IA.inversion_in_strm = (FileInputStream) new FileInputStream(FileDescriptor.in);
		try 
		  {
			// FIXME ! Can anyone think of a better way to detect if stdin is not a console
			// but has no bytes available. The usage() call below should happen if
			// the input is a console. Need java equivalent of POSIX's isatty()
			boolean haveStdIn=false;
			long tEnd = System.currentTimeMillis()+1000;
			for (;;)
			  {
				if (IA.inversion_in_strm.available()!=0)
				  {
					haveStdIn=true;
					break;
				  }
				if (System.currentTimeMillis()>tEnd) break;
			  }
			if (!haveStdIn && !IA.waitForStdin)
			  {
				IA.EH.warn("No su data found on stdin. Check usage");
				IA.inversion_in_strm.close();
				options.Usage();
			  }
		  }
		catch (IOException e) {throw new RuntimeException("couldn't check stdin bytes available");}
	  }

	// OK, there's data. Try to get the model format
	IA.setUpTraceInputs();
	IA.initialiseTraceCacheAndGetFormat();

	// Stuff commandline options into InversionAnalyser fields
	if (options.isOptionSet("-h"))
	  {
		options.Usage();
	  }

	
	if (options.isOptionSet("--realisation-header-word"))
	  {
		IA.MF.realisationNumberHeaderField = 
		  ((String)options.getOptionValues("--realisation-header-word")[0]);
	  }
	if (options.isOptionSet("-v")) 
	  IA.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();

	if (options.isOptionSet("-l") || options.isOptionSet("--layer") || 
		options.isOptionSet("-lf") || options.isOptionSet("-lNG") || 
		options.isOptionSet("-lZ_eff") || options.isOptionSet("-lprop"))
	  {
		if (options.isOptionSet("-l"))
		  {
			IA.layerCake= true;
		  }
		else if (options.isOptionSet("-lf")) 
		  {
			IA.layerCakeOfFluids = true;
		  }
		else if (options.isOptionSet("-lNG")) 
		  {
			IA.layerCakeOfNG = true;
		  }
		else if (options.isOptionSet("-lZ_eff")) 
		  {
			IA.layerCakeOfZ_eff = true;
		  }
		else if (options.isOptionSet("-lprop")) 
		  {
			IA.layerCakeOfProp = true;
			IA.layerCakeProp = options.getOptionString("-lprop",0);
		  }
		if (options.isOptionSet("-ldepth"))
		  {
			IA.layerCakesInDepth = true;
		  }
		if (options.isOptionSet("--layer"))
		  {
			IA.layerRangeSpecified = true;
			if (!IA.layerCakesInDepth)
			  {
				IA.t_start = options.getOptionFloat("--layer",0);
				IA.t_end = options.getOptionFloat("--layer",1);
			  }
			else
			  {
				IA.d_start = options.getOptionFloat("--layer",0);
				IA.d_end = options.getOptionFloat("--layer",1);
			  }
		  }
		if (options.isOptionSet("--layer-colors"))
		  {
			IA.layerColors= ((String)options.getOptionValues("--layer-colors")[0]);
			// colors must be comma or colon separated 
			StringTokenizer args = new  StringTokenizer(IA.layerColors,",:");
			if (args!=null)
			  {
				IA.layerColorArray = new ArrayList();
				while (args.hasMoreTokens())
				  {					
					IA.layerColorArray.add(new Integer(Integer.parseInt(args.nextToken())));
				  }
			  }
		  }
		if (options.isOptionSet("--layers-add-ensemble-boundaries"))
		  {
			IA.layersMakeFirstTraceInRealisationsBlack = true;
		  }
	  }
	if (options.isOptionSet("-f")) 
	  IA.fluidContentMap = true;

	if (options.isOptionSet("--full")) 
	  IA.makeFullSummaryTraces = true;

	if (options.isOptionSet("--full-means")) 
	  IA.makeFullMeanSummaryTraces = true;

	if (options.isOptionSet("--full-stddevs")) 
	  IA.makeFullStdDevSummaryTraces = true;

	if (options.isOptionSet("--full-quantile")) 
	  {
		IA.makeFullQuantileSummaryTraces = true;
		try 
		  {
			// convert percentile to quantile in (0,1)
			IA.quantile = 0.01*Double.parseDouble(((String)options.getOptionValues("--full-quantile")[0]));			
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("Quantile must be a number");
			options.Usage();
		  }			
		if (IA.quantile<0 || IA.quantile >1.0)
		  {
			IA.EH.die_err("Quantile must be in range [0,100]");
		  }
	  }
	if (options.isOptionSet("-p")) 
	  {
		IA.propertySUMap = true;
	  }
	else if (options.isOptionSet("-p-ascii")) 
	  {
		IA.propertyMap = true;
		IA.propertySUMap = false;  // switch this off. Can't have both
	  }
	else if (options.isOptionSet("-p-ascii2")) 
	  {
		IA.propertyMap2 = true;
		IA.propertyMap = false;    // switch this off. Can't have both
		IA.propertySUMap = false;  // switch this off. Can't have both
	  }
	else if (options.isOptionSet("-p-ascii2-gen")) 
	  {
		IA.propertyMap2gen = true;
		IA.propertyMap = false;    // switch this off. Can't have both
		IA.propertySUMap = false;  // switch this off. Can't have both
	  }
	else if (options.isOptionSet("-p-ascii-prop")) 
	  {
		IA.makePropertyAsciiSummary = true;
		IA.propertySUMap = false;  // switch this off. Can't have both
	  }
	else if (options.isOptionSet("-p-ascii-full")) 
	  {
		IA.makeFullAsciiSummary = true;
	  }
	else if (options.isOptionSet("-p-ascii-t")) 
	  {
		IA.makeTimesAsciiSummary = true;
	  }
	if (options.isOptionSet("--time-to-depth")) 
	  {
		IA.performTimeToDepthRemapping = true;
	  }
	if (options.isOptionSet("-s") || options.isOptionSet("-S")) 
	  {
		String flag = (options.isOptionSet("-S"))? "-S" : "-s";
		IA.makeSyntheticSeismic = true;
		IA.property_name1 = ((String)options.getOptionValues(flag)[2]);
		// if producing synthetic seismic for e.g. mean property, do
		if (options.isOptionSet("--mean") || options.isOptionSet("--quantile"))
		  {
			IA.property_name2 = ((String)options.getOptionValues(flag)[2]);
		  }
		if (options.isOptionSet("--stddev"))
		  {
			IA.EH.warn("-s and --stddev together make no sense. Try -s and --quantile or --mean");
			options.Usage();
		  }
		IA.wavelet_filename = ((String)options.getOptionValues(flag)[3]);
		if (options.isOptionSet("-S")) IA.seismic_filename = ((String)options.getOptionValues("-S")[4]);
		// tests
		try
		  {
			IA.tSyntheticStart = Double.parseDouble(((String)options.getOptionValues(flag)[0]));
			IA.tSyntheticEnd = Double.parseDouble(((String)options.getOptionValues(flag)[1]));
			if (IA.tSyntheticEnd <= IA.tSyntheticStart)
			  {
				IA.EH.warn("Synthetic end time must be greater than start time%n");
				options.Usage();
			  }
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("Synthetic end and start times must be numbers%n");
			options.Usage();
		  }
	  }
	if (options.isOptionSet("-Sl") || options.isOptionSet("-Slo")) 
	  {
		IA.showHorizonsAndSeismic = true;
		IA.graphics = true; // only sensible way to show this
		if (options.isOptionSet("-Sl")) IA.seismic_filename = options.getOptionString("-Sl",0);
		else 
		  {
			IA.seismic_filename = ((String)options.getOptionString("-Slo",0));
			IA.showOverLayerCake=true;
		  }
	  }
	if (options.isOptionSet("--fluid-prob") || options.isOptionSet("--mean") || 
		options.isOptionSet("--stddev") || options.isOptionSet("--quantile") || 
		options.isOptionSet("--histogram") || options.isOptionSet("--Raftery"))
	  {
		if (options.isOptionSet("--histogram"))
		  {
			IA.computeHistogram= true;
			IA.property_name1 = ((String)options.getOptionValues("--histogram")[0]);
		  }
		if (options.isOptionSet("--Raftery"))
		  {
			IA.computeRafteryThinning = true;
			IA.property_name1 = ((String)options.getOptionValues("--Raftery")[0]);
		  }
	    if (options.isOptionSet("--mean"))
		  {
			IA.computeMean = true;
			IA.property_name1 = ((String)options.getOptionValues("--mean")[0]);
		  }
	    if (options.isOptionSet("--stddev")) 
		  {
			IA.computeStddev = true;
			IA.property_name1 = ((String)options.getOptionValues("--stddev")[0]);
		  }
		if (options.isOptionSet("--fluid-prob")) 
		  {
			IA.computeFluidProb = true;
			IA.property_name1 = "fluid_type";
		  }
		if (options.isOptionSet("--quantile")) 
		  {
			IA.computeQuantile = true;
			IA.property_name1 = ((String)options.getOptionValues("--quantile")[0]);
			try 
			  {
				IA.quantile = 0.01*Double.parseDouble(((String)options.getOptionValues("--quantile")[2]));
			  }
			catch (NumberFormatException e)
			  {
				IA.EH.warn("Quantile must be a number");
				options.Usage();
			  }			
			if (IA.quantile<0 || IA.quantile >1.0)
			  {
				IA.EH.die_err("Quantile must be in range [0,100]%n");
			  }
		  }
		if (IA.completeWords.indexOf(IA.property_name1) == -1)
		  {
			IA.EH.warn("Property " + IA.property_name1 + " isn't available from inverted traces. Conceivable are:");
			IA.showAvailableProperties();
			IA.EH.die_err("%n");
		  }
		try
		  {
			if (options.isOptionSet("--histogram")) 
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--histogram")[1]));
			  }
			if (options.isOptionSet("--Raftery"))
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--Raftery")[1]));
			  }
			if (options.isOptionSet("--mean")) 
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--mean")[1]));
			  }
			if (options.isOptionSet("--stddev")) 
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--stddev")[1]));
			  }
			if (options.isOptionSet("--quantile")) 
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--quantile")[1]));
			  }
			if (options.isOptionSet("--fluid-prob")) 
			  {
				IA.layerNumber = 
				  Integer.parseInt(((String)options.getOptionValues("--fluid-prob")[1]));
				IA.desiredFluidType = 
				  Integer.parseInt(((String)options.getOptionValues("--fluid-prob")[0]));
			  }
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("Layer number must be a positive integer");
			options.Usage();
		  }
	  }

	if (options.isOptionSet("--filter") || options.isOptionSet("--filter-pad")) 
	  {
		String fopt;
		if (options.isOptionSet("--filter"))
		  {
			fopt = "--filter";
		  }
		else
		  {
			fopt = "--filter-pad";
			IA.filterPad = true;
		  }
		IA.applyFilter = true;
		IA.filter_property_name = ((String)options.getOptionValues(fopt)[0]);
		// parse, looking for --filter PROP [>|<|=] value layer_number structure
		// FIXME must be a better way to do this...
		String op = ((String)options.getOptionValues(fopt)[1]);
		IA.filter_operation = (op.equals("=")) ? '=' : 
		  ((op.equals("<")) ? '<' : '>');
		// filter value
		try 
		  {
			IA.filter_value = Double.parseDouble(((String)options.getOptionValues(fopt)[2]));
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("filter value must be a number");
			options.Usage();
		  }
		// now get layer number in which to apply filter
		try
		  {
			IA.filterLayerNumber = 
			  Integer.parseInt(((String)options.getOptionValues(fopt)[3]));
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("Filter layer number must be a positive integer");
			options.Usage();
		  }
		// checks
		if (IA.completeWords.indexOf(IA.filter_property_name) == -1)
		  {
			IA.EH.warn("Property " + IA.filter_property_name + " isn't available from inverted traces. Conceivable are:");
			IA.showAvailableProperties();
			IA.EH.die_err("%n");
		  }
	  }
	if (options.isOptionSet("--trace-filter")) 
	  {
		// need also a trace filter, even if its inactive. Ensure verbosity correct here
		IA.traceFilter.verbosity = IA.verbosity;
		IA.traceFilter.applyTraceFilter = true;
		IA.traceFilter.cmdLineTraceFilterString = (String)options.getOptionValues("--trace-filter")[0];
	  }


	// BHP commands
	if (options.isOptionSet("--BHPcommand"))
	  {
		IA.doBHPcommand = true;
		IA.BHPcommand_filename = (String)options.getOptionValues("--BHPcommand")[0];
		// IA.BHPcommandstring = (String)options.getOptionValues("--BHPcommand")[1];
	  }

	// Prior traces
	if (options.isOptionSet("--make-prior-traces"))
	  {
		IA.makePriors = true;

		// property words
		String s = ((String)options.getOptionValues("--make-prior-traces")[0]);
		IA.priorKeywords = s.split(",");
		if (IA.priorKeywords.length<1) 
		  {
			IA.EH.die_err("--make-prior-traces specifies no properties");
			options.Usage();
		  }
	  }


	// Massaging
	if (options.isOptionSet("--massage-analyse") || options.isOptionSet("--massage-analyse-ascii"))
	  {
		// keyword for parsing
		String optionKeyword;
		IA.massageAnalyse = true;
		if (options.isOptionSet("--massage-analyse-ascii")) 
		  {
			IA.massageAnalyseAscii = true;
			optionKeyword = "--massage-analyse-ascii";
		  } 
		else
		  {
			optionKeyword = "--massage-analyse";
		  }

		// header words
		String s = ((String)options.getOptionValues(optionKeyword)[0]);
		IA.massagingHeaderWords = s.split(",");

		// layer numbers
		s = ((String)options.getOptionValues(optionKeyword)[1]);
		String[] words = s.split(",");
		if (words.length<1) 
		  {
			IA.EH.die_err("--massage-analyse specifies no layer numbers");
			options.Usage();
		  }
		int[] analyse_layer_number = new int[words.length];
		int i;
		for (i=0;i<words.length;i++)
		  {
			try 
			  {
				analyse_layer_number[i] = Integer.parseInt(words[i]);
			  }
			catch (NumberFormatException e)
			  {
				IA.EH.warn("analysing layer numbers must be numbers only");
				options.Usage();
			  }
		  }

		// properties
		s = ((String)options.getOptionValues(optionKeyword)[2]);
		words = s.split(",");
		if (words.length<1) 
		  {
			IA.EH.die_err("--massage-analyse specifies no properties");
			options.Usage();
		  }
		IA.massagingKeywords = new ArrayList(words.length);
		for (int j=0;j<analyse_layer_number.length;j++)
		  {
			for (i=0;i<words.length;i++)
			  {
				IA.massagingKeywords.add(new AnalysingProperty(words[i],analyse_layer_number[j]));
			  }
		  }

		// output file
		IA.massagingOutputFileName = ((String)options.getOptionValues(optionKeyword)[3]);
	  }

	// test option validity
	if (IA.propertySUMap || IA.propertyMap || IA.propertyMap2 || IA.makePropertyAsciiSummary || IA.propertyMap2gen)
	  {
		if (IA.propertyMap) 
		  {
			IA.property_name1 = ((String)options.getOptionValues("-p-ascii")[0]);
		  }
		else if (IA.propertySUMap) 
		  {
			IA.property_name1 = ((String)options.getOptionValues("-p")[0]);
		  }
		else if (IA.makePropertyAsciiSummary) 
		  {
			IA.property_name1 = ((String)options.getOptionValues("-p-ascii-prop")[0]);
		  }
		else if (IA.propertyMap2)
		  {
			IA.property_name1 = ((String)options.getOptionValues("-p-ascii2")[0]);
			IA.property_name2 = ((String)options.getOptionValues("-p-ascii2")[1]);
		  }
		else if (IA.propertyMap2gen) 
		  {
			IA.property_name1 = ((String)options.getOptionValues("-p-ascii2-gen")[0]);
			IA.property_name2 = ((String)options.getOptionValues("-p-ascii2-gen")[1]);
		  }

		// checks
		if (IA.completeWords.indexOf(IA.property_name1) == -1)
		  {
			IA.EH.warn("property " + IA.property_name1 + " isn't a recognised keyword "+
					   (IA.gotFirstTrace ? " for this data set" : "") + ". Conceivable are:");
			IA.showAvailableProperties();
			IA.EH.die_err("");
		  }
		if ((IA.propertyMap2 || IA.propertyMap2gen) && IA.completeWords.indexOf(IA.property_name2) == -1)
		  {
			IA.EH.warn("property " + IA.property_name2 + " isn't a recognised keyword "+
					   (IA.gotFirstTrace ? " for this data set" : "") + ". Conceivable are:");
			IA.showAvailableProperties();
			IA.EH.die_err("");
		  }


		try
		  {
			if (IA.propertyMap) 
			  {
				IA.layerNumber = Integer.parseInt(((String)options.getOptionValues("-p-ascii")[1]));
			  }
			else if (IA.propertySUMap) 
			  {
				IA.layerNumber = Integer.parseInt(((String)options.getOptionValues("-p")[1]));
			  }
			else if (IA.propertyMap2) 
			  {
				IA.layerNumber = Integer.parseInt(((String)options.getOptionValues("-p-ascii2")[2]));
			  }
			if (IA.propertyMap2gen) 
			  {
				IA.layerNumber = Integer.parseInt(((String)options.getOptionValues("-p-ascii2-gen")[2]));
				IA.layerNumber2 = Integer.parseInt(((String)options.getOptionValues("-p-ascii2-gen")[3]));
			  }
		  }
		catch (NumberFormatException e)
		  {
			IA.EH.warn("Layer number must be a positive integer");
			options.Usage();
		  }
	  }

	// (2...) other checks...

	IA.runInversionAnalyser();

  }

  // make clean start and end range times: must be divisible by dt
  private void cleanSyntheticRange()
  {
	double EPS=1e-7; // floating point machine_epsilon, 
	// since tSyntheticStart might have been mangled in parsing
	tSyntheticStart = Math.floor(tSyntheticStart*(1+EPS)/delta_t)*delta_t;
	tSyntheticEnd = Math.floor(tSyntheticEnd*(1+EPS)/delta_t)*delta_t;
  }


  private void showAvailableProperties()
  {
	EH.warn(""+propertyWords);
  } 


  protected boolean runInversionAnalyser()
  {
	setUpTraceOutputs();
	// setUpTraceInputs();
	traceFilter.setUpTraceInputFilter(tr);
	checkHeaderWords();
	setUpAuxiliaryOutputs();
	assessIfFilterIsExclusiveOperation();
	initialiseGraphics();
	try
	  {
		loopOverTraces();
	  }
	catch (Exception ex)
	  {
		ex.printStackTrace();
		return false;
	  }
	finalizeTraceInputs();
	finalizeOutput();
	finalizeAuxiliaryOutputs();
	showGraphics();
	return true;
  }

  private void initialiseTraceCacheAndGetFormat()
  {
	// Open input file
	try 
	  {
		if (realisation_filename==null)
		  {
			// use stdin. This may have been opened in options parsing sanity checks
			if (inversion_in_strm==null)
			  {
				inversion_in_strm = (FileInputStream) new 
				  FileInputStream(FileDescriptor.in);
			  }
		  }
		else
		  {
			// use given file
			inversion_in_strm = (FileInputStream) new 
			  FileInputStream(realisation_filename);
		  } 
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open stdin or " + realisation_filename);
	  }

	try
	  {
		boolean gotTrace = (sureader.fgettr(inversion_in_strm,tr)!=0);
		if (gotTrace)
		  {
			gotFirstTrace = true;

			// get format stuff and initialise format
			int format_version  = (int)tr.get_short(ModelFormat.VERSION_NUMBER_HEADER_FIELD);
			int Nstacks  = (int)tr.get_short(ModelFormat.NUMBER_OF_STACKS_HEADER_FIELD);
			if (Nstacks!=0)
			  {
				EH.die_err("Realisations file "+realisation_filename
						   + " has (for nstacks) nhs!=0. Probably a delivery4D file."
						   + " can't parse these with legacy deliveryAnalyser");				
			  }
			setUpModelFormats(format_version);
		  }
		else
		  {
			EH.die_err("Failed to read format from first trace");
		  }
	  }
	catch (Exception ex)
	  {
		ex.printStackTrace();
	  }
  }

  private void loopOverTraces()
  {
	boolean gotTrace,skipTasks,detectedEndOfLastRealisationGroup;
	int nRealisationsPassedByFilters=0;
	int rN,last_rN=-1;

	// Loop over traces in each file
	for (trace_number=0;;trace_number++)
	  {
		// get the trace
		gotTrace = getTraceFromCache();

		if (doBHPcommand && trace_number==0) handleBHPcommand();

		// Filters & tracking of accepted realisations 
		skipTasks = false;
		detectedEndOfLastRealisationGroup = false;
		if (gotTrace)
		  {
			// filter by property if required
			if (applyFilter && !filterOutThisTrace()) skipTasks = true;
			// filter by header words if required
			if (!skipTasks && traceFilter.skipTraceIfRequired(tr)) skipTasks = true;
			rN = getRealisationNumber(tr);
			if (rN <= last_rN) detectedEndOfLastRealisationGroup = true;  // e.g. mark hdr word reset
			last_rN = rN;
			if (filterPad && trace_number==0) tr_out.copyHeader(tr); // store very first 1st header
		  }
		else
		  {
			detectedEndOfLastRealisationGroup = true; // ran out of traces
		  }
		if (detectedEndOfLastRealisationGroup)
		  {
			if (filterPad)
			  {
				if (nRealisationsPassedByFilters==0)  // NB filterPad will be exclusive program operation
				  {
					dumpNullTrace();
				  }
				tr_out.copyHeader(tr); // store header at start of group
			  }
			nRealisationsPassedByFilters = 0;
		  }

		if (skipTasks) continue;

		// here, an exclusive list of things might want to do on each trace
		if (massageAnalyse)
		  {
			generateMomentAnalysisforMassaging(!gotTrace);
		  }
		else if (makeFullSummaryTraces || makeFullAsciiSummary || makeTimesAsciiSummary || makePropertyAsciiSummary)
		  {
			if (gotTrace) makeFullSummaryTraces();
		  }
		else if (makeFullMeanSummaryTraces || makeFullStdDevSummaryTraces || makeFullQuantileSummaryTraces)
		  {
			makeSuiteOfStatisticalMomentMaps(completeWords,!gotTrace);
		  }
		else if ((layerCake || layerCakeOfFluids || layerCakeOfNG || layerCakeOfZ_eff || layerCakeOfProp) && gotTrace)
		  {
			makeLayerCakeTraces();
			// can also do these simultaneously!
			if (graphics && makeSyntheticSeismic) makeSyntheticSeismics();
			if (graphics && showHorizonsAndSeismic) showHorizonsAndSeismic();
		  }
		else if (showHorizonsAndSeismic && gotTrace)
		  {
			showHorizonsAndSeismic();
		  }
		else if (fluidContentMap && gotTrace)
		  {
			makeFluidContentTraces();
		  }
		else if (computeMean || computeStddev || computeQuantile || computeHistogram || 
				 computeRafteryThinning || computeFluidProb)
		  {
			makeStatisticalMomentMaps(computeMean,computeStddev, computeQuantile, 
									  computeHistogram, computeRafteryThinning, 
									  computeFluidProb, !gotTrace);
		  }
		else if ((propertySUMap || propertyMap || propertyMap2 || propertyMap2gen) && gotTrace)
		  {
			makePropertyMaps();
		  }
		else if (makePriors)
		  {
			makePriorTraces(!gotTrace);
		  }
		else if (makeSyntheticSeismic && gotTrace)
		  {
			makeSyntheticSeismics();
		  }
		else if (performTimeToDepthRemapping && gotTrace)
		  {
			performTimeToDepths();
		  }
		else if (gotTrace)
		  {
			// default action; reproduce traces that passed though filters...
			writeOutAcceptedTraces();
		  }
		if (!gotTrace) break;
		// Count total number of traces accepted
		accepted_trace_number++;
		// Count number traces per location
		nRealisationsPassedByFilters++;
	  }
  }

  // compute here whether --filter or --filter-pad is the only operation performed by the program
  private void assessIfFilterIsExclusiveOperation()
  {
	boolean isFilterExclusively = false;

	if (applyFilter) 
	  {
		isFilterExclusively = true;

		// list of functionality as above in loopOverTraces(). Add new functionality here as well in tests...
		if (massageAnalyse
			|| makeFullSummaryTraces || makeFullAsciiSummary || makeTimesAsciiSummary || makePropertyAsciiSummary
			|| makeFullMeanSummaryTraces || makeFullStdDevSummaryTraces || makeFullQuantileSummaryTraces
			|| layerCake || layerCakeOfFluids || layerCakeOfNG || layerCakeOfZ_eff || layerCakeOfProp
			|| fluidContentMap
			|| computeMean || computeStddev || computeQuantile || computeHistogram || computeRafteryThinning || computeFluidProb
			|| propertySUMap || propertyMap || propertyMap2 || propertyMap2gen
			|| makeSyntheticSeismic
			|| performTimeToDepthRemapping 
			)
		  {
			isFilterExclusively = false;
		  }
		// filter-pad must be exclusive
		if (filterPad && !isFilterExclusively)
		  {
			EH.die_err("--filter-pad option must be exclusive (i.e. only) operation of Analyser. "
					   +"Pipe into separate analyser operation after filter-pad!");
		  }
	  }
	return;
  }

  // use tr_out as wkspace to dump nulls
  private void dumpNullTrace()
  {
	int k,ns;

	// fill trace with -999s
	ns = (int)tr.get_short("ns");
	for (k=0; k<ns;k++) tr_out.set_data_value((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER,k);
	// set realisation number to 1
	setRealisationNumber(tr_out,1);
	// write trace out
	suwriter.fputtr(out_strm, tr_out);
  }


  private boolean isTraceNull(float[] values)
  {
	if ((values.length>0) && (values[0] <= (float)SeismicInverter.MODEL_TRACE_BAD_NUMBER))
	  {
		return true;
	  }
	return false;
  }

  private boolean getTraceFromCache()
  {
	boolean gotTrace;

	if (!gotFirstTrace) 
	  {
		return (sureader.fgettr(inversion_in_strm,tr)!=0);
	  }
	else
	  {
		// have already secured first trace. Reset this flag,
		// and continue on, so next time we get the next trace
		gotFirstTrace = false;
		return true;
	  }
  }

  /**
	 Hook to write external command script: typically a command that requires
	 the string "properties=A,B,C,.." appended to something like BHPwrite.
	 This method appends "properties=A,B,C,.." to BHPcommandstring, writes
	 it to BHPcommand_filename, and chmods this to "ugo+rx".
	 A,B,C etc are the properties formed from completeWords, but the
	 last, t_base, is omitted.
  */
  void handleBHPcommand()
  {
	int c=0;
	PrintWriter out = null;

	setNumberOfLayers();

	try
	  {
		out = new PrintWriter(new FileWriter(BHPcommand_filename));
	  }
	catch (IOException e) 
	  {
		EH.die_err("Can't open file "+BHPcommand_filename); 
	  }
	
	// Append the comma separated list of properties
	StringBuffer property_string = new StringBuffer("properties=");
	
	// Assemble the string. Omit the t and t_base entry at the end.
	for (int i=0;i<completeWords.size()-2;i++) 
	  {
		// Skip thinning if making summary traces
		if (((String)completeWords.get(i)).equals(THINNING_KEYWORD) && makeFullSummaryTraces) continue; 
		property_string.append((String)completeWords.get(i));
		// since we omit t, t_base, and THINNING, we are 3 short...
		if (i<completeWords.size()-3) property_string.append(",");  // FIXME is this ok when thinning omitted?
		c++; // count props
	  }
	// add "t" if required, i.e. when we run any of these three modes
	if (makeFullMeanSummaryTraces || makeFullStdDevSummaryTraces || makeFullQuantileSummaryTraces)
	  {
		// append "t",
		property_string.append(",t");
		c++;
	  }
	out.printf(BHPcommandstring+property_string+" %n");
	out.close();
	EH.warn("Dumped keywords for "+c+" properties");
	// Not needed
	// ScriptGenerator.chmod(BHPcommand_filename,"ugo+rx",EH);
  }
  
  void showHorizonsAndSeismic()
  {	
	int i,rN;
	float t,Xe;

	// Graphing stuff
	GraphingData gd=null;
	AnPlotInfo api=null;

	if (graphics) 
	  {
		gd = (GraphingData)GraphingDataMap.get(5);
	  }
	// read data  from trace
	float values[] = tr.get_data_array();

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		// times are last n_layer values. Compute n_layers
		setNumberOfLayers();
		refIndx_t = MF.reducedPropertyList.indexOf("t") * n_layers;
		// model trace times are in ms
		if (!layerRangeSpecified)
		  {
			// auto-guess these from first trace
			t_start = 0.001f*values[refIndx_t];
			t_end = 0.001f*values[refIndx_t + n_layers];
		  }
		EH.warn("start and end times = "+t_start+", "+t_end);
		if (showOverLayerCake) 
		  {
			// Set up index of NG
			propFetcher1 = new TracePropertyFetcher(this);
			setUpLayerColors(false,false);
		  }
		// read, save seismic, also stuff max/min values into SGI[0,1] for graphics
		// read an extra second of data above and below. Rely on graphics clipping to suppress it		
		SGI = gCaptureWindowedSeismicTraces(gd,t_start-1.0,t_end+1.0); 
	  }
	if (graphics) 
	  {
		// track where we are in realisation stream
		rN = getRealisationNumber(tr);
		if (rN <= lastRealisationNumberHeaderField) accepted_tracegroup_number++;
		lastRealisationNumberHeaderField = rN;

		// make little line segments at horizon picks
		Xe = max(-SGI.Xmin,SGI.Xmax);
		if (!showOverLayerCake)
		  {
			// make little line segments at horizon picks
			for(i=0;i<=n_layers;i++)
			  {
				api = new AnPlotInfo(false);
				DoubleArrayList X = api.x;
				DoubleArrayList Y = api.y;
				t = 0.001f*values[refIndx_t + i];
				// horizontal line segments for picks, across seismic span....
				X.add(-Xe);       // -ve value, largest amplitude
				Y.add(t);            // Horizon
				X.add(Xe);       // +ve value, largest amplitude
				Y.add(t);
				api.offsetNumber = accepted_tracegroup_number;
				api.rNumber = rN;
				api.trim();
				gd.plotList.add(api);
			  }
		  }
		else
		  {
			// Block plot underneath wiggles
			// Alloc a block-plot suite per column...
			api = new AnPlotInfo(true);
			DoubleArrayList X,X1,Y,Y1;
			X = api.x; X1 = api.x1; Y = api.y; Y1 = api.y1;
			api.offsetNumber = accepted_tracegroup_number;
			api.rNumber = rN;
			for(i=0;i<n_layers;i++)
			  {
				t = 0.001f*values[refIndx_t + i];            // layer top
				float t_b = 0.001f*values[refIndx_t + i+1];  // layer base
				// rectangles for layers, across seismic span....
				X.add(-Xe);       // -ve value, largest amplitude
				Y.add(t);         // Horizon
				X1.add(Xe);       // +ve value, largest amplitude
				Y1.add(t_b);      // base
				// colors
				//propFetcher1.setIndex("NG",i+1,false);			
				//float dcolor = (float) propFetcher1.getRequiredQuantityFromTrace(values,"Z_eff");
				float dcolor = ((float)i)/(n_layers);
				api.cList.add(dcolor);
				api.cgList.add(-1); // int gcolor=-1
			  }
			api.trim();
			gd.plotList.add(api);
		  }
	  }
  }

  void makeLayerCakeTraces()
  {
	boolean colourByHydroOnly = false;
	// NG colouring things
	boolean colourByNGOnly = false;
	boolean colourByZ_effOnly = false;
	boolean colourByPropOnly = false;

	// default t resolution (Got a better idea?);
	int i,ii,j,jmax,k,layer_indx,gcolor=-1;
	float t,t_layer_top,t_layer_bot,t_base,color,dcolor=0;
	float t_sublayer_top,t_sublayer_bot;
	float z_layer_top=0,z_layer_bot=0;
	float z_sublayer_top,z_sublayer_bot=0,depth=0,den;
	boolean haveContacts=false;

	// oil colouring things
	float oilPresent=0;
	boolean addOilColouring = !true;
	// contacts stuff
	int n_internal_contacts;
	short fluid_type=0;
	short sub_fluid_type=0;

	// Graphing stuff
	GraphingData gd=null;
	AnPlotInfo api=null;

	if (graphics) 
	  {
		gd = (GraphingData)GraphingDataMap.get(2);
	  }
	
	// read data  from trace
	float values[] = tr.get_data_array();

	if (layerCakeOfFluids) 
	  {
		// hydro colouring (layer independent)
		colourByHydroOnly = true;
	  }
	else if (layerCakeOfNG)
	  {
		// NG coloring
		colourByNGOnly = true;
	  }
	else if (layerCakeOfZ_eff)
	  {
		// Z_eff colouring
		colourByZ_effOnly = true;
	  }
	else if (layerCakeOfProp)
	  {
		colourByPropOnly = true;
	  }

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		if (graphics) 
		  {
			api = new AnPlotInfo(true);
			gd.plotList.add(api);
		  }
		// times are last n_layer values. Compute n_layers
		setNumberOfLayers();
		refIndx_t = MF.reducedPropertyList.indexOf("t") * n_layers;
		// view at quarter delta_t for detail
		delta_t = LAYER_CAKE_DEFAULT_DT; 
 		float dt_in = 1.0e-6f*tr.get_short("dt");
		if (dt_in <= (2.0*delta_t)) 
		  {
			delta_t = 0.25f*dt_in;
		  }


		if (colourByHydroOnly)
		  {
			// Set up index of NG
			propFetcher1 = new TracePropertyFetcher(this);
			propFetcher1.setIndex("NG",1,false);
		  }
		if (colourByNGOnly)
		  {
			refIndx0 = MF.reducedPropertyList.indexOf("NG") * n_layers;
		  } 
		else if (colourByZ_effOnly)
		  {
			refIndx0 = MF.reducedPropertyList.indexOf("Z_eff") * n_layers;
		  } 
		else if (layerCakeOfProp)
		  {
			propFetcher1 = new TracePropertyFetcher(this);
		  }
		else
		  {
			refIndx0 = MF.reducedPropertyList.indexOf("fluid_type") * n_layers;
		  }
		// for fluid contacts
		refIndx1 = MF.reducedPropertyList.indexOf("t_fc");
		if (refIndx1 >= 0) // have contacts
		  {
			refIndx1 *= n_layers;
			EH.warn("Showing contacts if present");
		  }
		if (layerCakesInDepth)
		  {
			// Set up index of depths
			propFetcher0 = new TracePropertyFetcher(this);
			propFetcher0.setIndex("d",1,false);
			// Set up index of vp_eff in last layer
			propFetcher2 = new TracePropertyFetcher(this);
			propFetcher2.setIndex("vp_eff",n_layers,false);
		  }

		// model trace times are in ms
		if (!layerRangeSpecified)
		  {
			// auto-guess these from first trace
			t_start = 0.001f*values[refIndx_t];
			t_end = 0.001f*values[refIndx_t + n_layers];
			if (layerCakesInDepth)
			  {
				// top
				depth = 0.001f*values[propFetcher0.refIndx0];
				d_start = depth;
				// do base via calculation depth_base =d_N+(1/2)0.001(t_base-t_N)*vp_{eff,N}
				// (times in ms and are two-way....)
				depth = values[propFetcher0.refIndx0+n_layers-1];
				depth += 0.0005*(values[refIndx_t+n_layers] - values[refIndx_t+n_layers-1])
				  * values[propFetcher2.refIndx0];
				depth *= 0.001;
				d_end = depth;
			  }
		  }
		EH.warn("start and end times = "+t_start+", "+t_end);
 		if (layerCakesInDepth)
		  {
			new_ns = (int)Math.ceil((d_end-d_start)/delta_t);
		  }
		else
		  {
			new_ns = (int)Math.ceil((t_end-t_start)/delta_t);
		  }
		EH.warn("fixed no. of vertical samples chosen = "+new_ns);
		// if a depth model, the ns can be all over the place. Use 1000 as a crisis default...
		if (new_ns>1000)
		  {
			if (layerCakesInDepth)
			  {
				delta_t = (float)((d_end-d_start)/1000.0);
				new_ns = (int)Math.ceil((d_end-d_start)/delta_t);
			  }
			else
			  {
				delta_t = (float)((t_end-t_start)/1000.0);
				new_ns = (int)Math.ceil((t_end-t_start)/delta_t);
			  }
			EH.warn("fixed no. of samples altered : "+new_ns);
		  }
		setUpLayerColors(colourByHydroOnly, colourByNGOnly);
	  } // first trace
	// ok, now make a fake trace by stepping across these times

	t_base = 0.001f*values[refIndx_t+n_layers];
	// set no. of samples correctly. make this the same for all traces
	if (!graphics)
	  {
		tr.set("ns",(short)new_ns);
		t = layerCakesInDepth ? ((float)d_start) : ((float)t_start);
		tr.set("f1",t);
		tr.set("dt",(short)(delta_t*1e6));
		tr.set("delrt",(short)Math.rint((t*1e3)));
		tr.set("d1",0f);
		tr.set("d2",0f);
		// clean data on output trace
		tr.zero_data_array();
	  }


	if (graphics) 
	  {
		api = gd.plotList.get(0);
		if (layerCakesInDepth) 
		  {
			// these are available from above
			api.Ymin = d_start;
			api.Ymax = d_end;
		  }
		else
		  {
			// take supplied times
			api.Ymin = t_start;
			api.Ymax = t_end;
		  }
	  }

	// here, i==-1 = top padding,
	// i==n_layers = lower padding
	for (i=-1;i<=n_layers;i++)
	  {
		t_layer_top =  (i==-1)? t_start : 0.001f*values[refIndx_t+i];
		t_layer_bot =  (i==n_layers)? t_end : 0.001f*values[refIndx_t+i+1];
		if (layerCakesInDepth)
		  {
			// Get layer top depth. Work in km or kilo-ft
			if (i==-1)
			  {
				depth = d_start;
			  }
			else if (i<n_layers) // depth of layer top
			  {
				depth = 0.001f*values[propFetcher0.refIndx0+i];
			  }
			else if (i==n_layers) // require here base depth of last layer
			  {
				// do base via calculation depth_base =d_N+(1/2)0.001(t_base-t_N)*vp_{eff,N}
				// (times in ms and are two-way....)
				depth = values[propFetcher0.refIndx0+n_layers-1];
				depth += 0.0005*(values[refIndx_t+n_layers] - values[refIndx_t+n_layers-1])
				  * values[propFetcher2.refIndx0];
				depth *= 0.001;
			  }
			z_layer_top = depth;

			// Get layer bottom depths. Work in km or kilo-ft
			if (i<(n_layers-1)) // top depth of layer below
			  {
				depth = 0.001f*values[propFetcher0.refIndx0+i+1];
			  }
			else if (i==(n_layers-1)) // base depth of last layer
			  {
				// do base via calculation depth_base =d_N+(1/2)0.001(t_base-t_N)*vp_{eff,N}
				// (times in ms and are two-way....)
				depth = values[propFetcher0.refIndx0+n_layers-1];
				depth += 0.0005*(values[refIndx_t+n_layers] - values[refIndx_t+n_layers-1])
				  * values[propFetcher2.refIndx0];
				depth *= 0.001;
			  }
			else // padding depth
			  {
				depth = t_end;
			  }
			z_layer_bot = depth;
		  }
		if ((!layerCakesInDepth && (t_layer_bot > t_layer_top)) || (layerCakesInDepth && (z_layer_bot > z_layer_top)))
		  {
			// layer indexed i has finite thickness. Try to depict it
			layer_indx = i;
			n_internal_contacts=0;
			if (i==-1) color = layer_color_step_values[0];
			else if (i==n_layers) color = layer_color_step_values[n_layers+1];
			else
			  {
				fluid_type = (short)values[refIndx0+layer_indx];
				n_internal_contacts = (fluid_type >> 10);
				// EH.warn("N contact = "+n_internal_contacts);
				color = layer_color_step_values[layer_indx+1];

				// decorate with oil present 
				if (!(refIndx1>=0) && addOilColouring && (values[refIndx0+layer_indx]>0)) color += 0.5;

				// colour by NG or Z_eff, or prop
				if (colourByNGOnly || colourByZ_effOnly)
				  {
					color = values[refIndx0+layer_indx];	
				  }
				else if (colourByPropOnly)
				  { 
					propFetcher1.setIndex(layerCakeProp,layer_indx+1,false);
					color = (float) propFetcher1.getRequiredQuantityFromTrace(values,layerCakeProp);
				  }
				else if (colourByHydroOnly)
				  {
					// layer-trend is 0 : detail added below...
					color=0;
				  }
			  }

			t_sublayer_top = t_layer_top;
			t_sublayer_bot = t_layer_bot;
			z_sublayer_top = z_layer_top;

			// internal model detail from contacts...
			for (k=0;k<n_layers;k++)
			  {
				dcolor=0;
				// get all possible sublayers
				// internal fluid contacts?
				if (refIndx1>=0)
				  {
					// scan for internal reflections
					for (ii=0;ii<n_layers;ii++)
					  {
						t = 0.001f*values[refIndx1+ii];
						if ((t > t_sublayer_top) && (t < t_layer_bot))
						  {
							// contact inside layer
							t_sublayer_bot = t;
							break;
						  }
					  }
					if (k>0 && t_sublayer_top==t_sublayer_bot)
					  {
						// fell through the last internal contact. set..
						t_sublayer_bot = t_layer_bot;
					  }
					if (n_internal_contacts>0)
					  {
						// grade internal colors by fluid
						sub_fluid_type = (short)(3 & ((fluid_type >> (2*(k+1)))));
						dcolor = sub_fluid_type; // =0,1,2,3 fluid type of kth sublayer by bit extraction 
						if (graphics && colourByHydroOnly)
						  {
							gcolor = canvasColor(sub_fluid_type);
						  }
						dcolor = (colourByHydroOnly ? 1.0f : 0.25f)*dcolor; // map to 0,0.25,0.5,0.75 or 0,1,2,3
					  }
					else
					  {
						if ((i==-1) || (i==n_layers))
						  {
							dcolor=0; // model padding	
						  }
						else
						  {
							sub_fluid_type = (short)(3 & fluid_type); // mask out all but lowest 2 bits
							if (graphics && colourByHydroOnly)
							  {
								gcolor = canvasColor(sub_fluid_type);
							  }
							// default color adjustment by majority fluid class
							dcolor = (colourByHydroOnly ? 1.0f : 0.25f) * sub_fluid_type; 
						  }
					  }
					if (colourByHydroOnly && (values[propFetcher1.refIndx0 + layer_indx]==0.0))
					  {
						// override hydrocarbon by "shale" indicator if NG==0 (propFetcher1 gets NG)
						dcolor=0;
						if (graphics) gcolor = 7; // brown, brine
					  }
					else if (colourByPropOnly)
					  {
						// No color morphing by fluid type....
						dcolor=0;
					  }
				  }
				dcolor += color;
				if (layerCakesInDepth)
				  {
					// Linearly interpolate depths from sublayer times:					
					den = (t_layer_bot - t_layer_top);
					den = (den>0.0)? den : 1.0f;
					z_sublayer_bot = z_layer_top + (z_layer_bot - z_layer_top)*(t_sublayer_bot - t_layer_top)/den;

					// fill trace with colour for samples in depth range
					j = Math.max(0,(int)Math.ceil((z_sublayer_top - d_start)/delta_t));
					t = j*delta_t + d_start;

					jmax = Math.min(new_ns-1, Math.max(0,(int)Math.floor((z_sublayer_bot - d_start)/delta_t)));
					// for (; t<=z_sublayer_bot && j<new_ns; t+=delta_t, j++) 
					if (!graphics)
					  {
						for (; j<=jmax; j++) tr.set_data_value(dcolor,j);
					  }
					else
					  {
						if (i>=0 && i<n_layers) // depict internals only
						  {
							api.x.add(accepted_trace_number);
							api.x1.add(accepted_trace_number+1);
							api.y.add(z_sublayer_top);
							api.y1.add(z_sublayer_bot);
							api.cList.add(dcolor);
							api.cgList.add(gcolor);
						  }
					  }
				  }
				else
				  {
					// fill trace with colour for samples in time range
					j = Math.max(0,(int)Math.ceil((t_sublayer_top - t_start)/delta_t));
					jmax = Math.min(new_ns-1, Math.max(0,(int)Math.floor((t_sublayer_bot - t_start)/delta_t)));
					t = j*delta_t + t_start;
					// EH.warn("DT: "+t_sublayer_top+" "+t_sublayer_bot+" col="+dcolor);
					// for (; t<=t_sublayer_bot && j<new_ns; t+=delta_t, j++) tr.set_data_value(dcolor,j);
					if (!graphics)
					  {
						for (; j<=jmax; j++) tr.set_data_value(dcolor,j);
					  }
					else
					  {
						if (i>=0 && i<n_layers) // depict internals only
						  {
							api.x.add(accepted_trace_number);
							api.x1.add(accepted_trace_number+1);
							api.y.add(t_sublayer_top);
							api.y1.add(t_sublayer_bot);
							api.cList.add(dcolor);
							api.cgList.add(gcolor);
						  }
					  }
				  }			  
				if (refIndx1>=0)
				  {
					t_sublayer_top = t_sublayer_bot;
					z_sublayer_top = z_sublayer_bot;
				  }
				else break;
			  }
		  }
	  }
	// add black trace if required
	if (layersMakeFirstTraceInRealisationsBlack)
	  {
		if ((int)tr.get_short(MF.realisationNumberHeaderField)==1)
		  {
			if (!graphics)
			  {
				for (j=0; j<new_ns; j++) tr.set_data_value(largest_layer_color+1,j);
			  }
		  }
	  }
	// mask nulls
	if (isTraceNull(values) && !graphics)
	  {
		tr.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);
	  }
	// now write trace out
	if (!graphics) suwriter.fputtr(out_strm, tr);
  }
  
  private void setUpLayerColors(boolean colourByHydroOnly,boolean colourByNGOnly)
  {
	layer_color_step_values = new float[n_layers+2];
	for (int j=0;j<=n_layers+1;j++)
	  {
		if (colourByHydroOnly || colourByNGOnly) // || colourByPropOnly) 
		  {
			// upper and lower padding = -1
			layer_color_step_values[j] = (j==0 || j==(n_layers+1)) ? -1f : 0f;
		  }
		else layer_color_step_values[j] = (float)j; // integer steps
		// optional colors
		if (layerColorArray!=null && j<layerColorArray.size() && layerColorArray.get(j)!=null) 
		  {
			layer_color_step_values[j]=(float)((Integer)layerColorArray.get(j)).intValue();
		  }
		largest_layer_color=Math.max(layer_color_step_values[j],largest_layer_color);
	  }
  }
  /**
	 fluid colouring index suitable for CanvasPlotter.defaultPalette[]
	 array to render fluid types
  */
  private int canvasColor(short fluid_type)
  {
	int gcolor=0;
	switch (fluid_type)
	  {
	  case 0:
		gcolor = 4; break; // blue, brine
	  case 1:
		gcolor = 2; break; // red, oil
	  case 2:
		gcolor = 5; break; // magenta, lsg
	  case 3:
		gcolor = 9; break; // orange, gas
	  }
	return gcolor;
  }

  void makeFluidContentTraces()
  {
	int j;

	// read data  from trace
	float values[] = tr.get_data_array();

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		// Compute n_layers
		setNumberOfLayers();
		propFetcher0 = new TracePropertyFetcher(this);
		propFetcher0.refIndx0 = MF.reducedPropertyList.indexOf("fluid_type") * n_layers;
		new_ns = n_layers;
	  }
	// ok, now make a fake trace by mapping fluid-type to the samples
	// set no. of samples correctly. make this the same for all traces
	tr.set("ns",(short)new_ns);
	tr.set("dt",(short)1e3);  // time scale is layer number, x 1e-3 (can't set to seconds)
	tr.set("delrt",(short)0);
	// clean data on output trace
	tr.zero_data_array();
	for (j=0; j<n_layers; j++)
	  {
		tr.set_data_value(values[propFetcher0.refIndx0+j],j);
	  }
	// mask nulls
	if (isTraceNull(values))
	  {
		tr.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);
	  }
	// now write trace out
	suwriter.fputtr(out_strm, tr);
  }

  /**
	 here we dump a master trace that has everything we could possibly want in it
  */
  private void makeFullSummaryTraces()
  {
	int i,j,k;
	double value;
	String word;
	boolean isNull;
	// Simply write out stuff of interest in a model-format trace
	// i.e. consecutive blocks corresponding to the quantities keyed by the
	// completeWords list.

	// read data from trace
	float values[] = tr.get_data_array();
	
	// detect nulls
	isNull = isTraceNull(values);

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		if (makePropertyAsciiSummary) 
		  {
			property_name2 = new String(property_name1);  // make safe copy -> property_name1 is abused below
		  }
		// Compute n_layers
		setNumberOfLayers();
		propFetcher0 = new TracePropertyFetcher(this);
		// count size omitting t_base
		Ns = (completeWords.size()-1)*n_layers+1;
		if (completeWords.contains(THINNING_KEYWORD)) Ns-=n_layers;  // no thinning block
		// headers for ascii output
		if (makeFullAsciiSummary || makeTimesAsciiSummary || makePropertyAsciiSummary)
		  {
			columnTitles = new ArrayList();
			for (i=0,j=0;i<completeWords.size();i++)
			  {
				// skip THINNING
				if (((String)completeWords.get(i)).equals(THINNING_KEYWORD)) continue;
		        // other words....
				for (k=1; k<=n_layers && j<Ns; k++,j++)
				  { // column titles are magic words, with _layer-number appended
					word = (String)completeWords.get(i);
					StringBuffer title = new StringBuffer(word).append("_").append(k);
					// ensure columns are at least 10 characters wide
					for (int l=title.length();l<=10;l++) title.insert(0,' ');
					columnTitles.add(j,title.toString());
					if (makePropertyAsciiSummary && !word.equals(property_name2)) continue; // skip non-matching props
					if (makeTimesAsciiSummary && (j<Ns-(n_layers+1))) continue;     // skip non-times
					EH.printf("%s ",title.toString());
				  }
			  }
			EH.println("");
		  }
	  }
	if (makeFullSummaryTraces) tr.set("ns",(short)Ns);
	// Loop over properties listed by completeWords
	// for interfacing with BHP stuff, use the --BHPcommand, which
	// will list the properties in completeWords in the script it creates
	for (i=0,j=0;i<completeWords.size();i++)
	  {
		property_name1 = (String) completeWords.get(i);
		// skip property if not required
		if (makePropertyAsciiSummary && !property_name1.equals(property_name2)) continue; 

		// skip thinning
		Integer mapInteger = (Integer)keyWordsMap.get(property_name1);
		if (mapInteger!=null && mapInteger.intValue()==InversionAnalyser.THINNING_KEY) continue;

		// and then loop over layers k. Bail out at j=Ns also (only one t_base entry)
		for (k=1; k<=n_layers && j<Ns; k++,j++)
		  {			
			if (makeTimesAsciiSummary && (j<Ns-(n_layers+1))) continue;
			// compute correct index of desired property
			propFetcher0.setIndex(property_name1,k,(k!=1));
			value = propFetcher0.getRequiredQuantityFromTrace(values,property_name1);
			if (isNull) value = (float)SeismicInverter.MODEL_TRACE_BAD_NUMBER;
			// output ascii if required
			if (makeFullAsciiSummary || makeTimesAsciiSummary || makePropertyAsciiSummary)
			  {
				// print columns in width of precision allowed by title: title is jth element in columnTitles
				EH.printf("%"+((String)columnTitles.get(j)).length()+"g ", value);
			  }
			else
			  {
				// stuff required property into the required float value in the output trace
				tr.set_data_value((float)value,j);
			  }
		  }
	  }
	if (makeFullAsciiSummary || makeTimesAsciiSummary || makePropertyAsciiSummary) EH.println("");

	if (makeFullSummaryTraces)
	  {
		// now write trace out
		suwriter.fputtr(out_strm, tr);
	  }
  }


  boolean filterOutThisTrace()
  {
	double value;

	// read data  from trace
	float values[] = tr.get_data_array();

	// on first trace set up the following....
	if (trace_number==0)
	  {
		// Compute n_layers
		setNumberOfLayers();
		if (filterLayerNumber<1) EH.die_err("Filter layer number must not be negative");
		if (filterLayerNumber>n_layers) EH.die_err("Filter layer number must be less than "+ (n_layers+1));
		// compute correct index of desired property
		filterPropertyFetcher0 = new TracePropertyFetcher(this);      // filter property
		filterPropertyFetcher0.setIndex(filter_property_name,filterLayerNumber,false);
	  }

	// the filter value
	value = filterPropertyFetcher0.getRequiredQuantityFromTrace(values,filter_property_name);
	switch (filter_operation) 
	  {
	  case '=':
		if (value == filter_value) return true;
		break;
	  case '<':
		if (value < filter_value) return true;	  
		break;  
	  case '>':
		if (value > filter_value) return true;	 
		break;
	  }
	// default no filtering
	return false;
  }


  void makePropertyMaps()
  {
	double value1,value2=0;
	AnPlotInfo api=null;
	GraphingData gd=null;

	if (graphics) 
	  {
		gd = (GraphingData)GraphingDataMap.get(3);
	  }
	// read data  from trace
	float values[] = tr.get_data_array();

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		if (graphics) 
		  {
			api = new AnPlotInfo();
			gd.plotList.add(api);
		  }

		// Compute n_layers
		setNumberOfLayers();
		if (layerNumber<1) EH.die_err("Layer number must not be negative");
		if (layerNumber>n_layers) EH.die_err("Layer number must be less than "+ (n_layers+1));

		// compute correct index of desired property
		propFetcher0 = new TracePropertyFetcher(this);      // first property
		propFetcher0.setIndex(property_name1,layerNumber,false);
		if (propertyMap2) 
		  {
			propFetcher1 = new TracePropertyFetcher(this);  // second property if reqd
			propFetcher1.setIndex(property_name2,layerNumber,false);
		  } 
		else if (propertyMap2gen) 
		  {
			if (layerNumber2<1) EH.die_err("Layer number (prop 2) must not be negative%n");
			if (layerNumber2>n_layers) EH.die_err("Layer number (prop 2) must be less than "+ (n_layers+1));

			propFetcher1 = new TracePropertyFetcher(this);  // second property if reqd
			propFetcher1.setIndex(property_name2,layerNumber2,false);
		  } 
	  }
	if (graphics) 
	  {
		api = gd.plotList.get(0);
	  }

	value1 = propFetcher0.getRequiredQuantityFromTrace(values,property_name1);
	if (propertySUMap)
	  {
		// Make a fake trace by mapping fluid-type to the samples
		// set no. of samples correctly. make this the same for all traces
		tr.set("ns",(short)1);

		// stuff required property in as the first float value
		tr.set_data_value((float)value1,0);

		// now write trace out
		suwriter.fputtr(out_strm, tr);
	  }
	else
	  {
		// simply write the ascii value(s) to stdout
		if (propertyMap)
		  {
			if (!graphics)
			  {
				EH.println("%10.8g",value1);
			  }
			else
			  {
				api.x.add(accepted_trace_number);
				api.y.add(value1);
			  }
		  }
		else if (propertyMap2 || propertyMap2gen)
		  {
			value2 = propFetcher1.getRequiredQuantityFromTrace(values,property_name2);
			if (!graphics)
			  {
				EH.println("%10.8g %10.8g",value1, value2);
			  }
			else
			  {
				api.x.add(value1);
				api.y.add(value2);
			  }
		  }
	  }
  } 
  
  /**
	 Gets realisation number from the trace tr, in the header word realisationNumberHeaderField.
	 The keyword realisationNumberHeaderField may point to an int, short or float.
	 It'd better be a short or int. Get this value and return it cast to an int.
  */
  private int getRealisationNumber(SuSegyFormat _tr)
  {
	int realisationNumber=1;

	switch (realisationNumberSuHeaderMapEntry.type) 
	  {
	  case SuSegyFormat.SuHeadTypes.short_type:
		realisationNumber = (int)_tr.get_short(MF.realisationNumberHeaderField);
		break;
	  case SuSegyFormat.SuHeadTypes.int_type:
		realisationNumber = _tr.get_int(MF.realisationNumberHeaderField);
		break;
	  }
	return realisationNumber;
  }

  /**
	 gets generic value from trace as a string
  */
  private String getStringEquivalent(SuSegyFormat _tr,String headerWord)
  {
	SuSegyFormat.SuHeaderEntry entry = (SuSegyFormat.SuHeaderEntry) _tr.headerMap.get(headerWord);
	if (entry.type == SuSegyFormat.SuHeadTypes.short_type)
	  {
		return Short.toString(_tr.get_short(headerWord));
	  }
	else if (entry.type == SuSegyFormat.SuHeadTypes.int_type)
	  {
		return Integer.toString(_tr.get_int(headerWord));
	  }
	else if (entry.type == SuSegyFormat.SuHeadTypes.float_type)
	  {
		return Float.toString(_tr.get_float(headerWord));
	  }
	return "";
  }


  /**
	 Sets realisation number on the trace tr, in the header word realisationNumberHeaderField.
	 The keyword realisationNumberHeaderField may point to an int, short or float.
	 It'd better be a short or int. Get this value and return it cast to an int.
  */
  private void setRealisationNumber(SuSegyFormat _tr, int value)
  {
	switch (realisationNumberSuHeaderMapEntry.type) 
	  {
	  case SuSegyFormat.SuHeadTypes.short_type:
		_tr.set(MF.realisationNumberHeaderField,(short)value);
		break;
	  case SuSegyFormat.SuHeadTypes.int_type:
		_tr.set(MF.realisationNumberHeaderField,value);
		break;
	  }
  }

  private void makeStatisticalMomentMaps(boolean computeMean, boolean computeStdDev, boolean computeQuantile, 
										 boolean computeHistogram, boolean computeRafteryThinning, 
										 boolean computeFluidProb, boolean finalise)
  {
	int i;
	double value=0.0;
	float[] values=null;
	SyntheticSeismicData requiredSyntheticSeismicData;
	GraphingData gd=null;
	AnPlotInfo api=null;

	// This new trace may be the start of a new sequence, or the end of file (finalise==true). 
	// If so complete current statistical ensemble
	if (finalise || (getRealisationNumber(tr) <= lastRealisationNumberHeaderField))
	  {
		if (dataList==null) 
		  {
			EH.println("Null data set. Probably no satisfactory traces...");
			return;
		  }
		// compute stats if not making hisotgrams
		if (!computeHistogram)
		  {
			// fluid prob of desiredFluidType in layer
			if (computeFluidProb)
			  {
				value = getFluidProb(dataList,desiredFluidType);
			  }
			// quantiles
			else if (computeQuantile)
			  {
				dataList.sort();
				value = Descriptive.quantile(dataList,quantile);
			  }
			else if (computeRafteryThinning)
			  {
				value = (double)getRafteryThinningEstimate(dataList);
			  }
			else if (computeMean)
			  {
				// means
				value = Descriptive.mean(dataList);
			  }
			else
			  {
				// Get the std dev. We need the mean
				double mean = Descriptive.mean(dataList);
				// then the variance
				value = Descriptive.sampleVariance(dataList,mean);
				// convert to std_dev
				if (value>0) value = Math.sqrt(value);
			  }
			if (realisationGroupIsNull) value = (float)SeismicInverter.MODEL_TRACE_BAD_NUMBER;
		  }

		// Form output.....
		if (!makeSyntheticSeismic)
		  {
			if (graphics) 
			  {
				gd = (GraphingData)GraphingDataMap.get(3);
			  }
			if (computeFluidProb || computeMean || computeQuantile || computeStdDev)
			  {
				// simple ascii indicator for mean, stddev or quantile
				EH.println("N(samples) = "+dataList.size()+" statistic= " +value);
			  }
			else if (computeRafteryThinning)
			  {
				// simple ascii indicator for mean, stddev or quantile
				EH.println("N(samples) = "+dataList.size()+" suggested thinning= " +value);
			  }
			else if (computeHistogram)
			  {
				// Histogram output
				dataList.sort();
				// Let's take sqrt(2N) bins (same default rule as in |stat), 
				// but no less than 4, and no more than 25.
				int Nbins = Math.min(25,Math.max(4,(int)Math.rint(2*Math.sqrt(dataList.size()))));
				double min = dataList.get(0);
				double max = dataList.get(dataList.size()-1);
				// In wierd cases, if data is a spike (all constant), 
				// we need to ensure we have a non-zero range
				if (max==min) 
				  {
					if (max>0.0) {max*=1.05; min/=1.05;}
					else if (max<0.0) {min*=1.05; max/=1.05;} 
					else
					  {
						max=1.0;
						min=-1.0;
					  }
				  }
				if (isFirstHistogram)
				  {
					histogram_max = max;
					histogram_min = min;
					isFirstHistogram = false;
				  }
				else
				  {
					max = histogram_max;
					min = histogram_min;
				  }


				// form histogram with some % extra range at either end (say 10%)
				// Ensure range is all positive too, *unless* we're examining reflection coeffs, or trace is null (-999s)
				IHistogram1D h1;
				if (realisationGroupIsNull)
				  {
					h1 = new Histogram1D(null,Nbins,SeismicInverter.MODEL_TRACE_BAD_NUMBER-1.0,
										 SeismicInverter.MODEL_TRACE_BAD_NUMBER+1.0);
				  }
				else if (property_name1.equals("R_near") || property_name1.equals("R_far")
						 || property_name1.equals("R_near_fc") || property_name1.equals("R_far_fc"))
				  {
					h1 = new Histogram1D(null,Nbins,min-0.1*(max-min),max+0.1*(max-min));
				  }
				else
				  {
					// all positive
					h1 = new Histogram1D(null,Nbins,Math.max(0,min-0.1*(max-min)),max+0.1*(max-min));
				  }
				// populate histogram
				if (!realisationGroupIsNull) 
				  {
					for (i=0; i < dataList.size(); i++) h1.fill(dataList.get(i));
				  }
				else
				  {
					for (i=0; i < dataList.size(); i++) h1.fill(SeismicInverter.MODEL_TRACE_BAD_NUMBER);
				  }
				// Generate 2 column output: bin_centre, count_in_bin.
				IAxis ax = h1.xAxis();
				int Nb = ax.bins();
				if (graphics) 
				  {
					gd = (GraphingData)GraphingDataMap.get(4);
					api = new AnPlotInfo();
					
					api.x.add(ax.binLowerEdge(0));
					api.y.add(0);
					for (i=0; i < Nb; i++) 
					  {
						api.x.add(ax.binLowerEdge(i));
						api.x.add(ax.binUpperEdge(i));
						api.y.add(h1.binHeight(i));
						api.y.add(h1.binHeight(i));
					  }
					api.x.add(ax.binUpperEdge(Nb-1));
					api.y.add(0);
					api.trim();
					gd.plotList.add(api);
				  }
				else
				  {
					// ascii stream to stdout
					EH.printf("%10.8g 0%n",ax.binLowerEdge(0));
					for (i=0; i < Nb; i++) 
					  {
						EH.printf("%10.8g %10.8g%n%10.8g %10.8g%n",ax.binLowerEdge(i), h1.binHeight(i),ax.binUpperEdge(i), h1.binHeight(i));
					  }
					EH.printf("%10.8g 0%n%n",ax.binUpperEdge(Nb-1));
				  }
			  }
		  }
		else
		  {
			// Must be computing the synthetic seismic of some `ensemble average'
			// if making a synthetic seismic corresponding to some average, ask for it here
			if (computeMean) 
			  {
				// data isn't sorted yet, we better do this for the sake of the rankInterpolation below
				// The mean is still sensible though
				dataList.sort();
			  }
			// index of element in sorted list corresponding to quantile
			int indx = (int)(Descriptive.rankInterpolated(dataList,value)-0.5);
			// EH.warn("retrieving element ...."+ indx + "array-size="+dataList.size());
			// now ensure value is an exact element extracted from the list. Rounding the quantile seems a fair idea
			value = dataList.get(indx);
			// and retrieve corresponding set of seismic data from the HashMap
			requiredSyntheticSeismicData = (SyntheticSeismicData)SyntheticSeismicDataMap.get(new Double(value));
			makeGivenSyntheticSeismic(requiredSyntheticSeismicData,realisationGroupIsNull);
			// and clean out HashMap
			SyntheticSeismicDataMap.clear();
		  }
		if (!finalise)
		  {
			// reset statistics summation and counter variables
			dataList.removeFromTo(0,dataList.size()-1);
			dataList.trimToSize();
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
		  }
		realisationGroupIsNull = false;
	  }

	// on an ordinary trace (not EOF), we must increment the running statistics
	if (!finalise)
	  {
		// read data from trace
		values = tr.get_data_array();
		// on first trace set up the following....
		if (accepted_trace_number==0)
		  {
			// Compute n_layers
			setNumberOfLayers();
			if (layerNumber<1) EH.die_err("Layer number must not be negative");
			if (layerNumber>n_layers) EH.die_err("Layer number must be less than "+ (n_layers+1));
		    
			propFetcher0 = new TracePropertyFetcher(this);      // first property
			propFetcher0.setIndex(property_name1,layerNumber,false);  // property_name1 is what we do stats on

			dataList = new DoubleArrayList();		
			realisationGroupIsNull = isTraceNull(values);

			if (makeSyntheticSeismic) 
			  {
				// we be needing to store reflection coeffs etc
				SyntheticSeismicDataMap = new HashMap();
				// need to store reflection coefficients etc for later retrieval, and
				// initialise seismic model
				SM = new OneDSeismicConvolutionalModel(this);
				SM.initialiseSeismicModel(n_layers);	
				// set time range of synthetic seismic
				SM.setSyntheticAcquisitionRange(tSyntheticStart, tSyntheticEnd);

				// read in wavelet
				SuReader suWaveletReader = new SuReader();
				SuSegyFormat tr_wavelet = new SuSegyFormat();
				wavelet = SM.getWavelet(suWaveletReader, tr_wavelet, wavelet_filename);
				tr_wavelet = null;

				propFetcher1 = new TracePropertyFetcher(this); // will be used to access R_near/R_far
				String R_word_fc=null;
				if (property_name2.equals("R_near"))
				  {
					propFetcher1.setIndex(property_name2,1,false); 	
					R_word_fc = "R_near_fc";
				  }
				else if (property_name2.equals("R_far"))
				  {
					propFetcher1.setIndex(property_name2,1,false);
					R_word_fc = "R_far_fc";
				  }
				else if (property_name2.equals("R_near_fc"))        // show only fluid contact synthetic...
				  {
					propFetcher1=null;
					R_word_fc = "R_near_fc";
				  }
				else if (property_name2.equals("R_far_fc"))    // show only fluid contact synthetic...
				  {
					propFetcher1=null;
					R_word_fc = "R_far_fc";
				  }
				else
				  {
					EH.die_err("Synthetic seismics must be computed against fields R_near[_fc] or R_far[_fc]");
				  }
				propFetcher2 = new TracePropertyFetcher(this); // will be used to access times
				propFetcher2.setIndex("t",1,false);

				// Set up access to fluid contact R's if available
				if (completeWords.indexOf(R_word_fc) != -1)
				  {
					// reflectivity
					propFetcher3 = new TracePropertyFetcher(this);
					propFetcher3.setIndex(R_word_fc,1,false);
					// fluid contact time
					propFetcher4 = new TracePropertyFetcher(this);
					propFetcher4.setIndex("t_fc",1,false);	
				  }
				

				// set sampling time same as wavelet
				delta_t = (float)((double)wavelet.dt/1.0e6);
				cleanSyntheticRange();
				// compute number of samples in synthetic seismic
				new_ns = 1 + (int)((tSyntheticEnd - tSyntheticStart)/delta_t);   
			  }
		  }

		if ((getRealisationNumber(tr) > lastRealisationNumberHeaderField) || accepted_trace_number==0)
		  {
			// test if dud trace
			realisationGroupIsNull |= isTraceNull(values);    // already dud, or is dud now...
			// get required quantity
			value = propFetcher0.getRequiredQuantityFromTrace(values,property_name1);
			// accumulate data in list
			dataList.add(value);
			lastRealisationNumberHeaderField = getRealisationNumber(tr);

			// store required stuff for seismic if desired
			if (makeSyntheticSeismic) 
			  {
				SyntheticSeismicData ssd = new SyntheticSeismicData(n_layers);
				// Store times and reflection coefficients in array for later retrieval
				for (i=0;i<n_layers;i++)
				  {
					ssd.times[i] = (double) values[propFetcher2.refIndx0 + i];
					ssd.refCoeffs[i] = (propFetcher1 != null) ? (double) values[propFetcher1.refIndx0 + i] : 0.0;
				  }
				// fluid contact info
				if (propFetcher4!=null)
				  {
					for (i=0;i<n_layers;i++)
					  {
						ssd.times_fc[i] = (double) values[propFetcher4.refIndx0 + i];
						ssd.refCoeffs_fc[i] = (double) values[propFetcher3.refIndx0 + i];
					  }
				  }
				// add to HashMap, as keyed by the value we are doing the stats on
				SyntheticSeismicDataMap.put(new Double(value), ssd);
				// must also keep track of header fields
				// => keep current copy of all header fields in the tr_out object, so they're in place
				// when we want to write out
				tr_out.copyHeader(tr);
			  }
		  }
	  }
  }

  void makeSuiteOfStatisticalMomentMaps(ArrayList words, boolean finalise)
  {
	int j,k,N=1,Navs=1,Nmodel,last_ppty_indx_start;
	double value;
	double estimate[];

	float[] values=null;

	// This new trace may be the start of a new sequence, or the end of file (finalise==true). 
	// If so complete current statistical ensemble
	if (finalise || (getRealisationNumber(tr) <= lastRealisationNumberHeaderField))
	  {
		N = words.size();
		Navs = (N-1)*n_layers+1;    // t_base is last word, and only gets one extra place, not another block
		Nmodel = Navs + n_layers;
		last_ppty_indx_start = (N-2)*n_layers; //OLD		last_ppty_indx_start = (N-2)*n_layers;		
		estimate = new double[Navs]; // last (n_layers+1) entries here will be average times
		// Form output. 
		// (1) Copy all header fields: NB!! This is done below on last (!finalise) trace
		// (2) set size of output trace suitably. t_base never has its own block; it's in the last block
		// extended by 1.
		tr_out.set("ns",(short)Nmodel);
		// and set realisation counter to zero
		setRealisationNumber(tr_out,0);   // will disappear in surange
		// check data is available (may be filtered to null)
	
		if (dataLists != null)
		  {
			// compute stats
			// Loop over words
			for (j=0,k=0;j<N;j++)
			  {
				property_name1 = (String)words.get(j);
				// loop over layers, also k<Navs to suppress extra t_base entries
				for (layerNumber=1; layerNumber<=n_layers && k<Navs; layerNumber++)
				  {
					// The last block should probably be the time 50% quantile (the median), since this
					// will show pinchouts nicely on the INT viewer display. If we're doing the
					// last block, get these quantiles...
					if (k>=last_ppty_indx_start)
					  {
						dataLists[k].sort();
						estimate[k] = Descriptive.quantile(dataLists[k],0.50); // median
						// here's the place to write these median times into this last block
						tr_out.set_data_value((float)estimate[k] , k+n_layers);  // ***
					  }
					// OK , now the other properties
					// One of the magic words is the Raftery thinning estimate
					if (property_name1.equals(THINNING_KEYWORD))
					  {
						// OK, this is the Raftery thinning estimate
						estimate[k] = (double)getRafteryThinningEstimate(dataLists[k]);
					  }
					else
					  {
						// Quantiles.
						if (makeFullQuantileSummaryTraces)
						  {
							// can skip sort if we're on the last block (done above)
							if (k<last_ppty_indx_start) dataLists[k].sort();
							estimate[k] = Descriptive.quantile(dataLists[k],quantile);
						  }
						else 
						  {
							// get mean
							estimate[k] = Descriptive.mean(dataLists[k]);
							// variance
							if (makeFullStdDevSummaryTraces)
							  {
								// now get variance (overwrite estimate[k])
								estimate[k] = Descriptive.sampleVariance(dataLists[k],estimate[k]);
								// convert to std_dev
								if (estimate[k]>0) estimate[k] = Math.sqrt(estimate[k]);
							  }
						  }
					  }
					/* ascii impression */
					/*
					  EH.warn("prop = "+words.get(j)+" layer = "+layerNumber+", N(samples) = "+dataLists[k].size()+
					  " statistic= " +estimate[k]);
					*/
					// stuff required property into the required float value. Don't overwrite element (N-1)n_layers.
					if (k<Navs-1)
					  {
						tr_out.set_data_value((float)estimate[k] , k);
					  }
					// ELSE: must be t_base. Don't overwrite element (N-1)n_layers. tbase has been grabbed above ****
					k++;
				  }
			  }
		  }
		else
		  {
			EH.die_err("No data surviving for sane statistics in makeSuiteOfStatisticalMomentMaps(). "
					   +"Too-stringent filters.");
		  }
		if (realisationGroupIsNull) tr_out.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);

		// write out trace here............
		suwriter.fputtr(out_strm, tr_out);

		if (!finalise)
		  {
			// reset statistics summation and counter variables
			// flush lists.
			if (dataLists != null)
			  {
				for (k=0;k<dataLists.length;k++) 
				  {
					dataLists[k].removeFromTo(0,dataLists[k].size()-1);
					dataLists[k].trimToSize();
				  }
				lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
			  }
			realisationGroupIsNull = false;
		  }
		estimate=null; // GC
	  }

	// on an ordinary trace (not EOF), we must increment the running statistics
	if (!finalise)
	  {
		// read data from trace
		values = tr.get_data_array();

		N=words.size();
		// on first trace set up the following....
		if (accepted_trace_number==0)
		  {
			// Compute n_layers
			setNumberOfLayers();
			propFetcher0 = new TracePropertyFetcher(this);
			Navs=N*n_layers;
			dataLists = new DoubleArrayList[Navs];
			for (j=0;j<dataLists.length;j++) dataLists[j] = new DoubleArrayList();
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
			realisationGroupIsNull = isTraceNull(values);
		  }
		if (getRealisationNumber(tr) > lastRealisationNumberHeaderField)
		  {
			// test if dud trace
			realisationGroupIsNull |= isTraceNull(values);    // already dud, or is dud now...
			// loop over words to acquire averages of, accumulating samples
			for (j=0,k=0;j<N;j++)
			  {
				property_name1 = (String) words.get(j);
				// loop over layers
				for (layerNumber=1; layerNumber<=n_layers; layerNumber++)
				  {
					// set correct index of desired property
					propFetcher0.setIndex(property_name1,layerNumber, (layerNumber!=1));
					value = propFetcher0.getRequiredQuantityFromTrace(values,property_name1);
					//accumulate data in lists
					dataLists[k].add(value);
					k++;
				  }
			  }
		  }
		lastRealisationNumberHeaderField = getRealisationNumber(tr);
		// and keep current copy of all header fields in the tr_out object, so they're in place
		// when we want to write out
		tr_out.copyHeader(tr);
	  }
	values = null; // GC
  }


  /**
	Here we form the full mean and covariance of the suite of properties
	listed in the XML file that drives the massing analysis, and write this out to a summary
	file.
  */
  void generateMomentAnalysisforMassaging(boolean finalise)
  {
	int i,j,k,N=1,n;
	double tmp;
	double value;
	double mean;
	boolean useBinaryFormat=false;
	float f;

	float[] values=null;
	
	useBinaryFormat = (WRITE_SCALE_COMPRESSED_BINARY_MOMENT_ANALYSIS && !massageAnalyseAscii);

	// This new trace may be the start of a new sequence, or the end of file (finalise==true). 
	// If so complete current statistical ensemble
	if (finalise || (getRealisationNumber(tr) <= lastRealisationNumberHeaderField))
	  {
		N = massagingKeywords.size();
		// compute stats using standard unbiased estimator of covariance
		n = dataLists[0].size();
		
		if (n <= 0)
		  {
			EH.warn("Too few realisations (maybe null traces?) at: "+hwb.toString());
		  }
		else
		  { 
			// use Y in transposed view compared to Trimm
			DenseDoubleMatrix2D Y = new DenseDoubleMatrix2D(N,n);
			DenseDoubleMatrix2D sampleCov;
			DenseDoubleMatrix2D mu = new DenseDoubleMatrix2D(1,N);
			for (j=0;j<N;j++)
			  {
				for (i=0;i<n;i++)
				  {
					Y.setQuick(j,i,dataLists[j].getQuick(i));
				  }
			  }

			// loop over massagingKeywords to get stats
			for (j=0;j<N;j++)
			  {
				// get mean
				mean = Descriptive.mean(dataLists[j]);
				mu.setQuick(0,j,mean);
			  }
			// form (-n/(n-1)) mu.mu^T
			double preFactor = (n>1) ? -n/(n-1.0) : 0.0;
			sampleCov = (DenseDoubleMatrix2D)mu.zMult(mu,null,preFactor,0,true,false); 
			// and add in Y.Y^T/(n-1)
			preFactor = (n>1) ? 1.0/(n-1.0) : 0.0;
			Y.zMult(Y,sampleCov,preFactor,1,false,true); 
		

			if (useBinaryFormat)
			  {
				if (false)
				  {
					// possible compression tests...
					double EPS_COMP = 1e-20, SVD_EPS=0.01;
					int BLOCK_SIZE=10;
					int blockNum,blOffset=0;
					DoubleMatrix1D sampleCovDiag = new DenseDoubleMatrix1D(N);
					for (i=0;i<N;i++) sampleCovDiag.setQuick(i,Math.sqrt(Math.max(EPS_COMP,sampleCov.getQuick(i,i))));
					DoubleMatrix2D sampleCovNormedBlock=null;
					blockNum = ((accepted_trace_number/lastRealisationNumberHeaderField)-1) % BLOCK_SIZE;
					blOffset = N*blockNum;
					if (sampleCovNormed == null)
					  {
						// allocate multi-trace blocked covariance every BLOCK_SIZE traces
						sampleCovNormed = new DenseDoubleMatrix2D(N,N);
						sampleCovNormedDiff = new DenseDoubleMatrix2D(N,N);
					  }
					// form view of block for this trace
					// sampleCovNormedBlock = sampleCovNormed.viewPart(blOffset,0,N,N);

					// Set block equal to scaled covariance (maybe? minus identity matrix:  sampleCovDiag^1.sampleCov - I)
					for (i=0;i<N;i++)
					  {
						for (j=i;j<N;j++)
						  {
							tmp = sampleCov.getQuick(i,j)/(sampleCovDiag.getQuick(i)*sampleCovDiag.getQuick(j));
							if (i==j) tmp=1.0;
							sampleCovNormedDiff.setQuick(i,j, tmp);
							sampleCovNormedDiff.setQuick(j,i, tmp);
						  }
					  }
					if (blockNum == 0) 
					  {
						sampleCovNormed.assign(sampleCovNormedDiff);
					  }
					else
					  {
						// form diffs
						sampleCovNormedDiff.assign(sampleCovNormed,cern.jet.math.PlusMult.minusMult(1.0));
					  } 
					// Perform SVD when diffed covariance is fully populated:
					Util.debugPrintMatrixToFile("Cres."+blockNum+".txt",sampleCovNormedDiff);
					SingularValueDecomposition svd = new SingularValueDecomposition(sampleCovNormedDiff);
					DoubleMatrix2D U = svd.getU();
					// DoubleMatrix2D V = svd.getV();
					double[] w = svd.getSingularValues();
			
					// Util.print(w,N);
					// SVD scale is unity, since this is the residual matrix after removing the identity matrix
					for (i=0;i<N;i++) if (w[i]<SVD_EPS) break;
					System.out.println("MAX_INT"+Short.MAX_VALUE+" SVD compression factor: "+i/(double)N);
				  }


				// could compute this from e.g. sqrt(N) type things....
				compressed_bytes_per_float = 2;


				// Set upper triangle equal to correlation coefficients
				// diagonal is still the correlation coeff.
				for (i=0;i<N;i++)
				  {
					// off diag turned into corr. coefficients
					// System.out.println("COVOUT: "+sampleCov.getQuick(i,i));	
					for (j=i+1;j<N;j++)
					  {
						if ( (sampleCov.getQuick(i,i) > 0.0) && (sampleCov.getQuick(j,j) > 0.0))
						  {
							tmp = Math.sqrt(sampleCov.getQuick(i,i)*sampleCov.getQuick(j,j));
							tmp = sampleCov.getQuick(i,j)/tmp;
						  }
						else
						  {
							tmp=0.0;
						  }

						// tmp is now a correlation coeff in [-1,1]
						// scale by Short.MAX_VALUE for later rounding
						sampleCov.setQuick(i,j, tmp * Short.MAX_VALUE);
						// System.out.println("COVOUT: "+(short)(Math.rint(sampleCov.getQuick(i,j))));				
					  }
				  }

			  }


			// Output

			if (!useBinaryFormat)
			  {
				// ASCII, unscaled
				// header words
				for (i=0;i<massagingHeaderWords.length;i++)
				  {
					aux_out_writer.print(massagingHeaderWordsValues[i]+" ");
				  }
				// mean
				for (i=0;i<N;i++)
				  {
					aux_out_writer.print(mu.getQuick(0,i)+" ");
				  }
				// covar
				for (i=0;i<N;i++)
				  {
					for (j=i;j<N;j++)
					  {
						aux_out_writer.print(sampleCov.getQuick(i,j)+" ");
					  }
				  }
				aux_out_writer.println("");
			  }
			else
			  {
				// binary
				try
				  {			
					// header words
					for (i=0;i<massagingHeaderWords.length;i++)
					  {
						f = Float.parseFloat(massagingHeaderWordsValues[i]);
						// System.out.println("HW "+f);					
						CastUtils.typecastFloatToBytes(f, work_bytes, 0);
						aux_out_binary_writer_buf.write(work_bytes);
					  }
					// Mean. Full big endian floats.
					for (i=0;i<N;i++)
					  {
						// aux_out_writer.print(mu.getQuick(0,i)+" ");
						f = (float)(mu.getQuick(0,i));
						CastUtils.typecastFloatToBytes(f, work_bytes, 0);
						aux_out_binary_writer_buf.write(work_bytes);
					  }
					// Covar. 
					// compressed_bytes_per_float-Byte floats (e.g. 2 most 
					// sig bytes from 4 byte float: omit least sig 2 bytes from mantissa)
					// off--diag elements are correl. coefficients scaled by MAX_SHORT.
					for (i=0;i<N;i++)
					  {
						// diagonal as float, with two least sig bytes zeroed: 
						// is a 1% approximation for scale, works for Big endian (java) floats
						CastUtils.typecastFloatToBytes((float)(sampleCov.getQuick(i,i)), work_bytes, 0);
						// implies work_bytes[2]=work_bytes[3]=0;
						aux_out_binary_writer_buf.write(work_bytes,0,compressed_bytes_per_float);
						// upper triangle correlation coeffs.
						for (j=i+1;j<N;j++)
						  {
							// write as signed shorts! sizeof(Short)=2
							CastUtils.typecastShortToBytes((short)(Math.rint(sampleCov.getQuick(i,j))), work_bytes, 0);
							aux_out_binary_writer_buf.write(work_bytes,0,2);
						  }
					  }
				  }
				catch (IOException  e) 
				  {
					EH.die_err("Failed to write binary output to "+massagingOutputFileName);
				  }
			  }
			// GC
			Y=null;
			sampleCov=null;
			mu=null;
		  }
		if (!finalise)
		  {
			// reset statistics summation and counter variables
			// flush lists.
			for (j=0;j<dataLists.length;j++) 
			  {
				dataLists[j].removeFromTo(0,dataLists[j].size()-1);
				dataLists[j].trimToSize();
			  }
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
		  }
	  }



	// on an ordinary trace (not EOF), we must increment the running statistics
	if (!finalise)
	  {
		// read data from trace
		values = tr.get_data_array();

		N = massagingKeywords.size();
		// on first trace set up the following....
		if (accepted_trace_number==0)
		  {
			// Compute n_layers
			setNumberOfLayers();
			propFetcher0 = new TracePropertyFetcher(this);
			dataLists = new DoubleArrayList[N];
			// alloc space for the averages to acquire
			for (i=0;i<N;i++) dataLists[i] = new DoubleArrayList();
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
			// check words
			for (i=0;i<N;i++)
			  {
				AnalysingProperty prop = ((AnalysingProperty)massagingKeywords.get(i));
				if (prop.layerNumber<1 || prop.layerNumber>n_layers)
				  {
					EH.die_err("layer number not legal (<= 0 or too great for model)");
				  }
			  }
			// check retained header words
			for (i=0;i<massagingHeaderWords.length;i++)
			  {
				if (tr.headerMap.get(massagingHeaderWords[i]) == null)
				  {
					EH.die_err(massagingHeaderWords[i] +" is an illegal su header word");
				  }
			  }
			massagingHeaderWordsValues = new String[massagingHeaderWords.length];

			// Header in summary file
			if (useBinaryFormat)
			  {
				aux_out_writer.println("Binary-compressed [corr_coeff_scale= "+Short.MAX_VALUE + " n_bytes_covariance= 2] Massaging analysis summary file");
			  }
			else
			  {
				aux_out_writer.println("Massaging analysis summary file");
			  }
			aux_out_writer.println(""+(massagingHeaderWords.length+N+(N*(N+1))/2));
			// header words
			for (i=0;i<massagingHeaderWords.length;i++)
			  {
				aux_out_writer.println(massagingHeaderWords[i]);
			  }
			// properties
			for (i=0;i<N;i++) 
			  {
				AnalysingProperty prop = ((AnalysingProperty)massagingKeywords.get(i));
				aux_out_writer.println(prop.name+"_"+prop.layerNumber);
			  }
			// covariances: top triangle only (is symmetric)
			for (i=0;i<N;i++) 
			  {
				for (j=i;j<N;j++) 
				  {
					aux_out_writer.println("C_"+(i+1)+""+(j+1));
				  }
			  }
			// if going to binary, close file here
			if (useBinaryFormat)
			  {
				// close ascii stream
				aux_out_writer.close();
				// append binary
				try
				  {
					aux_out_binary_writer = (FileOutputStream) new FileOutputStream(massagingOutputFileName,true);
				  }
				catch (IOException e) 
				  {
					throw new RuntimeException("couldn't append to " + massagingOutputFileName);
				  }
				aux_out_binary_writer_buf = new BufferedOutputStream(aux_out_binary_writer);
				work_bytes = new byte[4];  // needed work array
			  }
		  }
		boolean isNull = isTraceNull(values);
		if (getRealisationNumber(tr) > lastRealisationNumberHeaderField)
		  {
			if (!isNull)
			  {
				// loop over massagingKeywords to acquire averages of, accumulating samples
				for (j=0;j<N;j++)
				  {
					AnalysingProperty prop = ((AnalysingProperty)massagingKeywords.get(j));
					propFetcher0.setIndex(prop.name, prop.layerNumber,false);
					value = propFetcher0.getRequiredQuantityFromTrace(values,prop.name);
					//accumulate data in lists
					dataLists[j].add(value);
				  }			  
				// save some salient header words
				for (i=0;i<massagingHeaderWords.length;i++)
				  {
					massagingHeaderWordsValues[i] = getStringEquivalent(tr,massagingHeaderWords[i]);
				  }
			  }
		  }
		lastRealisationNumberHeaderField = getRealisationNumber(tr);
		// and keep current copy of all header fields in the tr_out object, so they're in place
		// when we want to write out
		tr_out.copyHeader(tr);

		if (isNull)
		  {
			// save some salient header words for err msg
			hwb.setLength(0);
			for (i=0;i<massagingHeaderWords.length;i++)
			  {
				hwb.append(massagingHeaderWords[i]).append("=").
				  append(getStringEquivalent(tr,massagingHeaderWords[i])).append(" ");
			  }
		  }
	  }
	values = null; // GC
  }



  /**
	Here we form the p50 stats of the named properties and generate this as 
	a SU property file, suitable for running the Inversion code.
  */
  void makePriorTraces(boolean finalise)
  {
	int i,j,k,N=1,n;
	double tmp;
	double value;
	double median;

	float[] values=null;
	

	// This new trace may be the start of a new sequence, or the end of file (finalise==true). 
	// If so complete current statistical ensemble
	if (finalise || (getRealisationNumber(tr) <= lastRealisationNumberHeaderField))
	  {
		// First trace never enters this block
		N = priorKeywordsArray.size();

		// compute stats using standard unbiased estimator of covariance
		n = dataLists[0].size();
		double[] mu = new double[N];

		// loop over keywords to get stats
		for (j=0;j<N;j++)
		  {
			// get median
			median = Descriptive.median(dataLists[j]);
			mu[j] = median;
			// AnalysingProperty prop = ((AnalysingProperty)priorKeywordsArray.get(j));
			// System.err.printf("prop = %10s ,layer=%5d, median=%10.6g%n", prop.name, prop.layerNumber, median);		
		  }
		// SU model traces output here!
		tr_out.set("ns",(short)N);
		// clean output header
		tr_out.set("delrt",(short)0);
		tr_out.set("f1",0f);
		// and set realisation counter to zero
		setRealisationNumber(tr_out,0);   // will disappear in surange
		for (j=0;j<N;j++)
		  {
			tr_out.set_data_value((float)mu[j],j);  // assumes t_base will be last
		  }
		// write out trace here............
		suwriter.fputtr(out_strm, tr_out);

		// GC
		mu=null;
		if (!finalise)
		  {
			// reset statistics summation and counter variables
			// flush lists.
			for (j=0;j<dataLists.length;j++) 
			  {
				dataLists[j].removeFromTo(0,dataLists[j].size()-1);
				dataLists[j].trimToSize();
			  }
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;
		  }
	  }

	// on an ordinary trace (not EOF), we must increment the running statistics
	if (!finalise)
	  {
		// read data from trace
		values = tr.get_data_array();

		// on first trace set up the following....
		if (accepted_trace_number==0)
		  {
			// Compute n_layers
			boolean mustAddTbase=false;
			setNumberOfLayers();
			propFetcher0 = new TracePropertyFetcher(this);
			lastRealisationNumberHeaderField = getRealisationNumber(tr)-1;

			// List wordsArray = Arrays.asList(priorKeywords).;
			ArrayList<String> wordsArray = new ArrayList();
			wordsArray.addAll(Arrays.asList(priorKeywords));
			wordsArray.remove("t_base");
			if (wordsArray.contains("t"))
			  {
			    mustAddTbase = true;
				if (wordsArray.indexOf("t") != wordsArray.size()-1)
				  {
					EH.die_err("property t must be last in list of properties for --make-prior-traces");
				  }
			  }

			for (i=0;i<wordsArray.size();i++)
			  {
				if (completeWords.indexOf(wordsArray.get(i)) == -1)
				  {
					EH.warn("property " + wordsArray.get(i) + " isn't a recognised keyword "+
							(gotFirstTrace ? " for this data set" : "") + ". Conceivable are:");
					showAvailableProperties();
					EH.die_err("");
				  }
			  }

			// gathers stats for all layers, all named properties
			priorKeywordsArray = new ArrayList();
			for (i=0;i<wordsArray.size();i++)
			  {
				for (j=1;j<=n_layers;j++)
				  {
					priorKeywordsArray.add(new AnalysingProperty(wordsArray.get(i),j));
				  }
			  }
			if (mustAddTbase) priorKeywordsArray.add(new AnalysingProperty("t_base",n_layers));

			priorKeywordsArray.trimToSize();

			N = priorKeywordsArray.size();

			dataLists = new DoubleArrayList[N];
			// alloc space for the averages to acquire
			for (i=0;i<N;i++) dataLists[i] = new DoubleArrayList();


		  }
		else
		  {
			N = priorKeywordsArray.size();
		  }

		if (getRealisationNumber(tr) > lastRealisationNumberHeaderField)
		  {
			// loop over keywords to acquire averages of, accumulating samples
			for (j=0;j<N;j++)
			  {
				AnalysingProperty prop = ((AnalysingProperty)priorKeywordsArray.get(j));
				propFetcher0.setIndex(prop.name, prop.layerNumber,false);
				value = propFetcher0.getRequiredQuantityFromTrace(values,prop.name);
				//accumulate data in lists
				dataLists[j].add(value);
			  }
		  }
		lastRealisationNumberHeaderField = getRealisationNumber(tr);
		// and keep current copy of all header fields in the tr_out object, so they're in place
		// when we want to write out
		tr_out.copyHeader(tr);
	  }
	values = null; // GC
  }

  private void makeSyntheticSeismics()
  {
	int i;
	float values[];
	String R_word_fc=null;
	GraphingData gd=null;
	AnPlotInfo api=null;

	// read data  from trace
	values = tr.get_data_array();

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		// Compute n_layers
		setNumberOfLayers();
		// set up Seismic model & initialise it
		SM = new OneDSeismicConvolutionalModel(this);
		SM.initialiseSeismicModel(n_layers);

		// read in wavelet into new trace
		SuReader suWaveletReader = new SuReader();
		SuSegyFormat tr_wavelet = new SuSegyFormat();
		wavelet = SM.getWavelet(suWaveletReader, tr_wavelet, wavelet_filename);
		tr_wavelet = null;
		// set sampling time same as wavelet
		delta_t = (float)((double)wavelet.dt/1.0e6);
		// Slightly bump up tSyntheticStart if needed
		tSyntheticStart = delta_t*Math.floor(tSyntheticStart/delta_t);

		// set time range of synthetic seismic
		SM.setSyntheticAcquisitionRange(tSyntheticStart, tSyntheticEnd);

		synthSeisFetchers = new FetcherSet();
		synthSeisFetchers.propFetcher0 = new TracePropertyFetcher(this);      // will be used to access R_near/R_far
		if (property_name1.equals("R_near"))
		  {
			synthSeisFetchers.propFetcher0.setIndex(property_name1,1,false); 	
			R_word_fc = "R_near_fc";
		  }

		// Set up index of reflection coefficients
		synthSeisFetchers.propFetcher0 = new TracePropertyFetcher(this);      // will be used to access R_near/R_far
		if (property_name1.equals("R_near"))
		  {
			synthSeisFetchers.propFetcher0.setIndex(property_name1,1,false); 	
			R_word_fc = "R_near_fc";
		  }
		else if (property_name1.equals("R_far"))
		  {
			synthSeisFetchers.propFetcher0.setIndex(property_name1,1,false);
			R_word_fc = "R_far_fc";
		  }
		else if (property_name1.equals("R_near_fc"))        // show only fluid contact synthetic...
		  {
			synthSeisFetchers.propFetcher0 = null;
			R_word_fc = "R_near_fc";
		  }
		else if (property_name1.equals("R_far_fc"))    // show only fluid contact synthetic...
		  {
			synthSeisFetchers.propFetcher0 = null;
			R_word_fc = "R_far_fc";
		  }
		else
		  {
			EH.die_err("Synthetic seismics must be computed against fields R_near[_fc] or R_far[_fc]");
		  }
		// Set up index of times
		synthSeisFetchers.propFetcher1 = new TracePropertyFetcher(this);
		synthSeisFetchers.propFetcher1.setIndex("t",1,false);
		
		// Set up access to fluid contact Rs if available
		if (completeWords.indexOf(R_word_fc) != -1)
		  {
			// reflectivity
			synthSeisFetchers.propFetcher2 = new TracePropertyFetcher(this);
			synthSeisFetchers.propFetcher2.setIndex(R_word_fc,1,false);
			// fluid contact time
			synthSeisFetchers.propFetcher3 = new TracePropertyFetcher(this);
			synthSeisFetchers.propFetcher3.setIndex("t_fc",1,false);	
		  }

		
		// compute number of samples in synthetic seismic
		new_ns = 1 + (int)((tSyntheticEnd - tSyntheticStart)/delta_t);

		if (seismic_filename != null)
		  {
			// prepend true seismic traces
			outputWindowedSeismicTraces(new_ns,tSyntheticStart,wavelet.dt);
		  }
	  }

	if (graphics) 
	  {
		gd = (GraphingData)GraphingDataMap.get(1);
		api = new AnPlotInfo();
		if (layerRangeSpecified)
		  {
			api.Ymin = t_start;
			api.Ymax = t_end;
		  }
	  }

	// now for all traces...
	// test if dud trace
	boolean isNull = isTraceNull(values);
	// get whatever other header fields are useful
	tr_out.copyHeader(tr);
	// Set times and reflection coefficients in OneDSeismicConvolutionalModel object
	// convert times to secs
	int offset0, offset1, offset2=0,offset3=0;
	offset0 = (synthSeisFetchers.propFetcher0 != null) ? synthSeisFetchers.propFetcher0.refIndx0 : 0;
	offset1 = synthSeisFetchers.propFetcher1.refIndx0;
	if (synthSeisFetchers.propFetcher3!=null)
	  {
		offset2 = synthSeisFetchers.propFetcher2.refIndx0;
		offset3 = synthSeisFetchers.propFetcher3.refIndx0;
	  }
	for (i=0;i<n_layers;i++)
	  {
		SM.layer_top_times[i] = (double) 0.001*values[offset1+i]; 
		// set layer reflectivities to 0 if only fluid contacts required...	
		SM.R_near[i] = (synthSeisFetchers.propFetcher0 != null) ? (double) values[offset0 + i] : 0;
		// fluid contact info
		if (synthSeisFetchers.propFetcher3!=null)
		  {
			SM.t_fc[i] = (double) 0.001*values[offset3+i];
			SM.R_near_fc[i] = (double) values[offset2+i];
		  }
	  }
	// set no. of samples correctly
	if (!graphics) tr_out.set("ns",(short)new_ns);
	if (!isNull)
	  {
		// endpoint warnings
		if (SM.layer_top_times[0] < tSyntheticStart)
		  {
			EH.die_err("Error: Synthetic seismic must be computed for range starting less than "+SM.layer_top_times[0]);
		  }
		if (SM.layer_top_times[n_layers-1] > tSyntheticEnd)
		  {
			EH.die_err("Error: Synthetic seismic must be computed for range at least as deep as "+SM.layer_top_times[n_layers-1]);
		  }
		
		// now form Synthetic
		
		SM.generateSyntheticTrace_kernel(wavelet);
		
		// copy Synthetic seismic back to trace
		if (!graphics)
		  {
			for (i=0; i<new_ns ; i++)
			  {
				tr_out.set_data_value((float)SM.padded_synth_near[i],i);
			  }
		  }
		else
		  {
			DoubleArrayList X = api.x;
			DoubleArrayList Y = api.y;
			for (i=0; i<new_ns ; i++)
			  {
				Y.add(tSyntheticStart + i*delta_t);
				X.add(SM.padded_synth_near[i]);
			  }
			api.rNumber = getRealisationNumber(tr);
			if (api.rNumber <= lastRealisationNumberHeaderField) groupOffsetNumber++;
			api.offsetNumber = groupOffsetNumber;
			api.trim();
			gd.plotList.add(api);
		  }
		lastRealisationNumberHeaderField = getRealisationNumber(tr);
	  }
	else
	  {
		tr_out.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);
	  }
	if (!graphics)
	  {
		// Set trace header fields suitably
		// set start time
		tr_out.set("f1",(float)tSyntheticStart);
		// set dt as for wavelet
		tr_out.set("dt",wavelet.dt);
		// and delay time
		tr_out.set("delrt",(short)Math.rint((tSyntheticStart*1e3)));
		
		// now write trace out
		suwriter.fputtr(out_strm, tr_out);	
	  }	
  }

  /**
	 Write out trace-filtered su traces in file seismic_filename
	 for supplied seismic file, windowed to match synthetic seismic with parameters
	 new_ns, tSyntheticStart, dt
  */
  void outputWindowedSeismicTraces(int new_ns, double tSyntheticStart, short dt)
  {
	int i,j,i_end,j_end,k,s_trace_number;
	FileInputStream in_strm;
	boolean gotTrace;
	double t,t_d;
	int i_seis_start,i_seis_end,i_synth_start,i_synth_end;

	// Open input file
	try 
	  {
		in_strm = (FileInputStream) new FileInputStream(seismic_filename);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + seismic_filename);
	  }
	SuSegyFormat trs = new SuSegyFormat();
	SuReader seisReader = new SuReader();

	i_synth_start = (int)Math.rint((tSyntheticStart/(dt * 1e-6)));
	i_synth_end = i_synth_start + new_ns - 1;

	// Loop over seismic traces 
	for (s_trace_number=0;;)
	  {
		// get the trace
		if (seisReader.fgettr(in_strm,trs) == 0) break;

		// filter by header words if required
		if (traceFilter.skipTraceIfRequired(trs)) continue;


		// now for all traces...
		// get whatever other header fields are useful
		tr_out.copyHeader(trs);
		if ((s_trace_number==0) && (wavelet.dt != trs.get_short("dt")))
		  {
			throw new RuntimeException("Supplied " + seismic_filename+" has different dt to requested synthetic. Match them!");
		  }
		// Set trace header fields suitably
		// set no. of samples correctly
		tr_out.set("ns",(short)new_ns);
		// set start time
		tr_out.set("f1",(float)tSyntheticStart);
		// set dt as for wavelet
		tr_out.set("dt",dt);
		// and delay time
		tr_out.set("delrt",(short)Math.rint((tSyntheticStart*1e3)));

		short delrt_s;

		if ((delrt_s = trs.get_short("delrt")) != 0)
		  {
			i_seis_start = (int)Math.rint(((1000*delrt_s)/dt));
		  }
		else if ((delrt_s = (short)(1000*trs.get_float("f1"))) != 0)
		  {
			i_seis_start = (int)Math.rint((1000*delrt_s)/dt);
		  } else i_seis_start=0;


		j_end = (int)trs.get_short("ns");
		i_seis_end = i_seis_start + j_end - 1;

		if (!graphics)
		  {
			// padding up to start of seismic data (if any)
			for (i=i_synth_start,k=0; i < i_seis_start; i++,k++) tr_out.set_data_value(0f,k);
			// seismic
			i_end = Math.min(i_seis_end,i_synth_end);
			j = i-i_seis_start;
			for (; ((i<=i_end) && (j<j_end)); i++,j++,k++) 
			  {
				tr_out.set_data_value(trs.get_data_value(j),k);
			  }
			// post padding if needed
			for (; i<=i_synth_end; i++,k++) tr_out.set_data_value(0f,k);
		  }
		else
		  {
			// same as above, stuffed into arrays
			GraphingData gd = (GraphingData)GraphingDataMap.get(1);
			AnPlotInfo api = new AnPlotInfo();
			
			
			DoubleArrayList X = api.x;
			DoubleArrayList Y = api.y;
			// padding up to start of seismic data (if any)
			for (i=i_synth_start,k=0; i < i_seis_start; i++,k++)
			  {
				Y.add(tSyntheticStart + k*delta_t);
				X.add(0);
			  }
			// seismic
			i_end = Math.min(i_seis_end,i_synth_end);
			j = i-i_seis_start;
			for (; ((i<=i_end) && (j<j_end)); i++,j++,k++) 
			  {
				Y.add(tSyntheticStart + k*delta_t);
				X.add(trs.get_data_value(j));
			  }
			// post padding if needed
			for (; i<=i_synth_end; i++,k++)
			  {
				Y.add(tSyntheticStart + k*delta_t);
				X.add(0);
			  }
			api.trim();
			// make thick? & red
			api.lineWidth=1.0f;
			api.c = Color.red;
			if (layerRangeSpecified) 
			  {
				api.Ymin = t_start;
				api.Ymax = t_end;
			  }
			api.offsetNumber = s_trace_number;
			gd.plotList.add(api);
		  }
		// now write trace out
		if (!graphics) suwriter.fputtr(out_strm, tr_out);	

		s_trace_number++;
	  }
	try 
	  {
		in_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + seismic_filename);
	  }
  }


  /**
	 Capture seismic to superpose on layer images in Graphics mode
  */
  SeismicGraphicInfo gCaptureWindowedSeismicTraces(GraphingData gd,double tStart, double tEnd)
  {
	int i,j,i_end,j_end,k,s_trace_number;
	FileInputStream in_strm;
	boolean gotTrace;
	double t,t_d;
	float x,Xmax,Xmin;
	int i_seis_start,i_seis_end,i_start;

	// Open input file
	try 
	  {
		in_strm = (FileInputStream) new FileInputStream(seismic_filename);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + seismic_filename);
	  }
	SuSegyFormat trs = new SuSegyFormat();
	SuReader seisReader = new SuReader();


	// Loop over seismic traces 
	Xmax =-Float.MAX_VALUE;
	Xmin = Float.MAX_VALUE;
	for (s_trace_number=0;;)
	  {
		// get the trace
		if (seisReader.fgettr(in_strm,trs) == 0) break;

		// filter by header words if required
		if (traceFilter.skipTraceIfRequired(trs)) continue;


		// now for all traces...
		float dt = 1.0e-6f*trs.get_short("dt");
		short delrt_s;
		if ((delrt_s = trs.get_short("delrt")) != 0)
		  {
			i_seis_start = (int)Math.rint(((1000*delrt_s)/(1e6*dt)));
		  }
		else if ((delrt_s = (short)(1000*trs.get_float("f1"))) != 0)
		  {
			i_seis_start = (int)Math.rint((1000*delrt_s)/(1e6*dt));
		  } else i_seis_start=0;

		i_start = (int)Math.floor(tStart/dt);
		i_end = (int)Math.ceil(tEnd/dt);

		j_end = (int)trs.get_short("ns");
		i_seis_end = i_seis_start + j_end - 1;

		// stuff into arrays
		AnPlotInfo api = new AnPlotInfo();			
			
		DoubleArrayList X = api.x;
		DoubleArrayList Y = api.y;
		// padding up to start of seismic data (if any)
		for (i=i_start,k=0; i < i_seis_start; i++,k++)
		  {
			Y.add(tStart + k*dt);
			X.add(0);
		  }
		// seismic
		i_end = Math.min(i_seis_end,i_end);
		j = i-i_seis_start;
		for (; ((i<=i_end) && (j<j_end)); i++,j++,k++) 
		  {
			Y.add(tStart + k*dt);
			x = trs.get_data_value(j);
			Xmax = max(Xmax,x);
			Xmin = min(Xmin,x);
			X.add(x);
		  }
		// post padding if needed
		for (; i<=i_end; i++,k++)
		  {
			Y.add(tStart + k*dt);
			X.add(0);
		  }
		api.trim();
		// make thick? & red
	    api.lineWidth=1.0f;
		api.c = Color.red;
		if (layerRangeSpecified) 
		  {
			api.Ymin = t_start;
			api.Ymax = t_end;
		  }
		api.offsetNumber = s_trace_number;
		// api.rNumber=0 (in constructor)
		gd.plotList.add(api);
		  
		s_trace_number++;
	  }
	try 
	  {
		in_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + seismic_filename);
	  }

	if (s_trace_number==0) throw new RuntimeException("seismic with filters allows no traces: " + seismic_filename);

	// return vector with ranges of seismic
	SGI = new SeismicGraphicInfo();
	SGI.Xmin = Xmin;
	SGI.Xmax = Xmax;
	SGI.nPlots = s_trace_number;
	return SGI;
  }



  /**
	 Duplicate input trace, but map last block of (n+1) times to depths
  */
  private void performTimeToDepths()
  {
	int i;
	float depth=0f;
	// read data  from trace
	float values[] = tr.get_data_array();

	// on first trace set up the following....
	if (accepted_trace_number==0)
	  {
		// Compute n_layers
		setNumberOfLayers();
		// Set up index of depths
		propFetcher0 = new TracePropertyFetcher(this);      // will be used to access depths
		propFetcher0.setIndex("d",1,false);   
		// Set up index of times
		propFetcher1 = new TracePropertyFetcher(this);
		propFetcher1.setIndex("t",1,false);
		// Set up index of vp_eff in last layer
		propFetcher2 = new TracePropertyFetcher(this);
		propFetcher2.setIndex("vp_eff",n_layers,false);
	  }
	// mask nulls
	if (isTraceNull(values))
	  {
		tr.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);
	  }
	else
	  {
		// replace depths with times for tops
		for (i=0;i<n_layers;i++)
		  {
			depth = values[propFetcher0.refIndx0+i];
			tr.set_data_value(depth,propFetcher1.refIndx0+i);
		  }
		// and do base via calculation depth_base =d_N+(1/2)0.001(t_base-t_N)*vp_{eff,N}
		// (times in ms and are two-way....)
		depth += 0.0005*(values[propFetcher1.refIndx0+n_layers] -values[propFetcher1.refIndx0+n_layers-1])
	      * values[propFetcher2.refIndx0];
		tr.set_data_value(depth,propFetcher1.refIndx0+n_layers);
	  }
	// now write trace out
	suwriter.fputtr(out_strm, tr);		
  }

 /**
	write accepted traces to stdout
  */
  private void writeOutAcceptedTraces()
  {
	// write trace out
	suwriter.fputtr(out_strm, tr);
  }

  /**
	 This is essentially the same as makeSyntheticSeismics(), but here
	 the times and reflection coefficients are provided: these are
	 usually procured by a statistical procedure (e.g. R_near/R_far
	 may correspond to the median seismic chisq error, or something
	 similar. Also, the calling program must have set up the 
	 OneDSeismicConvolutionalModel object.
  */
  private void makeGivenSyntheticSeismic(SyntheticSeismicData requiredSyntheticSeismicData, boolean isNull)
  {
	int i;

	// Set trace header fields suitably
	// set no. of samples correctly
	tr_out.set("ns",(short)new_ns);
	// set start time
	tr_out.set("f1",(float)tSyntheticStart);
	// set dt as for wavelet
	tr_out.set("dt",wavelet.dt);
	// and delay time
	tr_out.set("delrt",(short)Math.rint((tSyntheticStart*1e3)));
	
	if (!isNull)
	  {
		// Set times and reflection coefficients in OneDSeismicConvolutionalModel object
		// Convert times from ms also
		for (i=0;i<n_layers;i++)
		  {
			SM.layer_top_times[i] = 0.001 * requiredSyntheticSeismicData.times[i];
			SM.R_near[i] = requiredSyntheticSeismicData.refCoeffs[i];
			SM.t_fc[i] = 0.001 * requiredSyntheticSeismicData.times_fc[i];
			SM.R_near_fc[i] = requiredSyntheticSeismicData.refCoeffs_fc[i];
		  }
		// now form Synthetic
		SM.generateSyntheticTrace_kernel(wavelet);
		
		// Here we *could* set the realisation number to zero if desired, but
		// maybe having the number of realisations around is useful ??? 
		// If we want the latter, comment out the next line
		setRealisationNumber(tr_out,0);
		
		// copy Synthetic seismic back to trace (note tr.ns must be correct first)
		for (i=0; i<new_ns ; i++)
		  {
			tr_out.set_data_value((float)SM.padded_synth_near[i],i);
		  }
	  }
	else
	  {
		tr_out.fill_data_array((float)SeismicInverter.MODEL_TRACE_BAD_NUMBER);
	  }

	// now write trace out
	suwriter.fputtr(out_strm, tr_out);	
  }

  // Routines for Raftery style thinning estimates
  private int getRafteryThinningEstimate(DoubleArrayList data)
  {
	int val;
	MCMCConvergenceMonitor CM;
 
	CM = new  MCMCConvergenceMonitor(EH,verbosity,false,true);
 
	// set estimands..
	CM.initialiseRafteryStyleConvergenceDiagnostics(1,data.size());
	for (int k=0; k<data.size(); k++)
	  {
		CM.setEstimand(0,data.get(k));	
		CM.completeScalarEstimandsUpdate(k);
	  }
	CM.finishRafteryStyleConvergenceDiagnostics();
	val = CM.Raftery_kthin[0];
	CM =null;
	return val;
   }


  private void setNumberOfLayers()
  {
	MF.recomputeLayersAndMaps((int)tr.get_short("ns"),
							  (int)tr.get_int(ModelFormat.OUTPUT_FIELDS_SU_KEYWORD));
	n_layers = MF.getNumberOfLayers();

	// recompute words and maps
	//	assembleWordsAndMaps();

	EH.warn("number of layers deduced = "+n_layers);
  }

  private void setUpTraceInputs()
  {
	tr = new SuSegyFormat();
	sureader = new SuReader();
  }

  private void setUpAuxiliaryOutputs()
  {
	if (massageAnalyse)
	  {
		try
		  {
			aux_out_writer = new PrintWriter(new FileWriter(massagingOutputFileName));
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + massagingOutputFileName);
		  }		
	  }
  }

  private void finalizeAuxiliaryOutputs()
  {
	if (massageAnalyse)
	  {
		// close output streams
		aux_out_writer.close();

		if (WRITE_SCALE_COMPRESSED_BINARY_MOMENT_ANALYSIS && !massageAnalyseAscii)
		  {
			try
			  {
				aux_out_binary_writer_buf.flush();
				aux_out_binary_writer_buf.close();
				aux_out_binary_writer.close();
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close binary file" + massagingOutputFileName);
			  }	
		  }
	  }
  }

  /**
	 Elementary checks to see that the realisationNumberHeaderField
	 is a respectable header word
  */
  private void checkHeaderWords()
  {
	realisationNumberSuHeaderMapEntry = (SuSegyFormat.SuHeaderEntry) 
	  tr.headerMap.get(MF.realisationNumberHeaderField);

	if (realisationNumberSuHeaderMapEntry==null)
	  throw new RuntimeException(MF.realisationNumberHeaderField + " is an illegal header word");
	if (realisationNumberSuHeaderMapEntry.type == SuSegyFormat.SuHeadTypes.float_type)
	  {
		throw new RuntimeException(MF.realisationNumberHeaderField + " is a float header word: it must be a short or int.");
	  }
  }


  private void setUpTraceOutputs()
  {
	out_strm = (FileOutputStream) new FileOutputStream(FileDescriptor.out);
	suwriter = new SuWriter();
	// separate trace for outputting if required
	tr_out = new SuSegyFormat();
  }

  private void finalizeTraceInputs()
  {
	// close input streams
	try 
	  {
		inversion_in_strm.close();
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + realisation_filename);
	  }
  }

  private void finalizeOutput()
  {
	// close output streams
	try 
	  {
		suwriter.flush(out_strm);
		out_strm.close();
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close output stream");
	  }
  }

  /**
	 Algorithms for getting frequency of fluid type in data set
  */
 private double getFluidProb(DoubleArrayList dataList, int desiredFluidType)
  {
	int i,count=0;
	for (i=0; i < dataList.size(); i++)
	  {
		if (((int)dataList.get(i)) == desiredFluidType) count++;
	  }
	return (double)count/(double)dataList.size();
  }


  /** 
	  Master method for collecting graphics outputs and displaying
  */
  void initialiseGraphics()
  {
	if (!graphics) return;
	else
	  {
		// Master List of plots
		GraphingDataMap = new OpenIntObjectHashMap(3); // usually not more than 3
		// Synthetic seismic
		if (makeSyntheticSeismic)
		  {
			GraphingData gd = new GraphingData();
			GraphingDataMap.put(1,gd);  // code 1 for synthetic seismic
		  }
		// Layercake
		if (layerCake || layerCakeOfFluids || layerCakeOfNG || layerCakeOfZ_eff || layerCakeOfProp)
		  {
			GraphingData gd = new GraphingData();
			GraphingDataMap.put(2,gd);  // code 2 for layercakes
		  }
		// Ascii scatterplots
		if (propertySUMap || propertyMap || propertyMap2 || propertyMap2gen)
		  {
			GraphingData gd = new GraphingData();
			GraphingDataMap.put(3,gd);  // code 3 for scatterplots
		  }
		// histograms
		if (computeHistogram)
		  {
			GraphingData gd = new GraphingData();
			GraphingDataMap.put(4,gd);  // code 4 for histograms
		  }
		// Structure on seismic
		if (showHorizonsAndSeismic)
		  {
			GraphingData gd = new GraphingData();
			GraphingDataMap.put(5,gd);  // code 5 for structure on seismic
		  }		
	  }
  }

  void showGraphics()
  {
	if (!graphics) return;
	if (GraphingDataMap==null || GraphingDataMap.size()==0) return;

	// here, show stuff
	int i,j,k,nRows=1,nCols,indxPG=0,pSize,N;
	double a=0;
	CanvasPlotter CP;

	GraphingData gd;
	ArrayList<CanvasPlotter.TwoDPlotGroup>Pgroups = null;
	CanvasPlotter.TwoDPlotGroup pg = null;
	CanvasPlotter.TwoDPlotGroup pgSeis = null;
	CanvasPlotter.TwoDPlot plot = null;

	CP = new CanvasPlotter(true);

	CP.reverseVideo=false;

	// At present, line all plots up in a row
	nCols = GraphingDataMap.size();

	Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(nRows * nCols);

	// Synthetic seismic (1)
	gd = (GraphingData)GraphingDataMap.get(1) ;
	if (gd != null) 
	  {
		pgSeis = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pgSeis);
		pgSeis.GG.plotGroupNX = nCols; 
		pgSeis.GG.plotGroupNY = nRows;
		pgSeis.GG.plotGroupXindx = 0;
		indxPG++;
		pgSeis.GG.plotGroupYindx = 0;
		pgSeis.GG.reverseY = true;
		pgSeis.GG.Yword ="t";
		pgSeis.GG.Xword = "seismic amplitude";
		pgSeis.GG.setTickNumbersEngStyle(false,false);
		AnPlotInfo api;
		// stagger traces if needed
		pSize = gd.plotList.size();
		boolean stagger = true;
		if (stagger)
		  {
			// compute good spacing
			double dx, xMin=Double.MAX_VALUE, xMax=-Double.MAX_VALUE;
			int maxOff=0;
			for (i=0;i<pSize;i++) 
			  {
				api = gd.plotList.get(i);
				maxOff = max(maxOff,api.offsetNumber);
				double[] x = api.x.elements(); N = x.length;
				for (j=0;j<N;j++) 
				  {
					xMax = max(xMax,x[j]);
					xMin = min(xMin,x[j]);
				  }
			  }
			// stagger them magic factor here could be controlled....
			double scaleFactor = 2.2*max(xMax,-xMin);
			boolean changeXlabel = false;
			for (i=0;i<pSize;i++) 
			  {
				api = gd.plotList.get(i);
				if (api.offsetNumber != 0) changeXlabel = true;
				dx = scaleFactor*api.offsetNumber;
				double[] x = api.x.elements(); N = x.length;
				for (j=0;j<N;j++) x[j] += dx;
			  }
			if (changeXlabel)
			  {
				pgSeis.GG.Xword = "seismic amplitudes";
				pgSeis.GG.doXticks=false;
			  }
		  }
		// plot in reverse order, so true seismic is uppermost, if present
		for (i=pSize-1; i>=0; i--)
		  {
			api = gd.plotList.get(i);
			plot = new CanvasPlotter.TwoDPlot();
			plot.x = api.x.elements();
			plot.y = api.y.elements();
			if (api.c != null) plot.c = api.c; //color
			plot.lineWidth = api.lineWidth;
			pgSeis.plots.add(plot);
			pgSeis.GG.computeRanges(plot.x,plot.y);		  
			if (layerRangeSpecified) pgSeis.GG.setYRange(api.Ymax,api.Ymin);
		  } 
	  }

	// Layer cake images (2)
	gd = (GraphingData)GraphingDataMap.get(2) ;
	if (gd != null) 
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = nCols; 
		pg.GG.plotGroupNY = nRows;
		pg.GG.plotGroupXindx = indxPG;
		indxPG++;
		pg.GG.plotGroupYindx = 0;
		pg.GG.reverseY = true;
		pg.GG.Yword = layerCakesInDepth ? "depth" : "t";
		pg.GG.Xword = "sample/trace";
		if (gd.plotList.size()==1)
		  {
			// (1) Layer cake images
			AnPlotInfo api = gd.plotList.get(0);
			CanvasPlotter.TwoDShapePlot plotS = new CanvasPlotter.TwoDShapePlot();
			api.trim();
			plotS.x = api.x.elements();
			plotS.x1 = api.x1.elements();
			plotS.y = api.y.elements();
			plotS.y1 = api.y1.elements();
			plotS.col = api.cList.elements();
			plotS.gcol = api.cgList.elements();
			pg.GG.computeRanges(plotS.x,plotS.y);
			pg.GG.computeRanges(plotS.x1,plotS.y1);
			pg.GG.setTickNumbersEngStyle(false,false);
			if (layerRangeSpecified)
			  {
				pg.GG.setYRange(api.Ymax,api.Ymin);
			  }
			if (pgSeis!=null && !layerCakesInDepth)
			  {
				// try to align time scales with seismic
				pg.GG.setYRange(pgSeis.GG.yMax,pgSeis.GG.yMin);
			  }				  
			plotS.normalize();
			if (layerCakeOfFluids) 
			  {
				// black-blue hydro colouring (layer independent)
				plotS.setColorScheme(Color.black,Color.blue);
			  }
			else if (layerCakeOfNG)
			  {
				// Brown-yellow
				plotS.setColorScheme(new Color(101, 67, 33),Color.yellow);
			  }
			else if (layerCakeOfZ_eff)
			  {
				// Z_eff colouring
				plotS.setColorScheme(Color.black,Color.white);
			  }
			else
			  {
				plotS.setColorScheme(Color.black,Color.white);
			  }
			pg.plots.add(plotS);
		  }
	  }


	// Ascii scatterplots (3)
	gd = (GraphingData)GraphingDataMap.get(3) ;
	if (gd != null) 
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = nCols; 
		pg.GG.plotGroupNY = nRows;
		pg.GG.plotGroupXindx = indxPG;
		indxPG++;
		pg.GG.plotGroupYindx = 0;
		if (propertyMap2 || propertyMap2gen) 
		  {
			pg.GG.Xword = property_name1 + ","+layerNumber;
			pg.GG.Yword = property_name2 + ","+(propertyMap2gen ? layerNumber2 : layerNumber);
		  }
		else 
		  {
			pg.GG.Xword = "sample";
			pg.GG.Yword = property_name1+ ","+layerNumber;
		  }
		if (gd.plotList.size()==1)
		  {
			AnPlotInfo api = gd.plotList.get(0);
			api.trim();
			plot = new CanvasPlotter.TwoDPlot();
			plot.x = api.x.elements();
			plot.y = api.y.elements();
			if (api.c != null) plot.c = api.c; //color
			plot.drawline = false;
			plot.symbol='r';
			plot.symSize=0.005;
			pg.plots.add(plot);
			pg.GG.setTickNumbersEngStyle(false,false);
			pg.GG.computeRanges(plot.x,plot.y);
		  }
	  }

	// Histograms (4)
	gd = (GraphingData)GraphingDataMap.get(4) ;
	if (gd != null) 
	  {
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);
		pg.GG.plotGroupNX = nCols; 
		pg.GG.plotGroupNY = nRows;
		pg.GG.plotGroupXindx = indxPG;
		indxPG++;
		pg.GG.plotGroupYindx = 0;
		pg.GG.Xword = property_name1+ ","+layerNumber;
		pg.GG.Yword = "count";
		for (i=0; i<gd.plotList.size(); i++)
		  {
			AnPlotInfo api = gd.plotList.get(i);
			plot = new CanvasPlotter.TwoDPlot(); 
			plot.x = api.x.elements();
			plot.y = api.y.elements();
			// cycle colors from CanvasPlotter.defaultPalette. Skip black=0 & white=1
			plot.c = CanvasPlotter.defaultPalette[i%(CanvasPlotter.defaultPalette.length - 2)+2];
			plot.lineWidth=2.0f;
			pg.plots.add(plot);
			pg.GG.computeRanges(plot.x,plot.y);
		  }
	  }

	// Structure on seismic (5)
	gd = (GraphingData)GraphingDataMap.get(5) ;
	if (gd != null) 
	  {
		pgSeis = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pgSeis);
		pgSeis.GG.plotGroupNX = nCols; 
		pgSeis.GG.plotGroupNY = nRows;
		pgSeis.GG.plotGroupXindx = indxPG;
		indxPG++;
		pgSeis.GG.plotGroupYindx = 0;
		pgSeis.GG.reverseY = true;
		pgSeis.GG.Yword ="t";
		pgSeis.GG.Xword = "seismic amplitude";
		AnPlotInfo api;
		// stagger traces if needed
		pSize = gd.plotList.size();
		// check we got some data
		if (SGI==null) throw new RuntimeException("No data accepted from seismic! Filters too severe?");

		// first, rescale seismic by time range from models
		double yMax =-Double.MAX_VALUE, yMin = Double.MAX_VALUE, Y[];
		if (!layerRangeSpecified)
		  {
			for (i=SGI.nPlots; i<pSize; i++) // range of horizon plots 
			  {
				api = gd.plotList.get(i);
				Y = api.y.elements();
				yMin = min(yMin,Y[0]);
				yMax = max(yMax,Y[Y.length-1]);
			  	if (api.blockPlot) 
				  {
					Y = api.y1.elements();
					yMax = max(yMax,Y[Y.length-1]);
				  }
			  }
			// apply these to seismic plots below..
		  }
		// for blockplot color range
		double cMax = -Double.MAX_VALUE;
		double cMin = -cMax;
		boolean stagger = true;
		if (stagger)
		  {
			double[] x,x1;
			int maxRealisations=0,rNumberMin=Integer.MAX_VALUE;
			// compute good spacing
			double dx, xMin=Double.MAX_VALUE, xMax=-Double.MAX_VALUE;
			int maxOff=0;
			for (i=0;i<pSize;i++) 
			  {
				api = gd.plotList.get(i);
				maxOff = max(maxOff,api.offsetNumber);
				x = api.x.elements(); N = x.length;
				for (j=0;j<N;j++) 
				  {
					xMax = max(xMax,x[j]);
					xMin = min(xMin,x[j]);
				  }
				if (api.blockPlot)
				  {
					x = api.x1.elements(); N = x.length;
					for (j=0;j<N;j++) 
					  {
						xMax = max(xMax,x[j]);
						xMin = min(xMin,x[j]);
					  }
					maxRealisations = max(maxRealisations, api.rNumber);
					rNumberMin = min(rNumberMin,api.rNumber);
				  }
			  }
			// some files have "mark=0" (e.g. MAP). Look after this
			if (maxRealisations==0) maxRealisations=1;

			// Staggering: (a) wiggles superposed (b) blocks offset
			// magic factor here could be controlled....
			double scaleFactor = 2.2*max(xMax,-xMin);
			boolean changeXlabel = false;
			for (i=0;i<pSize;i++) 
			  {
				api = gd.plotList.get(i);
				if (api.offsetNumber != 0) changeXlabel = true;
				if (!api.blockPlot)
				  { // wiggles
					dx = scaleFactor*api.offsetNumber;
					x = api.x.elements(); N = x.length;
					for (j=0;j<N;j++) x[j] += dx;
				  }
				else
				  {
					// blocks
					// a = (x1-x)/nR
					// x1 -> x+(rN+1)*a
					// x ->x+rN*a
					k = api.rNumber - rNumberMin;
					dx = scaleFactor * api.offsetNumber;
					x = api.x.elements(); N = x.length;
					x1 = api.x1.elements(); // assume same length
					if (N>0) a = (x1[0]-x[0])/maxRealisations;
					for (j=0;j<N;j++)
					  {
						x1[j] = dx + (x[j]+(k+1)*a);
						x[j] = dx + (x[j]+k*a);
					  }
					// common color scale
					x = api.cList.elements(); 
					for (j=0;j<N;j++)
					  {
						cMax=max(cMax,x[j]);
						cMin=min(cMin,x[j]);
					  }
				  }
			  }
			if (changeXlabel)
			  {
				pgSeis.GG.Xword = "seismic amplitudes";
				pgSeis.GG.doXticks=false;
			  }
		  }

		// plot in reverse order, so true seismic is uppermost
		for (i=pSize-1; i>=0; i--)
		  {
			api = gd.plotList.get(i);

			if (!api.blockPlot)
			  {
				// pick dashes
				plot = new CanvasPlotter.TwoDPlot();
				plot.x = api.x.elements();
				plot.y = api.y.elements();
				if (api.c != null) plot.c = api.c; //color
				plot.lineWidth = api.lineWidth;
			  
				pgSeis.plots.add(plot);
				pgSeis.GG.computeRanges(plot.x,plot.y);		  
				if (layerRangeSpecified) pgSeis.GG.setYRange(api.Ymax,api.Ymin);
				else pgSeis.GG.setYRange(yMax,yMin);
			  } 
			else 
			  {
				// blocks
				CanvasPlotter.TwoDShapePlot plotS = new CanvasPlotter.TwoDShapePlot();
				plotS.x = api.x.elements();
				plotS.x1 = api.x1.elements();
				plotS.y = api.y.elements();
				plotS.y1 = api.y1.elements();
				plotS.col = api.cList.elements();
				plotS.gcol = api.cgList.elements();
				plotS.setColorRange(cMin,cMax);
				/*
				  if (layerCakeOfFluids) 
				  {
					// black-blue hydro colouring (layer independent)
					plotS.setColorScheme(Color.black,Color.blue);
				  }
				else if (layerCakeOfNG)
				  {
					// Brown-yellow
					plotS.setColorScheme(new Color(101, 67, 33),Color.yellow);
				  }
				else if (layerCakeOfZ_eff)
				  {
					// Z_eff colouring
					plotS.setColorScheme(Color.black,Color.white);
				  }
				else
				  {
					plotS.setColorScheme(Color.black,Color.white);
				  }
				*/
				plotS.setColorScheme(Color.black,Color.white);

				pgSeis.plots.add(plotS);
				pgSeis.GG.computeRanges(plotS.x,plotS.y);	
				
				if (layerRangeSpecified) pgSeis.GG.setYRange(api.Ymax,api.Ymin);
				else pgSeis.GG.setYRange(yMax,yMin);
			  }
		  } 
	  }

	CP.windowTitle = "DeliveryAnalyser";
	CP.run(Pgroups); 
  }

  // info for Horizon plotting
  private static class SeismicGraphicInfo
  {
	float Xmax,Xmin;
	int nPlots;
  }

  // Class for handing properties for Massaging or other things
  private static class AnalysingProperty
  {
	int layerNumber=0;
	String name=null;
	
	AnalysingProperty(String name, int layerNumber)
	{
	  this.layerNumber = layerNumber;
	  this.name = name;
	}
  }


  // Class for handing properties in trace and extracting them conveniently
  private class TraceProperty
  {
	double value=0.0;
	String name=null;
	
	TraceProperty(String name, double value)
	{
	  this.value = value;
	  this.name = name;
	}
  }

  class FetcherSet 
  {
	// Safe set of property fetchers 
	TracePropertyFetcher propFetcher0;   // for property 0
	TracePropertyFetcher propFetcher1;   // for property 1
	TracePropertyFetcher propFetcher2;   // for property 2
	TracePropertyFetcher propFetcher3;   // for property 3
	TracePropertyFetcher propFetcher4;   // for property 4
  }

  class TracePropertyFetcher
  {
	InversionAnalyser _IA;
	ModelFormat _MF;
	// useful for most modes
	int layerNumber;
	int refIndx0;     // generally useful indexes for positions in model trace
	int refIndx1;     // less useful indexes for positions in model trace
	int refIndx_t;
	int refIndx_NG;
	int refIndx_NG_hard;
	int refIndx_vp_eff;
	int refIndx_phi_s;
	int refIndx_S_h;
	int refIndx_fluid_type;

	private int lastIndexOfPropertyNameTimesNlayers=0;

	public TracePropertyFetcher(InversionAnalyser IA)
	{
	  _IA = IA;
	  _MF = IA.MF;
	}
	
	/**
	   here, we set the indexes for the named property in the layer number. 

	   As yet *****unimplemented****:
	   If the property hasn't changed name since the last call, we can gain efficiencies by setting
	   property_name_unchanged==true
	*/
	void setIndex(String property_name, int layerNumberForIndex, boolean property_name_unchanged)
	{
	  // hashMap based efficient method
	  int indx,nameKeyValue;
	  Integer mapInteger = (Integer)_IA.keyWordsMap.get(property_name);
	  if (mapInteger!=null)
		{
		  nameKeyValue = mapInteger.intValue();
		  switch (nameKeyValue)
			{
			case InversionAnalyser.THINNING_KEY:
			  {
				indx = _MF.nameMap.get("t").intValue();
				refIndx_t = (layerNumberForIndex-1) + indx * n_layers;
				break;
			  }
			case InversionAnalyser.THICKNESS_KEY:
			case InversionAnalyser.NET_SAND_KEY:
			case InversionAnalyser.NET_HYDROCARBON_KEY:
			case InversionAnalyser.NET_OIL_KEY:
			case InversionAnalyser.NET_LSG_KEY:
			case InversionAnalyser.NET_GAS_KEY:
			case InversionAnalyser.NET_HARD_KEY:
			  {
				// It's a special combination requiring some calculation. 
				// We need to set up this stuff.....
				// vp_eff indx
				indx = (_MF.nameMap.get("vp_eff")).intValue();
				refIndx_vp_eff = (layerNumberForIndex-1) + indx * n_layers;
				// bounding times: t_above; t_below if obviously at refIndx_t+1...
				indx = (_MF.nameMap.get("t")).intValue();
				refIndx_t = (layerNumberForIndex-1) + indx * n_layers;
				// NG indx (not needed for thickness only)	
				indx = (_MF.nameMap.get("NG")).intValue();
				refIndx_NG = (layerNumberForIndex-1) + indx * n_layers;
				if (nameKeyValue == InversionAnalyser.NET_HARD_KEY)
				  {
					indx = (_MF.nameMap.get("NG_hard")).intValue();
					refIndx_NG_hard = (layerNumberForIndex-1) + indx * n_layers;
				  }
				break;
			  }
			case InversionAnalyser.PHIE_KEY:
			  {
				// NG indx (not needed for thickness only)	
				indx = (_MF.nameMap.get("NG")).intValue();
				refIndx_NG = (layerNumberForIndex-1) + indx * n_layers;
				// porosity
				indx = (_MF.nameMap.get("phi_s")).intValue();
				refIndx_phi_s = (layerNumberForIndex-1) + indx * n_layers;
				break;
			  }
			case InversionAnalyser.BRINE_PROB_KEY:
			case InversionAnalyser.OIL_PROB_KEY:
			case InversionAnalyser.LSG_PROB_KEY:
			case InversionAnalyser.GAS_PROB_KEY:
			  {
				// we just need access to the fluid_type field
				indx = (_MF.nameMap.get("fluid_type")).intValue();
				refIndx_fluid_type = (layerNumberForIndex-1) + indx * n_layers;
				break;
			  }
			}
		  if ((nameKeyValue == InversionAnalyser.NET_HYDROCARBON_KEY) ||
			  (nameKeyValue == InversionAnalyser.NET_OIL_KEY) ||
			  (nameKeyValue == InversionAnalyser.NET_LSG_KEY) ||
			  (nameKeyValue == InversionAnalyser.NET_GAS_KEY))			  
			{
			  // porosity
			  indx = (_MF.nameMap.get("phi_s")).intValue();
			  refIndx_phi_s = (layerNumberForIndex-1) + indx * n_layers;
			  // fluid_type
			  indx = (_MF.nameMap.get("fluid_type")).intValue();
			  refIndx_fluid_type = (layerNumberForIndex-1) + indx * n_layers;

			  // saturation
			  switch (nameKeyValue)
				{
				case InversionAnalyser.NET_HYDROCARBON_KEY:
				  indx = (_MF.nameMap.get("S_h")).intValue();
				  break;
				case InversionAnalyser.NET_OIL_KEY:
				  indx = (_MF.nameMap.get("S_oil")).intValue();
				  break;
				case InversionAnalyser.NET_LSG_KEY:
				  indx = (_MF.nameMap.get("S_lsg")).intValue();
				  break;
				case InversionAnalyser.NET_GAS_KEY:
				  indx = (_MF.nameMap.get("S_gas")).intValue();
				  break;
				}
			  refIndx_S_h = (layerNumberForIndex-1) + indx * n_layers;
			}
		}
	  else
		{
		  // It's a naive property directly in the model trace
		  // compute correct index of desired property
		  mapInteger =  _MF.nameMap.get(property_name);
		  if (mapInteger!=null)
			{
			  nameKeyValue = mapInteger.intValue();
			  if (nameKeyValue == _MF.T_BASE)
				{
				  refIndx0 = nameKeyValue * n_layers;
				}
			  else
				{
				  refIndx0 = (layerNumberForIndex-1) + nameKeyValue * n_layers;
				}
			}
		  else
			{
			  // must be an illegal word
			  EH.warn("Property " + property_name + " isn't available from inverted traces. Available are:");
			  showAvailableProperties();
			  EH.die_err("%n");
			}
		}
	}
	
	// extract the required quantity from the traces using the indices set up in setIndices()
	double getRequiredQuantityFromTrace(float[] values, String property_name)
	{
	  double value=0.0;
	  short fluid_type;

	  // Set of special "hybrid" properties...
	  Integer mapInteger = (Integer)_IA.keyWordsMap.get(property_name);
	  if (mapInteger!=null)
		{
		  int nameKeyValue = mapInteger.intValue();
		  switch (nameKeyValue)
			{
			case InversionAnalyser.BRINE_PROB_KEY:
			  {
				// set unity if brine detected
				fluid_type = (short)values[refIndx_fluid_type];
				fluid_type &= 3; // use majority class
				if (fluid_type == FluidStates.FluidTypes.brine) value=1.0;
				break;
			  }
			case InversionAnalyser.OIL_PROB_KEY:
			  {
				// set unity if oil detected
				fluid_type = (short)values[refIndx_fluid_type];
				fluid_type &= 3; // use majority class
				if (fluid_type == FluidStates.FluidTypes.oil) value=1.0;
				break;
			  }
			case InversionAnalyser.LSG_PROB_KEY:
			  {
				fluid_type = (short)values[refIndx_fluid_type];
				fluid_type &= 3; // use majority class
				// set unity if lsg detected
				if (fluid_type == FluidStates.FluidTypes.lsg) value=1.0;
				break;
			  }
			case InversionAnalyser.GAS_PROB_KEY:
			  {
				// set unity if gas detected
				fluid_type = (short)values[refIndx_fluid_type];
				fluid_type &= 3; // use majority class
				if (fluid_type == FluidStates.FluidTypes.gas) value=1.0;
				break;
			  }
			case InversionAnalyser.THINNING_KEY:
			  {
				// thinning estimate performed on sequence of times (in ms <- scale doesn't matter)
				value = values[refIndx_t];
				break;
			  }
			case InversionAnalyser.THICKNESS_KEY:
			  {
				// thickness = 0.5*(t_below - t_above) x vp_eff
				// times in ms too x factor of 1/2 since times are two-way...
				value = 0.0005*(double)((values[refIndx_t+1] - values[refIndx_t]) * values[refIndx_vp_eff]);
				break;
			  }
			case InversionAnalyser.NET_SAND_KEY:
			  {
				// net sand = thickness x NG = = (t_below - t_above) x vp_eff x NG
				// times in ms too x factor of 1/2 since times are two-way...
				value = 0.0005*(double)((values[refIndx_t+1] - values[refIndx_t]) * values[refIndx_vp_eff] * values[refIndx_NG]);
				break;
			  }
			case InversionAnalyser.NET_HARD_KEY:
			  {
				// net-hard = thickness x NG_hard = = (t_below - t_above) x vp_eff x NG_hard
				// times in ms too x factor of 1/2 since times are two-way...
				value = 0.0005*(double)((values[refIndx_t+1] - values[refIndx_t]) * values[refIndx_vp_eff] * values[refIndx_NG_hard]);
				break;
			  }
			case InversionAnalyser.PHIE_KEY:
			  {
				// phi_e = porosity * NG
				value = values[refIndx_phi_s]*values[refIndx_NG];
				break;
			  }
			case InversionAnalyser.NET_HYDROCARBON_KEY:
			case InversionAnalyser.NET_OIL_KEY:
			case InversionAnalyser.NET_LSG_KEY:
			case InversionAnalyser.NET_GAS_KEY:
			  {
				// saturation * porosity x net-sand. Times in ms too x factor of 1/2 since times are two-way...
				value = 0.0005*(double)( values[refIndx_S_h] * values[refIndx_phi_s] * 
										 (values[refIndx_t+1] - values[refIndx_t]) * values[refIndx_vp_eff] * 
										 values[refIndx_NG]);
				// kill it if the fluid is brine
				fluid_type = (short)values[refIndx_fluid_type];
				fluid_type &= 3; // use majority class
				if (fluid_type == FluidStates.FluidTypes.brine) value=0.0;
				break;
			  }
			}
		}
	  else
		{
		  value = (double)values[refIndx0];
		}
	  return value;
	}
  }

  /**
	 Simple class used to store the coefficients we need when extracting a certain
	 realisation (e.g. the median) from a suite of realisations in order to compute
	 its synthetic
  */
  class SyntheticSeismicData
  {
	double[] times;
	double[] refCoeffs;

	// for fluid contacts
	double[] times_fc;
	double[] refCoeffs_fc;

	public SyntheticSeismicData(int n_layers)
	{
	  times = new double[n_layers];
	  refCoeffs = new double[n_layers];
	  // for fluid contacts
	  times_fc = new double[n_layers];
	  refCoeffs_fc = new double[n_layers];
	}
  }

  // Graphics classes
  /**
	 Necessary info for a 2D plot (within a set of plots)
  */
  class GraphingData
  {
	ArrayList<AnPlotInfo>plotList = new ArrayList<AnPlotInfo>();
  }
  
  /** 
	  Data class to collect numbers needed in diverse plots
  */
  class AnPlotInfo
  {
	boolean blockPlot=false;
	float lineWidth=1.0f;
	Color c = null;

	DoubleArrayList x = new DoubleArrayList();
	DoubleArrayList y = new DoubleArrayList();

	DoubleArrayList x1;
	DoubleArrayList y1;
	DoubleArrayList cList;   // primary colour list
	IntArrayList cgList;     // secondary colour-palette list
	String Xword=null;
	String Yword=null;
	// optional specified ranges
	double Ymin = Double.MAX_VALUE;
	double Ymax = Double.MAX_VALUE;

	// offset number (for traces)
	int offsetNumber=0;
	// realisation number
	int rNumber=0;

	public AnPlotInfo() {}

	public AnPlotInfo(boolean _blockPlot) 
	{
	  x1 = new DoubleArrayList();
	  y1 = new DoubleArrayList();
	  cList = new DoubleArrayList();
	  cgList = new IntArrayList();
	  blockPlot = _blockPlot;
	}

	void trim()
	{
	  x.trimToSize();
	  y.trimToSize();
	  if (x1!=null) 
		{
		  x1.trimToSize();
		  y1.trimToSize();
		  cList.trimToSize();
		  cgList.trimToSize();
		}
	}
  }
}
