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

import cern.colt.list.*;
import cern.jet.stat.*;
import cern.colt.matrix.linalg.*;


import static au.csiro.JamesGunning.lib.util.MathUtil.*;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.Numerical.*;
   

/**
   Basic class for performing wavelet extraction.
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/



public final class WaveletExtractor
{
  // basic objects we need
  Model m;
  // Geometry geom;
  // WellLogDefinitions wellLogDefinitions;
  // CheckshotAndMarkerDefinitions chkshotAndMarkerDefinitions;
  GeneralWellInfo[] wellInfo;
  LogData[] logData;
  SegmentedLogData[] sLogData;
  CheckshotAndMarkerData[] chkShotandMarkerData;
  SeismicStackInfo[][] seismicInfo;
  ConvolutionalModel convModel;
  MultiModelSetUp multiModelSetUp;
  UTMConversionInfo utmInfo;
  WaveletPrior waveletPrior;

  /**
	 filename of su output data
  */
  String su_output_filename = "properties.su";

  /**
	 Default filename for su wavelet samples
  */
  String wavelet_realisations_filename = "wavelet_realisations.su";

  /**
	 Default filename for wavelet prior
  */
  String wavelet_prior_filename = null;

  /**
	 Number of wavelet realisations to generate from posterior
  */
  int nSamples=0;

  /**
	 Begininning of window for wavelet estimation
  */
  double tStart;
  /**
	 End of window for wavelet estimation
  */
  double tEnd;
  /** 
	  least acceptable time for a wavlet precursor, seconds
	  Always negative.
   */
  double tPrecursorMin=-0.030;

  /** 
	  greatest acceptable time for a wavlet coda, seconds
   */
  double tCodaMax=0.030;

  /**
	 time discretisation of seismic work (2ms, 4ms etc)
  */
  
  /**
	 SU header word expected for seismic X position
  */
  String seismicXword;
  /**
	 SU header word expected for seismic Y position
  */
  String seismicYword;

  /**
	 SU header word expected for seismic Inline position
  */
  String seismicInlineWord;
  /**
	 SU header word expected for seismic Xline position
  */
  String seismicXlineWord;

  double dt=0.004;   // default for entire extraction. Overridden by seismic minicubes anyway

  ErrorHandler EH;

  /**
	 Parameter object for formatted I/O
  */
    

  // Some basic operating flags which will probably be fixed (#defines)
  boolean TIME_UNITS_AUTO_DETECTION = true;
  /**
	 This will detect velocity comparisons between sonic logs and checkshots that may have mixed ft/s or
	 m/s units
  */
  public static boolean AUTO_DETECT_VELOCITY_UNITS = true;

  boolean useSecantOptimisation=true;  // FIXME add option?????

  
  // Runtime variables

  /**
	 this is the string for the filename of sample model outputs in su form,
	 which take the form realisations.su etc. Traces for a given realisation
	 are tagged by the header word realisationNumberHeaderField.
  */
  public static final String SU_REALISATIONS_FILENAME="realisations.su";
  String su_realisations_filename=SU_REALISATIONS_FILENAME;
  /**
	 Suffix string for blocked logs
  */
  public static String BLOCKED_LOGS_SUFFIX="_blocked.txt";

  /**
	 Switch to allow incomplete log data in estimation interval (substitutes 
	 zero reflections outside, confines likelihood to prior-mean inferred time range)
  */
  public static boolean ALLOW_INCOMPLETE_LOGS = true;

  /**
	 Switch to use full numerical post BFGS Hessian, rather than linearised form
  */
  public static boolean USE_FULL_NUMERICAL_HESSIAN = false;

  /**
	Default maximum number of multi-start global optimiser starting points
  */
  public static int MAX_GLOBAL_OPT_STARTING_POINTS = 40;

  /**
	 minimum number of seismic cube samples for spectral estimation
  */
  public static int MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION = 50;
  /**
	 Default fractional interval velocity constraint
  */
  public static double DEFAULT_INTERVAL_VELOCITY_FRACTIONAL_CONSTRAINT = 0.2;

  String params_xml_filename;


  String commandline;
  int verbosity=0;
  boolean doDumpDebug = false;
  boolean fakeVsLog=false;
  boolean allowLogInterpolation=false;
  double  logInterpolationMaxDist=0;
  boolean eliminateSpikes=false;
  boolean drawMLModel = false;
  boolean monteCarlo = false;
  double monteCarloLateralScale = 10.0;
  boolean useUTMconversion = false;
  boolean fixMinicubeRegistrations = false;
  boolean fixCheckshots = false;
  boolean fixWellPositions = false;
  boolean fixNoise = false;
  double  fixNoiseValue = 0.0;
  boolean constrainPhase = false; 
  double  constantPhase=0;
  double  constantPhaseStdDev=0;
  boolean constrainPeakArrival = false;
  boolean constrainPeakArrivalSign = false;
  boolean peakArrivalPositive = true;
  double  peakArrivalTime=0;
  double  peakArrivalStdDev=0;
  boolean constrainIntervalVelocities = true;
  double  intervalVelocityConstraintStdDev = DEFAULT_INTERVAL_VELOCITY_FRACTIONAL_CONSTRAINT;
  boolean multiStart = false;
  boolean printXsections = false;
  boolean printPrincipalXsections = false;
  boolean printDodgyXsections = false;
  boolean printMinimiserXsections = false;
  boolean useSeed = false;
  String external_wavelet_filename = null;
  boolean BBsynthetic=false;
  boolean makeMLSynthetic = false;
  boolean makeSyntheticAVOgather = false;
  boolean makeMLTiltedSynthetic = false;
  boolean makeCheckshotsAndMarkers = false;
  boolean computeForLongestWaveletOnly = false;
  boolean computeForSelectedWaveletOnly = false;
  boolean showBB = true;
  boolean zeroConstrainPmoments = false;
  int     zeroConstrainPmoments_p=0;
  boolean bbTail = false;
  double  bbTail_T=20;
  double  bbTail_W=100;
  boolean assumeNoReflectionsOutsideLogs = false;
  int selectedWaveletModeNumber = 0;
  boolean computeForCenteredWaveletsOnly = false;
  boolean dumpMLsynthetics = false;
  boolean dumpMLwavelets = false;
  boolean dumpMLparameters = false;
  boolean dumpAnglesToFiles = false;
  boolean ignoreSeismic = false;
  boolean useBIC = false;
  boolean freeAVO = false;
  boolean doShuey = false;
  boolean doIsoZoeppritz = false;
  boolean useStretchedModel = false;
  boolean dumpBlockedLogs = false;
  boolean useMLblocking = false;
  double aggregativeBlockingCutoff = 1.0; 
  boolean doAnisotropicModel = false;
  boolean checkAnisotropicModel = false;
  boolean dumpPandSreflectivity = false;
  boolean synchroniseStackTregistrationErrors = false;
  boolean doBHPcommand=false;
  String BHPcommand_filename = null;
  String runCommand = null;
  boolean useWaveletPrior=false;
  boolean useRestrictedWellSet=false;
  String restrictedWellNumbers=null;
  boolean useRestrictedStackSet=false;
  String restrictedStackNumbers=null;

  boolean useSaturationModel=false;








  boolean useJNI = true;

  /**
	 Default synthetic wavelet bandwidth
  */
  float syntheticFreq = 40.0f;
  /**
	 Default synthetic sampling period (s). 4 ms data
  */
  double syntheticSampleTime=0.004;

  /**
	 Default synthetic noise level (zero)
  */
  double syntheticNoiseLevel=0;

  
  int seed = cern.jet.random.engine.MersenneTwister.DEFAULT_SEED ;

  // Member variables
  SuReader sureader;
  /**
	 seismic amplitude near-stack traces
  */
  SuSegyFormat tr;
  /**
	 wavelet trace
  */
  SuSegyFormat trW;
  /**
	 Array of recent wavelet filenames
  */
  String[] lastWaveletFilenames;
  /**
	Writer for various su output files
  */
  SuWriter suwriter;

  /** 
	output stream for realisations, ML wavelet etc
  */

  FileOutputStream out_strm;

  /**
	 Try to autodetect velocity units. Logs and stacking velocities  must be in same units
	 in AVO reflectivity equation. If true, will try to do this automatically. BEWARE if false
  */
  public static boolean   autoDetectVelocityUnits = true;


  WaveletExtractor()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();

	// set up blas
	int cpu_s = 1;  // FIXME can get this from runtime in 1.4 on (not 1.3 tho')
	SmpBlas.allocateBlas(cpu_s, SeqBlas.seqBlas);
  }


  /** 
	  Sample driver from input streams etc
  */
  public static void main(String argv[])
  {
	WaveletExtractor WE = new WaveletExtractor();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("WaveletExtractor",WE.EH);

	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber);

   	options.addMandatoryArg(1,"file.xml","xml file containing the full extraction specification");	
	options.addOption("-v",null,1,"int","level","Execute with stated level of verbosity");
	options.addOption("-N",null,1,"int","number","Generate N samples of the wavelet from the posterior (default, 0, i.e. no realisations)");
	options.addOption("-b","--dump-blocked-logs",0,null,null,"Dump blocked logs to ascii files",1);
	options.addOption(null,"--ML-blocking",0,null,null,"Use Maximum-likelihood blocking (expensive)",1);
	options.addOption(null,"--blocking-cutoff",1,"float","cutoff","Override aggregative-method blocking cutoff (suggested default 1.0)",1);
	options.addOption("-o",null,1,"String",WE.wavelet_realisations_filename,
					  "Output su file for samples of the wavelet from the posterior "
					  +"(default "+WE.wavelet_realisations_filename+", not written unless -N number provided also)",1);
	options.addOption("--fake-Vs",null,0,null,null,"Form fake Vs log from Vp/"+FundamentalConstants.FAKE_VS_FACTOR+
					  ", or use regression based on useful Vs logs to predict missing data (if only some logs missing)");
	options.addOption("--allow-log-interpolation",null,1,"float","depth-range","permit interpolation of missing logs over a gap less than the range specified",1);
	options.addOption("--eliminate-log-spikes",null,0,null,null,"eliminate spikes in log data",1);
	options.addOption("--anisotropy-test",null,0,null,null,"test anisotropic model for significance",1);
	options.addOption("--seed",null,1,"int","seed value","Seed all internal random generators from this value",1);
	options.addOption(null,"--synthetic-wavelet",1,"String","wavelet.su",
					  "Generate synthetic model(s) (--make-synthetic | --make-synthetic-AVO-gather | --make-tilted-synthetic) using supplied wavelet",1);
	options.addOption(null,"--make-synthetic",3,"float","freq(Hz) sample-interval(ms) noise level (RMS)",
					  "Generate synthetic model with stated seismic wavelet bandwidth, sampling time dt, and noise rms level",1);
	options.addOption(null,"--make-synthetic-AVO-gather",2,"float","freq(Hz) sample-interval(ms)",
					  "Generate synthetic AVO stack with stated seismic wavelet bandwidth, sampling time dt",1);
	options.addOption(null,"--make-tilted-synthetic",3,"float","freq(Hz) sample-interval(ms) noise level (RMS)",
					  "Generate synthetic model with stated seismic wavelet bandwidth, sampling time dt, and noise rms level. "
					  + "Synthetic seismic image has gentle dip, useful for testing purposes",1);
	options.addOption(null,"--BBsynthetic",0,null,null,"--make-synthetic* options use test broadband wavelet",1);
	options.addOption(null,"--make-checkshots-markers",0,null,null,
					  "Using the log and survey data only, make some fake checkshot and marker data",1);
	options.addOption("-l",null,0,null,null,"Perform computations for longest possible wavelet only");
	options.addOption("-c","--centered-only",0,null,null,"Perform computations for centered-span wavelets only");
	options.addOption("-cn",null,1,"int","centered-mode-number",
					  "Perform computations for wavelet with given mode number (as per centered list) only",1);

	if (WE.showBB) options.addOption("-z",null,1,"int","p","Constrain extractions so wavelet--coefficient moments [0..p] are zero",1);
	if (WE.showBB) options.addOption("-t",null,1,"String","T[,weight]","Impose strongly correlated (smooth) tail on wavelet coefficients outside [-T,T]. Default correlation length (samples)="+WE.bbTail_W+". Suggest, e.g. \"-t 20\"",1);
	options.addOption("--dump-ML-synthetics",null,0,null,null,"Dump max-likelihood (ML) synthetics to "+
					  "SU files for all wells + stacks. First trace is true, second is synthetic. Uses "+
					  "seismic minicube time-scale if possible");
	options.addOption("--dump-ML-wavelets",null,0,null,null,"Dump max-likelihood (ML) wavelets to SU files for all stacks");
	options.addOption("--dump-ML-parameters",null,0,null,null,"Dump max-likelihood (ML) parameters to files");
	options.addOption("--dump-PS-reflectivity",null,0,null,null,"Dump P, S and density reflectivities to SU files also, when using --dump-ML-synthetics",1);
	options.addOption("--dump-angles",null,0,null,null,"Dump \"ray-traced\" segment angles to files",1);
	options.addOption("--fix-time-to-depth",null,0,null,null,"Fix time-to-depth parameters at prior-suggested values "+
					  "(in checkshots and registration offsets)");
	options.addOption("--fix-checkshots",null,0,null,null,"Fix time-to-depth parameters in checkshots at prior-suggested values "+
					  "(registration offsets can still vary)");
	options.addOption("--fix-well-positions",null,0,null,null,"Fix well-positioning errors at 0, if these are modelled as unknown");
	options.addOption("--synch-stacks",null,0,null,null,"Synchronise stack registration errors at each well, if these are modelled as unknown");
	options.addOption("--fix-noise",null,1,"float","noise-rms-value",
					  "Fix the noise level at the value given in the estimation (i.e. omit from mode finding). "
					  +"A value of zero or less will default to a rough guess, "
					  +" which is about 1/4 of the rms signal amplitude");
	options.addOption("--constrain-phase",null,2,"float","angle angle-stddev",
					  "Constrain wavelet phase to a constant in pass-band. "
					  +"Angle specifies desired phase (use >360 if unknown - it will float). "
					  +"Angle-std-dev specifies phase-angle uncertainty (std. dev)");
	options.addOption("--constrain-peak-arrival",null,2,"float","time(ms) time-stddev(ms)",
					  "Constrain wavelet peak to arrive at time +/- time-stddev");
	options.addOption("--constrain-peak-arrival-sign",null,3,"float","time(ms) time-stddev(ms) sign",
					  "Constrain wavelet peak of specified sign (use 1 or -1) to arrive at time +/- time-stddev");
	options.addOption("--constrain-interval-velocities",null,1,"float","percentage variation",
					  "Constrain interval velocities: argument is %%-deviation from prior corresponding to one std. "+
					  "deviation (this option is on by default, with 20%% deviation). "+
					  "Setting to 100%% or more will remove velocity constraint altogether.");
	options.addOption("--Xsections",null,0,null,null,"Print log-posterior cross sections at mode points",1);
	options.addOption("--Xsections-principal",null,0,null,null,
					  "Print log-posterior cross sections along principal axes at mode points",1);
	options.addOption("--Xsections-dodgy",null,0,null,null,"Print log-posterior cross sections of 'dodgy' cross sections at mode points",1);
	options.addOption("--Xsections-minimiser",null,0,null,null,"Print log-posterior cross sections pre and post minimiser",1);
	options.addOption("-IS",null,0,null,null,"Ignore seismic data (effectively: optimise/sample from prior)",1);
	options.addOption("-NLR",null,0,null,null,"Assume no reflections outside log data scope");
	options.addOption("--BIC",null,0,null,null,"Use approximate Bayesian Information Criterion (BIC) instead of formal Bayes factors",1);
	options.addOption("--free-AVO",null,0,null,null,"Make AVO B coefficient free for each stack",1);
	options.addOption("--Shuey",null,0,null,null,"Use Shuey equation for AVO",1);
	options.addOption("--Zoeppritz",null,0,null,null,"Use isotropic Zoeppritz equation for AVO",1);
	options.addOption("--multi-start",null,0,null,null,"Try multi-start optimiser points for multimodal problems",1);
	options.addOption("--dump-debug-package",null,0,null,null,
					  "dump debug zip package containing all needed files for maintainers");	
	options.addOption("--monte-carlo",null,1,"float","correlation-length(traces)",
					  "Perform specialised parameter uncertainty study by MC generation of fake seismic data sets",1);
	options.addOption("--BHPcommand",null,1,"String","filename",
					  "write to filename the string \"properties=x,x..\" "+
					  "(describes the properties in the layered-model-like properties.su file)",1);
	options.addOption("--prior",null,1,"String","MLparameters-file",
					  "Use MLparameters file from a previous extraction to set the wavelet prior",1);
	options.addOption("--wells",null,1,"String","comma-separated-well-numbers",
					  "restrict extraction to subset of wells in XML file (numbered 1,2..). E.g. --wells 1,3,5");
	options.addOption("--stacks",null,1,"String","comma-separated-stack-numbers",
					  "restrict extraction to subset of stacks in XML file (numbered 1,2..). E.g. --stacks 1,2");
	options.addOption(null,"--no-jni",0,null,null,"Disable native libraries (use java fallbacks)",1);

	// FIXME RESTORE, see also occurrence below
	//options.addOption("-D",null,0,null,null,"Secret time-lapse mode ;-) ",1); 

	// leave this to last!
	options.addExtraHelpOption("-X","Show help with all options (there are lots more!)");

 

	options.parseOptions(argv,true);
 

	// Save commandline 
	StringBuffer cl = new StringBuffer();
	for (int i=0;i<argv.length;i++) cl.append(argv[i]+" ");
	WE.commandline = cl.toString();

	// stuff commandline options into Seismicinverter fields
	WE.params_xml_filename = options.getMandatoryArg(1);

	

	if (options.isOptionSet("-X")) options.Usage();

	if (options.isOptionSet("--dump-debug-package"))
	  { 
		WE.doDumpDebug = true;
	  }

	if (options.isOptionSet("-v")) 
	  WE.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();

	if (options.isOptionSet("-N")) 
	  WE.nSamples = ((Integer)options.getOptionValues("-N")[0]).intValue();

	if (options.isOptionSet("--dump-PS-reflectivity")) 
	  WE.dumpPandSreflectivity = true;

	if (options.isOptionSet("-b")) 
	  WE.dumpBlockedLogs=true;

	if (options.isOptionSet("--ML-blocking")) 
	  WE.useMLblocking=true;

	if (options.isOptionSet("--blocking-cutoff")) 
	  {
		WE.aggregativeBlockingCutoff = (double)((Float)options.getOptionValues("--blocking-cutoff")[0]).floatValue();
		if (WE.aggregativeBlockingCutoff <=0)
		  {	
			WE.EH.die_err("--blocking-cutoff cutoff < 0 makes no sense. Require cutoff > 0: cutoff around 1 is sensible");
		  }
	  }
	
	if (options.isOptionSet("-o")) 
	  WE.wavelet_realisations_filename = ((String)options.getOptionValues("-o")[0]);

	if (options.isOptionSet("--prior")) 
	  {
		WE.useWaveletPrior = true;
		WE.wavelet_prior_filename = ((String)options.getOptionValues("--prior")[0]);
	  }

	if (options.isOptionSet("--seed"))
	  {
		WE.useSeed = true;
		WE.seed = ((Integer)options.getOptionValues("--seed")[0]).intValue();
	  }
	if (options.isOptionSet("--allow-log-interpolation"))
	  {
		WE.allowLogInterpolation = true;
		WE.logInterpolationMaxDist = (double)((Float)
											  options.getOptionValues("--allow-log-interpolation")[0]).floatValue();
	  }
	if (options.isOptionSet("--eliminate-log-spikes"))
	  {
		WE.eliminateSpikes = true;
	  }
	if (options.isOptionSet("--anisotropy-test"))
	  {
		WE.checkAnisotropicModel = true;
	  }
	if (options.isOptionSet("--fake-Vs"))
	  {
		WE.fakeVsLog = true;
	  }
	if (options.isOptionSet("--make-checkshots-markers"))
	  {
		WE.makeCheckshotsAndMarkers = true;
	  }
	if (options.isOptionSet("--BBsynthetic"))
	  {
		WE.BBsynthetic = true;
	  }
	if (options.isOptionSet("--make-synthetic"))
	  {
		WE.makeMLSynthetic = true;
		// set synthetic freq
		WE.syntheticFreq = ((Float)(options.getOptionValues("--make-synthetic")[0])).floatValue();
		// and dt. Convert to seconds
		WE.syntheticSampleTime = 1.0e-3*((Float)(options.getOptionValues("--make-synthetic")[1])).floatValue();
		WE.syntheticNoiseLevel = ((Float)(options.getOptionValues("--make-synthetic")[2])).floatValue();
		WE.assumeNoReflectionsOutsideLogs = true; // good safe thing to do
	  }
	if (options.isOptionSet("--make-tilted-synthetic"))
	  {
		WE.makeMLTiltedSynthetic= true;
		WE.makeMLSynthetic = false;
		// set synthetic freq
		WE.syntheticFreq = ((Float)(options.getOptionValues("--make-tilted-synthetic")[0])).floatValue();
		// and dt. Convert to seconds
		WE.syntheticSampleTime = 1.0e-3*((Float)(options.getOptionValues("--make-tilted-synthetic")[1])).floatValue();
		WE.syntheticNoiseLevel = ((Float)(options.getOptionValues("--make-tilted-synthetic")[2])).floatValue();
		WE.assumeNoReflectionsOutsideLogs = true; // good safe thing to do
	  }
	if (options.isOptionSet("--make-synthetic-AVO-gather"))
	  {
		WE.makeMLSynthetic = false;
		WE.makeMLTiltedSynthetic= false;
		WE.makeSyntheticAVOgather = true;
		// set synthetic freq
		WE.syntheticFreq = ((Float)(options.getOptionValues("--make-synthetic-AVO-gather")[0])).floatValue();
		// and dt. Convert to seconds
		WE.syntheticSampleTime = 1.0e-3*((Float)(options.getOptionValues("--make-synthetic-AVO-gather")[1])).floatValue();
		WE.assumeNoReflectionsOutsideLogs = true; // good safe thing to do
	  }
	if (options.isOptionSet("--synthetic-wavelet"))
	  {
		WE.external_wavelet_filename = (String)options.getOptionValues("--synthetic-wavelet")[0];
	  }
	if (options.isOptionSet("-l"))
	  {
		WE.computeForLongestWaveletOnly = true;
	  }
	if (options.isOptionSet("-c"))
	  {
		WE.computeForCenteredWaveletsOnly = true;
	  }
	if (options.isOptionSet("-cn"))
	  {
		WE.computeForSelectedWaveletOnly = true;
		WE.selectedWaveletModeNumber = ((Integer)options.getOptionValues("-cn")[0]).intValue();
	  }
	if (WE.showBB && options.isOptionSet("-z"))
	  {
		WE.zeroConstrainPmoments = true;
		WE.zeroConstrainPmoments_p = ((Integer)options.getOptionValues("-z")[0]).intValue();
		if(WE.zeroConstrainPmoments_p<0 || WE.zeroConstrainPmoments_p>4) 
		  {
			WE.EH.warn("Constraining p-moments needs 0 <= p <= 4");
			options.Usage();
		  }
	  }
	if (WE.showBB && options.isOptionSet("-t"))
	  {
		WE.bbTail = true;
		String s = options.getOptionString("-t",0);
		String[] sa = s.split(",");     // split ,
		try
		  {
			WE.bbTail_T = Double.parseDouble(sa[0]);
			if(sa.length>1) WE.bbTail_W = Double.parseDouble(sa[1]);
		  }
		catch (NumberFormatException e)
		  {
			WE.EH.die_err("Option -t "+s+" must be -t number[,number]");
		  }
		if(	WE.bbTail_T <= 0) WE.EH.die_err("Option -t "+s+" must be a T[,weight], with T>0");
		if(	WE.bbTail_W <= 0) WE.EH.die_err("Option -t "+s+" must be a T[,weight], with weight>0");

	  }
	if (options.isOptionSet("--dump-ML-synthetics"))
	  {
		WE.dumpMLsynthetics = true;
	  }
	if (options.isOptionSet("--dump-ML-wavelets"))
	  {
		WE.dumpMLwavelets = true;
	  }
	if (options.isOptionSet("--dump-ML-parameters"))
	  {
		WE.dumpMLparameters = true;
	  }
	if (options.isOptionSet("--dump-angles"))
	  {
		WE.dumpAnglesToFiles = true;
	  }
	if (options.isOptionSet("--multi-start"))
	  {
		WE.multiStart = true;
	  }
	if (options.isOptionSet("--monte-carlo"))
	  {
		WE.monteCarlo = true;
		WE.monteCarloLateralScale = ((Float)(options.getOptionValues("--monte-carlo")[0])).floatValue();
		if (WE.monteCarloLateralScale < 0.1 || WE.monteCarloLateralScale > 100)
		  {
			WE.EH.warn("Monte-Carlo correlation length probably unrealistic. Try [0.1 ... 100]");
			options.Usage();
		  }
	  }
	if (options.isOptionSet("--fix-noise"))
	  {
		WE.fixNoise = true;
		WE.fixNoiseValue = ((Float)(options.getOptionValues("--fix-noise")[0])).floatValue();
	  }
	if (options.isOptionSet("--constrain-phase"))
	  {
		WE.constrainPhase = true;
		WE.constantPhase = ((Float)(options.getOptionValues("--constrain-phase")[0])).floatValue();
		WE.constantPhaseStdDev = Math.abs(((Float)(options.getOptionValues("--constrain-phase")[1])).floatValue());
	  }	
	if (options.isOptionSet("--constrain-peak-arrival"))
	  {
		WE.constrainPeakArrival = true;
		WE.peakArrivalTime = 1e-3*((Float)(options.getOptionValues("--constrain-peak-arrival")[0])).floatValue();
		WE.peakArrivalStdDev = 1e-3*Math.abs(((Float)(options.getOptionValues("--constrain-peak-arrival")[1])).floatValue());
		if (WE.peakArrivalStdDev  < 1e-4)
		  {
			WE.EH.warn("Peak arrival tolerances smaller than 0.1ms are illegal");
			options.Usage();
		  }
	  }	
	if (options.isOptionSet("--constrain-peak-arrival-sign"))
	  {
		WE.constrainPeakArrival = true;
		WE.constrainPeakArrivalSign = true;
		Float[] ops = (Float[])(options.getOptionValues("--constrain-peak-arrival-sign"));
		WE.peakArrivalTime = 1e-3*ops[0].floatValue();
		WE.peakArrivalStdDev = 1e-3*Math.abs(ops[1].floatValue());
		WE.peakArrivalPositive = (ops[2].floatValue() > 0 );

		if (WE.peakArrivalStdDev  < 1e-4)
		  {
			WE.EH.warn("Peak arrival tolerances smaller than 0.1ms are illegal");
			options.Usage();
		  }

	  }	
	if (options.isOptionSet("--constrain-interval-velocities"))
	  {
		// this is already set
		// WE.constrainIntervalVelocities = true;
		// convert % to true number
		WE.intervalVelocityConstraintStdDev 
		  = 0.01*((Float)(options.getOptionValues("--constrain-interval-velocities")[0])).floatValue();
		if (WE.intervalVelocityConstraintStdDev < 0)
		  {
			WE.EH.warn("Error: Interval velocity constraint must be > 0 ");
			options.Usage();
		  }
		if (WE.intervalVelocityConstraintStdDev == 0)
		  {
			// take this to mean  fixed t2d in checkshots
			WE.fixCheckshots = true;
			WE.constrainIntervalVelocities = false;
		  }
		if (WE.intervalVelocityConstraintStdDev < 0.02)
		  {
			WE.EH.warn("WARNING: constraining interval velocities to less than 2%% may be unstable");
		  }
		if ((WE.intervalVelocityConstraintStdDev > 0.3) && (WE.intervalVelocityConstraintStdDev < 1.0))
		  {
			WE.EH.warn("WARNING: constraining interval velocities by more than 30%% may be unphysical");
		  }
		if (WE.intervalVelocityConstraintStdDev >= 1.0)
		  {
			WE.EH.println("Switching off interval velocity constraints");
			WE.constrainIntervalVelocities = false;
		  }
	  }
	if (options.isOptionSet("--Xsections"))
	  {
		WE.printXsections= true;
	  }
	if (options.isOptionSet("--Xsections-principal"))
	  {
		WE.printPrincipalXsections= true;
	  }
	if (options.isOptionSet("--Xsections-dodgy"))
	  {
		WE.printDodgyXsections= true;
	  }
	if (options.isOptionSet("--Xsections-minimiser"))
	  {
		WE.printMinimiserXsections= true;
	  }
	if (options.isOptionSet("-IS"))
	  {
		WE.ignoreSeismic = true;
	  }
	if (options.isOptionSet("--BIC"))
	  {
		WE.useBIC = true;
	  }
	if (options.isOptionSet("--free-AVO"))
	  {
		WE.freeAVO = true;
	  }
	if (options.isOptionSet("--Shuey"))
	  {
		WE.doShuey = true;
	  }
	if (options.isOptionSet("--Zoeppritz"))
	  {
		WE.doShuey = false; // can't have Shuey
		WE.doIsoZoeppritz = true;
	  }
	if (options.isOptionSet("--fix-time-to-depth"))
	  {
		// fix both
		WE.fixMinicubeRegistrations = true;
		WE.fixCheckshots = true;
	  }

	if (options.isOptionSet("--fix-checkshots"))
	  {
		WE.fixCheckshots = true;
	  }

	if (options.isOptionSet("--fix-well-positions"))
	  {
		WE.fixWellPositions = true;
	  }
	if (options.isOptionSet("--synch-stacks"))
	  {
		WE.synchroniseStackTregistrationErrors = true;
	  }

	if (options.isOptionSet("-NLR"))
	  {
		WE.assumeNoReflectionsOutsideLogs = true;
	  }

	// BHP commands
	if (options.isOptionSet("--BHPcommand"))
	  {
		WE.doBHPcommand = true;
		WE.BHPcommand_filename = (String)options.getOptionValues("--BHPcommand")[0];
	  }


	if (!WE.computeForLongestWaveletOnly && WE.useWaveletPrior)
	  {
		WE.EH.warn("Error: Wavelet --prior must be used with -l (longest-only wavelet option)");
		options.Usage();		
	  }

	// Restricted wells
	if (options.isOptionSet("--wells"))
	  {
		WE.useRestrictedWellSet = true;
		WE.restrictedWellNumbers = (String)options.getOptionValues("--wells")[0];
	  }
	// Restricted stacks
	if (options.isOptionSet("--stacks"))
	  {
		WE.useRestrictedStackSet = true;
		WE.restrictedStackNumbers = (String)options.getOptionValues("--stacks")[0];
	  }
	if (options.isOptionSet("--no-jni")) WE.useJNI = false;
	//if (options.isOptionSet("-D")) WE.useSaturationModel = true;

	// Store the command string
	WE.runCommand = options.getFullArgsString();

	// Show command for posterity
	String classname = ((new Throwable()).getStackTrace())[0].getClassName();
	WE.EH.println("java "+classname+" "+WE.runCommand+"%n");

	// Do the extraction !
	if (!WE.runWaveletExtractor())
	  {
		options.Usage();
	  }
  }


  /**
	 Main sequence of calls to load data and perform extraction
  */
  boolean runWaveletExtractor()
  {
	boolean success=false;

	// allocate all the basic objects we need
	m = new Model(this);
	convModel = new ConvolutionalModel(this);

	// (1) Parse the XML file
	ParameterFileParser MDP = new ParameterFileParser(params_xml_filename, this);
	if (!MDP.readParameters())
	  {
		return false;
	  }
	if (verbosity>0) 
	  {
		MDP.summariseParameters();
	  }

	if (doDumpDebug) 
	  {
		// Dump debug package here, if required, then bail out
		formWavExDebugPackage(null);
		return true;
	  }
	try
	  {
		setUpInput();
		setUpOutput();
		loadData();
		setUpPropertiesOutput();
		setUpWellInterpolations();
		setUpWaveletModel();
		integrateCheckShotAndMarkerData();
		segmentLogsAndCastToTrueDepth();
		createInitialTimeDepthMap();
		setUpSaturationAdjustmentModel();
		getWaveletPriorIfNeeded();

		performExtraction();
		finalizeOutputs();
	  }
	catch (Exception ex)
	  {
		EH.println("Wavelet extraction failed to finish. Maybe do a --dump-debug-package for maintainers?");
		ex.printStackTrace();
	  }
	return true;
  }

  void formWavExDebugPackage(Exception ex)
  {
	int i,j;
	DebugPackageDumper dumper;


	ArrayList ordinaryFilesAL;
	ArrayList accessibleSuFilesAL;
	String[] ordinaryFiles=null;
	String[] accessibleSuFiles=null;
	DebugPackageDumper.SuTraceInfo[] traceInfo=null;

	// ordinary files
	ordinaryFilesAL = new ArrayList();

	// XML
	ordinaryFilesAL.add(params_xml_filename);

	// logs, surveys, checkshots, markers
	for (i=0; i<m.nWells; i++)
	  {
		if (wellInfo[i].surveyFilename != null) 
		  ordinaryFilesAL.add(wellInfo[i].surveyFilename);
		if (wellInfo[i].logFilename != null) 
		  ordinaryFilesAL.add(wellInfo[i].logFilename);
		if (wellInfo[i].checkshotFilename != null) 
		  ordinaryFilesAL.add(wellInfo[i].checkshotFilename);
		if (wellInfo[i].markersFilename != null) 
		  ordinaryFilesAL.add(wellInfo[i].markersFilename);
	  }
	ordinaryFilesAL.trimToSize();
	ordinaryFiles = (String[]) (ordinaryFilesAL.toArray(new String[ordinaryFilesAL.size()]));

	// su files
	accessibleSuFilesAL = new ArrayList();
	for (i=0; i<m.nWells; i++)
	  {
		for (j=0; j<m.nStacks; j++)
		  {
			if (seismicInfo[i][j].filename != null) 
			  accessibleSuFilesAL.add(seismicInfo[i][j].filename);
		  }
	  }
	accessibleSuFilesAL.trimToSize();

	accessibleSuFiles = (String[]) accessibleSuFilesAL.toArray(new String[accessibleSuFilesAL.size()]);



	dumper = new DebugPackageDumper(EH, ordinaryFiles,accessibleSuFiles,
									null,commandline,ex);


	dumper.dumpPackage();	
  }

 
  private void loadData()
  {
	int i,j;


	// loop over wells, reading in survey, log, checkshot, marker and seismic data
	logData = new LogData[m.nWells];
	chkShotandMarkerData = new CheckshotAndMarkerData[m.nWells];


	// logs and surveys
	if (verbosity>0) EH.println("Reading survey and log data...");
	for (i=0; i<m.nWells; i++)
	  {
		if (verbosity>0) EH.println("  For well: "+ wellInfo[i].well_name);
		logData[i] = new LogData(EH);
		if (verbosity>0) EH.println("    Reading survey data: "+wellInfo[i].surveyFilename.replaceAll("%","%%"));
		readSurveys(wellInfo[i],logData[i]);
		// Quick sanity checks on survey
		logData[i].survey.checkSurveyDataForSanity(EH);

		if (verbosity>0) EH.println("    Reading log data: "+wellInfo[i].logFilename.replaceAll("%","%%"));
		readLogs(wellInfo[i],logData[i]);
	  }
	// generate fake Shear logs if required
	if (fakeVsLog)
	  {
		LogPredictor LP = new LogPredictor(this);
		LP.formPredictedShearData();
	  }
	
		



	// checkshots 
	if (makeCheckshotsAndMarkers)
	  {
		// early call of this
		segmentLogsAndCastToTrueDepth();
	  }
	if (verbosity>0) EH.println("Reading checkshots and markers...");
	for (i=0; i<m.nWells; i++)
	  {
		if (verbosity>0) EH.println("  For well: "+ wellInfo[i].well_name);
		chkShotandMarkerData[i]= new CheckshotAndMarkerData(this,i);
		if (makeCheckshotsAndMarkers)
		  {
			// create and write fake data
			wellInfo[i].checkshotFilename = ("fake_checkshot."+(i+1)+".txt");
			chkShotandMarkerData[i].GenerateFakeCheckshotdata(sLogData[i],wellInfo[i]);
		  }

		if (verbosity>0) EH.println("    Reading checkshot data: "+wellInfo[i].checkshotFilename.replaceAll("%","%%"));
		readCheckShots(wellInfo[i],chkShotandMarkerData[i]);
		// automake checkshots if suggested
		if (!makeCheckshotsAndMarkers && wellInfo[i].checkshotAutoMake)
		  {
			// Take first point in checkshot-file as reference point only
			// checkshot "suspended" from this
			chkShotandMarkerData[i].autoGenerateCheckshots(logData[i],wellInfo[i]);
		  }
	  }
	// markers
	for (i=0; i<m.nWells; i++)
	  {
		if (verbosity>0) EH.println("  For well: "+ wellInfo[i].well_name);
		if (makeCheckshotsAndMarkers)
		  {
			// create and write fake data
			wellInfo[i].markersFilename = ("fake_markers."+(i+1)+".txt");
			chkShotandMarkerData[i].GenerateFakeMarkerdata(sLogData[i],wellInfo[i]);
		  }
		if (wellInfo[i].markersFilename != null)
		  {
			if (verbosity>0) EH.println("    Reading marker data: "+wellInfo[i].markersFilename);
			readMarkers(wellInfo[i],chkShotandMarkerData[i]);
		  }
		// Quick sanity checks on checkshots and markers
		chkShotandMarkerData[i].checkCheckshotAndMarkerDataForSanity(logData[i],wellInfo[i]);
	  }

	// seismic data will be loaded here.......
	if (!(makeMLSynthetic || makeMLTiltedSynthetic || makeSyntheticAVOgather))
	  {
		if (verbosity>0) EH.println("Reading seismic data...");
		for (i=0; i<m.nWells; i++)
		  {
			if (verbosity>0) EH.println("  For well: "+ wellInfo[i].well_name);
			for (j=0; j<m.nStacks; j++)
			  {
				// ......certain useful stuff is extracted (eg dt) here etc
				seismicInfo[i][j].loadSeismicData(this,i);
			  }
		  }
	  }
	else
	  {
		// set sample times etc
		dt = syntheticSampleTime;
		if (dt<=0) EH.die_err("sample time dt illegal"); 
	  }


  }


  /**
	 Method to set up properties.su file and handle BHPcommand.
	 This is a hook to enable external command scripts: typically a command that requires
	 the string "properties=A,B,C,.." is appended to something like BHPwrite.
	 This method writes the string "properties=A,B,C,.." to BHPcommand_filename, so it can be
	 used in scripts using e.g. backticking.
	 A,B,C etc are the properties that follow the order in 
	 Model.writeOutStochasticProperties()
  */
  private void setUpPropertiesOutput() 
  {
	PrintWriter out=null;

	// we always need this, even if doBHPcommand is false;
	String[] stringList = m.getPropertiesStringList();

	if (doBHPcommand) 
	  {
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(BHPcommand_filename))); 
		  }
		catch (IOException e) 
		  {
			EH.die_err("Can't open file "+BHPcommand_filename); 
		  }
		
		// Append the comma separated list of properties
		StringBuffer property_string = new StringBuffer("properties=");
		// Leave out t & t-base here, (always last entry in property_string)
		for (int i=0;i<stringList.length;i++) 
		  {
			property_string.append(stringList[i]);
			if (i<stringList.length-1) property_string.append(",");	
		  }
		out.printf("%s %n",property_string);
		out.close();
	  }
  }








  // we need a bunch of splines at each well
  private void setUpWellInterpolations()
  {
	double Zm_min, Zm_max;

	for (int i=0; i<m.nWells; i++)
	  {
		logData[i].survey.setUpXYSplines(this.EH,chkShotandMarkerData[i]);
	  }
  }

  private void setUpWaveletModel()
  {
	// then complete initialisation of convolutional model
	convModel.initialiseConvolutionalModel();
  }



  // integrateCheckShotAndMarkerData
  private void integrateCheckShotAndMarkerData()
  {
	if (verbosity>0) 
	  {
		EH.printf("%nAmalgamating checkshot and marker data...%n");
	  }
	for (int i=0; i<m.nWells; i++)
	  {
		chkShotandMarkerData[i].setUpCheckShotAndMarkerData(logData[i].survey);
	  }
  }

  // set up time to depth maps
  private void createInitialTimeDepthMap()
  {
	if (verbosity>0) EH.println("Setting up time to depth map...");
	m.timeToDepthMaps = new TimeToDepthMap[m.nWells];
	for (int i=0; i<m.nWells; i++)
	  {
		m.timeToDepthMaps[i] = new TimeToDepthMap(this, chkShotandMarkerData[i].t2d,i);
	  }
  }


 // Set up useful input streams
  private void setUpInput()
  {
	tr = new SuSegyFormat();
	trW = new SuSegyFormat();
	sureader = new SuReader();
  }

  // Set up useful output streams
  private void setUpOutput()
  {
	suwriter = new SuWriter();
	if (!(makeMLSynthetic || makeMLTiltedSynthetic || makeSyntheticAVOgather))
	  {
		try 
		  {
			out_strm = new FileOutputStream(su_output_filename);
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open" + su_output_filename);
		  }
	  }
  }

  // clean up all streams etc opened in setUpOutput()
  private void finalizeOutputs()
  {
	if (!(makeMLSynthetic || makeMLTiltedSynthetic || makeSyntheticAVOgather))
	  {
		try
		  {
			suwriter.flush(out_strm); // must do this!
			out_strm.close();
			EH.println("Wrote su data to: " + su_output_filename);
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close" + su_output_filename);
		  }
	  }
  }


  /**
	 Read info from wavelet prior if needed
  */
  void getWaveletPriorIfNeeded()
  {
	if (useWaveletPrior)
	  {
		waveletPrior = new WaveletPrior(this);
	  }
  }

  /**
	 Encoding string for wavelets in MLgraphing file
  */
  String waveletFilenamesEncodingString()
  {
	if (lastWaveletFilenames==null) return "";

	StringBuffer code = new StringBuffer(" Wavelets:[");
	int nWavelets = useStretchedModel ? m.nStacks : 1;
	for (int i=0; i<nWavelets; i++)
	  {
		code.append(lastWaveletFilenames[i]);
		if (i!=(nWavelets-1)) code.append(",");
	  }
	code.append("]");
	return code.toString();
  }

  /**
	 Do the required things
  */ 
  void performExtraction()
  {
	if (makeMLSynthetic || makeMLTiltedSynthetic || makeSyntheticAVOgather)
	  {
		makeSynthetic();
	  }
	else
	  {
		// Actual extractions
		multiModelSetUp = new MultiModelSetUp(this);

		// set up the prior model
		if (verbosity>0) EH.println("Setting up prior model...");
		multiModelSetUp.setUpPriorModel();	


		// (2) Find modes
		if (verbosity>0) EH.println("Finding modes...");
		multiModelSetUp.findModes();

		// (3) Summary output
		multiModelSetUp.outputPosteriorEntities();
		
		// (4) possible monteCarlo
		multiModelSetUp.performMCcalculations();

		// Due to possible monteCarlo, the in-situ ML wavelet 
		// can no longer be assumed to be the ML wavelet from here......

		// (5) Stochastic samples if required
		makeSamples();

		// instrumentation (usually silent)
		convModel.writeOptimiserTimings();
		
	  } 
  }
  /**
	Make a synthetic using Maximum likelihood values from the prior, but a standard
	Ricker wavelet
  */
  void makeSynthetic() 
  {
	if (verbosity>0) EH.println("Making synthetic...");
	m.checkTimeToDepthRanges(); // Will extrapolate OK if checkshot too narrow
	m.checklogDataRanges();     // I think this is OK if logs run out....

	// Anisotropy common to all fwd synthetic models
	m.petro.hasAnisotropy = doAnisotropicModel;

	if (makeMLSynthetic || makeMLTiltedSynthetic)
	  {
		if (useSaturationModel) m.recomputeLogProperties(false);
		m.computeReflectivities(false,false);
		if (external_wavelet_filename==null)
		  {
			m.wavelet = new Wavelet(this, syntheticFreq, syntheticSampleTime, !BBsynthetic /* default ricker2 */);
		  }
		else
		  {
			m.wavelet = new Wavelet(this,external_wavelet_filename);
		  }
		if (useSaturationModel) 
		  {
			// copy reflectivities
			m.makeReflectivityCopy();
			// compute effect of sat change,
			m.recomputeLogProperties(true);
			// new reflectivities
			m.computeReflectivities(false,false);
			// then difference
			m.computeReflectivityDifference();
		  }

		m.addNoiseToReflectivity();
		convModel.computeSyntheticSeismic();
		convModel.writeOutFullSyntheticSUset(makeMLTiltedSynthetic);
		if (dumpMLwavelets)
		  {
			m.wavelet.writeOutSUWavelet("MLwavelet.su",false,false,0,0,false,0,false);
		  }
	  }
	else if (makeSyntheticAVOgather)
	  {
		convModel.writeOutAVOSyntheticSUset();
		m.wavelet.writeOutSUWavelet("SynthWavelet.su",false,false,0,0,false,0,false);
	  }
  }


  void makeSamples()
  {
	// (3) sample...
	if (verbosity>0) EH.println("Generating wavelet samples...");
	multiModelSetUp.sampleWavelets();
  }

  /**
	 Here we will segment the logs on the basis of the impedance, and
	 generate Vp,Vs,rho data at the layer tops, with the equivalent true depth
	 filled in as well
  */
  private void segmentLogsAndCastToTrueDepth()
  {
	WellSurveyDepthConverter converter;
	// loop over wells, upscaling data

	sLogData = new SegmentedLogData[m.nWells];

	for (int i=0; i<m.nWells; i++)
	  {
		if (verbosity>0) EH.println("Checking log data integrity...");
		// simple checks
		if (chkShotandMarkerData[i].Zchm[0] >logData[i].Zm[0])
		  {
			EH.println("WARNING! Well log data starts before first checkshot point");
		  }
		if (chkShotandMarkerData[i].Zchm[chkShotandMarkerData[i].Zchm.length-1] 
			< logData[i].Zm[logData[i].Zm.length-1])
		  {
			EH.println("WARNING! Well log extends below last checkshot point");
		  }
		boolean skipSubselect; 
		skipSubselect = wellInfo[i].checkshotAutoMake && wellInfo[i].checkshotAutoMakeModelSelect;
		if (skipSubselect)
		  {
			EH.println("FIXME: do log reduction consistently with with segIndexValuesArray realignment");
		  }
		else
		  {
			logData[i].cleanUpDataAndCheckIntegrity(logData[i].ZmMin,logData[i].ZmMax);
		  }

		if (verbosity>0) EH.println("Segmenting log data...well "+wellInfo[i].logFilename.replaceAll("%","%%"));

		sLogData[i] = new SegmentedLogData(wellInfo[i],EH,verbosity);

		// pass in first stack to get (for debugging only) idea of angles
		sLogData[i].segmentLogSuite(this, logData[i],
									useMLblocking, dumpPandSreflectivity, seismicInfo[0][0].getVdenom());
		
		if (autoDetectVelocityUnits)
		  {
			// do this based on first stack. Ought to be good enough
			sLogData[i].autoDetectVunits(this,seismicInfo[i][0].V_stack, 
										 chkShotandMarkerData[i].meanNaturalVelocity);
		  }
		// make sure logs have true depth available:
		converter = new WellSurveyDepthConverter(logData[i].survey,EH);
		sLogData[i].Z = converter.getTrueVerticalDepths(sLogData[i].Zm);
		// check units
		
	  }
  }



  private void readSurveys(GeneralWellInfo info, LogData logData)
  {
	WellSurveyReader reader;
	
	// FIXME this is format dependent
	reader = new GEOEASWellSurveyReader();

	// set the header words we need
	reader.setSurveyDataKeyword(KeyWords.Zm, info.geom.Zm);
	reader.setSurveyDataKeyword(KeyWords.Z, info.geom.Z);
	reader.setSurveyDataKeyword(KeyWords.X, info.geom.X);
	reader.setSurveyDataKeyword(KeyWords.Y, info.geom.Y);

	// get the data
	reader.extractSurveyData(info.surveyFilename, logData, EH);
  }

  private void readLogs(GeneralWellInfo info, LogData logData)
  {
	WellLogReader reader;
	double z1,z2;	
	// FIXME this is format dependent
	reader = new GEOEASWellLogReader();

	// set the header words we need
	reader.setLogDataKeyword(KeyWords.Zm, info.wellLogDefinitions.MD);
	reader.setLogDataKeyword(KeyWords.DT, info.wellLogDefinitions.DT);
	if (info.wellLogDefinitions.shearTreatment == 1)
	  {
		// look for this, we need it
		reader.setLogDataKeyword(KeyWords.DTshear, info.wellLogDefinitions.DTshear);
	  }
	reader.setLogDataKeyword(KeyWords.rho, info.wellLogDefinitions.rho);
	// Optional gamma
	if (info.wellLogDefinitions.gamma != null)
	  {
		reader.setLogDataKeyword(KeyWords.gamma, info.wellLogDefinitions.gamma);
	  }

	// get the data
	if (!reader.extractData(info.logFilename, logData, EH))
	  {
		EH.die_err("Unsuccessful parse of file "+info.logFilename);
	  }

	// check for monotonically increasing depths?
	if (logData.Zm[1] < logData.Zm[0])
	  {
		EH.println("   LAS data is in reverse depth order. Reverse it to save time!");
		logData.reverseDepthOrderOfData();
	  }
	z1=info.wellLogDefinitions.MD_uppermute;
	z2=info.wellLogDefinitions.MD_lowermute;
	if(z1!= Double.NaN  || z2 != Double.NaN)
	  {
		if (z1==Double.NaN) z1 = -Double.MAX_VALUE;
		if (z2==Double.NaN) z2 =  Double.MAX_VALUE;
		EH.warn("Muting data in file %s outside [%f,%f]%n",info.logFilename,z1,z2);
		logData.muteData(z1,z2);
	  }

	// generate fake Shear log if required
	if (fakeVsLog)
	  {
		/*
		logData.DTshear = (double[])logData.DT.clone();
		if (true)
		  {
			EH.println("FIXME Vs REGRESSION");
			logData.DTshear = (double[])logData.DT.clone();
			for (int i =0; i< logData.DT.length ; i++) 
			  {
				logData.DTshear[i] = FundamentalConstants.FAKE_VS_FACTOR * logData.DT[i];
				logData.DTshear[i] =  3.43 * logData.DT[i] -136.108;
			  }
		  }
		*/
	  }
	// mark gamma status
	if (logData.gamma == null || (logData.gamma.length == 0))
	  {
		logData.haveGamma = false;
	  }
	// check if gamma presence is consistent with specified shale thresholds
	if (doAnisotropicModel) 
	  {		
		if (logData.haveGamma)
		  {
			if (!m.petro.GRThreshProvided)
			  {
				EH.println("   Anisotropy model warning: have gamma log but no XML petrophysics_info/shale_threshold/GR_threshold entry");
			  }
		  }
		else
		  {
			if (!m.petro.densityThreshProvided)
			  {
				EH.die_err("   Anisotropy model failure: have no gamma log and no petrophysics_info/shale_threshold/density_threshold");
			  }
		  }
	  }
  }

  private void readCheckShots(GeneralWellInfo info,CheckshotAndMarkerData chkShotandMarkerData)
  {
	CheckshotReader reader;
	
	// FIXME this is format dependent
	reader = new GEOEASCheckshotReader();

	// set the header words we need
	reader.setCheckshotDataKeyword(KeyWords.Zm, info.chkshotAndMarkerDefinitions.Zm_checkshot);
	reader.setCheckshotDataKeyword(KeyWords.T, info.chkshotAndMarkerDefinitions.T_checkshot);
	reader.setCheckshotDataKeyword(KeyWords.sigmaT, info.chkshotAndMarkerDefinitions.sigmaT_checkshot);

	// get the data
	reader.extractData(info.checkshotFilename, chkShotandMarkerData, EH);

  }

  private void readMarkers(GeneralWellInfo info,CheckshotAndMarkerData chkShotandMarkerData)
  {
	MarkerReader reader;
	
	if (info.markersFilename != null)
	  {
		// FIXME this is format dependent
		reader = new GEOEASMarkerReader();

		// set the header words we need
		reader.setMarkerDataKeyword(KeyWords.Zm, info.chkshotAndMarkerDefinitions.Zm_marker);
		reader.setMarkerDataKeyword(KeyWords.T,   info.chkshotAndMarkerDefinitions.T_marker);
		reader.setMarkerDataKeyword(KeyWords.sigmaT, info.chkshotAndMarkerDefinitions.sigmaT_marker);

		// get the data
		reader.extractData(info.markersFilename, chkShotandMarkerData, EH);
	  }
  }

  // we need a bunch of splines at each well
  private void setUpSaturationAdjustmentModel()
  {
	int i,j,k,n,i_surf,NsatSegs;
	double r,tmin=0,rmin,z,zA;
	double xMid,yMid;
	double[] X,Y,T;
	double[][] Zbr;
	String filename=null;

	if (!useSaturationModel) return;

	convModel.BVR = new ConvolutionalModel.BaseVintageReflectivities();

	// set up parameter structures etc (now in xml parser)
	//for (i=0; i<m.nWells; i++) wellInfo[i].setupSaturationInfo();

	// keep depth window

	ASCIIfileReader rdr = new ASCIIfileReader(EH,true);
	DoubleArrayList Xl = new DoubleArrayList();
	DoubleArrayList Yl = new DoubleArrayList();
	DoubleArrayList Tl = new DoubleArrayList();
	rdr.addToContainerList(null,Xl);
	rdr.addToContainerList(null,Yl);
	rdr.addToContainerList(null,Tl);

	// Sweep 1. Find saturation-shift vertical depth bracket, for each well
	TimeLapseInfo _TLI = m.timelapseinfo;
	NsatSegs=0;
	Zbr = new double[m.nWells][2];
	for (i=0; i<m.nWells; i++)
	  {
		// time to TVD mapping 
		TimeToDepthMap t2dMap  = m.timeToDepthMaps[i];
		LinearInterpolation zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,EH);
		// x,y's of midpoint of survey file
		UTM survey = logData[i].survey;
		xMid = (survey.X[0]+survey.X[survey.X.length-1])/2.0;
		yMid = (survey.Y[0]+survey.Y[survey.Y.length-1])/2.0;

		// read surface
		for(i_surf=0;i_surf<=1;i_surf++)
		  {
			// upper then lower surface of shifted saturation
			// "../horizons/OWCs/OWC_2002_02-fullstk_PPTIME.dat":  <- 02
			// filename = (i_surf==0) ? "../horizons/OWCs/OWC_2008_08-02_IPdiff_mergedw_08-fullstk_PPTIME.dat" : 
			//  "../horizons/OWCs/OOWC_1991_91-fullstk_PPTIME.dat";
			filename = (i_surf==0) ? _TLI.monitor_contacts_filename : _TLI.base_contacts_filename;
			EH.println("Getting contact surface estimation from "+filename);
			// These surfaces are 3column files: Xutm Yutm TWT
			Xl.clear(); Yl.clear(); Tl.clear();
			rdr.extractData(filename);
			X = Xl.elements();
			Y = Yl.elements();
			T = Tl.elements();
			n = X.length;
			for (j=0,rmin=Double.MAX_VALUE; j<n; j++)
			  {
				r= sqr(X[j]-xMid)+ sqr(Y[j]-yMid);
				if (r<rmin)
				  {
					rmin = r;
					tmin=T[j];
				  }
			  }			
			  double RADIUS_CUTOFF = 20.0;
			  // suite of pts within r-nearest+20m...
			  double R=sqr(sqrt(rmin) + RADIUS_CUTOFF);
			  boolean f;
			  zA = (i_surf==0)? Double.MAX_VALUE : -Double.MAX_VALUE;
			  for (j=0,f=false; j<n; j++)
				{
				  r= sqr(X[j]-xMid)+ sqr(Y[j]-yMid);
				  if (r<R) 
					{
					  z = zInterp.interpolate(T[j]/1000,false);  //tmin is in ms FIXME
					  zA = (i_surf==0)? min(z,zA) : max(z,zA);					  
					  //EH.printf("UU close-points TVD=%4f%n",z);				  
					  f=true;
					}
				}
			  // EH.printf("UU %n");
			  if(!f)
				{
				  EH.die_err("No surface info within "+ RADIUS_CUTOFF + " of well "+i); 
				}			
			  boolean inSecs = (i_surf==0) ? _TLI.monitor_timeunitssecs : _TLI.base_timeunitssecs;
			  if (!inSecs) tmin/=1000; //tmin is in ms 
			  // TVD corresponding to surface
			  Zbr[i][i_surf] = zA; // zInterp.interpolate(tmin,false);
			  // Under e.g. primary production expect Zbr[i][1] > Zbr[i][0]
			  EH.printf("Nearest surface point rmin=%4g T=%6f%n",sqrt(rmin),tmin);
		  }
		EH.printf("Well %d Saturation-change depth bracket: Z: [%6f, %6f]%n",i,Zbr[i][0],Zbr[i][1]);
	  }

	// 2. Set up segmented log properties to be morphable under a saturation parameter
	// over these depth brackets
	double[] segDepths=null;    // layer tops in true depth (from segmented data)
	double[] segGamma=null;    // gamma/NG
	double Zl,Zu,wt,NG,EPS=1e-4;
	boolean haveGamma = false;
	NsatSegs=0;
	for (i=0; i<m.nWells; i++)
	  {
		SegmentedLogData _sLogData = sLogData[i];
		_sLogData.SMI = new SegmentedLogData.SaturationModelInfo();
		SegmentedLogData.SaturationModelInfo _SMI = _sLogData.SMI;
		_sLogData.autoDetectRhounits(); // ensure this is correct
		segDepths = _sLogData.Z;
		segGamma = _sLogData.gamma;
	
	
		// store indicator of "swept" zone via facies_indicator
		_sLogData.facies_indicator = new byte[_sLogData.Zm.length]; // same as n of SegmentedLogData
		Arrays.fill(_sLogData.facies_indicator,(byte)1);

		Zu = min(Zbr[i][0],Zbr[i][1]);  // if primary expect Zbr[i][1] > Zbr[i][0]
		Zl = max(Zbr[i][0],Zbr[i][1]);

		if (segGamma != null) haveGamma=true;
		if(!haveGamma) EH.printf("    No VSH log supplied! - assuming net-gross=1 in saturation-change depth bracket%n");
		if(!haveGamma) EH.printf("    Check units of log data and Jacobian sensitivities wrt S below are correct!%n");
	
		// wellInfo[i].saturationInfo.satProp values set in xml parser etc
		// Ctual distribution & mapper properties overridden in ModelVariables (grep sig_Sa)
		
		double[] vp_vs_rho = new double[3];
		double[] sat_slopes= new double[3];
		for(j=0;j<segDepths.length-1;j++)
		  {
			if (segDepths[j+1] > Zu && segDepths[j] <= Zl)
			  {
				// segment partially intersects "sweep zone" [Zu,Zl]
				// 1: "weight" factor. Fraction of segment entirely "in zone"
				wt = (min(segDepths[j+1],Zl) - max(segDepths[j],Zu))/(segDepths[j+1] - segDepths[j]);
				NsatSegs++;
				_sLogData.facies_indicator[j] = 0; //"mark" it

				// Regression approximations to saturation effect
				// Assume we have available a Vsh log, plonked instead of gamma. If not, assume all sand
				NG = haveGamma ? (1-segGamma[j]) : 1.0;
				
				// Typical "stock" numbers from user: ensure unit consistency!				
				double vp_s = _TLI.vp_s;
				double vs_s =_TLI.vs_s;
				double phi_s =_TLI.phi_s;
				double vp_m = _TLI.vp_m;
				double vs_m = _TLI.vs_m;
				double rho_m = _TLI.rho_m;
				double vp_h = _TLI.vp_h;
				double rho_h = _TLI.rho_h;
				double vp_b = _TLI.vp_b;
				double rho_b = _TLI.rho_b;
				double rho_g =_TLI.rho_g;
				double K_g = _TLI.K_g;

				double S_h = wellInfo[i].saturationInfo.Sa; // estimate Jacobian at monitor level mean sat
				
				// compute slope of property responses wrt saturation Gassman, by finite diffs, at supplied
				// **** SI units point *****
				GeophysicsUtils.GassmannBackus(vp_s,vs_s,phi_s, vp_m,vs_m,rho_m, NG, S_h,vp_h, rho_h, 
											   vp_b,rho_b, K_g,rho_g,sat_slopes);
				GeophysicsUtils.GassmannBackus(vp_s,vs_s,phi_s, vp_m,vs_m,rho_m, NG, S_h+EPS,vp_h, rho_h, 
											   vp_b,rho_b, K_g,rho_g,vp_vs_rho);

				
				// Compute finite diff derivative, also discount by segment weight, in logData units
				for (k=0;k<3;k++) sat_slopes[k] = wt * (vp_vs_rho[k] - sat_slopes[k])/EPS; 
				// compensate if logdensity in gm/cc 
				sat_slopes[2] /= _sLogData.densityUnitsConversionFactor;

				// save lookup to "saturation slope", and this segment index, for this log
				_SMI.sensLookups.add(new SegmentedLogData.SegSensLookup(j,sat_slopes,wt,NG));
				EH.printf("Well %d, Gassmann Jacobian: NG= %4g wt= %4g {Vp,Vs,rho}_sat_slope={%6f %6f %6f}%n",(i+1),
						  NG,wt,sat_slopes[0],sat_slopes[1],sat_slopes[2]);
			  }
		  }
	  }
	if (NsatSegs==0)
	  {
		EH.println("No well zones in allowed depth widow with variable saturation for supplied wells.");		
		EH.println("Wavelet extraction will be nonsensical");
		EH.die_err("");
	  }
  }
}
