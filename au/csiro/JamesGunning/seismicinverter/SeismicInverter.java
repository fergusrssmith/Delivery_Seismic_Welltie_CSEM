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

import java.lang.*;
import java.text.*;
import java.io.*;
import java.util.*;

import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;



/**
   Basic class for performing trace by trace stochastic inversion
   based on principles in the paper:

   James Gunning and Michael Glinsky,
   Delivery: an open-source model-based Bayesian seismic inversion program,
   Computers and Geosciences,2004, v30, pp619-636.
   Updated pdf versions available in the distribution.

   @author James.Gunning@csiro.au
   @version 1.0, 9/10/2001 
*/

public final class SeismicInverter
{
  String classname=null;
  Model m;
  ErrorHandler EH;
  ModelFormat MF;
  FundamentalScalingConstants FSC;


  // runtime parameters
  String commandline;
  String rockandfluids_xml_filename;
  String modelDescription_xml_filename;
  String su_zero_offset_filename;
  String su_far_offset_filename;
  String su_model_filename=null;


  /**
	 this is the string for the filename of sample model outputs in su form,
	 which take the form realisations.su etc. Traces for a given realisation
	 are tagged by the header word realisationNumberHeaderField.
  */
  public static final String SU_REALISATIONS_FILENAME="realisations.su";
  String su_realisations_filename=SU_REALISATIONS_FILENAME;

  /**
	 this is the string for the filename of model MAP outputs in su form,
	 which take the form MAP_realisation.su etc. 
  */
  public static final String SU_MAP_REALISATION_FILENAME="realisation.MAP.su";
  String su_MAP_realisation_filename = SU_MAP_REALISATION_FILENAME;


  /**
	 debugOptions: turn on for debugging options to appear
  */
  static boolean debugOptions = true;

  int verbosity=0;
  boolean drawPriorPlusIsopachModel = false;
  boolean overwritePriorTimesFromDepths = false;
  boolean drawMAPmodels = false;
  boolean drawMinChisqModels = false;
  boolean drawPriorPlusIsopachModelForAllFluidCombos = false;
  boolean doInversion = false;
  boolean modelTracesAreProvided = false;
  boolean timesOnlyInOptimiser = false;
  boolean ignoreIsopachs = false;
  boolean ignoreDepthConsistency = false;
  boolean ignoreVerticalStackConstraint = false;
  boolean lockDepths = false;
  boolean outputParametricDepths = false;
  boolean ignoreSeismic = false;
  boolean writeHiddenSeismicError = false;
  boolean useFluidProbabilitiesAsPostOrdering = false;
  boolean ignoreAllButFirstFluidState = false;
// fluid state to confine interest to if ignoreAllButFirstFluidState=true
  int restrictedFluidStateIndex=0;
  boolean dumpSyntheticSeismic = false;
  boolean truncateSynthetics = false;
  boolean suppressModelDumping = false;
  int     density_ordering_mode = FluidStates.fullDensityOrdering;
  boolean outputDepthModels = false;
  boolean makeSyntheticAndWavelet = false;
  boolean usePinchoutDetection;
  boolean useMAPcheck = true;
  boolean useGeneticAlgorithmsIfNeeded=true;
  boolean useDescentTest=true;
  boolean useEventDetection = true;
  boolean UseRandomWalkSampler = false;
  boolean useApproximateMixedGaussianIndependenceSampler=false;
  boolean reParameteriseToAdditiveTimes = false;
  boolean useSmartModeFindingStartingPoints = true;
  boolean decimateRandomWalkSamplerOutput = true;
  boolean performGelmanConvergenceDiagnostics = false;
  boolean performRafteryConvergenceDiagnostics = false;
  boolean useSeed = false;
  int seed = cern.jet.random.engine.MersenneTwister.DEFAULT_SEED;
  boolean plotTrend = false;
  String  plotTrendRock,plotTrendType;
  boolean useJNI = true;
  double waveletTrimEpsilon = 0.0; // wavelet coda/precursor trimming epsilon

  /**
	 Set this true to always do a Raftery style diagnostic on runs, whether we choose to
	 do Gelman, auto--decimation or whatever else. Overhead is very small.
  */
  public static final boolean alwaysPerformRafteryConvergenceDiagnostics = true;
  boolean performRafteryAutoDecimation = false;
  TraceFilter traceFilter = null;
  boolean doDumpDebug = false;
  boolean dumpOptimalSeismic = false;
  boolean dumpRealisations = true;
  String  cmdLineNearWaveletFilename = null;
  String  cmdLineFarWaveletFilename = null;
  String  cmdLineNearStackFilename = null;
  String  cmdLineFarStackFilename = null;
  boolean dumpGassmanDiagnostics = false;
  boolean checkPriorSanity = false;
  boolean printXsections = false;
  double  XsectionsDeviation = 2.0;
  boolean printModeProfiles = false;
  private boolean doBHPcommand=false;
  boolean doStackedModel = false;
  boolean doSimpleFluidUnitsAndContacts = false;
  // ensure fluid combos are balanced ?? FIXME: compensate in sampler??????
  boolean equaliseModeCountAcrossFluidCombinations = !true;


  // if we switch this on, the -N flag yields number of sample before decimation, i.i
  // not taking decimation into account 
  boolean haveSpecifiedPreDecimationNsamples = false;
  int n_samples_per_trace_undecimated=10;

  /** 
	  Switch to enable depth self consistency. Disabled for now
  */
  public static final boolean ENABLE_DEPTH_SELF_CONSISTENCY = false;

  /**
	 A switch to detect use of ms in time samples, rather than seconds
  */
  final boolean MODEL_TRACE_TIME_UNITS_AUTO_DETECTION = true;

  /**
	 Time values in priors or model traces above this value signify
	 use of milliseconds in model traces, if
	 MODEL_TRACE_TIME_UNITS_AUTO_DETECTION is switched on.  Rationale:
	 no realistic times are going to be longer than 50secs, but will
	 certainly be longer than 50mS
  */
  static final double TIME_UNITS_AUTO_DETECTION_TIME_MAX = 50;
 
  /**
   Signing "Unknown" or "Bad Value" entries in model traces is performed
	 by setting them to this number:
  */
  public static final double MODEL_TRACE_BAD_NUMBER = -999.0;

  /**
	 flag in realisations file to signify failed inversion (optimiser defeated)
  */
  public static final double NO_MODES_FOUND_CHISQ_BIGMARKER = MODEL_TRACE_BAD_NUMBER;

  /**
	 Minimum acceptable depth uncertainty used as threshold on sigma_dz
  */
  public static final double MINIMUM_ACCEPTABLE_ISOPACH_UNCERTAINTY = 0.05;

  /**
	 minimum acceptable noise-signal ratio for inversions
  */
  public static final double MIN_ACCEPTABLE_NS_RATIO=1e-3;

  /**
	 This flag enforces all internal representations of depth to be positive
  */
  final boolean MAKE_MODEL_TRACE_PRIOR_DEPTHS_POSITIVE = true;

  /**
	 A switch to detect layer synchronisation
  */
  final boolean DETECT_LAYER_SYNCHING = true;

  


  /**
	 Default synthetic wavelet bandwidth
  */
  float syntheticFreq = 40.0f;
  /**
	 Default synthetic sampling period (s). 4 ms data
  */
  double syntheticSampleTime = 0.004;
  /**
	 Start time for synthetics if making of prior synthetic required
  */
  double tSyntheticStart = 0;

  /**
	 End time for synthetics if making of prior synthetic required
  */
  double tSyntheticEnd = 0;


  // Member variables
  

  SuReader sureader;
  /**
	 seismic amplitude near-stack traces
  */
  SuSegyFormat tr_near;

  /**
	 seismic amplitude far-stack traces
  */
  SuSegyFormat tr_far;

  /**
	 synthetic seismic amplitude near-stack traces
  */
  SuSegyFormat tr_near_synth;

  /*
	 synthetic amplitude far-stack traces
  */
  SuSegyFormat tr_far_synth;

  /**
	 model traces
  */
  SuSegyFormat tr_model;
  /** 
	  input stream for model traces
  */
  FileInputStream model_in_strm;

  /** 
	  input stream for near-stack data
  */
  FileInputStream near_stack_su_in_strm=null;

  /** 
	  input stream for far-stack data
  */
   FileInputStream far_stack_su_in_strm;
  

  /**
	Writer for various su output files
  */
  SuWriter suwriter;
  

  /** 
	  output stream for near-stack synthetic traces
	  if desired
  */ 
  FileOutputStream near_synth_out_strm=null;
  /** 
	  output stream for far-stack synthetic traces
	  if desired
  */
  FileOutputStream far_synth_out_strm=null;

  /** 
	output stream for realisation model traces
  */

  FileOutputStream realisations_out_strm;

  /** 
	output stream for realisation MAP traces
  */

  FileOutputStream realisations_MAP_out_strm;

  /**
	 prefix for synthetic data sets to attach to near/far stack filenames
	 (to avoid overwriting data)
  */
  public static final String SYNTH_NEAR_FILENAME="near_stack_synth.su";
  public static final String SYNTH_FAR_FILENAME="far_stack_synth.su";

  // Debug package dumper
  DebugPackageDumper dumper=null;

  // flag for catching when zero times flag trace-to-be-skipped
  private boolean haveZeroTimesFlaggedSkipTrace=false;

  // BHP commands to be forked as processes, formed partially from this string
  private String BHPcommandstring="";
  private String BHPcommand_filename;

  int n_traces=0;
  
  SeismicInverter()
  {
	classname = ((new Throwable()).getStackTrace())[0].getClassName();

	// get fundamental scaling constants. Will/may need this early on to check against
	// commandline args for sanity etc
	FSC = new FundamentalScalingConstants();
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();
	// Always need a trace filter object
	traceFilter = new TraceFilter(EH,0); // update verbosity later
  }


  /** 
	  Sample driver from input streams etc
  */
  public static void main(String argv[])
  {
	SeismicInverter SI = new SeismicInverter();


	

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("SeismicInverter",SI.EH);

	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber); 


   	options.addMandatoryArg(1,"model_file.xml","xml file containing the full Model description"+
							", plus rock physics and fluid properties, and details of AVO seismic model");
	
	options.addOption("-m",null,1,"String","model_traces.su",
					  "su file containing model traces, if model priors vary");	
	options.addOption("-o",null,1,"String",SU_REALISATIONS_FILENAME,
					  "su file for output realisation traces (default "+SU_REALISATIONS_FILENAME+")");
	options.addOption("-oMAP",null,1,"String",SU_MAP_REALISATION_FILENAME,
					  "su file for outputting MAP realisation traces (default "+SU_MAP_REALISATION_FILENAME+")",1);
	options.addOption("-v",null,1,"int","level","Execute with stated level of verbosity");
	options.addOption("-p",null,0,null,null,"Generate maximum likelihood model from prior+isopach(+depth_consistency) only, for most-likely prior fluid combination");
	options.addOption("--p-depths",null,0,null,null,"Set prior times by prior depths and prior-mean layer velocities, relative to master layer time, using most-likely prior fluid combination",1);
	options.addOption("--p-fluid",null,0,null,null,"Generate models from prior +isopach only, for all fluid combinations",1);
	options.addOption("--MAP",null,0,null,null,"Dump MAP model of most likely mode in realisation file (-oMAP)");
	options.addOption("--ML",null,0,null,null,"Dump minimum chisq model in realisation file (-oMAP). Might not be most--likely fluid state.");

	options.addOption("-N",null,1,"int","number of samples","Attempt this many samples per trace in MCMC method");
	options.addOption("--Npre-decimation",null,0,null,null,"signifies that -N flag refers to number of samples before decimation",1);
	options.addOption(null,"--make-synthetic",4,"float","freq(Hz) sample-interval(ms) t_start(sec) t_end(sec)",
					  "Generate synthetic model with stated seismic wavelet bandwidth and sampling time, "+
					  "over interval specified. Wavelet filenames (in XML or commandline) overwritten with "+
					  "new Ricker wavelet. Seismic written to "+SYNTH_NEAR_FILENAME+", "+SYNTH_FAR_FILENAME);
	options.addOption("-II",null,0,null,null,"Ignore isopach specifications");
	options.addOption("-IVS",null,0,null,null,"Ignore vertical-stacked-models constraints",1);
	options.addOption("-ID",null,0,null,null,"Ignore (skip) depth self-consistency "+ 
						  ((!SeismicInverter.ENABLE_DEPTH_SELF_CONSISTENCY)? "(N.B. depth self--consistency disabled since 13/4/04))":""),1);
	options.addOption("--lock-depths",null,0,null,null,"lock all depths in prior to that of master depth layer in XML model file",1);
	options.addOption("--parametric-depths",null,0,null,null,"Output parametric depths, not time*velocity hung depths",1);
	options.addOption("-IS",null,0,null,null,"Ignore seismic (sets S/N ratio to very low number)");
	options.addOption("--write-hidden-seismic-error",null,0,null,null,"Write \"hidden\" seismic error to output traces, even if -IS is specified",1);
	options.addOption("-IF",null,0,null,null,"Ignore all but first fluid state (testing tool for speed)",1);
	options.addOption("--IF",null,1,"int","fluid-combo","Examine only the specified fluid state (useful testing tool)");
	options.addOption("-S",null,2,"float","t_start(sec) t_end(sec)","Dump synthetic seismic as inversion proceeds",1);
	options.addOption("--truncate-synthetics",null,0,null,null, 
					  "Sets synthetic seismic to zero outside error-acquisition range (flavour of -S option)",1);
	options.addOption("--suppress-model-dump",null,0,null,null,"Suppress dumping of model su files",1);
	options.addOption("-RWS",null,0,null,null,"Use Random Walk (small step Metropolis) sampler");
	options.addOption("-AIS","--approx-indep-sampler",0,null,null,"Use approximate independence sampler",1);
	options.addOption("--Gelman-diagnostics",null,0,null,null,"Compute Gelman-style performance diagnostics on MCMC method (expensive)",1);
	options.addOption("--Raftery-diagnostics",null,0,null,null,"Compute Raftery-style performance diagnostics on MCMC method",1);
	options.addOption("--Raftery-auto-decimation",null,0,null,null,"Use Raftery diagnostics to perform auto-decimation (expensive)",1);
	options.addOption("--no-decimation",null,0,null,null,"Stop decimation of Hopping sampler output if used",1);
	options.addOption("-t",null,1,"float","epsilon","Trim wavelet coda/precursors below this relative threshold",1);
	options.addOption("--trace-filter",null,1,"String","Expression",
					  "invert only on  traces satisfying shell-protected, comma-separated header word filter "+
					  "expressions involving comparisons <,>,=,!=,>=,<= e.g. \"fldr>100,ep=234\"");

	
	options.addOption("--dump-debug-package",null,0,null,null,
					  "dump debug zip package on *first* trace satisfying --trace-filter, if any");	
	options.addOption("--times-only",null,0,null,null,"Use only layer-times as parameters when mode hunting",1);
	options.addOption("--additive-times",null,0,null,null,"Use additive layer-times reparameterisation when mode hunting",1);
	options.addOption("--no-genetics",null,0,null,null,"Never use genetic methods when mode hunting (can speed things up if prior model is close)",1);
	options.addOption("--no-descent-test",null,0,null,null,"Omit aggressive descent test (makes global optimization slower & more thorough)",1);
	options.addOption("--no-event-detection",null,0,null,null,"Skip event detection algorithms and consequent mode searches",1);
	options.addOption("--disable-MAP-check",null,0,null,null,"Don't do fancy stuff if check of mode points with time collocation fails",1);
	options.addOption("--no-smart-mode-finding-starting-points",null,0,null,null,
					  "Don't use smart starting points when mode hunting (slower, but good for tests)",1);
	options.addOption("--dump-optimal-seismic",null,0,null,null,"Dump optimal seismic traces "
					  +"when mode hunting (sampling suppressed)",1);
	options.addOption("--dont-dump-realisations",null,0,null,null,"Don't dump realisation traces",1);
	options.addOption("--seed",null,1,"int","seed value","Seed all internal random generators from this value. Use -ve value (quote-protect) "+
                                            " to seed from \"System.nanoTime()\", e.g. --seed '-1'",1);
	options.addOption("--Gassman-diagnostics",null,0,null,null,"Dump Gassman diagnostics if Gassman substitution is giving trouble",1);
	options.addOption("-PD",null,0,null,null,"Use smart pinchout detection in mode finding",1);
	options.addOption("--check-prior",null,0,null,null,"Check prior model: no sampling or mode searching. "+
					  "Simple sanity checks, especially on fluid substitution",1);
	/*
	  options.addOption("--BHPcommand",null,2,"String","filename shell-protected-command",
					  "write to filename the script \"shell-protected-command properties=x,x..\" and make it executable");
	*/
	options.addOption("--BHPcommand",null,1,"String","filename",
					  "write to filename the string \"properties=x,x..\" ",1);
	options.addOption("--near-stack-wavelet-filename",null,1,"String","near-wavelet.su","near-stack wavelet filename",1);
	options.addOption("--far-stack-wavelet-filename",null,1,"String","far-wavelet.su","far-stack wavelet filename",1);
	options.addOption("--near-stack-filename",null,1,"String","near-stack.su","near-stack seismic data file or named-pipe",1);
	options.addOption("--far-stack-filename",null,1,"String","far-stack.su","far-stack seismic data file or named-pipe",1);
	options.addOption(null,"--no-jni",0,null,null,"Disable native libraries (use java fallbacks)",1);

	// debugging options
	if (SI.debugOptions)
	  {
		options.addOption("-D1",null,0,null,null,"Turn on debug switch 1",1);
		options.addOption("--Dval",null,1,"float","value","use this value in debugging",1);
		options.addOption("--Xsections",null,1,"float","n(standard-deviations)","Print Xsections at modes to posteriorXsection* files (very extensive output!). Suggest n(standard-deviations)=2");
		options.addOption("--mode-profiles",null,0,null,null,"Print mode-joining profiles when mode checking (posterior prob along line joing ML points for same fluid-state modes)",1);
		options.addOption("--plot-trend",null,2,"String","type rock-name",
						  " plot to stdout a {x,y} set for named trend (vp|vs|rho|phi|imp) over model range");
	  }

	// leave this to last!
	options.addExtraHelpOption("-X","Show help with all options (there are lots more!)");

	options.parseOptions(argv,true);


	// stuff commandline options into Seismicinverter fields
	SI.modelDescription_xml_filename = options.getMandatoryArg(1);
	// now amalgated...
	SI.rockandfluids_xml_filename = SI.modelDescription_xml_filename;

	if (options.isOptionSet("-X")) options.Usage();

	if (options.isOptionSet("-v")) 
	  SI.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();

	if (options.isOptionSet("--seed"))
	  {
		SI.useSeed = true;
		SI.seed = ((Integer)options.getOptionValues("--seed")[0]).intValue();
	  }

	if (options.isOptionSet("-m")) 
	  {
		SI.su_model_filename = ((String)options.getOptionValues("-m")[0]);
		SI.modelTracesAreProvided = true;
	  }

	if (options.isOptionSet("-o")) 
	  SI.su_realisations_filename = ((String)options.getOptionValues("-o")[0]);

	if (options.isOptionSet("-oMAP")) 
	  SI.su_MAP_realisation_filename = ((String)options.getOptionValues("-oMAP")[0]);

	if (options.isOptionSet("--p-fluid")) 
	  SI.drawPriorPlusIsopachModelForAllFluidCombos = true;

	if (options.isOptionSet("-p")) 
	  SI.drawPriorPlusIsopachModel = true;

	if (options.isOptionSet("--p-depths")) 
	  SI.overwritePriorTimesFromDepths = true;

	if (options.isOptionSet("--MAP")) 
	  SI.drawMAPmodels = true;

	if (options.isOptionSet("--ML")) 
	  {
		SI.drawMinChisqModels = true;
		SI.drawMAPmodels = false;
	  }

	if (options.isOptionSet("-II")) 
	  SI.ignoreIsopachs = true;

	if (options.isOptionSet("-IS")) 
	  SI.ignoreSeismic = true;

	if (options.isOptionSet("-IVS")) 
	  SI.ignoreVerticalStackConstraint = true;

	if (options.isOptionSet("--write-hidden-seismic-error")) 
	  SI.writeHiddenSeismicError = true;

	if (options.isOptionSet("-ID")) SI.ignoreDepthConsistency = true;
	// override 
	if (!SeismicInverter.ENABLE_DEPTH_SELF_CONSISTENCY) SI.ignoreDepthConsistency = true;

	if (options.isOptionSet("--lock-depths")) 
	  SI.lockDepths = true;

	if (options.isOptionSet("--parametric-depths")) 
	  SI.outputParametricDepths = true;

	if (options.isOptionSet("-IF")) 
	  SI.ignoreAllButFirstFluidState = true;

	if (options.isOptionSet("--IF"))
	  {
		SI.ignoreAllButFirstFluidState = true;
		SI.restrictedFluidStateIndex = ((Integer)options.getOptionValues("--IF")[0]).intValue();
	  }

	if (options.isOptionSet("-S"))
	  {
		SI.dumpSyntheticSeismic = true;
		SI.tSyntheticStart = (double)((Float)(options.getOptionValues("-S")[0])).floatValue();
		SI.tSyntheticEnd = (double)((Float)(options.getOptionValues("-S")[1])).floatValue();
		SI.cleanSyntheticRange();
		if (SI.tSyntheticEnd <= SI.tSyntheticStart)
		  {
			SI.EH.warn("Synthetic end time must be greater than start time");
			options.Usage();
		  }		
	  }
	if (options.isOptionSet("--truncate-synthetics"))
	  {
		SI.truncateSynthetics = true;
	  }

	if (options.isOptionSet("--suppress-model-dump")) 
	  SI.suppressModelDumping = true;

	if (options.isOptionSet("-PD")) 
	  SI.usePinchoutDetection = true;

	if (options.isOptionSet("--no-genetics")) 
	  SI.useGeneticAlgorithmsIfNeeded = false;

	if (options.isOptionSet("--no-descent-test")) 
	  SI.useDescentTest = false;

	if (options.isOptionSet("--no-event-detection")) 
	  SI.useEventDetection = false;

	if (options.isOptionSet("--disable-MAP-check")) 
	  SI.useMAPcheck = false;

	if (options.isOptionSet("-RWS")) 
	  SI.UseRandomWalkSampler = true;

	if (options.isOptionSet("--approx-indep-sampler") || options.isOptionSet("-AIS"))
	  {
		SI.UseRandomWalkSampler = false; // this had better be off
		SI.useApproximateMixedGaussianIndependenceSampler = true;
	  }

	if (options.isOptionSet("--no-decimation"))
	  SI.decimateRandomWalkSamplerOutput = false;

	if (options.isOptionSet("-N")) 
	  {
		SI.n_samples_per_trace_undecimated = ((Integer)options.getOptionValues("-N")[0]).intValue();
		if (SI.n_samples_per_trace_undecimated > Short.MAX_VALUE) 
		  SI.EH.warn("N must be <= "+Short.MAX_VALUE); // reason: realisation # is stored in a short header field
	  }

	if (options.isOptionSet("--Npre-decimation"))
	  SI.haveSpecifiedPreDecimationNsamples = true;

	if (options.isOptionSet("--Gelman-diagnostics")) 
	  SI.performGelmanConvergenceDiagnostics = true;

	if (options.isOptionSet("--Raftery-diagnostics") || SI.alwaysPerformRafteryConvergenceDiagnostics) 
	  SI.performRafteryConvergenceDiagnostics = true;

	if (options.isOptionSet("--Raftery-auto-decimation"))
	  {
		SI.performRafteryConvergenceDiagnostics = true;
		SI.performGelmanConvergenceDiagnostics = false; // have to turn this off
		SI.decimateRandomWalkSamplerOutput = true;      // sanest choice here, to get correlations in ballpark
		SI.performRafteryAutoDecimation =true;
	  }
	if (options.isOptionSet("-t"))
	  {
		SI.waveletTrimEpsilon = (double)((Float)(options.getOptionValues("-t")[0])).floatValue();
		if (SI.waveletTrimEpsilon<0 || SI.waveletTrimEpsilon>=1)
		  {
			SI.EH.warn("Wavelet trimming epsilon must be in [0,1)");
			options.Usage();
		  }		
	  }

	if (options.isOptionSet("--times-only")) 
	  SI.timesOnlyInOptimiser = true;

	if (options.isOptionSet("--additive-times")) 
	  SI.reParameteriseToAdditiveTimes = true;

	if (options.isOptionSet("--no-smart-mode-finding-starting-points")) 
	  SI.useSmartModeFindingStartingPoints = false;

	if (options.isOptionSet("--trace-filter")) 
	  {
		// need also a trace filter, even if its inactive. Ensure verbosity correct here
		SI.traceFilter.verbosity = SI.verbosity;
		SI.traceFilter.applyTraceFilter = true;
		SI.traceFilter.cmdLineTraceFilterString = (String)options.getOptionValues("--trace-filter")[0];
	  }
	if (options.isOptionSet("--dump-debug-package"))
	  {
		SI.doDumpDebug = true;
	  }

	if (options.isOptionSet("--dump-optimal-seismic")) 
	  SI.dumpOptimalSeismic = true;

	if (options.isOptionSet("--dont-dump-realisations"))
	  SI.dumpRealisations = false;

	if (options.isOptionSet("--near-stack-filename"))
	  {
		SI.cmdLineNearStackFilename = (String)options.getOptionValues("--near-stack-filename")[0];
	  }
	if (options.isOptionSet("--far-stack-filename")) 
	  {
		SI.cmdLineFarStackFilename = (String)options.getOptionValues("--far-stack-filename")[0];
	  }
	if (options.isOptionSet("--no-jni")) SI.useJNI = false;

	if (options.isOptionSet("--near-stack-wavelet-filename")) 
	  {
		SI.cmdLineNearWaveletFilename = (String)options.getOptionValues("--near-stack-wavelet-filename")[0];
	  }

	if (options.isOptionSet("--far-stack-wavelet-filename")) 
	  {
		SI.cmdLineFarWaveletFilename = (String)options.getOptionValues("--far-stack-wavelet-filename")[0];
	  }

	if (options.isOptionSet("--Gassman-diagnostics")) 
	  SI.dumpGassmanDiagnostics = true;


	if (options.isOptionSet("--make-synthetic")) 
	  {
		SI.makeSyntheticAndWavelet = true; 
		SI.drawPriorPlusIsopachModel = true;
		// kill all other options...
		SI.drawPriorPlusIsopachModelForAllFluidCombos = false;
		// set synthetic freq
		SI.syntheticFreq = ((Float)(options.getOptionValues("--make-synthetic")[0])).floatValue();
		// and dt. Convert to seconds
		SI.syntheticSampleTime = 1.0e-3*((Float)(options.getOptionValues("--make-synthetic")[1])).floatValue();
		SI.tSyntheticStart = (double)((Float)(options.getOptionValues("--make-synthetic")[2])).floatValue();
		SI.tSyntheticEnd = (double)((Float)(options.getOptionValues("--make-synthetic")[3])).floatValue();
		SI.cleanSyntheticRange();
		if (SI.tSyntheticEnd <= SI.tSyntheticStart)
		  {
			SI.EH.warn("Synthetic end time must be greater than start time%n");
			options.Usage();
		  }
	  }

	if (options.isOptionSet("--check-prior"))
	  {
		SI.checkPriorSanity = true;
		// kill other options
		SI.makeSyntheticAndWavelet = false; 
		SI.drawPriorPlusIsopachModelForAllFluidCombos = false;
	  }

	// BHP commands
	if (options.isOptionSet("--BHPcommand"))
	  {
		SI.doBHPcommand = true;
		SI.BHPcommand_filename = (String)options.getOptionValues("--BHPcommand")[0];
		// SI.BHPcommandstring = (String)options.getOptionValues("--BHPcommand")[1];
	  }

	// debugging options
	if (SI.debugOptions)
	  {
		if (options.isOptionSet("--Xsections"))
		  {
			SI.printXsections= true;
			SI.XsectionsDeviation = (double)(((Float)(options.getOptionValues("--Xsections")[0])).floatValue());
			if (SI.XsectionsDeviation<=0)
			  {
				SI.EH.warn("--Xsections value must have value>0");
				options.Usage();
			  }
		  }
		if (options.isOptionSet("--mode-profiles"))
		  SI.printModeProfiles = true;
		if (options.isOptionSet("--plot-trend"))
		  {
			SI.plotTrend=true;
			SI.plotTrendType = (String)options.getOptionValues("--plot-trend")[0];
			SI.plotTrendRock = (String)options.getOptionValues("--plot-trend")[1];
		  }
	  }

	// Save commandline

	StringBuffer cl = new StringBuffer();
	for (int i=0;i<argv.length;i++) cl.append(argv[i]+" ");
	SI.commandline = cl.toString();

	// Show for posterity:
	SI.EH.println("java "+SI.classname+" "+SI.commandline+"%n");

	if (!SI.setUpSeismicInverter())
	  {
		options.Usage();
	  }
	// Now do it.
	SI.runSeismicInverter();
  }

  boolean setUpSeismicInverter()
  {
	boolean success=false;
	// basic objects needed for model
	RockandFluidData RFD = new RockandFluidData();

	// create model in which to stuff layer and rock/fluid data
	m = new Model(this);

	// Set up Output formatter
	MF = new ModelFormat(true,0);
	// handle BHP command (done below now, after dropping bloat)
	// as early as possible (MF must be set up)
	// if (doBHPcommand) handleBHPcommand();

	// Current implementation uses one-D convolutions
	m.SM = new OneDSeismicConvolutionalModel(true, m);

	// (1) get layer descriptions, seismic model, default properties etc	
	ModelDescriptionParser MDP = new ModelDescriptionParser(modelDescription_xml_filename, m,
															this.verbosity, this.EH);
	success = MDP.readInModelAndSeismicDescriptions(modelDescription_xml_filename,this.m);
	if (!success) return false;
	
	if (plotTrend)
	  {
		m.RFD.plotTrend(this.m,plotTrendType,plotTrendRock);
		System.exit(0);
	  }

	// initialise Seismic model
	if (!(m.SM.initialiseSeismicModel(m.n_layers))) return false;


	success = checkFilesRequired();
	if (!success) return false;



	// complete set up of output formatter
	recomputeLeanFormat();

	// dump out what we use if asked
	if (doBHPcommand) handleBHPcommand();

	MF.initialiseModelFormat(m.n_layers);

	return true;
  }

  // Try to figure out less wasteful output formats, on the XML basis
  void recomputeLeanFormat()
  {
	m.determineIfHardRockModelNeeded();
	m.determineIfFloatingGrainModelNeeded();

	if (!m.needFloatingGrainTerms)
	  {
		MF = MF.getNewModelFormatSansFloatingGrainFields();
	  }
	if (!m.needHardRockTerms)
	  {
		MF = MF.getNewModelFormatSansHardRockFields();
	  }

  }

  boolean runSeismicInverter()
  {
	setUpTraceInputs();
	traceFilter.setUpTraceInputFilter(tr_near);
	setUpRealisationOutput();
	initialiseSyntheticSeismicOutput();

	try
	  {
		if (!loopOverTraces()) return false;
		// loopOverTraces();
	  }
	catch (Exception ex)
	  {
		formDebugPackage(ex);
		ex.printStackTrace();
		return false;
	  }

	finalizeTraceInputs();
	finalizeRealisationOutput();
	finalizeSyntheticSeismicOutput();
	return true;
  }
  
  void formDebugPackage(Exception ex)
  {
	int i=0;

	String[] oFiles=null;
	String[] accessibleSuFiles=null;
	DebugPackageDumper.SuTraceInfo[] traceInfo=null;

	// XML
	oFiles = new String[1];
	oFiles[0] = modelDescription_xml_filename;


	// num stacks
	int nS = (m.SM.getFarstackData() != null) ? 2 : 1;

	// Wavelets
	accessibleSuFiles=new String[nS];
	accessibleSuFiles[0] = m.SM.getNearstackData().wavelet_filename;
	if (nS==2) accessibleSuFiles[1] = m.SM.getFarstackData().wavelet_filename;

	// Inplace traces
	int nTraces = 0;
	if (modelTracesAreProvided) nTraces++;
	// plus 1,2 stacks
	nTraces+=nS;
	traceInfo = new DebugPackageDumper.SuTraceInfo[nTraces];
	if (modelTracesAreProvided)
	  {
		traceInfo[0] = new DebugPackageDumper.SuTraceInfo(su_model_filename,
														  tr_model);
		i++;
	  }
	traceInfo[i] = new DebugPackageDumper.SuTraceInfo(m.SM.getNearstackData().filename,
													  tr_near);
	i++;
	if (nS==2)
	  {
		traceInfo[i] = new DebugPackageDumper.SuTraceInfo(m.SM.getFarstackData().filename,
														  tr_far);
	  }

	if (dumper==null)
	  {
		dumper = new DebugPackageDumper(EH, oFiles,accessibleSuFiles,
										traceInfo,commandline,ex);
	  }

	dumper.dumpPackage();	
  }

  private void setUpTraceInputs()
  {
	// model traces
	if (modelTracesAreProvided) 
	  {
		try 
		  {
			model_in_strm = (FileInputStream) new FileInputStream(su_model_filename);
		  }
		catch (IOException ex)
		  {
			throw new RuntimeException("couldn't open " + su_model_filename);
		  }
	  }

	// seismic traces
	tr_near = new SuSegyFormat();
	if (!makeSyntheticAndWavelet)
	  {
		// must be able to input seismic data
		// nears
		try 
		  {
		    near_stack_su_in_strm = (FileInputStream) new 
			  FileInputStream(m.SM.getNearstackData().filename);			
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + 
									   m.SM.getNearstackData().filename);
		  }
		// fars
		try 
		  {
			// there may be only near-stack data....
			StackData sd = m.SM.getFarstackData();
			
			if (!(sd==null))
			  {
				if (m.SM.getFarstackData().filename!=null)
				  {
					far_stack_su_in_strm = (FileInputStream) new 
					  FileInputStream(m.SM.getFarstackData().filename);
					// only need this if there are fars????
					tr_far = new SuSegyFormat();
				  }
			  }
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + 
									   m.SM.getFarstackData().filename);
		  }
	  }
	else
	  {
		tr_far = new SuSegyFormat();
	  }
  }

  private void finalizeTraceInputs()
  {
	if (!makeSyntheticAndWavelet)
	  {
		// close seismic data streams
		// nears
		try 
		  {
		    near_stack_su_in_strm.close();
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close " + 
									   m.SM.getNearstackData().filename);
		  }
		// fars
		if (m.SM.getFarstackData()!=null)
		  {
			try 
			  {
				far_stack_su_in_strm.close();
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close " + 
										   m.SM.getFarstackData().filename);
			  }
		  }
	  }
	if (modelTracesAreProvided) 
	  {
		try 
		  {
		    model_in_strm.close();
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close " + su_model_filename);
		  } 
	  }	
  }


  private void initialiseSyntheticSeismicOutput()
  {
	// set up synthetic seismic output
	if (dumpSyntheticSeismic || makeSyntheticAndWavelet || dumpOptimalSeismic)
	  {
		// nears
		if (m.SM.nearstack_data != null)
		  {
			tr_near_synth = new SuSegyFormat();		
			// Nears. Open stream.
			try 
			  {
				near_synth_out_strm = (FileOutputStream) new 
				  FileOutputStream(SYNTH_NEAR_FILENAME);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" 
										   + SYNTH_NEAR_FILENAME);
			  }
		  }
		// fars 
		if (m.SM.farstack_data != null)
		  {
			tr_far_synth = new SuSegyFormat();

			// fars. Open stream.
			try 
			  {
				far_synth_out_strm = (FileOutputStream) new 
				  FileOutputStream(SYNTH_FAR_FILENAME);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" 
										   + SYNTH_FAR_FILENAME);
			  }
		  }
	  }	  
  }


  // flush and close synthetic seismic output streams
  void finalizeSyntheticSeismicOutput()
  {
	if (dumpSyntheticSeismic || makeSyntheticAndWavelet || dumpOptimalSeismic)
	  {
		if (near_synth_out_strm!=null)
		  {
			try
			  {
				suwriter.flush(near_synth_out_strm); // must do this!
				near_synth_out_strm.close();
				if (verbosity>2) 
				  EH.println("%nWrote near-stack synthetic data to " + SYNTH_NEAR_FILENAME);
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close " + SYNTH_NEAR_FILENAME);
			  }
		  }
		if (far_synth_out_strm!=null)
		  {
			try
			  {
				suwriter.flush(far_synth_out_strm); // must do this!
				far_synth_out_strm.close();
				if (verbosity>2) 
				  EH.println("%nWrote far-stack synthetic data to " + SYNTH_FAR_FILENAME);
			  }
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't close " + SYNTH_FAR_FILENAME);
			  }
		  }
	  }
  }

  // opens the file for the realisations in su form
  void setUpRealisationOutput()
  {
	suwriter = new SuWriter();

	if (!suppressModelDumping)
	  {
		try 
		  {
			// open stream for model traces
			realisations_out_strm = (FileOutputStream) new 
			  FileOutputStream(su_realisations_filename);
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open" 
									   + su_realisations_filename);
		  }
		if (drawMAPmodels || drawMinChisqModels)
		  {
			try 
			  {
				// open stream for model traces
				realisations_MAP_out_strm = (FileOutputStream) new 
				  FileOutputStream(su_MAP_realisation_filename);
			  } 
			catch (IOException e) 
			  {
				throw new RuntimeException("couldn't open" 
										   + su_MAP_realisation_filename);
			  }	
		  }
	  }
  }

  // Closes the output file for the current realisation su stream
  void finalizeRealisationOutput()
  {
	if (!suppressModelDumping)
	  {
		try 
		  {
			suwriter.flush(realisations_out_strm); // must do this!
			realisations_out_strm.close();
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close realisations-output file " + su_realisations_filename);
		  }
		if (verbosity>1)
		  {
			EH.println("Wrote "+su_realisations_filename);
		  }

	  }
	if (drawMAPmodels || drawMinChisqModels)
	  {
		try 
		  {
			suwriter.flush(realisations_MAP_out_strm); // must do this!
			realisations_MAP_out_strm.close();
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close MAP-output file " + su_MAP_realisation_filename);
		  }
		if (verbosity>2)
		  {
			EH.println("Wrote "+su_MAP_realisation_filename);
		  }
	  }
  }

  // make clean start and end range times: must be divisible by dt
  private void cleanSyntheticRange() 
  {
	double EPS=1e-7; // floating point machine_epsilon, 
	// since tSyntheticStart might have been mangled in parsing
	tSyntheticStart = Math.floor(tSyntheticStart*(1+EPS)/syntheticSampleTime)*syntheticSampleTime;
	tSyntheticEnd = Math.floor(tSyntheticEnd*(1+EPS)/syntheticSampleTime)*syntheticSampleTime;
  }

  // Simple checks to see if model and seismic traces are present if required
  private boolean checkFilesRequired()
  {
	if (m.MPD.PriorIndexInfo.size()>0)
	  {
		// we have a prior mdel that varies with each trace (prior traces necessary)
		if (su_model_filename==null) 
		  {
			EH.println("ERROR: you must provide an SU " +
					   "model-traces file if properties vary%n");
			return false;
		  }
	  }
		
	return true;
  }

  private boolean loopOverTraces() throws Exception
  {
	int ntr;
	// loop over seismic and model traces 
	sureader = new SuReader();
	tr_model = new SuSegyFormat();
	boolean skip;

	// check the realisation_number header field is OK
	if (!tr_model.set(MF.realisationNumberHeaderField,(short)0))
		throw new RuntimeException("Header field " + MF.realisationNumberHeaderField +
								   " can't be used for the realisation number");

	// write the output format version
	if (!tr_model.set(ModelFormat.VERSION_NUMBER_HEADER_FIELD, (short)MF.format_version))
		throw new RuntimeException("Header field " + ModelFormat.VERSION_NUMBER_HEADER_FIELD +
								   " can't be used for the format version number");
	  

	n_traces=0;


	if (makeSyntheticAndWavelet) 
	  {
		// read new data traces from model : there is no seismic to compare to, or,
		// we don't care about the seismic for the purposes of generating synthetics
		while ((ntr = populateFromModelTrace(false,false))!=-1)
		  {
			n_traces++;
			// skip trace if required
			if (traceFilter.skipTraceIfRequired(tr_near)) continue;
			if (haveZeroTimesFlaggedSkipTrace) 
			  {
				if (verbosity>2) m._SI.EH.println("TRACE trapped zero times and skipped");
				continue;
			  }
			m.setUpReducedModel();
			// adjust prior times if required
			m.resetPriorTimesFromDepthsIfNeeded();
			m.generateForwardModelsFromPrior();
			if (ntr==0) break;  
			// this happens if model is constant => there are 
			// no prior traces, might as well do this only once
		  }		
	  }
	else
	  {
		// Drive loop reading as long as there are seismic traces
		// loop over nears for want of a better choice
		while (sureader.fgettr(near_stack_su_in_strm,tr_near) != 0) 
		{
		  n_traces++;

		  skip = traceFilter.skipTraceIfRequired(tr_near);

		  // read new model traces;
		  populateFromModelTrace(true,skip);

		  // read new far traces;
		  getFarStackTrace();
	
		  // setup stack data arrays
		  setUpStackDataArrays();

		  // dump package on first acceptable trace if required
		  if (!skip && doDumpDebug) 
			{
			  EH.println("Dumping debug package on first accepted trace in trace-filter.%n");
			  formDebugPackage(null); 
			  return false;
			}
		  // skip trace if required
		  if (skip) continue;

		  if (haveZeroTimesFlaggedSkipTrace) 
			{
			  if (verbosity>1) EH.println("TRACE trapped zero times and skipped...");
			  continue;
			}

		  // setting up the prior model
		  m.setUpReducedModel();
		  
		  // adjust prior times if required
		  m.resetPriorTimesFromDepthsIfNeeded();

		  
		  // do one of these....
		  if (checkPriorSanity)
			{ 
			  /* do no mode searching or sampling. Checks mainly fluid substiutions
				 and sanity of rock velocity relations, as per setUpReducedModel()
			  */
			}
		  else if (drawPriorPlusIsopachModel || drawPriorPlusIsopachModelForAllFluidCombos)
			{
			  // draw forward models from prior
			  m.generateForwardModelsFromPrior();
			}
		  else
			{
			  // the main inversion is here
			  m.drawSamples();
			}
		}
	  }
		
	return true;
  }

  void getFarStackTrace()
  {
	// if we've specified far-stack data, get the matching trace
	if (m.SM.getFarstackData()!=null)
	  {
		if (m.SM.getFarstackData().filename!=null)
		  {
			if (sureader.fgettr(far_stack_su_in_strm,tr_far)==0)
			  throw new RuntimeException("Ran out of far-stack traces in " + 
										 su_far_offset_filename + " at trace " +
										 n_traces);
		  }
	  }
  }

  private void setUpStackDataArrays()
  {
	StackData stackData;

	// set up ptr for near-stack data 
	if ((stackData = m.SM.getNearstackData())!=null)
	  {
		if (stackData.filename!=null)
		  {
			stackData.currentTrace = tr_near.get_data_array();
		  }
	  }

	// set up ptr for far-stack data 
	if ((stackData = m.SM.getFarstackData())!=null)
	  {
		if (stackData.filename!=null)
		  {
			stackData.currentTrace = tr_far.get_data_array();

		  }
	  }
  }


  /** 
	  Efficiency initialisations required for all models,
	  trace-varying or not
  */
  private void doSimpleInitialisations()
  {
	// Numerical efficiency initialisations
	m.computeEfficiencyInitialisations();
	// And synthetic window
	m.computeSyntheticWindowInitialisations();	
  }

  /**
   update model trace from prior trace
  */
  int populateFromModelTrace(boolean mustBeThere,boolean skip) throws RuntimeException
  {
	ModelPriorDescriptions.MetaInfo MI;


	haveZeroTimesFlaggedSkipTrace = false;

	// this survives as true if traces import zero times, and we skip
	// the inversion

	// if everything's constant, no need to do anything
	if ( !modelTracesAreProvided || m.MPD.PriorIndexInfo.size() == 0) 
	  {
		doSimpleInitialisations();
		return 0; 
	  }
	// otherwise, we need input model traces

	int n_model = sureader.fgettr(model_in_strm,tr_model);
	if (n_model == 0)
	  {
		if (mustBeThere)   // must have dipped out 
		  {
			throw new RuntimeException("too few model traces in file " 
									   + this.su_model_filename + 
									   ": missing for traces >=  " + this.n_traces);
		  }
		else 
		  return -1;
	  }
	// Rudimentary checks. Check trace size etc
	int ns = (int)tr_model.get_short("ns");
	
	// The largest index will be the block index corresponding to t_base
	// times the number of layers
	int index_max = m.MPD.t_base.index * m.n_layers;
	if ((index_max+1) > ns)	  
	  throw new RuntimeException("model traces in file " + su_model_filename + 
								 " have fewer properties than specified by " + 
								 modelDescription_xml_filename +
								 ". Variable prior specification requires " + (index_max+1) +
								 " properties in the model trace, only " + ns + " are supplied");
		
	// If we're skipping, bail out now (we've read the trace)
	if (skip) return 1;

	// load in needed values
	ModelPriorDescriptions.propertyInfo info = new ModelPriorDescriptions.propertyInfo();
	// read data  from trace
	float values[] = tr_model.get_data_array();


	int index;
	// We expect blocks of m.n_layers properties in the input traces.
	// The last block will be times, of size (m.n_layers+1) if present

	// assume we will skip this trace until non-zero times prove otherwise. Activate this 
	// only if prior_times_vary
	if (m.MPD.prior_times_vary) haveZeroTimesFlaggedSkipTrace = true;

	// loop over property - indices
	for (int i=0;i<m.MPD.PriorIndexInfo.size();i++)
	  {
		info = (ModelPriorDescriptions.propertyInfo) m.MPD.PriorIndexInfo.get(i);
		if (info.PD == m.MPD.t_base) 
		  // base times are special. This will be the last possible entry
		  {
			index = info.PD.index * m.n_layers;
		  }
		else
		  {
			// other properties (including sigma_t's). Normal square offset formula 
			index = (info.PD.index-1) * m.n_layers + info.layer_no;
		  }
		if (index+1 > ns)	  
		  throw new RuntimeException("model traces in file " + su_model_filename + 
									 " don't have properties for index " + index);
		double value = (double)values[index]; 
		//EH.println("From model-trace: block-index = "+info.PD.index+", pos. in trace = "+ index+", value = "+value);

		// update prior if the value read in is not a 'bad-number'
		if (value != MODEL_TRACE_BAD_NUMBER)
		  {
			// Set of possible patches to priors as they come in.......

			// (1) Time-units autodetection -> map ms data to seconds
			MI = m.MPD.MetaInfo[info.layer_no];
			if (MODEL_TRACE_TIME_UNITS_AUTO_DETECTION)
			  {
				if (info.PD == MI.t || info.PD == m.MPD.t_base)
				  {
					if (value>TIME_UNITS_AUTO_DETECTION_TIME_MAX)
					  {
						// input probably (as it ought to be) in ms. convert to seconds
						value*=1.0e-3;
					  }
				  }
				if (info.PD == MI.sigma_t || info.PD == m.MPD.sigma_t_base)
				  {
					// Disable test .***** sigma_t must be in ms. Conversion always applies ******
					// if (value>TIME_UNITS_AUTO_DETECTION_TIME_MAX) value*=1.0e-3;
					value*=1.0e-3;
				  }
			  }

			

			// (2) Skip trace detection
			if (info.PD == MI.t || info.PD == m.MPD.t_base)
			  {
				if (value!=0) 
				  {
					// a non zero was detected, so we mustn't be skipping
					haveZeroTimesFlaggedSkipTrace = false;
				  }
			  }


			// (3) Depth handling
			// (3a) Depths must be positive
			if (info.PD == MI.depth && MAKE_MODEL_TRACE_PRIOR_DEPTHS_POSITIVE)
			  {
				if (value<0) value = Math.abs(value);
			  }

			// (4) Isopach protection
			// if ((info.PD == MI.sigma_dz) && value<0) 
			if (info.PD == MI.sigma_dz)
			  {
				if (value!=0 && value<MINIMUM_ACCEPTABLE_ISOPACH_UNCERTAINTY)
				  {
					// a non physical value was detected. Zero is OK (forced pinchout). Threshold it
					EH.println("Warning: isopach uncertainty sigma_dz unphysical: corrected to threshold "
							   +MINIMUM_ACCEPTABLE_ISOPACH_UNCERTAINTY);
					value = MINIMUM_ACCEPTABLE_ISOPACH_UNCERTAINTY;
				  }
			  }

			// (5) Sigma_t protection
			if (((info.PD == MI.sigma_t) || (info.PD == m.MPD.sigma_t_base)) && value<=0) 
			  {
				// a non physical value was detected. Threshold it
				throw new RuntimeException("Time uncertainty sigma_t unphysical: must be >=0 ");
			  }

			// (6) Signal-noise ratio protection
			if (m.SM.doSNchecks) 
			  { 
				if (!m.SM.isSNratioOK())
				  {
					throw new RuntimeException("Signal to noise ratio is unrealistically large. "+
											   "Check specified noise level in XML or model traces!!");
				  }
			  }
			
			// update prior information in ModelPriorDescriptions
			info.PD.setValue(value);
		  }
	  }
	// 1b Depth locking...
	if (lockDepths)
	  {
		double depth_value = m.MPD.MetaInfo[m.master_depth_layer_number-1].depth.value;
		for (int i=0;i<m.n_layers;i++)
		  {
			m.MPD.MetaInfo[i].depth.setValue(depth_value);
		  }
	  }
	// Forced pinchouts
	m.initialiseLayerSynching();

	doSimpleInitialisations();

	if (verbosity>0)  EH.println("%nTRACE number " + n_traces);

	if (verbosity>2) showPriorInfo(true);

	return 1;
  }

  void showPriorInfo(boolean isUpdated)
  {
	// get location from seismic if we have it, or prior if we don't
	printLocation((near_stack_su_in_strm != null) ? tr_near : tr_model);
	m.MPD.printPriorInfo(isUpdated);
  }


  /**
	 Dump out to stdout the ep, cdp, and fldr locations in the seismic data, 
	 for help in debugging when needed. Gets location from trace suggested
	 (probably seismic if available, or model if no seismic)
  */
  private void printLocation(SuSegyFormat tr)
  {
	int val;
	if ((val=tr.get_int("fldr"))!=0) EH.printf("fldr = "+val+", ");
	if ((val=tr.get_int("ep"))!=0)   EH.printf("ep = "+val+", ");
	if ((val=tr.get_int("cdp"))!=0)  EH.println("cdp = "+val);
  }

  
  // Dump out model realisation
  // Reuses the model traces header and data parts
  // A new state will be dumped, but a repetition will dump the
  // old values in the trace
  // A new state requires the values to be in-situ in the 
  // Model.LayerProperties member fields

  public void dumpRealisation(boolean newState, int realisation_number, boolean write_state, boolean dumpToMAPstrm)
  {
	float value=0.0f,value2=0.0f;
	double depth_ref=0.0;
	double seisError=0;
	double depths[];
	int i,ii,j,jj,k;
	boolean useFG,useSC;
	String name;

	// suppress dumping if required 
	if (suppressModelDumping) return;

	// If it's a new state, fill the trace with the new properties, 
	// Otherwise repeat the dump of what's already there (must be in place already)
	if (newState)
	  {
		LayerProperties.Property prop;


		// Allow model trace to inherit full info from near-stack seismic
		// This is ok, unless we are making synthetic, inwhich case there is no decent near stack header info
		if (!makeSyntheticAndWavelet) tr_model.copyHeader(tr_near);
		// version number must also always be correct
		tr_model.set(ModelFormat.VERSION_NUMBER_HEADER_FIELD, (short)MF.format_version);

		// Write how many samples there will be in the model
		tr_model.set("ns",(short)MF.ns);

		// might as well inherit inline/xline fields from model
		// setInlineXlineHeaderFields(tr_near, tr_model);

		// Set trace identification code for missing fields if any
		tr_model.set(ModelFormat.OUTPUT_FIELDS_SU_KEYWORD,MF.getFieldBitsSuppressed());

		// precompute all suspended or supported depths
		depths = new double[m.n_layers+1];
		if (outputParametricDepths)
		  {
			// original parameteric depths
			for (j = 0 ; j<m.n_layers; j++) depths[j] = m.LP[j].d.value;
		  }
		else
		  {
			for (j = 0 ; j<m.n_layers; j++) depths[j] = m.LP[j].suspended_depth_top;
		  }

		// Extra features
		useFG = MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL;
		useSC = MF.ENABLE_SIMPLE_CONTACTS_MODEL;

		// and samples
		k=0;
		for (i=0;i<MF.complete_field_names.length;i++)
		  {
			if (MF.includeField[i])
			  {
				String MFfieldNameI =  MF.complete_field_names[i];
				if (MFfieldNameI.equals("fluid_type"))
				  {
					// Dump an indication of fluid type. 0,1,2,3 for brine, oil, lsg, gas
					if (!doSimpleFluidUnitsAndContacts) 
					  {
						for (j=0;j<m.n_layers;j++) 
						  {						
							tr_model.set_data_value((float)m.MS.fluid_state.types[j], k++); 
						  }
					  }
					else
					  {
						for (j=0;j<m.n_layers;j++) 
						  {
							if (useSC && m.LP[j].fcInfo.is_permeable) 
							  {
								value = ((LayerProperties.InternalPermeableFluidContactInfo)m.LP[j].fcInfo).fluid_type_eff;
							  }
							else value = (float)m.MS.fluid_state.types[j]; // will be brine I'm sure
							tr_model.set_data_value(value, k++); 
						  }
					  }
				  }
				else if (MFfieldNameI.equals("vp_eff"))
				  {
					// Dump mixed-rock, fluid-substituted effective p-wave velocity of layer 
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.LP[j].vp_eff, k++); 
					  }
				  }
				else if (MFfieldNameI.equals("rho_eff"))
				  {
					// Dump mixed-rock, fluid-substituted effective density of layer
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.LP[j].rho_eff, k++);  
					  }
				  }
				else if (MFfieldNameI.equals("Z_eff"))
				  {
					// Dump mixed-rock, fluid-substituted effective Impedance (zero-offset)
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)(m.LP[j].vp_eff * m.LP[j].rho_eff), k++); 
					  }
				  }
				else if (MFfieldNameI.equals("vs_eff"))
				  {
					// Dump mixed-rock, fluid-substituted effective s-wave velocity of layer 
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.LP[j].vs_eff, k++); 
					  }
				  }
				else if (MFfieldNameI.equals("Zs_eff"))
				  {
					// Dump mixed-rock, fluid-substituted effective s-wave impedance of layer 
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)(m.LP[j].vs_eff * m.LP[j].rho_eff), k++);
					  }
				  }
				else if (MFfieldNameI.equals("rho_h"))
				  {
					// Dump an indication of Substituted Fluid density
					if (!doSimpleFluidUnitsAndContacts) 
					 {
					   for (j=0;j<m.n_layers;j++) 
						 {
						   if (m.MS.fluid_state.types[j] == FluidStates.FluidTypes.brine)
							 value = (float)m.LP[j].rho_b.value;
						   else
							 value = (float)m.LP[j].rho_h;
						   tr_model.set_data_value(value, k++);  
						 }
					 }
					else
					  {
						// contacts model. set directly
						for (j=0;j<m.n_layers;j++) tr_model.set_data_value((float)m.LP[j].rho_h, k++);  	
					  }

				  }
				else if (MFfieldNameI.equals("S_h"))
				  {
					// Dump an indication of Substituted Fluid saturation
					if (!doSimpleFluidUnitsAndContacts) 
					  {
						for (j=0;j<m.n_layers;j++) 
						  {
							// set zero here if the fluid state is brine
							if (m.MS.fluid_state.types[j] == FluidStates.FluidTypes.brine) 
							  value=(float)0.0;
							else
							  value = (float)m.LP[j].S_h;
							tr_model.set_data_value(value, k++); 
						  }
					  }
					else
					  {
						// contacts model. set directly
						for (j=0;j<m.n_layers;j++) tr_model.set_data_value((float)m.LP[j].S_h, k++);  	
					  }
				  }
				else if (MFfieldNameI.equals("S_oil"))
				  {
					// Dump an indication of oil saturation
					for (j=0;j<m.n_layers;j++) tr_model.set_data_value((float)m.LP[j].S_oil_eff, k++);  	
				  }
				else if (MFfieldNameI.equals("S_lsg"))
				  {
					// Dump an indication of lsg saturation
					for (j=0;j<m.n_layers;j++) tr_model.set_data_value((float)m.LP[j].S_lsg_eff, k++);  	
				  }
				else if (MFfieldNameI.equals("S_gas"))
				  {
					// Dump an indication of gas saturation
					for (j=0;j<m.n_layers;j++) tr_model.set_data_value((float)m.LP[j].S_gas_eff, k++);  	
				  }
				else if (MFfieldNameI.equals("R_near"))
				  {
					// Dump near reflection coefficients
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.SM.R_near[j], k++); 
					  }
				  }
				else if (MFfieldNameI.equals("R_far"))
				  {
					// Dump near reflection coefficients
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.SM.R_far[j], k++); 
					  }
				  }
				else if (useSC && MFfieldNameI.equals("t_fc"))
				  // Suite for fluid-contacts
				  {
					// Dump near reflection coefficients times in ms
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)(1000*m.SM.t_fc[j]), k++); 
					  }
				  }
				else if (useSC && MFfieldNameI.equals("R_near_fc"))
				  {
					// Dump near reflection coefficients
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.SM.R_near_fc[j], k++); 
					  }
				  }
				else if (useSC && MFfieldNameI.equals("R_far_fc"))
				  {
					// Dump near reflection coefficients
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.SM.R_far_fc[j], k++); 
					  }
				  }
				else if (useSC && MFfieldNameI.equals("d_fc"))
				  {
					// Dump contact depths: total number guaranteed to be less than
					// num layers. Dump out in unrolled fashion
					j=0;
					if (m._SI.doSimpleFluidUnitsAndContacts)
					  {
						for (ii=0; ii<m.fluidUnitsAndContacts.fluidUnits.length; ii++)
						  {
							FluidUnitsAndContacts.FluidUnit fu = m.fluidUnitsAndContacts.fluidUnits[ii];
							for (jj=0; jj<fu.contactInfo.length; jj++)
							  {
								tr_model.set_data_value((float)fu.contactInfo[jj].depth.value, k++); 
								j++;
							  }
						  }
					  }
					// zero pad...
					for (;j<m.n_layers;j++) tr_model.set_data_value((float)0.0, k++); 
				  }
				else if (useFG && MFfieldNameI.equals("float_fraction"))
				  {
					// Terms here for floating fraction and adjusted porosity
					// Dump float fractions
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)m.LP[j].float_fraction.value, k++); 
					  }
				  }
				else if (useFG && MFfieldNameI.equals("phi_s"))
				  {					
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)Math.max(m.LP[j].phi_s.value, 0), k++); 
					  }
				  }
				else if (MFfieldNameI.equals("seismic_chisquare"))
				  {
					if (m.SM.seismic_error_vector!=null)
					  {
						// Equivalent number of measurements...
						double nu = m.SM.seismic_error_vector.length * m.SM.delta_t_ratio_correction_factor;
						// m.MS.negLogLikSeisOnly = 0.5*sum {errors^2}, so undo factor of 2
						seisError = (writeHiddenSeismicError ? m.MS.negLogLikHiddenSeisOnly : m.MS.negLogLikSeisOnly);
						value = (float)(2.0*seisError/nu);
						// value = (float)seisError;


						// Usable full chisquare for post run model comparisons...
						// chisq = negLogLikSeisOnly + negLogIsopachLikelihood 
						//         + negLogDepthLikelihood + negLogVStackLikelihood + negLogPrior;
						// full chisq of posterior0
						// value2 is chisq...
						value2 = (float)(m.MS.mode.lastAcceptedNegLogAllLikelihoods + m.MS.mode.lastAcceptedNegLogLikPriorOnly);
						if (writeHiddenSeismicError) value2 += m.MS.negLogLikHiddenSeisOnly - m.MS.negLogLikSeisOnly;
					    value2 += m.MS.mode.normalisationConstants; // normalisation from seismic, normalisation from isopachs
						                   // normalisation from vertical stacking *likelihood* (thicknesses)
						                   // note (1/2)log(2Pi) factors cancel here
					  }
					else
					  {
						value=0.0f;
					  }
					// emergency marker for failed models
					if (m.MS.noModesFound) value = value2 = (float)NO_MODES_FOUND_CHISQ_BIGMARKER; // -999s
					// seismic chisq in layer 1,2....
					for (j=0;j<(m.n_layers-1);j++) tr_model.set_data_value(value, k++);
					// full posterior chisq in last layer
					tr_model.set_data_value(value2, k++); 
				  }
				else if (MFfieldNameI.equals("diagnostics"))
				  {
					// Dump Raftery thinning parameter as indication of chain sanity
					if (m.MS.CM!=null && m.MS.CM.RafteryStyleDiagnostics)
					  {
						tr_model.set_data_value(0.0f, k++); // top layer has no Raftery diagnostic. Internal times only
						for (j=1;j<m.n_layers;j++) tr_model.set_data_value((float)m.MS.CM.Raftery_kthin[j-1], k++); 
					  }
					else
					  {
						for (j=0;j<m.n_layers;j++) tr_model.set_data_value(0.0f, k++); 
					  }
				  }
				else if (MFfieldNameI.equals("t"))
				  // Last block. May be depth or time.......
				  // Times must be in ms
				  {
					if (!outputDepthModels)
					  {
						for (j=0;j<m.n_layers;j++) 
						  {
							tr_model.set_data_value((float)(1000*m.LP[j].t.value), k++); 
						  }
					  }
					else
					  {
						// and write these
						for (j = 0; j<m.n_layers;j++) 
						  {
							tr_model.set_data_value((float)depths[j], k++);
						  }
					  }
				  }
				else if (MFfieldNameI.equals("t_base"))
				  // Supplement to last block. May be depth or time.......
				  {
					if (!outputDepthModels)
					  {
						// Only one entry for the base time of the bottom layer
						// time in ms
						tr_model.set_data_value((float)(1000f*m.LP[m.n_layers-1].t_base.value), k++);
					  }
					else
					  {
						// If we're dumping depth models, the last block has to be depth
						tr_model.set_data_value((float)depths[m.n_layers], k++);
					  }
				  }
				else if (MFfieldNameI.equals("d"))
				  // depth field has to be overwritten with hung depths computed from delta t * vp_effs
				  {
					for (j=0;j<m.n_layers;j++) 
					  {
						tr_model.set_data_value((float)depths[j], k++);
					  }
				  }
				else
				  {
					// All other ordinary blocks: vp_m, vp_s etc				
					name = MF.complete_field_names[i];
					for (j=0;j<m.n_layers;j++) 
					  {
						prop = (LayerProperties.Property) m.LP[j].completePropertyNameMap.get(name);
						if (prop==null)	throw new RuntimeException("No such field in dumpRealisation(). Fix this bug!");
						tr_model.set_data_value((float)prop.value, k++);
					  }
				  }
			  }
		  }

	  }



	// set trace identification to the realisation number. 
	// this is done whether to state is new or not
	// Assume we can't have realisation_number  > Short.MAX_VALUE
	tr_model.set(MF.realisationNumberHeaderField,(short)realisation_number);	

	// Now dump out trace if insisted on
	if (write_state) 
	  {
		if (!dumpToMAPstrm)
		  {
			suwriter.fputtr(realisations_out_strm,tr_model);
		  }
		else
		  {
			suwriter.fputtr(realisations_MAP_out_strm,tr_model);
		  }
	  }
  }


  /**
	 Hook to write external command script: typically a command that requires
	 the string "properties=A,B,C,.." appended to something like BHPwrite.
	 This method appends "properties=A,B,C,.." to BHPcommandstring, writes
	 it to BHPcommand_filename, and chmods this to "ugo+rx".
	 A,B,C etc are the properties formed from MF.complete_field_names, but the
	 last, t_base, is omitted.
  */

  void handleBHPcommand()
  {
	PrintWriter out = null;

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
	// Leave out t & t-base here, (always last entry in property_string)
	int N = MF.complete_field_names.length-2;
	for (int i=0;i<N;i++) 
	  {
		property_string.append((String)MF.complete_field_names[i]);
		if (i<N-1) property_string.append(",");		
	  }

	out.printf(BHPcommandstring+property_string+" %n");
	out.close();
	EH.warn("BHPcommand: Dumped keywords for "+N+" properties");
	// Unsure whether we need this
	// ScriptGenerator.chmod(BHPcommand_filename,"ugo+rx",EH);
  }

}
