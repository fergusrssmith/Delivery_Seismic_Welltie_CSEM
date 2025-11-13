/*************************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
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
//------------------------- EOCOPYRIGHT --------------------------------- 

package au.csiro.JamesGunning.CSEM.OneDinversion;

import java.io.*;
import java.nio.*;
import java.util.*;
import static java.lang.Math.*;

import cern.colt.map.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;  
import cern.colt.matrix.linalg.*;
import cern.colt.list.*; 
import cern.jet.math.*;
 

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.Numerical.ModifiedCholeskyLeftDecomposition;
import au.csiro.JamesGunning.lib.Numerical.ModifiedCholeskyRightDecomposition;
import au.csiro.JamesGunning.lib.Numerical.DenseFloatMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.DenseIntMatrix2D;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.util.Util.*;

import au.csiro.JamesGunning.lib.Minpack.optimization.*;

import static au.csiro.JamesGunning.lib.util.MathUtil.*;


/** 
   Java Bayesian 1D CSEM code, using forward computation engine from
   Scripps consortium fortran codes. Has variety of inversion features
   and methodologies.  Basic version uses startup, data and model
   files as per legacy Scripps code (the old f77 "OCCAM" code).  The
   --NS flag supports files in formats for new f90 Scripps kernel
   (Kerry Key's Dipole1D engine) See documentation, self docs.
   
   The bulk of the comments correspond to the mathematics in the "deliveryCSEM"
   paper.

   @author James.Gunning@csiro.au
   @version 1.0, 10/2/2008
*/


public class CSEM1DInverter
{
   /**
	Lumped transmitter height bucket size (in m)
  */
  public static double LTH_BUCKET_SIZE = 10.0;
  
  /**
	 Difference in mode max. posterior values less than this (for the same
	 fluid type) is taken as evidence the modes are in fact the same,
	 and they are amalgamated. This kind of energy barrier is easily
	 overcome in diffusion (easy tunneling)
  */
  public static double MODE_SEPARATION_EPSILON = 0.5;

  /**
	 Range of negative log posterior values over which distinct
	 local optima may be regarded as of "comparable significance". MCMC
	 sampling scales are computed from average covariances of the modes
	 which fall within this window
  */
  public static double NEG_LOG_POSTERIOR_PLAUSIBLE_MODE_GROUP_RANGE = 3.0;

  /**
	 Scaled separation of log resistivity solution indicating probable
	 dislocation in solution
  */
  public static double MULTISTART_DIST_THRESH = 1.0;

  /**
	 Do Bayes mu update only every BAYES_SMOOTHING_MU_UPDATE_SKIP
	 iterations (allows other parameters to relax a bit)
  */
  public int BAYES_SMOOTHING_MU_UPDATE_SKIP = 2;

  /**
	 Starting value of log_10(mu) in auto-smoothing mode.
	 Overruled by Hessian eigenvalue balance calculation if possible.
  */
  public float PMU_DEFAULT_START = 5.0f;
  /**
	 When setting of pmu = log_10(mu) is performed by 
	 eigenvalue balances, this is the ratio of the
	 smoothing term to the sensitivity Hessian. I.e. initial  
	 guess is to set <p>
	 smoothing_part_of_Hessian = INITIAL_SMOOTHING_EXCESS_WEIGHT * sensitivity_part_of_Hessian
  */
  float INITIAL_SMOOTHING_EXCESS_WEIGHT = 10000.0f;  
  /**
	 Safe prior std deviation for models where this is set to zero 
	 (used in prior covariance deteriminant only)
  */
  public static float DEFAULT_SAFE_INV_RHO_PRIOR_STD_DEV = 0.2f;  //  1/10.0
  
  /**
	 Power of 2 used to generate underlying grid for splitting methods
  */
  public static int LOG2_MAX_NUM_SPLIT_GRID_ELEMENTS = 5;  //  split-grids can have as many as 2^5 = 32

  /** 
	  Geometric factor in logarithmic gridding
  */
  public static float GEOMETRIC_LOGMODEL_FACTOR = 1.15f;
  
  /**
	 Useful constant
  */
  static float SQRT2_180_ON_PI = (float)sqrt(2.0)*180f/(float)Math.PI;
  

  /** 
	  Defaults constants for termination criteria
  */
  public static float EPS_CHISQ = 0.001f;
  public static float SQRT_EPS_CHISQ = (float)sqrt(EPS_CHISQ);
  public static float CUBE_ROOT_EPS_CHISQ = (float)pow(EPS_CHISQ,1.0/3.0);

  /** 
	  Here, the tau_f "number of sig. figures in f"
	  precision constant discussed by Gill, Murray & Wright, section 8.2.3.6
	  Very important constant for termination of Marquardt methods, effectively
	  setting the small-gradient threshold etc. See also 8.2.3.2 U1,2,3 conditions
  */
  public static float GMW_TAU_F= 5e-5f;

  float GMW_TAU_F_PWR_ONETHIRD = (float)pow(GMW_TAU_F,0.333333333);

  au.csiro.JamesGunning.lib.ErrorHandler EH;
  Date date = new Date();
  String classname=null;
  PrintWriter out_logger = null;
  String log_filename = "logfil";

  public static boolean debug = false;

  String altStartFile = null;
  boolean doBayesianAutoSmoothing=false;
  float muTrustSize = 0.5f;
  boolean bayesianAutoSmoothingExpStyle = true;
  // For extra term in optimisation of mu, compared to pure ML theory:
  boolean doEmpBayes = true;   

  boolean doNoiseCooling=false;
  boolean doNoiseEstimation=false;
  boolean doRelaxedSmoothing=false;



  boolean prefitLinearTrend=false;
  boolean doMCMC = false;
  boolean doModelSelectAutoSplitting = false;
  boolean splitSuggest = false;

  boolean makeSynthetic = false;
  String syntheticSpecFloor=null, syntheticSpecFrac = null;
  NoiseSpec nSpec = null;

  int N_MCMC = 0;
  String samplesFilename = null;
  boolean doMDstyleMCMC = false; 
  boolean doMultistart = false;
  boolean multistartFlipping,multistartRandom;
  int multistartRandomN=0;
  boolean doXsections = false;
  boolean dumpPrfl = false;
  boolean quietIO = false;
  boolean zeroSmoothingWeights = false;
  boolean THApproxTest = false;
  boolean approximateTxHeights = true;
  boolean estDataErrors = false;
  boolean graphDataErrors = false;
  boolean useBoundConstraints = true;
  float option_bounds_lower = -Float.MAX_VALUE;
  float option_bounds_upper = Float.MAX_VALUE;
  boolean dumpConvDiagnostics = false;
  boolean useSVD = false;  
  float eps_SVD = 1e-6f;
  boolean doPathIntegralTests = false;

  boolean overwriteLogModel=false;
  float logModelFirstThickness=100f;
  float logModelRatio = GEOMETRIC_LOGMODEL_FACTOR;
  boolean readProfile=false;
  String readProfileFilename= null;
  boolean readProfileSetMean = false;

  boolean resampleMeansForBayesianBootstrapping = false;
  boolean doBootstrapping = false;
  int N_PB = 0;
  boolean doBootstrappingRecenter=false;
  boolean checkMLmodelOnNewData = false;
  String checkNewDataFile = null;

  boolean doModelStudies=false;
  String modelStudy=null;
  ModelStudyList[] modelStudyList=null;

  boolean doResamplingBootstrap=false;
  int NBootstrap=0;

  boolean modelSystematics=false;
  String systematicsFlavour=null;
  SystematicErrorsModel SEM  = new SystematicErrorsModel();

  boolean dilute=false;
  int Ndilute=0;
  float diluteFactor=1.0f;
  
  boolean doWishart=false;
  float Wishart_nu=10f;

  boolean doShrinkage=false;
  float shrinkage_alpha=0.1f;

  int verbosity=1;
  boolean dumpDebug=false;

  // stuff for new Scripps kernel
  public static boolean newScripps = false;
  ScrippsModelInterface SMI;

  // Experimental Marquardt: not finished yet. 
  // 1) Terminate with lambda=0, so matrices unpolluted
  // 2) Ensure starting lambda guesses are right for all modes, especially --MS, -r, -s
  boolean useLMstrategy = false;



  // Collected optima
  ArrayList<LocalMode>optima;

  // Useful stuff for doModelSelectAutoSplitting
  ModelSelectInfo msInfo=null;


  // useful globals
  PriorTrendInfo PTI = new PriorTrendInfo(this);
  float global_error_scale = -1.0f;
  float[] depths;
  byte[] act_PN=null;  // projected Newton working actives
  float[] m_lower;     // bounds
  float[] m_upper;     // bounds
  float[] grad;
  float[] WJdeltaD;
  float minSkinDepth,freqMin;
  boolean haveLargestEigenvaluePTP=false;
  double largestEigenvaluePTP=0.0;
  boolean haveLargestEigenvalueWJTWJ = false;
  double largestEigenvalueWJTWJ = 0.0;
  boolean avoidTerminationBayesianAutoSmoothing = false;
  LevMarInfo LMinfo = null;

  ArrayList<cern.colt.list.FloatArrayList>convergenceResiduals;

  // Current filenames
  String current_profile_filename =null;
  String current_profile_err_filename = null;
  String current_resp_filename = null;
  String current_iter_filename = null;

  // JNI pointer management
  boolean nativeModelPointersSet = false;

  /* Global parameters across itrin, main loop, itout */
  String itform;
  String descr;
  String modfil;
  String datfil;
  int maxitr=0;
  float tolreq=0.0f;  // JG: unused except in relaxed smoothing
  int iruf;
  int nit=0;
  int nit_ms=0;
  float pmu;
  float rlast;
  int konv=0;
  float tobt_in;
  int ifftol;         // JG : unused
  /* End itrin, main loop, itout globals */

  DenseFloatMatrix2D WJ_trans;  // contiguous memory  version of float[][] WJ_trans, Jacobean
  DenseDoubleMatrix2D del; // will be (np,np)
  DenseDoubleMatrix2D del2=null; // will be (np,np) if iruf==2
  /* Material correponding to occam commons (revision 3.00, 26 may 2001) */
  // include numbers
  // int npp,ndd;
  long tNatives=0;
  long tOverall=0;

  // data
  int nd, ndp;
  float[] d, sd;         // size nd
  float[] sd_orig=null;  // optional copy of sd, possibly needed for PB
  // float[][] dp;       // size [6][nd]
  DenseFloatMatrix2D dp; // contiguous-memory  version of dp[6][nd]
  float[] origTransHeights; // keeps original data since we use lumping

  // model
  int np;
  float[] pm;       // size np
  float[] dm;       // size nd
  float[] dm_true;  // points to dm if no systematics, else fwd model
  // result
  int idebug, nfor, iout;
  DenseDoubleMatrix2D ptp,wjtwj; // size [np][np] 
  DenseDoubleMatrix2D H; // size [Nptot][Nptot], of which H_11 is wjtwj
  float[] pwk1, pwk2, pwk3, premod, prewts; // size [np]
  float[] dwk1, dwk3;  // size [nd]
  float frac;

  // fmodel
  // In Occam commons: thickf(0:npp),rho(0:npp),freqs(ndd),
  // More correctly (I think, in modules.f90):
  // thickf_ptr(0:nlayer), rho_ptr(0:nlayer), freqs(nd)
  // ipen1_ptr(0:nlayer), ipen2_ptr(0:nlayer)
  float[] thickf, rho, freqs, ipen1,ipen2;
  int[] iff;
  float thite,tangle;
  int nfreq,nlayer;

  // Useful fields interaction with new Scripps kernel
  // new: number of free sealayers (0 in Legacy, ? in New)
  int nFreeSeaLayers=0;
  int nFixedSurfaceLayers=0;
  
  // Sampler
  float ms_MCMC,msprej_MCMC;

  Random RD = new Random();

  /** 
	  main driver from command line
  */
  public static void main(String argv[]) 
  {
	CSEM1DInverter CI = new CSEM1DInverter();
	ErrorHandler _EH = CI.EH;
	String s;
	String[] opts;

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("CSEM1DInverter",CI.EH);
	options.addLicense(DeliveryLicense.csem_license_text);
	options.addVersion(VersionNumber.versionNumber);
	options.addSpielText("CSEM1DInverter runs CSEM 1D inversion in working directory with "
						 +"startup,model, and data files; formats currently in (default) Scripps legacy-code OCCAM style."
						 +"New Scripps formats are usable under the --NS option, but not heavily tested");

	options.addOption(null,"--startup",1,"String","startup_file","Use specified startup file, not default 'startup'");
	options.addOption(null,"--MS",1,"String","[d][:r,number]","Multistart global optimization: try the union of these strategies"
					  +" 1) d=default (layer-flipping) 2) r,N (N pure random start-points in bounds)"
					  +" e.g. --MS d:r,100");
	options.addOption(null,"--auto-split",0,null,null,"Model-select by auto-splitting");
	options.addOption(null,"--split-suggest",0,null,null,"Create suitable model file (model_suggestion.mod) based on auto--splitting results");
	options.addOption("-v",null,1,"int"," verbosity","level of verboseness");
	options.addOption("-s",null,0,null,null,"Bayesian auto-smoothing mode (\"Constrained Bayesian inversion...spatial correlation\" in paper)");
	options.addOption("-r",null,0,null,null,"Do Newman-like 'gradually relaxed smoothing' mode. Smoothing parameter geometrically relaxed until fit is adequate.");
	options.addOption("-c",null,0,null,null,"Do 'noise level cooling'. Noise level geometrically relaxed. A kind of controlled Marquardt approach",1);
	options.addOption("-n",null,0,null,null,"Do simultaneous noise level estimation. Noise variance an additonal parameter, not fixed by input files.");
	options.addOption("-q",null,0,null,null,"Quiet I/O: dump only final useful files (i.e. at optima)");
	options.addOption("-lt",null,0,null,null,"Roughly fit linear trend as prior model before inverting (overrides prior mean/prejudice in files)");
	options.addOption("-N",null,1,"int"," N samples","Draw N MCMC posterior samples");
	options.addOption("-oSamp",null,1,"String"," filename ","MCMC or Bootstrapping output filename");
	options.addOption(null,"--MD",0,null,null,"Use MD (molecular dynamics) method for MCMC proposals. Experimental.",1);
	options.addOption(null,"--bounds",2,"float","lower upper","Define new bounds on log_10(resistivity) (default [-0.1,4])");
	options.addOption(null,"--NB",0,null,null,"Disable bound constraints on log_10(resistivity) (default [-0.1,4])");
	options.addOption(null,"--make-synthetic",2,"String","errorE[,errorEphase_degrees,errorB,errorBphase_degrees] noise_floor[,noise_floorB]",
					  "Make synthetic data set based on model in startup file, "
					  +"geometry in data file. Add Gaussian error of errorX%% (relative) to 'measurements'. And "
					  +"no data error less than noise_floor. If using B field too, use e.g. --make-synthetic 5,5,10,10 1e-15,1e-18");
	options.addOption(null,"--noise-levels",2,"String","errorE[,errorEphase_degrees,errorB,errorBphase_degrees] noise_floor[,noise_floorB]",
					  "Set noise levels for e.g. uncontrolled bootstrap. Meanings as per --make-synthetic");

	options.addOption(null,"--Xsections",0,null,null,"Dump file X-sections of -log posterior at optima");
	options.addOption("-P",null,0,null,null,"Dump file depth-profiles of solutions while optimizing");
	options.addOption(null,"--no-smoothing-weights",0,null,null,"No smoothing. Zero all smoothing weights");
	options.addOption(null,"--conv",0,null,null,"Dump convergence diagnostics file",1);
	options.addOption(null,"--TH-approx-test",0,null,null,"Test accuracy of lumped-transmitter-height approximations");
	options.addOption(null,"--no-TH-approx",0,null,null,"Do not perform lumped-transmitter-height approximations "+
					  "(i.e. use of average transmitter height over observation-groups). Can slow things down markedly...");
	options.addOption("-e",null,1,"String","q|g",
					  "Overwrite large data errors with crude median-filter error estimate: q for quietly, g for show-graphs. "
					  +"A quick'n'dirty outlier-management tool");
	options.addOption(null,"--log-grid",1,"String","first-layer-thickness[,ratio]","Overrule input file model with standard logarithmic grid. Supply extra ratio if desired [default="+GEOMETRIC_LOGMODEL_FACTOR+"]");
	options.addOption(null,"--read-profile",1,"String","profile_filename","Read existing solution from profile file, regrid on to existing model grid as starting point, and then invert as usual (not usable for auto-splitting mode)");
	options.addOption(null,"--read-profile-set-mean",1,"String","profile_filename",
					  "As per --read-profile, but also set prior mean from profile");
	options.addOption("--PB",null,1,"int","N-samples","Try to do full parametric bootstrap calculation with N samples");
	options.addOption("--PB-recenter",null,0,null,null,"Do recentering with --PB bootstrap resampling",1);
	options.addOption("--BPB",null,0,null,null,"Resample means for Bayesian Parametric Bootstrapping (automatic with --PB)");
	options.addOption("--check-new-data",null,1,"String"," filename ",
					  "verify fit against new data filename. Useful for bootstrap checks",1);
	options.addOption("--SVD",null,1,"float","eps","specify relative cutoff and use SVD solver. Suggest eps=1e-6 ?",1);
	options.addOption(null,"--s-trust",1,"float","step","sometimes the trust region (allowable step in mu) in -s mode needs "
					  +"to be smaller than the default 0.5: E.g. big or heterogeneous data sets. Suggest ? 0.1",1);
	options.addOption(null,"--NS",0,null,null,"Use new Scripps (Kerry Key's) kernel");
	options.addOption("--LM",null,0,null,null,"Use Levenberg Marquardt trust-region method, not Gauss-Newton line search. Often superior.",1);
	options.addOption("--PI",null,0,null,null,"Try to connect distinct modes using path-integral tests (expensive, and experimental)",1);
	options.addOption("--model-study",null,1,"String","spec1:spec2:...",
					  "Perform model study on layers (writes ModelStudy.txt). Each layer spec is of form "+
                      "lay_num,start_thickness,inc_thickness,end_thickness. (Remember "
					  +"lay_num=1=first-subsea-layer (in legacy Scripps model)",1);	
	options.addOption("--systematics",null,1,"String","spec"," Model systematic effects. Spec=colon separated lists of styles. "
					  +"Current and only style is d,Np,sigma_sys: order Np polynomials on |E|,|B|, sigma_sys = ratio of tolerable "
					  +"systematics to i.i.d Gaussian noise. See Appendix E.",1);	
	options.addOption("--Wishart",null,1,"float","nu","Use inv-Wishart formulation, prior degree-of-belief nu. C_0=diag{input file variances}. Experimental.",1);	
	options.addOption("--shrinkage",null,1,"float","alpha","Use covariance shrinkage based on Leung formula. Experimental.",1);	
	options.addOption("--dilute",null,1,"int","N","Reduce data to equivalently N 'effective' data. See Appendix E.");	
	options.addOption("--bootstrap",null,1,"int","N","Perform classical bootstrap (resample errors) study for a fixed model.",1);	
	options.addOption("--default-prior-constants",null,1,"String","A_mean,B_mean,A_sigma,B_sigma",
					  "Override default constants ("+CI.PTI.getConsts()+") for depth trends (used in splitting etc) of Bayesian prior: "+
					  "%n[ log_10(rho) ~ N(A_mean+B_mean z,(A_sigma+B_sigma z)^2), z in m]",1);
	options.addOption("--dump-debug-package",null,0,null,null,"Zip up all files for debugging if a problem occurs",1);
	options.addOption("--gdb-debug",null,0,null,null,"Wait for user prompt to continue (useful only for debugging)",1);
	options.addOption("-h",null,0,null,null,"Show help");

	// leave this to last!
	options.addExtraHelpOption("-X","Show help with all options (there are lots more!)");


	// parse options
	options.parseOptions(argv,true);


	if (options.isOptionSet("--gdb-debug"))
	  {
		Util.writeGDBJavaScript();
		Util.getStdinUserPrompt();  // delay till user acknowledges (for debugging)
	  }

	if (options.isOptionSet("-v")) 
	  {
		CI.verbosity = ((Integer)options.getOptionValues("-v")[0]).intValue();
	  }
	if (options.isOptionSet("--startup")) 
	  {
		CI.altStartFile = ((String)options.getOptionValues("--startup")[0]);
	  }
	if (options.isOptionSet("-s")) 
	  {
		CI.doBayesianAutoSmoothing = true;
	  }
	if (options.isOptionSet("-lt")) 
	  {
		CI.prefitLinearTrend = true;
	  }
	if (options.isOptionSet("-c")) 
	  {
		CI.doNoiseCooling= true;
		CI.doBayesianAutoSmoothing = false;
	  }
	if (options.isOptionSet("-n")) 
	  {
		CI.doNoiseEstimation = true;
		CI.doNoiseCooling = false;
	  }
	if (options.isOptionSet("-r")) 
	  {
		CI.doRelaxedSmoothing = true;
		CI.doNoiseCooling = false;
		CI.doBayesianAutoSmoothing = false;
	  }		
	if (options.isOptionSet("--MS")) 
	  {
		CI.doMultistart= true;
		CI.setMultistartStrategies(((String)options.getOptionValues("--MS")[0]));
	  }
	if (options.isOptionSet("--auto-split"))
	  {
		CI.doModelSelectAutoSplitting = true;
		CI.doNoiseCooling = false;
		CI.doBayesianAutoSmoothing =false;
		CI.doRelaxedSmoothing =false;
	  }
	if (options.isOptionSet("--split-suggest"))
	  {
		_EH.warn("--split-suggest turns on --auto-split");
		CI.doModelSelectAutoSplitting = true;
		CI.doNoiseCooling = false;
		CI.doBayesianAutoSmoothing =false;
		CI.doRelaxedSmoothing =false;
		CI.splitSuggest=true;
	  }
	if (options.isOptionSet("--Xsections")) 
	  {
		CI.doXsections = true;
	  }
	if (options.isOptionSet("-P")) 
	  {
		CI.dumpPrfl = true;
	  }
	if (options.isOptionSet("-q")) 
	  {
		CI.quietIO = true;
	  }
	if (options.isOptionSet("-N")) 
	  {
		CI.doMCMC = true;
		CI.N_MCMC = ((Integer)options.getOptionValues("-N")[0]).intValue();
		if (CI.N_MCMC < 0) _EH.die_err("Require -N number>=0");
	  }
	if (options.isOptionSet("-oSamp"))
	  {
		CI.samplesFilename = ((String)options.getOptionValues("-oSamp")[0]);
	  }
	if (options.isOptionSet("--MD"))
	  {
		CI.doMDstyleMCMC = true;
	  }
	if (options.isOptionSet("--NB"))
	  {
		CI.useBoundConstraints = false;
	  }
	if (options.isOptionSet("--bounds"))
	  {
		CI.option_bounds_lower = ((Float)options.getOptionValues("--bounds")[0]).floatValue();
		CI.option_bounds_upper = ((Float)options.getOptionValues("--bounds")[1]).floatValue();
		if (CI.option_bounds_lower > CI.option_bounds_upper) _EH.die_err("lower bound must be less than upper bound");
		if (CI.option_bounds_lower <-10.0f ) _EH.die_err("lower bound must be > -10.0");
		if (CI.option_bounds_upper > 10.0f ) _EH.die_err("upper bound must be < 10.0");
	  }

	if (options.isOptionSet("--no-smoothing-weights")) 
	  {
		CI.zeroSmoothingWeights = true;
		if (CI.doBayesianAutoSmoothing)
		  {
			CI.EH.println("Makes no sense to zero smoothing weights in Bayesian smoothing mode. Switch one off!");
			System.exit(0);
		  }
	  }
	if (options.isOptionSet("--conv")) 
	  {
		CI.dumpConvDiagnostics = true;
	  }
	if (options.isOptionSet("--TH-approx-test")) 
	  {
		CI.THApproxTest = true;
	  }
	if (options.isOptionSet("--no-TH-approx")) 
	  {
		CI.approximateTxHeights = false;
	  }
	if (options.isOptionSet("-e")) 
	  {
		CI.estDataErrors = true;
        if (((String)options.getOptionValues("-e")[0]).equals("g")) 
		  {
			CI.graphDataErrors = true;
		  }
	  }
	if (options.isOptionSet("--make-synthetic"))
	  {
		CI.makeSynthetic = true;
		CI.syntheticSpecFrac = ((String)options.getOptionValues("--make-synthetic")[0]);
		CI.syntheticSpecFloor = ((String)options.getOptionValues("--make-synthetic")[1]);
		CI.parseSynthNoiseSpecs();

		// these other options must be disabled
		CI.doMultistart = false;
		CI.doModelSelectAutoSplitting =false; 
		CI.doBayesianAutoSmoothing = false;
		CI.doNoiseCooling= false;
	  }

	if (options.isOptionSet("--noise-levels"))
	  {
		CI.syntheticSpecFrac = ((String)options.getOptionValues("--noise-levels")[0]);
		CI.syntheticSpecFloor = ((String)options.getOptionValues("--noise-levels")[1]);
		CI.parseSynthNoiseSpecs();
	  }

	if (options.isOptionSet("--log-grid"))
	  {
		CI.overwriteLogModel = true;
		s = ((String)options.getOptionValues("--log-grid")[0]);
		opts = s.split(",");
		try {CI.logModelFirstThickness = Float.parseFloat(opts[0]);}
		catch (NumberFormatException e) {CI.EH.die_err("--log-grid thickness[,ratio] has thickness not a number");}
		if (CI.logModelFirstThickness<=0.0) CI.EH.die_err("--log-grid thickness[,ratio] must have thickness > 0 ");
		if (opts.length > 1) 
		  {
			try {CI.logModelRatio = Float.parseFloat(opts[1]);}
			catch (NumberFormatException e) {CI.EH.die_err("--log-grid thickness[,ratio] has ratio not a number");}
			if (CI.logModelRatio < 1.0) CI.EH.die_err("--log-grid thickness[,ratio] must have ratio >= 1");
		  }
	  }
	if (options.isOptionSet("--read-profile"))
	  {
		CI.readProfile = true;
		CI.readProfileFilename = ((String)options.getOptionValues("--read-profile")[0]);		
	  }
	if (options.isOptionSet("--read-profile-set-mean"))
	  {
		CI.readProfile = true;
		CI.readProfileFilename = ((String)options.getOptionValues("--read-profile-set-mean")[0]);		
		CI.readProfileSetMean = true;
	  }

	if (options.isOptionSet("--check-new-data"))
	  {
		CI.checkMLmodelOnNewData = true;
		CI.checkNewDataFile = (String)(options.getOptionValues("--check-new-data")[0]);
	  }
	if (options.isOptionSet("--SVD"))
	  {
		CI.useSVD=true;
		CI.eps_SVD = ((Float)options.getOptionValues("--SVD")[0]).floatValue();
	  }
	if (options.isOptionSet("--s-trust"))
	  {
		CI.muTrustSize = ((Float)options.getOptionValues("--s-trust")[0]).floatValue();
		if (CI.muTrustSize <0.01) CI.EH.die_err("--s-trust must supply size > 0.01");
	  }
	if (options.isOptionSet("--LM"))
	  {
		CI.useLMstrategy = true;
	  }
	if (options.isOptionSet("--PI"))
	  {
		CI.doPathIntegralTests = true;
	  }
	if (options.isOptionSet("--BPB"))
	  {
		CI.resampleMeansForBayesianBootstrapping = true;
		if (CI.doRelaxedSmoothing)
		  {
			CI.EH.die_err("Bootstrapping mean-resampling not suitable (yet) with relaxed (Newman) mode");
		  }
		if (CI.doNoiseCooling)
		  {
			CI.EH.die_err("Bootstrapping mean-resampling not suitable (yet) with noise-cooling mode");
		  }
		if (CI.doBayesianAutoSmoothing)
		  {
			CI.EH.die_err("Bootstrapping mean-resampling not suitable (yet) with Bayesian smoothing");
		  }
		if (CI.doModelSelectAutoSplitting)
		  {
			CI.EH.die_err("Bootstrapping mean-resampling not suitable for Bayesian splitting ");
		  }
	  }
	if (options.isOptionSet("--PB"))
	  {
		CI.doBootstrapping = true;
		CI.N_PB = ((Integer)options.getOptionValues("--PB")[0]).intValue();
		if (CI.N_PB < 0) _EH.die_err("Require --PB number>=0");

		if (CI.nSpec==null)
		  {
			CI.syntheticSpecFrac = "5.0,5,5.0,5";
			CI.syntheticSpecFloor = "1e-15,1e-18";
			_EH.warn("FIXME: Warning: using default noise specs: "+CI.syntheticSpecFrac+" "+CI.syntheticSpecFloor);
			CI.parseSynthNoiseSpecs();
		  }

		CI.checkMLmodelOnNewData = false; // not needed
		if (CI.doRelaxedSmoothing)
		  {
			CI.EH.die_err("Bootstrapping not suitable (yet) with relaxed (Newman) mode");
		  }
		if (CI.doNoiseCooling)
		  {
			CI.EH.die_err("Bootstrapping not suitable (yet) with noise-cooling mode");
		  }
		if (CI.doBayesianAutoSmoothing)
		  {
			CI.EH.die_err("Bootstrapping not suitable (yet) with Bayesian smoothing");
		  }
		if (CI.doModelSelectAutoSplitting)
		  {
			CI.EH.die_err("Bootstrapping not suitable for Bayesian splitting ");
		  }
		if (CI.doPathIntegralTests) CI.EH.die_err("Bootstrapping and Path-Integral tests can't be done together");
	  }
	if (options.isOptionSet("--PB-recenter"))
	  {
		CI.doBootstrappingRecenter = true;
	  }

	if (options.isOptionSet("--model-study")) 
	  {
	    CI.doModelStudies = true;
		CI.modelStudy =  ((String)options.getOptionValues("--model-study")[0]);	
		CI.doBootstrapping=false;	
		CI.doModelSelectAutoSplitting =false; 
	  }
	if (options.isOptionSet("--systematics")) 
	  {
	    CI.modelSystematics = true;
		CI.systematicsFlavour =  ((String)options.getOptionValues("--systematics")[0]);	
	  }
	if (options.isOptionSet("--dilute")) 
	  {
	    CI.dilute = true;
		CI.Ndilute = ((Integer)options.getOptionValues("--dilute")[0]).intValue();
		if (CI.Ndilute<1) CI.EH.die_err("--dilute N must have N>=1");
	  }	
	if (options.isOptionSet("--bootstrap")) 
	  {
	    CI.doResamplingBootstrap = true;
		CI.NBootstrap = options.getOptionInt("--bootstrap",0);
		if (CI.NBootstrap<1) CI.EH.die_err("--bootstrap N must have N>=1");
	  }
	if (options.isOptionSet("--Wishart")) 
	  {
		CI.doWishart = true;
		CI.Wishart_nu = options.getOptionFloat("--Wishart",0);
		if (CI.Wishart_nu <= 0) _EH.die_err("Require --Wishart nu >0");
	  }
	if (options.isOptionSet("--shrinkage")) 
	  {
		CI.doShrinkage = true;
		CI.shrinkage_alpha = options.getOptionFloat("--shrinkage",0);
		if (CI.shrinkage_alpha <= 0 || CI.shrinkage_alpha>=1) _EH.die_err("Require --Shrinkage alpha, 0<alpha<1");
	  }
	if (options.isOptionSet("--NS"))
	  {
		CI.newScripps = true;
	  }
	if (options.isOptionSet("--dump-debug-package")) CI.dumpDebug = true;
	if (options.isOptionSet("--default-prior-constants"))
	  {
		CI.PTI.overrideDefaults((String)options.getOptionValues("--default-prior-constants")[0]);	
	  }

	if (options.isOptionSet("-X") || options.isOptionSet("-h")) 
	  {
		options.Usage();
	  }
	
	if (CI.dumpDebug)
	  {
		CI.dumpDebugPackage(argv);
		System.exit(0);
	  }

	CI.EH.printf("java "+CI.classname+" ");
	for (int i=0;i<argv.length;i++) CI.EH.printf(argv[i]+" ");
	CI.EH.println("");

	CI.setComputationalKernel();
	CI.setSystematicsFlavour();
	CI.startTimer();

	// CI.testSplittingAlg();

	if (CI.doBootstrapping)
	  {
		CI.performPBsampling();
	  }
	else if (CI.doModelStudies)
	  {
		CI.performModelStudies();
	  }
	else if (CI.doResamplingBootstrap)
	  {
		CI.performResamplingBootstrapStudy();
	  }
	else if (CI.doModelSelectAutoSplitting)
	  {
		CI.runWithResolutionModelSelect(true);
	  }
	else
	  {
		CI.run(true);
	  }

	CI.showTiming();
  }

  public CSEM1DInverter()
  {
	// handle errors: can change how these are handled if desired
	EH = new au.csiro.JamesGunning.lib.ErrorHandler(); 
	classname = ((new Throwable()).getStackTrace())[0].getClassName();
  }

  void setComputationalKernel()
  {
	if (newScripps)
	  {
		SMI = new ScrippsModelInterface(this);
	  }
  }

  void setMultistartStrategies(String multistartOpts)
  {
	String[] strategies;  
	
	
	strategies = multistartOpts.split(":");
	for (int i=0;i<strategies.length;i++)
	  {
		if (strategies[i].equals("d")) 
		  {
			multistartFlipping = true;
		  }
		else if (strategies[i].startsWith("r,"))
		  {
			multistartRandom = true;
			try
			  {
				multistartRandomN = Integer.parseInt(strategies[i].substring(2));
			  }
			catch  (NumberFormatException e)
			  {
				EH.die_err("Number of random starts suggested in multistart option "
						   +multistartOpts+" not a positive integer");
			  }
			if ((multistartRandomN < 1) || (multistartRandomN > 10000))
			  {
				EH.die_err("Number of random starts suggested in multistart option "
						   +multistartOpts+" not in [1,10000]");
			  }
		  }
		else
		  {
			EH.die_err("Uninterpreted multistart style "+strategies[i]+" in "+multistartOpts);
		  }
	  }
  }

  void startTimer()
  {
	tOverall = System.currentTimeMillis();
  }

  void showTiming()
  {
	tOverall = System.currentTimeMillis() - tOverall; 
	if (verbosity > 1) EH.printf("Time in native methods=%6.2f s, %5.2f%% of total run-times%n",
								 (1e-9*tNatives),100*(1e-9*tNatives)/(1e-3*tOverall));  
  }
 
  void readModelAndInputs()
  {
	// Read & check against startup files...
	read_startup_and_model();

	// input the data
	if (verbosity > 0) EH.println("Reading data from file: "+datfil);
	inputData();
	if (estDataErrors) medianFilterReEstimateDataErrors();
	if (dilute) diluteErrorBars();

	if (overwriteLogModel && !readProfile)
	  {
		overwriteLogGridModel();
	  }
	if (readProfile && !newScripps)
	  {
		overwriteGridModelFromProfile(readProfileFilename);
	  }

	clean_smoothing(); 
	// Summarise runmode
	dumpRunmode();

	// make transmitter height approximations
	lumpTransmitterHeights();

  }

  void setUpSystematicModelPieces()
  {
	int i,j,k,ii,iS=0,iE=0,n=0,nB,type;
	int Nrows,Ncols;
	float f, fOld=-1f,range,dr,scale;
	DoubleMatrix2D block=null;
	SystematicGroup s=null;

	SEM.systematicGroups = new ArrayList<SystematicGroup>();

	// n = data length
	if (!newScripps)
	  {
		n = d.length;
		for (i=0;i<n;i++)
		  {		
			// dp[2] = freq, 
			// dp[1] = range,
			// dp[0] = type.
			// dp[5] = trans. height.
			f = dp.getQuick(2,i);
			range = dp.getQuick(1,i);
			type = (int)dp.getQuick(0,i);
			if (type==1)
			  {
				if (s==null || f!=fOld)
				  {
					fOld=f;
					s = new SystematicGroup();
					s.f = f;
					s.indexes.add(i);
					s.ranges.add(range);
					s.rangemin = s.rangemax = range;
					SEM.systematicGroups.add(s); // add to complete list
				  }
				else
				  {
					if (f==fOld)
					  {
						s.indexes.add(i);
						s.ranges.add(range);
						s.rangemin = min(s.rangemin ,range);
						s.rangemax = max(s.rangemax ,range);
					  }
				  }
			  }
		  }
	  }
	else
	  {
		// New Scripps format
		int iFreq,iFreqOld=0,typeOld=0,iTx,iTxOld=0,iRx;
		DenseIntMatrix2D _dp = SMI.OM.dp_trans; 
		ScrippsModelInterface.Csem1d_mod _CM = SMI.CM;

		// sort data into groups of common transmitter #, freq, type.
		n = nd;
		int[][] tA = new int[nd][5];		
		for (i=0;i<nd;i++)
		  {
			tA[i][0] = _dp.getQuick(0,i); // type
			tA[i][1] = _dp.getQuick(1,i); // freq
			tA[i][2] = _dp.getQuick(2,i); // transmitter #
			tA[i][3] = _dp.getQuick(3,i); // receiver #
			tA[i][4] = i;                 // original index		  
		  }
		Arrays.sort(tA, new NSDataComparator());

		// data now sorted by transmitter #, freq, type. Chop into groups as above
		// Group must have common freq, transmitter number, type
		for (i=0;i<nd;i++)
		  {		
			type = tA[i][0];
			iFreq = tA[i][1] - 1; // zero-offset indices - freq
			iTx = tA[i][2] - 1;   // zero-offset indices - trans
			iRx = tA[i][3] - 1;   // zero-offset indices - rec
			j = tA[i][4];         // original index
			if (type==21 || type==23 || type==25 ||  // |Ex|, |Ey|, |Ez|
				type==31 || type==33 || type==35)    // |Bx|, |By|, |Bz|
			  {
				// range -> use euclidean TX-Rx distance (got a better idea?)
				range =	(float)Math.hypot(_CM.xRxIn[iRx] - _CM.xTxIn[iTx], _CM.yRxIn[iRx] - _CM.yTxIn[iTx]);

				if (s==null || (iFreq != iFreqOld) || (type != typeOld) || (iTxOld != iTx))
				  {
					iFreqOld = iFreq;
					typeOld = type;
					iTxOld = iTx;
					s = new SystematicGroup();
					s.f = (float)_CM.fTxIn[iFreq];
					s.indexes.add(j); // use original index
					s.ranges.add(range);
					s.rangemin = s.rangemax = range;
					SEM.systematicGroups.add(s); // add to complete list
				  }
				else
				  {
					if ((iFreq==iFreqOld)  && (type == typeOld) && (iTxOld == iTx))
					  {
						s.indexes.add(j);
						s.ranges.add(range);
						s.rangemin = min(s.rangemin ,range);
						s.rangemax = max(s.rangemax ,range);
					  }
				  }
			  }
		  }
		tA = null; // GC
	  }

	// code hereon, independent of forward kernel

	SEM.systematicGroups.trimToSize();
	SEM.nPsystematics = 0;
	// Assemble systematics model-vector effective prior
	Ncols = SEM.systematicsRegressionOrder + 1;  // since we need terms 0,1..systematicsRegressionOrder
	Nrows = Ncols;
	// spacing of "test" points"
	dr = (Ncols>1) ? (1.0f/(Ncols-1.0f)) : 0.0f;
	for (i=0; i<SEM.systematicGroups.size(); i++)
	  {
		s = SEM.systematicGroups.get(i);
		s.indexes.trimToSize(); 
		s.X_NT = new DenseDoubleMatrix2D(Nrows,Ncols);
		// fill s.X_NT_ij =  u_{i,test}^{j-1}/sigma_sys
		for (j=0;j<Nrows;j++)
		  {	   
			// range_j = equispaced test points along offset span. Used normalised variable in [0,1], so
			// "distance" corresponds to an actual unitised length:
			// range = (j/ncols)*(s.rangemax-s.rangemin) + s.rangemin ; 
			// typically, e.g. test points are 0, 0.5, 1.0
			range = j * dr;	
			for (k=0;k<Ncols;k++)
			  {
				// regression is power series in range: Block parameter m_{N,k} is coeff of kth-power
				// Absorb also "test-point" prior std-deviation
				s.X_NT.setQuick(j,k,pow(range,k)/SEM.sigmaSystematicTestPoint);
			  }
		  }
		SEM.nPsystematics += Ncols;
	  }

	// Form C_BsysInv and its determinant

	SEM.C_BsysInv = new DenseDoubleMatrix2D(SEM.nPsystematics,SEM.nPsystematics);

	
	cern.colt.matrix.linalg.Algebra alg = new cern.colt.matrix.linalg.Algebra();
	SEM.logDetInvC_Bsys =0.0;
	for (i=0; i<SEM.systematicGroups.size(); i++)
	  {
		s = SEM.systematicGroups.get(i);
		// point at sub-block on diagonal
		block = SEM.C_BsysInv.viewPart(i*Nrows,i*Ncols,Nrows,Ncols);
		// fill block with X_NT^T.X_NT
		SeqBlas.seqBlas.dgemm(true,false,1.0,s.X_NT,s.X_NT,0,block);
		// log determinant is sum over blocks
		SEM.logDetInvC_Bsys += Math.log(alg.det(block));
	  }
	if (resampleMeansForBayesianBootstrapping || doBootstrapping)
	  {
		// Form The cholesky L factor of the inverse while we're here, for bootstrap
		// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("C.txt",SEM.C_BsysInv);
		SEM.L_Bsys = (DenseDoubleMatrix2D)(new CholeskyDecomposition(alg.inverse(SEM.C_BsysInv))).getL();
		// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("L.txt",SEM.L_Bsys);
	  }

	// Form also effective full design matrix for all parameters
	// Forward model correction is of form Cd^{1/2}.{polynomial-in-normalised-offset}
	// All of which is absorbed in SEM.X

	SEM.X = new DenseDoubleMatrix2D(n,SEM.nPsystematics);
	for (i=0; i<SEM.systematicGroups.size(); i++)
	  {
		// point at block for this frequency's parameters
		block = SEM.X.viewPart(0,i*Ncols,n,Ncols);
		s = SEM.systematicGroups.get(i);
		nB = s.indexes.size();
		scale = s.rangemax - s.rangemin;
		if (scale==0.0) scale=1.0f;
		for (j=0;j<nB;j++)
		  {
			ii = s.indexes.get(j); // index of data point in original data
			range = (s.ranges.get(j) - s.rangemin)/scale; 
			for (k=0;k<Ncols;k++)
			  {
				block.setQuick(ii,k,sd[ii] * pow(range,k));
			  }	
		  }
	  }
	// EH.println("Xfull: "+SEM.X);
	// all model parameters
	SEM.m_N = new double[SEM.nPsystematics];
	// extra systematics coefficients: mean = 0
	SEM.m_N_bar = new double[SEM.nPsystematics];
  }

  void run(boolean readNewData)
  {
	int i,j,k,kk,n_cool,itmax,multistart_cntr=0,k_multistart=0,multistart_phase=0;
	float r,m_tmp,tobt=0f,tobt_full=0,tobt_full_old=0,chisq_prej=0,tmp,noise_stddev=1.0f;
	boolean write_final_it=false,keep_iterating,seekNewMultistartSolution;
	boolean termTest;
	cern.colt.list.FloatArrayList convergenceResidualArray = null;
	ModeCovarianceAndStdDev modeCov;
	MultiStartStrategyManager MSSM=null;
	float gr,dist;

	float INITIAL_COOLING_MULTIPLIER = 50;
	double GUESS_COOLING_RATIO = 0.8;
	/*
	  threshold criterion to catch over-relaxed Newman-style style: if
	  chisq-data < (1-RELAXED_SMOOTHING_CAPTURE_EPS)* rms_data_misfit ( = required chisq_data/N)
	*/
	float RELAXED_SMOOTHING_CAPTURE_EPS = 0.05f;
 
	// stuff for relaxed smoothing
	float GUESS_RELAX_LOGMU_DECREMENT = 0.2f;
	float PMU_RELAX_MIN = 0.5f;  // terminate here if chisq never reaches tolreq
	float PMU_RELAX_EPS = 0.01f;  // terminate in relaxation if delta_pmu < 0.01
	float pmu_last=0f,tobt_last=0f;
	float[] pm_last=null;

	double cooling_ratio=GUESS_COOLING_RATIO;


	EH.println("Starting CSEM1DInverter");

	if (readNewData)
	  {
		readModelAndInputs();
		
		if (modelSystematics)
		  {
			setUpSystematicModelPieces();
			allocWorkspace(np);
		  }

		initialise();
		conformStartupToBounds();
		// Set JNI pointers to fortran module data
		setNativeModelPointers();	
	  }

	if (doModelStudies)
	  {
		releaseNativeModelPointers();
		overWriteWithModelStudyValues();
		setNativeModelPointers();	
	  }

	if (doBootstrapping)
	  {
		// must do runs from clean start
		nit = 0;
	  }

	itmax = maxitr + nit;
	if (nit == 0)
	  {
		rlast = Float.MAX_VALUE;
		tobt = Float.MAX_VALUE;
		konv=0;
		tobt_full_old = Float.MAX_VALUE;
		write_final_it = false;
		multistart_cntr = 1;
		k_multistart = 0;
		itrout(null,maxitr,rlast,tobt,0,multistart_cntr);
	  }

	if (resampleMeansForBayesianBootstrapping)
	  {
		resampleMean();
	  }

    EH.println("Starting major loop...");

	optima = new ArrayList();
	seekNewMultistartSolution = true;


	if (dumpConvDiagnostics) convergenceResiduals = new ArrayList();

	if (doNoiseCooling)
	  {
		global_error_scale = INITIAL_COOLING_MULTIPLIER;
		n_cool = (int) floor(-Math.log(global_error_scale)/Math.log(GUESS_COOLING_RATIO));
		cooling_ratio = Math.exp(-Math.log(global_error_scale)/n_cool);
		if (verbosity > 1) EH.println("Cooling schedule: Initial factor= %g, N=%d, ratio=%g%n",
									  global_error_scale,n_cool,cooling_ratio);
	  }
	if (doNoiseEstimation)
	  {
		global_error_scale = 1.0f;
		if (verbosity > 1) EH.println("Performing noise level estimation");
	  }

	if (doRelaxedSmoothing)
	  {
		pmu = PMU_DEFAULT_START + GUESS_RELAX_LOGMU_DECREMENT; // overwritten by eigenvalue guess anyway
		pmu_last = pmu;
		pm_last = (float[])pm.clone();
	  }

	if (doMultistart) 
	  {
		MSSM = new MultiStartStrategyManager(multistartFlipping,multistartRandom);
		if (!doNoiseCooling && !doRelaxedSmoothing) MSSM.initialiseStrategies();
	  }

	// set premod by approximately best linear trend
	if (prefitLinearTrend) findOptimalLinearBackground();


	nit_ms = 0;
	while (seekNewMultistartSolution)
	  { 
		if (dumpConvDiagnostics) 
		  {
			convergenceResidualArray = new cern.colt.list.FloatArrayList();
			convergenceResiduals.add(convergenceResidualArray);
		  }
		if (doMultistart)
		  {
			// in this mode, these numbers need to be clean
			rlast = Float.MAX_VALUE;
			nit = 0;
			itmax = maxitr;
		  }
		if (doNoiseCooling)
		  {
			global_error_scale *= cooling_ratio;
			rlast = Float.MAX_VALUE;
			nit = 0;
			itmax = maxitr;
			if (verbosity > 1) EH.println("Cooling: new noise multiplier %g%n",global_error_scale);
		  }
		if (doRelaxedSmoothing)
		  {
			rlast = Float.MAX_VALUE;
			nit = 0;
			itmax = maxitr;
			if (nit_ms==0) pmu_last = pmu;
			pmu -= GUESS_RELAX_LOGMU_DECREMENT;
			if ((verbosity > 1) && (nit_ms>0)) EH.printf("Relaxing: new smoothing parameter pmu=%4f%n",pmu);
		  }
		// start iterative loop for local optimization
		keep_iterating = true;
		while (keep_iterating)
		  {
			r = rlast;
			tobt = GaussNewtonUpdate(tobt,nit,doBayesianAutoSmoothing,multistart_cntr);
			// tobt = raw misfit
			nit++;

			
			// At this point, always,
			// ptp = partial^T.partial + W_B^T.W_B
			// wjwjt = (WJ)^T WJ

			// write an iteration file
			itrout(null,maxitr,rlast,tobt,0,multistart_cntr);
			// let the console know roughly what is going on
			if (verbosity > 1) 
			  {
  				EH.println("Iteration %d, Math.log(mu) = %10.6g, RMS misfit = %10.6g",nit,pmu,doNoiseEstimation ? global_error_scale : tobt);
			  }
			// test for convergence: continue inversion if 
			// 1) iteration number does not exceed maximum, and
			// 2) stepsize is still significant, and
			// 3) rms misfit larger than required, and 
			// 4) there are no irrecoverable convergence problems

			// JG...
			tobt_full = nd * tobt*tobt;
			noise_stddev = tobt;
			if (doNoiseEstimation)
			  {
				// tobt_full = (float)(nd,np) * (1.0+Math.log(tobt_full/nd))); // incorporates optimization on noise-scalar
				// EB version of above.
				tobt_full *= sqr(global_error_scale); // restores misfit wrt original data-file noise levels
				float dof = DOF(nd,np);
				noise_stddev = (float)sqrt(tobt_full/dof);
	  	        tobt_full = dof + nd*(float)Math.log(sqr(noise_stddev));   
			  }
			chisq_prej = chisqPrej(pmu,pm);
			tobt_full += fpow10(pmu) * chisqRoughness(pm) + chisq_prej;
			if (doBayesianAutoSmoothing && !doEmpBayes) tobt_full += negLogDetCpEffInv(pmu,true);

			// compute  norm of step
			float step_size = (float)sqrt(anorm(pwk3));

 			if (verbosity > 2) 
			  {
				EH.println(" -2log_posterior_mismatch= %-15.8g    -2log_posterior_full= %-15.8g",
						   (nd * tobt*tobt),tobt_full);
				// typical convergence termination criterion
				EH.println(" objective function decrease: "+(tobt_full-tobt_full_old));
				EH.printf("|model| relative_step_size= %-15.8g%n",step_size/SQRT_EPS_CHISQ);
				EH.printf("|grad|= %-15.8g |grad|_test= %-15.8g%n",(float)sqrt(anorm(grad)),CUBE_ROOT_EPS_CHISQ*(tobt_full+1.0));
			  }
			if (dumpConvDiagnostics) convergenceResidualArray.add(tobt_full);

			termTest = getTerminationCriteria(tobt_full,tobt_full_old,step_size);


			if (doBayesianAutoSmoothing && avoidTerminationBayesianAutoSmoothing)
			  {
				// do nothing - not allowed to terminate on a step with mu  frozen
			  }
			else if ((nit >= 1) && termTest)
			  {
				if (verbosity > 1) 
				  {
					EH.println("Satisfactory full convergence stop");
					showParams(pwk2,true); 
					if (verbosity > 2)
					  {
						EH.println("chisq_misfit= %-15.8g chisq_prior= %-15.8g",nd * tobt*tobt,chisq_prej);
						EH.println("Noise-adjusted chisq_tot (for no smoothing) = %-15.8g",(nd-np)+chisq_prej);
					  }
				  }
				H = computeFullHessian(pmu);
				write_final_it = true;
				keep_iterating = false;
			  }
			else if (doRelaxedSmoothing &&  // opportunities for early escape in undersmoothed Newman models
					 (nit >= 1) && (tobt <= tolreq*(1-RELAXED_SMOOTHING_CAPTURE_EPS)))
			  {
				if (verbosity > 2) EH.println("Premature optimisation break, undersmoothed Newman-mode");
				write_final_it = true;
				H = computeFullHessian(pmu);
				keep_iterating = false;
			  }

			tobt_full_old = tobt_full;
		

			if (keep_iterating) 
			  {
				if (konv != 0)
				  {
					// out_logger.println(" convergence problems of type "+konv);
					if (verbosity > 1) EH.println(" convergence problems of type "+ konv);
					keep_iterating = false;
				  }
				else if (nit >= itmax) 
				  {
					// out_logger.println(" stop on maximum iterations");
					if (verbosity > 1) EH.println(" stop on maximum iterations");
					write_final_it = true;
					H = computeFullHessian(pmu);
					keep_iterating = false;
				  } 
			  }
		  } // major local-optimisation loop: 

		if (write_final_it) 
		  {
			// write final iteration file
			itrout(null,maxitr,rlast,noise_stddev,2,multistart_cntr);


			// uncertainty files
			// use this to show uncertainties without smoothing (not true local minimum,
			// but estimates uncertainty scales if smoothing were not present)
			makjtj();
			float pmu_tmp = -10f;
			makptp(pmu_tmp);
			H = computeFullHessian(pmu_tmp);
			modeCov = computeUncertainty();
			writeUncertaintyFiles(null,true,false,multistart_cntr,modeCov);
			
			// To show cross sections through local optimum, 
			if (doXsections) Xsections(pmu,false,multistart_cntr);
 
			// collect optima in list
			if (verbosity > 1) EH.println("Multistart solution "+multistart_cntr+" obj="+tobt_full);
			if (doNoiseCooling || doRelaxedSmoothing)
			  {
				// overwrite mode
				LocalMode mode = new LocalMode(np,current_resp_filename,current_profile_filename,current_iter_filename,
											   pm,pmu,tobt_full/2.0f, chisq_prej, noise_stddev, 
											   // wjtwj,modeCov.stdv);
											   H, modeCov.stdv);
				optima.clear(); 
				optima.add(mode);
			  }
			else if (!foundSimilarExistingMode(tobt_full/2.0f))
			  {			
				LocalMode mode = new LocalMode(np,current_resp_filename,current_profile_filename,current_iter_filename,
											   pm,pmu,tobt_full/2.0f, chisq_prej, noise_stddev ,
											   // wjtwj,modeCov.stdv);
											   H, modeCov.stdv);
				optima.add(mode);
			  }

			if (!doMultistart && !doNoiseCooling && !doRelaxedSmoothing)
			  { 
				// single pass only. Bail out ASAP 
				seekNewMultistartSolution = false;
			  }
			else if (doNoiseCooling)
			  {
				itrout(null,maxitr,rlast, noise_stddev,2,multistart_cntr);
				if (global_error_scale > 1.01)
				  {
					seekNewMultistartSolution = true;
				  }
				else
				  {
					seekNewMultistartSolution = false;	
				  }
				multistart_cntr++;
			  }
			else if (doRelaxedSmoothing)
			  {
				itrout(null,maxitr,rlast,noise_stddev,2,multistart_cntr);
				if (verbosity > 1) EH.println("Relaxed solution "+multistart_cntr+" pmu= "+pmu+" obj="+tobt_full+
											  " rms="+noise_stddev+" tot= "+tobt_full);

				// terminate if 1 or 2 or 3 or 4
				if ((pmu < PMU_RELAX_MIN) ||     // (1) over-relaxed. No point relaxaing any more
					((tobt < tolreq) && (tobt > tolreq*(1.0f-RELAXED_SMOOTHING_CAPTURE_EPS))) ||// (2) answer in "capture window"
					(abs(pmu_last - pmu) < PMU_RELAX_EPS) ||      // (3) or pseudo interpolation on pmu contracted to within sufficient PMU_RELAX_EPS
					((tobt < tolreq) && (multistart_cntr==1))) // (4) first optimization achieved OK solution
				  {
					// termination at required chisq
					seekNewMultistartSolution = false;
				  }
				else 
				  {
					if ((tobt < tolreq) && (multistart_cntr>1))
					  {
						// went too far
						pmu = pmu + (pmu_last - pmu)*(tolreq - tobt)/(tobt_last - tobt);
						pmu += GUESS_RELAX_LOGMU_DECREMENT; // anticipate factor on loop return
						System.arraycopy(pm_last,0,pm,0,np); // overwrite last decent oversmoothed solution
						if (verbosity > 1) EH.println("Relaxed solution: resetting...");
					  }
					if (tobt >= tolreq)					  
					  {
						tobt_last = tobt;
						pmu_last = pmu;
						System.arraycopy(pm,0,pm_last,0,np);	// pm_last <- pm   
					  }
					seekNewMultistartSolution = true;
					multistart_cntr++;
				  }
			  }
			else if (doMultistart && (optima.size() > 0))
			  {
				// Multistart restarts only sensible if at least one optima found so far.
				// set new start point
				if (MSSM.tryPhaseFlipping && MSSM.phaseFlipping)
				  {
					// change start points suitably, and recommence
					if (verbosity > 1) EH.println("Multistart solution "+multistart_cntr+" obj="+tobt_full);
					itrout(null,maxitr,rlast,tobt,2,multistart_cntr); 
					// Restart value based on first optimum..
					LocalMode mode0 = optima.get(0);
					// search forward for break in profile likely to be ambiguous...
					seekNewMultistartSolution = false;
					for (k = MSSM.k_flips ; k< (np-1); k++)
					  {
						dist = abs(mode0.pmLoc[k] - mode0.pmLoc[k+1]);
						if (mode0.stdv[k] > 0) dist/= mode0.stdv[k];
						// large dist value indicated possible "ambiguous" dislocation in solution
						// between indices k,k+1. Start new optimization based on flipped pair
						if (dist > MULTISTART_DIST_THRESH)
						  {
							System.arraycopy(mode0.pmLoc,0,pm,0,np);	// pm <- pmLoc   
							// flip layers & fire off new optimization
							m_tmp = pm[k];
							pm[k] = pm[k+1];
							pm[k+1] = m_tmp;
							seekNewMultistartSolution = true;
							MSSM.k_flips  = k+1; // next search starts here
							if (verbosity > 1) EH.println("MultistartNewSolution "+multistart_cntr);
							break;
						  }	
					  }
					if ((k >= (np-1)) && !seekNewMultistartSolution) 
					  {
						MSSM.exitFlipping(); // opens up random...
					  }
				  }
				if (MSSM.tryRandom && MSSM.phaseRandom)
				  {
					// random starts
					if (MSSM.random_ctr < multistartRandomN)
					  {
						// next phase - randomness
						if ((optima.size() > 1) && MSSM.random_ctr==0) 
						  {
							if (verbosity > 1) EH.printf("Multistart: presort of local optima%n",optima.size());
							Collections.sort(optima);
						  }
						LocalMode mode0 = optima.get(0);
						for (kk=0;kk<np;kk++)
						  {
							pm[kk] = mode0.pmLoc[kk] + 2.0f*((float)Math.random() - 0.5f);								
						  }
						conformStartupToBounds();
						seekNewMultistartSolution = true;
						MSSM.random_ctr++;
						if (verbosity > 1) EH.println("MultistartNewSolution (random start) "+multistart_cntr);							
					  }
					else
					  {
						seekNewMultistartSolution = false;
						MSSM.exitRandom();
					  }
				  }
				
				if (seekNewMultistartSolution)
				  {
					// go back and try new local optimization
					multistart_cntr++;
				  }
				else
				  {
					// give up multistart. Take what we've got and push on
					seekNewMultistartSolution = false;
				  }
			  }
		  } // last optimization "successful"
		else
		  {
			// some kind of convergence problem. Bail out ASAP 
			seekNewMultistartSolution = false;			
		  }
		if (seekNewMultistartSolution) nit_ms++;
	  }  // multistart search


	cleanUpMultistartSolutions(null);

	if (dumpConvDiagnostics) 
	  {
		dumpConvergenceDiagnostics();
	  }
	
	if (THApproxTest)
	  {
		doTHApproxTest();
	  }

	if (checkMLmodelOnNewData)
	  {
		checkNewData();
		itrout(null,maxitr,rlast,noise_stddev,2,multistart_cntr+1);
	  }
	
	if (doMCMC) 
	  {
		GenerateMCMCPosterior();
	  }
	// out_logger.close();
  }




  void runWithResolutionModelSelect(boolean readNewData)
  {
	int i,j,k,kk,itmax,multistart_cntr=0;
	float r,m_tmp,tobt=0f,tobt_full=0,tobt_full_old=0,chisq_prej=0,noise_stddev=1f,dof;
	boolean write_final_it=false,keep_iterating,seekNewMultistartSolution,termTest;
	cern.colt.list.FloatArrayList convergenceResidualArray = null;
	ModeCovarianceAndStdDev modeCov;
	MultiStartStrategyManager MSSM=null;
	float dist;

	float LOGMML_SAFETY_WINDOW = 4.0f;

	System.out.println("Starting CSEM1DInverter");

	if (readNewData)
	  {
		readModelAndInputs();
	  }
	if (newScripps) SMI.copyFixedUpperLayers();// Copy upper fixed layers
					

	itmax = maxitr + nit;

	if (doMultistart) MSSM = new MultiStartStrategyManager(multistartFlipping,multistartRandom);

	// model selection bits...
	int i_grid,NlayerMax;
	boolean firstRun=true;
	float max_thickness,rhoSea,thickSea,pm_startval;
	float best_split_min,best_split_rms,best_MML;
	float zBML,dz,rho_mean,rho_sigma;
	int n_extra_safety_splits;
	int best_split_indx =0; 
	float[] pm_best_fine;
	float[] pm_best_last_level_fine;
	float[] thickf_best=null;
	int NresMax = LOG2_MAX_NUM_SPLIT_GRID_ELEMENTS;    // 2^NresMax possible subdivisions: probably 16 or 32 max
	int grid_size,grid_start_size,grid_max_splits;
	int[] grid;
	int[] old_grid;

	msInfo = new ModelSelectInfo();
	NlayerMax = 1<<NresMax;  // 2^NresMax 
	// initial total subsea model thickness
	max_thickness = getMaxThickness();
	rhoSea = rho[0]; // save these
	thickSea = thickf[0];
	pm_startval = 0.1f;// pm[0]; // top layer start resistivity
	pm_best_fine = new float[NlayerMax];
	pm_best_last_level_fine = new float[NlayerMax];

	// Loop over family of models. Start at 1-layer model
	grid_start_size = 1;
	// Don't attempt any more splits than this
	grid_max_splits = NlayerMax; // 2^NresMax
	n_extra_safety_splits=0;

	old_grid = new int[grid_start_size];
	for (i=0; i<grid_start_size; i++) old_grid[i] = NlayerMax/grid_start_size;

	best_MML = Float.MAX_VALUE;
	for (grid_size = grid_start_size; grid_size < (grid_start_size + grid_max_splits); grid_size++)
	  {
		grid = new int[old_grid.length+1];

		best_split_min = best_split_rms = Float.MAX_VALUE;
		// examine all possible splits of this grid
		for (i_grid=0; i_grid<old_grid.length; i_grid++)
		  {
			String extraText = ".level_"+grid_size+".split_"+i_grid; 
			if (old_grid[i_grid]==1) continue; // can't split this layer

			// model selection....impose constant agnostic model
			// Set JNI pointers to fortran module data
			if (!firstRun) // Release JNI pointers of last model before new allocs...
			  {
				// OLD ScrippsLegacyModelInterface.releaseModelPointersNative(thickf,rho,freqs,ipen1,ipen2,iff);
				releaseNativeModelPointers();
			  }
			np = grid.length;
			nlayer = np;
			allocateModel();
			allocateGeoModel();
			if (newScripps)
			  {
				// expect, by here, DM.nlay1D is set consistently (as per allocateModel())
				SMI.CM.allocateL(SMI.DM.nlay1D);
				SMI.DM.allocateL(SMI.DM.nlay1D);
				SMI.restoreFixedUpperLayers();
			  }
			allocWorkspace(np);
			// repopulate
			rho[0] = rhoSea;
			thickf[0] = thickSea;

			for (i=0,j=0; i<old_grid.length; i++,j++)
			  {
				if (i==i_grid) // split layer i_grid
				  {
					grid[j] = old_grid[i]/2;
					grid[++j] = old_grid[i]/2;
				  }
				else grid[j] = old_grid[i];
			  }

			if (newScripps) j = SMI.CM.Nfixed;
			for (i=1,dz=0; i<=nlayer; i++)
			  {
				rho[i] = 0.0f;
				if (firstRun) pm[i-1] = pm_startval;
				thickf[i] = grid[i-1]*max_thickness/NlayerMax;
				dz += thickf[i];
				zBML = dz + thickf[i]/2.0f;
				rho_mean = PTI.getDefaultPriorMean(zBML);
				rho_sigma = PTI.getDefaultPriorSigma(zBML);
				ipen1[i] = 0.0f;   // roughness penalty
				premod[i-1] = rho_mean; // prior mean
				ipen2[i] = 1.0f/rho_sigma;   // prior weight
				if (newScripps)
				  {
					// [top_depth, rho, penalty, Preference, prej_penalty] 
					SMI.CM.layer_params_trans.setQuick(0,j, SMI.CM.ZtopFree + (dz -thickf[i]) );
					SMI.CM.layer_params_trans.setQuick(1,j,-1.0); // free
					// no need to set the rest
					j++;
				  }
			  }
			if (newScripps)
			  {
				SMI.freeJacobianAllocs();
				SMI.finishModelSetup();
			  }



			for (i = 0; i<nlayer;i++) iff[i] = i+1;

			if (verbosity > 1)
			  {
				EH.println("Equivalent new model file contents:");
				EH.println("Thickness rho smooth-wt <log10(rho)>_prior prior_wt (prior-std-dev)");
				for (i=1; i<=nlayer; i++)
				  {
					EH.println("%10.6g -1 0       %10.6g         %10.6g    (%-10.6g)",
							   thickf[i],premod[i-1], ipen2[i],1/ipen2[i]);
				  }
			  }
			// map in new starts from best solution on previous grid
			if (!firstRun) 
			  {
				for (i=0,k=0; i<grid.length; i++)
				  {
					pm[i] = pm_best_last_level_fine[k];
					k+= grid[i];
				  }
			  }
			nit = 0;


			initialise();
			conformStartupToBounds();
			// Set JNI pointers to fortran module data
			setNativeModelPointers();
			
			firstRun=false;

			if (nit == 0)
			  { 
				pmu = 0.0f; // splitting mode should have no smoothing anyway (ipen1[i] = 0.0f = roughness penalty)
				rlast = Float.MAX_VALUE;
				tobt = Float.MAX_VALUE;
				konv=0;
				tobt_full_old = Float.MAX_VALUE;
				write_final_it = false;
				multistart_cntr = 1;
				itrout(extraText,maxitr,rlast,tobt,0,multistart_cntr);
			  }
			if (verbosity > 1) EH.println("Major i/o accomplished:");
			EH.println("Starting major loop...");

			optima = new ArrayList();
			seekNewMultistartSolution = true;
			if (doMultistart && !doNoiseCooling && !doRelaxedSmoothing) MSSM.initialiseStrategies();

			if (dumpConvDiagnostics) convergenceResiduals = new ArrayList();

			while (seekNewMultistartSolution)
			  {
				if (dumpConvDiagnostics) 
				  {
					convergenceResidualArray = new cern.colt.list.FloatArrayList();
					convergenceResiduals.add(convergenceResidualArray);
				  }
				if (doMultistart)
				  {
					// in this mode, these numbers need to be clean
					rlast = Float.MAX_VALUE;
					nit = 0;
					itmax = maxitr;
				  }
				// start iterative loop for local optimization
				keep_iterating = true;
				while (keep_iterating)
				  {
					r = rlast;
					tobt = GaussNewtonUpdate(tobt,nit,doBayesianAutoSmoothing,multistart_cntr);
					nit++;

					// write an iteration file
					itrout(extraText,maxitr,rlast,tobt,0,multistart_cntr);
					// let the console know roughly what is going on
					if (verbosity > 2) 
					  {
						EH.println("Iteration "+nit+" RMS misfit= "+tobt);
					  }

					// test for convergence: continue inversion if 
					// 1) iteration number does not exceed maximum, and
					// 2) stepsize is still significant, and
					// 3) rms misfit larger than required, and 
					// 4) there are no irrecoverable convergence problems

					if (verbosity > 2) EH.println("Chisq_mismatch= "+tobt);
					// JG...
					tobt_full = nd * tobt*tobt;
					noise_stddev = tobt;
					if (doNoiseEstimation)
					  {
						// tobt_full = (float)(nd * (1.0+Math.log(tobt_full/nd))); // incorporates optimization on noise-scalar
						// EB version of above
						tobt_full *= sqr(global_error_scale); // restores misfit wrt original data-file noise levels
						dof = DOF(nd,np);
						noise_stddev = (float)sqrt(tobt_full/dof);
						tobt_full = dof + nd*(float)Math.log(sqr(noise_stddev));  
					  }
					chisq_prej = chisqPrej(pmu,pm);
					tobt_full += fpow10(pmu) * chisqRoughness(pm) + chisq_prej;

					if (doBayesianAutoSmoothing && !doEmpBayes) tobt_full += negLogDetCpEffInv(pmu,true);

					// compute  norm of step
					float step_size = (float)sqrt(anorm(pwk3));

					if (verbosity >2 ) 
					  {
						EH.println(" -2log_posterior_mismatch= %-15.8g    -2log_posterior_full= %-15.8g",
								   (nd * tobt*tobt),tobt_full);
						// typical convergence termination criterion
						EH.println(" objective function decrease: "+(tobt_full-tobt_full_old));
						EH.printf("|model| step size= %-15.8g%n",step_size);
						EH.printf("|grad|= %-15.8g |grad|_test= %-15.8g%n",(float)sqrt(anorm(grad)),CUBE_ROOT_EPS_CHISQ*(tobt_full+1.0));
					  }
					if (dumpConvDiagnostics) convergenceResidualArray.add(tobt_full);


					termTest = getTerminationCriteria(tobt_full,tobt_full_old,step_size);


					if (doBayesianAutoSmoothing && avoidTerminationBayesianAutoSmoothing)
					  {
						// do nothing - not allowed to terminate on a step with mu  frozen
					  }
					else if ((nit >= 1) && termTest)
					  /*
						if ((nit >= 1) && 
						(((abs(tobt_full-tobt_full_old) <= EPS_CHISQ) &&  // conventional Newton stop
						(tobt_full<=tobt_full_old) &&           // ... on an *improvement*
						(step_size < SQRT_EPS_CHISQ))     // ... having taken a small step..U2 in  GMW
						|| (konv == 2))) // line search degenerate
					  */
					  {
						if (verbosity > 1) 
						  {
							EH.println("Satisfactory full convergence stop");
							showParams(pwk2,true);
						  }
						H = computeFullHessian(pmu);
						write_final_it = true;
						keep_iterating = false;
					  }
					tobt_full_old = tobt_full;
		
					if (keep_iterating) 
					  {
						if (konv != 0)
						  {
							// out_logger.println(" convergence problems of type "+konv);
							if (verbosity > 1) EH.println(" convergence problems of type "+ konv);
							keep_iterating = false;
						  }
						else if (nit >= itmax) 
						  {
							// out_logger.println(" stop on maximum iterations");
							if (verbosity > 1) EH.println(" stop on maximum iterations");
							write_final_it = true;
							keep_iterating = false;
							// JGHACK in full Bayes these criterion no longer apply
							// else if ((tobt <= 1.01*tolreq)) then
							// write (iof1,*) ' stop on converged model'
							// EH.println(" stop on converged model'
							// close(iof1)
							// stop
							// else if ((tobt<=1.01*tolreq) && (abs(r-rlast)<=0.00001)) then
							// write (iof1,*) ' stop on converged roughness'
							// EH.println(" stop on converged roughness'
							// close(iof1)
						  }
					  }
				  } // major local-optimisation loop: 

				if (write_final_it) 
				  {
					// write final iteration file
					itrout(extraText,maxitr,rlast,noise_stddev,2,multistart_cntr);
					// uncertainty files
					modeCov = computeUncertainty();
					writeUncertaintyFiles(extraText,true,false,multistart_cntr,modeCov);
			
					// To show cross sections through local optimum, 
					if (doXsections) Xsections(pmu,false,multistart_cntr);
 
					// collect optima in list
					if (verbosity > 1) EH.println("Multistart solution "+multistart_cntr+" obj="+tobt_full);
					if (!foundSimilarExistingMode(tobt_full/2.0f))
					  {			
						LocalMode mode = new LocalMode(np,current_resp_filename,current_profile_filename,current_iter_filename,
													   pm,pmu,tobt_full/2.0f,chisq_prej, noise_stddev, 
													   // wjtwj,modeCov.stdv);
													   H,modeCov.stdv);
						optima.add(mode);
					  }
					if (!doMultistart)
					  { 
						// single pass only. Bail out ASAP 
						seekNewMultistartSolution = false;
					  }
					else if (optima.size() > 0)
					  {
						// Multistart restarts only sensible if at least one optima found so far.
						// set new start point
						if (MSSM.tryPhaseFlipping && MSSM.phaseFlipping)
						  {
							// change start points suitably, and recommence
							itrout(extraText,maxitr,rlast,noise_stddev,2,multistart_cntr); 
							// Restart value based on first optimum..
							LocalMode mode0 = optima.get(0);
							// search forward for break in profile likely to be ambiguous...
							seekNewMultistartSolution = false;
							for (k = MSSM.k_flips ; k< (np-1); k++)
							  {
								dist = abs(mode0.pmLoc[k] - mode0.pmLoc[k+1]);
								if (mode0.stdv[k] > 0) dist/= mode0.stdv[k];
								// large dist value indicated possible "ambiguous" dislocation in solution
								// between indices k,k+1. Start new optimization based on flipped pair
								if (dist > MULTISTART_DIST_THRESH)
								  {
									System.arraycopy(mode0.pmLoc,0,pm,0,np);	// pm <- pmLoc   
									// flip layers & fire off new optimization
									m_tmp = pm[k];
									pm[k] = pm[k+1];
									pm[k+1] = m_tmp;
									seekNewMultistartSolution = true;
									MSSM.k_flips = k+1; // next search starts here
									if (verbosity > 1) EH.println("MultistartNewSolution "+multistart_cntr);
									break;
								  }
							  }
							if ((k >= (np-1)) && !seekNewMultistartSolution) 
							  {
								MSSM.exitFlipping(); // opens up random...
							  }
						  }
						if (MSSM.tryRandom && MSSM.phaseRandom)
						  {
							// random starts
							if (MSSM.random_ctr < multistartRandomN)
							  {
								LocalMode mode0 = optima.get(0);
								for (kk=0;kk<np;kk++)
								  {
									pm[kk] = mode0.pmLoc[kk] + 2.0f*((float)Math.random() - 0.5f);	
								  }
								conformStartupToBounds();
								seekNewMultistartSolution = true;
								if (verbosity > 1) EH.println("MultistartNewSolution  (random) " + MSSM.random_ctr);
								MSSM.random_ctr++;
							  }
							else
							  {
								seekNewMultistartSolution = false;
								MSSM.exitRandom(); 
							  }
						  }
						if (seekNewMultistartSolution)
						  {
							multistart_cntr++;
							// go back and try new local optimization
						  }
						else
						  {
							// give up multistart. Take what we've got and push on
							seekNewMultistartSolution = false;
						  }
					  }
				  } // last optimization "successful"
				else
				  {
					// some kind of convergence problem. Bail out ASAP 
					seekNewMultistartSolution = false;			
				  }
			  }  // multistart search

			cleanUpMultistartSolutions(".level_"+grid_size);
			LocalMode mi = optima.get(0);
			if (mi.logMMLapprox < best_split_min) 
			  {
				best_split_min = mi.logMMLapprox;
				best_split_rms = mi.rmsMisfit;
				best_split_indx = i_grid;
				if (best_split_min < best_MML)
				  {
					// also better than solution on previous level
					// store corresponding best solution on fine mesh
					for (i=0,k=0;i<grid.length;i++) 
					  {
						for (j=0;j<grid[i];j++,k++) pm_best_fine[k] = pm[i];
					  }
					msInfo.best_resp_filename = mi.resp_filename;
					msInfo.best_profile_filename = mi.profile_filename;
					msInfo.best_iter_filename = mi.iter_filename;
				  }
			  }
			if (verbosity > 1) 
			  {
				EH.printf("Model solution: level %2d split %2d Gridding: ",grid_size,i_grid);
				for (i=0;i<grid.length;i++) 
				  {
					for(j=0;j<grid[i];j++) EH.printf("%2d ",i);
				  }
				EH.printf(" %-10.6g %-10.6g %-10.6g ",mi.negLogPost, mi.rmsMisfit, mi.logMMLapprox);
				EH.println("");
				if (verbosity > 1) 
				  {
					if (i_grid==0)
					  {
						if (grid_size == grid_start_size)
						  {
							EH.printf("GriddingValues: #DEPTHS ");
							for (j=0;j<NlayerMax;j++) EH.printf("%5g ",j*max_thickness/(double)NlayerMax);
							EH.println("");
						  }
						EH.printf("GriddingValues: ");
						for (i=0,k=0;i<grid.length;i++) for (j=0;j<grid[i];j++,k++) EH.printf("-1 ");;
						EH.println("");
					  }
					EH.printf("GriddingValues: ");
					for (i=0,k=0;i<grid.length;i++) 
					  {
						for (j=0;j<grid[i];j++,k++) EH.printf("%8g ",pm[i]);;
					  }
					EH.println("");
				  }
			  }
		  }  // split loop
		if (best_split_min < best_MML)
		  {
			best_MML = best_split_min;
		
			// and embed best split
			i_grid = best_split_indx;
			for (i=0,j=0; i<old_grid.length; i++,j++)
			  {
				if (i==i_grid) // split layer i_grid
				  {
					grid[j] = old_grid[i]/2;
					grid[++j] = old_grid[i]/2;
				  }
				else grid[j] = old_grid[i];
			  }
			// Save to template for new level
			System.arraycopy(pm_best_fine,0,pm_best_last_level_fine,0,pm_best_fine.length);
			// save best thicknesses & write in correct values
			thickf_best = (float[])thickf.clone();
			for (i=1; i<=nlayer; i++) thickf_best[i] = grid[i-1]*max_thickness/NlayerMax;
			// Now old_grid replace by new one..
			old_grid = grid;
		  }
		else if ((best_split_min < (best_MML + LOGMML_SAFETY_WINDOW))
				 && (best_split_rms > tolreq) && (n_extra_safety_splits<2))
		  {
			if (verbosity > 1) EH.println("Gridding: Inadequate fit to data - attempting further split of first non-minimal layer");
			// best_MML = best_split_min;		
			n_extra_safety_splits++;
			// Enforce split of first non-minimal layer
			best_split_indx = -1;
			for (i=0;i<grid.length;i++) 
			  {
				if (grid[i]>1) 
				  {
					best_split_indx=i;
					break;
				  }
			  }
			if (best_split_indx == -1)
			  {
				if (verbosity > 1) EH.println("Gridding: no further split possible. Giving up");
				break;
			  }
			i_grid = best_split_indx;
			for (i=0,j=0; i<old_grid.length; i++,j++)
			  {
				if (i==i_grid) // split layer i_grid
				  {
					grid[j] = old_grid[i]/2;
					grid[++j] = old_grid[i]/2;
				  }
				else grid[j] = old_grid[i];
			  }
			// Don't save to template for new level - solution not an improvement
			// System.arraycopy(pm_best_fine,0,pm_best_last_level_fine,0,pm_best_fine.length);
			// Now old_grid replace by new one..
			old_grid = grid;
		  }
		else if ((best_split_min < (best_MML + LOGMML_SAFETY_WINDOW))
				 && (thickf[1] > min(thickf[0],minSkinDepth)))
		  {
			if (verbosity > 1) EH.println("Gridding safety - splitting first layer based on sea depth & minimal skin depth");
			// can't expect to see anything unless first layer is no deeper than sea depth
			// best_MML = best_split_min;		
			// Enforce split of layer 1
			best_split_indx = 0;
			i_grid = best_split_indx;
			for (i=0,j=0; i<old_grid.length; i++,j++)
			  {
				if (i==i_grid) // split layer i_grid
				  {
					grid[j] = old_grid[i]/2;
					grid[++j] = old_grid[i]/2;
				  }
				else grid[j] = old_grid[i];
			  }
			// Don't save to template for new level - solution not an improvement
			// System.arraycopy(pm_best_fine,0,pm_best_last_level_fine,0,pm_best_fine.length);
			// Now old_grid replace by new one..
			old_grid = grid;
		  }
		else
		  {
			// No useful split found.
			// Best split was that of level above
			System.arraycopy(pm_best_last_level_fine,0,pm_best_fine,0,pm_best_fine.length);
			grid = old_grid;
			if (verbosity > 1) 
			  {
				EH.printf("Grid splitting stopped. Best%nModel solution: Gridding: ");
				for (i=0;i<grid.length;i++) {for(j=0;j<grid[i];j++) EH.printf("%2d ",i);}
				EH.printf(" minimum: %-10.6g%n",best_MML);
				EH.println("Best solution: "+msInfo.best_resp_filename +" "+msInfo.best_profile_filename);
			  }
			String destFile = "profile_final.txt";
			if (copyFile(msInfo.best_profile_filename,destFile))
			  {
				if (verbosity > 1) EH.println("File "+msInfo.best_profile_filename +" copied to "+destFile);	
			  }
			destFile = "resp_final";
			if (copyFile(msInfo.best_resp_filename,destFile))
			  {
				if (verbosity > 1) EH.println("File "+msInfo.best_resp_filename +" copied to "+destFile);	
			  }
			destFile = "iter_final";
			if (copyFile(msInfo.best_iter_filename,destFile))
			  {
				if (verbosity > 1) EH.println("File "+msInfo.best_iter_filename +" copied to "+destFile);	
			  }
			break;
		  }
		// N.B: pm_best_fine now has best solution
	  }

	if (splitSuggest)
	  {
		outputModel(thickf_best);
	  }

	if (dumpConvDiagnostics) 
	  {
		dumpConvergenceDiagnostics();
	  }
	
	if (THApproxTest)
	  {
		doTHApproxTest();
	  }

	if (doMCMC) 
	  {
		GenerateMCMCPosterior();
	  }
	// out_logger.close();
  }

  // termination criterion, shared by all run() methods

  boolean getTerminationCriteria(float tobt_full,float tobt_full_old,float step_size)
  {
	boolean termTest;

	if (!useLMstrategy)
	  {
		if (true)
		  {
			// Termination criteria of Gill/Murray/Wright Ch.8 8.2.3.2 at this point
			float eps_grad;
			termTest = false;
			eps_grad = nd * GMW_TAU_F_PWR_ONETHIRD; 
			termTest = (abs(tobt_full-tobt_full_old) < GMW_TAU_F*(1.0f+abs(tobt_full)));
			//termTest = termTest && (sqrt(anorm(pwk3)) < sqrt(GMW_TAU_F));
			termTest =  termTest && (getInfInactiveNormGrad() < eps_grad);
			if (verbosity > 2) EH.printf("Gauss-Newton: |grad_inactive|_inf== %-15.8g%n",
										 getInfInactiveNormGrad());		
		  }
		else
		  { // OLD
			termTest = (((abs(tobt_full-tobt_full_old) <= EPS_CHISQ) &&  // conventional Newton stop
						 (tobt_full<=tobt_full_old) &&           // ... on an *improvement*
						 (step_size < SQRT_EPS_CHISQ))           // ... having taken a small step..U2 in  GMW
						|| (konv == 2));                         // line search degenerate
		  }		
	  }
	else
	  {
		// Termination criteria of Gill/Murray/Wright Ch.8 8.2.3.2 at this point
		if (false)
		  {
			float Tau_F = 1e-5f;
			float gr = getInfInactiveNormGrad();
			termTest = (abs(tobt_full-tobt_full_old) < Tau_F*(1.0f+abs(tobt_full)));
			//termTest &= (sqrt(anorm(pwk3)) < sqrt(Tau_F));
			//				termTest &= (getInfInactiveNormGrad() < (pow(Tau_F,0.33)*(1.0f+abs(tobt_full))));
			termTest &= (gr < (pow(Tau_F,0.33)*nd));
			termTest |= (konv==2); 
		  }
		// Termination criteria now captured within Marquardt update as LMinfo.terminateLM
		termTest = LMinfo.terminateLM || (konv==2);                      // Encapsulates above
		if (verbosity > 2) EH.printf("LM: |grad_inactive|_inf== %-15.8g, step_size== %-15.8g, lambda=%-10.5g%n",
									 getInfInactiveNormGrad(),step_size,LMinfo.lambda_LM);			  
	  }

	return termTest;
  }


  float DOF(int Nd,int Np)
  {
	int dof = Nd-Np;

	return (dof>0) ? (float)dof : 1f;
  }

  void performPBsampling()
  {
	// Explicit loops to try and do bootstrap calculations
	int i,j,Nptot;
	float pmu_best;
	double y,rms;
	float[] pm_best;
	float[] d_true;
	float[] premod_orig;
	LocalMode m0;
	float LOGMU=-10.0f; // no smoothing at present
	float chisqPrej,chsqNLPmin=0,chsqNLPminActual=0;
	double errorScale=1.0;
	PrintWriter out;



	String PBFilename = (samplesFilename == null) ? "PBsamples.txt" : samplesFilename;


	// get best model. Use actual mean
	resampleMeansForBayesianBootstrapping = false;

	// Inversion with real data
	run(true);

	

	// Fix according to best LocalMode solution.
	Nptot = np + SEM.nPsystematics;
	m0 = optima.get(0);
	pmu_best = m0.pmuLoc;
	pm_best = (float[])(m0.pmLoc).clone();

	System.arraycopy(m0.pmLoc,0,pm_best,0,np);

	if (doNoiseEstimation)
	  {
		// for now, freeze noise at estimated level, and switch off
		// adjust also for guessed DOF factor
		errorScale = m0.rmsMisfit; // this is the EB estimated noise level
		global_error_scale  = m0.rmsMisfit;  // set this so the inversion routines in run() see the right noise level
		if (verbosity>1) EH.printf("Freezing noise level adjustment at %6g%n",errorScale);
		doNoiseEstimation = false;
	  }

	// this is true, if no smoothing etc
	chsqNLPmin = (float) ( (nd-np) + m0.negLogNoSmoothingPrior); 
	chsqNLPminActual = nd * m0.rmsMisfit*m0.rmsMisfit + m0.negLogNoSmoothingPrior; 

	if (verbosity>1) EH.printf("Expected minimum chisq_total of samples = %6g, actual = %6g%n",chsqNLPmin, chsqNLPminActual);


	// loop over synthetic data generation, inversion, etc

	if (verbosity>1) EH.println("Starting parametric bootstrap major loop..");


	// copy real data & errors
	d_true = (float[])d.clone();
	sd_orig = (float[])sd.clone();

	// use resampled means: but copy true one
	premod_orig = (float[])premod.clone();
	resampleMeansForBayesianBootstrapping = true;

	// Open PBsamples file, write header
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(PBFilename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + PBFilename + " for output");
	  }
	out.printf("#DEPTHS ");
	// top depths of free layers
	for (i=0; i<(nlayer+1); i++) if (rho[i] <= 0.0f) out.printf("%10.5g ",depths[i]);
	out.println("");
	// bottom model depth
	out.printf("#BASE_DEPTH %10.5g%n",depths[nlayer+1]);
	out.printf("#nd= %d%n",nd);
	out.printf("#np= %d%n",np);
	out.printf("#Nptot= %d%n",Nptot);
	out.printf("#error-rescaling-factor= %g%n",errorScale);
	out.printf("#expected-min-chisq_total= %g%n",chsqNLPmin);
	out.printf("#actual-min-chisq_total= %g%n",chsqNLPminActual);
	out.printf("# ");
	for (i=1;i<=Nptot;i++) out.printf("m_%d ",i);
	out.println("rms_misfit  chisq_data chisq_prior chisq_total rms_misfit_orig_noise rms_misfit_orig_noise_pre_dilution");

	
	
	for (i=1; i<= N_PB;i++)
	  {
		if (verbosity>1) EH.println("Making synthetic data set "+i);
		// restore MAP model
		System.arraycopy(pm_best,0,pm,0,Nptot);
		// Original mean of earth model
		System.arraycopy(premod_orig,0,premod,0,np);
		pmu = pmu_best;
		doForwardModel(pm);
		createSyntheticFromCurrentModel(true,errorScale, quietIO ? null : ("synthetic_data."+i+".dat"));

		// invert - no need to reread data etc - is all in memory
		// If you want to recenter mean..
		if (doBootstrappingRecenter) System.arraycopy(pm_best,0,premod,0,np);
		run(false);
		// set best solution in-situ
		System.arraycopy(optima.get(0).pmLoc,0,pm,0,np);

		y = optima.get(0).rmsMisfit;

		// test fit against original true data, but with re-estimated noise level

		// Overwrite d, so d_true not corrupted
		if (errorScale != 1.0) System.arraycopy(sd_orig,0,sd,0,nd);
		System.arraycopy(d_true,0,d,0,nd);
		System.arraycopy(premod_orig,0,premod,0,np);

		if (false)
		  {
			EH.printf("TESTING:\n");
			for (j=0; j<np; j++)
			  {
				EH.printf("%g %g %g\n",premod[j],prewts[j],pm[j]);
			  }
			for (j=0; j<10; j++)
			  {
				EH.printf("DD %g %g\n",d[j],sd[j]);
			  }
		  }
		doForwardModel(pm);

		// Dump best model & its characteristic to sample file
		y = chisqMisFit()/sqr(errorScale);

		rms = sqrt(y/nd);
		showParams(pm,false); 

		chisqPrej = chisqPrej(LOGMU,pm);

		// samples
		for (j=0;j<Nptot;j++) out.printf("%-12.7g ",pm[j]);
		// fit characteristics: rms_misfit, chisq_data, chisq_prior, chisq_total, rms_misfit_orig_noise rms_misfit_orig_noise_pre_dilution
		out.printf("%-12.7g %-12.7g %-12.7g %-12.7g %-12.7g %-12.7g%n",rms, y , chisqPrej, (y+chisqPrej), rms*errorScale, rms*errorScale*diluteFactor);
		out.flush();
	  }

	out.close();

	if (verbosity>1) EH.println("Wrote Bootstrap samples to "+PBFilename);
  }

  void performModelStudies()
  {
	int i,j,i_free,Nstud;
	float x,tmp;
	String s;
	String[] layer_spec,spec;  
	LocalMode m0;
	ModelStudyList[] _sL;
	float[] pm_best=null;

	// Parse commandline

	layer_spec = modelStudy.split(":");
	Nstud = layer_spec.length;
	if (Nstud>3) 
	  {
		EH.warn("Model study "+layer_spec+" will take a long time!");
	  }

	modelStudyList = new ModelStudyList[Nstud];
	_sL = modelStudyList;
	for (i=0;i<layer_spec.length;i++)
	  {
		s = layer_spec[i];
		spec = s.split(",");
		if (spec.length!=4) 
		  {
			EH.die_err("Model study layer spec "+s+" needs to be layer_num,start,spacing,end");
		  }
		else
		  {
			_sL[i]= new ModelStudyList();
			// layer
			try {_sL[i].layer = Integer.parseInt(spec[0]);}
			catch  (NumberFormatException e) {
			  EH.die_err(spec[0] +" not an integer in "+s+ ", inside "+modelStudy);}

			// start
			try {_sL[i].start = Float.parseFloat(spec[1]);}
			catch  (NumberFormatException e) {
			  EH.die_err(spec[1] +" not a number in "+s+ ", inside "+modelStudy);}
			// increment
			try {_sL[i].inc = Float.parseFloat(spec[2]);}
			catch  (NumberFormatException e) {
			  EH.die_err(spec[2] +" not a number in "+s+ ", inside "+modelStudy);}
			// end
			try {_sL[i].end= Float.parseFloat(spec[3]);}
			catch  (NumberFormatException e) {
			  EH.die_err(spec[3] +" not a number in "+s+ ", inside "+modelStudy);}
		  }
	  }

	// Now do study
	// First, get model size
	read_startup_and_model(); // gets nlayer for header

	// Open output file, write header
	PrintWriter out;
	String filename = "ModelStudy.txt";
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }
	out.printf("# ");
	for (j=0;j<Nstud;j++) out.printf("m_%d ",_sL[j].layer);
	for (j=0;j<Nstud;j++) out.printf("rho_%d ",_sL[j].layer);
	out.printf("logMMLapprox rmsMisfit ");
	for (j = 0; j<(nlayer+2); j++) out.printf("d_%d ",j); // depths, including base
	// finally, complete resistivities
	for (j = 0; j<(nlayer+1); j++) out.printf("log10rho_%d ",j); // depths
	out.printf("%n");

	int[] modSamplesUnitCube = new int[Nstud];
	for(i=0;i<Nstud;i++)
	  {
		tmp = _sL[i].inc;
		modSamplesUnitCube[i] = (int)rint((tmp>0)? ((_sL[i].end - _sL[i].start)/tmp + 1) : 1);
	  }

	float[][] hc = au.csiro.JamesGunning.lib.Numerical.Hypercube.createHypercubeList(modSamplesUnitCube);

	// loop over models in hypercube
	double[] mmls = new double[hc.length]; 
	_sL = modelStudyList;
	for (i=0;i<hc.length;i++)
	  {
		// set values
		for(j=0;j<Nstud;j++)
		  {
			_sL[j].current_value = hc[i][j] * (_sL[j].end - _sL[j].start)+ _sL[j].start; 
		  }

		// Set startup, conform it
		if ((i>0) && (i % (modSamplesUnitCube[Nstud-1])) > 0)
		  {
			// within raster scan of fastest dimension => should be close to last solution: use it to start with...
			for(j=0;j<np;j++) pm[j] = pm_best[j];
		  }
		else
		  {
			for(j=0;j<np;j++) pm[j]=0;
		  }

		if (i!=0) conformStartupToBounds(); 
		// reset nit. Resets LM apparatus etc
		nit=0;

		// run, reading data first time only
		run(i==0);

		// Fix according to best LocalMode solution.
		m0 = optima.get(0);
		pm_best = m0.pmLoc;
		// save MML
		mmls[i] = m0.logMMLapprox;

		// IO: first Nstud columns are model-parm thicknesses
		for(j=0;j<Nstud;j++) out.printf("%-12.7g ",_sL[j].current_value);
		// then resistivities
		for(j=0;j<Nstud;j++) out.printf("%-12.7g ",fpow10(pm_best[_sL[j].layer-1]));  // Assumes Legacy scripps indexing!
		// then MML, RMS
		out.printf("%-12.7g %-12.7g ",m0.logMMLapprox,m0.rmsMisfit);
		// complete depths
		for (j = 0; j<=(nlayer+1); j++) out.printf("%12.7g ",depths[j]);
		// finally, complete resistivities
		for (j = 0,i_free=0; j<(nlayer+1); j++) 
		  {
			out.printf("%12.7g ",((rho[j] <= 0.0f) ? pm[i_free++] : ((j==0) ? flog10(rho[j]) : rho[j])));
		  }	  
		out.printf("%n");
		// for(j=0;j<=nlayer;j++) out.printf("%-12.7g ",thickf[j]);
		//pmu_best = m0.pmuLoc;
		//pm_best = (float[])(m0.pmLoc).clone();
		
		//System.arraycopy(m0.pmLoc,0,pm_best,0,np);
		optima.clear(); 
	  }
	out.close();
  }

  /* 
	 Study fixed input model using resampling of data using classical
	 bootstrap (sampling with replacement
  */
  void performResamplingBootstrapStudy()
  {
	int i,j,k,Nb,iB;
	float mml0=0;
	DataComparator cmp = new DataComparator();
	PrintWriter out;

	// First, get model size
	read_startup_and_model(); // gets nlayer for header

	// Open output file, write header
	String filename = "BootstrapStudy.txt";
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }
	out.printf("# Bootstrap values: logMML rmsMisfit%n");

	// run, reading data first time only
	run(true);

	// Copy original data
	ArrayList<SimpleEMData> data_list_cpy = new ArrayList<SimpleEMData>(nd);

	for (i=0;i<nd;i++) 
	  {
		SimpleEMData EM = new SimpleEMData(6);
		EM.d = d[i];
		EM.sd = sd[i];
		for (j=0;j<6;j++) EM.dp[j] = dp.getQuick(j,i);
		data_list_cpy.add(EM);
	  }



	// Bootstrap
	Nb=NBootstrap;
	ArrayList<SimpleEMData> data_list_BS = new ArrayList<SimpleEMData>(nd);
	for (iB=0;iB<Nb;iB++)
	  {
		// sample with replacement 
		data_list_BS.clear();
		for(j=0;j<nd;j++)
		  {
			k = min((int)(Math.random()*nd),nd-1); // random index in [0...nd-1]
			data_list_BS.add(data_list_cpy.get(k));
		  }
		// sort data by freq, range, type...
		Collections.sort(data_list_BS,cmp);
		for (i=0;i<nd;i++)
		  {
			SimpleEMData data = data_list_BS.get(i);
			d[i] = data.d;
			sd[i] = data.sd;
			for (j=0;j<6;j++) dp.setQuick(j,i,data.dp[j]);
		  }

		conformStartupToBounds(); 
		// reset nit. Resets LM apparatus etc
		nit=0;

		// Alternate between these 2 models...
		for(j=0;j<2;j++)
		  {
			releaseNativeModelPointers();
			if (j==0)
			  {
				thickf[1] = 950.0f;
				thickf[2] = 110f;
			  }
			else
			  {
				thickf[1] = 850.0f;
				thickf[2] = 250f;
			  }
		
			// re-compute depths
			depths[0] = 0f;
			for (i=1 ; i < nlayer+2 ; i++)
			  {
				depths[i] = depths[i-1]+thickf[i-1];
			  }

			setNativeModelPointers();

			// run, not rereading data 
			run(false);

			// Fix according to best LocalMode solution.
			LocalMode m0 = optima.get(0);
			float[] pm_best = m0.pmLoc;
			if (j==0)
			  {
				mml0 = m0.logMMLapprox;
			  }
			else
			  {
				// save MML
				out.printf("%-12.7g %-12.7g %-12.7g%n",m0.logMMLapprox, m0.rmsMisfit, m0.logMMLapprox - mml0);
				// mmls[i] = m0.logMMLapprox;
			  }
			optima.clear(); 
		  }
	  }
	out.close();
	EH.printf("Wrote %s%n",filename);
  }

  private float GaussNewtonUpdate(float tobt,int nit,
								  boolean doBayesianAutoSmoothing, 
								  int i_multistart)
  {
	float rho_Goldstein = 0.25f;
	float tolm= 0.1f, toli = 0.001f;
  
	int i,j,Nptot;
	int rescale_count;
	float tobt_full0 = 0f;
	float[] dmhat,pm0; // npp
	float[] dmhat_norm = new float[1];
	float slope,q,alpha_new;

	// JG additions

	float tol0, tmin=0f, tint=0;
	float ruf=0f, pmu2,tobt_full=0f,tobt_full_pre_mu_opt, tobt_new;
	float chiSqData,chiSqData0,chiSqRough,chiSqRoughFull,chiSqPrej,negLogDet;
	float rhs_Goldstein1=0f,rhs_Goldstein2;
	boolean Goldstein_test_fails;

   	float tmp;

	float pmuopt=0;
	boolean optimise_mu_Newton;

	Nptot = np + SEM.nPsystematics;

	dmhat= new float[Nptot];
	pm0 = new float[Nptot];
  
  	// Grad is global in scope, but allocate here....
	grad = new float[Nptot];

	optimise_mu_Newton = false; // FIXME not finished yet - leave false.....;
  
	// pmu controlled external to here on input except initial eigenvalue guesses, or doBayesianAutoSmoothing mode
	// if (!doBayesianAutoSmoothing && !doRelaxedSmoothing && !zeroSmoothingWeights) pmu = 0.0f;
  
	// frac controls the step size; normally will remain at 1.0 unless we have
	// convergence problems
	frac = 1.0f;
	

	// construct wjtwj, and find initial misfit (JG: with model currently in dm)
	if (verbosity > 2) EH.println(" constructing derivative dependent matrices..");
	
	tol0 = makjtj();


	if (nit == 0) 
	  {
		// assume model is feasible
		ifftol = 1;
	  }
	if (verbosity > 2) EH.println("Starting RMS = "+tol0);

	// 
	// construct the penalty matrix multiplied into its transpose
	if (verbosity > 2) EH.println(" constructing penalty matrix..");
	// Very first iteration: guess from eigenvalue balance
	if ((doBayesianAutoSmoothing || doRelaxedSmoothing) && !zeroSmoothingWeights) 
	  {
		if (nit == 0 && (doBayesianAutoSmoothing || (doRelaxedSmoothing && nit_ms==0)))
		  {
			if (nit_ms==0) computeLargestEigenValuePTP();
			if ((largestEigenvaluePTP > 0) && (largestEigenvalueWJTWJ>0))
			  {
				pmu = flog10(largestEigenvalueWJTWJ/largestEigenvaluePTP);
				// start Newman more cautiously...				
				if (doRelaxedSmoothing) pmu += flog10(INITIAL_SMOOTHING_EXCESS_WEIGHT);
				if (verbosity > 2) EH.printf("Eigenvalue-balance guess initial pmu = %4g%n",pmu);
			  }
		  }
	  }
	if (useLMstrategy)
	  {
		if (nit == 0)
		  {
			if (verbosity > 2) EH.println("LM: Scale resetting");
			LMinfo.setLMscalesByDiag(wjtwj);
		  }
		if (LMinfo.restartLM)
		  {
			LMinfo.lambda_LM = (float)getLargestEigenvalue(wjtwj);
			LMinfo.nu = 2.0f;
			// guess of scales based on Hessian diagonal
			LMinfo.setLMscalesByDiag(wjtwj);
			EH.println("LM: Reinitialize-emergency");
			LMinfo.restartLM = false;
		  }
	  }

	makptp(pmu);   // forms ptp = del^T.del + W_B.W_B
	// Trick to get model writeout for *initialisation* run
	if (verbosity > 2) showParams(pm,false);
	if (nit==0) 
	  { 
		if (makeSynthetic)   
		  {			
			createSyntheticFromCurrentModel(false,-1.0,"synthetic_data.dat");   // recomputes error..
			rlast = tobt = 0;
			itrout(null,maxitr,rlast,tobt,2,i_multistart);
			writeProfileSnapshot("profile_synthetic.txt");			
			System.exit(0);
		  }
		itrout(null,maxitr,rlast,tobt,1,i_multistart);
	  }

	// JG need full objective function in determination of step size
	chiSqData0 = nd * tol0*tol0; // chisq_data, data misfit
	if (doNoiseEstimation)
	  {
		// tobt_full0 = (float)(nd * (1.0+Math.log(chiSqData0/nd)));     // incorporates optimization on noise-scalar
	  	tobt_full0 = DOF(nd,np) + nd*(float)Math.log(sqr(global_error_scale));  // EB version of above, since global_error_scale = sqrt(Chisq_misfit/DOF)
	  } 
	else
	  {
		tobt_full0 = chiSqData0;
	  }
	chiSqRough = chisqRoughness(pm); // roughness component without mu multiplier

	// the next block of code controls the selection of the smoothing parameter mu or alpha
	// 
	if (doBayesianAutoSmoothing && (!optimise_mu_Newton))
	  {
		// OLDMUOPT tobt_full_pre_mu_opt = tobt_full0 + fpow10(pmu)*chiSqRough;
		tobt_full_pre_mu_opt = tobt_full0 + fpow10(pmu)*chiSqRough + chisqPrej(pmu,pm);
		if (!doEmpBayes) tobt_full_pre_mu_opt += (float)negLogDetCpEffInv(pmu,true);


		if (doEmpBayes || (tol0 < 10.0))  // maybe too aggressive early on here, in simple ML method, so delay
		  {
			// Optimise mu based on current point in pm(), 1D optimisation
			if ((nit % BAYES_SMOOTHING_MU_UPDATE_SKIP) == 0) 
			  {
				// perform a mu update every BAYES_SMOOTHING_MU_UPDATE_SKIP (=5, say)
				pmuopt = minimiseMuOverDetCpAndOtherMuPieces(pmu);
				avoidTerminationBayesianAutoSmoothing = false; // can only terminate on one of these steps
			  }
			else
			  {
				pmuopt = pmu; // stay with current mu
				avoidTerminationBayesianAutoSmoothing = true;
			  }
		  }
		else
		  {
			pmuopt = pmu; // stay with current mu
			avoidTerminationBayesianAutoSmoothing = false; 
			if (verbosity>2) EH.println("Skipping muopt: chisq_data too large");
		  }
		chiSqRoughFull = fpow10(pmuopt)*chiSqRough;
		chiSqPrej = chisqPrej(pmuopt,pm);
		tobt_new = tobt_full0 + chiSqRoughFull + chiSqPrej;
		if (!doEmpBayes)
		  {
			negLogDet = (float)negLogDetCpEffInv(pmuopt,true);
			tobt_new += negLogDet;
		  }
		float pmu_balanced = flog10(tobt_full0/((chiSqRough>0) ? chiSqRough : 1.0f));
		pmu_balanced = min(pmu_balanced , PMU_DEFAULT_START);
		if (pmu_balanced > pmuopt)
		  {
			if (true)
			  {
				// Do nothing!

				// Heuristic: sometimes works OK
				// pmuopt = 0.5f*(pmu_balanced+pmuopt);
			  }
			else
			  {
				// set as a ratio ????
				pmuopt = flog10(0.1 * pow(10.0,pmu_balanced));
			  }

			chiSqRoughFull = fpow10(pmuopt)*chiSqRough;			
			// OLDMUOPT tobt_new = tobt_full0 + chiSqRoughFull;
			chiSqPrej = chisqPrej(pmuopt,pm);
			tobt_new = tobt_full0 + chiSqRoughFull + chiSqPrej;
			if (!doEmpBayes)
			  {
				negLogDet = (float)negLogDetCpEffInv(pmuopt,true);
				tobt_new += negLogDet;
			  }

			// EH.printf("NEW   TERM BALANCE %g %g %g %g%n",tobt_full0,chiSqRoughFull,negLogDetCpEffInv(pmu,true),pmuopt);
			// force acceptance in test branch below
			// tobt_full_pre_mu_opt = tobt_new+1;
		  }

		if ((tobt_full_pre_mu_opt > tobt_new) || doEmpBayes)
		  {
			// Genuine improvement, or with EmpBayes, function is not monotonic
			pmu = pmuopt;
			// Set optimisation tobt here omitting det() term, as this is not present below
			tobt_new = tobt_full0 + chiSqRoughFull + chiSqPrej;
		  }
		else
 		  {
			// Algorithm is starting to cycle. Freeze at old pmu
			pmuopt = pmu;
			// OLDMUOPT tobt_new = tobt_full0 + fpow10(pmu)*chiSqRough;
			// Set optimisation tobt here omitting det() term, as this is not present below
			tobt_new = tobt_full0 + fpow10(pmu)*chiSqRough + chisqPrej(pmu,pm);
		  }

		tobt_full0 = tobt_new;
	  }
	else
	  {
		// OLDMUOPT tobt_full0 += pow(10.0,pmu)*chiSqRough;
		tobt_full0 += (pow(10.0,pmu)*chiSqRough + chisqPrej(pmu,pm));
	  }
	if (verbosity > 2) EH.println("Search start: pmu="+pmu+" obj. function = "+tobt_full0);

   	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  {
		// need to recompute ptp since the "prejudice weights" Wb change with mu. 
		// Ensure we recompute prewts
		makptp(pmu);
	  }





	// copy of starting point
	System.arraycopy(pm,0,pm0,0,Nptot);	// pm0 <- pm	   
	rescale_count=0;

	boolean try_again = true;
	while (try_again)
	  {
		for (;;)
		  {
			// look for successful backtrack
			nfor = 0;

			// produce the misfit function if required (this will be used to bracket min)
			// if (abs(idebug) == 1) EH.println( " ..bracketing minimum.."
			if (doBayesianAutoSmoothing) 
			  {
				pmu = pmuopt;
				if (optimise_mu_Newton) 
				  {
					// Do joint Newton step on {m,mu} system              
					// FIXME finish this function....
					EH.println("optimise_mu_Newton not functional yet. Finish this function!");
					System.exit(0);
				  }
				else
				  {
					// Use existing 1D mu optimum found well above.
					// based on current point in pm()
					// This optimisation need only be done for rescale_count==0 
					// - all backtracks are the same...
					tmin = tofmu(pmu,grad);
				  }
			  }
			else
			  {
				// mu fixed. Look for Newton point in np-dimensional problem
				if (verbosity > 2) EH.println( "Skipping mu-optimisation (fixed pmu)");
				tmin = tofmu(pmu,grad);
			  }
			// New point now in pwk2

			// model pwk1 must contain current model....
			// copy of starting point
			System.arraycopy(pwk2,0,pwk1,0,Nptot);	// pwk1 <- pwk2

			// JG need full objective function in determination of step size

			chiSqData = nd * tmin*tmin;
			if (doNoiseEstimation)
			  {
				// tobt_full = (float)(nd * (1.0+Math.log(chiSqData/nd)));   // incorporates optimization on noise-scalar
				// EB version of above, since global_error_scale = sqrt(Chisq_misfit/DOF)
				chiSqData *= sqr(global_error_scale); // restores misfit wrt original data-file noise levels
				tobt_full = DOF(nd,np) + nd*(float)Math.log(chiSqData/DOF(nd,np));  
			  }
			else
			  {
				tobt_full = chiSqData;
			  }
			ruf = chisqRoughness(pwk1);
			tobt_full += pow(10.0,pmu)*ruf;
			tobt_full += chisqPrej(pmu,pwk1);

			if (useLMstrategy)
			  {
				LMinfo.terminateLM = false;
		
				// Test if predicted improvement is about correct based on quadratic model, adjust
				// lambda using Madsen or Fletcher's strategy
				// Observed:
				float delta_f, delta_fP = 0f, qd=0f,ratio;
				delta_f = tobt_full - tobt_full0;

				// Predicted function change from quadratic model.
				// First, model step Delta m
				for (i=0;i<Nptot;i++) pwk3[i] = pwk2[i] - pm[i];				
				// Hessian piece
				for (i=0;i<Nptot;i++)
				  {
					for (qd=0f, j=0; j<Nptot; j++) qd += LMinfo.H_modified.getQuick(i,j)*pwk3[j];
					delta_fP += pwk3[i] * qd;
				  }
				delta_fP *= 0.5;
				// grad piece
				for (i=0;i<Nptot;i++) delta_fP += LMinfo.grad_modified[i] * pwk3[i];
				
				// "gain" ratio of observed decrease to predicted, under quadratic model: 
				// same as rho in  Madsen's notes (e.g. 3.16)
				ratio = delta_f/((delta_fP!=0)?delta_fP:1.0f);
				if (verbosity > 2)
				  {
					EH.println("LM: fN=%12.7g fO=%12.7g lambda=%12.7g df=%10.7g dfP=%10.7g, ratio = %12.7g,",
							   tobt_full,tobt_full0,LMinfo.lambda_LM,delta_f,delta_fP,ratio);
				  }
				// Magic constants below as as per Madsen's notes
				if (true)
				  {
					// Madsen
					float lambda_max=1e20f;
					// Attempt to encode termination criteria of GMW at this point
					float eps_grad;
					boolean GMWtest = false,GMW_df_test;
					eps_grad = nd * GMW_TAU_F_PWR_ONETHIRD; 
					GMW_df_test = (abs(delta_f) < GMW_TAU_F*(1.0f+abs(tobt_full)));
					//GMWtest = GMW_df_test && (sqrt(anorm(pwk3)) < sqrt(GMW_TAU_F));
					GMWtest =  GMW_df_test && (getInfInactiveNormGrad() < eps_grad);
					// EH.println("GRAD_INF = "+getInfInactiveNormGrad()+" GRAD_TEST_EPS= "+eps_grad);
					if (GMWtest)
					  {
						LMinfo.terminateLM = true;
						if (verbosity > 2) printGradient();
						break;
					  }
					// detect also roundoff limit in formation of  delta_f (1e-7 about smallest sensible diff. in SP)
					if (false)
					  {
						if (abs(delta_f) < 5e-7*(1.0f+abs(tobt_full)))
						  {
							EH.println( "BREAK PRECISION LM1");
						    LMinfo.restartLM=true;
							//LMinfo.terminateLM = true; 
							break;
						  }
					  }
					
				  

					if (ratio >= 0)
					  {
						if (tobt_full <= tobt_full0)
						  {
							tmp = LMinfo.lambda_LM * max(0.33333333f, 1 - cube(2*ratio-1));
							if (GMW_df_test && tmp>1)
							  {
								// An increase of lambda. Will only reduce delta_f
								LMinfo.terminateLM = true;
								if (verbosity > 2) EH.println("LM: emergency overflow termination 1");
								//LMinfo.restartLM=true;
								break;
							  }
							if ((tmp >1e-20)  // avoid underflow !
								&& (tmp < lambda_max)) // avoid overflow!
							  {
								// change lambda -- likely an improvement, but may increase...
								if (verbosity > 2) EH.println( "Successful LM point at -2log_posterior_full= "+tobt_full);
								LMinfo.lambda_LM =tmp;
								LMinfo.nu=2.0f;
							  }
							else
							  {
								// desperation...
								if (verbosity > 2) EH.println("LM: over/underflow issues "+LMinfo.lambda_LM );
							  }
							// break only on improvement
							if (GMWtest)
							  {
								LMinfo.terminateLM = true;
								// EH.println( "BREAK LM2");
								printGradient();
								break;
							  }
						  }
						else
						  {
							// increase lambda
							if (verbosity > 2) EH.println( "Unsuccessful LM update: point is uphill!");						  
							if (GMW_df_test)
							  {
								// last step was close, but gradient not quite zero. Probably precision limited????
								LMinfo.terminateLM = true;
								if (verbosity > 2) EH.println("LM: emergency overflow termination 3");
								
								break;
							  }

							if ((tmp = (LMinfo.lambda_LM * LMinfo.nu)) < lambda_max) // avoid overflow!
							  {
								LMinfo.lambda_LM = tmp;
								LMinfo.nu *= 2.0;
							  }
							else
							  {
								// desperation...
								if (verbosity > 2) EH.println("LM : BIG LAMBDA STAGNATION ISSUES");
								LMinfo.restartLM=true;
							  }
						  }
					  }
					else
					  {
						if (GMW_df_test)
						  {
							// last step was close, but gradient not quite zero. Probably precision limited????
							LMinfo.terminateLM = true;
							if (verbosity > 2) EH.println("LM: emergency overflow termination 4");
							break;
						  }
						if ((tmp = (LMinfo.lambda_LM * LMinfo.nu)) < lambda_max) // avoid overflow!
						  {
							// increase lambda
							if (verbosity > 2) EH.println( "Unsuccessful LM update: point is uphill!");						  
							LMinfo.lambda_LM = tmp;
							LMinfo.nu *= 2.0;
						  }
						else
						  {		
							// desperation...
							if (verbosity > 2) EH.println("LM ISSUES: lambda overflow??");
							LMinfo.restartLM=true;
						  }						  
					  }
				  }
				else
				  {	
					// Fletcher/Marquardt
					if (tobt_full > tobt_full0)
					  {
						if (verbosity > 2) EH.println( "Unsuccessful LM update: point is uphill!");
						LMinfo.lambda_LM *= 4.0; 
					  } 
					else if ((ratio < 0.2) && (LMinfo.lambda_LM < 1e20)) // avoid overflow!  
					  {
						LMinfo.lambda_LM *= 4.0; 
					  }
					else if (ratio > 0.75) 
					  {
						if (verbosity > 2) EH.println( "Successful LM point at -2log_posterior_full= "+tobt_full);
						LMinfo.lambda_LM /= 2.0;
						break;
					  }
					else break;
				  }
			  }

			if (true) 
			  {
				// JG: get slope, as part of Goldstein or Wolfe conditions
				slope = getDirectionalSlope(pwk1,pm0,grad,dmhat,dmhat_norm);
				// note here, dmhat_norm reduces in prop. to 1/frac, so
				// dmhat_norm*slope = alpha * f"(0) of Fletcher eqn 2.5.1
				if (verbosity > 2) EH.println("Goldstein slope f'(0) = "+slope);
				if (slope > 0)
				  {
					if (verbosity > 2) EH.println("BUG: line search in ascent direction. FIXME");
				  }
				rho_Goldstein=1e-2f;
				rhs_Goldstein1 = tobt_full0 + rho_Goldstein*dmhat_norm[0]*slope;
				// rhs_Goldstein2 = tobt_full0+(1.0-rho_Goldstein)*dmhat_norm*slope
				Goldstein_test_fails = (tobt_full > rhs_Goldstein1);
				// *                || (tobt_full<rhs_Goldstein2)
			  }
			else
			  {
				// old very approximate test
				rhs_Goldstein1 = tobt_full0;
				// rhs_Goldstein2 = 0.0
				Goldstein_test_fails = (tobt_full>rhs_Goldstein1);
			  }

			if (verbosity > 2) 
			  {
				EH.println("Line search: [frac=%7.4g, old obj. func=%-12.7g], new obj. func=%-12.7g",
						   1.0/frac,tobt_full0,tobt_full);
			  }


			// JG: simple replacement
			if (useLMstrategy)
			  {
				if (LMinfo.terminateLM)
				  {
					if (verbosity > 2) EH.println("Degenerate (successful?) LM search at -2log_posterior_full= "+tobt_full);
					// return
					break; // jump out of backtrack loop
				  }
				else
				  {
					break; 
				  }
			  }
			else
			  {
				// Line-search methods
				if (Goldstein_test_fails)
				  {
					// out_logger.println( "divergence problems, cutting step size");
					if (true)
					  {
						// standard factor-of-2 backtrack...
						frac *= 2.0;
					  }
					else if ((rescale_count>0) && ((tobt_full - tobt_full0)<1000.0)) 
					  { 
						// JG: experimental simple halving algorithm above works better..(hence the .false.)
						// When close enough, try parabolic search?
						// new point in pwk1, starting point in dm0
						q = (tobt_full - tobt_full0 - slope/frac)/(1.0f/(frac*frac));
						alpha_new = -slope/(2.0f*q);
						frac = min(1.0f/alpha_new,1.0e+05f);
						if (verbosity > 2) EH.println( "Using parabolic-fitted frac= "+frac);
					  }
					
					if (frac > 1.0e+05) 
					  {
						// we have cut step size a lot to no effect: give up
						// *  JG this can also happen if Newton step if virtually perfect....
						// so it shouldn't look like something's wrong...;
						if (verbosity > 2) EH.println("Degenerate (successful?) line search at -2log_posterior_full= "+tobt_full);
						// return
						break; // jump out of backtrack loop
					  }
					else
					  {
						// backtrack...
						rescale_count++;   // back to backtrack loop
					  }
				  } // Goldstein test
				else
				  {
					// Satisfied Goldstein
					if (verbosity > 2) EH.println( "Successful line search at -2log_posterior_full= "+tobt_full);
					break; // jump out of backtrack loop
				  }
			  }
		  
		  } // .. backtrack loop

		// tolerance is no yet small enough.  we will keep the minimum.
		System.arraycopy(pwk1,0,pwk2,0,Nptot);

		try_again = false;
	  }

	// save new model and compute step size
	if (useLMstrategy && (tobt_full >= tobt_full0))
	  {
		// don't go anywhere!
		System.arraycopy(pm,0,pwk2,0,Nptot);
	  }
	for (i=0;i<Nptot;i++)
	  {
		pwk3[i] = pwk2[i] - pm[i];
		pm[i] = pwk2[i];
	  }
	if (modelSystematics)
	  {
		for (i=np,j=0; i<Nptot; i++,j++) SEM.m_N[j] = pwk2[i];
	  }

	// the stepsize is the actual change in the model, normalised by Nptot:
	// out_logger.println("roughness is = "+ruf);
	rlast = ruf; 
	if (verbosity > 2) EH.println( "leaving GaussNewtonUpdate");

	return tmin;
  }

  void setNativeModelPointers()
  {
	if (newScripps)
	  {
		SMI.CM.setCSEMModelPointers(SMI.OM.nd, SMI.DM.nlay1D);
		SMI.DM.setDipoleModelPointers();
		SMI.OM.setOccamModelPointers(SMI.CM.cnDataParams);
	  }
	else
	  {
		ScrippsLegacyModelInterface.setModelPointersNative(np, nd, nlayer, nfreq, thickf,rho,freqs,ipen1,ipen2,iff);
	  }
	nativeModelPointersSet = true;
  } 

  void releaseNativeModelPointers()
  {
	if (nativeModelPointersSet)
	  {
		if (newScripps)
		  {
			SMI.CM.releaseCSEMModelPointers(SMI.OM.nd, SMI.DM.nlay1D);
			SMI.DM.releaseDipoleModelPointers();
			SMI.OM.releaseOccamModelPointers();
		  }
		else
		  {
			ScrippsLegacyModelInterface.releaseModelPointersNative(thickf,rho,freqs,ipen1,ipen2,iff);
		  }

		nativeModelPointersSet = false;
	  } 
  }




  /** 
	  Resamples the means to be used in the inversion from their prior, as required
	  for bootstrap with Bayesian priors
  */
  private void resampleMean()
  {
	int i,j,NTRYMAX = 10000;
	float sample_stddev,sample_mean=0f;

	for (i=0; i<np; i++)
	  {
		prewts[i] = ipen2[iff[i]];
		sample_stddev = (prewts[i]>0) ? 1.0f/prewts[i] : 0.0f;
		if (!useBoundConstraints) 
		  {
			sample_mean = (float)(premod[i] + RD.nextGaussian() * sample_stddev);
		  }
		else
		  {
			// sample a mean inside the bounds
			for (j=0;j<NTRYMAX;j++)
			  {
				sample_mean = (float)(premod[i] + RD.nextGaussian() * sample_stddev);
				if ((sample_mean >= m_lower[i]) && (sample_mean <= m_upper[i])) break;
			  }
			if (j>=NTRYMAX) 
			  {
				EH.printf("WARNING! Could not resample mean inside bounds. Setting bootstrap mean at lower bound!");
				sample_mean = m_lower[i];
			  }
		  }
		premod[i] = sample_mean;
	  }
	if (modelSystematics)
	  {
		// add sample SEM.m_N_bar = L.w from prior
		double[] w_samp = new double[SEM.nPsystematics];
		for (i=0; i<SEM.nPsystematics; i++) 
		  {
			w_samp[i] = RD.nextGaussian();
			for (j=0,sample_mean=0;j<=i;j++) sample_mean += SEM.L_Bsys.getQuick(i,j) * w_samp[j];
			SEM.m_N_bar[i] = sample_mean;
		  }
	  }
	if (verbosity > 1) 
	  {
		EH.printf("Bayesian bootstrap new mean: ");
		for (i=0; i<np; i++) EH.printf("%-10.3g ",premod[i]);
		if (modelSystematics) for (i=0; i<SEM.nPsystematics; i++) EH.printf("%-10.3g ",SEM.m_N_bar[i]);
		EH.println("");
	  }
  }

  void dumpRunmode()
  {
	int i;

	if (verbosity > 0) 
	  {
		EH.println("BayesianAutoSmoothing = "+doBayesianAutoSmoothing);
		EH.println("MCMC = "+doMCMC);
		EH.println("Multistart = "+doMultistart);
		EH.println("Zeroed smoothing weights = "+zeroSmoothingWeights);
		EH.println("Dumping Xsections = "+doXsections);
		EH.println("Dumping profiles  = "+dumpPrfl);
		EH.println("Using Lev. Marq.= "+useLMstrategy);
		if (modelSystematics)
		  {
			EH.println("Modelling systematics: default (|E| or |B| effects), order=%d, ratio=%4g",
					   SEM.systematicsRegressionOrder,SEM.sigmaSystematicTestPoint);
			if (doBayesianAutoSmoothing)
			  {
				EH.die_err("Bayesian smoothing (-s) and systematics (--systematics) not yet jointly supported! A to-do.");
			  }
		  }
	  }
  }

  /**
	 Assemble Hessian matrix
	 A = (WJ)^T WJ + mu * partial^T.partial + W_B^T.W_B
	 On input,
	 ptp   = partial^T.partial + W_B^T.W_B
	 wjtwj = (WJ)^T WJ

	 Sets A = (WJ)^T WJ + mu * partial^T.partial + W_B^T.W_B
  */
  void assembleHessian(float mu, DenseDoubleMatrix2D A)
  {
	int i,j;
	float s;

	// add mu.trans(partial).partial to trans(w.j).w.j
	// By here,
	// wjtwj = (WJ)^T WJ
 	// ptp   = partial^T.partial + W_B^T.W_B

	// copy ptp into a
	A.viewPart(0,0,np,np).assign(ptp);
	// restore pure smoothing operator (undo weights for the model prejudice)
	for (i=0;i<np;i++) A.setQuick(i,i, A.getQuick(i,i) - prewts[i]*prewts[i]);
  
	// Now A = partial^T.partial

	// Add (mu partial^T.partial + W_B^T.W_B). (m_prior - pm) to -GRAD
	int jmax;
	for (i=0; i<np; i++)
	  {
		jmax = min(i+2,np-1);
		s=0.0f;
		for (j=max(0,i-2); j <= jmax; j++) 
		  {
			s += A.getQuick(i,j) * (premod[j] - pm[j]);
		  }
		s = s*mu + prewts[i]*prewts[i] * (premod[i] - pm[i]);
		grad[i] = -(WJdeltaD[i] + s);
	  }
	
	// Get lower sub-block of grad[] from systematics...
	if (modelSystematics)
	  {
		int ii;
		for (i=np,ii=0; ii<SEM.nPsystematics; i++,ii++)
		  {
			// This piece is s_i = C_BsysInv.m_N:
			for (j=0,s=0f; j<SEM.nPsystematics; j++) 
			  {
				s += SEM.C_BsysInv.getQuick(ii,j) * (SEM.m_N[j] - SEM.m_N_bar[j]);
			  }
			grad[i] = -SEM.X_NdeltaD[ii] + s; 
		  }
	  }


	for (i=0;i<np;i++)
	  {
		for (j=0;j<np;j++) 
		  {
			// smoothing plus Hessian
			A.setQuick(i,j, mu*A.getQuick(i,j) + wjtwj.getQuick(i,j)); 
		  }
		// Add prejudice/Bayesian prior covariance W_B^T W_B
		A.setQuick(i,i, A.getQuick(i,i) + prewts[i]*prewts[i]);
	  }
	// Now, (1,1) subblock of A matrix is: 
	// A_11 = (WJ)^T WJ + mu * partial^T.partial + W_B^T.W_B	

	// Add other sub-blocks (12,21,22) from systematics if requested...
	if (modelSystematics && (SEM.nPsystematics>0))
	  {
		// stuff in H12 block
		A.viewPart(0,np,np,SEM.nPsystematics).assign(SEM.H12);
		// H21 by symmetry
		A.viewPart(np,0,SEM.nPsystematics,np).assign(SEM.H12.viewDice()); 
		// Bottom H22 block
		DoubleMatrix2D _H22 = A.viewPart(np,np,SEM.nPsystematics,SEM.nPsystematics);
		_H22.assign(SEM.H22);

		// and add H22 prior precision: H22 += C_BsysInv
		_H22.assign(SEM.C_BsysInv, cern.jet.math.Functions.plus);
		//EH.println("Hessian+H+C_BsysInv: "+A);
	  }
	// EH.println("Hessian: "+A);
  }

  /**
	 Compute effective "gradient" grad[], and effective Hessian A
	 of the Newton system
	 Requires space: (Np>=np)
	 A = DenseDoubleMatrix2D(Np,Np);
	 grad = float[Np]
	 On entry, expects
	 ptp   = partial^T.partial + W_B^T.W_B
	 i.e. mu = 1
  */
  void computeGradandHessian(float alogmu,DenseDoubleMatrix2D A)
  {
	// on input:
	// alogmu = log10 of the lagrange multiplier 


	int i, j,jmax,npm1;
	float s;
	float amu; 


	amu = fpow10(alogmu);

	assembleHessian(amu,A);  // <--  computes also grad[] here

  }

  /**
	 tofmu: Main objective function
  */
  float tofmu(float alogmu,float[] grad)
  {
	// function tofmu returnd the rms misfit of the response of a model constructed
	// using the given value of lagrange multiplier
	// 
	// on input:
	// alogmu = log10 of the lagrange multiplier 
	// (-99. if a zero lagrange multiplier needed)
	// on output
	// tofmu = r.m.s. misfit or 1.0e+10 if cholesky decomposition failed
	// pwk2[] contains the model which produces the misfit tofmu
	// dm[] contians the response of the model 
	// 

	boolean Lagrange_smoothing_only;
	DenseDoubleMatrix2D A;
	// float a(npp,npp+1);
	float[] grad_modified;
	float s, amu, part,diff;
	float alphaMax,alphaMin;
	double dM,dMP,ratio;
	float amax,tofmu;
	int i, j, Nptot,pass;
	float max_s,alpha,step_jump,MAX_STEP;
	// projected Newton variables
	boolean stepLimited;
	boolean flipToSD= false;
	float[] m_PN;
	float ep_PN;
	ModifiedCholeskyLeftDecomposition chol;
	float EP_MAX_ACTIVE_THRESHOLD = 0.1f;
	float MACH_EPS_F = 1e-7f;  // about macheps float (usually about 6e-8)

	nfor++;

	Nptot = np + SEM.nPsystematics;

	A = new DenseDoubleMatrix2D(Nptot,Nptot);

	computeGradandHessian(alogmu,A);
	grad_modified = grad.clone();


	// EH.printf("Lambda span: %s%n",showEigenvalueSpan(A));

	// Projected Newton method
	if (useBoundConstraints)
	  {
		m_PN = new float[np];

		float[] proj_test_dirn;		
		proj_test_dirn = grad;

		ep_PN =0.0f;
		float ep_PN_min = Float.MAX_VALUE;
		for (i=0;i<np;i++) 
		  {
			// Literal Bertsekas: m_PN[i] = min(max(pm[i] - grad[i], m_lower[i]),m_upper[i]); // truncate
			m_PN[i] = min(max(pm[i] - proj_test_dirn[i], m_lower[i]),m_upper[i]); // Based on full Hessian
			ep_PN_min =  min(abs(pm[i] - m_PN[i]),ep_PN_min);
			ep_PN += sqr(pm[i] - m_PN[i]);
		  }
		if (false)
		  {
			// Bertsekas idea: use distance to PN point to update epsilon_k
			ep_PN = (float)sqrt(ep_PN);
		  }
		else
		  {
			// My severe alternative: use shortest distance to constraint surface - will avoid premature "activation"
			// ep_PN -> 0 if method converges to a constraint plane
			ep_PN = max(1.00001f*ep_PN_min, 1e-5f);  // allow for some FP error
			ep_PN = min(ep_PN, EP_MAX_ACTIVE_THRESHOLD);
		  }
		// Accumulate active set
		for (i=0;i<np;i++) 
		  {
			act_PN[i] = 0;
			// Bertsekas inspired version
			if (((pm[i] - m_lower[i]) <= ep_PN) && (grad[i] >0)) act_PN[i] = -1;  // lower active 
			if (((m_upper[i] - pm[i]) <= ep_PN) && (grad[i] <0)) act_PN[i] = 1;   // upper active
		  }
		if (debug)
		  {
			EH.printf("A, pre-PN-mangling:%n");
			for (i=0;i<np;i++)
			  {
				for (j=0;j<np;j++) EH.printf("%15.8g ",A.getQuick(i,j));
				EH.println("");
			  }
			EH.printf("RHS-modified:%n");
			for (i=0;i<np;i++) EH.printf("%15.8g ",grad_modified[i]);
		    EH.println("");
		  }

		// Compute projected Hessian. Zrero off-diagonals
		for (i=0;i<np;i++) 
		  {
			if (act_PN[i] != 0)
			  {
				// active component. Diagonalise wrt these
				for (j=0;j<np;j++)  
				  {
					if (j==i) continue;  // don't touch diagonal!
					A.setQuick(i,j,0.0);
					A.setQuick(j,i,0.0);
				  }
				// A.setQuick(i,i,(grad[i] != 0f) ? abs(grad[i]) : 1.0); // could also be, say 1.0				
				// A.setQuick(i,i,1);				
			  }
		  } 
		if (verbosity > 2) 
		  {
			EH.printf("Active set (ep=%4g): ",ep_PN);
			for (i=0;i<np;i++) EH.printf("%2d ",act_PN[i]);
			EH.println("");
		  }
		// Create new appropriate effective gradient (technique for Marquardt is in later chunk)
		if (!useLMstrategy)
		  {
			for (i=0;i<np;i++) 
			  {								   			
				if (act_PN[i] == -1)
				  {
					// jamming up against lower bound
					grad_modified[i] = (pm[i] - m_lower[i]) * (float)A.getQuick(i,i); 
				  }
				else if (act_PN[i] == 1)
				  {
					// jamming up against upper bound
					grad_modified[i] = (pm[i] - m_upper[i]) * (float)A.getQuick(i,i); 
				  }
				else // head to the interior..									  
				  {
					grad_modified[i] = grad[i];
					// for (j=0;j<np;j++) RHS[i] += A.getQuick(i,j)*pm[j];
				  }
			  }
			if (modelSystematics) // systematic components unmodified
			  {
				for (i=np,j=0; j<SEM.nPsystematics; i++,j++) grad_modified[i] = grad[i];
			  }
		  }
	  }

	// Marquardtize if desired
	if (useLMstrategy)
	  {
		if (useBoundConstraints)
		  {
			// experiment. Reset lambda_LM if active set changes, to multiple of gershgorin bound
			// on inactive submatrix of A
			boolean changed = false;
			if (LMinfo._act_PN != null)
			  {
				for (i=0;i<np;i++) if ( act_PN[i] != LMinfo._act_PN[i]) {changed=true;break;} 
			  }
			else changed = true;
			if (changed)
			  {
				float maxA = 0f;
				for (i=0;i<np;i++) 
				  {
					if (act_PN[i]==0) 
					  {
						for(j=0;j<=i;j++) maxA = (float)max(maxA,abs(A.getQuick(i,j)));
					  }
				  }
				LMinfo.lambda_LM = 0.01f*maxA;
				if (verbosity > 2) EH.printf("LM: changed active set, changed lambda_LM: %-6g%n",LMinfo.lambda_LM);
			  }
		  }

		// Now the Hessian modification
		if (!useBoundConstraints)
		  {
			// works sort of OK also in bound case, but chunk below is better 
			for (i=0;i<Nptot;i++) 
			  {
				// leave quadratic problem for active constraints unchanged ?
				// if (act_PN[i] == 0) A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
				A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
				// don't adjust RHS?
				grad_modified[i] = grad[i];
			  }
		  }
		else
		  {
			// first, estimate good diagonal scale of inactive submatrix
			boolean rescaleActiveDiags = false;
			double Amin = Double.MAX_VALUE;
			for (i=0;i<np;i++) 
			  {
				if (act_PN[i]==0) Amin = min(Amin,A.getQuick(i,i));
			  }
			if (Amin!=Double.MAX_VALUE && Amin>0.0) {rescaleActiveDiags = true;}
			if (true)
			  {
				// Here: rhs adjusted so active solutions point at bounds with current Hessian
				//       then Hessian adjusted for all components using LM machinery
				// Marquardt trajectory will shrink in inactive plus active directions

				for (i=0;i<np;i++) 
				  {
					// scale down active components in Hessian....
					if ((act_PN[i] != 0) && rescaleActiveDiags) A.setQuick(i,i,Amin);

					// adjust also RHS to aim for bounds, given current modified Hessian
					if (act_PN[i] == 0) 
					  {
						grad_modified[i] = grad[i];
					  }
					else if (act_PN[i] == -1)
					  {
						// jamming up against lower bound: independent of lambda_LM
						grad_modified[i] = (pm[i] - m_lower[i]) * (float)A.getQuick(i,i); 
					  }
					else if (act_PN[i] == 1)
					  {
						// jamming up against upper bound: independent of lambda_LM
						grad_modified[i] = (pm[i] - m_upper[i]) * (float)A.getQuick(i,i); 
					  }
					// but then shrink accordingly, all components
					A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
				  }
				if (modelSystematics) // no bounds on systematics
				  { 
					for (i=np;i<Nptot;i++) A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
				  }
			  }
			else
			  {
				// operate LM machinery only on inactive subspace 
				for (i=0;i<np;i++) 
				  {
					// scale down active components in Hessian....
					if ((act_PN[i] != 0) && rescaleActiveDiags) A.setQuick(i,i,Amin);

					// adjust also RHS
					if (act_PN[i] == 0) 
					  {
						// here, we shrink/grow in inactive subspace
						A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
						grad_modified[i] = grad[i];
					  }
					else if (act_PN[i] == -1)
					  {
						// jamming up against lower bound: independent of lambda_LM
						grad_modified[i] = (pm[i] - m_lower[i]) * (float)A.getQuick(i,i); 
					  }
					else if (act_PN[i] == 1)
					  {
						// jamming up against upper bound: independent of lambda_LM
						grad_modified[i] = (pm[i] - m_upper[i]) * (float)A.getQuick(i,i); 
					  }
				  }
				if (modelSystematics) // no bounds on systematics
				  { 
					for (i=np;i<Nptot;i++) A.setQuick(i,i,A.getQuick(i,i) + LMinfo.lambda_LM * LMinfo.scales[i]);
				  }
			  }			  
		  }
		LMinfo.grad_modified = grad_modified;
		if ((LMinfo.H_modified==null) || (LMinfo.H_modified.rows() != Nptot ))
		  {
			LMinfo.H_modified = new DenseDoubleMatrix2D(Nptot,Nptot);
		  }
		LMinfo.H_modified.assign(A);
		LMinfo._act_PN = act_PN;
	  }

	// Gauss-Newton line search...
	// solve the linear system

    if (verbosity > 2) printGradient();

	if (debug)
	  {
		EH.printf("A:%n");
		for (i=0;i<np;i++)
		  {
			for (j=0;j<np;j++) EH.printf("%15.8g ",A.getQuick(i,j));
			EH.println("");
		  }
		EH.printf("grad_modified:%n");
		for (i=0;i<np;i++) EH.printf("%15.8g ",grad_modified[i]);
		EH.println("");
	  }

	if (!useSVD)
	  {
		// Standard cholesky
		if (false)
		  {
			// interesting check of Hessian condition number
			EigenvalueDecomposition eig = new EigenvalueDecomposition(A); 
			DoubleMatrix1D eR = eig.getRealEigenvalues();
			double minEig = abs(eR.getQuick(0));
			double maxEig = minEig;
			for (i=1;i<np;i++)
			  {
				minEig = min(minEig, abs(eR.getQuick(i)));
				maxEig = max(maxEig, abs(eR.getQuick(i)));
			  }
			EH.println("Total Bayes-Hessian condition number "+ maxEig/minEig);
		  }
		// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("Abefore.txt",A);
		chol = new ModifiedCholeskyLeftDecomposition(A); 
		// diagnoseMatrix(A); // FIXME remove!!!
		// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("Aafter.txt",A);
		// Solve for model change -delta_m = H^{-1}.grad
		chol.solve(pwk2,grad_modified);

		// Active block seems to get inaccurate solutions, maybe due to poor condition of inactive block
		// Correct -delta_m for active block:
		if (useBoundConstraints) correctActiveComponents(A);

		// stuff into absolute new model pwk2 = pm - (-delta_m)
		for (i=0;i<Nptot;i++) pwk2[i] = pm[i] - pwk2[i];
	  }
	else
	  {
		// Linear algebra check of y = A^{-1}.RHS
		DenseDoubleMatrix2D Ainv;
		DenseDoubleMatrix1D rhscpy = new DenseDoubleMatrix1D(Nptot);
		DoubleMatrix1D y = rhscpy.like();
		DenseDoubleMatrix2D Acpy = new DenseDoubleMatrix2D(Nptot,Nptot);
		for (i=0;i<Nptot;i++)
		  {
			for (j=0;j<Nptot;j++) Acpy.setQuick(i,j,A.getQuick(i,j));
		  }
		if (true)
		  {
			// SVD version with cutoff
			SingularValueDecomposition svd = new SingularValueDecomposition(Acpy);
			DoubleMatrix2D S = svd.getS();
			DoubleMatrix2D U = svd.getU();
			double sThresh,sMax;
			int nSVD = Nptot;
			sMax = S.getQuick(0,0);    // Largest eigenvalue
			sThresh = sMax  * (double)eps_SVD;
			// must have at least 1 eigenvalue..
			for (i=1;i<Nptot;i++)
			  {
				if (S.getQuick(i,i) < sThresh)
				  {
					nSVD=i;
					break;
				  }
			  }
			if (verbosity > 2) EH.println("SVD cutoff uses: "+nSVD+ " dominant eigenvalues");
			// form S^{-1}, adjusted for SVD cutoff
			for (i=0;i<nSVD;i++) S.setQuick(i,i,1.0/S.getQuick(i,i));
			for (i=nSVD;i<Nptot;i++) S.setQuick(i,i,0.0);
			// Linear algebraic matrix-vector multiplication; z = alpha * A * y + beta*z.
			// solve for increment only...
			for (i=0;i<Nptot;i++) rhscpy.setQuick(i, -grad_modified[i]);
			// zMult(y,z,alpha,beta,transposeA)
			// update is U^T.rhs
			DoubleMatrix1D z = U.zMult(rhscpy,null, 1.0 ,0.0, true);
			// S^{-1}.(U^T.rhs)
			for (i=0;i<Nptot;i++) z.setQuick(i,z.getQuick(i) * S.getQuick(i,i));
			DoubleMatrix1D ans = U.zMult(z,null, 1.0 ,0.0,false);
			for (i=0;i<Nptot;i++) pwk2[i] = -(float)ans.getQuick(i);
		  }
		else
		  {
			// Vanilla-style inverse, as check of ModifiedCholeskyLeftDecomposition style
			for (i=0;i<Nptot;i++) rhscpy.setQuick(i,grad_modified[i]);
			Algebra alg = new Algebra();
			Ainv = (cern.colt.matrix.impl.DenseDoubleMatrix2D)alg.inverse(Acpy);
			// Linear algebraic matrix-vector multiplication; z = alpha * A * y + beta*z.
			Ainv.zMult(rhscpy,y,1,0,false);
			for (i=0;i<Nptot;i++) pwk2[i] = (float)y.getQuick(i);
		  }

		// Active block seems to get inaccurate solutions, maybe due to poor condition of inactive block
		// Correct -delta_m for active block:
		if (useBoundConstraints) correctActiveComponents(A);

		// add back to current point 
		for (i=0;i<Nptot;i++) pwk2[i] = pm[i] - pwk2[i];
	  }

	// Newton point is now in pwk2().
	// Old model in pm(i). 

    if (debug)
      {
        EH.printf("Solution: ");
        for (i=0;i<np;i++) EH.printf("%15.8g ",pwk2[i]);
        EH.println("");
		EH.printf("Delta-pm: ");
        for (i=0;i<np;i++) EH.printf("%15.8g ",pwk2[i]-pm[i]);
        EH.println("");
      }

	// EH.warn("RESTORE MAX_STEP=0.3");	
	MAX_STEP = 0.3f; // 0.3
    if (useLMstrategy)
	  {
		// maximal step reduction, should occur only early on, and rarely
		reduceStep(MAX_STEP,pwk2,pm);
		s = getDirectionalSlope(pwk2,pm,grad,null,null); // NB overwrites m_prime[0]!!!
	  }
	else
	  {
		if (useBoundConstraints)
		  {
			// Reduce Newton step to smallest that leaves projection/truncation invariant
			alphaMax=0f;
			alphaMin=Float.MAX_VALUE;

			for (pass=0; pass<2; pass++)
			  {
				for (i=0;i<np;i++) 
				  {
					dMP = min(max(pwk2[i], m_lower[i]),m_upper[i]) - pm[i];
					dM = abs(pwk2[i]-pm[i]);
					ratio = (dM > 0) ? abs(dMP)/dM : 0.0;
					alphaMax = (float)max(alphaMax, ratio);
					if (ratio > MACH_EPS_F)  //about macheps float
					  {
						// parameter not-on-boundary change
						alphaMin = (float)min(alphaMin, ratio); 
					  }
				  }
				if (verbosity > 2) EH.println("alphaMax = "+alphaMax+" alphaMin = "+alphaMin);

				// alphaMax is now maximum reduction of Newton step that leaves Box-projection unchanged
				// if point is not on box surface, alphaMin>0
				if (pass==0) 
				  {
					alpha = (alphaMin != Float.MAX_VALUE) ? alphaMin : alphaMax; // Newton
				  }
				else alpha = alphaMin;         // SD
				// reduce Newton point to here
				for (i=0;i<Nptot;i++) 
				  {
					pwk2[i] = alpha * (pwk2[i]-pm[i]) + pm[i];
				  }
				if (verbosity > 3)
				  {
					EH.println("Confining PN step to box using scaling factor alpha = "+alpha);
					EH.printf("Rescaled delta_pm: ");
					for (i=0; i<Nptot; i++) EH.printf("%15.8g ",pwk2[i] - pm[i]); 
					EH.printf("%n");
				  }
				// Reduce step size further if Projected Newton point is still very far away
				for (alpha=1.0f,stepLimited=false;;)
				  {
					for (i=0,max_s=0f;i<np;i++)
					  {
						dMP = min(max(pwk2[i], m_lower[i]),m_upper[i]) - pm[i];
						step_jump = (float)abs(dMP);
						max_s = (float)max(max_s, step_jump);
					  }
					if (max_s > MAX_STEP)
					  {
						if (verbosity > 2) EH.println("Reducing Projected Newton step...max_s= "+max_s);
						alpha *= 0.5;
						// Reduced change
						for (i=0;i<Nptot;i++) pwk2[i] = alpha * (pwk2[i] - pm[i]) + pm[i];
						stepLimited = true;
					  }
					else break;
				  }
				if ((verbosity > 3) && stepLimited)
				  {
					EH.printf("Final step-limited delta_pm: ");
					for (i=0; i<Nptot; i++) EH.printf("%15.8g ",pwk2[i] - pm[i]); 
					EH.printf("%n");
				  }
				// Check directional derivative is still a descent
				s = getDirectionalSlope(pwk2,pm,grad,null,null);
				if (verbosity > 2) EH.println("Newton direction slope f'(0) = "+s);
				flipToSD = false;
				if (true)
				  {
					if (verbosity > 2) 
					  {
						EH.printf("Proposed change, pass=%d SD=%s: ",pass,flipToSD);
						for (i=0; i<np; i++) 
						  {
							EH.printf("%15.8g ",min(max((pwk2[i] - pm[i])/frac + pm[i], m_lower[i]),m_upper[i]) - pm[i]);
						  }
						EH.printf("%n");
					  }
				  }
				if (pass==0) 
				  {
					if (s > 0) 
					  {
						flipToSD = true; // Rotation to Newton direction failed for some (Bizarre?) reason. 
						if (verbosity > 1) EH.println("Newton direction BUG: uphill?? Switching to steepest descent.");
					  }
					else
					  {
						float[] m_prime = new float[Nptot];
						for (i=0; i<Nptot; i++) m_prime[i] = (pwk2[i] - pm[i])/frac + pm[i];
						for (i=0; i<np; i++) m_prime[i] = min(max(m_prime[i], m_lower[i]),m_upper[i]);
					
						s = getDirectionalSlope(m_prime,pm,grad,null,m_prime); // NB overwrites m_prime[0]!!!
						if (verbosity > 2) EH.printf("Projected-Newton (frac=%8g) direction slope f'(0) = %10.6g%n",frac,s);
						if (s > 0) 
						  {
							flipToSD = true; // Rotation to Projected Newton direction failed.
							if (verbosity > 1) EH.printf("Projected Newton direction BUG: uphill?? (Note |Delta_m| = %8.4g). " 
														 +"Switching to steepest descent.%n",m_prime[0]);
						  }
					  }
				  }
				if (flipToSD)
				  {
					// Flip to steepest descent & continue
					for (i=0;i<Nptot;i++) pwk2[i] = pm[i] - grad[i];
					continue;
				  }
				else
				  {
					break;  // out of pass loop
				  }
			  }
		  }

		if (!useBoundConstraints || (useBoundConstraints && flipToSD))
		  {
			// unbounded method, or bounded method with desperate flip to SD needs step limiting....
			// Reduce step size  if proposed point is very far away
			reduceStep(MAX_STEP,pwk2,pm);
		  }
	  }
	// Newton point is now in pwk2(). Arjimo backtrack by frac (won't happen in Marquardt)
	if (frac > 1.0f)
	  {
		part = 1.0f/frac;
		for (i=0;i<Nptot;i++)
		  {
			pwk2[i] = pm[i] + part*(pwk2[i] - pm[i]);
		  }
	  }

	if (useBoundConstraints)
	  {		
		// Projected line search target point 
		for (i=0;i<np;i++) pwk2[i] = min(max(pwk2[i], m_lower[i]),m_upper[i]);	
	  }

	if (debug)
	  {
		EH.printf("Current-point : ");
		for (i=0;i<Nptot;i++) EH.printf("%15.8g ",pm[i]);
		EH.println("");
		EH.printf("Proposed-new-point: ");
		for (i=0;i<Nptot;i++) EH.printf("%15.8g ",pwk2[i]);
		EH.println("");
	  }

	// calculate the model response and misfit
	long t = System.nanoTime();
	doForwardModel(pwk2);
	tNatives += (System.nanoTime()-t); 
	
	if (doNoiseEstimation) 
	  {
		s = (chisqMisFit()/sqr(global_error_scale));
	  }
	else if (global_error_scale > 0f)  // either doNoiseLevelCooling, or PBsampling, say
	  {	
		s = (chisqMisFit()/sqr(global_error_scale));
	  }
	else
	  {
		// here, sd[] as per files
		s = chisqMisFit(); // raw misfit using supplied sd[] values
	  }
	tofmu = (float)sqrt(s/nd);

	return tofmu;
  }

  // In the bound-constrained methods, the active block seems to get 
  // inaccurate solutions, maybe due to poor condition of inactive block
  // Write in the answers we know analytically 
  void correctActiveComponents(DenseDoubleMatrix2D A)
  {
	int i;
	double Aii;
	
	for (i=0;i<np;i++) 
	  {
		if (act_PN[i] == -1)
		  {
			// jamming up against lower bound
			Aii = A.getQuick(i,i);
			pwk2[i] = (float)((pm[i] - m_lower[i]) * (useLMstrategy ? (Aii/(Aii+LMinfo.lambda_LM * LMinfo.scales[i])) : 1.0));
		  }
		else if (act_PN[i] == 1)
		  {
			// jamming up against upper bound
			Aii = A.getQuick(i,i);
			pwk2[i] = (float)((pm[i] - m_upper[i]) * (useLMstrategy ? (Aii/(Aii+LMinfo.lambda_LM * LMinfo.scales[i])) : 1.0)); 
		  }
	  }
  }

  void printGradient()
  {
	int i,Nptot;

	Nptot = grad.length;
	EH.printf("-grad: ");
	for (i=0;i<Nptot;i++) EH.printf("%15.8g ",-grad[i]);
	EH.println("");
  }

  // ensure jump from m_old to m_new has no component jump larger than delta_m_max
  void reduceStep(float delta_m_max,float[] m_new,float[] m_old)
  {
	int i,Nptot;
	float max_s, step_jump;

	Nptot = m_new.length;
	for (i=0,max_s=0f;i<Nptot;i++)
	  {
		step_jump = (float)abs(m_new[i] - m_old[i])/(float)max(abs(m_old[i]),1.0);
		max_s = (float)max(max_s, step_jump);
	  }
	if (max_s > delta_m_max)
	  {
		if (verbosity > 2) EH.println("Reducing step...");
		max_s = delta_m_max/max_s;
		// Reduced change
		for (i=0;i<Nptot;i++) m_new[i] = max_s * (m_new[i] - m_old[i]) + m_old[i];
	  }
  }

  /* 
	 Computes directional derivative of vector pointing to m_new from m0, 
	 with respect to gradient g. Can return delta_m[0...np-1] as
	 delta_m = m_new - m0
	 and delta_m_norm[0] = ||delta_m||
	 if space is provided. 
	 If not, supply delta_m=null, delta_m_norm =null
	 and these results are discarded.
  */
  float getDirectionalSlope(float[] m_new,float[] m0, float[] g,float[] delta_m,float delta_m_norm[])
  {
	int i;
	float slope,norm;
	float[] dM;

	if (delta_m == null) 
	  {
		dM = new float[m0.length]; // allocate new space
	  }
	else 
	  {
		dM = delta_m;
	  }
	// assume here m0.length = m_new.length
	for (i=0; i< m0.length; i++) dM[i] = m_new[i] - m0[i];
	slope = 0.0f;
	for (i=0; i< m0.length; i++) slope += g[i]* dM[i];
	norm = (float)sqrt(anorm(dM));
	if (norm>0) slope /= norm;

	// return useful norm too
	if (delta_m_norm != null) delta_m_norm[0] = norm;

	return slope;
  }

  /**
	 Returns the ||_infty norm of the active part of the 
	 gradient applicable at the last tofmu() call. 
  */
  float getInfInactiveNormGrad()
  {
	int i=0;
	float normInf=0f;
	float[] x;

	if (useLMstrategy)
	  {
		x = LMinfo.grad_modified;  // Use modified form of gradient!
	  }
	else
	  {
		x = grad; // Use full gradient
	  }
	if (useBoundConstraints)
	  {
		for (i=0;i<np;i++) 
		  {								   			
			if (act_PN[i] == 0) normInf = max(normInf,abs(x[i]));
		  }
	  }
	else 
	  {
		for (i=0;i<np;i++) normInf = max(normInf,abs(x[i])); 
	  }
	return normInf;
  }


  float anorm(float[] x)
  {
	int i;
	float s = 0.0f;
	for (i=0;i<x.length;i++) s+= x[i]*x[i];
	return s;
  }



  /**
	 chisqMisFitfor fixed noise. Useful utility.
  */
  float chisqMisFit()
  {
	int k;
	float chisq, v;
	
	chisq=0.0f;
	for (k = 0; k < nd; k++) 
	  {
		v= (d[k] - dm[k])/sd[k];
		chisq += v*v;
	  }
	return chisq;
  }


  /**
	 chisqRoughness: computes additional objective terms from l2 roughness terms
  */
  float chisqRoughness(float[] p_mod)
  {
	// local vars
	int  i,npM1,npM2;
	float  y,chisq_ruf;
	boolean roughnessRelativeToPrior = true;

	// Roughness bits:
	// iruf = 0 for no roughnening (i.e. min. norm models) (added in v1.5)
	// = 1 for 1st derivative roughness
	// 2 for 2nd derivative roughness
	// /model/ np = number of model parameters
	// NB objmat() sets npen = np

	objmat();

	// del has at worst a diagonal and 2 upper diagonals. 
	npM1 = np-1;
	npM2 = np-2;
	chisq_ruf=0.0f;
	if (roughnessRelativeToPrior)
	  {
		// Roughness is ||del.(m-m_prior)||^2
		for (i=0; i<np; i++)
		  {
			// y = ith component of del.m ...
			y = (float)del.getQuick(i,i) * (p_mod[i] - premod[i]);
			if (i<(npM1)) 
			  {
				y += del.getQuick(i,i+1)* (p_mod[i+1]- premod[i+1]);
				if (i<(npM2)) 
				  {
					y += del.getQuick(i,i+2) * (p_mod[i+2] - premod[i+2]);
				  }
			  }
			// here, the y norm of ||del.(m-m_prior)|| is assembled...
			chisq_ruf += y*y;
		  }
	  }
	else
	  {
		// Roughness is ||del.m||^2
		for (i=0; i<np; i++)
		  {
			// y = ith component of del.m ...
			y = (float)del.getQuick(i,i)*p_mod[i];
			if (i<(npM1)) 
			  {
				y += del.getQuick(i,i+1)*p_mod[i+1];
				if (i<(npM2)) 
				  {
					y += del.getQuick(i,i+2)*p_mod[i+2];
				  }
			  }
			// here, the y norm of ||del.m|| is assembled...
			chisq_ruf += y*y;
		  }
	  }
	return chisq_ruf;
  }

  /**
	 chisqPrej: computes additional prejudice 
	 objective terms, based on parameters in p_mod[0...Nptot-1]
  */
  float chisqPrej(float log10_mu,float[] p_mod)
  {
	int i;
	float x,sum;

	// Prejudice bits
	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  {
		updatePrewtsForExpStyleSmoothing(pow(10.0,log10_mu));
	  }

	sum = 0.0f;
	for(i=0; i<np; i++)
	  {
		x = prewts[i] * (p_mod[i]-premod[i]);
		sum += x*x;
	  }
	if (modelSystematics)
	  {
		int ii,j,jj;
		float sum_sys=0f;
		DoubleMatrix2D _C_BsysInv = SEM.C_BsysInv;
		for (i=np,ii=0; ii<SEM.nPsystematics; i++,ii++)
		  {
			for (jj=0,x=0,j=np; jj<SEM.nPsystematics; j++,jj++)
			  {
				x += _C_BsysInv.getQuick(ii,jj) * p_mod[j];
			  }
			sum_sys += x * p_mod[i];
		  }
		sum += sum_sys;
		// EH.printf("Prior sys loglik: %10g%n",sum_sys);
	  }
	return sum;
  }


  float makjtj()
  {
	// Constructs matrices that depend on the Jacobian (wjtwj).
	
	float sum,tol0,s,noise_est_nd=1f;
	double S;
	int i, j,k;
	
	//-----------------------------------------------------------------------
	// calculate the matrix of partials and model response
	if (debug)
	  {
		if (!newScripps)
		  {
			for (i=0;i<6;i++)
			  {
				EH.printf("dp (java): ");
				for(j=0; j<min(nd,6); j++) EH.printf("%10.6g ",dp.getQuick(i,j));
				EH.printf("\n"); 
			  }		 
		  }
		EH.printf("pm: "); 
		for (j=0; j<np; j++) EH.printf("%10.6g ",pm[j]);
		EH.printf("%n");
	  }

	long t = System.nanoTime();
	doForwardModelPlusJac();
	tNatives += (System.nanoTime()-t); 

	if (false)
	  {
		EH.printf("dm: "); 
		for (j=0; j<nd; j++) EH.printf("%10.6g ",dm[j]);
		EH.printf("%n");

		PrintWriter out;
		String filename="X.txt";
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		EH.printf("X: (slow x fast) %d x %d%n",np,nd); 
		for (i=0;i<np;i++)
		  {
			for (j=0; j<nd; j++) out.printf("%11.5g ",WJ_trans.getQuick(i,j));
			out.printf("%n");
		  }
		out.close();
		EH.println("Wrote "+filename);
		System.exit(0);
	  }


	// calc misfit vector and misfit
	// form dwk1 = W.(d-F(m))
	sum=0.0f;
	if (doNoiseEstimation)
	  {
		// Need replacement as if W -> W/(sqrt(chisq_data)/n_d)
		for (i=0;i<nd;i++)
		  {
			s = (d[i] - dm[i])/sd[i];
			dwk1[i] = s;
			sum += s*s;
		  }
		noise_est_nd = (float)sqrt(sum/nd);
		// global_error_scale = (float)sqrt(sum/nd);
		// this is the EB estimate
		global_error_scale = (float)sqrt(sum/DOF(nd,np));
		if (verbosity > 1) EH.printf("Noise estimation: noise_scale=%10.6g%n",global_error_scale);
		// and linear theory requires (even for EB) adjustment of W by W -> W/(sqrt(chisq_data)/n_d)
		for (i=0; i<nd; i++)
		  {
			dwk1[i] /= noise_est_nd; // global_error_scale; // first rescale
		  }
		//tol0 = global_error_scale;
		tol0 = (float)sqrt(DOF(nd,np)/nd); // = sqrt(sum/(nd*global_error_scale^2))
	  }
	else if (global_error_scale > 0f)
	  {
		// here, global_error_scale is "fixed & known". eg. Noise cooling, PB sampling
		for (i=0; i<nd; i++)
		  {
			dwk1[i] = (d[i] - dm[i])/(global_error_scale * sd[i]);
			sum += dwk1[i]*dwk1[i];
		  }
		tol0 = (float)sqrt(sum/nd);
	  }
	else
	  {
		for (i=0; i<nd; i++)
		  {
			dwk1[i] = (d[i] - dm[i])/sd[i];
			sum += dwk1[i]*dwk1[i];
		  }
		tol0 = (float)sqrt(sum/nd);
	  }

	if (debug)
	  {
		EH.printf("WJ:%n");
		for (i=0; i<nd; i++)
		  {
			for (j=0;j<np;j++)
			  {
				EH.printf("%15.8g ",WJ_trans.getQuick(j,i));
			  }
			EH.println("");
		  }
	  }	

	// weight the jacobian matrix: form WJ_trans = J^T.C_d^{-1/2}/sigma_n
	if (doNoiseEstimation)
	  {
		for (j=0;j<np;j++)
		  {
			for (i=0; i<nd; i++)
			  {
				// WJ_trans.mult(j,i,1.0f/(global_error_scale * sd[i]));
				WJ_trans.mult(j,i,1.0f/(noise_est_nd * sd[i]));
			  }
		  }
	  }
	else if (global_error_scale > 0f)
	  {
		for (j=0;j<np;j++)
		  {
			for (i=0; i<nd; i++)
			  {
				WJ_trans.mult(j,i,1.0f/(global_error_scale * sd[i]));
			  }
		  }
	  }
	else
	  {
		for (j=0;j<np;j++)
		  {
			for (i=0; i<nd; i++)
			  {
				WJ_trans.mult(j,i,1.0f/sd[i]);
			  }
		  }
	  }
	
	// form wjtwj = (WJ)^T.(WJ)= WJ_trans. WJ_trans^T

	for (i=0; i<np; i++)
	  {
		for (j=i;j<np;j++)   // use symmetry
		  {
			for (k=0,S=0.0; k<nd; k++) S += WJ_trans.getQuick(i,k) * WJ_trans.getQuick(j,k);
			wjtwj.setQuick(i,j,S);
			wjtwj.setQuick(j,i,S);
		  }
	  }

	if (false)
 	  {
		EH.printf("wjtwj:%n");
		for (i=0; i<np; i++)
		  {
			for (j=0;j<np;j++)
			  {
				EH.printf("%15.8g ",wjtwj.getQuick(i,j));
			  }
			EH.println("");
		  }
	  }

	if ((haveLargestEigenvalueWJTWJ == false) && (!makeSynthetic))
	  {
		// estimate largest eigenvalue of (WJ)^T. WJ (can skip if just making synthetic)
		largestEigenvalueWJTWJ = getLargestEigenvalue(wjtwj);
		haveLargestEigenvalueWJTWJ = true;
	  }


	// Compute WJdeltaD = (WJ)^T. W(d-F(m))
	mult(np,nd,WJ_trans,dwk1,WJdeltaD);

	if (modelSystematics)
	  {
		// Form X_N^T.W^T . W(d-F(m)). By here, dwk1 = W(d-F(m)) = C_d^{-1/2}.(y-f(m))/sigma_n
		// (bottom block-vector of gradient component from data-mismatch)
		DenseDoubleMatrix2D _X = SEM.X;
		
		for (i=0;i<SEM.nPsystematics; i++)
		  {
			for (j=0,S=0.0; j<nd; j++) S += _X.getQuick(j,i) * dwk1[j]/sd[j];
			SEM.X_NdeltaD[i] = (float)S;
		  }
		if (doNoiseEstimation)
		  {
			for (i=0;i<SEM.nPsystematics; i++) SEM.X_NdeltaD[i] /= noise_est_nd;
		  }
		else if (global_error_scale > 0f)
		  {
			for (i=0;i<SEM.nPsystematics; i++) SEM.X_NdeltaD[i] /= global_error_scale;
		  }


		// Also H12 and H22 blocks:
		// H12:
		for (i=0;i<np;i++)
		  {
			for (j=0; j<SEM.nPsystematics; j++)
			  {
				for (k=0,S=0.0;k<nd;k++) S += WJ_trans.getQuick(i,k) * _X.getQuick(k,j)/sd[k];
				SEM.H12.setQuick(i,j,S);
			  }
		  }	
		if (doNoiseEstimation)
		  {
			SEM.H12.assign(Functions.div(noise_est_nd));
		  }
		else if (global_error_scale > 0f)
		  {
			SEM.H12.assign(Functions.div(global_error_scale));
		  }	

		// Form & weight the H22 block.
		// Form X_N^T.C_d^{-1}.X_N/sigma_n^2. Use symmetry.
		for (i=0;i<SEM.nPsystematics;i++)
		  {
			for (j=i; j<SEM.nPsystematics; j++)
			  {
				for (k=0,S=0.0;k<nd;k++) S += _X.getQuick(k,i) * _X.getQuick(k,j)/(sd[k]*sd[k]);
				SEM.H22.setQuick(i,j,S);
				if (i!=j) SEM.H22.setQuick(j,i,S);
			  }
		  }
		if (doNoiseEstimation)
		  {
			SEM.H22.assign(Functions.div(noise_est_nd * noise_est_nd));
		  }		
		else if (global_error_scale > 0f)
		  {
			SEM.H22.assign(Functions.div(global_error_scale * global_error_scale));
		  }
	  }

	return tol0;
  }



  /**
	 Scan solution to see if it resembles an existing one 
  */
  boolean foundSimilarExistingMode(float negLogPost)
  {
	int i,j;
	float[] m1;
	float pmu1;
	boolean isSimilar;
	LocalMode mi;

	if (optima.size() > 0)
	  {
		float MU_OBJ = 0.1f;   // epsilon for "matching objective values"
		float MU_EPS = 0.1f;   // epsilon for "matching mu values"
		float MP_EPS = 0.1f;   // epsilon for "matching paramater values"
		for (i=0;i<optima.size();i++)
		  {
			mi= optima.get(i);
			m1 = mi.pmLoc;
			pmu1 = mi.pmuLoc;
			isSimilar = true;
			// test objective
			if (abs(negLogPost - mi.negLogPost) > MU_OBJ) 
			  {
				isSimilar = false;   // sufficiently different
				break;
			  }
			// test mu
			if (abs(pmu1 - pmu) > MU_EPS) 
			  {
				isSimilar = false;  // sufficiently different 
				break;
			  }
			// test parameters
			for (j=0;j<np;j++)
			  {
				if (abs(pm[j] - m1[j]) > MP_EPS) 
				  {
					isSimilar = false;   // sufficiently different
					break;
				  }
			  }
			// if we got here, it looks like a duplicate
			if (isSimilar) 
			  {
				// found a very similar looking mode.
				return true;
			  }
		  }
	  }
	return false;
  }

  /**
	 wrappers for each supported forward model
  */
  void doForwardModelPlusJac()
  {
	if (newScripps)
	  {
		SMI.OM.forwardModelPlusJac(pm);
		if (false)
		  { 
			PrintWriter out;
			String filename="dm.txt";
			try
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			  } 
			catch (IOException ioe) 
			  {
				throw new RuntimeException("couldn't open " + filename+" for output");
			  }
			int i;
			for (i=0;i<nd;i++) out.printf("%15.8g   %15.8g%n",d[i],dm[i]);
			out.close();			
		  }
	  }
	else
	  {
		ScrippsLegacyModelInterface.fordiv(np,nd,pm,dp.elements,dm,WJ_trans,d);
	  }
	addSystematicsIfRequired(pm);
  }

  void doForwardModel(float[] model)
  {
	if (newScripps)
	  {
		SMI.OM.forwardModel(model);
	  }
	else
	  {
		ScrippsLegacyModelInterface.formod(np,nd,model,dp.elements,dm,d);
	  }
	addSystematicsIfRequired(model);
  }

  // Add Linear regression term to fwd model for systematic effects if needed
  // Note Jacobian associated with this is kept separately!
  void addSystematicsIfRequired(float[] model)
  {
	if (modelSystematics)
	  {
		// take copy of actual CSEM model in dm_true[]
		System.arraycopy(dm,0,dm_true,0,dm.length);

		// Add on general effect of systematics into dm[]
		int i,j,jj;
		double u;
		DenseDoubleMatrix2D _X = SEM.X;
		for (i=0;i<nd;i++)
		  {
			// for(j=0,jj=np,u=0.0; j<SEM.nPsystematics; j++,jj++) u += _X.getQuick(i,j) * SEM.m_N[j];
			for(j=0,jj=np,u=0.0; j<SEM.nPsystematics; j++,jj++) u += _X.getQuick(i,j) * model[jj];
			dm[i] += u;
		  }
	  }
  }

  /**
	 Clean up list of multistart solutions. Eliminate duplicates etc
  */
  void cleanUpMultistartSolutions(String extraText)
  {
	int i,j,ii,jj,k,Np_eff,Nopt,Nptot;
	int NSTEPS = 25;  
	double frac,step,val,logDetInvPosterior,logDetInvPrior;
	LocalMode mi,mj;
	float m1[];
	float m2[];
	float[] xvals;
	float[] fvals;
	float f,tmp,dist,max,Cpii;
	float pmu1,pmu2,mueff;
	boolean remove=false;
	boolean printProfiles = true;

	if (verbosity > 1) EH.printf("Sorting %d local optima%n",optima.size());

	Nptot = np + SEM.nPsystematics;

	Nopt = optima.size();
	if (Nopt > 1)
	  {
		Collections.sort(optima);
		// loop over solutions, eliminating duplicates
		// form array of mode-pairs, with estimated "distances" between them
		// use, for distance, std-dev-normalised Euclidean distance between components of mode MAP points
		ArrayList<ModePair>Mp = new ArrayList();
		for (i=0; i<Nopt-1; i++)
		  {
			mi= optima.get(i);
			m1 = mi.pmLoc;
			for (j=i+1; j<Nopt; j++)
			  {
				mj = optima.get(j);
				m2 = mj.pmLoc;
				for (k=0,f=0.0f; k<Nptot; k++) 
				  {
					tmp = sqr(mi.stdv[k]) + sqr(mj.stdv[k]);
					f += sqr(m1[k]-m2[k])/tmp;
				  }
				Mp.add(new ModePair(i,j,f));
			  }
		  }
		// sort so "close--seeming" mode-pairs are early on in list of pairs
		Collections.sort(Mp);
		// now run over pairs, marking duplicate modes
		xvals = new float[NSTEPS];
		fvals = new float[NSTEPS];
		step=1.0/(NSTEPS-1);
		boolean[] removeMode = new boolean[Nopt];
		for (k=0; k < Mp.size(); k++)
		  {
			ModePair mp = Mp.get(k);
			if (removeMode[mp.i] || removeMode[mp.j]) continue;
			i = mp.i; 
			mi= optima.get(i);
			m1 = mi.pmLoc;
			pmu1 = mi.pmuLoc;
			j = mp.j;
			mj = optima.get(j);
			m2 = mj.pmLoc;
			pmu2 = mj.pmuLoc;
			if (verbosity>2) EH.printf("Examining mode-pair %d %d, mode distance = %g%n",i,j,mp.dist);
			// Form discrete profile of NSTEPS
			max = -Float.MAX_VALUE;
			for (ii=0,frac=0.0;ii<NSTEPS;ii++,frac+=step)
			  {
				// interpolate model vector
				for (jj=0;jj < Nptot; jj++)
				  {
					//val = (1.0-frac)*m1[jj] + frac*m2[jj]; // linear in Math.log(rho) space
					val = log10( (1.0-frac)*pow(10.0,m1[jj]) + frac * pow(10.0,m2[jj]));  // linear in rho space
					pm[jj] = (float)val;
				  }
				mueff = (float)((1.0-frac)*pmu1 + frac*pmu2);

				// get function				
				f = twiceNegLogPosterior(mueff);
				// store result
				xvals[ii] = (float)frac;
				fvals[ii] = f;
				max = max(max,f);
				// early breakout if desired ?
				// if (((f-fvals[0]) > MODE_SEPARATION_EPSILON) && ((max-f) > MODE_SEPARATION_EPSILON)) break; 
				// would need to clean up rest of fvals[] in above break.....
			  }

			// Possibilities for failure
				
			// (1) We can have instances of two modes that ought to be the same, differing by only
			// floating point errors.
			if (((max-fvals[0] < MODE_SEPARATION_EPSILON) && (max-fvals[NSTEPS-1] < MODE_SEPARATION_EPSILON))
				|| (max < (fvals[NSTEPS-1]+MODE_SEPARATION_EPSILON)))
			  {
				// One of the modes looks spurious. j is worse than i by sorting, so 
				// eliminate mode j
				removeMode[j] = true;
				if (verbosity > 2) EH.println("Removed duplicate mode "+j);
			  }
			// else... some kind of significant intermediate energy barrier. Retain both.
			if (printProfiles)
			  {
				String filename;
				PrintWriter out;
				filename = "ModeRemovalTestProfile."+ (removeMode[j]  ? "removed.":"retained.")+k+"."+i+"-"+j;
				if (extraText != null) filename = filename.concat(extraText);
				filename = filename.concat(".txt");
				try
				  {
					out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
				  } 
				catch (IOException ioe) 
				  {
					throw new RuntimeException("couldn't open " + filename+" for output");
				  }
				for (ii=0;ii<xvals.length;ii++) out.printf("%15.8g %15.8g%n",xvals[ii],fvals[ii]);
				out.close();
			  }
		  }
		// discard duplicate modes from optima (reverse order)
		for (k=Nopt-1; k>=0; k--) if (removeMode[k]) optima.remove(k);
		removeMode = null;
		Mp = null;
	  }
	
		

	// Get logMML of surviving modes..
	
	cern.colt.matrix.linalg.Algebra alg = new cern.colt.matrix.linalg.Algebra();
	Np_eff = doNoiseEstimation ? (np+1) : np;

	if (modelSystematics) Np_eff += SEM.nPsystematics;

	for (i=0; i<optima.size(); i++)
	  {
		mi = optima.get(i);
		// wjtwj = (cern.colt.matrix.impl.DenseDoubleMatrix2D) mi.wjtwj;		  
		// wjtwj.assign(mi.H.viewPart(0,0,np,np));
		H = mi.H;		
		wjtwj = mi.wjtwj;		
		CholeskyDecomposition chol = new CholeskyDecomposition(mi.H);
		DoubleMatrix2D L = chol.getL();
		// determinant of inverse posterior covariance
		for (j=0,logDetInvPosterior=0.0; j<Nptot; j++) 
		  {
			logDetInvPosterior += Math.log(L.getQuick(j,j));
		  }
		// determinant of inverse prior CP^{-1} = W_B^T.W_B
		// Assumption here is that smoothing is *zero* (del=0)
		if (doBayesianAutoSmoothing && !doEmpBayes)
		  {
			logDetInvPrior=0.0; // it's not really, but the term is already absorbed into negLogPost
		  }
		else
		  {
			for (j=0,logDetInvPrior=0.0; j<np; j++) 
			  {
				Cpii= prewts[j];
				if (Cpii <= 0f) Cpii = DEFAULT_SAFE_INV_RHO_PRIOR_STD_DEV;
				logDetInvPrior += Math.log(Cpii * Cpii);
			  }
			if (modelSystematics) logDetInvPrior += SEM.getLogDetC_BsysInv();
		  }
		if (doWishart)
		  {
			// rough'n'ready approx
			double wishartBit = 0.5*(Wishart_nu+1.0-np)*Math.log(1.0 + 2.0*mi.negLogPost/Wishart_nu);
			mi.logMMLapprox = (float) (wishartBit
									   + logDetInvPosterior/2.0f - logDetInvPrior/2.0f
									   + Np_eff*MathUtil.HALF_LOG_TWOPI);
			// Hass empirical Bayes...
			// mi.logMMLapprox = (float) (0.5*np*(Math.log(2.0*mi.negLogPost)));
		  }
		else if (doShrinkage)
		  {
			// recompute regression based on shrunken covariance
			computeShrinkageEstimates(i);
			// adjust using prior bits too
			mi.logMMLapprox += logDetInvPrior/2.0f + Np_eff*MathUtil.HALF_LOG_TWOPI;
		  }
		else
		  {
			mi.logMMLapprox = (float) (mi.negLogPost + logDetInvPosterior/2.0f - logDetInvPrior/2.0f
									   + Np_eff*MathUtil.HALF_LOG_TWOPI);
		  }
	  }

	if (verbosity > 1) 
	  {
		EH.println("Sorted list of local optima and parameters:");
		EH.printf("Mode NegLogPost RMS_misfit -logMML Parameters: %n");
		for (i=0;i<optima.size();i++)
		  {
			mi = optima.get(i);
			EH.printf("%4d %-10.6g %-10.6g %-10.6g ",(i+1),mi.negLogPost, mi.rmsMisfit, mi.logMMLapprox);
			for(j=0;j<Nptot;j++) EH.printf("%-10.4g ",mi.pmLoc[j]);
			EH.println("");
		  }
	  }
	// dump uncertainty files
	for (i=0; i<optima.size(); i++)
	  {
		ModeCovarianceAndStdDev modeCov;

		mi = optima.get(i);
		// wjtwj = (cern.colt.matrix.impl.DenseDoubleMatrix2D) mi.wjtwj;
		H = mi.H;
		wjtwj = mi.wjtwj;
		System.arraycopy(mi.pmLoc,0,pm,0,Nptot);
		modeCov = computeUncertainty();
		writeUncertaintyFiles(extraText,true,true,i+1,modeCov);
		if (i==0)
		  {
			// copy best files
			if (copyFile(current_profile_filename,"profile_final.txt"))
			  {
				if (verbosity > 2) EH.println("File "+current_profile_filename +" copied to profile_final.txt");
			  }
			if (copyFile(mi.resp_filename,"resp_final"))
			  {
				if (verbosity > 2) EH.println("File "+mi.resp_filename +" copied to resp_final");	
			  }
			if (copyFile(mi.iter_filename,"iter_final"))
			  {
				if (verbosity > 2) EH.println("File "+mi.iter_filename +" copied to iter_final");	
			  }
		  }
	  }
	if (doPathIntegralTests)
	  {
		// Try to connect optima using path-integral tricks....
		ModeLinkFinder MLF = new ModeLinkFinder(this);
		for (i=0; i<optima.size(); i++)
		  {
			mi = optima.get(i);
			for (j=i+1; j<optima.size(); j++)
			  {
				mj = optima.get(j);
				if (verbosity>1) EH.println("Finding mode-links, modes %d..%d",i+1,j+1);
				MLF.findModePairConnection(mi.pmLoc,mj.pmLoc,"ModePair."+(i+1)+"."+(j+1)); 
			  }
		  }
	  }


	// point at first mode for safety hereon
	mi = optima.get(0);
	System.arraycopy(mi.pmLoc,0,pm,0,np);
	// wjtwj = (cern.colt.matrix.impl.DenseDoubleMatrix2D) mi.wjtwj;
	H = mi.H;
	wjtwj = mi.wjtwj;
  }

  void computeShrinkageEstimates(int optima_index)
  {
	int i,j,k;
	double S_ij,s,logDetChat,logDetH,logDetHsh,sumSQ;
	DoubleMatrix2D Hsh;
	LocalMode mi;

	mi = optima.get(optima_index);

	// this is a selection of the required things. Makes no sense if 
	// doNoiseEstimation=true, or global_error_scale>0, or doSystematics=true

	// set model
	System.arraycopy(mi.pmLoc,0,pm,0,np);	// pm <- pmLoc 
	// get Jacobean etc
	long t = System.nanoTime();
	doForwardModelPlusJac();
	tNatives += (System.nanoTime()-t); 

	// scaled residuals by initial diagonal covariance guess
	for (i=0,sumSQ=0.0; i<nd; i++) 
	  {
		dwk1[i] = (d[i] - dm[i])/sd[i];
		// test robustness to noise: if alpha->1, logMML below insenstitive to noise :-)
		// dwk1[i] += .7*RD.nextGaussian();
		sumSQ += dwk1[i] * dwk1[i];
	  }
	// form Shrinkage:
	DenseDoubleMatrix2D Chat = new DenseDoubleMatrix2D(nd,nd);
	for(i=0;i<nd;i++)
	  {
		for(j=0;j<=i;j++)
		  {
			S_ij = dwk1[i] * dwk1[j];
			// My formula; with dropped nd
			s = shrinkage_alpha * (nd/sumSQ)*(S_ij) + ((1.0-shrinkage_alpha) * ((i==j) ? 1.0 : 0));
			// Leung
			// s = shrinkage_alpha * S_ij+ ((1.0-shrinkage_alpha) * (sumSQ/nd) * ((i==j) ? 1.0 : 0));
			Chat.setQuick(i,j,s);
			if (i!=j) Chat.setQuick(j,i,s);  
		  }
	  }
	// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("Chat.txt",Chat);

	// weight the Jacobian matrix: form WJ_trans = J^T.C_d^{-1/2}
	for (j=0;j<np;j++)
	  {
		for (i=0; i<nd; i++)
		  {
			WJ_trans.mult(j,i,1.0f/sd[i]);
		  }
	  }
	// fill scaled sensitivity matrix X = WJ
	DenseDoubleMatrix2D X = new DenseDoubleMatrix2D(nd,np);
	for (i=0; i<nd; i++)
	  {
		for (j=0;j<np;j++) X.setQuick(i,j,WJ_trans.getQuick(j,i));
	  }
	// Get quantities needed for MML based on linearisation
	// (1) log|Chat|
	CholeskyDecomposition chol = new CholeskyDecomposition(Chat); 
	DoubleMatrix2D L = chol.getL();

	// determinant of "prior" covariance
	for (i=0,logDetChat=0.0; i<nd; i++) logDetChat += Math.log(L.getQuick(i,i));
	
	// (2) residuals SS term from new covariance: sumSQ = dy^T.Chat^{-1}.dy
	DenseDoubleMatrix2D dy = new DenseDoubleMatrix2D(nd,1);
	for (i=0; i<nd; i++) dy.setQuick(i,0,dwk1[i]);
	DoubleMatrix2D ChatInvDy = chol.solve(dy);
	for (i=0,sumSQ=0.0; i<nd; i++) sumSQ += dwk1[i] * ChatInvDy.getQuick(i,0);

	// (3) Hessian contribution from   log|X^T.Chat^{-1}.X|
	DoubleMatrix2D Xcpy = X.like().assign(X);
	DoubleMatrix2D ChatInvX = chol.solve(X);
	// form Hsh = X^T.Chat^{-1}.X = X^T.ChatInvX
	// Linear algebraic matrix-matrix multiplication; C = alpha * A x B + beta*C.
	Hsh = Xcpy.zMult(ChatInvX,null ,1.0,0.0,true,false);

	// au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("H_sh.txt",Hsh);

	CholeskyDecomposition chol_Hsh = new CholeskyDecomposition(Hsh); 
	DoubleMatrix2D L_Hsh = chol_Hsh.getL();
	for (i=0,logDetHsh=0.0; i<np; i++) logDetHsh += Math.log(L_Hsh.getQuick(i,i));

	// set overall mode evidence, or MML (approx)
	mi.negLogPost = (float) sumSQ/2.0f;
	mi.logMMLapprox = (float) (mi.negLogPost + logDetHsh/2.0 + logDetChat/2.0);
	// prior contributions added outside this method
  }
 


  /**
	 Computes c = a^T.b
	 for the system of sizes: a[m][n],   b[m], c[n]
  */
  void trmult(int m,int n,DenseFloatMatrix2D a,float[] b,float[] c)
  {
	int i, k;
	float ci;

	for (i=0;i<n;i++)
	  {
		ci = 0.0f;
		for (k=0; k<m; k++) ci += a.getQuick(k,i)*b[k];
		c[i] = ci;
	  }
  }

  /**
	 Computes c = a.b
	 with c[n], a[n][m], b[m]
  */
  void mult(int n,int m,DenseFloatMatrix2D a,float[] b,float[] c)
  {
	// Computes c = a.b
	// with c[n], a[n][m], b[m]
	int i, j;
	float ci;
	for (i=0;i<n;i++)
	  {
		ci = 0.0f;
		for (j=0;j<m;j++) ci += a.getQuick(i,j)*b[j];
		c[i] = ci;
	  }
  }
  double getLargestEigenvalue(DenseDoubleMatrix2D A)
  {
	int i,N;
	double lMax = -Double.MAX_VALUE;

	EigenvalueDecomposition eig = new EigenvalueDecomposition(A); 
	DoubleMatrix1D eR = eig.getRealEigenvalues();
	N = A.rows();
	for (i=0;i<N;i++) lMax = max(lMax, abs(eR.getQuick(i)));
	eig = null;
	eR=null;
	return lMax;
  }

  String showEigenvalueSpan(DenseDoubleMatrix2D A)
  {
	int i,N;
	double lMax = -Double.MAX_VALUE;
	double lMin = Double.MAX_VALUE;

	EigenvalueDecomposition eig = new EigenvalueDecomposition(A); 
	DoubleMatrix1D eR = eig.getRealEigenvalues();
	N = A.rows();
	for (i=0;i<N;i++) 
	  {
		lMax = max(lMax, abs(eR.getQuick(i)));
		lMin = min(lMin, abs(eR.getQuick(i)));
	  }
	
	eig = null;
	eR=null;
	return new String(lMax+" "+lMin);
  }

  /* initialisation call: estimates largest eigenvalue of PTP */

  void computeLargestEigenValuePTP()
  {
	int i,npen;

	if (haveLargestEigenvaluePTP == false)
	  {
		// make partials matrix:
		npen = np; // as per objmat...
		objmat();

		// form ptp = del^T.del
		// A.zmult: C = alpha * A x B + beta*C.

		del.zMult(del,ptp,1.0,0.0,true,false);

		// estimate largest eigenvalue of del^T.del
		largestEigenvaluePTP = -Double.MAX_VALUE;
		EigenvalueDecomposition eig = new EigenvalueDecomposition(ptp); 
		DoubleMatrix1D eR = eig.getRealEigenvalues();
		for (i=0;i<np;i++)
		  {
			largestEigenvaluePTP = max(largestEigenvaluePTP,abs(eR.getQuick(i)));
		  }
		// emergency care..
		if (largestEigenvaluePTP <= 0)
		  {
			// compare to prior...
			for (i=0;i<np;i++) 
			  {
				largestEigenvaluePTP = max(largestEigenvaluePTP,prewts[i]*prewts[i]);
			  }
			// if still in trouble...
			if (largestEigenvaluePTP <= 0) largestEigenvaluePTP = 1.0;
		  }
		haveLargestEigenvaluePTP = true;
		eig = null;
		eR = null;		
	  }
  }

  void makptp(float log10_mu)
  {
	// JG: corrected comments:
	// constructs penalty matrix (del) multiplied by its transpose. Uses:
	// iruf = 0 for no roughening (i.e. min. norm models) (added in v1.5)
	// iruf = 1 for 1st derivative roughness
	// iruf = 2 for 2nd derivative roughness
	// On exit ptp() = del^T * del + W_B^T W_B

	// if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle) W_B are functions of (log10_mu)
	// & these are updated in this call		


	int npen,i;
	// float[][] del = new float[np][np];
	// DenseDoubleMatrix2D del = new DenseDoubleMatrix2D(np,np);

	// make partials matrix:
	npen = np; // as per objmat...
	objmat();

	// form ptp = del^T.del
	// A.zmult: C = alpha * A x B + beta*C.

	del.zMult(del,ptp,1.0,0.0,true,false);

	if (false)
	  {
		if (haveLargestEigenvaluePTP == false)
		  {
			// estimate largest eigenvalue of del^T.del
			largestEigenvaluePTP = -Double.MAX_VALUE;
			EigenvalueDecomposition eig = new EigenvalueDecomposition(ptp); 
			DoubleMatrix1D eR = eig.getRealEigenvalues();
			for (i=0;i<np;i++)
			  {
				largestEigenvaluePTP = max(largestEigenvaluePTP,abs(eR.getQuick(i)));
			  }
			// emergency care..
			if (largestEigenvaluePTP <= 0)
			  {
				// compare to prior...
				for (i=0;i<np;i++) 
				  {
					largestEigenvaluePTP = max(largestEigenvaluePTP,prewts[i]*prewts[i]);
				  }
				// if still in trouble...
				if (largestEigenvaluePTP <= 0) largestEigenvaluePTP = 1.0;
			  }
			haveLargestEigenvaluePTP = true;
			eig = null;
			eR = null;		
		  }
	  }

	//	EH.println("MAKPTP: PTP="+ptp.viewPart(0,0,8,8));

	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  {		
		updatePrewtsForExpStyleSmoothing(pow(10.0,log10_mu));
	  }
	// add the weights for the model prejudice
	for (i=0;i<np;i++) 
	  {
		ptp.setQuick(i,i,ptp.getQuick(i,i) + prewts[i]*prewts[i]);
	  }
	//au.csiro.JamesGunning.lib.util.Util.debugPrintMatrixToFile("PTP.txt",ptp);
	return;
  }

  void objmat()
  {
	// iruf = 0 for no roughnening (i.e. min. norm models) 
	// = 1 for 1st derivative roughness
	// 2 for 2nd derivative roughness
	// input del[np][np]
	// on output:
	// del() = penalty matrix
	// npen = number of rows in penalty matrix

	int i, j,npen;
	double x;
	float[] weights;

	//
	// for the 1d case the row count of penalty matrix is just np
	npen = np;
	del.assign(0.0);

	if (!zeroSmoothingWeights)
	  {
		if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
		  {
			// This mode corresponds to first derivative smoothing (iruf=1)
			// Weights are simple here, to achieve exp-covariance inverse of tridiagonal form
			// corrections are absorbed into prejudice on top and bottom corners of matrix
			for (i=0; i<np-1; i++)
			  {
				del.setQuick(i,i,-1.0);
				del.setQuick(i,i+1,1.0);
			  }
		  }
		else
		  {
			weights = new float[nlayer];
			// now turn this into a depth to thickness array
			// JG weights by depth from seafloor is equivalent to occam weighting
			// (SC remark: multiplying by 100 adds stability)
			weights[0] = thickf[1];
			for (i = 2; i<= (nlayer-1);i++) 
			  {
				weights[i-1] = weights[i-2] + thickf[i];
			  }
			weights[nlayer-1] = weights[nlayer-2] + thickf[nlayer-1];
			// now turn this into a depth to thickness array
			// (mulitplying by 100 adds stability)
			for (i = 1; i<=(nlayer-1); i++)
			  {
				weights[i-1] = 100.0f*thickf[i]/weights[i-1];
			  }
			weights[nlayer-1] = 100.0f*thickf[nlayer-1]/weights[nlayer-1];

			// min. norm models just pick off the size of the model parameters
			if (iruf == 0)
			  {
				for (i = 0; i<np; i++) del.setQuick(i,i,(double)(ipen1[iff[i]]));
			  }
			else
			  {
				// (iruf ==1 or 2)
				// min. first derivative take differences of adjacent layers
				for (i = 0;i<np-1;i++)
				  {
					// form a weight based on thickness to depth ratio
					// del(i,i) = -1.*float(ipen1(iff(i)))/depth(iff(i))
					// del(i,i+1) = 1.*float(ipen1(iff(i)))/depth(iff(i))
					x = (double)(ipen1[iff[i]]/weights[iff[i]-1]);
					del.setQuick(i,i,-x);
					del.setQuick(i,i+1,x);
				  }
			  }
			
			// square matrix and copy back into del if second deriv. smoothing required
			if (iruf == 2)
			  {
				// call mult(nnp,npen,nnp,npen,npen,del,del,del2);
				// DenseDoubleMatrix2D del2 = new DenseDoubleMatrix2D(np,np);
				del.zMult(del,del2,1.0,0.0,false,false); // A.zMult(..) : C = alpha * A x B + beta*C. 
				del.assign(del2);
			  }
		  }
		// EH.println("del="+del);
	  }
	
	// Set the weights for the prejudical model
	if (!(doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle))
	  {
		for (i=0;i<np;i++)
		  {
			prewts[i] = ipen2[iff[i]];
		  }
	  }
  }
      

  float minimiseMuOverDetCpAndOtherMuPieces(float pmu)
  {
	//-----------------------------------------------------------------------
	// JG.
	// Very crudely finds an approximate minimum of the sum of the objective terms 
	// mu ||(partial (pm-mbar))||^2 + ||Wp(mu)(pm-mbar)||^2 - log(|Cinv_p_eff(mu)|) 
	// + (in case of emprical Bayes) log|Ctilde(mu)|, where
	// Ctilde^{-1} = (WJ)^TWJ + mu \partial^T\partial+W_B^T(mu).W_B(mu)

	// with respect to the parameter pmu = log_10(mu), and sets
	// pmu to this minimum. Here Cinv_p_eff = mu partial^T partial + Wp^T Wp
	// i.e. the weighted sum of the smoothing and the prior prejudice.
	// In the exponential-covariance model, we have also Wp=Wp(mu)
	// Notice |partial^T partial| = 0, so we must always have |Wp^T Wp| != 0
	// i.e. *every* parameter must have a prior prejudice.
	// Also, approximate minimum is confined to local interval
	// max(pmumin,pmu-1/2) <= pmu <= min(pmumax,pmu+1/2)
	// so smoothing optimium is bounded in [-5,5]


	// local variables
	int i,imin; 
	float amu1,amu2,amu3,dmu;
	float fL,fR,f,mL,mR;
	float rough,muopt,muopt_min=0f;
	float[] min_mu={1.0f};
	int nsteps;
	float pmumax,pmumin;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D Ctilde=null;


	// Empirical Bayes bits
	if (doEmpBayes)
	  {
		// By here, 
		// wjtwj = (WJ)^T.(WJ)= WJ_trans. WJ_trans^T
		// Need additional term
		// log |(WJ)^T.(WJ)+ ptp(mu)|
		// Allocate space here, and use in called routines
		Ctilde = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D(np,np);
	  }

	if (zeroSmoothingWeights)
	  {
		// optimum is not well defined, and Lagrange mu has dropped out of problem
		// return simple solution
		min_mu[0] = 1.0f;
		return min_mu[0];
	  }


	// Algorithm parameters (make static?)....
	pmumax = 15.0f;
	pmumin = -1.0f;
	nsteps = 200;

	// muTrustSize = 0.5 (default)
	amu1 = (float)max(pmumin, pmu - muTrustSize);
	amu2 = (float)min(pmumax, pmu + muTrustSize);
	rough = chisqRoughness(pm);


	// initial scan for bisection...;
	dmu = (amu2-amu1)/(nsteps-1);
	for (i=1;i<=nsteps; i++)
	  {
		amu3 = amu1+(i-1)*dmu;
		muopt = mu_objective(amu3,rough,Ctilde);
		if ((i==1) || (muopt < muopt_min)) 
		  {
            muopt_min = muopt;
            min_mu[0] = amu3;
            imin = i;
		  }
		if (verbosity > 3) EH.printf("muopt: %15.8g  %15.8g%n",amu3,muopt);
      }
	// left ;
	mL = max(pmumin, min_mu[0] - dmu);
	fL = mu_objective(mL,rough,Ctilde);
	// right;
	mR = min(pmumax, min_mu[0] + dmu);
	fR = mu_objective(mR,rough,Ctilde);
	// middle;
	amu3 = 0.5f*(mL+mR);
	f = mu_objective(amu3,rough,Ctilde);
	if ((f<fL) && (f<fR))
	  {
		// 1D optimising search function ;
		muopt_min = findmin_det(mL,amu3,mR,rough,1.0e-5f,min_mu,Ctilde);
	  }
	if (verbosity > 2) EH.println("Lagrange multiplier scan: mu_min: "+min_mu[0]+" "+muopt_min);
	if (!true)
	  {
		PrintWriter out;	
		String filename = "mu_section."+nit+".txt";
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		for(amu1=0f;amu1<=pmumax;amu1+=0.1)
		  {
			f = mu_objective(amu1,rough,Ctilde);
			out.printf("%8.3g %15.8g%n",amu1,f);
		  }
		out.close();
		EH.println("Wrote: "+ filename);
	  }
	return min_mu[0];
  }



  double negLogDetCpEffInv(float alogmu,boolean alogmu_is_log10)
  {
	//-----------------------------------------------------------------------
	// 
	// Returnes log determinant of A = Cp_eff^{-1}, (combined smoothing 
	// and prejudice matrix), for current value of mu = 10^alogmu 
	// Code exploits (and assumes!) matrix A having 
	// centered bandwidth maximum of 5 columns (diag +/- 2)
	
	int i,j;
	double negLogDetCpEffInv;
	double amu,v,scale_guess;
	double c_im2,b_im2,beta_im2,beta_im3 ;
	double[] d;
	double[] a;          // a[i] = a(i,i), diagonal  
	double[] a_d1=null;  // a_d1[i] = a(i,i+1) next-upper diagonal
	double[] a_d2=null;  // a_d2[i] = a(i,i+2) next-next upper diagonal

	boolean debugDet=false;
	boolean debugMarg=false;
	
	if (alogmu == -99.0f)
	  {
		amu = 0.0;
	  }
	else
	  {
		if (alogmu_is_log10) amu = pow(10.0,(double)alogmu);
		else amu = (double)alogmu;
	  }
	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  {
		// use known analytical form
		double exp_m2alpha = getExpMinus2Alpha(amu);
		double sill = getExpPriorCovarianceSill();
		negLogDetCpEffInv = np*Math.log(sill) + (np-1)*Math.log(1 - exp_m2alpha);
		
		if (debugMarg || debugDet)
		  {
			DenseDoubleMatrix2D Q = new DenseDoubleMatrix2D(np,np);
			for (i=0;i<np;i++)
			  {
				double e = sqr(prewts[i]);
				e += ((i==0 || i==(np-1))? amu : 2*amu);  
				Q.setQuick(i,i,e);
				if (i<(np-1)) Q.setQuick(i,i+1,-amu);
				if (i>0) Q.setQuick(i,i-1,-amu);
			  }
			Algebra alg = new Algebra();
			// EH.println("Q:="+Q);
			if (debugDet) EH.println("CERN: -Math.log( |Cinvp_eff| ) = "+(-Math.log(alg.det(Q))));		
			DenseDoubleMatrix2D Cp_eff;
			Cp_eff = (DenseDoubleMatrix2D)alg.inverse(Q);
			
			EH.println("Marginal prior std deviations: ");
			for (i=0;i<np;i++) EH.printf("%12.6g ",Cp_eff.getQuick(i,i));
			EH.printf("%n");
		  }
	  }
	else
	  {
		d = new double[np+2];
		a = new double[np]; // diag of DenseDoubleMatrix2D(np,np);

		// ptp   = partial^T.partial + W_B^T.W_B

		// reconstruct upper triangle matrix representation of
		// a = amu*partial^T.partial + W_B^T.W_B
		// (1) Diagonal.
		for (i=0;i<np;i++)
		  {
			a[i]= amu*ptp.getQuick(i,i) + (1.0 - amu)*prewts[i]*prewts[i];
		  }

		// (2) upper triangle off-diagonals
		if (iruf>0)
		  {
			a_d1 = new double[np-1];
			for (i=0;i<np-1;i++)
			  {
				// a.setQuick(i,i+1, amu*ptp.getQuick(i,i+1));
				a_d1[i] = amu*ptp.getQuick(i,i+1);
			  }
		  }
		if (iruf==2)
		  {
			a_d2 = new double[np-2];
			for (i=0;i<np-2;i++)
			  {
				// a.setQuick(i,i+2, amu*ptp.getQuick(i,i+2));
				a_d2[i] = amu*ptp.getQuick(i,i+2);
			  }
		  }

		// compute average multiplier scale: scale_guess = 1/<diag{a}>
		scale_guess=0.0;
		for (i=0;i<np;i++) scale_guess += a[i];
		scale_guess = np/scale_guess;

		if (debugDet || debugMarg) 
		  {
			debugPrintCpInv(a,a_d1,a_d2);		

			DenseDoubleMatrix2D Q = new DenseDoubleMatrix2D(np,np);
			for (i=0;i<np;i++)
			  {
				for (j=0;j<np;j++)
				  {
					if (j==i) {Q.setQuick(i,i,a[i]);}
					else if ((j==i+1) && (a_d1!=null)) Q.setQuick(i,j,a_d1[i]);
					else if ((j==i-1) && (a_d1!=null)) Q.setQuick(i,j,a_d1[i-1]);
					else if ((j==i+2) && (a_d2!=null)) Q.setQuick(i,j,a_d2[i]);
					else if ((j==i-2) && (a_d2!=null)) Q.setQuick(i,j,a_d2[i-2]);
				  }
			  }
			Algebra alg = new Algebra();
			//		EH.println("CERN: det(Cinvp_eff) = "+alg.det(Q));		
			if (debugDet)	
			  {
				EH.println("Q:="+Q);
				EH.println("CERN: -log( |Cinvp_eff| ) = "+(-Math.log(alg.det(Q))));		
			  }
			if (debugMarg)
			  {
				DenseDoubleMatrix2D Cp_eff;
				Cp_eff = (DenseDoubleMatrix2D)alg.inverse(Q);
			
				EH.println("Marginal prior std deviations: ");
				for (i=0;i<np;i++)
				  {
					EH.printf("%12.6g ",Cp_eff.getQuick(i,i));
				  }
				EH.printf("%n");
			  }
		  }



		// Form determinant of Cp^{-1}: pull out constant factor of scale_guess 
		// to try to keep scales sane and avoid overflows etc
		if (iruf==0) 
		  {			
			// diagonal form
			d[1] = 1.0;
			for(i=0;i<np;i++)
			  {
				d[1] *= a[i] * scale_guess;
			  }
		  }
		else if (iruf==1)
		  {
			// tridiagonal form
			d[np+1] = 1.0;
			d[np] = a[np-1]*scale_guess;
			for (i=np-1; i>=1; i--)
			  {
				v = a_d1[i-1]*scale_guess;
				d[i] = (a[i-1]*scale_guess)*d[i+1] - (v*v)*d[i+2];
			  }
		  }
		else if (iruf==2)
		  {
			// Symmetric pentadiagonal form
			// From "A recursive relation for the determinant of a pentadiagonal matrix"
			// Roland Sweet, Comm. ACM, 12,6,1969
			// Definitions...
			// a_i = a(i,i)
			// b_i = a(i,i+1)**2                      // i < np
			// beta_i = a(i,i+2)**2                   // i < np-1
			// c_i = a(i,i+1) * a(i+1,i+2) * a(i,i+2) // i < np-1
			// Initialisations
			// Perform recursions with all elements scaled by scale_guess
			// scale_guess=1.0;
			for (i=0;i<np;i++) a[i] *= scale_guess;
			for (i=0;i<np-1;i++) a_d1[i] *= scale_guess;
			for (i=0;i<np-2;i++) a_d2[i] *= scale_guess;

			double tmp;
			d[0] = 1.0;
			d[1] = a[0];
			if (np > 1) d[2] = a[0]*a[1] - a_d1[0]*a_d1[0];
			// 6 term recursion, Eqn 11 of Sweet
			for (i = 3; i<=np; i++)
			  {
				// Legible version of equations....
				// d(i)=(a(i)-c(i-2)/b(i-2))*d(i-1)
				// -(b(i-1)-a(i-1)*c(i-2)/b(i-2))*d(i-2)
				// -(beta(i-2)*a(i-1)-c(i-2))*d(i-3)
				// if (i > 4) then
				// d(i)=d(i)+beta(i-3)*(beta(i-2)-a(i-2)*c(i-2)/b(i-2))*d(i-4)
				// if (i > 5) then
				// d(i) = d(i)+ (beta(i-3)*beta(i-4)*c(i-2)/b(i-2))*d(i-5)
				// end if
				// end if

				c_im2 = a_d1[i-3] * a_d1[i-2]* a_d2[i-3];
				b_im2 = a_d1[i-3] * a_d1[i-3];
				beta_im2 = a_d2[i-3] * a_d2[i-3];
				tmp = c_im2/((b_im2!=0.0) ? b_im2 : 1.0);
				d[i] = (a[i-1] - tmp)*d[i-1]
				  -(a_d1[i-2]*a_d1[i-2] - a[i-2]*tmp) * d[i-2]
				  -(beta_im2*a[i-2] - c_im2)*d[i-3];
				if (i>3) 
				  {
					beta_im3 = a_d2[i-4] * a_d2[i-4];
					d[i] += (beta_im3*(beta_im2 - a[i-3]*tmp)*d[i-4]);
					if (i>4)
					  {
						d[i] += (beta_im3 * a_d2[i-5] * a_d2[i-5] * tmp)*d[i-5];
					  }
				  }
			  }

			d[1]=d[np];
		  }
		// EH.printf("Function returns: %g%n",d[1]);

		negLogDetCpEffInv = -(Math.log(d[1]) - np*Math.log(scale_guess));

	  }

   	if (debugDet)  EH.println("Function returns: "+negLogDetCpEffInv);

	// negLogDetCpEffInv -=log(amu); // from Jeffreys prior?

	return negLogDetCpEffInv;
  }


  double negLogDetCtildeEffInv(float alogmu, DoubleMatrix2D Ctilde)
  {
	// Returns log determinant of A = Ctilde_eff^{-1} = J^W^TWJ + mu \partial^T\partial+W_B^T(mu).W_B(mu)
	// and prejudice matrix), for current value of mu = 10^alogmu 

	int i,Np;
	double amu,logDet;

	Np = Ctilde.rows();
	logDet = 1.0;

	amu = pow(10.0,(double)alogmu);

	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  { 
		// need to recompute W_b since the "prejudice weights" Wb change with mu. 
		// Ensure we recompute prewts
		// The matrix del is already in-place from objmat():
		// Util.debugPrintMatrixToFile("del.txt",del);
		// form Ctilde = amu * del^T.del
		if (false)
		  {
			// works ok
			del.zMult(del,Ctilde,amu ,0.0,true,false);
		  }
		else
		  {
			// faster, band multiply
			double s;
			int j,k,ks,ke,b=1; // b = bandwidth
			for(i=0;i<Np;i++)
			  {
				ks = max(0,i-b);
				ke = min(Np-1,i+b);
				for(j=0;j<=i;j++)
				  {
					for(k = ks,s=0.0; k<= ke; k++) s += del.getQuick(k,i) * del.getQuick(k,j);
					Ctilde.setQuick(i,j,amu * s);
					if (j!=i) Ctilde.setQuick(j,i,amu * s);
				  }
			  }
		  }
		//Util.debugPrintMatrixToFile("Ct1.txt",Ctilde);
		updatePrewtsForExpStyleSmoothing(amu);
		// add the weights for the model prejudice (\cal{W}_p(\mu))
		for (i=0;i<Np;i++) 
		  {
			Ctilde.setQuick(i,i,Ctilde.getQuick(i,i) + prewts[i]*prewts[i]);
		  }
		
		// and add Jacobean piece: Ctilde += wjtwj
		Ctilde.assign(wjtwj, Functions.plus);

		// now get log|Ctilde|. I think Cholesky is cheapest way
		CholeskyDecomposition chol = new CholeskyDecomposition(Ctilde);
		DoubleMatrix2D L = chol.getL();
		// determinant of inverse posterior covariance
		for (i=0,logDet=0.0; i<Np; i++) logDet += Math.log(L.getQuick(i,i)); 
		if (modelSystematics)
		  {
			EH.println("FIXME Update 10 lines above or systematics!!!!!!!!!");
		  }
	  }
	return logDet;
  }

  // Objective for use in univariate minimization on mu in findmin_det
  // requires current model to be in pm[]

  float mu_objective(float log10_mu, float rough, DoubleMatrix2D Ctilde)
  {
	float f;


	f = chisqPrej(log10_mu,pm); // evaluate first to set prewts...
	f += fpow10(log10_mu)*rough;
	f += (float)negLogDetCpEffInv(log10_mu,true);
	if (Ctilde != null)
	  {
		// extra piece for empirical Bayes
		f += (float)negLogDetCtildeEffInv(log10_mu,Ctilde);
	  }

	return f;
  }

  float findmin_det(float ax,float bx,float cx,float rough,float tol,float[] xmin,DoubleMatrix2D Ctilde)
  {
	//-----------------------------------------------------------------------
	// JG: Optimising search method for negLogDetCpEffInv() function
	// on input:
	// ax,bx,cx = independent variable which bracket the minimum
	// on output:
	// xmin[0] = abscissa of minimum

    // Algorithm control parameters...
	float r=0.61803399f;
	float c=1.0f-r;
	int MAX_ITS = 100;

	int i;
	float f_min, f;
	float f1,f2,x0,x1,x2,x3;

	if (zeroSmoothingWeights)
	  {
		// optimum is not well defined, and Lagrange mu has dropped out of problem
		// return simple solution
		xmin[0]=ax;
		return 0.0f;
	  }

	x0=ax;
	x3=cx;
	if (abs(cx-bx) > abs(bx-ax))
	  {
		x1=bx;
		x2=bx+c*(cx-bx);
	  }
	else
	  {
		x2=bx;
		x1=bx-c*(bx-ax);
	  }

	// f1=f(x1)
	f1 = mu_objective(x1,rough,Ctilde); // (float)negLogDetCpEffInv(x1,true) + fpow10(x1)*rough;
	// f2=f(x2)
	f2 = mu_objective(x2,rough,Ctilde); // (float)negLogDetCpEffInv(x2,true) + fpow10(x2)*rough;
	for(i=1;i<=MAX_ITS;i++)
	  {
		if (abs(x3-x0) > tol*(abs(x1)+abs(x2)))
		  {
            if (f2 < f1)
			  {
				x0=x1;
				x1=x2;
				x2=r*x1+c*x3;
				f1=f2;
				// f2=f(x2)
				f2 = mu_objective(x2,rough,Ctilde);// (float)negLogDetCpEffInv(x2,true) + fpow10(x2)*rough;
			  }
            else
			  {
				x3=x2;
				x2=x1;
				x1=r*x2+c*x0;
				f2=f1;
				// f1=f(x1)
				f1 = mu_objective(x1,rough,Ctilde); // (float)negLogDetCpEffInv(x1,true) + fpow10(x1)*rough;
			  }
		  }
		else
          {
			break;
		  }
	  }
	if (i > MAX_ITS)
	  { 
		if (verbosity > 2) EH.println("Warning: golden search failed (MAX_ITS) in findmin_det()");
	  }

	if (f1 < f2) 
	  {
		f_min=f1;
		xmin[0]=x1;
	  }
	else
	  {
		f_min=f2;
		xmin[0]=x2;
	  }

	return f_min;
  }

  /**
	 Computation of inverse covariance matrix full-Hessian 
	 at optimisation termination.

	 On calling, expects
	 wjtwj = (WJ)^T WJ 

	 N.B.: OVERWRITES wjtwj with
	 wjtwj = (WJ)^T WJ + mu partial^T.partial + W_B^T.W_B
	 for leading block.

	 Other blocks of Hessian filled in as per current parameter values also.

	 returns current full Hessian.
  */
  DenseDoubleMatrix2D computeFullHessian(float log10mu)
  {
	int i,j;
	cern.colt.matrix.impl.DenseDoubleMatrix2D H;
	double amu;

	int Nptot = np + SEM.nPsystematics;

	H = new DenseDoubleMatrix2D(Nptot,Nptot);
	assembleHessian(fpow10(log10mu),H);

	// Overwrite inverse-covariance back into wjtwj for subsequent usage...
	wjtwj.assign(H.viewPart(0,0,np,np));
	return H;
  }


  /**
	 Computation of of covariance matrix & stddev
	 Expects, on calling:

	 (1) wjtwj = (WJ)^T WJ + mu partial^T.partial + W_B^T.W_B
	 (2) H = correct full Hessian, with H11 = wjtwj

	 returns ModeCovarianceAndStdDev object
	 containing modeCov.C and sqrt of its diagonal
	 in modeCov.stdv
  */
  ModeCovarianceAndStdDev computeUncertainty()
  {
	int i,n,Nptot;
	float[] stdv;
	cern.colt.matrix.impl.DenseDoubleMatrix2D Cpost;
	ModeCovarianceAndStdDev modeCov;
	DoubleMatrix2D Ainv;

	// By here, we expect wjtwj = (WJ)^T WJ + mu partial^T.partial + W_B^T.W_B, and
	// H11 = wjtwj, also other elements of H in place
	
	modeCov = new ModeCovarianceAndStdDev();

	Nptot = H.rows();
	
	// Bulletproof inverse
	// modeCov.C  = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(wjtwj); 
	modeCov.C  = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(H); 

	modeCov.stdv = new float[Nptot];
	for (i=0; i<Nptot; i++)
	  {
		modeCov.stdv[i] = (float)sqrt(modeCov.C.getQuick(i,i));
	  }

	return modeCov;
  }


  
  float getMaxThickness()
  {
	int i,N;
	float max_thickness=0f;

	N = newScripps ? (nlayer-1) : nlayer;
	for (i=1; i <= N; i++) max_thickness += thickf[i]; 
	return max_thickness;
  }

  void overwriteLogGridModel()
  {
	int i,j,k;
	float max_thickness,dz,ztot,zSeafloor,rhoSea,thickSea,z,pm_startval;

	// initial total subsea model thickness
	max_thickness = getMaxThickness();

	if (!newScripps)
	  {
		rhoSea = rho[0]; // save these
		thickSea = thickf[0];
		pm_startval = pm[0];
		dz = logModelFirstThickness;
		ztot=dz;
		for (i=1; ztot<max_thickness; i++) {dz *= logModelRatio; ztot+=dz;}		
		np = nlayer = i; // OLD i-1;
		if (verbosity > 1) 
		  {
			//EH.printf("New log-grid (total-original-thickness =%3g, total-new-thickness= %3g): N logarithmic layers = %d%n",max_thickness,ztot,np);
			EH.printf("New log-grid: N logarithmic layers = %d%n",np);
		  }
		allocateModel();
		allocateGeoModel();
		// repopulate
		rho[0] = rhoSea;
		thickf[0] = thickSea;
		dz = logModelFirstThickness;
		for (i=1,z=0; i<=nlayer; i++)
		  {
			rho[i] = 0.0f;
			pm[i-1] = pm_startval;
			thickf[i] = dz;
			ipen1[i] = 1.0f;   // roughness penalty
			premod[i-1] = PTI.getDefaultPriorMean(z+dz/2f);   // prior mean
			ipen2[i] = 1.0f/PTI.getDefaultPriorSigma(z+dz/2f);      // prior weight
			z += thickf[i];
			dz *= logModelRatio;
		  }
		for (i = 0; i<nlayer;i++) iff[i] = i+1;
	  }
	else
	  {
		// New Scripps

		// ALL THE CHUNK BELOW NEEDS CAREFUL SCRUTINY. Probably way too long & messy. Fix it!!.

		if (nFixedSurfaceLayers!=2)
		  {
			EH.die_err("New scripps kernel and --log-grid option works only for fixed upper 2 layers at present. Fix your grid manually!");
		  }
		int i_start=1;
		SMI.copyFixedUpperLayers();   // Copy upper fixed layers
		zSeafloor  =(float)SMI.CM.layer_params_trans.getQuick(0,2); // depth of seafloor
		float[] rho_copy = (float[])rho.clone();
		float[] thickf_copy = (float[])thickf.clone();
		pm_startval = pm[0];
		ztot = dz = logModelFirstThickness;
		for (i=1; ztot<max_thickness; i++) {dz *= logModelRatio; ztot+=dz;}		
		np = i; 
		nlayer = np; // number subsea layers 
		if (verbosity > 1) EH.printf("New log-grid: N (free) logarithmic layers = %d%n",np);

		allocateModel();
		allocateGeoModel();
		// expect, by here, DM.nlay1D is set consistently (as per allocateModel())
		SMI.CM.allocateL(SMI.DM.nlay1D);
		SMI.DM.allocateL(SMI.DM.nlay1D);
		SMI.restoreFixedUpperLayers();

		// repopulate
		for (i=0; i < i_start; i++) 
		  {
			rho[i] =  rho_copy[i];
			thickf[i] = thickf_copy[i];
		  }
		rho_copy=null;
		thickf_copy=null;

		dz = logModelFirstThickness;
		k = 2;
		ztot = zSeafloor;
		for (i=i_start,z=0,j=0; i<=nlayer; i++,j++,k++)
		  {
			// here, in loop ztot= absolute depth
			//               z = depth below seafloor
			rho[i] = 0.0f;
			pm[j] = pm_startval;
			thickf[i] = dz;
			ipen1[i] = 1.0f;   // roughness penalty
			premod[j] = PTI.getDefaultPriorMean(z+dz/2f);        // prior mean
			ipen2[i] = 1.0f/PTI.getDefaultPriorSigma(z+dz/2f);   // prior weight

			// SMI interface bits: [top_depth, rho, penalty, Preference, prej_penalty] 
			SMI.CM.layer_params_trans.setQuick(0,k,ztot);
			SMI.CM.layer_params_trans.setQuick(1,k,-1.0); // free
			// no need to set the rest
			z += thickf[i];
			ztot += z;
			dz *= logModelRatio;
		  }
		// update free layers
		for (i = 1,j=0; i<=nlayer;i++)
		  {
			if (rho[i] <= 0.0) iff[j++] = i;
		  }
		SMI.freeJacobianAllocs();
		SMI.finishModelSetup();
	  }
  }


  /* 
	 A method to read in a "profile_final.txt" file, & optionally 
	 map the model therein to a logarithmically gridded model (using conservation of
     rho*thickness based upscaling), resets the model variables, etc, so as
     to solve a regularly-gridded model "inspired" by the supplied
     profile. Prior variances are set to unity, and the prejudice to
     the supplied re-gridded model also
  */
  void overwriteGridModelFromProfile(String filename)
  {
	int NresMax = LOG2_MAX_NUM_SPLIT_GRID_ELEMENTS;    // 2^NresMax possible subdivisions: probably 16 or 32 max
	int i,j,js,je,nz,NlayerMax;
	float z=0f,r=0f,max_thickness,rhoSea,thickSea,dz,zs=0f,ze=0f,rhoH;
	float[] zL,rL;
	boolean fS,fE;
	String s;
	String[] st;
	BufferedReader in;
	FloatArrayList zList = null,rList = null;

	if (verbosity > 1) EH.println("Reading profile file for regridding: "+filename); 
	try 
	  {
		in = new BufferedReader(new FileReader(filename));
		zList = new FloatArrayList();
		rList = new FloatArrayList();

		j=0;
		while ((s = in.readLine()) != null)
		  {
			s = s.trim();
			if (s.startsWith("#")) continue; // skip comments
			st = s.split("\\s+");
			if (st.length==2) 
			  {
				try
				  {
					z = Float.parseFloat(st[0]);
					r = Float.parseFloat(st[1]);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Line containing "+s+" not two numbers: depth log10(rho)");
				  }
				// save every second line
				if (j%2 == 0)
				  {
					zList.add(z);
					rList.add(r);
				  }
				j++;
			  }
		  }
		in.close();
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open or close " + filename);
	  }
	zList.add(z);   // add also thickness of last layer
	zList.trimToSize();
	nz = zList.size();

	rList.trimToSize();

	thickSea = zList.get(1);
	// total subsea model thickness. z has last read-in depth
	max_thickness = z - thickSea;
	rhoSea = rList.get(0); // save this

	if (overwriteLogModel)
	  {
		dz = logModelFirstThickness;
		float ztot = dz;
		for (i=1; ztot<max_thickness; i++) {dz *= logModelRatio; ztot+=dz;}		
		NlayerMax = i-1;
		if (verbosity > 1) EH.println("New log-grid: N logarithmic layers = "+NlayerMax);
	  }
	else
	  {
		// copy profile to target model exactly as-is
		NlayerMax = rList.size()-1; // omit sea
	  }
	np = nlayer = NlayerMax; // # subsea layers
	if (verbosity > 1) EH.println("New log-grid mapped from profile file: N layers = "+NlayerMax);
	allocateModel();
	allocateGeoModel();

	// repopulate
	rho[0] = (float)pow(10,rhoSea);
	thickf[0] = thickSea;

	zL = zList.elements();
	rL = rList.elements();
	for (j=0; j<rL.length; j++) rL[j] = (float)pow(10,rL[j]);
	js=0; 
	je=0;
	ze = thickSea;
	dz = logModelFirstThickness;
	z=0f;
	for (i=0 ; i< nlayer; i++)
	  {
		if (overwriteLogModel)
		  {
			zs = ze;
			ze = zs + dz;
		  }
		else
		  {
			zs = zL[i+1];
			ze = zL[i+2];
			dz = ze - zs;
		  }
		// search for start index js of profile segment that straddles start of new-grid segment i
		fS=fE=false;
		for (j=max(js-1,0); j<(nz-1); j++)
		  {
			if ((zL[j]<zs) && (zL[j+1]>=zs))
			  {
				js = j; fS=true; break;
			  }
		  }
		// search for start index je of profile segment that straddles end of new-grid segment i
		for (j=max(js-1,0),je=js;j<(nz-1);j++)
		  {
			if ((zL[j]<ze) && (zL[j+1]>=ze))
			  {
				je = j; fE=true; break;
			  }
		  }
		if (!fS || !fE) EH.warn("Regridding suspect!");

		// accumulate resistivity-thickness integral over start & end fragments, plus internal pieces
		rhoH = (min(ze,zL[js+1]) - zs) * rL[js];
		if (je > js) 
		  {
			rhoH += (ze - max(zs,zL[je])) * rL[je];
		  }
		for (j=js+1; j<je; j++)
		  {
			rhoH += (zL[j+1] - zL[j]) * rL[j];
		  }
		// scale back to equivalent resistance required for segment of depth range zs<z<ze
		rhoH /= dz;
		// and use, as starting point:		
		rhoH = flog10(rhoH);
		pm[i] = rhoH;
		thickf[i+1] = dz;
		z += dz; // TVDBML
		// But prior distribution informed by...
		ipen1[i+1] = 1.0f;   // roughness penalty
		if (readProfileSetMean)
		  {
			premod[i] = pm[i];   // prior mean
		  }
		else
		  {
			premod[i] = 0;       // prior mean
		  }
		ipen2[i+1] = 1f/PTI.getDefaultPriorSigma(z+dz/2f);       // prior weight

		if (overwriteLogModel)
		  {
			dz *= logModelRatio;
		  }
	  }
	// but ensure resistivities still inverted for
	for (i = 0; i<nlayer;i++) 
	  {
		rho[i+1] = -1.0f;
		iff[i] = i+1;
	  }

	if (!false)
	  {
		// syserr printout of original model & upscaled/regridded one		
		filename = "regridded_model.txt";
		try
		  {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			for(j=0; j<nz-1; j++)
			  {
				out.println(zL[j]+" "+log10(rL[j]));
				out.println(zL[j+1]+" "+log10(rL[j]));
			  }
			out.println("  ");
			zs = thickSea;
			for (i=0 ; i< nlayer; i++)
			  {
				ze = zs + thickf[i+1];
				out.println(zs+" "+pm[i]);
				out.println(ze+" "+pm[i]);
				zs = ze;
			  }
			out.close();
			EH.println("Wrote graphs of regridded model to "+filename);
		  }
		catch (IOException e)
		  {
			throw new RuntimeException("couldn't open or close " + filename);	  
		  }
	  }
  }

  void diagnoseMatrix(DoubleMatrix2D a)
  {
	System.out.printf("Eigenvalues: ");
	EigenvalueDecomposition eig = new EigenvalueDecomposition(a); 
	DoubleMatrix1D eR = eig.getRealEigenvalues();
	for (int i=0;i<eR.size();i++) System.out.printf("%g ",eR.getQuick(i));
	System.out.println("");
	double eMax=-Double.MAX_VALUE,eMin=Double.MAX_VALUE;
	for (int i=0;i<eR.size();i++)
	  {
		double x = abs(eR.getQuick(i));
		eMax = max(eMax,x);
		eMin = min(eMin,x);
	  }	
	System.out.printf("Condition number = %g%n",eMax/eMin);
  }

  String getFileText(boolean postSorting,int indx,String extraText)
  {
	String text;
	if (postSorting)
	  {
		text= String.format(".mode_%04d",indx);
	  }
	else if (doMultistart || doNoiseCooling || doRelaxedSmoothing)
	  {
		text= String.format(".start_%04d",indx);
	  }
	else
	  {
		text = "";
	  }
	if (extraText != null) text = text.concat(extraText);
	return text;
  }

  /**
	 Dumps of covariance matrix & other statistical entities at optimisation 
	 termination. Ascii dumps if dumpAscii=true.
	 Expects, on calling:

	 (1) modeCov to contain correct posterior matrix and stddv at local optimum
	 (2) pm to contain correct model at local optimum

	 when called, i.e. from computeFullHessian();

  */
  void writeUncertaintyFiles(String extraText, boolean dumpAscii, boolean postSorting, 
							 int indx, ModeCovarianceAndStdDev modeCov)
  {
	int i,j,Nptot,i_free;
	float mv,mean,err,y,z;
	float[] stdv;
	String filename;
	PrintWriter out,out_simple;

	Nptot = modeCov.C.rows();

	if (dumpAscii)
	  {
		String text;

		text = getFileText(postSorting,indx,extraText);

		// EH.println("Cinv = "+a);
		if (verbosity > 1) EH.printf("Param std. devs: ");
		for (i=0;i<Nptot;i++)
		  {
			EH.printf("%-15.6g ",modeCov.stdv[i]);
		  }
		EH.printf("%n");
		
		if (verbosity > 1) EH.println("Writing Covariance matrix...");
		filename = String.format("Covariance"+text+".txt");
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		for (i=0;i<Nptot;i++)
		  {
			for (j=0;j<Nptot;j++) out.printf("%15.8g ",modeCov.C.getQuick(i,j));
			out.printf("%n");
		  }
		out.close();
		
		if (verbosity > 1) EH.println("Writing Correlation. coeff. matrix...");
		filename = String.format("CorrCoeff"+text+".txt");
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		for (i=0;i<Nptot;i++)
		  {
			for (j=0;j<Nptot;j++) 
			  {
				out.printf("%15.8g ",modeCov.C.getQuick(i,j)/(modeCov.stdv[i]*modeCov.stdv[j]));
			  }
			out.printf("%n");
		  }
		out.close();


		filename = String.format("SummaryStats"+text+".txt");
		if (verbosity > 1) EH.println("Writing SummaryStats file: "+filename);
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		// write params & std deviations: 3 rows: 
		for (i = 0; i< Nptot; i++) 
		  {
			y = pm[i]-modeCov.stdv[i];
			if (useBoundConstraints && i<np) y = min(max(y, m_lower[i]),m_upper[i]); // truncate
			out.printf("%15.7g ",y);   // mean -1 std dev
		  }
		out.println();
		for (i = 0; i< Nptot; i++) out.printf("%15.7g ",pm[i]);           // mean 
		out.println();
		for (i = 0; i< Nptot; i++) 
		  {
			y = pm[i]+modeCov.stdv[i];
			if (useBoundConstraints && i<np) y = min(max(y, m_lower[i]),m_upper[i]); // truncate
			out.printf("%15.7g ",y);   // mean +1 std dev
		  }
		out.println();
		out.close();



		if (verbosity > 1) EH.println("Writing Profile files");
		filename = String.format("profile_final_err"+text+".txt");
		current_profile_err_filename = filename;
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		filename = String.format("profile_final"+text+".txt");
		current_profile_filename = filename;
		try
		  {
			out_simple = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		// write depths versus resistivities
		for (j=-1; j<=1; j++)
		  {
			i_free = 0;
			for (i = 0; i<= nlayer; i++)
			  {
				if (rho[i] <= 0.0)
				  {
					err = modeCov.stdv[i_free];
					mean = pm[i_free];
					mv = mean + ((j!=0)? (j*err) : 0);
					if (useBoundConstraints) 
					  {
						mv = min(max(mv, m_lower[i_free]),m_upper[i_free]); // truncate
					  }
					i_free++;
				  }
				else 
				  {
					err = 0.0f;
					mean = (i==0) ? flog10(rho[0]) : rho[i]; // CSEM case has rho_sea in real units for sea-layer
					mv = mean;
				  }
				// error file
				out.printf("%15.7g %15.7g%n",depths[i],mv);
				out.printf("%15.7g %15.7g%n",depths[i+1],mv);

				// simple file
				if (j==0)
				  {
					out_simple.printf("%15.7g %15.7g%n",depths[i],mv);
					out_simple.printf("%15.7g %15.7g%n",depths[i+1],mv);
				  }
			  }
			if (j!=1) out.printf("%n");
		  }
		out_simple.close();
		out.close();
	  }
  }

  void Xsections(float pmu,boolean postSorting, int indx)
  {
	int i,j,k,npts,Nptot;
	float tmp,mu,tobt_full,v,scale;
	String filename,text;
	PrintWriter out=null;

	// number of points in Xsection 
	npts = 6;

	mu = fpow10(pmu);
	if (verbosity > 1) EH.println("Writing Xsections:");
	if (postSorting)
	  {
		text= String.format("mode_%04d",indx);
	  }
	else if (doMultistart) 
	  {
		text= String.format("start_%04d",indx);
	  }
	else 
	  {
		text = "";
	  }
	Nptot = np + (modelSystematics ? SEM.nPsystematics : 0);
	for (i=0;i<Nptot;i++)
	  {
		filename = String.format("Xsection.%s.%d.txt",text,i);
		try
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		scale = (float)sqrt(H.getQuick(i,i));		
		scale = 1.0f/scale;
		tmp = pm[i];
		// scale=1.0f;
		for (j=-npts;j<=npts;j++) 
		  {
			pm[i] = tmp + (scale*j)/(float)npts;
			doForwardModel(pm);
			// ScrippsLegacyModelInterface.formod(np,nd,pm,dp.elements,dm,d);
			tobt_full = chisqMisFit();
			if (doNoiseEstimation)
			  {
				//tobt_full = (float)(nd * (1.0+Math.log(tobt_full/nd))); // incorporates optimization on noise-scalar
				tobt_full = DOF(nd,np) + nd * (float)Math.log(tobt_full/DOF(nd,np)); // EB version 
			  }
			tobt_full += mu * chisqRoughness(pm);
			tobt_full += chisqPrej(pmu,pm);
			out.printf("%15.8g %15.8g%n",pm[i],tobt_full);
		  }
		out.close();
		if (verbosity > 1) EH.println("Wrote "+filename);
		pm[i] = tmp;
	  }
	if (verbosity > 1) EH.println("Finished Xsections");
  }


  void findOptimalLinearBackground()
  {
	int i,j,k,jmin=0,i_free,NTRY=20;
	float rhoM,v,tobt_full,tobt_full_min,jmin_opt;

	float RHO_MAX=3.0f;
	float RHO_MIN=0.0f;

	if (verbosity>1) EH.println("Estimating optimal linear prior profile");

	// save initial guess read from files
	float[] pm_cpy = (float[])pm.clone();
	float[] f = new float[NTRY];

	tobt_full_min = Float.MAX_VALUE;
	for (j=0; j<NTRY; j++)
	  {
		rhoM = j*(float)RHO_MAX/(float)NTRY + RHO_MIN;

		// triangular profile from RHO_MIN to rhoM
		for (i = 0,i_free=0; i<= nlayer; i++)
		  {
			if (rho[i] <= 0.0)
			  {
				pm[i_free] = rhoM * (depths[i+1]- depths[1])/(depths[np]- depths[1]) + RHO_MIN;
				i_free++;
			  }
		  }
		doForwardModel(pm);
		// ScrippsLegacyModelInterface.formod(np,nd,pm,dp.elements,dm,d);
		tobt_full = chisqMisFit();
		f[j] = (float)sqrt(tobt_full/nd);
		if (verbosity>2) EH.printf("linear profile: base log10(rho)=%-5.2g,  RMS_misfit=%-5.2g%n",rhoM,f[j]);
		if (tobt_full < tobt_full_min)
		  {
			jmin=j;
			tobt_full_min = tobt_full;
		  }
	  }
	if (jmin>0 && jmin<(NTRY-1))
	  {
		// maxima is internal. Use effective jmin from minima of parabolic fit
		jmin_opt = 2 * (f[jmin-1] - 2 * f[jmin] + f[jmin+1]);
		jmin_opt = (jmin_opt!=0.0) ? (jmin + (f[jmin-1] - f[jmin+1])/jmin_opt) : (float) jmin;
	  }
	else jmin_opt = (float) jmin;

	// restore profile
	rhoM = jmin_opt*(float)RHO_MAX/(float)NTRY + RHO_MIN;	
	// triangular profile from RHO_MIN to rhoM
	// And set prejudice trend
	for (i = 0,i_free=0; i<= nlayer; i++)
	  {
		if (rho[i] <= 0.0)
		  {
			pm[i_free] = rhoM * (depths[i+1]- depths[1])/(depths[np]- depths[1]) + RHO_MIN;
			premod[i_free] = pm[i_free];
			i_free++;
		  }
	  }
	if (verbosity>1) 
	  {
		EH.printf("Estimated profile: base log10(rho)=%-5.2g,  with RMS_misfit=%-7.3g%n",rhoM,f[jmin]);
	  }

	PrintWriter out;
	String filename = "profile.estimated-trend.txt";
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }
	// write depths verses log10 resistivities as a plottable line-graph
	for (i = 0, i_free=0; i<(nlayer+1); i++)
	  {
		if (rho[i] <= 0.0f)
		  {
			v = pm[i_free];
			i_free++;
		  }
		else 
		  {
			v = (i==0) ? flog10(rho[i]) : rho[i]; 
		  }
		out.printf("%15.7g %15.7g%n",depths[i],v);
		out.printf("%15.7g %15.7g%n",depths[i+1],v);
	  }
	out.close();
	EH.println("Wrote "+filename);

	// Restore initial guess from files ?
	// System.arraycopy(pm_cpy,0,pm,0,np);
  }



  private void showParams(float[] x,boolean showPrej)
  {
	int i,j;
	
	if (doModelSelectAutoSplitting)
	  {
		EH.printf("Thicknesses: ");
		for(i=1;i<=np;i++) EH.printf("%-15.8g ",thickf[i]);
		EH.println("");
	  }
	EH.printf("Parameters: ");
	for(i=0;i<np;i++) EH.printf("%-15.8g ",x[i]);
	EH.println("");
	if (SEM.nPsystematics > 0)
	  {
		EH.printf("Systematic Parameters: ");
		for (i=np,j=0; j<SEM.nPsystematics; i++,j++) EH.printf("%-15.8g ",x[i]);
		EH.println("");
	  }
	if (showPrej)
	  {
		EH.printf("%nParam prejudices: ");
		for(i=0;i<np;i++) EH.printf("%-15.8g ",premod[i]);
		EH.println("");
	  }
  }

  void allocateGeoModel()
  {
	// Legacy Scripps: there's only one sea layer, fixed. All indexing goes from sea-surface down
	// nlayer = # subsea layers.
	// from input files...
	depths = new float[nlayer+2];  // depths[0]= sea surface
	thickf = new float[nlayer+1];  // thickf[0] = sea thickness
	rho = new float[nlayer+1];
	ipen1 = new float[nlayer+1];
	ipen2 = new float[nlayer+1];
  }

  void allocateModel()
  {
	// and associated space for model. Will have np <= nlayers
	pm = new float[np];
	iff = new int[np];
	pwk1 = new float[np];
	pwk2 = new float[np];
	pwk3 = new float[np];
	prewts = new float[np]; 
	premod = new float[np];
	del = new DenseDoubleMatrix2D(np,np);
	ptp = new DenseDoubleMatrix2D(np,np);
	wjtwj = new DenseDoubleMatrix2D(np,np);
	WJdeltaD = new float[np];
	if (iruf==2)
	  {
		del2  = new DenseDoubleMatrix2D(np,np);
	  }
	if (newScripps)
	  {
		// set pointers of new Scripps modules to matching bits
		SMI.DM.nlay1D = nlayer+2;  // nlayer in Legacy = # subsea layers. Assume here NewScripps has 1 atmos + 1 sea layer
		//  LAYER_CHECK

		SMI.OM.nParams = np;
		SMI.OM.premod = premod;
		SMI.OM.prewts = prewts;
	  }
  }

  public void allocWorkspace(int np)
  {
	if (np>0 && SEM.nPsystematics>0)
	  {
		SEM.X_NdeltaD = new float[SEM.nPsystematics];
		SEM.H12 = new DenseDoubleMatrix2D(np,SEM.nPsystematics);
		SEM.H22 = new DenseDoubleMatrix2D(SEM.nPsystematics,SEM.nPsystematics);
	  }

	if (SEM.nPsystematics>0)
	  {
		int Nptot = np + SEM.nPsystematics;
		// inherit any pm[] components
		float[] new_pm = new float[Nptot];
		System.arraycopy(pm,0,new_pm,0,min(np,min(pm.length,new_pm.length)));
		pm = new_pm;
		pwk1 = new float[Nptot];
		pwk2 = new float[Nptot];
		pwk3 = new float[Nptot];
	  }
  }

  void allocateData()
  {
	if (newScripps)
	  {		
		if (modelSystematics)
		  {
			// model computed by fwd engine
			dm_true = SMI.OM.dm.clone();
		  }
		else dm_true = SMI.OM.dm;
	  }
	else
	  {
		d = new float[nd];
		sd = new float[nd];
		dm = new float[nd];	
		if (modelSystematics)
		  {
			// model computed by fwd engine
			dm_true = new float[nd];	
		  }
		else dm_true = dm;

		freqs = new float[nd];
		dp = new DenseFloatMatrix2D(6,nd);   // ala float[6][nd] ... note transposed form!
	  }
	// need wkspace in both modes
	dwk1 = new float[nd];
	dwk3 = new float[nd];
  }

  void inputModel()
  {
	if (newScripps)
	  {
		SMI.readModel();
		// Allocate (np dependent) model
		SMI.reconcileGeoModelVariables();
		allocateModel();
		SMI.reconcileModelVariables();
		imposeOrCheckPriorStability();
	  }
	else
	  {
		// Legacy model
		inputm();
	  }
  }

  void inputm()
  {
	int i,j;
	BufferedReader in=null;
	String s;
	String[] st;
	float[] ppmod;
	LineSplitter LS = new LineSplitter();
	boolean isLegacy = false;;

    if (verbosity > 1) EH.println("reading model file..."+modfil);
	try 
	  {
		in = new BufferedReader(new FileReader(modfil));
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + modfil);
	  }
	try 
	  {
		// read character descriptors      
		for (;;) 
		  {
			s = in.readLine();
			LS.parse(s);
			if (LS.bComment) continue;
			if (LS.sValue.startsWith("OCCAM1DMODEL_1.0")) // String mform = stringPastCol18(in.readLine());
			  {
				isLegacy=true; 
				break;
			  }
			else if (LS.sValue.toLowerCase().startsWith("occam1dmod_1.0"))
			  {
				isLegacy=false;
				if (!newScripps)
				  {
					EH.die_err("Looks like you want to use new Scripps format: use the --NS flag!");
				  }
				break;
			  }
			else
			  {
				EH.die_err("wrong model file type: not OCCAM1DMODEL_1.0 or Occam1DMod_1.0");	
			  }
		  }
		if (isLegacy)
		  {
			String mdescr =  stringPastCol18(in.readLine());
			// read in the number of subsea layers
			nlayer = intPastCol18(in.readLine());
			if (nlayer<=0) 
			  {
				EH.die_err("nlayers must be > 0");
			  }
		  }
		else
		  {
			for (;;) 
			  {
				s = in.readLine();
				LS.parse(s);
				if (LS.bComment) continue;
				if (LS.sCode.equals("#layers"))
				  {
					try
					  {
						nlayer = Integer.parseInt(LS.sValue);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Number of layers not an integer: "+s);
					  }
				  }
			  }
		  }
		if (nlayer<=0) 
		  {
			EH.die_err("nlayers must be > 0");
		  }

		allocateGeoModel();			

		// headers
		for (;;) 
		  {
			s = in.readLine();
			LS.parse(s);
			if (!LS.bComment) break;
		  }

		// tmp storage
		ppmod = new float[nlayer+1];

		// read in the thicknesses etc
		// the first row will be the seawater thickness and resistivity
		for (i=0; i<=nlayer; i++) //i = 0, nlayer
		  {
			if ((s = in.readLine()) != null)
			  {
				st = s.trim().split("\\s+");
				if (st.length >= 5)
				  {
					thickf[i] = Float.parseFloat(st[0]);
					rho[i] = Float.parseFloat(st[1]);
					ipen1[i] = Float.parseFloat(st[2]);
					ppmod[i] = Float.parseFloat(st[3]);
					ipen2[i] = Float.parseFloat(st[4]);	
				  }
				else
				  {
					EH.die_err("Insufficient no. of columns in line "+s+" : of file "+modfil);
				  }	
			  }
			else
			  {
				EH.die_err("Can't read model-file line of file "+modfil+". Expect "+nlayer+" layers???");
			  }
		  }
		in.close();
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + modfil);
	  }

	// note that in this version we will not consider the sea to be free.
	// now count free parameters...
	// if all layers free, np = nlayer
	np = 0;
	for (i = 1; i <= nlayer; i++)
	  {
		if (rho[i] <= 0.0) np++;
	  }
	if (verbosity > 1) EH.println("Expect "+np+" free parameters");

	// Allocate model
	allocateModel();

	// sort out fixed versus free resistivities:
	int nfree = 0;
	for (i = 1; i<=nlayer;i++)
	  {
		if (rho[i] <= 0.0)
		  {
			iff[nfree] = i;
			nfree++;
		  }
	  }

	// c set the prejudical model
	for (i = 0; i<np; i++)
	  {
		premod[i] = ppmod[iff[i]];
	  }

	imposeOrCheckPriorStability();

    ppmod = null;
  }

  void imposeOrCheckPriorStability()
  {
	int i;
	
	if (doBayesianAutoSmoothing)
	  {
		int warnCount=0;
		for (i = 0; i<np; i++)
		  {
			if (ipen2[iff[i]] <= 0.0)
			  {
				// Bayesian prior is improper: required determinant won't exist.
				if (verbosity > 1)
				  {
					if (warnCount<5) 
					  {
						EH.printf("WARNING: Bayesian auto-smoothing with parameter_weight <= 0 is "
								  +"nonsensical: free-layer %d. => Setting weight to broad value %3g%n",
								  (i+1),DEFAULT_SAFE_INV_RHO_PRIOR_STD_DEV);
					  }
					else if (warnCount==5)
					  {
						EH.println("WARNING: Suppressing further Bayesian auto-smoothings warnings...");
					  }
				  }
				ipen2[iff[i]] = DEFAULT_SAFE_INV_RHO_PRIOR_STD_DEV;
				warnCount++;
			  }
		  }
		// Check if auto-smoothing makes any sense to request
		boolean smoothingApplies = false;
		for (i = 0; i<np; i++)
		  {
			if (ipen1[iff[i]] > 0.0) smoothingApplies = true;
		  }
		if (!smoothingApplies)
		  {
			// no coefficients of smoothing, so determinant is independent of mu.
			// Warn and override
			EH.warn("WARNING: Bayesian auto-smoothing requested with no smoothing weights specified. Setting all to 1");
			for (i = 1; i<=nlayer; i++) ipen1[i] = 1.0f;
		  }
	  }
	else
	  {
		boolean needToStabiliseDet = false;
		for (i = 0; i<np; i++)
		  {
			if (ipen2[iff[i]] <= 0.0) {needToStabiliseDet=true; break;}
		  }
		if (needToStabiliseDet)
		  {
			EH.warn("WARNING: Some prior weights are zero. Stabilization applied to generate "
					+"safe |C_prior| in logMML values (does not change prior)");	
		  }
	  } 
  }

  void clean_smoothing()
  {
	// clean smoothing to imply Exponential covariance model if required
	// this implies constant superdiagonal weights.
	int i;

	if (newScripps)
	  {
		EH.warn("FIXME CLEAN THIS FOR NEW SCRIPPS");
		return;
	  }

	if (doBayesianAutoSmoothing && bayesianAutoSmoothingExpStyle)
	  {
		float firstWpEntry=1.0f, firstWEntry=1.0f;

		firstWpEntry = ipen2[iff[0]];
		if (verbosity > 1) EH.printf("WARNING: Bayesian auto-smoothing - setting all prior marginal "
									 +"std-deviations from first weight as %6g%n",1.0/firstWpEntry);
		for (i = 1; i<np; i++)
		  {
			ipen2[iff[i]] = firstWpEntry;
		  }

		firstWEntry = ipen1[iff[0]];
		if (firstWEntry <= 0) firstWEntry=1.0f;
		for (i = 0; i<np; i++) ipen1[iff[i]] = firstWEntry;
		if (verbosity > 1) EH.printf("WARNING: Bayesian auto-smoothing - setting all smoothing weights"
									 +" to constant %6g for exponential covariance%n",firstWEntry);
	  }
	if (zeroSmoothingWeights)
	  {
		if (verbosity > 1) EH.printf("WARNING: Zeroing all smoothing weights");
		for (i = 0; i<np; i++) ipen1[iff[i]] = 0.0f;
	  }
  }

  void diluteErrorBars()
  {
	int i;

	// dilute data-errors by suitable ratio sqrt(nd/Neffective) 
	diluteFactor = (float)sqrt(nd/(float)Ndilute);
	for (i=0;i<nd;i++) sd[i] *= diluteFactor;	  
  }

  void medianFilterReEstimateDataErrors()
  {
	int MED_FILT_SIZE = 2;
	float SD_SAFETY_FACTOR = 2.0f;
	int i,j,k,jmin,jmax,iS=0,type,type_i,n,nOutliers;
	float f,f_i,r,median,y,sigma_y,new_sd;
	cern.colt.list.IntArrayList indxlist;
	int indx[];
	float xmed[];
	ArrayList<DataQCchecker.DataArray>dArrs = null;
	DataQCchecker.DataArray dArr = null;

	PrintWriter out=null;
	String filename=null;
	boolean dumpMedianResults = false;

	if (newScripps)
	  {
		EH.println("Median-filtering not yet implemented for new Scripps kernel/formats");
		return;
	  }


	indxlist = new cern.colt.list.IntArrayList();
	xmed = new float[2*MED_FILT_SIZE+1];

	if (graphDataErrors) dArrs = new ArrayList<DataQCchecker.DataArray>();

	if (verbosity>1) EH.println("Median-filtering data");
	
	f = -Float.MAX_VALUE;
	type = -1;
	nOutliers=0;
	boolean suspect;
	boolean used[] = new boolean[nd];
	for (;;)
	  {
		// find start of new group
	
		// dp[2] = freq, then
		// dp[1] = range, then
		// dp[0] = type.
		// dp[5] = trans. height.
		for (i=0;i<nd;i++)
		  {
			f_i = dp.getQuick(2,i);
			type_i = (int)dp.getQuick(0,i);
			if (((f_i != f) || (type_i!=type)) && !used[i]) 
			  {
				iS = i;
				f = f_i;
				type = type_i;
				break;
			  }
		  }
		indxlist.clear();
		indxlist.trimToSize();
		for (i=iS;i<nd;i++)
		  {
			f_i = dp.getQuick(2,i);
			type_i = (int)dp.getQuick(0,i);
			if ((f_i == f) && (type_i==type) && !used[i])  
			  {
				indxlist.add(i);
				used[i]=true;
			  }
		  }
		if (indxlist.size()==0) break;
		indxlist.trimToSize();
		// now x is a list of data suitable for median filtering etc
		if (graphDataErrors) 
		  {
			dArr = new DataQCchecker.DataArray();
			dArrs.add(dArr);
			dArr.f = f;
			// plot on log scale if useful
			if ((type==1)||(type==3)||((type>=14) && (type<=17)))
			  {
				dArr.logView = true;
			  }
			dArr.code = type;
		  }
		if (dumpMedianResults)
		  {
			filename = String.format("Median_filtered_data.type=%d.f=%.2f.txt",type,f);
			try
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			  } 
			catch (IOException ioe) 
			  {
				throw new RuntimeException("couldn't open " + filename+" for output");
			  }
			out.printf("#no. d median sd_revised d-sd_revised d+sd_revised%n");
		  }
		n = indxlist.size();
		indx = indxlist.elements();
		for (i=0;i<n;i++)
		  {
			jmin=max(0,i-MED_FILT_SIZE);
			jmax=min(n-1,i+MED_FILT_SIZE);

			for(j=jmin,k=0; j<=jmax;j++,k++) xmed[k] = d[indx[j]];
			Arrays.sort(xmed,0,jmax-jmin+1);
			median = xmed[(k)/2]; // central element in sorted array
			r = dp.getQuick(1,indx[i]);
			y = d[indx[i]];
			sigma_y = sd[indx[i]];
			suspect=false;
			if (abs(y - median) > SD_SAFETY_FACTOR*sigma_y)
			  {
				// suspected outlier
				new_sd =  abs(y - median);
				suspect=true;
				nOutliers++;
			  }
			else 
			  {
				// retain existing error bar
				new_sd = sigma_y;
			  }
			if (dumpMedianResults) 
			  {
				// out.printf("%d %g %g %g %g %g%n",i,y,median,new_sd,y-new_sd,y+new_sd);
				out.printf("%d %g 0%n",i,median);
			  }
			if (graphDataErrors) 
			  {
				if (!suspect)
				  {
					dArr.X.add(r);
					dArr.Y.add(y);
					dArr.E.add(sigma_y);
				  }
				else
				  {
					dArr.Xo.add(r);
					dArr.Yo.add(y);
					dArr.Eo.add(new_sd);
				  }
			  }
			// overwrite existing sd
			sd[indx[i]] = new_sd;
		  }
		if (graphDataErrors) 
		  {
			dArr.trim(); 
		  }
		if (dumpMedianResults) 
		  {
			out.println("");
			for (i=0;i<n;i++) out.printf("%d %g %g%n",i,d[indx[i]],sd[indx[i]]);
		  }
		if (dumpMedianResults)
		  {
			out.close();
			EH.println("Wrote Median filtering summary file %s",filename);
		  }
	  }
	if (verbosity>1) EH.println("Median-filtering: detected %d outliers",nOutliers);
	if (graphDataErrors) 
	  {
		DataQCchecker DQG = new DataQCchecker(); 
		DQG.plotData(dArrs);
	  }
  }

  /**
	 Here we define the sill of the prior exponential covariance used ifbayesianAutoSmoothingExpStyle
	 applies. for layers i=0...np-1, the correlation structure is defined as
	 C_ij = sigma_p_sq exp(-alpha |i-j|)
	 We use the first prejudice weight in the input file to define the inverse sill
	 sigma_p_sq = 1/ipen2[iff[0]] ^ 2
  */
  double getExpPriorCovarianceSill()
  {
	double firstWpEntry;

	firstWpEntry  = ipen2[iff[0]];

	return 1/(firstWpEntry*firstWpEntry);

  }

  /** 
	  map smoothing weight mu to effective exponent alpha
  */
  double getExpMinus2Alpha(double mu)
  {
	double sigma_p_sq, sinh_alpha, alpha, exp_m2alpha;

	sigma_p_sq = getExpPriorCovarianceSill();
	sinh_alpha = 1/(2*mu*sigma_p_sq);
	alpha = Math.log(sinh_alpha + sqrt(sinh_alpha*sinh_alpha+1)); // arsinh()
	exp_m2alpha = Math.exp(-2*alpha);

	return exp_m2alpha;
  }

  void updatePrewtsForExpStyleSmoothing(double mu)
  {
	double sigma_p_sq, sinh_alpha, alpha, exp_m2alpha;
	double e1,e2,e3;
	double Wp1,Wp2;

	sigma_p_sq = getExpPriorCovarianceSill();
	exp_m2alpha = getExpMinus2Alpha(mu);
	// entries in tridiagonal inverse covariance
	e1 = 1/(1-exp_m2alpha);
	e2 = -mu*sigma_p_sq;
	e3 = (1+exp_m2alpha)/(1-exp_m2alpha);

	// diagonal entries of additional weights, these must be in the order
	// {Wp1,Wp2,Wp2,.....Wp2,Wp1}
	// to give standard exponential covariance tridiagonal inverse form
	if (exp_m2alpha > 0.99999)
	  {
		alpha = -0.5*Math.log(exp_m2alpha);
		Wp1 = sqrt(1.0/(sigma_p_sq *(2.0-alpha)));
		Wp2 = sqrt(alpha/(sigma_p_sq *(2.0-alpha)));
	  }
	else
	  {
		Wp1 = e1/sigma_p_sq - mu;
		Wp1 = (Wp1>0) ? sqrt(Wp1) : 0.0;
		Wp2 = e3/sigma_p_sq - 2*mu;
		// e3/sigma_p_sq - 2*mu should theoretically be positive, but to guard against fp error...
		Wp2 = (Wp2>0) ? sqrt(Wp2) : 0;
	  }
	// update prewts vector to reflect this: {Wp1,Wp2,Wp2,.....Wp2,Wp1}
	prewts[0] = (float)Wp1;
	prewts[np-1] = (float)Wp1;
	for (int i=1; i<(np-1); i++) prewts[i] = (float)Wp2;
  }

  void outputModel(float[] thickf_best)
  {
	int i;
	PrintWriter out;

	String filename="model_suggestion.mod";
	try
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }

	int Nl = thickf_best.length;

	// Print best model file
	out.printf("FORMAT:           OCCAM1DMODEL_1.0%n"
			   +"MODEL NAME:       Autosplit best model file%n"
			   +"NUM LAYERS:       %d%n"
			   +"THICKNESS     RESISTIVITY  PENALTY  PREJUDICE  PENALTY%n",Nl-1);
	
	out.printf("%12.6g %12.2g %12.2g %12.2g %12.2g%n",
			   thickf_best[0],rho[0],
			   0.0,0.0,1.0);

	for (i=1; i<Nl; i++)
	  {
		out.printf("%12.6g %12.2g %12.2g %12.2g %12.2g%n",
				   thickf_best[i],-1.0,
				   1.0,0.0,1.0);
	  }
	out.close();
	if (verbosity > 1) EH.printf("Wrote split-suggest model %s%n",filename);
  } 

  /*
	Parse cmdline options to make noise
  */
  void parseSynthNoiseSpecs()
  {
	int i,N;
	String spec;

	nSpec = new NoiseSpec(EH,syntheticSpecFrac,syntheticSpecFloor);
  }

  /**
	 Compute "synthetic data set" based on existing model. 

	 Overwrite existing data d[0...nd-1] if overwriteData == true.

	 Write to file synthFilename if synthFilename!=null

	 If (scaleExistingErrors > 0) will scale the existing sd_orig[] array.
	 If sd_orig is null, this method will create it from existing sd.
	 Values are scaled by supplied scalar
	 scaleExistingErrors. I.e. this method corrupts sd[]

	 If (scaleExistingErrors <=0), will use distributions specified by
	 globals synthetic_error and synthetic_error_noise_floor
  */
  void createSyntheticFromCurrentModel(boolean overwriteData, double scaleExistingErrors, String synthFilename)
  {
	int i,ip,j,k,type,typep;
	float frac_error,y,freq,s;
	float yR,yI,syR,syI,sA;
	float sSQR;
	boolean writePair = false;
	boolean rescaleInSituErrors;
	float[] d_new, _dm;
	PrintWriter out;

	//if (newScripps) EH.die_err("Synthetics from new Scripps kernel not yet implemented.");

	if (overwriteData) 
	  {
		d_new = d;
	  }
	else
	  {
		d_new = new float[nd];
	  }

	// frac_error = 0.01f * synthetic_error;

	rescaleInSituErrors = (scaleExistingErrors > 0);
	if (rescaleInSituErrors)
	  {
		// we will work with original sd data, make sure we have a copy, if not already
		if (sd_orig==null) sd_orig = (float[])sd.clone();
	  }

	// Write out actual fields, pre systematic adding, if invoked (not likely)
	_dm = modelSystematics ? dm_true : dm;

	if (!newScripps)
	  {
		// Legacy
		sSQR = (float)MathUtil.sqr(Math.toRadians(nSpec.phase_error_E));
		for (i=0;i<nd;i++)
		  {
			writePair = false;

			type = (int)dp.getQuick(0,i);
			freq = dp.getQuick(2,i);
			ip=(i+1);
			if (ip<nd)
			  {
				typep = (int)dp.getQuick(0,ip);
				if (((type==14 && typep==15) || (type==16 && typep==17))  // possible Re-Im pair
					&& (dp.getQuick(1,i) == dp.getQuick(1,ip))            // same range
					&& (freq == dp.getQuick(2,ip))                        // freq
					&& (dp.getQuick(3,i) == dp.getQuick(3,ip))            // az
					&& (dp.getQuick(4,i) == dp.getQuick(4,ip))            // rot
					&& (dp.getQuick(5,i) == dp.getQuick(5,ip))            // theight
					)
				  {
					// A Re-Im pair
					writePair = true;
					yR = _dm[i];
					yI = _dm[ip];
					frac_error = 0.01f * nSpec.synthetic_error_E;   // Legacy: these are E fields
					syR = (float)sqrt(frac_error*frac_error*yR*yR + sSQR*yI*yI);
					syI = (float)sqrt(frac_error*frac_error*yI*yI + sSQR*yR*yR);
					if (!rescaleInSituErrors)
					  {
						sd[i] = syR;  // Real
						sd[ip] = syI; // Im
					  }
				  }
			  } 
			if (!writePair && !rescaleInSituErrors)
			  {
				if (type==2 || type==4)
				  {
					// an E phase, degrees. 
					// OLD sd[i] = frac_error * SQRT2_180_ON_PI;
					sd[i] = nSpec.phase_error_E;
				  }
				else
				  {
					frac_error = 0.01f * nSpec.synthetic_error_E;   // Legacy: these are E fields
					sd[i] = frac_error * abs(_dm[i]);
				  }			
			  }
			for (j=i; j<=(writePair ? ip : i); j++)       // one or 2 points
			  {
				if (rescaleInSituErrors) 
				  {
					// scale original copy
					sd[j] =  (float)(scaleExistingErrors * sd_orig[j]);
				  }
				else
				  {
					// apply requested recipe: Legacy E fields
					if (sd[j] < nSpec.synthetic_error_noise_floor_E) sd[j] = nSpec.synthetic_error_noise_floor_E;
				  }
				for(;;)
				  {
					// correct version
					s = (float)(RD.nextGaussian());				
					d_new[j] = _dm[j] + s * sd[j];  // resample data
					// don't allow -ve |E| fields
					if ((d_new[j]>=0) || !(type==1 || type==3)) break;
				  }
			  }
			if (writePair) i++;
		  }


		if (synthFilename != null)
		  {
			// Dump synthetic to new file
			try 
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(synthFilename)));
			  } 
			catch (IOException ioe)  
			  {
				throw new RuntimeException("couldn't open " + synthFilename + " for output");
			  }
			out.printf("FORMAT:           OCCAM1DDATA_1.0%n"
					   +"DESCRIPTION:      synthetic output%n"
					   +"NUM DATA:         %d%n"
					   +"Type Datum Error range/period frq az rot theight%n",nd);

			for (i=0; i<nd; i++)
			  {
				out.printf("%5d ",(int)dp.getQuick(0,i)); // type
				out.printf("%15.5g ", d_new[i]);           // New datum
				out.printf("%15.5g ",sd[i]);              // Error
				for (k=1;k<=5;k++) out.printf("%15.5g ",dp.getQuick(k,i)); // range||period frq az rot theight
				out.printf("%n");
			  }
			out.close();
			if (verbosity > 1) EH.println("Wrote synthetic file " + synthFilename);
		  }
	  }	
	else
	  {
		// New Scripps format
		DenseIntMatrix2D _dp = SMI.OM.dp_trans; 
		int iFreq,iTx,iRx,ipFreq,ipTx,ipRx;
		float freqp;
		boolean isAnAbsField=false;
		boolean isAnEField=false;
		for (i=0;i<nd;i++)
		  {
			writePair = false;

			type = _dp.getQuick(0,i);
			iFreq = _dp.getQuick(1,i);
			iTx = _dp.getQuick(2,i);
			iRx = _dp.getQuick(3,i);
			freq = (float)SMI.CM.fTxIn[iFreq-1];

			isAnEField = ((type>=1 && type<=6) || type==21 || type==23 || type==25);

			ip=(i+1);  // next datum
			if (ip<nd)
			  {
				typep = _dp.getQuick(0,ip);
				ipFreq = _dp.getQuick(1,ip);
				ipTx = _dp.getQuick(2,ip);
				ipRx = _dp.getQuick(3,ip);
				freqp =(float)SMI.CM.fTxIn[ipFreq-1];

				if (   (iRx == ipRx)                                        // same Rx
					&& (iTx == ipTx)                                        // same Tx					   
					&& (freq == freqp)                                      // same freq
					&& (  (type==1 && typep==2) || (type==3 && typep==4) || (type==5 && typep==6)             // E: possible Re-Im pair
						  ||(type==11 && typep==12) || (type==13 && typep==14) || (type==15 && typep==16)))   // B: possible Re-Im pair
				  {
					// A Re-Im pair
					writePair = true;
					yR = _dm[i];
					yI = _dm[ip];
					if (isAnEField)
					  {
						frac_error = 0.01f * nSpec.synthetic_error_E;
						sSQR = (float)MathUtil.sqr(Math.toRadians(nSpec.phase_error_E));
					  }
					else
					  {
						frac_error = 0.01f * nSpec.synthetic_error_B;
						sSQR = (float)MathUtil.sqr(Math.toRadians(nSpec.phase_error_B));
					  }
					syR = (float)sqrt(frac_error*frac_error*yR*yR + sSQR*yI*yI);
					syI = (float)sqrt(frac_error*frac_error*yI*yI + sSQR*yR*yR);
					if (!rescaleInSituErrors)
					  {
						sd[i] = syR;  // Real
						sd[ip] = syI; // Im
					  }
				  }
			  } 
			if (!writePair && !rescaleInSituErrors)
			  {
				if (type==22 || type==24 || type==26)         // E phase
				  {
					// an E phase, degrees.
					// OLD sd[i] = frac_error * SQRT2_180_ON_PI;
					sd[i] = nSpec.phase_error_E;
				  }
				else if (type==32 || type==34 || type==36)    // B phase
				  {
					sd[i] = nSpec.phase_error_B;
				  }
				else
				  {
					// An E or B, not a phase
					frac_error = 0.01f * ( isAnEField ? nSpec.synthetic_error_E : nSpec.synthetic_error_B);
					sd[i] = frac_error * abs(_dm[i]);
				  }			
			  }
			for (j=i; j<=(writePair ? ip : i); j++)       // one or 2 points
			  {
				if (rescaleInSituErrors) 
				  {
					// scale original copy
					sd[j] =  (float)(scaleExistingErrors * sd_orig[j]);
				  }
				else
				  {
					// apply requested recipe
					if (isAnEField) {if (sd[j] < nSpec.synthetic_error_noise_floor_E) sd[j] = nSpec.synthetic_error_noise_floor_E;}
					else {if (sd[j] < nSpec.synthetic_error_noise_floor_B) sd[j] = nSpec.synthetic_error_noise_floor_B;}
				  }
				isAnAbsField = (type==21 || type==23 || type==25         // |E| types
								|| type==31 || type==33 || type==35);    // |B| types
				for(;;)
				  {
					s = (float)(RD.nextGaussian());				
					d_new[j] = _dm[j] + s * sd[j];  // resample data
					// don't allow -ve |E| or |B| fields
					if ((d_new[j]>=0) || !isAnAbsField) break;
				  }
			  }
			if (writePair) i++;
		  }

 		if (synthFilename != null)
		  {
			// Dump synthetic to new file
			try 
			  {
				out = new PrintWriter(new BufferedWriter(new FileWriter(synthFilename)));
			  } 
			catch (IOException ioe)  
			  {
				throw new RuntimeException("couldn't open " + synthFilename + " for output");
			  }
			
			// New scripps "response format": clone of KKey's usage - EMDATA_1.1"
			out.printf("Format:    EMDATA_1.1%n");
			out.printf("# Transmitters:      %d%n! X Y Z Azimuth Dip%n",   SMI.CM.nTxIn);
			for (i = 0;i<SMI.CM.nTxIn;i++)
			  {
				out.printf("%15.5g %15.5g %15.5g %15.5g %15.5g%n",
						   SMI.CM.xTxIn[i], SMI.CM.yTxIn[i], SMI.CM.zTxIn[i], SMI.CM.azimuthTxIn[i], SMI.CM.dipTxIn[i]);
			  }
			out.printf("# Frequencies: %d%n",SMI.CM.nFreqIn);
			for (i = 0; i<SMI.CM.nFreqIn; i++) out.printf("%15.5g%n",SMI.CM.fTxIn[i]);			
			out.printf("# Receivers: %d%n",   SMI.CM.nRxIn);
			out.printf("! X Y Z Theta Alpha Beta%n");
			for (i = 0;i<SMI.CM.nRxIn;i++)
			  {
				out.printf("%15.5g %15.5g %15.5g %15.5g %15.5g %15.5g%n",SMI.CM.xRxIn[i],SMI.CM.yRxIn[i],SMI.CM.zRxIn[i], 
						   SMI.CM.ThetaRxIn[i],SMI.CM.AlphaRxIn[i],SMI.CM.BetaRxIn[i]);
			  }    
			out.printf("# Data: %d%n", nd);
			out.printf("%8s %8s %8s %8s %15s %15s %15s %15s%n","! Type","Freq#","Tx#","Rx#","Data","StdError","Response","Residual");
			for (i=0;i<nd;i++)
			  {
				for (j=0; j<4; j++) out.printf("%8d ",SMI.OM.dp_trans.getQuick(j,i));
				s = (sd[i] <=0)? 0 : ((d_new[i] - _dm[i])/sd[i]);
				out.printf("%15.5g %15.5g %15.5g %15.5g%n",d_new[i],sd[i],_dm[i],s);
			  } 
			out.close();
			if (verbosity > 1) EH.println("Wrote synthetic file " + synthFilename);
		  }
	  }
  }

  void itrout(String extraText,int maxitr,float rlast,float tobt,int task_code, int i_multistart)      
  {
	// task_code = 0 for usual dump during optimisation
	// task_code = 1 for startup dump
	// task_code = 2 for final dump with special file name _final

	// local variables
	int i, j, i_free;
	float v=0.0f;
	String filename,num=null,start=null;
	PrintWriter out;
 
	//-----------------------------------------------------------------------
	// Quiet unless final dumps
	if (quietIO && (task_code != 2)) return;


	// create file name suffix
	num = String.format("%04d",nit); 
	// open iteration file
	if (doMultistart) 
	  {
        start= String.format("start_%04d",i_multistart);
	  } 
	else if (doNoiseCooling || doRelaxedSmoothing)
	  {
		start= String.format("cool_%04d",i_multistart);
	  }
	  
	if (task_code == 2) 
	  {
		filename = doMultistart ? ("iter_final."+start) : "iter_final";
	  }
	else
	  {
		if (doMultistart) 
		  {
            filename = "iter."+start+"."+num;
		  }
		else
		  {
            filename = "iter"+num;
		  }
      }
	if (extraText != null) filename = filename.concat(extraText);
	current_iter_filename = filename;
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		// we could be in the middle of an iteration, so don't die on error, but
		// send a message and forge on regardless
        if (verbosity > 0) EH.println("Error opening iteration file..");
		return;
		// throw new RuntimeException("couldn't open " + filename+" for output");
	  }
	String datetm = date.toString(); 
	// 
	// write stuff      
	out.println("format:           " + itform);
	out.println("description:      " + descr);
	out.println("model file:       " + modfil);
	out.println("data file:        " + datfil);
	out.println("date/time:        " + datetm);
	out.printf("max iter:         %4d%n",maxitr);
	out.printf("req tol:          %-6.3f%n",tolreq);
	out.printf("iruf:             %-4d%n",iruf);
	out.printf("debug level:      %-4d%n",idebug);
	out.printf("iteration:        %-6d%n",nit);
	out.printf("pmu:              %-15.7g%n",pmu);
	out.printf("rlast:            %-15.7g%n",rlast);
	out.printf("tlast:            %-15.7g%n",tobt);
	out.printf("ifftol:           %-3d%n",ifftol);
	out.printf("no. parms:        %-4d%n",np);
	for (i=0;i<np;i++) out.printf("%-17.10g ",pm[i]);
	out.printf("%n");
	out.close();

	// JG graph files for resistivity profile
	if (dumpPrfl) 
	  {	
		// Dump only iterate-profiles
		if (task_code != 2) 
		  {
			if (doMultistart) 
			  {
				filename = "profile."+start+"."+num;
			  }
			else
			  {
				filename = "profile"+num;
			  }
			if (extraText != null) filename = filename.concat(extraText);
			writeProfileSnapshot(filename);
		  }
	  }


	// open response file
	// JG allow for init call
	if ((nit!=0) || (task_code == 1)) 
	  { 
		if (task_code == 2) 
		  {
			filename = doMultistart ? ("resp_final."+start) : "resp_final";
		  }
		else
		  {
            filename = "resp"+num;
		  }
		if (extraText != null) filename = filename.concat(extraText);
		try 
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename+" for output");
		  }
		// write out actual raw CSEM fields (before systematics, if used)
		float[] _dm = modelSystematics ? dm_true : dm;
		if (!newScripps)
		  {
			for (i=0;i<nd;i++)
			  {
				out.printf("%5d ",(int)dp.getQuick(0,i));
				for (j=1;j<ndp;j++) out.printf("%15.5g ",dp.getQuick(j,i));
				out.printf("%15.5g %15.5g %15.5g%n",d[i],_dm[i],(d[i] - _dm[i])/sd[i]);
			  }
		  }
		else
		  {
			// New scripps "response format": clone of KKey's usage -  EMResp_1.1"
			out.printf("Format:    EMResp_1.1%n");
			out.printf("# Transmitters:      %d%n! X Y Z Azimuth Dip%n",   SMI.CM.nTxIn);
			for (i = 0;i<SMI.CM.nTxIn;i++)
			  {
				out.printf("%15.5g %15.5g %15.5g %15.5g %15.5g%n",
						   SMI.CM.xTxIn[i], SMI.CM.yTxIn[i], SMI.CM.zTxIn[i], SMI.CM.azimuthTxIn[i], SMI.CM.dipTxIn[i]);
			  }
			out.printf("# Frequencies: %d%n",SMI.CM.nFreqIn);
			for (i = 0; i<SMI.CM.nFreqIn; i++) out.printf("%15.5g%n",SMI.CM.fTxIn[i]);			
			out.printf("# Receivers: %d%n",   SMI.CM.nRxIn);
			out.printf("! X Y Z Theta Alpha Beta%n");
			for (i = 0;i<SMI.CM.nRxIn;i++)
			  {
				out.printf("%15.5g %15.5g %15.5g %15.5g %15.5g %15.5g%n",SMI.CM.xRxIn[i],SMI.CM.yRxIn[i],SMI.CM.zRxIn[i], 
						   SMI.CM.ThetaRxIn[i],SMI.CM.AlphaRxIn[i],SMI.CM.BetaRxIn[i]);
			  }    
			out.printf("# Data: %d%n", nd);
			out.printf("%8s %8s %8s %8s %15s %15s %15s %15s%n","! Type","Freq#","Tx#","Rx#","Data","StdError","Response","Residual");
			for (i=0;i<nd;i++)
			  {
				for (j=0; j<4; j++) out.printf("%8d ",SMI.OM.dp_trans.getQuick(j,i));
				out.printf("%15.5g %15.5g %15.5g %15.5g%n",d[i],sd[i],_dm[i],(d[i] - _dm[i])/sd[i]);
			  }
		  }
		out.close();
		if (verbosity > 1) EH.println("Wrote response file "+filename);
		current_resp_filename = filename;
	  }
  }

  
  private void writeProfileSnapshot(String filename)
  {
	int i,i_free;
	float v;

	PrintWriter out;

	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }
	// write depths verses log10 resistivities as a plottable line-graph
	for (i = 0, i_free=0; i<(nlayer+1); i++)
	  {
		if (rho[i] <= 0.0f)
		  {
			v = pm[i_free];
			i_free++;
		  }
		else 
		  {
			// v = (i==0) ? flog10(rho[i]) : rho[i]; 
			v = flog10(rho[i]); // +ve numbers are in absolute units, bith new and legacy Scripps I think
		  }
		out.printf("%15.7g %15.7g%n",depths[i],v);
		out.printf("%15.7g %15.7g%n",depths[i+1],v);
	  }
	out.close();

	if (verbosity > 1) EH.println("Wrote profile file "+filename);
  }

  private void initialise()
  {
	int i,Nptot;

	Nptot = np;
	if (modelSystematics)
	  {
		Nptot += SEM.nPsystematics;
	  } 
	
	// gradient space
	grad = new float[Nptot];


	// Sensitivity allocation. Transposed form since this comes from fortran
    // WJ_trans = new float[np][nd];  // WJ = float[nd][np]; 
	WJ_trans = new DenseFloatMatrix2D(np,nd);
	if (newScripps)
	  {
		SMI.OM.wj_trans = WJ_trans;
	  }

	if (useLMstrategy)
	  {
		LMinfo = new LevMarInfo(Nptot);
	  }


	// precompute depths
    depths[0] = 0f;
	for (i=1 ; i < nlayer+2 ; i++)
	  {
		depths[i] = depths[i-1]+thickf[i-1];
	  }
	// Bounded constraint variables
	// TODO: we should be able to read these in on a layer basis etc. 
	// Merely requires suitable I/O, formats etc.
	// For now, set same bounds for all free layers.
	if (useBoundConstraints) 
	  {
		float lower_bound,upper_bound;
		act_PN = new byte[np];
		m_lower = new float[np];  // fill with lower bounds
		m_upper = new float[np];  // fill with upper bounds

		lower_bound = -0.1f; // this is log_{10}(0.8 Ohm-m)..a pretty good lwr bound based on Hashin-Shtrikman
		upper_bound = 4.0f;  // reasonable for most CSEM cases
		if (option_bounds_lower != -Float.MAX_VALUE) lower_bound = option_bounds_lower;
		if (option_bounds_upper != Float.MAX_VALUE) upper_bound = option_bounds_upper;
		for (i=0;i<np;i++)
		  {
			m_lower[i] = lower_bound;
			m_upper[i] = upper_bound;
		  }
		// skin depth
		minSkinDepth = (float)sqrt(pow(10,lower_bound)/(PI*freqMin*(PI*4e-7)));
	  }
	else
	  {
		// estimate skin depth by very-porous shale (in SI) value: 0.8 Ohm.m
		minSkinDepth = (float)sqrt(0.8/(PI*freqMin * PI*4.0e-7));
	  }
	if (verbosity > 1) EH.printf("Estimated min skin-depth = %4g (freq_min = %3f Hz)%n",minSkinDepth,freqMin);
  }

  private void conformStartupToBounds()
  {
	if (useBoundConstraints) 
	  {
		int i;
		for (i=0;i<np;i++)
		  {
			pm[i] = min(max(pm[i],m_lower[i]),m_upper[i]);
		  }
	  }
  }
  
  private void overWriteWithModelStudyValues()
  {
	if (modelStudyList!=null)
	  {
		int i,j;
		for (i=0; i<modelStudyList.length; i++)
		  {
			j = modelStudyList[i].layer;
			if (j<=0 || j >nlayer)
			  {
				EH.printf("Model study specifies layer not present in model: %d%n",j);
				EH.die();
			  }
			thickf[j] = modelStudyList[i].current_value;
		  }
		// re-compute depths
		depths[0] = 0f;
		for (i=1 ; i < nlayer+2 ; i++)
		  {
			depths[i] = depths[i-1]+thickf[i-1];
		  }
		if (verbosity>1)
		  {
			EH.printf("DEPTHS: ");
			for (i=0;i<nlayer+2;i++) EH.printf("%-10.2f ",depths[i]);
			EH.println("");
			EH.printf("THICKF: ");
			for (i=0;i<nlayer+1;i++) EH.printf("%-10.2f ",thickf[i]);
			EH.println("");
			if (false)
			  {
				EH.println("Equivalent new model file contents:");
				EH.println("Thickness rho smooth-wt <log10(rho)>_prior prior_wt (prior-std-dev)");
				for (i=1; i<=nlayer; i++)
				  {
					EH.println("%10.6g -1 0       %10.6g         %10.6g    (%-10.6g)",
							   thickf[i],premod[i-1], ipen2[i],1/ipen2[i]);
				  }
			  }
		  }
	  }
  }

  void read_startup_and_model()
  {
	int i,j, np_startup=0;
	String s;
	String st[];
    BufferedReader in;
	String filename;
	boolean isOCCAMITER_FLEX=false;

	filename = (altStartFile != null) ? altStartFile : "startup";

	if (verbosity > 1) EH.println("reading startup file: "+filename); 
	try 
	  {
		in = new BufferedReader(new FileReader(filename));
		in.mark(10000);
		// check if it's a "new scripps" file
		try 
		  {
			for(;;)
			  {
				if ((s = in.readLine()) == null) 
				  {
					break;
				  }
				else 
				  {
					s = s.trim().toLowerCase();
					if (s.startsWith("format") && s.contains("occamiter_flex"))
					  {
						isOCCAMITER_FLEX = true;
						if (verbosity > 1) EH.println("Using occamiter_flex format.."); 
					  }
					break;
				  }
			  }
		  }
		catch (IOException e)
		  {
			EH.die_err("Can't read from file "+filename);
		  }
		in.reset();
		if (!isOCCAMITER_FLEX)
		  {
			s = getNextNonCommentLine(in, filename);
			itform = stringPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			descr = stringPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			modfil = stringPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			datfil = stringPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			s = stringPastCol18(s); // ignore date
			s = getNextNonCommentLine(in, filename);
			maxitr = intPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			tolreq = floatPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			iruf = intPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			idebug = intPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			nit = intPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			pmu = floatPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			rlast = floatPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			tobt_in = floatPastCol18(s);	   
			s = getNextNonCommentLine(in, filename);
			ifftol = intPastCol18(s);	
		  }
		else
		  {
			// isOCCAMITER_FLEX format
			// code should probbaly be run with --NS
			if (!newScripps)
			  {
				EH.die_err("Startup file indicates new Scripps formats. I think you need to run with --NS");
			  }

			String  sLine, sCode, sValue; // These are for reading the lines of the file
			LineSplitter PT = new LineSplitter();
			for(;;)
			  {
				// Get the next code/value pair
				// ParseCode forces the code portion to be all lowercase with no
				//   padding and ending colon stripped off.  User comments are 
				//   stripped from the value portion.
				try
				  {
					if ((s = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					break; //skip reading
				  }
				PT.parse(s);
				if (PT.bComment) continue;
				sCode = PT.sCode.toLowerCase();
				sValue = PT.sValue;
        
				// What do we have?
				if (sCode.equals("format"))
				  {
					itform = sValue;
				  }
				else if (sCode.equals("description"))
				  {
					descr = sValue; 
				  }
				else if (sCode.equals("model file"))
				  {
					modfil = sValue; 
				  }
				else if (sCode.equals("data file"))
				  {
					datfil = sValue; 
				  }
				else if (Util.matches(sCode,"max iter","iterations to run"))
				  {
					maxitr = parseInt(sValue); 
				  }
				else if (Util.matches(sCode,"req tol","target misfit"))
				  {
					tolreq = parseFloat(sValue); 
				  }
				else if (Util.matches(sCode,"iruf","roughness type"))
				  {
					iruf= parseInt(sValue); 
				  }
				else if (sCode.equals("debug level"))
				  {
					idebug = parseInt(sValue); 
				  }
				else if (Util.matches(sCode,"pmu","lagrange value"))
				  {
					pmu = parseFloat(sValue); 
				  }
				else if (Util.matches(sCode,"rlast","roughness value"))
				  {
					rlast = parseFloat(sValue); 
				  }
				else if (Util.matches(sCode,"tlast","misfit value"))
				  {
					tobt_in = parseFloat(sValue); 
				  }
				else if (Util.matches(sCode,"no. parms", "param count"))
				  {
					np_startup = parseInt(sValue); 
					break; // next stuff must be parameter startup values
				  }
				// skipping: (these make no sense in deliveryCSEM)
				// 'stepsize cut count'
				// 'delta misfit limit'
				// 'model limits'
				// 'model value steps'
				// 'modify roughness'
			  }
		  }
		// Now we know filenames etc, 
		// input the fixed model characteristics, get np, allocate storage etc
		if (abs(idebug) == 1) if (verbosity > 1) EH.println("reading model..");   
		inputModel();
	
		if (!isOCCAMITER_FLEX)
		  {		
			// continue parsing startup file....
			s = getNextNonCommentLine(in, filename);
			np_startup = intPastCol18(s);
		  }

		// Check implied model size
		// npp = np;
		if (np_startup != np) 
		  {
			EH.die_err("inconsistent number of free parameters in startup file c.f. model file");
		  }

		// Allocations should be in place for this...

		j=0;
		LOOP:
		while ((s = in.readLine()) != null)
		  {
			s = s.trim();
			if (s.startsWith("#")) continue; // skip comments
			st = s.split("\\s+");
			for (i=0; i<st.length; i++)
			  {
				try
				  {
					pm[j++] = Float.parseFloat(st[i]);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Line containing "+s+" contains non-number. "
							   +"Or insufficient startup parameters");
				  }
				if (j==np) break LOOP;
			  }
		  }
		if (j<(np)) EH.die_err("Insufficient startup parameters");

		// JG: Extra options
		if ((s = in.readLine()) !=null)
		  {
			st = s.trim().split("\\s+");
			if (st.length>=2)
			  {
				if (st[0].equals("XSECTIONS:") && st[1].equals("1")) doXsections = true;
			  }
		  }
		if ((s = in.readLine()) !=null) 
		  {
			st = s.trim().split("\\s+");
			if (st.length>=2)
			  {
				if (st[0].equals("PROFILES:") && st[1].equals("1")) dumpPrfl = true;				
			  }
		  }
		in.close();
	  }
	catch (FileNotFoundException e)
	  {
		EH.die_err("Can't find file "+filename);
	  }
	catch (IOException e)
	  {
		EH.die_err("Can't read from file "+filename);
	  }
  }


  /**
	 Speculative MCMC sampler of posterior
  */
  void GenerateMCMCPosterior()
  {
	int i,i_r,nModes,nr,j=0,k,i_acc,i_acc_loc,i_acc_sw,i_acc_sw_try,i_acc_hop,i_acc_hop_try,info;
	int N_burn, N_loop,Nptot;
	int TEST_WIN_SIZE;
	int Nmodes,jmode,jmode_old;
	double acc_rate;
	float  tmp,rms=1f,ms=1f,obj,obj_old,ms_old=1f,msprej_old=0f;
	String filename,filename_prof;
	boolean  try_swap=false,add_pen,set_lwr_bound;
	boolean  write_MC_resp,write_MC_profiles;
	float  mu,log10mu, ex,m_tmp,R,LWR_BOUND_M,LWR_BOUND;
	float  BAD_TRY_PENALTY = 100;
	float  PERMUTE_PROB;
	float  modeExchangeProb, MODE_EXCHANGE_PROB;
	float  u,rho_jp1,rho_j,m_new_jp1,m_new_j;
	float JFAC=0.0f;
	DoubleMatrix2D  a,Rtri;
	float[] delta_m=null,mopt;
	float[]  m_last_accepted;
	float[] stdv; 
	int flippableLayers[]=null;
	float  Rtot,Rtot_new,Rmin,mv,err,y,z,scale;
	PrintWriter out=null,out_prof=null;
	LocalMode m0,mi;
	boolean proposalLegal;
	boolean decimate = true;
	int decimation_rate = 50;
	boolean calibrating = true;
	double TARGET_MAX_ACC_RATE,TARGET_MIN_ACC_RATE;
	double min,prob,pmin;
	float chsqNLPmin=0f,errorScale=1.0f;
	double[] p_mode_prop = null, p_mode_CDF=null;

	// MD variables...
	int N_MD=50;
	double[]  p_MD=null,m_MD=null;
	float[] grad_MD = null;
	DenseDoubleMatrix2D a_MD = null;
	float deltaT_MD=0f;

	// debugging
	int   ii,jj,nrs,i_xs;
	float  mj_old,mjp1_old;
	float  mmin,mmax,mCav,mRav;

	Nptot = np + SEM.nPsystematics;

	boolean adjustNoiseLevelToEBestimate, doNoiseEstimation_old  = doNoiseEstimation;

	filename = (samplesFilename == null) ? "MCMCsamples.txt" : samplesFilename;
	filename_prof = "MCMCprofiles.txt";

	// Default setting of minimum chisq (data and prior bits, no determinants)
	m0 = optima.get(0);
	chsqNLPmin = nd * m0.rmsMisfit*m0.rmsMisfit + m0.negLogNoSmoothingPrior; 

	adjustNoiseLevelToEBestimate = true;
	if (adjustNoiseLevelToEBestimate)
	  {
		// chunk to adjust noise level to Emp Bayes estimate
		if (doNoiseEstimation)
		  {
			// reset global error scale to EB estimate, fix the noise, and proceed
			global_error_scale = m0.rmsMisfit;
			errorScale = global_error_scale;
			doNoiseEstimation = false;
			// get lower bound on twiceNegLogPosterior
			// this is true, if no smoothing etc
			chsqNLPmin = (float) ( (nd-Nptot) + m0.negLogNoSmoothingPrior); 
			if (verbosity > 1) 
			  {
				EH.printf("MCMC: Noise level adjusted to %-6g, minimum_chisq_total= %-6g%n",
						  global_error_scale,chsqNLPmin);
			  }
		  }
	  }

	TEST_WIN_SIZE = 100;
	TARGET_MIN_ACC_RATE=0.1;
	TARGET_MAX_ACC_RATE=0.25;

	MODE_EXCHANGE_PROB = 0.5f;
	Nmodes = optima.size();
	modeExchangeProb = Nmodes>1 ? MODE_EXCHANGE_PROB : 0;
	jmode = jmode_old = 0;
	if (Nmodes>1)
	  {
		// set up probabilities of mode jumps. watch out for over/underflow etc
		p_mode_prop = new double[Nmodes];
		p_mode_CDF = new double[Nmodes];
		min = Double.MAX_VALUE;
		for (i=0;i<Nmodes;i++) min = min((optima.get(i)).logMMLapprox, min);
		// allow heavier tail on weak modes: no mode will have probability less than pmin
		pmin = 0.05/Nmodes;
		for (i=0,prob=0.0; i<Nmodes; i++) 
		  {
			prob += max(Math.exp(-((optima.get(i)).logMMLapprox - min)), pmin);
		  }
		// form PDF
		for (i=0;i<Nmodes;i++) 
		  {
			p_mode_prop[i] = max(Math.exp(-((optima.get(i)).logMMLapprox - min)),pmin)/prob;
		  }
		// turn into CDF
		p_mode_CDF[0] = p_mode_prop[0]; 
		for (i=1;i<Nmodes;i++) p_mode_CDF[i] = p_mode_CDF[i-1] + p_mode_prop[i]; 
	  }
	// Do sampling for pure unsmoothed model at present
	mu = 0.0f;
	log10mu = -10.0f;

	Algebra alg = new Algebra();

	// Start at best LocalMode solution.
	// copy mode point, inverse covariance etc.
	m0 = optima.get(0);
	pmu = m0.pmuLoc;
	System.arraycopy(m0.pmLoc,0,pm,0,Nptot);
	// set scaling
	if (!true)
	  {
		if (verbosity > 1) EH.println("Scaling set by Born-like approx. chunk");
		// set scaling by sensitivity to flat model (Born-like approximation)....
		for (i=0;i<Nptot;i++) pm[i] = 0.0f;
		makjtj();			
		// by here, wjtwj = (WJ)^T.(WJ)
		H = computeFullHessian(pmu);
		// now wjtwj = (WJ)^T WJ + partial^T.partial + W_B^T.W_B
		// but evaluated for the flat (m_i = 0) model 
	  }
	else
	  {
		if (verbosity > 1) EH.println("Scaling set by covariance of best solution");
		// set from localMode covariance of best solution so far
		H.assign(m0.H);

		// weight with other mode solutions if they seem reasonably plausible too
		nModes=1;
		for (i=1; i<optima.size(); i++)
		  {
			mi = optima.get(i);
			if (mi.negLogPost < m0.negLogPost + NEG_LOG_POSTERIOR_PLAUSIBLE_MODE_GROUP_RANGE)
			  {
				// accumulate sum of matrices in H
				H.assign(mi.H, Functions.plus);
				if (verbosity > 1) EH.println("...averaged with mode "+(i+1));
				nModes++;
			  }
		  }
		// average covariance
		H.assign(Mult.div((double)nModes));
		// Compute approx std deviations
		// a = (cern.colt.matrix.impl.DenseDoubleMatrix2D)alg.inverse(H);
		a = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(H);
		stdv = new float[Nptot];
		for (i=0;i<Nptot;i++) stdv[i] = (float)sqrt(a.getQuick(i,i));

		if (verbosity > 1) 
		  {
			EH.printf("Approx (for sampler) guessed param std. devs:%n");
			for (i=0;i<Nptot;i++) EH.printf("%-15.8g ",stdv[i]);
			EH.println("");
		  }    
	  }

	// Get Scaling covariance
	// At this point, after computeFullHessian(), 
	// H = (WJ)^T WJ + partial^T.partial + W_B^T.W_B
	a = new DenseDoubleMatrix2D(Nptot,Nptot);
	a.assign(H);
	
	// if (verbosity > ) EH.println("Cinv = "+a);
	if (!true)
	  {
		// Set scaling to be independent, based on diagonal of covariance of best solution
		// a = (cern.colt.matrix.impl.DenseDoubleMatrix2D)alg.inverse(a);
		a = ModifiedCholeskyLeftDecomposition.getInverseViaCholesky(a);
		stdv = new float[Nptot];
		for (i=0;i<Nptot;i++) stdv[i] = (float)sqrt(a.getQuick(i,i));
		if (verbosity > 1) EH.printf("Using Param std. devs: ");
		for (i=0;i<Nptot;i++) EH.printf("%-15.8g ",stdv[i]);
		EH.printf("%n");
		// EH.println("Cinv2 = "+a);
		a.assign(0.0);
		for (i=0;i<Nptot;i++) a.setQuick(i,i,stdv[i]);
		Rtri = a;
	  }
	else if (true)
	  {
		// Set scaling by Cholesky factor of local approximate covariance
		ModifiedCholeskyRightDecomposition cholDR = new ModifiedCholeskyRightDecomposition(a); 
		Rtri =  cholDR.getR();
		cholDR.getCholeskyRightInverse(Rtri);
		// Rtri is now an upper triangular Cholesky factor of C
		// viz C = {\tilde C} = Rtri^T . Rtri
		//EH.println("R="+Rtri+"\n");

		// Diagonal only????
		min = Rtri.getQuick(0,0);
		Rtri.assign(0);
		for (i=0;i<Nptot;i++) Rtri.setQuick(i,i,min);
	  }
	else
	  {
		// set scale on ith parameter to minimum error bar over supported modes
		EH.printf("FIXME - new MCMC scaling experiment!!!!!");
		Rtri = new DenseDoubleMatrix2D(Nptot,Nptot);
		for(j=0;j<Nptot;j++) Rtri.setQuick(j,j,m0.stdv[j]);
		for (i=1; i<optima.size(); i++)
		  {
			mi = optima.get(i);
			if (mi.negLogPost < m0.negLogPost + NEG_LOG_POSTERIOR_PLAUSIBLE_MODE_GROUP_RANGE)
			  {
				// find minimum error bar over supported modes
				for (j=0;j<Nptot;j++) Rtri.setQuick(j,j,min(mi.stdv[j], Rtri.getQuick(j,j)));
			  }
		  }
	  }
	if (verbosity > 0) EH.printf("FIXME GET SCALING SIMILAR HERE IN DEPTH???....R:"+Rtri);
	EH.println("");

	// Add ? small diagonal noise
	if (true)   
	  {
		float w_noise = 1.5f;
		for (i=0;i<Nptot;i++) Rtri.setQuick(i,i,Rtri.getQuick(i,i) * w_noise);
	  }

	// setups for layer-flipping. Only flip layers of reasonably comparable thickness, or thickness<1/5 depth
	flippableLayers = getSuggestedFlippableLayerList();
	
	// Setups for MD sampling
	if (doMDstyleMCMC)
	  {
		m_MD = new double[Nptot]; // masses
		p_MD = new double[Nptot]; // leapfrog momenta
		a_MD = new DenseDoubleMatrix2D(Nptot,Nptot);
		deltaT_MD = 0.1f;   // time step
		// Simplest case. Masses  ~ 1/stdv
		if (verbosity > 1) EH.printf("MD masses: ");
		for (i=0;i<Nptot;i++) 
		  {
			m_MD[i] = 1.0/stdv[i];
			EH.printf("%7.2g ",m_MD[i]);
		  }
		EH.println("");
	  }

	mopt = new float[Nptot];
	m_last_accepted = new float[Nptot];
	delta_m = new float[Nptot];

	if (verbosity > 1) EH.println("Writing MCMC samples");
	EH.println("(NOTE NEGATIVE EXCLUSION CODE....)");   

	write_MC_resp = !true;
	write_MC_profiles = !true;

	try
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }


	if (write_MC_profiles)
	  {
		try
		  {
			out_prof = new PrintWriter(new BufferedWriter(new FileWriter(filename_prof)));
		  } 
		catch (IOException ioe) 
		  {
			throw new RuntimeException("couldn't open " + filename_prof+" for output");
		  }
	  }

	// guessed scale factor
	scale = 0.3f/(float)sqrt(Nptot);  // 2.4/sqrt(Nptot) is magic according to Gelman etc
	obj = obj_old = 0f;
	BAD_TRY_PENALTY = 10000;

	PERMUTE_PROB = 0.5f;

	if (useBoundConstraints)
	  {
		for (i=0;i<np;i++)
		  {
			if (pm[i] <= m_lower[i])
			  {
				// start sampler a little away from bounds for safety
				pm[i] = m_lower[i] + 0.1f*stdv[i];
			  }
		  }
	  }

	System.arraycopy(pm,0,mopt,0,Nptot);
	System.arraycopy(pm,0,m_last_accepted,0,Nptot);

	i_acc=0;
	i_acc_sw =0;
	i_acc_sw_try =0;
	i_acc_hop =0;
	i_acc_hop_try =0;
	i_acc_loc=0;
	i_xs=0;

	// main chain loop. First 10% for equilibration, remaining chunk for sampling

	// MCMC Output file header..
	out.printf("#DEPTHS ");
	// top depths of free layers
	for (i=0; i<(nlayer+1); i++) if (rho[i] <= 0.0f) out.printf("%10.5g ",depths[i]);
	out.printf("%n");
	out.printf("#BASE_DEPTH %10.5g%n",depths[nlayer+1]);
	out.printf("#nd= %d%n#np= %d%n#Nptot= %d%n",nd,np,Nptot);
	out.printf("#error-rescaling-factor= %g%n",errorScale);
	out.printf("#expected-min-chisq_total= %g%n",chsqNLPmin);
	out.printf("# ");
	for (i=1;i<=Nptot;i++) out.printf("m_%d ",i);
	out.println("rms_misfit  chisq_data chisq_prior chisq_total");

	N_burn =  (int)min(0.1*N_MCMC,1000);
	N_loop = N_MCMC + N_burn;
	EH.println("FIX: PERMUTING scale calibration...");
	for (i_r = 0; i_r<N_loop; i_r++)
	  {
		if (i_r%TEST_WIN_SIZE == 0)
		  {
			// OLD acc_rate = i_acc/(double)max(1,i_r);
			acc_rate = i_acc_hop/(double)max(1,i_acc_hop_try);
			if (verbosity > 1) EH.println("Realisation %d, hop-acc.-rate = %-7.4g, %d-sample-acc.-rate = %-7.4g, swap-proposal-acc.-rate = %-7.4g",
										  i_r, acc_rate,TEST_WIN_SIZE,
										  (i_acc_loc/(double)TEST_WIN_SIZE), 
										  i_acc_sw/(double)max(1,i_acc_sw_try));
			i_acc_loc = 0;
			if (calibrating && (i_r>0))
			  {
				if (acc_rate < TARGET_MIN_ACC_RATE)
				  {
					scale /=2.0;  // half scale
					if (verbosity > 1) EH.println("Reducing scale...");
					i_acc = i_r = i_acc_hop = i_acc_hop_try = 0; 
				  }
				else if (acc_rate > TARGET_MAX_ACC_RATE)
				  {
					scale *= 2.0;  // double scale
					if (verbosity > 1) EH.println("Increasing scale...");
					i_acc = i_r = i_acc_hop = i_acc_hop_try = 0;
				  }
				else
				  {
					// assume we're OK. freeze scale
					calibrating = false;
					if (verbosity > 1) EH.println("Freezing scale for subsequent chain...");
				  }
			  }
		  }
		// Profiles......
		if (!quietIO && write_MC_profiles && (i_r%10==0)) 
		  {
			// write depths verses resistivities
			// FIXME loop over layers, use iff[]....
			for(i=0;i<np;i++)
			  {
				out_prof.printf("%15.7g %15.7g%n",depths[i+1],m_last_accepted[i]);
				out_prof.printf("%g15.7g %15.7g%n",depths[i+2],m_last_accepted[i]);
			  }
		  }
		// Responses
		if (!quietIO && write_MC_resp && (i_r%10 == 0)) writeMCMCresp(i_r/10);
		add_pen = false;
		if (!doMDstyleMCMC)
		  {
			// Standard pure random style proposals
			if ((Math.random() < PERMUTE_PROB) && (flippableLayers.length != 0))
			  {
				// Adjacent-layer "resistivity flip".
				// Choose random layer j in [0...np-2], from flippableLayers
				try_swap = true;
				i_acc_sw_try++;
				// OLD: all layers: j = (int) min(np-2,max(0,floor((np-1)*Math.random()))); 
				// Choose from amenable layers (thinnish ones which are almost a pair)
				j = (int)min(flippableLayers.length -1, floor(flippableLayers.length*Math.random())); 	
				j = flippableLayers[j];
				// Resample weighted j & j+1 values of adjacent layer-pair
				// Aim here is to generate sample with similar thickness * resistivity sum
				// Track Jacobean for reversibility
				rho_j   = fpow10(pm[j]);
				rho_jp1 = fpow10(pm[j+1]);

				float flipScale=0.4f, flipScale2=0.02f;
				float T1,T2,tmp2;
				T1 = thickf[j+1];     // thickness associated with pm[j] ..remember thicknesses start at 0 for sealayer
				T2 = thickf[j+2];

				// This algorithmn checked to be OK for bound problems, when sampling from prior only
				m_new_j = pm[j] + flipScale * (float)RD.nextGaussian(); // jump for layer j
				// assume illegal jump till proven otherwise (simplifies code)...
				add_pen = true; // kicks in rejection machinery 
				if (useBoundConstraints) pm[j] = m_lower[j]-1; // kicks in rejection machinery 
				if (m_new_j > m_lower[j])
				  {
					tmp = (rho_j * T1+  rho_jp1 * T2 - T1*fpow10(m_new_j))/T2;
					if (tmp > 0)
					  {
						tmp2 = (float)RD.nextGaussian();
						m_new_jp1 = flog10(tmp) + tmp2 * flipScale2;
						if (m_new_jp1 > m_lower[j+1])
						  {
							Rtot_new = (fpow10(m_new_j) * T1  + fpow10(m_new_jp1) * T2 - T1 * rho_j)/T2;
							if (Rtot_new > 0)
							  {
								tmp = flog10(Rtot_new);
								JFAC = 0.5f * (sqr(tmp2) -  sqr((tmp - pm[j+1])/flipScale2));
								//     assign new samples 
								pm[j] = m_new_j;
								pm[j+1] = m_new_jp1;  
								add_pen = false;
							  }
						  }
					  }
				  }
			  }  // Layer flip
			else
			  {
				if (!calibrating && (Math.random() < modeExchangeProb))
				  {
					// (Don't do mode flips when calibrating scale on random jumps)
					// a mode flip: choose another mode at random
					prob = Math.random();
					for (jmode=0; jmode<Nmodes; jmode++) if (prob <= p_mode_CDF[jmode]) break;

					float[] mbar_new = (optima.get(jmode)).pmLoc;
					float[] mbar_old = (optima.get(jmode_old)).pmLoc;
					for (i = 0;i<Nptot;i++) pm[i] += (mbar_new[i]-mbar_old[i]);
					JFAC = (jmode != jmode_old) ? (float)Math.log(p_mode_prop[jmode_old]/p_mode_prop[jmode]): 0f;
				  }
				else
				  {
					// do Metropolis Gaussian jump
					try_swap=false;
					jmode = jmode_old;
					JFAC = 0.0f;
					i_acc_hop_try++;
					for (i = 0;i<Nptot;i++)
					  {
						// uncorrelated deviates
						delta_m[i] = scale * (float)RD.nextGaussian();
					  }
					// blas: delta_m <- Rtri^T.delta_m
					trimult(Rtri,delta_m);
					// pm <- pm + delta_m
					for (i=0;i<Nptot;i++) pm[i] += delta_m[i];
				  }
			  }
		  }
		else
		  {
			// Hybrid MD method
			float Hp,Hp_old=0f;   // twice Kinetic energy for before/after MD trajectory
			for (i=0;i<Nptot;i++)
			  {
				// Random Gaussian momenta for initial momenta p(t-dt/2)
				p_MD[i] = RD.nextGaussian();
				Hp_old +=  p_MD[i]*p_MD[i];
			  }
			// Leapfrog forward stepping. Current state is pm[]
			deltaT_MD = .001f;
			N_MD = 50;
			if (verbosity > 2) EH.printf("MD:%n");
			for (j=0;j<N_MD;j++)
			  {
				// Get gradient & full objective
				rms = makjtj();
				ms = rms*rms;  
				makptp(log10mu);
				computeGradandHessian(log10mu,a_MD);
				// Full (potential) objective function
				obj = nd * ms;             // data misfit
				if (doNoiseEstimation)
				  {
					obj = (float)(nd * (1.0+Math.log(obj/nd))); // incorporates optimization on noise-scalar
				  }
				obj += mu * chisqRoughness(pm);         // prior roughness component
				obj += chisqPrej(log10mu,pm);           // prior prejudice

			  
				if (j==0) 
				  {
					obj_old = obj + Hp_old;  // initial 2xHamiltonian
				  }
				// leapfrog equations - see Jun Liu, Ch. 9
				if (verbosity > 2) EH.printf("MD: ");
				for (i=0;i<Nptot;i++)
				  {
					p_MD[i] +=  deltaT_MD * grad[i];        // momentum half-time update; p(t+dt/2)
					pm[i] += deltaT_MD * p_MD[i]/m_MD[i];      // parameter update
					EH.printf("%15.8g ",pm[i]);
				  }
				EH.printf("%n");
			  }
			// final twice kinetic energy
			for (i=0, Hp=0f; i<Nptot; i++) Hp +=  p_MD[i]*p_MD[i];
			// Final 2xHamiltonian
			obj += Hp;
			// No Jacobian piece
			JFAC = 0.0f;
		  }



		// exclude negative zones.....
		proposalLegal = true;
		if (useBoundConstraints)
		  {
			for (i=0;i<np;i++)
			  {
				if ((pm[i] < m_lower[i]) || (pm[i] > m_upper[i]))
				  {
					proposalLegal = false;
					break;
				  }
			  }
		  }
		// EH.println("It "+i_r+" legal = "+proposalLegal);
		// Evaluate proposal...
		if (!doMDstyleMCMC)
		  {
			if (proposalLegal)
			  {
				obj = twiceNegLogPosterior(mu);
				ms = ms_MCMC;
			  }
			else 
			  {
				// Illegal proposal. Use approximate cheap call to avoid fwd  evaluation...
				obj = BAD_TRY_PENALTY;
				if (i_r==0) 
				  {
					System.arraycopy(m_last_accepted,0,pm,0,Nptot); // recover from bad first proposal
				  }
			  }
		  }
		else
		  {
			// MD method. obj already evaluated...
			// Check Illegal proposal. Use approximate cheap call to avoid fwd  evaluation...
			if (!proposalLegal) 
			  {
				obj = BAD_TRY_PENALTY;
				if (verbosity > 2) EH.println("MD Illegal");
			  }
		  }
			
		// Kill ludicrous proposals (from try_swap branch...)
		if (add_pen) obj += BAD_TRY_PENALTY;
	

		//     Metropolis test
		ex = -JFAC + (obj - obj_old)/2.0f;

		if (proposalLegal && ((i_r==0) || (ex < 0) || (Math.random() <= Math.exp(-ex))))
		  {
			// Accept & save
			System.arraycopy(pm,0,m_last_accepted,0,Nptot);
			obj_old = obj;
			ms_old = ms;
			msprej_old = msprej_MCMC;
			i_acc++;
			i_acc_loc++;
			if (try_swap) 
			  {
				i_acc_sw++;
			  }
			else
			  {
				i_acc_hop++;
			  }
			if (jmode != jmode_old) 
			  {
				jmode_old = jmode;
				EH.println("Mode-flip accepted");
			  }
		  }
		else
		  {
			// Reject. Restore old state
			System.arraycopy(m_last_accepted,0,pm,0,Nptot);
		  }

		// MCMC Sampling output
		if (i_r >= N_burn)
		  {
			if (!decimate || (decimate && (i_r % decimation_rate == 0)))
			  { 
				// format: Nptot columns of model parameters
				for (i=0;i<Nptot;i++) out.printf("%15.8g ",pm[i]);
				
				// RMS (wrt new noise level if EB)
				out.printf("%15.8g",sqrt(ms_old));
				// Chisq data misfit
				out.printf("%15.8g",ms_old*nd);
				// Chisq prior
				out.printf("%15.8g",msprej_old);
				// then overall 2*NLP  = chisq_total
				out.printf("%15.8g%n",obj_old);
			  }
		  }
	  } // realisations loop
	
	// Clean up files etc
	out.close();

	if (verbosity > 1) EH.println("Wrote "+filename);

	if (write_MC_profiles) 
	  {
		out_prof.close();
	  }

	if (adjustNoiseLevelToEBestimate)
	  {
		// reset noise to defaults etc...
		doNoiseEstimation_old = doNoiseEstimation_old;
		global_error_scale = -1.0f;
	  }

	if (verbosity > 1) EH.println("Finished MCMC sampler, acceptance rate = %8.4g%n",i_acc/(double)(N_loop));
  }


  int[] getSuggestedFlippableLayerList()
  {
	int i,j;
	float tmp,y,z;
	float FLIPPABLE_THICKDEPTH_THRESHOLD= 0.25f;

	IntArrayList flL = new cern.colt.list.IntArrayList();
	for (i = 0, j=0; i<nlayer; i++)  // scan only to second-last layer
	  {
		if (rho[i] <= 0.0f)
		  {
			if (rho[i+1] <= 0.0f)
			  {
				tmp = depths[i]-thickf[0]; // depth below seawater
				y = (float) ((tmp>0) ? (thickf[i]/tmp) : Float.MAX_VALUE);
				tmp = depths[i+1]-thickf[0]; // depth below seawater
				z = (float) ((tmp>0) ? (thickf[i+1]/tmp) : Float.MAX_VALUE);
				if (max(y,z) < FLIPPABLE_THICKDEPTH_THRESHOLD) flL.add(j);  // FLIPPING THRESHOLD
			  }
			j++;
		  }
	  }
	flL.trimToSize();
	int[] flippableLayers = flL.elements();
	if (verbosity>1)
	  {
		EH.printf("Flippable layers indexes: ");
		for (i=0;i<flippableLayers.length;i++) EH.printf("%d ",flippableLayers[i]); 
		EH.println("");
	  }
	return flippableLayers;
  }
  /**
	 In-place Triangular multiplication. R is upper triangular. Forms:
	 delta_m <- Rtri^T.delta_m
  */
  void trimult(DoubleMatrix2D R,float[] delta_m)
  {
	// blas: delta_m <- R.delta_m
	int i,j;
	float s;

	// float[] tmp = delta_m.clone();

	for (i = delta_m.length-1 ; i>=0 ;i--)
	  {
		s = 0.0f;
		for (j=0; j<=i; j++) s += R.getQuick(j,i) * delta_m[j];
		delta_m[i] = s;  // is not used again...
		// tmp[i] = s;
	  }
	// for (i = 0; i<delta_m.length ;i++) delta_m[i]=tmp[i];
	//	tmp=null;
  }


  /**
	 Writes an MCMC resp file
  */
  void writeMCMCresp(int n)
  {
	int i,j;
	String filename;
	PrintWriter out;

	// Make suitable file
	filename = String.format("MCMCresp_%06d",n);
	try
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }
	for(i=0;i<nd;i++)
	  {
		out.printf("%5d ",(int)dp.getQuick(0,i));
		for(j=1;j<ndp;j++) out.printf("%15.5g ",dp.getQuick(j,i));
		out.printf("%15.5g %15.5g %15.5g",d[i], dm[i], (d[i]-dm[i])/sd[i]);
	  }
	out.close();

	if (verbosity > 1) EH.println("Wrote: "+filename);
  }

  /**
	 Returns total -ve log posterior at current model parameters pm[], parameter mu
  */

  float twiceNegLogPosterior(float mu)
  {
	// locals
	float x, y;
	int i;
	boolean samplePriorOnly = false;


	long t = System.nanoTime();
	if (!samplePriorOnly) 
	  {
		doForwardModel(pm);
		// ScrippsLegacyModelInterface.formod(np,nd,pm,dp.elements,dm,d);
	  }
	tNatives += (System.nanoTime()-t); 

	y = 0.0f;
	if (!samplePriorOnly) 
	  {
		if (doNoiseEstimation)
		  {
			// Need replacement C_D -> chisq^2_D/n_p C_D
			y += chisqMisFit();
			ms_MCMC = y/nd;
			y = (float)(nd * (1.0+Math.log(y/nd))); // incorporates optimization on noise-scalar
		  }
		else if (global_error_scale > 0f)
		  {
			y += (chisqMisFit()/sqr(global_error_scale));
			ms_MCMC = y/nd;
		  }
		else
		  {
			y += chisqMisFit();
			ms_MCMC = y/nd;
		  }
	  }
	y += mu*chisqRoughness(pm);
	msprej_MCMC = chisqPrej( (mu>0) ? flog10(mu) : -10f, pm);
	y += msprej_MCMC;

	return y;
  }


  void inputData()
  {
	if (newScripps)
	  {
		SMI.readData();
		nd = SMI.OM.nd;
		SMI.finishModelSetup();
		SMI.reconcileDataVariables();
		allocateData();
	  }
	else
	  {
		inputd();
		// Debug fake data set & fwd model
		if (ScrippsLegacyModelInterface.TEST_OPTIMISER) 
		  {
			ScrippsLegacyModelInterface.makeTestOptimiserCalibrationDataSet(this);
		  }
	  }
	getMinimumFrequency();
  }
  
  // Compute lowest freq of data
  void getMinimumFrequency()
  {
	if (newScripps)
	  {
		int i;
		for (i=1, freqMin=(float)(SMI.CM.fTxIn[0]);i<SMI.CM.fTxIn.length; i++) 
		  {
			freqMin = (float)min(freqMin,SMI.CM.fTxIn[i]);
		  }
	  }
	else
	  {
		freqMin = dp.getQuick(2,0);
	  }
  }


  /**
	 Read old-style scripps input data files
  */
  void inputd()
  {
	String dform, ddescr, string;
	int lerr, i, j, k, type=0,l_count,iform;
	boolean ifcsem, ifmt;
	BufferedReader in=null;
	SimpleEMData data;
	String s,colheaders;
	String[] st;

	ifcsem = false;
	ifmt = false;
	  
	if (verbosity > 1) EH.println(" reading data file..."+datfil);
	// open data file 
	try 
	  {
		in = new BufferedReader(new FileReader(datfil));
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't open file "+datfil);
	  }
	try
	  {
		// c read character descriptors      
		dform = stringPastCol18(in.readLine());
		if (!dform.equals("OCCAM1DDATA_1.0"))
		  {
			EH.die_err("wrong data file type: not OCCAM1DDATA_1.0");
		  }
		// description
		ddescr = stringPastCol18(in.readLine());
		// read in the number of data
		nd = intPastCol18(in.readLine());
		// headers
		colheaders = stringPastCol18(in.readLine());
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't parse file "+datfil);
	  }
	// read & grow array...
	ArrayList<SimpleEMData> data_list = new ArrayList<SimpleEMData>(nd);

	l_count=0;
	for (i=0;i<nd;i++)
	  {
        try
		  {
			if ((s = in.readLine()) != null)
			  {
				st = s.trim().split("\\s+");
				try
				  {
					type = Integer.parseInt(st[0]);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Line containing "+s+" has unintelligible type-integer in "+datfil);
				  }
				if (type < 20) 
				  {
					ndp = 6;
					ifcsem = true;
				  } 
				else 
				  {
					ndp = 2;
					ifmt = true;
				  }
				data = new SimpleEMData(ndp);

				data.dp[0] = (float)type;
				if (st.length == (ndp+2))
				  {
					try
					  {					
						data.d = Float.parseFloat(st[1]);
						data.sd = Float.parseFloat(st[2]);
						if (data.sd<=0)
						  {
							EH.die_err("Line containing "+s+" has zero error in file "+datfil+" Can't invert!");
						  }
						for (j=1,k=3; j<ndp; j++,k++)
						  {
							data.dp[j] = Float.parseFloat(st[k]);
						  }
						//dp(i,1), d(i), sd(i), (dp(i,j),j=2,ndp)
						data_list.add(data);
						l_count++;
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Line containing "+s+" has unintelligible number in "+datfil);
					  }
				  }
				else
				  {
					EH.die_err("Incorrect number of columns in line: "+s+" of "+datfil);
				  }
			  }
		  }
		catch (Exception e)
		  {
			EH.die_err("Can't read line in file "+datfil);
		  }
	  }
	if (l_count < nd)
	  {
		EH.die_err("Number of data read doesn't match specified number in "+datfil);
	  }


	try
	  {	
		in.close();
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't close file "+datfil);
	  }

	if (verbosity > 1) EH.println("Read "+nd+" data points");

	// sort data 
	DataComparator cmp = new DataComparator();
	// sort data by freq, range, type...
	Collections.sort(data_list,cmp);

	// Now move into standard fortran-friendly 
	// structures (ensure contiguousness here)
	// Allocations!
    allocateData();

	for (i=0;i<data_list.size();i++)
	  {
		data = data_list.get(i);
		d[i] = data.d;
		sd[i] = data.sd;

		// dp[2] = freq, 
		// dp[1] = range,
		// dp[0] = type.
		// dp[5] = trans. height.

		for (j=0;j<6;j++) dp.setQuick(j,i,data.dp[j]);
	  }
	// data_list = null; // GC	
  }

  void lumpTransmitterHeights()
  {
	int i,j,imin,imax,type;
	float zxmitt,freq1,freq,txAv,hmin;

	float dZBucketSize = (float)LTH_BUCKET_SIZE;

	if (newScripps)
	  {
		EH.warn("FIXME check lumpTransmitterHeights() for new Scripps kernel (already done internally I think)");
		return;
	  }

	origTransHeights = new float[nd];

	// copy heights;
	for (i=0;i<nd;i++) origTransHeights[i] = dp.getQuick(5,i);

	if (!approximateTxHeights) return;

	// seek groups. By here, data sorted by freq, range, Tx height

	for (i=0;i<nd;i++)
	  {
		imin = i;
		imax = nd-1;
		freq1 = dp.getQuick(2,i);
		hmin = dp.getQuick(5,i);        // transmitter height - dp[5,j]
		for (j = imin+1; j< nd; j++)
		  { 
			// transmitter height - dp[5,j]
			// transmitter freq - dp[2,j]
			// data type dp[0,j]
			zxmitt = dp.getQuick(5,j);	
			type = (int)dp.getQuick(0,j);
			freq = dp.getQuick(2,j); 
			if ((abs(freq - freq1) > 1e-6) || (type >= 20) || (abs(zxmitt -hmin) > dZBucketSize) ) 
			  {
				// different freq or not CSEM data
				imax = j-1;
				break;
			  }
		  }
		// Groupable data spans [imin..imax]
		// EH.printf("THL group %d %d%n",imin,imax);
		// Average it
		for (j = imin,txAv = 0.0f; j<= imax ; j++) txAv += origTransHeights[j];
		txAv /= (imax-imin+1);
		// And overwrite
		for (j = imin; j<= imax ; j++) 
		  {
			dp.setQuick(5,j,txAv);
		  }
		i=imax; // makes next increment good  
	  }

  }


  void checkNewData()
  {
	LocalMode m0;
	String old_datfil;
	int i;
	double x,y,rms;

	EH.println("Checking MAP model on new data in "+checkNewDataFile);

	// Fix according to best LocalMode solution.
	m0 = optima.get(0);
	pmu = m0.pmuLoc;
	System.arraycopy(m0.pmLoc,0,pm,0,np);

	// read in new data to test with in-situ model
	old_datfil = new String(datfil);
	datfil = checkNewDataFile;
	if (newScripps)
	  { 
		SMI.readData();
	  }
	else
	  {
		inputd();
	  }
	lumpTransmitterHeights();
	doForwardModel(pm);

	y = chisqMisFit();

	rms = sqrt(y/nd);
	showParams(pm,false); 

	float LOGMU=-10.0f;
	float chisqPrej = chisqPrej(LOGMU,pm);

	EH.printf("New_data_file_check(%s): ",checkNewDataFile);
	for(i=0;i<np;i++) EH.printf("%-12.7g ",pm[i]);
	EH.printf("%-12.7g %-12.7g %-12.7g%n",rms,chisqPrej,(y+chisqPrej)/(nd+np));

	datfil = old_datfil; 
  }


  void doTHApproxTest()
  {
	int i,j,k;
	float norm,err,errN,mu=0f;
	PrintWriter out;
	LocalMode m0;

	if (newScripps)
	  {
		EH.println("Transmitter height approximations test not yet implemented for new Scripps kernel");
		return;
	  }

	// Fix according to best LocalMode solution.
	m0 = optima.get(0);
	pmu = m0.pmuLoc;
	System.arraycopy(m0.pmLoc,0,pm,0,np);

	mu = fpow10(pmu);
	if (verbosity > 1) EH.println("Performing transmitter height approximations test:");
	
	String filename = String.format("Trans_height_approximation_errors.txt");
	try
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }
	// evaluate in-situ model
	doForwardModel(pm);
	// ScrippsLegacyModelInterface.formod(np,nd,pm,dp.elements,dm,d);
	float[] dm_approx = (float[])dm.clone();
	float[] h_approx = new float[nd];
	System.arraycopy(dm,0,dm_approx,0,dm.length);

	// restore old depths, and get fwd model
	for (j = 0;j<nd; j++)
	  {
		h_approx[j] = dp.getQuick(5,j);
		dp.setQuick(5,j,origTransHeights[j]);
	  }
	doForwardModel(pm);
	// ScrippsLegacyModelInterface.formod(np,nd,pm,dp.elements,dm,d);

	// scan errors
	float errMax=0.0f,errMaxStd=0.0f;
	out.printf("#Type range freq. azimuth rot theight Approx-t-height True-model Approximate-model Relative-error Noise-relative-error%n");
	for (j = 0;j<nd; j++) 
	  {
		norm = (dm[j] != 0.0f) ? abs(dm[j]) : 1.0f;
		err = (dm[j] - dm_approx[j])/norm;
		errN = (dm[j] - dm_approx[j])/sd[j];
		out.printf("%4d ",(int)dp.getQuick(0,j));
		for (k=1;k<6;k++) 
		  {
			out.printf("%11.4g ",dp.getQuick(k,j));
		  }
		out.printf("%11.4g ",h_approx[j]);
		out.printf("%15.5g %15.5g %15.5g %15.5g%n",dm[j],dm_approx[j],err,errN);
		errMax = max(errMax,abs(err));
		errMaxStd = max(errMaxStd,abs(errN));
	  }
	out.close();
	if (verbosity > 1) 
	  {
		EH.println("Wrote "+filename);
		EH.printf("Maximum approximation-error= %7.3g, Noise-relative-error=%7.3g%n",errMax,errMaxStd);
	  }

	// restore approximations
	for (j = 0;j<nd; j++) dp.setQuick(5,j,h_approx[j]);
	h_approx=null;
	dm_approx=null;
  }

  /**
	 Dump residuals of optimization to file, for all optimization runs
  */
  private void dumpConvergenceDiagnostics()
  {
	int i,j,n;
	PrintWriter out;
	cern.colt.list.FloatArrayList res;
	String filename = String.format("ConvergenceDiagnostics.txt");

	try
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename+" for output");
	  }	
	for (i=0; i< convergenceResiduals.size(); i++)
	  {
		res = convergenceResiduals.get(i);  
		float[] y = res.elements();
		n = res.size();
		float ylast = y[n-1];
		for (j=0; j<(n-1); j++) out.printf("%15.8g%n",y[j] - ylast);
		out.printf("%n");
	  }
	out.close();
	if (verbosity > 1) EH.println("Wrote "+filename);
  }

  private String getNextNonCommentLine(BufferedReader in,  String filename)
  {
	String s = null;
	try 
	  {
		while ((s = in.readLine()) != null)
		  {
			if (!(s.startsWith("#"))) break;
		  }
	  }
	catch (IOException e)
	  {
		EH.die_err("Can't read from file "+filename);
	  }
	return s;
  }

 private String stringPastCol18(String s)
  {
	if (s!=null && s.length()>18) 
	  {
		return s.substring(18).trim(); 
	  }
	else
	  {
		EH.die_err("startup file ended prematurely");
	  }
	return null;
  }
  
  private int intPastCol18(String s)
  {
	try
	  {
		return Integer.parseInt(stringPastCol18(s));
	  }
	catch (NumberFormatException e)
	  {
		EH.die_err("Line containing "+s+" not an integer");
	  }
	return 0;
  } 
 private float floatPastCol18(String s)
  {
	try
	  {
		return Float.parseFloat(stringPastCol18(s));
	  }
	catch (NumberFormatException e)
	  {
		EH.die_err("Line containing "+s+" not a floating point number");
	  }
	return 0f;
  }

  private int parseInt(String s)
  {
	try
	  {
		return Integer.parseInt(s);
	  }
	catch (NumberFormatException e)
	  {
		EH.die_err("Line containing "+s+" not an integer");
	  }
	return 0;
  } 

 private float parseFloat(String s)
  {
	try
	  {
		return Float.parseFloat(s);
	  }
	catch (NumberFormatException e)
	  {
		EH.die_err("Line containing "+s+" not a floating point number");
	  }
	return 0f;
  }


  private boolean copyFile(String fromFile, String toFile)
  {
    try
	  {
		if (toFile.equals(fromFile)) return true;

		File f1 = new File(fromFile);
		File f2 = new File(toFile);
		InputStream in = new FileInputStream(f1);
		OutputStream out = new FileOutputStream(f2);

		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0)
		  {
			out.write(buf, 0, len);
		  }
		in.close();
		out.close();
	  }
    catch(FileNotFoundException ex)
	  {
		EH.warn(ex.getMessage() + " in the specified directory.");
		return false;
	  }
    catch(IOException e)
	  {
		EH.warn(e.getMessage());      
		return false;
	  }
	return true;
  }

  private void debugPrintCpInv(double[] a,double[] a_d1,double[] a_d2)
  {
	int i,j;
	PrintWriter out=null;
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter("Cp_inv.txt")));
	  } 
	catch (IOException ioe) {}
	for (i=0;i<np;i++)
	  {
		for(j=0;j<np;j++)
		  {
			if (j==i) out.printf("%16.9g ",a[i]);
			else if ((j==i+1) && (a_d1!=null)) out.printf("%16.9g ",a_d1[i]);
			else if ((j==i-1) && (a_d1!=null)) out.printf("%16.9g ",a_d1[i-1]);
			else if ((j==i+2) && (a_d2!=null)) out.printf("%16.9g ",a_d2[i]);
			else if ((j==i-2) && (a_d2!=null)) out.printf("%16.9g ",a_d2[i-2]);
			else  out.printf("%8.3g ",0.0);
		  }
		out.println(""); 
	  }
	out.close();
  }

  /**
	 Assemble problematic files into a zip file
  */
  void dumpDebugPackage(String[] argv)
  {
	int i;
	String s,filename,commandline;
	BufferedReader in  = null;
	DebugPackageDumper dumper;
	String[] txtFiles=null;

	// Save commandline
	StringBuffer cl = new StringBuffer();
	cl.append("java "+classname+" ");
	for (i=0; i<argv.length; i++) cl.append(argv[i]+" ");
	commandline = cl.toString();

	if (!newScripps)
	  {
		// read startup file to get files....
		filename = (altStartFile != null) ? altStartFile : "startup";
		try 
		  {
			in = new BufferedReader(new FileReader(filename));			
			for (i=0;i<2;i++) getNextNonCommentLine(in, filename); // read format, description
			s = getNextNonCommentLine(in, filename);
			modfil = stringPastCol18(s);
			s = getNextNonCommentLine(in, filename);
			datfil = stringPastCol18(s);
			in.close();
		  }
		catch (FileNotFoundException e)
		  {
			EH.die_err("Can't find file "+filename);
		  }
		catch (IOException e)
		  {
			EH.die_err("Can't read from file "+filename);
		  }
		txtFiles = new String[] {filename,modfil,datfil};
	  }
	else
	  {
		EH.println("FIXME!: dump-debug not implemented yet for new Scripps kernel...");
		System.exit(0);
	  }
	dumper = new DebugPackageDumper(EH, txtFiles, null, null, commandline, null);
	dumper.dumpPackage();
  }


  /**
	 Info class for making synthetics with prescribed noise
  */
  public static class NoiseSpec
  {
	float synthetic_error_E = 0f;
	float synthetic_error_noise_floor_E = 1e-15f;
	float synthetic_error_B = 0f;
	float synthetic_error_noise_floor_B = 1e-18f;
	float phase_error_E = 5;
	float phase_error_B = 5;

	/*
	  Set up noise levels & floors from class arguments: inherited from commandline
	*/
	public NoiseSpec(ErrorHandler EH, String syntheticSpecFrac, String syntheticSpecFloor)
	{
	  int N;
	  String s[];

	  // Fractions & Angles
	  s = syntheticSpecFrac.split(",");
	  N = s.length;
	  if (N==1)
		{
		  try 
			{
			  synthetic_error_E = Float.parseFloat(s[0]);
			  synthetic_error_B = synthetic_error_E;
			  // default rough guess from uncorrelated errors in Re(E), Im(E)
			  // Rough formula assuming frac_error relative error
			  // in Er,Ei, phase_deg = (180/pi) atan(Ei/Er)
			  phase_error_E = phase_error_B = 0.01f * synthetic_error_E * SQRT2_180_ON_PI;
		  }
		  catch  (NumberFormatException e) 
			{
			  EH.die_err("Not a number in "+syntheticSpecFrac);
			}
		}
	  else if (N==4)
		{		  
		  try 
			{
			  synthetic_error_E = Float.parseFloat(s[0]);
			  phase_error_E = Float.parseFloat(s[1]);
			  synthetic_error_B = Float.parseFloat(s[2]);
			  phase_error_B = Float.parseFloat(s[3]);
			}
		  catch  (NumberFormatException e) 
			{
			  EH.die_err("Not a comma-separated list of 4 numbers in "+syntheticSpecFrac);
			}
		}
	  else
		{
		  EH.die_err("Noise fractional specification "+syntheticSpecFrac+
					 " must be of form: percent_E[,degrees_Ephase,percent_B,degrees_Bphase]");
		}


	  // Floors
	  s = syntheticSpecFloor.split(",");
	  N = s.length;
	  if (N==1)
		{
		  try 
			{
			  synthetic_error_noise_floor_E = Float.parseFloat(s[0]);
			  synthetic_error_noise_floor_B = 0.001f * synthetic_error_noise_floor_E; // Default guess
			}
		  catch  (NumberFormatException e) 
			{
			  EH.die_err("Not a number in "+syntheticSpecFloor);
			}
		}
	  else if (N==2)
		{		  
		  try 
			{
			  synthetic_error_noise_floor_E = Float.parseFloat(s[0]);
			  synthetic_error_noise_floor_B = Float.parseFloat(s[1]);
			}
		  catch  (NumberFormatException e) 
			{
			  EH.die_err("Not a comma-separated list of 2 numbers in "+syntheticSpecFloor);
			}
		}
	  else
		{
		  EH.die_err("Noise floor specification "+syntheticSpecFloor+
				 " must be of form    Efloor[,Bfloor]");
		}		

	  // Checks
	  if ((synthetic_error_E < 0) || (phase_error_E < 0) || (synthetic_error_B < 0) || (phase_error_B < 0))
		{
		  EH.die_err("Need numbers >= 0 in "+syntheticSpecFrac);
		}
	  if (synthetic_error_noise_floor_E < 0) EH.die_err("Need E_floor number>=0 in "+syntheticSpecFloor);
	  if (synthetic_error_noise_floor_B < 0) EH.die_err("Need B_floor number>=0 in "+syntheticSpecFloor);
	}
  }


  /**
	 Prior trend info for resistivity
  */
  public static class PriorTrendInfo
  {
	CSEM1DInverter _CI;

	/**
	   Linear regression for mean log_10(rho):   <log_10(rho)> = A + B * zTVDBML
	*/
	float A_mean = 0f;         // 1  Ohm.m at seafloor
	float B_mean = 0.25e-3f;   // 10 Ohm.m. at 4km

	/**
	   Linear regression for mean log_10(rho):   std_dev ( log_10(rho) )= A + B * zTVDBML
	*/
	float A_sigma = 1f;         // +/- 1  Ohm.m at seafloor
	float B_sigma = 0.25e-3f;   // +- 2 Ohm.m. at 4km

	public PriorTrendInfo(CSEM1DInverter CI)
	{
	  _CI = CI;
	}

	float getDefaultPriorMean(float zBML)
	{
	  return A_mean + B_mean * zBML;
	}

	float getDefaultPriorSigma(float zBML)
	{
	  return A_sigma + B_sigma * zBML;
	}

	String getConsts()
	{
	  return A_mean+","+B_mean+","+A_sigma+","+B_sigma;
	}

	void overrideDefaults(String s)
	  {
		String spec[];
		spec = s.split(",");
		if (spec.length != 4) 
		  {
			_CI.EH.die_err("prior-depth-trend spec "+s+" needs to be mean_intercept,mean_slope(1/m),stddev_interecpt,stddev_slope(1/m)");
		  }
		else
		  {
			try {A_mean = Float.parseFloat(spec[0]);}
			catch  (NumberFormatException e) { _CI.EH.die_err(spec[0] +" not a number in "+s);}
			try {B_mean = Float.parseFloat(spec[1]);}
			catch  (NumberFormatException e) { _CI.EH.die_err(spec[1] +" not a number in "+s);}
			try {A_sigma = Float.parseFloat(spec[2]);}
			catch  (NumberFormatException e) { _CI.EH.die_err(spec[2] +" not a number in "+s);}
			try {B_sigma = Float.parseFloat(spec[3]);}
			catch  (NumberFormatException e) { _CI.EH.die_err(spec[3] +" not a number in "+s);}
		  }
	  }
  }


  /** 
	  Info for Marquardt strategies
  */
  public static class LevMarInfo
  {
	boolean terminateLM;
	boolean restartLM;
	float lambda_LM,nu;
	float[] scales;
	float[] grad_modified = null;

	DenseDoubleMatrix2D H_modified = null;
	byte[] _act_PN=null;                    // projected Newton working actives (reference)

	public LevMarInfo(int n)
	{
	  scales = new float[n];
	  Arrays.fill(scales,1.0f); // Default choice
	}
	/** 
		set scales[] by diag of supplied matrix (A = np x np), scaled s.t.
		largest entry is 1, i.e. scales[i]=A_ii/(max{A_jj})
	*/
	void setLMscalesByDiag(DenseDoubleMatrix2D A)
	{
	  int i,Nptot;

	  Nptot = scales.length;

	  if (false)
		{
		  float max=0f;
		  for (i=0;i<Nptot;i++) max = max(abs((float)A.getQuick(i,i)),max);
		  if (max==0f) max=1.0f;
		  for (i=0;i<Nptot;i++) scales[i] = abs((float)A.getQuick(i,i))/max;
		}
	  else
		{
		  Arrays.fill(scales,1.0f); 
		}
	}
  }

  /**
	 Local covariance and std dev storage class
  */
  public static class ModeCovarianceAndStdDev
  {
	float[] stdv;
	DoubleMatrix2D C; 
  }

  /**
	 Class for finding optimal link connecting modes, using Uncmin's BFGS etc
   */
  public static class ModeLinkFinder implements Uncmin_methods  
  {
	CSEM1DInverter _CI;

	/**
	   Number of internal nodes along path integral
	*/
	public static int N_PTS = 30; //30
	public static boolean TRY_TO_EQUISPACE = true;
	public static boolean EQUISPACE_TEST = false;
	public static float EQUISPACE_A_DEFAULT = 1.0f;
	float EQUISPACE_A = EQUISPACE_A_DEFAULT;

	// useful shared info across methods
	int np=0; // copy of _CI.np (dim of CSEM problem)
	int N=0;  // number of free intermediate points in np-space along trajectory
	int n=0;  // N*np
	int nIts=0;
	float[][] m;
	float[]   f,dm2,gradL;
	float[][] gradf;
	float msbar;
	String _label=null;
	boolean savePaths=true;
	ArrayList<float[][]>path = null;
	ArrayList<float[]>path_f=null;

	// Control variables for UncMin optimizer.
	// see dfault_f77 for example of how this is done in Uncmin....

	double[] f_best = new double[2];             // f[1] has the function minimum at exit
	int[] info = new int[2];                     // status information about the code
	double[] fscale = new double[2];
	int[] method = new int[2];
	int[] is_exp = new int[2];
	int[] msg = new int[2];
	int[] n_good_digits = new int[2];
	int[] max_its = new int[2];
	int[] is_gradient_supplied = new int[2];
	int[] is_hessian_supplied = new int[2];
	double[] trust_radius = new double[2];
	double[] gradient_tol = new double[2];
	double[] max_step_size = new double[2];
	double[] step_tol = new double[2];

	public ModeLinkFinder(CSEM1DInverter CI)
	{
	  _CI = CI;

	  // Number of points collocating the minimum path
	  N = N_PTS;  // 30

	  // dimension of optimisation problem
	  np = _CI.np;
	  n = np*N;
	  // suite of points
	  m = new float[N+2][np];
	  dm2 = new float[N+1];
	  gradf = new float[N+2][np];
	  f = new float[N+2];
	  gradL = new float[n+1];
	}

	/**
	   Two mode estimates are passed in in mA[0...np-1] and mB[0...np-1]
	   We attempt to find by collocation and discretizing a path-integral, a connecting pathway
	   between the two modes using BFGS with line-search optimization.
	*/
	private void findModePairConnection(float[] mA,float[] mB,String label)
	{
	  int i,j,k,l,iM=0,d=2;
	  float s,frac,L,Lo,lambda;
	  
	  PrintWriter out=null;
	  String filename=null;

	  _label=label;

	  // linearly interpolate starting state
	  for (i=0;i<=(N+1);i++)
		{
		  frac = i/(float)(N+1);
		  for(j=0;j<np;j++)
			{
			  // linear in log(rho) space
			  // m[i][j] = frac*mB[j] + (1.0f-frac)*mA[j];    // linear in log(rho) space
			  //linear in rho space
			  m[i][j] = flog10(frac * pow(10.0,mB[j]) +  (1.0-frac)*pow(10.0,mA[j]));
			  if (EQUISPACE_TEST) 
				{
				  // "tricky" starting configuration to test spacing part of objective
				  if (frac<0.5) frac/=2;
				  m[i][j] = frac*mB[j] + (1.0f-frac)*mA[j];
				}		
			}
		}
	  if (TRY_TO_EQUISPACE)
		{
		  getSegmentLengths();
		  EQUISPACE_A /= sqr(msbar); // set scaling		 
		  if (_CI.verbosity>2)
			{
			  _CI.EH.println("Initial node spacings: ");
			  for (i=0; i<=N; i++) _CI.EH.printf("%g%n",dm2[i]);  
			}
		}

	  // BFGS search

	  // optimiser setup
	  // Tolerances...
	  double epsilon  = 1e-8; // default 1.12e-16;
	  gradient_tol[1] = pow(epsilon,1.0/3.0);
	  max_step_size[1] = 1.0;     // 0.0 default
	  step_tol[1] = sqrt(epsilon);
	  trust_radius[1]=-1;   // -1 is default

	  // Flags for methods/outputs...
	  method[1]=1; // 1 = line-search
	  is_exp[1]=1; // expensive
	  msg[1] = 10;  // sum of: 2 = don't check supplied gradient by FDs; 8 = don't print input state, final state etc
		  
	  n_good_digits[1]=7;  // -1 gives defaults
	  max_its[1]=1000;
	  is_gradient_supplied[1]=1;
	  is_hessian_supplied[1]=0;

	  // workspace
	  double[] fullParameters,fullParameters_best,grad_sub_best,typical_values,Hessian_diag_wk;
	  double[][] Hessian_wk;
	  fullParameters = new double[n+1];
	  fullParameters_best = new double[n+1];
	  grad_sub_best = new double[n+1];
	  typical_values = new double[n+1];
	  Hessian_wk = new double[n+1][n+1];
	  Hessian_diag_wk = new double[n+1];

	  for (i=1;i<=n;i++) typical_values[i] = 1;

	  fscale[1] = EQUISPACE_TEST ? N : _CI.nd;
	  // set start point
	  offsetCopyMtoX(fullParameters);

	  if (false)
		{
		  // test FD gradient: verifies that coding below has overall gradient correct (and, 
		  // thus, by implication, grad[] in CI at each path node)
		  float gn,fn,mo,EPS = 0.0001f;
		  for (i=0; i<=N+1; i++) 
			{
			  f[i] = f_obj_and_grad(m[i],gradf[i],true);
			  for(j=0;j<np;j++)
				{
				  mo = m[i][j];
				  m[i][j] += EPS;
				  fn = f_obj_and_grad(m[i],null,false);
				  gn = (fn-f[i])/EPS;
				  _CI.EH.println("GRADTEST %d %d %g %g %g%n",i,j,gradf[i][j],gn,abs(gradf[i][j] - gn)/abs(gn));
				  // restore
				  m[i][j] = mo;
				}
			}
		}

	  // compute f at two endpoints (always fixed)
	  f[0] = f_obj_and_grad(m[0],gradf[0],false);
	  f[N+1] = f_obj_and_grad(m[N+1],gradf[N+1],false);
	  nIts=0;

	  if (savePaths) 
		{
		  path = new ArrayList<float[][]>();
		  path_f = new ArrayList<float[]>();
		}

	  // Go!
	  Uncmin_f77.optif9_f77(n, fullParameters, this,
							typical_values,fscale,method,is_exp,msg,n_good_digits,
							max_its,is_gradient_supplied,is_hessian_supplied,trust_radius,
							gradient_tol,max_step_size,step_tol,
							fullParameters_best,f_best,grad_sub_best,
							info,Hessian_wk,Hessian_diag_wk);
	  // take copy of final answer
	  offsetCopyMtoX(fullParameters_best);

	  dumpOptimizationPathtoFiles();
	  path=null;
	  path_f=null;

	  // Dump result to a file
	  filename = "ModeLinks."+label+".txt";
	  try
		{
		  out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		} 
	  catch (IOException ioe) 
		{
		  throw new RuntimeException("couldn't open " + filename+" for output");
		}
	  out.printf("#DEPTHS ");
	  // top depths of free layers
	  for (i=0; i<(_CI.nlayer+1); i++) if (_CI.rho[i] <= 0.0f) out.printf("%10.5g ",_CI.depths[i]);
	  out.printf("%n# ");
	  // top depths of free layers
	  for (i=0,j=0; i<(_CI.nlayer+1); i++) if (_CI.rho[i] <= 0.0f) {out.printf("m_%d ",j); j++;}
	  out.printf("Lambda_AB%n");
	  // now mode-link points
	  for (i=0;i<=N+1;i++)
		{
		  for (j=0;j<np;j++) out.printf("%g ",m[i][j]);
		  out.printf("%g%n",f[i]);
		}
	  out.close(); 
	  _CI.EH.println("Wrote "+filename);
	  if (TRY_TO_EQUISPACE && _CI.verbosity>2)
		{
		  getSegmentLengths();
		  _CI.EH.println("Node spacings: ");
		  for (i=0; i<=N; i++) _CI.EH.printf("%g%n",dm2[i]); 
		}
	}

	/** 
		Objective function for Uncmin. x[] is unrolled over node points & dimensions
		to the parameters m[0..N+1][0..np-1] through offsetCopyXtoM()
	*/
	public double f_to_minimize(double[] x)  
	{
	  int i;
	  float s,Lambda_AB=0f;
	  
	  offsetCopyXtoM(x);
	  getSegmentLengths();

	  // get f's and grad components on internal points (assume i=0,N+1 precomputed)
	  if (!EQUISPACE_TEST)
		{
		  for (i=1; i<=N; i++) f[i] = f_obj_and_grad(m[i],gradf[i],false);
		  
		  // Path-integral quadrature approximation
		  for (i=0,Lambda_AB=0; i<=N; i++) Lambda_AB += (f[i]+f[i+1])*dm2[i];
		}
	  if (TRY_TO_EQUISPACE)
		{
		  // additional term to encourage nodes to remain equispaced
		  // penalise deviation from mean-spacing
		  float chisqEq;
		  for (i=0,chisqEq=0.0f; i<=N; i++) chisqEq += sqr(dm2[i] - msbar);
		  chisqEq *= EQUISPACE_A;
		  Lambda_AB = EQUISPACE_TEST ? chisqEq : (Lambda_AB+chisqEq) ;
		}
	  return (double)Lambda_AB;
	}

	// get unit--offset gradient gradL[1....]
	public void gradient(double x[], double g[])
	{
	  int i,j,k;
	  float s,Lambda_AB;
	  
	  offsetCopyXtoM(x);
	  getSegmentLengths();
	  if (!EQUISPACE_TEST)
		{	  
		  // get f's and grad components on internal points (assume i=0,N+1 precomputed)
		  for (i=1; i<=N; i++) f[i] = f_obj_and_grad(m[i],gradf[i],true);
		  
		  // Here the effective function is Lambda_AB
		  for (i=0,Lambda_AB=0; i<=N; i++) Lambda_AB += (f[i]+f[i+1])*dm2[i];
		
		  // Get gradL
		  for (k=1,i=1; i<=N; i++)
			{
			  for(j=0;j<np;j++)
				{
				  g[k] = gradf[i][j] * (dm2[i-1]+dm2[i])
					- (f[i]+f[i+1]) * (m[i+1][j] - m[i][j])/((dm2[i]>0)?dm2[i]:1.0)
					+ (f[i]+f[i-1]) * (m[i][j] - m[i-1][j])/((dm2[i-1]>0)?dm2[i-1]:1);

				  k++;
				}
			}
		}
	  if (TRY_TO_EQUISPACE)
		{
		  boolean useAnalyticGradient = true;
		  float gr=0,coeff;
		  coeff = 2*EQUISPACE_A;

		  if (useAnalyticGradient) 
			{
			  for (k=1,i=1; i<=N; i++)
				{
				  for(j=0;j<np;j++)
					{
					  // analytic - OK I think.
					  gr =  (m[i][j] - m[i+1][j])*(1f - msbar/( (dm2[i]>0) ? dm2[i] : 1f))
					  	  - (m[i-1][j] - m[i][j])*(1f - msbar/((dm2[i-1]>0) ? dm2[i-1] : 1f));
					  gr *= coeff;
					
					  if (EQUISPACE_TEST) g[k] = gr; //	EQUISPACE TEST
					  else g[k] += gr;
					  k++;
					}
				}
			}
		  else
			{
			  // Fwd diffs for gradient
			  int ii;
			  float chisqEqO=0f,chisqEq,EPS=3.16e-4f;
			  for (ii=0,chisqEqO=0.0f; ii<=N; ii++) chisqEqO += sqr(dm2[ii] - msbar);
			  for (k=1,i=1; i<=N; i++)
				{
				  for(j=0;j<np;j++)
					{
					  m[i][j] += EPS;
					  getSegmentLengths();
					  for (ii=0,msbar=0f; ii<=N; ii++) msbar += dm2[ii]; 
					  msbar/=(N+1);

					  for (ii=0,chisqEq=0.0f; ii<=N; ii++) chisqEq += sqr(dm2[ii] - msbar);
					  gr = EQUISPACE_A*(chisqEq - chisqEqO)/EPS;
					  // restore
					  m[i][j] = (float)x[k];
					
					  if (EQUISPACE_TEST) g[k] = gr; //	EQUISPACE TEST
					  else g[k] += gr;
					  k++;
					}
				}
			  getSegmentLengths(); // restore
			}
		}
	}


	// Optimisation function for findModePairConnection()

	private float f_obj_and_grad(float[] m,float[] gradient,boolean doGradient)
	{
	  int i;
	  long t;
	  double chisq;

	  // FIXME JOIN SMOOTH MODELS HOW????
	  float alogmu = -10.0f; 

	  // copy input params into global pm[]
	  System.arraycopy(m,0,_CI.pm,0,_CI.np);
	
	  if (doGradient)
		{
		  // construct wjtwj and find misfit...
		  _CI.makjtj();
		  _CI.makptp(alogmu); 
		  DenseDoubleMatrix2D A = new DenseDoubleMatrix2D(_CI.np,_CI.np);
		  _CI.computeGradandHessian(alogmu,A);
		  
		  // gradient is now in global _CI.grad[]. Feed back to caller...
		  System.arraycopy(_CI.grad,0,gradient,0,_CI.np);
		}
	  // calculate the model response and misfit
	  t = System.nanoTime();
	  _CI.doForwardModel(_CI.pm);
	  _CI.tNatives += (System.nanoTime()-t); 	
	  chisq = _CI.chisqMisFit();
	  chisq += fpow10(alogmu) * _CI.chisqRoughness(_CI.pm) + _CI.chisqPrej(alogmu,_CI.pm);
	  chisq /= 2.0;

	  

	  return (float)chisq;
	}


	private void offsetCopyXtoM(double[] x)
	{
	  int i,j,k;
	  float[] mptr;
	  for (i=1,k=1;i<=N;i++)
		{
		  mptr = m[i];
		  for(j=0;j<np;j++,k++) mptr[j] = (float)(x[k]);
		}
	}

	private void offsetCopyMtoX(double[] x)
	{
	  int i,j,k;
	  float[] mptr;
	  for (i=1,k=1;i<=N;i++)
		{
		  mptr = m[i];
		  for(j=0;j<np;j++,k++) x[k] = (double)(mptr[j]);
		}
	}

	private void getSegmentLengths()
	{
	  int i,j;
	  float s;
	  for (i=0;i<=N;i++)
		{
		  for (s=0,j=0;j<np;j++) s += sqr(m[i+1][j]-m[i][j]);
		  dm2[i] = (float)sqrt(s);
		}
	  if (TRY_TO_EQUISPACE)
		{
		  for (i=0,msbar=0f; i<=N; i++) msbar += dm2[i]; 
		  msbar/=(N+1);
		}
	}

	void dumpOptimizationPathtoFiles()
	{

	  if (!savePaths) return;

	  int i,j,l,k;
	  float [][] _m;
	  PrintWriter out=null;

	  // Parameter crossplots
	  for(i=0;i<np-1;i++)
		{
		  for(j=i+1;j<np;j++)
			{				  
			  String filename = "ModeLinks."+_label+".params="+i+"."+j+".steps.txt";
			  try
				{
				  out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
				} 
			  catch (IOException ioe) 
				{
				  throw new RuntimeException("couldn't open " + filename+" for output");
				}
			  // now mode-link points
			  for (l=0;l<path.size();l++)
				{
				  _m = (float[][])path.get(l);
				  for(k=0;k<=N+1;k++)
					{
					  out.println(_m[k][i]+" "+_m[k][j]);
					}
				  out.println("");
				}
			  out.close();
			  _CI.EH.println("Wrote "+filename);
			  
			}
		} 

	  // Function profile
	  String filename = "ModeLinks."+_label+".func.steps.txt";
	  try
		{
		  out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		} 
	  catch (IOException ioe) 
		{
		  throw new RuntimeException("couldn't open " + filename+" for output");
		}
	  for (l=0;l<path.size();l++)
		{
		  float[] _f = (float[])path_f.get(l);
		  for(k=0;k<=N+1;k++) out.println(_f[k]);
		  out.println("");
		}
	  out.close();
	  _CI.EH.println("Wrote "+filename);	  
	}
  
	
	/**
	   The solver throws execution to here everytime it gets a major new
	   point. This is where we can monitor  things.
	*/
	public void updateNewPoint(double fnew, double[] xnew)
	{
	  if (_CI.verbosity > 1) _CI.EH.println("ModeLink optimizer: f = "+fnew);
	  
	  if (savePaths)
		{
		  int i,j;
		  float[][] mcpy = new float[N+2][np];
		  for (i=0;i<N+2;i++) System.arraycopy(m[i],0,mcpy[i],0,np);
		  path.add(mcpy);
		  path_f.add((float[])f.clone());
		}
	  nIts++;
	}

	public boolean abortNow() {return false;};

	public boolean wasAborted() {return false;};

	public void setFlags(boolean setOrReset) {};

	// null methods required since we can't figure these out analytically
	public void hessian(double times[], double h[][]) 
	{};
  }
  

/**
	 Local modes from optimizer stored here
  */
  public static class LocalMode implements Comparable
  {
	float[] pmLoc;
	float pmuLoc;
	float negLogPost;
	float negLogNoSmoothingPrior;
	float rmsMisfit;
	float logMMLapprox;
	DenseDoubleMatrix2D wjtwj;   // copy of inverse covariance at optimum
	DenseDoubleMatrix2D H;   // copy of inverse covariance (hessian) at optimum
	float[] stdv;
	String resp_filename;
	String profile_filename;
	String iter_filename;

	public LocalMode(int np,String resp_filename,String profile_filename,String iter_filename,
					 float[] pm,float pmuLoc,float negLogPost,float negLogNoSmoothingPrior,
					 float rmsMisfit, 
					 // DenseDoubleMatrix2D wjtwj,
					 DenseDoubleMatrix2D H,
					 float[] stdv)
	{
	  pmLoc = (float[]) pm.clone();
	  System.arraycopy(pm,0,pmLoc,0,pm.length);
	  this.resp_filename = new String(resp_filename);
	  this.profile_filename = new String(profile_filename);
	  this.iter_filename = new String(iter_filename);
	  this.pmuLoc = pmuLoc;
	  this.negLogPost = negLogPost;
	  this.negLogNoSmoothingPrior = negLogNoSmoothingPrior;
	  this.rmsMisfit = rmsMisfit;
	  this.H = (DenseDoubleMatrix2D)(H.like().assign((DoubleMatrix2D)H));
	  this.wjtwj =(DenseDoubleMatrix2D) this.H.viewPart(0,0,np,np);
	  this.stdv = (float[]) stdv.clone();
	}
	
	/** sort modes by chisq, then fields */
	public int compareTo(Object o) 
	{ 
	  int i,l;

	  LocalMode m1 = this;
	  LocalMode m2 = (LocalMode)o;

	  if (m1.negLogPost != m2.negLogPost)  // freq
		{
		  return (m1.negLogPost < m2.negLogPost) ? -1 : 1;
		}
	  else 
		{
		  // sort by values
		  l = min(m1.pmLoc.length, m2.pmLoc.length);
		  for (i=0;i<l;i++)
			{
			  if (m1.pmLoc[i] != m2.pmLoc[i])
				return (m1.pmLoc[i] < m2.pmLoc[i]) ? -1 : 1;
			}
		}
	  // at this point, don't care...
	  return 1;
	}
  }
  
  /**
	ModePair class. Stores info about mode i,j, and a measure of "distance" between them
  */
  public static class ModePair implements Comparable
  {
	int i,j;
	float dist;

	public ModePair(int i,int j,float dist)
	{
	  this.i=i;
	  this.j=j;
	  this.dist=dist;
	}

  	/** sort modes by dist */
	public int compareTo(Object o) 
	{ 
	  int i,l;

	  ModePair mp1 = this;
	  ModePair mp2 = (ModePair)o;

	  if (mp1.dist != mp2.dist) return (mp1.dist < mp2.dist) ? -1 : 1;
	  else return  (mp1.i < mp2.i) ? -1 : 1;
	}
  }
  
  /** 
	  simple class to store commandline model study spec
  */
  public static class ModelStudyList
  {
	int layer;
	float start,inc,end;
	float current_value;
  }

  /**
	 Simple class for storing contiguous 2D arrays (COLT-like) for communications with
	 fortran. Underlying storage is a NIO float buffer. Can be used to
	 pass memory through JNI via pointers without copying. Experimental. Not used.
  */
  public static class DenseBufferedFloatMatrix2D
  {
	int rows;
	int columns;
	java.nio.FloatBuffer buf = null;

	public DenseBufferedFloatMatrix2D(int rows, int columns)
	{
	  this.rows = rows;
	  this.columns = columns;

	  buf = ByteBuffer.allocateDirect(4*rows*columns).asFloatBuffer();
	  System.out.println("FloatBuffer direct = "+buf.isDirect());
	}

	public void setQuick(int row, int column, float value) 
	{
	  // manually inlined:
	  buf.put(row*columns + column, value);
	}

	public float getQuick(int row, int column) 
	{
	  // manually inlined:
	  return buf.get(row*columns + column); 
	}

	public void mult(int row, int column, float value) 
	{
	  // manually inlined:
	  int i = row*columns + column;
	  buf.put(i, buf.get(i) * value);
	}
	
  }

  
  // Parse commandline for systematic flavours
  void setSystematicsFlavour()
  {
	int i;
	  
	if (!modelSystematics) return;

	String[] types = systematicsFlavour.split(":");
	  
	for (i=0;i<types.length;i++)
	  {
		String s = types[i];
		String[] spec = s.split(",");
		if (spec.length==3 && spec[0].equals("d"))
		  {
			SEM.systematicsRegressionOrder = Integer.parseInt(spec[1]);
			if (SEM.systematicsRegressionOrder < 0 || SEM.systematicsRegressionOrder>3) // up to cubic only for sanity
			  {
				EH.die_err("Systematic flavour d,order,ratio: accepts only 0 <= order <= 3 for safety");
			  }

			SEM.sigmaSystematicTestPoint = Double.parseDouble(spec[2]);
			if (SEM.sigmaSystematicTestPoint < 0.01) // about as small as is useful...
			  {
				EH.die_err("Systematic flavour d,order,ratio: currently, require ratio>0.01");
			  }
		  }
		else
		  {
			EH.die_err("Systematic flavours: currently, understands only: d,order,ratio. Can't parse: "+s);		  
		  }
	  }
  }
  /**
	 Simple class for keeping highest level generic info on systematic error modelling
  */
  public static class SystematicErrorsModel
  {
	public static int DEFAULT_SYSTEMATICS_REGRESSION_ORDER = 2;             // = highest power of offset in regressions
	public double sigmaSystematicTestPoint = 0.5; // relative fraction of systematic component (use overwrites)

	int systematicsRegressionOrder = DEFAULT_SYSTEMATICS_REGRESSION_ORDER;  // Parabolic(?) trend across offsets
	// extra parameters from systematics
	int nPsystematics=0;

	ArrayList<SystematicGroup> systematicGroups = null;
	DenseDoubleMatrix2D C_BsysInv=null;
	DenseDoubleMatrix2D L_Bsys=null;
	DenseDoubleMatrix2D X=null;
	double logDetInvC_Bsys=0; // logDet of C_BsysInv

	// systematic error model. Working variables
	double[] m_N = null;
	double[] m_N_bar = null;
	float[] X_NdeltaD = null;
	DenseDoubleMatrix2D H12 = null;
	DenseDoubleMatrix2D H22 = null;
  
	double getLogDetC_BsysInv()
	{
	  return logDetInvC_Bsys;
	}
  }


  /**
	 Simple class for keeping info on data points with systematic errors
  */
  public static class SystematicGroup  
  {
	float f; // freq;
	IntArrayList indexes = new IntArrayList();
	FloatArrayList ranges = new FloatArrayList();
	float rangemin,rangemax;
	DenseDoubleMatrix2D X_NT = null;
  }

  /**
	 Simple class for model select info
  */
  public static class  ModelSelectInfo 
  {
	String best_profile_filename;
	String best_resp_filename;
	String best_iter_filename;
  }

  /**
	 Simple class for storing Scripps old-style CSEM data
  */
  public static class SimpleEMData 
  {
	float d;
	float sd;
	float[] dp;

	public SimpleEMData(int ndp)
	{
	  dp = new float[ndp];
	}
  }


  /**
	 Simple class for managing multistart optimization strategies and book-keeping
  */
  static class MultiStartStrategyManager
  {
	// what to try
	boolean tryPhaseFlipping = false;
	boolean tryRandom = false;


	// bookkeeping flags
	boolean phaseFlipping = false;
	int k_flips=0;
	boolean phaseRandom = false;
	int random_ctr=0;

	public 	MultiStartStrategyManager(boolean tryPhaseFlipping, boolean tryRandom)
	{
	  this.tryPhaseFlipping = tryPhaseFlipping;
	  this.tryRandom  = tryRandom;
	}

	void initialiseStrategies()
	{
	  phaseFlipping = false;
	  phaseRandom = false;

	  k_flips = 0;
	  random_ctr = 0;

	  if (tryPhaseFlipping)
		{
		  phaseFlipping = true;
		}
	  else if (tryRandom)
		{
		  phaseRandom = true;
		}
	  // else ... add more strategies...
	}

	void exitFlipping()
	{
	  phaseFlipping = false;
	  k_flips = 0;
	  if (tryRandom)
		{
		  phaseRandom = true;
		}
	  else
		{
		  // init other strategies..
		}
	}
	
	void exitRandom()
	{
	  phaseRandom = false;
	  random_ctr = 0;
	  /* etc, for next strategy	 
	  if (someNewStrategy)
		{
		  phaseNewStrategy = true;
		}
	  */
	}
}

 

  /**
	 Simple class for sorting  Scripps old-style CSEM data
  */
  class DataComparator implements Comparator<SimpleEMData> 
  {
	public static final boolean MAKE_TYPE_SORT_LEAST_PRIORITY = false;

	// sort by 
	// dp[2] = frequency, then
	// dp[1] = range, then
	// if (MAKE_TYPE_SORT_LEAST_PRIORITY)
	//     dp[5] = trans. height.
	//     dp[0] = type.
	// else
	//     dp[0] = type.
	//     dp[5] = trans. height. 

	/**
	   Comparator to sort by frequency, then range, then type, then transmitter height
	*/
	public int compare(SimpleEMData data1, SimpleEMData data2) 
	{ 
	  if (data1.dp[2] != data2.dp[2])  // freq
		{
		  return (data1.dp[2] < data2.dp[2]) ? -1 : 1;
		}
	  else if (data1.dp[1] != data2.dp[1]) // range
		{
		  return (data1.dp[1] < data2.dp[1]) ? -1 : 1;
		} 
	  else
		{
		  if (MAKE_TYPE_SORT_LEAST_PRIORITY)
			{
			  if (((data1.dp.length>5) && (data2.dp.length >5)) 
				  && (data1.dp[5] != data2.dp[5])) // trans height
				{
				  return (data1.dp[5] < data2.dp[5]) ? -1 : 1;
				}
			  else if (data1.dp[0] != data2.dp[0]) // type
				{
				  return (data1.dp[0] < data2.dp[0]) ? -1 : 1;
				} 
			  else return 1; // don't care
			}
		  else
			{
			  if (data1.dp[0] != data2.dp[0]) // type
				{
				  return (data1.dp[0] < data2.dp[0]) ? -1 : 1;
				}
			  else if (((data1.dp.length>5) && (data2.dp.length >5)) 
					   && (data1.dp[5] != data2.dp[5])) // trans height
				{
				  return (data1.dp[5] < data2.dp[5]) ? -1 : 1;
				}
			  else return 1; // don't care
			}
		}
	}
  }
  
  /**
	 Simple class for sorting  new Scripps data by iTx, then iFreq, then type
  */
  class NSDataComparator implements Comparator<int[]> 
  {
	public int compare(int[] X, int[] Y) 
	{ 
	  if (X[2] != Y[2])                   // ascending transmitter #
		{
		  return (X[2] < Y[2]) ? -1 : 1;
		}
	  else if (X[1] != Y[1])              // freq
		{
		  return (X[1] < Y[1]) ? -1 : 1;
		} 
	  else if (X[0] != Y[0])              // type
		{
		  return (X[0] < Y[0]) ? -1 : 1;
		} 
	  else if (X[3] != Y[3])              // sort by ascending receiver #
		{
		  return (X[3] > Y[3]) ? -1 : 1;
		}
	  else return 1; // don't care
	}
  }
}