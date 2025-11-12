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



package au.csiro.JamesGunning.lib.stat.util;

import java.util.*;



import cern.colt.matrix.*;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;
import cern.jet.random.ChiSquare;

import au.csiro.JamesGunning.lib.util.MathUtil;
import au.csiro.JamesGunning.lib.util.Util;
import au.csiro.JamesGunning.lib.ErrorHandler;

/**
   <p>  A class designed to monitor convergence in two possible ways.<br>
   <p>
   (1) as per suggestions on p 331
   of Gelman et al "Bayesian Statistics". Read these pages to make
   sense of this Class. The discussion is very similar to Ch. 8 of
   "Markov Chain Monte Carlo in Practice", Gilks et al; entitled
   "Inference and monitoring convergence".  For the time being, we use
   internal layer times as the scalar estimands.  These seem to be
   most susceptible to non-stationarity, or slow convergence.<p>
   <p>
   (2) There are also some utility routines for computing the most optimal
   decimation of a single sequence, which we call Raftery-style decimation,
   as per Raftery et al, Ch 7, "Markov Chain Monte Carlo in Practice", Gilks et al. 
   See also Raftery, A.E. and Lewis, S.M. (1992) "How many
   iterations in the Gibbs sampler?"  In Bayesian Statistics, Vol. 4
   (J.M. Bernardo, J.O.  Berger, A.P. Dawid and A.F.M. Smith,
   eds.). Oxford, U.K.: Oxford University Press, 763-773.  This paper
   is available via the World Wide Web by linking to URL
   http://www.stat.washington.edu/tech.reports/pub/tech.reports and
   then selecting the "How Many Iterations in the Gibbs Sampler" link.
   A useful reference is also
   Bishop, Feinberg & Holland, "Discrete Multivariate Analysis" MIT Press.
   <br>
   Here we have freely translated a few salient routines from the gibbsit
   code - mainly the calls required to compute the thinning constant. 
   <p>
*/
public final class MCMCConvergenceMonitor
{
  // Error Handler
  private ErrorHandler _EH;
  

  private int verbosity=0;

  /** 
	State variable, set to true for Gelman style diagnostics
  */
  public boolean GelmanStyleDiagnostics = false;
  /** 
	State variable, set to true for Raftery style diagnostics
  */
  public boolean RafteryStyleDiagnostics = true;


  /**
	 Array of scalar estimands for convergence checking. Used for both Gelman and Raftery styles
	 Typically the layer times.
  */
  private double[] scalar_estimands;

  /**
	 Number of independent chains to run in Gelman style diagnostics (Gelmans' J, p 331)
  */
  public static final int N_INDEPENDENT_CHAINS = 5;
  
  // useful stuff for acquiring these stats
  private double[] sum_psi_sq;
  private DenseDoubleMatrix2D psiBar_dotj;

  private int ensemble_cntr;
  
  /*********************************************************
  Variables for the Raftery style gibbsit routines
  *********************************************************/

  /**
	 this will store the requisite continuous sequences we wish to examine (typically
	 times, which are thresholded...)
  */
  DoubleMatrix2D estimands_data;
  public int[] Raftery_kthin;
  private double g2,bic;
  private double alpha, beta;
  private int kthin_max=1;
  private boolean isRafteryCalibrationRun = false;

  /**
	 The method getRafteryWorstCaseThinningRecommendation()
	 recommending extra decimation will be forced to recommend nothing
	 greater than this constant
  */
  public static final int MAXIMUM_RAFTERY_DECIMATION_FACTOR = 50;

  /**
	 Make a new MCMCConvergenceMonitor object.
	 @param doGelmanStyleDiagnostics If true, do Gelman style multiple random walk monitoring with
	 scalar estimands, or if false, do Raftery style monitoring like the gibbsit program
  */
  public MCMCConvergenceMonitor(ErrorHandler EH, int verbosity,
						 boolean doGelmanStyleDiagnostics, boolean doRafteryStyleDiagnostics)
  {
	this._EH = EH;
	this.verbosity = verbosity;
	// style of diagnostics
	this.GelmanStyleDiagnostics = doGelmanStyleDiagnostics;
	this.RafteryStyleDiagnostics = doRafteryStyleDiagnostics;
  }


  /**
	 Compute some convergence diagnostics as per suggestions on p 331 of Gelman et al 
	 "Bayesian Statistics". Read these pages to make sense of this Class. The discussion is
	 very similar to Ch. 8 of "Markov Chain Monte Carlo in Practice", Gilks et al;
	 entitled "Inference and monitoring convergence".
	 We use i,j, in Gelman's sense, and k to run over the scalar estimands of choice.
	 For the time being, we use internal layer times as the scalar estimands.
	 These seem to be most susceptible to non-stationarity, or slow convergence. Here we initilise things.
  */
  public void initialiseGelmanStyleConvergenceDiagnostics(int n_estimands, int n_samples)
  {
	
	// Let's choose as estimands the set of internal layer times.
	scalar_estimands = new double[n_estimands];
	if (GelmanStyleDiagnostics)
	  {
		psiBar_dotj = new DenseDoubleMatrix2D(N_INDEPENDENT_CHAINS,n_estimands);
		sum_psi_sq = new double[n_estimands];
	  }
	// there *may* be some Raftery style stuff going on as well, so we'll need this
	if (RafteryStyleDiagnostics)
	  {
		// save the estimands data here
		// calibration run is the same size as the -N ... flag 
		// (_m._SI.n_samples_per_trace_undecimated)
		estimands_data = new DenseDoubleMatrix2D(scalar_estimands.length,n_samples);
		Raftery_kthin = new int[scalar_estimands.length];
	  }
  }

  /**
	 Compute some convergence diagnostics as per suggestions on p 331 of Gelman et al 
	 "Bayesian Statistics". Read these pages to make sense of this Class. The discussion is
	 very similar to Ch. 8 of "Markov Chain Monte Carlo in Practice", Gilks et al;
	 entitled "Inference and monitoring convergence".
	 We use i,j, in Gelman's sense, and k to run over the scalar estimands of choice.
	 For the time being, we use internal layer times as the scalar estimands.
	 These seem to be most susceptible to non-stationarity, or slow convergence. <p>
	 Here we finalise the calculation.
  */
  public void finishGelmanStyleConvergenceDiagnostics(int n_samples_per_chain)
  {
	int j,k,n;
	
	double sum;
	double[] psiBar_dotdot;
	double[] B,W;        // Will hold the between and within sequence variances as per Gelman
	double[] varPlus;    // will hold eqn 11.2 in Gelman for k=0...scalar_estimands.length-1 estimands

	// between and within-sequence stats
	W =  new double[scalar_estimands.length];
	B =  new double[scalar_estimands.length];
	varPlus =  new double[scalar_estimands.length];

	// psiBar_dotj = new DenseDoubleMatrix2D(N_INDEPENDENT_CHAINS,scalar_estimands.length);
	// sum_psi_sq = new double[scalar_estimands.length];
	psiBar_dotdot = new double[scalar_estimands.length];
	

	// Do any Raftery stuff we may require. This will act on the last run in the Gelman ensemble.
	if (RafteryStyleDiagnostics) finishRafteryStyleConvergenceDiagnostics();

	// Now attend to Gelman stuff
	// finish off psiBar_dotj  & psiBar_dotdot calculation
	// number of samples n is from second half of chain
	n = n_samples_per_chain;
	_EH.warn("Scalar estimands based on "+ n_samples_per_chain +" samples, "+ N_INDEPENDENT_CHAINS+" chains");
	for (k=0; k<scalar_estimands.length ; k++)
	  {
		for (j=0, psiBar_dotdot[k]=0.0 ; j<N_INDEPENDENT_CHAINS ; j++)
		  {
			psiBar_dotj.setQuick(j,k, psiBar_dotj.getQuick(j,k)/(double)n);
			psiBar_dotdot[k] += psiBar_dotj.getQuick(j,k);
		  }
		psiBar_dotdot[k] /= (double)N_INDEPENDENT_CHAINS;
	  }
	// now get B's and W's
	for (k=0; k<scalar_estimands.length ; k++)
	  {
		sum = B[k] = 0.0;
		for (j=0; j<N_INDEPENDENT_CHAINS ; j++)
		  {
			sum += MathUtil.sqr(psiBar_dotj.getQuick(j,k)); 
			B[k] += MathUtil.sqr(psiBar_dotj.getQuick(j,k) - psiBar_dotdot[k]);
		  }
		W[k] = (sum_psi_sq[k] - n*sum)/((double)(N_INDEPENDENT_CHAINS*(n-1)));
		B[k] *= n/((double)(N_INDEPENDENT_CHAINS-1));

		varPlus[k] = (n-1.0)*W[k]/n + B[k]/n;
		double sqrtR = Math.sqrt(Math.max(varPlus[k]/W[k], 0));
		_EH.warn("Estimand "+k+" = "+sqrtR+", ens average = "+psiBar_dotdot[k]);
		// things you might want also.........
		// _EH.warn("                               B stddev = "+ Math.sqrt(B[k]));
		// _EH.warn("                               W stddev = "+ Math.sqrt(W[k]));
		for (j=0; j<N_INDEPENDENT_CHAINS ; j++)
		  {
			_EH.warn("Chain"+j+" average = "+psiBar_dotj.getQuick(j,k));
		  }
	  }

	// clean up
	scalar_estimands=null;
	B=null;
	W=null;
	varPlus=null;
	psiBar_dotdot=null;
	psiBar_dotj=null;
	
	// raftery
	estimands_data = null;
  }


  /**
	 Here we can set the psi scalar estimands of index k to the value
  */
  public void setEstimand(int indx, double value)
  {
	if (indx>=scalar_estimands.length) 
	  throw new RuntimeException("indx greater than scaler estimands.length");
	scalar_estimands[indx] = value;
  }

  /**
	 Here we complete operations on the psi scalar estimands. The
	 scalar_estimands[] array is expected to be filled with the
	 desired values by this stage.  <p> If doing Gelman monitoring, as
	 per Gelman p 331, the indices i (current position in the chain)
	 and j (current member of ensemble (j=1...N_INDEPENDENT_CHAINS,
	 (we use ensemble_cntr)) are taken as implicit.  <p> If doing
	 Raftery-style, the desired scalars are stored for later use
	 @param sample_number is the sample in the current sequence, must
	 be 0....(_m._SI.n_samples_per_trace_undecimated-1)
  */
  public void completeScalarEstimandsUpdate(int sample_number)
  {
	int k;

	if (GelmanStyleDiagnostics)
	  {
		// Gelman style....
		// Acquire/compute the stats on p 331 of Gelman:
		for (k=0; k<scalar_estimands.length ; k++)
		  {
			psiBar_dotj.setQuick(ensemble_cntr,k, psiBar_dotj.getQuick(ensemble_cntr,k) + scalar_estimands[k]);
			sum_psi_sq[k] += scalar_estimands[k] * scalar_estimands[k];
		  }
	  }
	if (RafteryStyleDiagnostics)
	  {
		// Raftery style. Store the data for now.
		for (k=0; k<scalar_estimands.length ; k++) 
		  {
			estimands_data.setQuick(k,sample_number,scalar_estimands[k]);
		  }
	  }
  }






  /******************************************************************************/
  //                                    Raftery routines
  /******************************************************************************/

  /**
	 After a successful call to
	 performConvergenceDiagnosticsWhileSamplingFromLocalModes(), if
	 RafteryStyleDiagnostics==true (we are doing Raftery style), this
	 method reveals the recommended thinning parameter for the Markov
	 chain.  No recommendation greater than
	 MAXIMUM_RAFTERY_DECIMATION_FACTOR will be returned.
  */
  public int getRafteryWorstCaseThinningRecommendation()
  {
	// retrn the worst case thinning parameter required
	return Math.min(kthin_max,MAXIMUM_RAFTERY_DECIMATION_FACTOR);
  }

  /**
	 <p>
	 <b>Raftery:</b>
	 <p>
	 call this method to switch the calibration flag: typically set to true
	 before a calibration run, then false for the real run, with
	 thinning parameters computed from the calibration
  */
  public void setRafteryCalibrationRun(boolean state)
  {
	isRafteryCalibrationRun = state;
  }

  /**
	 Call this method to see the calibration flag. if false, we
	 may want to suppress output etc
  */
  public boolean isRafteryCalibrationRun()
  {
	return isRafteryCalibrationRun;
  }




  /**
	 Compute some convergence diagnostics as per Raftery-style
	 suggestions in "Markov Chain Monte Carlo in Practice", Gilks et
	 al; Chapter 7. See also Raftery, A.E. and Lewis, S.M. (1992).
	 How many iterations in the Gibbs sampler?  In Bayesian
	 Statistics, Vol. 4 (J.M. Bernardo, J.O.  Berger, A.P. Dawid and
	 A.F.M. Smith, eds.).  We use internal layer times as the scalar
	 estimands - these are thresholded at suitable values to produce
	 the 0,1 sequences required.  These seem to be most susceptible to
	 non-stationarity, or slow convergence<p>
	 Initialisation step.
  */
  public void initialiseRafteryStyleConvergenceDiagnostics(int n_estimands, int n_samples_in_chain)
  {
	// Let's choose as estimands the set of internal layer times.
	scalar_estimands = new double[n_estimands];

	// save the estimands data here
	// calibration run is the same size as the -N ... flag 
	// (_m._SI.n_samples_per_trace_undecimated)
	estimands_data = new DenseDoubleMatrix2D(scalar_estimands.length,n_samples_in_chain);
	Raftery_kthin = new int[scalar_estimands.length];
  }




  /**
	 Compute some convergence diagnostics as per Raftery-style
	 suggestions in "Markov Chain Monte Carlo in Practice", Gilks et
	 al; Chapter 7. See also Raftery, A.E. and Lewis, S.M. (1992).
	 How many iterations in the Gibbs sampler?  In Bayesian
	 Statistics, Vol. 4 (J.M. Bernardo, J.O.  Berger, A.P. Dawid and
	 A.F.M. Smith, eds.).  We use internal layer times as the scalar
	 estimands - these are thresholded at suitable values to produce
	 the 0,1 sequences required.  These seem to be most susceptible to
	 non-stationarity, or slow convergence <p>
	 Finish statistics calculations.
  */
  public void finishRafteryStyleConvergenceDiagnostics()
  {
	int i,j;
	double min,max,thresh;
	int[] indicator_data;

	// Now analyse the data
	// set thresholds at the middle of the observed range for all
	// estimands
	indicator_data = new int[estimands_data.columns()];
	for (i=0;i<scalar_estimands.length;i++)
	  {
		// compute data max & min
		min=Double.MAX_VALUE;
		max=-min;
		DoubleMatrix1D data = estimands_data.viewRow(i);
		for (j=0;j<estimands_data.columns();j++)
		  {
			min = Math.min(min,data.getQuick(j));
			max = Math.max(max,data.getQuick(j));
		  }
		// compute threshold as midpoint of data range
		thresh = 0.5*(max+min);
		if (verbosity>3) _EH.warn("Estimand "+(i+1)+", test threshold = "+thresh);
		// form indicator sequence
		for (j=0;j<estimands_data.columns();j++)
		  {
			indicator_data[j] = (data.getQuick(j) <= thresh) ? 1:0;
		  }
		Raftery_kthin[i] = getThinningParameter(indicator_data);
		// keep track of worst thinning required
		kthin_max = Math.max(Raftery_kthin[i], kthin_max);
	  }
	indicator_data = null;
	estimands_data = null;
  }


  /**
	 Find kthin, the degree of thinning at which the indicator series is
	 first-order Markov. data[] must be an 0,1 indicator sequence.
	 N.B. root searching algorithms look sequentially from left since the bic
	 has many zero crossings, so e.g. bisection isn't good enough to find the
	 *first* zero-crossing.
  */
  private int getThinningParameter(int data[])
  {
	int thinned_data[];
	int kmind,kthin;

	for (kthin=1;;kthin++)
	  {
		thinned_data = thin(data,kthin);
		mctest(thinned_data);
		if (bic<=0.0 || kthin>=data.length) break;
	  }
	/*
	  Calculate both the alpha and beta transition probabilities (in the  
	  Cox & Miller parametrization) of the two state first-order Markov   
	  chain determined above. This computes the alpha, beta.                                             
	*/                                                                
	mcest(thinned_data);
	/*
	  Now compute just how big the spacing needs to be so that a thinned  
	  chain would no longer be a Markov chain, but rather would be an     
	  independence chain.  This thinning parameter must be at least as    
	  large as the thinning parameter required for a first-order Markov   
	  chain.                                                              
	*/
	for (kmind=kthin;;kmind++)
	  {
		if (kmind>kthin) thinned_data = thin(data,kmind);
		indtest(thinned_data);
		if (bic<=0 || kmind>=thinned_data.length) break;
	  }

	if (verbosity>2)
	  {
		_EH.warn("Raftery auto-decimation: estimated thinning parameter = "+kmind);
	  }
	//Util.printfileZO("series.txt",data,data.length);
	//Util.printfileZO("thinseries.txt",thinned_data,thinned_data.length);
	return kmind;
  }

  /**
	 From the indicator sequence data[]=0,1, form & return the sequence obtained
	 by subsampling every kthin'th sample
  */
  private int[] thin(int data[],int kthin)
  {
	int i,j,sz;
	sz = (data.length-1)/kthin+1;
	int[] thinned_data = new int[sz];
	for (i=0,j=0 ; i<data.length && j<sz ; i+=kthin,j++) thinned_data[j]=data[i];
	return thinned_data;
  }



  /**
	 mctest function from gibsit. Computes bic and g2.
	 @param data is an indicator array of 0,1's
  */
  private void mctest(int[] data)
  {
	double fitted;
	double focus;
	int   i,j,k;

	//  Initialize the transition counts array to all zeroes.              
	int[][][] tran = new int[2][2][2]; // will be preset to zeros
	/*
	  Count up the number of occurrences of each possible type of        
	  transition.  Keep these counts in the transition counts array.     
	*/
	for (i=2;i<data.length;i++) tran[data[i-2]][data[i-1]][data[i]]++;
	/*
	  Compute the log likelihood ratio statistic for testing a second-   
	  order Markov chain against only a first-order Markov chain.  This  
	  is Bishop, Feinberg and Holland's G2 statistic.                    
	*/
	g2 = 0.0;

	for(i=0;i<=1;i++)
	  {
		for(j=0;j<=1;j++) 
		  {
			for(k=0;k<=1;k++) 
			  {
				if (tran[i][j][k]==0) continue;
				fitted = (double)(tran[i][j][0] + tran[i][j][1]) * (tran[0][j][k] + tran[1][j][k]) 
				  /(double)(tran[0][j][0] + tran[0][j][1] + tran[1][j][0] + tran[1][j][1]);
				focus = (double)tran[i][j][k];
				g2 += Math.log(focus/fitted) * focus;
			  }
		  }
	  }
	g2 *= 2.0;

	// Finally calculate the associated bic statistic
	bic = g2 - 2.0*Math.log(data.length-2);
  }
  /**
	 <p>
	 <b>Raftery:</b>
	 <p>
	 mcest function from gibsit. Computes bic and g2.
	 @param data is an indicator array of 0,1's
  */

  private void mcest(int[] data)
  {
	int   i;

	//  Initialize the transition counts array to all zeroes.              
	int[][] tran = new int[2][2]; // will be set to zeros
	// Arrays.fill(tran,0);
	/*
	  Count up the number of occurrences of each possible type of        
	  transition.  Keep these counts in the transition counts array.     
	*/
	for(i=1;i<data.length;i++) tran[data[i-1]][data[i]]++;
	/*
	  Calculate the empirical transition probabilities between 0's and
	  1's in the input (returned in alpha) and between 1's and 0's in the
	  input (returned in beta).
	*/
	alpha =(double)tran[0][1]/(double)(tran[0][0] + tran[0][1]);
	beta = (double)(tran[1][0])/(double)(tran[1][0] + tran[1][1]);
  }




  /**
	 <p>This subroutine tests for an independence chain against a first-
	 order Markov chain using the log-linear modeling formulation.  In
	 our case the independence model is the [1][2][3] model, while the
	 first-order model is the [12][23] model.  Both the [1][2][3] and
	 the [12][23] models have closed form estimates - see Bishop,
	 Feinberg and Holland (1975).
	 <p>
	 @param data an integer vector containing the series of 0's and 1's
	 for which this subroutine is to determine whether an
	 independence chain is sufficient or whether a first-
	 order Markov chain is needed to model the data.<br>

	 <p>	Computes:<br>
	 g2       = the log likelihood ratio statistic for
	 testing a first-order Markov chain against simply an
	 independence chain.  Bishop, Feinberg and Holland
	 denote this statistic as G2.<br>

	 bic      = the BIC value corresponding to the log
	 likelihood ratio statistic, g2.
	 <p>
  */
  private void indtest(int[] data)
  {
	double fitted, focus,dcm1;
	int[][] tran = new int[2][2];
	int i,j;

	/*
	  Count up the number of occurrences of each possible type of
	  transition.  Keep these counts in the transition counts array.
	*/
	for(i=1;i<data.length;i++) tran[data[i-1]][data[i]]++;
	/*
	  Compute the log likelihood ratio statistic for testing a first-
	  order Markov chain against simply an independence chain.  This is
	  Bishop, Feinberg and Holland's G2 statistic.
	*/
	dcm1 = data.length-1.0;
	g2 = 0.0;
    for(i=0;i<=1;i++)
	  {
		for(j=0;j<=1;j++)
		  {
            if (tran[i][j]==0) continue;
            fitted = ((tran[i][0] + tran[i][1])*(tran[0][j]+tran[1][j]))/dcm1;
			focus = (double)tran[i][j];
			g2 += Math.log(focus/fitted)*focus;
		  }
	  }      
	g2 *= 2.0;
	/*
	  Finally calculate the associated bic statistic and return to the
	  caller.  Note that the first-order Markov chain model contains just
	  one more parameter than does the independence chain model, so p=1.
	*/
	bic = g2 - Math.log(dcm1);
	tran=null; // GC
  }
}
