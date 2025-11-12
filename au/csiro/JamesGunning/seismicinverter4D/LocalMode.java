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



package au.csiro.JamesGunning.seismicinverter4D;

import java.util.*;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.util.ArrayUtils;
import cern.colt.matrix.*;

/**
	 Class to encapsulate properties of a local mode. A mode is a
	 characterised by a mean vector, and a local covariance (or
	 'scale' matrix) whose Cholesky decomposition is of form: R.R^T =
	 C.  Use of the mean is subtle. The mode may only have meaning in
	 some subset of the space where the layer times satisfy a
	 particular ordering (`pinchout configuration'), and in the high
	 probability part of the space satisfying this constraint, there
	 may not be a maximum away from the boundaries defined by the
	 pinchout configuration. (i.e. the probability `mountain' has
	 meaning only on it's shoulders say - the top region may be cut
	 off or illegal.). A Gaussian approximation to the mode may then
	 have its centre outside the legal region, and this point may be
	 quite different from the maximum likelihood point on the region
	 boundary. The mode stores information about both the hypothetical
	 centre of the quadratic approximation to the mode
	 neg-log-likelihood, and also the maximum likelihood point in the
	 legal region.
  */
public final class LocalMode implements Comparable
{
  /**
	 Approximate posterior weight of mode, perhaps incorporating posterior fits
	 etc. Initially stored as a negative log-likelihood, then sorted
	 and normalised to a conventional absolute +ve probability. Sum
	 over all modes = 1
  */
  protected double weight;

  /**
	CDF of this mode from posterior weights (to be computed later). Assembled from the
	normalised cdf constructed from "weight"
  */
  protected double CDF=0.0;
  /**
	 Prior weight of mode, computed *only* from fluid-combination prior probabilities etc.
	 Used always as a conventional absolute +ve probability. Sum over all modes = 1
  */
  protected double prior_weight;
  /**
	 Hybrid weight of mode, equal to 0.5(prior_weight+weight)
  */
  protected double hybrid_weight;
 /**
	CDF of this mode from prior weights (to be computed later)
  */
  protected double prior_CDF=0.0;
  /**
	 Chisq of mode, i.e. value of neg-log-likelihood omitting
	 determinant terms and prior fluid combination contribution
	 (the function that is optimised in mode hunting)
  */
  protected double chisq;
  /**
	 Seismic chisq of mode, i.e. value of neg-log-likelihood 
	 coming entirely from quality of fit to seismic. Used as an indicator of
	 being trapped in local minima and => in criteria
	 to initiate global searches
  */
  protected double seismic_chisq;

  /**
	   Fluid state of this mode
	*/
  protected int fs;

  /**
	 Ptr to fluid state of this mode
  */
  protected FluidStates.FluidState FS;

  /**
	   "Mean" vector at mode. If the mode is not pinched out, this
	   will be the true minimum of the quadratic approximation the the
	   smooth posterior. If the minima exists at a pinchout, This will
	   be the optimal point on the {t_i} boundary (and therefore not
	   the centre of the best quadratic approximation to the posterior
	   in the acceptable time subspace - i.e. that set of times
	   yielding the correct layer-code) 
  */
  protected DoubleMatrix1D QP_optimum;

  /**
	   See the blurb for "QP_optimum". This point is the true quadratic centre of the 
	   best quadratic approximation to the posterior. It will likely exist in an illegal
	   time state though, so evaluation at this point will produce the wrong layer-code.
	*/
  protected DoubleMatrix1D quadratic_mean;

  /**
	   Cholesky right decompositon of local-covariance or scale matrix
	*/
  protected DoubleMatrix2D R;
  /** 
	  Diagonal entries in the Cholesky Right factor of the inverse posterior covariance, stored
	  for all parameters
  */
  protected double[] CtildeInvDiag;

  /**
	 viewDice backed transpose of R (Cholesky right decompositon of local-covariance or scale matrix)
	 This is in fact the left cholesky factor one would normally use to draw samples
  */
  protected DoubleMatrix2D Rtransp;

 /**
	inverse of Rtransp (Linv) used in some samplers (MALA)
  */
  protected DoubleMatrix2D RtranspInv=null;
  /**
	   log(Determinant of local_covariance)
  */
  protected double logDetCtilde;
  
  /**
	 log(Determinant of prior covariance)
  */
  protected double logDetCprior;

  /**
	 -log(marginal model likelihood), used 
	 for between-delivery comparisons. Not normalised, and includes
	 all contributions from data, priors etc.
  */
  protected double logMargModelLikelihood;

  /**
	   vector of scaled master-layer time-thickness uncertainties
	   (std dev/mean)
	*/
  protected double[] time_uncertainty;

  /**
	 Information about the properties that are good candidates for big steps
  */
  protected ArrayList bigStepInfo;

  /**
	   layer_code is a simple binary encoding representation of the configuration of
	   existing layers (i.e. not pinched out) where the derivates etc
	   used to form the mode are valid
	*/
  protected int layer_code;

  /**
	 Flag to signify initialised state in Hopping sampler
  */

  protected boolean initialised=false;

  /**
	   Some kinds of samplers need to remember the last accepted state of this mode, so
	   this is a reference to it
  */
  protected double[] last_accepted_state;

  /**
	   Some kinds of samplers need to remember the last accepted values of the
	   reflection coefficients for this mode, so these are stored here
  */
  protected double[][] last_accepted_R_pp;
  /**
	 Ditto for contacts here [stack][indx]
  */
  protected double[][] last_accepted_R_pp_fc;
  /**
	 Ditto for times. Here [vintage][indx]
  */
  protected double[][] last_accepted_t_fc;

  /**
	   Some kinds of samplers need to remember the last accepted negative Log-likelihood
	   for isopach+seismic+depth+VS:
  */
  protected double lastAcceptedNegLogAllLikelihoods;
  /**
	   Some kinds of samplers need to remember the last accepted negative Log-Prior density
  */
  protected double lastAcceptedNegLogLikPriorOnly;
  /**
	   Some kinds of samplers need to remember these normailisations constants
  */
  protected double normalisationConstants;

  /**
	 possible list of big jump directions
  */
  protected ArrayList bigJumpDirections=null;


  protected LocalMode(double weight, double chisq, int fs, FluidStates.FluidState FS, int layer_code, DoubleMatrix1D QP_optimum, DoubleMatrix2D R)
  {
	this.weight = weight;
	this.chisq = chisq;
	this.fs = fs;
	this.QP_optimum = QP_optimum;
	this.R = R;
	this.layer_code = layer_code;
	this.FS = FS;
  }
	
  /**
	 Instantiator amounting to rough deep copy. Watch out
	 for uninitialised fields
  */
  protected LocalMode(LocalMode LM)
  {
	this.weight = LM.weight;
	this.CDF = LM.CDF;
	this.prior_weight = LM.prior_weight;
	this.hybrid_weight = LM.hybrid_weight;
	this.prior_CDF = LM.prior_CDF;
	this.chisq = LM.chisq;
	this.seismic_chisq = LM.seismic_chisq;
	this.fs = LM.fs;
	this.FS = LM.FS;
	this.QP_optimum = LM.QP_optimum;
	this.quadratic_mean = LM.quadratic_mean;
	this.R = LM.R;

	this.Rtransp = LM.Rtransp;
	this.RtranspInv = LM.RtranspInv;
	this.logDetCtilde = LM.logDetCtilde;
	this.logDetCprior = LM.logDetCprior;
	this.logMargModelLikelihood = LM.logMargModelLikelihood;

	this.bigStepInfo = null;
	this.layer_code = LM.layer_code;
	this.initialised = LM.initialised;

	this.lastAcceptedNegLogAllLikelihoods = LM.lastAcceptedNegLogAllLikelihoods;
	this.lastAcceptedNegLogLikPriorOnly = LM.lastAcceptedNegLogLikPriorOnly;
	this.normalisationConstants = LM.normalisationConstants;
	this.bigJumpDirections = null;

	this.CtildeInvDiag = (LM.CtildeInvDiag != null) ? (double[])LM.CtildeInvDiag.clone() : null;
	this.time_uncertainty = (LM.time_uncertainty != null) ? (double[])LM.time_uncertainty.clone() : null;
	this.last_accepted_state = (LM.last_accepted_state != null) ? (double[])LM.last_accepted_state.clone() : null;

	this.last_accepted_R_pp = (LM.last_accepted_R_pp != null) ? ArrayUtils.deepCopy(LM.last_accepted_R_pp) : null;
	this.last_accepted_R_pp_fc = (LM.last_accepted_R_pp_fc != null) ? ArrayUtils.deepCopy(LM.last_accepted_R_pp_fc) : null;
	this.last_accepted_t_fc = (LM.last_accepted_t_fc != null) ? ArrayUtils.deepCopy(LM.last_accepted_t_fc) : null;
  }


  public double getPriorFSAdjustedWeight()
  {
	return logMargModelLikelihood - Math.log(FS.prob);
  }


  /** 
	 Comparator to sort modes by prior_weight
  */
  public int compareTo(Object o)
  {
	LocalMode n = (LocalMode)o;
	int cmp = (prior_weight<n.prior_weight)? -1 : ((prior_weight>n.prior_weight)? 1 : 0 );
	return cmp;
  }


  public final static class LocalModeSorter implements java.util.Comparator
  {
	/** 
		Class for sorting modes on basis of weight. Also possible to
		sort on fs first by calling the appropriate flag in the constructor
	*/
	private boolean sortBy_fs_First = false;
	private boolean sortBy_prior_weight_not_weight = true;

	public LocalModeSorter(boolean sortBy_fs_First, boolean sortBy_prior_weight_not_weight)
	{
	  this.sortBy_fs_First = sortBy_fs_First;
	  this.sortBy_prior_weight_not_weight = sortBy_prior_weight_not_weight;
	}

	/**
	   Sort on fs first, then weight (assume weights are -log(weight))
	*/
	public int compare(Object o1, Object o2)
	{
	  LocalMode n1 = (LocalMode)o1;
	  LocalMode n2 = (LocalMode)o2;
	  int cmp;
	  double w1,w2;

	  if (sortBy_prior_weight_not_weight)
		{
		  w1 = n1.prior_weight;
		  w2 = n2.prior_weight;
		}
	  else
		{
		  w1 = n1.weight;
		  w2 = n2.weight;
		} 
	  if (sortBy_fs_First)
		{
		  cmp = (n2.fs<n1.fs)? 1 : ((n2.fs>n1.fs)? -1 : 
									((w2<w1)? 1 : ((w2>w1)? -1 : 0)));
		}
	  else
		{
		  cmp = ((w2<w1)? 1 : ((w2>w1)? -1 : 0));
		}
	  return cmp;
	}
  }
}
