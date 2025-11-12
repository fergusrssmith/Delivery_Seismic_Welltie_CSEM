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

import au.csiro.JamesGunning.lib.Minpack.optimization.*;
import au.csiro.JamesGunning.lib.util.MathUtil;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.DifferentialEvolution.*;


import java.util.*;

import cern.colt.matrix.*;
import cern.jet.random.engine.MersenneTwister;

import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix1D;
import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.*;

/**
	 Class to find optimal (maximum-likelihood) layer times via genetic
	 algorithms (Differential Evolution)

	 The function to minimise is the sum of the negative log-likelihood error
	 terms from the prior, the seismic mismatch, and the isopach requirement.
  */	 
public final class GeneticOptimalTimesFinder implements DE_methods 
{
  Model _m;
  ModelSampler _MS;                          // reference to _m.MS

  private int n_sub;                         // dimension of current times sub-vector

  private double negLogIsopachLikelihood;    // terms available via getters() to get pieces of overall posterior
  private double negLogDepthLikelihood;      // terms available via getters() to get pieces of overall posterior
  private double negLogSeisLikelihood;       // these are set by each f_to_minimise call
  private double negLogPriorLikelihood;     
  private double negLogTruncationLikelihood; // term to force solutions away from problem time boundaries

  private double f_min;                      // snapshot of current best f
  private DoubleMatrix1D times;              // contains snapshot of current times

  /**
	 Constructor.
	 nonTimeProperties must be an ArrayList of PropertyMapper objects, pointing to
	 suitable layer properties to optimise over
  */

  public GeneticOptimalTimesFinder(Model m)
  {
	this._m = m;
	this._MS = _m.MS;

	// Set up size of time sub problem immediately
	// n_sub = _m.n_layers-1;
	n_sub = _m.n_active_layers-1;
	times = new DenseDoubleMatrix1D(n_sub);
  }

  /**
	 Find optimal solution using DifferentialEvolution class with non-time properties fixed
	 at their means after isopach incorporation
  */
  
  public double findOptimalSolution()
  {
	int k=0;

	// set outer times at prior+isopach posterior means
	// top
	_MS.property_sample[0] = _MS.mu.getQuick(0);
	_MS.setValuesFromIndex(0,_MS.property_sample[0]);
	// base
	_MS.property_sample[_m.n_active_layers] = _MS.mu.getQuick(_m.n_active_layers);
	_MS.setValuesFromIndex(_m.n_active_layers,_MS.property_sample[_m.n_active_layers]);

	// set non-time parameters at their prior+isopach posterior means
	for (k=_m.n_active_layers+1;k<_MS.mu.size();k++)
	  {
		_MS.property_sample[k] = _MS.mu.getQuick(k);
		_MS.setValuesFromIndex(k,_MS.property_sample[k]);
	  }

	// order the times
	_MS.orderTimes();
	// generate effective impedance etc for all layers
	_m.computeEffectivePropertiesForAllLayers(true,false,false,false,false); // FCOK

	// instantiate Genetic solver and solve!

	DifferentialEvolution DE 
	  = new DifferentialEvolution(n_sub, (_m._SI.useSeed ? _m._SI.seed : MersenneTwister.DEFAULT_SEED) , _m._SI.EH);

	// major control knobs...
	// Negative numbers signify internal defaults
	DE.setControlParameters(-1.0,-1.0,-1,-1.0,-1);
	// DE.setControlParameters(-1.0,-1.0,-1,-1.0,100*n_sub);

	DE.minimize(this);

	// ensure best solution is in place
	f_to_minimize(times);
	
	// feedback
	if (_m._SI.verbosity>3)
	  {
		_m._SI.EH.printf("Genetic method: Best match at overall f = %10.8g, prior = %10.8g, isopach = %10.8g, "+
						"depth-self-cons = %10.8g, seismic NLL =  %10.8g%n",
						f_min, negLogPriorLikelihood,
						negLogIsopachLikelihood,
						negLogDepthLikelihood,
						negLogSeisLikelihood);
		_m._SI.EH.print("Genetic method: best match times: ");
		for (k=1;k<=_m.n_active_layers-1;k++) _m._SI.EH.printf("%9.6g ",_MS.property_sample[k]);
		_m._SI.EH.print("%n");
		_m._SI.EH.printf("Genetic method: population variance %9.6g, %d function evaluations%n",
						DE.getPopulationCostVariance(), DE.getNumberOfFunctionEvaluations());
	  }
   
	return f_min; // FIXME make this f_best
  }

  public boolean abortNow(){return false;}
  public boolean wasAborted(){return false;}
  public void setFlags(boolean setOrReset){}

  /**
	 Method to get feedback from the solver...
  */
  public void updateNewPoint(double f, DoubleMatrix1D x, int itCount, int  nfeval)
  {
	if (false && itCount%1==0)   /* display after every refresh=10 generations */
	  {
		_m._SI.EH.println("%n%n%n Best-so-far cost funct. value= %-15.10g, num. func. evaluations=%d%n",
						  f, nfeval);			
		for (int j=0;j<x.size();j++)
		  {
			_m._SI.EH.printf("best[%d]=%-15.10g%n",j, x.getQuick(j));
		  }
	  }
	// take a copy
	f_min=f;
	times.assign(x);
  }

  /**
	 Make initial population. Must have c.rows()>=NP
  */
  public void create_initial_generation(int NP, DoubleMatrix2D c)
  {
	int i,ii,j,jj,Npop,dim;
	double mean;

	double[] inibound_l = new double[n_sub];
	double[] inibound_h = new double[n_sub];

	// layer time uncertainties estimated from prior+isopach
	DoubleMatrix1D t_uncert = _MS.t_layer_uncertainty.viewRow(_MS.fs);

	// set array bounds for each time
	// try mean +- three std deviations......
	for (j=0,jj=1; j<n_sub; j++,jj++)
	  {
		mean = _MS.mu.getQuick(jj);
		inibound_l[j] = mean - 3.0 * t_uncert.getQuick(jj);
		inibound_h[j] = mean + 3.0 * t_uncert.getQuick(jj);
	  }

	// Generate population....
	i=0; // running counter over population...


	// Case A, small dimension: try an experimental-design like spanning set
	Npop = (1<<n_sub);
	dim = c.columns();
	if (Npop <= NP)
	  {
		if (true)
		  {
			float[][] hypPts;
			hypPts = au.csiro.JamesGunning.lib.Numerical.Hypercube.createHypercubeList(dim,2);
			for (ii=0; ii<hypPts.length; i++,ii++) /* spread initial population members */
			  {
				for (j=0;j<dim; j++) /* spread initial population members */
				  {
					c.setQuick(i,j, inibound_l[j] + hypPts[ii][j]*(inibound_h[j]-inibound_l[j]));
				  }
			  }
			hypPts=null; // GC
			// rest of population can be random
			for (; i<NP; i++)
			  {
				for (j=0;j<dim; j++) /* spread initial population members */
				  {
					c.setQuick(i,j, inibound_l[j] + _MS.Ugen.nextDouble()*(inibound_h[j] - inibound_l[j]));
				  }
			  }
		  }
		else
		  {
			for (j=0,jj=0; j<dim; j++) /* spread initial population members */
			  {
				//			for (i=0; i<(1<<n_sub) ; i++,jj++)  
				for (; i<Npop ; i++,jj++)  
				  // loop over initial population, values alternating hi:lo like a Box 2-level 
				  // experimental design. Values sit at all corners in a d-dimensional hypercube
				  {
					c.setQuick(i,j, ((jj & (1<<j))>0) ? inibound_h[j] : inibound_l[j]);
					// _m._SI.EH.println("initial value: "+ (((jj & (1<<j))>0) ? 1 : -1));
				  }
				// rest of population can be random
				for (; i<NP; i++)
				  {
					c.setQuick(i,j, inibound_l[j] + _MS.Ugen.nextDouble()*(inibound_h[j] - inibound_l[j]));
				  }
			  }
		  }
	  }
	else
	  {
		// All initial population is random
		// for (i=0; i<NP; i++)
		for (; i<NP; i++)
		  {
			for (j=0; j<c.columns(); j++) /* spread initial population members */
			  {
				c.setQuick(i,j, inibound_l[j] + _MS.Ugen.nextDouble()*(inibound_h[j] - inibound_l[j]));
			  }
		  }
	  }
	if (_m._SI.verbosity>5)	_m._SI.EH.println("Finished creating genetic population");
	// What do these look like?
	if (!true)
	  {
		for (i=0; i<NP; i++)
		  {
			_m._SI.EH.print("GEP: ");
			for (j=0; j<c.columns(); j++)
			  {
				_m._SI.EH.printf("%10.8g ",c.getQuick(i,j));
			  }
			_m._SI.EH.print("%n");
		  }
	  }
  }

  /**
	 Here parameters[] is the union of the times and non-times properties we wish to
	 optimise over. The times are the components parameters[1....n_sub], and the
	 non-time variables are parameters[n_sub+1...n_tot]. The time component is the 
	 set of independent `cluster times' used in smart pinchout detection.
  */

  public double f_to_minimize(DoubleMatrix1D parameters)
  {
	double negLogLik;

	// (1) Constraint on times to fit into total reservoir interval
	testAndImposeTimeTruncation(parameters);

	// map the sub-vector parameters[0...n_sub-1] into required fields (internal times)
	mapSubTimesVectorToCompleteTimes(parameters);

	// order the times
	_MS.orderTimes();

	// For fluid contacts - time shifting may alter internal fluids and reflections....
	if (_m._SI.doSimpleFluidUnitsAndContacts)
	  {
		// time perturbation will affect depths and => internal layer FC's
		_m.computeEffectivePropertiesForAllLayers(false,false,false,false,false);
	  }


	// (1) Acquire logLikelihood contribution from seismic mismatch
	negLogLik =  _m.SM.computeSyntheticAndSeismicError();
	negLogSeisLikelihood = negLogLik;

	// (2) Acquire logLikelihood contribution from prior and isopachs

	// USE EXACT ISOPACH LIKELIHOOD WITH PRIOR
	// Acquire logLikelihood contribution from isopachs if we're not ignoring them
	// and there are some isopachs to attend to
	negLogIsopachLikelihood = _MS.computeIsopachNegLogLikelihood(false);
	negLogLik += negLogIsopachLikelihood;

	// (3) Contribution from depth self-consistency
	negLogDepthLikelihood = _MS.computeConsistentDepthNegLogLikelihood(false); 
	negLogLik += negLogDepthLikelihood;

	// (3) Contribution from vertical stacking
	if (_m._SI.doStackedModel)
	  {
		double negLogVStackLikelihood = _m.stackedModelInfo.computeVerticalStackingNegLogLikelihood(false); 
		negLogLik += negLogVStackLikelihood; 
	  }

	//  (5) Also the contribution from the raw prior 
	// (i.e. before isopachs were used in the linearised update)
	// Acquire logLikelihood contribution from prior.
	// The times are the variables indexed 0...(_m.n_active_layers-1)
	// top and base are fixed at their modes. 	

	// Compute prior contribution. There may be efficiencies to be gained here (only a 
	// subvector varies, and this component is almost certainly constant), 
	// but for now, be slow but safe.
	negLogPriorLikelihood = _MS.computePriorNegLogLikelihood();
	negLogLik += negLogPriorLikelihood;

	// Add truncation error too
	negLogLik+= negLogTruncationLikelihood;

	return negLogLik;
  }

  /** 
	  Tests whether optimiser has times at unacceptable values.
	  All times must lie in the interval <br>
	  ( _m.MPD.MetaInfo[0].t.value < t <  _m.MPD.t_base.value ) <br>
	  I.e. prior mean value of model-region top, and base
  */
  protected boolean testAndImposeTimeTruncation(DoubleMatrix1D parameters)
  {
	int i;
	double t_i=0.0,t_iminus1=0.0;
	negLogTruncationLikelihood=0.0;

	boolean truncated = false;
	for (i=0;i<n_sub;i++) 
	  {
		t_i = parameters.getQuick(i);
		if (t_i<=_m.MPD.MetaInfo[0].t.value) 
		  {
			// form penalty
			negLogTruncationLikelihood += MathUtil.sqr(_m.MPD.MetaInfo[0].t.value - t_i);
			// truncate??
			parameters.setQuick(i,_m.MPD.MetaInfo[0].t.value);	
			truncated = true;
			if (_m._SI.verbosity>4) _m._SI.EH.println("Upper T truncation imposed!");
		  }
		if (t_i>=_m.MPD.t_base.value)
		  {
			// form penalty
			negLogTruncationLikelihood += MathUtil.sqr(_m.MPD.t_base.value - t_i);
			// truncate??
			parameters.setQuick(i,_m.MPD.t_base.value); 	
			truncated = true;
			if (_m._SI.verbosity>4) _m._SI.EH.println("Bottom T truncation imposed!");
		  }
	  }
	return truncated;
  }

  /* 
	 map the times parameters[0....nsub-1] 
	 back to true internal layer times in the full model parameter vector
   */

  protected void mapSubTimesVectorToCompleteTimes(DoubleMatrix1D parameters)
  {
	int i,j;
	double t_i;

	for (i=0,j=1; i<n_sub; i++,j++)
	  {
		// direct times are parameters for internal layers...
		t_i = parameters.getQuick(i);
		// and write these values into the in-situ fields
		_MS.setValuesFromIndex(j,t_i);
		// update current snapshot too
		_MS.property_sample[j] = t_i;
	  }
  }
}
