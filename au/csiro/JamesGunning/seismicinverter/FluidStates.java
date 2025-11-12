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
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;

import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import java.util.ArrayList;
import static java.lang.Math.*;

public final class FluidStates
{
  /**
	 Fluids are assembled into a set of fluid combinations. The prior probability of
	 a given fluid combination is proportional to the product of the prior probabilities
	 of each fluid in each layer, and normalised so the sum over all possible fluid
	 combinations is unity.<br>

	 Some combinations may end up with negligible prior probability, so these are eliminated
	 from the list of interesting prior combinations. This variable is the threshold for
	 eliminating combinations from the list.
	 Combos with normalised probability less than this are eliminated
  */
  public static final double PriorFluidCombinationThreshold = 0.01;

  /**
	 Defines for labelling the varieties of possible density ordering:<br>
	 noDensityOrdering = any fluids allowed in any permeable layer<br>
  */
  public static final int noDensityOrdering = 0;
   /**
	 Defines for labelling the varieties of possible density ordering:<br>
	 partialDensityOrdering = fluids are density ordered for all adjacent permeable layers not separated by a pure shale<br>
  */
  public static final int partialDensityOrdering = 1;
  /**
	 Defines for labelling the varieties of possible density ordering:<br>
	 fullDensityOrdering = fluids are density ordered across all the reservoir model, regardless of whether there are
	 impermeable layers in the middle<br>
  */
  public static final int fullDensityOrdering = 2;
  
  /** 
	  array of fluid states for current fluid configuration
   */



  short[] fluid_type;
  double[] prob_state;
  double[] P_brine; 

  ArrayList<FluidState> possibleFluidStates;

  /** 
	  reference to underlying model
  */
  Model _m;

  /**
	 convenient reference to most likely fluid state index
  */
  int fs_max=0;

  /**
	 Cache for rechecking fluid state probabilities if targets used
	 rather than default method
  */
  ArrayList<IntShortFloatTriple> fs_probs_cache=null; 
  double[] fs_probs_cache_values=null;

  FluidStates(Model m)
  {
	this._m = m;
	fluid_type = new short[_m.n_layers];
	prob_state = new double[_m.n_layers];
	P_brine = new double[_m.n_layers];
  }

  protected void  setUpFluidStates()
  {
	// Map the number of possible fluid combinations
	
	ModelPriorDescriptions.MetaInfo MI;
	possibleFluidStates = new ArrayList<FluidState>();

	// St up matrix of possible fluids
	boolean C[][] = new boolean[_m.n_layers][FluidTypes.Nfluids]; //layers x {brine,oil,lsg,gas}
	for (int i=0;i<_m.n_layers;i++)
	  {
		MI = _m.MPD.MetaInfo[i];
		P_brine[i]=1.0;
		C[i][0]=true; // always possible to have brine
		if (!_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			// usual explicit accounting
			if (MI.oil_substitution_possible && MI.P_oil.value > 0) 
			  {
				C[i][1]=true; // might be oil
				P_brine[i]=Math.max(P_brine[i]-MI.P_oil.value,0.0); // update brine prob
			  }
			if (MI.lsg_substitution_possible && MI.P_lsg.value > 0) 
			  {
				C[i][2]=true; // might be lsg
				P_brine[i]=Math.max(P_brine[i]-MI.P_lsg.value,0.0); // update brine prob
			  }
			if (MI.gas_substitution_possible && MI.P_gas.value > 0) 
			  {
				C[i][3]=true; // might be gas
				P_brine[i]=Math.max(P_brine[i]-MI.P_gas.value,0.0); // update brine prob
			  }
		  }
		else
		  {
			// Under _m._SI.doSimpleFluidUnitsAndContacts, permit only one fluid mode (all "brine", as a 
			// placeholder for what the contacts prescribe)
			// and augment model as required later....
		  }
	  }
	// enumerate fluid combos as seen from layer "-1", with fictitious fluid 3 (gas) in it
	enumerateFluidCombos(-1,(short)(FluidTypes.Nfluids-1),C);

	if (possibleFluidStates==null || possibleFluidStates.size()==0)
	  {
		throw new RuntimeException("Density ordering mode and current sequence of fluid probabilities "
								   +"do not allow any legal fluid states. Rethink density ordering or "
								   +"your fluid probabilities!");
	  }
	normaliseFluidProbabilities();

	// if (_m._SI.verbosity>2) printFluidCombos();

	// apply --IF flags....
	trimToSubsetOfPossibleCombinationsIfRequired();
	
	if (_m._SI.useFluidProbabilitiesAsPostOrdering) 
	  {
		recalibrateStateProbabilitiesToOverallTargets();
	  }
	// and thresholding on priors
	trimPossibleCombinationsFromThresholds(true);
	
	possibleFluidStates.trimToSize();

	createFluidStateStrings();

	if (_m._SI.verbosity>2) printFluidCombos();


  }
  

  
  private void enumerateFluidCombos(int l,short j,boolean C[][])
  {
	int end_fluid_type,next_layer=0,k;
	double prob=0.0;
	if (l==_m.n_layers-1)
	  {
		// We have hit the bottom layer
		// product of layer probabilities for this fluid state
		prob=1.0;
		for (k=0;k<_m.n_layers;k++)
		  {
			prob*=prob_state[k];
		  }
		// add complete fluid config and prior prob to list.
		possibleFluidStates.add(new FluidState(fluid_type,prob));
		return;
	  }
	else
	  {
		if (_m._SI.density_ordering_mode == fullDensityOrdering)	  
		  {
			// here we skip impermeable barriers so fluids 'see' across them
			// find the next permeable layer
			for (k=l+1;k<_m.n_layers;k++)
			  {
				if (!_m.MPD.MetaInfo[k].is_permeable)
				  {
					// a shale
					fluid_type[k]=FluidTypes.brine; // has to be brine
					prob_state[k]=P_brine[k];
				  }
				else
				  {
					next_layer=k;
					break;
				  }
			  }
			if (k==_m.n_layers)  // we had a shale sequence and ran out
			  {
				enumerateFluidCombos(_m.n_layers-1,FluidTypes.brine,C);
				return;
			  }
		  }
		else
		  {
			next_layer=l+1;
		  }

		// To completely ignore density, the end_fluid type is always FluidTypes.Nfluids-1, but
		// if we want ordering, it must be chosen carefully
		if (_m._SI.density_ordering_mode != noDensityOrdering)
		  {
			if (l<0 || !(_m.MPD.MetaInfo[l].is_permeable))// is a shale above??
			  {
				end_fluid_type = FluidTypes.Nfluids-1;
			  }
			else
			  {
				end_fluid_type = j; // layer above is permeable => density ordering has to apply
			  }
		  }
		else
		  {
			// noDensityOrdering. Pretend layer above is impermeable, so this can have anything.....
			end_fluid_type = FluidTypes.Nfluids-1;
		  }
		

		if (true)
		  {
			for (fluid_type[next_layer]=0; fluid_type[next_layer]<=FluidTypes.Nfluids-1; fluid_type[next_layer]++) 
			  {
				if (FluidTypes.density_priority[fluid_type[next_layer]] > FluidTypes.density_priority[end_fluid_type])
				  {
					// fluid_type[next_layer] is lighter, so skip...
					continue;
				  }
				if (C[next_layer][fluid_type[next_layer]]) 
				  {
					// this fluid is possible at layer next_layer
					switch (fluid_type[next_layer])
					  {
					  case FluidTypes.brine: 
						prob = P_brine[next_layer];break;
					  case FluidTypes.oil: 
						prob = _m.MPD.MetaInfo[next_layer].P_oil.value;break;
					  case FluidTypes.lsg:
						prob = _m.MPD.MetaInfo[next_layer].P_lsg.value;break;
					  case FluidTypes.gas:
						prob = _m.MPD.MetaInfo[next_layer].P_gas.value;break;
					  }
					prob_state[next_layer]=prob;
					if (prob>0) enumerateFluidCombos(next_layer,fluid_type[next_layer],C);
				  }
			  }
		  }
		else
		  {
			// Old version, with br,oli,lsg,gas ordering
			for (fluid_type[next_layer]=0; fluid_type[next_layer]<=end_fluid_type; fluid_type[next_layer]++) 
			  {
				if (C[next_layer][fluid_type[next_layer]]) 
				  {
					// this fluid is possible at layer next_layer
					switch (fluid_type[next_layer])
					  {
					  case FluidTypes.brine: 
						prob = P_brine[next_layer];break;
					  case FluidTypes.oil: 
						prob = _m.MPD.MetaInfo[next_layer].P_oil.value;break;
					  case FluidTypes.lsg:
						prob = _m.MPD.MetaInfo[next_layer].P_lsg.value;break;
					  case FluidTypes.gas:
						prob = _m.MPD.MetaInfo[next_layer].P_gas.value;break;
					  }
					prob_state[next_layer]=prob;
					if (prob>0) enumerateFluidCombos(next_layer,fluid_type[next_layer],C);
				  }
			  }
		  }
		return;
	  }
  }

  
  protected void normaliseFluidProbabilities()
  {
	int i;
	double norm=0.0;
	FluidState fs;
	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		fs = possibleFluidStates.get(i);
		norm += fs.prob;
	  }
	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		fs = possibleFluidStates.get(i);
		fs.prob /= norm;
		fs.log_prob = ( fs.prob > 0 ) ? Math.log(fs.prob) : -Double.MAX_VALUE;
	  }
  }

  protected void trimPossibleCombinationsFromThresholds(boolean force)
  {
	int i,count=0;
	double thresh,probMax=0.0;

	// double pass to make sure this is sensible
	// get highest fluid state prior probability (probMax)
	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		FluidState fs = possibleFluidStates.get(i);
		probMax = Math.max(probMax,fs.prob);
	  }
	// use elimination threshold relative to most probable prior model
	if (!_m._SI.useFluidProbabilitiesAsPostOrdering || force)
	  {
		thresh = PriorFluidCombinationThreshold*probMax;
	  }
	else thresh=1e-5;
	// now drop prior states that look too relatively improbable
	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		FluidState fs = possibleFluidStates.get(i);
		if (fs.prob < thresh)
		  {
			possibleFluidStates.remove(i); // causes left shift in successors
			i--;
		  }
	  }
	if (possibleFluidStates.size() > 100)
	  {
		_m._SI.EH.println("WARNING: You probably have too many states in a flat prior. Strongly suggest simplifying the model");
	  }
  }

 protected void recalibrateStateProbabilitiesToOverallTargets()
  {
	int i,j,i_f,pass,Np,Nfs,Nfsm1;
	short k;
	double prob=0,s,sp,f,fn,df,dx,a,eNH,EPS;
	FluidState fs;
	ModelPriorDescriptions.MetaInfo MI;
	DenseDoubleMatrix2D XF,XR,H,rhs,ans;


	if (_m._SI.doSimpleFluidUnitsAndContacts)
	  {
		_m._SI.EH.println("WARNING: recalibration of fluid probabilities not possible under contact model!");
		return;
	  }
	else if (_m._SI.verbosity>2) 
	  {
		_m._SI.EH.println("%nRecalibrating fluid-state probabilities to overall targets");
	  }

	Nfs = possibleFluidStates.size();

	// count dimension of target probabilities Np & assemble equations	
	ArrayList<IntShortFloatTriple> fs_probs = new ArrayList<IntShortFloatTriple>();
	for(j=0;j<_m.n_layers;j++)
	  {
		MI = _m.MPD.MetaInfo[j];
		for(k=0;k<FluidTypes.Nfluids;k++)
		  {
			prob=0;
			switch (k)
			  {
			  case FluidTypes.brine: 
				continue;
			  case FluidTypes.oil: 
				if (MI.oil_substitution_possible) prob = MI.P_oil.value; 
				break;
			  case FluidTypes.lsg:				
				if (MI.lsg_substitution_possible) prob = MI.P_lsg.value; 
				break; 
			  case FluidTypes.gas:
				if (MI.gas_substitution_possible) prob = MI.P_gas.value; 
				break; 
			  }
			if (prob>0) fs_probs.add(new IntShortFloatTriple(j,k,prob));
		  }
	  }
	fs_probs.trimToSize();
	Np = fs_probs.size();

	// check against cache to see if any of the calculations need be done 
	// (i.e. is trace same as last time)
	boolean doCompute=false;
	if (fs_probs_cache != null)
	  {
		// compare to existing cache

		if (fs_probs_cache.size() == Np)
		  {
			IntShortFloatTriple triple,triple_c;
			for(i=0;i<Np;i++)
			  {
				triple = fs_probs.get(i);
				triple_c = fs_probs_cache.get(i);
				if (triple.A != triple_c.A || triple.B != triple_c.B || triple.p !=  triple_c.p)
				  {
					doCompute = true;
					break;	
				  }
			  }
		  } 
		else 
		  {
			doCompute=true;  // have to resolve
			if (fs_probs_cache_values.length != Nfs)
			  fs_probs_cache_values = new double[Nfs];
			fs_probs_cache = fs_probs;
		  }
	  }
	else
	  {
		// alloc cache
		fs_probs_cache  = fs_probs;
		fs_probs_cache_values = new double[Nfs];
		doCompute=true;
	  }
	if (!doCompute)
	  {
		// can reuse cache, avoid solving system
		for(i=0;i<Nfs;i++) 
		  {
			fs = possibleFluidStates.get(i);
			fs.prob = fs_probs_cache_values[i];
		  }
		return;  // early cheap exit
	  }
	
	// Assemble polytope-bounded OLS system.
	XF = new DenseDoubleMatrix2D(Np,Nfs);
	for (i_f=0 ;i_f<Nfs;i_f++)
	  {
		fs = possibleFluidStates.get(i_f);
		for(i=0;i<Np;i++) // loop over target probabilities
		  {
			IntShortFloatTriple triple = fs_probs.get(i);
			j = triple.A; // layer index
			k = triple.B; // fluid type
			if (fs.types[j] == k) XF.s(i,i_f,1.0); // fill Xmatrix
		  }
	  }
	//_m._SI.EH.println("Matrix X:\n"+XF);

	// Create reduced OLS problem incorporating sum-constraint on fluid-state probabilities
	// i.e. variables are state probabilities p_1...p_{Nfs-1}, last one automatically filled by sum constraint
	DenseDoubleMatrix1D x,x_new,g,P;
	Nfsm1 = Nfs-1;
	XR = new DenseDoubleMatrix2D(Np,Nfsm1);  // reduced design matrix incorporating sum constraint
	for(i=0;i<Np;i++) {
	  for(j=0;j<Nfsm1;j++) XR.s(i,j,XF.g(i,j)-XF.g(i,Nfsm1));
	}
	// Form reduced data PR
	P = new DenseDoubleMatrix1D(Np);
	// stuff target probabilities into vector & adjust for reduced system
	for (i=0;i<Np;i++) P.s(i,fs_probs.get(i).p - XF.g(i,Nfsm1));  

	// (1) OLS solution via Cholesky. Try first, as fastest. Not polytope-bounded
	DenseDoubleMatrix2D Pvars = new DenseDoubleMatrix2D(Np,1);
	H = new DenseDoubleMatrix2D(Nfsm1,Nfsm1);
	for (i=0;i<Np;i++) Pvars.s(i,0,P.g(i));  // stuff target probabilities into vector
	rhs = new DenseDoubleMatrix2D(Nfsm1,1);
	SeqBlas.seqBlas.dgemm(true,false,1.0,XR,XR,0,H);           // form X^T.X = H
	// stabilise, since sometimes XR is short-rank
	EPS=1e-8;
	for(i=0;i<Nfsm1;i++) H.s(i,i,H.g(i,i)+EPS); 
 	CholeskyDecomposition chol = new CholeskyDecomposition(H); 

	if (chol.isSymmetricPositiveDefinite())
	  {
		SeqBlas.seqBlas.dgemm(true,false,1.0,XR,Pvars,0,rhs);  // form X^T.P = rhs
		ans = (DenseDoubleMatrix2D)chol.solve(rhs);
		doCompute = false;
		EPS = 1.0e-5;  // epsilon for fp error
		for(i=0,s=0;i<Nfsm1;i++) 
		  {
			prob = ans.g(i,0);
			if (prob<-EPS || prob>(1+EPS)) {doCompute=true; break;} // needs to be bounded
			s+=prob;
		  }
		if (s>1.0) doCompute=true;
		if (!doCompute)
		  {
			// Cholesky solution got through. Set it and return early
			for(i=0;i<(Nfsm1);i++) 
			  {
				fs = possibleFluidStates.get(i);
				fs.prob = min(max(ans.g(i,0),0),1);
				fs_probs_cache_values[i] = fs.prob;
			  }
			// last state from sum-constraint. 
			fs = possibleFluidStates.get(Nfsm1);
			fs.prob = max(1.0-s,0);
			fs_probs_cache_values[Nfsm1] = fs.prob;	
			if (_m._SI.verbosity>2) 
			  _m._SI.EH.println("Satisfactory recalibrated state probabilities from unbounded Cholesky");
			return;
		  }
	  }
	
	// (2) Fallback Crude steepest descent polytope-bound optimization for new fluid-state probabilities
	class OLScompute
	{
	  DenseDoubleMatrix2D X;
	  DenseDoubleMatrix1D P,resid;
	  double EPS=1e-14;
	  public OLScompute(DenseDoubleMatrix2D X,DenseDoubleMatrix1D P)
	  {
		this.X=X; 
		this.P=P;
		resid = new DenseDoubleMatrix1D(X.rows());
	  }
	  // Computes, given point m, f=(1/2) ||X.m-P||^2. And, if computeGrad==true, g = X^T(X.m-rhs)
	  // to stablise against possibly rank-short X, add:
	  // f=(1/2) (||X.m-P||^2 + EPS||m||^2), so g = X^T(X.m-rhs) + EPS*m
	  public double f(DenseDoubleMatrix1D m,boolean computeGrad,DenseDoubleMatrix1D g)
	  {
		int j,Np;
		double s,sum,f;
		SeqBlas.seqBlas.dgemv(false,1.0,X,m,0,resid);  // resid <- X.m
		Np=resid.size();
		for(j=0,f=0; j<Np; j++) 
		  {
			s = resid.g(j)-P.g(j);   // X.m - P
			resid.s(j,s);   // resid = X.m-rhs
			f += s*s;
		  }
		for(j=0,sum=0;j<m.size();j++) 
		  {
			s = m.g(j);
			sum+=s*s;
		  }
		f += EPS * sum;
		if (computeGrad) 
		  {
			SeqBlas.seqBlas.dgemv(true,1.0,X,resid,0,g);      // g = X^T(X.m-rhs)
			for(j=0;j<m.size();j++) g.s(j,g.g(j)+EPS*m.g(j)); //    + EPS.m
		  }
		return 0.5*f;
	  }
	}
		
	OLScompute OLS = new OLScompute(XR,P);

	// Working space for SD method
	x = new DenseDoubleMatrix1D(Nfsm1);
	x_new = new DenseDoubleMatrix1D(Nfsm1);
	g = new DenseDoubleMatrix1D(Nfsm1);
	double EPS_TERM;

	// safe starting guess for x
	x.assign(1.0/Nfs);
	int it,MAX_ITS=500;
	EPS = 1.0e-5;     // backtracking epsilon
	EPS_TERM = 5e-5;  // termination ||dx|| epsilon
	eNH = (Nfsm1>0) ? 1.0/sqrt((Nfsm1)): 1.0; 
	double[] t = new double[Nfsm1];
	boolean terminate=false;
	for(df=dx=f=fn=0,it=0;it<MAX_ITS;it++)
	  {
		f = OLS.f(x,true,g);  // compute func and gradient
		for(a=1.0; a>EPS; a*=0.5)
		  {
			for(j=0,sp=0;j<Nfsm1;j++) 
			  {
				s = x.g(j) - a * g.g(j); // steepest descent
				// x_new = new point from x - a.g
				x_new.s(j,s);
				sp+=s;
			  }
			// test if point fails inequality sum_i x(i) >=1
			if (sp>1)
			  {
				double a0,sX,sum;
				// OK have to project along hyperplane sum_i x_i=1. 
				// Outward normal to this plane is n_hat={eNH,eNH,...eNH}
				// s = (x'-x).n_hat = -g.n_hat, (with a=1.0)
				for(j=0,s=0;j<Nfsm1;j++) s -= g.g(j)*eNH; 
				for(j=0,sum=0;j<Nfsm1;j++) 
				  {
					t[j] = -g.g(j) - s * eNH;  // plane tangent vector t = (xn-x) -((xn-x).n_hat)n_hat
					sum -= g.g(j);
				  }
				//  value of a at which descent vector intersects inequality plane
				for(j=0,sX=1.0; j<Nfsm1; j++) sX -= x.g(j);
				a0 = (abs(sum)>1e-12) ? (sX/sum) : 0;  // avoid small division
				// project onto constraint plane when a > a0					  
				// x(a) = x + a0(xn-x)+(a-a0)t
				for(j=0;j<Nfsm1;j++) x_new.s(j, x.g(j) - a0*g.g(j)+(a-a0)*t[j]);
				// a check
				//for(j=0,s=0;j<Nfsm1;j++) s+= x_new.g(j);
			 	//_m._SI.EH.println("ONPLANE: "+s);
			  }

			for(j=0;j<Nfsm1;j++) x_new.s(j,min(max(x_new.g(j),0),1.0));
			fn = OLS.f(x_new,false,null);
			if (_m._SI.verbosity>3) 
			  {
				_m._SI.EH.println("recalibration fluid state optimiser points: "
								  +"a=%-15.8g  fO=%-15.8g f=%-15.8g",a,f,fn);
			  }
			// distance to new point
			for(j=0,dx=0;j<Nfsm1;j++) dx += MathUtil.sqr(x_new.g(j)-x.g(j));
			df = fn-f;
			if (df<(-1e-4*dx/a)) // like 5.13, Kelley, sufficient decrease
			  {
				f=fn;
				x.assign(x_new);
				break;
			  }
		  }

		if (sqrt(dx)<EPS_TERM) break;  // Same as Kelley 5.18, since r_0 = O(1)
	  }

	// set new probabilities
	for(i=0,s=0;i<(Nfsm1);i++) 
	  {
		fs = possibleFluidStates.get(i);
		fs.prob = x.g(i);
		s += fs.prob;
		// save in cache
		fs_probs_cache_values[i] = fs.prob;
	  }
	// last state filled in from sum-constraint. Thresholding for fp error
	fs = possibleFluidStates.get(Nfsm1);
	fs.prob = Math.max(0,1.0-s); 
	fs_probs_cache_values[Nfsm1] = fs.prob;	// save in cache 
	if (_m._SI.verbosity>3) 
	  {
		for(i=0,s=0;i<Nfs;i++) 
		  {
			fs = possibleFluidStates.get(i);
			_m._SI.EH.printf("SD fluid state solution fs=%d prob=%8.4g%n",i,fs.prob);
		  }
	  }
	if (false)
	  {
		// test profiles
		int Npoints=10;
		double delta=0.01;
		double[] fvals=new double[2*Npoints+1];
		double[] xvals=(double[])fvals.clone();
		for(i=0;i<Nfsm1;i++)
		  {
			x_new.assign(x);
			for(j=-Npoints,k=0;j<=Npoints;j++,k++) 
			  {
				xvals[k] = x.g(i)+j*delta;
				x_new.s(i,xvals[k]);
				fvals[k] = OLS.f(x_new,false,null);
			  }
			Util.printfileZO("ProbXsection."+i+".txt",xvals,fvals,fvals.length);
		  }
	  }
  }


 /**
	 Certain options require the confining of attention to only
	 one fluid state. This function will strip all other
	 states out of the list of fluid states
  */
  protected void trimToSubsetOfPossibleCombinationsIfRequired()
  {
	int i;

	// for testing purposes we may want to confine ourselves to one fluid state
	// without altering the prior model. Do this via the global flag
	if (_m._SI.ignoreAllButFirstFluidState)
	  {
		// remove trailing ones 
		if (_m._SI.restrictedFluidStateIndex<0)
		  {
			throw new RuntimeException("Desired restricted fluid state (--IF "+
									   _m._SI.restrictedFluidStateIndex+") is less than 0");
		  }
		if (_m._SI.restrictedFluidStateIndex >= possibleFluidStates.size())
		  {
			throw new RuntimeException("Desired restricted fluid state (--IF "+
									   _m._SI.restrictedFluidStateIndex+") is too large");
		  }
		int fullSize = possibleFluidStates.size();
		// remove trailing ones 
		for (i=_m._SI.restrictedFluidStateIndex+1; i<fullSize; i++)
		  {
			possibleFluidStates.remove(_m._SI.restrictedFluidStateIndex+1); // causes left shift in successors
		  }
		for (i=0; i<_m._SI.restrictedFluidStateIndex; i++)
		  {
			possibleFluidStates.remove(0); // causes left shift in successors
		  }
	  }

	possibleFluidStates.trimToSize();
  }

  

  // Create useful strings to represent state
  void createFluidStateStrings()
  {	
	int i,j;
	StringBuilder fsstring;

	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		fsstring = new StringBuilder(_m.n_layers);
		FluidState fs = possibleFluidStates.get(i);
		if (!_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			for (j=0 ;j<_m.n_layers;j++) fsstring.append(FluidTypes.getFluidChar(fs.types[j]));
		  }
		else
		  {
			for (j=0 ;j<_m.n_layers;j++) fsstring.append((_m.MPD.MetaInfo[j].is_permeable)? 'C' : 'b');
		  }
		fs.fState = fsstring.toString();
		fsstring = null;
	  }
  }

  protected void printFluidCombos()
  {
	int i,j;
	char fl='b';
	
	_m._SI.EH.printf("%nFluid states %s:%n", (_m._SI.doSimpleFluidUnitsAndContacts) ? "(C = contact-determined)" : "");

	for (i=0 ;i<possibleFluidStates.size();i++)
	  {
		FluidState fs = possibleFluidStates.get(i);
		_m._SI.EH.printf("Fluid state %2d: prior prob = %1.4f :   ",i,fs.prob);
		if (!_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			for (j=0 ;j<_m.n_layers;j++)
			  {
				fl = FluidTypes.getFluidChar(fs.types[j]);
				_m._SI.EH.print(fl + " ");
			  }
		  }
		else
		  {
			// more informative output for contacts mode
			for (j=0 ;j<_m.n_layers;j++)
			  {
				fl = (_m.MPD.MetaInfo[j].is_permeable)? 'C' : 'b';
				_m._SI.EH.print(fl + " ");
			  }
		  }
		_m._SI.EH.println("");
	  }
  }

  /**
	 simple check to see if prior fluid probs are sane for a given layer i
  */

  void checkPriorFluidProbabilities(int i)
  {
	if (_m.MPD.MetaInfo[i].is_permeable)
	  {
		double P_hydro=0;
		if (_m.MPD.MetaInfo[i].oil_name != null)
		  P_hydro += _m.MPD.MetaInfo[i].P_oil.value;
		if (_m.MPD.MetaInfo[i].gas_name != null)
		  P_hydro +=_m.MPD.MetaInfo[i].P_gas.value;
		if (_m.MPD.MetaInfo[i].lsg_name != null)
		  P_hydro +=_m.MPD.MetaInfo[i].P_lsg.value;

		if (P_hydro>1.0)
		  {
			_m._SI.EH.println("%nWARNING!!!!!: total hydrocarbon probability > 1%n");
		  }		  
	  }
  }
			
			
  /**
	 Each possible fluid configuration has an array of types in each
	 layer and the associated prior probability
  */

  public final class FluidState
  {
	short[] types;
	double prob;
	double log_prob; // efficiency calculation; to be kept equal to log(prob)
	String fState;  // String concatenation representation of types[0...m.n_layers-1]

	public FluidState(short[] types,double prob)
	{
	  types =(short[])types.clone();
	  this.types=types;
	  this.prob=prob;
	  this.log_prob = ( this.prob > 0 ) ? Math.log(this.prob) : -Double.MAX_VALUE;
	}
  }

  /**
	 enum type fluid - labelling class
  */
  public static final class FluidTypes
  {
	public static final short Nfluids=4;
	
	/**
	   Basic fluid codes
	*/
	public static final short brine=0;
	public static final short oil=1;
	public static final short lsg=2;
	public static final short gas=3;


	/**
	   This sequence defines the density ordering priority of the
	   fluid types above.
	   This can be changed to alter density ordering without
	   needing to alter the fluid codes.

	   Currently the new (8/Aug 07) default ordering is
	   {brine, lsg, oil, gas}
	   The idea is that lsg may be trapped by 
	   residual saturation mechanisms etc.
	*/
	public static final short[] density_priority = {0,2,1,3};

	static char getFluidChar(short type)
	{
	  char fl='\0';
	  switch (type)
		{
		case FluidTypes.brine: 
		  fl= 'b';break;
		case FluidTypes.oil: 
		  fl= 'o';break;
		case FluidTypes.lsg:
		  fl= 'l';break;
		case FluidTypes.gas:
		  fl= 'g';break;
		default:
		  fl='?'; break;
		}
	  return fl;
	}
  }

  /**
	 Utility class for fluid state recomputation from target probabilities
  */
  public final class IntShortFloatTriple 
  {
	int A; short B; double p;
	public IntShortFloatTriple(int A,short B,double p) 
	{
	  this.A=A; this.B=B; this.p=p;
	} 	   
  }
}
