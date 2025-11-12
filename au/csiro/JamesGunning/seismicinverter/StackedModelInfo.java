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
import cern.colt.matrix.*;
import au.csiro.JamesGunning.lib.util.MathUtil;

public final class StackedModelInfo
{
  /*
	pointer to underlying model
  */
  Model _m;
  /*
	shortcut to seismic inverter too
  */
  SeismicInverter _SI;

  /*
	Possible stacked model information
  */
  List<StackedPropertyMappingDetails> layer_maps = new ArrayList<StackedPropertyMappingDetails>();

  StackedPropertyMappingDetails[] l_maps; // array backing of above

  /**
	 Array of seismic exclusion zones: for each zone, the exclusion
	 zone starts 3*sigma_t+coda above layer
	 seismic_exclusion_start_layer top mean time, and ends
	 3*sigma_t+precusor+t below named layer
  */
  SeismicExclusionLayers[]  seismicExclusionLayers = null;



  StackedPropertyInfoDetails vp_m;
  StackedPropertyInfoDetails vp_s;
  StackedPropertyInfoDetails vs_m;
  StackedPropertyInfoDetails vs_s;
  StackedPropertyInfoDetails rho_m;
  StackedPropertyInfoDetails phi_s;
  StackedPropertyInfoDetails NG;
  StackedPropertyInfoDetails thickness;
  
  /* for later use */
  int n_constraints=0;

  double logDetCvertstackAnd2PiPieces=0.0;


  StackedModelInfo(Model m)
  {
	this._m = m;
	this._SI = m._SI;

	vp_m = new StackedPropertyInfoDetails(this);
	vp_s = new StackedPropertyInfoDetails(this);
	vs_m = new StackedPropertyInfoDetails(this);
	vs_s = new StackedPropertyInfoDetails(this);
	rho_m = new StackedPropertyInfoDetails(this);
	phi_s = new StackedPropertyInfoDetails(this);
	NG = new StackedPropertyInfoDetails(this);
	thickness = new StackedPropertyInfoDetails(this);
  }

  void setUpMapsArray()
  {
	l_maps = layer_maps.toArray(new StackedPropertyMappingDetails[layer_maps.size()]);
  }

  void setUpConstraints()
  {
	n_constraints =  computeVerticalStackingErrorsSize();
  }

  /** 
	  Count number of "effective data" points in extra constraints from stacking.
	  Keep internals synchronized (in terms of functionality) with 
	  computeVerticalStackingNegLogLikelihood()! 
	  May possibly have dependencies on fluid state, but I don't think so.

	  Updates also normalisation term logDetCvertstackAnd2PiPieces associated with this likelihood.

	  Last revision confines extra constraints *only* to thickness (nonlinear) term, as this cannot
	  be explicitly absorbed in the Gaussian prior.
  */
  private int computeVerticalStackingErrorsSize()
  {	
	int i;

	i=0;
	logDetCvertstackAnd2PiPieces=0.0;

	if (_SI.doStackedModel)
	  {	
		int i_b,j,jp,k;
		Integer I;
		double sigma=0.0;

		// set useful maps...
		_m.MS.map_to = _m.MS.map_to_index[0];



		// sum over target layer range
		for (i_b=0; i_b<l_maps.length; i_b++)
		  {
			StackedPropertyMappingDetails map = l_maps[i_b];
			j = map.base_layer-1; // indx of base layer

			LayerProperties LPI = _m.LP[j];

			for (k=0; k< map.target_layers.length; k++)
			  {
				jp = map.target_layers[k] - 1;  // index of comparison layer

				if (false)
				  {
					// All the linear properties below subsumed into prior
					
					// vp_m
					I = (Integer)(_m.MS.map_to.get(LPI.vp_m));
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vp_m)))!=null && vp_m.doConstraint && vp_m.isApplicable[j])  // the vp_m is uncertain
					  {
						i++;
						sigma = vp_m.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }   
					// vp_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vp_s)))!=null && vp_s.doConstraint && vp_s.isApplicable[j])  // the vp_s is uncertain
					  {
						i++;
						sigma = vp_s.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
					// vs_m
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vs_m)))!=null && vs_m.doConstraint && vs_m.isApplicable[j])  // the vs_m is uncertain
					  {
						i++;
						sigma = vs_m.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
					// vs_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vs_s)))!=null && vs_s.doConstraint && vs_s.isApplicable[j])  // the vs_s is uncertain
					  {
						i++;
						sigma = vs_s.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
					// rho_m
					if ((I = (Integer)(_m.MS.map_to.get(LPI.rho_m)))!=null && rho_m.doConstraint && rho_m.isApplicable[j])  // the rho_m is uncertain
					  {
						i++;
						sigma = rho_m.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
					// phi_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.phi_s)))!=null && phi_s.doConstraint && phi_s.isApplicable[j])  // the phi_s is uncertain
					  {
						i++;
						sigma = phi_s.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
					// NG
					if ((I = (Integer)(_m.MS.map_to.get(LPI.NG)))!=null && NG.doConstraint && NG.isApplicable[j])  // the NG is uncertain
					  {
						i++;
						sigma = NG.sigma;
						logDetCvertstackAnd2PiPieces += Math.log(sigma);
					  }
				  }

				// thickness
				if (thickness.doConstraint && thickness.isApplicable[j])  // the thickness is uncertain
				  {
					i++;
					sigma = thickness.sigma;
					logDetCvertstackAnd2PiPieces += Math.log(sigma);
				  }
			  }
		  }
		// every term contributes one of these
		logDetCvertstackAnd2PiPieces += (double)i * MathUtil.HALF_LOG_TWOPI;
	  }

	return i;
  }



  /** 
	  Likelihood function for imposing property "similarity". Computed
	  after times are truncated and ordered, so depths are non-negative.
	  LayerProperty effective properties must be current and valid.
  */

  double computeVerticalStackingNegLogLikelihood(boolean showErrors)
  {	
	double negLogLik=0.0;

	if (_SI.doStackedModel)
	  {	
		int i=0,i_b,j,jp,k,indx;
		double chisq=0.0, y1=0,y2=0,mean,sigma=0;
		Integer I;
		double[] vertical_stacking_errors = _m.MS.vertical_stacking_errors;

		_m.MS.logDetCvsAnd2PiPieces=0.0;

		// sum over target layer range
		for (i_b=0; i_b<l_maps.length; i_b++)
		  {
			StackedPropertyMappingDetails map = l_maps[i_b];
			j = map.base_layer-1; // indx of base layer

			LayerProperties LPI = _m.LP[j];

			for (k=0; k< map.target_layers.length; k++)
			  {
				jp = map.target_layers[k] - 1;  // index of comparison layer

				LayerProperties LPI_compare = _m.LP[jp];
			
				if (false)
				  {
					// All linear terms below now absorbed in prior !!!!

				  
					// vp_m
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vp_m)))!=null && vp_m.doConstraint && vp_m.isApplicable[j])  // the vp_m is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = vp_m.sigma;
						y1 = LPI.vp_m.value;
						y2 = LPI_compare.vp_m.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }   
					// vp_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vp_s)))!=null && vp_s.doConstraint && vp_s.isApplicable[j])  // the vp_s is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = vp_s.sigma;
						y1 = LPI.vp_s.value;
						y2 = LPI_compare.vp_s.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
					// vs_m
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vs_m)))!=null && vs_m.doConstraint && vs_m.isApplicable[j])  // the vs_m is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = vs_m.sigma;
						y1 = LPI.vs_m.value;
						y2 = LPI_compare.vs_m.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
					// vs_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.vs_s)))!=null && vs_s.doConstraint && vs_s.isApplicable[j])  // the vs_s is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = vs_s.sigma;
						y1 = LPI.vs_s.value;
						y2 = LPI_compare.vs_s.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
					// rho_m
					if ((I = (Integer)(_m.MS.map_to.get(LPI.rho_m)))!=null && rho_m.doConstraint && rho_m.isApplicable[j])  // the rho_m is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = rho_m.sigma;
						y1 = LPI.rho_m.value;
						y2 = LPI_compare.rho_m.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
					// phi_s
					if ((I = (Integer)(_m.MS.map_to.get(LPI.phi_s)))!=null && phi_s.doConstraint && phi_s.isApplicable[j])  // the phi_s is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = phi_s.sigma;
						y1 = LPI.phi_s.value;
						y2 = LPI_compare.phi_s.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
					// NG
					if ((I = (Integer)(_m.MS.map_to.get(LPI.NG)))!=null && NG.doConstraint && NG.isApplicable[j])  // the NG is uncertain
					  {
						//indx = I.intValue();
						//mean = mu.getQuick(indx);
						sigma = NG.sigma;
						y1 = LPI.NG.value;
						y2 = LPI_compare.NG.value;
						vertical_stacking_errors[i] = (y1-y2)/sigma;
						chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
						if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
						i++;
					  }
				  }

				// Nonlinear thickness term!

				// thickness
				if (thickness.doConstraint && thickness.isApplicable[j])  // the thickness is uncertain
				  {
					double t_base;
					
					// layer time thickness...
					t_base = (j==_m.n_layers-1) ? LPI.t_base.value : _m.LP[j+1].t.value;
					// and => DZ = thickness. Factor of 1/2 since times are two-way
					y1 = LPI.vp_eff * 0.5 * (t_base - LPI.t.value);
					// comparison layer....
					t_base = (jp==_m.n_layers-1) ? LPI_compare.t_base.value : _m.LP[jp+1].t.value;
					// thickness of comparison layer..
					y2 = LPI_compare.vp_eff * 0.5 * (t_base - LPI_compare.t.value);
					// => error
					sigma = thickness.sigma;
					vertical_stacking_errors[i] = (y1-y2)/sigma;
					chisq += vertical_stacking_errors[i] * vertical_stacking_errors[i];
					
					_m.MS.logDetCvsAnd2PiPieces += 0.5 * Math.log(sigma*sigma);
					_m.MS.logDetCvsAnd2PiPieces += MathUtil.HALF_LOG_TWOPI;


					if (showErrors) _SI.EH.printf("%10.4g ",vertical_stacking_errors[i]);
					i++;
				  }
			  }
		  }
		negLogLik += chisq/2.0;
	  }
	return negLogLik;
  }

  /**
	 data class for keeping list of exclusion-zone layer numbers
  */
  static class SeismicExclusionLayers
  {
	int seismic_exclusion_start_layer;
	int seismic_exclusion_end_layer;
  }

  static class StackedPropertyMappingDetails
  {
	/**
	   this layer is used as the basis for comparisons over delta_properties
	*/
	int base_layer;
	/**
	   and these are the "other layers" over whom a sum of "delta-property" squared type
	   errors contribute to the likelihood, is the property is stochastic and the constraint is applicable
	*/
	int[] target_layers;
	
	StackedPropertyMappingDetails() {}
  }
  
  static class StackedPropertyInfoDetails
  {
	StackedModelInfo _sMI;
	boolean doConstraint=false;
	double sigma = Double.NaN;
	boolean[] isApplicable;

	StackedPropertyInfoDetails(StackedModelInfo sMI)
	{
	  _sMI = sMI;
	  isApplicable = new boolean[sMI._m.n_layers];
	  setAllApplicable(true);  // switched on by default
	}

	void checkApplicableConstraints()
	{
	  if (!Double.isNaN(sigma)) doConstraint = true;
	}


	void setAllApplicable(boolean yes)
	{
	  Arrays.fill(isApplicable,yes);
	}
  }

}

