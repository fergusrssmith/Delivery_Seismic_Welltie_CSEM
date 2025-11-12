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
//********************************************************************************//
//                                                                                //
//    LICENCE FOR BHP SEISMIC INVERSION SUITE OF PROGRAMS                         //
//                                                                                //
//    The following is the license that applies to this copy of the BHP Seismic   //
//    Inversion Suite of Programs software (hereinafter the "Software"). BHP's    //
//    Software Manager may be contacted at the following address:                 //
//                                                                                //
//          CSIRO Petroleum                                                       //
//          P.O. Box 3000                                                         //
//          Glen Waverley, Victoria                                               //
//          Australia 3150                                                        //
//          Telephone Number:  613 9545 8396                                      //
//                                                                                //
//    Copyright 2001                                                              //
//                                                                                //
//    BHP Software.  All rights reserved.                                         //
//                                                                                //
//    Permission is hereby granted, free of charge, to any person obtaining       //
//    a copy of this software and associated documentation files (the             //
//    "Software") to deal in the Software, without restriction, except as         //
//    hereinafter provided, including without limitation the rights to use,       //
//    copy, modify merge, publish, and distribute the Software and to permit      //
//    persons to whom the Software is furnished to do so, provided that the       //
//    above copyright notice and this permission notice appear in all copies      //
//    of the Software and that both the above copyright notice and this           //
//    permission notice appear in supporting documentation.  No charge may        //
//    be made for any redistribution of the Software, including modified or       //
//    merged versions of the Software.  The complete source code must be          //
//    included in any distribution. For an executable file, complete source       //
//    code means the source code for all modules it contains                      //
//                                                                                //
//    Modified or merged versions of the Software must be provided to the         //
//    Software Manager, regardless of whether such modified or merged             //
//    versions are distributed to others.                                         //
//                                                                                //
//    THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS      //
//    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                  //
//    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGMENT        //
//    OF THIRD PARTY RIGHTS.  IN NO EVENT SHALL THE COPYRIGHT HOLDER              //
//    INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM OR ANY SPECIAL              //
//    INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING      //
//    FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,        //
//    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION        //
//    WITH THE USE OF PERFORMANCE OF THIS SOFTWARE.                               //
//                                                                                //
//    The name of the copyright holder shall not be used in advertising or        //
//    otherwise to promote the use or other dealings in this Software,            //
//    without prior written consent of the copyright holder.                      //
//                                                                                //
//********************************************************************************//


package au.csiro.JamesGunning.waveletextraction;


import java.util.*;
import java.io.*;
import static java.lang.Math.*;


import cern.colt.matrix.*;
//import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import cern.colt.map.*;

import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import static au.csiro.JamesGunning.lib.util.MathUtil.*;
import au.csiro.JamesGunning.lib.util.*;
/**
   Basic model parameters needed for each distinct wavelet span and position. Will have HashMaps etc
   for mapping of index to properties etc
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class ModelVariables
{
  WaveletExtractor _WE;

  Wavelet wavelet;

  OpenIntObjectHashMap mapFromIndex;
  HashMap mapToIndex;
  HashMap mapToName;

  // caches
  ArrayList<Property> amplitudeCachingTestProperties;
  double[] amplitudeCachingOldValues;
  ArrayList<Property> reflectivityCachingTestProperties;
  double[] reflectivityCachingOldValues;

  /**
	 Typical coefficient used to estimate noise from typical seismic amplitude
  */
  public final static double TYP_NOISE_TO_AMP_RATIO = 0.4;
  /**
	 Useful boolean flag for turning off noise estimation if required. It gets fixed at the
	 prior mean value
  */
  boolean fixNoise = false;

  // Auxiliary (non-synthetic) constraint variables
  /**
	 Dimension of Auxiliary constraints
  */
  int dimConstraints = 0;
  /**
	 Snapshot arrays of constraint errors at current model configuration
  */
  double[] constraintErrors;
  double chiSqPhase=0,chisqVelMismatch=0,chisqPeakMismatch=0;

  // full number of parameters in this model, including noise
  int nParams;
  // full number of Gaussian (untransformed) parameters in this model
  int nGaussian;
  
  int firstTDparamIndex=-1;
  int lastTDparamIndex=-1;
  int firstTregistrationParamIndex=-1;
  int lastTregistrationParamIndex=-1;
  int firstXYpositioningParamIndex=-1;
  int lastXYpositioningParamIndex=-1;

  // log of Laplace weight of mode
  double logModeWeight;
  // 0.5 log(det(Ceff)) of mode
  double halfLogDetCeff;

  // actual posterior probability weight after normalisation etc
  double weight;

  // CDF
  double CDF;

  // prior covariance of subvector of parameters with Gaussian prior
  DoubleMatrix2D C_prior;
  // inverse of this
  DoubleMatrix2D Cinv_prior;
  // determinant of it 
  double logDetC_prior_only;
  // and prior mean of subvector of parameters with Gaussian prior
  DoubleMatrix1D mu_prior;

  // Cholesky L factor of posterior covariance of mode
  double[] CtildeInvDiag;
  // Cholesky L factor of posterior covariance of mode
  DoubleMatrix2D L_post;
  // Permutation proper to the Cholesky factor L_post
  int[] L_Permutation;
  // Diagonal of posterior covariance of mode
  double[] CtildeDiag;

  // posterior mean of mode
  DoubleMatrix1D mu_post;

  // posterior noiseCorrCoeffMatrix
  double[][] noiseCorrCoeffMatrix=null;

  // Typical scales for the full suite of parameters
  double[] typicalScale;
  // Typical initial guess for the full suite of parameters
  double[] initialGuess;

  //some useful efficiency flags?
  boolean haveXYseismicErrors=false;
  boolean haveTseismicErrors=false;

  /**
	 Allowable start time of this wavelet extraction
  */
  double tStart;

  /**
	 Allowable end time of this wavelet extraction
  */
  double tEnd;

  int phaseConstraintIndxStart;
  int phaseConstraintIndxEnd;

  // whether model has anisotriopy
  boolean hasAnisotropy = false;

  // checkshot model index if used
  int checkshotModelIndex=0;  

  // Useful stuff for wavelet peak handling
  boolean modeHasPeakAmbiguity=false;
  boolean useNearestPeakArrivalConstraint;
  boolean constrainPeakArrivalSign;
  boolean peakArrivalPositive;

  // If set, constrain the phase via ratios to avoid branch discontinuities etc
  public static boolean CONSTRAIN_PHASE_VIA_PHASE_LINE=true;


  /**
	 Constructor instantiates necessary wavelet
  */
  public ModelVariables(WaveletExtractor WE, int waveletCoarseIndexMin, int waveletCoarseIndexMax, 
						double fine_scale_dt, double coarse_scale_dt, boolean includeAnisotropy,
						int checkshotModelIndex)
  {
	this._WE = WE;
	wavelet = new Wavelet(_WE, waveletCoarseIndexMin,waveletCoarseIndexMax,fine_scale_dt,coarse_scale_dt);

	constrainPeakArrivalSign = _WE.constrainPeakArrivalSign;
	peakArrivalPositive = _WE.peakArrivalPositive;

	// Set up tStart, tEnd
	if (_WE.assumeNoReflectionsOutsideLogs)
	  {
		this.tStart = _WE.tStart;
		this.tEnd = _WE.tEnd;
	  }
	else 
	  {
		if (_WE.convModel.RIGOROUSLY_EXCLUDE_EDGE_REGIONS)
		  {
			// The wavelet.*Duration are time distances to the last 
			// samples at the edges with 
			// non-zero values, so to rigorously preclude 
			// discontinuities arising from reflections wandering into
			// the interval [_WE.tStart, _WE,tEnd], we have to add
			// dt padding (since the Lagrange 4pt formula propagates things  an extra dt
			this.tStart = _WE.tStart + wavelet.codaDuration + _WE.dt;
			this.tEnd = _WE.tEnd - wavelet.precursorDuration - _WE.dt;
		  }
		else
		  {
			this.tStart = _WE.tStart;
			this.tEnd = _WE.tEnd;
		  }
	  }


	fixNoise = _WE.fixNoise;

	// anisotropy flag
	if (includeAnisotropy && _WE.doAnisotropicModel) 	hasAnisotropy = true;

	// model index of checkshotModelSelect
	this.checkshotModelIndex = checkshotModelIndex;
  }


  /**
	 Set checkshot for this model according to 
	 checkshotModelIndex established at instantiation. 
	 Does nothing if checkshotModelSelect not in use
  */
  void updateCheckshotForModelSelection()
  {
	if (_WE.multiModelSetUp.checkshotModelSelect)
	  {
		// Set up new checkshot if varying checkshot model		
		int[] model_index = _WE.multiModelSetUp.checkshot_model_select_index_array.get(checkshotModelIndex);
		for (int i=0; i<_WE.m.nWells; i++)
		  { 
			if (_WE.multiModelSetUp.checkshotModelSelectWells[i])
			  {
				// set up new checkshot for model checkshotModelIndex, at well i
				_WE.chkShotandMarkerData[i].autoGenSetCheckshotModel(model_index[i]);
				// reconstruct t2d maps etc
				_WE.chkShotandMarkerData[i].setUpCheckShotAndMarkerData(_WE.logData[i].survey);
				// Update maps
				_WE.m.timeToDepthMaps[i] = new TimeToDepthMap(_WE, _WE.chkShotandMarkerData[i].t2d,i);
			  }
		  }
	  }
  }

  /**
	 Assembly of basic model for a given wavelet <p>
	 The model is the list of all parameters that vary. These parameters have distinct 
	 Property objects associated with them that control the getting and setting. For the
	 linear algebra, the parameters have to be assembled into a suitable vector, so we need
	 ways of managing this vector, getting indices, setting properties, etc. This is usually done
	 with HashMaps.

	 The model vector is the union of<p>
	 (1) wavelet coarse coefficients<p>
	 (2) Time to depth knots<p>
	 (3) Reflectivity errors<p>
	 (4) Noise levels per stack<p>
	 (5) Registration errors in t per stack<p>
	 (6) AVO coefficient B (and Bfar if free)<p> 
	 (6.5) Anisotropy delta of shales if possible<p>
	 (7) Stack scale for far stack if free
	 (8) x,y location errors per stack<p>


	 Also, since the prior model is Gaussian, this is the place to assemble the prior covariances as the 
	 parameter vector is assembled. The model vector is assembled with The Guassian parameters first, and
	 any non-Gaussian parameters later.
 
  */
  void initialiseModel(boolean createWaveletCoefficientModelOnly)
  {
	int i,ii,j,k,p;
	double var,val,tmp,t_mean,typicalAmp,typicalWaveletCoeff;
	int synchIndex;
	PropertyMapper pM=null;
	Integer I;
	int verb = _WE.verbosity;
	DoubleMatrix2D Cp_invWavCoeffs=null;


	// OK run over all variables to get indexing in place
	mapFromIndex = new OpenIntObjectHashMap();
	mapToIndex = new HashMap();
	mapToName = new HashMap();
	
	// caches for speed
	if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
	  {
		amplitudeCachingTestProperties = new ArrayList<Property>();
	  }
	if (Model.USE_CACHED_REFLECTIVITY)
	  {
		reflectivityCachingTestProperties = new ArrayList<Property>();
	  }

	// Set up Gaussian priors
	if (verb>3) _WE.EH.println("Model variables & prior means"); 
	// (1) Wavelet coefficients
	k=0;  // runs over model parameters
    for (i=0; i<wavelet.valuesCoarseProp.length; i++,k++)
	  {
		mapFromIndex.put(k, wavelet.valuesCoarseProp[i]);
		mapToIndex.put(wavelet.valuesCoarseProp[i], new Integer(k));
		mapToName.put(wavelet.valuesCoarseProp[i],"wavelet_coeff_"+k);
		if (verb>3) _WE.EH.println("%3d wavelet coeffs, %-10.8g",k,wavelet.valuesCoarseProp[i].mean); 
	  }
	nParams = k;

	// (2) Time to depth knots. 
	// Skip if t2d is fixed
	if (!_WE.fixCheckshots && !createWaveletCoefficientModelOnly)
	  {
		firstTDparamIndex=k;
		for (i=0; i<_WE.m.nWells; i++)
		  {
			TimeToDepthMap t2dmap = _WE.m.timeToDepthMaps[i];
			for (j=0;j<t2dmap.tProp.length;j++)
			  {
				if (t2dmap.sigma_t[j] > 0)
				  {
					mapFromIndex.put(k,t2dmap.tProp[j]);
					mapToIndex.put(t2dmap.tProp[j],new Integer(k));
					mapToName.put(t2dmap.tProp[j],"time_to_depth_knot_"+j+"_well_"+i);
					if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
					  {
						amplitudeCachingTestProperties.add(t2dmap.tProp[j]); // CACHE
					  }
					if (Model.USE_CACHED_REFLECTIVITY)
					  {
						reflectivityCachingTestProperties.add(t2dmap.tProp[j]); // CACHE
					  }
					nParams++;
					if (verb>3) 
					  {
						// t_mean = _WE.chkShotandMarkerData[i].t2d[j].T_mean;
						_WE.EH.println("%3d time to depth knots, well %s, %-10.8g",
									   k, _WE.wellInfo[i].well_name, t2dmap.tProp[j].mean); 
					  }
					k++;
				  }
			  }
		  }
		lastTDparamIndex=k-1;
	  }

	// (3) Reflectivity or Impedance errors
	// ....not used

	// (4) time registration error for different stacks. Currently modelled independently at this stage
	haveTseismicErrors = false;

	if (!createWaveletCoefficientModelOnly && !_WE.fixMinicubeRegistrations)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			for (j=0; j<_WE.m.nStacks; j++)
			  {
				if (_WE.wellInfo[i].synchronizeWellName == null)
				  {
					if (j==0 || !_WE.synchroniseStackTregistrationErrors)
					  {
						// A "master" well. We can add genuine entries
						if ((_WE.seismicInfo[i][j]).systematicTProp.sigma > 0)
						  {
							mapFromIndex.put(k, (_WE.seismicInfo[i][j]).systematicTProp);
							mapToIndex.put((_WE.seismicInfo[i][j]).systematicTProp, new Integer(k));
							mapToName.put((_WE.seismicInfo[i][j]).systematicTProp,
										  "t_shift_"+_WE.seismicInfo[i][j].filename);
							
							if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
							  {
								amplitudeCachingTestProperties.add((_WE.seismicInfo[i][j]).systematicTProp); // CACHE
							  }

							if (verb>3) _WE.EH.println("%3d t_shift, stack %3d, %s %-10.8g",
													   k, j+1, _WE.seismicInfo[i][j].filename,
													   (_WE.seismicInfo[i][j]).systematicTProp.mean);
							haveTseismicErrors = true;
							// keep track of block of params controlling these;
							if (firstTregistrationParamIndex == -1) firstTregistrationParamIndex = k;
							nParams++;
							k++;
						  }
					  }
					else
					  {
						// a synchronised image...						
						// synchIndex = ((Integer)mapToIndex.get((_WE.seismicInfo[i][0]).systematicTProp)).intValue();
						if (verb>3) _WE.EH.println(" .. t_shift, stack %3d, %s, synchronised to first stack",
												   j+1, _WE.seismicInfo[i][j].filename);
					  }
				  }
				else
				  {
					// a synchronised image...
					synchIndex = _WE.wellInfo[i].synchronizeWellIndex;
					if ((_WE.seismicInfo[synchIndex][_WE.synchroniseStackTregistrationErrors ? 0 : j]).systematicTProp.sigma !=0 )
					  {
						if (verb>3) _WE.EH.println(" .. t_shift, %s, synchronised to well %s",
												   _WE.wellInfo[i].well_name,
												   _WE.wellInfo[synchIndex].well_name);
					  }
				  }
			  }
		  }
	  }
	if (haveTseismicErrors) lastTregistrationParamIndex = k-1;


	// (5) positioning error dx dy for each well
	haveXYseismicErrors = false;
	if (!createWaveletCoefficientModelOnly && !_WE.fixWellPositions)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			if (_WE.wellInfo[i].synchronizeWellName == null)
			  {
				// A "master" well. We can add genuine entries
				if (_WE.wellInfo[i].systematicXProp.sigma !=0 )
				  {
					mapFromIndex.put(k,_WE.wellInfo[i].systematicXProp);
					mapToIndex.put(_WE.wellInfo[i].systematicXProp,new Integer(k));
					mapToName.put(_WE.wellInfo[i].systematicXProp,"X_shift_well_"+i);
					if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
					  {
						amplitudeCachingTestProperties.add(_WE.wellInfo[i].systematicXProp);
					  }
					if (verb>3) _WE.EH.println("%3d X_shift, well %s, %-10.8g",
											   k, _WE.wellInfo[i].well_name,
											   _WE.wellInfo[i].systematicXProp.mean); 
					haveXYseismicErrors = true;
					if (firstXYpositioningParamIndex == -1) firstXYpositioningParamIndex = k;
					nParams++;
					k++;
				  }
				if (_WE.wellInfo[i].systematicYProp.sigma !=0 )
				  {
					mapFromIndex.put(k,_WE.wellInfo[i].systematicYProp);
					mapToIndex.put(_WE.wellInfo[i].systematicYProp,new Integer(k));
					mapToName.put(_WE.wellInfo[i].systematicYProp,"Y_shift_well_"+i);
					if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
					  {
						amplitudeCachingTestProperties.add(_WE.wellInfo[i].systematicYProp);
					  }
					if (verb>3) _WE.EH.println("%3d Y_shift, well %s, %-10.8g",
											   k, _WE.wellInfo[i].well_name,
											   _WE.wellInfo[i].systematicYProp.mean); 
					haveXYseismicErrors = true;
					if (firstXYpositioningParamIndex == -1) firstXYpositioningParamIndex = k;
					nParams++;
					k++;
				  }
			  }
			else
			  {
				// a synchronised image...
				synchIndex = _WE.wellInfo[i].synchronizeWellIndex;
				if (_WE.wellInfo[synchIndex].systematicXProp.sigma !=0 )
				  {
					if (verb>3) _WE.EH.println(" .. X_shift, well %s, synchronised to well %s",
											   _WE.wellInfo[i].well_name,
											   _WE.wellInfo[synchIndex].well_name);
				  }
				if (_WE.wellInfo[synchIndex].systematicYProp.sigma !=0 )
				  {
					if (verb>3) _WE.EH.println(" .. Y_shift, well %s, synchronised to well %s",
											   _WE.wellInfo[i].well_name,
											   _WE.wellInfo[synchIndex].well_name);
				  }
			  }
		  }
	  }
	if (haveXYseismicErrors) lastXYpositioningParamIndex = k-1;



	// (6) B coefficient in AVO equations if required
	if (!createWaveletCoefficientModelOnly && (_WE.m.B_AVO[0].sigma>0))
	  {
		for (j=0; j<_WE.m.nStacks; j++)
		  {
			mapFromIndex.put(k,_WE.m.B_AVO[j]);
			mapToIndex.put(_WE.m.B_AVO[j],new Integer(k));
			mapToName.put(_WE.m.B_AVO[j],"AVO_B_coeff_"+(j+1));
			if (Model.USE_CACHED_REFLECTIVITY && !_WE.doIsoZoeppritz)
			  {
				// Zoeppritz doesn't use a "B" factor
				reflectivityCachingTestProperties.add(_WE.m.B_AVO[j]);
			  }
			if (verb>3) _WE.EH.println("%3d AVO B term_%d, %-10.8g",k, (j+1), _WE.m.B_AVO[j].mean); 
			nParams++;
			k++;
			// only one B_AVO unless B's are all free
			if (!_WE.freeAVO) break;
		  }
		// extra factor if stack B's free
		// if (_WE.freeAVO)
		// 		  {
		// 			mapFromIndex.put(k,_WE.m.B_AVO_far);
		// 			mapToIndex.put(_WE.m.B_AVO_far,new Integer(k));
		// 			mapToName.put(_WE.m.B_AVO_far,"AVO_B_far_coeff");
		// 			if (Model.USE_CACHED_REFLECTIVITY)
		// 			  {
		// 				reflectivityCachingTestProperties.add(_WE.m.B_AVO_far);
		// 			  }
		// 			if (verb>3) _WE.EH.println("%3d AVO B term (far stack), %-10.8g",k, _WE.m.B_AVO_far.mean); 
		// 			nParams++;
		// 			k++;
		// 		  }
	  }

	// (6.5) Anisotropy
	if (hasAnisotropy) // _WE.doAnisotropicModel)
	  {
		if (!createWaveletCoefficientModelOnly && (_WE.m.deltaThomsen_shale.sigma > 0))
		  {
			mapFromIndex.put(k,_WE.m.deltaThomsen_shale);
			mapToIndex.put(_WE.m.deltaThomsen_shale,new Integer(k));
			mapToName.put(_WE.m.deltaThomsen_shale,"Thomsen_shale_delta");
			if (Model.USE_CACHED_REFLECTIVITY)
			  {
				reflectivityCachingTestProperties.add(_WE.m.deltaThomsen_shale);
			  }
			if (verb>3) _WE.EH.println("%3d Anisotropy delta term, %-10.8g",k, _WE.m.deltaThomsen_shale.mean); 
			nParams++;
			k++;
		  }
	  } 

	// (7) Far stack wavelet stretching. Can't cache this I think
	if (!createWaveletCoefficientModelOnly && (_WE.m.nStacks>1) && _WE.useStretchedModel)
	  {
		// use a common remapper-tool for all stretched wavelets
		if (_WE.m.farStackWaveletTimeStretchSigma>0)
		  {
			_WE.m.farStackWaveletTimeStretchMapper = new PropertyMapper(PropertyMapper.LOGIT, 1, Wavelet.MAX_WAVELET_STRETCH);
			pM = _WE.m.farStackWaveletTimeStretchMapper; 
		  }
		
		for (j=0; j<(_WE.m.nStacks-1); j++)
		  {
			if (_WE.m.farStackWaveletTimeStretchSigma>0)
			  {
				// Use broad uninformative prior: use *likelihood* to enforce true prior
				// (nonlinearity requires this kind of treatment)
				// Treat the property _WE.m.farStackWaveletTimeStretch as a logitized 
				// variable. Real values will have to be recovered from it every time it is used.

				// Logit transform goes silly at endpoint 1.0, so protect against this...
				double protectedMean = _WE.m.farStackWaveletTimeStretchMean;
				if (_WE.m.farStackWaveletTimeStretchMean <= 1.01) protectedMean = 1.1;
				_WE.m.farStackWaveletTimeStretch[j].setMean(pM.getTransformedValue(protectedMean));
				_WE.m.farStackWaveletTimeStretch[j].setSigma(10.0); // so large it makes the prior here uninformative
				mapFromIndex.put(k,_WE.m.farStackWaveletTimeStretch[j]);
				mapToIndex.put(_WE.m.farStackWaveletTimeStretch[j],new Integer(k));
				mapToName.put(_WE.m.farStackWaveletTimeStretch[j],"far_wavelet_stretch_stack_"+(j+2));
				if (verb>3) _WE.EH.println("%3d far_wavelet_stretch_stack_%d, %-10.8g",k,(j+2), _WE.m.farStackWaveletTimeStretchMean); 
				nParams++;
				k++;
			  }
			if (_WE.m.farStackWaveletAmplitude[j].sigma > 0)
			  {
				mapFromIndex.put(k,_WE.m.farStackWaveletAmplitude[j]);
				mapToIndex.put(_WE.m.farStackWaveletAmplitude[j],new Integer(k));
				mapToName.put(_WE.m.farStackWaveletAmplitude[j],"far_wavelet_amp_stack_"+(j+2));
				if (verb>3) _WE.EH.println("%3d far_wavelet_amp_stack_%d term, %-10.8g",k,(j+2),_WE.m.farStackWaveletAmplitude[j].mean); 
				nParams++;
				k++;
			  }
		  }
	  }


	// (7.5) Saturation parameters for time-lapse well adjustments
	if (!createWaveletCoefficientModelOnly && _WE.useSaturationModel)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			SaturationInfo _sInfo = _WE.wellInfo[i].saturationInfo;
			if (_sInfo==null) continue;
			double mu_S = min(max(_sInfo.Sa,0.01),0.99);  // "after" saturation. 
			if(_sInfo.sig_Sa <= 0)
			  {
				_sInfo.satProp.setMean(mu_S); // set default mean & value if fixed
			  }
			else // variable: (_sInfo.satProp.sigma > 0 )
			  {
				pM = _sInfo.satPropMapper = new PropertyMapper(PropertyMapper.LOGIT,0.0,1.0); // FIXME USER SPEC ENDPOINTS

				// Logit transform goes silly at endpoints - protect against this...
				double sigma_L = sqrt(sqr(_sInfo.sig_Sa)/(mu_S*(1-mu_S)));  // approximate equivalent std deviation of logitised variable

				//_sInfo.satProp.setSigma(PropertyMapper.LOGIT_NORM_TO_UNIFORM_APPROX_STDDEV); // FIXME???(1.4) makes prior here uninformative
				// sigma_L is the width required from rough analytical expansions based on small sig_Sa
				// if sig_Sa is large, this makes no sense. The remapped prior also goes multimodal for very wide
				// sigma_L, so apply cutoff:
				sigma_L = min(sigma_L, PropertyMapper.LOGIT_NORM_TO_UNIFORM_APPROX_STDDEV);
				_sInfo.satProp.setMean(pM.getTransformedValue(mu_S));
				_sInfo.satProp.setSigma(sigma_L);


				//_WE.EH.println("satProp DISTRIBUTION: m=%4g stddev=%5g%n",_sInfo.satProp.mean,_sInfo.satProp.sigma);
				mapFromIndex.put(k,_sInfo.satProp);
				mapToIndex.put(_sInfo.satProp,new Integer(k));
				mapToName.put(_sInfo.satProp,"saturation_after_well_"+i);
				if (verb>3) _WE.EH.println("%3d saturation_after, well %s, %-10.8g",
										   k, _WE.wellInfo[i].well_name,mu_S);										  
				nParams++;
				k++;
			  }
		  }
	  }


	nGaussian = nParams;

	// (8) Stack noise levels....these do not have Gaussian priors!!
	// These must be the last parameters!!!!!
	// first, simple sanity check...
	if (fixNoise)
	  {
		if ((_WE.fixNoiseValue > 0.0) && 
			_WE.fixNoiseValue < ModeFinder.MINIMUM_TOLERABLE_SN_RATIO * _WE.multiModelSetUp.typicalRMSSeismicAmplitude)
		  {
			_WE.EH.die_err("ERROR: Fixed value of noise provided ("+_WE.fixNoiseValue +
						   ") is impractically small. Try a larger value");
		  }
	  }
	if (!fixNoise && !createWaveletCoefficientModelOnly)
	  {	
		for (i=0; i<_WE.m.nStacks; i++,k++)
		  {
			mapFromIndex.put(k, _WE.m.seismic_noise[i]);
			mapToIndex.put(_WE.m.seismic_noise[i],new Integer(k));
			mapToName.put(_WE.m.seismic_noise[i],"noise_stack_"+i);
			nParams++;  // but not nGaussian!
			if (_WE.ignoreSeismic)
			  {
				var = 1.0;
			  }
			else if (_WE.fixNoiseValue>0)
			  {
				var = _WE.fixNoiseValue;
			  }
			else
			  {
				// default: set scale to approx 1/4 of peak amplitude
				var = TYP_NOISE_TO_AMP_RATIO*_WE.multiModelSetUp.typicalRMSSeismicAmplitude;
			  }
			if (verb>3) 
			  {
				_WE.EH.println("%3d noise level: (typical amp:) %-10.8g",k,var); 
			  }
		  }
	  }

	// Allocate caching stuff
	if (ConvolutionalModel.USE_CACHED_AMPLITUDES) 
	  {
		amplitudeCachingTestProperties.trimToSize();
		amplitudeCachingOldValues = new double[amplitudeCachingTestProperties.size()];
		_WE.EH.println("Amplitude caching: "+amplitudeCachingOldValues.length+" parameters to check");
	  }
	if (Model.USE_CACHED_REFLECTIVITY)
	  {
		reflectivityCachingTestProperties.trimToSize();
		reflectivityCachingOldValues = new double[reflectivityCachingTestProperties.size()];
		_WE.EH.println("Reflectivity caching: "+reflectivityCachingOldValues.length+" parameters to check");
	  }



	// DON'T ADD ANY MORE PARAMETERS HERE!!!! DO IT BEFORE THE NOISE LEVELS!





	// Now only some of the model parameters have Gaussian priors.
	// This is every thing but the Seismic noise levels.

	// Now construct prior covariance, its inverse & mean for these parameters
	C_prior = new DenseDoubleMatrix2D(nGaussian,nGaussian);
	Cinv_prior = new DenseDoubleMatrix2D(nGaussian,nGaussian);
	int nW = wavelet.valuesCoarseProp.length;
	mu_prior = new DenseDoubleMatrix1D(nGaussian);

	// typical scales and initial Guess for all parameters
	typicalScale = new double[nParams];
	initialGuess = new double[nParams];

	//(1) Wavelet coefficients
	// assuming typical seismic amplitudes scale like reflectivity*wavelet,
	// => wavelet amp ~ typicalRMSSeismicAmplitude/(typical Reflection coeff),
	// so a loose prior is about...
    typicalWaveletCoeff = _WE.multiModelSetUp.typicalRMSSeismicAmplitude/
	             FundamentalConstants.TYPICAL_REFLECTION_COEFF;
	// My guess is that peak events are often about 3-4 times the RMS value, so
	// a peak wavelet coefficient would be about 3*typicalWaveletCoeff.
	// We should arrange the prior so this is about 1 std deviation. Hence 10 ~ 3^2, giving
	double waveletCoeffVariance;
	waveletCoeffVariance = MathUtil.sqr(FundamentalConstants.TYPICAL_PEAK_TO_RMS_SEISMIC_RATIO 
										* typicalWaveletCoeff);
	if (verb>3) _WE.EH.println("Prior-wavelet variance=%-10.8g",waveletCoeffVariance); 

	// determinant to be accrued
	logDetC_prior_only=0.0;

	k=0;  // runs over model parameters
	// Coarse wavelet coefficients

	for (i=0; i<nW; i++,k++)
	  {
		if (_WE.useWaveletPrior)
		  {
			// override prior from previous derivation
			waveletCoeffVariance = MathUtil.sqr(_WE.waveletPrior.amp_stddev.get(i));
		  }
		// Default case: wavelet coeffs are uncorrelated and sloppy
		C_prior.setQuick(k,k,waveletCoeffVariance);
		Cinv_prior.setQuick(k,k,1.0/waveletCoeffVariance);
		logDetC_prior_only += Math.log(waveletCoeffVariance);

		// means
		if (_WE.useWaveletPrior)
		  {
			val = _WE.waveletPrior.amp_means.get(i);
		  }
		else
		  {
			// with mean zero, 
			val = 0.0;
		  
			// FIXME PHASE: except if there is a phase constraint? (The likelihood has 
			// a singularity at 0: avoid this by starting with a non-zero wavelet)
			if (_WE.constrainPhase) val = 0.1*typicalWaveletCoeff;
		  }
		mu_prior.setQuick(k,val);
		// and this scale (about half peak response?)
		typicalScale[k] = 0.5*typicalWaveletCoeff;
		initialGuess[k] = val;
	  }

	if (_WE.bbTail) // EXPERIMENTAL bband ideas
	  { 
		double lam,v,vt; 
		int jj,kk,iTail;
		v = waveletCoeffVariance;      // large variance, inner region

		vt = sqr(typicalWaveletCoeff); // smaller variance for tail region
		
		// bbTail_T in ms, dt_sparse in s. Get coarse sample iTail from where heavy correlations kick in
		iTail = (int)Math.ceil(_WE.bbTail_T/(1000*_WE.multiModelSetUp.dt_sparse)); // usually about eg 3 
		lam = _WE.bbTail_W;                                                 // about 200
		int N = max(-wavelet.firstCoarseSample, wavelet.finalCoarseSample);

		// kriging covariance so fixed 0-data point at coefficient wavelet.firstCoarseSample
		N = N-iTail+1;
		DoubleMatrix2D C_ptail=null;
		CholeskyDecomposition CD;
		DoubleMatrix2D L;
		DoubleMatrix2D C_ptail_kriged=null;
		if(N>1)
		  {
			C_ptail = new DenseDoubleMatrix2D(N,N);
			for(j=0;j<N;j++) 
			  {
				for(k=j;k<N;k++) 
				  {
					C_ptail.setQuick(j,k,tmp = Math.exp(-abs(j-k)/lam)); 
					if(j!=k) C_ptail.setQuick(k,j,tmp); 
				  }
			  }
			CD = new CholeskyDecomposition(C_ptail);
			L = CD.getL();
			DoubleMatrix2D P_tail = CD.solve(DenseDoubleMatrix2D.identity(N));
			CD = new CholeskyDecomposition(P_tail.viewPart(1,1,N-1,N-1));
			L = CD.getL();
			C_ptail_kriged = CD.solve(DenseDoubleMatrix2D.identity(N-1));
		  }
		for(j=0,jj=wavelet.firstCoarseSample+1; j<nW; j++,jj++) 
		  {
			for(k=0,kk=wavelet.firstCoarseSample+1; k<nW; k++,kk++) 
			  {
				if ((jj>=iTail) && (kk>=iTail))
				  {
					// dx = wavelet.finalCoarseSample-jj
					// j_k = dx-1 
					val = vt * C_ptail_kriged.getQuick(wavelet.finalCoarseSample-jj-1,
												  wavelet.finalCoarseSample-kk-1);
				  }
				else if ((jj<=-iTail) && (kk<=-iTail))
				  {
					// dx = jj-wavelet.firstCoarseSample
					// j_k = dx-1 
					val = vt * C_ptail_kriged.getQuick(jj-wavelet.firstCoarseSample-1,
												  kk-wavelet.firstCoarseSample-1);
				  }
				else
				  {
					val = (j==k) ? v : 0;
				  }				  
				C_prior.setQuick(j,k, val);
			  }
		  }
		CD = new CholeskyDecomposition(C_prior.viewPart(0,0,nW,nW));
		L = CD.getL();
		// recompute logDetC_prior_only while we're here
		for(k=0,val=0; k<nW; k++) val += Math.log(L.getQuick(k,k));
		logDetC_prior_only = 2.0*val;
		Cinv_prior.viewPart(0,0,nW,nW).assign(CD.solve(DenseDoubleMatrix2D.identity(nW)));

		//Util.debugPrintMatrixToFile("CP.txt",C_prior);
		//Util.debugPrintMatrixToFile("CPI.txt",Cinv_prior);
	  }
	if (_WE.zeroConstrainPmoments)
	  {
		// Enforces constraints 
		// \sum_istart^iend i^p wavelet_coeffs[i]=0 (zero mean)
		// via a pseudo-likelihood trick Cp_eff =(Cp^{-1}+Xp^T.Xp)^{-1}, X_p = moment constraints
		// In below, 0.1 = relative-stddev of  \sum wavelet_coeffs[i]
		double ep = sqr(1.0/(0.1 * typicalWaveletCoeff)); 
		val = 1.0/waveletCoeffVariance;
		DoubleMatrix2D Xp  = new DenseDoubleMatrix2D(_WE.zeroConstrainPmoments_p+1, nW);	
		DoubleMatrix2D C_prior_eff  = new DenseDoubleMatrix2D(nW,nW);	

		for(p=0; p<= _WE.zeroConstrainPmoments_p; p++)
		  {
			//firstCoarseSample is at fixed endpoint entry 0, first free index at firstCoarseSample+1
			for(j=0,k=wavelet.firstCoarseSample+1; j<nW; j++,k++) 
			  {
				Xp.setQuick(p,j,ipow(k,p)); // X_{p,j}= k^p. typically k = -5...5, say
			  }
		  }
		if (!_WE.bbTail)
		  {
			for(j=0;j<nW;j++) C_prior_eff.setQuick(j,j,val);  // Basic precision
		  }
		else
		  {
			C_prior_eff.assign(Cinv_prior.viewPart(0,0,nW,nW));// from above
		  }
		// Util.debugPrintMatrixToFile("CP.txt",C_prior_eff);
		// Util.debugPrintMatrixToFile("X.txt",Xp);
		// Cinv_prior =(Cp^{-1}+Xp^T.Xp)
		SeqBlas.seqBlas.dgemm(true,false,ep,Xp,Xp,1.0,C_prior_eff);
		// fill leading sub-block of Cinv_prior
		for(j=0;j<nW;j++)
		  {
			for(k=0;k<nW;k++)Cinv_prior.setQuick(j,k,C_prior_eff.getQuick(j,k));
		  }
		// Cp_eff =(Cp^{-1}+Xp^T.Xp)^{-1},
		CholeskyDecomposition CD = new CholeskyDecomposition(C_prior_eff);
		DoubleMatrix2D L = CD.getL();
		// recompute logDetC_prior_only while we're here
		for(k=0,val=0; k<nW; k++) val += Math.log(L.getQuick(k,k));
		logDetC_prior_only = -2*val;
		C_prior_eff = CD.solve(DenseDoubleMatrix2D.identity(nW));
		// fill leading sub-block of C_prior with C_prior_eff
		for(j=0;j<nW;j++)
		  {
			for(k=0;k<nW;k++) C_prior.setQuick(j,k,C_prior_eff.getQuick(j,k));
		  }
		CD = null;
		//Util.debugPrintMatrixToFile("CP.txt",C_prior);
	  }
	
	// (2) Time to depth knots. One per well
	// omit these if time to depth is fixed
	if (!_WE.fixCheckshots && !createWaveletCoefficientModelOnly)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			TimeToDepthMap t2dmap = _WE.m.timeToDepthMaps[i];
			for (j=0;j<t2dmap.tProp.length;j++)
			  {
				if (t2dmap.sigma_t[j]>0)
				  {
					// variance set by  knot sloppiness
					var = t2dmap.sigma_t[j] * t2dmap.sigma_t[j];
					C_prior.setQuick(k,k,var);
					Cinv_prior.setQuick(k,k, 1.0/var);
					logDetC_prior_only += Math.log(var);
					// and mean by the actual value
					t_mean = t2dmap.tProp[j].mean;
					mu_prior.setQuick(k,t_mean);

					// and this scale. Make it the expected uncertainty
					// but if this is wide, settle for dt
					typicalScale[k] = Math.min(t2dmap.sigma_t[j],_WE.dt); 
					initialGuess[k] = t_mean;
					k++;
				  }
			  }
		  }
	  }
	else
	  {
		// ensure properties are reset, even if not using them, as previous models may have tweaked them
		for (i=0; i<_WE.m.nWells; i++)
		  {
			TimeToDepthMap t2dmap = _WE.m.timeToDepthMaps[i];
			for (j=0;j<t2dmap.tProp.length;j++)
			  {
				t2dmap.tProp[j].resetToMean();
			  }
		  }
	  }

	// Must add these Gaussian pieces.........FIXME!!!
  	// (3) Reflectivity or Impedance errors


	// (4) registration error for different stacks
	if (haveTseismicErrors)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			// only add these if a master...
			if (_WE.wellInfo[i].synchronizeWellName == null)
			  {
				for (j=0; j<_WE.m.nStacks; j++)
				  {
					if (j==0 || !_WE.synchroniseStackTregistrationErrors)
					  {
						// only add master params if stacks are synchronised
						if (_WE.seismicInfo[i][j].systematicTProp.sigma > 0)
						  {
							// variance set by  systematic sloppiness
							var = MathUtil.sqr(_WE.seismicInfo[i][j].systematicTProp.sigma);
							C_prior.setQuick(k,k,var);
							Cinv_prior.setQuick(k,k, 1.0/var);
							logDetC_prior_only += Math.log(var);
							// and mean is the value (will be zero)
							mu_prior.setQuick(k,_WE.seismicInfo[i][j].systematicTProp.mean);
							// and this scale. Make it the expected uncertainty
							// these get bound tightly, so take scale down here.... by 0.1
							typicalScale[k] = 0.1*Math.min(_WE.seismicInfo[i][j].systematicTProp.sigma , _WE.dt);
							initialGuess[k] = _WE.seismicInfo[i][j].systematicTProp.mean;
							k++;
						  }
					  }
				  }
			  }
		  }
	  }
	else
	  {
		// ensure properties are reset, even if not using them, as previous models may have tweaked them
		for (i=0; i<_WE.m.nWells; i++)
		  {
			// only add these if a master...
			if (_WE.wellInfo[i].synchronizeWellName == null)
			  {
				for (j=0; j<_WE.m.nStacks; j++)
				  {
					if (j==0 || !_WE.synchroniseStackTregistrationErrors)
					  {
						// only add master params if stacks are synchronised
						if (_WE.seismicInfo[i][j].systematicTProp.sigma > 0)
						  {
							// variance set by  systematic sloppiness
							_WE.seismicInfo[i][j].systematicTProp.resetToMean();
						  }
					  }
				  }
			  }
		  }
	  }


	// (5) positioning error dx dy for each well
	if (haveXYseismicErrors)
	  {
		for (i=0; i<_WE.m.nWells; i++)
		  {
			// only add these if a master...
			if (_WE.wellInfo[i].synchronizeWellName == null)
			  {
				if (_WE.wellInfo[i].systematicXProp.sigma !=0 )
				  {
					// variance set by  systematic sloppiness
					var = MathUtil.sqr(_WE.wellInfo[i].systematicXProp.sigma);
					C_prior.setQuick(k,k,var);
					Cinv_prior.setQuick(k,k, 1.0/var);
					logDetC_prior_only += Math.log(var);
					// and mean is the value (will usually bezero)
					mu_prior.setQuick(k,_WE.wellInfo[i].systematicXProp.mean);
					// and this scale. Make it the expected uncertainty
					typicalScale[k] = _WE.wellInfo[i].systematicXProp.sigma;
					initialGuess[k] = _WE.wellInfo[i].systematicXProp.mean;
					k++;
				  }
				if (_WE.wellInfo[i].systematicYProp.sigma !=0 )
				  {
					// variance set by  systematic sloppiness
					var = MathUtil.sqr(_WE.wellInfo[i].systematicYProp.sigma);
					C_prior.setQuick(k,k,var);
					Cinv_prior.setQuick(k,k, 1.0/var);
					logDetC_prior_only += Math.log(var);
					// and mean is the value (will usually be zero)
					mu_prior.setQuick(k,_WE.wellInfo[i].systematicYProp.mean);
					// and this scale. Make it the expected uncertainty
					typicalScale[k] = _WE.wellInfo[i].systematicYProp.sigma;
					initialGuess[k] = _WE.wellInfo[i].systematicYProp.mean;
					k++;
				  }				
			  }
		  }
	  }
	else
	  {
		// ensure properties are reset, even if not using them, as previous models may have tweaked them
		for (i=0; i<_WE.m.nWells; i++)
		  {
			// only add these if a master...
			if (_WE.wellInfo[i].synchronizeWellName == null)
			  {
				if (_WE.wellInfo[i].systematicXProp.sigma !=0 )
				  {
					_WE.wellInfo[i].systematicXProp.resetToMean();
				  }
				if (_WE.wellInfo[i].systematicYProp.sigma !=0 )
				  {
					// variance set by  systematic sloppiness
					_WE.wellInfo[i].systematicYProp.resetToMean();
				  }				
			  }
		  }
	  }



	// (6) B coefficient in AVO equations
	if (!createWaveletCoefficientModelOnly && (_WE.m.B_AVO[0].sigma>0))
	  {
		for (j=0; j<_WE.m.nStacks; j++)
		  {
			C_prior.setQuick(k,k, var = (_WE.m.B_AVO[j].sigma * _WE.m.B_AVO[j].sigma));
			Cinv_prior.setQuick(k,k, 1.0/var);
			logDetC_prior_only += Math.log(var);
			// and mean by the actual value
			mu_prior.setQuick(k,_WE.m.B_AVO[j].mean);
			typicalScale[k] = _WE.m.B_AVO[j].sigma;
			initialGuess[k] = mu_prior.getQuick(k);
			k++;
			// only one B_AVO unless B's are all free
			if (!_WE.freeAVO) break;
		  }

		// 		if (_WE.freeAVO)
		// 		  {
		// 			C_prior.setQuick(k,k, var = (_WE.m.B_AVO_far.sigma * _WE.m.B_AVO_far.sigma));
		// 			Cinv_prior.setQuick(k,k, 1.0/var);
		// 			logDetC_prior_only += Math.log(var);
		// 			// and mean by the actual value
		// 			mu_prior.setQuick(k,_WE.m.B_AVO_far.mean);
		// 			typicalScale[k] = _WE.m.B_AVO_far.sigma;
		// 			initialGuess[k] = mu_prior.getQuick(k);
		// 			k++;
		// 		  }
	  }
	else
	  {
		// ensure properties are reset, even if not using them, as previous models may have tweaked them
		for (j=0; j<_WE.m.nStacks; j++)
		  {
			_WE.m.B_AVO[j].resetToMean();
			if (!_WE.freeAVO) break;
		  }
	  }


	// (6.5) Anisotropy delta term
	if (hasAnisotropy) // _WE.doAnisotropicModel)
	  {
		if (!createWaveletCoefficientModelOnly && (_WE.m.deltaThomsen_shale.sigma>0))
		  {
			C_prior.setQuick(k,k, var = (_WE.m.deltaThomsen_shale.sigma * _WE.m.deltaThomsen_shale.sigma));
			Cinv_prior.setQuick(k,k, 1.0/var);
			logDetC_prior_only += Math.log(var);
			// and mean by the actual value
			mu_prior.setQuick(k,_WE.m.deltaThomsen_shale.mean);
			typicalScale[k] = _WE.m.deltaThomsen_shale.sigma;
			initialGuess[k] = mu_prior.getQuick(k);
			k++;
		  }
		else
		  {
			// ensure properties are reset, even if not using them, as previous models may have tweaked them
			_WE.m.deltaThomsen_shale.resetToMean();
		  }
	  }

	// (7) Far wavelet stretching
	if (!createWaveletCoefficientModelOnly && (_WE.m.nStacks>1) && _WE.useStretchedModel)
	  {
		for (j=0; j<(_WE.m.nStacks-1); j++)
		  {
			if (_WE.m.farStackWaveletTimeStretch[j].sigma>0)
			  {
				C_prior.setQuick(k,k, var = sqr(_WE.m.farStackWaveletTimeStretch[j].sigma));
				Cinv_prior.setQuick(k,k, 1.0/var);
				logDetC_prior_only += Math.log(var);
				// and mean by the actual value
				mu_prior.setQuick(k,_WE.m.farStackWaveletTimeStretch[j].mean);
				typicalScale[k] = 1.0; // _WE.m.farStackWaveletTimeStretch[j].sigma;
				initialGuess[k] = mu_prior.getQuick(k);
				k++;
			  }
			if (_WE.m.farStackWaveletAmplitude[j].sigma>0)
			  {
				C_prior.setQuick(k,k, var = sqr(_WE.m.farStackWaveletAmplitude[j].sigma));
				Cinv_prior.setQuick(k,k, 1.0/var);
				logDetC_prior_only += Math.log(var);
				// and mean by the actual value
				mu_prior.setQuick(k,_WE.m.farStackWaveletAmplitude[j].mean);
				typicalScale[k] = _WE.m.farStackWaveletAmplitude[j].sigma;
				initialGuess[k] = mu_prior.getQuick(k);
				k++;
			  }
		  }
	  }
	else
	  {
		if (_WE.useStretchedModel)
		  {
			for (j=0; j<(_WE.m.nStacks-1); j++)
			  {
				// ensure properties are reset, even if not using them, as previous models may have tweaked them
				if (_WE.m.farStackWaveletTimeStretch[j].sigma>0)
				  {
					_WE.m.farStackWaveletTimeStretch[j].resetToMean();
				  }
				if (_WE.m.farStackWaveletAmplitude[j].sigma>0)
				  {
					_WE.m.farStackWaveletAmplitude[j].resetToMean();
				  }
			  }
		  }
	  }

	// (7.5) Saturation effects
	if (_WE.useSaturationModel)
	  {
		if (!createWaveletCoefficientModelOnly)
		  {
			for (i=0; i<_WE.m.nWells; i++)
			  {
				SaturationInfo _sInfo = _WE.wellInfo[i].saturationInfo;
				if (_sInfo.satProp.sigma > 0 )
				  {
					C_prior.setQuick(k,k, var = sqr(_sInfo.satProp.sigma));
					Cinv_prior.setQuick(k,k, 1.0/var);
					logDetC_prior_only += Math.log(var);
					// and mean by the actual value
					mu_prior.setQuick(k,tmp = _sInfo.satProp.mean);
					typicalScale[k] = 1.0;
					initialGuess[k] = tmp;
					k++;
				  }
			  }
		  }
		else
		  {		
			// ensure properties are reset, even if not using them, as previous models may have tweaked them	
			for (i=0; i<_WE.m.nWells; i++)
			  {
				SaturationInfo _sInfo = _WE.wellInfo[i].saturationInfo;
				if (_sInfo.satProp.sigma > 0 ) _sInfo.satProp.resetToMean();
			  }
		  }
	  }

	// (8) Stack noise levels....
	// these have conjugate prior (say, Jeffrey's), so are not part of Gaussian model
	for (i=0; i<_WE.m.nStacks; i++,k++)
	  {
		// set scale to 1/4 of peak amplitude
		typicalAmp = _WE.multiModelSetUp.typicalRMSSeismicAmplitude;
		if (fixNoise && (_WE.fixNoiseValue > 0))
		  {
			var = _WE.fixNoiseValue;
		  }
		else
		  {
			// start initial guess too large to save optimiser going ballistic in some cases
			//var = 3*TYP_NOISE_TO_AMP_RATIO*typicalAmp;
			var = 2*typicalAmp;
		  }
		// Useful gedanken test of model comparison formalism..
		// set var to be much greater than a typical amplitude 
		if (createWaveletCoefficientModelOnly)
		  {
			// Here we are starting a wavelet-only optimisation. 
			// there are no initialGuess or typicalScale fields for noise
			// But need to set in-situ noise to a decent value...
			_WE.m.seismic_noise[i].setValue(_WE.ignoreSeismic ? 1.0 : var);
			// set the mean too, as this is useful
			_WE.m.seismic_noise[i].setMean(_WE.ignoreSeismic ? 1.0 : var);
		  }
		else
		  {
			if (!fixNoise)
			  {
				typicalScale[k] = var;
				initialGuess[k] = _WE.ignoreSeismic ? 1.0 : var;
				// set the mean too, as this is useful
				_WE.m.seismic_noise[i].setMean(_WE.ignoreSeismic ? 1.0 : var);
			  }
		  }
	  }


	// Set up NonSynthetic constraints
	initialiseNonSyntheticConstraints();
  }


  /**
	 useful method to reset the initialGuess vector to the values currently in the fields
  */
  void recomputeInitialGuessFromInSituFields()
  {
	int i;
	double tOffset=0.0;
	Property prop;

	if (_WE.constrainPeakArrival && haveTseismicErrors && _WE.peakArrivalStdDev < 0.1) 
	  {
		// We have time registration variables available. Set these 
		// from the wavelet peak. Do this only 
		// if the likelihood tolerance seems stringent enough (we use 100ms)
		tOffset = wavelet.getPeakTime(false,constrainPeakArrivalSign,peakArrivalPositive);
		// Shift the wavelet so the initial guess is close to peak time, 
		// Can possibly restrict shift to case if there
		// seems to be little ambiguity about which peak is the correct one, e.g.
		// if (!wavelet.isSecondPeakStrong())
		if (true)
		  {
			if (_WE.verbosity>3)
			  {
				_WE.EH.println("Peak arrival time from wavelet applied to registration "
							   +"shifts before full optimisation = %10.4g",tOffset);
			  }
			wavelet.embedTimeShift(tOffset);
		  }
		else 
		  {
			// make offset embedding benign below
			tOffset=0;
		  }
	  }

	
	for (i=0; i<nParams; i++)
	  {
		prop = (Property)(mapFromIndex.get(i));
		// String name = (String)mapToName.get(prop);
		if ((i>=firstTregistrationParamIndex) && (i<=lastTregistrationParamIndex))
		  {
			// Set time registration variables available
			// registration properties
			initialGuess[i] = prop.value + tOffset;
		  }
		else
		  {
			initialGuess[i] = prop.value;
		  }
	  }
  }




  // computes prior negative log likelihood of currrent model
  double computePriorNegLogLikelihood(boolean verbose)
  {
	int i,j;
	Property prop;
	double JeffreyTerm=0.0,chisq;
	DoubleMatrix1D delta_m,ans;

	// (1) First , the Gaussian pieces

 	delta_m = new DenseDoubleMatrix1D(mu_prior.size());
	for (j=0;j<mu_prior.size();j++) 
	  {
		prop = (Property)mapFromIndex.get(j);
		delta_m.setQuick(j,prop.value - mu_prior.getQuick(j));
	  }

	ans = delta_m.like();

	// form C^1.delta_m -> ends up in-place in "ans"
	Cinv_prior.zMult(delta_m,ans);

	// then delta_m^T.C^-1.delta_m/2.0
	chisq = delta_m.zDotProduct(ans)/2.0;

	ans = null; // GC
	delta_m = null; //GC

	if (verbose) _WE.EH.println("Prior: contrib. from exponent= "+chisq);

    // add determinant piece too
	chisq += logDetC_prior_only/2.0;
	if (verbose) _WE.EH.println("Prior: contrib. from det= "+logDetC_prior_only/2.0);

	// (2) now the non-Gaussian pieces
	// First Guess. Jeffrey's prior
	if (!_WE.fixNoise && !ModeFinder.ABSORB_JEFFREYS_PRIOR_INTO_LIKELIHOOD)
	  {
		for (i=0; i<_WE.m.nStacks; i++)
		  {
			// make sure -ve sigma's do not wreak havoc
			JeffreyTerm += Math.log(_WE.m.seismic_noise[i].value*_WE.m.seismic_noise[i].value);
		  }
		chisq += JeffreyTerm;
		if (verbose) _WE.EH.println("Prior: contrib. from Jeffrey term= "+JeffreyTerm);
	  }

	return chisq;
  }

  void initialiseNonSyntheticConstraints()
  {
	int k;

	// dimension count

	dimConstraints = 0;

	if (_WE.constrainIntervalVelocities)
	  {
		if (!_WE.fixCheckshots)
		  {
			setUpSonicLogIntervalVelocityConstraints();
			// (1) checkshot velocities
			for (k=0;k<_WE.m.nWells;k++)
			  {
				dimConstraints += (_WE.m.timeToDepthMaps[k].t.length-1);
			  }
		  }
	  }
	if (_WE.constrainPhase)
	  {
		// Constant phase constraints
		phaseConstraintIndxStart = (int)Math.rint((0.3*_WE.multiModelSetUp.f_peak/wavelet.phaseVectorFrequencySpacing));
		phaseConstraintIndxEnd = (int)Math.rint((1.9*_WE.multiModelSetUp.f_peak/wavelet.phaseVectorFrequencySpacing));
		phaseConstraintIndxEnd = Math.min(phaseConstraintIndxEnd, wavelet.phase.length-1);
		if (_WE.verbosity>3)
		  {
			_WE.EH.println("Phase constraint index range: "+phaseConstraintIndxStart+" "+phaseConstraintIndxEnd);
		  }
		if (!CONSTRAIN_PHASE_VIA_PHASE_LINE)
		  {
			// twice as many extra constraints for phase via spectral ratios
			dimConstraints += 2*((phaseConstraintIndxEnd - phaseConstraintIndxStart)+1) ;
		  }
		else
		  {
			dimConstraints += (phaseConstraintIndxEnd - phaseConstraintIndxStart)+1 ;
		  }
	  }
	if (_WE.constrainPeakArrival || _WE.constrainPeakArrivalSign) 
	  {
		// Peak arrival constraint. Only one term here
		dimConstraints ++;
	  }

	// Add more as you please....
	// ...

	// finally, allocate error space vector
	constraintErrors = new double[dimConstraints];
  }

  void setUpSonicLogIntervalVelocityConstraints()
  {
	int i=0,k,j,iZStart,maxSegIndx;
	TimeToDepthMap timeDepth=null;
	double[] segDepths=null;    // layer tops in true depth (from segmented data)
	double[] segVelocity=null;  // layer velocities (from segmented data)
	double sum=0.0,norm=0.0,wt,dz;
	boolean intervalIsComplete = true;

	for (k=0; k<_WE.m.nWells; k++)
	  {
		timeDepth = _WE.m.timeToDepthMaps[k];
		segDepths = _WE.sLogData[k].Z;
		segVelocity = _WE.sLogData[k].VpHarm; // use slowness average here
		maxSegIndx = segDepths.length-1;

		if (_WE.verbosity>4)
		  {
			_WE.EH.println("VC: Sonic log velocity constraints, well " +_WE.wellInfo[k].logFilename);
			_WE.EH.println("VC: %-20s %-20s %-20s","Depth_interval_start",
						   "Vp_harmonic_mean",
						   "initial_checkshot_scaling_velocity");
		  }
		// loop over intervals in t2d map for this well
		for (i=0; i<timeDepth.t.length-1; i++)
		  {
			// Try to find range of segmented indices that straddle this interval
			// we assume we can compute an average velocity from logs until proven otherwise
			intervalIsComplete = true;
			iZStart = Model.getStraddlingLowerIndex(timeDepth.Z[i],segDepths);
			if (iZStart>=0 && iZStart<(maxSegIndx-1))
			  {
				// the left edge is bounded OK && there is at least one segmented boundary at larger Z
				wt = Math.min(segDepths[iZStart+1],timeDepth.Z[i+1]) - timeDepth.Z[i];
				// start of harmonic mean
				sum = wt/segVelocity[iZStart];
				norm = wt;
				if (segDepths[iZStart+1] < timeDepth.Z[i+1])
				  {
					// more than one segment in this checkshot interval
					for (j=iZStart+1;;j++)
					  {
						// loop over segments within interval
						if (j>=maxSegIndx)
						  {
							intervalIsComplete = false;
							break;
						  }
						else
						  {
							if (segDepths[j+1]<=timeDepth.Z[i+1])
							  {
								// segment is entirely within interval
								sum += (segDepths[j+1] - segDepths[j])/segVelocity[j];
								norm += (segDepths[j+1] - segDepths[j]);
							  }
							else
							  {
								// interval finishes in this segment
								dz = timeDepth.Z[i+1] - Math.max(segDepths[j],timeDepth.Z[i]);
								sum += dz/segVelocity[j];
								norm += dz;
								break;
							  }
						  }
					  }
				  }
			  }
			else
			  {
				intervalIsComplete = false;
			  }
			if (intervalIsComplete)
			  {
				// harmonic mean
				double vBar = 1.0/(sum/norm);
				timeDepth.upscaledLogVelocity[i] = vBar;
				if (WaveletExtractor.AUTO_DETECT_VELOCITY_UNITS)
				  {
					// suspect vBar and checkshot velocities in different units 
					// (e.g former in in ft/s, checkshot in m/s)
					timeDepth.upscaledLogVelocity[i] *= _WE.sLogData[k].lengthUnitsConversionFactorCheckshot;
				  }
			  }

			if (_WE.verbosity>4)
			  {
				_WE.EH.printf("VC: %-20g ",timeDepth.Z[i]);
				if (intervalIsComplete) 
				  {
					_WE.EH.printf("%-20g ",timeDepth.upscaledLogVelocity[i]);
				  }
				else _WE.EH.printf("%-20s ","unformable");
				_WE.EH.println("%-20g",timeDepth.scalingVelocity[i]);
			  }
		  }
	  }
  }




  double computeAuxiliaryConstraints()
  {
	int i,j,k;
	double chisq=0.0;

	j=0; // runs over all constraints......
	Arrays.fill(constraintErrors,0.0);

	// (1) Checkshot velocity penalties
	if (_WE.constrainIntervalVelocities && !_WE.fixCheckshots) 
	  {
		double vP;
		TimeToDepthMap timeDepth=null;
		chisqVelMismatch = 0;
		for (k=0; k<_WE.m.nWells; k++)
		  {
			timeDepth = _WE.m.timeToDepthMaps[k];
			for (i=1; i<timeDepth.t.length; i++,j++)
			  {
				// current model interval velocity (remember timeDepth.t's are two-way times)
				vP = 2.0*(timeDepth.Z[i]-timeDepth.Z[i-1])/(timeDepth.t[i]-timeDepth.t[i-1]);
				// recall: scalingVelocity[i-1] = 2*(Z[i]-Z[i-1])/(t[i]-t[i-1]);
				if (vP>0)
				  {
					// define errors as (model - data)/sigma, same as seismic error
					constraintErrors[j] = (vP - timeDepth.upscaledLogVelocity[i-1])/
					  (_WE.intervalVelocityConstraintStdDev * timeDepth.upscaledLogVelocity[i-1]);
					// A test: constrain slownesses..
					//constraintErrors[j] = (1.0/vP - 1.0/timeDepth.upscaledLogVelocity[i-1])/
					//        (_WE.intervalVelocityConstraintStdDev * (1.0/timeDepth.upscaledLogVelocity[i-1]));
					// _WE.EH.println("VpInt error = "+constraintErrors[j]);
				  }
				else
				  {
					constraintErrors[j]=10000; // FIXME PROMOTE ???
				  }
				// accumulate error
				chisqVelMismatch += constraintErrors[j]*constraintErrors[j];
			  }
			chisq += chisqVelMismatch;
		  }
	  }

	// add other stuff here if desired.......

	// (2) Constant phase constraints	
	if (_WE.constrainPhase && !_WE.ignoreSeismic)
	  {
		// if ignoring seismic, wavelet will be zeros, so no contribution here

		chiSqPhase=0.0;
		wavelet.computePhaseSpectrum(null,0,0);
		// constrain phase to desired number
		if (CONSTRAIN_PHASE_VIA_PHASE_LINE)
		  {
			// Here we constrain the phase spectrum points at each 
			// frequency of interest to lie "along"
			// a line of the desired slope in the complex plane. 
			// This avoid the branch cut involved when working directly with angles
			// Phase flips of 180 degrees associated with these cuts are just 
			// sign flips anyway.
			double targetSlope=0,targetAngle;
			if (Math.abs(_WE.constantPhase) > 360)
			  {
				// constant phase angle is free and "floating"
				for (i=phaseConstraintIndxStart,targetAngle=0.0; i<=phaseConstraintIndxEnd; i++) 
				  {
					targetAngle +=wavelet.phase[i];
				  }
				targetAngle /= (phaseConstraintIndxEnd - phaseConstraintIndxStart +1);
			  }
			else
			  {
				// user-specified angle
				targetAngle = _WE.constantPhase;
			  }
			// don't constrain to better than 2 degrees (unstable)
			// typical line uncertainty must be sigma_theta*typical_wavelet_coeff.
			double coeffStdDev = Math.toRadians(Math.max(2,_WE.constantPhaseStdDev))*_WE.multiModelSetUp.typicalRMSSeismicAmplitude/FundamentalConstants.TYPICAL_REFLECTION_COEFF;
			if (!true)
			  {
				if (Math.abs(Math.toRadians(targetAngle)) < Math.PI/4.0)
				  {
					targetSlope = Math.tan(Math.toRadians(targetAngle));
					for (i=phaseConstraintIndxStart; i<=phaseConstraintIndxEnd; i++,j++)
					  {
						// define errors as (model - data)/sigma, same as seismic error
						constraintErrors[j] = (wavelet.phaseIm[i] - wavelet.phaseRe[i] * targetSlope)/coeffStdDev;
						chiSqPhase += constraintErrors[j]*constraintErrors[j];
					  } 
				  }
				else
				  {
					// targetSlope = cot(angle)
					targetSlope = Math.tan(Math.toRadians(targetAngle));
					targetSlope = (targetSlope != 0.0) ? 1.0/targetSlope : Double.MAX_VALUE; 
					for (i=phaseConstraintIndxStart; i<=phaseConstraintIndxEnd; i++,j++)
					  {
						// define errors as (model - data)/sigma, same as seismic error
						// define errors as (model - data)/sigma, same as seismic error
						constraintErrors[j] = (wavelet.phaseRe[i] - wavelet.phaseIm[i] * targetSlope)/coeffStdDev;
						chiSqPhase += constraintErrors[j]*constraintErrors[j];
					  }
				  }
			  }
			else
			  {
				// form contours about target ray with penalties for polarity flips.
				// Think hippodrome-like contours about the target ray
				double dotProd,xProd,cosTheta,sinTheta;
				cosTheta = Math.cos(Math.toRadians(targetAngle));
				sinTheta = Math.sin(Math.toRadians(targetAngle));
				for (i=phaseConstraintIndxStart; i<=phaseConstraintIndxEnd; i++,j++)
				  {
					dotProd = cosTheta * wavelet.phaseRe[i] + sinTheta * wavelet.phaseIm[i];
					if (dotProd < 0) // point is "behind" target ray
					  {
						constraintErrors[j] = Math.sqrt((wavelet.phaseRe[i]*wavelet.phaseRe[i] + 
														 wavelet.phaseIm[i]*wavelet.phaseIm[i]))/coeffStdDev;
					  }
					else
					  {
						xProd = sinTheta * wavelet.phaseRe[i] - cosTheta * wavelet.phaseIm[i];
						// define errors as (model - data)/sigma, same as seismic error
						constraintErrors[j] = Math.abs(xProd)/coeffStdDev;
					  }

					chiSqPhase += constraintErrors[j]*constraintErrors[j];
				  } 
			  }
		  }
		else
		  {
			// !CONSTRAIN_PHASE_VIA_PHASE_LINE
			// Work directly on angles. Try to avoid branch discontinuities this way...
			double targetCos,targetSin;
			double[] cosPhase,sinPhase;
			int nPC = phaseConstraintIndxEnd - phaseConstraintIndxStart+1;
			cosPhase = new double[nPC];
			sinPhase = cosPhase.clone();
			targetCos = targetSin = 0;
			for (i=phaseConstraintIndxStart,k=0; i<=phaseConstraintIndxEnd; i++,k++) 
			  {
				cosPhase[k] = Math.cos(Math.toRadians(wavelet.phase[i]));
				sinPhase[k] = Math.sin(Math.toRadians(wavelet.phase[i]));
				targetCos += cosPhase[k];
				targetSin += sinPhase[k];
			  }
			targetCos /= nPC;
			targetSin /= nPC;
			// override these targets if angle is specified
			if (Math.abs(_WE.constantPhase) < 360)
			  {
				// user-specified angle
				targetCos = Math.cos(Math.toRadians(_WE.constantPhase));
				targetSin = Math.sin(Math.toRadians(_WE.constantPhase));
			  }
			// don't constrain to better than 2 degrees (unstable) (0.035 = sin(2degrees))
			double phaseStdDev = Math.max(0.035,Math.toRadians(_WE.constantPhaseStdDev));
			for (k=0;k<nPC;j++,k++)
			  {
				// define errors as (model - data)/sigma, same as seismic error
				constraintErrors[j] = (cosPhase[k] - targetCos)/phaseStdDev;
				chiSqPhase += constraintErrors[j]*constraintErrors[j];
				j++;
				constraintErrors[j] = (sinPhase[k] - targetSin)/phaseStdDev;
				chiSqPhase += constraintErrors[j]*constraintErrors[j];
			  }
			sinPhase = cosPhase = null;
		  }
		chisq += chiSqPhase;
	  }

	// (3) Peak arrival
	if (_WE.constrainPeakArrival && !_WE.ignoreSeismic)
	  {
		// if ignoring seismic, wavelet will be zeros, so no contribution here
	    chisqPeakMismatch=0;
		// switch on peak arrivals only for wavelets of appreciable magnitude
		double ampThresh = 0.1*_WE.multiModelSetUp.typicalRMSSeismicAmplitude/
		  FundamentalConstants.TYPICAL_REFLECTION_COEFF;
		double peak = wavelet.getCoarseAbsPeak();
		if (true)
		  {

			if (!ModeFinder.HANDLE_AMBIGUOUS_PEAK_DETECTION_TERMINATION_USING_PENALTY
				&& useNearestPeakArrivalConstraint)
			  {
				// (1) Use Greedy peak detection to latch onto nearest of the 2 largest peaks, 
				// not purely the largest: Nearest of the 2 largest
				constraintErrors[j] = (wavelet.getNearestPeakTime(_WE.peakArrivalTime, constrainPeakArrivalSign, peakArrivalPositive) 
									   - _WE.peakArrivalTime)/_WE.peakArrivalStdDev;
			  }
			else
			  {
				// Acts if HANDLE_AMBIGUOUS_PEAK_DETECTION_TERMINATION_USING_PENALTY is on (default)
				// (2) Snap to largest, with weak arbitration for competing peaks...
				constraintErrors[j] = (wavelet.getPeakTime(false, constrainPeakArrivalSign, peakArrivalPositive) 
									   - _WE.peakArrivalTime)/_WE.peakArrivalStdDev;
			  }
		  }
		else
		  {
			constraintErrors[j]=0;
		  }
		// Ramp this up to threshold level
		if (wavelet.getCoarseAbsPeak() < ampThresh)
		  {
			constraintErrors[j] *= (peak/ampThresh); 
		  }


		chisqPeakMismatch += constraintErrors[j]*constraintErrors[j];
		j++;

		chisq += chisqPeakMismatch;
	  }

	// (4) Far wavelet stretch 
	if (_WE.m.farStackWaveletTimeStretchMapper!=null) 
	  {
		double stretch;
		for (j=0;j<(_WE.m.nStacks-1);j++)
		  {
			stretch = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(_WE.m.farStackWaveletTimeStretch[j].value);
			chisq += MathUtil.sqr((stretch - _WE.m.farStackWaveletTimeStretchMean)/_WE.m.farStackWaveletTimeStretchSigma);
		  }
	  }
	
	// make chisq a conventional 0.5-normalised value
	chisq/=2.0;


	return chisq;
  }




  /**
	 ASCII dump of most likely parameters	 
  */
  void dumpMLparameters(String filename)
  {
	int i,ii,j,k,nG;
	double MLval, priorVal, priorStdDev, deviation, postStdDev;
	String name;
	SaturationInfo _sInfo;
	Property prop;

	PrintWriter out;

	// open ascii summary file
	_WE.EH.println("Writing summary of ML parameters to "+filename);
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }

	// Number of Posterior Gaussian parameters
	nG = C_prior.rows();


	out.printf("List of parameters and deviations. ( Wavelet coarse dt= %g f1= %g )%n",
			   _WE.m.wavelet.dtCoarse,_WE.m.wavelet.tCoarse[1]);
	out.printf("%-45s %-10s %-10s %-12s %-30s%n","property","prior-mean",
			   "ML-value","sigma_post","(<post>-<prior>)/sigma_prior");
		
	for (i=0; i<nG; i++)
	  {
		MLval = mu_post.getQuick(i);
		priorVal = mu_prior.getQuick(i);
		prop = (Property)mapFromIndex.get(i);
		priorStdDev =  Math.sqrt(C_prior.getQuick(i,i));
		name = 	(String)mapToName.get(mapFromIndex.get(i));
		postStdDev = Math.sqrt(CtildeDiag[i]);
		// some properties may have been remapped. Retrieve true values here...
		// (1) Wavelet stretch factors
		if (_WE.m.farStackWaveletTimeStretchList.contains(prop) &&
			_WE.m.farStackWaveletTimeStretchMapper!=null) 
		  {
			MLval = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(MLval);
			// old -> _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(priorVal);
			// untransformed std dev is common to all stretch parameters
			priorVal = _WE.m.farStackWaveletTimeStretchMean; 
			double slope = _WE.m.farStackWaveletTimeStretchMapper.getFDSlope(MLval);
			// approximate posterior std deviation
			postStdDev = postStdDev/slope;
			// new
			priorStdDev = _WE.m.farStackWaveletTimeStretchSigma;
		  }
		// (2) Saturations
		if (_WE.useSaturationModel)
		  {
			// lookup for saturations. Only variable ones appear, so they will have a mapper			
			if ((_sInfo = getParentSatInfo(prop)) != null)
			  {
				MLval = _sInfo.satPropMapper.getTrueValue(prop.value); // will be in (0,1)
				priorVal = _sInfo.Sa;
				priorStdDev = _sInfo.sig_Sa;
				// approximate equivalent std deviation of logitised variable
				postStdDev = sqrt(sqr(postStdDev)*(MLval*(1-MLval)));  
			  }
		  }	

		deviation = (priorStdDev>0) ? (MLval-priorVal)/priorStdDev : 0; 

		out.printf("%-45s %10.6g %10.6g %12.6g %30g%n",
				   name, priorVal, MLval,
				   postStdDev ,deviation);
	  }
	// add the noise bits in here too
	for (ii=0;ii<_WE.m.nStacks;ii++)
	  {
		if (!_WE.fixNoise)
		  {
			i = ((Integer)mapToIndex.get(_WE.m.seismic_noise[ii])).intValue();
			MLval = Math.exp(mu_post.getQuick(i));
			// approximate std deviation from linearisation
			postStdDev = MLval*Math.sqrt(CtildeDiag[i]);
			name = 	(String)mapToName.get(mapFromIndex.get(i));
		  } 
		else 
		  {
			name = ("noise_stack_"+ii);
			MLval = _WE.m.seismic_noise[ii].value;
			postStdDev=0.0;
		  }
		out.printf("%-45s %10s %10.6g %12.6g %30s%n",
				   name,"-", MLval,
				   postStdDev,"-");
	  }
	// And the noise corr. coeff matrix
	if (noiseCorrCoeffMatrix != null)
	  {
		for (i=0;i<_WE.m.nStacks;i++)
		  {
			for (ii=i+1;ii<_WE.m.nStacks;ii++)
			  {
				out.printf("%-45s %10s %10.6g %12s %30s%n",
						   "noise_corr_coeff_"+(i)+"_"+(ii),"-",noiseCorrCoeffMatrix[i][ii],
						   "-","-");
			  }
		  }
	  }
	
	out.close();
	_WE.EH.println("Wrote " + filename);
  }



  SaturationInfo getParentSatInfo(Property prop)
  {
	// lookup for saturations. Only variable ones appear, so they will have a mapper			
	for(int k=0;k<_WE.m.nWells;k++)
	  {
		SaturationInfo _sInfo = _WE.wellInfo[k].saturationInfo;
		if(_sInfo==null || _sInfo.satProp==null || _sInfo.satPropMapper==null) continue;
		if (prop == _sInfo.satProp) 
		  {
			return _sInfo;
		  }
	  }
	return null;
  }

  /**
	 ASCII dump of most likely time to depth map 
  */
  void dumpMLTimeToDepthMap(String filename)
  {
	int i,j,k;
	double zMD, MLval, MLoffset_val, priorVal, priorStdDev, deviation, postStdDev;
	String name;
	Property prop;
	boolean writeSeismicMinicubeTime = true;

	PrintWriter out;
	

	for (i=0; i<_WE.m.nWells; i++)
	  {

		// open ascii summary file
		String fullFilename = filename + "."+Util.fileBasename(_WE.wellInfo[i].logFilename);
		fullFilename = fullFilename.replaceAll("(?i)\\.las|\\.txt","") + ".txt";

		_WE.EH.println("Writing summary of ML parameters to "+fullFilename.replaceAll("%","%%"));
		try 
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(fullFilename)));
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + fullFilename);
		  }
		out.printf("Most likely time to depth map, with errors. ");

		// Interpolation for TVD and MD columns
		// TVD to MD mapping from survey. MUST BE MONOTONIC
		LinearInterpolation zToMDInterp = 
		  new LinearInterpolation(_WE.logData[i].survey.Z, _WE.logData[i].survey.Zm, _WE.logData[i].survey.Z.length, _WE.EH);

		if (!_WE.m.registrationErrorsAreCommon)
		  {
			// 2 different offsets. Can't absorb both into an "effective" t2d conversion
			writeSeismicMinicubeTime = false;
			out.printf("No seismic-minicube-time column as stacks have different registration offsets!");
		  }
		out.printf("%n");

		if (writeSeismicMinicubeTime)
		  {
			out.printf("%-10s %-10s %-22s %-18s %-10s%n","TVD","MD","t_seismic","t_checkshot","sigma-time");		
		  }
		else
		  {
			out.printf("%-10s %-10s %-18s %-10s%n","TVD","MD","t_checkshot","sigma-time"); // 
		  }

		TimeToDepthMap t2dmap = _WE.m.timeToDepthMaps[i];
		for (j=0;j<t2dmap.tProp.length;j++)
		  {
			// was probabilistic
			Integer Indx = ((Integer)mapToIndex.get(t2dmap.tProp[j]));
			if (Indx != null)
			  {
				k = Indx.intValue();
				MLval = mu_post.getQuick(k);
				postStdDev = Math.sqrt(CtildeDiag[k]);
			  }
			else
			  {
				// fixed
				MLval = t2dmap.tProp[j].value;
				postStdDev=0;
			  }

			// get MD
			zMD = zToMDInterp.interpolate(t2dmap.Z[j],j!=0);

			// priorStdDev =  Math.sqrt(C_prior.getQuick(k,k));
			// name = 	(String)mapToName.get(mapFromIndex.get(k));
			// deviation = (priorStdDev>0) ? (MLval-priorVal)/priorStdDev : 0;					
			if (writeSeismicMinicubeTime)
			  {
				Indx = ((Integer)mapToIndex.get(t2dmap.systematicTProp));
				if (Indx != null)
				  {
					MLoffset_val = mu_post.getQuick(Indx.intValue());
				  }
				else
				  {
					// fixed
					MLoffset_val = t2dmap.systematicTProp.value;
				  }
				// GET ML minicube offset value
				// MLoffset_val = t2dmap.systematicTProp.value;
				out.printf("%-10.6g %-10.6g %-22.6g %-18.6g %-10.6g%n",
						   t2dmap.Z[j], zMD, MLval + MLoffset_val, MLval , postStdDev);
			  }
			else
			  {
				out.printf("%10.6g %-10.6g %18.6g %10.6g%n",
						   t2dmap.Z[j], zMD, MLval, postStdDev);
			  }
		  }
	
		out.close();
		_WE.EH.println("Wrote data to: " + fullFilename.replaceAll("%","%%"));
	  }
  }

  // experimental guff
  void getInversionPredictionWavelet_exptl()
  {
	if (_WE.m.nWells!=1) return; // can only do 1 well, near, Fatti :-)

	int i,j,k,kk,kTD,n,index,L,nd,nk,nd3=0,nm1,its,MAXITS;
	double t,tOffset=0, dt,dtSub,zTVD,zMD,tFinish,rho,DT_scaling,vp_scaling,Sscale=1;
	double DT_Vp_sonic,lam,a,b,c,x,v,varS,varA,E=0,E0,Ec=0,u,w,grad,lam_NP,SS=0,S2=0;
	double stdDev;
	boolean haveLogData;
	double logDataZmin, logDataZmax;
	LinearInterpolation[] sTrueInterp;
	double[] _Zm;

	_WE.EH.printf("%nInversion prediction experiments....%n");

	// compute the fwd model before extracting AI values
	computeMLForwardModel();

	kTD = firstTDparamIndex;
	for (i=0; i<_WE.m.nWells; i++)
	  {
		LogData LDI = _WE.logData[i];

		DT_scaling = _WE.logData[i].getSonicLogScaling();
		vp_scaling = _WE.sLogData[i].lengthUnitsConversionFactor;
		// arrays for seismic interpolation
		n = _WE.convModel.nZ;
		double[] times = new double[n];
		dt = _WE.m.wavelet.delta_t;
		for (t=_WE.convModel.tStartReff,k=0; k<n; k++,t+=dt) times[k] = t;
		tFinish = times[n-1];
		// array of interpolations for stacks
		sTrueInterp = new LinearInterpolation[_WE.m.nStacks];

		// time sampling: use standard dt
		dtSub = _WE.m.wavelet.delta_t;
			
		TimeToDepthMap t2dMap = _WE.m.timeToDepthMaps[i];
		// time to TVD mapping 
		LinearInterpolation zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);

		// TVD to MD mapping from survey. MUST BE MONOTONIC
		LinearInterpolation zToMDInterp = 
		  new LinearInterpolation(_WE.logData[i].survey.Z, _WE.logData[i].survey.Zm, _WE.logData[i].survey.Z.length, _WE.EH);
		// but don't ask for log data outside this range....
		_Zm = _WE.logData[i].Zm; // shortcut

		logDataZmin = _Zm[0];
		logDataZmax = _Zm[_Zm.length - 1];

		// density
		// raw..from measured depths			
		LinearInterpolation rhoInterp = new LinearInterpolation(_Zm,
																_WE.logData[i].rho,_Zm.length,_WE.EH);
		// DT log
		LinearInterpolation sonicInterp = new LinearInterpolation(_Zm,
																  _WE.logData[i].DT,_Zm.length,_WE.EH);

		for (j=0; j<_WE.m.nStacks; j++) 
		  {
			// seismic interp for jth stack
			sTrueInterp[j] = new LinearInterpolation(times, _WE.convModel.S_true[i][j], n ,_WE.EH);
		  }

		n = (int)floor((_WE.tEnd -_WE.tStart)/dtSub);
		nm1 = n-1;
		double[] AI = new double[n];
		double[] S = new double[n];
		double[] mu = new double[n];
		// run over fine-scale times
		boolean missing=false;
		kk=0;
		for (t=_WE.tStart,k=0,varS=0; k<n ; k++,t+=dtSub)
		  {
			// no-reflection-outside logs option can extend the range of t  a little  so the
			// depth interpolations don't work. Skip these times
			tOffset = t;
			if (_WE.m.registrationErrorsAreCommon) tOffset += t2dMap.systematicTProp.value;

			if (t<_WE.tStart || t>_WE.tEnd) continue;
			// get TVD
			zTVD = zInterp.interpolate(t,k!=0);
			// get MD
			zMD = zToMDInterp.interpolate(zTVD,k!=0);
			// see if log data is there
			haveLogData = (zMD >=logDataZmin && zMD <= logDataZmax);

			// get log data
			if (haveLogData)
			  {
				if(false)
				  {
					// no upscaling...
					rho = rhoInterp.interpolate(zMD,k!=0);
					// sonic
					DT_Vp_sonic = sonicInterp.interpolate(zMD,k!=0);
					AI[k] = Math.log(rho* vp_scaling * 1e6/DT_Vp_sonic);
				  }
				else
				  {
					// full Backus upscaling
					double tp,rhop,DT_Vp_sonicp,dtf=dtSub/2.0;
					for(tp=max(t-dtf,_WE.tStart),rho=0,u=0,j=0; tp<=min(t+dtf,_WE.tEnd); tp+=dtSub/20.0)
					  {
						// get TVD
						zTVD = zInterp.interpolate(tp,false);
						// get MD
						zMD = zToMDInterp.interpolate(zTVD,k!=0);
						if(zMD >=logDataZmin && zMD <= logDataZmax)
						  {
							rhop = rhoInterp.interpolate(zMD,k!=0);
							// sonic
							DT_Vp_sonicp = sonicInterp.interpolate(zMD,k!=0);
							rho += rhop;
							u += 1.0/(rhop*sqr(1e6/DT_Vp_sonicp));	// harmonic av of Mp = rho.Vp^2
							j++;
						  }
					  }
					rho /=j;
					u = j/u;	       // u <- harmonic average modulus Mp			  
					u = sqrt(u/rho);   // effective Vp, in log units
					u *= vp_scaling;   // to m/s
					AI[k] = Math.log(rho * u); // log( rho Vp_Backus )
				  }
			  }  
			else
			  {
				_WE.EH.printf("Missing log data at MD=%5f%n",zMD); 
				missing=true; // mark
			  }
			S[k] = sTrueInterp[0].interpolate(t,k!=0);
			varS += sqr(S[k]);
		  } 
		if (missing) {_WE.EH.printf("Can't do inversionPredictionWavelet(): log data incomplete over interval!%n"); return;}

		// variance estimate
		varS/=n;		
		// heron, scale S by variance
		Sscale = sqrt(varS);
		for (k=0; k<n ; k++) S[k]/=Sscale;
		S2 = n ; // == ||S||^2, used later on
		
		_WE.EH.printf("Seismic scaled by %4f%n",Sscale);
		if (false)
		  {
			// linear background
			// dopey linear model in sample-number for mu = <log(AI)>
			double AIbar=0;
			t = nm1/2.0;
			for (k=0; k<n ; k++) AIbar += AI[k];
			AIbar/=n;
			for (k=0,a=0,b=0; k<n ; k++) {a += (k-t)*(AI[k]-AIbar); b += sqr(k-t);}
			// slope
			b = a/b;
			// intercept
			a = AIbar-b*t;
			for (k=0; k<n ; k++) mu[k] = a + b*k;
		  }
		else
		  {
			int Lma=30;
			double[] wt = new double[Lma+1];
			for(i=0;i<=Lma;i++) wt[i] = Math.exp(-sqr(2.0*i/(double)Lma));
			// MA filter
			for (k=0; k<n ; k++)
			  {
				for(j=max(0,k-Lma),u=0,b=0; j<=min(n-1,k+Lma); j++) 
				  {					
					if ((i=abs(k-j)) <= Lma) 
					  {
						a = wt[i];
						u += AI[j]*a;
						b += a;
					  }
				  }
				mu[k] = u/b;
			  }
		  }
		for (k=0,varA=0; k<n ; k++) varA += sqr(AI[k] - mu[k]);
		varA/=n;

		Util.printfileZO("AI.txt",AI,n);
		Util.printfileZO("S.txt",S,n);
		Util.printfileZO("mu.txt",mu,n);

		if (true)
		  {
			// A "Basic" inversion
			_WE.EH.printf("Vanilla near-stack inversion test...%n");
			// wavelet half-length in samples
			Wavelet W = _WE.m.wavelet;
			nk = W.getLength();
			// inverse prediction problem...		
			DoubleMatrix2D mw,Sy,ym,X,XtS,H;
			CholeskyDecomposition cholDec;
			mw = new DenseDoubleMatrix2D(nk,1); // has starting vector 0

			double[] _values = W.values[0];
			// precursor: 
			for (j=W.finalFineSample+1,k=0; j<_values.length; j++,k++) mw.setQuick(k,0,_values[j]/Sscale);
			// coda
			for (j=0; j<=W.finalFineSample; j++,k++) mw.setQuick(k,0,_values[j]/Sscale);

			Util.debugPrintMatrixToFile("mw.txt",mw);
			H = new DenseDoubleMatrix2D(n,n);
			// X Jacobian
			X = new DenseDoubleMatrix2D(n,n);
			Sy = new DenseDoubleMatrix2D(n,1);
			for(i=0;i<n;i++) Sy.setQuick(i,0,S[i]);
			
			fillX_exptl(n,nk,W.finalFineSample,mw,X);
			
			// noise level a, rescaled to new seismic scale
			if (!_WE.fixNoise)
			  { // take as stack 0
				j = ((Integer)mapToIndex.get(_WE.m.seismic_noise[0])).intValue();
				a = Math.exp(mu_post.getQuick(j))/Sscale;
			  } 
			else 
			  {
				a = _WE.m.seismic_noise[0].value/Sscale;
			  }
			_WE.EH.printf("Noise scale in rescaled units = %6f%n",a);
			// a = sigma_noise. Cd = diag(a^2)
			// H = X'.C_d^{-1}.X + Cp^{-1}			
			lam = 1.0/(a*a);
			X.zMult(X,H,lam,0,true,false); 
			for(i=0;i<n;i++) H.setQuick(i,i,H.getQuick(i,i) + 1.0/varA); // += Cp^{-1}
			
			// rhs  = X'.Cd^{-1}.S + Cp^{-1}.mu				
			XtS = X.zMult(Sy,null,lam,0,true,false);
			for(i=0;i<n;i++) XtS.setQuick(i,0,XtS.getQuick(i,0) + mu[i]/varA);
			ym =(new CholeskyDecomposition(H)).solve(XtS);
			
			// ym = prediction for log(AI)
			_WE.EH.printf("Basic near-stack inversion in ym_basic.txt%n");
			Util.debugPrintMatrixToFile("ym_basic.txt",ym);
			X.zMult(ym,Sy,1.0,0,false,false); // Sy has ML synthetic
			// rms misprediction error, various PEP diagnostics...Xcc = corr. coefficient for AI prediction
			double sxy=0,sx=0,sy=0,sxx=0,syy=0;
			for(i=0,b=0,u=0,x=0;i<n;i++)
			  {
				v = Math.exp(AI[i]);
				w = Math.exp(ym.getQuick(i,0));
				x += sqr(w-v);				
				b += sqr(S[i]-Sy.getQuick(i,0));
				u += sqr(S[i]);
				sx+=v; sy+= w; sxy+=w*v; sxx+= v*v; syy+= w*w; 
			  }
			double PEP = 1-b/u, PEP_AI = 1-x/sxx, xcc = (n*sxy-sx*sy)/sqrt((n*sxx-sqr(sx))*(n*syy-sqr(sy))); 
			_WE.EH.printf("Basic near-stack inversion RMSE(AI) = %6g, RMS-seismic-misfit=%6g, PEP(seismic)=%6g, "
                          +"PEP(AI)=%6g XCC(AI)=%5g%n",
						  sqrt(x/n),Sscale*Math.sqrt(b/n),PEP,PEP_AI,xcc);
			if (!true)// FIXME 
			  {	
				_WE.EH.printf(".....Inversion prediction experiments%n%n");
				return;
			  }
		  }
	

		// std inference
		Wavelet W = _WE.m.wavelet;
		DoubleMatrix2D J,H;
		boolean uselaplacepenalty=true;
		boolean usez0 = true;
		boolean optimiseSigmaNoise = false; 
		boolean addFwdAmplitudeConstraints = false;

		// joint amplitude-fidelity weight (used if addFwdAmplitudeConstraints=true)
		double beta_f= 1; 
		
		// Cd^{-1} diag
		double lambda_start = 10.0;

		// smoothness constraints
		double lambdaE = 10;  // 100. ends, first moment
		double lambdaS = 2;    // 1.0 smoothness

		// wavelet half-length in samples
		L=13;
		//for(L=2;L<60;L+=2)
		//for(L=2;L<30;L+=2)
		  {
			_WE.EH.printf("Runnng L=%d%n",L);
			nk = 2*L+1;  // # wav samples
		
			nd = n + nExtraConstraints_exptl(nk,uselaplacepenalty,usez0);
			if (addFwdAmplitudeConstraints)
			  {
				nd3 = n;
				if (optimiseSigmaNoise) nd3++;
				nd += nd3;
			  }
			// inverse prediction problem...		
			DoubleMatrix2D mw,Sy,y_tmp =null;;
			CholeskyDecomposition cholDec;
			int nkpn = nk;  // wavelet coeffs

			if (optimiseSigmaNoise) nkpn++; // lambda = 1/sigSSqr

			mw = new DenseDoubleMatrix2D(nkpn,1); // has starting vector 0
				
			if (true)
			  {
				// Vanilla extraction as "starting guess"
				J = new DenseDoubleMatrix2D(nd,nk);					
				//true Jacobian
				a = beta_f;
				for(i=0,k=0;i<n;i++,k++)
				  {
					for(j=0;j<nk;j++)
					  {
						kk=i-j+L;
						if(kk>=0 && kk<(nm1)) J.setQuick(i,j,a * 0.5*(AI[kk+1]-AI[kk]));
					  }
				  }
				addConstraints_exptl(n,nk,nd,usez0,uselaplacepenalty,lambdaE,lambdaS,J);
				H = new DenseDoubleMatrix2D(nk,nk);				
				Sy = new DenseDoubleMatrix2D(nd,1);
				for(i=0;i<n;i++) Sy.setQuick(i,0,a*S[i]);
				// J'.J
				J.zMult(J,H,1.0,0,true,false); 
				// J'.S
				DoubleMatrix2D Jty = J.zMult(Sy,null,1,0,true,false); 
				// solve (J'.J). w = J'S.
				cholDec  = new CholeskyDecomposition(H);
				DoubleMatrix2D mwf = cholDec.solve(Jty); 
				// mw now has correct wavelet coeffs
				dumpSymmWaveletAscii_exptl("mw_basic.txt",mwf,dt);
				mw.viewPart(0,0,nk,1).assign(mwf);
			  }
			else
			  {
				mw = new DenseDoubleMatrix2D(nkpn,1); // has starting vector 0
			  }

			J = new DenseDoubleMatrix2D(nd,nkpn);
				
			MAXITS = 50;

			// noise level lam = 1/Cd in inversion ( S is scaled to rms=1)...
			lam = lambda_start; // = 1/sigS^2
			if (optimiseSigmaNoise) mw.setQuick(nk,0,lam);

			boolean conv=false;
			DoubleMatrix2D ym,X,Wa,D,XtS,dS;
			H = new DenseDoubleMatrix2D(n,n);
			// X Jacobian
			X = new DenseDoubleMatrix2D(n,n);
			Sy = new DenseDoubleMatrix2D(n,1);
			dS = new DenseDoubleMatrix2D(n,1);
			for(i=0;i<n;i++) Sy.setQuick(i,0,S[i]);

			for(its=0; its<MAXITS ;its++)
			  {
				// form X
				fillX_exptl(n,nk,L,mw,X);

				//Util.debugPrintMatrixToFile("X.txt",X);

				// H = X'.X + Cp^{-1}
				X.zMult(X,H,lam,0,true,false); 
				for(i=0;i<n;i++) H.setQuick(i,i,H.getQuick(i,i) + 1.0/varA);

				// rhs  = X'S.lam+Cp^{-1}mu				
				XtS = X.zMult(Sy,null,lam,0,true,false);
				for(i=0;i<n;i++) XtS.setQuick(i,0,XtS.getQuick(i,0) + mu[i]/varA);
				cholDec  = new CholeskyDecomposition(H);
				ym = cholDec.solve(XtS); 

				// Correct here.
				// ym = prediction for log(AI)
				Util.debugPrintMatrixToFile("ym.txt",ym);

				// dS = S-X.ym
				dS.assign(Sy);
				X.zMult(ym,dS,-1.0,1.0,false,false);
				//Util.debugPrintMatrixToFile("dS.txt",dS);			
				// A
				DoubleMatrix2D A = new DenseDoubleMatrix2D(n,nk);
				for(i=0;i<n;i++)
				  {
					for(j=0;j<nm1;j++)
					  {
						kk=i-j+L;
						if (kk<0 || kk>=nk) continue;
						a = ym.getQuick(j+1,0)-ym.getQuick(j,0);
						A.setQuick(i,kk,0.5*a);
					  }
				  }
				//Util.debugPrintMatrixToFile("A.txt",A);
				// R 
				DoubleMatrix2D R = new DenseDoubleMatrix2D(n,nk);
				for(i=0;i<n;i++)
				  {
					for(j=0;j<nk;j++)
					  {
						kk = i+j-L;
						if(kk<0)continue;
						a=0;
						if (i<nm1 && (kk<n)) a = -dS.getQuick(kk,0);
						if(i>0 && (kk>0) && (kk<=n)) a += dS.getQuick(kk-1,0);
						R.setQuick(i,j,0.5*a);
					  }
				  }
				//Util.debugPrintMatrixToFile("R.txt",R);

				// Total effective Jacobian by adjoint-state trick
				// R <- -X^T.A lam + R lam
				X.zMult(A,R,-lam, lam,true,false);

				//Util.debugPrintMatrixToFile("Rfinal.txt",R);

				J.assign(0);
				DoubleMatrix2D JJ = cholDec.solve(R);
				DoubleMatrix2D J1 = J.viewPart(0,0,n,nk);
				DoubleMatrix2D J1n = J.viewPart(0,0,n,nkpn);
				J1.assign(JJ); // copy JJ to top block of J

				//Util.debugPrintMatrixToFile("J.txt",J);
				if (optimiseSigmaNoise)
				  {
					DoubleMatrix2D XTDS = new DenseDoubleMatrix2D(n,1);
					// XTDS = X^T.(S-X.ym)
					X.zMult(dS,XTDS,1.0,0,true,false);
					DoubleMatrix2D Jn = cholDec.solve(XTDS);
					DoubleMatrix2D Jlc = J.viewPart(0,nk,n,1); // last column of J
					//Util.debugPrintMatrixToFile("Jlc.txt",Jlc);
					Jlc.assign(Jn);
					//Util.debugPrintMatrixToFile("J2.txt",J);
				  }
					
				//Util.debugPrintMatrixToFile("JJ.txt",JJ);
				// here JJ = n x nk
				
				if (false)
				  {
					// Compute JJ hard way using finite diffs & check....
					DoubleMatrix2D JJcheck=new DenseDoubleMatrix2D(n,nkpn);
					double eps = 1e-5;
					DoubleMatrix2D Ht = new DenseDoubleMatrix2D(n,n);
					for(kk=0,b=0;kk<nk;kk++) b = max(b,abs(mw.getQuick(kk,0)));
					for(kk=0;kk<nkpn;kk++)
					  {
						a = mw.getQuick(kk,0);
						double da = (a!=0) ? (((kk!=nk) ? b : abs(a)) * eps) : eps;
						mw.setQuick(kk,0,a+da);
						if (optimiseSigmaNoise) lam = mw.getQuick(nk,0); 
						// form X
						fillX_exptl(n,nk,L,mw,X);
						// X'.X lam + invCp
						X.zMult(X,Ht,lam,0,true,false); 
						for(i=0;i<n;i++) Ht.setQuick(i,i,Ht.getQuick(i,i) + 1.0/varA);
						// rhs X'S lam+Cp^{-1}mu				
						XtS = X.zMult(Sy,null,lam,0,true,false);
						for(i=0;i<n;i++) XtS.setQuick(i,0,XtS.getQuick(i,0) + mu[i]/varA);
						cholDec  = new CholeskyDecomposition(Ht);
						DoubleMatrix2D ymn = cholDec.solve(XtS);
						// FD Jacobian
						for(i=0;i<n;i++) JJcheck.setQuick(i,kk,(ymn.getQuick(i,0) - ym.getQuick(i,0))/da); 
						mw.setQuick(kk,0,a);
						if (optimiseSigmaNoise) lam = mw.getQuick(nk,0); 
					  }
					Util.debugPrintMatrixToFile("JJcheck.txt",JJcheck);
					if (true)
					  {
						// check J diffs
						for(i=0,a=0;i<n;i++)
						  {
							for(j=0;j<nkpn;j++) a= max(a,J.getQuick(i,j));
						  }
						DoubleMatrix2D ERR = J.like();
						for(i=0,u=0;i<n;i++)
						  {
							for(j=0;j<nkpn;j++)
							  {
								ERR.setQuick(i,j,b = abs(J.getQuick(i,j)
														 -JJcheck.getQuick(i,j))/a);
								u=max(u,b);
							  }
						  }
						_WE.EH.printf("Errmax %f%n",u);
						Util.debugPrintMatrixToFile("Jerr.txt",ERR);	
						//System.exit(0);
					  }
				  }

				// now Gauss-Newton				
				//  typical wav coeff ~ 10, with S scaled @ S_rms=1;

				// Rescale J1 for chisq measured in terms of log(AI) variability scale: ||log(AI) - f(m)||_sigAI
				// J1 -> J1/sigAI
				double sigAI = sqrt(varA); // target precision for AI
				J1n.assign(cern.jet.math.Functions.functions.div(sigAI));

				// add lower block J2 constraints
				addConstraints_exptl(n,nk,nd,usez0,uselaplacepenalty,lambdaE,lambdaS,J);
				//Util.debugPrintMatrixToFile("Jfull.txt",J);	
					
				// effective data vector
				// yA - ym(mw_0) + JJ.wm_0
				DoubleMatrix2D y_eff = new DenseDoubleMatrix2D(nd,1);
				DoubleMatrix2D y_eff1 = y_eff.viewPart(0,0,n,1);
				for(i=0,E0=0;i<n;i++)
				  {
					a = (AI[i] - ym.getQuick(i,0))/sigAI;
					y_eff.setQuick(i,0,a);                // yeff = yA - ym(m_0)
					E0 += sqr(a);
				  }
				// last blocks for constraints have effective-data yA=0 (constraints are perfectly linear)
				E0 += getEconstraints_exptl(n,nk,nd,usez0,uselaplacepenalty,lambdaE, lambdaS, mw);
				if (addFwdAmplitudeConstraints)
				  {
					DoubleMatrix2D S3 = X.zMult(ym,null,1.0,0,false,false);
					for(i=0,SS=0;i<n;i++) SS += sqr(S[i]-S3.getQuick(i,0));					
					E0 += sqr(beta_f) * lam * SS;
					if (optimiseSigmaNoise) E0 -= sqr(beta_f)* n* Math.log(lam);
				  }
				// top-block adjustment for += J.m0
				J1n.zMult(mw,y_eff1,1.0,1.0,false,false); 

				if (addFwdAmplitudeConstraints)
				  {
					// Block-3 adjustment for += J.m0
					DoubleMatrix2D J3 = J.viewPart(nd-nd3,0,nd3,nkpn);
					//true Jacobian
					u = beta_f * sqrt(lam);
					for(i=0,k=0;i<n;i++,k++)
					  {
						for(j=0;j<nk;j++)
						  {
							kk=i-j+L;
							if(kk>=0 && kk<(nm1)) J3.setQuick(i,j,u*(AI[kk+1]-AI[kk])/2); // std X = 0.5*(AI[kk+1]-AI[kk])
						  }						
					  }
					if (optimiseSigmaNoise)
					  {
						// rightmost column
						for(i=0;i<n;i++) J3.setQuick(i,n,u*S[i]);
						// nn element
						J3.setQuick(n,n,beta_f*sqrt(n/2.0)/lam);
					  }
					// y_eff elements. dS = S-Xe.m_w already
					DoubleMatrix2D y_eff3 = y_eff.viewPart(nd-nd3,0,nd3,1);
					b = beta_f * sqrt(lam);
					for(i=0;i<n;i++) y_eff3.setQuick(i,0,b*S[i]);					
					if (optimiseSigmaNoise)
					  {
						a = beta_f * ( ((n-lam*SS)/2-lam*S2)/sqrt(n/2.0) + sqrt(n/2.0));
						y_eff3.setQuick(n,0,a);
					  }
				  }

				//Util.debugPrintMatrixToFile("y_eff.txt",y_eff);
				//Util.debugPrintMatrixToFile("Jfull.txt",J);
				//System.exit(0);

				H = new DenseDoubleMatrix2D(nkpn,nkpn);
				// J'.J
				J.zMult(J,H,1.0,0,true,false);

				// J'.y_eff
				XtS = J.zMult(y_eff,null,1.0,0,true,false);
				// OLS solve
				cholDec  = new CholeskyDecomposition(H);
			
				DoubleMatrix2D mw_new = cholDec.solve(XtS);
				//mw_new = (new QRDecomposition(H)).solve(XtS);
				Util.debugPrintMatrixToFile("mw_new."+its+".txt",mw_new.viewPart(0,0,nk,1));
				dumpSymmWaveletAscii_exptl("mw_final.txt",mw_new.viewPart(0,0,nk,1),dt);
					
				// Linesearch backtrack
				H = new DenseDoubleMatrix2D(n,n);
				DoubleMatrix2D mw0 = new DenseDoubleMatrix2D(nkpn,1);
				mw0.assign(mw);
				if (true)
				  {
					// gradient check. (yA-ym).J.(m'-m0).
					DoubleMatrix2D delta_m = mw.like();
					for(i=0;i<nkpn;i++) delta_m.setQuick(i,0,mw_new.getQuick(i,0)-mw0.getQuick(i,0));
					DoubleMatrix2D dy = ym.like();
					for(i=0;i<n;i++) dy.setQuick(i,0,AI[i] - ym.getQuick(i,0));
					DoubleMatrix2D gradm = J1n.zMult(dy,null,-1.0/sigAI,0,true,false);

					DoubleMatrix2D Jc = J.viewPart(n,0,nd-n,nkpn);
					DoubleMatrix2D y2 = Jc.zMult(mw0,null,1.0,0,false,false);
					// grad += Jc'.(Jc.mw0-yeff)
					if (addFwdAmplitudeConstraints)
					  {							
						DoubleMatrix2D y_eff23 = y_eff.viewPart(n,0,nd-n,1);
						for(i=n,j=0;i<nd;i++,j++) y2.setQuick(j,0,y2.getQuick(j,0)-y_eff23.getQuick(j,0));
					  }
					Jc.zMult(y2,gradm,1.0,1.0,true,false);
						
					Util.debugPrintMatrixToFile("grad."+its+".txt",gradm);
					for(i=0,u=0;i<nkpn;i++) u += delta_m.getQuick(i,0) * gradm.getQuick(i,0);

					_WE.EH.printf("|grad|= %15.8f%n",u);						
					grad=u;
				  }
				lam_NP = lam;
				if(false)
				  {
					// test profile of E vs lambda
					for(a=.01;a<20;a*=1.25)
					  {
						// backtrack
						lam = a;
						// form X
						fillX_exptl(n,nk,L,mw,X);
						// X'.X*lam + invCp
						X.zMult(X,H,lam,0,true,false); 
						for(i=0;i<n;i++) H.setQuick(i,i,H.getQuick(i,i) + 1.0/varA);
						// rhs X'S lam+Cp^{-1}mu				
						XtS = X.zMult(Sy,null,lam,0,true,false);
						for(i=0;i<n;i++) XtS.setQuick(i,0,XtS.getQuick(i,0) + mu[i]/varA);
						CholeskyDecomposition cholDecLS  = new CholeskyDecomposition(H);
						DoubleMatrix2D ymn = cholDecLS.solve(XtS);
						for(i=0,E=0;i<n;i++) E += sqr(AI[i] - ymn.getQuick(i,0)); // RSS
						E /= sqr(sigAI);
						Ec = getEconstraints_exptl(n,nk,nd,usez0,uselaplacepenalty,lambdaE, lambdaS, mw);
						E += Ec;
						//E += Math.log(lam);
						_WE.EH.printf("LAM %15.8f %15.8f%n",a,E); 
					  }
				  }
				{
				  // step limit.
				  if (optimiseSigmaNoise)
					{
					  u = abs(mw_new.getQuick(nk,0) - mw.getQuick(nk,0))/mw.getQuick(nk,0);
					  if(u>1.1) u = (1.1/L);
					}
				  else u=1;
				}
				for(a=1.0*u; abs(a)>1e-5; a*=0.5)
				  {
					// backtrack
					for(i=0,b=1.0-a; i<nkpn; i++) mw.setQuick(i,0,b*mw0.getQuick(i,0) + a*mw_new.getQuick(i,0));
					if (optimiseSigmaNoise)
					  {
							
						lam_NP = mw.getQuick(nk,0);
						if (lam_NP<=0) continue; //contract newt pt,, noise must be +ve
					  }
					// form X
					fillX_exptl(n,nk,L,mw,X);
					// X'.X*lam_NP + invCp
					X.zMult(X,H,lam_NP,0,true,false); 
					for(i=0;i<n;i++) H.setQuick(i,i,H.getQuick(i,i) + 1.0/varA);
					// rhs X'S lam_NP + Cp^{-1}mu				
					XtS = X.zMult(Sy,null,lam_NP,0,true,false);
					for(i=0;i<n;i++) XtS.setQuick(i,0,XtS.getQuick(i,0) + mu[i]/varA);
					CholeskyDecomposition cholDecLS  = new CholeskyDecomposition(H);
					DoubleMatrix2D ymn = cholDecLS.solve(XtS);
					for(i=0,E=0;i<n;i++) E += sqr(AI[i] - ymn.getQuick(i,0)); // RSS
					E /= sqr(sigAI);
					double Ed = E;
					Ec = getEconstraints_exptl(n,nk,nd,usez0,uselaplacepenalty,lambdaE, lambdaS, mw);

					if (addFwdAmplitudeConstraints)
					  {
						DoubleMatrix2D S3 = X.zMult(ymn,null,1.0,0,false,false);
						for(i=0,b=0;i<n;i++) b += sqr(S[i]-S3.getQuick(i,0));
						E += sqr(beta_f)* lam_NP * b;
						if (optimiseSigmaNoise) E -= sqr(beta_f) * n * Math.log(lam_NP);
						Util.debugPrintMatrixToFile("Sy.txt",S3);
					  }
						
					E += Ec;
					//_WE.EH.printf("IT %d %10f %15.8f %15.8f%n",its,a,E0,E); 
					if (E < E0)
					  {
						if ((abs(E-E0) < 0.0001*E) && (abs(grad)<0.1)) conv=true;

						E0=E;
						if (optimiseSigmaNoise) lam = lam_NP;
						_WE.EH.printf("(S) alpha=%5f Ed=%15.8f En=%15.8f lambdaN= %f%n",a,Ed,E,lam_NP);
						break; // a--  backtrack
					  }
					else
					  {
						_WE.EH.printf("....Backtrack alpha=%5f, E0=%10f  E=%10f lambdaN= %f%n",a,E0,E,lam_NP);
					  }
				  } // backtrack
				if(abs(a)<=1e-5) conv=true;
				if (conv)
				  {
					// discount factor for coarse spacing
					double Npdiscount = 1.0;//dtSub/_WE.multiModelSetUp.dt_sparse;

					// BIC =(n/2)log(RSS/n)+(1/2)(logn)(p+1)
					double RSS = (E-Ec)*sqr(sigAI); // data
					RSS += Ec; 
					_WE.EH.printf("BIC DIAGNOSTIC L= %d %15f %15f %15f%n",L,E/2,Ec/2, 0.5*nd*Math.log(RSS/nd)+(nk*Npdiscount/2)*Math.log(nd));
					if(false)
					  {
						// tr{Hatmatrix} estimate of p - ends up yielding p=nk
						// H=J'.J
						DoubleMatrix2D Jt= J.viewDice();
						DoubleMatrix2D invHJt = cholDec.solve(Jt);
						for(i=0,a=0;i<nd;i++)
						  {
							for(j=0,u=0;j<nk;j++) u += J.getQuick(i,j)*invHJt.getQuick(j,i);
							a+=u;
						  }
						_WE.EH.printf("P-EFF %d tr{HAT}= %f%n",nk,a);
					  }
					break; // Gauss Newt
				  }
			  }
		  }
	  }// wavelet span
	//_WE.EH.printf(".....Inversion prediction experiments%n");
	System.exit(0);  
  }

  private void dumpSymmWaveletAscii_exptl(String filename,DoubleMatrix2D m,double dt)
  {
	PrintWriter out;
	// open ascii symm wavelet file
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }
	int i,n;
	n = m.rows();
	double t = -dt*(n-1)/2.0;
	for(i=0;i<n;i++,t+=dt) out.printf("%f %f%n",t,m.getQuick(i,0));
	out.close();
	_WE.EH.printf("Wrote: %s%n",filename);
  }

  private void fillX_exptl(int n,int nk,int L,DoubleMatrix2D mw,DoubleMatrix2D X)
  {
	int i,j,k,nm1=n-1;
	double b;
	for(i=0;i<n;i++)
	  {
		for(k=-1; k<nk; k++)
		  {
			j = i-k+L;
			if(j<0 || j>=n) continue;
			b = (j<nm1 && k>=0) ? -mw.getQuick(k,0) : 0;
			if (j!=0 && j<n && (k+1)<nk) b += mw.getQuick(k+1,0);
			X.setQuick(i,j,0.5*b);
		  }
	  }
  }
  
  private int nExtraConstraints_exptl(int nk,boolean uselaplacepenalty,boolean usez0)
  {
	int n;
	// (1) Fix ends X.w=0
	n=1;
	if (nk>2) n++;

	// (2) fix z0 constraint \sum w_j = 0
	if(usez0) n++;
	// (3) smoothness bit, roughness operator
	if(!uselaplacepenalty)
	  { 
		if(nk>1) n += (nk-1);
	  }
	else
	  { // laplace penalty 
		if(nk>0) n += nk;
	  }
	return n;
  }
  
  // To experimental sample-based wavelet, add end-fixing , smoothness, and z0 constraints to Jacobian
  private void addConstraints_exptl(int n,int nk,int nd,
									boolean usez0,  // set first moment=0
									boolean uselaplacepenalty, // smoothness via ||Dw|| or ||D^2 w||
									double lambdaE, //end-weight
									double lambdaS,  // smooth weight
									DoubleMatrix2D J)
  {
	int i,j,k,nkm1=nk-1;
	double t;

	// at this point, natural" J elements insitu are rows 0...n-1, cols 0...nk-1 (wrt samples)
	k=n;

	// (1) Fix ends X.w=0
	J.setQuick(k++,0,lambdaE);
	if(nk>2) J.setQuick(k++,nk-1,lambdaE);
	//J.setQuick(k++,1,lambdaE);
	//J.setQuick(k++,nk-2,lambdaE);
	
	// (2) fix z0 constraint \sum w_j = 0
	if(usez0) 
	  {
		for(j=0;j<nk;j++) J.setQuick(k,j,lambdaE);
		k++;
	  }
	
	// (3) smoothness bit E += lambda^2||R.w||^2, some roughness operator
	// =3.0 mild
	if(!uselaplacepenalty)
	  { // penalty ||D.m||^2
		for(i=0;i<nkm1;i++,k++)
		  {
			J.setQuick(k,i,-lambdaS);
			J.setQuick(k,i+1,lambdaS);
		  }
	  }
	else
	  { // laplace penalty ||D^2.m||^2
		t=2*lambdaS;
		for(i=0; i<nk;i++,k++)
		  {
			if(i>0) J.setQuick(k,i-1,-lambdaS);
			J.setQuick(k,i,t);
			if (i < nkm1)J.setQuick(k,i+1,-lambdaS);
		  }	
	  }
  }

  private double getEconstraints_exptl(int n,int nk,int nd,
									   boolean usez0,  // set first moment=0
									   boolean uselaplacepenalty, // smoothness via ||Dw|| or ||D^2 w||
									   double lambdaE, //end-weight
									   double lambdaS,  // smooth weight
									   DoubleMatrix2D ws)
  {
	int i;
	double E,u,ul;
	// at this point, natural" J elements insitu are rows 0...n-1, cols 0...nk-1 (wrt samples)

	// (1) Fix ends X.w=0
	u = sqr(ws.getQuick(0,0));
	if (nk>2) u += sqr(ws.getQuick(nk-1,0));
	E = sqr(lambdaE)*u;

	// (2) fix z0 constraint \sum w_j = 0
	if (usez0) 
	  {
		for(i=0,u=0;i<nk;i++) u += ws.getQuick(i,0);
		E += sqr(lambdaE *u);
	  }
	
	// (3) smoothness bit E += lambda^2||R.w||^2, some roughness operator
	int nkm1=nk-1;
	if (!uselaplacepenalty)
	  { // penalty ||D.m||^2
		for(i=0,u=0;i<nkm1;i++)
		  {
			u += sqr(ws.getQuick(i+1,0) - ws.getQuick(i,0));
		  }
		E += sqr(lambdaS)*u;
	  }
	else
	  {
		// laplace penalty ||D^2.m||^2 	   
		ul=0;
		for(i=0; i<nk;i++)
		  {			
			u = (i>0) ? -ws.getQuick(i-1,0) : 0;
			u += 2* ws.getQuick(i,0);
			if (i < nkm1) u -= ws.getQuick(i+1,0);
			ul += sqr(u);
		  }	
		 E += sqr(lambdaS)*ul;
	  }
	return E;
  }

  
  /**
	 ASCII dump of useful file for multiple graphing ala commercial packages.
	 Plots, on common time scale, salient logs, TVT, MD, seismics etc.
  */
  void dumpMLgraphingfile(String filenamePrefix)
  {
 	int i,j,k,kk,kTD,n,index;
	byte facies_indicator;
	double t,tOffset=0, dt,dtSub,zTVD,zMD,sTrue,sSynth,tFinish,rho,DT_scaling;
	double DT_Vp_sonic,DT_Vs_sonic,vApprox,DT_blocked=0,gamma=0;
	boolean haveLogData, haveGamma, haveFacies;
	double logDataZmin, logDataZmax;
	LinearInterpolation[] sTrueInterp,sSynthInterp;
	LinearInterpolation gammaInterp=null;
	String filename=null;
	Property prop;
	double[] _Zm;
	double T_EPS = FundamentalConstants.BASE_SCALING_TIME * 1.0e-3;
	double TIME_SUBSAMPLING_FACTOR=16.0;
	double[] segTimes=null;     // mapping of layer tops to times, using current realisation (from segmented data)
	double[] segDepths=null;    // layer tops in true depth (from segmented data)
	double[] segMDepths=null;    // layer tops in measured depth (from segmented data)
	double[] segVelocity=null;  // layer velocities (from segmented data)
	double[] segDensity=null;  // layer density (from segmented data)


	PrintWriter out;
	


	// we need to compute the fwd model before dumping the graphing file
	computeMLForwardModel();

	try
	  {	
		kTD = firstTDparamIndex;
		for (i=0; i<_WE.m.nWells; i++)
		  {
			LogData LDI = _WE.logData[i];
			filename = filenamePrefix + "."+Util.fileBasename(_WE.wellInfo[i].logFilename);
			filename = filename.replaceAll("(?i)\\.las|\\.txt","") + ".txt";
			filename = filename.replaceAll("%","%%");
			
			haveGamma = (_WE.logData[i].gamma != null);
			haveFacies = (_WE.doAnisotropicModel || _WE.useSaturationModel);

			DT_scaling = _WE.logData[i].getSonicLogScaling();

			// open ascii summary file
			_WE.EH.println("Writing summary of ML generic graphs to "+filename);
			// out = new PrintWriter(new FileWriter(filename));
			out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));

			out.printf("Generic ML summary file for plotting purposes, well %s",_WE.wellInfo[i].logFilename);

			if (_WE.m.petro.petrophysicsProvided)
			  {
				out.printf(_WE.m.petro.encodingString());
			  } 
			// units
			out.printf(" Units:["+ LDI.getDTandRhoUnits()+"]");

			// add wavelet filenames
			out.printf(_WE.waveletFilenamesEncodingString());

			out.println("");
			out.printf("%-15s ","t");
			if (_WE.m.registrationErrorsAreCommon) out.printf("%-15s ","t_seismic");
			out.printf("%-15s %-15s %-15s %-15s %-15s %-15s","TVD","MD",
						   "rho","DT","V-interval","DT_blocked");
			out.printf("%-15s ","DTS");  // Shear log too
			if (haveGamma) out.printf("%-15s ","GR");
			if (haveFacies) out.printf("%-17s ","facies_indicator");
			// seismic
			out.printf("%-15s %-15s ","S_synth","S_true");
			for (j=2; j<=_WE.m.nStacks; j++) 
			  {
				out.printf("%-15s %-15s ",("S_synth_"+j),("S_true_"+j));
			  }
			out.println("");


			// arrays for seismic interpolation
			n = _WE.convModel.nZ;
			double[] times = new double[n];
			dt = _WE.m.wavelet.delta_t;
			for (t=_WE.convModel.tStartReff,k=0; k<n; k++,t+=dt) times[k] = t;
			tFinish = times[n-1];
			// array of interpolations for stacks
			sTrueInterp = new LinearInterpolation[_WE.m.nStacks];
			sSynthInterp = new LinearInterpolation[_WE.m.nStacks];


			// time sampling for graphs..
			dtSub = _WE.m.wavelet.delta_t/TIME_SUBSAMPLING_FACTOR;
			
			TimeToDepthMap t2dMap = _WE.m.timeToDepthMaps[i];
			// time to TVD mapping 
			LinearInterpolation zInterp = new LinearInterpolation(t2dMap.t, t2dMap.Z, t2dMap.t.length,_WE.EH);

			// TVD to MD mapping from survey. MUST BE MONOTONIC
			LinearInterpolation zToMDInterp = 
			  new LinearInterpolation(_WE.logData[i].survey.Z, _WE.logData[i].survey.Zm, _WE.logData[i].survey.Z.length, _WE.EH);
			// but don't ask for log data outside this range....
			_Zm = _WE.logData[i].Zm; // shortcut

			logDataZmin = _Zm[0];
			logDataZmax = _Zm[_Zm.length - 1];

			// density
			// raw..from measured depths			
			LinearInterpolation rhoInterp = new LinearInterpolation(_Zm,
																	_WE.logData[i].rho,_Zm.length,_WE.EH);
			// DT log
			LinearInterpolation sonicInterp = new LinearInterpolation(_Zm,
																	_WE.logData[i].DT,_Zm.length,_WE.EH);
			// DTS log
			LinearInterpolation shearInterp = new LinearInterpolation(_Zm,_WE.logData[i].DTshear,_Zm.length,_WE.EH);

			// possible gamma
			if (haveGamma)
			  {
				gammaInterp = new LinearInterpolation(_Zm,
													  _WE.logData[i].gamma,_Zm.length,_WE.EH);
			  }
			for (j=0; j<_WE.m.nStacks; j++) 
			  {
				// seismic interp for jth stack
				sTrueInterp[j] = new LinearInterpolation(times, _WE.convModel.S_true[i][j], n ,_WE.EH);
				sSynthInterp[j] = new LinearInterpolation(times, _WE.convModel.S_synth[i][j], n ,_WE.EH);
			  }

			// Blocked data
			SegmentedLogData sLogData = _WE.sLogData[i];

			// run over fine-scale times

			//for (t=_WE.convModel.tStartReff,k=0; t<tEnd; k++,t+=dtSub)
			kk=0;
			for (t=_WE.convModel.tStartReff,k=0; t<times[n-1]; k++,t+=dtSub)
			  {
				// no-reflection-outside logs option can extend the range of t  a little  so the
				// depth interpolations don't work. Skip these times
				tOffset = t;
				if (_WE.m.registrationErrorsAreCommon)
				  {
					tOffset += t2dMap.systematicTProp.value;
				  }

				if (t<_WE.tStart || t>_WE.tEnd) continue;
				// get TVD
				zTVD = zInterp.interpolate(t,k!=0);
				// get MD
				zMD = zToMDInterp.interpolate(zTVD,k!=0);
				// see if log data is there
				haveLogData = (zMD >=logDataZmin && zMD <= logDataZmax);

				// get log data

				rho = haveLogData ? rhoInterp.interpolate(zMD,k!=0) : 0;
				// sonic
				DT_Vp_sonic = haveLogData ? sonicInterp.interpolate(zMD,k!=0) : 0;
				//DT_Vp_sonic *= DT_scaling;

				// shear sonic
				DT_Vs_sonic = haveLogData ? shearInterp.interpolate(zMD,k!=0) : 0;
				//DT_Vs_sonic *= DT_scaling;

				// approximate velocity from finite differences: ensure two-way time is taken into account... 
				vApprox = 2.0*(zInterp.interpolate(t + T_EPS, k!=0) - 
							   zInterp.interpolate(t,k!=0))/T_EPS;

				// Blocked DT
				if (haveLogData)
				  {
					// 			 Get corresponding index in segmented suite
					// The depth lData.Zm[0,1...] has segments for the range of indices
					// [segIndexValues.getQuick(k),segIndexValues.getQuick(k+1)-1]
					// last entry in sLogData.segIndexValues is index or end of data (end segment)
					// search only over range of segments-start-indices
					int maxSegIndexValueToCheck = sLogData.segIndexValues.size()-2;
					for (;kk<=maxSegIndexValueToCheck; kk++)
					  {
						index = sLogData.segIndexValues.getQuick(kk);
						if (_Zm[index] > zMD) 
						  {
							if(kk>0) kk--;
							// expect zMD to be in segment 
							// sLogData.segIndexValues.getQuick(kk) <= index <= sLogData.segIndexValues.getQuick(k+1)-1
							break;
						  }
						// the blocked data for index 
					  }
					if (kk==maxSegIndexValueToCheck+1) kk--; // overran; step back 1

					//DT_blocked = (sLogData.Vp[kk]!=0) ? FundamentalConstants.SONIC_LOG_TIME_SCALING/sLogData.Vp[kk] : 0;					
					DT_blocked = (sLogData.Vp[kk]!=0) ? DT_scaling/sLogData.Vp[kk] : 0;					
					facies_indicator =  (haveFacies) ? sLogData.facies_indicator[kk] : 0;
				  }
				else 
				  {
					facies_indicator = 1; // default "shale"
				  }

				// gamma
				if (haveGamma)
				  {
					gamma = haveLogData ? gammaInterp.interpolate(zMD,k!=0) : 0;
				  }

				// print out log data & times & depths
				out.printf("%15.7g ",t);
				if (_WE.m.registrationErrorsAreCommon) out.printf("%15.7g ",tOffset);
				out.printf("%15.7g %15.7g %15.7g %15.7g %15.7g %15.7g ",
							   zTVD, zMD, rho, DT_Vp_sonic, vApprox, DT_blocked);
				out.printf("%15.7g ",DT_Vs_sonic);  // shear sonic
				if (haveGamma) out.printf("%15.7g ",gamma);
				if (haveFacies) out.printf("%17d ",facies_indicator);

				// append seismic stacks: synth, then true, in pairs, for each stack
				for (j=0; j<_WE.m.nStacks; j++) 
				  {
					// get amplitudes
					sSynth = sSynthInterp[j].interpolate(t,k!=0);
					sTrue = sTrueInterp[j].interpolate(t,k!=0);
					out.printf("%15.7g %15.7g",sSynth,sTrue);
				  }
				out.println("");
			  }
			out.close();
			_WE.EH.println("Wrote data to: " + filename);
		  }
	  }
	catch (java.io.IOException e) 
	  {
		_WE.EH.handleIOException(e);
	  }
  }


  /**
	 Compute the reflectivities, synthetic etc, as posted by values in
	 mu_post, assuming this field now contains the maximum -likelihood estimation
  */
  boolean computeMLForwardModel()
  {
	int i,n;
	double stretch,val;
	n = mapFromIndex.size();
	// stuff values from mu_post into fields
	for (i=0;i<n;i++)
	  {
		Property prop = (Property)mapFromIndex.get(i);
		// restore log(noise) to actual noise
		val = mu_post.getQuick(i);
		if (i>= nGaussian)
		  {
			// a log-transformed noise parameter
			val = Math.exp(val);
		  }
		prop.setValue(val);
	  }
	_WE.m.updateTimeToDepthMaps();
	if (!_WE.m.timesAreDecentlyOrdered)
	  {
		// time to depth map was illegal
		// return false;
	  }
	// update Wavelet
	if (_WE.useStretchedModel)
	  {
		// remap wav. stretch parameters
		for (i=0; i< (_WE.m.nStacks-1); i++)
		  {
			if (_WE.m.farStackWaveletTimeStretchMapper!=null) 
			  {
				_WE.m.farWavStretch[i] = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(_WE.m.farStackWaveletTimeStretch[i].value);
			  }
			else
			  {
				_WE.m.farWavStretch[i] = _WE.m.farStackWaveletTimeStretch[i].value;
			  }
			_WE.m.farWavAmp[i] = _WE.m.farStackWaveletAmplitude[i].value;
		  }
	  }
	// update Wavelet
	_WE.m.wavelet.updateWavelet();

	// compute new reflectivity arrays. Do not use cache, for safety
	if (_WE.useSaturationModel) _WE.m.recomputeLogProperties(false);
	_WE.m.computeReflectivities(false,false);

	if (_WE.useSaturationModel) 
	  {
		// copy reflectivities
		_WE.m.makeReflectivityCopy();
		// compute effect of sat change,
		_WE.m.recomputeLogProperties(true);
		// new reflectivities
		_WE.m.computeReflectivities(false,false);
		// then difference
		_WE.m.computeReflectivityDifference();
	  }

	// compute synth seismic 
	_WE.convModel.computeSyntheticSeismic();
	// and punch in the true amplitude via this call. 
	// we don't need the chisq, but we do need the "true seis" to have the right values
	// Avoid cache and logNoise parameterisation
	_WE.convModel.computeSyntheticNegLogLikelihood(false,false);

	return true;
  }

  /**
	 Return snapshot class of useful posterior entities in this class, at current state
  */
  ModelVariablesPosteriorSnapshot takePosteriorSnapshot()
  {
	return new ModelVariablesPosteriorSnapshot(this);
  }

  /**
	 Snapshot class to keep track of posterior quantities for intermediate solutions
  */
  public class ModelVariablesPosteriorSnapshot
  {
	/** 
	  Backup elements for the parts of ModelVariables that this corrupts, and then
	  must restore on completion
	*/
	private DoubleMatrix2D L_post;
	private int[] L_Permutation;
	private double logModeWeight,halfLogDetCeff,CDF,weight;
	private DoubleMatrix1D mu_post;
	private double[] CtildeInvDiag;

	/**
	   Instantiator makes deep copies of these things
	*/
	public ModelVariablesPosteriorSnapshot(ModelVariables modelVariables)
	{
	  L_post = modelVariables.L_post.like().assign(modelVariables.L_post);
	  L_Permutation = modelVariables.L_Permutation.clone();
	  logModeWeight = modelVariables.logModeWeight;
	  halfLogDetCeff = modelVariables.halfLogDetCeff;
	  CDF = modelVariables.CDF;
	  weight = modelVariables.weight;
	  mu_post = modelVariables.mu_post.like().assign(modelVariables.mu_post);
	  CtildeInvDiag = modelVariables.CtildeInvDiag.clone();
	}

	/** 
		Copy snapshot fields back into supplied modelVariables
	*/
	void restoreFromSnapshot(ModelVariables modelVariables)
	{
	  modelVariables.C_prior = C_prior;
	  modelVariables.L_post =	L_post;
	  modelVariables.L_Permutation = L_Permutation;
	  modelVariables.logModeWeight = logModeWeight;
	  modelVariables.halfLogDetCeff =	halfLogDetCeff;
	  modelVariables.CDF = CDF;
	  modelVariables.weight = weight;
	  modelVariables.mu_post = mu_post;
	  modelVariables.CtildeInvDiag = 	CtildeInvDiag;
	}
  }
}
