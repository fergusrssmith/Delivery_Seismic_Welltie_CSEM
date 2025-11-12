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

import java.util.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;

import cern.jet.random.*;
import cern.colt.matrix.linalg.*;
import cern.jet.random.engine.MersenneTwister;

import au.csiro.JamesGunning.lib.util.*;

/**
   Data class for computing basic setup of multi model required for sampling etc
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/


public class MultiModelSetUp
{
  // Link to calling extractor
  WaveletExtractor _WE;

  ModelVariables[] modelVariables;
  int mostLikelyWaveletIndex;
  
  int waveletCoarseIndexMin;
  int waveletCoarseIndexMax;

  /**
	 least safety margin required over set of all wavelets (+ve)
  */
  double maxPrecursorMargin;
  /**
	 max safety margin required over set of all wavelets (+ve)
  */
  double maxCodaMargin;

  /** 
	  Useful quantity for scaling, computed in setUpPriorModel()
  */
  double typicalRMSSeismicAmplitude;
  /** 
	  Useful quantity for scaling, computed in setUpPriorModel()
  */
  double typicalPeakSeismicAmplitude;

  /**
	 Useful estimate of data's peak frequency
  */
  double f_peak;

  /**
	 shared dt_sparse coarse spacing in time
  */
  double dt_sparse;
  /**
	 mode number indexed as we loop through all the various wavelet span possibilities
  */
  int modeNumber;

  boolean checkAnisModel=!true; // do both aniso and isotropic models

  boolean checkshotModelSelect = false;
  boolean[] checkshotModelSelectWells = null;
  ArrayList<int[]> checkshot_model_select_index_array=null;

  // times for longest wavelet (used in multi-Wavelet i/o routines)
  double tMinLongestWavelet=0,tMaxLongestWavelet=0;

  // standard Gaussian and Uniform generators
  MersenneTwister ranEngine;
  AbstractDistribution Ngen;
  AbstractDistribution Ugen;

  /**
	 The sampler can sometimes generate illegal states. 
	 Keep trying to sample this many times before 
	 bombing out
  */
  public static final int MAX_NO_SAMPLE_ATTEMPTS = 10000;

  /**
  Constructor.
  */
  public MultiModelSetUp( WaveletExtractor WE)
  {
	this._WE=WE;
	// deviate generators
	ranEngine = new MersenneTwister(_WE.useSeed ? _WE.seed : MersenneTwister.DEFAULT_SEED);
	Ngen = new Normal(0,1,ranEngine);
	Ugen = new Uniform(0,1,ranEngine);
  }

  /**
	 Here we create arrays of basic models, one for each different 
	 wavelet span. Sets of useful proeprties are also instantiated
  */

  void setUpPriorModel()
  {
	int i,j,k,i_c,Imax=0,Nd=0,iS, Nwa;
	int iStart, iEnd,winSize,halfWinPad;
	int size,halfLoopCoarseIndex=0;
	double Arms=0,u;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D threeDdata = null;

	// Average bandwidth estimates etc over all wells, stack 0.
	// else estimate this from 1st well, 1st stack, 1st trace.
	boolean useWell1Only=false;

	double data[];
	

	// get basic wavelet time discretisation from pwr spectrum
	typicalPeakSeismicAmplitude = 0;
	dt_sparse=0;
	Nwa = useWell1Only ? 1 : _WE.m.nWells;
	for (iS=0; iS < Nwa; iS++) 
	  {
		SeismicStackInfo ssInfo = _WE.seismicInfo[iS][0];
		threeDdata = ssInfo.data;
		// use mid-cube trace as example
		DoubleMatrix1D longdata = threeDdata.viewSlice(threeDdata.slices()/2).viewRow(threeDdata.rows()/2);

		iStart = Math.max(0, (int)((_WE.tStart - ssInfo.tStartData)/_WE.dt));
		iEnd = Math.min(longdata.size(), (int)((_WE.tEnd - ssInfo.tStartData)/_WE.dt)) -1;
		winSize = iEnd-iStart+1;
		if (winSize < WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION)
		  {
			halfWinPad = (int)Math.ceil(WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION/2.0);
			if (((iStart - halfWinPad) >=0)
				&& ((iEnd + halfWinPad) < longdata.size()))
			  {
				// (1) try to expand up and down	
				iStart -= halfWinPad;
				iEnd += halfWinPad;
			  }
			else if ((iStart - WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION) >=0)
			  {
				// (2) try to expand upwards
				iStart -= WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION;
			  } 
			else if ((iEnd + WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION) < longdata.size())
			  {
				// (3) try to expand down
				iEnd += WaveletExtractor.MIN_WINDOW_SAMPLES_FOR_SPEC_ESTIMATION;
			  }
			else
			  {
				// (4) desperate. Use what we can
				iStart = 0;
				iEnd = longdata.size()-1;
				_WE.EH.warn("WARNING! spectrum of seismic data %s is very unreliable from short time gate!!!", ssInfo.filename);
			  }
		  }

		// extract chunk in window of interest
		data = longdata.viewPart(iStart,iEnd-iStart+1).toArray();

		// get typical amplitude in window of interest (useful in other stuff)
		// make this the peak value in window required
		// use RMS amplitude 		
		for (i = 0; i < data.length; i++)
		  {
			Arms += data[i]*data[i];
			typicalPeakSeismicAmplitude = Math.max(typicalPeakSeismicAmplitude,Math.abs(data[i]));
		  }
		Nd += data.length;
		
		// find suggested sparse sampling rate		
		u = Wavelet.estimateWaveletIndependentDT(data,_WE.dt,_WE.verbosity,_WE.EH);
		dt_sparse += u;
		//_WE.EH.println("DT estimate data %s %f%n",ssInfo.filename,u);
	  }
	typicalRMSSeismicAmplitude = Math.sqrt(Arms/Nd);
	if (typicalRMSSeismicAmplitude == 0.0)
	  {
		_WE.EH.warn("Seismic minicube appears to have zero amplitudes! Extraction will be nonsense!");
	  }
	if (_WE.verbosity>2) _WE.EH.println("Typical amplitude taken as: %3.4g", typicalRMSSeismicAmplitude);

	// arithmetic mean for this
	dt_sparse/=Nwa;
  
	// compare to prior wavelet if supplied
	if (_WE.useWaveletPrior)
	  {
		double ratio = Math.abs(dt_sparse/_WE.waveletPrior.dtCoarse_prior);
		if (ratio < 0.5 || ratio > 2.0)
		  {
			_WE.EH.warn("Wavelet prior seems to have been derived from seismic with very different bandwidth!");
			_WE.EH.warn("Check this! Proceeding only blindly anyway...");
		  }
		dt_sparse = _WE.waveletPrior.dtCoarse_prior;
		_WE.EH.println("Wavelet prior dt_coarse ratio: "+ratio);	
	  }

	if (dt_sparse<_WE.dt)
	  {
		_WE.EH.warn("WARNING! spectrum of seismic data doesn't yields estimated sampling dt < dt of data");
		_WE.EH.warn("Check that your seismic data is not undersampled!");
		_WE.EH.warn("Setting likelihood-acquisition dt to that of data and continuing on blindly...");
		// Fix this and hope XXXX
		dt_sparse = _WE.dt;
	  }
	if (_WE.verbosity>2) 
	  {
		_WE.EH.println("Underlying coarse sample spacing for wavelet = %3.4g ms",1e3*dt_sparse);
		if(_WE.bbTail)
		  {
			_WE.EH.println("Applying broadband correlated-tail model for coefficients outside [%.2f,%.2f]ms",-_WE.bbTail_T,_WE.bbTail_T);
		  }
	  }

	// Set up f_peak for this class, which is used in other contexts
	f_peak = Wavelet.SAMPLE_TIME_MAGIC_CONSTANT/dt_sparse;

	// so within the time window allowed, these are the possible wavelets
	if (!_WE.useWaveletPrior)
	  {
		waveletCoarseIndexMin = -(int)Math.abs(_WE.tPrecursorMin/dt_sparse);
		waveletCoarseIndexMax = (int)Math.abs(_WE.tCodaMax/dt_sparse);
	  }
	else
	  {
		double firstCoarseSample = _WE.waveletPrior.f1_prior/_WE.waveletPrior.dtCoarse_prior-1 ;

		// override wavelet range thus....
		waveletCoarseIndexMin = -(int)Math.floor(-firstCoarseSample);
		waveletCoarseIndexMax = waveletCoarseIndexMin + _WE.waveletPrior.amp_means.size() + 1;

		if (_WE.verbosity>2) 
		  {
			_WE.EH.println("Use of wavelet prior imposes new wavelet span: %8.4g, %8.4g\n",
						   dt_sparse * waveletCoarseIndexMin,
						   dt_sparse * waveletCoarseIndexMax);
		  }
		if (waveletCoarseIndexMin*dt_sparse < _WE.tPrecursorMin)
		  {
			_WE.EH.die_err("Wavelet precursor time "+_WE.tPrecursorMin+
						   " not large enough to cope with supplied wavelet using --prior");
		  }
		if (waveletCoarseIndexMax*dt_sparse > _WE.tCodaMax)
		  {
			_WE.EH.die_err("Wavelet coda time "+_WE.tCodaMax+
						   " not large enough to cope with supplied wavelet using --prior");
		  }

		// c.f.
		// waveletCoarseIndexMin = -(int)Math.abs(_WE.tPrecursorMin/dt_sparse);
		// waveletCoarseIndexMax = (int)Math.abs(_WE.tCodaMax/dt_sparse);
	  }
 
	// ok, so here's the wavelet array

	// but we may just do the longest wavelet
	if (_WE.computeForLongestWaveletOnly)
	  {
		size = 1;
	  }
	else if (_WE.computeForCenteredWaveletsOnly)
	  {
		// count size over all possible combinations...
		// old...
		Imax = Math.min(Math.abs(waveletCoarseIndexMin),waveletCoarseIndexMax);
		size = Imax - halfLoopCoarseIndex + 1;
		// ...old
		// new
		for (i = waveletCoarseIndexMin, j=waveletCoarseIndexMax,size=0;; i++,j--)
		  {
			if (((j-i) >= (2*halfLoopCoarseIndex+1)) && (i<0) && (j>0)) size++;
			else break;
		  }
		if (size<1) _WE.EH.die_err("Specified wavelet span is too short to be reasonable");
		// ...new
	  }
	else if (_WE.computeForSelectedWaveletOnly) 
	  {
		// count size over all possible combinations, as per centered mode
		Imax = Math.min(Math.abs(waveletCoarseIndexMin),waveletCoarseIndexMax);
		// but just use named mode...
		size = 1;
	  }
	else
	  {
		// comprehensive mode (all possible spans)
		// size = (-halfLoopCoarseIndex - waveletCoarseIndexMin + 1)  *(waveletCoarseIndexMax - halfLoopCoarseIndex+1);
		size=0;
		// comprehensive mode : same loop as below, just count.....(all possible spans)
		for (i = waveletCoarseIndexMin; i<= -halfLoopCoarseIndex; i++)
		  {
			// ...and loop over wavelet end index
			for (j=waveletCoarseIndexMax; j>= halfLoopCoarseIndex; j--)
			  {
				if (j-i <= 1) continue;
				size++;
			  }
		  }
	  }

	// Anisotropy checks
	if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel) 
	  {
		// try also each model with anisotropy both on and off
		size*=2;
	  }

	// Auto-checkshot model-selection candidates

	int checkshotModelSelectSize = 1;
	for (i=0; i<_WE.m.nWells;i++) 
	  {	  
		if (_WE.wellInfo[i].checkshotAutoMake && _WE.wellInfo[i].checkshotAutoMakeModelSelect) 
		  {
			checkshotModelSelect = true; 
			if (checkshotModelSelectWells==null) checkshotModelSelectWells = new boolean[_WE.m.nWells];
			checkshotModelSelectWells[i] = true;
		  }
	  }
	if (checkshotModelSelect)
	  {
		int[] curDim = {-1};
		checkshot_model_select_index_array = new ArrayList<int[]>();
		int[] model_index = new int[_WE.m.nWells];
		populateCheckShotModels(model_index,checkshot_model_select_index_array,curDim);
		checkshot_model_select_index_array.trimToSize();
		// checkshot_model_select_index_array now has "hypercube" listing of checkshot models, so
		// k = (checkshot_model_select_index_array.get(i))[j] is the ith model, setting checkshot model k for well j
		checkshotModelSelectSize = checkshot_model_select_index_array.size();
		size *= checkshotModelSelectSize;
	  }


	modelVariables = new ModelVariables[size];

	// loop over possible wavelets, creating complete model for each.
	// but first we need a few other pieces
	// complete setup of common model properties
	_WE.m.setUpCommonProperties(typicalRMSSeismicAmplitude);

	// ...loop over wavelet start index
	if (_WE.computeForLongestWaveletOnly) 
	  {
		k=0;
		for (i_c = 0; i_c < checkshotModelSelectSize; i_c++)
		  {
			modelVariables[k++] = new ModelVariables(_WE,waveletCoarseIndexMin,waveletCoarseIndexMax,
													 _WE.dt,dt_sparse,_WE.doAnisotropicModel, i_c);
			if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel)
			  {
				// comparison model sans anisotropy
				modelVariables[k++] = new ModelVariables(_WE,waveletCoarseIndexMin,
														 waveletCoarseIndexMax,_WE.dt,dt_sparse, false, i_c);
			  }		
		  }
	  }
	else if (_WE.computeForCenteredWaveletsOnly)
	  {
		k=0;
		for (i_c = 0; i_c < checkshotModelSelectSize; i_c++)
		  {
			for (i = waveletCoarseIndexMin, j=waveletCoarseIndexMax;; i++,j--,k++)
			  {
				if (((j-i) < (2*halfLoopCoarseIndex+1)) || (i==0) || (j==0)) break;
				// _WE.EH.println("span "+i+" "+j);
				modelVariables[k] = new ModelVariables(_WE,i,j,_WE.dt,dt_sparse,_WE.doAnisotropicModel, i_c);
				if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel)
				  {
					k++;
					// comparison model sans anisotropy
					modelVariables[k] = new ModelVariables(_WE,i,j,_WE.dt,dt_sparse,false, i_c);
				  }	
			  }
		  }
	  }
	else if (_WE.computeForSelectedWaveletOnly)
	  {
		int k_w;
		k=0;
		for (i_c = 0; i_c < checkshotModelSelectSize; i_c++)
		  {
			for (i = -Imax,k_w=0; i<= -halfLoopCoarseIndex; i++,k++,k_w++)
			  {
				if (k_w == _WE.selectedWaveletModeNumber)
				  {
					modelVariables[0] = new ModelVariables(_WE,i,-i,_WE.dt,dt_sparse,_WE.doAnisotropicModel, i_c);
					if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel)
					  {
						// comparison model sans anisotropy
						modelVariables[1] = new ModelVariables(_WE,i,-i,_WE.dt,dt_sparse,false, i_c);
					  }	
				  }
			  }
		  }
		// check..
		if (modelVariables[0] == null)
		  {
			_WE.EH.die_err("Selected mode number "+_WE.selectedWaveletModeNumber+" doesn't exist");
		  }
	  }
	else
	  {
		// comprehensive mode (all possible spans)
		k=0;
		for (i_c = 0; i_c < checkshotModelSelectSize; i_c++)
		  {
			for (i = waveletCoarseIndexMin; i<= -halfLoopCoarseIndex; i++)
			  {
				// ...and loop over wavelet end index
				for (j=waveletCoarseIndexMax; j>= halfLoopCoarseIndex; j--)
				  {
					if (j-i<=1) continue;
					
					// Create distinct model for each wavelet configuration
					modelVariables[k] = new ModelVariables(_WE,i,j,_WE.dt,dt_sparse,_WE.doAnisotropicModel, i_c);
					k++;
					if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel)
					  {
						k++;
						// comparison model sans anisotropy
						modelVariables[k] = new ModelVariables(_WE,i,j,_WE.dt,dt_sparse,false, i_c);
					  }	
				  }
			  }
		  }
	  }

	// set up appropriate estimation interval
	setUpEstimationIntervalForAllModels();

  }

  /*
	Recursive method to enumerate checkshot models in autoGenerate from sonic-logs mode.
	curDim is effectively an integer ptr: initialise by setting to -1.
  */
  private void populateCheckShotModels(int[] model_index, ArrayList checkshot_model_select_index_array, int[] curDim)
  {
	int i_well,nC=0;
	curDim[0]++;


	i_well = curDim[0];

	if (_WE.wellInfo[i_well].checkshotAutoMake && _WE.wellInfo[i_well].checkshotAutoMakeModelSelect) 
	  {
		nC = _WE.chkShotandMarkerData[i_well].autoCheckshotInfo.segIndexValuesArray.length;
	  }
	for (int i=0;i<nC;i++)
	  {
		model_index[i_well] = i;
		if (i_well < (_WE.m.nWells-1)) 
		  {
			populateCheckShotModels(model_index,checkshot_model_select_index_array,curDim);
		  }
		else
		  {
			// store in array
			int[] index_array = model_index.clone();
			for (int k=0; k< _WE.m.nWells; k++)	index_array[k] = model_index[k];
			checkshot_model_select_index_array.add(index_array);
			//System.arraycopy(model_index,0,checkshot_model_select_index_array[checkshot_model_select_index_array_ctr],0,dim);
			//checkshot_model_select_index_array_ctr++;
		  }
	  }
	curDim[0]--;
  } 


  /**
	 loop over all possible wavelets to find modes in posterior
  */
  void findModes()
  {
	int i;
	ModeFinder finder;
	double f_opt;

	setLongestWaveletSpanVariables();
	ModelVariables _modelVariables;

	// loop over possible wavelet spans, locating modes for each wavelet
	for (modeNumber=0; modeNumber<modelVariables.length;modeNumber++)
	  {
		_modelVariables = modelVariables[modeNumber];

		// set shared objects to see current model
		_WE.m.setSharedVariables(_modelVariables);

		// update checkshot if needed
		_modelVariables.updateCheckshotForModelSelection();
		
		// create a mode finder for the current model
		finder = new ModeFinder(_WE,modelVariables[modeNumber]);

		// find local optimum starting from defaults
		finder.findOptimalSolution();

		if(false)
		  { 
			// experimental
			finder.modelVariables.getInversionPredictionWavelet_exptl();
		  }

		// Show some X sections (post-minimiser) if desired
		if (_WE.printXsections)
		  {		
			finder.debugPrintModeXsections("posteriorXsection.span."+modeNumber+".",false,false);
		  }
		// Principal axes X-sections (solution surface usually has oblique narrow valleys)
		// Maybe do this only on dodgy ones (false ->true)
		if (_WE.printPrincipalXsections)
		  {
			finder.debugPrintModePrincipalXsections("posteriorPrincipalXsection.span."+modeNumber+".",false);
		  }
		// Show some dodgy X sections (post-minimiser) if desired
		if (_WE.printDodgyXsections)
		  {		
			finder.debugPrintModeXsections("posteriorDodgyXsection.span."+modeNumber+".",true,false);
		  }

		// dump some synthetics if required
		if (_WE.dumpMLsynthetics)
		  {
			_WE.convModel.writeOutTrueAndSyntheticSUset("MLsynth_plus_true_seis.mode."+modeNumber);
		  }
		if (_WE.dumpMLwavelets || _WE.dumpMLparameters)
		  {
			// need wavelet if dumpMLgraphingfile() called
			_WE.m.wavelet.writeOutSUWavelet("MLwavelet.mode."+modeNumber,true,true,tMinLongestWavelet,tMaxLongestWavelet,
											false,0,_modelVariables.haveTseismicErrors);
		  }
		if (_WE.dumpMLparameters)
		  {
			finder.modelVariables.dumpMLparameters("MLparameters.mode."+modeNumber+".txt");
			finder.modelVariables.dumpMLTimeToDepthMap("MLtime-to-depth-map."+modeNumber);
			finder.modelVariables.dumpMLgraphingfile("MLgraphing-file.mode."+modeNumber);
		  }
		_WE.EH.println("Wavelet centroid time: %10.5f, peak arrival: %10.5f",
					   _WE.m.wavelet.getMeanWaveletTime(false),
					   _WE.m.wavelet.getPeakTime(false,finder.modelVariables.constrainPeakArrivalSign,
													 finder.modelVariables.peakArrivalPositive));
		//		_WE.m.wavelet.getPeakTime(true);
		if (finder.problemHasDegenerateNoise) 
		  {
			_WE.EH.println("WARNING: Problem has inferred a perfect fit (zero noise): "
						   +"you have too little data to make a meaningful well tie");
			_WE.EH.println("       : Zero noise estimates make uncertainty parameters meaningless");
		  }
	  }


	// finally, order the CDF etc
	computeModelProbabilities();

  }

  private void setUpEstimationIntervalForAllModels()
  {
	int i;
	double tStartMax = -Double.MAX_VALUE;
	double tEndMin = Double.MAX_VALUE;

	if (_WE.convModel.RIGOROUSLY_EXCLUDE_EDGE_REGIONS)
	  {
		// get smallest interval
		maxCodaMargin = maxPrecursorMargin =0.0;
		for (i=0; i<modelVariables.length;i++)
		  {
			tStartMax = Math.max(tStartMax, modelVariables[i].tStart);

			tEndMin = Math.min(tEndMin, modelVariables[i].tEnd);
			maxCodaMargin = Math.max(maxCodaMargin, 
										 modelVariables[i].wavelet.codaDuration + _WE.dt);
			maxPrecursorMargin = Math.max(maxPrecursorMargin, 
											modelVariables[i].wavelet.precursorDuration + _WE.dt);
		  }
		// Check this still leaves a useful estimation interval
		// I think we need at least one data point...
		if (tEndMin - tStartMax < _WE.dt)
		  {
			_WE.EH.println("Error: Wavelet estimation interval is too narrow");
			if (_WE.convModel.RIGOROUSLY_EXCLUDE_EDGE_REGIONS)
			  {
				_WE.EH.println("Remember the window is narrowed by the largest coda and precursor length "+
							   "of the wavelets you examine");
			  }
			_WE.EH.die_err("");
		  }
		// Now hard write these into the modelVariables objects
		for (i=0; i<modelVariables.length;i++)
		  {
			_WE.convModel.tStartGlobal = tStartMax;
			_WE.convModel.tEndGlobal = tEndMin;
		  }
	  }
  }











  void setLongestWaveletSpanVariables()
  {
	// Longest wavelet is first model in array
	_WE.m.setSharedVariables(modelVariables[0]);

	tMinLongestWavelet = -_WE.m.wavelet.getDTprecursor();
	tMaxLongestWavelet = _WE.m.wavelet.getDTcoda();
	
	// but far wavelets can be even longer if stretching permitted
	// here we have an upper bound
	if (_WE.useStretchedModel)
 	  {
 		tMinLongestWavelet = -_WE.m.wavelet.precursorStretchedDuration;
 		tMaxLongestWavelet = _WE.m.wavelet.codaStretchedDuration;
	  }
  }




  void sampleWavelets()
  {
	// There are several ways one might want to do this.
	// here we implement the simplest: an approximate independence sampler
	sampleWaveletsUsingIndependenceSampler();
  }

  private void computeModelProbabilities()
  {
	int i,j,k;
	double mostLikelyWaveletWeight;
	double norm,sum;


	// loop over modes to establish normalisation
	double weight_norm = -Double.MAX_VALUE;
	for (k=0; k<modelVariables.length; k++)
	  {
		weight_norm = Math.max(modelVariables[k].logModeWeight, weight_norm);
	  }

	
	// normalise relative 'posterior' probabilities to be less than 1
	double max_LL = -Double.MAX_VALUE;
	for (k=0; k< modelVariables.length; k++) 
	  {
		_WE.m.setSharedVariables(modelVariables[k]);
		// normalise weight so most probable state has probability of 1
		modelVariables[k].logModeWeight -= weight_norm;
		// retain record of most likely mode in posterior
		max_LL = Math.max(max_LL,modelVariables[k].logModeWeight);
		if (_WE.verbosity>3) _WE.EH.println("Initial mode list entry: "+k+" weight = "+modelVariables[k].logModeWeight );
	  }

	// drop modes that seem impossible wrt best mode, where
	// the best mode has the negloglikelihood max_LL
	/*
	  for (k=kk=0; k<modelVariables.length; k++,kk++) 
	  {
		if (((LocalMode)localModes.get(k)).weight > least_nLL + MODE_NLL_ELIMINATION_MARGIN)
		  {
			if (_verbosity>3) _WE.EH.println("Dropping mode "+kk);
			localModes.remove(k--);
		  }
	  }
	*/

	for (k=0,norm=0.0; k< modelVariables.length; k++) 
	  {
		norm += Math.exp(modelVariables[k].logModeWeight);
	  }
	// convert weight probabilities to real ones now
	sum=0.0;
	if (_WE.verbosity>2) _WE.EH.printf("%nNormalised local modes summary:%n");
	mostLikelyWaveletIndex=0;
	mostLikelyWaveletWeight=0;
	boolean isAmbiguous = false;
	for (k=0; k< modelVariables.length; k++) 
	  {
		// mode.weight's are henceforth absolute +ve probabilities
		modelVariables[k].weight = Math.exp(modelVariables[k].logModeWeight)/norm;
		sum += modelVariables[k].weight;
 		// compute both CDF's
 		modelVariables[k].CDF = sum;
		// keep track of ML wavelet
		if (modelVariables[k].weight > mostLikelyWaveletWeight) 
		  {
			mostLikelyWaveletWeight = modelVariables[k].weight; 
			mostLikelyWaveletIndex = k;
		  }
		if (modelVariables[k].modeHasPeakAmbiguity) isAmbiguous = true;
		Wavelet wavelet = modelVariables[k].wavelet;
		_WE.EH.println("Normalised model entry: mode number %2d "
					   +" weight= %-14.6g"
					   +" CDF= %-14.6g"
					   + (modelVariables[k].modeHasPeakAmbiguity ? "(ambiguous peak!)" : "                 ")
					   +" wavelet span: %-10.5g %-10.5g",
					   k, modelVariables[k].weight, modelVariables[k].CDF,
					   wavelet.firstCoarseSample*wavelet.dtCoarse,
					   wavelet.finalCoarseSample*wavelet.dtCoarse);
	  }
	_WE.EH.println("");
	if (isAmbiguous)
	  {
		_WE.EH.println("Modes markes with (ambiguous peak!) are non-interpretable. "
					   +"If these appear, ");
		_WE.EH.println("I suggest trying again, using --constrain-peak-arrival-sign with specified sign.");
		_WE.EH.println("");
	  }
	// set last mode CDF to be 1.0 exactly (precaution against roundoff)
	modelVariables[k-1].CDF = 1.0;
  }

  /**
	 Print out everything required about the wavelet extraction modes etc
  */
  void outputPosteriorEntities()
  {
	int i, indx;
	double mean,stdDev;

	// DUMP stuff to do with most likely model
	ModelVariables mostLikelyMV = modelVariables[mostLikelyWaveletIndex];
	_WE.EH.println("Most likely mode = "+mostLikelyWaveletIndex);

	// update checkshot if needed
	mostLikelyMV.updateCheckshotForModelSelection();

	// DUMP MOST likely overall wavelet
	if (_WE.dumpMLwavelets || _WE.dumpMLparameters)
	  {
		// wavelets needed as dumpMLparameters below calls dumpMLgraphingfile()
		_WE.m.setSharedVariables(mostLikelyMV);
		_WE.m.wavelet.writeOutSUWavelet("MLwavelet.MOST_LIKELY",true,false,0,0,false,0, mostLikelyMV.haveTseismicErrors);
		_WE.m.wavelet.computePhaseSpectrum("MLwavelet-phase-spectrum.MOST_LIKELY.txt",
										   mostLikelyMV.phaseConstraintIndxStart,
										   mostLikelyMV.phaseConstraintIndxEnd);

		// FIXME
		//_WE.m.wavelet.computeSpectrum("MLwavelet-pwr-spectrum.MOST_LIKELY.txt");

		// _WE.m.wavelet.debugPrintFineWavelet();
		// _WE.m.wavelet.getNearestPeakTime(_WE.peakArrivalTime);
		// _WE.m.wavelet.debugPrintCoarseWavelet();
	  }

	// Dump blocked logs
	if (_WE.dumpBlockedLogs)
	  {
		// must set insitu variables correctly first, for time converison
		mostLikelyMV.computeMLForwardModel();
		for (i=0; i<_WE.m.nWells; i++)
		  {
			String blockedLogFilename = (Util.fileBasename(_WE.wellInfo[i].logFilename)).replaceAll("(?i)\\.las|\\.txt","")
			                            +WaveletExtractor.BLOCKED_LOGS_SUFFIX;
			_WE.sLogData[i].writeOutBlockedLogFile(blockedLogFilename, _WE.m.timeToDepthMaps[i],
												   _WE.m.registrationErrorsAreCommon);
		  }
	  }
	// and most-likely overall parameters
	if (_WE.dumpMLparameters)
	  {
		mostLikelyMV.dumpMLparameters("MLparameters.MOST_LIKELY.txt");
		mostLikelyMV.dumpMLTimeToDepthMap("MLtime-to-depth-map.MOST_LIKELY");
		mostLikelyMV.dumpMLgraphingfile("MLgraphing-file.MOST_LIKELY");
	  }

	// Set and summarise most likely noise values
	if (_WE.verbosity>2) _WE.EH.println("Most likely noise characteristics:");
	for (i=0;i<_WE.m.nStacks;i++)
	  {
		if (!_WE.fixNoise)
		  {
			indx = ((Integer)mostLikelyMV.mapToIndex.get(_WE.m.seismic_noise[i])).intValue();
			mean = Math.exp(mostLikelyMV.mu_post.getQuick(indx));
			// set it, for convenience
			_WE.m.seismic_noise[i].setValue(mean);
			// approximate std deviation from linearisation
			stdDev = mean*Math.sqrt(mostLikelyMV.CtildeDiag[indx]);
		  } 
		else 
		  {
			mean = _WE.m.seismic_noise[i].value;
			stdDev=0.0;
		  }
		if (_WE.verbosity>2) 
		  {
			_WE.EH.println("Stack %d most-likely noise = %-10g, "
						   +"std-deviation %-10g, fraction of RMS amplitude = %-10g",
						   i, mean, stdDev, mean/_WE.multiModelSetUp.typicalRMSSeismicAmplitude);
		  }
	  }

  }
	
  /**
	 Perform Monte-Carlo studies here if required
  */
  void performMCcalculations()
  {
	// Bootstrapping 
	if (_WE.monteCarlo)
	  {
		DataBootstrapper bootstrapper = new DataBootstrapper(_WE);
		bootstrapper.performMCBootstrapOnSpecifiedMode(mostLikelyWaveletIndex);
		// Due to possible bootstrapping, the in-situ ML wavelet 
		// can no longer be assumed to be the ML wavelet from here......
	  }
  }



  private void sampleWaveletsUsingIndependenceSampler()
  {
	int i,j,k,n,nG,indx;
	double cdf;
	double mean,stdDev,val,stretch;
	ModelVariables sampleMV;

	// sample from approximate CDF

	setLongestWaveletSpanVariables();


	for (i=0; i< _WE.nSamples; i++) 
	  {
		// see writeOutStochasticProperties() below for interpretation (function not yet debugged)
		if (i==0) System.out.println("TODO: patch up writeOutStochasticProperties() for SU velocity realisations");
		for (k=0; k<MAX_NO_SAMPLE_ATTEMPTS; k++)
		  {
			cdf = Ugen.nextDouble();
			// find associated local mode
			for (k=0;k<modelVariables.length;k++)
			  {
				if (cdf <= modelVariables[k].CDF) break;
			  }
			// ok now we have chosen mode k
			_WE.EH.println("sampler, chose mode:: "+k);
			// set the wavelet etc
			sampleMV = modelVariables[k];
			_WE.m.setSharedVariables(sampleMV);
			nG = sampleMV.C_prior.rows();
			n = sampleMV.L_post.rows();

			DoubleMatrix1D w = new DenseDoubleMatrix1D(n);
			DoubleMatrix1D tmp = w.like();
			for (j=0;j<n;j++) w.setQuick(j,Ngen.nextDouble());
			// a realisation is given by P.L_post^T.w + mean
			// form w = L_post^T.w
			SeqBlas.seqBlas.dtrmv(false, true, false, sampleMV.L_post, w);
			for (j=0;j<n;j++) 
			  {
				// undo permutation & add back mean
				indx = sampleMV.L_Permutation[j];
				tmp.setQuick(indx, w.getQuick(j)+sampleMV.mu_post.getQuick(indx));
				// realisation is now in tmp. Write it into in-situ fields
				Property prop = (Property)sampleMV.mapFromIndex.get(indx);
				// restore log(noise) to actual noise
				val = tmp.getQuick(indx);
				if (indx >= sampleMV.nGaussian)
				  {
					// a log-transformed noise parameter
					val = Math.exp(val);
				  }
				prop.setValue(val);
			  }
			// noise ?
			// _WE.EH.println("SIGMA: "+_WE.m.seismic_noise[0].value);
			// update time-depth map
			_WE.m.updateTimeToDepthMaps();
			if (_WE.m.timesAreDecentlyOrdered)
			  {
				// time to depth map was legal, so continue on merrily
				break;
			  }
		  }
	  
		if (k>=MAX_NO_SAMPLE_ATTEMPTS)
		  {
			_WE.EH.warn("bad time ordering");
			_WE.EH.die_err("Couldn't get a sane t2d map in the sampling: posterior is wierd. Giving up");
		  }
		// update Wavelet
		if (_WE.useStretchedModel)
		  {
			// remap wav. stretch parameters
			for (k=0; k< (_WE.m.nStacks-1); k++)
			  {
				if (_WE.m.farStackWaveletTimeStretchMapper!=null) 
				  {
					_WE.m.farWavStretch[k] = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(_WE.m.farStackWaveletTimeStretch[k].value);
				  }
				else
				  {
					_WE.m.farWavStretch[k] = _WE.m.farStackWaveletTimeStretch[k].value;
				  }
				_WE.m.farWavAmp[k] = _WE.m.farStackWaveletAmplitude[k].value;
			  }
		  }
		// update Wavelet
		_WE.m.wavelet.updateWavelet();

		// compute new reflectivity arrays
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

		// make the mark header field go from 1...N
		_WE.m.wavelet.writeOutSUWavelet(_WE.wavelet_realisations_filename,false,true,
										tMinLongestWavelet,tMaxLongestWavelet,(i>0),i+1,false);


		// _WE.m.computeCheckshotIntervalVelocities();
		_WE.m.writeOutStochasticProperties(i);
	  }

  }
}
