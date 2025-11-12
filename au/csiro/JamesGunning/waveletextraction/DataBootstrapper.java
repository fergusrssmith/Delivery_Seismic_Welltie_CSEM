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

import cern.colt.matrix.*;
//import cern.colt.matrix.impl.*;
import cern.jet.random.*;
import cern.colt.matrix.linalg.*;
import cern.jet.random.engine.MersenneTwister;
import cern.colt.list.*;
import cern.jet.stat.*;
import hep.aida.*;
import hep.aida.ref.*;

import au.csiro.JamesGunning.lib.Numerical.*;
import au.csiro.JamesGunning.lib.util.*;



/**
   Class to perform an ensemble of extractions, for a specified wavelet mode, 
   using seismic data sets synthetically generated using the ML parameters extracted
   in the usual way. Statistics on the parameter ensemble can then be be formed as usual.
   (A MC way of exploring parameter uncertainty). On instantiation, assumes the in-situ
   seismic noise levels are the ML estimates.

   @author James.Gunning@csiro.au
   @version 1.0, 12/11/2004
*/

public class DataBootstrapper
{
  WaveletExtractor _WE;

  public static final int N_MC_REALISATIONS = 30;

  ModelVariables modelVariables;
  AbstractDistribution _Ngen;
  ModeFinder finder;

  double[] seismic_noise; // noise levels per stack. ML values on instantiation stored here
  double[] sigToNoise; 

  /**
	 Members describing the filter function
  */
  private double filterXfactor;
  private double filterYfactor;
  private double filterTfactor;
  private double filterNorm;
  private double filterScale;

  /**
	 MonteCarlo stats members
  */
  private int i_MC=0;
  private DoubleArrayList[] dataLists;
  private String[] propNames;
  private DoubleArrayList[] seismicQuantiles;
  private DoubleArrayList centiles;
  /** 
	  Number of quantiles used in CDF matching the fake seismic minicubes
   */
  public static int N_QUANTILES = 100;


  /** 
	  Backup elements for the parts of ModelVariables that this corrupts, and then
	  must restore on completion
  */
  private DoubleMatrix2D C_prior;
  private DoubleMatrix2D L_post;
  private int[] L_Permutation;
  private double logModeWeight,halfLogDetCeff,CDF,weight;
  private DoubleMatrix1D mu_post;
  private double[] CtildeInvDiag;


  /**
	 Constructor. Must be called when maximum a posteriori model is in-place
	 in the relevant objects, modelVariables[], multiModelSetUp etc.
  */
  public DataBootstrapper(WaveletExtractor WE)
  {
	this._WE = WE;
	_Ngen = _WE.multiModelSetUp.Ngen;

	seismic_noise = new double[_WE.m.nStacks];
	sigToNoise = new double[_WE.m.nStacks];

	for (int i=0; i<_WE.m.nStacks; i++) seismic_noise[i] = _WE.m.seismic_noise[i].value;
  }


  void performMCBootstrapOnSpecifiedMode(int modeNumber)
  {
	if (_WE.verbosity>2) 
	  {
		_WE.EH.printf("%nPerforming MC study for fake seismic minicubes...%n");
	  }
	modelVariables = _WE.multiModelSetUp.modelVariables[modeNumber];

	saveModelVariables(modelVariables);

	_WE.multiModelSetUp.setLongestWaveletSpanVariables();

	// set shared objects to see current model
	_WE.m.setSharedVariables(modelVariables);
		
	// create a mode finder for the current model
	finder = new ModeFinder(_WE,modelVariables);
	
	initialiseStats();
	
	initialiseFakeSeismicGeneration();

	for (i_MC=0; i_MC<N_MC_REALISATIONS; i_MC++)
	  {
		// generate Fake seismic data
		generateAllFakeSeismicStacks();

		// find local optimum starting from defaults
		// reset it_count so amplitudeCache resets
		finder.resetCache();
		finder.findOptimalSolution();
		
		

		// Accumulate requisite statistics from MAP solution
		accumulateStats();

		// Dump some useful output
		if (_WE.dumpMLsynthetics)
		  {
			_WE.convModel.writeOutTrueAndSyntheticSUset("MC_MLsynth_plus_true_seis."+i_MC);
		  }

		if (_WE.dumpMLwavelets || _WE.dumpMLparameters)
		  {
			// need wavelet if dumpMLgraphingfile() called
			_WE.m.wavelet.writeOutSUWavelet("MC_MLwavelet."+i_MC,true,false,0,0,false,0,false);
		  }
		if (_WE.dumpMLparameters)
		  {			
			finder.modelVariables.dumpMLparameters("MC_MLparameters."+i_MC+".txt");
			finder.modelVariables.dumpMLgraphingfile("MC_MLgraphing-file."+i_MC);
		  }
	  }
	finaliseFakeSeismicGeneration();
	computeMCstatsAndDumpToFile("MC_ParameterSummary.txt");

	restoreModelVariables(modelVariables);
  }



  /**
	 Makes copies of the seismic minicubes & obtain histogram data for
	 preservation of well-amplitude univariate stats
  */
  void initialiseFakeSeismicGeneration()
  {
	int i,j,k;
	double[] _data;
	DoubleArrayList dataList;

	seismicQuantiles = new DoubleArrayList[_WE.m.nStacks];

	// standard list of centiles 0, 0.01,....1
    centiles = new DoubleArrayList(N_QUANTILES);
	for (i=0;i<N_QUANTILES;i++) centiles.add(Math.min(1.0,i/(double)(N_QUANTILES-1)));

	for (j=0; j<_WE.m.nStacks; j++)
	  {
		// rough guess of size needed
		dataList = new DoubleArrayList(_WE.m.nWells * _WE.convModel.S_true[0][0].length);
		for (i=0;i<_WE.m.nWells;i++)
		  {
			_WE.seismicInfo[i][j].dataBk = new 
			  au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D(
																	  _WE.seismicInfo[i][j].data);
			// Form histogram from True traces at wells. We assume last 
			// run was with well position in ML point.
			// Ignore zeros as they are almost certainly from the padding at missing log data
			_data = _WE.convModel.S_true[i][j];
			for (k=0;k<_data.length;k++) if (_data[k]!=0.0) dataList.add(_data[k]);
		  }
		// complete quantile calculation from actual amplitude data at wells
		dataList.trimToSize();
		if (_WE.verbosity>2) 
		  {
			// printHistogram(dataList.elements(),"MC_SeisWellHistogram."+j+".txt",20);
			_WE.EH.println("Stack %d well amplitudes RMS = %5g%n",j,
						   Descriptive.rms(dataList.size(),
										   Descriptive.sumOfSquares(dataList)));
		  }
		dataList.sort();
		// now the quantiles of the seismic amps at the wells
		seismicQuantiles[j] = Descriptive.quantiles(dataList, centiles);
		seismicQuantiles[j].trimToSize();
		dataList = null;
	  }
  }

  /**
	 restore initial seismic minicubes
  */
  void finaliseFakeSeismicGeneration()
  {
	int i,j;
	for (i=0;i<_WE.m.nWells;i++)
	  {
		for (j=0;j<_WE.m.nStacks; j++)
		  {
			_WE.seismicInfo[i][j].data = _WE.seismicInfo[i][j].dataBk;		  
		  }
	  }
  }


  void generateAllFakeSeismicStacks()
  {
	int i,j,k,ii,jj;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D data;

	if (_WE.verbosity>2) _WE.EH.println("Generating fake seismic stacks...");
	for (i=0;i<_WE.m.nWells;i++)
	  {
		for (j=0;j<_WE.m.nStacks; j++)
		  {
			data = _WE.seismicInfo[i][j].data;
			generateFakeSeismic(data,j);
			if (false)
			  {
				_WE.EH.println("Field = %d x %d x %d%n",
							   data.slices(), data.rows(), data.columns());
				Util.debugPrintMatrixToFile("Field.txt",data);
			  }
		  }
	  }
  }

  /**
	 Routine to generate correlated seismic-ish data cube embedded in
	 dense Cern-style matrix

	 Uses FFT method with white frequencies filtered by 
	 filt(i,j,k,Filt), a filter function     
	 operating in frequency space.
  */
  void generateFakeSeismic(au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D data,int stackNumber)
  {
	double w1,w2;
	double vert_res = 5.0;
	double hor_res = 10.0;
	int ii,jj,kk,k2;
	int Nx,Ny,Nz,Nxpad,Nypad,Nzpad; 
	double[][] specf;
	double amp;
	au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D f=null;

	Nx = data.slices();
	Ny = data.rows();
	Nz = data.columns();
	
	double nTracesCorrLength= _WE.monteCarloLateralScale;


	/* pad to powers of 2 as required */
	/* Be fussy about xy padding, but not Z */
	/* Ricker (1-x^2) exp(-x^2/2) has to get to span x=4 to have died out enough */	
	Nxpad = FFT.getNearestPowerOfTwo(Nx+(int)(nTracesCorrLength));
	Nypad = FFT.getNearestPowerOfTwo(Ny+(int)(nTracesCorrLength));
	Nzpad = Math.max(2,FFT.getNearestPowerOfTwo(Nz)); // Nzpad must be at least 2
	
	/* fix filter parameters */

	filterXfactor = MathUtil.sqr(2*Math.PI*nTracesCorrLength/Nxpad);
	filterYfactor = MathUtil.sqr(2*Math.PI*nTracesCorrLength/Nypad);
	filterTfactor = MathUtil.sqr(2/(Nzpad*_WE.multiModelSetUp.f_peak * _WE.dt));

  
	filterNorm=2.0/(Nxpad*Nypad*Nzpad);
	
	try
	  {
		f = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D(Nxpad,Nypad,Nzpad);
	  }
	catch (OutOfMemoryError e)
	  {
		_WE.EH.die_err("Data minicubes are too big to perform bootstrapping on. Sorry :-(");
	  }
	specf = new double[Nxpad][2*Nypad];
	
	/* Set f = white noise. This loop treats all but the Nyquist frequencies */
	for (ii=0;ii<Math.max(1,Nxpad/2);ii++) {
	  for (jj=0;jj<Math.max(1,Nypad/2);jj++) {
		for (kk=0,k2=1;kk<Nzpad/2;kk++,k2+=2) {   /* k2=2*k */
		  /* ii=jj=kk==0 corresponds to DC */
		  w1 = _Ngen.nextDouble();
		  w2 = _Ngen.nextDouble();
		  f.setQuick(ii,jj,k2-1,w1);                     // Re
		  f.setQuick(ii,jj,k2,(ii==0 && jj==0)? 0 : w2); // Im (=0 for DC)
		  /* now for -ve non-zero frequencies */
		  if (jj>0) {
			f.setQuick(ii,Nypad-jj,k2-1,w1);
			f.setQuick(ii,Nypad-jj,k2,w2);
		  }
		  if (ii>0) {
			f.setQuick(Nxpad-ii,jj,k2-1,w1);
			f.setQuick(Nxpad-ii,jj,k2,w2);
			if (jj>0) {
			  f.setQuick(Nxpad-ii,Nypad-jj,k2-1,w1);
			  f.setQuick(Nxpad-ii,Nypad-jj,k2,-w2);
			}
		  }
		}
	  }
	}

	/* now do filtering. This loop treats all but the Nyquist frequencies */
	for (ii=0;ii<Math.max(1,Nxpad/2);ii++) {
	  for (jj=0;jj<Math.max(1,Nypad/2);jj++) {
		for (kk=0,k2=1;kk<Nzpad/2;kk++,k2+=2) {   /* k2=2*k */
		  /* i==1 or j==1 or k==1 corresponds to DC frequencies */
		  f.setQuick(ii,jj,k2,f.getQuick(ii,jj,k2) * filt(ii,jj,kk));                  /* Re +ve freqs i,j*/
		  f.setQuick(ii,jj,k2-1,f.getQuick(ii,jj,k2-1) * filt(ii,jj,kk));				 /* Im +ve freqs i,j*/
		  /* now for -ve non-zero frequencies */
		  if (jj>0) {
			f.setQuick(ii,Nypad-jj,k2,f.getQuick(ii,Nypad-jj,k2) * filt(ii,jj,kk));       /* Re -ve freqs j*/
			f.setQuick(ii,Nypad-jj,k2-1,f.getQuick(ii,Nypad-jj,k2-1) * filt(ii,jj,kk));       /* Im -ve freqs j*/
		  }
		  if (ii>0) {
			f.setQuick(Nxpad-ii,jj,k2,f.getQuick(Nxpad-ii,jj,k2) * filt(ii,jj,kk));          /* Re -ve freqs i*/
			f.setQuick(Nxpad-ii,jj,k2-1,f.getQuick(Nxpad-ii,jj,k2-1) * filt(ii,jj,kk));          /* Re -ve freqs i*/
			if (jj>0) {
			  f.setQuick(Nxpad-ii,Nypad-jj,k2,f.getQuick(Nxpad-ii,Nypad-jj,k2) * filt(ii,jj,kk));  /* Re -ve freqs */
			  f.setQuick(Nxpad-ii,Nypad-jj,k2-1,f.getQuick(Nxpad-ii,Nypad-jj,k2-1) * filt(ii,jj,kk));  /* Re -ve freqs */
			}
		  }
		}
	  }
	}
	/* now filter Nyquist components in the same way: Hopefully they are completely passed
	   by the radial filtering, and largely suppressed by the vertical filtering... */
	/* Nyquist i */
	if (Nxpad>=2) {
	  for (jj=0;jj<Math.max(1,Nypad/2);jj++) {
		for (kk=0,k2=1;kk<Nzpad/2;kk++,k2+=2) {   /* k2=2*k */
		  /* j==1,k==1 corresponds to DC frequencies */
		  f.setQuick(Nxpad/2 ,jj,k2,f.getQuick(Nxpad/2 ,jj,k2) * filt(Nxpad/2,jj,kk));				  /* Re +ve freqs i,j*/
		  f.setQuick(Nxpad/2 ,jj,k2-1,f.getQuick(Nxpad/2 ,jj,k2-1) * filt(Nxpad/2,jj,kk));				 /* Im +ve freqs i,j*/
		  /* now for -ve non-zero frequencies j */
		  if (jj>0) {
			f.setQuick(Nxpad/2 ,Nypad-jj,k2,f.getQuick(Nxpad/2 ,Nypad-jj,k2) * filt(Nxpad/2,jj,kk));	   /* Re -ve freqs j*/
			f.setQuick(Nxpad/2 ,Nypad-jj,k2-1,f.getQuick(Nxpad/2 ,Nypad-jj,k2-1) * filt(Nxpad/2,jj,kk));	   /* Im -ve freqs j*/
		  }
																						}
	  }
	}
	/* Nyquist j */
	if (Nypad>=2) {
	  for (ii=0;ii<Math.max(1,Nxpad/2);ii++) {
		for (kk=0,k2=1;kk<Nzpad/2;kk++,k2+=2) {   /* k2=2*k */
		  /* i==1,k==1 corresponds to DC frequencies */
		  f.setQuick(ii,Nypad/2 ,k2,f.getQuick(ii,Nypad/2 ,k2) * filt(ii,Nypad/2,kk));                /* Re +ve freqs i,j*/
		  f.setQuick(ii,Nypad/2 ,k2-1,f.getQuick(ii,Nypad/2 ,k2-1) * filt(ii,Nypad/2,kk));				/* Im +ve freqs i,j*/
		  /* now for -ve non-zero i frequencies */
		  if (ii>0) {
			f.setQuick(Nxpad-ii,Nypad/2 ,k2,f.getQuick(Nxpad-ii,Nypad/2 ,k2) * filt(ii,Nypad/2,kk));    /* Re -ve freqs i*/
			f.setQuick(Nxpad-ii,Nypad/2 ,k2-1,f.getQuick(Nxpad-ii,Nypad/2 ,k2-1) * filt(ii,Nypad/2,kk));  /* Re -ve freqs i*/
		  }
		}
	  }
	}
	/* Nyquist k */
	for (ii=0;ii<Math.max(1,Nxpad/2);ii++) {
	  for (jj=0;jj<Math.max(1,Nypad/2);jj++) {
		/* i==1,j==1, corresponds to DC frequencies */
		specf[ii][2*jj]*=filt(ii,jj,Nzpad/2);                  /* Re +ve freqs i,j*/
		specf[ii][2*jj+1]*=filt(ii,jj,Nzpad/2);				      /* Im +ve freqs i,j*/
		/* now for -ve non-zero frequencies */
		if (jj>0) {
		  specf[ii][2*(Nypad+1-jj)-2]*=filt(ii,jj,Nzpad/2);       /* Re -ve freqs j*/
		  specf[ii][2*(Nypad+1-jj)-1]*=filt(ii,jj,Nzpad/2);         /* Im -ve freqs j*/
		}
		if (ii>0) {
		  specf[Nxpad-ii][2*jj]*=filt(ii,jj,Nzpad/2);          /* Re -ve freqs i*/
		  specf[Nxpad-ii][2*jj+1]*=filt(ii,jj,Nzpad/2);            /* Re -ve freqs i*/
		  if (jj>0) {
			specf[Nxpad-ii][2*(Nypad+1-jj)-2]*=filt(ii,jj,Nzpad/2);/* Re -ve freqs */
			specf[Nxpad-ii][2*(Nypad+1-jj)-1]*=filt(ii,jj,Nzpad/2);  /* Re -ve freqs */
		  }
		}
	  }
	}
	
	/* now normalise the output */
	for(ii=0; ii<Nxpad; ii++) {
	  for(jj=0; jj<Nypad; jj++) {
		for(kk=0; kk<Nzpad; kk++) {
		  f.setQuick(ii,jj,kk,f.getQuick(ii,jj,kk) * filterNorm);   /* robust filter Normalisation */
		}
	  }
	  for(jj=0; jj<2*Nypad; jj++) {
		specf[ii][jj]*=filterNorm;      /* robust normalisation */
	  }
	}
	/* invert the FFT */
	FFT.rlft3(f,specf,-1); /* reverse ft */
	
	
	// Now map histogram of FFT simulation onto that of true data, for this stack

	double[] amps = f.viewElements().clone();
	DoubleArrayList fList = new DoubleArrayList(amps);
	fList.trimToSize();
	fList.sort();
	DoubleArrayList quantiles = Descriptive.quantiles(fList, centiles);
	quantiles.trimToSize();
	// now loop to restore elements
	LinearInterpolation interp = new LinearInterpolation(quantiles.elements(),seismicQuantiles[stackNumber].elements(),
														 N_QUANTILES,_WE.EH);
	  
	// set cube histogram to that of wells
	for(ii=0; ii<Nx; ii++) {
	  for(jj=0; jj<Ny; jj++) {
		for(kk=0; kk<Nz; kk++)
		  {
			amp = f.getQuick(ii,jj,kk);
			amp = interp.interpolate(amp,kk!=0);
			data.setQuick(ii,jj,kk,amp);
		  }
	  }
	} 

	if (_WE.verbosity>2) 
	  {
		// printHistogram(data,"MC_SynthWellHistogram."+i_MC+".txt",20);
		_WE.EH.println("Synthetic stack %d, run %d, well amplitudes RMS = %5g%n",stackNumber,i_MC,
					   Descriptive.rms(Nx*Ny*Nz,
									   Descriptive.sumOfSquares(new DoubleArrayList(data.viewElements()))));
	  }

	quantiles = null;
	amps = null;
	fList = null;
	interp=null;

	// clean up
	specf=null;
	f=null;
  }
  
  /**
	 Power filter H(k), where i,j,k are discretised frequencies
	 on the DFT lattice. Will generate a correlated image which looks
	 vaguely seismic like when scale factors are precomputed sensibly
  */
  
  private double filt(int i,int j,int k)
  {
	double f,arg1,arg2;

	arg1 = (filterTfactor*k*k);                      // FT of (1-t^2) exp[-t^2/2] is w^2 exp(-w^2/2) (Ricker)
	arg2 = (filterYfactor*j*j + filterXfactor*i*i);  // FT of exp[-(3 x)^2] is exp(-w^2/36)
	f = arg1*Math.exp(-arg1/2.0 - arg2/36.0);
	return f;
  }

  /**
	 Initialise data containers for properties of interest
  */
  void initialiseStats()
  {
	int i,ii,indx;
	Property prop;
	String name;

	dataLists = new DoubleArrayList[finder.modelVariables.nGaussian + ((!_WE.fixNoise)?_WE.m.nStacks:0)];
	for (i=0;i<dataLists.length;i++) dataLists[i] = new DoubleArrayList(N_MC_REALISATIONS);
	propNames = new String[dataLists.length];
	// include standard list of Gausian variables (wavelet coefficients etc)
	for (i=0; i<finder.modelVariables.nGaussian; i++)
	  {
		prop = (Property)finder.modelVariables.mapFromIndex.get(i);
		name = 	(String)finder.modelVariables.mapToName.get(prop);
		propNames[i] = name;
	  }
	if (!_WE.fixNoise)
	  {
		// add the noise bits in here too.
		double[] rms = _WE.convModel.getTrueSeisRMSatWells();
		for (ii=0;ii<_WE.m.nStacks;ii++,i++)
		  {
			indx = ((Integer)finder.modelVariables.mapToIndex.get(_WE.m.seismic_noise[ii])).intValue();
			name = (String)finder.modelVariables.mapToName.get(finder.modelVariables.mapFromIndex.get(indx));
			propNames[i] = new String("noise-to-signal-stack_"+ii);
			sigToNoise[ii] = seismic_noise[ii]/((rms[ii]>0)?rms[ii]:1);
			if (_WE.verbosity>2) 
			  {
				_WE.EH.println("Noise-to-RMS ratio for well %d, is %g",ii, sigToNoise[ii]);
			  }
		  }
	  }
  }

  /**
	 Accumulate stats for properties of interest
  */
  void accumulateStats()
  {
	int i,j,ii,indx;
	double MLval;
	Property prop;

	
	for (i=0; i<finder.modelVariables.nGaussian; i++)
	  {
		MLval = finder.modelVariables.mu_post.getQuick(i);
		prop = (Property)finder.modelVariables.mapFromIndex.get(i);
		// some properties may have been remapped. Retrieve true values here...
		if (_WE.m.farStackWaveletTimeStretchList.contains(prop) &&
			(_WE.m.farStackWaveletTimeStretchMapper != null)) 
		  {
			MLval = _WE.m.farStackWaveletTimeStretchMapper.getTrueValue(MLval);
		  }
		dataLists[i].add(MLval);
	  }

	if (!_WE.fixNoise)
	  {
		// get rms fake seismic along wells. Treat noise independently per stack
		double[] rms = _WE.convModel.getTrueSeisRMSatWells();
		double ratio;
		for (j=0;j<_WE.m.nStacks; j++,i++)
		  {
			// add the noise ratio
			indx = ((Integer)finder.modelVariables.mapToIndex.get(_WE.m.seismic_noise[j])).intValue();
			MLval = Math.exp(finder.modelVariables.mu_post.getQuick(indx));
			ratio = MLval/((rms[j]>0)?rms[j]:1);
			dataLists[i].add(ratio);
			if (_WE.verbosity>2) 
			  {
				_WE.EH.println("MC run: Stack %d; ML noise-amp = %6g, RMS-well-amp = %6g, Noise-to-RMS ratio = %6g",j, MLval,rms[j], ratio);
			  }
		  } 
		rms=null;
	  }
  } 


  /** 
	  Standard analysis of uncertainties in the parameter list. Mean and covariance
	  are computed
  */
  void computeMCstatsAndDumpToFile(String filename)
  {
	int i,j,n,N;
	double mean,percentile=0;
	PrintWriter out;

	// no. of properties
	N = dataLists.length;
	
	// compute stats using standard unbiased estimator of covariance
	// no. of samples
	dataLists[0].trimToSize();
	n = dataLists[0].size();
	// use Y in transposed view compared to Trimm
	DoubleMatrix2D Y = new cern.colt.matrix.impl.DenseDoubleMatrix2D(N,n);
	DoubleMatrix2D sampleCov;
	DoubleMatrix2D mu = new cern.colt.matrix.impl.DenseDoubleMatrix2D(1,N);
	for (j=0;j<N;j++)
	  {
		for (i=0;i<n;i++)
		  {
			Y.setQuick(j,i,dataLists[j].getQuick(i));
		  }
	  }
	// loop over parameters to get stats
	for (j=0;j<N;j++)
	  {
		// get mean
		mean = Descriptive.mean(dataLists[j]);
		mu.setQuick(0,j,mean);
	  }
	// form mu.mu^T
	sampleCov = (cern.colt.matrix.impl.DenseDoubleMatrix2D)mu.zMult(mu,null,-n/(n-1.0),0,true,false); 
	// and add in Y.Y^T/(n-1)
	Y.zMult(Y,sampleCov,1.0/(n-1),1,false,true); 

	// output
	// open ascii summary file
	_WE.EH.println("Writing summary of Monte Carlo parameters to "+filename);
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }

	out.printf("Reference: stack no., noise-to-signal level, percentile of MC stats%n");
	// noise-to-rms are the last members of the lists
	for (i=0,j=N-_WE.m.nStacks; i<_WE.m.nStacks; i++,j++)
	  {
		dataLists[j].sort();
		percentile = cern.jet.stat.Descriptive.quantileInverse(dataLists[j], sigToNoise[i]); 
		out.printf("Stack %d, %10.6g %10.2g%n",i, sigToNoise[i], percentile);
	  }
	out.printf("Monte-carlo run parameter means and std-deviations%n");
	for (i=0;i<N;i++)
	  {
		out.printf("%-30s %10.6g %10.6g%n",
				   propNames[i],mu.getQuick(0,i), Math.sqrt(sampleCov.getQuick(i,i)));
	  }
	
	out.close();
	_WE.EH.println("Wrote data to: " + filename);
  }


  // check histogram
  // If Nbins>0, use suggested number, else auto compute
  private void printHistogram(au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D data, String filename,int NBins)
  {
	printHistogram(data.viewElements(), filename,NBins);
  }
  
  



  // check histogram
  // If Nbins>0, use suggested number, else auto compute
  // Counts are normalised
  private void printHistogram(double[] data, String filename, int NBins)
  {
	int i,Nbins;
	DoubleArrayList fList;
	PrintWriter out;

	double[] amps = data.clone();
	fList = new DoubleArrayList(amps);
	// Histogram output
	fList.sort();
	if (NBins<=0)
	  {
		// let's take sqrt(2N) bins (same default rule as in |stat), but no less than 4
		Nbins = Math.max(4,(int)Math.rint(2*Math.sqrt(fList.size())));
	  }
	else Nbins = NBins;
	double min = fList.get(0);
	double max = fList.get(fList.size()-1);
	// form histogram with some% extra range at either end (say 10%)
	IHistogram1D h1 = new Histogram1D(null,Nbins,min-0.1*(max-min),max+0.1*(max-min));

	// populate histogram
	for (i=0; i < fList.size(); i++) h1.fill(fList.get(i));

	// Generate 2 column output: bin_centre, count_in_bin.
	IAxis ax = h1.xAxis();

	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + filename);
	  }

	double norm = h1.sumBinHeights();
	for (i=0; i < ax.bins(); i++) 
	  {
		out.printf("%10.8g %10.8g%n",ax.binCentre(i), h1.binHeight(i)/norm);
	  } 

	out.close();
	_WE.EH.println("Wrote data to: " + filename);	
  }

  private void saveModelVariables(ModelVariables modelVariables)
  {
	C_prior = modelVariables.C_prior.like().assign(modelVariables.C_prior);
	L_post = modelVariables.L_post.like().assign(modelVariables.L_post);
	L_Permutation = modelVariables.L_Permutation.clone();
	logModeWeight = modelVariables.logModeWeight;
	halfLogDetCeff = modelVariables.halfLogDetCeff;
	CDF = modelVariables.CDF;
	weight = modelVariables.weight;
	mu_post = modelVariables.mu_post.like().assign(modelVariables.mu_post);
	CtildeInvDiag = modelVariables.CtildeInvDiag.clone();
  }

  private void restoreModelVariables(ModelVariables modelVariables)
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
