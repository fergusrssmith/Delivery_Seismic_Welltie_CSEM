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
import java.util.*;

import cern.colt.matrix.linalg.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;

/**
   Basic class for predicting missing log values, especially shear
   
   @author James.Gunning@csiro.au
   @version 1.0, 8/3/07
*/

public final class LogPredictor
{
  WaveletExtractor _WE;
  Model _m;

  LogPredictor(WaveletExtractor WE)
  {
	_WE = WE;
	_m = _WE.m;
  }

  void formPredictedShearData()
  {
	int i,j,k,n=0,nIllegal=0;
	int Npredicted=0;
	GeneralWellInfo info;
	LogData logData;
	double f,SQRT2 = 1.4142136;
	boolean dumpStats=false;

	dumpStats = (_WE.verbosity>3);

	if (_WE.verbosity>0) _WE.EH.printf("%nRegression: fake shear log prediction...%n");

	// Count shear data
	for (i=0; i<_m.nWells; i++)
	  {
		info = _WE.wellInfo[i];
		logData = _WE.logData[i];
		if (info.wellLogDefinitions.shearTreatment != 0)
		  {
			for (j =0; j< logData.DT.length ; j++) 
			  {
				if ((logData.DTshear[j]>0) && logData.DT[j]>0) 
				  {
					// possible data. Count only if physical
					if (logData.DTshear[j] > SQRT2 * logData.DT[j])	
					  {
						n++;
					  }
					else nIllegal++;
				  } 
			  }
		  }
	  }



	if (n==0)
	  {
		if (_WE.verbosity>1)
		  {
			_WE.EH.println("Regression: no shear data available for regression. Using crude DTS = "
						   +FundamentalConstants.FAKE_VS_FACTOR+"x DT");
		  }
		// No shear data. crude approx
		for (i=0; i<_m.nWells; i++)
		  {
			info = _WE.wellInfo[i];
			logData = _WE.logData[i];
			if (info.wellLogDefinitions.shearTreatment == 0)
			  {	
				logData.DTshear = logData.DT.clone();
				Arrays.fill(logData.DTshear, FundamentalConstants.BAD_NUMBER);
				for (j =0; j< logData.DT.length ; j++) 
				  {
					if (logData.DT[j]>0) 
					  {				
						logData.DTshear[j] = FundamentalConstants.FAKE_VS_FACTOR * logData.DT[j];
					  }
					if (dumpStats)
					  {
						// temporarily mark by -ve value; restored after dump
						logData.DTshear[j] = -logData.DTshear[j];
					  }
				  }
			  }
		  }
	  }
	else
	  {
		if (nIllegal>10)
		  {
			_WE.EH.println("Regression: some unphysical shear sonic log observations - warning!");
		  }
		if (_WE.verbosity>1)
		  {
			_WE.EH.println("Regression: using available data"+
						   " in simple regression A+BxDT to infill missing shear values");
		  }
		// deal with possibly different sonic units
		int iRef=-1;
		double DTbar;
		double[] meanDT = new double[_m.nWells];
		double[] unitsConversion = new double[_m.nWells];
		Arrays.fill(unitsConversion,1.0);
		for (i=0; i<_m.nWells; i++)
		  {
			info = _WE.wellInfo[i];
			logData = _WE.logData[i];
			if (info.wellLogDefinitions.shearTreatment != 0)
			  {
				if (iRef==-1) iRef=i;  // index of first reference well
				DTbar=0;
				for (j=0,k=0; j< logData.DT.length ; j++) 
				  {
					if (logData.DT[j]>0) 
					  {
						DTbar += logData.DT[j];
						k++;
					  }
				  }
				if (k>0) DTbar/=k;
				meanDT[i] = DTbar;
			  }
		  }
		// compare to other wells
		for (i=iRef+1; i<_m.nWells; i++)
		  {
			info = _WE.wellInfo[i];
			logData = _WE.logData[i];
			if (info.wellLogDefinitions.shearTreatment != 0)
			  {
				// Compare mean DT statistics.
				// Should be within a factor of two if units are commensurable
				// There are two mixed units cases we ought to consider....
				if (((meanDT[i]*FundamentalConstants.METRES_PER_FOOT) < 2.0* meanDT[iRef]) && 
					((meanDT[i]*FundamentalConstants.METRES_PER_FOOT) > 0.5* meanDT[iRef]))
				  {
					// Case: well meanDT[i] in uS/m, meanDT[iRef] in uS/ft:
					unitsConversion[i] = FundamentalConstants.METRES_PER_FOOT;
					if (_WE.verbosity>1)
					  {
						_WE.EH.println("Regression: adjusting sonic units for wells %s (uS/m) and %s (uS/ft)%n",
									   _WE.wellInfo[i].well_name, _WE.wellInfo[iRef].well_name);
					  }
				  }
				else if (((meanDT[i]/FundamentalConstants.METRES_PER_FOOT) < 2.0* meanDT[iRef]) && 
						 ((meanDT[i]/FundamentalConstants.METRES_PER_FOOT) > 0.5* meanDT[iRef]))
				  {
					// Case: well meanDT[i] in uS/ft, meanDT[iRef] in uS/m:
					unitsConversion[i] = 1.0/FundamentalConstants.METRES_PER_FOOT;
					if (_WE.verbosity>1)
					  {
						_WE.EH.println("Regression: adjusting sonic units for wells %s (uS/ft) and %s (uS/m)%n",
									   _WE.wellInfo[i].well_name, _WE.wellInfo[iRef].well_name);
					  }
				  }	
			  }
		  }



		// Set up regression system	

		DoubleMatrix2D X = new DenseDoubleMatrix2D(n,2);
		DoubleMatrix1D y = new DenseDoubleMatrix1D(n);


		for (i=0,k=0; i<_m.nWells; i++)
		  {
			info = _WE.wellInfo[i];
			logData = _WE.logData[i];
			if (info.wellLogDefinitions.shearTreatment != 0)
			  {
				f = unitsConversion[i];
				for (j =0; j< logData.DT.length ; j++) 
				  {
					if ((logData.DTshear[j]>0) && logData.DT[j]>0) 
					  {
						// possible data. Count only if physical
						if (logData.DTshear[j] > SQRT2 * logData.DT[j])
						  {
							X.setQuick(k,0,1.0);
							X.setQuick(k,1,f * logData.DT[j]);
							y.setQuick(k,f * logData.DTshear[j]);
							k++;
						  }
					  }
				  }
			  }
		  }
		// Get LS estimate
		// (1) Form C = X^T.X
		DoubleMatrix2D C  = new DenseDoubleMatrix2D(2,2);
		DoubleMatrix1D rhs = new DenseDoubleMatrix1D(2);
		DoubleMatrix1D ans  = rhs.like();
		Blas blas = SeqBlas.seqBlas;

		blas.dgemm(true,false,1.0,X,X,1.0,C);
		// Form rhs = X^T.y 
		blas.dgemv(true,1.0,X,y,0,rhs);
		Algebra alg = new Algebra();
		C = alg.inverse(C); 
		blas.dgemv(false,1.0,C,rhs,0,ans);
		double A,B;
		A = ans.getQuick(0);
		B = ans.getQuick(1);
		if (_WE.verbosity>2)
		  {
			_WE.EH.printf("Fake seismic regression: DTS = %5g + %5g * DT%n",A,B);
			_WE.EH.printf("N.B. Set verbosity>3 to see regression dump, and XML attribute "+
						  "<shear_sonic_DT_keyword use=\"false\">%n" +
						  "to omit a particular well's shear data from regression if desired%n");
		  }



		// repopulate missing data



		for (i=0; i<_m.nWells; i++)
		  {
			info = _WE.wellInfo[i];
			logData = _WE.logData[i];
			f = unitsConversion[i];
			if (info.wellLogDefinitions.shearTreatment == 0)
			  {
				// allocate new data
				logData.DTshear = logData.DT.clone();
			  }
			// populate...
			for (j =0; j< logData.DT.length ; j++) 
			  {
				// re-estimate if either (1) data not used in regression, or
				// (2) log-data used, but this was a missing chunk (-ve values)
				if ((info.wellLogDefinitions.shearTreatment == 0)
					|| (logData.DTshear[j]<=0))
				  {
					logData.DTshear[j] = FundamentalConstants.BAD_NUMBER;
					if (logData.DT[j]>0) 
					  {
						// missing shear, Predict !!!
						Npredicted++;
						logData.DTshear[j] = A/f + B*logData.DT[j];
						if (logData.DTshear[j] <= SQRT2 * logData.DT[j])
						  {
							// safe bound if regression is wierd
							logData.DTshear[j] = (SQRT2+.01)*logData.DT[j];
						  }					
						if (dumpStats)
						  {
							// temporarily mark by -ve value; restored after dump
							logData.DTshear[j] = -logData.DTshear[j];
						  }
					  }
				  }				   	
			  }
		  }
	  }

	if (dumpStats)
	  {
		PrintWriter out=null;
		String regressionFilename = "shear_regression.txt";
		if (_WE.verbosity>1)
		  {
			_WE.EH.println("Regression: number of predicted data=%d%n",Npredicted);
		  }
		_WE.EH.println("Printing lumped shear regression data to "+regressionFilename);
		try 
		  {
			out = new PrintWriter(new BufferedWriter(new FileWriter(regressionFilename)));
		  } 
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open " + regressionFilename);
		  }
		out.println("DT DTShear");


		// pass over regression input data...
		
		//		int mode=-1;
		// out.printf("#m="+mode+",S=4%n");
		for (i=0; i<_m.nWells; i++)
		  {
			logData = _WE.logData[i];
			info = _WE.wellInfo[i];
			if (info.wellLogDefinitions.shearTreatment != 0)
			  {
				for (j =0; j< logData.DT.length ; j++)
				  {
					if ((logData.DTshear[j]>0) && logData.DT[j]>0) 
					  {
						out.printf("%9.6g %9.6g%n",logData.DTshear[j],logData.DT[j]);
					  }
				  }
				// mode--; out.printf("#m="+mode+",S=4%n");
				out.printf("%n");
			  }
		  }
		
		
		// then regressed data
		for (i=0; i<_m.nWells; i++)
		  {
			logData = _WE.logData[i];
			for (j =0; j< logData.DT.length ; j++)
			  {
				if ((logData.DTshear[j]<0) && logData.DT[j]>0) 
				  {
					if (logData.DTshear[j] != FundamentalConstants.BAD_NUMBER)
					  {
						// one of the regressed guys. Restore + print
						logData.DTshear[j] = -logData.DTshear[j];
						out.printf("%9.6g %9.6g%n",logData.DTshear[j],logData.DT[j]);
					  }
				  }
			  }
			// mode--; out.printf("#m="+mode+",S=4%n");
			out.printf("%n");
		  }		
		out.close();
	  }
	if (_WE.verbosity>2)
	  {
		_WE.EH.println("");
	  }
  }
}
 
