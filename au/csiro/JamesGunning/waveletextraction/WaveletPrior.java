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
import java.io.*;

import cern.colt.list.*;
import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;




/**
   Class to read in a wavelet prior from a previous wavelet extraction
   
   @author James.Gunning@csiro.au
   @version 1.0, 6/12/2006
*/

public class WaveletPrior
{
  /**
	Calling Extractor
  */
  private WaveletExtractor _WE;

  /**
	 Calling Extractor's error handler
  */
  private ErrorHandler _EH;

  double dtCoarse_prior=-1;
  double f1_prior=0;

  DoubleArrayList amp_means;
  DoubleArrayList amp_stddev;

  boolean test=false;

  public WaveletPrior(WaveletExtractor WE)
  {
	this._WE = WE;
	this._EH = _WE.EH;

	parseMLparameterFile(_WE.wavelet_prior_filename);
  }

  // test version with main() only
  public WaveletPrior(ErrorHandler EH, String filename)
  {
	this._EH = EH;

	test=true;
	parseMLparameterFile(filename);
  }

  /** 
	  Sample driver to show first few lines of a named file after parsing
  */
  public static void main(String argv[])
  {
	ErrorHandler EH = new ErrorHandler();
	String filename = null;

	OptionsParser options = new OptionsParser("WaveletPrior",EH);

   	options.addMandatoryArg(1,"filename","MLparameters file to parse");	

	options.parseOptions(argv,true);

	filename = options.getMandatoryArg(1);
	
	new WaveletPrior(EH,filename);
  }

  void parseMLparameterFile(String filename)
  {

	String s,word;
	StringTokenizer st;
	double val=0;
	BufferedReader f;
	int n_coeffs=0;

	
	try 
	  {
        f = new BufferedReader(new FileReader(filename));			
		int lineNo=0;

		amp_means = new DoubleArrayList();
		amp_stddev = new DoubleArrayList();

		_EH.println("Parsing wavelet prior from %s",filename);
		while ((s=f.readLine()) != null)
		  {
			// strip leading & trailing white space
			s = s.trim();
			lineNo++;

			// lineNo==1 is title
			if (lineNo==1)
			  {
				// expect format:
				// "List of parameters and deviations. ( Wavelet coarse dt= 0.00252776 f1= -0.0151666 )"
				st = new StringTokenizer(s);
				while (st.hasMoreTokens())
				  {
					word = st.nextToken();
					if (word.equals("dt="))
					  {
						if (st.hasMoreTokens())
						  {
							word = st.nextToken();
							try
							  {
								dtCoarse_prior = Double.parseDouble(word);
							  }
							catch (NumberFormatException e)
							  {
								_EH.die_err("dt not a number in first line of "+filename);
							  }
							if (dtCoarse_prior < 0)
							  {
								_EH.die_err("dt negative in first line of "+filename);
							  }
						  }
					  }
					if (word.equals("f1="))
					  {
						if (st.hasMoreTokens())
						  {
							word = st.nextToken();
							try
							  {
								f1_prior = Double.parseDouble(word);
							  }
							catch (NumberFormatException e)
							  {
								_EH.die_err("f1 not a number in first line of "+filename);
							  }
						  }
					  }
				  }
				// now finished line 1. Do checks
				if (dtCoarse_prior==-1)
				  {
					// probably missing: old file
					_EH.warn("dt and f1 missing from first line of "+filename+". Use a recent extraction which writes these entries!");
					_EH.die_err("Example first line should be:\nList of parameters and deviations. ( Wavelet coarse dt= 0.00252776 f1= -0.0151");
				  }
			  }
			else
			  {
				// Expect format (e.g.)
				// wavelet_coeff_3   0.00000    66764.6      38654.9                       0.102990
				st = new StringTokenizer(s);
				word = st.nextToken();
				if (word.startsWith("wavelet_coeff"))
				  {
					if (st.countTokens() < 4)
					  {
						_EH.die_err("wavelet_coeff lines in "+filename+" have insufficient entries");
					  }
					word = st.nextToken(); // prior mean
					word = st.nextToken(); // ML value
					try
					  {
						val = Double.parseDouble(word);
					  }
					catch (NumberFormatException e)
					  {
						_EH.die_err("wavelet_coeff MLvalue not a number in line "+lineNo+" of "+filename);
					  }
					amp_means.add(val);

					word = st.nextToken(); // stddev
					try
					  {
						val = Double.parseDouble(word);
					  }
					catch (NumberFormatException e)
					  {
						_EH.die_err("wavelet_coeff stddev not a number in line "+lineNo+" of "+filename);
					  }
					amp_stddev.add(val);
					
					n_coeffs++;
				  }
				// finish line
				while (st.hasMoreTokens()) {st.nextToken();}				
			  }
		  }
		f.close();
		
		amp_means.trimToSize();
		amp_stddev.trimToSize();

		if (amp_means.size()<1)
		  {
			_EH.die_err(filename+" does not appear to contain useful wavelet_coeff lines. "
						+"Check it - it should be a MLparameters type file.");
		  }
		if (test)
		  {
			_EH.println("Wavelet prior: dt = "+dtCoarse_prior+" f1="+f1_prior);
			_EH.println("Coefficients:");
			for (int i=0; i<amp_means.size(); i++)
			  {
				_EH.printf("%10.5g %10.5g\n",amp_means.get(i),amp_stddev.get(i));
			  }
			System.exit(0);
		  }	  
	  }
	catch (Exception e) 
	  {
		_EH.die_err("Can't open file "+filename);
	  }
  }
}
