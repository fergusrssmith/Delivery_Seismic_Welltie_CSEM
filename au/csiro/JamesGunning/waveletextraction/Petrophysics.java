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
import au.csiro.JamesGunning.lib.*;

/**
   Basic petrophysics entries, constants etc
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/11/2006
*/

public class Petrophysics
{
  boolean hasAnisotropy = false;  // set for each model

  boolean petrophysicsProvided = false;

  boolean GRThreshProvided = false;
  boolean densityThreshProvided = false;

  float GR_intercept=70f;
  float GR_slope=0f;
  float density_intercept=2.3f;
  float density_slope = 0f;

  int version=1;

  /**
	 Default constructor
  */
  Petrophysics()  {}

  /**
	 Constructor for ExtractionViewer
  */
  Petrophysics(String encodingString, ErrorHandler EH)
  {
	if (encodingString==null) return;

	String[] fields = encodingString.split("[=,:]");
	if (fields!=null && fields.length>1)
	  {
		if (!fields[0].equals("version"))
		  {
			EH.warn("Petrophysics encodingString cannot be deciphered for version: ignoring petrophysics info!");
			return;
		  }
		try
		  {
		    version= Integer.parseInt(fields[1]);
		  }
		catch (NumberFormatException e)
		  {
			EH.die_err("Petrophysics encodingString "+encodingString+" has non-integral version number! ");
		  }
	  }
	else
	  {
		EH.warn("Petrophysics encodingString cannot be parsed. Ignoring petrophysics info! ");
		return;
	  }

	// Get to here only if encodingString looks sensible enough for version-specifics
	// set ups for various versions
	if (version==1)
	  {
		// Typically......"version=1:true,true,70.0,0.0,true,2.2,0.0"
					   
		if (fields.length != 9)
		  {
			EH.warn("Version 1 Petrophysics encodingString "+encodingString+" cannot be parsed. Ignoring petrophysics info! ");
			return;
		  }
		else
		  {
			try
			  {
				hasAnisotropy  = Boolean.parseBoolean(fields[2]);
				GRThreshProvided = Boolean.parseBoolean(fields[3]);
				GR_intercept = Float.parseFloat(fields[4]);
				GR_slope = Float.parseFloat(fields[5]);

				densityThreshProvided = Boolean.parseBoolean(fields[6]);
				density_intercept = Float.parseFloat(fields[7]);
				density_slope = Float.parseFloat(fields[8]);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Petrophysics encodingString "+encodingString+" does not have format "+
						   "version=1:boolean,boolean,float,float,boolean,float,float");
			  }
		  }
	  }
	// else
	// other versions....


	petrophysicsProvided = true;
  }

  void setVersion(int v)
  {
	version = v;
  }

  boolean isShale(boolean haveGamma, double... args)
  {
	if (version==1)
	  {
		// varargs for version 1 are: zTVD, gamma, rho
		double zTVD, gamma, rho;
		zTVD = args[0];
		gamma = args[1];
		rho = args[2];
		
		if (haveGamma && GRThreshProvided)
		  {			
			double gammaShaleThreshold = GR_intercept + GR_slope * zTVD;
			if (densityThreshProvided)
			  {
				// joint criterion
				double densityShaleThreshold = density_intercept + density_slope * zTVD;
				return ((gamma > gammaShaleThreshold) && (rho > densityShaleThreshold));
			  }
			else
			  {
				// GR threshold only
				return (gamma > gammaShaleThreshold);
			  }
		  }
		else
		  {
			double densityShaleThreshold = density_intercept + density_slope * zTVD;
			return (rho > densityShaleThreshold);
		  }
	  }
	// else 
	// later versions.....fill in..

	// never gets here!!!
	return false;
  }
  
  String encodingString()
  {
	StringBuffer code = new StringBuffer(" Petrophysics:[version="+version+":");
	if (version==1)
	  {
		code.append(hasAnisotropy+",");
		code.append(GRThreshProvided+","+GR_intercept+","+GR_slope+",");
		code.append(densityThreshProvided+","+density_intercept+","+density_slope);
	  }
	// else .... 
	// other versions...




	code.append("]");

	return code.toString();
  }
}
