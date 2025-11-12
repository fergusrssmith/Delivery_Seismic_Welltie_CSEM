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

import cern.colt.list.*;
import au.csiro.JamesGunning.lib.ErrorHandler;


import java.io.*;
import java.util.*;

/**
   Basic GEOEAS style ascii log parser, as example implementation
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class GEOEASWellLogReader implements WellLogReader
{
  private String ZmName = null;
  private String DTName = null;
  private String DTshearName = null;
  private String rhoName = null;
  private String gammaName = null;

  private int DTCol,DTshearCol,ZmCol,rhoCol,gammaCol;

  /** 
	  Test harness: makes up a fictitious well to stdout.
   */
  public static void main(String argv[])
  {
	int i,state=1;
	double DT,rho,TRANSITION_EPS=0.01;
	ErrorHandler EH = new ErrorHandler();
	

	EH.println("well data%n4%nMD%nDT%nDTTS%nRHOB");
	// make Log from binary markov chain
	for (i=50;i<15000;i++)
	  {
		if (Math.random()<TRANSITION_EPS) state = -state;
		DT = 80 +10*state+10.0*Math.random();
		rho = 2.2 +.15*state+.15*Math.random();
		EH.printf("%10.8g %10.8g %10.8g %10.8g%n",i, DT, DT/1.7, rho);
	  }
  }

 /**
	 Set column -header words to appropriate meanings
  */
  public void setLogDataKeyword(String keyword,String  word)
  {
	if (keyword.equals(KeyWords.Zm)) ZmName = word;
	if (keyword.equals(KeyWords.DT)) DTName = word;
	if (keyword.equals(KeyWords.DTshear)) DTshearName = word;
	if (keyword.equals(KeyWords.rho)) rhoName = word;
	if (keyword.equals(KeyWords.gamma)) gammaName = word;
  }

  /**
	 Open the named file and stuff all the required data (assuming the header words have been correctly
	 initialised using the setLogDataKeyword_*() methods into the correct pieces of logData
  */
  public boolean extractData(String filename, LogData logData, ErrorHandler EH)
  {
	// data buckets
	DoubleArrayList Zm = new DoubleArrayList();
	DoubleArrayList DT = new DoubleArrayList();
	DoubleArrayList DTshear = new DoubleArrayList();
	DoubleArrayList rho = new DoubleArrayList();
	DoubleArrayList gamma = new DoubleArrayList();

	// get a general geoEAS reader
	ASCIIfileReader gReader = new ASCIIfileReader(EH,true,true);
 
	// form list of required things
	if (ZmName != null ) gReader.addToContainerList(ZmName,Zm);
	if (DTName != null ) gReader.addToContainerList(DTName,DT);
	if (DTshearName != null ) gReader.addToContainerList(DTshearName,DTshear);
	if (rhoName != null ) gReader.addToContainerList(rhoName,rho);
	if (gammaName != null ) gReader.addToContainerList(gammaName,gamma);

	// now get the goodies
	gReader.extractData(filename);
	  
	// Populate relevent fields of logdata class
	logData.Zm = Zm.elements();
	logData.DT = DT.elements();
 	logData.DTshear = DTshear.elements();
	logData.rho = rho.elements();
	if (gamma != null && (gamma.size()!=0))
	  { 
		logData.gamma = gamma.elements();
	  }
	else logData.gamma = null;

	// check units

	if (gReader.parsingInfo.isLAS2)
	  {
		if (DTName!=null && DTshearName!=null)
		  {
			String DTunit,DTSunit;
			DTunit = gReader.parsingInfo.unitsMap.get(DTName);
			DTSunit = gReader.parsingInfo.unitsMap.get(DTshearName);
			if (DTunit != null && DTSunit != null)
			  {
				if (!DTunit.equals(DTSunit)) 
				  {
					EH.println("\t\tMismatching sonic log units ("+DTunit+","+DTSunit+") in file "+filename);
					EH.println("\t\tCorrecting mismatching units:");

					int i;
					double factor=1.0;
					if (DTunit.contains("m") && DTSunit.contains("ft"))
					  {
						factor = 1.0/FundamentalConstants.METRES_PER_FOOT;
						EH.println("\t\tShear log: */ft -> */m");
					  }
					else  if (DTunit.contains("ft") && DTSunit.contains("m"))
					  {
						factor = FundamentalConstants.METRES_PER_FOOT;
						EH.println("\t\tShear log: */m -> */ft");
					  }
					else
					  {
						EH.println("\t\tUnable to correct mismatching sonic log units. Correct log file manually!");
						return false;
					  }
					double[] DTS = logData.DTshear;
					for (i=0;i<DTS.length;i++)
					  {
						if (DTS[i]>0) DTS[i] *= factor;
					  }
				  }
			  }
		  }

		// check of sonic units for auto checkshots
		if (gReader.parsingInfo.isLAS2)
		  {
			String DTunit,Zmunit;

			Zmunit = gReader.parsingInfo.unitsMap.get(ZmName);
			DTunit = gReader.parsingInfo.unitsMap.get(DTName);
			if (!WellLogDefinitions.MD_allowableunits.contains(Zmunit))
			  {
				EH.println("WARNING: code does not recognize MD unit "+Zmunit+" of log file "+filename);
			  }
			if (!WellLogDefinitions.DT_allowableunits.contains(DTunit))
			  {
				EH.println("WARNING: code does not recognize DT unit "+DTunit+" of log file "+filename);
			  }
			if (Zmunit != null)
			  {
				// from definition of WellLogDefinitions.MD_allowableunits we contrive:
				if (Zmunit.startsWith("m")) logData.Zm_unit = WellLogDefinitions.MD_unit.M;
				if (Zmunit.startsWith("f")) logData.Zm_unit = WellLogDefinitions.MD_unit.FT;
			  }
			if (DTunit != null)
			  {
				// from definition of WellLogDefinitions.DT_allowable units we contrive:
				if (DTunit.contains("us/m")) 
				  {
					logData.DT_unit = WellLogDefinitions.DT_unit.US_PER_M;
				  }
				else if (DTunit.contains("s/m")) 
				  {
					logData.DT_unit = WellLogDefinitions.DT_unit.S_PER_M;
				  }
				else if (DTunit.contains("us/ft")) 
				  {
					logData.DT_unit = WellLogDefinitions.DT_unit.US_PER_FT;
				  }
			  }
		  }

		// density
		if (rhoName!=null)
		  {
			logData.RHO_unit = gReader.parsingInfo.unitsMap.get(rhoName);
		  }
	  }
	return true;
  }
}
