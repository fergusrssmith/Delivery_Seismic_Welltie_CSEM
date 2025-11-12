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
Basic GEOEAS style ascii well survey reader, as example implementation

@author James.Gunning@csiro.au
@version 1.0, 12/2/2003
*/

public class GEOEASWellSurveyReader implements WellSurveyReader
{
private String XName = "X";    // Easting?
private String YName = "Y";  // northing?
private String ZName = "Z";  // true vertical depth (with kelly bushing taken into account?)
private String ZmName = "Zm"; // measured depth (as per the log)

/**
Set column -header words to appropriate meanings
*/
public void setSurveyDataKeyword(String keyword,String  word)
  {
	if (keyword.equals(KeyWords.Zm)) ZmName = word;
	if (keyword.equals(KeyWords.Z)) ZName = word;
	if (keyword.equals(KeyWords.X)) XName = word;
	if (keyword.equals(KeyWords.Y)) YName = word;
  }

  /**
	 Open the named file and stuff all the required data (assuming the header words have been correctly
	 initialised using the setLogDataKeyword_*() methods into the correct pieces of logData
  */
  public boolean extractSurveyData(String filename, LogData logData, ErrorHandler EH)
  {
	// data buckets
	DoubleArrayList X = new DoubleArrayList();
	DoubleArrayList Y = new DoubleArrayList();
	DoubleArrayList Z = new DoubleArrayList();
	DoubleArrayList Zm = new DoubleArrayList();

	// get a general geoEAS reader
	ASCIIfileReader gReader = new ASCIIfileReader(EH,true,false);

	// form list of required things
	gReader.addToContainerList(XName,X);
	gReader.addToContainerList(YName,Y);
	gReader.addToContainerList(ZName,Z);
	gReader.addToContainerList(ZmName,Zm);

	// now get the goodies
	gReader.extractData(filename);
	  
	// Populate relevent fields of logdata class
	logData.survey.X = X.elements();
	logData.survey.Y = Y.elements();
 	logData.survey.Z = Z.elements();
	logData.survey.Zm = Zm.elements();
	logData.survey.N = logData.survey.X.length;

	return true;
  }

}
