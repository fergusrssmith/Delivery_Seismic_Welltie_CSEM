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
/**
   Data class for storing info for each well
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class GeneralWellInfo
{
  
  /**
	 Name of well data group
  */
  String well_name;
  /**
	 Filename of log data
  */
  String logFilename;
  /**
	 Filename of survey data
  */
  String surveyFilename;
  /**
	  Filename of seismic min-cube with adequate traces for this well
  */
  String seismicFilename;
  /**
	 Filename of checkshots
  */
  String checkshotFilename;
  /**
	 automake checkshots
  */
  boolean checkshotAutoMake = false;
  /**
	 model-select knots in automake checkshots
  */
  boolean checkshotAutoMakeModelSelect = false;

  boolean needPSreflectivities=false;
  /**
	 Filename of markers;
  */
  String markersFilename=null;

  /**
	 Name of well which we might synchronize to...
  */
  String synchronizeWellName;
  /**
	 Index of well which we are synchronized to, if synchronisation operates
  */
  int synchronizeWellIndex = 0;

  /**
	 Geometry info for the well data
  */
  Geometry geom;
  /**
	 Well log definitions
  */
  WellLogDefinitions wellLogDefinitions;
  /**
	 Checkshot and marker definitions
  */
  CheckshotAndMarkerDefinitions chkshotAndMarkerDefinitions;

  /**
	 Global uncertainty parameters on well position and systematic time error
  */
  double systematicSigmaT = 0.0;
  double systematicSigmaX = 0.0;
  double systematicSigmaY = 0.0;

  
  /**
	 Property for the systematic X uncertainty
  */
  Property systematicXProp;

  /**
	 Property for the systematic X uncertainty
  */
  Property systematicYProp;

  /**
	 Collected info for any 4D / saturation modelling
  */
  SaturationInfo saturationInfo = null;

  public GeneralWellInfo()
  {
	geom = new Geometry();
	wellLogDefinitions = new WellLogDefinitions();
	chkshotAndMarkerDefinitions = new CheckshotAndMarkerDefinitions();
	systematicXProp = new Property(0.0, false, false);
	systematicYProp = new Property(0.0, false, false);
  }

  // set up saturation class
  void setupSaturationInfo()
  {
	saturationInfo = new SaturationInfo();
  }  
}
