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
   Basic magic words used for parsing purposes, independent of the
   file format of well logs
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class WellLogDefinitions
{
  /**
	 word for DT data (P-wave DT)
  */
  String DT="DT";
	
  /**
	 word for DTshear data (S-wave DT)
  */
  String DTshear="DTshear";
  /**
	 treatment of shear log:
	 0 = ignore
	 1 = use
  */
  char shearTreatment=1;  // defaults to use until contradicted in parsing
  /**
	 word for rho data
  */
  String rho="rho";
  /**
	 word for measured depth
  */
  String MD="MD";
  /** 
	  word for (opt) gamma
  */
  String gamma="gamma";
  
  /** 
	  Optional MD mute thresholds: data will be set to BAD values outside [MD_uppermute,MD_lowermute]
  */
  double MD_uppermute=Double.NaN;
  double MD_lowermute=Double.NaN;
  /** 
	  Array of acceptable units (lower case) for measured depth
  */
  public static final String MD_allowableunits = "m,meter,meters,ft,feet";
  public enum MD_unit {UNKNOWN,FT,M};

  /** 
	  Array of acceptable units (lower case) for DT
  */
  public static final String DT_allowableunits = "us/ft,us/m,s/m"; // need unicode mu here...
  public enum DT_unit {UNKNOWN,US_PER_FT,US_PER_M,S_PER_M};
}
 
