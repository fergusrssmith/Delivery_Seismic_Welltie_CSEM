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

public class KeyWords 
{
  /**
	  String for use in parsers to mean measured depth
  */
  public static final String Zm = "Zm";
  /** 
	  String for use in parsers to mean the Pwave DT
  */
  public static final String DT = "DT";
  /** 
	  String for use in parsers to mean the shear wave DT
  */
  public static final String DTshear = "DTshear";
  /** 
	  String for use in parsers to mean density
  */
  public static final String rho = "rho";
  /** 
	  String for use in parsers to mean gamma
  */
  public static final String gamma = "gamma";
  /** 
	  String for use in parsers to mean "X"
  */
  public static final String X="X";
 /** 
	  String for use in parsers to mean "Y"
  */
  public static final String Y="Y";
 /** 
	  String for use in parsers to mean "Z"
  */
  public static final String Z="Z";
  /**
	  String for use in parsers to mean measured one-way time
  */
  public static final String T = "T";
  /**
	  String for use in parsers to mean measured one-way time errors
  */
  public static final String sigmaT = "sigmaT";
}
