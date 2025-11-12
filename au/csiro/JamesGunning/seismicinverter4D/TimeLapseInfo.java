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
package au.csiro.JamesGunning.seismicinverter4D;

import java.util.*;

import au.csiro.JamesGunning.lib.ErrorHandler;

/**
   Basic class of extra Model class stuff needed
   under time lapse
   
   @author James.Gunning@csiro.au
   @version 1.0, 15Apr/2011
*/
public final class TimeLapseInfo
{
  ModelPriorDescriptions.PriorDetails[] DeltaP4Dinfo=null;
  ModelPriorDescriptions.PriorDetails[] SigmaDeltaP4Dinfo=null;
  ModelPriorDescriptions.PriorDetails[] SigmaDeltaTinfo=null;
  LayerProperties.Property[] DeltaP4D = null; // pressure change properties, 4D case
  LayerProperties.Property[] DeltaT = null; // time registration shift properties, 4D case

  // for use in the design matrix
  List<LayerProperties.Property> DeltaProp4Dlist=null;

  // Experimental: for use in ordered samplers (not used at present - implemented another way)  
  List<LayerProperties.Property[]> monotonicDeltaProp4Dlist=null;

  public TimeLapseInfo(int N)
  {
	DeltaP4Dinfo      = new ModelPriorDescriptions.PriorDetails[N];
	SigmaDeltaP4Dinfo = new ModelPriorDescriptions.PriorDetails[N];
	SigmaDeltaTinfo = new ModelPriorDescriptions.PriorDetails[N];
	DeltaP4D =  new LayerProperties.Property[N];
	DeltaT =  new LayerProperties.Property[N];
  }
}
