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



package au.csiro.JamesGunning.seismicinverter;

/**
   To get co-located seismic traces from the layer model, we need
   a seismic model that does the job of simulating seismic response from 
   a medium with a source. The simplest model is convolution, but other
   flavours are of interest too. So let's make the SeismicModel an
   interface, with Convolution one possible implementation
*/

public interface SeismicModel
{
  /**
	 Any seismic model must have parameters that 
	 can be defined via an XML input.
	 These can be parsed via a suitable parser.
	 Different models may have wildy different parameters, so they need
	 a unique parser
  */
  public ModelPriorDescriptions.PriorDetails getInternalStochasticPriorVariable(int index);
  public ModelPriorDescriptions.PriorDetails[] getInternalStochasticPriorVariablesArray();
  public boolean initialiseSeismicModel(int n_layers);
  public void generateSyntheticTrace(); 
  public void dumpSynthetics(boolean newState, int realisation_number);
  public double computeSyntheticAndSeismicError();
  /**
	 Returns error vector of model mismatches at all relevant data points,
	 scaled by their respective error standard deviations
  */
  public double[] getErrorVector();

  public void setNearstackData(StackData sd);
  public StackData getNearstackData();
  public void setFarstackData(StackData sd);
  public StackData getFarstackData();
}


