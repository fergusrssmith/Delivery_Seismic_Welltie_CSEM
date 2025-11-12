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
import static au.csiro.JamesGunning.seismicinverter4D.FluidStates.FluidTypes.*;

/**
   Basic class for data storage of the concept of fluid units and
   contacts within a fluid unit.

   A fluid unit is a contiguous sequence of layers within which fluids
   are density ordered. Each fluid unit can have a only one each of
   the possible set of {OB, GB, GO...} fluid contacts. Users specify
   priors for the known contacts in terms of depth.
   
   @author James.Gunning@csiro.au
   @version 1.0, 13 Sept/2006
*/

public final class FluidUnitsAndContacts
{
  Model _m;
  FluidUnit[] fluidUnits;
  boolean[][] fluidRequired;
  ErrorHandler _EH;
  boolean contactsVariable=false;

  // working utility array in sensitivity calculations
  ArrayList<LayerProperties.Property> fc_props = null;

  public FluidUnitsAndContacts(Model m)
  {
	this._m = m;
	_EH = m._SI.EH;
  }

 /**
	Builds boolean matrix of whether layer i
	requires fluid type j (j=0,1,2,3 for b,o,l,g)
 */
  void buildRequiredFluidsTable()
  {
	int i,j,v,k,i_h;
	FluidUnit fu;
	ContactInfo CI;

	if (!_m._SI.doSimpleFluidUnitsAndContacts) return;

	fluidRequired = new boolean[_m.n_layers][FluidStates.FluidTypes.Nfluids];
	
	for (i=0;i<fluidUnits.length; i++)
	  {
		fu = fluidUnits[i];
		for (j=fu.firstLayerIndx-1;
			 j <= fu.finalLayerIndx-1; j++)
		  {
			// j can be expected to span all layers
			// in this fluid unit
			// Mark brine as required by default
			fluidRequired[i][FluidStates.FluidTypes.brine] = true;
			// examine all contacts
			for(v=0;v<_m.Nvintages;v++)
			  {
				for (k=0;k<fu.contactInfo[v].length; k++)
				  {
					CI = fu.contactInfo[v][k];
					fluidRequired[j][CI.UF] = true;
					fluidRequired[j][CI.LF] = true;
				  }
			  }
		  }
		if (fu.fluids_synched)
		  {
			fu.masterFluidLayer = new int[FluidStates.FluidTypes.Nfluids];
			// mark master layer of each fluid type for synchronisation purposes
			for (i_h=0; i_h<FluidStates.FluidTypes.Nfluids; i_h++) // brine, oil, lsg, gas
			  {	
				// default marker
				fu.masterFluidLayer[i_h] = -1;
			  }
		  }
	  }
  }

  void updateMasterFluidLayerTable()
  {
	if (!_m._SI.doSimpleFluidUnitsAndContacts) return;

	int i,j,i_h;
	ModelPriorDescriptions.MetaInfo MI;

	// loop on layers
	for (i=0;i<fluidUnits.length; i++)
	  {
		FluidUnitsAndContacts.FluidUnit fu = fluidUnits[i];
		if (fu.fluids_synched)
		  {
			// mark master layer of each fluid type for synchronisation purposes
			for (i_h=0; i_h<FluidStates.FluidTypes.Nfluids; i_h++) // brine, oil, lsg, gas
			  {	
				for (j=fu.firstLayerIndx-1;
					 j <= fu.finalLayerIndx-1; j++)
				  {
					// j can be expected to span all layers in FU
					MI = _m.MPD.MetaInfo[j];
					if (MI.is_permeable && MI.isActive && fluidRequired[j][i_h])
					  {						
						// This is first layer in this FU that really controls
						// fluid properties under synchronisation. Mark and break
						fu.masterFluidLayer[i_h] = j;
						break;
					  }
				  }
			  }
		  }
	  }
  }

  /**
	 returns the fluid unit responsible for layer i_layer (zero offset)
  */
  FluidUnit getFUofThisLayer(int i_layer)
  {
	FluidUnit fu = fluidUnits[0];

	for (int i=0;i<fluidUnits.length; i++)
	  {
		fu = fluidUnits[i];
		if ((i_layer >= (fu.firstLayerIndx-1)) && (i_layer <= (fu.finalLayerIndx-1))) return fu;
	  }
	// never gets here
	return fu;
  }

  /**
	 return true if fluid unit owning (zero offset) layer i_layer has synching applied
  */
  boolean fluidsSynchedForFUowningLayer(int i_layer)
  {
	FluidUnit fu = getFUofThisLayer(i_layer);
	return fu.fluids_synched;
  }

  void updateRequiredFluidsTable(int i, boolean is_permeable)
  {
	if (!_m._SI.doSimpleFluidUnitsAndContacts) return;

	if (!is_permeable)
	  {
		// layer index i is impermeable, => requires nothing
		for (int j= FluidStates.FluidTypes.oil; j < FluidStates.FluidTypes.Nfluids; j++) 
		  {
			fluidRequired[i][j]=false;
		  }
	  }
  }
  /**
	Performs two possibly fatal checks  here:
	// (1) within a fluid unit, all instances of e.g. "oil" must be the same actual fluid.
	// (2) interface between fluid units must have at least one impermeable layer to make sense
	*/
  void checkFluidConsistencyWithinUnits()
  {
	int i,j;
	String g,l,o,b;
	FluidUnit fu;
	ModelPriorDescriptions.MetaInfo MI;

	if (!_m._SI.doSimpleFluidUnitsAndContacts) return;

	for (i=0;i<fluidUnits.length; i++)
	  {
		fu = fluidUnits[i];
		g=null;
		l=null;
		o=null;
		b=null;
		for (j=fu.firstLayerIndx-1; j <= fu.finalLayerIndx-1; j++)
		  {
			MI = _m.MPD.MetaInfo[j];
			if (MI.gas_name!=null)
			  {
				if (g==null) g =  MI.gas_name;
				if (!MI.gas_name.equals(g))
				  {
					_m._SI.EH.die_err("Gas names ["+g+","+MI.gas_name + 
									  "] within fluid unit spanning layers "+fu.firstLayerIndx+".."
									  +fu.finalLayerIndx+" are not identical");
				  }
			  }
			if (MI.lsg_name!=null)
			  {
				if (l==null) l =  MI.lsg_name;
				if (!MI.lsg_name.equals(l))
				  {
					_m._SI.EH.die_err("Lsg names ["+l+","+MI.lsg_name + 
									  "] within fluid unit spanning layers "+fu.firstLayerIndx+".."
									  +fu.finalLayerIndx+" are not identical");
				  }
			  }
			if (MI.oil_name!=null)
			  {
				if (o==null) o =  MI.oil_name;
				if (!MI.oil_name.equals(o))
				  {
					_m._SI.EH.die_err("Oil names ["+o+","+MI.oil_name + 
									  "] within fluid unit spanning layers "+fu.firstLayerIndx+".."
									  +fu.finalLayerIndx+" are not identical");
				  }
			  }
			if (MI.brine_name!=null)
			  {
				if (b==null) b =  MI.brine_name;
				if (!MI.brine_name.equals(b))
				  {
					_m._SI.EH.die_err("Brine names ["+b+","+MI.brine_name + 
									  "] within fluid unit spanning layers "+fu.firstLayerIndx+".."
									  +fu.finalLayerIndx+" are not identical");
				  }
			  }
		  }
	  }
	for (i=0;i<fluidUnits.length-1; i++)
	  {
		// at least one of final layer of fluid unit above or first layer of
		// fluid unit below must be impermeable
		if (_m.MPD.MetaInfo[fluidUnits[i].finalLayerIndx-1].is_permeable
			&& _m.MPD.MetaInfo[fluidUnits[i+1].firstLayerIndx-1].is_permeable)
		  {
			_m._SI.EH.die_err("Fluid units of layer spans ["+fluidUnits[i].firstLayerIndx+".."
							  +fluidUnits[i].finalLayerIndx+"] and ["+fluidUnits[i+1].firstLayerIndx+".."
							  +fluidUnits[i+1].finalLayerIndx+"] are not bounded by an "
							  +" impermeable layer in at least one fluid unit. Contacts don't make sense");
		  }
	  }
  }

  boolean isFluidContactMeanSequenceSane()
  {
	if (!_m._SI.doSimpleFluidUnitsAndContacts) return true;

	int i,j,k;
	double z,zlast=0;
	for (i=0;i<	_m.fluidUnitsAndContacts.fluidUnits.length; i++)
	  {
		FluidUnitsAndContacts.FluidUnit fu = _m.fluidUnitsAndContacts.fluidUnits[i];
		for(j=0;j<_m.Nvintages;j++)
		  {
			for (k=0; k< fu.contactInfo[j].length; k++)
			  {
				FluidUnitsAndContacts.ContactInfo CI = 	fu.contactInfo[j][k];
				z = CI.depthPD.value;  // prior mean of this contact
				if (k>0 && z<zlast)
				  {
					_m._SI.EH.println("Strange sequence of mean contact depths for specified "
									  +"contacts: vintage "+(j+1)
									  + " fluid unit "+(i+1)+" .Depths: "+z+", "+zlast);
					return false;
				  }
				zlast = z;
			  }			  
		  }
	  }  
	return true;
  }


 /**
	 Interrogate contact info to see if fluid f
	 is needed in layer indexed i_layer (zero offset sense)
  */
  boolean contactsRequireThisFluid(int i_layer,short f)
  {
	int i,i_test,k;
	FluidUnit fu;
	ContactInfo CI;
	boolean need=false;
	
	if (!_m._SI.doSimpleFluidUnitsAndContacts) return false;

	return fluidRequired[i_layer][f];
  }

  boolean checkFUConfigurationIsPossible()
  {
	int i,k;
	boolean checkPassed=true;
	ContactInfo CI;
	FluidUnit fu,fu2;
	short fl;

	// (1) layer-span check - should be non-overlapping and a complete covering after sort

	if (fluidUnits[0].firstLayerIndx != 1)
	  {
		_EH.println("Uppermost fluid_unit must have top_layer=1");
		_EH.println("(Fluid units must span entire model layer range)");
		_EH.die_err("XML fluid contact entries not sensible");
	  }
	if (fluidUnits[fluidUnits.length-1].finalLayerIndx != _m.n_layers)
	  {
		_EH.println("Lowermost fluid_unit must have bottom_layer="+_m.n_layers);
		_EH.println("(Fluid units must span entire model layer range)");
		_EH.die_err("XML fluid contact entries not sensible");
	  }

	layerchck:
	for (i=0;i<fluidUnits.length-1; i++)
	  {
		fu = fluidUnits[i];
		if (fu.finalLayerIndx >= fluidUnits[i+1].firstLayerIndx)
		  {
			_EH.println("Overlapping fluid_unit layer spans (top_layer..bottom_layer) of in fluid units "+(i+1)+" and "+(i+2));
			checkPassed=false;
			break layerchck;
		  }
		if (fu.finalLayerIndx != (fluidUnits[i+1].firstLayerIndx - 1))
		  {
			_EH.println("Fluid_units have missing 'orphaned' layers (owned by no fluid unit) between "+
						" layers "+fu.finalLayerIndx + " and "+(fluidUnits[i+1].firstLayerIndx));
			checkPassed=false;
			break layerchck;
		  }
	  }
	if (!checkPassed) _EH.die_err("XML fluid contact entries not sensible");

	return true;
  }


  public static class FluidUnit implements Comparable
  {
	int firstLayerIndx;
	int finalLayerIndx;
	// this is set to indicate that stochastic fluid properties within this fluid unit
	// should be synchronised
   	boolean fluids_synched = false;
	// masterFluidLayer[i_h] is master fluid layer for fluid type i_h in this fluid unit, if required
	int[] masterFluidLayer =null;
	/**
	   array of contacts for this fluid unit
	*/
	FluidUnitsAndContacts.ContactInfo[][] contactInfo;
  
	/**
	   Comparator to sort by first then second layer index
	*/
	public int compareTo(Object o)
	{
	  FluidUnit n = (FluidUnit)o;
	  int cmp = (firstLayerIndx<n.firstLayerIndx)? -1 : ((firstLayerIndx>n.firstLayerIndx)? 1 : 0 );
	  if (cmp==0)
		{
		  cmp = (finalLayerIndx<n.finalLayerIndx)? -1 : ((finalLayerIndx>n.finalLayerIndx)? 1 : 0 );
		}
	  return cmp;
	}
  }

  /**
	 basic info of a possible contact within a fluid unit
  */
  public static class ContactInfo
  {
	short UF;              // upper fluid type
	short LF;              // lower fluid type

	ModelPriorDescriptions.PriorDetails depthPD;
	ModelPriorDescriptions.PriorDetails sigma_depthPD;

	LayerProperties.Property depth;         // depth of this contact
	LayerProperties.Property sigma_depth;   // sigma_depth of this contact

	/**
	   return true if supplied key (eg "GO") matches this contact fluid type
	*/
	boolean fcTypeMatches(String key)
	{
	  char U,L;
	  if (key==null || key.length() != 2) return false;
	  U = key.charAt(0);
	  L = key.charAt(1);

	  if (UF != getFluidShort(U)) return false;
	  if (LF != getFluidShort(L)) return false;
	  
	  // must match
	  return true;
	}

	/**
	   returns error message if something fails, null otherwise
	*/
	String setContactTypeFromXMLKey(String key)
	{
	  short i,c=0;
	  char U,L;

	  if (key.length() != 2) 
		{
		  return new String("key has more than 2 characters");
		}

	  U = key.charAt(0);
	  L = key.charAt(1);

	  UF = getFluidShort(U);
	  LF = getFluidShort(L);

	  if (UF==-1 || LF==-1) return new String("Contact fluid key "+key+" is illegal");
	  if (UF==LF) return new String("Contact fluid key "+key+" has matching fluids"); 
	  // if (UF<LF) 
	  if (FluidStates.FluidTypes.density_priority[UF] < FluidStates.FluidTypes.density_priority[LF])
		{
		  // return new String("Contact fluid key "+key+" is density-inverted: must be an left-right ordered pair from {G,L,O,B} ");
		  return new String("Contact fluid key "+key+" is density-inverted: must be an left-right ordered pair from {G,O,L,B} ");
		}
	  return null;
	}  
  }
}


