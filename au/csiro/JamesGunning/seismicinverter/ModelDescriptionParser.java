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

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.XmlUtils;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;



import cern.colt.map.*;
import cern.colt.list.*;

import javax.xml.parsers.*;
import javax.xml.xpath.*;
import org.xml.sax.*;
import org.w3c.dom.*;
import java.io.*;
import java.lang.reflect.*;
import java.util.*;


public final class ModelDescriptionParser
{
  /** 
	  Error Handler for instance of this class
  */
  ErrorHandler EH;

  Document document;

  String xml_filename=null;

  int verbosity = 0;

  /**
	 String to form link tags with, for variable prior data
  */
  public static final String LINK_SUFFIX="_link";

  /**
	 true String int attributes
  */
  public static final String TRUE_TAG="true";

  /**
	 delimits block describing model
  */
  public static final String MODEL_DESCRIPTION_TAG="model_description";

 /**
	 delimits top layer element
 */
  public static final String TOP_LAYER_TAG="top_layer";
  /**
	 delimits middle layer element
  */
  public static final String MIDDLE_LAYER_TAG="middle_layer";
  /**
	 delimits layer element
  */
  public static final String BOTTOM_LAYER_TAG="bottom_layer";
  /**
	 delimits sand_member name (permeable rock)
  */
  public static final String SAND_MEMBER_TAG="reservoir_endmember";
  
  /**
	 delimits shale_member name
  */
  public static final String SHALE_MEMBER_TAG="nonreservoir_endmember";

  /**
	 delimits shale2 name
  */
  public static final String HARD_SHALE_MEMBER_TAG="nonreservoir_hard_endmember";
  /**
	 is a layer synchronised with another
	 any stochastic rock properties are identical
  **/
  public static final String SYNCHRONISE_TAG="synchronise_rock_matrix";

  /**
	 flags presence of variable quantity in prior-input traces
  */
  public static final String VARIES_TAG="varies_areally";

  /**
	 delimits fluid names
  */
  public static final String NAME_TAG="name";

  /**
	 delimits brine name
  */
  public static final String BRINE_TAG="brine";

  /**
	 delimits oil name
  */
  public static final String OIL_TAG="oil";

  /**
	 delimits gas name
  */
  public static final String GAS_TAG="gas";

  /**
	 delimits low-sat-gas name
  */ 
  public static final String LSG_TAG="lsg";
  /**
	 delimits prior oil probability
  */
  public static final String P_OIL_TAG="probability_of_oil";
  /**
	 delimits prior gas probability
  */
  public static final String P_GAS_TAG="probability_of_gas";
  /**
   delimits prior lsg probability
  */
  public static final String P_LSG_TAG="probability_of_lsg";

  /**
	 delimits prior oil saturation
  */
  public static final String S_OIL_TAG="saturation_of_oil";
  /**
	 delimits prior gas saturation
  */
  public static final String S_GAS_TAG="saturation_of_gas";
  /**
   delimits prior lsg saturation
  */
  public static final String S_LSG_TAG="saturation_of_lsg";

  /**
	 delimits prior oil saturation std dev
  */
  public static final String SIGMA_S_OIL_TAG="sigma_saturation_of_oil";
  /**
	 delimits prior gas saturation std dev
  */
  public static final String SIGMA_S_GAS_TAG="sigma_saturation_of_gas";
  /**
   delimits prior lsg saturation std dev
  */
  public static final String SIGMA_S_LSG_TAG="sigma_saturation_of_lsg";

  /**
	Tags for brine properties
  */
  public static final String RHO_BRINE_TAG="rho_brine";
  public static final String SIGMA_RHO_BRINE_TAG="sigma_rho_brine";
  public static final String VP_BRINE_TAG="vp_brine";
  public static final String SIGMA_VP_BRINE_TAG="sigma_vp_brine";

 /**
	Tags for oil properties
  */
  public static final String RHO_OIL_TAG="rho_oil";
  public static final String SIGMA_RHO_OIL_TAG="sigma_rho_oil";
  public static final String VP_OIL_TAG="vp_oil";                 // JGCO2
  public static final String SIGMA_VP_OIL_TAG="sigma_vp_oil";
 /**
	Tags for gas properties
  */
  public static final String RHO_GAS_TAG="rho_gas";
  public static final String SIGMA_RHO_GAS_TAG="sigma_rho_gas";
  public static final String VP_GAS_TAG="vp_gas";
  public static final String SIGMA_VP_GAS_TAG="sigma_vp_gas";
  /**
	Tags for lsg properties
  */
  public static final String RHO_LSG_TAG="rho_lsg";
  public static final String SIGMA_RHO_LSG_TAG="sigma_rho_lsg";
  public static final String VP_LSG_TAG="vp_lsg";
  public static final String SIGMA_VP_LSG_TAG="sigma_vp_lsg";


  /**
   delimits prior NG (net-to-gross)
  */
  public static final String NG_TAG="net_to_gross";
  
  /**
   delimits prior NG std dev
  */
  public static final String SIGMA_NG_TAG="sigma_net_to_gross";

  /**
   delimits prior NG_hard (fraction of hard-rock in impermeable part)
  */
  public static final String NG_HARD_TAG="net_to_gross_hard";
  
  /**
   delimits prior NG_hard std dev
  */
  public static final String SIGMA_NG_HARD_TAG="sigma_net_to_gross_hard";

  /**
   delimits prior float_fraction
  */
  public static final String FLOAT_FRACTION_TAG="float_fraction";
  
  /**
   delimits prior float_fraction std dev
  */
  public static final String SIGMA_FLOAT_FRACTION_TAG="sigma_float_fraction";

  /**
   delimits prior time of layer top
  */
  public static final String T_TAG="time";
  /**
   delimits prior std dev in time of layer top
  */
  public static final String SIGMA_T_TAG="sigma_time";

  /**
   delimits prior time of layer top
  */
  public static final String DZ_TAG="thickness";
  /**
   delimits prior std dev in time of layer top
  */
  public static final String SIGMA_DZ_TAG="sigma_thickness";


  /**
   delimits prior time of last-layer bottom
  */
  public static final String T_BASE_TAG="t_base";
  /**
   delimits prior std dev in time of last-layer bottom
  */
  public static final String SIGMA_T_BASE_TAG="sigma_t_base";

  /**
   delimits model depth of layer, also rock-curve regression depth if 
   "depth_rock_curves" is not provided
  */
  public static final String DEPTH_TAG="depth";
  /**
	 delimits rock-curve regression depth of layer
  */
  public static final String DEPTH_ROCK_CURVES_TAG="depth_rock_curves";
  /**
   delimits low-frequency interval velocity
  */ 
  public static final String LFIV_TAG="LFIV";

  /**
	 Elements such as LFIV and depth can carry this tag,
	 meaning the deviates are to be synchronised with the layers above
  */
  public static final String SAME_AS_LAYER_ABOVE_TAG="same_as_layer_above";

  /**
   delimits std dev of depth of layer
  */
  public static final String SIGMA_DEPTH_TAG="sigma_depth";

  /**
   delimits std dev of low-frequency interval velocity
  */
  public static final String SIGMA_LFIV_TAG="sigma_LFIV";




  /**
   delimits group of variable prior indices in model traces
   eg, P_oil, and NG may vary. This requires an

   <PROPERTY_INDICES_TAG>
   <P_OIL_TAG>1</P_OIL_TAG>
   <NG_TAG>2</NG_TAG>
   </PROPERTY_INDICES_TAG>
   
   entry in the xml model file
  */
  public static final String PROPERTY_INDICES_TAG="property_indices";
 
  /**
	 This tag wraps the keyword defining the desired fluid density ordering
	 flavour. There are 3:<br>
	 <b>"full"</b>=fluids are density-ordered in permeable layers for the full 
	 extent of the model, even if two permeable layers are separated byan impermeable layer<br>
	 <b>"partial"</b>= fluids are density-ordered in adjacent permeable layers not separated by an impermeable layer<br>
	 <b>"none"</b>= all fluids are allowed in all layers, irrespective of density<br>
  */
  public static final String DENSITY_ORDERING_TAG="density_ordering";

  /**
	 keyword string for no density ordering in PROPERTY_INDICES_TAG field
  */
  public static final String NO_DENSITY_ORDERING_TAG="none";
  /**
	 keyword string for partial density ordering in PROPERTY_INDICES_TAG field
  */
  public static final String PARTIAL_DENSITY_ORDERING_TAG="partial";
  /**
	 keyword string for full density ordering in PROPERTY_INDICES_TAG field
  */
  public static final String FULL_DENSITY_ORDERING_TAG="full";

  /**
	 Tag to force models to be written in depth (last block in model traces
	 will become depth, not time)
  */
  public static final String OUTPUT_DEPTH_MODELS_TAG="depth_models";

  /**
	 Tag to indicate which layer number to hang depths from
  */
  public static final String MASTER_DEPTH_LAYER_NUMBER_TAG="master_depth_layer_number";

  /**
	 This tag delimits a string for the header field into which we will
	 write the realisation number in both model realisations and synthetic seismics
  */
  public static final String REALISATION_NUMBER_HEADER_FIELD_TAG="realisation_number_header_field";
  
  /**
	 This tag delimits the section controlling output
  */
  public static final String OUTPUT_SECTION_TAG="output";
  /**
	 This tag delimits the section controlling properties output
  */
  public static final String PROPERTIES_OUTPUT_TAG="properties";






  //********************** All the seismic stuff now ******************************

  /**
	 delimits block of seismic data
  */
  public static final String SEISMIC_DATA_TAG="seismic_data";

  /**
	 Tag for single trace minicube in parameter file
  */
  public static String CONSTANT_ANGLE_STACK_KEYWORD_TAG = "is-constant-angle-stack-with-angle-degrees";
  /**
	 Tag for PS data stack
  */
  public static String PS_STACK_KEYWORD_TAG = "is-PS-data";

  /*
	 delimits set of AVO descriptions
  */
  public static final String AVO_TERMS_TAG="AVO_terms";

  /**
	 Zero offset AVO fudge coefficient A
  */
  public static final String A_TAG="A";

  /**
	 Zero offset AVO fudge coefficient A std dev
  */
  public static final String SIGMA_A_TAG="sigma_A";

  /**
	 Quadratic-in-theta AVO fudge coefficient B
  */
  public static final String B_TAG="B";

  /**
	 Quadratic-in-theta AVO fudge coefficient B std dec
  */
  public static final String SIGMA_B_TAG="sigma_B";

  /**
	Stack tag
  */
  public static final String STACK_TAG="stack";
  /**
	Tag for the filename of the seismic stack
  */
  public static final String FILENAME_TAG="filename";

  /**
	 tag to identify the wavelet block
  */
  public static final String WAVELET_BLOCK_TAG="wavelet";
  /**
	 Tag to identify the user's override of the seismic error-acquisition range time-distance above
	 the top reflector (normally defaults to the wavelet coda)
  */
  public static final String WAVELET_T_ABOVE_TAG="time_above";
  /**
	 Tag to identify the user's override of the seismic error-acquisition range time-distance below
	 the bottom reflector (normally defaults to the wavelet precursor)
  */
  public static final String WAVELET_T_BELOW_TAG="time_below";
  /**
	 tag to identify the wavelet filename
  */
  public static final String WAVELET_FILENAME_TAG="filename";
  /**
	 tag to identify the RMS noise amplitude in a stack-data entry.
	 Second entry is legacy style
  */
  public static final String NOISE_RMS_TAG="noise_rms";
  public static final String LEGACY_NOISE_RMS_TAG="signal_to_noise";

  /**
	 Tag to identify the RMS noise amplitude of seismic stacks
	 if this property varies and is supplied in the model traces.
	 First is a legacy value, second the recommend usage.
  */
  public static final String NEAR_STACK_NOISE_RMS_TAG="near_stack_noise_rms";
  public static final String LEGACY_NEAR_STACK_NOISE_RMS_TAG="near_stack_signal_to_noise";
  /**
	 tag to identify the RMS noise amplitude of far stacks
	 if this property varies and is supplied in the model traces
	 First is a legacy value, second the recommend usage.
  */
  public static final String FAR_STACK_NOISE_RMS_TAG="far_stack_signal_to_noise";
  public static final String LEGACY_FAR_STACK_NOISE_RMS_TAG="far_stack_noise_rms";

  /**
	tag for near offset distance
  */
  public static final String NEAR_OFFSET_TAG="min_offset";

  /**
	 tag for far offset distance
  */
  public static final String FAR_OFFSET_TAG="max_offset";

  /**	
		tag for event time
  */
  public static final String REFLECTOR_TIME_TAG="reflector_time";

  /**	
		tag for event stacking velocity
  */
  public static final String STACK_VELOCITY_TAG="stack_velocity";

  /** 
	  Model into which seismic constants ought to be stuffed
  */
  OneDSeismicConvolutionalModel SM;

 

  // WORKING variables

  /**
	 reference to current set of ModelPriorDescriptions we are parsing for
  */
  ModelPriorDescriptions _MPD;

  /**
	 reference to current Model
  */
  Model _m;

  /**
	 Useful hashmap of properties and their name & values
  */
  HashMap propMap;

  private MutableInteger mutI;

  // XML parsing for xmlPaths
  XPath xPath = XPathFactory.newInstance().newXPath();

  /**
	 Make a parser for the named file, which will populate certain
	 fields of the Model m, with the desired certain level of
	 verbosity and handle errors using the ErrorHandler
  */
  public ModelDescriptionParser(String filename, Model m, 
								int verbosity, ErrorHandler EH)
  {

	this.verbosity = verbosity;
	this.EH = EH;
	this._m = m;
	this.SM = (OneDSeismicConvolutionalModel) m.SM;

	// working Integer for map retrieval etc
	this.mutI=new MutableInteger(0);

	propMap = new HashMap();

	propMap.put(T_TAG,new propertyTagInfo("t", T_TAG, false));
	propMap.put(SIGMA_T_TAG,new propertyTagInfo("sigma_t",SIGMA_T_TAG,false));
	propMap.put(DZ_TAG,new propertyTagInfo("dz",DZ_TAG,false));
	propMap.put(SIGMA_DZ_TAG,new propertyTagInfo("sigma_dz",SIGMA_DZ_TAG,false));
	propMap.put(NG_TAG,new propertyTagInfo("NG",NG_TAG,true));
	propMap.put(SIGMA_NG_TAG,new propertyTagInfo("sigma_NG",SIGMA_NG_TAG,true));
	if (_m._SI.MF.ENABLE_HARD_ROCK_MODEL)
	  {
		propMap.put(NG_HARD_TAG,new propertyTagInfo("NG_hard",NG_HARD_TAG,true));
		propMap.put(SIGMA_NG_HARD_TAG,new propertyTagInfo("sigma_NG_hard",SIGMA_NG_HARD_TAG,true));	
	  }
	if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		propMap.put(FLOAT_FRACTION_TAG,new propertyTagInfo("float_fraction",FLOAT_FRACTION_TAG,true));
		propMap.put(SIGMA_FLOAT_FRACTION_TAG,new propertyTagInfo("sigma_float_fraction",SIGMA_FLOAT_FRACTION_TAG,true));
	  }
	propMap.put(P_OIL_TAG,new propertyTagInfo("P_oil",P_OIL_TAG,true));
	propMap.put(P_GAS_TAG,new propertyTagInfo("P_gas",P_GAS_TAG,true));
	propMap.put(P_LSG_TAG,new propertyTagInfo("P_lsg",P_LSG_TAG,true));
	propMap.put(S_OIL_TAG,new propertyTagInfo("S_oil",S_OIL_TAG,true));
	propMap.put(S_GAS_TAG,new propertyTagInfo("S_gas",S_GAS_TAG,true));
	propMap.put(S_LSG_TAG,new propertyTagInfo("S_lsg",S_LSG_TAG,true));
	propMap.put(SIGMA_S_OIL_TAG,new propertyTagInfo("sigma_S_oil",SIGMA_S_OIL_TAG,false));
	// Useful variable fluids...(JGCO2)
	propMap.put(VP_OIL_TAG,new propertyTagInfo("vp_oil",VP_OIL_TAG,false));                   
	propMap.put(SIGMA_VP_OIL_TAG,new propertyTagInfo("sigma_vp_oil",SIGMA_VP_OIL_TAG,false)); 
	propMap.put(VP_GAS_TAG,new propertyTagInfo("vp_gas",VP_GAS_TAG,false));                   
	propMap.put(SIGMA_VP_GAS_TAG,new propertyTagInfo("sigma_vp_gas",SIGMA_VP_GAS_TAG,false)); 
	propMap.put(VP_LSG_TAG,new propertyTagInfo("vp_lsg",VP_LSG_TAG,false));                   
	propMap.put(SIGMA_VP_LSG_TAG,new propertyTagInfo("sigma_vp_lsg",SIGMA_VP_LSG_TAG,false));
	
	propMap.put(RHO_OIL_TAG,new propertyTagInfo("vp_oil",RHO_OIL_TAG,false));                   
	propMap.put(SIGMA_RHO_OIL_TAG,new propertyTagInfo("sigma_vp_oil",SIGMA_RHO_OIL_TAG,false)); 
	propMap.put(RHO_GAS_TAG,new propertyTagInfo("vp_gas",RHO_GAS_TAG,false));                   
	propMap.put(SIGMA_RHO_GAS_TAG,new propertyTagInfo("sigma_vp_gas",SIGMA_RHO_GAS_TAG,false)); 
	propMap.put(RHO_LSG_TAG,new propertyTagInfo("vp_lsg",RHO_LSG_TAG,false));                   
	propMap.put(SIGMA_RHO_LSG_TAG,new propertyTagInfo("sigma_vp_lsg",SIGMA_RHO_LSG_TAG,false));  
	// ...variable fluids

	propMap.put(SIGMA_S_GAS_TAG,new propertyTagInfo("sigma_S_gas",SIGMA_S_GAS_TAG,false));
	propMap.put(SIGMA_S_LSG_TAG,new propertyTagInfo("sigma_S_lsg",SIGMA_S_LSG_TAG,false));
	propMap.put(DEPTH_TAG,new propertyTagInfo("depth",DEPTH_TAG,false));
	propMap.put(DEPTH_ROCK_CURVES_TAG,new propertyTagInfo("depth_rock_curves",DEPTH_ROCK_CURVES_TAG,false));
	propMap.put(LFIV_TAG,new propertyTagInfo("LFIV",LFIV_TAG,false));
	propMap.put(SIGMA_DEPTH_TAG,new propertyTagInfo("sigma_depth",SIGMA_DEPTH_TAG,false));
	propMap.put(SIGMA_LFIV_TAG,new propertyTagInfo("sigma_LFIV",SIGMA_LFIV_TAG,false));
	propMap.put(T_BASE_TAG,new propertyTagInfo("t_base",T_BASE_TAG,false));
	propMap.put(SIGMA_T_BASE_TAG,new propertyTagInfo("sigma_t_base",SIGMA_T_BASE_TAG,false));


	// 200104
	propMap.put(RHO_BRINE_TAG,new propertyTagInfo("rho_brine",RHO_BRINE_TAG,false));
	propMap.put(SIGMA_RHO_BRINE_TAG,new propertyTagInfo("sigma_rho_brine",SIGMA_RHO_BRINE_TAG,false));
	propMap.put(VP_BRINE_TAG,new propertyTagInfo("vp_brine",VP_BRINE_TAG,false));
	propMap.put(SIGMA_VP_BRINE_TAG,new propertyTagInfo("sigma_vp_brine",SIGMA_VP_BRINE_TAG,false));

	// noise levels in seismic stacks

	// Legacy (backcompatbility with old XML files only)
	propMap.put("near_stack_signal_to_noise",new propertyTagInfo("near_stack_signal_to_noise",LEGACY_NOISE_RMS_TAG,false));
	propMap.put("far_stack_signal_to_noise",new propertyTagInfo("far_stack_signal_to_noise",LEGACY_NOISE_RMS_TAG,false));
	// new
	propMap.put("near_stack_noise_rms",new propertyTagInfo("near_stack_noise_rms",NOISE_RMS_TAG,false));
	propMap.put("far_stack_noise_rms",new propertyTagInfo("far_stack_noise_rms",NOISE_RMS_TAG,false));

	// AVO constants
	propMap.put(A_TAG,new propertyTagInfo("A_avo",A_TAG,false));
	propMap.put(SIGMA_A_TAG,new propertyTagInfo("sigma_A_avo",SIGMA_A_TAG,false));
	propMap.put(B_TAG,new propertyTagInfo("B_avo",B_TAG,false));
	propMap.put(SIGMA_B_TAG,new propertyTagInfo("sigma_B_avo",SIGMA_B_TAG,false));

	// Simple fluid contacts info
	propMap.put("contact_depth",new propertyTagInfo("contact_depth","contact_depth",false));
	propMap.put("sigma_contact_depth",new propertyTagInfo("sigma_contact_depth","sigma_contact_depth",false));
  }

  /** 
	  general call to read layer descriptions into relevant fields in model.
	  Model must exist
  */
  
   boolean readInModelAndSeismicDescriptions(String xml_filename, Model m)
  {
	boolean success;
	this.xml_filename = xml_filename;

	 
	DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
	 
	try 
	  {
		DocumentBuilder builder = factory.newDocumentBuilder();
		File fin = new File(xml_filename);
		document = builder.parse(fin);
		
		// (1) Get layer data for properties 
		success = getLayerData(m);
		if (!success) return false;
		
		// (1.1) Possible vertical stacking
		success = getStackedModelsInfo(m);
		if (!success) return false;

		// (1.2) Check possible fluid units
		success = checkFluidUnitsAndContacts(m);
		if (!success) return false;
		
		// (2) read seismic stuff
		success = getOneDSeismicParameters();
		if (!success) return false;

		// (2) get rock_properties and fluid data
		RockandFluidDataParser RDP = new RockandFluidDataParser(xml_filename,verbosity,EH);
		success = RDP.getRockandFluidData(m);
		if (!success) return false;

		// Now we have stored in each _MPD.MetaInfo[i].xxxx.varies
		// whether the property in question (xxxxx)  needs to vary
		// This call finds the index of these varying properties in another
		// section (wrapped in PROPERTY_INDICES_TAG) of the xml file
		success = obtainVariablePriorData(m);
		if (!success) return false;

		// (3) get a range of generic parameters
		if (!obtainFluidDensityAndProb(m)) return false;
		obtainRealisationNumberHeaderField(m);
		obtainDepthModelFlags(m);
		obtainOutputPropertiesList(m);

		// build some useful things for handling Layer synchronisation
		setupLayerSynchronisation(m);
		
		if (this.verbosity > 0) 
		  {
			m.initialiseLayerSynching();
			m.MPD.printPriorInfo(false);
		  }
	  }
	catch (Exception ex) 
	  {
		ex.printStackTrace();
		return false;
	  };
	return true;
  }
  
  boolean getLayerData(Model m) throws SAXException
  {
	int i,j,N;
	int LFIV_master_layer_no = 0;
	int depth_master_layer_no = 0;
	int float_fraction_master_layer_no = 0;
	int depth_rock_curves_master_layer_no = 0;
	int bottom_layer_length=0,top_layer_length=0,layer_length=0,total_length;
	boolean need=false;
	String msg_frag,s;
	Node n;
	NodeList fl_list;
	Element reservoir_elem=null;
	Element fluid_elem=null;
	ModelPriorDescriptions.MetaInfo MI;
	FluidUnitsAndContacts _FUC;
	NodeList layerlist=null;
	try 
	  {
		// ok, look for the MODEL_DESCRIPTION_TAG chunk
		try
		  {
			layerlist = (NodeList)xPath.evaluate("/inversion/"+MODEL_DESCRIPTION_TAG,document,XPathConstants.NODESET);
		  }
		catch (Exception e) 
		  {
			_m._SI.EH.println("Parsing exception: "+e);
		  }

		if (layerlist.getLength()<1) 
		  throw new SAXException("No " + MODEL_DESCRIPTION_TAG + " entry in: " + this.xml_filename);

		// There should be only one entry. lets take the first
		// & look for the *LAYER_TAG chunks

		// In the new format there will also be TOP_LAYER_TAG entries. Augment the nodelist by these
		NodeList nlist_top = ((Element)layerlist.item(0)).getElementsByTagName(TOP_LAYER_TAG);
		top_layer_length = nlist_top.getLength();

		if (top_layer_length>1) throw new SAXException("Unworkable model: must be only 1 "+TOP_LAYER_TAG+
													   " entry in: " + this.xml_filename);

		NodeList nlist_middle;
		if (top_layer_length == 0)
		  {
			// must be old-style XML, with layer,layer...bottom_layer sequence.
			nlist_middle = ((Element)layerlist.item(0)).getElementsByTagName("layer");
		  }
		else
		  {
			// new XML, with layer,middle_layer...bottom_layer sequence.
			nlist_middle = ((Element)layerlist.item(0)).getElementsByTagName(MIDDLE_LAYER_TAG);	
		  }
		layer_length = nlist_middle.getLength();

		// Last layer has BOTTOM_LAYER_TAG entries. Augment the nodelist by these
		NodeList nlist_bottom = ((Element)layerlist.item(0)).getElementsByTagName(BOTTOM_LAYER_TAG);
		bottom_layer_length = nlist_bottom.getLength();
		
		if (bottom_layer_length != 1) throw new SAXException("Unworkable model: must be 1 "+BOTTOM_LAYER_TAG+
												   " entry in: " + this.xml_filename);

		total_length = top_layer_length + bottom_layer_length + layer_length;

		if (total_length<2) throw new SAXException("Unworkable model: less than total of 2 layer-type entries in: " 
												   + this.xml_filename);


		// Now we know how big the model has to be, set up the model fully
		m.setUpFully(total_length,EH,verbosity);		

		// look for fluid-contacts info now, as if changes how we intrpret fluid probabilities & identifications
		if (!getFluidUnitsAndContacts(m)) return false;

		m.MS.completeSetUpUsingFluidContactInfo();

		_FUC = m.fluidUnitsAndContacts;
		_FUC.buildRequiredFluidsTable();

		_MPD = m.MPD;
 
		// loop over layers
		for (i=0;i<total_length;i++)
		  {
			// Messy hack for back compatibility...
			Element elem;
			if (i == (total_length-1)) 
			  {
				// bottom layer
				elem = (Element)nlist_bottom.item(0); 
			  }           
			else if (i>0) 
			  {
				// middle layers 
				elem = (Element)nlist_middle.item((top_layer_length==0)? i : i-1); 
			  }
			else // i==0
			  {
				if (top_layer_length > 0)
				  {
					elem = (Element)nlist_top.item(0); 
				  }
				else
				  {
					elem = (Element)nlist_middle.item(0); 
				  }
			  }


			// allocate
			_MPD.MetaInfo[i] = new ModelPriorDescriptions.MetaInfo();
			MI = _MPD.MetaInfo[i];

			msg_frag=new String(" entry for layer " + (i+1));

			// Get layer name, if any
			MI.layer_name = extractValue(elem,NAME_TAG,NAME_TAG + msg_frag);

			// First look for SHALE_MEMBER_TAG ENTRIES
			// in this layer:                 ((ELEMENT)nlist.item(i))
			// look for                       SHALE_MEMBER
			// type entries. Take the first:  item(0)
			Element el = (Element) (elem.getElementsByTagName(SHALE_MEMBER_TAG).item(0));
			if (el != null)
			  {
				n = el.getElementsByTagName(NAME_TAG).item(0);
				if (n == null)
				  {
					throw new SAXException("In: " + xml_filename + ", entry " + 
										   SHALE_MEMBER_TAG +" is missing " + NAME_TAG);
				  }
				String value = XmlUtils.getTextNodeValue(n);
				if (value != null)
				  {
					MI.mud_rock_member = value;				  
				  }
				else
				  {
					throw new SAXException("In: " + xml_filename + ", must have intelligible " 
										   + SHALE_MEMBER_TAG + msg_frag);
				  }
			  }
			// is this layer  synchronised ???
			NamedNodeMap nnmap = elem.getAttributes();
			if (nnmap != null)
			  {
				Attr nv = (Attr)(nnmap.getNamedItem(SYNCHRONISE_TAG));
				if (nv != null)
				  {
					String attr_val = nv.getValue();
					if (attr_val != null && attr_val.equals(TRUE_TAG))
					  {
						MI.synchronised = true;
					  }
				  }
			  }

			// hard rocks
			if (m._SI.MF.ENABLE_HARD_ROCK_MODEL)
			  {
				// ok, look for the hard rocks...
				// CASE 1. Fixed, fluid-insenstive hard rock
				NodeList list_hard = elem.getElementsByTagName(HARD_SHALE_MEMBER_TAG);
				if (list_hard != null && list_hard.getLength()>0)
				  {
					// looks like there's a second impermeable rock
					Element elem_hard = (Element)(list_hard.item(0));				
					MI.hard_shale_rock_member 
					  = this.extractValue(elem_hard,NAME_TAG,NAME_TAG + msg_frag);
				  }
			
				// CASE 2: list of fluid-sensitive rocks
				list_hard = elem.getElementsByTagName("fluid_dependent_nonreservoir_hard_endmembers");
				if (list_hard != null && (list_hard.getLength() > 0))
				  {
					if (MI.hard_shale_rock_member != null) 
					  {
						throw new SAXException("In: " + xml_filename + ", layer i="+(i+1)+", "+MI.layer_name+
											   ", entry "+HARD_SHALE_MEMBER_TAG+
											   " cannot coexist with fluid_dependent_nonreservoir_hard_endmembers" + 
											   " Use one or the other!");
					  }
					N=0;
					MI.fluid_dep_hard_shale_rock_members = new String[FluidStates.FluidTypes.Nfluids]; 
					el = (Element) (elem.getElementsByTagName("brine_case").item(0));
					if ((el != null) && ((s = XmlUtils.getTextNodeValue(el)) != null))
					  {
						MI.fluid_dep_hard_shale_rock_members[FluidStates.FluidTypes.brine] = s; N++;
					  }
					el = (Element) (elem.getElementsByTagName("lsg_case").item(0));
					if ((el != null) && ((s = XmlUtils.getTextNodeValue(el)) != null))
					  {
						MI.fluid_dep_hard_shale_rock_members[FluidStates.FluidTypes.lsg] = s; N++;
					  }
					el = (Element) (elem.getElementsByTagName("oil_case").item(0));
					if ((el != null) && ((s = XmlUtils.getTextNodeValue(el)) != null))
					  {
						MI.fluid_dep_hard_shale_rock_members[FluidStates.FluidTypes.oil] = s; N++;
					  }
					el = (Element) (elem.getElementsByTagName("gas_case").item(0));
					if ((el != null) && ((s = XmlUtils.getTextNodeValue(el)) != null))
					  {
						MI.fluid_dep_hard_shale_rock_members[FluidStates.FluidTypes.gas] = s; N++;
					  }
					  
					if (N==0)
					  {
						throw new SAXException("In: " + xml_filename + ", layer i="+(i+1)+", "+MI.layer_name+
											   ", entry fluid_dependent_nonreservoir_hard_endmembers" + 
											   "require at least one fluid-case named rock type!");
					  }
				  }
	
				// Get NG_hard if the mixture is non-trivial, either CASE 1 or CASE 2 above
				if ((MI.hard_shale_rock_member != null) || (MI.fluid_dep_hard_shale_rock_members !=null))
				  {
					// We must have a  net-to-gross_hard
					MI.NG_hard = extractValueAndIndex(elem,NG_HARD_TAG,msg_frag);
					
					// And maybe net-to-gross uncertainty (defaults to zero)
					MI.sigma_NG_hard = extractValueAndIndex(elem,SIGMA_NG_HARD_TAG,msg_frag);
					// default null value of prior (sigma_NG=0)
					if (MI.sigma_NG_hard == null) MI.sigma_NG_hard = new ModelPriorDescriptions.PriorDetails(false);

					if (MI.NG_hard.varies || (MI.NG_hard.value > 0) || MI.sigma_NG_hard.varies || (MI.sigma_NG_hard.value > 0)) 
					  {
						// it's a non-trivial mixture that may have hard-rock
						MI.has_hard_rock_mixture = true;
					  }				  			  
				  }
			  }

			// ok, look for the SAND_MEMBER
			NodeList reservoir_list = elem.getElementsByTagName(SAND_MEMBER_TAG);
			if (reservoir_list != null && reservoir_list.getLength()>0)
			  {
				// looks like there's a permeable rock
				reservoir_elem = (Element)(reservoir_list.item(0));				
				MI.sand_rock_member 
				  = this.extractValue(reservoir_elem,NAME_TAG,
										   NAME_TAG + msg_frag);
			  }

			// If rock could be a mixture (sand name is provided), we will need fluid properties
			// if the mixture is non-trivial
			if (MI.sand_rock_member != null) 
			  {
				// There seems to be sand. We must have a  net-to-gross
				MI.NG = extractValueAndIndex(elem,NG_TAG,msg_frag);
				
				// And maybe net-to-gross uncertainty (defaults to zero)
				MI.sigma_NG = extractValueAndIndex(elem,SIGMA_NG_TAG,msg_frag);
				// default null value of prior (sigma_NG=0)
				if (MI.sigma_NG == null)
				  {
					MI.sigma_NG = new ModelPriorDescriptions.PriorDetails(false);
				  }

				if (MI.NG.varies || (MI.NG.value > 0) || MI.sigma_NG.varies || (MI.sigma_NG.value > 0)) 
				  {
					// it's a non-trivial mixture that may have sand
					MI.is_permeable=true;
				  }

				// Check also how this affects hard rocks....
				if (m._SI.MF.ENABLE_HARD_ROCK_MODEL)
				  {
					if ((MI.hard_shale_rock_member != null) || (MI.fluid_dep_hard_shale_rock_members !=null))  
					  {
						if (!MI.NG.varies && (MI.NG.value==1.0) && (MI.sigma_NG.value==0.0))
						  {
							// Fixed, pure sand. Hard rocks won't matter
							MI.has_hard_rock_mixture = false;
						  }
					  }
				  }
			  }

			// If it's non-trivial, we'll need the fluid properties
			// of (1) *named* possible fluids in the layer descriptions
			// that have (2) non-zero default prior probabilities

			// update contacts table based on is_permeable flag
			_FUC.updateRequiredFluidsTable(i,MI.is_permeable);

			if (MI.is_permeable)
			  {
				// (1) OIL. Get the oil name and default prior probability
				// look for an oil entry
				need = _FUC.contactsRequireThisFluid(i,FluidStates.FluidTypes.oil);
				MI.oil_name
				  = this.extractNamedValue(reservoir_elem,OIL_TAG,
										   OIL_TAG + msg_frag);
				if (MI.oil_name != null)
				  // I've named an oil, so oil is one possible fluid. What's its prior probabilities? 
				  {
					fluid_elem = (Element)(reservoir_elem.getElementsByTagName(OIL_TAG).item(0));
					MI.P_oil = extractValueAndIndex(fluid_elem,P_OIL_TAG,msg_frag);
					if (MI.P_oil.varies || MI.P_oil.value > 0 || need)
					  {		
						// we really must have the properties of this oil
						MI.oil_substitution_possible = true;
					  }
				  } 
				else if (need)
				  {
					throw new SAXException("In: " + xml_filename + 
										   ", fluid contacts require oil specification in layer "+(i+1));
				  }
			
				// (2) GAS. Get the gas name and default prior probability
				// looks for a gas entry
				need = _FUC.contactsRequireThisFluid(i,FluidStates.FluidTypes.gas);
				MI.gas_name
				  = this.extractNamedValue(reservoir_elem,GAS_TAG,
										   GAS_TAG + msg_frag);
				if (MI.gas_name != null)
				  // I've named an gas, so gas is one possible fluid. What's its prior probabilities? 
				  {
					fluid_elem = (Element)(reservoir_elem.getElementsByTagName(GAS_TAG).item(0));
					MI.P_gas = extractValueAndIndex(fluid_elem,P_GAS_TAG,msg_frag);
					if (MI.P_gas.varies|| MI.P_gas.value > 0 || need)
					  {		
						MI.gas_substitution_possible = true;
					  }
				  }
				else if (need)
				  {
					throw new SAXException("In: " + xml_filename + 
										   ", fluid contacts require gas specification in layer "+(i+1));
				  }

				// (3) LSG. Get the lsg name and default prior probability
				// look for an an lsg entry
				need = _FUC.contactsRequireThisFluid(i,FluidStates.FluidTypes.lsg);
				MI.lsg_name
				  = this.extractNamedValue(reservoir_elem,LSG_TAG,
										   LSG_TAG + msg_frag);
				if (MI.lsg_name != null)
				  // I've named an lsg, so lsg is one possible fluid. What's its prior probabilities? 
				  {
					fluid_elem = (Element)(reservoir_elem.getElementsByTagName(LSG_TAG).item(0));
					MI.P_lsg = extractValueAndIndex(fluid_elem,P_LSG_TAG,msg_frag);
					if (MI.P_lsg.varies|| MI.P_lsg.value > 0 || need)
					  {		
						MI.lsg_substitution_possible = true;
					  }
				  }
				else if (need)
				  {
					throw new SAXException("In: " + xml_filename + 
										   ", fluid contacts require lsg specification in layer "+(i+1));
				  }

				// (4) BRINE. Get the brine name and default prior probability
				// N.B. we *always* need the brine to get density from porosity in permeable rocks
				MI.brine_name
				  = this.extractNamedValue(reservoir_elem,BRINE_TAG, BRINE_TAG + msg_frag);
				if (MI.brine_name == null)
				  {
					throw new SAXException("In: " + xml_filename + ", must have " 
										   + NAME_TAG + msg_frag);
				  }

				// Brine fluid properties
				// density
				MI.rho_b
				  = extractValueAndIndex(reservoir_elem,RHO_BRINE_TAG,msg_frag);
				MI.sigma_rho_b
				  = extractValueAndIndex(reservoir_elem,SIGMA_RHO_BRINE_TAG,msg_frag);
				// vp
				MI.vp_b
				  = extractValueAndIndex(reservoir_elem,VP_BRINE_TAG,msg_frag);
				MI.sigma_vp_b
				  = extractValueAndIndex(reservoir_elem,SIGMA_VP_BRINE_TAG,msg_frag);

				// Also, any non-trivial floating grain behaviour:
				if (m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
				  {
					// We have reservoir rock => we might have floating grains
					MI.float_fraction = extractValueAndIndex(reservoir_elem,FLOAT_FRACTION_TAG,msg_frag);
					// And maybe float_fraction uncertainty
					MI.sigma_float_fraction = extractValueAndIndex(reservoir_elem,SIGMA_FLOAT_FRACTION_TAG,msg_frag);
					if ((MI.float_fraction!=null && (MI.float_fraction.varies || (MI.float_fraction.value > 0)))
						|| (MI.sigma_float_fraction!=null && (MI.sigma_float_fraction.varies || (MI.sigma_float_fraction.value > 0))))
					  {
						// we have to handle float stuff
						MI.needs_float_capability = true;
					  }
				  }
			  }

			if (m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
			  {
				// We *always* need float_fraction elements for sanity, even if the feature is not used, so:
				if (MI.float_fraction == null)
				  {
					// Assume zero for default.
					MI.float_fraction = new ModelPriorDescriptions.PriorDetails(true);
				  }
				// default null value of prior (sigma_float_fraction=0)
				if (MI.sigma_float_fraction == null)
				  {
					MI.sigma_float_fraction = new ModelPriorDescriptions.PriorDetails(false);
				  }
			  }
			
			// Now all the routine stuff
			// (1) mandatory entries			
			MI.t 
			  = extractMandatoryValueAndIndex(elem,T_TAG,msg_frag);
			// inputs are in ms. convert to sec.
			MI.t.value = 0.001 * MI.t.value;

			MI.sigma_t 
			  = extractMandatoryValueAndIndex(elem,SIGMA_T_TAG,msg_frag);
			// inputs are in ms. convert to secs.
			MI.sigma_t.value = 0.001 * MI.sigma_t.value;

			if (MI.sigma_t.value <= 0)
			  {
				throw new SAXException("In: " + xml_filename + ", layer "
								 + i + ", must have positive sigma_t: " 
								 + SIGMA_T_TAG);
			  }
			
			MI.depth 
			  = extractMandatoryValueAndIndex(elem,DEPTH_TAG,msg_frag);

			if (MI.depth.synched_above)
			  {
				MI.depth.master_layer_number = depth_master_layer_no;
			  }
			else 
			  {
				// next layer's depth will mirror this layer if it's synched_above
				depth_master_layer_no = i; 
			  }

			// Optional (back-compatibility) loading depth
			MI.depth_rock_curves
			  = extractValueAndIndex(elem,DEPTH_ROCK_CURVES_TAG,msg_frag);
			if (MI.depth_rock_curves == null)
			  {
				// Not provided: use the model depth as a surrogate:
				MI.depth_rock_curves = MI.depth;
			  }
			else
			  {
				if (MI.depth_rock_curves.synched_above)
				  {
					MI.depth_rock_curves.master_layer_number = depth_rock_curves_master_layer_no;
				  }
				// next layer's depth_rock_curves will mirror this layer if it's synched_above
				else 
				  {
					depth_rock_curves_master_layer_no = i;
				  }
			  }
	
			MI.LFIV 
			  = extractMandatoryValueAndIndex(elem,LFIV_TAG,msg_frag);

			if (MI.LFIV.synched_above)
			  {
				MI.LFIV.master_layer_number = LFIV_master_layer_no;
			  }
			// next layer's LFIV will mirror this layer if it's synched_above
			else 
			  {
				LFIV_master_layer_no = i;
			  }

			// (2) some optional things

			MI.sigma_depth 
			  = extractValueAndIndex(elem,SIGMA_DEPTH_TAG,msg_frag);
			if (MI.sigma_depth==null)
			  { 
				// default prior with zero sigma
				MI.sigma_depth = new ModelPriorDescriptions.PriorDetails(false);
			  }
			// uncertainty in trend depth
			if (MI.depth_rock_curves == MI.depth)
			  {
				// duplicate depths, make sigmas duplicate too
				MI.sigma_depth_rock_curves = MI.sigma_depth;
			  }
			else
			  {
				// Distinct loading depth. Set it's sigma to zero (required for later code).
				MI.sigma_depth_rock_curves = new ModelPriorDescriptions.PriorDetails(false);
			  }


			MI.sigma_LFIV 
			  = extractValueAndIndex(elem,SIGMA_LFIV_TAG,msg_frag);
			if (MI.sigma_LFIV==null)
			  {
				// default prior with zero sigma
				MI.sigma_LFIV = new ModelPriorDescriptions.PriorDetails(false);
			  }


			


			// (2) Saturations & uncertainties
			if (MI.oil_substitution_possible)
			  {
				fluid_elem = (Element)(reservoir_elem.getElementsByTagName(OIL_TAG).item(0));
				MI.S_oil
				  = extractMandatoryValueAndIndex(fluid_elem,S_OIL_TAG,msg_frag);
				MI.sigma_S_oil
				  = extractValueAndIndex(fluid_elem,SIGMA_S_OIL_TAG,msg_frag);
				// Useful variable fluids (JGCO2)
				MI.vp_oil
				  = extractValueAndIndex(fluid_elem,VP_OIL_TAG,msg_frag);
				MI.sigma_vp_oil
				  = extractValueAndIndex(fluid_elem,SIGMA_VP_OIL_TAG,msg_frag);
				// variable priors should override fluid defaults 
				MI.vp_oil_VA = (MI.vp_oil!=null) || (MI.sigma_vp_oil!=null);   
				MI.rho_oil
				  = extractValueAndIndex(fluid_elem,RHO_OIL_TAG,msg_frag);
				MI.sigma_rho_oil
				  = extractValueAndIndex(fluid_elem,SIGMA_RHO_OIL_TAG,msg_frag);
				// variable priors should override fluid defaults 
				MI.rho_oil_VA = (MI.rho_oil!=null) || (MI.sigma_rho_oil!=null);        
			  }
			if (MI.gas_substitution_possible)
			  {
				fluid_elem = (Element)(reservoir_elem.getElementsByTagName(GAS_TAG).item(0));
				MI.S_gas
				  = extractMandatoryValueAndIndex(fluid_elem,S_GAS_TAG,msg_frag);
				MI.sigma_S_gas
				  = extractValueAndIndex(fluid_elem,SIGMA_S_GAS_TAG,msg_frag);
				// Useful variable fluids (JGCO2)
				MI.vp_gas
				  = extractValueAndIndex(fluid_elem,VP_GAS_TAG,msg_frag);
				MI.sigma_vp_gas
				  = extractValueAndIndex(fluid_elem,SIGMA_VP_GAS_TAG,msg_frag);
				// variable priors should override fluid defaults JGCO2
				MI.vp_gas_VA = (MI.vp_gas!=null) || (MI.sigma_vp_gas!=null);
				MI.rho_gas
				  = extractValueAndIndex(fluid_elem,RHO_GAS_TAG,msg_frag);
				MI.sigma_rho_gas
				  = extractValueAndIndex(fluid_elem,SIGMA_RHO_GAS_TAG,msg_frag);
				// variable priors should override fluid defaults JGCO2
				MI.rho_gas_VA = (MI.rho_gas!=null) || (MI.sigma_rho_gas!=null);     
			  }
			if (MI.lsg_substitution_possible)
			  {
				fluid_elem = (Element)(reservoir_elem.getElementsByTagName(LSG_TAG).item(0));
				MI.S_lsg
				  = extractMandatoryValueAndIndex(fluid_elem,S_LSG_TAG,msg_frag);
				MI.sigma_S_lsg
				  = extractValueAndIndex(fluid_elem,SIGMA_S_LSG_TAG,msg_frag);
				// variable priors should override fluid defaults JGCO2
				MI.vp_lsg
				  = extractValueAndIndex(fluid_elem,VP_LSG_TAG,msg_frag);
				MI.sigma_vp_lsg
				  = extractValueAndIndex(fluid_elem,SIGMA_VP_LSG_TAG,msg_frag);
				// variable priors should override fluid defaults JGCO2
				MI.vp_lsg_VA = (MI.vp_lsg!=null) || (MI.sigma_vp_lsg!=null);
				MI.rho_lsg
				  = extractValueAndIndex(fluid_elem,RHO_LSG_TAG,msg_frag);
				MI.sigma_rho_lsg
				  = extractValueAndIndex(fluid_elem,SIGMA_RHO_LSG_TAG,msg_frag);
				// variable priors should override fluid defaults JGCO2
				MI.rho_lsg_VA = (MI.rho_lsg!=null) || (MI.sigma_rho_lsg!=null);       
			  }

			// (3) isopachs
			MI.dz
			  = extractValueAndIndex(elem,DZ_TAG,msg_frag);
			// must have +ve uncertainties on isopachs if they are present
			if (MI.dz != null)
			  {
				MI.sigma_dz
				  = extractMandatoryValueAndIndex(elem,SIGMA_DZ_TAG,msg_frag);

				// precompute for efficiency
				if (MI.sigma_dz.value > 0) MI.log_sigma_dz = Math.log(MI.sigma_dz.value);

				// Count genuine isopachs only if they are not trivial layer-synchronisation flags
				if (!((MI.dz.value==0) && (MI.sigma_dz.value==0))) 
				  {
					_MPD.n_isopachs++;
				  }
				if (MI.sigma_dz.value < 0)
				  {
					new SAXException("In: " + xml_filename + ", layer "
									 + i + ", must have non-negative isopach uncertainty: " 
									 + SIGMA_DZ_TAG);
				  }
			  }

			// (4) last layer *must* have bottom
			if (i == total_length-1)
			  {
				_MPD.t_base
				  = extractMandatoryValueAndIndex(elem,T_BASE_TAG,msg_frag);
				// inputs are in ms. convert to secs.
				_MPD.t_base.value = 0.001 * _MPD.t_base.value;
				_MPD.sigma_t_base
				  = extractMandatoryValueAndIndex(elem,SIGMA_T_BASE_TAG,msg_frag);
				// inputs are in ms. convert to secs.
				_MPD.sigma_t_base.value = 0.001 * _MPD.sigma_t_base.value;
				if (_MPD.sigma_t_base.value <= 0)
				  {
					new SAXException("In: " + xml_filename + ", layer "
									 + i + ", must have positive sigma_t: " 
									 + SIGMA_T_BASE_TAG);
				  }
			  }

			// fix units etc
			performUnitsAndSignCorrections(i,MI);

			// a few fluid checks
			m.FS.checkPriorFluidProbabilities(i);
		  }

		_FUC.checkFluidConsistencyWithinUnits();

	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }

	

	return true;
  }

  



  /**
	 Bunch of operations required to deal with unusual units or signs
	 conventions in XML file defining the priors
  */

  private void performUnitsAndSignCorrections(int layer_no, ModelPriorDescriptions.MetaInfo MI)
  {
	// (1) times may need to be forced into seconds....
	if (_m._SI.MODEL_TRACE_TIME_UNITS_AUTO_DETECTION)
	  {
		if (MI.t.value > _m._SI.TIME_UNITS_AUTO_DETECTION_TIME_MAX)
		  {
			MI.t.value *= 1.0e-3;
		  }
		if (MI.sigma_t.value > _m._SI.TIME_UNITS_AUTO_DETECTION_TIME_MAX)
		  {
			MI.sigma_t.value *= 1.0e-3;	
		  }
	  }
	// (2) depths may need to be forced to positive values
	if (_m._SI.MAKE_MODEL_TRACE_PRIOR_DEPTHS_POSITIVE)
	  {
		MI.depth.value = Math.abs(MI.depth.value);
		MI.depth_rock_curves.value = Math.abs(MI.depth_rock_curves.value);
	  }
  }





  
  /** 
	  Obtain and check indices of prior elemnts in model that vary with each trace
	  from master xml file
  */
  boolean obtainVariablePriorData(Model m)
  {
	String msg,name=null;
	int i,j,k,MaxIndx=0;
	ModelPriorDescriptions.MetaInfo MI;
	ModelPriorDescriptions.PriorDetails PD=null;
	RockandFluidData.MaterialProperties _mp;
	RockandFluidData.RegressionSet RFDregSet=null;

	try 
	  {
		// ok, look for the PROPERTY_INDICES_TAG chunk
		NodeList nlist=null;
		try
		  {
			nlist = (NodeList)xPath.evaluate("/inversion/"+PROPERTY_INDICES_TAG,document,XPathConstants.NODESET);
		  }
		catch (Exception e) 
		  {
			_m._SI.EH.println("Parsing exception: "+e);
		  }
		if (nlist.getLength()<1) 
		  throw new SAXException("No " + PROPERTY_INDICES_TAG + " entry in: " + this.xml_filename);

		Element e = (Element)nlist.item(0);
		
		msg=new String(PROPERTY_INDICES_TAG);

		_MPD.PriorIndexInfo = new ArrayList();

		// Find things that vary. They are of several types:

		// (1) Actual layer properties; need to loop over layers
		ModelFormat _MF = m._SI.MF;
		for (i=0;i<_MPD.MetaInfo.length;i++)
		  {
			MI = _MPD.MetaInfo[i];
			
			extractVariableIndex(e,i,MI, MI.t,T_TAG+LINK_SUFFIX, msg,false);
			if (MI.t.varies) m.MPD.prior_times_vary = true;
			extractVariableIndex(e,i,MI,MI.sigma_t,SIGMA_T_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.dz,DZ_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_dz,SIGMA_DZ_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.NG,NG_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_NG,SIGMA_NG_TAG+LINK_SUFFIX,msg,true);
			if (_MF.ENABLE_HARD_ROCK_MODEL)
			  {
				extractVariableIndex(e,i,MI,MI.NG_hard,NG_HARD_TAG+LINK_SUFFIX,msg,true);
				extractVariableIndex(e,i,MI,MI.sigma_NG_hard,SIGMA_NG_HARD_TAG+LINK_SUFFIX,msg,true);
			  }
			if (_MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
			  {
				extractVariableIndex(e,i,MI,MI.float_fraction,
									 FLOAT_FRACTION_TAG+LINK_SUFFIX,msg,true);
				extractVariableIndex(e,i,MI,MI.sigma_float_fraction,
									 SIGMA_FLOAT_FRACTION_TAG+LINK_SUFFIX,msg,true);
			  }
			extractVariableIndex(e,i,MI,MI.P_oil,P_OIL_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.P_gas,P_GAS_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.P_lsg,P_LSG_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.S_oil,S_OIL_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.S_gas,S_GAS_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.S_lsg,S_LSG_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_S_oil,SIGMA_S_OIL_TAG+LINK_SUFFIX,msg,true);
			// variables fluids for eg CO2...
			extractVariableIndex(e,i,MI,MI.vp_oil,VP_OIL_TAG+LINK_SUFFIX,msg,true); 
			extractVariableIndex(e,i,MI,MI.rho_oil,RHO_OIL_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.vp_gas,VP_GAS_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.rho_gas,RHO_GAS_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.vp_lsg,VP_LSG_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.rho_lsg,RHO_LSG_TAG+LINK_SUFFIX,msg,true);
			//....variable fluids

			extractVariableIndex(e,i,MI,MI.sigma_S_gas,SIGMA_S_GAS_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_S_lsg,SIGMA_S_LSG_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.depth,DEPTH_TAG+LINK_SUFFIX,msg,true);
			if (MI.depth_rock_curves != MI.depth)
			  {
				extractVariableIndex(e,i,MI,MI.depth_rock_curves,DEPTH_ROCK_CURVES_TAG+LINK_SUFFIX,msg,true);
			  }
			extractVariableIndex(e,i,MI,MI.LFIV,LFIV_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_depth,SIGMA_DEPTH_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_LFIV,SIGMA_LFIV_TAG+LINK_SUFFIX,msg,true);

			// Variable fluid properties 
			extractVariableIndex(e,i,MI,MI.rho_b,RHO_BRINE_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_rho_b,SIGMA_RHO_BRINE_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.vp_b,VP_BRINE_TAG+LINK_SUFFIX,msg,true);
			extractVariableIndex(e,i,MI,MI.sigma_vp_b,SIGMA_VP_BRINE_TAG+LINK_SUFFIX,msg,true);

			// these should only occur once, in the last layer
			// base times are tagged by t (an extra entry) so use the t-tag
			if (i == _MPD.MetaInfo.length-1)
			  {
				extractVariableIndex(e,i,MI,_MPD.t_base,T_TAG+LINK_SUFFIX,msg,false);
				if (_MPD.t_base.varies) m.MPD.prior_times_vary = true;
				// safety checks
				if (MI.t.varies && !_MPD.t_base.varies)
				  throw new SAXException("if "+T_TAG+" varies, so must the "+T_BASE_TAG+
										 " entry in: " + this.xml_filename);
				if (!MI.t.varies && _MPD.t_base.varies)
				  throw new SAXException("if "+T_BASE_TAG+" varies, so must the "+T_TAG+
										 " entries in: " + this.xml_filename);
			  }

			
		  }

		MI = _MPD.MetaInfo[0];

		// (2) See if the seismic noise levels vary
		// There will be an entire block of the
		// input model traces devoted to the stack rms-noise-level (deliberate redundancy)
		// The HashMap _MPD.PriorIndexInfo needs to have an entry only for the first of these
		
		// near stack noise rms
		StackData st;
		String noise_link_tag;
		st = _m.SM.getNearstackData();
		if (st !=null) 
		  {
			PD = _m.SM.getNearstackData().noise_rms;
			if (PD!=null) 
			  {   
				if (PD.varies)
				  {
					// we must find the variable - index of the property MI
					noise_link_tag = getUsedAlternateTag(e, LEGACY_NEAR_STACK_NOISE_RMS_TAG+LINK_SUFFIX,
														    NEAR_STACK_NOISE_RMS_TAG+LINK_SUFFIX);
					if (noise_link_tag==null) throw new SAXException("No link specified for near-stack varying noise rms in: " 
																	 + this.xml_filename);
					extractVariableIndex(e,0,MI,PD,noise_link_tag, msg,true);				
					_m.SM.doSNchecks = true;
				  }
			  }
		  }
		// far stack noise rms
		st = _m.SM.getFarstackData();
		if (st !=null)
		  {
			PD = _m.SM.getFarstackData().noise_rms;
			if (PD!=null) 
			  {
				if (PD.varies)
				  {
					// we must find the variable - index of the property MI
					noise_link_tag = getUsedAlternateTag(e, LEGACY_FAR_STACK_NOISE_RMS_TAG+LINK_SUFFIX,
														    FAR_STACK_NOISE_RMS_TAG+LINK_SUFFIX);
					if (noise_link_tag==null) throw new SAXException("No link specified for far-stack varying noise rms in: " 
																	 + this.xml_filename);
					extractVariableIndex(e,0,MI,PD,noise_link_tag, msg,true);
					_m.SM.doSNchecks = true;
				  }
			  }
		  }

		// (3) See if AVO constants are need in model traces
		if (_m.SM.A_avo != null) extractVariableIndex(e,0,MI,_m.SM.A_avo, A_TAG+LINK_SUFFIX, msg,true);
		if (_m.SM.sigma_A_avo != null) extractVariableIndex(e,0,MI,_m.SM.sigma_A_avo, SIGMA_A_TAG+LINK_SUFFIX, msg,true);
		if (_m.SM.B_avo != null) extractVariableIndex(e,0,MI,_m.SM.B_avo, B_TAG+LINK_SUFFIX, msg,true);
		if (_m.SM.sigma_B_avo != null) extractVariableIndex(e,0,MI,_m.SM.sigma_B_avo, SIGMA_B_TAG+LINK_SUFFIX, msg,true);
		

		// (4) Simple contacts
		if (_m._SI.doSimpleFluidUnitsAndContacts)
		  {
			// read these from unrolled layer_no fields
			for (i=0,k=0;i<	_m.fluidUnitsAndContacts.fluidUnits.length; i++)
			  {
				FluidUnitsAndContacts.FluidUnit fu = _m.fluidUnitsAndContacts.fluidUnits[i];
				for (j=0; j< fu.contactInfo.length; j++)
				  {
					FluidUnitsAndContacts.ContactInfo cI = 	fu.contactInfo[j];
					extractVariableIndex(e,k,MI,cI.depthPD,"contact_depth"+LINK_SUFFIX, msg,true);
					extractVariableIndex(e,k,MI,cI.sigma_depthPD,"sigma_contact_depth"+LINK_SUFFIX, msg,true);
					k++;
				  }
			  }
		  }
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }

	// _MPD.PriorIndexInfo.trimToSize();
  
	return true;
  }


  /** 
	Obtain flag from xml file controlling the mode of fluid density ordering
  */
  boolean obtainFluidDensityAndProb(Model m)
  {
	try
	  {
		Node n = (Node)xPath.evaluate("/inversion/output",document,XPathConstants.NODE);
		String s = (String)xPath.evaluate(DENSITY_ORDERING_TAG,n,XPathConstants.STRING);
		s=s.trim();

		// ok, look for the DENSITY_ORDERING_TAG chunk
		// this tag may not be there
		if (s.equals("")) 
		  {
			// assume full density ordering
			m._SI.density_ordering_mode=FluidStates.fullDensityOrdering;
			return true;
		  }
		else
		  {
			// if it's there, impose the desired ordering
			if (s.equals(NO_DENSITY_ORDERING_TAG))
			  {
				m._SI.density_ordering_mode = FluidStates.noDensityOrdering;
			  }
			else if (s.equals(PARTIAL_DENSITY_ORDERING_TAG))
			  {
				m._SI.density_ordering_mode = FluidStates.partialDensityOrdering;
			  }
			else if (s.equals(FULL_DENSITY_ORDERING_TAG))
			  {
				m._SI.density_ordering_mode = FluidStates.fullDensityOrdering;
			  }
			else
			  {
				try 
				  {
					throw new SAXException("Density ordering must be one of (" 
										   + NO_DENSITY_ORDERING_TAG +", "
										   + PARTIAL_DENSITY_ORDERING_TAG + ", "
										   + FULL_DENSITY_ORDERING_TAG+ ") for tag "
										   + DENSITY_ORDERING_TAG
										   +" entries in: " + this.xml_filename);
				  }
				catch (SAXException ex)
				  {
					ex.printStackTrace();
					return false;
				  }
			  }
		  }

		// Now fluid probability flavour
		s = (String)xPath.evaluate("take_fluid_probabilities_as_post_ordering",n,XPathConstants.STRING);
		if (s.trim().equals(TRUE_TAG)) 
		  {
			_m._SI.useFluidProbabilitiesAsPostOrdering = true;
		  }
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }
	return true;
  }

  /** 
	Obtain flag from xml file controlling whether models are dumped in depth, and what the
	master depth layer number is
  */
  void obtainDepthModelFlags(Model m) throws SAXException
  {
	try 
	  {
		// (1) whether to dump in depths

		// ok, look for the OUTPUT_DEPTH_MODELS_TAG chunk
		Node n = (Node)xPath.evaluate("/inversion/output",document,XPathConstants.NODE);
		String s = (String)xPath.evaluate(OUTPUT_DEPTH_MODELS_TAG,n,XPathConstants.STRING);
		s = s.trim();
		if (s==TRUE_TAG) m._SI.outputDepthModels = true;

		// (2) What layer to hang the depths from
		s = (String)xPath.evaluate(MASTER_DEPTH_LAYER_NUMBER_TAG,n,XPathConstants.STRING);
		if (s!="")
		  {
			try
			  {
				m.master_depth_layer_number = Integer.parseInt(s);
				if (m.master_depth_layer_number<1 || m.master_depth_layer_number > _m.n_layers)
				  {
					throw new SAXException("In: " + xml_filename + ", Master depth layer "
										   +m.master_depth_layer_number+" is too high or low");
				  }
			  }
			catch (NumberFormatException e)
			  {
				throw new SAXException("In: " + xml_filename + 
									   ", Not an integer: tag " + MASTER_DEPTH_LAYER_NUMBER_TAG);
			  }
		  }
		else
		  {
			throw new SAXException("In: " + xml_filename + 
								   ", No entry for tag: " + MASTER_DEPTH_LAYER_NUMBER_TAG);
		  }
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }
	return;
  }

  /**
	 Obtain the realisation number header field
  */
  void obtainRealisationNumberHeaderField(Model m)
  {
	try
	  {
		// ok, look for the REALISATION_NUMBER_HEADER_FIELD_TAG
		Node n = (Node)xPath.evaluate("/inversion/output",document,XPathConstants.NODE);
		String s = (String)xPath.evaluate(REALISATION_NUMBER_HEADER_FIELD_TAG,n,XPathConstants.STRING);
		s = s.trim();
		if (s.equals("")) return; // this tag may not be there. Default will then apply

		// else set the headerField variable to what we find and check it's sensible
		_m._SI.MF.realisationNumberHeaderField = s;
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }
	return;
  }

  /**
	 Obtain the properties subselection if present
  */
  void obtainOutputPropertiesList(Model m) throws SAXException
  {
	// ok, look for the OUTPUT_SECTION_TAG
	NodeList nlist = null;
	try
	  {
		nlist = (NodeList)xPath.evaluate("/inversion/"+OUTPUT_SECTION_TAG,document,XPathConstants.NODESET);
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }


	// this tag may not be there. Default will then apply
	if (nlist.getLength()<1) 
	  {
		throw new SAXException("In: " + xml_filename + 
							   ", No output section tagged as "+OUTPUT_SECTION_TAG);
	  }
	else
	  {
		// There should be only one entry. Lets take the first
		// & look for the PROPERTIES_OUTPUT_TAG
		NodeList nlist2 = ((Element)nlist.item(0)).getElementsByTagName(PROPERTIES_OUTPUT_TAG);
		if (nlist2.getLength()<1) 
		  {
			return; // take defaults
		  }
		else
		  {
			Element elem = (Element)nlist2.item(0); 
			// scan attributes to switch on fields
			NamedNodeMap nnmap = elem.getAttributes();
			if (nnmap != null)
			  {
				if (m._SI.MF.complete_field_names.length > 32)
				  {
					throw new SAXException("Model format has too many fields to flag in trid header word%n");
				  }
				// check all fields except t & t_base to possibly suppress them
				for (int i=0;i<m._SI.MF.complete_field_names.length-2;i++)
				  {
					Attr nv = (Attr)(nnmap.getNamedItem(m._SI.MF.complete_field_names[i]));
					if (nv != null)
					  {
						String attr_val = nv.getValue();
						if (attr_val != null)
						  {
							if (!attr_val.equals(TRUE_TAG))
							  {
								m._SI.MF.suppressField(i);
							  }
						  }
					  }
				  }
			  }
		  }
	  }
	return;
  }


  /**
	 Set up synchronisations for rock matrices between adjacent
	 layers.  Synchronised rock matrices can occur iff the layers have
	 matching end-member rocks. This method builds up the HashMap
	 model.MPD.synchMap, which points to the lowest indexed layer
	 above which the current layer is effectively synchronised to.
  */
  private void setupLayerSynchronisation(Model model)
  {
	int i;
	// nearly always we don't use this. Default to save room
	model.MPD.synchMap = new OpenIntIntHashMap(1);
	ModelPriorDescriptions.MetaInfo MI,MJ;

	// loop over layers in reverse order to get minimal, rooted HashMap of synchronisations
	// (synched layers will point at the smallest-numbered layer which is their image)
	for (i=_MPD.MetaInfo.length-1;i>=0;i--)
	  {
		MI = _MPD.MetaInfo[i];
		if (MI.synchronised)
		  {
			// find matching layer searching forwards
			for (int j=i+1;j<_MPD.MetaInfo.length;j++)
			  {
				MJ = _MPD.MetaInfo[j];
				if (MJ.synchronised)
				  {
					// shales must match
					if (MI.mud_rock_member.equals(MJ.mud_rock_member))
					  {
						// FIXME HARD Rocks Too!...add condition

						// so must sands if it's permeable
						if ((MI.is_permeable && MJ.is_permeable && 
							 MI.sand_rock_member.equals(MJ.sand_rock_member))
							|| (!MI.is_permeable && !MJ.is_permeable))
						  {
							// add new mapping to HashMap of synchronised layers
							model.MPD.synchMap.put(j,i);
						  }
					  }
				  }
			  }
		  }
	  }
	// get rid of this HashMap if nothing's required
	if (model.MPD.synchMap.size()==0) 
	  {
		model.MPD.synchMap = null;
	  }
	else
	  {
		if (verbosity>1) EH.print("%nLayer synchronisations summary:%n%n");
		model.MPD.synch_copies=new int[_MPD.MetaInfo.length];
		for (i=0;i<=_MPD.MetaInfo.length;i++)
		  {
			// Integer J = ((Integer)model.MPD.synchMap.get(mutI.set(i)));
			if (model.MPD.synchMap.containsKey(i))
			  {
				int j = model.MPD.synchMap.get(i);
				model.MPD.synch_copies[j]++;
				if (verbosity>1) 
				  {
					EH.printf("  Layer %d synchronised with %d%n%n",i, j);
				  }
			  } 
		  }
	  }
	
	// Finally, ensure the synchronisation of LFIV and depth (and possibly float_fraction) terms is
	// consistent with overall layer synchronisations
	// Ensure any _MPD.MetaInfo[i].(LFIV|depth|depth_rock_curves|float_fraction).master_layer_number values point at the
	// *lowest* master_layer_number that is consistent with the two flavours
	// of Synchronisation (1) overall layer synchs, (2) individual LFIV, depth, depth_rock_curves, layer synchs
	// Combinations of (1) and (2) can create dependencies that spread across multiple layers
	// (draw a diagram)
	int nStart = _MPD.MetaInfo.length-1;

	for (i=nStart; i>=1; i--)
	  {
		MI = _MPD.MetaInfo[i];
		MI.depth.master_layer_number = getLowestMirrorIndex("depth",i);
		MI.LFIV.master_layer_number = getLowestMirrorIndex("LFIV",i);

		if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			MI.float_fraction.master_layer_number = getLowestMirrorIndex("float_fraction",i);
			if (MI.float_fraction.master_layer_number != i)
			  {
				_MPD.MetaInfo[MI.float_fraction.master_layer_number].hasFloat_fractionDependents = true;
			  }
		  }

		if (MI.depth_rock_curves != MI.depth)
		  {
			MI.depth_rock_curves.master_layer_number = getLowestMirrorIndex("depth_rock_curves",i);
			if (MI.depth_rock_curves.master_layer_number != i)
			  {
				_MPD.MetaInfo[MI.depth_rock_curves.master_layer_number].hasDepthDependents = true;
			  }
		  }
		if (MI.depth.master_layer_number != i)
		  {
			_MPD.MetaInfo[MI.depth.master_layer_number].hasDepthDependents = true;
		  }
		if (MI.LFIV.master_layer_number != i)
		  {
			_MPD.MetaInfo[MI.LFIV.master_layer_number].hasLFIVDependents = true;
		  }
	  }
	if (verbosity>3) 
	  {
		for (i=nStart; i>=1; i--)
		  {
			MI = _MPD.MetaInfo[i];
			if (_MPD.synchMap!=null)
			  {
				if (MI.LFIV.synched_above || _MPD.synchMap.containsKey(i))
				  {
					EH.println("   Layer " + i + " LFIV mirrors layer "+ MI.LFIV.master_layer_number);
				  }
				else
				  EH.println("   Layer " + i + " LFIV is free");
			  }
		  }
		for (i=nStart; i>=1; i--)
		  {
			MI = _MPD.MetaInfo[i];
			if (MI.depth_rock_curves != MI.depth)
			  {
				if (_MPD.synchMap!=null)
				  {
					if (MI.depth_rock_curves.synched_above || _MPD.synchMap.containsKey(i))
					  EH.println("   Layer " + i + " depth_rock_curves mirrors layer "+ MI.depth_rock_curves.master_layer_number);
					else
					  EH.println("   Layer " + i + " depth_rock_curves is free");
				  }				
			  }
		  }
		for (i=nStart; i>=1; i--)
		  {
			MI = _MPD.MetaInfo[i];
			if (_MPD.synchMap!=null)
			  {
				if (MI.depth.synched_above || _MPD.synchMap.containsKey(i))
				  EH.println("   Layer " + i + " depth mirrors layer "+ MI.depth.master_layer_number);
				else
				  EH.println("   Layer " + i + " depth is free");
			  }
		  }
		if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
		  {
			for (i=nStart; i>=1; i--)
			  {
				MI = _MPD.MetaInfo[i];
				if (_MPD.synchMap!=null)
				  {
					if (MI.float_fraction.synched_above || _MPD.synchMap.containsKey(i))
					  {
						EH.println("   Layer " + i + " float_fraction mirrors layer "+ MI.float_fraction.master_layer_number);
					  }
					else
					  EH.println("   Layer " + i + " float_fraction is free");
				  }
			  }
		  }
	  }
	return;

  }

  // Recursive routines to get lowest mirror indices for LFIV, depth, float_fraction

  private int getLowestMirrorIndex(String priorDetailsName, int layer_number)
  {
	int j = layer_number,j_s;
	Field f=null;
	ModelPriorDescriptions.PriorDetails PD=null,PD_above=null;

	Class modelPriorDescriptionsClass = ((Object)_MPD.MetaInfo[layer_number]).getClass();

	try
	  {
		f = modelPriorDescriptionsClass.getField(priorDetailsName);
		PD = (ModelPriorDescriptions.PriorDetails) 
		  f.get((Object)_MPD.MetaInfo[layer_number]);
	  }
	// these exceptions never happen
	catch (NoSuchFieldException e) {}
	catch (IllegalAccessException e) {};

	if (PD.synched_above)
	  {
		j = layer_number - 1;
	  }
	if (_MPD.synchMap != null)
	  {
		if (_MPD.synchMap.containsKey(layer_number))
		  {
			j_s = _MPD.synchMap.get(layer_number);

			j = Math.min(j,j_s);
			// set this in case it wasn't set
			// PD.synched_above = true;

			// point any forks to the deeper prong
			if (PD.synched_above && (j_s < layer_number-1))
			  {
				try
				  {
					PD_above =  (ModelPriorDescriptions.PriorDetails) 
					  f.get((Object)_MPD.MetaInfo[layer_number-1]);
				  }
				// these exceptions never happen
				catch (IllegalAccessException e) {};
				PD_above.synched_above = true;
				PD_above.master_layer_number = j_s;
			  }
		  }
	  }

	if (j==layer_number) 
	  return j;
	else 
	  return getLowestMirrorIndex(priorDetailsName, j);
  }










  /** 
	Set the variable index of property PD , tagged by tag, within the xml node
	e. If needed, it must be found. If we fail, throw an exception. Must be
	kept identical to similar method in RockandFluidDataParser
  */
  void extractVariableIndex(Element e, int layer_no, ModelPriorDescriptions.MetaInfo MI,
									ModelPriorDescriptions.PriorDetails PD, 
									String tag,String element_name, boolean test) throws SAXException
  
  {
	int index;
	
	if (PD != null && PD.varies)
	  {
		// we must find the variable - index of the property MI
		Node n = e.getElementsByTagName(tag).item(0);
		if (n!= null)
		  {
			String str = XmlUtils.getTextNodeValue(n);
			if (str != null)
			  {
				try
				  {
					index = Integer.parseInt(str);
				  }
				catch (NumberFormatException ex)
				  {
					throw new SAXException("In: " + xml_filename + 
										   ", Not an integer: tag " + tag + " element " + element_name);
				  }
				// we seem to have got a reasonable index for tag_name in this Element
				
				if (index>0) 
				  {		
					// add to useful ArrayList of changing things
					PD.index = index;
					_MPD.PriorIndexInfo.add(new ModelPriorDescriptions.propertyInfo(layer_no,PD));

					// EH.println("tag = "+ tag + ", PriorIndexInfo SIZE = "+_MPD.PriorIndexInfo.size()+ " block = "+PD.index);
					if (MI.t.varies)
					  {
						// t, t_base have the largest index
						if (test)
						  {
							// is a property other than time, so must have smaller index
							if (index >= MI.t.index)
							  throw new SAXException("if " + tag + " varies, its index must be less than "
													 + MI.t.index + ", "+ this.xml_filename);
						  }
						else
						  {
							// must be t or t_base
							if (index != MI.t.index)
							  throw new SAXException("since " + tag + " varies, its index *must* be "
													 + MI.t.index + ", "+ this.xml_filename);
						  }
					  }
					return;
				  }
				else
				  throw new SAXException("In: " + xml_filename + 
										 ", negative index for " 
										 + tag + " in element " + element_name);
			  }
		  }
		// if we got here, we're in trouble
		throw new SAXException("In: " + xml_filename + 
							   ", can't find an variable index for " 
							   + tag + " in element " + element_name);
	  }
  }

  /** 
	  Return the text value from the named Element e,  as tagged by
	  tag_name. If we find it but don't understand it
	  send failure message to the error handler.
	  If not found, return null
  */
  private String extractValue(Element e,String tag_name,String message) throws SAXException
  {
	Node n= e.getElementsByTagName(tag_name).item(0);
	if (n!= null)
	  {
		String value = XmlUtils.getTextNodeValue(n);
		if (value != null)
		  {
			// we seem to have got a reasonable value for tag_name in this Element
			return value;
		  }
		else
		  {
			// we found the <tag_name> entry, but couldn't understand it
			throw new SAXException("In: " + xml_filename + ", unintelligible entry: " + message);
		  }
	  }
	// we didn't find the <tag_name> entry, so return null
	return null;
  }

  /** 
	  Return the text value from the first NAME_TAG subelement of
	  Element e, where e is tagged by tag_name. If we find it but
	  don't understand it send failure message to the error handler.
	  If not found, return null
  */
  private String extractNamedValue(Element e,String tag_name,String message) throws SAXException
  {
	Element el = (Element) (e.getElementsByTagName(tag_name).item(0));
	if (el!= null)
	  {
		Node n = (el.getElementsByTagName(NAME_TAG)).item(0);
		if (n == null)
		  {
			throw new SAXException("In: " + xml_filename + ", entry " + 
								   tag_name +" is missing " + NAME_TAG+ ": "+ message);
		  }
		String value = XmlUtils.getTextNodeValue(n);
		if (value != null)
		  {
			// we seem to have got a reasonable value for tag_name in this Element
			return value;
		  }
		else
		  {
			// we found the <tag_name> entry, but couldn't understand it
			throw new SAXException("In: " + xml_filename + ", unintelligible entry: " + message);
		  }
	  }
	// we didn't find the <tag_name> entry, so return null
	return null;
  }

  /** 
	  Return which of two alternative names 
	  in the XML file are available from element e in the tree, where
	  s1 and s2 anre tags expected to be found under e.
	  Tries s1 first, then s2, or returns null if neither found.
  */
  String getUsedAlternateTag(Element e,String s1, String s2)
  {
	NodeList n;

	if ((e.getElementsByTagName(s1).getLength() > 0))
	  {
		return s1;
	  }
	else
	  {
		// try alternative
		if ((e.getElementsByTagName(s2).getLength() > 0))
		  {
			return s2;
		  }
	  }
	return null;
  }

  /** 
	  Return which of two alternative names 
	  in the XML file are available from element e in the tree, where
	  s1 s2 are keys to propertyTagInfo propMap entries
	  Tries s1 first, then s2, or returns null if neither found.
  */
  String getUsedPDAlternateTag(Element e,String s1, String s2)
  {
	NodeList n;

	propertyTagInfo info = (propertyTagInfo) propMap.get(s1);
	if (info != null && (e.getElementsByTagName(info.tag).getLength() > 0))
	  {
		return s1;
	  }
	else
	  {
		// try alternative
		info = (propertyTagInfo) propMap.get(s2);
		if (info != null && (e.getElementsByTagName(info.tag).getLength() > 0))
		  {
			return s2;
		  }
	  }
	return null;
  }



  /** 
	  Return the text value from the named Element e,  as tagged by
	  tag_name. If we find it but don't understand it
	  send failure message to the error handler.
	  If not found, return null
	  Also see if the VARIES tag seems to be enclosed, and return this.
	  This looks for xml entries of the form

	  <tag_name VARIES_TAG="true|false">value</tag_name>

	  within the named Element e. Only the first occurrence is recognised
  */
  private ModelPriorDescriptions.PriorDetails 
	extractValueAndIndex(Element e,String var_name, String msg) throws SAXException
  {
	Node n_sub;
	Attr nv;
	NamedNodeMap nnmap;
	String attr_val;
	ModelPriorDescriptions.PriorDetails PD=null;
	double value;

	propertyTagInfo info = (propertyTagInfo) propMap.get(var_name);
	String tag_name = info.tag;
	boolean NormalisedValue = info.normalised;
	String message = new String(tag_name + " " + msg);

	Node n= e.getElementsByTagName(tag_name).item(0);
	

	if (n!= null)
	  {
		String str = XmlUtils.getTextNodeValue(n);
		if (str != null)
		  {
			try
			  {
				value = Double.parseDouble(str);
			  }
			catch (NumberFormatException ex)
			  {
				throw new SAXException("In: " + xml_filename + ", Not a number : " + message);
			  }
			// we seem to have got a reasonable value for tag_name in this Element
			PD = new ModelPriorDescriptions.PriorDetails(info.normalised);
			if (NormalisedValue)
			  {
				if (value<0) throw new SAXException("In: " + xml_filename + ", value < 0 : " + message);
				if (value>1) throw new SAXException("In: " + xml_filename + ", value > 1 : " + message);
			  }
			PD.value = value;
			nnmap = n.getAttributes();
			if (nnmap != null)
			  {
				nv = (Attr)(nnmap.getNamedItem(VARIES_TAG));
				if (nv != null)
				  {
					attr_val = nv.getValue();
					if (attr_val != null && attr_val.equals(TRUE_TAG))
					  {
						PD.varies=true;
					  }
				  }
			  }
			if (nnmap != null)
			  {
				nv = (Attr)(nnmap.getNamedItem(SAME_AS_LAYER_ABOVE_TAG));
				if (nv != null)
				  {
					attr_val = nv.getValue();
					if (attr_val != null && attr_val.equals(TRUE_TAG))
					  {
						PD.synched_above = true; 
					  }
				  }
			  }
		  }
		else
		  {
			// we found the <tag_name> entry, but couldn't understand it
			throw new SAXException("In: " + xml_filename + ", unintelligible entry: " + message);
		  }
		return PD;
	  }
	else
	  {
		// we didn't find the <tag_name> entry, so return null
		return null;
	  }
  }


  /** 
	  Return the text value from the named Element e,  as tagged by
	  tag_name. If we find it but don't understand it
	  send failure message to the error handler.
	  If not found, throw an error
	  Also see if the VARIES tag seems to be enclosed, and return this.
	  This looks for xml entries of the form

	  <tag_name>value<VARIES_TAG/></tag_name>

	  within the named Element e. Only the first occurrence is recognised
  */
  private ModelPriorDescriptions.PriorDetails 
	extractMandatoryValueAndIndex(Element e,String var_name,String message) throws SAXException
  {
	ModelPriorDescriptions.PriorDetails pd;
	
	pd = extractValueAndIndex(e,var_name,message);

	if (pd == null)
	  {
		propertyTagInfo info = (propertyTagInfo) propMap.get(var_name);
		String tag_name = info.tag;
		throw new SAXException("In: " + xml_filename + 
							   ", cannot find " + tag_name + " entry: " + message);
	  }

	return pd;
  }

  /********************************************************/
  /*****                Seismic parsing methods     *******/
  /********************************************************/

  /**
   * Function to read all seismic one-D terms from xml filename,
   with given level of verbosity
   * @return a <code>boolean</code> value
   */
  boolean getOneDSeismicParameters()
  {
	NodeList nlist_wavelet; 
	NamedNodeMap nnmap=null;
	Attr nv=null;
	String attr_val=null;

	try 
	  {
		// (1) Extract parameters 
		
		// get the chunk of seismic data
		NodeList nlist = null;
		try
		  {
			nlist = (NodeList)xPath.evaluate("/inversion/"+SEISMIC_DATA_TAG,document,XPathConstants.NODESET);
		  }
		catch (Exception e) 
		  {
			_m._SI.EH.println("Parsing exception: "+e);
		  }

		if (nlist.getLength()!=1) throw new SAXException("No " + SEISMIC_DATA_TAG + " entry in this file");
		
		// there should be only one occurence of <SEISMIC_DATA_TAG>, so the first is what we want
		Element seismic_data_element = ((Element)nlist.item(0));

		
		// ok, look for the AVO_TERMS_TAG chunk
		nlist = seismic_data_element.getElementsByTagName(AVO_TERMS_TAG);
		if (nlist.getLength()!=1) throw new SAXException("No " + AVO_TERMS_TAG + " entry in this file");

		// is the synchronised model A=B switched on
		Element	elem = ((Element)nlist.item(0));
		nnmap = elem.getAttributes();
		if (nnmap != null)
		  {
			nv = (Attr)(nnmap.getNamedItem("synchronise_B_from_A"));
			if (nv != null)
			  {
				attr_val = nv.getValue();
				if (attr_val != null && attr_val.equals(TRUE_TAG))
				  {
					SM.A_B_avo_synchronised = true;
				  }
			  }
		  }


		// there should be only one occurence of <AVO_TERMS_TAG>, so the first is what we want
		Element entry = ((Element)nlist.item(0));
	
		// get A
		SM.A_avo = extractValueAndIndex(entry,A_TAG," (AVO section)");
		if (SM.A_avo==null) SM.A_avo = new ModelPriorDescriptions.PriorDetails((double)1.0);


		// get sigma_AS
		SM.sigma_A_avo = extractValueAndIndex(entry,SIGMA_A_TAG," (AVO section)");
		if (SM.sigma_A_avo==null) SM.sigma_A_avo = new ModelPriorDescriptions.PriorDetails((double)0);

		// B
		SM.B_avo = extractValueAndIndex(entry,B_TAG," (AVO section)");
		if (SM.B_avo==null) SM.B_avo = new ModelPriorDescriptions.PriorDetails((double)1.0);

		// sigma_B
		SM.sigma_B_avo = extractValueAndIndex(entry,SIGMA_B_TAG," (AVO section)");
		if (SM.sigma_B_avo==null) SM.sigma_B_avo = new ModelPriorDescriptions.PriorDetails((double)0);

		// Create the actual properties. these should be set to the means of the
		// distributions for A & B defined above
		// add these to the map of propertynames too
		SM.A = new LayerProperties.Property(SM.A_avo.value,false);
		SM.B = new LayerProperties.Property(SM.B_avo.value,false);
		if (SM.A_B_avo_synchronised)
		  {
			SM.B_avo.setValue(SM.A_avo.value);
			SM.sigma_B_avo.setValue(SM.sigma_A_avo.value);
			SM.B.setPosValue(SM.A_avo.value);
		  }
		_m.completePropertyNameMap.put(SM.A,"A"); 
		_m.completePropertyNameMap.put(SM.B,"B");

		if (verbosity > 0)
		  {
			EH.print("%nSEISMIC MODEL DETAILS:%nAVO constants:%n");
			EH.println("  A,B synchronisation = "+SM.A_B_avo_synchronised);
			EH.println("  A = "+ SM.A_avo.value);
			EH.println("  sigma_A = "+ SM.sigma_A_avo.value);
			EH.println("  B = "+ SM.B_avo.value);
			EH.println("  sigma_B = "+ SM.sigma_B_avo.value);
		  }

		// (2) get stacks
		// look for the STACK_TAG chunk
		// We must have either one or two stacks.
		// If there's only one, we'll call it `near' (whatever the angle be)
		nlist = seismic_data_element.getElementsByTagName(STACK_TAG);
		if (nlist.getLength()==0) throw new SAXException("No " + STACK_TAG + " entry in " +xml_filename);
		if (nlist.getLength()>2) throw new SAXException("Can only handle 2 stacks: " +xml_filename);
		// loop over stacks (1 or 2)
		for (int i=0 ; i<nlist.getLength() ; i++)
		  {
			entry = ((Element)nlist.item(i));
			StackData stackdata =  new StackData();

			// New piece: look for constant angle attribute
			nnmap = entry.getAttributes();
			nv = (Attr)(nnmap.getNamedItem(CONSTANT_ANGLE_STACK_KEYWORD_TAG));
			if (nv != null)
			  { 
				attr_val = nv.getValue();
				if (attr_val != null)
				  {
					stackdata.theta_specified = Double.parseDouble(attr_val);
					if (stackdata.theta_specified >= 0)
					  {
						if (stackdata.theta_specified >= 90)
						  {
							throw new SAXException("Specified constant-angle attribute for seismic stack "+ xml_filename + 
												   " must be in range [0,90] (realistically, probably [0-50])");
						  } 
						// convert to radians
						stackdata.theta_specified *= Math.PI/180.0;
						// & set up angles...
						stackdata.computeStackParams();
					  }
				  }
			  }
			
			// is PS data ?
			nv = (Attr)(nnmap.getNamedItem(PS_STACK_KEYWORD_TAG));
			if (nv != null)
			  {
				attr_val = nv.getValue();
				if (attr_val != null && attr_val.equals(TRUE_TAG))
				  {
					stackdata.isPSdata = true;
				  }
			  }

			stackdata.filename = extractString(entry,FILENAME_TAG);

			// override stack filenames if runtime arguments insist
			if (i==0 && _m._SI.cmdLineNearStackFilename!=null) stackdata.filename = _m._SI.cmdLineNearStackFilename;
			if (i==1 && _m._SI.cmdLineFarStackFilename!=null) stackdata.filename = _m._SI.cmdLineFarStackFilename;

			stackdata.setStackRange(extractNumber(entry,FAR_OFFSET_TAG,false,true), 
									extractNumber(entry,NEAR_OFFSET_TAG,false,true));
			// reflector time comes in in ms ...Convert to seconds in stackdata
			stackdata.setTime(0.001*extractNumber(entry,REFLECTOR_TIME_TAG,true,true));
			stackdata.setVstack(extractNumber(entry,STACK_VELOCITY_TAG,true,true));

			// obtain any time_above/below entries. Inputs in ms...Convert to seconds in stackdata
			// (missing entries passed over quietly - wavelet defaults apply)
			stackdata.setTabove(0.001*extractNumber(entry,WAVELET_T_ABOVE_TAG,true,false));
			stackdata.setTbelow(0.001*extractNumber(entry,WAVELET_T_BELOW_TAG,true,false));

			// get wavelet entries
			nlist_wavelet = entry.getElementsByTagName(WAVELET_BLOCK_TAG);
			if (nlist_wavelet.getLength()==0) throw new SAXException("No " + WAVELET_BLOCK_TAG + " entry in " +xml_filename);
			entry = ((Element)nlist_wavelet.item(0));
			stackdata.wavelet_filename = extractString(entry,WAVELET_FILENAME_TAG);

			// override wavelet filenames if runtime arguments insist
			if (i==0 && _m._SI.cmdLineNearWaveletFilename!=null) 
			  stackdata.wavelet_filename = _m._SI.cmdLineNearWaveletFilename;
			if (i==1 && _m._SI.cmdLineFarWaveletFilename!=null) 
			  stackdata.wavelet_filename = _m._SI.cmdLineFarWaveletFilename;


			String noise_tag;
			if (i==0) 
			  {
				noise_tag = getUsedPDAlternateTag(entry,LEGACY_NEAR_STACK_NOISE_RMS_TAG,NEAR_STACK_NOISE_RMS_TAG);
				stackdata.noise_rms = extractValueAndIndex(entry,noise_tag,null);
				// first stack. call this 'near' for now
				SM.setNearstackData(stackdata);
			  }
			if (i==1)
			  {
				// there's another. Make sure these are ordered correctly
				noise_tag = getUsedPDAlternateTag(entry,LEGACY_FAR_STACK_NOISE_RMS_TAG,FAR_STACK_NOISE_RMS_TAG);
				stackdata.noise_rms = extractValueAndIndex(entry, noise_tag, null);
				if (SM.getNearstackData().getOffsetMeanSquare() > stackdata.getOffsetMeanSquare())
				  {
					// swap since these are backwards
					SM.setFarstackData(SM.getNearstackData());
					SM.setNearstackData(stackdata);
				  }
				else
				  {
					SM.setFarstackData(stackdata);
				  }
			  }
			// set up log noises..
			SM.precomputeLogNoise();

			// Give rough idea of noise & angles.....
			if (verbosity > 0)
			  {
				EH.println("  stack noise amplitude = %3.2g",stackdata.getNoise());
				EH.printf( "  %s",(stackdata.isPSdata) ? "PSdata, " : "");
				if (stackdata.theta_specified >= 0)
				  {
					EH.printf("stack RMS angle (fixed-incidence) = %3.1f degrees%n",
							   Math.toDegrees(stackdata.theta_specified));
				  }
				else
				  {
					EH.printf("stack RMS angle (approx) = %3.1f degrees%n",
						   (180/Math.PI)*Math.sqrt(stackdata.getOffsetMeanSquare()/
												   MathUtil.sqr(stackdata.getTime()*stackdata.getVstack())));
				  }
			  }
		  }		
	  }
	catch (Exception ex) 
	  {
		ex.printStackTrace();
		return false;
	  }
	return true;
  }

  private double extractAVOConstants(Element entry,String tag)
	throws SAXException
  {
	// main property
	double value;
	Node n  = entry.getElementsByTagName(tag).item(0);
	if (n!=null)
	  {
		try
		  {
			value = Double.parseDouble(XmlUtils.getTextNodeValue(n));
			// if (verbosity > 3)   EH.println(tag + " = " + value);
			return value;
		  }
		catch (NumberFormatException nfe)
		  {
			throw new SAXException(tag + " value for AVO entry not a number in " + xml_filename);
		  }
	  }
	else
	  {
		throw new SAXException(tag + " value for AVO entry not found in " + xml_filename);
	  }
  }

  private String extractString(Element entry,String tag)
	throws SAXException
  {
	String value;
	// main property
	Node n  = entry.getElementsByTagName(tag).item(0);
	if (n!=null)
	  {
		value = XmlUtils.getTextNodeValue(n);
		if (value==null)
		  throw new SAXException(tag + " value not found in " + xml_filename);			  
		return value;
	  }
	else
	  {
		throw new SAXException(tag + " value for AVO entry not found in " + xml_filename);
	  }
  }

  /** 
	 Finds number tagged by tag within the entry. Throws error if tag found but
	 number is not intelligible, or mustBePositive is set true and number is -ve.
	 Will operate silently if mustBePresent=true and no tag entry found.
  */
  private double extractNumber(Element entry,String tag, boolean mustBePositive, boolean mustBePresent)
	throws SAXException
  {
	// main property
	double value = Double.NaN;
	Node n  = entry.getElementsByTagName(tag).item(0);
	Node n_sub;
	if (n!=null)
	  {
		try
		  {
			value = Double.parseDouble(XmlUtils.getTextNodeValue(n));
			// if (verbosity > 3)   EH.println(tag + " = " + value);
			if (mustBePositive && value<0)
			  throw new SAXException(tag + " value must be >= 0 in " + xml_filename);

			return value;
		  }
		catch (NumberFormatException nfe)
		  {
			throw new SAXException(tag + " value not a number in " + xml_filename);
		  }
	  }
	else
	  {
		if (mustBePresent)
		  {
			throw new SAXException(tag + " value not found in " + xml_filename);
		  }
	  }
	return value;
  }


  /**
	 Finds number tagged by tag within the entry. Throws error if tag found but
	 number is not intelligible, or mustBePositive is set true and number is -ve.
	 Will operate silently if mustBePresent=true and no tag entry found.
  */
  private void extractStackedModelInfo(Element entry,String tag, boolean mustBePositive, StackedModelInfo.StackedPropertyInfoDetails sPID)
	throws SAXException
  {
	// main property
	int i;
	double value = Double.NaN;
	Node n  = entry.getElementsByTagName(tag).item(0);
	if (n!=null)
	  {
		try
		  {
			String content;
			content = ((XmlUtils.getTextNodeValue(n)).split(" "))[0];
			value = Double.parseDouble(content);
			// if (verbosity > 3)   EH.println(tag + " = " + value);
			if (mustBePositive && value<0)
			  throw new SAXException(tag + " value must be >= 0 in " + xml_filename);


			// get attrs for sub_layer applicability. Strings like "1,3,5"
			NamedNodeMap nnmap = n.getAttributes();
			if (nnmap != null)
			  {
				Attr nv = (Attr)(nnmap.getNamedItem("applicable_layers"));
				if (nv != null)
				  {
					String attr_val = nv.getValue();
					if (attr_val != null)
					  {
						StringTokenizer args = new  StringTokenizer(attr_val,",:");
						if (args!=null)
						  {
							// a sublist - override defaults
							sPID.setAllApplicable(false); // switch all off except the named ones
							while (args.hasMoreTokens())
							  {
								try
								  {
									i = Integer.parseInt(args.nextToken());
								  }
								catch (NumberFormatException e)
								  {
									throw new SAXException("applicable layers "+attr_val+" for " 
														   +tag+" must be comma separated numbers in "+xml_filename);
								  }
								if (i<=0 || i> sPID.isApplicable.length)
								  {
									throw new SAXException("applicable layers "+attr_val+" for " 
														   +tag+" must be in model range in "+xml_filename);
								  }
								sPID.isApplicable[i-1] = true;
							  }
						  }
					  }
				  }
			  }
			// set value
			sPID.sigma = value;
		  }
		catch (NumberFormatException nfe)
		  {
			throw new SAXException(tag + " value not a number in " + xml_filename);
		  }
	  }
	sPID.sigma = value;
	sPID.checkApplicableConstraints();
  }


  /**
	 Parses <mapping>n:m,p,q</mapping type entries>
  */
  private void extractStackedModelMapInfo(Element entry)
	throws SAXException
  {
	int i;
	StringTokenizer args, target_args;
	String content;
	StackedModelInfo.StackedPropertyMappingDetails details;
	IntArrayList l;

	content = XmlUtils.getTextNodeValue(entry);
	// split into orig_layer : other_args
	args = new StringTokenizer(content,":");
	if (args.countTokens() != 2) 
	  {
		throw new SAXException("<mappings> entry must be of form orig_layer:mapped_layer1,mapped_layer_2... in "+xml_filename);
	  }
	else
	  {
		while (args.hasMoreTokens())
		  {
			try
			  {
				i = Integer.parseInt(args.nextToken());
			  }
			catch (NumberFormatException e)
			  {
				throw new SAXException("<mappings> entry must be of form orig_layer:mapped_layer1,mapped_layer_2... in "
									   +xml_filename+": orig_layer entry is not a number.");
			  } 
			details = new StackedModelInfo.StackedPropertyMappingDetails();
			details.base_layer = i;
			_m.stackedModelInfo.layer_maps.add(details);

			content = args.nextToken();
			// split other_args into n1,n2,n3....
			l = new IntArrayList();
			target_args = new StringTokenizer(content,",");
			if (target_args.countTokens() == 0) 
			  {
				EH.warn("No mapped_layers for <mappings> entry in "+xml_filename);
			  }
			else
			  {
				while (target_args.hasMoreTokens())
				  {
					try
					  {
						i = Integer.parseInt(target_args.nextToken());
					  }
					catch (NumberFormatException e)
					  {
					throw new SAXException("<mappings> entry must be of form orig_layer:mapped_layer1,mapped_layer_2... in "
										   +xml_filename+": one of the mapped_layer entries is not a number.");
					  }
					l.add(i);
				  }
			  }
			l.trimToSize();
			details.target_layers = l.elements();
		  }
	  }
  }

  /**
	 Finds integer number tagged by tag within the entry. Throws error if tag found but
	 number is not intelligible, or mustBePositive is set true and number is -ve.
	 Will operate silently if mustBePresent=true and no tag entry found.
  */
  private int extractIntNumber(Element entry,String tag, boolean mustBePositive, boolean mustBePresent)
	throws SAXException
  {
	// main property
	int value = Integer.MAX_VALUE; 
	Node n  = entry.getElementsByTagName(tag).item(0);
	Node n_sub;
	if (n!=null)
	  {
		try
		  {
			value = Integer.parseInt(XmlUtils.getTextNodeValue(n));
			if (mustBePositive && value<0)
			  throw new SAXException(tag + " value must be >= 0 in " + xml_filename);

			return value;
		  }
		catch (NumberFormatException nfe)
		  {
			throw new SAXException(tag + " value not an integer number in " + xml_filename);
		  }
	  }
	else
	  {
		if (mustBePresent)
		  {
			throw new SAXException(tag + " value not found in " + xml_filename);
		  }
	  }
	return value;
  }



  boolean getStackedModelsInfo(Model m) throws SAXException
  {
	int i;

	// skip parsing is we want none of this guff
	if (m._SI.ignoreVerticalStackConstraint) return true;

	boolean windowSpecified=false;

	NodeList nlist = null;
	try 
	  {
		// ok, look for the MODEL_DESCRIPTION_TAG chunk
		try
		  {
			nlist = (NodeList)xPath.evaluate("/inversion/"+MODEL_DESCRIPTION_TAG+
											 "/stacked_models_info",document,XPathConstants.NODESET);
		  }
		catch (Exception e) 
		  {
			_m._SI.EH.println("Parsing exception: "+e);
		  }

		if (nlist.getLength() != 1) return true;
		
		// info is there. Allocate and parse
		m.stackedModelInfo = new StackedModelInfo(m);
		m._SI.doStackedModel=true;

		NodeList lmlist = ((Element)nlist.item(0)).getElementsByTagName("layer_mappings");
		if (lmlist.getLength() != 1) 
		  {
			throw new SAXException("In: " + xml_filename + ", must have one layer_mappings entry in stacked_models_info");
		  }
		else
		  {
			Element elem = ((Element)lmlist.item(0));
			NodeList maplist = elem.getElementsByTagName("mapping");
			if (maplist.getLength() == 0) 
			  {
				throw new SAXException("In: " + xml_filename + ", must have at least one <mapping> entry "+
									   "in stacked_models_info/layer_mappings");
			  }
			else
			  {
				for (i=0;i< maplist.getLength();i++)
				  {
					Element map_elem = ((Element)maplist.item(i));
					extractStackedModelMapInfo(map_elem);
				  }
				m.stackedModelInfo.setUpMapsArray();
			  }
		  }
		lmlist = ((Element)nlist.item(0)).getElementsByTagName("property_difference_std_deviations");
		if (lmlist.getLength() != 1) 
		  {
			throw new SAXException("In: " + xml_filename + ", must have one <property_relative_std_deviations> entry in stacked_models_info"); 
		  }
		else
		  {
			Element elem = ((Element)lmlist.item(0));
			StackedModelInfo mI = m.stackedModelInfo;
			extractStackedModelInfo(elem ,"vp_m",true, mI.vp_m); 
			extractStackedModelInfo(elem ,"vp_s",true,mI.vp_s); 
			extractStackedModelInfo(elem ,"vs_m",true,mI.vs_m); 
			extractStackedModelInfo(elem ,"vs_s",true,mI.vs_s); 
			extractStackedModelInfo(elem ,"rho_m",true,mI.rho_m); 
			extractStackedModelInfo(elem ,"phi_s",true,mI.phi_s); 
			extractStackedModelInfo(elem ,"NG",true,mI.NG); 
			extractStackedModelInfo(elem ,"thickness",true,mI.thickness); 
		  }
		lmlist = ((Element)nlist.item(0)).getElementsByTagName("seismic_likelihood_exclusion_range"); 
		if (lmlist.getLength() > 0)
		  {
			windowSpecified = true;
			m.stackedModelInfo.seismicExclusionLayers = new StackedModelInfo.SeismicExclusionLayers[lmlist.getLength()];
			for (i=0;i< m.stackedModelInfo.seismicExclusionLayers.length; i++)
			  {
				m.stackedModelInfo.seismicExclusionLayers[i] = new StackedModelInfo.SeismicExclusionLayers();
				StackedModelInfo.SeismicExclusionLayers sel = m.stackedModelInfo.seismicExclusionLayers[i];
				Element elem = ((Element)lmlist.item(i));
				sel.seismic_exclusion_start_layer = extractIntNumber(elem ,"start_layer",true,true);
				sel.seismic_exclusion_end_layer = extractIntNumber(elem ,"end_layer",true,true);		  
			  }
		  }
		
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }
	// checks

	for (i=0;i<m.stackedModelInfo.layer_maps.size();i++)	  
	  {
		if ((m.stackedModelInfo.layer_maps.get(i)).base_layer > m.n_layers)
		  {
			throw new SAXException("In: " + xml_filename + ", base_layer too high in stacked_models_info/layer_mappings/mapping"); 
		  }
		if ((m.stackedModelInfo.layer_maps.get(i)).base_layer < 1)
		  {
			throw new SAXException("In: " + xml_filename + ", base_layer too low (<1) in stacked_models_info/layer_mappings/mapping"); 
		  }
		int[] l = m.stackedModelInfo.layer_maps.get(i).target_layers;
		for (int j=0; j<l.length; j++)
		  {
			if (l[j] > m.n_layers)
			  {
				throw new SAXException("In: " + xml_filename + ", target_layer too high in stacked_models_info/layer_mappings/mapping"); 
			  }
			if (l[j] < 1)
			  {
				throw new SAXException("In: " + xml_filename + ", target_layer too low (<1) in stacked_models_info/layer_mappings/mapping"); 
			  }
		  }
	  }

	if (windowSpecified)
	  {
		for (i=0; i< m.stackedModelInfo.seismicExclusionLayers.length; i++)
		  {
			if (m.stackedModelInfo.seismicExclusionLayers[i].seismic_exclusion_start_layer > m.n_layers)
			  {
				throw new SAXException("In: " + xml_filename + ", seismic_exclusion_start_layer > "+ m.n_layers +" in stacked_models_info"); 
			  }
			if (m.stackedModelInfo.seismicExclusionLayers[i].seismic_exclusion_end_layer <1)
			  {
				throw new SAXException("In: " + xml_filename + ", seismic_exclusion_end_layer < 1 in stacked_models_info"); 
			  } 
		  }
	  }

	return true;
  }

  /**
	 Parse fluid units type XML entries. These look like

     <fluid_units_and_contacts>
        <fluid_unit fluid_properties_synchronised="false">
          <top_layer>1</top_layer>
          <bottom_layer>4</bottom_layer>
          <contact>
                <type>OB</type>
                <depth>2000</depth>
                <sigma_depth>10</sigma_depth>
          </contact>
		  <contact>
                <type>GO</type>
                <depth>1800</depth>
                <sigma_depth>10</sigma_depth>
          </contact>
        </fluid_unit>
     </fluid_units_and_contacts>
  */

  boolean getFluidUnitsAndContacts(Model m) throws SAXException
  {
	int i,j;
	String value=null,err=null;
	int n_contacts,n_fluid_units;
	int tot_n_contacts = 0;			
	NodeList flist=null;

	m.fluidUnitsAndContacts = new FluidUnitsAndContacts(m);

	try 
	  {
		// Look for the fluid_units_and_contacts chunk
		try
		  {
			flist = (NodeList)xPath.evaluate("/inversion/fluid_units_and_contacts",document,XPathConstants.NODESET);
		  }
		catch (Exception e) 
		  {
			_m._SI.EH.println("Parsing exception: "+e);
		  }
		if (flist.getLength()>1) 
		  throw new SAXException("Too many fluid_units_and_contacts entries in: " + this.xml_filename);
		if (flist.getLength()==0) return true;

		// info is there. Allocate and parse
		m._SI.doSimpleFluidUnitsAndContacts = true;

		NodeList fu_list = ((Element)flist.item(0)).getElementsByTagName("fluid_unit");
		n_fluid_units = fu_list.getLength();
		if (n_fluid_units == 0) 
		  {
			throw new SAXException("In: " + xml_filename + ", must have at least one fluid_units_and_contacts/fluid_unit entry");
		  }
		else
		  {
			m.fluidUnitsAndContacts.fluidUnits = new FluidUnitsAndContacts.FluidUnit[fu_list.getLength()];
			for (i=0; i< fu_list.getLength(); i++)
			  {
				m.fluidUnitsAndContacts.fluidUnits[i] = new FluidUnitsAndContacts.FluidUnit();
				FluidUnitsAndContacts.FluidUnit fu = 	m.fluidUnitsAndContacts.fluidUnits[i];

				Element elem = ((Element)fu_list.item(i));
				// Parse for synchronised fluids attr
				NamedNodeMap nnmap = elem.getAttributes();
				if (nnmap != null)
				  {
					Attr nv = (Attr)(nnmap.getNamedItem("fluid_properties_synchronised"));
					if (nv != null)
					  {
						String attr_val = nv.getValue();
						if (attr_val != null && attr_val.equals(TRUE_TAG))
						  {
							fu.fluids_synched = true;
						  }
					  }
				  }
				// Extract ingredients
				NodeList numlist = elem.getElementsByTagName("top_layer");
				if (numlist.getLength() != 1 ) 
				  {
					throw new SAXException("In: " + xml_filename + ", must have 1 fluid_units_and_contacts/fluid_unit/top_layer entry");
				  }
				else
				  {
					fu.firstLayerIndx = extractIntNumber(elem ,"top_layer",true,true);
					if (fu.firstLayerIndx <1 || fu.firstLayerIndx >_m.n_layers)
					  {
						throw new SAXException("In: " + xml_filename + ", fluid_units_and_contacts/fluid_unit/top_layer" 
											   +" number must be in [1,"+_m.n_layers+"]");
					  }
				  }
				numlist = elem.getElementsByTagName("bottom_layer");
				if (numlist.getLength() != 1 ) 
				  {
					throw new SAXException("In: " + xml_filename + ", must have 1 fluid_units_and_contacts/fluid_unit/bottom_layer entry");
				  }
				else
				  {
					fu.finalLayerIndx = extractIntNumber(elem ,"bottom_layer",true,true);
					if (fu.finalLayerIndx <1 || fu.finalLayerIndx >_m.n_layers)
					  {
						throw new SAXException("In: " + xml_filename + ", fluid_units_and_contacts/fluid_unit/bottom_layer" 
											   +" number must be in [1,"+_m.n_layers+"]");
					  }
					if (fu.finalLayerIndx < fu.firstLayerIndx)
					  {
						throw new SAXException("In: " + xml_filename + ", fluid_units_and_contacts/fluid_unit/bottom_layer" 
											   +" number must be greater than top_layer");
					  }
				  }
				// Get list of internal contacts
				NodeList c_list = elem.getElementsByTagName("contact");
				if (c_list.getLength() == 0) 
				  {
					throw new SAXException("In: " + xml_filename + ", must have at least one <contact> entry for each fluid_unit");
				  }
				else
				  {
					fu.contactInfo = new FluidUnitsAndContacts.ContactInfo[c_list.getLength()];
					n_contacts = c_list.getLength();
					tot_n_contacts += n_contacts;

					for (j=0; j< n_contacts; j++)
					  {
						fu.contactInfo[j] = new FluidUnitsAndContacts.ContactInfo();
						FluidUnitsAndContacts.ContactInfo cI = 	fu.contactInfo[j];

						Element c_elem = ((Element)c_list.item(j));
						// type
						NodeList ctlist = c_elem.getElementsByTagName("type");
						if (ctlist.getLength() != 1 ) 
						  {
							throw new SAXException("In: " + xml_filename + ", must have 1 <type> entry in every fluid_unit");
						  }
						else
						  {
							value = XmlUtils.getTextNodeValue(ctlist.item(0));
							if (value == null || (err=cI.setContactTypeFromXMLKey(value))!=null)				  
							  {
								throw new SAXException("In: " + xml_filename + ", must have sensible " 
													   + " fluid_units_and_contacts/fluid_unit/contact/type: "+err);
							  }
						  }

						// depth
						cI.depthPD = extractMandatoryValueAndIndex(c_elem,"depth"," (fluid units and contacts section)");

						// sigma-depth

						cI.sigma_depthPD = extractValueAndIndex(c_elem,"sigma_depth"," (fluid units and contacts section)");
						if (cI.sigma_depthPD==null) cI.sigma_depthPD = new ModelPriorDescriptions.PriorDetails((double)0.0);


						// Create the actual properties. these should be set to the means of the
						// distributions defined above
						cI.depth = new LayerProperties.Property(cI.depthPD.value,false);

						// add depth to the map of property names too add _fu_contact number suffixes
						_m.completePropertyNameMap.put(cI.depth,"contact_depth_"+(i+1)+"_"+(j+1)); 
					  } 
				  }
			  }
		  }
		// ensure there are not too many contacts
		if (tot_n_contacts > _m.n_layers)
		  {
			throw new SAXException("In: " + xml_filename + ", there are a total of "+tot_n_contacts+" fluid contacts."
								   +" There must be no more contacts than the total number of model layers");
		  }
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }


	return true;
  }
  
  boolean checkFluidUnitsAndContacts(Model m)
  {
	// checks
	if (m._SI.doSimpleFluidUnitsAndContacts)
	  {	
		if (!m.fluidUnitsAndContacts.checkConfigurationIsPossible())
		  {
			return false;
		  }
	  }
	// allocate required extra fields in LayerProperties...
	m.setUpFluidContactInfoInLayers();
	
	return true;
  }
  
  /**
	Small utility data class for parsing tags for physical properties (NG, v_p etc)
  */
  public final class propertyTagInfo
  {
	/**
	   Normalised properties must occupy the interval [0,1]
	*/
	boolean normalised =false;
	/**
	   property name
	*/
	String name=null;
	/**
	   tag to find this property
	*/
	String tag;
	
	public propertyTagInfo(String name, String tag, boolean normalised)
	{
	  this.name = name;
	  this.tag = tag;
	  this.normalised = normalised;
	}
  }
}
