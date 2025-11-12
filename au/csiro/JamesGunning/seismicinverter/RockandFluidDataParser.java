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
import au.csiro.JamesGunning.lib.util.MathUtil;


import javax.xml.parsers.*;
import javax.xml.xpath.*;
import org.xml.sax.*;
import org.w3c.dom.*; 
import java.io.*;
import java.util.*;

public final class RockandFluidDataParser 
{
  /** 
	  Error Handler for instance of this class
  */
  ErrorHandler EH;


  /**
	 level of noise to make while parsing
  */
  int verbosity=0;

  /**
	 file to open and parse
  */
  String filename;

  /** 
	  underlying model reference
   */
  Model _m;

 /**
	 If this is true, input density regressions are power laws of the 
	 Gardner-Gardner-Gregory type: rho = A vp^B
  */
  public static final boolean ASSUME_DENSITY_REGRESSIONS_ARE_POWER_LAWS = true;

  /**
	 delimits set of rock & fluid properties
  */
  public static final String ROCK_FLUID_PROPERTIES_TAG="rock_fluid_properties";

  /**
	 delimits set of fluid properties
  */
  public static final String FLUID_PROPERTIES_TAG="fluid_properties";

  /**
	 delimits reference fluid properties
  */
  public static final String REFERENCE_FLUID_TAG="reference_fluid";

  /**
	 delimits a fluid properties block
  */
  public static final String FLUID_TAG="fluid";

  /**
	 delimits single p-wave velocity entry
  */
  public static final String PWAVE_TAG="vp";
  /**
	 delimits single p-wave velocity entry
  */
  public static final String SIGMA_PWAVE_TAG="sigma_vp";

  /**
	 delimits single s-wave velocity entry
  */
  public static final String SWAVE_TAG="vs";
  
  /**
	 delimits single density entry
  */
  public static final String DENSITY_TAG="rho";
  /**
	 delimits single density entry
  */
  public static final String SIGMA_DENSITY_TAG="sigma_rho";

  /**
	 delimits single grain p-wave velocity entry
  */
  public static final String GRAIN_PWAVE_TAG="vp_grain";

  /**
	 delimits single grain s-wave velocity entry
  */
  public static final String GRAIN_SWAVE_TAG="vs_grain";
/**
	 delimits single grain density entry
  */
  public static final String GRAIN_DENSITY_TAG="rho_grain";

  /**
	 delimits set of rock properties, for all facies
  */
  public static final String ROCK_PROPERTIES_TAG="rock_properties";
  /**
	 delimits tage of rock name
  */
  public static final String NAME_TAG="name";
  /**
	 delimits set of rock curves, for reservoir facies
  */
  public static final String RESERVOIR_FACIES_TAG="reservoir_endmember";
  /**
	 delimits set of rock curves, for non-reservoir facies
  */
  public static final String NONRESERVOIR_FACIES_TAG="nonreservoir_endmember";
  /**
	 reservoir prefix in variable trends tags
  */
  public static final String RESERVOIR_PREFIX="reservoir_";
  /**
	 nonreservoir prefix in variable trends tags
  */
  public static final String NONRESERVOIR_PREFIX="nonreservoir_";

  /**
	 delimits tag for pwave curve in XML file
  */
  public static final String PWAVE_CURVE_TAG="vp_curve";
  /**
	 delimits tag for swave curve in XML file
  */
  public static final String SWAVE_CURVE_TAG="vs_curve";
  /**
	 delimits tag for density curve (shales) in XML file
  */
  public static final String DENSITY_CURVE_TAG="density_curve";
  /**
	 delimits tag for prosity curve (sands) in XML file
  */
  public static final String POROSITY_CURVE_TAG="porosity_curve";
  /**
	 delimits intercept in a regression curve
  */
  public static final String INTERCEPT_TAG="intercept";
  /**
	 delimits leading coefficient in a regression curve
  */
  public static final String FACTOR_TAG="factor";
  /**
	 delimits intercept in a regression curve
  */
  public static final String SLOPE_TAG="slope";
  /**
	 delimits power exponent in a regression curve
  */
  public static final String EXPONENT_TAG="exponent";
  /**
   delimits depth coefficient in vp regression
  */
  public static final String DEPTH_COEFF_TAG="depth_coefficient";
  /**
   delimits LFIV coefficient in vp regression
  */
  public static final String LFIV_COEFF_TAG="LFIV_coefficient";

  // ENABLE_CAPTURED_FLOATING_GRAIN_MODEL...
		 
  /**
   delimits float_fraction coefficient in vp regression
  */
  public static final String FLOAT_FRACTION_COEFF_TAG="phi_float_coefficient";
  /**
	 string with which to construct link tag for  float_fraction coefficient in vp regressions
  */
  public static final String VP_FLOAT_FRACTION_COEFF="vp_phi_float_coefficient";
  /**
	 string with which to construct link tag for  float_fraction coefficient in phi regressions
  */
  public static final String PHI_FLOAT_FRACTION_COEFF="phi_phi_float_coefficient";

  // ....ENABLE_CAPTURED_FLOATING_GRAIN_MODEL

  /**
	 delimits  density regression error
  */
  public static final String SIGMA_TAG="sigma";

  /**
	 string with which to construct link tag for  intercept in a vp regression curve
  */
  public static final String VP_INTERCEPT="vp_intercept";
  /**
	 string with which to construct link tag for  intercept in a vs regression curve
  */
  public static final String VS_INTERCEPT="vs_intercept";
  /**
	 string with which to construct link tag for  intercept in a vp regression curve
  */
  public static final String PHI_INTERCEPT="phi_intercept";
/**
	 string with which to construct link tag for sigma in a vp regression curve
  */
  public static final String VP_SIGMA="vp_sigma";
  /**
	 string with which to construct link tag for  sigma in a vs regression curve
  */
  public static final String VS_SIGMA="vs_sigma";
  /**
	 string with which to construct link tag for  sigma in a phi regression curve
  */
  public static final String PHI_SIGMA="phi_sigma";
 /**
	 string with which to construct link tag for  sigma in a rho regression curve
  */
  public static final String RHO_SIGMA="rho_sigma";



  /**
	 string with which to construct link tag for  leading coefficient in a regression curve
  */
  public static final String FACTOR="factor";
  /**
	 string with which to construct link tag for  intercept in a regression curve
  */
  public static final String SLOPE="slope";
  /**
	 string with which to construct link tag for  power exponent in a regression curve
  */
  public static final String EXPONENT="exponent";
  /**
   string with which to construct link tag for  depth coefficient in vp regression
  */
  public static final String DEPTH_COEFF="vp_depth_coefficient";

  /**
   string with which to construct link tag for  LFIV coefficient in vp regression
  */
  public static final String LFIV_COEFF="vp_LFIV_coefficient";
  /**
   string with which to construct link tag for vp coefficient in vs regression
  */
  public static final String VS_VP_COEFF="vs_vp_coefficient";
  /**
   string with which to construct link tag for  phi coefficient in vp regression
  */
  public static final String PHI_VP_COEFF="phi_vp_coefficient";




  /**
	 string with which to construct link tag for   density regression error
  */
  public static final String SIGMA="sigma";

  // XML parsing for xmlPaths
  XPath xPath = XPathFactory.newInstance().newXPath();

  Document document;

  /**
	 Make a parser with this level of chatter (0 = silent)
  */
  public RockandFluidDataParser(String filename,int verbosity,ErrorHandler EH)
  {
	this.filename = filename;
	this.verbosity = verbosity;
	this.EH = EH;

  }

  /** 
	  Function to read all rock and fluid data into a general Model class
	  The model_filename and rocks_filename may be the same, but we'll assume 
	  they can be different
  */


  protected boolean getRockandFluidData(Model m)
  {
	// set up internal reference to model
	_m = m;

	// First, parsing is dependent on having layer meta-data in place is this true?

	// Create HashMap of end_member rocks, indexed by rock name as per the layers xml file.
	m.RFD.RockTrends = new HashMap();

	// xml file parser
	DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();

	try 
	  {
		DocumentBuilder builder = factory.newDocumentBuilder();
		this.document = builder.parse( new File(filename) );
		
		// stuff data into RockandFluidData
		
		// (1) Extract rock curves: loop over facies names etc 
		extractRockProperties(_m);
		
		// (2) Extract fluid properties
		extractFluidProperties(_m);
	  } 
	catch (Exception ex) 
	  {
		ex.printStackTrace();
		return false;
	  }
	return true;
  }

  private void extractRockProperties(Model m) throws SAXException
  {
	String rock_name = null;
	boolean isRockPermeable;
	int i,j;
	ModelPriorDescriptions.MetaInfo MI;

	RockandFluidData.MaterialProperties grain_props;
	// loop over all rock types we seem to need from the meta-Layer data
	m.RFD.GrainTypes = new HashMap();

	for (i=0;i<m.MPD.MetaInfo.length;i++)
	  // loop over layers
	  {
		MI = m.MPD.MetaInfo[i];

		ROCK_TYPES:
		for (j=0 ; j<=6 ; j++)
		  {
			// do sand-end-members, then mud, then hard_rock members
			switch (j)
			  {
			  case 0:
				rock_name = MI.sand_rock_member ; break;
			  case 1:
				rock_name = MI.mud_rock_member ; break;
			  case 2: 
				rock_name = MI.hard_shale_rock_member ; break;
			  case 3:
			  case 4:
			  case 5:
			  case 6:
				if (MI.fluid_dep_hard_shale_rock_members == null) 
				  {
					break ROCK_TYPES;
				  }
				else
				  {
					rock_name = MI.fluid_dep_hard_shale_rock_members[j-3]; // maps to indexes 0,1,2,3
				  }
			  default:
			  }

			isRockPermeable = (j==0);
			// All rocks with non--null names need a rock-curve
			// A rock may have a null name if it's the sand end-member in a pure shale
			// and therefore contributes nothing to the physics
			if (rock_name != null)
			  {				
				if (m.RFD.RockTrends.get(rock_name) == null)
				  {
					if (verbosity>0) 
					  {
						EH.println("%nSeeking trend curves for: " + rock_name);
					  }

					// Creat new Regression set for the hashmap database
					m.RFD.curveSet = new RockandFluidData.RegressionSet();

					// must have grain properties of sands
					if (isRockPermeable) 
					  {
						extractGrainPropertiesAndRefFluidName(rock_name,m,i);
						// Each permeable layer will  need to see its grain properties
						grain_props = (RockandFluidData.MaterialProperties) m.RFD.GrainTypes.get(rock_name);

						m.RFD.curveSet.rho_g = grain_props.rho;
						m.RFD.curveSet.K_g = m.RFD.curveSet.rho_g * (MathUtil.sqr(grain_props.vp) 
															  - (4.0/3.0) * MathUtil.sqr(grain_props.vs));
						// Layer reference to suite of regressions. This carries
						// ref fluid props and grain props
						m.LP[i].RFDregSet = m.RFD.curveSet;

						if (m.RFD.curveSet.K_g < 0)
						  throw new SAXException("Grain properties lead to -ve K_g for " 
												 + rock_name + " in file " + filename);
					  }
						

					// pwave curve

					if (getCurveData(rock_name,PWAVE_CURVE_TAG,m.RFD.curveSet.Vp_versus_depth,false,isRockPermeable))
					  {
						if (verbosity>0) EH.println("      Found : " + rock_name + ", " + PWAVE_CURVE_TAG);
					  }
					else
					  {
						throw new SAXException("Couldn't get "+ PWAVE_CURVE_TAG + 
											   " for " + rock_name + " in file " + filename);
					  }

					// swave curve

					if (getCurveData(rock_name,SWAVE_CURVE_TAG,m.RFD.curveSet.Vs_versus_Vp,false,isRockPermeable))
					  {
						if (verbosity>0) EH.println("      Found : " + rock_name + ", " + SWAVE_CURVE_TAG);
					  }
					else
					  {
						throw new SAXException("Couldn't get "+ SWAVE_CURVE_TAG + 
											   " for " + rock_name + " in file " + filename);
					  }


					// now for sands we need porosity and ref fluid
					if (isRockPermeable)
					  {
						if (getCurveData(rock_name,POROSITY_CURVE_TAG,m.RFD.curveSet.Rho_or_Phi_versus_Vp, false,isRockPermeable))
						 {
						   if (verbosity>0) EH.println("      Found : " + rock_name + ", " + POROSITY_CURVE_TAG);
						 }
						else
						  {
							throw new SAXException("Couldn't get "+ POROSITY_CURVE_TAG +
												   " for " + rock_name + " in file " + filename);
						  }

						// Get reference. fluid
						getCurveRefFluid(rock_name,m.RFD.curveSet);						
					  }
					else
					// for muds we need density
					  {
						if (getCurveData(rock_name,DENSITY_CURVE_TAG,m.RFD.curveSet.Rho_or_Phi_versus_Vp, true, isRockPermeable))
						  {
							if (verbosity>0) EH.println("      Found : " + rock_name + ", " + DENSITY_CURVE_TAG);
						  }
						else
						  {
							throw new SAXException("Couldn't get "+ DENSITY_CURVE_TAG + 
												 " for " + rock_name + " in file " + filename);
						  }
					  }

					// put into hashmap for easy retrieval
					m.RFD.RockTrends.put(rock_name,m.RFD.curveSet);
				  }
				else
				  {
					if (isRockPermeable) 
					  {
						// we have this sand-rock already but must shove a reference to its Regression set into the LayerProperties						
						RockandFluidData.RegressionSet RFDregSet = (RockandFluidData.RegressionSet)m.RFD.RockTrends.get(rock_name);
						m.LP[i].RFDregSet = RFDregSet;
					  }
				  }
			  } // rock_name!=null
		  } //sand, mud, hard?
	  } // loop over layers


	// get any areally-variable links associated with trends
	obtainVariablePriorRockData(m);

	// cement in parsed values, forcing the mapping on this first occasion
	m.RFD.resetTrendData(m,true);


	// final sanity test on defaults from XML
	m.RFD.checkRockSpeedRegressions(m,true,filename);
  }

  /** 
	  extract coefficients of rock physics regressions. If isPowerLaw==true,
	  some tags are modified (INTERCEPT_TAG -> FACTOR_TAG, SLOPE_TAG -> EXPONENT_TAG)
  */

  protected boolean getCurveData(String facies_name, String tagname, 
								 RockandFluidData.Regression rockdata_curve,
								 boolean isPowerLaw, boolean isRockPermeable) throws SAXException
  {
	Node n,nn,entry;
	boolean found = false;
	String tag=null;

	// ok, look for the ROCK_FLUID_PROPERTIES_TAG chunk
	NodeList nlist = null;
	try
	  {
		nlist = (NodeList)xPath.evaluate("/inversion/"+ROCK_FLUID_PROPERTIES_TAG,document,XPathConstants.NODESET);
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }

	if (nlist.getLength()!=1) throw new SAXException("No " + ROCK_FLUID_PROPERTIES_TAG + " entry in this file");

	// ok, look for the ROCK_PROPERTIES_TAG chunk <rock_properties>
	nlist = ((Element) nlist.item(0)).getElementsByTagName(ROCK_PROPERTIES_TAG);
	if (nlist.getLength()!=1) throw new SAXException("less or more than 1" + 
													 ROCK_PROPERTIES_TAG + " entry in this file");
	// there should be only one occurence of <rock_properties_tag>, so the first is what we want
	Element flist = ((Element)nlist.item(0));

	// loop over non_reservoir, then reservoir members
	for (int end_type=0; end_type<=1 ; end_type++)
	  {
		// get the list of nodes with RESERVOIR_FACIES_TAG or NONRESERVOIR_FACIES_TAG, as appropriate
		NodeList facies_list  = flist.getElementsByTagName((end_type==0) ? 
														   NONRESERVOIR_FACIES_TAG : RESERVOIR_FACIES_TAG);
		// loop over facies nodes
		for (int i=0;i<facies_list.getLength();i++)
		  {
			// Does the <name> entry in this <facies> block match the desired name
			nn = ((Element)facies_list.item(i)).getElementsByTagName(NAME_TAG).item(0);
			String name = XmlUtils.getTextNodeValue(nn);

			if (name.equals(facies_name))
			  {
				// Ok, facies_list.item(i) has the subtree for the facies we required within ROCK_PROPERTIES_TAG
				//	find the Node in facies_list.item(i) that is tagged by <tagname>. (eg <vp_curve>)
				//  It'll be the first of these
				if (this.verbosity > 0) 
				  EH.println("    Parsing facies: " + facies_name + ", curve: " + tagname);
			
				n = ((Element)facies_list.item(i)).getElementsByTagName(tagname).item(0);
				if (n!= null)
				  {
					// Ok there is an entry for this rock type, for the nominated regression relation
					// Eg facies_name = "sand A", tagname='vp_curve"

					// Get the regression constants
				
					// (1) All regressions have an intercept
					tag = (isPowerLaw)? FACTOR_TAG : INTERCEPT_TAG;
					entry = ((Element) n).getElementsByTagName(tag).item(0);
					if (entry==null) throw new SAXException("Couldn't get " + tag + " in " + tagname + 
															" in " + facies_name + " in file " + filename);
					// Does it vary?
					setVariabilityFlag(entry,rockdata_curve,rockdata_curve.intercept_PD);
					rockdata_curve.intercept_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));

					// Do each regression on a case by case basis
					if (tagname == PWAVE_CURVE_TAG)
					  {
						// (2) depth-coeff
						entry = ((Element) n).getElementsByTagName(DEPTH_COEFF_TAG).item(0);
						if (entry==null) 
						  throw new SAXException("Couldn't get " + DEPTH_COEFF_TAG+ " in " + tagname + 
												 " in " + facies_name + " in file " + filename);
						setVariabilityFlag(entry,rockdata_curve,rockdata_curve.depth_coeff_PD); // Does it vary???

						rockdata_curve.depth_coeff_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));

						// (3) LFIV-coeff
						entry = ((Element) n).getElementsByTagName(LFIV_COEFF_TAG).item(0);
						if (entry==null) throw new SAXException("Couldn't get " + LFIV_COEFF_TAG+ " in " + tagname + 
																" in " + facies_name + " in file " + filename);
						rockdata_curve.LFIV_coeff_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
						setVariabilityFlag(entry,rockdata_curve,rockdata_curve.LFIV_coeff_PD); // Does it vary???
						
						if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
						  {
							if (isRockPermeable)
							  {
								// (3) phi_float-coeff
								entry = ((Element) n).getElementsByTagName(FLOAT_FRACTION_COEFF_TAG).item(0);
								if (entry==null) 
								  {
									// Accept as benign, with warning. Set coefficient to 0, fixed.
									EH.println("      Warning: no " + FLOAT_FRACTION_COEFF_TAG+ " in " + tagname + 
											   " in " + facies_name + " in file " + filename + ". Assuming 0.");
									rockdata_curve.float_fraction_coeff_PD.value = 0.0;
								  }
								else
								  {
									rockdata_curve.float_fraction_coeff_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
									setVariabilityFlag(entry,rockdata_curve,rockdata_curve.float_fraction_coeff_PD); // Does it vary???							
								  }
							  }
						  }		
					  }
					else if ((tagname == DENSITY_CURVE_TAG) || (tagname == POROSITY_CURVE_TAG)) 
					  {
						// Regression is against vp and maybe float_fraction

						// (2) slope, or exponent
						tag = (isPowerLaw)? EXPONENT_TAG : SLOPE_TAG;
						entry = ((Element) n).getElementsByTagName(tag).item(0);
						if (entry==null) throw new SAXException("Couldn't get " + tag + " in " + tagname + 
																" in " + facies_name + " in file " + filename);
						rockdata_curve.slope_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
						setVariabilityFlag(entry,rockdata_curve,rockdata_curve.slope_PD); // Does it vary???

						if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
						  {
							// extra term for porosity trend of float dependence
							if (tagname == POROSITY_CURVE_TAG && isRockPermeable)
							  {
								// (3) phi_float-coeff
								entry = ((Element) n).getElementsByTagName(FLOAT_FRACTION_COEFF_TAG).item(0);
								if (entry==null)
								  {
									EH.println("      Warning: no " + FLOAT_FRACTION_COEFF_TAG+ " in " + tagname + 
											   " in " + facies_name + " in file " + filename + ". Assuming 0.");
									rockdata_curve.float_fraction_coeff_PD.value = 0.0;
								  }
								else
								  {
									rockdata_curve.float_fraction_coeff_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
									setVariabilityFlag(entry,rockdata_curve,rockdata_curve.float_fraction_coeff_PD); // Does it vary???
								  }
							  }
						  }

					  }
					else if (tagname == SWAVE_CURVE_TAG) // Regression is only against vp
					  {
						// (2) slope, or exponent
						tag = (isPowerLaw)? EXPONENT_TAG : SLOPE_TAG;
						entry = ((Element) n).getElementsByTagName(tag).item(0);
						if (entry==null) throw new SAXException("Couldn't get " + tag + " in " + tagname + 
																" in " + facies_name + " in file " + filename);
						rockdata_curve.slope_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
						setVariabilityFlag(entry,rockdata_curve,rockdata_curve.slope_PD); // Does it vary???

					  }
					else
					  {
						throw new SAXException("Couldn't understand " + tagname + 
											   " entry for rock " + facies_name + " in file " + filename);
					  }
				
					// (3) error
					entry = ((Element) n).getElementsByTagName(SIGMA_TAG).item(0);
					if (entry==null) throw new SAXException("Couldn't get " + SIGMA_TAG+ " in " + tagname + 
															" in " + facies_name + " in file " + filename);
					rockdata_curve.error_PD.value = Float.parseFloat(XmlUtils.getTextNodeValue(entry));
					setVariabilityFlag(entry,rockdata_curve,rockdata_curve.error_PD); // Does it vary???

					// for now have to have non-negative errors (redundant cases not acceptable)
					if (rockdata_curve.error_PD.value<0)
					  throw new SAXException("Regression error must be +ve: " + SIGMA_TAG+ " in " + tagname + 
											 " in " + facies_name + " in file " + filename);

					// there should only be one rock-match. Flag success and break out
					found=true;
					break;
				  }
				else
				  {
					throw new SAXException("Can't find a " + tagname + 
										   " entry for rock-type " + facies_name + " in " + filename);
				  }
			  }
		  }
	  }
	return found;
  }

  
  /** 
	  Extract reference fluid for rock with name facies_name, stuff it into
	  RFDregSet
  */

  void getCurveRefFluid(String facies_name, RockandFluidData.RegressionSet RFDregSet) throws SAXException
  {
	Node n,nn,entry;
	NodeList nlist=null;

	// ok, look for the ROCK_FLUID_PROPERTIES_TAG chunk
	try
	  {
		nlist = (NodeList)xPath.evaluate("/inversion/"+ROCK_FLUID_PROPERTIES_TAG,document,XPathConstants.NODESET);
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }

	if (nlist.getLength()!=1) throw new SAXException("No " + ROCK_FLUID_PROPERTIES_TAG + " entry in this file");

	// ok, look for the ROCK_PROPERTIES_TAG chunk <rock_properties>
	nlist = ((Element) nlist.item(0)).getElementsByTagName(ROCK_PROPERTIES_TAG);
	if (nlist.getLength()!=1) throw new SAXException("less or more than 1" + 
													 ROCK_PROPERTIES_TAG + " entry in this file");
	// there should be only one occurence of <rock_properties_tag>, so the first is what we want
	Element flist = ((Element)nlist.item(0));

	// get the list of nodes with RESERVOIR_FACIES_TAG
	NodeList facies_list  = flist.getElementsByTagName(RESERVOIR_FACIES_TAG);
	// loop over facies nodes
	for (int i=0;i<facies_list.getLength();i++)
	  {
		// Does the <name> entry in this <facies> block match the desired name
		nn = ((Element)facies_list.item(i)).getElementsByTagName(NAME_TAG).item(0);
		String name = XmlUtils.getTextNodeValue(nn);
		
		if (name.equals(facies_name))
		  {
			// Ok, facies_list.item(i) has the subtree for the facies we required within ROCK_PROPERTIES_TAG
			//	find the Node in facies_list.item(i) that is tagged by REFERENCE_FLUID_TAG
			//  It'll be the first of these
	
			n = ((Element)facies_list.item(i)).getElementsByTagName(REFERENCE_FLUID_TAG).item(0);
			if (n!= null)
			  {
				RFDregSet.ref_brine_name = XmlUtils.getTextNodeValue(n);
			  }
			else
			  {
				throw new SAXException("Can't find a reference fluid entry for rock-type " + facies_name + " in " + filename);
			  }
			break;
		  }
	  }
  }



  private void setVariabilityFlag(Node entry, RockandFluidData.Regression curve, ModelPriorDescriptions.PriorDetails pd)
  {
	NamedNodeMap nnmap;
	Attr nv;
	String attr_val;


	nnmap = entry.getAttributes();
	if (nnmap != null)
	  {
		nv = (Attr)(nnmap.getNamedItem(ModelDescriptionParser.VARIES_TAG));
		if (nv != null)
		  {
			attr_val = nv.getValue();
			if (attr_val != null && attr_val.equals(ModelDescriptionParser.TRUE_TAG))
			  {
				curve.regressionsVary = true;   // so later on we know if a rock has variable regressions
				pd.varies = true;
			  }
		  }
	  }
  }

  protected void extractGrainPropertiesAndRefFluidName(String facies_name,Model m, int layerIndx) throws SAXException
  {
	Node n,entry;
	boolean found = false;

	// ok, look for the ROCK_PROPERTIES_TAG chunk <rock_properties>
	NodeList nlist = null; 
	try
	  {
		nlist = (NodeList)xPath.evaluate("/inversion/"+ROCK_FLUID_PROPERTIES_TAG+"/"+
										 ROCK_PROPERTIES_TAG,document,XPathConstants.NODESET);
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }

	if (nlist.getLength()!=1) throw new SAXException("less or more than 1" + 
													 ROCK_PROPERTIES_TAG + " entry in this file");
	// there should be only one occurence of <rock_properties_tag>, so the first is what we want
	Element flist = ((Element)nlist.item(0));

	// get the list of nodes with RESERVOIR_FACIES_TAG <reservoir_endmember>
	NodeList facies_list  = flist.getElementsByTagName(RESERVOIR_FACIES_TAG);
	// loop over facies nodes
	for (int i=0;i<facies_list.getLength();i++)
	  {
		entry = ((Element)facies_list.item(i));
		// Does the <name> entry in this <facies> block match the desired name
		String name = XmlUtils.getTextNodeValue(((Element)entry).getElementsByTagName(NAME_TAG).item(0));

		if (name.equals(facies_name))
		  {
			// Ok, facies_list.item(i) has the subtree for the facies we required within ROCK_PROPERTIES_TAG
			//  It'll be the first of these
			if (this.verbosity > 0) EH.println("    Parsing facies: " + facies_name + ", grain properties");

			// Ok there is an entry for this rock type, 
			// From hereon, we don't already have these grain types and ref_flui name in our hashmap database
			// make a new one
			m.RFD.material_props = new RockandFluidData.MaterialProperties();	


			// (1) look for vp
			n = ((Element) entry).getElementsByTagName(GRAIN_PWAVE_TAG).item(0);
			if (n!= null)
			  {
				try 
				  {
					m.RFD.material_props.vp = Float.parseFloat(XmlUtils.getTextNodeValue(n));
					if (verbosity > 1)   EH.println("      vp = " + m.RFD.material_props.vp);
				  }
				catch (NumberFormatException nfe)
				  {
					throw new SAXException("Property vp for grains of " + facies_name + " not a number: " + filename);
				  }
			  }
							
			// (1) look for vs
			n = ((Element) entry).getElementsByTagName(GRAIN_SWAVE_TAG).item(0);
			if (n!= null)
			  {
				try 
				  {
					m.RFD.material_props.vs = Float.parseFloat(XmlUtils.getTextNodeValue(n));
					if (verbosity > 1)   EH.println("      vs = " + m.RFD.material_props.vs);
				  }
				catch (NumberFormatException nfe)
				  {
					throw new SAXException("Property vs for grains of " + facies_name + " not a number: " + filename);
				  }
			  }
							
			// (1) look for rho
			n = ((Element) entry).getElementsByTagName(GRAIN_DENSITY_TAG).item(0);
			if (n!= null)
			  {
				try 
				  {
					m.RFD.material_props.rho = Float.parseFloat(XmlUtils.getTextNodeValue(n));
					if (verbosity > 1)   EH.println("      rho = " + m.RFD.material_props.rho);
				  }
				catch (NumberFormatException nfe)
				  {
					throw new SAXException("Property density for grains of " + facies_name + " not a number: " + filename);
				  }
			  }
			// Add to hashmap
			m.RFD.material_props.name = facies_name;
			m.RFD.GrainTypes.put(facies_name,m.RFD.material_props);
			found=true;
			break; // we've got what we want. 
		  }
	  }
	if (!found)
	  {
		throw new SAXException("Can't find grain moduli for " + facies_name + " in database ");			
	  }
  }





  private void extractFluidProperties(Model m) throws SAXException
  {
	NodeList sublist;
	Node n;
	String value,tag=null,fluid_name=null,name=null;
	boolean found_desired_entry=false;
	RockandFluidData.RegressionSet RFDregSet=null;
	NodeList nlist = null;

	//  Look for /inversion/rock_fluid_properties/fluid_properties
	try
	  {
		nlist = (NodeList)xPath.evaluate("/inversion/"+ROCK_FLUID_PROPERTIES_TAG
										 +"/"+FLUID_PROPERTIES_TAG,document,XPathConstants.NODESET);
	  }
	catch (Exception e) 
	  {
		_m._SI.EH.println("Parsing exception: "+e);
	  }
	if (nlist.getLength()!=1) 
	  throw new SAXException("No /inversion/rock_fluid_properties/fluid_properties entry in this file");


	// there should be only one occurence of <FLUID_PROPERTIES_TAG>, so the first is what we want
	Element flist = ((Element)nlist.item(0));

	// hash map of required fluids
	m.RFD.FluidTypes = new HashMap();
	  

	// loop over all layers
	if (verbosity > 0)   EH.println("");
	for (int i=0;i<m.MPD.MetaInfo.length;i++)
	  {
		// don't need any fluid properties for pure shales
		if (!m.MPD.MetaInfo[i].is_permeable) continue;

		// but in mixtures, maybe.
		// loop over required fluids in this layer
		// o = ref fluid, 1,2,3,4 =b,o,g,lsg
		for (int fluid_cntr=0;fluid_cntr<=4;fluid_cntr++)
		  {
			switch (fluid_cntr)
			  {
			  case 0:  
				// Reference brine. Pull its name out of RockandFluidData.RegressionSet object in HashMap 
				// indexed by reservoir rock name
				RFDregSet = (RockandFluidData.RegressionSet)
				  m.RFD.RockTrends.get(m.MPD.MetaInfo[i].sand_rock_member);
				if (RFDregSet == null)
				  {
					throw new SAXException("No reservoir rock regression set for rock "+
										   m.MPD.MetaInfo[i].sand_rock_member);
				  }
				else
				  {
					name = RFDregSet.ref_brine_name;
				  }					
				break;
			  case 1: 
				// We always need the new brine, since we need density 
				name = m.MPD.MetaInfo[i].brine_name;					
				break;
				
			  case 2:
				if (m.MPD.MetaInfo[i].oil_substitution_possible)
				  {
					name = m.MPD.MetaInfo[i].oil_name;					
				  }
				else continue;			
				break;

			  case 3:
				if (m.MPD.MetaInfo[i].gas_substitution_possible)
				  {
					name = m.MPD.MetaInfo[i].gas_name;					
				  }
				else continue;			
				break;

			  case 4:
				if (m.MPD.MetaInfo[i].lsg_substitution_possible)
				  {
					name = m.MPD.MetaInfo[i].lsg_name;					
				  }
				else continue;			
				break;
			  }
			tag=FLUID_TAG;
			
			// OK at this point ,we will look for a fluid with matching name.
			
			m.RFD.material_props = (RockandFluidData.MaterialProperties)m.RFD.FluidTypes.get(name);
			if (m.RFD.material_props != null) 
			  {
				continue; 
			  }
			else
			  {
				// nothing in the data base. Create a new fluid material MaterialProperties object
				if (verbosity > 0)   EH.println("Parsing fluid properties of : " + name);
	
				// get the list of nodes for the named fluid-type
				sublist  = flist.getElementsByTagName(tag);
				// loop over, all named fluids
				for (int j=0;j<sublist.getLength();j++)
				  {
					Element entry = (Element) sublist.item(j);
					found_desired_entry=false;
					if (entry!=null)
					  {
						// ok we might have a set of properties to add..
						// does its name match the fluid we want?
						fluid_name = XmlUtils.getTextNodeValue((Element)entry.getElementsByTagName(NAME_TAG).item(0));
						if (!fluid_name.equals(name)) continue;
						// we must have found it
						found_desired_entry=true;

						m.RFD.material_props = new RockandFluidData.MaterialProperties();	
					

						// (1) Vp - will be a list, where we want the first
						
						extractFluidConstants(m.RFD.material_props, 0 ,entry,fluid_name,PWAVE_TAG,true);
						extractFluidConstants(m.RFD.material_props, 1, entry,fluid_name,
																			  SIGMA_PWAVE_TAG,false);
					
						// (2) Then Density - will be a list, where we want the first //

						extractFluidConstants(m.RFD.material_props, 2, entry,fluid_name,DENSITY_TAG,true);
						extractFluidConstants(m.RFD.material_props, 3, entry,fluid_name,
																			   SIGMA_DENSITY_TAG,false);
						// Form a property
						// m.RFD.material_props.vpProp = new LayerProperties.Property(m.RFD.material_props.vp,false);// 190104Revs
						// m.RFD.material_props.vpProp.setDistr(m.RFD.material_props.vp ,m.RFD.material_props.sigma_vp);// 190104Revs
						// m.RFD.material_props.rhoProp = new LayerProperties.Property(m.RFD.material_props.rho,false);// 190104Revs
						// m.RFD.material_props.rhoProp.setDistr(m.RFD.material_props.rho ,m.RFD.material_props.sigma_rho);// 190104Revs

						// Add to hash map of valid fluid types			
						m.RFD.material_props.name = fluid_name;
						m.RFD.FluidTypes.put(fluid_name,m.RFD.material_props);

						// Add to reference fluid properties in Rock regressions on first encounter
						if (fluid_cntr ==0)
						  {
							RFDregSet.rho_b_ref = m.RFD.material_props.rho;
							RFDregSet.vp_b_ref = m.RFD.material_props.vp;
							RFDregSet.K_b_ref = RFDregSet.rho_b_ref * RFDregSet.vp_b_ref * RFDregSet.vp_b_ref;
						  }
						break;   // no point looping over all other eg "brines": we've got what we need
					  }
					
				  }
				if (!found_desired_entry)
				  {
					throw new SAXException("Properties for fluid " + name + " not found in " + filename);
				  }
			  }
			
			// set layer reference to required MaterialProperties object for this fluid
			/*
			  switch (fluid_cntr)
			  {
			  case 0:  
			  // ref brine is fixed
			  break;
			  case 1: 
			  m.LP[i].brineProps = m.RFD.material_props;				
				break;
				case 2:
				m.LP[i].oilProps = m.RFD.material_props;
				break;
				case 3:
				m.LP[i].gasProps = m.RFD.material_props;
				break;
				case 4:
				m.LP[i].lsgProps = m.RFD.material_props;
				break;
			  }
			*/
		  } // loop over possible fluids in layer i
		// Complete set up of LayerProperties hashmaps
		// m.LP[i].addFluidPropertiesToHashMap();
	  } // loop over layers
  }

  private void extractFluidConstants(RockandFluidData.MaterialProperties mp, int type, Element entry, String fluid_name, String tag, boolean mandatory)
	throws SAXException
  {
	// main property
	double value;
	NamedNodeMap nnmap;
	Attr nv;
	String attr_val;
	Node n  = entry.getElementsByTagName(tag).item(0);   
	boolean varies = false;

	if (n!=null)
	  {
		try 
		  {
			value = Double.parseDouble(XmlUtils.getTextNodeValue(n));
			if (verbosity > 1)   EH.println("    " + tag + " = " + value);
		  }
		catch (NumberFormatException nfe)
		  {
			throw new SAXException(tag + " value for fluid " + fluid_name + " not a number " + filename);
		  }
		nnmap = n.getAttributes();
		if (nnmap != null)
		  {
			nv = (Attr)(nnmap.getNamedItem(ModelDescriptionParser.VARIES_TAG));
			if (nv != null)
			  {
				attr_val = nv.getValue();
				if (attr_val != null && attr_val.equals(ModelDescriptionParser.TRUE_TAG)) varies = true;
			  }
		  }
	  }
	else
	  {
		if (mandatory) 
		  throw new SAXException(tag + " value for fluid " + fluid_name + " not found in " + filename);
		value=0.0;
	  }
	// Set in-situ fields
	switch (type)
	  {
	  case 0:
		mp.vp = value;
		mp.vp_varies = varies;
		break;
	  case 1:
		mp.sigma_vp = value;
		mp.sigma_vp_varies = varies;
		break;
	  case 2:
		mp.rho = value;
		mp.rho_varies = varies;
		break;
	  case 3:
		mp.sigma_rho = value;
		mp.sigma_rho_varies = varies;
		break;
	  }
	return;
  }

  /**
	 Obtain link indices for variable rock trends if any
  */
  private void obtainVariablePriorRockData(Model m)
  {
	String rock_name;
	int i;
	ModelPriorDescriptions.MetaInfo MI;
	ArrayList rockList;
	RockandFluidData.RegressionSet curveSet;
	RockandFluidData.Regression curve;

	String msg, prefix, LINK_SUFFIX = ModelDescriptionParser.LINK_SUFFIX;
	Element e;


	try 
	  {
		// Variable links part of XML
		// ok, look for the PROPERTY_INDICES_TAG chunk
		NodeList nlist = null;
		try
		  {
			nlist = (NodeList)xPath.evaluate("/inversion/"+ModelDescriptionParser.PROPERTY_INDICES_TAG,
											 document,XPathConstants.NODESET);
		  }
		catch (Exception ex) 
		  {
			_m._SI.EH.println("Parsing exception: "+ex);
		  }



		if (nlist.getLength()<1) 
		  throw new SAXException("No " + ModelDescriptionParser.PROPERTY_INDICES_TAG + 
								 " entry in: " + filename);
		e = (Element)nlist.item(0);		
		msg=new String(ModelDescriptionParser.PROPERTY_INDICES_TAG);

		// loop over layers
		for (i=0;i<m.MPD.MetaInfo.length;i++)
		  {
			MI = m.MPD.MetaInfo[i];
			//loop over end members
			for (int end_type=0; end_type<=1 ; end_type++)
			  {
				rock_name = (end_type==0) ? MI.mud_rock_member : MI.sand_rock_member;
				prefix    = (end_type==0) ? NONRESERVOIR_PREFIX : RESERVOIR_PREFIX;

				if (rock_name==null) continue; // break out if there's nothing to do
				
				// get curve set from hashmap
				curveSet = (RockandFluidData.RegressionSet) m.RFD.RockTrends.get(rock_name);
				if (curveSet==null) throw new SAXException("No regression curves in layer "+i+" for rock "
												 + rock_name + " in file " + filename);
				// Explicitly code three trend curves
				
				// (1) Vp versus depth
				curve = curveSet.Vp_versus_depth;			
				if (curve.intercept_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.intercept_PD, 
															   prefix+VP_INTERCEPT+LINK_SUFFIX, msg,true);
				  }
				if (curve.depth_coeff_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.depth_coeff_PD, 
															   prefix+DEPTH_COEFF+LINK_SUFFIX, msg,true);
				  }
				if (curve.LFIV_coeff_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.LFIV_coeff_PD, 
															   prefix+LFIV_COEFF+LINK_SUFFIX, msg,true);
				  }
				if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
				  {
					if (curve.float_fraction_coeff_PD.varies)
					  {
						extractVariableIndex(e,0,MI,curve.float_fraction_coeff_PD, 
											 prefix+VP_FLOAT_FRACTION_COEFF+LINK_SUFFIX, msg,true);
					  }
				  }
				if (curve.error_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.error_PD, 
															   prefix+VP_SIGMA+LINK_SUFFIX, msg,true);
				  }

				// (2) Vs versus Vp
				curve = curveSet.Vs_versus_Vp;			
				if (curve.intercept_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.intercept_PD, 
															   prefix+VS_INTERCEPT+LINK_SUFFIX, msg,true);
				  }
				if (curve.slope_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.slope_PD, 
															   prefix+VS_VP_COEFF+LINK_SUFFIX, msg,true);
				  }
				if (curve.error_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.error_PD, 
															   prefix+VS_SIGMA+LINK_SUFFIX, msg,true);
				  }

				// (3) rho or phi versus Vp
				curve = curveSet.Rho_or_Phi_versus_Vp;			
				if (curve.intercept_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.intercept_PD, prefix+
																((end_type==0) ? FACTOR : PHI_INTERCEPT)
																+LINK_SUFFIX, msg,true);
				  }
				if (curve.slope_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.slope_PD, prefix+
																((end_type==0) ? EXPONENT : PHI_VP_COEFF)
																+LINK_SUFFIX, msg,true);
				  }
				if (_m._SI.MF.ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
				  {
					if (curve.float_fraction_coeff_PD.varies && (end_type==1))
					  {
						// for reservoir member only
						extractVariableIndex(e,0,MI,curve.float_fraction_coeff_PD, prefix+ PHI_FLOAT_FRACTION_COEFF+LINK_SUFFIX, msg,true);
					  }
				  }
				if (curve.error_PD.varies)
				  {
					extractVariableIndex(e,0,MI,curve.error_PD, prefix+
																((end_type==0) ? RHO_SIGMA : PHI_SIGMA)
																+LINK_SUFFIX, msg,true);
				  }
			  } // end members
		  } // layers
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
	  }
  }

  /** 
	Set the variable index of property PD , tagged by tag, within the xml node
	e. If needed, it must be found. If we fail, throw an exception. Must be
	kept identical to similar method in ModelDescriptionParser
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
					throw new SAXException("In: " + filename + 
										   ", Not an integer: tag " + tag + " element " + element_name);
				  }
				// we seem to have got a reasonable index for tag_name in this Element
				
				if (index>0) 
				  {		
					// add to useful ArrayList of changing things
					PD.index = index;
					_m.MPD.PriorIndexInfo.add(new ModelPriorDescriptions.propertyInfo(layer_no,PD));
					EH.println("tag = "+ tag + ", PriorIndexInfo SIZE = "+_m.MPD.PriorIndexInfo.size()+ " block = "+PD.index);
					if (MI.t.varies)
					  {
						// t, t_base have the largest index
						if (test)
						  {
							// is a property other than time, so must have smaller index
							if (index >= MI.t.index)
							  throw new SAXException("if " + tag + " varies, its index must be less than "
													 + MI.t.index + ", "+ this.filename);
						  }
						else
						  {
							// must be t or t_base
							if (index != MI.t.index)
							  throw new SAXException("since " + tag + " varies, its index *must* be "
													 + MI.t.index + ", "+ this.filename);
						  }
					  }
					return;
				  }
				else
				  throw new SAXException("In: " + filename + 
										 ", negative index for " 
										 + tag + " in element " + element_name);
			  }
		  }
		// if we got here, we're in trouble
		throw new SAXException("In: " + filename + 
							   ", can't find an variable index for " 
							   + tag + " in element " + element_name);
	  }
  }
}// RockandFluidDataParser


