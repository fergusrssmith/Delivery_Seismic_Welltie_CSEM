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

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.XmlUtils;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;



import javax.xml.parsers.*;
import org.xml.sax.*;
import org.w3c.dom.*;
import java.io.*;

import java.util.*;

import cern.colt.map.*;

public final class ParameterFileParser
{
  /** 
	  Error Handler for instance of this class
  */
  ErrorHandler EH;

  Document document;

  String xml_filename=null;

  int verbosity = 0;

  /**
	 Tag for true variable
  */
  public static String TRUE_TAG = "true";

  /**
	 Tag for well_info in parameter file
  */
  public static String WELL_INFO_TAG = "well_info";

  /**
	 Tag for well in parameter file
  */
  public static String WELL_TAG = "well";

  /**
	 Tag for name in parameter file
  */
  public static String NAME_TAG = "name";
/**
	 Tag for well_log_info in parameter file
  */
  public static String WELL_LOG_INFO_TAG = "well_log_info";
  /**
	 Tag for well_log_filename in parameter file
  */
  public static String WELL_LOG_FILENAME_TAG = "well_log_filename";

  /**
	 Tag for measured_depth_keyword in parameter file
  */
  public static String MEASURED_DEPTH_KEYWORD_TAG = "measured_depth_keyword";

  /**
	 Tag for sonic_DT_keyword in parameter file
  */
  public static String DT_KEYWORD_TAG = "sonic_DT_keyword";

  /**
	 Tag for shear_sonic_DT_keyword in parameter file
  */
  public static String SHEAR_DT_KEYWORD_TAG = "shear_sonic_DT_keyword";
  /**
	 Tag for density_keyword in parameter file
  */
  public static String DENSITY_KEYWORD_TAG = "density_keyword";
  /**
	 Tag for gamma_keyword in parameter file
  */
  public static String GAMMA_KEYWORD_TAG = "gamma_keyword";
  /**
	 Tag for MD upper mute in parameter file
  */
  public static String MD_UPPER_MUTE_KEYWORD_TAG = "MD_uppermute";

  /**
	 Tag for MD lower mute in parameter file
  */ 
  public static String MD_LOWER_MUTE_KEYWORD_TAG = "MD_lowermute";

  /**
	 Tag for well_survey_info in parameter file
  */
  public static String WELL_SURVEY_INFO_TAG = "well_survey_info";
  /**
	 Tag for well_survey_filename in parameter file
  */
  public static String WELL_SURVEY_FILENAME_TAG = "well_survey_filename";
  /**
	 Tag for total_vertical_depth_keyword in parameter file
  */
  public static String TOTAL_VERTICAL_DEPTH_KEYWORD_TAG = "total_vertical_depth_keyword";
  /**
	 Tag for x_position_keyword in parameter file
  */
  public static String X_POSITION_KEYWORD_TAG = "x_position_keyword";
  /**
	 Tag for y_position_keyword in parameter file
  */
  public static String Y_POSITION_KEYWORD_TAG = "y_position_keyword";
  /**
	 Tag for checkshot_info in parameter file
  */
  public static String CHECKSHOT_INFO_TAG = "checkshot_info";
  /**
	 Tag for checkshot_filename in parameter file
  */
  public static String CHECKSHOT_FILENAME_TAG = "checkshot_filename";
  /**
	 Tag for checkshot_time_keyword in parameter file
  */
  public static String CHECKSHOT_TIME_KEYWORD_TAG = "checkshot_time_keyword";
  /**
	 Tag for sigma_checkshot_time_keyword in parameter file
  */
  public static String SIGMA_CHECKSHOT_TIME_KEYWORD_TAG = "sigma_checkshot_time_keyword";
  /**
	 Tag for checkshot_MD_keyword in parameter file
  */
  public static String CHECKSHOT_MD_KEYWORD_TAG = "measured_depth_keyword";
  /**
	 Tag for markers_info in parameter file
  */
  public static String MARKERS_INFO_TAG = "markers_info";
  /**
	 Tag for markers_filename in parameter file
  */
  public static String MARKERS_FILENAME_TAG = "markers_filename";
  /**
	 Tag for marker_time_keyword in parameter file
  */
  public static String MARKERS_TIME_KEYWORD_TAG = "marker_time_keyword";
  /**
	 Tag for sigma_marker_time_keyword in parameter file
  */
  public static String SIGMA_MARKERS_TIME_KEYWORD_TAG = "sigma_marker_time_keyword";
  /**
	 Tag for markers_MD_keyword in parameter file
  */
  public static String MARKERS_MD_KEYWORD_TAG = "measured_depth_keyword";
  /**
	 Tag for seismic_stack_info in parameter file
  */
  public static String SEISMIC_STACK_INFO_TAG = "seismic_stack_info";
  /**
	 Tag for seismic_minicube_filename in parameter file
  */
  public static String SEISMIC_MINICUBE_FILENAME_TAG = "seismic_minicube_filename";
  /**
	 Tag for min_offset in parameter file
  */
  public static String MIN_OFFSET_TAG = "min_offset";
  /**
	 Tag for max_offset in parameter file
  */
  public static String MAX_OFFSET_TAG = "max_offset";
  /**
	 Tag for reflector_time in parameter file
  */
  public static String REFLECTOR_TIME_TAG = "reflector_time_ms";
  /**
	 Tag for stack_velocity in parameter file
  */
  public static String STACK_VELOCITY_TAG = "stack_velocity";
  /**
	 Tag for extraction_info in parameter file
  */
  public static String EXTRACTION_INFO_TAG = "extraction_info";
  /**
	 Tag for extraction_window_start_time in parameter file
  */
  public static String EXTRACTION_WINDOW_START_TIME_TAG = "extraction_window_start_time";
  /**
	 Tag for extraction_window_end_time in parameter file
  */
  public static String EXTRACTION_WINDOW_END_TIME_TAG = "extraction_window_end_time";
  /**
	 Tag for wavelet_parameters in parameter file
  */
  public static String WAVELET_PARAMETERS_TAG = "wavelet_parameters";
  /**
	 Tag for maximum_coda_time_millis in parameter file
  */
  public static String MAXIMUM_CODA_TIME_MILLIS_TAG = "maximum_coda_time_millis";
  /**
	 Tag for maximum_precursor_time_millis in parameter file
  */
  public static String MAXIMUM_PRECURSOR_TIME_MILLIS_TAG = "maximum_precursor_time_millis";
  /**
	 Tag for seismic_minicube_header_words in parameter file
  */
  public static String SEISMIC_MINICUBE_HEADER_WORDS_TAG = "seismic_minicube_header_words";
 
  /**
	 Tag for X_header_keyword in parameter file
  */
  public static String X_HEADER_KEYWORD_TAG = "X_header_keyword";
  /**
	 Tag for Y_header_keyword in parameter file
  */
  public static String Y_HEADER_KEYWORD_TAG = "Y_header_keyword";
  /**
	 Tag for inline_header_keyword in parameter file
  */
  public static String INLINE_HEADER_KEYWORD_TAG = "inline_header_keyword";
  /**
	 Tag for xline_header_keyword in parameter file
  */
  public static String XLINE_HEADER_KEYWORD_TAG = "xline_header_keyword";


  /**
	 Tag for single trace minicube in parameter file
  */
  public static String IS_SINGLE_TRACE_ONLY_KEYWORD_TAG = "is-single-trace-only";
  /**
	 Tag for constant angle stack
  */
  public static String CONSTANT_ANGLE_STACK_KEYWORD_TAG = "is-constant-angle-stack-with-angle-degrees";
  /**
	 Tag for PS data stack
  */
  public static String PS_STACK_KEYWORD_TAG = "is-PS-data";
  /**
	 Tag for utm section in parameter file
  */
  public static String TWOD_UTM_CONVERSIONS_KEYWORD_TAG = "two-D-utm-conversions";
  /**
	 Tag for utm geometry-setting section in parameter file
  */
  public static String SET_2D_GEOM_FROM_HEADERS_TAG = "set-two-D-geometry-from-XY-header-words";
  /**
	 Tag for utm section in parameter file
  */
  public static String UTM_CONVERSIONS_KEYWORD_TAG = "utm-conversions";
  /**
	 Tag for utm point in parameter file
  */
  public static String UTM_POINT_KEYWORD_TAG = "point";
  /**
	 Tag for utm x in parameter file
  */
  public static String UTM_X_KEYWORD_TAG = "x";
  /**
	 Tag for utm y in parameter file
  */
  public static String UTM_Y_KEYWORD_TAG = "y";
  /**
	 Tag for utm inline in parameter file
  */
  public static String UTM_INLINE_KEYWORD_TAG = "inline";
  /**
	 Tag for utm xline in parameter file
  */
  public static String UTM_XLINE_KEYWORD_TAG = "xline";

  /**
	 Tag for AVO_terms_prior_distribution in parameter file
  */
  public static String AVO_TERMS_PRIOR_DISTRIBUTION_TAG = "AVO_terms_prior_distribution";
  /**
	 Tag for B_coefficient_prior_mean in parameter file
  */
  public static String B_COEFFICIENT_PRIOR_MEAN_TAG = "B_coefficient_prior_mean";
  /**
	 Tag for B_coefficient_prior_std_deviation in parameter file
  */
  public static String B_COEFFICIENT_PRIOR_STD_DEVIATION_TAG = "B_coefficient_prior_std_deviation";


  /**
	 Tag for far wavelet stuff in parameter file
  */
  public static String FAR_WAVELET_STRETCH_PRIOR_DISTRIBUTION_TAG = "far_wavelet_stretch_prior_distribution";
  /**
	 Tag for far_wavelet_time_stretch_prior_mean in parameter file
  */
  public static String FAR_WAVELET_TIME_STRETCH_COEFFICIENT_PRIOR_MEAN_TAG = "far_wavelet_time_stretch_prior_mean";
  /**
	 Tag for far_wavelet_time_stretch_prior_std_deviation in parameter file
  */
  public static String FAR_WAVELET_TIME_STRETCH_COEFFICIENT_PRIOR_STD_DEVIATION_TAG = "far_wavelet_time_stretch_prior_std_deviation";
/**
	 Tag for far_wavelet_scale_prior_mean in parameter file
  */
  public static String FAR_WAVELET_SCALE_PRIOR_MEAN_TAG = "far_wavelet_scale_prior_mean";
  /**
	 Tag for far_wavelet_scale_prior_std_deviation in parameter file
  */
  public static String FAR_WAVELET_SCALE_PRIOR_STD_DEVIATION_TAG = "far_wavelet_scale_prior_std_deviation";



  /**
	 Tag for systematic time uncertainty mean info-block in seismic minicube block parameter file (in ms)
  */
  public static String SYSTEMATIC_TIME_SHIFT_MEAN_TAG = "time_to_depth_local_systematic_time_shift_mean_ms";
  /**
	 Tag for systematic time uncertainty info-block in seismic minicube block parameter file (in ms)
  */
  public static String SYSTEMATIC_TIME_SHIFT_TAG = "time_to_depth_local_systematic_time_shift_uncertainty_ms";
  /**
	 Tag for systematic time uncertainty info-block in parameter file
  */
  public static String SYSTEMATIC_WELL_POSITION_SHIFT_BLOCK_TAG = "well_positioning_systematic_shift_uncertainty";
  /**
	 Tag for systematic X position offset bias in parameter file
  */
  public static String SYSTEMATIC_X_SHIFT_OFFSET_TAG = "x_position_offset";
  /**
	 Tag for systematic Y position offset bias in parameter file
  */
  public static String SYSTEMATIC_Y_SHIFT_OFFSET_TAG = "y_position_offset";
  /**
	 Tag for systematic X position uncertainty info in parameter file
  */
  public static String SYSTEMATIC_X_SHIFT_SIGMA_TAG = "sigma_x_position";
  /**
	 Tag for systematic Y position uncertainty info in parameter file
  */
  public static String SYSTEMATIC_Y_SHIFT_SIGMA_TAG = "sigma_y_position";
  /**
	 Tag for synchronising well positions and times with other well
  */
  public static String SYNCHRONIZE_SHIFTS_WITH_OTHER_WELL_TAG = "synchronize_shifts_with_other_well";

  /**
	 Tag for output in parameter file
  */
  public static String OUTPUT_TAG = "output";
  /**
	 Tag for maximum_likelihood_wavelet_filename in parameter file
  */
  public static String MAXIMUM_LIKELIHOOD_WAVELET_FILENAME_TAG = "maximum_likelihood_wavelet_filename";
  /**
	 Tag for maximum_likelihood_parameters_filename in parameter file
  */
  public static String MAXIMUM_LIKELIHOOD_PARAMETERS_FILENAME_TAG = "maximum_likelihood_parameters_filename";
  /**
	 Tag for wavelet_realisations_filename in parameter file
  */
  public static String WAVELET_REALISATIONS_FILENAME_TAG = "wavelet_realisations_filename";
  /**
	 Tag for parameter_realisations_filename in parameter file
  */
  public static String PARAMETER_REALISATIONS_FILENAME_TAG = "parameter_realisations_filename";
  /**
	 Tag for realisation_number_header_field in parameter file
  */
  public static String REALISATION_NUMBER_HEADER_FIELD_TAG = "realisation_number_header_field";

  /**
	 Tag for XXX in parameter file
  */
  public static String XXX_TAG = "XXX";



  private Model _m;
  private WaveletExtractor _WE;

  /**
	 Make a parser for the named file, which will populate certain
	 fields of the Model m, with the desired certain level of
	 verbosity and handle errors using the ErrorHandler
  */
  public ParameterFileParser(String filename, WaveletExtractor WE)
  {
	this.verbosity = WE.verbosity;
	this._WE = WE;
	this.EH = WE.EH;
	this.xml_filename = filename;
  }

  /** 
	  General call to read layer descriptions into relevant fields in model.
  */  
  public boolean readParameters()
  {
	boolean success;
	 
	DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
	 
	try 
	  {
		DocumentBuilder builder = factory.newDocumentBuilder();
		document = builder.parse( new File(xml_filename));
		
		// (1) Get well and seismic data
		success = getWellandSeismicInfoData();
		if (!success) return false;

		// (2) Get other data (wavelet, extraction window etc)
		success = getOtherInfo();
		if (!success) return false;
	  }
	catch (Exception ex) 
	  {
		ex.printStackTrace();
		return false;
	  };
	return true;
  }
  
  private boolean getWellandSeismicInfoData() throws SAXException
  {
	int i,i_w,j,j_s,k,nWellsXML,nStacksXML=0;
	String msg_frag;
	Node n;
	NodeList nListA, nListB, nListC,listB, listC;
	Element elem=null,elemA=null,elemC=null,elemS;
	NamedNodeMap nnmap=null;
	Attr nv;
	String attr_val;
	OpenIntIntHashMap wMap = null;  // --wells
	OpenIntIntHashMap sMap = null;  // --stacks


	try 
	  {
		// Preliminary parse of allowed stacks
		if (_WE.useRestrictedStackSet)
		  {
			String[] st = _WE.restrictedStackNumbers.split("[,;:]");
			sMap = new OpenIntIntHashMap();
			for (i=0; i< st.length; i++)
			  {
				try
				  {
					j = Integer.parseInt(st[i]);
				  }
				catch (NumberFormatException e)
				  {
					 throw new RuntimeException("--stacks "+_WE.restrictedStackNumbers+" contains a non-integer");
				  }
				if (j<1)
				  {
					throw new RuntimeException("--stacks  "+_WE.restrictedStackNumbers+" contains number "+j+"<1");
				  }
				sMap.put(j,j);
			  }
			// OK, by here, restricted list has at least sane entries
			sMap.trimToSize();
			if (sMap.size()==0)
			  {
				throw new RuntimeException("--stacks "+_WE.restrictedStackNumbers+" permits no stacks at all!");				
			  }
		  }




		// ok, look for the WELL_INFO_TAG chunk
		nListA = document.getElementsByTagName(WELL_INFO_TAG);
		if (nListA.getLength()<1) 
		  throw new SAXException("No " + WELL_INFO_TAG + " entry in: " + xml_filename);

		// There should be only one entry. lets take the first
		// & look for the WELL_TAG chunk
		NodeList nlistB = ((Element)nListA.item(0)).getElementsByTagName(WELL_TAG);
		_WE.m.nWells = nlistB.getLength();
		nWellsXML = _WE.m.nWells;

		// reduce this if restrictions apply
		if (_WE.useRestrictedWellSet)
		  {
			String[] st = _WE.restrictedWellNumbers.split("[,;:]");
			if (st.length > nWellsXML)
			  {
			   throw new SAXException("Fewer wells in " + xml_filename+" than available to subselect using --wells "+_WE.restrictedWellNumbers);
			  }
			wMap = new OpenIntIntHashMap();
			for (i=0; i< st.length; i++)
			  {
				try
				  {
					j = Integer.parseInt(st[i]);
				  }
				catch (NumberFormatException e)
				  {
					 throw new SAXException("--wells "+_WE.restrictedWellNumbers+" contains a non-integer");
				  }
				if (j<1 || j>nWellsXML)
				  {
					throw new SAXException("--wells "+_WE.restrictedWellNumbers+" contains number "+j+
										   ", not in range [1.."+nWellsXML+"] (number of wells in XML)");
				  }
				wMap.put(j,j);
			  }
			// OK, by here, restricted list should be OK
			_WE.m.nWells = st.length;
			wMap.trimToSize();				
		  }


		if (_WE.m.nWells <1) throw new SAXException("Cannot proceed: no "+WELL_TAG+
											   " entry in: " + xml_filename);
		// allocate room for everything
		_WE.wellInfo = new GeneralWellInfo[_WE.m.nWells];
		
		// loop over wells
		for (i=0,i_w=0; i<nWellsXML; i++)
		  {
			if (_WE.useRestrictedWellSet)
			  {
				if (!wMap.containsKey(i+1)) 
				  {
					continue; // skip this well - it's not in the list
					// will skip i_w++ at end of brace too
				  }
			  }

			_WE.wellInfo[i_w] = new GeneralWellInfo();

			msg_frag = " for well "+(i+1);
			
			elemA = ((Element)nlistB.item(i));
			GeneralWellInfo _WI = _WE.wellInfo[i_w]; // shortcut

			_WI.well_name  = extractValue(elemA,NAME_TAG , msg_frag); // "fake_well";

			// Well log info

			listB = elemA.getElementsByTagName(WELL_LOG_INFO_TAG);
			if (listB.getLength()<1) 
			  throw new SAXException("No " + WELL_LOG_INFO_TAG + " entry in: " + xml_filename);
			elem = (Element)listB.item(0); 
			_WI.logFilename  = extractValue(elem,WELL_LOG_FILENAME_TAG , msg_frag); // "fake_well.txt";

			_WI.wellLogDefinitions.MD  = extractValue(elem,MEASURED_DEPTH_KEYWORD_TAG , msg_frag); // "MD";
			if (_WI.wellLogDefinitions.MD == null) 
			  throw new SAXException("No " + MEASURED_DEPTH_KEYWORD_TAG + " entry in: " + xml_filename);

			_WI.wellLogDefinitions.DT  = extractValue(elem,DT_KEYWORD_TAG , msg_frag); // "DT";
			if (_WI.wellLogDefinitions.DT  == null) 
			  throw new SAXException("No " + DT_KEYWORD_TAG + " entry in: " + xml_filename);

			_WI.wellLogDefinitions.DTshear  = extractValue(elem,SHEAR_DT_KEYWORD_TAG , msg_frag); // "DT_shear";
			if (_WI.wellLogDefinitions.DTshear == null)
			  {
				if (!_WE.fakeVsLog)
				  {
					throw new SAXException("No " + SHEAR_DT_KEYWORD_TAG + " entry in: " + xml_filename);
				  }
				// Shear keyword not supplied, ignore shear
				_WI.wellLogDefinitions.shearTreatment = 0; // ignore
			  }
			
			listC = elem.getElementsByTagName(SHEAR_DT_KEYWORD_TAG);
			if (listC.getLength()>0)
			  { 
				elemS = (Element)listC.item(0);
				// Seek attributes controlling shear estimation
				nnmap = elemS.getAttributes();
				if (nnmap != null)
				  {
					nv = (Attr)(nnmap.getNamedItem("use"));
					if (nv != null)
					  {
						attr_val = nv.getValue();
						if (attr_val != null)
						  {
							if (attr_val.equals("false"))
							  {
								_WI.wellLogDefinitions.shearTreatment = 0; // ignore 
							  }
						  } 
					  }
				  }
			  }

			_WI.wellLogDefinitions.rho  = extractValue(elem,DENSITY_KEYWORD_TAG , msg_frag); // "RHOB";
			if (_WI.wellLogDefinitions.rho == null) 
			  throw new SAXException("No " + DENSITY_KEYWORD_TAG + " entry in: " + xml_filename);

			// optional gamma
			_WI.wellLogDefinitions.gamma  = extractValue(elem,GAMMA_KEYWORD_TAG , msg_frag); // "gamma";

			// optional mutes 
			_WI.wellLogDefinitions.MD_uppermute  = extractNumber(elem,MD_UPPER_MUTE_KEYWORD_TAG,false,false);
			_WI.wellLogDefinitions.MD_lowermute  = extractNumber(elem,MD_LOWER_MUTE_KEYWORD_TAG,false,false);
  
			// Well survey info

			listB = elemA.getElementsByTagName(WELL_SURVEY_INFO_TAG);
			if (listB.getLength()<1) 
			  throw new SAXException("No " + WELL_SURVEY_INFO_TAG + " entry in: " + xml_filename);
			elem = (Element)listB.item(0); 
			_WI.surveyFilename  = extractValue(elem,WELL_SURVEY_FILENAME_TAG , msg_frag);

			_WI.geom.Zm  = extractValue(elem,MEASURED_DEPTH_KEYWORD_TAG , msg_frag);
			if (_WI.geom.Zm  == null) 
			  throw new SAXException("No " + MEASURED_DEPTH_KEYWORD_TAG + " entry in: " + xml_filename);

			_WI.geom.Z  = extractValue(elem,TOTAL_VERTICAL_DEPTH_KEYWORD_TAG, msg_frag);
			if (_WI.geom.Z == null) 
			  throw new SAXException("No " + TOTAL_VERTICAL_DEPTH_KEYWORD_TAG + " entry in: " + xml_filename);

			_WI.geom.X  = extractValue(elem,X_POSITION_KEYWORD_TAG , msg_frag);
			if (_WI.geom.X == null) 
			  throw new SAXException("No " + X_POSITION_KEYWORD_TAG + " entry in: " + xml_filename);

			_WI.geom.Y  = extractValue(elem,Y_POSITION_KEYWORD_TAG , msg_frag);
			if (_WI.geom.Y == null) 
			  throw new SAXException("No " + Y_POSITION_KEYWORD_TAG + " entry in: " + xml_filename);

			// Checkshot info

			listB = elemA.getElementsByTagName(CHECKSHOT_INFO_TAG);
			if (listB.getLength()<1) 
			  throw new SAXException("No " + CHECKSHOT_INFO_TAG + " entry in: " + xml_filename);
			elem = (Element)listB.item(0); 
			_WI.checkshotFilename  = extractValue(elem,CHECKSHOT_FILENAME_TAG , msg_frag);
			_WI.chkshotAndMarkerDefinitions.T_checkshot  = extractValue(elem,CHECKSHOT_TIME_KEYWORD_TAG , msg_frag);
			_WI.chkshotAndMarkerDefinitions.sigmaT_checkshot  = extractValue(elem,SIGMA_CHECKSHOT_TIME_KEYWORD_TAG, msg_frag);
			_WI.chkshotAndMarkerDefinitions.Zm_checkshot  = extractValue(elem,CHECKSHOT_MD_KEYWORD_TAG, msg_frag);
			// Possible auto-checkshot entry
			nnmap = elem.getAttributes();
			if (nnmap != null)
			  {
				nv = (Attr)(nnmap.getNamedItem("automake-from-sonic-log"));
				if (nv != null)
				  {
					attr_val = nv.getValue();
					if (attr_val != null && attr_val.equals("true")) _WI.checkshotAutoMake = true;
				  }
				if (_WI.checkshotAutoMake)
				  {
					nv = (Attr)(nnmap.getNamedItem("model-select-knots"));
					if (nv != null)
					  {
						attr_val = nv.getValue();
						if (attr_val != null && attr_val.equals("true")) _WI.checkshotAutoMakeModelSelect = true;
					  }
				  }
			  }

			// Marker info. This may be omitted

			listB = elemA.getElementsByTagName(MARKERS_INFO_TAG);
			if (listB.getLength()>0)
			  { 
				elem = (Element)listB.item(0); 
				_WI.markersFilename  = extractValue(elem,MARKERS_FILENAME_TAG , msg_frag);
				_WI.chkshotAndMarkerDefinitions.T_marker  = extractValue(elem,MARKERS_TIME_KEYWORD_TAG , msg_frag);
				_WI.chkshotAndMarkerDefinitions.sigmaT_marker  = extractValue(elem,SIGMA_MARKERS_TIME_KEYWORD_TAG, msg_frag);
				_WI.chkshotAndMarkerDefinitions.Zm_marker  = extractValue(elem,MARKERS_MD_KEYWORD_TAG, msg_frag);
			  }
			else
			  {
				// This may be omitted
				// throw new SAXException("No " + MARKERS_INFO_TAG + " entry in: " + xml_filename);
			  }
			// (a) well position
			listB = elemA.getElementsByTagName(SYSTEMATIC_WELL_POSITION_SHIFT_BLOCK_TAG);
			if (listB.getLength()>0)
			  { 
				elem = (Element)listB.item(0);
				_WI.systematicXProp.setMean(extractNumber(elem, SYSTEMATIC_X_SHIFT_OFFSET_TAG,false,true));
				_WI.systematicYProp.setMean(extractNumber(elem, SYSTEMATIC_Y_SHIFT_OFFSET_TAG,false,true));
				_WI.systematicXProp.setSigma(extractNumber(elem, SYSTEMATIC_X_SHIFT_SIGMA_TAG,true,true));
				_WI.systematicYProp.setSigma(extractNumber(elem, SYSTEMATIC_Y_SHIFT_SIGMA_TAG,true,true));
			  }

			// (b) possible synchronization
			listB = elemA.getElementsByTagName(SYNCHRONIZE_SHIFTS_WITH_OTHER_WELL_TAG);
			if (listB.getLength()>0)
			  { 
				elem = (Element)listB.item(0);
				_WI.synchronizeWellName  = extractValue(elem, NAME_TAG, msg_frag);
			  }


			// Seismic stacks/minicubes info

			listB = elemA.getElementsByTagName(SEISMIC_STACK_INFO_TAG);
			if (listB.getLength()<1) 
			  throw new SAXException("No " + SEISMIC_STACK_INFO_TAG + " entry in: " + xml_filename);
			if (_WE.m.nStacks==0)
			  {
				// first pass (i.e. first well)
				nStacksXML = listB.getLength();

				if (_WE.useRestrictedStackSet)
				  {
					// now check restricted stack list for sanity
					for (j=0;j<nStacksXML; j++) if (sMap.containsKey(j+1)) _WE.m.nStacks++;
					if (_WE.m.nStacks != sMap.size())
					  {
						throw new RuntimeException("--stacks "+_WE.restrictedStackNumbers+
											   " contains unreasonable stack number. Must be in [1.."+nStacksXML+"]");
					  }
				  }
				else
				  {
					_WE.m.nStacks = nStacksXML;
				  }
				_WE.m.initialiseStretchParameters();
				_WE.seismicInfo = new SeismicStackInfo[nWellsXML][_WE.m.nStacks];
				_WE.lastWaveletFilenames = new String[_WE.m.nStacks];
			  }
			else if (nStacksXML != listB.getLength())
			  {
				throw new SAXException("Number of seismic stacks available for well" + 	_WI.well_name + " does not "
									   +"match that of other wells, in file " +xml_filename);
			  }
			if (_WE.m.nStacks==1 || _WE.synchroniseStackTregistrationErrors)
			  {
				_WE.m.registrationErrorsAreCommon = true;
			  }
			for (j=0,j_s=0;j<nStacksXML; j++)
			  {
				if (_WE.useRestrictedStackSet)
				  {
					if (!sMap.containsKey(j+1)) 
					  {
						continue; // skip this stack - it's not in the list
					  }
				  }
				// seismic parsing
				elem = (Element)listB.item(j);

				String filename  = extractValue(elem,SEISMIC_MINICUBE_FILENAME_TAG, msg_frag);
				double min_offset = extractNumber(elem,MIN_OFFSET_TAG,  true, true);
				double max_offset = extractNumber(elem,MAX_OFFSET_TAG,  true, true);
				// stack velocity in seconds
				double reflector_time = 
				  FundamentalConstants.autoTimeConvertMsToS(extractNumber(elem,REFLECTOR_TIME_TAG,  true, true));
				double stack_velocity = extractNumber(elem,STACK_VELOCITY_TAG,  true, true);

				_WE.seismicInfo[i_w][j_s] = new SeismicStackInfo(filename,max_offset,min_offset,reflector_time,stack_velocity);
				SeismicStackInfo _SSI = _WE.seismicInfo[i_w][j_s];

				// Is it a single trace only???
				// is the utm info obtainable from XY header words? 
				nnmap = elem.getAttributes();
				if (nnmap != null)
				  {
					nv = (Attr)(nnmap.getNamedItem(IS_SINGLE_TRACE_ONLY_KEYWORD_TAG));
					if (nv != null)
					  {
						attr_val = nv.getValue();
						if (attr_val != null && attr_val.equals(TRUE_TAG))
						  {
							_SSI.is1DdataAndUseSnapping = true;
						  }
					  }
					// is constant angle?
					nv = (Attr)(nnmap.getNamedItem(CONSTANT_ANGLE_STACK_KEYWORD_TAG));
					if (nv != null)
					  {
						attr_val = nv.getValue();
						if (attr_val != null)
						  {
							_SSI.theta_specified = Double.parseDouble(attr_val);
							if (_SSI.theta_specified >= 0)
							  {
								if (_SSI.theta_specified >= 90)
								  {
									throw new SAXException("Specified constant-angle attribute for seismic stack "+ filename + 
														   " must be in range [0,90] (realistically, probably [0-50])");
								  } 
								// convert to radians
								_SSI.theta_specified *= Math.PI/180.0;
								// & set up angles...
								_SSI.computeStackParams();
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
							_SSI.isPSdata = true;
							_WI.needPSreflectivities=true; // mark this to assist segmentation
						  }
					  }
				  } 

				// Systematic time shifts (inputs in ms): These may be omitted, and thus not modelled
				_SSI.systematicTProp.setMean(0.001*extractNumber(elem, SYSTEMATIC_TIME_SHIFT_MEAN_TAG ,false, false));
				if (Math.abs(_SSI.systematicTProp.mean) > SeismicStackInfo.MAXIMUM_SYSTEMATIC_TIME_SCALE_SECS)
				  {
					throw new SAXException("Systematic time mean for seismic stack "+ filename + 
										   " cannot be reasonably more than "
										   + SeismicStackInfo.MAXIMUM_SYSTEMATIC_TIME_SCALE_SECS+ "(s)");
				  }

				// Systematic time shifts uncertainty
				_SSI.systematicTProp.setSigma(0.001*
																 extractNumber(elem, SYSTEMATIC_TIME_SHIFT_TAG ,true,false));
				if (_SSI.systematicTProp.sigma > SeismicStackInfo.MAXIMUM_SYSTEMATIC_TIME_SCALE_SECS)
				  {
					throw new SAXException("Systematic time uncertainty for seismic stack "+ filename + 
										   " cannot be reasonably more than "
										   + SeismicStackInfo.MAXIMUM_SYSTEMATIC_TIME_SCALE_SECS+ "(s)");
				  }

				// 2D UTM conversions if any
				// look for the 2DUTM_CONVERSIONS_KEYWORD_TAG chunk
				listC = elem.getElementsByTagName(TWOD_UTM_CONVERSIONS_KEYWORD_TAG);
				if (listC.getLength()>0)
				  {
					// OK, its a 2D line.
					if (_SSI.is1DdataAndUseSnapping)
					  {
						throw new SAXException("Minicube data has been flagged as both single trace and a line in "+ filename + 
											   " Please correct XML suitably.");
					  }
					_SSI.is2DdataAndUseSnapping = true;

					elem = (Element)listC.item(0);

					// is the utm info obtainable from XY header words? 
					nnmap = elem.getAttributes();
					if (nnmap != null)
					  {
						nv = (Attr)(nnmap.getNamedItem(SET_2D_GEOM_FROM_HEADERS_TAG));
						if (nv != null)
						  {
							attr_val = nv.getValue();
							if (attr_val != null && attr_val.equals(TRUE_TAG))
							  {
								_SSI.is2DdataAndUseXYwords  = true;
							  }
						  }
					  }
					else
					  {
						// OK, can't trust XY header words. Use 2 UTM points following in the XML file
						_SSI.is2DdataAndUseXYwords  = false;

						nListC = elem.getElementsByTagName(UTM_POINT_KEYWORD_TAG);
						if (nListC.getLength()!=2)
						  {
							throw new SAXException("Fewer or more than 2 2D utm geometry-setting points in: " + xml_filename);
						  }
						// read the 2 utm points. They must all be "inline" points (xLine=0)
						_SSI.utmInfo = new UTMConversionInfo();
						UTMConversionInfo utmInfo = _SSI.utmInfo;
						for (k=0;k<2;k++)
						  {
							elem = (Element)nListC.item(k); 
							utmInfo.x[k] = extractNumber(elem, UTM_X_KEYWORD_TAG ,false, true);
							utmInfo.y[k] = extractNumber(elem, UTM_Y_KEYWORD_TAG ,false, true);
							utmInfo.inLine[k] = extractNumber(elem, UTM_INLINE_KEYWORD_TAG ,false, true);
							utmInfo.xLine[k] = extractNumber(elem, UTM_XLINE_KEYWORD_TAG ,false, true);
						  }
						if (utmInfo.xLine[0] != utmInfo.xLine[1])
						  {
							throw new SAXException("2D utm geometry-setting points must have common xLine in: " 
												   + xml_filename);
						  }
						// set the third point by trivial 90deg rotation, preserving squareness of geometry
						// Points then end up like this:          2    1
						//                                          0
						utmInfo.x[2] = utmInfo.x[0]+utmInfo.y[0]-utmInfo.y[1];
						utmInfo.y[2] = utmInfo.y[0]+utmInfo.x[1]-utmInfo.x[0];
						utmInfo.inLine[2] = utmInfo.inLine[0]+utmInfo.xLine[0]-utmInfo.xLine[1];
						utmInfo.xLine[2]  = utmInfo.xLine[0]+utmInfo.inLine[1]-utmInfo.inLine[0];

						// setup UTMS and return false if they're bogus
						if (!utmInfo.setUpUTMmapping(_WE.EH)) return false;
					  }
				  }
				j_s++; // to manage stack subselection
			  }

			// Per-well timelapse info
			listB = elemA.getElementsByTagName("time_lapse_well_info");
			if (listB.getLength()>0) 
			  {
				_WE.useSaturationModel = true;
				_WI.setupSaturationInfo();
				elem = (Element)listB.item(0);
				_WI.saturationInfo.Sb = extractNumber(elem,"base_saturation",true,true);
				_WI.saturationInfo.Sa = extractNumber(elem,"mean_monitor_saturation",true,true);
				_WI.saturationInfo.sig_Sa = extractNumber(elem,"sigma_monitor_saturation",true,true);
			  }

			i_w++;	// to manage well subselection
		  }

		// End of well parsing



		// Check: if synchronising, check we have suitable targets
		// loop over wells
		for (i=0;i<_WE.m.nWells;i++)
		  {
			// we will use this index whether synchronised or not
			_WE.wellInfo[i].synchronizeWellIndex = i;
			if (_WE.wellInfo[i].synchronizeWellName != null)
			  {
				// This well thinks it's synchronised. Find its master among the other wells
				boolean foundTarget=false;
				for (j=0;j<nWellsXML ;j++)
				  {
					// cannot synchronize to oneself (skip)
					if (j==i) continue;
					// look at different wells...
					if (_WE.wellInfo[j].well_name.equals(_WE.wellInfo[i].synchronizeWellName))
					  {
						// found a possible master. Check master itself is not synchronised
						if (_WE.wellInfo[j].synchronizeWellName != null)
						  {
							// target is not an independent parameterised master-well
							// Reference is nonsensical
							throw new SAXException("Synchronization target well "+
												   _WE.wellInfo[j].well_name+" specified for well "
												   +_WE.wellInfo[i].well_name+ " is synchronised to another well. "
												   +" (one can only synchronise to a master well). In: " 
												   + xml_filename);
						  }
						else
						  {
							// master looks good
							foundTarget = true;
							// override index
							_WE.wellInfo[i].synchronizeWellIndex = j;
							break;
						  }
					  }
				  }
				if (!foundTarget)
				  {
					throw new SAXException("No synchronization target well "+
										   _WE.wellInfo[i].synchronizeWellName+" found, as specified for well "+
										   _WE.wellInfo[i].well_name+", in: " + xml_filename);
				  }
			  }
		  }

		// seismic header words

		// look for the SEISMIC_MINICUBE_HEADER_WORDS_TAG chunk
		nListA = document.getElementsByTagName(SEISMIC_MINICUBE_HEADER_WORDS_TAG);
		if (nListA.getLength()<1) 
		  throw new SAXException("No " + SEISMIC_MINICUBE_HEADER_WORDS_TAG + " entry in: " + xml_filename);


		elem = (Element)nListA.item(0); 
		
		_WE.seismicXword = extractValue(elem,X_HEADER_KEYWORD_TAG,null);
		_WE.seismicYword = extractValue(elem,Y_HEADER_KEYWORD_TAG,null);
		if (_WE.seismicXword.equals(_WE.seismicYword))
		  {
			throw new SAXException("Can't have "+X_HEADER_KEYWORD_TAG+" == "
								   +Y_HEADER_KEYWORD_TAG+" , in: " + xml_filename);
		  }
		_WE.seismicInlineWord = extractValue(elem,INLINE_HEADER_KEYWORD_TAG,null);
		_WE.seismicXlineWord= extractValue(elem,XLINE_HEADER_KEYWORD_TAG,null);
		if (_WE.seismicInlineWord.equals(_WE.seismicXlineWord))
		  {
			throw new SAXException("Can't have "+INLINE_HEADER_KEYWORD_TAG+" == "+
								   XLINE_HEADER_KEYWORD_TAG+" , in: " + xml_filename);
		  }
		// and any (optional) utm data inside this
		// look for the UTM_CONVERSIONS_KEYWORD_TAG chunk
		nListA = document.getElementsByTagName(UTM_CONVERSIONS_KEYWORD_TAG);
		if (nListA.getLength()>0)
		  {
			_WE.useUTMconversion = true;
			elem = (Element)nListA.item(0);
			nListB = elem.getElementsByTagName(UTM_POINT_KEYWORD_TAG);
			if (nListB.getLength()<3)
			  {
				throw new SAXException("Fewer or more than 3 utm gemotry-setting points in: " + xml_filename);
			  }
			// read the 3 utm points
			_WE.utmInfo = new UTMConversionInfo();
			for (i=0;i<3;i++)
			  {
				elem = (Element)nListB.item(i); 
				_WE.utmInfo.x[i] = extractNumber(elem, UTM_X_KEYWORD_TAG ,false, true);
				_WE.utmInfo.y[i] = extractNumber(elem, UTM_Y_KEYWORD_TAG ,false, true);
				_WE.utmInfo.inLine[i] = extractNumber(elem, UTM_INLINE_KEYWORD_TAG ,false, true);
				_WE.utmInfo.xLine[i] = extractNumber(elem, UTM_XLINE_KEYWORD_TAG ,false, true);
			  }
			// setup UTMS and return false if they're bogus
			if (!_WE.utmInfo.setUpUTMmapping(_WE.EH)) return false;
		  }
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }
	return true;
  }

  
  private boolean getOtherInfo() throws SAXException
  {
	NodeList nListA, nListB, nListC;
	Element elem=null, elemB, elemC;

	try  
	  {
		//  Possible time-lapse surfaces
		String TIMELAPSE_INFO_TAG = "time_lapse_info";
		nListA = document.getElementsByTagName(TIMELAPSE_INFO_TAG);
		if (nListA.getLength()>0) 
		  {
			_WE.useSaturationModel = true;
			
			TimeLapseInfo _TL = _WE.m.timelapseinfo = new TimeLapseInfo();
			// check other well-specific xml info is OK
			for (int i=0;i<_WE.m.nWells ;i++)
			  {
				if(_WE.wellInfo[i].saturationInfo ==null)
				  throw new SAXException("No well specific time_lapse_well_info entries in: " + xml_filename);
			  }
			// parse
			elem = (Element)nListA.item(0);
			for(int i=0;i<2;i++)
			  {
				String s = (i==0) ? "base_contacts_filename" : "monitor_contacts_filename";
				String a = extractValue(elem,s,null);
				if(i==0) _TL.base_contacts_filename = a;
				else _TL.monitor_contacts_filename = a;

				nListA = elem.getElementsByTagName(s);
				if (nListA.getLength()>0)
				  {
					elemB = (Element)nListA.item(0); 
					NamedNodeMap nnmap = elemB.getAttributes();
					if (nnmap != null)
					  {
						Attr nv = (Attr)(nnmap.getNamedItem("time_unit"));
						if (nv != null)
						  {
							String attr_val = nv.getValue();
							if (attr_val != null)
							  {
								boolean t = attr_val.equals("s");
								if(i==0) _TL.base_timeunitssecs = t;
								else _TL.monitor_timeunitssecs = t;
							  } 
						  }
					  }
				  }
			  }
			_TL.vp_s = extractNumber(elem,"vp_sand_representative_SI",true,true);
			_TL.vs_s = extractNumber(elem,"vs_sand_representative_SI",true,true);
			_TL.phi_s = extractNumber(elem,"phi_sand_representative",true,true);
			_TL.vp_m = extractNumber(elem,"vp_shale_representative_SI",true,true);
			_TL.vs_m = extractNumber(elem,"vs_shale_representative_SI",true,true);
			_TL.rho_m = extractNumber(elem,"rho_shale_representative_SI",true,true);
			_TL.vp_h = extractNumber(elem,"vp_h_representative_SI",true,true);
			_TL.rho_h = extractNumber(elem,"rho_h_representative_SI",true,true); 
			_TL.vp_b = extractNumber(elem,"vp_brine_representative_SI",true,true);
			_TL.rho_b = extractNumber(elem,"rho_brine_representative_SI",true,true);
			_TL.rho_g = extractNumber(elem,"rho_grain_representative_SI",true,true);
			_TL.K_g = extractNumber(elem,"K_grain_representative_SI",true,true);
		  }

		//  Wavelet extraction window
		nListA = document.getElementsByTagName(EXTRACTION_INFO_TAG);
		if (nListA.getLength()<1) 
		  throw new SAXException("No " + EXTRACTION_INFO_TAG + " entry in: " + xml_filename);

		elem = (Element)nListA.item(0); 
		_WE.tStart = extractNumber(elem,EXTRACTION_WINDOW_START_TIME_TAG,true,true);
		_WE.tStart = FundamentalConstants.autoTimeConvertMsToS(_WE.tStart);
		_WE.tEnd = extractNumber(elem,EXTRACTION_WINDOW_END_TIME_TAG,true,true);
		_WE.tEnd = FundamentalConstants.autoTimeConvertMsToS(_WE.tEnd);


		//  Wavelet parameters
		// wavelet min and max precursor/coda times be sure to QC these
		nListA = document.getElementsByTagName(WAVELET_PARAMETERS_TAG);
		if (nListA.getLength()<1) 
		  throw new SAXException("No " + WAVELET_PARAMETERS_TAG + " entry in: " + xml_filename);

		elem = (Element)nListA.item(0); 
		// convert both of these to seconds 
		_WE.tCodaMax = 0.001*extractNumber(elem,MAXIMUM_CODA_TIME_MILLIS_TAG,true,true);
		// this must be negative, so insist on a positive input????
		_WE.tPrecursorMin = -0.001*extractNumber(elem,MAXIMUM_PRECURSOR_TIME_MILLIS_TAG,true,true);


		//  AVO terms: these are (are optional: absence will leave defaults as B ~ N(1,0) )
		nListA = document.getElementsByTagName(AVO_TERMS_PRIOR_DISTRIBUTION_TAG);
		if (nListA.getLength()>0) 
		  {
			elem = (Element)nListA.item(0); 
			_WE.m.B_AVO[0].setMean(extractNumber(elem,B_COEFFICIENT_PRIOR_MEAN_TAG,true,false));
			_WE.m.B_AVO[0].setSigma(extractNumber(elem,B_COEFFICIENT_PRIOR_STD_DEVIATION_TAG, true,false));
			// FIXME SHOULD BE ABLE TO SET INDEP FOR EACH STACK!!!
			if (_WE.freeAVO)
			  {
				for (int j=1; j<_WE.m.nStacks; j++)
				  {
					_WE.m.B_AVO[j].setMean(_WE.m.B_AVO[0].mean);
					_WE.m.B_AVO[j].setSigma(_WE.m.B_AVO[0].sigma);
				  }
			  }
		  }

		// far stack stretch/amplitude terms
		nListA = document.getElementsByTagName(FAR_WAVELET_STRETCH_PRIOR_DISTRIBUTION_TAG);
		if (nListA.getLength()>0) 
		  {
			_WE.useStretchedModel = true;
			elem = (Element)nListA.item(0);
			_WE.m.farStackWaveletTimeStretchMean = extractNumber(elem,FAR_WAVELET_TIME_STRETCH_COEFFICIENT_PRIOR_MEAN_TAG,true,false);

			// if mean stretch is illegal, bail out
			if (_WE.m.farStackWaveletTimeStretchMean <1 || _WE.m.farStackWaveletTimeStretchMean>=Wavelet.MAX_WAVELET_STRETCH)
			  {
				throw new SAXException("Illegal wavelet stretch mean: must be in range (1,"+ Wavelet.MAX_WAVELET_STRETCH+") in "
									   + xml_filename);
			  }
			_WE.m.farStackWaveletTimeStretchSigma = extractNumber(elem,FAR_WAVELET_TIME_STRETCH_COEFFICIENT_PRIOR_STD_DEVIATION_TAG,true,false);

			_WE.m.farStackWaveletAmplitudeMean = extractNumber(elem,FAR_WAVELET_SCALE_PRIOR_MEAN_TAG,true,false); 
			_WE.m.farStackWaveletAmplitudeSigma = extractNumber(elem,FAR_WAVELET_SCALE_PRIOR_STD_DEVIATION_TAG,true,false);

			// if the entries are trivial we can switch stretching off
			if (_WE.m.farStackWaveletTimeStretchMean==1 && _WE.m.farStackWaveletTimeStretchSigma == 0
				&& _WE.m.farStackWaveletAmplitudeMean==1 && _WE.m.farStackWaveletAmplitudeSigma == 0)
			  {
				_WE.useStretchedModel = false;
			  }
			if (_WE.m.nStacks==1) _WE.useStretchedModel = false;
			if (_WE.useStretchedModel)
			  {
				for (int i=0;i< (_WE.m.nStacks-1); i++)
				  {
					_WE.m.farStackWaveletTimeStretch[i].setMean(_WE.m.farStackWaveletTimeStretchMean); 
					_WE.m.farStackWaveletTimeStretch[i].setSigma(_WE.m.farStackWaveletTimeStretchSigma); 
					_WE.m.farStackWaveletAmplitude[i].setMean(_WE.m.farStackWaveletAmplitudeMean);
					_WE.m.farStackWaveletAmplitude[i].setSigma(_WE.m.farStackWaveletAmplitudeSigma);
				  }
			  }
		  }

		// Petrophysics data
		nListA = document.getElementsByTagName("petrophysics_info");
		if (nListA.getLength()>0)
		  {
			_WE.m.petro.petrophysicsProvided = true;
			elem = (Element)nListA.item(0);

			// thresholds
			nListB = elem.getElementsByTagName("shale_thresholds");
			if (nListB.getLength()>0)
			  {
				elemB = (Element)nListB.item(0);
				nListC = elemB.getElementsByTagName("GR_threshold");
				if (nListC.getLength()>0)
				  {
					elemC = (Element)nListC.item(0);
					_WE.m.petro.GRThreshProvided = true;
					_WE.m.petro.GR_intercept = (float)extractNumber(elemC,"GR_threshold_intercept",true,true);
					_WE.m.petro.GR_slope =  (float)extractNumber(elemC,"GR_threshold_TVD_slope",true,true);
				  }
				nListC = elemB.getElementsByTagName("density_threshold");
				if (nListC.getLength()>0)
				  {
					elemC = (Element)nListC.item(0);
					_WE.m.petro.densityThreshProvided = true;
					_WE.m.petro.density_intercept = (float)extractNumber(elemC,"density_threshold_intercept",true,true);
					_WE.m.petro.density_slope = (float) extractNumber(elemC,"density_threshold_TVD_slope",true,true);
				  }
				if (!(_WE.m.petro.GRThreshProvided || 	_WE.m.petro.densityThreshProvided))
				  {
					throw new SAXException("Must provide either a GR or density shale threshold in petrophysics_info/shale_thresholds in"
										   + xml_filename);
				  }
			  }

			// lumped shale Thomsen delta prior
			nListB = elem.getElementsByTagName("lumped_shale_Thomsen_delta_prior_distribution");
			if (nListB.getLength()>0)
			  {
				_WE.doAnisotropicModel = true;

				_WE.m.deltaThomsen_shale = new Property(0.05,false,true);
				elemB = (Element)nListB.item(0);				
				_WE.m.deltaThomsen_shale.setMean(extractNumber(elemB,"Thomsen_delta_prior_mean",false,true));
				_WE.m.deltaThomsen_shale.setSigma(extractNumber(elemB,"Thomsen_delta_prior_std_deviation",true,true));
			  }
		  }
	  }
	catch (SAXException ex)
	  {
		ex.printStackTrace();
		return false;
	  }
	return true;
  }

  public void summariseParameters()
  {
	int i,j;
	String s;
	EH.printf("Runtime arguments:%n");
	EH.println(_WE.runCommand);


	EH.printf("%nParameter summary:%n");

	// loop over wells
	for (i=0;i<_WE.m.nWells ;i++)
	  {
		GeneralWellInfo wellInfo = _WE.wellInfo[i];
		
		s = wellInfo.well_name.replaceAll("%","%%");  
		EH.println("For well: "+ s);
		EH.println("           " + " log filename = " +wellInfo.logFilename.replaceAll("%","%%"));
		EH.println("           " + " meas. depth keyword = "+wellInfo.wellLogDefinitions.MD);
		EH.println("           " + " sonic log keyword   = "+wellInfo.wellLogDefinitions.DT);
		EH.println("           " + " shear-log keyword   = "+wellInfo.wellLogDefinitions.DTshear);
		EH.println("           " + " density keyword     = "+wellInfo.wellLogDefinitions.rho);
		if (wellInfo.wellLogDefinitions.gamma != null)
		  {
			EH.println("           " + " gamma keyword   = "+wellInfo.wellLogDefinitions.gamma);
		  }
		if(!Double.isNaN(wellInfo.wellLogDefinitions.MD_uppermute))
		  {
			EH.println("           " + " data muted above MD = "+wellInfo.wellLogDefinitions.MD_uppermute);
		  }
		if(!Double.isNaN(wellInfo.wellLogDefinitions.MD_lowermute))
		  {
			EH.println("           " + " data muted below MD = "+wellInfo.wellLogDefinitions.MD_lowermute);
		  }
		EH.println("");
		EH.println("           " + " survey filename = " + wellInfo.surveyFilename.replaceAll("%","%%"));
		EH.println("           " + " survey MD  keyword = "	+ wellInfo.geom.Zm);
		EH.println("           " + " survey TVD keyword = "+ wellInfo.geom.Z);
		EH.println("           " + " survey X   keyword = "	+ wellInfo.geom.X); 
		EH.println("           " + " survey Y   keyword = "	+ wellInfo.geom.Y);
		EH.println("");
		EH.println("           " + " checkshot filename = " + wellInfo.checkshotFilename.replaceAll("%","%%"));
		EH.println("           " + " checkshot T keyword = " + wellInfo.chkshotAndMarkerDefinitions.T_checkshot);
		EH.println("           " + " checkshot sigmaT keyword = " + wellInfo.chkshotAndMarkerDefinitions.sigmaT_checkshot);
		EH.println("           " + " checkshot MD keyword = "+	wellInfo.chkshotAndMarkerDefinitions.Zm_checkshot);
		EH.println("");
		if (wellInfo.markersFilename==null)
		  {
			EH.println("           " + " No markers supplied!");
		  }
		else
		  {
			EH.println("           " + " marker filename = " + wellInfo.markersFilename.replaceAll("%","%%"));
			EH.println("           " + " marker T keyword = " + wellInfo.chkshotAndMarkerDefinitions.T_marker);
			EH.println("           " + " marker sigmaT keyword = " + wellInfo.chkshotAndMarkerDefinitions.sigmaT_marker);
			EH.println("           " + " marker MD keyword = "+	wellInfo.chkshotAndMarkerDefinitions.Zm_marker);
		  }
		if (wellInfo.synchronizeWellName == null)
		  {
			EH.println("           " + " systematic X shift mean = %5.2g",wellInfo.systematicXProp.mean);
			if (wellInfo.systematicXProp.sigma>0)
			  {
				EH.println("           " + " systematic X shift prior uncertainty = %5.2g",wellInfo.systematicXProp.sigma);
			  }
			EH.println("           " + " systematic Y shift mean = %5.2g",wellInfo.systematicYProp.mean);
			if (wellInfo.systematicYProp.sigma>0)
			  {
				EH.println("           " + " systematic Y shift prior uncertainty = %5.2g",wellInfo.systematicYProp.sigma);
			  }
		  }
		else
		  {
			EH.println("           " + " systematic XY offsets synchronized with "+wellInfo.synchronizeWellName);
		  }
		
		
		
		for (j=0;j<_WE.m.nStacks; j++)
		  {
			EH.println("");
			EH.println("           " + " minicube(stack) filename = " + _WE.seismicInfo[i][j].filename);
			if (_WE.seismicInfo[i][j].isPSdata) EH.println("            stack modelled as PS data");
			EH.println("           " + " stack Xmin, Xmax = " + _WE.seismicInfo[i][j].offsetXmin + ", "+ _WE.seismicInfo[i][j].offsetXmax);
			EH.println("           " + " stack Tevent (secs), Vstack = %-6.3g, %-6.3g",
					   _WE.seismicInfo[i][j].T_event,_WE.seismicInfo[i][j].V_stack);
			if (_WE.seismicInfo[i][j].theta_specified >= 0)
			  {
				EH.println("           " + " stack angle (specified, degrees) = %5.2g",_WE.seismicInfo[i][j].theta_specified*180.0/Math.PI);				
			  }
			else if (_WE.seismicInfo[i][j].theta_specified == -Double.MAX_VALUE)
			  {
				EH.println("           " + " stack angle (approx., degrees) = %5.2g",_WE.seismicInfo[i][j].thetaApprox*180.0/Math.PI);
			  }
			else
			  {
				EH.println("            ..taking negative specified constant-angle attribute to mean non-constant-angle (i.e. constant ray-parameter, or ray-traced) stack");
				EH.println("            stack angle (approx., degrees) = %5.2g",_WE.seismicInfo[i][j].thetaApprox*180.0/Math.PI);
			  }
			if (wellInfo.synchronizeWellName == null)
			  {
				if ((_WE.seismicInfo[i][j]).systematicTProp.mean!=0)
				  {
					EH.println("           " + " systematic time shift mean (s) = %5.2g",(_WE.seismicInfo[i][j]).systematicTProp.mean);
				  }
				if ((_WE.seismicInfo[i][j]).systematicTProp.sigma>0)
				  {
					EH.println("           " + " systematic time shift uncertainty (s) = %5.2g",(_WE.seismicInfo[i][j]).systematicTProp.sigma);
				  }
			  }
			else
			  {
				EH.println("           " + " systematic time offsets synchronized with "+wellInfo.synchronizeWellName);
			  }
		  }
	  }

	EH.println("");
	EH.println("  Inline, Xline seismic keywords = " + _WE.seismicInlineWord + ", " +_WE.seismicXlineWord); 
	EH.println("  X, Y seismic keywords = " + _WE.seismicXword + ", " +_WE.seismicYword); 
	
	if (_WE.useUTMconversion)
	  {
		EH.println("  Setting geometry from UTM/XY triples");
	  }
	EH.println("");
	EH.println("  Wavelet time extraction window (secs) = [" + _WE.tStart+ "," +_WE.tEnd+"]"); 

	EH.println("");
	EH.println("  Wavelet maximum coda time (secs)     = "+_WE.tCodaMax); 
	EH.println("  Wavelet maximum precursor time (secs) = "+(-_WE.tPrecursorMin)); 

	EH.println("");
	EH.println("  AVO B coefficient (of theta^2) prior mean          = "+_WE.m.B_AVO[0].mean);
	EH.println("  AVO B coefficient (of theta^2) prior std deviation = "+_WE.m.B_AVO[0].sigma);

	EH.println("");

	if (_WE.doAnisotropicModel)
	  {
		EH.println("  Shale Thomsen epsilon prior mean          = "+_WE.m.deltaThomsen_shale.mean);
		EH.println("  Shale Thomsen epsilon prior std deviation = "+_WE.m.deltaThomsen_shale.sigma);
		
		if (_WE.m.petro.petrophysicsProvided)
		  {
			if (_WE.m.petro.GRThreshProvided)
			  {
				EH.println("  Anisotropic model: shale GR threshold: "+_WE.m.petro.GR_intercept+" + "+_WE.m.petro.GR_slope+" x TVD");
			  }
			if (_WE.m.petro.densityThreshProvided)
			  {
				EH.println("  Anisotropic model: shale density threshold: "+_WE.m.petro.density_intercept+" + "+_WE.m.petro.density_slope+" x TVD");
			  }
		  }
		EH.println("");
	  }

	if (_WE.useStretchedModel)
	  {
		if (_WE.m.farStackWaveletTimeStretchMean != 1 || _WE.m.farStackWaveletTimeStretchSigma!=0)
		  {
			EH.println("  Far stack stretch coefficient prior mean          = "+_WE.m.farStackWaveletTimeStretchMean);
			EH.println("  Far stack t-stretch coefficient prior std deviation = "+_WE.m.farStackWaveletTimeStretchSigma);
		  }
		if (_WE.m.farStackWaveletAmplitudeMean != 1 || _WE.m.farStackWaveletAmplitudeSigma!=0)
		  {
			EH.println("  Far stack wavelet scaling coefficient prior mean          = "+_WE.m.farStackWaveletAmplitudeMean);
			EH.println("  Far stack wavelet scaling prior std deviation = "+_WE.m.farStackWaveletAmplitudeSigma);
			EH.println("");
		  }
	  }

	if(_WE.useSaturationModel)
	  {
		EH.println("  Time-lapse mode: treating seismic data as difference, adding saturation/Gassman machinery");
	 	EH.println("");
	  }
	
	EH.println("Runtime modifications:");
	if (_WE.doShuey)
	  {
		EH.println("Using Shuey equation for reflectivity");
	  }
	else if (_WE.doIsoZoeppritz)
	  {
		EH.println("Using isotropic Zoeppritz equation for reflectivity");
	  }
	if (_WE.freeAVO)
	  {
		EH.println("AVO B coefficients free and independent per stack");
	  }
	if (_WE.ignoreSeismic)
	  {
		EH.println("Seismic likelihood switched off");
	  }
	if (_WE.computeForLongestWaveletOnly)
	  {
		EH.println("Longest wavelet only");
	  }
	if (_WE.fixNoise)
	  {
		EH.println("Noise level fixed at "+_WE.fixNoiseValue);
	  }
	if (!_WE.constrainIntervalVelocities)
	  {
		EH.println("Interval velocity constraints disabled");
	  }
	else
	  {
		EH.println("Interval velocity constrained to (normalised) std. dev. of %10.8g",
				   _WE.intervalVelocityConstraintStdDev);
	  }
	if (_WE.constrainPhase)
	  {
		EH.println("Phase constraint: target = %10.8g std-dev = %10.8g",
				   _WE.constantPhase,_WE.constantPhaseStdDev);
	  }
	if (_WE.constrainPeakArrival)
	  {
		EH.println("Peak arrival constraint: target = %10.8g std-dev = %10.8g",
				   _WE.peakArrivalTime,_WE.peakArrivalStdDev);
	  }
	if (_WE.useBIC)
	  {
		EH.println("Mode probabilities (wavelet length) chosen using BIC");
	  }
	if (_WE.fixCheckshots)
	  {
		EH.println("Time-to-depth (checkshots) fixed");
	  }
	if (_WE.fixMinicubeRegistrations)
	  {
		EH.println("Time-to-depth (minicube registrations) fixed");
	  }
	if (_WE.fixWellPositions)
	  {
		EH.println("Well positions fixed");
	  }
	if (_WE.assumeNoReflectionsOutsideLogs)
	  {
		EH.println("Assuming no reflections outside log data");
	  }
	if (_WE.allowLogInterpolation)
	  {
		EH.println("Allowing log interpolation of missing log data over maximum distance "+
				   _WE.logInterpolationMaxDist);
	  }
	if (_WE.eliminateSpikes)
	  {
		EH.println("Eliminating log spikes based on DT log");
	  }
	if (_WE.synchroniseStackTregistrationErrors)
	  {
		EH.println("Synchronizing stack registration errors if any");
	  }

	if (_WE.doAnisotropicModel && _WE.checkAnisotropicModel)
	  {
		EH.println("Anisotropic model: checking models without anisotropy also (Thomsen delta=0)");
	  }
	if (_WE.useWaveletPrior)
	  {
		EH.println("Using prior wavelet specs gleaned from "+_WE.wavelet_prior_filename);
	  }
	//....

	EH.println("");
  }

  /** 
	  Return the text value from the named Element e, as tagged by
	  tag_name (returning the first entry found in the tree below the
	  current point). If we find it but don't understand it send
	  failure message to the error handler.  If not found, return null
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
			throw new SAXException("In: " + xml_filename + ", unintelligible entry: for " 
								   + tag_name + " " +message);
		  }
	  }
	// we didn't find the <tag_name> entry, so return null
	return null;
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
}
