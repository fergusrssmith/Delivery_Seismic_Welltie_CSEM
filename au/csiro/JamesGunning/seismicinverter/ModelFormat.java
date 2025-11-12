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
import java.util.*;

/**
   Class containing simple fields describing the output format of the
   SeismicInverter code. The output format has to be intelligible to other
   binaries, so it is encapsulated in this small class and shared amongst
   the various executables.
*/
public final class ModelFormat
{
  /**
	 A switch to enable captured floating grain model
  */
  final boolean ENABLE_CAPTURED_FLOATING_GRAIN_MODEL = true;

  /**
	 A switch to enable extra fields for simple contacts model
  */
  final boolean ENABLE_SIMPLE_CONTACTS_MODEL = true;

  /**
	 A switch to enable hard-rock model
  */
  final boolean ENABLE_HARD_ROCK_MODEL = true; 

  /**
	 A switch to enable extra shear fields
  */
  final boolean ENABLE_EXTRA_SHEAR_FIELDS = true;

  /**
	 complete_field_names is a list of field-names from the LayerProperties class that are output, and
	 other "code-words" that indicate particular fields or computed quantities
	 to output. It is assembled by the Class constructor.
	 NB, ensure the map LayerProperties.completePropertyNameMap
	 is kept current with this, and also the SeismicInverter.dumpRealisation() method
	 Code words indicated by CODE, otherwise LayerProperties fields
  */
  String[] complete_field_names =null;

  /**
	Basic list of fields available in early versions. Will always be present.
	Later versions augment this list, and update the format_version field.
	NB, ensure the map LayerProperties.completePropertyNameMap
	is kept current with this, and also the SeismicInverter.dumpRealisation() method
	Code words indicated by CODE, otherwise LayerProperties fields
  */
  String[] fundamental_field_names = {"d",          // depth
								   "LFIV",       // low freq interval velocity
								   "NG",         // net to gross
								   "vp_m",       // impermeable rock p-wave velocity
								   "vs_m",       // "" s-wave velocity
								   "rho_m",      // "" density
								   "vp_s",       // permeable rock p-wave velocity
								   "vs_s",       // "" s-wave velocity
								   "phi_s",      // "" porosity
								   "S_h",        // hydrocarbon sat (or brine if layer is brine-filled (CODE)
								   "rho_h",      // hydrocarbon density (or brine "") (CODE)
								   "fluid_type", // 0,1,2,3 for brine, oil, lsg, gas (CODE)
								   "vp_eff",     // rock-mixed and fluid-substituted effective layer velocity  (CODE)
								   "rho_eff",    // rock-mixed and fluid-substituted effective layer density (CODE)
								   "Z_eff",      // rock-mixed and fluid-substituted effective layer impedance (CODE)
								   "R_near",     // near-stack reflection coefficient (CODE)
								   "R_far",      // far stack reflection coefficient (CODE)
								   "seismic_chisquare",  // normalised chisq/n for seismic (CODE)
								   "t",          // layer-top times. MUST LEAVE AS SECOND LAST ENTRY
								   // this can also be depth 
								   "t_base"      // bottom layer base-time. MUST LEAVE AS LAST ENTRY
								   // this can also be last layer bottom depth 
  };

  /**
	 Start of items that may be added to complete_field_names, depending on version
  */
  static String FLOAT_FIELD_NAME = "float_fraction"; // fraction of float in reservoir rock (CODE)
  
  // Fluid contacts extra fields
  static String T_FC_FIELD_NAME = "t_fc";                // fluid contact time 
  static String D_FC_FIELD_NAME = "d_fc";                // fluid contact depth 
  static String R_NEAR_FC_FIELD_NAME = "R_near_fc";      // fluid contacts near refl coeff.
  static String R_FAR_FC_FIELD_NAME = "R_far_fc";        // fluid contacts far  refl coeff.
  static String S_OIL_EFF_FIELD_NAME = "S_oil";          // Effective Oil saturation
  static String S_LSG_EFF_FIELD_NAME = "S_lsg";          // Effective Lsg saturation
  static String S_GAS_EFF_FIELD_NAME = "S_gas";          // Effective Gas saturation


  // Hard rocks - mixing in of second impermeable facies
  static String NG_HARD_FIELD_NAME = "NG_hard";          // fraction of hard rock in non-reservoir volume-space
  static String VP_HARD_FIELD_NAME = "vp_hard";          // fraction of hard rock in non-reservoir volume-space
  static String VS_HARD_FIELD_NAME = "vs_hard";          // fraction of hard rock in non-reservoir volume-space
  static String RHO_HARD_FIELD_NAME = "rho_hard";        // fraction of hard rock in non-reservoir volume-space

  // shear info
  static String VS_EFF_FIELD_NAME = "vs_eff";             // effective shear speed
  static String ZS_EFF_FIELD_NAME = "Zs_eff";             // effective shear modulus
  
  
  /**
	 SU header word for the field which labels the version number.
	 Must correspond to a SU short (2 bytes).
  */
  public static final String VERSION_NUMBER_HEADER_FIELD = "duse";
  /**
	 SU header word for the field which contains the number of stacks.
	 Must correspond to a SU short (2 bytes).
  */
  public static final String NUMBER_OF_STACKS_HEADER_FIELD = "nhs";

  /** 
	  Encoding of output format version
  */
  int format_version = 0;
  
  /**
	 Su header-word used to save model file info for partially suppressed outputs.
	 Must correspond to a SU int (4 bytes).
  */
  public static final String OUTPUT_FIELDS_SU_KEYWORD = "sdepth";

  private int fieldBitsSuppressed;
  /**
	 boolean representaion of fieldBitsSuppressed
  */
  boolean[] includeField;

  private int n_layers;

  /**
	 SU header word for the field which labels the realisation number
  */
  String realisationNumberHeaderField="mark";

  /**
	 Resulting number of values in the output trace (computed at object creation from complete_field_names)
  */
  int ns;
  /**
	 useful List of field names and their values
  */
  List reducedPropertyList;
  /**
	 useful hashMap of field names and their indices
  */
  Map<String,Integer> nameMap;

  /**
	 Quick reference value to determine whether key used in HashMap is t_base
  */
  int T_BASE=-1;

  /**
	 Simple constructor, usable to set up basic fields. Complete with a
	 initialiseModelFormat(n_layers) call as soon as n_layers is known.
	 Set setDefaults=true for all available features (in Inverter), or use
	 setDefaults=false && format_version=code from header field for older formats (Analyser).
  */
  public ModelFormat(boolean setDefaults, int required_format_version)
  {
	int i;
	ArrayList<String> cfl;

	// Assemble complete list
	cfl = new ArrayList<String>(Arrays.asList(fundamental_field_names));

	// ----------------Addition of extra features---------------------
	// Feature #1: Augment model format by extra fields if required
	if (ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		// Add floating grain piece if allowed
		// use bitwise comparisons if required_format_version is to be followed
		// ---- Floating grain sets bit 1 ----
		if (setDefaults || ((required_format_version & 1)==1))
		  {
			cfl.add(cfl.size() - 2, FLOAT_FIELD_NAME);
			format_version += 1;    // if extending to new features , use pwrs of 2 for binary encoding
		  }
	  }
	
	// Feature #2: Augment model format by extra fields if required
	if (ENABLE_SIMPLE_CONTACTS_MODEL)
	  {
		// Add simple-contacts fields
		// use bitwise comparisons if required_format_version is to be followed
		// ---- Simple contacts sets bit 2 ----
		if (setDefaults || ((required_format_version & 2)==2))
		  {
			cfl.add(cfl.size() - 2, T_FC_FIELD_NAME);
			cfl.add(cfl.size() - 2, D_FC_FIELD_NAME);
			cfl.add(cfl.size() - 2, R_NEAR_FC_FIELD_NAME);
			// Bug-fix of missing R_FAR_FC_FIELD_NAME
			// ---- Set bit 3 ----
			if (setDefaults || ((required_format_version & 4)==4))
			  {
				cfl.add(cfl.size() - 2, R_FAR_FC_FIELD_NAME);
				format_version += 4;    // if extending to new features , use pwrs of 2 for binary encoding
			  }
			cfl.add(cfl.size() - 2, S_OIL_EFF_FIELD_NAME);  // (CODE)
			cfl.add(cfl.size() - 2, S_LSG_EFF_FIELD_NAME);  // (CODE)
			cfl.add(cfl.size() - 2, S_GAS_EFF_FIELD_NAME);  // (CODE)
			format_version += 2;    // if extending to new features , use pwrs of 2 for binary encoding
		  }
	  }

	// Feature #3: Hard rocks
	if (ENABLE_HARD_ROCK_MODEL)
	  {
		// Hard rock pieces
		// use bitwise comparisons if required_format_version is to be followed
		// ---- sets bit 4: 1248 ----
		//                     ^
		if (setDefaults || ((required_format_version & 8)==8))
		  {
			cfl.add(cfl.size() - 2, NG_HARD_FIELD_NAME);
			cfl.add(cfl.size() - 2, VP_HARD_FIELD_NAME);
			cfl.add(cfl.size() - 2, VS_HARD_FIELD_NAME);
			cfl.add(cfl.size() - 2, RHO_HARD_FIELD_NAME);
			format_version += 8;    // if extending to new features , use pwrs of 2 for binary encoding
		  }
	  }

	// Feature #4: Export shear info
	if (ENABLE_EXTRA_SHEAR_FIELDS)
	  {
		// use bitwise comparisons if required_format_version is to be followed
		// ---- sets bit 5: 1248,16 ----
		//                       ^
		if (setDefaults || ((required_format_version & 16)==16))
		  {
			cfl.add(cfl.size() - 2, VS_EFF_FIELD_NAME);
			cfl.add(cfl.size() - 2, ZS_EFF_FIELD_NAME);
			format_version += 16;    // if extending to new features , use pwrs of 2 for binary encoding
		  }
	  }

	// ---------------- Clean Up ------------
	cfl.trimToSize();
	complete_field_names = cfl.toArray(new String[0]);

	// ---- Supplementary fields we need ---
	reducedPropertyList = (List) Arrays.asList(complete_field_names);
	includeField = new boolean[complete_field_names.length];	
	nameMap = new HashMap<String,Integer>();
	for (i=0;i<complete_field_names.length;i++) 
	  {
		nameMap.put(complete_field_names[i],new Integer(i));
		includeField[i]=true;
	  }
	fieldBitsSuppressed = 0;
	// where t-base occurs in map...
	T_BASE = ((Integer)nameMap.get("t_base")).intValue();
  }

  /** 
	  Return a new modelFormat class, dropping the floating pieces
  */
  public ModelFormat getNewModelFormatSansFloatingGrainFields()
  {
	// Feature #1
	if (ENABLE_CAPTURED_FLOATING_GRAIN_MODEL)
	  {
		// mask off bit 1
		if ((format_version & 1)==1)
		  {
			format_version -= 1;
			return new ModelFormat(false,format_version);
		  }
		else return this;
	  }
	else return this;
  }

  /** 
	  Return a new modelFormat class, dropping the hard rock pieces
  */
  public ModelFormat getNewModelFormatSansHardRockFields()
  {
	// Feature #3: Hard rocks
	if (ENABLE_HARD_ROCK_MODEL)
	  {
		// mask off bit 4
		if ((format_version & 8)==8)
		  {
			format_version -= 8;
			return new ModelFormat(false,format_version);
		  }
		else return this;
	  }
	else return this;
  }


  /**
	 Method to set a flag encoding the suppression of field i in the master list
	 complete_field_names[]
  */
  public void suppressField(int i)
  {
	if (i>=0 && i<complete_field_names.length)
	  {
		fieldBitsSuppressed = (fieldBitsSuppressed | (1 << i));
		includeField[i]=false;
	  }
  }

  /**
	 Method to complete construction from incomplete ModelFormat()
	 constructor. Call as soon as n_layers is known.
  */
  public void initialiseModelFormat(int n_layers)
  {
	reducedPropertyList = new ArrayList();
	for (int i=0;i<complete_field_names.length - 1;i++) 
	  {
		if (includeField[i])
		  {
			// field is present
			reducedPropertyList.add(complete_field_names[i]);
			ns += n_layers;
		  }
	  }
	// make sure tbase is on
	includeField[complete_field_names.length-1] = true;

	ns++; // for time of t_base

  }



  /**
	 Method to compute internal fields etc from ns and the fieldMaskingbit. Usage e.g.:
	 recomputeLayersAndMaps((int)tr.get_short("ns"),
	 (int)tr.get_int(ModelFormat.OUTPUT_FIELDS_SU_KEYWORD));
	 A full model has fieldBitsSuppressed=0 naturally,
	 and ns will be read from the su header field.

	 The masking field should always ensure the last two fields (t,tbase) are not masked off
  */
  void recomputeLayersAndMaps(int ns, int fieldBitsSuppressed)
  {
	int i;
	int nFields;

	this.ns=ns;
	this.fieldBitsSuppressed = fieldBitsSuppressed;

	reducedPropertyList = new ArrayList();
	nameMap = new HashMap<String,Integer>();
	nFields=0;
	for (i=0; i<complete_field_names.length;i++) 
	  {
		if (((fieldBitsSuppressed >> i) & 1) == 0)
		  {
			// field is present
			includeField[i] = true;
			reducedPropertyList.add(complete_field_names[i]);
			nameMap.put(complete_field_names[i],new Integer(nFields));
			nFields++;
		  }
		else
		  {
			includeField[i] = false;
		  }
	  }

	// Since nFields includes tbase,
	// ns = 1 (from tbase)  + (nFields-1)*nLayers;
	n_layers = (ns-1)/(nFields-1);
  }

  /**
	 returns integer coding of suppressed fields
  */
  int getFieldBitsSuppressed()
  {
	return fieldBitsSuppressed;
  }

  /**
	 Returns number of layers, best called after recomputeLayersAndMaps(.,.)
	 reinitialises the internal maps and arrays.
  */
  int getNumberOfLayers()
  {
	return n_layers;
  }
}
