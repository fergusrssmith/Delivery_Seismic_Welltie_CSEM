/**********************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion
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
package au.csiro.JamesGunning.CSEM.OneDinversion;


import java.io.*;
import java.nio.*;
import java.util.HashMap;
import java.util.Arrays;

import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.util.*; 

import au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.DenseFloatMatrix2D;
import au.csiro.JamesGunning.lib.Numerical.DenseIntMatrix2D;

/**
   Class to provide JNI interface to new f90 
   1D CSEM modelling routines from the 
   Scripps; CSEM1D.f90

   Simple test for loading...

   java au.csiro.JamesGunning.CSEM.OneDinversion.ScrippsModelInterface
*/

public final class ScrippsModelInterface 
{
  /* Static members */

  private static boolean haveNativeLibrary = true;

  // public static native int testNative(int i,float[] y);


  /**
	 Location of native library builds wrt root of java packages.
	 This location contains builds in the directories<p>

	 ScrippsOneDRelativeDir/Linux-i686-unknown/libScrippsModelInterface.so,<p>
	 ScrippsOneDRelativeDir/Linux-x86_64-unknown/libScrippsModelInterface.so,<p>
	 ScrippsOneDRelativeDir/Darwin-Power-Macintosh-unknown/libScrippsModelInterface.jnilib<p>

	 etc etc
  */
  public static final String ScrippsOneDRelativeDir = "au/csiro/JamesGunning/CSEM/OneDinversion";
  /*
	name of library
  */
  public static final String libName = "libScrippsModelInterface";

  /* JNI interface members */

  //  public static native void setModelPointersNative(int np, rho_ptr);
  public static native void doEBFieldAllocs(int n);
  public static native void doJacobianAllocs(int n,int m);
  public static native void freeJacobianAllocs();
  public static native void setSimpleTxRxVariables(int nTxIn, int nRxIn);
  public static native void initDefaultsDipole1D();

  public static native void setOccamModPointers(int np, int n, int cnDataParams,int[] dp_ptr, 
												float[] d_ptr, float[] sd_ptr,float[] dm_ptr);
  public static native void releaseOccamModPointers(int[] dp_ptr, float[] d_ptr, float[] sd_ptr,float[] dm_ptr);
	
  public static native void setCsemModPointers 
	(
	 int nTxIn_v, double[] azimuthTxIn, double[] dipTxIn, double[] xTxIn, double[] yTxIn, double[] zTxIn,                    // allocateTX
	 int nFreqIn_v,double[] fTxIn,                                                                                           // allocateF  
	 int nRxIn_v, double[] xRxIn, double[] yRxIn, double[] zRxIn, double[] ThetaRxIn, double[] AlphaRxIn, double[] BetaRxIn, // allocateR,
	 int nRxTfm_v,double[]  xRxTfm,double[]  yRxTfm,double[]  zRxTfm, int[] iRxTfm,                                          // allocateRfm
	 int nd, int[] dpTfm,                                                                                                     // allocate_dpTfm,
	 int nTxTfm_v, double[] xTxTfm, double[] yTxTfm, double[] zTxTfm, double[] azimuthTxTfm, double[] dipTxTfm,               // allocateTfm 
	 float[] d_copy, float[] sd_copy,                                                                                         // allocateDcopies 
	 double[] sigsite, int [] RxInd, int[] rx1DtoRxIn,                                                                        // allocateDM
	 int nlay1D_v, double[] layer_params,                                                                                     // allocateL
	 int nRxRotSolve_v, int[] IndRxRotSolve                                                                                   // allocateIndRxRotSolve
	);
	

  public static native void releaseCsemModPointers 
	(
	 int nTxIn_v, double[] azimuthTxIn, double[] dipTxIn, double[] xTxIn, double[] yTxIn, double[] zTxIn,                    // allocateTX
	 int nFreqIn_v,double[] fTxIn,                                                                                           // allocateF  
	 int nRxIn_v, double[] xRxIn, double[] yRxIn, double[] zRxIn, double[] ThetaRxIn, double[] AlphaRxIn, double[] BetaRxIn, // allocateR,
	 int nRxTfm_v,double[]  xRxTfm,double[]  yRxTfm,double[]  zRxTfm, int[] iRxTfm,                                          // allocateRfm
	 int nd, int[] dpTfm,                                                                                                     // allocate_dpTfm,
	 int nTxTfm_v, double[] xTxTfm, double[] yTxTfm, double[] zTxTfm, double[] azimuthTxTfm, double[] dipTxTfm,               // allocateTfm 
	 float[] d_copy, float[] sd_copy,                                                                                         // allocateDcopies 
	 double[] sigsite, int [] RxInd, int[] rx1DtoRxIn,                                                                        // allocateDM
	 int nlay1D_v, double[] layer_params,                                                                                     // allocateL
	 int nRxRotSolve_v, int[] IndRxRotSolve                                                                                   // allocateIndRxRotSolve
	);

  public static native void setDipoleModPointers(int nlay1D, double[] zlay1D, double[] sig1D);  // allocateL
  public static native void  releaseDipoleModPointers(double[] zlay1D, double[] sig1D);
	

  public static native void computeFwd(boolean bDoPartials, float[] pm, float[] dm, float[] g);


  /** Load JNI library, using explicit full path name for library */
  static 
  {
	String path;
	
	path = getScrippsOneDLibraryPath();
	try
	  {
		// System.loadLibrary("ScrippsModelInterface");
		System.load(path);
	  }
	catch (SecurityException se) 
	  {
		haveNativeLibrary=false;
		System.err.println("Native fortran interface not accessible (security exception)");
		System.exit(1);
	  }
	catch (UnsatisfiedLinkError le) 
	  {
		haveNativeLibrary=false;
		System.err.println("\nError: "+le.toString()+"\n");
		//System.err.println("Message: " +le.getMessage());
		//le.printStackTrace();
		System.err.println("\nNative fortran interface not available. Try rebuilding.");
		System.exit(2);
	  }
  }

  /*
	Method to return full native library path of the correct
	platform-specific libName.* library.  On 64bit Linux platforms
	("Linux-x86_64-unknown") running 32bit jvms, the path will be the
	32bit build in "Linux-i686-unknown"

	Typical full library names are:
	ROOT/au/csiro/JamesGunning/CSEM/OneDinversion/Linux-i686-unknown/libScrippsModelInterface.so
  */
  public static String getScrippsOneDLibraryPath()
  {
	String scrippsPath;

	// Resolve absolute name to this class
	try
	  {
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		scrippsPath = loader.getResource(ScrippsOneDRelativeDir).toURI().getPath();
	  }
	catch (Exception e) 
	  {
		System.err.println("ClassLoader-error: "+e);
		scrippsPath="";
	  }
	
	scrippsPath = scrippsPath + File.separatorChar + JNIPathUtil.getJNILibraryPath(libName,null);

	return scrippsPath;
  }


  /**
	 test routine
  */
  public static final void test(int i)
  {
	if (haveNativeLibrary)
	  {
		int n=20000;
		float[] y = new float[n];
		y[0]=2.0f; y[1]=3.0f;
		System.out.println("Java inserted values "+y[0]+" "+y[1]);
		//testNative(i,y);
		System.out.println("New values seen in java "+y[0]+" "+y[1]);
	  }
  }

  /** 
	means to interrogate library type
  */
  public static boolean haveNativeLib()
  {
	return haveNativeLibrary;
  }

  /** 
	  Dummy testbed for method.
  */
  public static void main(String[] args) 
  {
	  System.err.println("Have native library = "+haveNativeLib());

	  ScrippsModelInterface.test(5);
  }


  /* Instance members */


  // ptr to calling inverter class
  CSEM1DInverter _CI;
  
  // References to module objects in calling inverter class
  ScrippsModelInterface.Occam_interface_mod OM;
  ScrippsModelInterface.Csem1d_mod CM;
  ScrippsModelInterface.Dipole1d_io DM;
  //  ScrippsModelInterface.Mt1d_mod MM;


  public ScrippsModelInterface(CSEM1DInverter CI)
  {
	this._CI = CI;		
	OM = new ScrippsModelInterface.Occam_interface_mod();
	CM = new ScrippsModelInterface.Csem1d_mod();
	DM = new ScrippsModelInterface.Dipole1d_io();
	//	MM = new ScrippsModelInterface.Mt1d_mod();
  }

 
	

  /**
	 Translations of Kerry Key's parsing code.
	 Reading of input files. Allocations etc.
  */
  void readData()
  {
	// Version 2.1  September 24, 2008.  Fixed a bug where iRxTfm was undefined.
	// Version 2.0  September 2008.    Added some error checks. Created Tx-Rx transformations 
	//                                 to speed up computations when simple transmitters are modeled.
	//                                 The speed of Dipole1D depends on the number of transmitters, so 
	//                                 transforming many simple transmitters into a single transmitter 
	//                                 can greatly speed up the 1D inversion.
	//                                 If only the x,y position of the transmitter varies,
	//                                 but not the azimuth, dip or z position, then this can be modeled as a single
	//                                 transmitter with many pseudo-receivers added to the Rx array.  There are even
	//                                 further transformations that can be added for variable azimuth Tx's, but 
	//                                 I'll skip that for now.
	//                                   
	// Version 1.0  Spring 2008.     
	//
	//
	  
    int          i, err, iTxRead, iRxRead, iFreqRead, iDataRead, nAllocErr, nRxRots=0, iRx, icnt, iTx,code=0;
	double       rxtheta,rxalpha,rxbeta;
	String sLine, sCode, sValue,s; // These are for reading the lines of the file
    boolean      bComment, lerror, lTxSame;
	Integer I;
    
    BufferedReader in=null;
	LineSplitter PT = new LineSplitter();
	LineSplitter PT1 = new LineSplitter();
	String[] st;
	
	ErrorHandler EH = _CI.EH;

	//
	// First step is to initialize Dipole1D defaults.  These can then be overridden
	// by values in the model file:
	// 
    DM.init_defaults_Dipole1D();

    CM.nTxIn  = 0;
    CM.nRxIn  = 0;
    CM.nRxRotSolve = 0;

 	try 
	  {
		in = new BufferedReader(new FileReader(_CI.datfil));
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't open file "+_CI.datfil);
	  }
	
	//
	// Loop through the model file looking for the Format and Number of Layers flags
	// and skipping over any comments
	//
    for(;;)
	  {
        // Get the next code/value pair
        // ParseCode forces the code portion to be all lowercase with no
        //   padding and ending colon stripped off.  User comments are 
        //   stripped from the value portion.        
        try
		  {
			if ((sLine = in.readLine()) == null) break;
		  }
		catch (IOException e)
		  {
			break; //skip reading
		  }
		PT.parse(sLine);
		if (PT.bComment) continue;
		sCode = PT.sCode;
		sValue = PT.sValue;

		// What do we have?
        if (sCode.equals("format"))
		  {
			// Set cnDataParams for the input data format        
            sValue = sValue.toLowerCase(); 
            if (sValue.equals("emdata_1.1"))
			  { 
				CM.cnDataParams = 4;
			  } 
			else
			  { 
				EH.die_err("Error: data format unsupported: "+sValue+" Try using EMDATA_1.1");
			  }
		  }
		else if (sCode.equals("# transmitters"))
		  {
			try
			  {
				CM.nTxIn = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Number of transmitters not a positive integer");
			  }
            //  Allocate arrays dependent on nTx
            CM.allocateTX(CM.nTxIn);
                
            //
            // Now read in block of transmitters, skipping any comment lines:
            //
            for (iTxRead=0; iTxRead < CM.nTxIn;) 
			  { 
                // Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient transmitters");
				  }
        
                // Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
                // Else it is not a comment, so read in the Tx parameters from sLine
                st = sLine.trim().split("\\s+");            
                if (st.length >= 5)
				  {
					try
					  {
						CM.xTxIn[iTxRead] = Util.parseDouble(st[0]);
						CM.yTxIn[iTxRead] = Util.parseDouble(st[1]);
						CM.zTxIn[iTxRead] = Util.parseDouble(st[2]);
						CM.azimuthTxIn[iTxRead] = Util.parseDouble(st[3]);
						CM.dipTxIn[iTxRead] = Util.parseDouble(st[4]);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("X Y Z Azimuth Dip contains a non-number:"+sLine);
					  }
				  }
				else
				  {
					EH.die_err("Not 5 number entries on transmitter line "+sLine);
				  }
				iTxRead++;   
			  }
            //
            // Check and make sure all the transmitters were read in:
            //
            if (iTxRead != CM.nTxIn) 
			  {
				EH.die_err("Error reading transmitters: iTxRead!=nTxIn: "+iTxRead+" "+CM.nTxIn);
			  }
		  }        
		else if (Util.matches(sCode,"# frequencies","nfreq","#freq","#freqs", "# freq"))
		  {
			// Get the next code/value pair
			try
			  {
				CM.nFreqIn = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Indecipherable # frequencies - not an integer");
			  }
            CM.allocateF(CM.nFreqIn);

            //
            // Now read in block of frequencies, skipping any comment lines:
            //
            for (iFreqRead = 0; iFreqRead < CM.nFreqIn;) 
			  {
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient frequencies");
				  }
				// Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;

                // Else it is not a comment, so read in the Tx parameters from sLine
				try   
				  {
					CM.fTxIn[iFreqRead] = Util.parseDouble(sLine);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Frequency not a number");
				  }
                iFreqRead++;
			  }        
            //
            // Check and make sure all the frequencies were read in:
            //
            if (iFreqRead != CM.nFreqIn) 
			  {
				EH.die_err(" Error reading frequencies: iFreqRead!=nFreqIn: "+iFreqRead+" "+CM.nFreqIn);
			  }
		  }
		else if (Util.matches(sCode,"# receivers","nrx","#receivers"))
		  {
			try  
			  {
				CM.nRxIn = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("# receivers not an integer");
			  }

            // Allocate the field arrays: 
            //
            CM.allocateR(CM.nRxIn);
            // initialize rotations to 0 degrees (auto in java)
            // Now read in block of receiver locations, skipping any comment lines:
            //
			for (iRxRead = 0;iRxRead < CM.nRxIn;) 
			  {
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient receivers");
				  }
        
                // Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
                // Else it is not a comment, so read in the Tx parameters from sLine
                st = sLine.trim().split("\\s+");            
                if (st.length >= 3)
				  {
					try
					  {
						CM.xRxIn[iRxRead] = Util.parseDouble(st[0]);
						CM.yRxIn[iRxRead] = Util.parseDouble(st[1]);
						CM.zRxIn[iRxRead] = Util.parseDouble(st[2]);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err(sLine+": Not numbers for receivers");
					  }
				  }
				else
				  {
					EH.die_err("Not 3 columns of numbers for receivers");
				  }
				iRxRead++;   
			  }
            //
            // Check and make sure all the layers were read in:
            //
            if (iRxRead!=CM.nRxIn)
			  {
				EH.die_err("Error reading receivers: iRxRead != nRxIn: "+iRxRead+" "+CM.nRxIn);
			  }
			
		  }
		else if (sCode.equals("# receiver orientations"))
		  {
			//
			// Optional parameter for receiver rotations.  This allows for modeling data
			// at the orientation of the receiver, which can be different from x,y,z.
			//  So this is actually the angle the model responses are rotated by. Degrees clockwise from x toward y.
			//
			try
			  {
				nRxRots = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Indecipherable # receiver orientations - not an integer");
			  }
            for(iRxRead = 0; iRxRead < nRxRots;) 
			  {
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient receiver orientations");
				  }        
                // Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
                // Else it is not a comment, so read in the Tx parameters from sLine
                st = sLine.trim().split("\\s+");            
                if (st.length >= 4)
				  {
					try
					  {
						int irx = Integer.parseInt(st[0]);
						CM.ThetaRxIn[irx-1] = Util.parseDouble(st[1]);
						CM.AlphaRxIn[irx-1] = Util.parseDouble(st[2]);
						CM.BetaRxIn[irx-1] = Util.parseDouble(st[3]);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Not numbers for receiver orientations");
					  }
				  }
				else
				  {
					EH.die_err("Not 4 columns of numbers for receiver orientations");
				  }
				iRxRead++;   
			  }
		  }
		else if (Util.matches(sCode,"# solve orientations","# unknown orientations","# unknown rotations"))
		  {
			//
			// Should we solve for Rx orientation from x?
			//
			try  
			  {
				CM.nRxRotSolve = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("# Unknown Rotations: not an integer");
			  }
            // Allocate;
            CM.allocateIndRxRotSolve(CM.nRxRotSolve);

			for (iRxRead = 0; iRxRead < CM.nRxIn;) 
			  {
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient unknown rotations");
				  }
				
                // Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
				try  
				  {
					CM.IndRxRotSolve[iRxRead] = Integer.parseInt(sLine);
				  }
				catch (NumberFormatException e)
				  {
					EH.die_err("Unknown rotation not an integer");
				  }
				iRxRead++;
			  }

            // These are stored in pm_assoc  (associated parameters)
            OM.npm_assoc = CM.nRxRotSolve*3; // solving for 3 angles to get full 3D orientation (tilts set to zero if only horiz data).
            OM.pm_assoc = new float[OM.npm_assoc]; 
		  }
		else if (Util.matches(sCode,"# data","#data","ndata"))
		  {
			//
			// Lastly we need to read in the DATA block.
			//
			try  
			  {
				OM.nd = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Number of data not an integer");	
			  }
			
            // Allocate the nd data arrays: 
            
            OM.allocateD(OM.nd,CM.cnDataParams);
            CM.allocateDcopies(OM.nd);
         
			// read data. Types may be string codes, so get HashMap which may be needed
			CM.setUpCodeMap();

            for(iDataRead = 0; iDataRead < OM.nd;) 
			  {
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient data");
				  }        
                // Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
                // Else it is not a comment, so read in the Tx parameters from sLine
                st = sLine.trim().split("\\s+");            
                if (st.length >= 6)
				  {
					try
					  {
						s = st[0];
						if (Util.isInteger(s))
						  {
							code = Integer.parseInt(s);
						  }
						else
						  {
							// a string code
							if ((I = CM.codeMap.get(s.toLowerCase()))==null) 
							  {
								EH.die_err("Unrecognizable type-code "+s+" in data section, data-line "+(iDataRead+1));
							  }
							else code = I.intValue();
						  }
						OM.dp_trans.setQuick(0, iDataRead, code);
						OM.dp_trans.setQuick(1, iDataRead, Integer.parseInt(st[1]));
						OM.dp_trans.setQuick(2, iDataRead, Integer.parseInt(st[2]));
						OM.dp_trans.setQuick(3, iDataRead, Integer.parseInt(st[3]));
						OM.d[iDataRead] = Util.parseFloat(st[4]);
						OM.sd[iDataRead] = Util.parseFloat(st[5]);
					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Not numbers in data section, data-line "+(iDataRead+1));
					  }
				  }
				else
				  {
					EH.die_err("Not 6 columns of numbers for data");
				  }
				iDataRead++;
			  }
			
            // Check and make sure all the data was read in:
            if (iDataRead != OM.nd)
			  {
                EH.die_err("Error reading receivers: iDataRead!=nd: "+iDataRead+ " "+OM.nd);
			  }
            //
            // Copy data array for OPRA:
            //
            CM.d_copy = (float[])OM.d.clone();
            CM.sd_copy = (float[])OM.sd.clone();
		  }
	  } // The main data reading loop
	try
	  {
		in.close();
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't close " + _CI.datfil);
	  }

	//
	// Check that parameters have been defined:
	//
	lerror = false;
	if (CM.nTxIn < 1)
	  {
    	lerror = true;
    	EH.println(" Error in Data File, insufficient number of transmitters: nTx = "+CM.nTxIn);
	  }
   	if (CM.nRxIn < 1)
	  {
    	lerror = true;
    	EH.println(" Error in Data File, insufficient number of receivers: nRx = "+CM.nRxIn);
	}
   	if (CM.nFreqIn < 1)
	  {
    	lerror = true;
    	EH.println(" Error in Data File, insufficient number of frequencies: nFreq = "+CM.nFreqIn);
	  }
   	if (OM.nd < 1)
	  {
    	lerror = true;
    	EH.println(" Error in Data File, insufficient number of data: # Data = "+OM.nd);
	  }	
	
	if (lerror) EH.die_err("Giving up. Fix data input!");
	
	//------------------------------------------------------------------------------
	// Create Rx and Tx arrays for fast computations.  This means reducing the 
	// number of discrete transmitters since Dipole1D's runtime is proportional 
	// to the number of transmitters.  So, we will look through the array of 
	// transmitters. If all transmitters have the same zTx,AzimuthTx and DipTx 
	// then we can simply reduce them to a single transmitter and 
	// many sites (in x and y).
	//------------------------------------------------------------------------------
	
	//
	// Allocate a dpTfm array that contains transformed data parameter indices:
	//
	
	CM.allocate_dpTfm(OM.nd);

	// Loop through nTx  and check for identical transmitters parameters:
	// 
	
	lTxSame = true;
	double MATCH_EPS = 0.01;
	for (iTx = 1; iTx<CM.nTxIn;iTx++)
	  {
	 	if ( Math.abs( CM.zTxIn[iTx]       - CM.zTxIn[0]       ) > MATCH_EPS ) lTxSame = false;  
		if ( Math.abs( CM.azimuthTxIn[iTx] - CM.azimuthTxIn[0] ) > MATCH_EPS ) lTxSame = false;
	 	if ( Math.abs( CM.dipTxIn[iTx]     - CM.dipTxIn[0]     ) > MATCH_EPS ) lTxSame = false;
	  }

	//
	// If not similar, then just copy default Rx and Tx params:
	
	// COMMENT -> restore this when translation done
	//  if (!lTxSame) 
	if (true) 
	  {
		if (!lTxSame)
		  {
			EH.println("FIXME ---finish the translations in the code at this remark to get maximum forward kernel speed...");
		  }
		else
		  {
			EH.println("Warning: not using Tx-Rx transformations yet: code translation *incomplete*. Accuracy is OK, but speed"
					   +" slower for multiple transmitters. A to-do.");
		  }
		CM.allocateTfm(CM.nTxIn);
		
		CM.nTxTfm       = CM.nTxIn;
		CM.xTxTfm       = (double[])CM.xTxIn.clone();

		CM.yTxTfm       = (double[])CM.yTxIn.clone();
		CM.zTxTfm       = (double[])CM.zTxIn.clone();
		CM.azimuthTxTfm = (double[])CM.azimuthTxIn.clone();
		CM.dipTxTfm     = (double[])CM.dipTxIn.clone();
		
		CM.allocateRfm(CM.nRxIn);
		
		CM.nRxTfm = CM.nRxIn;
		CM.xRxTfm = CM.xRxIn;
		CM.yRxTfm = CM.yRxIn;
		CM.zRxTfm = CM.zRxIn;
		// JG next line has to be fortran indexing....
		for (i=0; i<CM.nRxIn; i++) CM.iRxTfm[i] = (i+1);

		for (i=0; i<OM.nd; i++)
		  {
			CM.dpTfm_trans.setQuick(0,i, OM.dp_trans.getQuick(2,i)); // [i][2]);
			CM.dpTfm_trans.setQuick(1,i, OM.dp_trans.getQuick(3,i)); // [i][3]);
		  }
 	  }
	//
	// If Tx's are similar, then create transformation arrays:
	//	
	else
	  {
		EH.die_err("FIXME ---finish this translation of the else bit...");
		/*
	
       // Define a single transmitter:
		nTxTfm = 1
		allocate ( xTxTfm(1), yTxTfm(1), zTxTfm(1), azimuthTxTfm(1), dipTxTfm(1), &
		         & stat = nAllocErr)
		if (nAllocErr .ne. 0) then
			EH.println("Out of memory.  Too many transmitters  (', 1, ')'
			stop 
		endif   		
	   
		xTxTfm       = 0d0
		yTxTfm       = 0d0
		zTxTfm       = zTxIn(1)
		azimuthTxTfm = azimuthTxIn(1)
		dipTxTfm     = dipTxIn(1)
		
	  // Define array of sites:
	  	nRxTfm = nTxIn*nRxIn
		allocate ( xRxTfm(nRxTfm), yRxTfm(nRxTfm), zRxTfm(nRxTfm),iRxTfm(nRxTfm), stat = nAllocErr)
		
		if (nAllocErr .ne. 0) then
			EH.println("Out of memory.  Too many receivers  (', nRxIn, ')'
			stop 
		endif   	  
	
		icnt = 0	
		do iRx = 1,nRxIn
		   do iTx = 1,nTxIn
		       icnt = icnt + 1
				xRxTfm(icnt) = xRxIn(iRx) - xTxIn(iTx)
				yRxTfm(icnt) = yRxIn(iRx) - yTxIn(iTx)
				zRxTfm(icnt) = zRxIn(iRx)
				iRxTfm(icnt) = iRx // store the original site number for the Tfm site! 
				// JG UNDER TRANSLATION check offsetting of iRxTfm!
				
		   enddo
		enddo
	    
	    // Now update the dpTfmarray to point to transformed transmitter and sites:
	    dpTfm(1:nd,1) = 1  // only 1 transformed transmitter 
	    
		do i = 1,nd  // 
		   iTx = dp(i,3)		
		   iRx = dp(i,4)
		   dpTfm(i,2) = iTx + (iRx-1)*nTxIn  // in other words transformed receiver index is now in blocks of iTx
		enddo
		        
		write(*,*)  ' '
		EH.println(" Transmitters are similar, so I transformed the transmitters into '
		EH.println(" many receivers to speed up the computations. '
		EH.println(" Now the number of effective receivers is: ',icnt 
		*/
	  }
	
	
	//   
	// Allocate the arrays passed in dipole1D_io module.  Note that derivative arrays are allocated in readModel.
	//
	DM.allocateXYZDM(CM.nRxTfm);
	CM.allocateDM(CM.nRxTfm);

  }     

  /**
	 Translation of KKs model reader for OCCAM
  */
  void readModel()
  {              
	ErrorHandler EH = _CI.EH;

	// References to module objects in calling inverter class
    // ScrippsModelInterface.Csem1d_mod CM; // stores all the variables needed to go between Occam and Dipole1D 
    // ScrippsModelInterface.Occam_interface_mod OM = CI.OM;   
    // ScrippsModelInterface.Mt1d_mod MM = CI.MM;   // used for passing 1D MT parameters and responses

	// ScrippsModelInterface.Dipole1d_io DM = CI.DM;    
    // Model nLay1D and zlay1D set, allocates sig1D(nlay1D), zlay1D(nlay1D)
	// note that sig1D is not set since it needs to be updated to 
	// the current model on each call to Dipole1D (done in subr. Fwd_Calc)
                     
	//
	// Layer block in file is
	// 
	// [top_depth    resistivity  penalty    Preference   prej_penalty]
	// 
	// Use 0 for Preference to not have a Preferenced parameter
	// Penalty is penalty to apply at the layer top boundary.  Use 0 for no roughness penalty at layer top
	//

    int     i=0, err, iof, ilayread, ichecknParams, nAllocErr;
    String  sLine, sCode, sValue; // These are for reading the lines of the file
    boolean bComment, lhasPref;
	boolean foundFirstFree=false;

	BufferedReader in=null;
	LineSplitter PT = new LineSplitter();
	LineSplitter PT1 = new LineSplitter();
	String[] st;

	// Initialize number of layers in model:

    DM.nlay1D = 0;
	try 
	  {
		in = new BufferedReader(new FileReader(_CI.modfil));
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't open file "+_CI.modfil);
	  }
    EH.println ("Reading model from file: "+_CI.modfil);
	//
	// Loop through the model file looking for the Format and Number of Layers flags
	// and skipping over any comments
	//
    for(;;)
	  {
        // Get the next code/value pair
        // ParseCode forces the code portion to be all lowercase with no
        //   padding and ending colon stripped off.  User comments are 
        //   stripped from the value portion.
        try
		  {
			if ((sLine = in.readLine()) == null) break;
		  }
		catch (IOException e)
		  {
			break; //skip reading
		  }

        PT.parse(sLine);
		if (PT.bComment) continue;
		sCode = PT.sCode;
		sValue = PT.sValue;
        
        // What do we have?
		if (sCode.equals("format"))
		  {
			sValue = sValue.toLowerCase();
            if (!(sValue.equals("occam1dmod_1.0")))
			  {
                EH.die_err("Error: model format unsupported: "+sValue+
						   "Try using Occam1DMod_1.0");
			  }
		  }
		else if (Util.matches(sCode,"# layers","#layers", "number of layers","nlayers"))   
		  {
            // Get the number of layers:
			try
			  {
				DM.nlay1D = Integer.parseInt(sValue);
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Number of layers not a positive integer");
			  }                                    
			CM.allocateL(DM.nlay1D);
			DM.allocateL(DM.nlay1D);
            // read in the layering blocks (skipping any comment lines)
            for (ilayread = 0;ilayread < DM.nlay1D;) 
			  {
                // read in  a line:
				// Get the next code/value pair
				try
				  {
					if ((sLine = in.readLine()) == null) break;
				  }
				catch (IOException e)
				  {
					EH.die_err("Insufficient layer specs, line "+(ilayread+1));
				  }
				// Check if line is a comment, if not strip off comments at end of line:
                PT1.parse(sLine);
                if (PT1.bComment ) continue;
                // Not a comment, so read in the values from sLine:  
				// Expected format
				// [top_depth,rho,penalty,Preference,prej_penalty]    
				st = sLine.trim().split("\\s+");            
                if (st.length >= 5)
				  {
					try
					  {
						for (i=0;i<5;i++) CM.layer_params_trans.setQuick(i,ilayread, Util.parseDouble(st[i]));

						if (CM.layer_params_trans.getQuick(1,ilayread) < 0)
						  {
							if (!foundFirstFree) CM.ZtopFree = CM.layer_params_trans.getQuick(0,ilayread); // save ZtopFree
							foundFirstFree = true;
						  }

					  }
					catch (NumberFormatException e)
					  {
						EH.die_err("Layer-parameter "+st[i]+" not a number");
					  }
				  }
				else
				  {
					EH.die_err("Layer-parameters line does not have 5 numbers");
				  }
				ilayread++;
			  }
            //
            // Check and make sure all the layers were read in:
            //
            if (ilayread != DM.nlay1D)
			  {
                EH.die_err("Error: model ilayread != DM.nlay1D: "+ilayread+" "+DM.nlay1D);
			  }
		  }
		else if (sCode.equals("ht filters"))
		  // The rest of the cases are for overriding the Dipole1D default setup
		  { 
            DM.HTmethod1D = sValue.toCharArray();
			EH.println("Using Hankel TF Filters: "+DM.HTmethod1D);
		  }
		else
		  {            
            EH.println("Error reading Model file, unknown or unsupported code: "+sCode);
		  }
	  }
	try
	  {
		in.close();
	  }
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't close " + _CI.modfil);
	  }
	//
	// Check to make sure layers are defined: 
	//
    if (DM.nlay1D <= 0)
	  {
        EH.die_err("Error reading model file: number of layers undefined");
	  }
	// 
	//  JG: Count the number of free parameters to see if it agrees with the startup iteration file
	//  LAYER_CHECK
    ichecknParams = 0;

	for (i = 0; i<DM.nlay1D; i++)
	  {
		if (CM.layer_params_trans.getQuick(1,i) < 0) 
		  {
			// negative resistivity is free parameter
            ichecknParams++;
		  }
	  }
	setNFreeSeaLayers();
	_CI.nFixedSurfaceLayers = getNFixedShallowLayers();
	if (_CI.nFixedSurfaceLayers != 2)
	  {
        EH.warn("Number of fixed shallow layers != 2. Code not fully debugged for this case yet => not trustworthy/faulty?");
        EH.warn("See also code occurences of LAYER_CHECK");
	  }

	_CI.np = ichecknParams;
	OM.nParams = _CI.np;
  }

  int getNFixedShallowLayers()
  {
	int nFixedShallowLayers=0;

	for (int i = 0; i<DM.nlay1D; i++)
	  {
		if (CM.layer_params_trans.getQuick(1,i) < 0) 
		  {
			// negative resistivity is free parameter. Break out
			break;
		  }
		else 
		  {
			// fixed
			nFixedShallowLayers++;
		  }
	  }
	return nFixedShallowLayers; 
  }

  void setNFreeSeaLayers()
  {
	boolean inSea=false,inAir=true;
	double flag;

	_CI.nFreeSeaLayers=0;
	for (int i = 0; i<DM.nlay1D; i++)
	  {
		// flag these in input files by number <-1 (e.g. -2)
		flag = CM.layer_params_trans.getQuick(1,i);
		if (flag < -1) 
		  {
			if (inAir)
			  {
				// entered sea layer with a -2 flag
				inAir = false;
				inSea = true;
			  }
			if (inSea) _CI.nFreeSeaLayers++;
		  }
		else if (flag >= -1f)
		  {
			// entered earth
			if (inSea) inSea = false;
		  }
	  }
  }

  void copyFixedUpperLayers()
  {
	CM.Nfixed = getNFixedShallowLayers();
	CM.fixed_upper_layer_params_trans = new double[5][CM.Nfixed];
	for(int i=0;i<CM.Nfixed;i++)
	  {
		for(int j=0;j<5;j++) CM.fixed_upper_layer_params_trans[j][i] = CM.layer_params_trans.getQuick(j,i);
	  }
  }

  void restoreFixedUpperLayers()
  {
	for(int i=0;i<CM.Nfixed;i++)
	  {
		for(int j=0;j<5;j++) CM.layer_params_trans.setQuick(j,i,CM.fixed_upper_layer_params_trans[j][i]);
	  }
  }

  void finishModelSetup()
  {
	int i;

	//
	// Assign layer top depths:
	//
	for (i=0; i<DM.nlay1D; i++) DM.zlay1D[i] = CM.layer_params_trans.getQuick(0,i); // [i][0];
    
	//
	// Notes: upon each call to dipole1d, sig1d will be filled using layer_params(i,2) or 
	// pm where layer_params(i,2) < 0 (ie free parameter)
	//
   
	//
	// Set Preference and weights:
	//
    int ichecknParams = 0;

    Arrays.fill(OM.premod,0); // initialize these to 0 (unused)
    Arrays.fill(OM.prewts,0);
    boolean lhasPref = false;
	double x;
	for (i=0;i<DM.nlay1D;i++)
	  {
        if (CM.layer_params_trans.getQuick(1,i) < 0)
		  {
			// negative resistivity is free parameter
            if (Math.abs(x = CM.layer_params_trans.getQuick(3,i)) != 0) 
			  {
				// convert to log10 scale resistivity (only iteration file has log10 resistivity)
                OM.premod[ichecknParams] = (float)Math.log10(x);
			  }
			else
			  {
				OM.premod[ichecknParams] = 0;
			  }
            
            OM.prewts[ichecknParams] = (float)CM.layer_params_trans.getQuick(4,i); 
            if ( OM.prewts[ichecknParams] != 0 ) 
			  {
				lhasPref = true;
			  }
			ichecknParams++;
		  }
	  }
    if (!lhasPref) _CI.EH.println("No model preference values used.");
	
	//
	// Allocate the jacobian arrays io from dipole1d:
	// (native as complex)
	doJacobianAllocs(CM.nRxTfm,DM.nlay1D);

	/*
	allocate( dexdsig(nRxTfm,DM.nlay1D),deydsig(nRxTfm,DM.nlay1D),djzdsig(nRxTfm,DM.nlay1D), &
			& dbxdsig(nRxTfm,DM.nlay1D),dbydsig(nRxTfm,DM.nlay1D),dbzdsig(nRxTfm,DM.nlay1D), &
			&  stat=nAllocErr )             
	*/
  }
    

  void reconcileDataVariables()
  {
	_CI.d = OM.d;
	_CI.sd = OM.sd;
	_CI.dm = OM.dm;
  }

  void reconcileGeoModelVariables()
  {
	int i,j,nfree;

	_CI.nlayer = DM.nlay1D - 2;  
	// Legacy Scripps: nlayer= #subsea layers 
	// New Scripps: even if third layer is not start of seabed, 
	// need space for _CI.nlayer+1 = (DM.nlay1D -1) variables anyway - do not use/store air
	//	_CI.nlayer = DM.nlay1D - _CI.nFixedSurfaceLayers - _CI.nFreeSeaLayers;


	_CI.allocateGeoModel();

	// legacy format, from water down : thick      rho ipen1    ppmod      ipen2
	// new format, from air down      : top_depth  rho penalty  Preference Prej_penalty  

	// this gives array of thickness from second layer down. Often, full sea-thickness
	for (i=0; i < _CI.nlayer ; i++) 
	  {
		_CI.thickf[i] = (float)(CM.layer_params_trans.getQuick(0,i+2)
								- CM.layer_params_trans.getQuick(0,i+1));
	  }
	// final thickness is assumed infinity anyway, so
	_CI.thickf[_CI.nlayer] = 500;



	// skip uppermost layer (probably air)
	// j= 0...nlayer,   i=1....(nlayer+1)=1...DM.nlay1D-1
	for (j=0,i=1; j<=_CI.nlayer; i++,j++)
	  {
		_CI.rho[j] = (float)CM.layer_params_trans.getQuick(1,i);
		_CI.ipen1[j] = (float)CM.layer_params_trans.getQuick(2,i);
		_CI.ipen2[j] = (float)CM.layer_params_trans.getQuick(4,i);
	  }
  }

  void reconcileModelVariables()
  {
	int i,nfree;

	// set up matching iff array
	// i.e. sort out fixed versus free resistivities:
	// This loop skips air layer (i=0 starts at second layer of input files) 
	// as I'll presume this is always fixed
	nfree=0;
	for (i=0; i<=_CI.nlayer; i++)
	  {
        if (_CI.rho[i] < 0)
		  {
			// negative resistivity is free parameter
			_CI.iff[nfree] = i;
			_CI.premod[nfree] =  (float)CM.layer_params_trans.getQuick(3,i+1);
			nfree++;
		  }
	  }
  }


  // Module members represented as classes
  public static class Occam_interface_mod
  {
	// (1) module Occam_interface_mod 
	//
	// This module contains objects that are used in the Forward
	// and sensitivity computations, and within Occam.
	//
	// *This is the only Occam module that Forward routines need to use.
    
	int nCurrentIter;     // Current Occam iteration number
	char cRootName;    // [50] The Occam iteration file root name. 
	// Defaults to 'ITER' if not given
	// as a command line argument
  
	// Data and model file names:
	char[] dataFileName, modelFileName; // [50]
  
	// Data to be inverted:
	int nd;
	float[] d, sd;
	DenseIntMatrix2D dp_trans;    // wrapper for int dp[][];
	// nd       = number of data.  + 1 or 2 if have statics
	// d(nd)    = data 
	// sd(nd)   = data errors (one standard error) - usually read from data file
	// dp(nd,:) = array of parameters associated with each datum (freq., etc)
	//            Note that dp() is not used in Occam, only in the forward routines
	//
	//-----------------------------------------------------------------------
	// model parameters for the current iteration and its response
	int nParams;
	float[] pm;
	float[] dm; 
	DenseFloatMatrix2D wj_trans;     // wrapper for float[][] wj;
	//   nParams        = number of free model parameters
	//   pm(nParams)    = array of model parameters
	//   dm(nd)         = array of model response (data domain)
	//   wj(nd,nParams) = model response Jacobian matrix
	//-----------------------------------------------------------------------
	// associated data/model parameters for the current iteration 
	// these are things that are derived from the current model and/or data
	// but are not used for the model roughness or misfit computations.
	// These are only along for the ride in the Occam routines and are entirely
	// dealt with in the Forward code.  They need to be included here
	// since Occam can keep track of pm_assoc for the best fitting model
	// of each iteration.
	int npm_assoc;
	float[] pm_assoc;  // associated model responses (not inverted for but important to Fwd responses)
    
	//   prewts(nParams) = weights for model prejudice
	//   premod(nParams) = model prejudice
	float[] prewts, premod;

	void setOccamModelPointers(int cnDataParams)
	{  
	  ScrippsModelInterface.setOccamModPointers(nParams, nd, cnDataParams, dp_trans.elements, d, sd, dm);
	}
	void releaseOccamModelPointers()
	{
	  ScrippsModelInterface.releaseOccamModPointers(dp_trans.elements, d, sd, dm);
	}

	void allocateD(int nd, int cnDataParams)
	{
	  this.nd = nd;
	  dp_trans = new DenseIntMatrix2D(cnDataParams,nd); // int[nd][cnDataParams];
	  d = new float[nd];
	  sd = new float[nd];
	  dm = new float[nd];
	}



	void forwardModel(float pm[])
	{
	  ScrippsModelInterface.computeFwd(false,pm,dm,null);
	}

	void forwardModelPlusJac(float pm[])
	{
	  ScrippsModelInterface.computeFwd(true,pm,dm,wj_trans.elements);
	}
  }
 

  public static class Csem1d_mod
  {
	// (2) module csem1d_mod

	// Module to store various parameters for 1D CSEM inversion
	// This is stuff that is not used in Occam.f90 and not Dipole1D.f90
	// but is required for interfacing between them.
	//
	// Data Type Index Dictionary
	//
	// These are the codes used in the data file parameter table.
	// 
	// Format:  EMData_1.0
	//
	// CSEM Data:
	// 
	// Note that this format will be compatible with 2D and 3D CSEM modeling.

	//   CSEM Electric Data, Real/Imag:
    public static final int indRealEx = 1;
    public static final int indImagEx = 2;
    public static final int indRealEy = 3;
    public static final int indImagEy = 4;
    public static final int indRealEz = 5;
    public static final int indImagEz = 6;

    
	//   CSEM Magnetic Data, Real/Imag:
    public static final int indRealBx = 11;
    public static final int indImagBx = 12;
    public static final int indRealBy = 13;
    public static final int indImagBy = 14;
    public static final int indRealBz = 15;
    public static final int indImagBz = 16;
    
	//   CSEM Electric Data, Amp/Phs:
    public static final int indAmpEx = 21;
    public static final int indPhsEx = 22;
    public static final int indAmpEy = 23;
    public static final int indPhsEy = 24;
    public static final int indAmpEz = 25;
    public static final int indPhsEz = 26;

    
	//   CSEM Magnetic Data, Amp/Phs:
    public static final int indAmpBx = 31;
    public static final int indPhsBx = 32;
    public static final int indAmpBy = 33;
    public static final int indPhsBy = 34;
    public static final int indAmpBz = 35;
    public static final int indPhsBz = 36;
    



	//   CSEM Horizontal Polarization Ellipses:
    public static final int iPEmax = 41;
    public static final int iPEmin = 42;
    
    public static final int iPBmax = 43;
    public static final int iPBmin = 44;
    

	
	// To do: E/E, E/H transfer functions, etc

	//   
	// Magnetotelluric Data:  
	//
	// This format is compatible with 2D and 3D MT modeling codes.
	//
	//   Apparent Resistivity and Phase:
    public static final int indRhoZXX  = 101;   // not used for 1D;
    public static final int indPhsZXX  = 102;   // not used for 1D;
    public static final int indRhoZXY  = 103;   // ;
    public static final int indPhsZXY  = 104;   //  ;
    public static final int indRhoZYX  = 105;   //;
    public static final int indPhsZYX  = 106;   // ;
    public static final int indRhoZYY  = 107;   // not used for 1D;
    public static final int indPhsZYY  = 108;   // not used for 1D;
    
	//   Impedances: 
    public static final int indRealZXX  = 111;   // not used for 1D;
    public static final int indImagZXX  = 112;   // not used for 1D;
    public static final int indRealZXY  = 113;   // ;
    public static final int indImagZXY  = 114;   // ;
    public static final int indRealZYX  = 115;   //;
    public static final int indImagZYX  = 116;   // ;
    public static final int indRealZYY  = 117;   // not used for 1D;
    public static final int indImagZYY  = 118;   // not used for 1D 

	// HashMap of above types for range 1-36 at least , for parsing
	HashMap<String,Integer> codeMap=null;

	public static final double pi = Math.PI;
	public static final double rad2deg = 57.295779513082323;
	public static final double deg2rad = 0.017453292519943;

	// While the data parameter array dp is passed in/out of Occam, it is only used
	// in the forward routines.  So we set cnDataParams in the forward routines
	// and then allocate the dp array in readData.
	//        
	int cnDataParams;
 
	DenseDoubleMatrix2D layer_params_trans;	// wrapper for double[][] layer_params ;
	int Nfixed;
	double[][] fixed_upper_layer_params_trans=null; // copy of fixed part of above used in splitting etc
	double ZtopFree = Double.MAX_VALUE; // top depth of uppermost free layer

	// layer_params is [ top_depth, rho, penalty, Preference, prej_penalty ]          
	//
	// top_depth  Depth to the top of the layer (ignores first value since top of 1D model is infinity)
	// rho        Either a fixed value or else -1, the flag for an inversion free parameter.
	//              If the resistivity is fixed, then none of the remaining 3 layer parameters are used 

	// penalty        The roughness penalty to apply across the layer top boundary. Use zero to allow for a jump
	// Preference      A linear resistivity Preference value that you can use to push the inversion in a certain direction.  
	//                  If you don't want a Preference, set this and the prej_penalty to 0.
	// prej_penalty   The roughness penalty to apply to the difference of the current model parameter value and the Preference.
	//                  If you really want to force the inversion towards the Preference value, make prej_penalty >> penalty for the layer
	//

	//
	// Store a copy of the d and sd vectors to use for the OPRA method:
	//
	float [] d_copy, sd_copy;
       
	//
	// Data dependent Tx-Rx variables read from data file in subroutine readData:
	//
	int nTxIn, nRxIn, nFreqIn;
	double[]   azimuthTxIn,dipTxIn;
	double[]   xTxIn,yTxIn,zTxIn,xRxIn,yRxIn,zRxIn;
	double[]   ThetaRxIn,AlphaRxIn, BetaRxIn;
	double[]   fTxIn, sigsite;


	//
	// Transformed Tx-Rx arrays used to speed up computations (reduces nTx when dip, azimuth and depth are constant): 
	//
	int nTxTfm;
	double[]     xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm;
	int          nRxTfm;
	double[]     xRxTfm, yRxTfm, zRxTfm;
	int[]    iRxTfm; // index pointing from Tfm site number to input site number (RxIn)
	int[]    RxInd;  // index pointing from Tfm site number to site sent to Dipole1D
	int[]    rx1DtoRxIn; // index pointing from Dipole1D site number to RxIn number, used for rotations
	DenseIntMatrix2D dpTfm_trans;	// int[][]  dpTfm;  // data parameter array for transformed Tx to Rx's
    
	//
	// Rx orientation solver variables:
	//
	int nRxRotSolve; // Number of Rx to solve for orientations
	int[] IndRxRotSolve; // [nRxRotSolve] array of Rx indexes to solve for rotations
	
	
	void allocateTX(int nTxIn)
	{
	  azimuthTxIn = new double[nTxIn];
	  dipTxIn = new double[nTxIn];
	  xTxIn= new double[nTxIn];
	  yTxIn= new double[nTxIn];
	  zTxIn= new double[nTxIn];
	}

	void allocateF(int nFreqIn)
	{
	  fTxIn = new double[nFreqIn];
	}

	void allocateR(int nRxIn)
	{
	  xRxIn = new double[nRxIn];
	  yRxIn = new double[nRxIn];
	  zRxIn = new double[nRxIn];
	  ThetaRxIn = new double[nRxIn];
	  AlphaRxIn= new double[nRxIn];
	  BetaRxIn  = new double[nRxIn];	  
	}

    void allocateRfm(int n)
	{
	  xRxTfm = new double[n];
	  yRxTfm = new double[n];
	  zRxTfm = new double[n];
	  iRxTfm = new int[n];
	}

	void allocate_dpTfm(int nd)
	{
	  dpTfm_trans = new DenseIntMatrix2D(2,nd); // new int[nd][2];
	}

	void allocateTfm(int n)
	{
	  xTxTfm = new double[n]; 
	  yTxTfm = new double[n]; 
	  zTxTfm = new double[n]; 
	  azimuthTxTfm = new double[n]; 
	  dipTxTfm = new double[n];
	}

	void allocateDcopies(int n)
	{
	  d_copy = new float[n];
	  sd_copy = new float[n];
	}
	
	void allocateDM(int n)
	{
	  sigsite = new double[n];
	  RxInd = new int[n];
	  rx1DtoRxIn = new int[n];
	}

    void allocateL(int nlay1D)
	{
	  layer_params_trans = new DenseDoubleMatrix2D(5,nlay1D); // double[nlay1D][5];	  
	}

	void allocateIndRxRotSolve(int nRxRotSolve)
	{
	  IndRxRotSolve = new int[nRxRotSolve]; 
	} 


	void setCSEMModelPointers(int nd, int nlay1D)
	{
	  ScrippsModelInterface.setCsemModPointers(
							  nTxIn, azimuthTxIn, dipTxIn, xTxIn, yTxIn, zTxIn,             // allocateTX
							  nFreqIn  ,fTxIn,                                              // allocateF  
							  nRxIn  , xRxIn, yRxIn, zRxIn, ThetaRxIn, AlphaRxIn, BetaRxIn, // allocateR,
							  nRxTfm  , xRxTfm, yRxTfm, zRxTfm, iRxTfm,                     // allocateRfm
							  nd, dpTfm_trans.elements,                                     // allocate_dpTfm,
							  nTxTfm  , xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm,     // allocateTfm 
							  d_copy, sd_copy,                                              // allocateDcopies 
							  sigsite, RxInd, rx1DtoRxIn,                                   // allocateDM
							  nlay1D, layer_params_trans.getElements(),                     // allocateL
							  nRxRotSolve  , IndRxRotSolve                                  // allocateIndRxRotSolve
							  );
	}
	
	void releaseCSEMModelPointers(int nd, int nlay1D)
	{
	  ScrippsModelInterface.releaseCsemModPointers(
							  nTxIn, azimuthTxIn, dipTxIn, xTxIn, yTxIn, zTxIn,             // allocateTX
							  nFreqIn  ,fTxIn,                                              // allocateF  
							  nRxIn  , xRxIn, yRxIn, zRxIn, ThetaRxIn, AlphaRxIn, BetaRxIn, // allocateR,
							  nRxTfm  , xRxTfm, yRxTfm, zRxTfm, iRxTfm,                     // allocateRfm
							  nd, dpTfm_trans.elements,                                     // allocate_dpTfm,
							  nTxTfm  , xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm,     // allocateTfm 
							  d_copy, sd_copy,                                              // allocateDcopies 
							  sigsite, RxInd, rx1DtoRxIn,                                   // allocateDM
							  nlay1D, layer_params_trans.getElements(),                     // allocateL
							  nRxRotSolve  , IndRxRotSolve                                  // allocateIndRxRotSolve
							  );
	}

	// HashMap of subset of EMDATA_1.1 types, from 1 to 36
	HashMap setUpCodeMap()
	{
	  codeMap = new HashMap<String,Integer>();

	  codeMap.put("realex" ,Integer.valueOf(1)); 
	  codeMap.put("imagex" ,Integer.valueOf(2)); 
	  codeMap.put("realey" ,Integer.valueOf(3)); 
	  codeMap.put("imagey" ,Integer.valueOf(4)); 
	  codeMap.put("realez" ,Integer.valueOf(5)); 
	  codeMap.put("imagez",  Integer.valueOf(6)); 
	  codeMap.put("realbx" ,Integer.valueOf(11));
	  codeMap.put("imagbx" ,Integer.valueOf(12));
	  codeMap.put("realby" ,Integer.valueOf(13));
	  codeMap.put("imagby" ,Integer.valueOf(14));
	  codeMap.put("realbz" ,Integer.valueOf(15));
	  codeMap.put("imagbz" ,Integer.valueOf(16));
	  codeMap.put("ampex" ,Integer.valueOf(21)); 
	  codeMap.put("phsex" ,Integer.valueOf(22)); 
	  codeMap.put("ampey" ,Integer.valueOf(23)); 
	  codeMap.put("phsey" ,Integer.valueOf(24)); 
	  codeMap.put("ampez" ,Integer.valueOf(25)); 
	  codeMap.put("phsez" ,Integer.valueOf(26)); 
	  codeMap.put("ampbx" ,Integer.valueOf(31)); 
	  codeMap.put("phsbx" ,Integer.valueOf(32)); 
	  codeMap.put("ampby" ,Integer.valueOf(33)); 
	  codeMap.put("phsby" ,Integer.valueOf(34)); 
	  codeMap.put("ampbz" ,Integer.valueOf(35)); 
	  codeMap.put("phsbz" ,Integer.valueOf(36)); 

	  return codeMap;
	}
  }

  public static class Dipole1d_io
  {
	// Transmitter paramters
	//
    double     sdm1D;             // (Am)      Source dipole moment 
    double     ftx1D;             // (Hz)      Source frequency 
    double     azimuthTx1D  ;     // (degrees) Source azimuth from x axis (positive clockwise)
    double     dipTx1D;           // (degrees) Vertical dip angle of source along azimuthTx1D, positive down 
    double     xTx1D,yTx1D,zTx1D; // (m)   Transmitter x, y and z location   
	//
	// Model parameters
	//
    int       nlay1D;    // Number of layers
    double[]  sig1D;     // (S/m) Layer conductivities 
    double[]  zlay1D;    // (m)   Depth to top of each layer, first layer ignored 
	//
	// Site locations (where E and B are computed)
	//
    int       n1D;          // number of receiver sites
    double[]  x1D, y1D, z1D;
	//
	// Option to only compute electric fields.  lbcomp = .false. or .true. 
	//
    boolean lbcomp;  
	//
	// Fast Hankel Transform Filter method.  Choose the filter coefficients:
	//
    char[] HTmethod1D = new char[32];   // Filter options listed in FilterModules.f90

	//
	// Output domain (either spatial (x,y,z) or (kx,y,z) wavenumber domain for 2.5D modeling primary fields)
	//
    char[] outputdomain1D = new char[32];     // Use 'spatial' or 'kx'
	//
	// If using outputdomain1D='kx', then set kx1d, kxmode1D and CTmethod1D:
	//
	// kxmode1D is used to define the principle axis for the Tx
	// 
    double kx1D;     // If outputdomain1D='kx', then the fields are 
	// computed at at this wavenumber.  
                              
    int kxmode1D; // 1,2,3 for x,y,z dipoles.  Ignored if spatial computation.
	//
	// (Co)Sine transform using digital filter method. 
	//
    char[] CTmethod1D = new char[32];  // Filter options listed in FilterModules.f90

	//
	// Use Spline interpolation of coefficients for massive speed up when there are many receivers
	//
	boolean lUseSpline1D;
    
	//
	//  Electric and magnetic field vectors.  Note that I output vertical current:
	//
	// complex(8), dimension(:), allocatable  ex1D,ey1D,jz1D,bx1D,by1D,bz1D  // units are V/m and T, e(-iwt) convention used

	//
	// Inversion sensitivity computations:
	//

	boolean linversion;  // set =.true. to output the field derivatives with respect to conductivity
	//
	//  Electric and magnetic field derivatives with respect to layer conductivities
	//  This are computed when linversion = .true.
	// 
	// complex(8),dimension(:,:), allocatable  dexdsig,deydsig,djzdsig,dbxdsig,dbydsig,dbzdsig   // df_i / dsig_j array

	

	void setDipoleModelPointers()
	{
	  ScrippsModelInterface.setDipoleModPointers(nlay1D,zlay1D,sig1D);
	}

	void releaseDipoleModelPointers()
	{
	  ScrippsModelInterface.releaseDipoleModPointers(zlay1D,sig1D);
	}



	void allocateXYZDM(int n)
	{
	  /*
		JG - move these to fortran only needed internally
		x1D = new double[n];
		y1D = new double[n];
		z1D = new double[n];
	  */
	  ScrippsModelInterface.doEBFieldAllocs(n);
	}
	
	void allocateL(int nlay1D)
	{
	  zlay1D = new double[nlay1D];
	  sig1D = new double[nlay1D];
	}

	void init_defaults_Dipole1D()
	{
	  // Specify some defaults for parameters required by Dipole1D:
	  /* default to fortran kernel for now */
	  /*
		HTmethod1D = (new String("kk_ht_201")).toCharArray();    // Use 201 point HT digital filters.
		outputdomain1D  = (new String("spatial")).toCharArray(); // Assume spatial domain comps
		lbcomp          = false;        // This is changed to true if magnetics in data file
		sdm1D           = 1.0;            // (Am), dipole moment. Normalize to unit source moment
		lUseSpline1D    = true;         // Use spline interpolation for faster 1D computations
		linversion      = false;        // Compute derivatives with respect to sigma(layers)
	  */
	  ScrippsModelInterface.initDefaultsDipole1D();
	}
  }

}
