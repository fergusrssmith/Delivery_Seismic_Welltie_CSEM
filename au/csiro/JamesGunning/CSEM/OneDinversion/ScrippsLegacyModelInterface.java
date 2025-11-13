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

import au.csiro.JamesGunning.lib.util.*; 

import java.io.*;
import java.nio.*;
import java.util.Arrays;

import au.csiro.JamesGunning.lib.Numerical.DenseFloatMatrix2D;

/**
   Class to provide JNI interface to legacy fortran
   1D CSEM modelling routines
   originally from the Scripps "Occam" fortran codes.

   Simple test for loading...

   java au.csiro.JamesGunning.CSEM.OneDinversion.ScrippsLegacyModelInterface
*/

public final class ScrippsLegacyModelInterface 
{
  private static boolean haveNativeLibrary = true;

  public static native void setModelPointersNative(int np, int nd, int nlayer, int nfreq,
												   float[] thickf_ptr, float[] rho_ptr, 
												   float[] freqs_ptr, float[] ipen1_ptr,
												   float[] ipen2_ptr, int[] iff_ptr);

  public static native void releaseModelPointersNative(float[] thickf_ptr, float[] rho_ptr, 
													   float[] freqs_ptr, float[] ipen1_ptr,
													   float[] ipen2_ptr, int[] iff_ptr);
  // short versions of above
  public static native void setSimpleModelPointersNative(int np_iff, int nlayer, float[] thickf_ptr, float[] rho_ptr,int[] iff_ptr);

  public static native void releaseSimpleModelPointersNative(float[] thickf_ptr, float[] rho_ptr,int[] iff_ptr);


  public static native void formodNative(int np, int nd, float[] pm, float[] dp,float[] dm,float[] d);

  public static native void fordivNative(int np, int nd, float[] pm, float[] dp, float[] dm, float[] wj_trans,float[] d);

  public static native int testmeNative(int i,float[] y);

  public static boolean TEST_OPTIMISER = false;

  /**
	 Location of native library builds wrt root of java packages.
	 This location contains builds in the directories<p>

	 ScrippsOneDRelativeDir/Linux-i686-unknown/libScrippsLegacyModelInterface.so,<p>
	 ScrippsOneDRelativeDir/Linux-x86_64-unknown/libScrippsLegacyModelInterface.so,<p>
	 ScrippsOneDRelativeDir/Darwin-Power-Macintosh-unknown/libScrippsLegacyModelInterface.jnilib<p>

	 etc etc
  */
  public static final String ScrippsOneDRelativeDir = "au/csiro/JamesGunning/CSEM/OneDinversion";
  /*
	name of library
  */
  public static final String libName = "libScrippsLegacyModelInterface";


  /** Load JNI library, using explicit full path name for library */
  static 
  {
	String path;
	
	path = getScrippsOneDLibraryPath();
	try
	  {
		// System.loadLibrary("ScrippsLegacyModelInterface");
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
	ROOT/au/csiro/JamesGunning/CSEM/OneDinversion/Linux-i686-unknown/libScrippsLegacyModelInterface.so
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

  public static final void setModelPointers(int np, int nd, int nlayer, int nfreq,
											float[] thickf_ptr, float[] rho_ptr, 
											float[] freqs_ptr, float[] ipen1_ptr,float[] ipen2_ptr, 
											int[] iff_ptr)
  {
	setModelPointersNative(np, nd, nlayer, nfreq, thickf_ptr,rho_ptr,freqs_ptr,ipen1_ptr,ipen2_ptr,iff_ptr);
  }

  public static final void releaseModelPointers(float[] thickf_ptr, float[] rho_ptr, 
												float[] freqs_ptr, float[] ipen1_ptr,float[] ipen2_ptr, 
												int[] iff_ptr)
  {
	releaseModelPointersNative(thickf_ptr,rho_ptr,freqs_ptr,ipen1_ptr,ipen2_ptr,iff_ptr);
  }

  public static final void setSimpleModelPointers(int np_iff, int nlayer,
												  float[] thickf_ptr, float[] rho_ptr,int[] iff_ptr)
  {
	setSimpleModelPointersNative(np_iff,nlayer,thickf_ptr,rho_ptr,iff_ptr);
  }

  public static final void releaseSimpleModelPointers(float[] thickf_ptr, float[] rho_ptr,int[] iff_ptr)
  {
	releaseSimpleModelPointersNative(thickf_ptr,rho_ptr,iff_ptr);
  }


  public static final void formod(int np, int nd, float[] pm, float[] dp,float[] dm,float[] d)
  {
	if (!TEST_OPTIMISER) formodNative(np,nd,pm,dp,dm,d); 
	else formodTestOptimiser(np,nd,pm,dp,dm,d);
  }

  public static final void fordiv(int np, int nd, float[] pm, float[] dp, float[] dm, DenseFloatMatrix2D wj_trans,float[] d)
  {
  	if (!TEST_OPTIMISER) fordivNative(np, nd, pm, dp, dm, wj_trans.elements,d);
	else fordivTestOptimiser(np, nd, pm, dp, dm, wj_trans,d);
  }


  /**
	 test routine
  */
  public static final void testme(int i)
  {
	if (haveNativeLibrary)
	  {
		int n=20000;
		float[] y = new float[n];
		y[0]=2.0f; y[1]=3.0f;
		System.out.println("Java inserted values "+y[0]+" "+y[1]);
		testmeNative(i,y);
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

	  ScrippsLegacyModelInterface.testme(5);
  }


  /** 
	  Dummy methods to test optimization framework
  */
  public static void makeTestOptimiserCalibrationDataSet(CSEM1DInverter inverter)
  {
	int i,np_g;

	inverter.EH.println("WARNING! TEST LINEAR DATA AND MODEL!");

	np_g = Math.min(inverter.np,inverter.nd); // almost certainly np
	for (i=0; i<np_g; i++)
	  {
		inverter.d[i] = (i<=15)? -0.5f+0.1f*i : 0f;
		inverter.sd[i] = 0.05f;
	  }
	for (i = np_g; i<inverter.nd ; i++)
	  {
		inverter.d[i] = 0.05f;
		inverter.sd[i] = 0.05f;
	  }
  }

  public static void formodTestOptimiser(int np, int nd, float[] pm, float[] dp,float[] dm,float[] d)
  {
	int i,np_g;
	np_g = Math.min(np,nd);

	for (i=0; i<np_g; i++)
	  {
		dm[i] = pm[i];
	  }
	for (i = np_g; i<nd ; i++) dm[i] = 0f; //
  }

  public static final void fordivTestOptimiser(int np, int nd, float[] pm, float[] dp, float[] dm, DenseFloatMatrix2D wj_trans,float[] d)
  {
	int i,j,np_g;

	formodTestOptimiser(np,nd,pm,dp,dm,d);

	np_g = Math.min(np,nd);
	for (i=0; i<np; i++)
	  {
		for (j=0;j<nd;j++) wj_trans.setQuick(i,j,(i==j && i<np_g) ? 1f : 0f);
	  }
  }
}

