/*************************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 

#include "au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface.h"
#include "au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_Csem1d_mod.h"
#include "au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_Dipole1d_io.h"
#include "au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_Occam_interface_mod.h"


#include "platform_specifics.h"

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>

/**
USE_PRIMITIVE_ARRAY_CRITICAL: If true, we will use
GetPrimitiveArrayCritical() type calls, which try hard not to copy
data. For this to be safe, native methods must not force deadlocks etc
by doing thread stuff. See

http://java.sun.com/j2se/1.5.0/docs/guide/jni/spec/functions.html#wp17314

A safe method is invoked if this macro is disabled (USE_PRIMITIVE_ARRAY_CRITICAL=0)
- then we use GetFloatArrayElements() type methods instead.
*/

#define USE_PRIMITIVE_ARRAY_CRITICAL 1

#if USE_PRIMITIVE_ARRAY_CRITICAL
  #define GET_FLOATS     GetPrimitiveArrayCritical
  #define RELEASE_FLOATS ReleasePrimitiveArrayCritical
#else
  #define GET_FLOATS     GetFloatArrayElements
  #define RELEASE_FLOATS ReleaseFloatArrayElements
#endif

#define DEBUG 0




typedef struct 
{
// SH_SP
   double*  azimuthTxIn_ptr;
   double*  dipTxIn_ptr;
   double*  xTxIn_ptr;
   double*  yTxIn_ptr;
   double*  zTxIn_ptr;
   double*  fTxIn_ptr;
   double*    xRxIn_ptr;
   double*    yRxIn_ptr;
   double*    zRxIn_ptr;
   double*    ThetaRxIn_ptr;
   double*    AlphaRxIn_ptr;
   double*    BetaRxIn_ptr;
   double*  xRxTfm_ptr;
   double*  yRxTfm_ptr;
   double*  zRxTfm_ptr;
   jint*    iRxTfm_ptr;
   jint*  dpTfm_ptr;
   double*    xTxTfm_ptr;
   double*    yTxTfm_ptr;
   double*    zTxTfm_ptr;
   double*    azimuthTxTfm_ptr;
   double*    dipTxTfm_ptr;
   float*    d_copy_ptr;
   float*    sd_copy_ptr;
   double*   sigsite_ptr;
   jint* RxInd_ptr;
   jint* rx1DtoRxIn_ptr;
   double*  layer_params_ptr; //  double[nlay1D][5];
   jint* IndRxRotSolve_ptr;
  // SH_EP
} NS_MODEL_POINTERS;
static THREAD NS_MODEL_POINTERS cptrs;

typedef struct 
{
  double* zlay1D_ptr; //  double[nlay1D]
  double* sig1D_ptr; //  double[nlay1D]
} D_MODEL_POINTERS;
static THREAD D_MODEL_POINTERS dptrs;

typedef struct 
{
  jint *dp_ptr;
  jfloat *d_ptr;
  jfloat *sd_ptr;
  jfloat *dm_ptr;
} CS_MODEL_POINTERS;
static THREAD CS_MODEL_POINTERS csptrs;


// Here, we have wrappers to the fortran library members. Check underscores via nm if needed */

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_doEBFieldAllocs
(JNIEnv *env, jclass obj, jint n)
{
  do_eb_field_allocs_(&n);
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_doJacobianAllocs
(JNIEnv *env, jclass obj, jint n, jint m)
{
  do_jacobian_allocs_(&n,&m);
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_freeJacobianAllocs
(JNIEnv *env, jclass obj)
{
  free_jacobian_allocs_();
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_initDefaultsDipole1D
(JNIEnv *env, jclass obj)
{
  set_defaults_dipole1d_();
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_setOccamModPointers
(JNIEnv *env, jclass obj, jint np, jint n, jint cnDataParams, jintArray dp, jfloatArray d, jfloatArray sd, jfloatArray dm)
{
  jboolean isCopy=JNI_FALSE;

  csptrs.dp_ptr=NULL;
  csptrs.d_ptr=NULL;
  csptrs.sd_ptr=NULL;
  csptrs.dm_ptr=NULL;

  csptrs.dp_ptr = (*env)->GetIntArrayElements(env,dp, &isCopy);
  csptrs.d_ptr = (*env)->GetFloatArrayElements(env,d, &isCopy);
  csptrs.sd_ptr = (*env)->GetFloatArrayElements(env,sd, &isCopy);
  csptrs.dm_ptr = (*env)->GetFloatArrayElements(env,dm, &isCopy);

  // Fortran call
  set_occam_mod_pointers_(&np,&n,&cnDataParams,csptrs.dp_ptr, csptrs.d_ptr, csptrs.sd_ptr, csptrs.dm_ptr);
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_releaseOccamModPointers
(JNIEnv *env, jclass obj, jintArray dp, jfloatArray d, jfloatArray sd, jfloatArray dm)
{
 if (csptrs.dp_ptr)  (*env)->ReleaseIntArrayElements(env,dp, csptrs.dp_ptr,0);
 if (csptrs.d_ptr)  (*env)->ReleaseFloatArrayElements(env,d, csptrs.d_ptr,0);
 if (csptrs.sd_ptr)  (*env)->ReleaseFloatArrayElements(env,sd, csptrs.sd_ptr,0);
 if (csptrs.dm_ptr)  (*env)->ReleaseFloatArrayElements(env,dm, csptrs.dm_ptr, 0);
}



// SH_ST
JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_setCsemModPointers
(env,obj,
 nTxIn_v, azimuthTxIn, dipTxIn, xTxIn, yTxIn, zTxIn,           // allocateTX
 nFreqIn_v,fTxIn,                                              // allocateF  
 nRxIn_v, xRxIn, yRxIn, zRxIn, ThetaRxIn, AlphaRxIn, BetaRxIn, // allocateR,
 nRxTfm_v, xRxTfm, yRxTfm, zRxTfm, iRxTfm,                     // allocateRfm
 nd, dpTfm,                                                    // allocate_dpTfm,
 nTxTfm_v, xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm,     // allocateTfm 
 d_copy, sd_copy,                                              // allocateDcopies 
 sigsite, RxInd, rx1DtoRxIn,                                   // allocateDM
 nlay1D_v, layer_params,                                       // allocateL
 nRxRotSolve_v,IndRxRotSolve                                   // allocateIndRxRotSolve
 )
JNIEnv *env; 
jclass obj;

// allocateTX
jint  nTxIn_v;
jdoubleArray  azimuthTxIn;
jdoubleArray  dipTxIn;
jdoubleArray  xTxIn;
jdoubleArray  yTxIn;
jdoubleArray  zTxIn;
// allocateF  
jint  nFreqIn_v;
jdoubleArray  fTxIn;
// allocateR
jint  nRxIn_v;
jdoubleArray    xRxIn;
jdoubleArray    yRxIn;
jdoubleArray    zRxIn;
jdoubleArray    ThetaRxIn;
jdoubleArray    AlphaRxIn;
jdoubleArray    BetaRxIn;
// allocateRfm
jint  nRxTfm_v;
jdoubleArray  xRxTfm;
jdoubleArray  yRxTfm;
jdoubleArray  zRxTfm;
jintArray  iRxTfm;
//allocate_dpTfm
jint  nd;
jintArray  dpTfm;
// allocateTfm;
jint  nTxTfm_v;
jdoubleArray    xTxTfm;
jdoubleArray    yTxTfm;
jdoubleArray    zTxTfm;
jdoubleArray    azimuthTxTfm;
jdoubleArray    dipTxTfm;
// allocateDcopies
jfloatArray    d_copy;
jfloatArray    sd_copy;
// allocateDM
jdoubleArray   sigsite;
jintArray   RxInd;
jintArray  rx1DtoRxIn;
//allocateL
jint  nlay1D_v;
jdoubleArray  layer_params; // double[nlay1D][5];
//allocateIndRxRotSolve
jint  nRxRotSolve_v;
jintArray  IndRxRotSolve;
{
  // SH_V
  jboolean isCopy=JNI_FALSE;

  // allocateTX
  //jint  nTxIn_v;
  cptrs.azimuthTxIn_ptr = NULL;
  cptrs.dipTxIn_ptr = NULL;
  cptrs.xTxIn_ptr = NULL;
  cptrs.yTxIn_ptr = NULL;
  cptrs.zTxIn_ptr = NULL;
  // allocateF  
  //jint  nFreqIn_v;
  cptrs.fTxIn_ptr = NULL;
  // allocateR
  //jint  nRxIn_v;
  cptrs.xRxIn_ptr = NULL;
  cptrs.yRxIn_ptr = NULL;
  cptrs.zRxIn_ptr = NULL;
  cptrs.ThetaRxIn_ptr = NULL;
  cptrs.AlphaRxIn_ptr = NULL;
  cptrs.BetaRxIn_ptr = NULL;
  // allocateRfm
  //jint  nRxTfm_v;
  cptrs.xRxTfm_ptr = NULL;
  cptrs.yRxTfm_ptr = NULL;
  cptrs.zRxTfm_ptr = NULL;
  cptrs.iRxTfm_ptr = NULL;
  //allocate_dpTfm
  //jint  nd;
  cptrs.dpTfm_ptr = NULL;
  // allocateTfm;
  //jint  nTxTfm_v;
  cptrs.xTxTfm_ptr = NULL;
  cptrs.yTxTfm_ptr = NULL;
  cptrs.zTxTfm_ptr = NULL;
  cptrs.azimuthTxTfm_ptr = NULL;
  cptrs.dipTxTfm_ptr = NULL;
  // allocateDcopies
  cptrs.d_copy_ptr = NULL;
  cptrs.sd_copy_ptr = NULL;
  // allocateDM
  cptrs.sigsite_ptr = NULL;
  cptrs.RxInd_ptr = NULL;
  cptrs.rx1DtoRxIn_ptr = NULL;
  //allocateL
  //jint  nlay1D_v;
  cptrs.layer_params_ptr = NULL; //  double[nlay1D][5];
  //double* zlay1D_ptr = NULL; //  double[nlay1D]
  //double* sig1D_ptr  = NULL; //  double[nlay1D]
  //allocateIndRxRotSolve
  //jint  nRxRotSolve_v;
  cptrs.IndRxRotSolve_ptr = NULL;


  // SH_EV

  // allocateTX
  if (nTxIn_v > 0)
	{
	  cptrs.azimuthTxIn_ptr = (*env)->GetDoubleArrayElements(env,  azimuthTxIn, &isCopy);
	  cptrs.dipTxIn_ptr = (*env)->GetDoubleArrayElements(env,  dipTxIn, &isCopy);
	  cptrs.xTxIn_ptr = (*env)->GetDoubleArrayElements(env,  xTxIn, &isCopy);
	  cptrs.yTxIn_ptr = (*env)->GetDoubleArrayElements(env,  yTxIn, &isCopy);
	  cptrs.zTxIn_ptr = (*env)->GetDoubleArrayElements(env,  zTxIn, &isCopy);
	}


  // allocateF  
  //jint  
  if (nFreqIn_v > 0)
	{
	  cptrs.fTxIn_ptr = (*env)->GetDoubleArrayElements(env,  fTxIn, &isCopy);
	}


  // allocateR
  //jint  nRxIn_v;
  if (nRxIn_v > 0)
	{
	  cptrs.xRxIn_ptr = (*env)->GetDoubleArrayElements(env,    xRxIn, &isCopy);
	  cptrs.yRxIn_ptr = (*env)->GetDoubleArrayElements(env,    yRxIn, &isCopy);
	  cptrs.zRxIn_ptr = (*env)->GetDoubleArrayElements(env,    zRxIn, &isCopy);
	  cptrs.ThetaRxIn_ptr = (*env)->GetDoubleArrayElements(env,    ThetaRxIn, &isCopy);
	  cptrs.AlphaRxIn_ptr = (*env)->GetDoubleArrayElements(env,    AlphaRxIn, &isCopy);
	  cptrs.BetaRxIn_ptr = (*env)->GetDoubleArrayElements(env,    BetaRxIn, &isCopy);
	}


  // allocateRfm
  //jint  nRxTfm_v;
  if (nRxTfm_v > 0)
	{
	  cptrs.xRxTfm_ptr = (*env)->GetDoubleArrayElements(env,  xRxTfm, &isCopy);
	  cptrs.yRxTfm_ptr = (*env)->GetDoubleArrayElements(env,  yRxTfm, &isCopy);
	  cptrs.zRxTfm_ptr = (*env)->GetDoubleArrayElements(env,  zRxTfm, &isCopy);
	  cptrs.iRxTfm_ptr = (*env)->GetIntArrayElements(env,    iRxTfm, &isCopy);
	}


  // allocate_dpTfm
  //jint  nd;
  if (nd > 0)
	{
	  cptrs.dpTfm_ptr = (*env)->GetIntArrayElements(env,  dpTfm, &isCopy);
	}
 

  // allocateTfm;
  //jint  nTxTfm_v;
  if (nTxTfm_v > 0)
	{
	  cptrs.xTxTfm_ptr = (*env)->GetDoubleArrayElements(env,    xTxTfm, &isCopy);
	  cptrs.yTxTfm_ptr = (*env)->GetDoubleArrayElements(env,    yTxTfm, &isCopy);
	  cptrs.zTxTfm_ptr = (*env)->GetDoubleArrayElements(env,    zTxTfm, &isCopy);
	  cptrs.azimuthTxTfm_ptr = (*env)->GetDoubleArrayElements(env,    azimuthTxTfm, &isCopy);
	  cptrs.dipTxTfm_ptr = (*env)->GetDoubleArrayElements(env,    dipTxTfm, &isCopy);
	}

  // allocateDcopies
  if (nd > 0)
	{
	  cptrs.d_copy_ptr = (*env)->GetFloatArrayElements(env,    d_copy, &isCopy);
	  cptrs.sd_copy_ptr = (*env)->GetFloatArrayElements(env,    sd_copy, &isCopy);
	}


  // allocateDM
  if (nRxTfm_v > 0)
	{
	  cptrs.sigsite_ptr = (*env)->GetDoubleArrayElements(env,   sigsite, &isCopy);
	  cptrs.RxInd_ptr = (*env)->GetIntArrayElements(env, RxInd, &isCopy);
	  cptrs.rx1DtoRxIn_ptr = (*env)->GetIntArrayElements(env, rx1DtoRxIn, &isCopy);
	}

  // allocateL
  //jint  nlay1D_v;
  if (nlay1D_v > 0)
	{
	  cptrs.layer_params_ptr = (*env)->GetDoubleArrayElements(env,  layer_params, &isCopy); //  double[nlay1D][5];
	}

  // allocateIndRxRotSolve
  //jint  nRxRotSolve_v;
  if (nRxRotSolve_v>0) 
	{
	  cptrs.IndRxRotSolve_ptr = (*env)->GetIntArrayElements(env,  IndRxRotSolve, &isCopy);
	}

  // fortran routine
  set_csem1d_mod_pointers_( 
	   &nTxIn_v, cptrs.azimuthTxIn_ptr, cptrs.dipTxIn_ptr, cptrs.xTxIn_ptr, cptrs.yTxIn_ptr, cptrs.zTxIn_ptr,               // allocateTX
	   &nFreqIn_v,cptrs.fTxIn_ptr,                                                                  // allocateF  
	   &nRxIn_v, cptrs.xRxIn_ptr, cptrs.yRxIn_ptr, cptrs.zRxIn_ptr, cptrs.ThetaRxIn_ptr, cptrs.AlphaRxIn_ptr, cptrs.BetaRxIn_ptr, // allocateR,
	   &nRxTfm_v, cptrs.xRxTfm_ptr, cptrs.yRxTfm_ptr, cptrs.zRxTfm_ptr, cptrs.iRxTfm_ptr,                             // allocateRfm
	   &nd, cptrs.dpTfm_ptr,                                                                        // allocate_dpTfm,
	   &nTxTfm_v, cptrs.xTxTfm_ptr, cptrs.yTxTfm_ptr, cptrs.zTxTfm_ptr, cptrs.azimuthTxTfm_ptr, cptrs.dipTxTfm_ptr,         // allocateTfm 
	   cptrs.d_copy_ptr, cptrs.sd_copy_ptr,                                                               // allocateDcopies 
	   cptrs.sigsite_ptr, cptrs.RxInd_ptr, cptrs.rx1DtoRxIn_ptr,                                                // allocateDM
	   &nlay1D_v, cptrs.layer_params_ptr,                                                           // allocateL
	   &nRxRotSolve_v,cptrs.IndRxRotSolve_ptr                                                       // allocateIndRxRotSolve
						   );
}
// SH_EF

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_releaseCsemModPointers
(env,obj,
 nTxIn_v, azimuthTxIn, dipTxIn, xTxIn, yTxIn, zTxIn,           // allocateTX
 nFreqIn_v,fTxIn,                                              // allocateF  
 nRxIn_v, xRxIn, yRxIn, zRxIn, ThetaRxIn, AlphaRxIn, BetaRxIn, // allocateR,
 nRxTfm_v, xRxTfm, yRxTfm, zRxTfm, iRxTfm,                     // allocateRfm
 nd, dpTfm,                                                    // allocate_dpTfm,
 nTxTfm_v, xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm,     // allocateTfm 
 d_copy, sd_copy,                                              // allocateDcopies 
 sigsite, RxInd, rx1DtoRxIn,                                   // allocateDM
 nlay1D_v, layer_params,                                       // allocateL
 nRxRotSolve_v,IndRxRotSolve                                   // allocateIndRxRotSolve
 )
JNIEnv *env; 
jclass obj;

// allocateTX
jint  nTxIn_v;
jdoubleArray  azimuthTxIn;
jdoubleArray  dipTxIn;
jdoubleArray  xTxIn;
jdoubleArray  yTxIn;
jdoubleArray  zTxIn;
// allocateF  
jint  nFreqIn_v;
jdoubleArray  fTxIn;
// allocateR
jint  nRxIn_v;
jdoubleArray    xRxIn;
jdoubleArray    yRxIn;
jdoubleArray    zRxIn;
jdoubleArray    ThetaRxIn;
jdoubleArray    AlphaRxIn;
jdoubleArray    BetaRxIn;
// allocateRfm
jint  nRxTfm_v;
jdoubleArray  xRxTfm;
jdoubleArray  yRxTfm;
jdoubleArray  zRxTfm;
jintArray  iRxTfm;
//allocate_dpTfm
jint  nd;
jintArray  dpTfm;
// allocateTfm;
jint  nTxTfm_v;
jdoubleArray    xTxTfm;
jdoubleArray    yTxTfm;
jdoubleArray    zTxTfm;
jdoubleArray    azimuthTxTfm;
jdoubleArray    dipTxTfm;
// allocateDcopies
jfloatArray    d_copy;
jfloatArray    sd_copy;
// allocateDM
jdoubleArray   sigsite;
jintArray   RxInd;
jintArray  rx1DtoRxIn;
//allocateL
jint  nlay1D_v;
jdoubleArray  layer_params; // double[nlay1D][5];
//allocateIndRxRotSolve
jint  nRxRotSolve_v;
jintArray  IndRxRotSolve;
{
  /* release pointers */
  if (cptrs.azimuthTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, azimuthTxIn, cptrs.azimuthTxIn_ptr ,0);
  if (cptrs.dipTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, dipTxIn, cptrs.dipTxIn_ptr ,0);
  if (cptrs.xTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, xTxIn, cptrs.xTxIn_ptr ,0);
  if (cptrs.yTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, yTxIn, cptrs.yTxIn_ptr ,0);
  if (cptrs.zTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, zTxIn, cptrs.zTxIn_ptr ,0);
  if (cptrs.fTxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, fTxIn, cptrs.fTxIn_ptr ,0);
  if (cptrs.xRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, xRxIn, cptrs.xRxIn_ptr ,0);
  if (cptrs.yRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, yRxIn, cptrs.yRxIn_ptr ,0);
  if (cptrs.zRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, zRxIn, cptrs.zRxIn_ptr ,0);
  if (cptrs.ThetaRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, ThetaRxIn, cptrs.ThetaRxIn_ptr ,0);
  if (cptrs.AlphaRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, AlphaRxIn, cptrs.AlphaRxIn_ptr ,0);
  if (cptrs.BetaRxIn_ptr) (*env)->ReleaseDoubleArrayElements(env, BetaRxIn, cptrs.BetaRxIn_ptr ,0);
  if (cptrs.xRxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, xRxTfm, cptrs.xRxTfm_ptr ,0);
  if (cptrs.yRxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, yRxTfm, cptrs.yRxTfm_ptr ,0);
  if (cptrs.zRxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, zRxTfm, cptrs.zRxTfm_ptr ,0);
  if (cptrs.iRxTfm_ptr) (*env)->ReleaseIntArrayElements(env, iRxTfm, cptrs.iRxTfm_ptr ,0);
  if (cptrs.dpTfm_ptr) (*env)->ReleaseIntArrayElements(env, dpTfm, cptrs.dpTfm_ptr ,0);
  if (cptrs.xTxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, xTxTfm, cptrs.xTxTfm_ptr ,0);
  if (cptrs.yTxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, yTxTfm, cptrs.yTxTfm_ptr ,0);
  if (cptrs.zTxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, zTxTfm, cptrs.zTxTfm_ptr ,0);
  if (cptrs.azimuthTxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, azimuthTxTfm, cptrs.azimuthTxTfm_ptr ,0);
  if (cptrs.dipTxTfm_ptr) (*env)->ReleaseDoubleArrayElements(env, dipTxTfm, cptrs.dipTxTfm_ptr ,0);
  if (cptrs.d_copy_ptr) (*env)->ReleaseFloatArrayElements(env, d_copy, cptrs.d_copy_ptr ,0);
  if (cptrs.sd_copy_ptr) (*env)->ReleaseFloatArrayElements(env, sd_copy, cptrs.sd_copy_ptr ,0);
  if (cptrs.sigsite_ptr) (*env)->ReleaseDoubleArrayElements(env, sigsite, cptrs.sigsite_ptr ,0);
  if (cptrs.RxInd_ptr) (*env)->ReleaseIntArrayElements(env, RxInd, cptrs.RxInd_ptr ,0);
  if (cptrs.rx1DtoRxIn_ptr) (*env)->ReleaseIntArrayElements(env, rx1DtoRxIn, cptrs.rx1DtoRxIn_ptr ,0);
  if (cptrs.layer_params_ptr) (*env)->ReleaseDoubleArrayElements(env, layer_params, cptrs.layer_params_ptr ,0);
  if (cptrs.IndRxRotSolve_ptr) (*env)->ReleaseIntArrayElements(env, IndRxRotSolve, cptrs.IndRxRotSolve_ptr ,0);
}



JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_setDipoleModPointers
(
 JNIEnv *env, jclass obj, 
 jint nlay1D,
 jdoubleArray  zlay1D,        // double[nlay1D]
 jdoubleArray  sig1D          // double[nlay1D]
 )
{
  jboolean isCopy=JNI_FALSE;
  dptrs.zlay1D_ptr = NULL; //  double[nlay1D]
  dptrs.sig1D_ptr  = NULL; //  double[nlay1D]

  //allocateL
  if (nlay1D > 0)
	{
	  dptrs.zlay1D_ptr = (*env)->GetDoubleArrayElements(env,  zlay1D, &isCopy);             //  double[nlay1D];
	  dptrs.sig1D_ptr = (*env)->GetDoubleArrayElements(env, sig1D, &isCopy);                //  double[nlay1D];
	}
  set_dipole_mod_pointers_(dptrs.zlay1D_ptr, dptrs.sig1D_ptr);
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_releaseDipoleModPointers
(
 JNIEnv *env, jclass obj, 
 jdoubleArray  zlay1D,        // double[nlay1D]
 jdoubleArray  sig1D          // double[nlay1D]
)
{
  //allocateL
  if (dptrs.zlay1D_ptr) (*env)->ReleaseDoubleArrayElements(env, zlay1D, dptrs.zlay1D_ptr ,0);
  if (dptrs.sig1D_ptr) (*env)->ReleaseDoubleArrayElements(env, sig1D, dptrs.sig1D_ptr ,0);
}



JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_computeFwd
(JNIEnv *env, jclass obj, jboolean bDoPartials, jfloatArray pm, jfloatArray dm, jfloatArray g)
{
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;
  int getJac;
  float *pm_ptr=NULL, *dp_ptr=NULL, *dm_ptr=NULL,*d_ptr=NULL,*g_ptr=NULL;

  pm_ptr = (*env)->GET_FLOATS(env,pm, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif

  if (bDoPartials && (g!=NULL))
	{
	  g_ptr = (*env)->GET_FLOATS(env,g, &isCopy);
#if DEBUG
	  madeCopy |= isCopy;
#endif
	}

  dm_ptr = (*env)->GET_FLOATS(env,dm, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif

  // set dm[] and g[] (Jacobian) pointers

  if (g_ptr != NULL) 
	{
	  set_dm_wj_pointer_(dm_ptr, g_ptr);
	}
  else 
	set_dm_pointer_(dm_ptr);

  getJac = bDoPartials ? 1 : 0;
  computefwdwrapper_(&getJac,pm_ptr);


  /* release pointers */
  (*env)->RELEASE_FLOATS(env,pm,pm_ptr,0);

  if (g_ptr != NULL)
	{
	  (*env)->RELEASE_FLOATS(env,g,g_ptr,0);
	  /*
		printf("GRAD ");
		int i;
		for(i=0;i<100;i++) printf("%g ",g_ptr[i]);
	  */
	}
  (*env)->RELEASE_FLOATS(env,dm,dm_ptr,0);
}


#if 0

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_fordivNative
(JNIEnv *env, jclass obj, jint np, jint nd, jfloatArray pm, jfloatArray dp, jfloatArray dm, jfloatArray g, jfloatArray d)
{
  float *pm_ptr, *dp_ptr, *dm_ptr, *g_ptr,*d_ptr;
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;

  pm_ptr = (*env)->GET_FLOATS(env,pm, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif

  dm_ptr = (*env)->GET_FLOATS(env,dm, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif


  /* point at start of 2d arrays...contiguity assured by unrolled DenseFloatMatrix2D element[] */
  dp_ptr =  (*env)->GET_FLOATS(env,dp, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif

  g_ptr = (*env)->GET_FLOATS(env,g, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
  printf("GetPrimitiveArrayCriticalFOR: copy = %d\n",(madeCopy==JNI_TRUE)? 1:0);
#endif

  d_ptr = (*env)->GET_FLOATS(env,d, &isCopy);
#if DEBUG
  madeCopy |= isCopy;
#endif

  /* call fortran routine */
  fordiv_(&np, &nd, pm_ptr, dp_ptr, dm_ptr, g_ptr,d_ptr);

  /* release pointers */
  (*env)->RELEASE_FLOATS(env,pm,pm_ptr,0);
  (*env)->RELEASE_FLOATS(env,dm,dm_ptr,0);
  (*env)->RELEASE_FLOATS(env,dp,dp_ptr,0);
  (*env)->RELEASE_FLOATS(env,g,g_ptr,0);
  (*env)->RELEASE_FLOATS(env,g,d_ptr,0);
}
#endif


/*
  test method
*/

JNIEXPORT jint JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsModelInterface_testmeNative
(JNIEnv *env, jclass obj, jint x,jfloatArray y)
{
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;
  jfloat *yn=NULL;
  int i; 

  yn = (jfloat*) ((*env)->GetPrimitiveArrayCritical(env,y,&isCopy));

  yn[0]=4.0; 
  yn[1]=5.0;

  printf("Tests: should have %g=4.0 %g=5.0 isCopy=%d\n",yn[0],yn[1],isCopy);
  (*env)->ReleasePrimitiveArrayCritical(env,y,yn,0);

}
