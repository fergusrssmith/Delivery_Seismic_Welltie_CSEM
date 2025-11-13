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

#include <jni.h>
#include "au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface.h"

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
  jfloat *thickf_ptr;
  jfloat *rho_ptr;
  jfloat *freqs_ptr;
  jfloat *ipen1_ptr;
  jfloat *ipen2_ptr;
  jint* iff_ptr; 
} MODEL_POINTERS;

static THREAD MODEL_POINTERS model_ptrs;

// Here, we have wrappers to the fortran library members. Check underscores via nm if needed */

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_formodNative
(JNIEnv *env, jclass obj, jint np, jint nd, jfloatArray pm, jfloatArray dp, jfloatArray dm, jfloatArray d)
{
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;
  float *pm_ptr, *dp_ptr, *dm_ptr,*d_ptr;

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
  printf("GetPrimitiveArrayCritical: copy = %d\n",(madeCopy==JNI_TRUE)? 1:0);
#endif

  d_ptr = (*env)->GET_FLOATS(env,d, &isCopy);
#if DEBUG
   madeCopy |= isCopy;
#endif

  // call fortran routine
  formod_(&np,&nd,pm_ptr,dp_ptr,dm_ptr,d_ptr);

  /* release pointers */
  (*env)->RELEASE_FLOATS(env,pm,pm_ptr,0);
  (*env)->RELEASE_FLOATS(env,dm,dm_ptr,0);
  (*env)->RELEASE_FLOATS(env,dp,dp_ptr,0);
  (*env)->RELEASE_FLOATS(env,dp,d_ptr,0);
}


JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_fordivNative
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
  (*env)->RELEASE_FLOATS(env,d,d_ptr,0);
}


JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_setModelPointersNative
  (JNIEnv *env, jclass obj, jint np, jint nd, jint nlayer, jint nfreq, 
   jfloatArray thickf, jfloatArray rho, jfloatArray freqs, 
   jfloatArray ipen1, jfloatArray ipen2, jintArray iff)
{
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;

  model_ptrs.thickf_ptr = (*env)->GetFloatArrayElements(env,thickf,&isCopy);
  model_ptrs.rho_ptr = (*env)->GetFloatArrayElements(env,rho, &isCopy);
  model_ptrs.freqs_ptr = (*env)->GetFloatArrayElements(env,freqs, &isCopy);
  model_ptrs.ipen1_ptr = (*env)->GetFloatArrayElements(env,ipen1, &isCopy);
  model_ptrs.ipen2_ptr = (*env)->GetFloatArrayElements(env,ipen2, &isCopy);
  model_ptrs.iff_ptr = (*env)->GetIntArrayElements(env,iff, &isCopy);


  // To work on *nix, 'doze: pass as standard ptrs in c, do offsetting in f90 set_module_pointers
  set_module_pointers_(&np, &nd, &nlayer,&nfreq, model_ptrs.thickf_ptr, 
					   model_ptrs.rho_ptr, model_ptrs.freqs_ptr, model_ptrs.ipen1_ptr,
					   model_ptrs.ipen2_ptr, model_ptrs.iff_ptr);
  
  /* release pointers by matching call as per below */
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_releaseModelPointersNative
 (JNIEnv *env, jclass obj,
   jfloatArray thickf, jfloatArray rho, jfloatArray freqs, 
   jfloatArray ipen1, jfloatArray ipen2, jintArray iff)
{
  /* release pointers */

  (*env)->ReleaseFloatArrayElements(env, thickf, model_ptrs.thickf_ptr, 0);
  (*env)->ReleaseFloatArrayElements(env, rho, model_ptrs.rho_ptr, 0);
  (*env)->ReleaseFloatArrayElements(env, freqs, model_ptrs.freqs_ptr, 0);
  (*env)->ReleaseFloatArrayElements(env, ipen1, model_ptrs.ipen1_ptr, 0);
  (*env)->ReleaseFloatArrayElements(env, ipen2, model_ptrs.ipen2_ptr, 0);
  (*env)->ReleaseIntArrayElements(env, iff, model_ptrs.iff_ptr, 0);
}


JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_setSimpleModelPointersNative
(JNIEnv *env, jclass obj, jint np, jint nlayer, 
   jfloatArray thickf, jfloatArray rho, jintArray iff)
{
  jboolean isCopy=JNI_FALSE, madeCopy=JNI_FALSE;

  model_ptrs.thickf_ptr = (*env)->GetFloatArrayElements(env,thickf,&isCopy);
  model_ptrs.rho_ptr = (*env)->GetFloatArrayElements(env,rho, &isCopy);
  model_ptrs.iff_ptr = (*env)->GetIntArrayElements(env,iff, &isCopy);

  // To work on *nix, 'doze: pass as standard ptrs in c, do offsetting in f90 set_module_pointers
  set_simple_module_pointers_(&np, &nlayer, model_ptrs.thickf_ptr, model_ptrs.rho_ptr, model_ptrs.iff_ptr);


  /* release pointers by matching call as per below */
}

JNIEXPORT void JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_releaseSimpleModelPointersNative
 (JNIEnv *env, jclass obj,
  jfloatArray thickf, jfloatArray rho, jintArray iff)
{
  /* release pointers */

  (*env)->ReleaseFloatArrayElements(env, thickf, model_ptrs.thickf_ptr, 0);
  (*env)->ReleaseFloatArrayElements(env, rho, model_ptrs.rho_ptr, 0);
  (*env)->ReleaseIntArrayElements(env, iff, model_ptrs.iff_ptr, 0);
  //   showtotaltiming_();
}


/*
  test method
*/

JNIEXPORT jint JNICALL Java_au_csiro_JamesGunning_CSEM_OneDinversion_ScrippsLegacyModelInterface_testmeNative
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
