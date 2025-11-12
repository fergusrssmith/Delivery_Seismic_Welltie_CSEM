#include <jni.h>
#include <stdio.h>
#include <math.h>

#include "au_csiro_JamesGunning_lib_util_MathUtil.h"

JNIEXPORT jdouble JNICALL Java_au_csiro_JamesGunning_lib_util_MathUtil_expNative(JNIEnv *e, jclass o, jdouble x)
{
  return exp(x);
}

JNIEXPORT jdouble JNICALL Java_au_csiro_JamesGunning_lib_util_MathUtil_logNative(JNIEnv *e, jclass o, jdouble x)
{
  return log(x);
}

JNIEXPORT jdouble JNICALL Java_au_csiro_JamesGunning_lib_util_MathUtil_sinNative(JNIEnv *e, jclass o, jdouble x)
{
  return sin(x);
}
JNIEXPORT jdouble JNICALL Java_au_csiro_JamesGunning_lib_util_MathUtil_cosNative(JNIEnv *e, jclass o, jdouble x)
{
  return cos(x);
}
JNIEXPORT jdouble JNICALL Java_au_csiro_JamesGunning_lib_util_MathUtil_tanNative(JNIEnv *e, jclass o, jdouble x)
{
  return tan(x);
}
