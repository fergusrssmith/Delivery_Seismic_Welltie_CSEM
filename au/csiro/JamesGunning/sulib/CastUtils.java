/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
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
// USE_NEW_COPYRIGHT


package au.csiro.JamesGunning.sulib;

import java.lang.*;

/**
   Utility functions for casting byte arrays to int, shorts, floats.   

   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/10/2001 
*/

public class CastUtils
{
  /**
   * Utility function to convert 2 bytes into an short.
   *
   * @param bytes[offset], bytes[offset+1]
   *        The bytes.
   * @return
   *        The equivalent short.
   */

  public static short typecastToShort(byte[] bytes,int offset)
  {
	short value = 0;

	for (int i=offset; i < offset+2; ++i )
	  {
		value = (short) ( ( value << 8 ) | (bytes[i] & 0xFF)); 
	  }
	return value;
  }
  
  /**
   * Utility function to convert short into 2 bytes at marked offset.
   *
   * @param bytes[offset], bytes[offset+1]
   *        The bytes.
   * @param x, the value
   * @param offset, where in bytes array to stick it
   */
  public static void typecastShortToBytes(short x, byte[] bytes,int offset)
  {
	bytes[offset]   = (byte) ((x & 0x0000ff00) >> 8);
	bytes[offset+1] = (byte) ((x & 0x000000ff));
  }

  /**
   * Utility function to convert 4 bytes into an int.
   *
   * @param bytes[offset], bytes[offset+1], bytes[offset+2], bytes[offset+3]
   *        The bytes.
   * @return
   *        The equivalent int.
   */
  public static int typecastToInt(byte[] bytes,int offset)
  {
	int value = 0;
	for (int i=offset; i < offset+4; ++i ) 
	  {
		value = ( value << 8 ) | ( bytes[i] & 0xFF ); 
	  }
	return value;
	// See http://java.sun.com/people/linden/faq_b2.html#I/O
  }

  /**
   * Utility function to convert int into 4 bytes at marked offset.
   *
   * @param bytes[offset], bytes[offset+1]
   *        The bytes.
   * @param x, the value
   * @param offset, where in bytes array to stick it
   *
   */
  public static void typecastIntToBytes(int x, byte[] bytes,int offset)
  {
	bytes[offset]  = (byte) ((x & 0xff000000) >> 24);
	bytes[offset+1]= (byte) ((x & 0x00ff0000) >> 16);
	bytes[offset+2]= (byte) ((x & 0x0000ff00) >> 8);
	bytes[offset+3]= (byte) ((x & 0x000000ff));
  }


  /**
   * Utility function to convert 4 bytes into an float.
   *
   * @param bytes[offset], bytes[offset+1], bytes[offset+2], bytes[offset+3]
   *        The bytes.
   * @return
   *        The equivalent float.
   */
  
  public static float typecastToFloat(byte[] bytes,int offset)
  {
	int value  = ( ( (bytes[offset] & 0xff) << 24) |
                   ( (bytes[offset+1] & 0xff) << 16) |
                   ( (bytes[offset+2] & 0xff) <<  8) |
                     (bytes[offset+3] & 0xff));

	return Float.intBitsToFloat(value);
  }
  
  /**
   * Utility function to convert float into 4 bytes at marked offset.
   *
   * @param bytes[offset], bytes[offset+1]
   *        The bytes.
   * @param x, the value
   * @param offset, where in bytes array to stick it
   *
   */
  public static void typecastFloatToBytes(float x,byte[] bytes,int offset)
  {
	int xint = Float.floatToIntBits(x);
	bytes[offset]  = (byte) ((xint & 0xff000000) >> 24);
	bytes[offset+1]= (byte) ((xint & 0x00ff0000) >> 16);
	bytes[offset+2]= (byte) ((xint & 0x0000ff00) >> 8);
	bytes[offset+3]= (byte) ((xint & 0x000000ff));
  }
  /**
   * Utility function to convert a byte into an unsigned byte.
   *
   * @param b
   *        The byte.
   * @return
   *        The unsigned byte.
   */
  protected static int toUByte(byte b)
  {
	return b < 0 ? b + 0x100 : b;
  }
}


