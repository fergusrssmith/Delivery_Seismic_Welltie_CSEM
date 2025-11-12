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

import java.util.*;
import java.lang.*;
/**
   <p>Class containing storage, getters, and setters for seismic data
   in Big endian (java-native) su (clean segy) format. 
   The header format follows the header
   file segy.h in the SeismicUnix distribution from the Colorado
   School of Mines. The XDR (platform-independent format) option is
   not currently supported</p>

   <p>Objects in this class have data in su format, stored natively as
   one copy in byte[] arrays for efficiency. The byte arrays will
   correspond to big endian seismic data. </p>

   <p>The get_xxx(String key) and put_xxx(String key, value) are getters and
   setters of the named header fields</p>

   <p> Looping over header fields is best accomplished using a
   Collection formed from the HashMap of the header fields. See
   SuRange for an example</p>

   <p> data from the data part of the traces is available from the
   calls get_data_value(indx) and get_data_array()</p>

   <p> data from the data part of the traces may be set by the calls
   set_data_value(float datum, int indx) and set_data_array(float
   data[])</p>
  
   <p><b>Header fields and their meanings</b></p>  

   <p> <b>tracl</b><br>
   int tracl; trace sequence number within line<br>
   Use get_int("tracl") and put("tracl",value) 
   </p>

   <p> <b>tracr</b><br>
   int tracr; trace sequence number within reel<br>
   Use get_int("tracr") and put("tracr",value)
   </p> 

   <p> <b>fldr</b><br>
   int fldr;	 field record number<br>
   Use get_int("fldr") and put("fldr",value)
   </p> 
  
   <p> <b>tracf</b><br>
   int tracf; trace number within field record<br>
   Use get_int("tracf") and put("tracf",value)
   </p> 

   <p> <b>ep</b><br>
   int ep;	energy source point number<br>
   Use get_int("ep") and put("ep",value)
   </p> 

   <p> <b>cdp</b><br>
   int cdp;	CDP ensemble number 
   Use get_int("cdp") and put("cdp",value)
   </p> 

   <p> <b>cdpt</b><br>
   int cdpt;	trace number within CDP ensemble<br>
   Use get_int("cdpt") and put("cdpt",value)
   </p> 

   <p> <b>trid</b><br>
   short trid;   trace identification code:<br>
   1 = seismic data<br>
   2 = dead<br>
   3 = dummy<br>
   4 = time break<br>
   5 = uphole<br>
   6 = sweep<br>
   7 = timing<br>
   8 = water break<br>
   9---, N = optional use (N = 32,767)<br>

   Following are CWP id flags:<br>

   9 = autocorrelation<br>

   10 = Fourier transformed - no packing<br>
   xr[0],xi[0], ..., xr[N-1],xi[N-1]<br>

   11 = Fourier transformed - unpacked Nyquist<br>
   xr[0],xi[0],...,xr[N/2],xi[N/2]<br>

   12 = Fourier transformed - packed Nyquist<br>
   even N:<br>
   xr[0],xr[N/2],xr[1],xi[1], ...,<br>
   xr[N/2 -1],xi[N/2 -1]<br>
   (note the exceptional second entry)<br>
   odd N:<br>
   xr[0],xr[(N-1)/2],xr[1],xi[1], ...,<br>
   xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]<br>
   (note the exceptional second & last entries)<br>

   13 = Complex signal in the time domain<br>
   xr[0],xi[0], ..., xr[N-1],xi[N-1]<br>

   14 = Fourier transformed - amplitude/phase<br>
   a[0],p[0], ..., a[N-1],p[N-1]<br>

   15 = Complex time signal - amplitude/phase<br>
   a[0],p[0], ..., a[N-1],p[N-1]<br>

   16 = Real part of complex trace from 0 to Nyquist<br>

   17 = Imag part of complex trace from 0 to Nyquist<br>

   18 = Amplitude of complex trace from 0 to Nyquist<br>

   19 = Phase of complex trace from 0 to Nyquist<br>

   21 = Wavenumber time domain (k-t)<br>

   22 = Wavenumber frequency (k-omega)<br>

   23 = Envelope of the complex time trace<br>

   24 = Phase of the complex time trace<br>

   25 = Frequency of the complex time trace<br>

   30 = Depth-Range (z-x) traces<br>

   43 = Seismic Data, Vertical Component<br>

   44 = Seismic Data, Horizontal Component 1<br>

   45 = Seismic Data, Horizontal Component 2<br>

   46 = Seismic Data, Radial Component<br>

   47 = Seismic Data, Transverse Component<br>

   101 = Seismic data packed to bytes (by supack1)<br>

   102 = Seismic data packed to 2 bytes (by supack2)<br>

   Use get_short("trid") and put("trid",value)
   </p> 

   <p> <b>nvs</b><br>
   short nvs;	 number of vertically summed traces (see vscode in bhed structure)<br>
   Use get_short("nvs") and put("nvs",value)
   </p> 

   <p> <b>nhs</b><br>
   short nhs; number of horizontally summed traces (see vscode in bhed structure)<br>
   Use get_short("nhs") and put("nhs",value)
   </p> 

   <p> <b>duse</b><br>
   short duse  data use: 1 = production, 2 = test<br>
   Use get_short("duse") and put("duse",value) 
   </p>

   <p> <b>offset</b><br> 
   int offset;	distance from source point to receiver group<br>
   (negative if opposite to direction
   in which the line was shot)<br>
   Use get_int("offset") and put("offset",value)
   </p> 

   <p> <b>gelev</b><br>
   int gelev;  receiver group elevation from sea level<br>
   (above sea level is positive)<br>
   Use get_int("gelev") and put("gelev",value)
   </p> 

   <p> <b>selev</b><br>
   int selev;	source elevation from sea level<br>
   (above sea level is positive)<br>
   Use get_int("selev") and put("selev",value)
   </p>

   <p> <b>sdepth</b><br>
   int sdepth;	 source depth (positive)<br>
   Use get_int("sdepth") and put("sdepth",value)
   </p> 

   <p> <b>gdel</b><br>
   int gdel; datum elevation at receiver group<br>
   Use get_int("gdel") and put("gdel",value)
   </p> 

   <p> <b>sdel</b><br>
   datum elevation at source<br>
   Use get_int("sdel") and put("sdel",value)  	
   </p>

   <p> <b>swdep</b><br>
   water depth at source<br>
   Use get_int("swdep") and put("swdep",value)
   </p>

   <p> <b>gwdep</b><br>
   water depth at receiver group<br>
   Use get_int("gwdep") and put("gwdep",value)
   </p>

   <p> <b>scalel</b><br>
   scale factor for previous 7 entries<br>
   with value plus or minus 10 to the
   power 0, 1, 2, 3, or 4 (if positive,
   multiply, if negative divide)<br>
   Use get_short("scalel") and put("scalel",value)
   </p>

   <p> <b>scalco</b><br>
   scale factor for next 4 entries
   with value plus or minus 10 to the
   power 0, 1, 2, 3, or 4 (if positive,
   multiply, if negative divide)<br>
   Use get_short("scalco") and put("scalco",value)
   </p>

   <p> <b>sx</b><br>
   X source coordinate<br>
   Use get_int("sx") and put("sx",value)
   </p>

   <p> <b>sy</b><br>
   Y source coordinate<br>
   Use get_int("sy") and put("sy",value)
   </p>

   <p> <b>gx</b><br>
   X group coordinate<br>
   Use get_int("gx") and put("gx",value)
   </p>

   <p> <b>gy</b><br>
   Y group coordinate<br>
   Use get_int("gy") and put("gy",value)
   </p>

   <p> <b>counit</b><br>
   coordinate units code:
   for previous four entries<br>
   1 = length (meters or feet)<br>
   2 = seconds of arc (in this case, the<br>
   X values are longitude and the Y values
   are latitude, a positive value designates
   the number of seconds east of Greenwich
   or north of the equator<br>
   Use get_short("counit") and put("counit",value)
   </p>

   <p> <b>wevel</b><br>
   weathering velocity<br>
   Use get_short("wevel") and put("wevel",value)
   </p>

   <p> <b>swevel</b><br>
   subweathering velocity<br>
   Use get_short("swevel") and put("swevel",value)
   </p>

   <p> <b>sut</b><br>
   uphole time at source<br>
   Use get_short("sut") and put("sut",value)
   </p>

   <p> <b>gut</b><br>
   uphole time at receiver group<br>
   Use get_short("gut") and put("gut",value)
   </p>

   <p> <b>sstat</b><br>
   source static correction<br>
   Use get_short("sstat") and put("sstat",value)
   </p>

   <p> <b>gstat</b><br>
   group static correction<br>
   Use get_short("gstat") and put("gstat",value)
   </p>

   <p> <b>tstat</b><br>
   total static applied<br>
   Use get_short("tstat") and put("tstat",value)
   </p>

   <p> <b>laga</b><br>
   lag time A, time in ms between end of 240-
   byte trace identification header and time
   break, positive if time break occurs after
   end of header, time break is defined as
   the initiation pulse which maybe recorded
   on an auxiliary trace or as otherwise
   specified by the recording system<br>
   Use get_short("laga") and put("laga",value)
   </p>

   <p> <b>lagb</b><br>
   lag time B, time in ms between the time break
   and the initiation time of the energy source,
   may be positive or negative<br>
   Use get_short("lagb") and put("lagb",value)
   </p>

   <p> <b>delrt</b><br>
   delay recording time, time in ms between
   initiation time of energy source and time
   when recording of data samples begins
   (for deep water work if recording does not
   start at zero time)<br>
   Use get_short("delrt") and put("delrt",value)
   </p>

   <p> <b>muts</b><br>
   mute time--start<br>
   Use get_short("muts") and put("muts",value)
   </p>

   <p> <b>mute</b><br>
   mute time--end<br>
   Use get_short("mute") and put("mute",value)
   </p>

   <p> <b>ns</b><br>
   number of samples in this trace<br>
   Use get_short("ns") and put("ns",value)
   </p>

   <p> <b>dt</b><br>
   sample interval in micro-seconds<br>
   get_short("dt");
   </p>

   <p> <b>gain</b><br>
   gain type of field instruments code:<br>
   1 = fixed<br>
   2 = binary<br>
   3 = floating point<br>
   4 ---- N = optional use<br>
   Use get_short("gain") and put("gain",value)
   </p>

   <p> <b>igc</b><br>
   instrument gain constant<br>
   Use get_short("igc") and put("igc",value)
   </p>

   <p> <b>igi</b><br>
   instrument early or initial gain<br>
   Use get_short("igi") and put("igi",value)
   </p>

   <p> <b>corr</b><br>
   correlated:<br>
   1 = no<br>
   2 = yes<br>
   Use get_short("corr") and put("corr",value)
   </p>

   <p> <b>sfs</b><br>
   sweep frequency at start<br>
   Use get_short("sfs") and put("sfs",value)
   </p>

   <p> <b>sfe</b><br>
   sweep frequency at end<br>
   Use get_short("sfe") and put("sfe",value)
   </p>

   <p> <b>slen</b><br>
   sweep length in ms<br>
   Use get_short("slen") and put("slen",value)
   </p>

   <p> <b>styp</b><br>
   sweep type code:<br>
   1 = linear<br>
   2 = cos-squared<br>
   3 = other<br>
   Use get_short("styp") and put("styp",value)
   </p>

   <p> <b>stas</b><br>
   sweep trace length at start in ms<br>
   Use get_short("stas") and put("stas",value)
   </p>

   <p> <b>stae</b><br>
   sweep trace length at end in ms<br>
   Use get_short("stae") and put("stae",value)
   </p>

   <p> <b>tatyp</b><br>
   taper type: 1=linear, 2=cos^2, 3=other<br>
   Use get_short("tatyp") and put("tatyp",value)
   </p>

   <p> <b>afilf</b><br>
   alias filter frequency if used<br>
   Use get_short("afilf") and put("afilf",value)
   </p>

   <p> <b>afils</b><br>
   alias filter slope<br>
   Use get_short("afils") and put("afils",value)
   </p>

   <p> <b>nofilf</b><br>
   notch filter frequency if used<br>
   Use get_short("nofilf") and put("nofilf",value)
   </p>

   <p> <b>nofils</b><br>
   notch filter slope<br>
   Use get_short("nofils") and put("nofils",value)
   </p>

   <p> <b>lcf</b><br>
   low cut frequency if used<br>
   Use get_short("lcf") and put("lcf",value)
   </p>

   <p> <b>hcf</b><br>
   high cut frequncy if used<br>
   Use get_short("hcf") and put("hcf",value)
   </p>

   <p> <b>lcs</b><br>
   low cut slope<br>
   Use get_short("lcs") and put("lcs",value)
   </p>

   <p> <b>hcs</b><br>
   high cut slope<br>
   Use get_short("hcs") and put("hcs",value)
   </p>

   <p> <b>year</b><br>
   year data recorded<br>
   Use get_short("year") and put("year",value)
   </p>

   <p> <b>day</b><br>
   day of year<br>
   Use get_short("day") and put("day",value)
   </p>

   <p> <b>hour</b><br>
   hour of day (24 hour clock)<br>
   Use get_short("hour") and put("hour",value)
   </p>

   <p> <b>minute</b><br>
   minute of hour<br>
   Use get_short("minute") and put("minute",value)
   </p>

   <p> <b>sec</b><br>
   second of minute<br>
   Use get_short("sec") and put("sec",value)
   </p>

   <p> <b>timbas</b><br>
   time basis code:<br>
   1 = local<br>
   2 = GMT<br>
   3 = other<br>
   Use get_short("timbas") and put("timbas",value)
   </p>

   <p> <b>trwf</b><br>
   trace weighting factor, defined as 1/2^N volts for the least sigificant bit<br>
   Use get_short("trwf") and put("trwf",value)
   </p>

   <p> <b>grnors</b><br>
   geophone group number of roll switch position one<br>
   Use get_short("grnors") and put("grnors",value)
   </p>

   <p> <b>grnofr</b><br>
   geophone group number of trace one within original field record<br>
   Use get_short("grnofr") and put("grnofr",value)
   </p>

   <p> <b>grnlof</b><br>
   geophone group number of last trace within original field record<br>
   Use get_short("grnlof") and put("grnlof",value)
   </p>

   <p> <b>gaps</b><br>
   gap size (total number of groups dropped)<br>
   Use get_short("gaps") and put("gaps",value)
   </p>

   <p> <b>otrav</b><br>
   overtravel taper code:<br>
   1 = down (or behind)<br>
   2 = up (or ahead)<br>
   Use get_short("otrav") and put("otrav",value)
   </p>

  <p><b>local assignments:</b></p>

  <p> <b>d1</b><br>
  sample spacing for non-seismic data<br>
  Use get_float("d1") and put("d1",value)
  </p>

  <p> <b>f1</b><br>
  first sample location for non-seismic data<br>
  Use get_float("f1") and put("f1",value)
  </p>

  <p> <b>d2</b><br>
  sample spacing between traces<br>
  Use get_float("d2") and put("d2",value)
  </p>

  <p> <b>f2</b><br>
  first trace location<br>
  Use get_float("f2") and put("f2",value)
  </p>

  <p> <b>ungpow</b><br>
  negative of power used for dynamic range compression<br>
  Use get_float("ungpow") and put("ungpow",value)
  </p>

  <p> <b>unscale</b><br>
  reciprocal of scaling factor to normalize range<br>
  Use get_float("unscale") and put("unscale",value)
  </p>

  <p> <b>ntr</b><br>
  number of traces<br>
  Use get_int("ntr") and put("ntr",value)
  </p>

  <p> <b>mark</b><br>
  mark selected traces<br>
  Use get_short("mark") and put("mark",value)
  </p>

  <p> <b>shortpad</b><br>
  alignment padding<br>
  Use get_short("shortpad") and put("shortpad",value)
  </p>


  <p>
  Bits in Segy.h which aren't needed...<br>

  short unass[14];	<br>
  unassigned--NOTE: last entry causes 
  a break in the word alignment, if we REALLY
  want to maintain 240 bytes, the following
  entry should be an odd number of short/UINT2
  OR do the insertion above the "mark" keyword
  entry
  </p>



   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/8/2001
*/

public class SuSegyFormat
{
  protected byte[] hdr_bytes;     // HDRBYTES  (240)
  protected byte[] data_bytes;    // SU_NFLTS  (32768); Arbitrary limit on data array size
  
  /**
	 A hashmap of all header fields created with
	 <em>setupHeaderMap()</em> (private) at each object instantiation.
	 It is always available for checking validity of keywords etc.
	 The <b>key</b> is the header field name (String), and the
	 <b>values</b> are instances of <b>SuHeaderEntry</b> (header name
	 (String), SuHeadTypes, byte_location) 
  */
  public HashMap headerMap=null;  // if you want to do smart stuff with this
  
  /**
	 Collection used for iterating over header fields. A slave to headerMap
  */
  Collection headerMapColl = null;

  /**
	 Constructor allocs room for header and data using constants as
	 per seismic constants 
  */
  public SuSegyFormat()
  {
	hdr_bytes = new byte[SuSeismicConstants.HDRBYTES];
	data_bytes = new byte[4 * SuSeismicConstants.SU_NFLTS];
	setupHeaderMap();
  }

  /**
	 Copying all existing header and data fields from the SuSegyFormat object tr
	 into <em>this</em>.
	 @param tr The source SuSegyFormat trace
	 @return
  */
  public void copy(SuSegyFormat tr) 
  {
	System.arraycopy(tr.hdr_bytes,0,hdr_bytes,0,SuSeismicConstants.HDRBYTES);
	System.arraycopy(tr.data_bytes,0,data_bytes,0,4*SuSeismicConstants.SU_NFLTS);
  } 	
  /**
	 Copying header fields from the SuSegyFormat object tr into <em>this</em>.
	 @param
	 tr The source SuSegyFormat trace
	 @return
  */
  public void copyHeader(SuSegyFormat tr) 
  {
	System.arraycopy(tr.hdr_bytes,0,hdr_bytes,0,SuSeismicConstants.HDRBYTES);
  } 
  /**
	 Copying header fields from <em>this</em> into user-supplied dest[0...SuSeismicConstants.HDRBYTES-1]
	 @param
	 dest the target byte array
	 @return
  */
  public void copyHeader(byte[] dest) 
  {
	System.arraycopy(hdr_bytes,0, dest,0,SuSeismicConstants.HDRBYTES);
  } 

  /**
	 Copying header fields into <em>this</em> from user-supplied dest[0...SuSeismicConstants.HDRBYTES-1]
	 @param
	 dest the src byte array
	 @return
  */
  public void setHeader(byte[] dest) 
  {
	System.arraycopy(dest,0,hdr_bytes,0,SuSeismicConstants.HDRBYTES);
  } 

  /**
	 Zero the trace header
  */
  public void zeroHeader() 
  {
	Arrays.fill(hdr_bytes,(byte)0);
  } 

  // functions to interrogate data

  /**
	 Get the ith data entry as a float from the raw data array in bytes
	 @param
	 i, the position in the float array
	 @return
	 A float corresponding to position i int the data array, cast from the 
	 4 corresponding  big-endian bytes
  */
  public float get_data_value(int i) 
  {
	if (i<0 || i>= SuSeismicConstants.SU_NFLTS) 
	  throw new IllegalArgumentException(
										 "indx i must be less than SU_NFLTS in SuSegyFormat.get_data_value");
	else
	  return CastUtils.typecastToFloat(data_bytes,4*i);
  } 	

  /**
	 Convert the bytes read into from the data part of the traces to a
	 float array if required
	  
	 @return
	 An array of floats of size this.get_short("ns") representing the data */
	 
	 
  public float[] get_data_array() 
  {
	int i,ib,n;

	n = this.get_short("ns");
	if (n<1) 
	  {
		throw new IllegalArgumentException(
										   "ns zero or less in SuSegyFormat.get_data_array");		
	  }	
	float[] data_floats = new float[n];
	for (i=0,ib=0;i<n;i++,ib+=4) 
	  {
		data_floats[i] = CastUtils.typecastToFloat(data_bytes,ib);
	  }
	return data_floats;
  }


  /**
	 Writes the supplied value into the data part of the trace.
	 The data is dumped for indices i=0...n-1
	 
	 @param data = array of floats.
  */
 
  public void fill_data_array(float value) 
  {
	int i,ib,n;

	n = (int)this.get_short("ns");
	if (n<1) 
	  {
		throw new IllegalArgumentException("ns zero or less in SuSegyFormat.set_data_array");		
	  }	
	for (i=0,ib=0; i<n; i++,ib+=4) 
	  {
		CastUtils.typecastFloatToBytes(value,data_bytes,ib);
	  }
  }

  /**
	 Convert the float array into bytes in the data part of the trace,
	 with no offsets. Virtually a data dump. The data is dumped
	 for indices i=0...n, with n = min(ns,data.length) (ns = no. samples
	 in trace)
	 
	 @param data = array of floats.
  */
	 
	 
  public void set_data_array(float data[]) 
  {
	int i,ib,n;

	n = Math.min((int)this.get_short("ns"),data.length);
	if (n<1) 
	  {
		throw new IllegalArgumentException("ns zero or less in SuSegyFormat.set_data_array");		
	  }	
	for (i=0,ib=0; i<n; i++,ib+=4) 
	  {
		CastUtils.typecastFloatToBytes(data[i],data_bytes,ib);
	  }
  }

  /**
	 Zero all the trace data values.
  */
	 
	 
  public void zero_data_array() 
  {
	java.util.Arrays.fill(data_bytes,(byte)0);
  }

  /**
	 Dump the the float into 4 bytes in the data part of the trace,
	 at the stated float index i (i.e. byte index 4i)
	  
	 @param f float value to dump into trace
	 @param i index at which to write it

  */
	 
	 
  public void set_data_value(float f, int i) 
  {
	int ib=4*i;

	if (i<0 || (i >= this.get_short("ns")))
	  {
		throw new IllegalArgumentException("index too low or high in SuSegyFormat.set_data_value");		
	  }
	
	CastUtils.typecastFloatToBytes(f,data_bytes,ib);
  }


  void swapHeaderValues()
  {
	SuSegyFormat.SuHeaderEntry fields;
	short s;
	int i;
	float f;
	Iterator I;

	// loop over header fields (screwed up order from HashMap)
	// Get range of fields. We use the 
	// Hashmaps
	// in order to loop over header fields
	// first call
	if (headerMapColl == null) headerMapColl = (Collection) headerMap.values();

	for (I = headerMapColl.iterator(); I.hasNext();) 
	  {
		fields = (SuSegyFormat.SuHeaderEntry)I.next();
		switch (fields.type) 
		  {
		  case SuSegyFormat.SuHeadTypes.short_type:
			s = CastUtils.typecastToShort(hdr_bytes,fields.byte_location);
			// flip
			s = (short)((((s>>8)&0xff) | ((s&0xff)<<8))); 
			// write back
			CastUtils.typecastShortToBytes(s,hdr_bytes,fields.byte_location);
			break;

		  case SuSegyFormat.SuHeadTypes.int_type:
			i = CastUtils.typecastToInt(hdr_bytes,fields.byte_location);
			// flip
			i = (((i>>24)&0xff) | ((i&0xff)<<24) | ((i>>8)&0xff00) | ((i&0xff00)<<8));  
			// write back
			CastUtils.typecastIntToBytes(i,hdr_bytes,fields.byte_location);
			break;

		  case SuSegyFormat.SuHeadTypes.float_type:
			
			/*
			  f = CastUtils.typecastToFloat(hdr_bytes,fields.byte_location);
			  // flip
			  i = Float.floatToIntBits(f);
			*/
			i = CastUtils.typecastToInt(hdr_bytes,fields.byte_location);
			i = (((i>>24)&0xff) | ((i&0xff)<<24) | ((i>>8)&0xff00) | ((i&0xff00)<<8));  
			f = Float.intBitsToFloat(i);
			// write back
			CastUtils.typecastFloatToBytes(f,hdr_bytes,fields.byte_location);
			break;
		  }
	  }
	return;
  }

  void swapDataValues()
  {
	int i,j,j4,_ns;

	_ns = get_short("ns");
	for(j=0,j4=0; j<_ns; j++,j4+=4)
	  {
		i = CastUtils.typecastToInt(data_bytes,j4);
		i = (((i>>24)&0xff) | ((i&0xff)<<24) | ((i>>8)&0xff00) | ((i&0xff00)<<8));  
		// write back
		CastUtils.typecastFloatToBytes(Float.intBitsToFloat(i),data_bytes,j4);
	  }
  }

  /** 
	  
	  <p> Function to set up a hashMap for all header fields, which can be
	  used to simply iterate over them in programs like SuRange. Note
	  that the hashMap screws up the natural order of the header fields
	  so a collection generated by headerMap.values() will produce
	  a sensible Iterator, but not in natural order.</p>


  */

  private void setupHeaderMap()
  {
	// Set up a hashMap for all header fields 
	// ensure initialCapacity*loadFactor > SuSeismicConstants.SU_NKEYS
	// default initialCapacity = 0.75
	headerMap = new HashMap(1+(int)((float)SuSeismicConstants.SU_NKEYS/0.75));


	headerMap.put("tracl", new SuHeaderEntry("tracl",SuHeadTypes.int_type,0));
	headerMap.put("tracr", new SuHeaderEntry("tracr",SuHeadTypes.int_type,4));
	headerMap.put("fldr", new SuHeaderEntry("fldr",SuHeadTypes.int_type,8));
	headerMap.put("tracf", new SuHeaderEntry("tracf",SuHeadTypes.int_type,12));
	headerMap.put("ep", new SuHeaderEntry("ep",SuHeadTypes.int_type,16));
	headerMap.put("cdp", new SuHeaderEntry("cdp",SuHeadTypes.int_type,20));
	headerMap.put("cdpt", new SuHeaderEntry("cdpt",SuHeadTypes.int_type,24));
	headerMap.put("trid", new SuHeaderEntry("trid",SuHeadTypes.short_type,28));
	headerMap.put("nvs", new SuHeaderEntry("nvs",SuHeadTypes.short_type,30));
	headerMap.put("nhs", new SuHeaderEntry("nhs",SuHeadTypes.short_type,32));
	headerMap.put("duse", new SuHeaderEntry("duse",SuHeadTypes.short_type,34));
	headerMap.put("offset", new SuHeaderEntry("offset",SuHeadTypes.int_type,36));
	headerMap.put("gelev", new SuHeaderEntry("gelev",SuHeadTypes.int_type,40));
	headerMap.put("selev", new SuHeaderEntry("selev",SuHeadTypes.int_type,44));
	headerMap.put("sdepth", new SuHeaderEntry("sdepth",SuHeadTypes.int_type,48));
	headerMap.put("gdel", new SuHeaderEntry("gdel",SuHeadTypes.int_type,52));
	headerMap.put("sdel", new SuHeaderEntry("sdel",SuHeadTypes.int_type,56));
	headerMap.put("swdep", new SuHeaderEntry("swdep",SuHeadTypes.int_type,60));
	headerMap.put("gwdep", new SuHeaderEntry("gwdep",SuHeadTypes.int_type,64));
	headerMap.put("scalel", new SuHeaderEntry("scalel",SuHeadTypes.short_type,68));
	headerMap.put("scalco", new SuHeaderEntry("scalco",SuHeadTypes.short_type,70));
	headerMap.put("sx", new SuHeaderEntry("sx",SuHeadTypes.int_type,72));
	headerMap.put("sy", new SuHeaderEntry("sy",SuHeadTypes.int_type,76));
	headerMap.put("gx", new SuHeaderEntry("gx",SuHeadTypes.int_type,80));
	headerMap.put("gy", new SuHeaderEntry("gy",SuHeadTypes.int_type,84));
	headerMap.put("counit", new SuHeaderEntry("counit",SuHeadTypes.short_type,88));
	headerMap.put("wevel", new SuHeaderEntry("wevel",SuHeadTypes.short_type,90));
	headerMap.put("swevel", new SuHeaderEntry("swevel",SuHeadTypes.short_type,92));
	headerMap.put("sut", new SuHeaderEntry("sut",SuHeadTypes.short_type,94));
	headerMap.put("gut", new SuHeaderEntry("gut",SuHeadTypes.short_type,96));
	headerMap.put("sstat", new SuHeaderEntry("sstat",SuHeadTypes.short_type,98));
	headerMap.put("gstat", new SuHeaderEntry("gstat",SuHeadTypes.short_type,100));
	headerMap.put("tstat", new SuHeaderEntry("tstat",SuHeadTypes.short_type,102));
	headerMap.put("laga", new SuHeaderEntry("laga",SuHeadTypes.short_type,104));
	headerMap.put("lagb", new SuHeaderEntry("lagb",SuHeadTypes.short_type,106));
	headerMap.put("delrt", new SuHeaderEntry("delrt",SuHeadTypes.short_type,108));
	headerMap.put("muts", new SuHeaderEntry("muts",SuHeadTypes.short_type,110));
	headerMap.put("mute", new SuHeaderEntry("mute",SuHeadTypes.short_type,112));
	headerMap.put("ns", new SuHeaderEntry("ns",SuHeadTypes.short_type,114));
	headerMap.put("dt", new SuHeaderEntry("dt",SuHeadTypes.short_type,116));
	headerMap.put("gain", new SuHeaderEntry("gain",SuHeadTypes.short_type,118));
	headerMap.put("igc", new SuHeaderEntry("igc",SuHeadTypes.short_type,120));
	headerMap.put("igi", new SuHeaderEntry("igi",SuHeadTypes.short_type,122));
	headerMap.put("corr", new SuHeaderEntry("corr",SuHeadTypes.short_type,124));
	headerMap.put("sfs", new SuHeaderEntry("sfs",SuHeadTypes.short_type,126));
	headerMap.put("sfe", new SuHeaderEntry("sfe",SuHeadTypes.short_type,128));
	headerMap.put("slen", new SuHeaderEntry("slen",SuHeadTypes.short_type,130));
	headerMap.put("styp", new SuHeaderEntry("styp",SuHeadTypes.short_type,132));
	headerMap.put("stas", new SuHeaderEntry("stas",SuHeadTypes.short_type,134));
	headerMap.put("stae", new SuHeaderEntry("stae",SuHeadTypes.short_type,136));
	headerMap.put("tatyp", new SuHeaderEntry("tatyp",SuHeadTypes.short_type,138));
	headerMap.put("afilf", new SuHeaderEntry("afilf",SuHeadTypes.short_type,140));
	headerMap.put("afils", new SuHeaderEntry("afils",SuHeadTypes.short_type,142));
	headerMap.put("nofilf", new SuHeaderEntry("nofilf",SuHeadTypes.short_type,144));
	headerMap.put("nofils", new SuHeaderEntry("nofils",SuHeadTypes.short_type,146));
	headerMap.put("lcf", new SuHeaderEntry("lcf",SuHeadTypes.short_type,148));
	headerMap.put("hcf", new SuHeaderEntry("hcf",SuHeadTypes.short_type,150));
	headerMap.put("lcs", new SuHeaderEntry("lcs",SuHeadTypes.short_type,152));
	headerMap.put("hcs", new SuHeaderEntry("hcs",SuHeadTypes.short_type,154));
	headerMap.put("year", new SuHeaderEntry("year",SuHeadTypes.short_type,156));
	headerMap.put("day", new SuHeaderEntry("day",SuHeadTypes.short_type,158));
	headerMap.put("hour", new SuHeaderEntry("hour",SuHeadTypes.short_type,160));
	headerMap.put("minute", new SuHeaderEntry("minute",SuHeadTypes.short_type,162));
	headerMap.put("sec", new SuHeaderEntry("sec",SuHeadTypes.short_type,164));
	headerMap.put("timbas", new SuHeaderEntry("timbas",SuHeadTypes.short_type,166));
	headerMap.put("trwf", new SuHeaderEntry("trwf",SuHeadTypes.short_type,168));
	headerMap.put("grnors", new SuHeaderEntry("grnors",SuHeadTypes.short_type,170));
	headerMap.put("grnofr", new SuHeaderEntry("grnofr",SuHeadTypes.short_type,172));
	headerMap.put("grnlof", new SuHeaderEntry("grnlof",SuHeadTypes.short_type,174));
	headerMap.put("gaps", new SuHeaderEntry("gaps",SuHeadTypes.short_type,176));
	headerMap.put("otrav", new SuHeaderEntry("otrav",SuHeadTypes.short_type,178));
	headerMap.put("d1", new SuHeaderEntry("d1",SuHeadTypes.float_type,180));
	headerMap.put("f1", new SuHeaderEntry("f1",SuHeadTypes.float_type,184));
	headerMap.put("d2", new SuHeaderEntry("d2",SuHeadTypes.float_type,188));
	headerMap.put("f2", new SuHeaderEntry("f2",SuHeadTypes.float_type,192));
	headerMap.put("ungpow", new SuHeaderEntry("ungpow",SuHeadTypes.float_type,196));
	headerMap.put("unscale", new SuHeaderEntry("unscale",SuHeadTypes.float_type,200));
	headerMap.put("ntr", new SuHeaderEntry("ntr",SuHeadTypes.int_type,204));
	headerMap.put("mark", new SuHeaderEntry("mark",SuHeadTypes.short_type,208));
	headerMap.put("shortpad", new SuHeaderEntry("shortpad",SuHeadTypes.short_type,210));

  }

  /** 
	  Class simulating an enum of storage types in su header (shorts/ints/floats)
  */
  public class SuHeadTypes
  {
	public static final short short_type=0;
	public static final short int_type=1;
	public static final short float_type=2;	
  }

  /** 
	  Class describing header entries in Hashmap values
  */
  public class SuHeaderEntry
  {
	/** 
		keyword name
	*/
	public String key;

	/** 
		Keyword storage type
	*/
	public short type;
	/** 
		Byte location of this fields in the su header
	*/
	public int byte_location;

	SuHeaderEntry(String key, short type, int byte_location)
	{
	  this.key = key;
	  this.type = type;
	  this.byte_location = byte_location;
	}
  }

  // generic header field reader (probably inefficient)

  /** 
	  Handy get function to retrieve a Short su header word by name.

	  @param key   name of header key word, e.g. "ns"
	  
	  @return value of field if found, 0 otherwise
  */

  public short get_short(String key)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return 0;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null) 
	  {
		return CastUtils.typecastToShort(hdr_bytes,fields.byte_location);
	  }
	else 
	  {
		System.err.println("get_short: unknown header key " + key);
		System.exit(1);
	  }
	return 0;
  }
  
  /** 
	  Handy get function to retrieve a Integer su header word by name.

	  @param key   name of header key word, e.g. "ns"
	  
	  @return value of field if found, 0 otherwise
  */
  public int get_int(String key)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return 0;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null) 
	  {
		return CastUtils.typecastToInt(hdr_bytes,fields.byte_location);
	  }
	else 
	  {
		System.err.println("get_int; unknown header key " + key);
		System.exit(1);
	  }
	return 0;
  }

  /** 
	  Handy get function to retrieve a Float su header word by name.

	  @param key   name of header key word, e.g. "ns"
	  
	  @return value of field if found, 0 otherwise
  */
  public float get_float(String key)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return 0;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null) 
	  {
		return CastUtils.typecastToFloat(hdr_bytes,fields.byte_location);
	  }
	else 
	  {
		System.err.println("get_float: unknown header key " + key);
		System.exit(1);
	  }
	return 0;
  }

/** 
	  Handy get function to retrieve a general su header word by name, typecast to a float.

	  @param key   name of header key word, e.g. "ns"
	  
	  @return value of field if found, 0 otherwise
  */

  public float getAsFloat(String key)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return 0;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null) 
	  {
		switch (fields.type)
		  {
		  case (SuHeadTypes.short_type):
			return (float)CastUtils.typecastToShort(hdr_bytes,fields.byte_location);
		  case (SuHeadTypes.int_type):
			return (float)CastUtils.typecastToInt(hdr_bytes,fields.byte_location);
		  case (SuHeadTypes.float_type):
			return CastUtils.typecastToFloat(hdr_bytes,fields.byte_location);
		  }
	  }
	else 
	  {
		System.err.println("getAsFloat(): unknown header key " + key);
		System.exit(1);
	  }
	return 0;
  }

  /** 
	  Handy setter function to set a Short su header word by name.

	  @param key   name of header key word, e.g. "ns"
	  @param value  value to set header key word to.
	  
	  @return false if failed, true if success
  */
  public boolean set(String key, short value)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return false;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null && fields.type == SuHeadTypes.short_type) 
	  {
		CastUtils.typecastShortToBytes(value,hdr_bytes,fields.byte_location);
		return true;
	  }
	else 
	  {
		System.err.println("unknown or wrong type header key " + key);
		System.exit(1);
		return false;
	  }
  }

  /** 
	  Handy setter function to set a Integer su header word by name.

	  @param key   name of header key word, e.g. "tracl"
	  @param value  value to set header key word to.
	  
	  @return false if failed, true if success
  */

  public boolean set(String key, int value)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return false;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null && fields.type == SuHeadTypes.int_type) 
	  {
		CastUtils.typecastIntToBytes(value,hdr_bytes,fields.byte_location);
		return true;
	  }
	else 
	  {
		System.err.println("unknown or wrong type header key " + key);
		System.exit(1);
		return false;
	  }
  }

  /** 
	  Handy setter function to set a Float su header word by name.
	  @param key   name of header key word, e.g. "d1"
	  @param value  value to set header key word to.
	
	  @return false if failed, true if success
  */
  
  public boolean set(String key, float value)
  {
	if (headerMap==null)
	  {
		System.err.println("headerMap has not been initialised!!");
		return false;
	  }
	SuHeaderEntry fields = (SuHeaderEntry) headerMap.get(key);
	if (fields != null && fields.type == SuHeadTypes.float_type) 
	  {
		CastUtils.typecastFloatToBytes(value,hdr_bytes,fields.byte_location);
		return true;
	  }
	else 
	  {
		System.err.println("unknown or wrong type header key "+key);
		System.exit(1);
		return false;
	  }
  }


}


