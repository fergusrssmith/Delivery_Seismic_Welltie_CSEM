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
import java.io.*;
import java.lang.*;

/**
   Code to clone the behaviour of the surange program from the 
   CWP SeismicUnix (su) distribution
   
   Reads big-endian traces from stdin and prints a record of the range of header fields

   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/8/2001 
*/

public class SuRange
{ 
  public static void main(String[] args)
  {
	SuReader sureader = new SuReader();
	SuSegyFormat tr = new SuSegyFormat();
	SuSegyFormat trmin = new SuSegyFormat();
	SuSegyFormat trmax = new SuSegyFormat();
    FileInputStream stdin=null;
	int n_traces=0;
	short svalue, svalue_min, svalue_max;
	int ivalue, ivalue_min, ivalue_max;
	float fvalue, fvalue_min, fvalue_max;
	SuSegyFormat.SuHeaderEntry fields;
	SuSegyFormat.SuHeaderEntry min_fields;
	SuSegyFormat.SuHeaderEntry max_fields;

	// flip bytes from little endian to BE option
	boolean fromLE=false;

	if (args.length>=1)
	  {
		if (args[0].equals("--LE"))
		  {
			// flip input bytes
			fromLE=true;
			System.err.println("flipping bytes...");
		  }
		else
		  {
			System.err.println("Usage:"); 
			System.err.println("   SuRange [--LE] < stdin");
			System.err.println("            --LE for little endian input seismic (byte flip)");
			System.exit(1);
		  }
	  }


	// get input stream (stdin)
	stdin = (FileInputStream) new FileInputStream(FileDescriptor.in);

	// Get range of fields. We use the 
	// Hashmaps for tr, trmin and trmax,
	// in order to loop over header fields
	Collection trvalues = (Collection) tr.headerMap.values();
	

	// Find the min and max across all traces
	// this loop will get entries in the (screwed up) order they appear in the HashMap
	// We will undo this at the end using a sort


	// loop over traces 
	while (sureader.fgettr(stdin,tr,fromLE) != 0)
	  {
		n_traces++;

		// System.out.println("Trace: nl="+tr.get_int("tracl"));		
		// loop over header fields (screwed up order from HashMap)
		for (Iterator i = trvalues.iterator(); i.hasNext();) 
		  {
			fields = (SuSegyFormat.SuHeaderEntry)i.next();
			switch (fields.type) 
			  {
			  case SuSegyFormat.SuHeadTypes.short_type:
				svalue = tr.get_short(fields.key);
				trmin.set(fields.key,(n_traces==1)? svalue : ShMin(svalue,trmin.get_short(fields.key)));
				trmax.set(fields.key,(n_traces==1)? svalue : ShMax(svalue,trmax.get_short(fields.key)));
				break;
			  case SuSegyFormat.SuHeadTypes.int_type:
				ivalue = tr.get_int(fields.key);
				trmin.set(fields.key,(n_traces==1)? ivalue : Math.min(ivalue,trmin.get_int(fields.key)));
				trmax.set(fields.key,(n_traces==1)? ivalue : Math.max(ivalue,trmax.get_int(fields.key)));
				break;
			  case SuSegyFormat.SuHeadTypes.float_type:
				fvalue = tr.get_float(fields.key);
				trmin.set(fields.key,(n_traces==1)? fvalue : Math.min(fvalue,trmin.get_float(fields.key)));
				trmax.set(fields.key,(n_traces==1)? fvalue : Math.max(fvalue,trmax.get_float(fields.key)));
				break;
			  }
		  }

	  }

	// Now acquire output for the interesting ranges we found, 
	// again, in screwed up order, as found in trmin & trmax

	SuRangeOutputs[] suoutputs = new SuRangeOutputs[SuSeismicConstants.SU_NKEYS];

	int j=0;
	String text=null;
	boolean make_entry;
	for (Iterator i = trvalues.iterator(); i.hasNext();) 
	  {
		fields = (SuSegyFormat.SuHeaderEntry)i.next();
		int loc = fields.byte_location;
		make_entry = false;   // assume uninteresting until shown otherwise
		switch (fields.type) 
		  {
		  case SuSegyFormat.SuHeadTypes.short_type:
			svalue_min = trmin.get_short(fields.key);
			svalue_max = trmax.get_short(fields.key);
			if (svalue_min !=0 || svalue_max !=0)
			  {
				make_entry = true;
				if (svalue_min < svalue_max)
				  text = fields.key + "=(" + svalue_min + "," + svalue_max + ") ";
				else 
				  text = fields.key + "=" + svalue_min + " ";
			  }
			break;

		  case SuSegyFormat.SuHeadTypes.int_type:
			ivalue_min = trmin.get_int(fields.key);
			ivalue_max = trmax.get_int(fields.key);
			if (ivalue_min !=0 || ivalue_max !=0)
			  {
				make_entry = true;
				if (ivalue_min < ivalue_max)
				  text = fields.key + "=(" + ivalue_min + "," + ivalue_max + ") ";
				else 
				  text = fields.key + "=" + ivalue_min + " ";
			  }
			break;

		  case SuSegyFormat.SuHeadTypes.float_type:
			fvalue_min = trmin.get_float(fields.key);
			fvalue_max = trmax.get_float(fields.key);
			if (fvalue_min !=0 || fvalue_max !=0)
			  {
				make_entry = true;
				if (fvalue_min < fvalue_max)
				  text = fields.key + "=(" + fvalue_min + "," + fvalue_max + ") ";
				else 
				  text = fields.key + "=" + fvalue_min + " ";
			  }
			break;
		  }
		if (make_entry ) 
		  {
			// the text fragment we want...
			suoutputs[j] = new SuRange.SuRangeOutputs(loc,text);	
			j++;
		  }
	  }
	
	// Now sort outputs into the original header (pre-HashMap) order
	Arrays.sort(suoutputs,0,j);

	// And lets see it all in sane format

	System.out.println(n_traces + " traces:");System.out.print(" ");
	for (int k = 0; k<j ; k++) 
	  {
		System.out.print(suoutputs[k].text);
		if (((k+1) % 5) == 0) 
		  {System.out.println("");System.out.print(" ");}
	  }	
	System.out.println("");
	
  }


  private static short ShMin(short a, short b)
  {
	return (a<b)? a : b;
  }

  private static short ShMax(short a, short b)
  {
	return (a>b)? a : b;
  }

  /**
	 Class for sorting header outputs back into order of byte_location, since HashMap screws these up
  */
  private static class SuRangeOutputs implements Comparable 
  {
	int byte_location;
	String text;

	SuRangeOutputs(int byte_location, String text)
	{
	  this.byte_location=byte_location;
	  this.text=text;
	}

	/**
	   Comparator to sort by byte_location number (will be the byte location int the header
	*/
	public int compareTo(Object o)
	{
	  SuRangeOutputs n = (SuRangeOutputs)o;
	  int cmp = (byte_location<n.byte_location)? -1 : ((byte_location>n.byte_location)? 1 : 0 );
	  return cmp;
	}
  }

}
