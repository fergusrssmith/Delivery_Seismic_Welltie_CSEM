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

import java.io.*;
import java.lang.*;
import java.util.*;
import java.text.*;

/**
   Test harness for SuReader. Use either stdin or
   a filename.su as arguments. 

   Simple mode: echos the traces from stdin.

   With arguments   "--dup filename.su",
   will duplicate filename.su -> dup.filename.su for as
   many traces as are in stdin. (This tests
   the stream/file hashmaps in the Readers/Writers)
   
   
   @author james.gunning@dpr.csiro.au
   @version 1.0, 9/10/2001 
*/

public class TestSuReader
{
  public static void main(String[] args)
  {
	SuReader sureader = new SuReader();
	SuWriter suwriter = new SuWriter();
	SuSegyFormat tr = new SuSegyFormat();
	SuSegyFormat tr2 = new SuSegyFormat();
	Object in_strm,in_strm2=null,out_strm,out_strm2=null;
	float[] data;
	String dup_name=null;
	boolean from_stdin,dup=false;

	if (args.length>0 && args[0].equals("--dup")) 
	  {
		if (args.length<2) SuErrorHandler.die_err("Need --dup filename");
		dup = true;  // inputs from both stdin and file args[1]
		dup_name = "dup." +args[1];
	  }

	// where's it coming from ? /
	from_stdin = (args.length==0 || dup) ? true : false;

	if (from_stdin) 
	  { // no file specified, try stdin
		in_strm = (FileInputStream) new FileInputStream(FileDescriptor.in); // FIXME undo
	  }
	else
	  {
		try {
		  in_strm = (RandomAccessFile) new RandomAccessFile(args[0],"r");
		}
		catch (IOException e) {System.out.println("couldn't open " + args[0]);return;}
	  }

	if (dup) // second input stream from file args[1]
	  {
		try {
		  in_strm2 = (RandomAccessFile) new RandomAccessFile(args[1],"r");
		}
		catch (IOException e) {System.out.println("couldn't open " + args[1]);return;}
	  }


	out_strm = (FileOutputStream) new FileOutputStream(FileDescriptor.out);

	if (dup)
	  {
		try {
		  out_strm2 = (FileOutputStream) new FileOutputStream(dup_name);
		} catch (IOException e) {SuErrorHandler.die_err("couldn't open" + dup_name);}
	  }
	int ntr=1;
	
	while ((from_stdin ? sureader.fgettr((FileInputStream)in_strm,tr) : 
			sureader.fgettr((RandomAccessFile)in_strm,tr)) != 0)
		{
		// Simple stuff.......
		//		data = tr.get_data_array_as_floats();
		//  		for (int i=0;i<5;i++ ) 
		//  		  {
		//  			System.out.println(data[i]);
		//  		  }
		// System.out.println("");

		// demanding stuff....
		if (dup)
		  {
			// duplicate file while stdin yields traces
			if (sureader.fgettr((RandomAccessFile)in_strm2,tr2) != 0)
			  {
				suwriter.fputtr((FileOutputStream)out_strm2,tr2);
			  }
		  }


		
		// and echo stdin
		suwriter.fputtr((FileOutputStream)out_strm,tr);
		
		ntr++;
	  }

	if (!from_stdin) 
	  {
		try {
		  ((RandomAccessFile)in_strm).close(); 
		}
		catch (IOException e) {System.out.println("couldn't close" + args[0]);}
	  }

	if (dup) 
	  {
		try {
		  ((RandomAccessFile)in_strm2).close(); 
		}
		catch (IOException e) {System.out.println("couldn't close" + args[1]);}
		
		try {
		  ((FileOutputStream)out_strm2).close(); 
		}
		catch (IOException e) {System.out.println("couldn't close" + dup_name);}

	  }
  }
}
