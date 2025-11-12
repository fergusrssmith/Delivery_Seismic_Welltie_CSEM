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
import java.io.*;

/**
   Basic class to put seismic traces in su format, rather
   like the fputtr.c suite of functions in SeismicUnix.

   Given a SuSegyFormat object to hold the trace, and an
   opened RandomAccessFile (for disk files) or
   FileOutputStream (for stdin), traces are written
   from the SuSegyFormat by the 

   fputtr(RandomAccessFile fp, SuSegyFormat tr)

   or

   fputtr(FileOutputStream fp, SuSegyFormat tr)

   calls. An example usage fragment is<br><br>



   FileInputStream stdin;<br>
   FileOutputStream stdout;<br>

   SuReader sureader = new SuReader();<br>
   SuWriter suwriter = new SuWriter();<br>
   SuSegyFormat tr = new SuSegyFormat();<br>

   // get input stream (stdin)<br>
   stdin = (FileInputStream) new FileInputStream(FileDescriptor.in);<br>

   // get out stream (stdout)<br>
   stdout = (FileOutputStream) new FileOutputStream(FileDescriptor.out);<br>
   fileout = (FileOutputStream) new FileOutputStream(filename);<br>

   while (sureader.fgettr(stdin,tr) != 0) // get traces<br>
   {<br>
     // do stuff with traces here....<br>
     suwriter.fputtr((FileOutputStream)stdout,tr);  // to stdout<br>
	 suwriter.fputtr((FileOutputStream)fileout,tr);  // to a file<br>
   }<br>
   ...<br>
   // To close, if the Output stream is attached to a file...<br>
   flush(fileout);<br>
   try 
     {<br>
       fileout.close();<br>
   etc....

*/

public class SuWriter
{
  HashMap outputStreams;

  /**
	 Switch on? trapping of broken pipes (not elegant, but seems to work).
	 Anyone got a better idea?
  */	 
  public static final boolean catchBrokenPipes = false;

  public SuWriter()
  {
	outputStreams = new java.util.HashMap();
  }

  
  // Functions for random access file output
  public void fputtr(RandomAccessFile fp,SuSegyFormat tp)
  {
	fputtr(fp,tp,false);
  }

  public void fputtr(RandomAccessFile fp,SuSegyFormat tp,boolean allowVariableNS)
  {
	int nwritten=0;	/* bytes seen by fwrite calls	*/

	// get last input stream from map 
	FileOutSegyInfo file_infoptr = (FileOutSegyInfo) outputStreams.get(fp);
        
	if ( file_infoptr == null )
	  {		
		/* initialize new segy input stream */
		file_infoptr = new FileOutSegyInfo(fp);
		file_infoptr.itr = 0;
		// add to map
		outputStreams.put(fp,file_infoptr);

		/* Sanity check the segy header */
		file_infoptr.nsfirst = tp.get_short("ns");
		if (file_infoptr.nsfirst > SuSeismicConstants.SU_NFLTS)
		  SuErrorHandler.die_err("unable to handle " + file_infoptr.nsfirst 
			  + " > " + SuSeismicConstants.SU_NFLTS + " samples per trace");

		switch(tp.get_short("trid")) 
		  {
		  case SuSeismicConstants.CHARPACK:
			file_infoptr.bytesper = 2; break;
		  case SuSeismicConstants.SHORTPACK:
			file_infoptr.bytesper = 4; break;
		  default:
			file_infoptr.bytesper = 4; break;
		  }
	  } 

	int ns = tp.get_short("ns");
	if (!allowVariableNS)
	  { 
		// impose check
		if (ns != file_infoptr.nsfirst)
		  {
			SuErrorHandler.die_err("On trace " + ((file_infoptr.itr)+1) + 
								   ", number of samples in header " + tp.get_short("ns") + 
								   " differs from number for first trace " + file_infoptr.nsfirst);
		  }
	  }

	try 
	  {
		file_infoptr._fp.write(tp.hdr_bytes, 0, SuSeismicConstants.HDRBYTES);
	  } catch (IOException  e) {SuErrorHandler.die_err("I/O header-write exception on non-first header write");}		  
		
	datawrite(tp, file_infoptr,ns);

	file_infoptr.itr++;
  } 
  


  // Overloaded function for writing to OutputStream, as opposed to RandomAccessFile
  public void fputtr(FileOutputStream fp,SuSegyFormat tp)
  {
	fputtr(fp,tp,false);
  }

  public void fputtr(FileOutputStream fp,SuSegyFormat tp,boolean allowVariableNS)
  {
	int nwritten=0;	/* bytes seen by fwrite calls	*/
	
	// get last input stream from map 
	StreamOutSegyInfo stream_infoptr = (StreamOutSegyInfo) outputStreams.get(fp);
        
	if ( stream_infoptr == null )
	  {		
		/* initialize new segy input stream */
		stream_infoptr = new StreamOutSegyInfo(fp);
		stream_infoptr.itr = 0;
		// add to map
		outputStreams.put(fp,stream_infoptr);

		/* Sanity check the segy header */
		stream_infoptr.nsfirst = tp.get_short("ns");
		if (stream_infoptr.nsfirst > SuSeismicConstants.SU_NFLTS)
		  SuErrorHandler.die_err("unable to handle " + stream_infoptr.nsfirst 
			  + " > " + SuSeismicConstants.SU_NFLTS + " samples per trace");

		switch(tp.get_short("trid")) 
		  {
		  case SuSeismicConstants.CHARPACK:
			stream_infoptr.bytesper = 2; break;
		  case SuSeismicConstants.SHORTPACK:
			stream_infoptr.bytesper = 4; break;
		  default:
			stream_infoptr.bytesper = 4; break;
		  }
	  } 

	int ns = tp.get_short("ns");
	if (!allowVariableNS) 
	  {
		if (ns != stream_infoptr.nsfirst)
		  {
			SuErrorHandler.die_err("On trace " + ((stream_infoptr.itr)+1) + 
								   ", number of samples in header " + tp.get_short("ns") + 
								   " differs from number for first trace " + stream_infoptr.nsfirst);
		  }
	  }

	try 
	  {
		stream_infoptr.buf_fp.write(tp.hdr_bytes, 0, SuSeismicConstants.HDRBYTES);
	  } 
	catch (IOException  e) 
	  {
		if (e.getMessage().equals("Broken pipe") && catchBrokenPipes)
		  {
			SuErrorHandler.dieSilently("SuWriter: Broken pipe");
		  }
		else 
		  {
			SuErrorHandler.die_err("I/O stream header-write exception on non-first header write");		  		  
		  }
	  }

	datawrite(tp, stream_infoptr,ns);

	stream_infoptr.itr++;
  } 

  /**
	 Write data into appropriate fields of SuSegyFormat tp, using file
	 info embedded in FileOutSegyInfo
  */
  private void datawrite(SuSegyFormat tp, FileOutSegyInfo iptr)
  {
	datawrite(tp,iptr,iptr.nsfirst);
  }
  /**
	 Write data into appropriate fields of SuSegyFormat tp, using file
	 info embedded in FileOutSegyInfo
  */
  private void datawrite(SuSegyFormat tp, FileOutSegyInfo iptr,int ns)
  {
	int databytes = iptr.bytesper * ns;

	try 
	  {
		iptr._fp.write(tp.data_bytes,0,databytes);
	  } 
	catch (IOException  e) 
	  {
		SuErrorHandler.die_err("On trace " + ((iptr.itr)+1) + ", tried to write " + databytes + 
			" bytes, and failed.");
	  }	
  }

  /**
	 Write data into appropriate fields of SuSegyFormat tp, using stream
	 info embedded in StreamOutSegyInfo
  */
  private void datawrite(SuSegyFormat tp, StreamOutSegyInfo iptr)
  {
	datawrite(tp,iptr,iptr.nsfirst);
  }

  private void datawrite(SuSegyFormat tp, StreamOutSegyInfo iptr,int ns)
  {
	int databytes = iptr.bytesper * ns;
	
	try 
	  {
		iptr.buf_fp.write(tp.data_bytes,0,databytes);
	  }
	catch (IOException  e) 
	  {
		// Broken pipes occur routinely if piping into, say, SU's surange. These
		// are not fatal errors, so things ought to die quietly.
		// Can't figure out a less clunky way to trap this......
		if (e.getMessage().equals("Broken pipe") && catchBrokenPipes)
		  {
			SuErrorHandler.dieSilently("SuWriter: Broken pipe");
		  }
		else
		  {
			SuErrorHandler.die_err("On trace " + ((iptr.itr)+1) + ", tried to write " +  databytes + 
								   " bytes, and failed.");
		  }	
	  }
  }

  /**
	 Flush the buffered file output stream before closing
	 (necessary for buffered output)
  */
  public void flush(FileOutputStream fp)
  {
	// get last input stream from map 
	StreamOutSegyInfo stream_infoptr = (StreamOutSegyInfo) outputStreams.get(fp);
        
	if ( stream_infoptr != null )
	  {	 
		try 
		  {
			stream_infoptr.buf_fp.flush();
		  } 
		catch (IOException  e) 
		  {
			SuErrorHandler.die_err("I/O exception on flush() of stream");
		  }
	  }
  }
  
  /**
	 Class of useful stuff to keep in a HashMap if we are opening 
	 multiple files, for efficient alternation among files
  */

  private class FileOutSegyInfo
  {
	RandomAccessFile _fp;         /* RA file from which traces are coming */
	long itr;	                  /* number of traces write	 */
	int nsfirst;        	      /* samples from 1st header	 */
	short bytesper;               /* bytes per datum		 */
	int nsegy; 		              /* segy bytes from nsfirst	 */
	int ntr;                      /* traces in input,if known */
	
	/** 
		Constructor 
	*/
	public FileOutSegyInfo(java.io.RandomAccessFile fp)
	{
	  _fp = fp;
	}
  }
  
  /**
	 Class of useful stuff to keep in a HashMap if we are opening 
	 multiple input streams, for efficient alternation among streams
  */
   
  private class StreamOutSegyInfo
  {
	FileOutputStream _fp;          /* stream from which traces are coming */
	BufferedOutputStream buf_fp;   /* buffered wrapper of above */
	long itr;	                  /* number of traces write	 */
	int nsfirst;        	      /* samples from 1st header	 */
	short bytesper;               /* bytes per datum		 */
	int nsegy; 		              /* segy bytes from nsfirst	 */
	int ntr;                      /* traces in input,if known */

	/** 
		Constructor
		
		On initialisation, wraps input stream into a buffered writer
	*/
	public StreamOutSegyInfo(java.io.FileOutputStream fp)
	{
	  _fp = fp;
	  buf_fp = new BufferedOutputStream(_fp);
	}
  }
}
