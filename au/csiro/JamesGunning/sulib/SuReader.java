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
   Basic class to get seismic traces in su format, rather
   like the fgettr.c suite of functions in SeismicUnix.

   Given a SuSegyFormat object to hold the trace, and an
   opened RandomAccessFile (for disk files) or
   FileInputStream (for stdin), traces are stuffed
   into the SuSegyFormat by the <br><br>

   fgettr(RandomAccessFile fp, SuSegyFormat tr)<br><br>

   or<br><br>

   fgettr(FileInputStream fp, SuSegyFormat tr)<br><br>

   calls. An example usage fragment is<br><br>

   <code>
   FileInputStream stdin = null;<br>
   SuSegyFormat tr = new SuSegyFormat();<br>
   // get input stream (stdin)<br>
   stdin = (FileInputStream) new FileInputStream(FileDescriptor.in);<br>
   while (sureader.fgettr(stdin,tr) != 0)<br>
   {<br>
     <ul>// do stuff with traces...</ul>
   }<br>
   </code>
*/

public class SuReader
{
  HashMap inputStreams;

  public SuReader()
  {
	inputStreams = new java.util.HashMap();
  }

  // Functions for random access file io

  public int fgettr(RandomAccessFile fp, SuSegyFormat tp, boolean fromLE)
  {
	return(fgettr_internal(fp,tp,true,fromLE));
  }

  public int fgettr(RandomAccessFile fp, SuSegyFormat tp)
  {
	return(fgettr_internal(fp,tp,true,false));
  }

  public int fvgettr(RandomAccessFile fp, SuSegyFormat tp)
  {
	return(fgettr_internal(fp,tp,false,false));
  }

  private int fgettr_internal(RandomAccessFile fp,SuSegyFormat tp, boolean fixed_length, boolean fromLE)
  {
	int nread=0;	/* bytes seen by fread calls	*/

	// get last input stream from map 
	FileInSegyInfo file_infoptr = (FileInSegyInfo) inputStreams.get(fp);
        
	if ( file_infoptr == null )
	  {		
		/* initialize new segy input stream */
		int databytes;               /* bytes from nsfirst	*/
		
		file_infoptr = new FileInSegyInfo(fp);
		file_infoptr.itr = 0;
		file_infoptr.ntr = -1;
		// file_infoptr.ftype = FileType.DISK; 
		// add to map
		inputStreams.put(fp,file_infoptr);
		
		/* Get the header */
		
		// switch (nread = (int) efread(tp, 1, HDRBYTES, file_infoptr->buf)) 
		try 
		  {
			nread = fp.read(tp.hdr_bytes,0,SuSeismicConstants.HDRBYTES);
		  } catch (IOException  e) {SuErrorHandler.die_err("I/O exception on first header read");}    

		if (nread<=0)
		  {
			return 0; /* no traces; trap in mains */
		  }
		else 
		  {
			if (nread != SuSeismicConstants.HDRBYTES)
			  SuErrorHandler.die_err("Bad first header");
		  }		
		
		/* Have the header, now for the data */
		// flip if need be
		if (fromLE) tp.swapHeaderValues();

		file_infoptr.nsfirst = tp.get_short("ns");
		if (file_infoptr.nsfirst > SuSeismicConstants.SU_NFLTS)
		  SuErrorHandler.die_err("Unable to handle " + file_infoptr.nsfirst + " > " + 
			  SuSeismicConstants.SU_NFLTS + " samples per trace");

		switch (tp.get_short("trid")) 
		  {
		  case SuSeismicConstants.CHARPACK:
			file_infoptr.bytesper = 2; break;
		  case SuSeismicConstants.SHORTPACK:
			file_infoptr.bytesper = 4; break;
		  default:
			file_infoptr.bytesper = 4; break;
		  }

		databytes = file_infoptr.bytesper * tp.get_short("ns");

		file_infoptr.nsegy = SuSeismicConstants.HDRBYTES + databytes;

		/* Inconvenient to bump nread here; do it in the switch */
		nread = dataread(tp, file_infoptr, fixed_length);

		switch (nread) 
		  {
		  case 0:   SuErrorHandler.die_err("No data on first trace");
		  default:  
			{
			  if (nread != databytes)
				{
				  SuErrorHandler.die_err("fgettr_internal: First trace: read only " + nread + " bytes of " +  databytes);
				}
			  else nread += SuSeismicConstants.HDRBYTES;
			}
		  // since this is a file, assuming constant ns, compute ntr
		  try 
			{
			  long filelength = fp.length();
			  file_infoptr.ntr = (int)((fp.length())/((long) file_infoptr.nsegy));							   
			} catch (IOException  e) {SuErrorHandler.die_err("fgettr_internal: cannot get file length");}
		  }
		if (fromLE) tp.swapDataValues();
	  } 
	else 
	  {		/* Not first entry */
		try 
		  {
			nread = file_infoptr._fp.read(tp.hdr_bytes, 0, SuSeismicConstants.HDRBYTES);
		  } catch (IOException  e) {SuErrorHandler.die_err("fgettr_internal: I/O exception on non-first header read");}		  
		switch (nread)
		  {
		  case -1: 
			return 0; /* finished */
		  default:  
			if (nread != SuSeismicConstants.HDRBYTES)
			  SuErrorHandler.die_err("fgettr_internal: On trace " + ((file_infoptr.itr)+1) + " read " + nread + " bytes, expected "
				  + SuSeismicConstants.HDRBYTES + " bytes");
		  }
		
		nread += dataread(tp, file_infoptr, fixed_length);

		if (fromLE) tp.swapHeaderValues();

		if (fixed_length && (tp.get_short("ns") != file_infoptr.nsfirst))
		  {
			SuErrorHandler.die_err("fgettr_internal: On trace " +  ((file_infoptr.itr)+1) +
								   ", number of samples in header " + tp.get_short("ns") +
								   " differs from number for first trace " + file_infoptr.nsfirst); 
		  }
		if (fromLE) tp.swapDataValues();		
	  }
	file_infoptr.itr++;
	return (nread);
  } 
  







  // Overloaded function for reading from InputStream, as opposed to RandomAccessFile

  public int fgettr(FileInputStream fp, SuSegyFormat tp,boolean fromLE)
  {
	return(fgettr_internal(fp,tp,true,fromLE));
  }

 public int fgettr(FileInputStream fp, SuSegyFormat tp)
  {
	return(fgettr_internal(fp,tp,true,false));
  }
  
  public int fvgettr(FileInputStream fp, SuSegyFormat tp)
  {
	return(fgettr_internal(fp,tp,false,false));
  }


  /** 
	  Overloaded function for reading from InputStream, as opposed to RandomAccessFile
  */
  private int fgettr_internal(FileInputStream fp,SuSegyFormat tp, boolean fixed_length, boolean fromLE)
  {
	int nread=0,off;	/* bytes seen by fread calls	*/
	int available;

	// get last input stream from map 
	StreamInSegyInfo stream_infoptr = (StreamInSegyInfo) inputStreams.get(fp);
        
	if ( stream_infoptr == null )
	  {		
		/* initialize new segy input stream */
		int databytes;               /* bytes from nsfirst	*/
		
		stream_infoptr = new StreamInSegyInfo(fp);
		stream_infoptr.itr = 0;
		stream_infoptr.ntr = -1;

		// add to map
		inputStreams.put(fp,stream_infoptr);
		
		/* Get the header */
		try 
		  {
			off=0;
			do
			  {
				nread =  stream_infoptr.buf_fp.read(tp.hdr_bytes,off,SuSeismicConstants.HDRBYTES-off);
				off += nread;
			  } while (off < (SuSeismicConstants.HDRBYTES) && nread != -1);
			nread = off;
		  } catch (IOException  e) {SuErrorHandler.die_err("Stream: I/O exception on first header read");}
		if (nread<=0)
		  {
			return 0; /* no traces; trap in mains */
		  }
		else 
		  {
			if (nread != SuSeismicConstants.HDRBYTES)
			  SuErrorHandler.die_err("Stream: Bad first header");
		  }		
		// flip if need be
		if (fromLE) tp.swapHeaderValues();

		/* Have the header, now for the data */
		stream_infoptr.nsfirst = tp.get_short("ns");

		if (stream_infoptr.nsfirst > SuSeismicConstants.SU_NFLTS)
		  SuErrorHandler.die_err("Stream: Unable to handle " + stream_infoptr.nsfirst + " > " + 
			  SuSeismicConstants.SU_NFLTS + " samples per trace");

		switch (tp.get_short("trid")) 
		  {
		  case SuSeismicConstants.CHARPACK:
			stream_infoptr.bytesper = 2; break;
		  case SuSeismicConstants.SHORTPACK:
			stream_infoptr.bytesper = 4; break;
		  default:
			stream_infoptr.bytesper = 4; break;
		  }

		databytes = stream_infoptr.bytesper * tp.get_short("ns");

		stream_infoptr.nsegy = SuSeismicConstants.HDRBYTES + databytes;

		/* Inconvenient to bump nread here; do it in the switch */
		nread = dataread(tp, stream_infoptr, fixed_length);

		switch (nread) 
		  {
		  case 0:   SuErrorHandler.die_err("Stream: No data on first trace");
		  default:  
			{
			  if (nread != databytes)
				{
				  SuErrorHandler.die_err("Stream: First trace: read only " + nread + " bytes of " +  databytes);
				}
			  else nread += SuSeismicConstants.HDRBYTES;
			}
		  // CAN'T get ntr (no random access)....skip this bit
  		  }

		// flip if need be
		if (fromLE) tp.swapDataValues();
	  } 
	else 
	  {		/* Not first entry */
		try 
		  {
			off=0;
			do
			  {
				nread =  stream_infoptr.buf_fp.read(tp.hdr_bytes,off,SuSeismicConstants.HDRBYTES-off);
				off += nread;
			  } while (off < (SuSeismicConstants.HDRBYTES) && nread != -1);   /// FIXME! CHECKME!!!!!!!!!!
			nread = off;

		  } catch (IOException  e) {SuErrorHandler.die_err("Stream: I/O exception on non-first header read");}		  
		switch (nread)
		  {
		  case -1: 
			return 0; /* finished */
		  default:  
			if (nread != SuSeismicConstants.HDRBYTES)
			  SuErrorHandler.die_err("Stream: On trace " + ((stream_infoptr.itr)+1) + " read " + nread + " bytes, expected "
				  + SuSeismicConstants.HDRBYTES + " bytes");
		  }
		// flip if need be
		if (fromLE) tp.swapHeaderValues();

		nread += dataread(tp, stream_infoptr, fixed_length);
		
		if (fixed_length && (tp.get_short("ns") != stream_infoptr.nsfirst))
		  SuErrorHandler.die_err("Stream: On trace " +  ((stream_infoptr.itr)+1) +
			  ", number of samples in header " + tp.get_short("ns") +
			  " differs from number for first trace " + stream_infoptr.nsfirst); 
		
		// flip if need be
		if (fromLE) tp.swapDataValues();
	  }		
	stream_infoptr.itr++;
	return (nread);
  }
  
  /**
	 Read data into appropriate fields of SuSegyFormat tp, using file
	 info embedded in FileInSegyInfo
  */
  private int dataread(SuSegyFormat tp, FileInSegyInfo iptr, boolean fixed_length)
  {
	int nsread = fixed_length ? iptr.nsfirst : tp.get_short("ns");
	int databytes = iptr.bytesper*nsread;
	int nread=0;

	try 
	  {
		nread =  iptr._fp.read(tp.data_bytes,0,databytes);
	  } catch (IOException  e) {SuErrorHandler.die_err("I/O exception on first header read");}	
	
	if ((nread > 0) && (nread != databytes))
	  {
		SuErrorHandler.die_err("On trace " + ((iptr.itr)+1) + ", tried to read " + databytes + 
			" bytes, read " + nread + " bytes ");
	  }
	return(nread);
  }

  /**
	 Read data into appropriate fields of SuSegyFormat tp, using stream
	 info embedded in StreamInSegyInfo
  */
  private int dataread(SuSegyFormat tp, StreamInSegyInfo iptr, boolean fixed_length)
  {
	int nsread = fixed_length ? iptr.nsfirst : tp.get_short("ns");
	int databytes = iptr.bytesper*nsread;
	int nread=0,off=0;

	
	//if (nsread<0)	SuErrorHandler.die_err("nsread = "+nsread+ " cf "+tp.get_short("ns"));
	try 
	  {
		do
		  {
			nread =  iptr.buf_fp.read(tp.data_bytes,off,databytes-off);
			off += nread;
		  } while (off < (databytes) && nread != -1);
		nread = off;
	  } 
	catch (IOException  e) 
	  {
		SuErrorHandler.die_err("Stream: I/O exception on first header read");
	  }	
	catch (IndexOutOfBoundsException  e) 
	  {
		SuErrorHandler.die_err("Stream: IndexOutOfBoundsException, off="+off+", nread="+
							   nread+" available="+tp.data_bytes.length+" no. data desired="+(databytes-off)+
							   "iptr.bytesper ="+iptr.nsfirst+" nsread="+nsread);
	  }

	if ((nread > 0) && (nread != databytes)) 
	  {
		SuErrorHandler.die_err("Stream: On trace " + ((iptr.itr)+1) + ", tried to read " + databytes + 
			" bytes, read " + nread + " bytes ");
	  }
	return(nread);
  }

  /**
	 Class of useful stuff to keep in a HashMap if we are opening 
	 multiple files, for efficient alternation among files
  */

  private class FileInSegyInfo
  {
	RandomAccessFile _fp;         /* RA file from which traces are coming */
	long itr;	                  /* number of traces read	 */
	int nsfirst;        	      /* ns = no. samples from 1st header	 */
	short bytesper;               /* bytes per datum		 */
	int nsegy; 		              /* segy bytes from nsfirst	 */
	int ntr;                      /* traces in input,if known */
	
	/** 
		Constructor 
	*/
	public FileInSegyInfo(java.io.RandomAccessFile fp)
	{
	  _fp = fp;
	}
  }
  
  /**
	 Class of useful stuff to keep in a HashMap if we are opening 
	 multiple input streams, for efficient alternation among streams
  */
   
  private class StreamInSegyInfo
  {
	FileInputStream _fp;          /* stream from which traces are coming */
	BufferedInputStream buf_fp;   /* buffered wrapper of above */
	long itr;	                  /* number of traces read	 */
	int nsfirst;        	      /* samples from 1st header	 */
	short bytesper;               /* bytes per datum		 */
	int nsegy; 		              /* segy bytes from nsfirst	 */
	int ntr;                      /* traces in input,if known */

	/** 
		Constructor
		
		On initialisation, wraps input stream into a buffered reader
	*/
	public StreamInSegyInfo(FileInputStream fp)
	{
	  _fp = fp;
	  buf_fp = new BufferedInputStream(_fp);
	}
  }
}
