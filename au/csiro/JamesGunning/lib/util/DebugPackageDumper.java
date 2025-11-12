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


package au.csiro.JamesGunning.lib.util;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.sulib.SuSegyFormat;

import java.io.*;
import java.util.*;
import java.util.zip.*;

/**
   Basic class for assembling files into a zip package for mailing to 
   developers, for reproducing bugs etc. 

   @author James.Gunning@csiro.au
   @version 1.0, 28/1/2006 
*/

public class DebugPackageDumper
{
  /**
	 Buffer size for zip creation
  */
  public static final int BUFFER = 2048;

  private String[] accessibleOrdinaryFiles=null;
  private String[] accessibleSuFiles=null;
  private SuTraceInfo[] traceInfo=null;
  private String commandline=null, filename=null;
  private Exception ex=null;
  private ErrorHandler _EH;
  /**
	 Constructor allocs room for header and data using constants as
	 per seismic constants 

	 @param EH ErrorHandler for i/o feedback
	 @param accessibleOrdinaryFiles . List of absolute filenames of ordinary files to be stuffed into
	 debug package. May be null.
	 @param accessibleSuFiles List of absolute filenames of Su files to be stuffed into
	 debug package. May be null.
	 @param traceInfo Array of in-place pairs of {String name,SuSegyFormat traces} traces to be
	 dumped. The trace has the actual data, the name is taken as a guide for naming (prefixes stripped etc)
	 May be null.
	 @param commandline String of invocation of program.
	 @param ex  Exception leading to the call of this class. May also be null 
	 if a snapshot if required.
  */
  public DebugPackageDumper(ErrorHandler EH,
							String[] accessibleOrdinaryFiles, String[] accessibleSuFiles,
							SuTraceInfo[] traceInfo, String commandline, Exception ex)
  {
	this._EH = EH;
	this.ex = ex;
	this.accessibleOrdinaryFiles = accessibleOrdinaryFiles;
	this.accessibleSuFiles = accessibleSuFiles;
	this.traceInfo = traceInfo;
	this.commandline = commandline;
  }

  /**
	 After construction, assemble all the files into a zip files
	 with a magic name <it>debug_package.suffix.zip</it> 
	 with suffix derived from the current time.
  */

  public void dumpPackage()
  {
	int i=0;
	File dir;
	String dirname,cwd,sufilename=null;
	InputStream s_in;
	OutputStream s_out;
	boolean abort=false;

	// make up some unique looking directory from current time
	Calendar cal = new GregorianCalendar();
	String suffix = cal.getTime().toString();
	// spaces and : entries are bad news for filenames
	suffix = suffix.replaceAll("[ :]", "_");

	_EH.println("");
	_EH.println("Making debug package...");
	dirname = new String("debug_package."+suffix);
	dir = new File(dirname);

	// make directory
	if (!dir.mkdir())
	  {
		_EH.println("Could not create debug package. Giving up...");
		return;
	  }


	// (1) General commandline and stack trace
	PrintWriter out;
	String cmdlinefilename = new String(dirname+File.separatorChar+"cmdline_and_stacktrace.txt");

	try 
	  {
		out = new PrintWriter(new FileWriter(cmdlinefilename));
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " + cmdlinefilename);
	  }
	if (commandline != null)
	  {
		out.println("Comandline args:");	
		out.println("  "+commandline);	
		out.println();	
	  }
	if (ex != null)
	  {
		out.println("Stacktrace:");	
		out.println(ex.toString());
		StackTraceElement[] exElements = ex.getStackTrace();
		for (i=0;i<exElements.length;i++)
		  {
			out.println(exElements[i].toString());
		  }
	  }
	out.close();
             
	// (2) Ordinary files                   
	cwd = dir.getPath();

	if (accessibleOrdinaryFiles != null)
	  {
		try 
		  {
			for (i=0; i<accessibleOrdinaryFiles.length; i++)
			  {
				s_in = new FileInputStream(accessibleOrdinaryFiles[i]);
				// strip any leading path to get target filename
				String[] spl = accessibleOrdinaryFiles[i].split(""+File.separatorChar);
				filename = spl[spl.length-1];
				s_out = new FileOutputStream(cwd+File.separatorChar+filename);
				_EH.println("Copying "+accessibleOrdinaryFiles[i]+" to "+ cwd+File.separatorChar+filename);
				FileUtil.copyStream(s_in,s_out);
			  }
		  }
		catch (FileNotFoundException e)
		  {
			_EH.println("Could not find file "+accessibleOrdinaryFiles[i]);
			abort=true;
		  }
		catch (SecurityException e)
		  {
			_EH.println("Couldn't (security access) read "+accessibleOrdinaryFiles[i]);
			abort=true;
		  }
		catch (Exception e)
		  {
			_EH.println("Could not copy "+accessibleOrdinaryFiles[i]+" to "+
						filename+". Aborting dumpPackage()");
			abort=true;
		  }
	  }

	// (3) Accessible su files
	if (!abort && accessibleSuFiles != null)
	  {
		for (i=0; i<accessibleSuFiles.length && !abort; i++)
		  {
			try 
			  {
				s_in = new FileInputStream(accessibleSuFiles[i]);
				// strip any leading path to get target filename
				String[] spl = accessibleSuFiles[i].split(""+File.separatorChar);
				filename = spl[spl.length-1];
				s_out = new FileOutputStream(cwd+File.separatorChar+filename);
				_EH.println("Copying "+accessibleSuFiles[i]+" to "+ cwd+File.separatorChar+filename);
				FileUtil.copyStream(s_in,s_out);
			  }
			catch (FileNotFoundException e)
			  {
				_EH.println("Could not find file "+accessibleSuFiles[i]);
				abort=true;
			  }
			catch (SecurityException e)
			  {
				_EH.println("Couldn't (security access) read "+accessibleSuFiles[i]);
				abort=true;
			  }
			catch (Exception e)
			  {
				_EH.println("Could not copy "+accessibleSuFiles[i]+" to "+
							filename+". Aborting dumpPackage()");
				abort=true;
			  }

		  }
	  }

	// (4) Inplace su traces. Dump into directory with supplied names

	if (!abort && traceInfo != null)
	  {
		FileOutputStream fout=null;
		SuWriter suwriter = new SuWriter();

		for (i=0; i<traceInfo.length && !abort; i++)
		  {
			// strip any leading path to get target filename
			String[] spl = traceInfo[i].name.split(""+File.separatorChar);
			sufilename = new String(cwd+File.separatorChar+spl[spl.length-1]);
			_EH.println("Dumping trace from "+traceInfo[i].name+" to "+ sufilename);

			try 
			  {
				fout = (FileOutputStream) new FileOutputStream(sufilename);
				// write out trace
				suwriter.fputtr(fout, traceInfo[i].tr);
				suwriter.flush(fout);
				fout.close();
			  } 
			catch (IOException e) 
			  {
				// throw new RuntimeException("couldn't open and write to " + sufilename);
				_EH.println("couldn't open and write to " + sufilename);
				abort = true;
			  }
		  }
	  }

	// (5) final zip , and clean up
	// Hacked from http://java.sun.com/developer/technicalArticles/Programming/compression/
	try 
	  {
		// get a list of files from current directory
		File f = new File(dirname);
		if (!abort)
		  {
			BufferedInputStream origin = null;
			FileOutputStream dest = new FileOutputStream(dirname+".zip");
			ZipOutputStream z_out = new ZipOutputStream(new BufferedOutputStream(dest));
			// z_out.setMethod(ZipOutputStream.DEFLATED);
			byte data[] = new byte[BUFFER];
			String files[] = f.list();
		
			for (i=0; i<files.length; i++) 
			  {
				filename = cwd+File.separatorChar+files[i];
				_EH.println("Zip file: adding: "+filename);
				FileInputStream fi = new FileInputStream(filename);
				origin = new BufferedInputStream(fi, BUFFER);
				ZipEntry entry = new ZipEntry(files[i]);
				z_out.putNextEntry(entry);
				int count;
				while((count = origin.read(data, 0, 
										   BUFFER)) != -1) 
				  {
					z_out.write(data, 0, count);
				  }
				origin.close();
			  }
			z_out.close();
			filename = dirname+".zip";
			_EH.println("Wrote zip archive "+filename);
		  }
		else
		  {
			_EH.println("Aborting debug package formation. If it's a file missing/access issue, "
						+ "fix your problem and try again!");
		  }
		// Cleanup
		File[] fFiles = f.listFiles();
		for (i=0; i<fFiles.length; i++) 
		  {
			if (!abort) _EH.println("Deleting "+fFiles[i].toString());
			fFiles[i].delete();
		  }
		fFiles = f.listFiles();
		_EH.println("Deleting "+f.toString());
		if (!f.delete())
		  {
			_EH.println("Warning!, could not clean up and delete directory "+f.toString());
		  } 
	  }
	catch(Exception e) 
	  {
		_EH.println("Couldn't write zip archive "+dirname+".zip");
		e.printStackTrace();
      }
	if (!abort)
	  {
		_EH.println("");		
		_EH.println("\nPlease email debug zip package " +filename+ " to developers!\n");	
	  }
  }

  

  /** 
	  Simple {name,trace} data class for su trace dumps
  */
  public static class SuTraceInfo
  {
	String name;
	SuSegyFormat tr;
	
	/**
	Creates pair with given name and trace.   
	@param name suggested name of trace
	@param tr Actual trace data in as a SuSegyFormat instance
	*/
	public SuTraceInfo(String name, SuSegyFormat tr)
	{
	  this.name = name;
	  this.tr = tr;
	}
  }
}




