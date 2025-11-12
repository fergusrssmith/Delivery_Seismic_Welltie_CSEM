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

package au.csiro.JamesGunning.lib.util;

import java.io.File;

public final class JNIPathUtil
{
 /*
	For library libname, returns, prefix/[add_lib/]libname.suffix

	e.g.

	Linux-i686-unknown/libname.so

	where the prefix, suffix .so and separator are platform determined.

	If add_lib==null, intermediate string is omitted, e.g.
	
	Linux-i686-unknown/libname.so

	Current choices of the platforms and suffixes are

	OS                JVMMODEL   suffix    prefix                            Remarks
	-----------------------------------------------------------------------------------------------------------
	Linux-i* or -x*   32         .so       Linux-i686-unknown
	Linux-x86_64      64         .so       Linux-x86_64-unknown
	OSX               32,64      .dylib    Darwin-Power-Macintosh-unknown   ought to be built as universal libs
	OSX               32,64      .dylib    Darwin-i386-unknown              ""
	windows/cygwin    32         .dll      cygwin                           mingw builds are useful

  */
  public static String getJNILibraryPath(String libName,String add_lib)
  {
	String model,fullLibName;
	String fullPath=null;
	String OS,arch,gcj;

	OS = System.getProperty("os.name", "?").toLowerCase();   // Usually e.g. "linux" or "mac os x"
	model = System.getProperty("sun.arch.data.model", "?");  // usually 32 or 64
	arch = System.getProperty("os.arch", "?").toLowerCase(); // e.g. "i386" or "ppc*"
	gcj = System.getProperty("gnu.classpath.vm.shortname","?").toLowerCase(); // returns libgcj if on gcj

	boolean found = true;

	if (gcj.equals("libgcj"))
	  {
		if (model.equals("?") && arch.startsWith("i")) // catches i*86 I think on libgcj
		  {
			model = "32";
		  }
	  }

	if (add_lib==null)
	  {
		fullLibName = libName;
	  }
	else
	  {
		fullLibName = add_lib + File.separatorChar + libName;
	  }
	if (OS.contains("linux")) 
	  {
		if (model.equals("32"))
		  {
			fullPath = "Linux-i686-unknown/"+fullLibName+".so";
		  }
		else
		  {
			fullPath = "Linux-x86_64-unknown/"+fullLibName+".so";
		  }
	  }
	else if (OS.startsWith("mac os x") || OS.startsWith("darwin")) 
	  {
		if (arch.startsWith("ppc"))
		  {
			fullPath = "Darwin-Power-Macintosh-unknown/"+fullLibName+".jnilib";
		  }
		else if (model.equals("32"))
		  {
			fullPath = "Darwin-i386-unknown/"+fullLibName+".jnilib";
		  }
		else if (model.equals("64"))
		  {
			fullPath = "Darwin-x86_64-unknown/"+fullLibName+".jnilib";
		  }
		else {found=false;}
		/*
		  if (model.equals("64")) System.err.println("No 64 bit libs available on OSX yet. Try setting your java to 32 bit. "+
		  "os.name="+OS+" os.arch="+arch+" sun.arch.data.model="+model);		
		*/
	  }
	else if (OS.startsWith("windows"))
	  {
		fullPath = "cygwin/"+fullLibName+".dll";
	  }
	else
	  {
		found=false;
	  }
	if (!found) System.err.println("Failed to locate native jni library: check uname and rebuild? "+
								   "os.name="+OS+" os.arch="+arch+" sun.arch.data.model="+model);		

	//System.err.println("USING: os.name="+OS+" os.arch="+arch+" sun.arch.data.model="+model);
	//System.err.println("Trying PATH="+fullPath);	

	return fullPath;
  }
}