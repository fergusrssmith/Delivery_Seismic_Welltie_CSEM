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


package au.csiro.JamesGunning.lib;
import java.io.*;
 
/**						
   Simple debugging oriented error handler
   Various flavours of error handler can be put here 
   
   @author james.gunning@csiro.au
   @version 2.0, 9/7/2006
 */
public class ErrorHandler
{
  private PrintWriter out_w,err_w;

  /**
	 ErrorHandler object with customised output and error streams
	 Autoflush on by default
  */
	
  public ErrorHandler(FileDescriptor fd_out,FileDescriptor fd_err)
  {
	out_w = new PrintWriter(new FileOutputStream(fd_out),true); 
	err_w = new PrintWriter(new FileOutputStream(fd_err),true); 
  }

  /**
	 ErrorHandler object with default output and error streams
	 Autoflush switched on.
  */
	
  public ErrorHandler()
  {
	out_w = new PrintWriter(new FileOutputStream(FileDescriptor.out),true);
	err_w = new PrintWriter(new FileOutputStream(FileDescriptor.err),true);
  }
	
  /**
	 Abort wrapper
  */
  public void die()
  {
	// .printStackTrace();
	System.exit(0);
  }

  /**
	 Error handler that calls system exit()

	 @param mesg
	 Spiel to user on exit
  */
  public void die_err(String mesg)
  {
	err_w.println(mesg);

	// printStackTrace();
	System.exit(0);
  }

  public void die_err(String mesg, Object... args)
  {
	err_w.printf(mesg,args);
	System.exit(0);
  }

  public void dieWithStackTrace(String mesg)
  {
	err_w.println(mesg);
	System.exit(0);
  }

  /** 
	  Generic warning message
  */
  public void warn(String mesg)
  {
	err_w.println(mesg);		
  }

  public void warn(String mesg, Object... args)
  {
	err_w.printf(mesg,args);
  }

  /**
	 General benign message, for general usage
  */
  public void print(String mesg)
  {
	out_w.printf(mesg);
  }

  /**
	 General benign varargs message, for general usage
  */
  public void printf(String mesg, Object... args)
  {
	out_w.printf(mesg,args);
  }


  /**
	 General benign message with linefeed, for general usage
  */
  public void println(String mesg)
  {
	out_w.printf(mesg+"%n");
  }
  /**
	 General benign varargs message with linefeed, for general usage
  */
  public void println(String mesg, Object... args)
  {
	out_w.printf(mesg+"%n",args);
  }

  public void handleIOException(java.io.IOException ex)
  {
	ex.printStackTrace();
	System.exit(0);
  }
	
  public void handleGeneralException(Exception ex)
  {
	ex.printStackTrace();
	System.exit(0);
  }
}
