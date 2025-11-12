package au.csiro.JamesGunning.lib.util;

import java.io.*;

/**
   Basic class of file utility functions

   @author James.Gunning@csiro.au
   @version 1.0, 1/2/2006 
*/

public class FileUtil
{
  public static int BUFFER_SIZE = 4096; //4KB buffer

  /**
	 Method for safely copying an input stream to an output stream
  */
  public static void copyStream(InputStream in, OutputStream out) throws IOException
  {	
    // Do not allow other threads to read from input or write to output while copying	
    synchronized (in)  
	  {
      synchronized (out) 
		{		  
		  byte[] buffer = new byte[BUFFER_SIZE];  //4KB buffer
		  while(true)
			{
			  int bytesRead = in.read(buffer);

			  if (bytesRead==-1)
				break;
			  out.write(buffer, 0, bytesRead);			  
			}
		}
	  }
  }
}
