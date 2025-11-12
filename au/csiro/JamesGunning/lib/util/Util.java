package au.csiro.JamesGunning.lib.util;

import au.csiro.JamesGunning.lib.ErrorHandler;
import java.io.*;
import java.util.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import au.csiro.JamesGunning.lib.Numerical.DenseFloatMatrix2D;

import java.lang.Math.*;

/**
   Basic class of static public function for use in debugging arrays, some 
   sundry String utils etc etc.

   @author james.gunning@dpr.csiro.au
   @version 2.0, 9/8/2008 
*/

public class Util
{
  /* 
	 Utility for halting execution till user responds..
  */
  public static void getStdinUserPrompt()
  {
	System.out.print("Press return to continue...");

	//  open up standard input
	BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
	//  read the from the command-line
	try 
	  {
		String s = br.readLine();
	  } 
	catch (IOException ioe) 
	  {
	  System.err.println("IO error reading from stdin");
	  System.exit(1);
	}
  }

  /**
	 GetPID. Returns -1 if exception occurs. Works on *nix flavours I believe
	 Tested on OSX, Linux.
  */
  public static int getPID()
  {
    int pid;
	try 
	  {
		String pids = java.lang.management.ManagementFactory.getRuntimeMXBean().getName();
		String[]Ids = pids.split("@");		
		pid = Integer.valueOf(Ids[0]);
		return pid;
	  }
	catch (Exception e1) 
	  {
		System.err.println("Can't get PID");
	  }
	return -1;
  }


  /**
	 Write script file gdbjava for gdb debugging of current process (with user execute permissions) 
  */
  public static void writeGDBJavaScript()
  {
	String filename="gdbjava";
	int pid=0;
	try 
	  {
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
		pid = Util.getPID();
		out.printf("gdb -q -p "+pid);
		out.close();
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't write " + filename + " script" );
	  }
	try
	  {
		Runtime.getRuntime().exec("chmod u+x "+filename);
	  }
	catch (IOException ioe) 
	  {
		System.err.println("chmod: couldn't set u+x " + filename);
	  }
	System.err.println("Wrote: "+filename+" (PID="+pid+")");
  }


  /** 
	  cheap method for checking if function is integer
  */
  public static boolean isInteger(String s)
  {
	if (s.matches("((-|\\+)?[0-9])+")) return true;
	else return false;
  }

  /** 
	  Parse doubles in c or fortran form
  */
  public static double parseDouble(String st)
  {
	if (st.contains("d")) 
	  {
		st = st.replace('d','e');
	  }
	else if (st.contains("D")) 
	  {
		st = st.replace('D','e');
	  }
	return Double.parseDouble(st);
  }

  /** 
	  Parse floats in c or fortran form
  */
  public static float parseFloat(String st)
  {
	if (st.contains("d")) 
	  {
		st = st.replace('d','e');
	  }
	else if (st.contains("D")) 
	  {
		st = st.replace('D','e');
	  }
	return Float.parseFloat(st);
  }

  /** 
	  Returns trailing basename of a full file path name, e.g.
	  fileBasename("/home/mystuff/files/foo")
	  returns foo
  */
  public static String fileBasename(String filename)
  {
	int c = filename.lastIndexOf(File.separatorChar); 
	if (c != -1)
	  {
		return filename.substring(c+1);
	  }
	else return new String(filename);
  }

  /** 
	  Returns leading part of filename sans suffix: suffix defined as fragment after '.' , e.g.
	  fileSansSuffix("foo.eps")
	  returns foo
  */
  public static String fileSansSuffix(String filename)
  {
	int c = filename.lastIndexOf('.');
	if (c != -1)
	  {
		return filename.substring(0,c);
	  }
	else return new String(filename);
  }

  /** 
	  Returns leading part of filename sans suffix and leading part: suffix defined as fragment after '.' , e.g.
	  fileSansSuffix("../../foo.eps")
	  returns foo
  */
  public static String fileBasenameSansSuffix(String filename)
  {
	String s = fileBasename(filename);
	return fileSansSuffix(s);
  }

  /**
	 return true if word is in the varargs array words
  */
  public static boolean matches(String word, String... words)
  {
	for (String w : words)
	  {
		if (word.equals(w)) return true;
	  }
	return false;
  }
  

  public static void print(double[] data, int n)
  {
	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		System.out.print("{");
		for (i=0;i<n;i++)
		  {  
			System.out.print(data[i]);
			if (i<n) System.out.print(", ");
		  }
		System.out.println("}");
	  }
  }

 public static void print(float[] data, int n)
  {
	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		System.out.print("{");
		for (i=0;i<n;i++)
		  {  
			System.out.print(data[i]);
			if (i<n) System.out.print(", ");
		  }
		System.out.println("}");
	  }
  }
 public static void print(int[] data, int n)
  {
	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		System.out.print("{");
		for (i=0;i<n;i++)
		  {  
			System.out.print(data[i]);
			if (i<n) System.out.print(", ");
		  }
		System.out.println("}");
	  }
  }
  /** 
	  Print range of data[] inclusive
  */
  public static void print(double[] data, int nFirst, int nLast)
  {
	int i,n;
	System.out.print("{");
	n = Math.min(data.length-1,nLast);
	for (i=Math.max(0,nFirst);i<=n;i++)
	  {  
		System.out.print(data[i]);
		if (i<n) System.out.print(", ");
	  }
	System.out.println("}");
  }


  public static void print(boolean[] data, int n)
  {
	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		System.out.print("{");
		for (i=0;i<n;i++)
		  {  
			System.out.print(data[i]);
			if (i<n) System.out.print(", ");
		  }
		System.out.println("}");
	  }
  }

  public static void printfileZO(double[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File("tmp.txt");
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }
  /**
	 print sequence of values x_i y_i\n , i=0...(n-1), to named file
  */
 public static void printfileZO(String filename, double[] dataX, double[] dataY, int n)
  {

	int i;
	if (dataX.length<n) 
	  {
		System.out.println("dataX has length " + dataX.length + " < " + n);
		return;
	  }
	else if (dataY.length<n) 
	  {
		System.out.println("dataY has length " + dataY.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(dataX[i] + " " + dataY[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  public static void printfileZO(String filename,double[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) 
		  {
			System.out.println("Couldn't open "+ filename);
		  }
	  }
  }

 public static void printfileZO(String filename,short[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

public static void printfileZO(String filename,int[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  public static void printfileZO(String filename,float[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=0;i<n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  public static void printfile(double[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File("tmp.txt");
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=1;i<=n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  public static void printfile(int[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File("tmp.txt");
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=1;i<=n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  public static void printfile(String filename,double[] data, int n)
  {

	int i;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			for (i=1;i<=n;i++)
			  {  
				out.println(data[i]);
			  }
			out.close();
		  }
		catch (IOException e) {System.out.println("Couldn't open tmp.txt");};
	  }
  }

  /** 
	  print matrix data[1...n][1...m] to file
  */
  public static void debugPrintf77Matrix(String filename,double[][] data, int n,int m)
  {

	int i,j;
	if (data.length<n) 
	  {
		System.out.println("data has length " + data.length + " < " + n);
		return;
	  }
	else
	  {
		try
		  {
			File f=new File(filename);
			ErrorHandler EH = new ErrorHandler();
			EH.println("Printing to : "+filename);
			PrintWriter out = new PrintWriter(new FileWriter(f));
			for (i=1;i<=n;i++)
			  {  
				for (j=1;j<=m;j++)
				  { 
					out.println(data[i][j]);
				  }
			  }
			out.close();
		  }
		catch (IOException e) 
		  {
			System.out.println("Couldn't open "+ filename);
		  }
	  }
  }

  public static void debugPrintMatrix(String name,DoubleMatrix2D A)
  {
	int i,j,n;
	ErrorHandler EH;

	EH = new ErrorHandler();

	n = A.rows();
	EH.printf("%s%n",name);
	for (i = 0;i<n;i++)
	  {
		for (j = 0;j < A.columns();j++)
		  {			
			EH.printf("%12.6g ",A.getQuick(i,j));
		  }
		EH.printf("%n");
	  }
  }
  
  public static void debugPrintMatrix(String name,DoubleMatrix1D A)
  {
	int i,j,n; 
	ErrorHandler EH;
	
	EH = new ErrorHandler();
	
	n = A.size();
	EH.printf("%s%n",name);
	
	for (i = 0;i<n;i++)
	  {
		EH.printf("%12.6g%n",A.getQuick(i));
	  }
  }

 public static void debugPrintMatrix(String name,DoubleMatrix1D A, int nn)
  {
	int i,j,n;
	ErrorHandler EH;
	
	EH = new ErrorHandler();
	
	n = A.size();
	EH.printf("%s%n",name);
	
	for (i = 0;i<Math.min(n,nn);i++)
	  {
		EH.printf("%16.10g%n",A.getQuick(i));
	  }
  } 

  public static void debugPrintMatrixToFile(String filename,DoubleMatrix2D A)
  {
	int i,j,n;
	ErrorHandler EH;

	try
	  {
		File f=new File(filename);
		FileOutputStream strm =  new FileOutputStream(f);
		EH = new ErrorHandler(strm.getFD(),FileDescriptor.err);
		ErrorHandler EHstdout = new ErrorHandler();
		EHstdout.println("Printing to : "+filename);
		n = A.rows();

		for (i = 0;i<n;i++)
		  {
			for (j = 0;j < A.columns();j++)
			  {			
				// EH.printf("%18.12g ",A.getQuick(i,j));
				EH.printf("%16g ",A.getQuick(i,j));
			  }
			EH.printf("%n");
		  }

		strm.close();
	  }
	catch (IOException e) {System.out.println("Couldn't open "+filename);};
  }

  public static void debugPrintMatrixToFile(String filename,DenseFloatMatrix2D A)
  {
	int i,j,n;
	ErrorHandler EH;

	try
	  {
		File f=new File(filename);
		FileOutputStream strm =  new FileOutputStream(f);
		EH = new ErrorHandler(strm.getFD(),FileDescriptor.err);
		ErrorHandler EHstdout = new ErrorHandler();
		EHstdout.println("Printing to : "+filename);
		n = A.rows();

		for (i = 0;i<n;i++)
		  {
			for (j = 0;j < A.columns();j++)
			  {			
				// EH.printf("%18.12g ",A.getQuick(i,j));
				EH.printf("%16g ",A.getQuick(i,j));
			  }
			EH.printf("%n");
		  }

		strm.close();
	  }
	catch (IOException e) {System.out.println("Couldn't open "+filename);};
  }

  public static void debugPrintMatrixToFile(String filename,DoubleMatrix3D A)
  {
	int i,j,k;
	ErrorHandler EH;

	try
	  {
		File f=new File(filename);
		FileOutputStream strm =  new FileOutputStream(f);
		EH = new ErrorHandler(strm.getFD(),FileDescriptor.err);
		ErrorHandler EHstdout = new ErrorHandler();
		EHstdout.println("Printing to : "+filename);

		for (i = 0;i<A.slices();i++)
		  {
			for (j = 0;j < A.rows();j++)
			  {			
				for (k=0; k < A.columns(); k++)
				  {
					EH.printf("%12.6g%n",A.getQuick(i,j,k));
				  }
			  }
		  }

		strm.close();
	  }
	catch (IOException e) {System.out.println("Couldn't open "+filename);};
  }

  public static void debugPrintMatrixToFile(String filename,DoubleMatrix1D A)
  {
	int i,j,n;
	ErrorHandler EH;

	try
	  {
		File f=new File(filename);
		FileOutputStream strm =  new FileOutputStream(f);
		EH = new ErrorHandler(strm.getFD(),FileDescriptor.err);		
		ErrorHandler EHstdout = new ErrorHandler();
		EHstdout.println("Printing to : "+filename);

		n = A.size();

		for (i = 0;i<n;i++)
		  {
			EH.printf("%12.6g%n",A.getQuick(i));
		  }

		strm.close();
	  }
	catch (IOException e) {System.out.println("Couldn't open "+filename);};
  }
}




