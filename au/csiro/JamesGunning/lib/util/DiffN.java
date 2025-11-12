package au.csiro.JamesGunning.lib.util;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.OptionsParser;

public class DiffN
{

  String filename1 = null;
  String filename2 = null;
  ErrorHandler EH = new ErrorHandler();
  double threshold=1e-14,floor=0;
  int n=1;
  boolean colNorm=false;
  int colMin,colMax;

  private int n_count;

  public static void main(String argv[]) 
  {
	int i,j;
	DiffN D = new DiffN();

	OptionsParser options = new OptionsParser("DiffN",D.EH);
	options.addSpielText("Numerical diff code: checks if file1 and file2 match line-by-line, with all "
						 +"numbers which are tokens matching up to relative precision epsilon (default "+D.threshold+")."
						 +"Stops at first difference, unless -n supplied");
	options.addMandatoryArg(1,"file1.txt","file1");
	options.addMandatoryArg(2,"file2.txt","file2");

	options.addOption("-e",null,1,"float","eps","specify relative epsilon (default "+D.threshold+")");
	options.addOption("-f",null,1,"float","floor","(default "+D.floor+") specify noise floor (if a relative difference > epsilon, but |numbers| < floor, we don't care)");
	options.addOption("-n",null,1,"int","n","show first n difference instances (default 1)");
	options.addOption("-c",null,1,"String","colMin,colMax","compare relative differences using scale from maximum |entry| in named column range, e.g. -c 1,3");

	// parse options
	options.parseOptions(argv,true);
	if (options.isOptionSet("-e"))
	  {
		D.threshold  = (double) abs(((Float)options.getOptionValues("-e")[0]).floatValue());
	  }	
	if (options.isOptionSet("-f"))
	  {
		D.floor = (double) abs(((Float)options.getOptionValues("-f")[0]).floatValue());
	  }	
	if (options.isOptionSet("-n"))
	  {
		D.n  = ((Integer)options.getOptionValues("-n")[0]).intValue();
	  }

	if (options.isOptionSet("-c")) 
	  {
		D.colNorm = true;
		String arg = (String)options.getOptionValues("-c")[0];
		String[] st = arg.split("[,;:]");
		int[] columns = new int[st.length];
		for (i=0; i< st.length; i++)
		  {
			try
			  {
				j = Integer.parseInt(st[i]);
			  }
			catch (NumberFormatException e)
			  {
				throw new RuntimeException("-t "+arg+" contains a non-number");
			  }
			if (j <= 0)
			  {
				throw new RuntimeException("-t "+columns+" contains no-positive type");
			  }
			columns[i] = j;
		  }
		D.colMin =	columns[0];
		D.colMax = columns[columns.length-1];
	  }

	D.filename1 = options.getMandatoryArg(1);
	D.filename2 = options.getMandatoryArg(2);

	D.run();
  }

  void run()
  {
	int i,j,k;
	double x1,x2,xMax,norm=0;
	BufferedReader in1=null,in2=null;
	boolean differs=false;
	String s1,s2;
	String[] st1,st2;

	try
	  {
		in1 = new BufferedReader(new FileReader(filename1));
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't find or use file "+filename1);
	  }	
	try
	  {
		in2 = new BufferedReader(new FileReader(filename2));
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't find or use file "+filename2);
	  }
	i=1;
	n_count=0;
	try
	  {
		while ((s1 = in1.readLine()) != null)
		  {
			s1 = s1.trim();
			if ((s2 = in2.readLine()) == null)
			  {
				EH.die_err("File "+filename2+" ends at line "+i+" of file "
						   +filename1+", whose line is:\n "+s1);
			  }
			s2 =s2.trim();
			//		if (s1.startsWith("#DEPTHS")) depths_header = s; 
			st1 = s1.split("\\s+");			
			st2 = s2.split("\\s+");			
			if (st1.length!=st2.length)
			  {
				EH.println("Line "+i+" lines unequal # of tokens:");
				showlines(s1,s2);
			  }
			else
			  {
				if (colNorm)
				  {
					norm=0;
					for(j=0,k=1;j<st1.length;j++,k++)
					  {
						try 
						  {
							x1 = Double.parseDouble(st1[j]);
							if (k>=colMin && k<=colMax) norm = max(norm,abs(x1));
						  }
						catch (NumberFormatException ne1) {}
					  }
				  }			   
				for(j=0,k=1;j<st1.length;j++,k++)
				  {
					try 
					  {
						x1 = Double.parseDouble(st1[j]);
						try 
						  {
							x2 = Double.parseDouble(st2[j]);
						  
							// OK, here, both are numbers
							if (x1!=0.0)
							  {
								if (x2!=0.0)
								  {
									xMax = max(abs(x1),abs(x2));
									differs=false;
									if (colNorm)
									  {
										// comparing "in-group" column
										if ((k>=colMin && k<=colMax) && (norm>0) && ((x1-x2)/norm > threshold))// comparing "in-group" column
										  {
											differs=true;
										  }
									  }
									else if ((x1-x2)/xMax > threshold) differs=true;
									if (differs && xMax>floor)
									  {
										EH.println("Tokens "+(j+1)+" disagree, ("+x1+", "+x2+"), line "+i);
										showlines(s1,s2);
									  }
								  }
							  }
						  }
						catch (NumberFormatException ne1)
						  {
							EH.println("Numbers not present on matching tokens, line "+i);
							showlines(st1[j],st2[j]);
						  }
					  }
					catch (NumberFormatException ne1) 
					  {
						try 
						  {
							x2 = Double.parseDouble(st2[j]);
							EH.println("Numbers not present on matching tokens, line "+i);
							showlines(st1[j],st2[j]);
						  }
						catch (NumberFormatException ne2) {}
					  }
				  }
			  }		  
			i++;
		  }
		in1.close();
		in2.close();
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't read files???");
	  }
  }

  void showlines(String s1,String s2)
	{
	  EH.println("< "+s1);
	  EH.println("> "+s2);
	  n_count++;
	  if (n_count==n) System.exit(0);
	}
}
