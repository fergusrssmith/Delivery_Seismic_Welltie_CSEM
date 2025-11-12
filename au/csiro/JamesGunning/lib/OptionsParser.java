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


package au.csiro.JamesGunning.lib;

import java.util.*;

/**
   Class for parsing set of commandline options to java executable
*/
public class OptionsParser
{
  /** 
	  Error Handler for instance of this class
  */
  ErrorHandler EH;

  //  private TreeMap options = null;
  private LinkedHashMap options = null; // preserves insertion order
  private TreeMap mandatoryargs = null;
  private String callingprog;
  private String extraArgs="";
  private String spiel="";
  private String version=null;
  private String license_text=null;
  private StringBuffer argb;
  private int hP=0; // help priority. Zero by default (highest)
  private boolean showHelpExtras=false;
  private String extraHelpOptionString=null;

  public OptionsParser(String callingprog, ErrorHandler EH)
  {
	mandatoryargs = new TreeMap();
	options = new LinkedHashMap();
	this.callingprog = callingprog;
	this.EH = EH;
  }
  
  public boolean parseOptions(String[] args, boolean matchExactly)
  {
	int i,j,k=0;
	Collection opts = options.values();

	/* first, restore any comma protected arguments that argv has split, eg 
	   the phrase "green rock" 
	*/
	if (args.length > 1)
	  {
		boolean clean=false;
		String[] nargs;
		ArrayList nargv;
		nargv = new ArrayList();
		for (i=0; i< args.length; i++)
		  {
			if (args[i].startsWith("\""))
			  {
				String s = args[i].substring(1);
				for (j=i+1 ; j< args.length ; j++)
				  {
					if (args[j].endsWith("\"")) 
					  {
						clean=true;
						i=j;
						s = s.concat(" "+args[j].substring(0,args[j].length()-1));
						break;
					  }
					else
					  {
						s = s.concat(" "+args[j]);
					  }
				  }
				nargv.add(s);
			  }
			else
			  {
				nargv.add(args[i]);
			  }
		  }
		if (clean)
		  {
			args = (String[])nargv.toArray(new String[nargv.size()]);
		  }
		else
		  {
			nargv=null;
		  }
	  }


	// See if args list requested extra help somewhere
	if (extraHelpOptionString!=null)
	  {
		for (i=0 ; i< args.length ; i++)
		  {
			if (args[i].equals(extraHelpOptionString)) showHelpExtras = true;
		  }
	  }
	
			

	argb = new StringBuffer();
	for (i=0 ; i< args.length ; i++) argb.append(args[i]).append(" ");

	for (i=0 ; i< args.length ; i++)
	  {
		if (options.containsKey(args[i]))
		  {
		   OptionDetails optvalues = (OptionDetails) options.get(args[i]);
		   optvalues.is_set=true;
		   
		   if (i + optvalues.n_args >= args.length)
			 {
			   EH.warn("Too few arguments for option " + args[i]);
			   Usage();
			 }
		   for (j=0; j<optvalues.n_args; j++)
			 {
			   try
				 {
				   if (optvalues.type.equals("float"))
					 {
					   optvalues.values[j]=new Float(Float.parseFloat(args[1+i+j]));
					 }
				   else if (optvalues.type.equals("int"))
					 {
					   optvalues.values[j]=new Integer(Integer.parseInt(args[1+i+j]));
					 }
				   else if (optvalues.type.equals("String"))
					 {
					   optvalues.values[j]=new String(args[1+i+j]);
					 }
				   else if (optvalues.type.equals(null))
					 {
					   optvalues.values[j]=null;
					 }
				   else
					 {
					   // should never get here...
					   EH.die_err("Illegal type in : parseOptions()" + optvalues.type);
					 }
				 }
			   catch (NumberFormatException nfe)
				 {
				   EH.warn("Incorrect arguments for option " + args[i]);
				   Usage();
				 }
			 }
		   i += optvalues.n_args;
		  }
		else
		  {
			// looks like a mandatory arg...
			// "-" alone is commonly used for stdin, but "-something" could be dodgy
			if (!args[i].equals("-") && args[i].startsWith("-") && matchExactly) 
			  {
				EH.println("Unrecognized (likely) option " + args[i]+" See usage:%n");
				Usage();
				return false;
			  }
			if (matchExactly)
			  {
				if (k >= mandatoryargs.size())
				  {
					EH.warn("Too many args: " + args[i]);
					Usage();
				  }
				MandatoryDetails mdets = (MandatoryDetails) mandatoryargs.get(new Integer(k+1));
				mdets.value=args[i];
			  }
			else
			  {
				// accept here any extra arguments without objection.
				MandatoryDetails mdets = (MandatoryDetails) mandatoryargs.get(new Integer(k+1));
				// Expand HashMap if necessary 
				if (mdets==null)
				  {
					mdets = new MandatoryDetails(k+1,null,null);
					mandatoryargs.put(new Integer(k+1),mdets);
				  }
				mdets.value = args[i];
			  }
			// count up through number of mandatory args we have seen
			k++;
		  }
	  }		
	if (k != mandatoryargs.size()) 
	  {
		Usage();
		return false;
	  }
	return true;
  }

  /**
	 Return a space-separated string representatiion of the commandline; available after parsing
	 with parseOptions
  */
  public String getFullArgsString()
  {
	return argb.toString();
  }

/**
	 Return a String array of all non-optional arguments
  */
  public String[] getAllMandatoryArgs()
  {
	String[] SA = new String[mandatoryargs.size()];
	for (int i=0 ; i< mandatoryargs.size() ; i++)
	  {
		SA[i] = ((MandatoryDetails) mandatoryargs.get(new Integer(i+1))).value;
	  }
	return SA;
  }

  /**
	 The Usage() method will print a line of form<br>
	 Program_name [options] extraArgs<br>
	 where the extraArgs text can be set by this method.<br>
	 It's default is ""
  */
  public void addMainCommandlineText(String text)
  {
	extraArgs = text;
  }

  /**
	 The Usage() method will print a line of form<br>
	 Program_name [options] extraArgs<br>
	    options summary<br>
	    spiel<br>
	 where the spiel text can be set by this method.<br>
	 It's default is ""
  */
  public void addSpielText(String text)
  {
	spiel = text;
  }

  /**
	 Add version string to the self doc usage call
  */
  public void addVersion(String text)
  {
	version = text;
  }
  /**
	 Add license string to the self doc usage call
  */
  public void addLicense(String text)
  {
	license_text = text;
  }

 public void addMandatoryArg(int arg_no, String name, String description)
  {
   	mandatoryargs.put(new Integer(arg_no),new MandatoryDetails(arg_no,name,description));
  }

  public void addExtraHelpOption(String short_key, String description)
  {
	extraHelpOptionString = short_key;
	addOption(short_key,null,0,null,null,description);
  }
  /**
	 Set a bunch of these to specify program options...
	 e.g.<br>
	 addOption("-a","--addgnus",2,"int","add two gnus to the collection");<br>

	 this creates space (an int array, int[2]) for a commandline option of form

	 --addgnus 3 4<br><br>
	 or<br><br>
	 -a 3 4<br>
  */

  public void addOption(String short_key, String long_key, int n_args, String type, String description)
  {
	if (short_key==null && long_key==null)
	  EH.die_err("addOption requires at least one short or long option");
	else
	  {
		if (short_key!=null && long_key!=null)
		  {
			options.put(short_key,new OptionDetails(short_key,long_key,n_args,type,description,hP));
		  }
		else if (short_key!=null)
		  {
			options.put(short_key,new OptionDetails(short_key,null,n_args,type,description,hP));
		  }
		else if (long_key!=null)
		  {
			options.put(short_key,new OptionDetails(long_key,null,n_args,type,description,hP));
		  }
	  }
  }
  
  /**
	 As above, but with helpPriority set as well. Use helpPriority>=0
  */
  public void addOption(String short_key, String long_key, int n_args, String type, String description, int helpPriority)
  {
	hP = helpPriority;
	addOption(short_key, long_key, n_args, type, description);
	hP=0; // restore highest level
  }


  public void addOption(String short_key, String long_key, int n_args, String type, String args_list, String description) 
  {
	if (short_key==null && long_key==null)
	  EH.die_err("addOption requires at least one short or long option");
	else
	  {
		if (short_key!=null && long_key!=null)
		  {
			OptionDetails details = new OptionDetails(short_key,long_key,n_args,type,args_list,description,hP);
			options.put(short_key,details);
			options.put(long_key,details);
		  }
		else if (short_key!=null)
		  {
			options.put(short_key,new OptionDetails(short_key,null,n_args,type,args_list,description,hP));
		  }
		else if (long_key!=null)
		  {
			options.put(long_key,new OptionDetails(long_key,null,n_args,type,args_list,description,hP));
		  }
	  }
  }

  /**
	 As above, but with possible help priority. Use helpPriority>=0
  */
  public void addOption(String short_key, String long_key, int n_args, String type, String args_list, String description, int helpPriority) 
  {
	hP = helpPriority;
	addOption(short_key, long_key, n_args, type, args_list, description);
	hP=0; // restore highest level
  }

  // function to retrieve stuff about options

  public boolean isOptionSet(String key)
  {
	if (((OptionDetails)options.get(key)).is_set==true)
	  return true;
	else
	  return false;
  }
  

  // return array of floats, ints, strings etc corresponding to the named option...
  public Object[] getOptionValues(String key)
  {
	if (options.get(key)==null)
	  {
		EH.die_err("No such option in getOptionValues(); " + key);
		return null;
	  }
	else
	  return ((OptionDetails) options.get(key)).values;
  }

  // return int for option argument n
  public int getOptionInt(String key,int n)
  {
	return ((Integer)getOptionValues(key)[n]).intValue();
  }

  // return float for option argument n
  public float getOptionFloat(String key,int n)
  {
	return ((Float)getOptionValues(key)[n]).floatValue();
  }
  
  // return String for option argument n
  public String getOptionString(String key,int n)
  {
	return ((String)getOptionValues(key)[n]);
  }

  // return the numbered mandatory argument
  public String getMandatoryArg(int arg_no)
  {
	MandatoryDetails details = (MandatoryDetails) mandatoryargs.get(new Integer(arg_no));
	if (details==null)
	  {
		EH.die_err("No such mandatory arg in getMandatoryArg(); " + (arg_no));
		return null;
	  }
	else
	  return (details.value);
  }
  
  // return the number of mandatory args given
  public int getNumberMandatoryArgsSupplied()
  {
	return mandatoryargs.size();
  }
 
  public void Usage()
  {
	int pass;
	StringTokenizer args=null;

	if (version!=null) 
	  {
		EH.println(callingprog+", version "+version);
		EH.println("");
	  }
	if (license_text != null)
	  {
		EH.printf(license_text+"%n");
	  }
	EH.println("Usage:");
	EH.println("");


	// mandatory bits;
	Collection c = mandatoryargs.values();
	EH.print(callingprog + " ");
	for (Iterator i = c.iterator(); i.hasNext();) 
	  {
		MandatoryDetails values = (MandatoryDetails)i.next();
		if(values.name!=null)
		  {
			EH.print(values.name);
			EH.print(" ");
		  }
	  }
	EH.println("[options] "+extraArgs);
	for (Iterator i = c.iterator(); i.hasNext();) 
	  {
		MandatoryDetails values = (MandatoryDetails)i.next();
		if(values.name!=null) EH.print("  " + values.name);
		if(values.description!=null)
		  {
			EH.print(": ");
			EH.print(makeIndentedString(values.description,true));
		  }
		EH.println("");
	  }

	c = null;

	// list of options
	EH.println("");
	EH.println("  Options:");
	Collection o = options.values();

	// If showHelpExtras, list extras at the end, as per "menu supplements"

	int duplicateCount=0;
	for (pass=0; pass<=(showHelpExtras?1:0); pass++)
	  {
		if (showHelpExtras && pass==1)
		  {
			EH.println("");
			EH.println("  Extra Options:");
		  }
		for (Iterator i = o.iterator(); i.hasNext();) 
		  {
			OptionDetails values = (OptionDetails)i.next();

			if (showHelpExtras)
			  {
				if (pass==0 && values.helpPriority!=0) continue; // skip nonessential options on pass 0
				if (pass==1 && values.helpPriority==0) continue; // skip essential options on pass 1
			  }
			else
			  {
				if (values.helpPriority!=0) continue; // skip nonessential options on pass 0
			  }
			
			// Don't print info for both long and short keys...
			if (values.key!=null && values.alt_key!=null) duplicateCount++; // has both long & short forms. Mark.
			else duplicateCount=0;
			if (duplicateCount==2) continue;

			
			// Start formatting...
			EH.print("    ");
			EH.print(values.key);
			if (values.alt_key!=null) EH.print(", " + values.alt_key);
			if (values.arg_string!=null)
			  {
				args = new  StringTokenizer(values.arg_string);
			  }
			for (int j=1;j<= values.n_args;j++) 
			  {
				if (values.arg_string==null)
				  {
					EH.print(" arg" + j);
				  }
				else
				  {				
					if (args.hasMoreTokens())
					  {					
						EH.print(" " + args.nextToken());
					  }
				  }
			  }
				
			if (values.type!=null) EH.print(" (" + values.type + ((values.n_args>1)?"s)":")"));
			EH.println("");EH.print("      ");
			EH.println(makeIndentedString(values.description,true));
		  }
	  }
	o = null;

	// finally the spiel
	EH.println("");
	EH.println(makeIndentedString(spiel,false));

	// and die
	EH.die_err("");
  }

  // method to make sure values.description descriptions do not wrap around 
  // badly when self-docs are generated
  // Words should start with %n to force a new line in description string

  private String makeIndentedString(String text, boolean useIndent)
  {
	int NCOLUMNS = 80;         // width to format inside
	int col,indentCol;
	StringTokenizer st = new StringTokenizer(text);
	StringBuffer sb = new StringBuffer(0);
	String word, iString="";
	String linespace = "%n      ";   // this is 6 column indentation
	String space = " ";
	boolean forcedLF=false;
	indentCol = space.length();
	col = indentCol;
	while (st.hasMoreTokens()) 
	   {
         word = st.nextToken();
		 forcedLF = word.startsWith("%n");
		 if (((col+word.length())< NCOLUMNS) && !forcedLF)
		   {
			 // concatenate tokens and space while within 80 cols
			 sb.append(word).append(space);
			 col += word.length()+1;
			 if (word.endsWith("%n"))
			   {
				 col = indentCol;
			   }
		   }
		 else
		   {
			 // add in an extra indentation 
			 if (forcedLF) word = word.substring(2);
			 if (useIndent)
			   {
				 sb.append(linespace).append(word).append(space);
			   }
			 else
			   {
				 sb.append("%n").append(word).append(space);
			   }
			 col=indentCol;
		   }
	   }
	return sb.toString();
  }



  private class OptionDetails
  {
	boolean is_set=false;
	int n_args=0;
	String key=null;
	String alt_key=null;
	String arg_string = null;
	String type=null;
	String description=null;
	Object[] values=null;
	int helpPriority=0;

	/**
	   Simple constructor with default argument names
	*/
	private  OptionDetails(String key, String alt_key, int n_args, String type,String description,int hP)
	{
	  this.n_args = n_args;
	  this.description = description;
	  this.type = type;
	  this.key = key;
	  this.alt_key = alt_key;
	  this.helpPriority = hP;

	  if (type==null)
		{
		  this.values = null;
		}
	  else
		{
		  if (type.equals("float"))
			{
			  this.values = (Object []) new Float[n_args];
			}
		  else if (type.equals("int"))
			{
			  this.values = (Object []) new Integer[n_args];
			}
		  else if (type.equals("String"))
			{
			  this.values = (Object []) new String[n_args];
			}
		  else
			{
			  EH.die_err("Illegal type in OptionDetails: " + type);
			}
		}
	}

	/**
	   Simple constructor with explicit argument names
	*/

	private  OptionDetails(String key, String alt_key, int n_args, String type,String arg_string, String description,int hP)
	{
	  this.n_args = n_args;
	  this.description = description;
	  this.type = type;
	  this.key = key;
	  this.alt_key = alt_key;
	  this.arg_string = arg_string;
	  this.helpPriority = hP;

	  if (type==null)
		{
		  this.values = null;
		}
	  else
		{
		  if (type.equals("float"))
			{
			  this.values = (Object []) new Float[n_args];
			}
		  else if (type.equals("int"))
			{
			  this.values = (Object []) new Integer[n_args];
			}
		  else if (type.equals("String"))
			{
			  this.values = (Object []) new String[n_args];
			}
		  else
			{
			  EH.die_err("Illegal type in OptionDetails: " + type);
			}
		}
	}
  }


  private class MandatoryDetails implements Comparable 
  {
	int arg_no;
	String name;
	String description;
	String value;

	private  MandatoryDetails(int arg_no,String name, String description)
	{
	  this.arg_no = arg_no;
	  this.name = name;
	  this.description = description;
	}

	/**
	   Comparator to sort by argument number
	*/
	public int  compareTo(Object o)
	{
	  MandatoryDetails n = (MandatoryDetails)o;
	  int cmp = (arg_no<n.arg_no)? -1 : ((arg_no>n.arg_no)? 1 : 0 );
	  return cmp;
	}
  }

  public static class OptionsComparator implements Comparator<String>
  {

	public int compare(String o1, String o2) 
	{
	  int order;
	  order = o1.compareTo(o2); 	  // natural String order
	  // make exception for ???
	  return order;
	}
	
	public boolean equals(Object o) 
	{
	  return false;
	}
  }
}		 

  
