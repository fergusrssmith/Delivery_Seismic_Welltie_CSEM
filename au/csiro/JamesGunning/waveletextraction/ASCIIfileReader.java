/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 



package au.csiro.JamesGunning.waveletextraction;

import cern.colt.list.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;



import java.io.*;
import java.util.*;

/**
   Basic ASCII style ascii log parser, as example implementation.
   Contains a naive attempt to auto-detect LAS 2 files and read them.
   Unwrapped LAS mode is only supported at this stage.
   
   
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/

public class ASCIIfileReader
{
  private ErrorHandler EH;
  private ArrayList containers;
  private boolean useXYZformat = false;
  private boolean useGEOEASformat = false;
  private boolean autodetectLAS2 = false;
  private boolean tryToParseAsLAS2 = false;
  ParsingInfo parsingInfo;

  /** 
	  a general ASCII reader with error handler as provided.
	  <p>
	  (1) <b>If useGEOEASformat == true && autodetectLAS2==false</b>, 
	  will expect the format described in GSLIB:
	  <p>
	  <i>
	  title<br>
	  N<br>
	  column_title_1<br>
	  column_title_2<br>
	  column_title_3<br>
	  ...<br>
	  column_title_N<br>
	  ascii data in N columns<br>
	  </i>
	  <p>

	  (2) <b>If useGEOEASformat == false && autodetectLAS2==false</b>, 
	  a generic ASCII format is sought; column_titles can appear on one line, but they must
	  contain no white space. N is omitted, and computed from the
	  number of column titles we see. Thus we expect either 
	  <p>
	  <i>title<br> 
	  column_title_1<br> 
	  column_title_2<br>
	  column_title_3<br> 
	  ...<br> 
	  column_title_N<br> 
	  ascii data in N columns<br> 
	  </i>
	  <p> 
	  OR 
	  <p>
	  <i> title<br>
	  column_title_1 column_title_2 column_title_3 ... column_title_N<br> 
	  ascii data in N columns<br> 
	  </i>
	  <p>
	  
	  (3) <b>If useXYZformat == true</b>, 
	  Simple 3-column header-less ascii format
	  <p>
	  <i>
	  ascii data in 3 columns<br> 
	  </i>
	  <p> 

	  (4) <b>If autodetectLAS2==true</b>, a naive attempt to autodetect LAS format is made. The first line
	  will contain ~V, ~VERSION, or ~version, and if this is detected, we switch to LAS format:
	  see http://www.cwls.org/las_info.htm
	  <p>
	  ~CURVE defines a section with column titles as the first words, possibly followed by other stuff
	  ~A or ~ASCII then starts the line preceding the actual columns of unwrapped data.

	  It looks roughly like: <br>
	  <pre>
	  ~VERSION INFORMATION
	  VERS.                 2.0  : CWLS LOG ASCII STANDARD -VERSION 2.0
	  WRAP.                 NO   : ONE LINE PER DEPTH STEP
	  ~WELL INFORMATION
	  #MNEM.UNIT              DATA        DESCRIPTION OF MNEMONIC
	  #----- -----          ----------    --------------------------
	  STRT .M           1670.0000          :START DEPTH
	  STOP .M           1669.7500          :STOP DEPTH
	  STEP .M           -0.1250            :STEP
	  NULL .            -999.25            :NULL VALUE
	  COMP .    ANY OIL COMPANY INC.       :COMPANY
	  WELL .    ANY ET AL 12-34-12-34      :WELL
	  FLD  .    WILDCAT                    :FIELD
	  LOC  .    12-34-12-34W5              :LOCATION
	  PROV .    ALBERTA                    :PROVINCE
	  SRVC .    ANY LOGGING COMPANY INC.   :SERVICE COMPANY
	  DATE .    13-DEC-86                  :LOG DATE
	  UWI  .    100123401234W500           :UNIQUE WELL ID
	  ~CURVE INFORMATION
	  #MNEM.UNIT          API CODES      CURVE DESCRIPTION
	  #-----------      ---------       -------------------------
	  DEPT .M                           :  1  DEPTH
	  DT   .US/M         60 520 32 00   :  2  SONIC TRANSIT TIME
	  RHOB .K/M3         45 350 01 00   :  3  BULK DENSITY
	  NPHI .V/V          42 890 00 00   :  4   NEUTRON POROSITY
	  SFLU .OHMM         07 220 04 00   :  5  RXORESISTIVITY
	  SFLA .OHMM         07 222 01 00   :  6  SHALLOW RESISTIVITY
	  ILM  .OHMM         07 120 44 00   :  7  MEDIUM RESISTIVITY
	  ILD  .OHMM         07 120 46 00   :  8  DEEP RESISTIVITY
	  ~PARAMETER INFORMATION
	  #MNEM.UNIT          VALUE            DESCRIPTION
	  #---------       ---------          ----------------------
	  MUD  .            GEL CHEM            :   MUD TYPE
	  BHT  .DEGC        35.5000             :   BOTTOM HOLE TEMPERATURE
	  BS   .MM          200.0000            :   BIT SIZE
	  FD   .K/M3        1000.0000           :   FLUID DENSITY
	  MATR .            SAND                :   NEUTRON MATRIX
	  MDEN .            2710.0000           :   LOGGING MATRIX DENSITY
	  RMF  .OHMM        0.2160              :   MUD FILTRATE RESISTIVITY
	  DFD  .K/M3        1525.0000           :   DRILL FLUID DENSITY
	  ~OTHER
	  Note: The logging tools became stuck at 625 metres
	  causing the data between 625 metres and 615 metres to be
	  invalid.
	  ~A  DEPTH    DT    RHOB       NPHI     SFLU    SFLA     ILM      ILD
	  1670.000  123.450 2550.000   0.450  123.450  123.450  110.200  105.600
	  1669.875  123.450 2550.000   0.450  123.450  123.450  110.200  105.600
	  1669.750  123.450 2550.000   0.450  123.450  123.450  110.200  105.600
	  </pre>
  */


  /** 
	  Sample driver to show first few lines of a named file after parsing
  */
  public static void main(String argv[])
  {
	ErrorHandler EH = new ErrorHandler();
	String filename = null;
	ASCIIfileReader rdr = new ASCIIfileReader(EH,true,true);

	OptionsParser options = new OptionsParser("ASCIIfileReader",EH);

   	options.addMandatoryArg(1,"filename","ascii file to parse");	
	options.addMandatoryArg(2,"property1,property2,...","supply properties to extract as a comma-separated list");

	options.parseOptions(argv,true);

	filename = options.getMandatoryArg(1);
	String words = options.getMandatoryArg(2);


	// Set up reader receptacles
	ArrayList props = new ArrayList();
	ArrayList numbers = new ArrayList();
	StringTokenizer st = new StringTokenizer(words,",:;\n\r");
	DoubleArrayList list=null;
	while (st.hasMoreTokens())
	  {
		String word = st.nextToken();
		props.add(word);
		list = new DoubleArrayList();
		numbers.add(list);
		rdr.addToContainerList(word,list);
	  }
	// now get the data
	rdr.extractData(filename);

	// Show some sample output
	int i,j;
	int length = Math.min(4,list.size());
	// int length = list.size();
	
	// column titles
	for (i=0;i<props.size();i++)
	  {
		EH.printf("%-14s ",(String)props.get(i));
	  }
	EH.println("");
	// data
	for (j=0;j<length;j++)
	  {
		for (i=0;i<props.size();i++)
		  {
			EH.printf("%-14g ",((DoubleArrayList)numbers.get(i)).get(j));
		  }
		EH.println("");
	  }
  }

  public ASCIIfileReader(ErrorHandler EH, boolean useGEOEASformat, boolean autodetectLAS2)
  {
	this.EH = EH;
	this.useGEOEASformat = useGEOEASformat;
	this.autodetectLAS2 = autodetectLAS2;
	containers = new ArrayList();
  }

  public ASCIIfileReader(ErrorHandler EH, boolean useXYZformat)
  {
	this.EH = EH;
	this.useGEOEASformat = false;
	this.autodetectLAS2 = false;
	this.useXYZformat = true;
	containers = new ArrayList();
  }

  /**
	 Method to help assemble the list of things required. This will 
	 need a column header word, and an array of double for the data
  */
  void addToContainerList(String word, DoubleArrayList receptable)
  {
	containers.add(new dataContainer(word,receptable));
  }

  /**
	 Method to help assemble the list of things required. This will 
	 need a column header word, and an array of double for the data
  */
  void addToContainerList(String word, DoubleArrayList receptable,boolean mandatory)
  {
	containers.add(new dataContainer(word,receptable,mandatory));
  }

  /**
	 Open the named file and stuff all the required data (assuming the header words have been correctly
	 added into the containers) methods into the container receptacles.
  */
  public boolean extractData(String filename)
  {
	boolean success=true;
	boolean seekVersHeaderLAS=true;
	boolean seekCurveHeaderLAS=false;
	boolean seekCurveInfoLAS=false;
	boolean seekAsciiHeaderLAS=false;
	boolean seekAsciiInfoLAS=false;
	int i,j,nColumns=0;
	boolean haveContinued = false;
	double val=0;
	DoubleArrayList data;
	dataContainer dC;
	String s=null,word,unitWord=null;
	StringTokenizer st;
	// array of column numbers for desired properties
	int[] columnNumber = new int[containers.size()];
	boolean[] fieldUsed=null;
	parsingInfo = new ParsingInfo();

	try 
	  {
        BufferedReader in = new BufferedReader(new FileReader(filename));
		//BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(filename), "UTF-8"));
		int lineNo=0;
		int field=1;
		int startData=3; // can't be earlier than this in either format

		while ((s=in.readLine()) != null)
		  {
			// strip leading & trailing white space
			s = s.trim();
			lineNo++;

			// lineNo==1 is title: look for LAS keywords
			// if (lineNo==1 && autodetectLAS2). But there might be a # comment
			if (seekVersHeaderLAS && autodetectLAS2)
			  {
				if (s.startsWith("#")) 
				  {
					// skip comments
					continue;
				  }
				// look for magic first word ~V, ~VERSION, ~Version
				st = new StringTokenizer(s);
				if (st.hasMoreTokens())
				  {
					word = st.nextToken();
					if (word.equals("~v") || word.equalsIgnoreCase("~version"))
					  {
						// yes it looks like LAS
						tryToParseAsLAS2 = true;
						parsingInfo.isLAS2 = true;
						// place to store units if found
						parsingInfo.unitsMap = new HashMap<java.lang.String,java.lang.String>();

						// now seek for curve header
						seekVersHeaderLAS=false;
						seekCurveHeaderLAS = true;
						haveContinued = true;
					  }
				  }
			  }
			if (tryToParseAsLAS2)
			  {
				// first, look for ~CURVE stuff
				if (seekCurveHeaderLAS)
				  {
					// look for words
					st = new StringTokenizer(s);
					if (st.hasMoreTokens())
					  {
						word = st.nextToken();
						if (word.equals("~C") || word.equalsIgnoreCase("~curve"))
						  {
							// OK we found the curve section
							seekCurveHeaderLAS = false;
							// now seek for curve info
							seekCurveInfoLAS = true;
							nColumns=0;
							haveContinued = false;
							continue;
						  }
					  }
				  }
				if (seekCurveInfoLAS)
				  {
					// now curve info
					// look for words. Some LAS files do funky things like WORD.UNIT,so we'll split at dots too
					// st = new StringTokenizer(s,". \t\r\n");
					String[] sa = s.split("\\s+");     // split just at white space
					//					if (st.hasMoreTokens())
					if (sa.length>0)
					  {
						// extract first word
						// word = st.nextToken();
						// Some LAS files do funky things like WORD.UNIT,so we'll split at dots here, to detect units
						String[] word1 = sa[0].split("\\."); 
						if (word1.length==0) EH.die_err("LAS file " + filename+" has wierd curve word starting with .");

						word = word1[0];
						if (word.startsWith("#") || word.equals("")) 
						  {
							// skip comments && cope with wierd blank line if required..
							haveContinued = false;
							continue;
						  }
						else if (word.startsWith("~"))
						  {
							// ~ Signals the end of curve section: ~ starts something new....
							seekCurveInfoLAS = false;
							// now look for ASCII header line
							seekAsciiHeaderLAS = true;
							// But before going on, check whether everything we need is here
							// Loop over column titles
							for (i=0;i<containers.size();i++)
							  {
								if (columnNumber[i]==0)
								  {
									// no column for this word
									dC = (dataContainer)(containers.get(i));
									if (dC.mandatory)
									  {
										EH.die_err("LAS file " + filename+" has no "+ dC.word + " data");
									  }
									else columnNumber[i]=-1; // flags non-necessity
								  }
							  }

							// Boolean array for parsing speedups
							fieldUsed = new boolean[nColumns+1];
							Arrays.fill(fieldUsed,false);
							for (i=0;i<containers.size();i++)
							  {
								if (columnNumber[i]>0) fieldUsed[columnNumber[i]] = true;
							  }
							// EH.print("Fields: ");
							// for (i=1;i<=nColumns;i++) EH.print(" "+fieldUsed[i]);
							// EH.println("");

							haveContinued = true;
						  }
						else
						  {
							// match column titles against first word of a line that doesn't start with ~ or #		
							nColumns++;
							for (i=0;i<containers.size();i++)
							  {
								dC = (dataContainer)(containers.get(i));
								if (word.equals(dC.word)) 
								  {
									columnNumber[i] = nColumns;
									// try to get units if they follow	
									char[] uarr;
									if (word1.length>1)
									  {
										// line started with "label.unit"
										unitWord = word1[1];
									  }	
									else if (word1.length==1 && (sa.length>1))
									  { 
										// line started with "label .unit"
										// so sa[1] = ".unit"
										unitWord = sa[1].substring(1);
									  }
									//unitWord.replaceAll("\u00b5", "u"); // Greek mu -> u. Rather harder than this, sadly...
									unitWord = unitWord.toLowerCase();
									parsingInfo.unitsMap.put(word,unitWord);
								  }
							  }							
						  }
					  }
				  }
				if (seekAsciiHeaderLAS)
				  {
					// now look for ascii header
					// look for words
					st = new StringTokenizer(s);
					if (st.hasMoreTokens())
					  {
						word = st.nextToken();
						if (word.equalsIgnoreCase("~a") || word.equalsIgnoreCase("~ascii")) 
						  {
							// ok start of ascii section
							seekAsciiHeaderLAS = false;
							// seek for ascii data
							seekAsciiInfoLAS = true;
							haveContinued = false;
							continue;
						  }
					  }
				  }
				if (seekAsciiInfoLAS)
				  {
					// now look for ascii data
					// actual data
					st = new StringTokenizer(s);
					int col=1;
					//EH.println("col_field: s="+s);
					while (st.hasMoreTokens()) 
					  {
						//EH.println("col_field "+col+" "+field);
						word = st.nextToken();
						if (col==1 && word.charAt(0) == '~')
						  {
							// Must have skipped into next section. jump out
							seekAsciiInfoLAS = false;
							break;
						  }
 						if (fieldUsed[field])
						  {
							// i.e. don't bother trying to understand columns we don't need
							try
							  {
								val = Double.parseDouble(word);
							  }
							catch (NumberFormatException e)
							  {
								EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
							  }
							// add data to lists
							for (i=0;i<containers.size();i++)
							  {
								// if (col == columnNumber[i]) // will skip optional quantities with columnNumber[i]=-1
								if (field == columnNumber[i]) // will skip optional quantities with columnNumber[i]=-1
								  {
									dC = (dataContainer)(containers.get(i));
									dC.data.add(val);
								  }
							  }
						  }
						col++;
						field++;
					  }
					if (field>nColumns) field=1;
				  }
				haveContinued = false;
			  } // EO LAS PARSING
			else
			  {
				if (useGEOEASformat)
				  {
					if (lineNo==2)
					  {
						// no of columns
						try
						  {
							nColumns= Integer.parseInt(s);
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Second line of geoEAS file " + filename+" not a number");
						  }
						startData = nColumns+3;
					  }
					else if ((lineNo > 2) && (lineNo < startData))		 
					  {
						// column titles
						for (i=0;i<containers.size();i++)
						  {
							dC = (dataContainer)(containers.get(i));
							if (s.equals(dC.word)) columnNumber[i] = lineNo-2;
						  }
					  }
					else if (lineNo == startData)
					  {
						// test all columns needed
						// column titles
						for (i=0;i<containers.size();i++)
						  {
							if (columnNumber[i]==0)
							  {
								// no column for this word
								dC = (dataContainer)(containers.get(i));
								if (dC.mandatory)
								  {
									EH.die_err("geoEAS file " + filename+" has no "+ dC.word + " data");
								  }
								else columnNumber[i]=-1; // flags non-necessity
							  }
						  }
					  }
				  }
				else if (!useXYZformat)
				  { 
					// Plain ascii with headers
					if (lineNo==1) nColumns=0;
					if ((lineNo > 1) && (lineNo < startData))
					  {
						// look for words
						st = new StringTokenizer(s);
						while (st.hasMoreTokens()) 
						  {
							boolean isANumber = true;
							word = st.nextToken();
							try
							  {
								val = Double.parseDouble(word);
							  }
							catch (NumberFormatException e)
							  {
								// looks like a word
								nColumns++;
								// match it up
								for (i=0;i<containers.size();i++)
								  {
									dC = (dataContainer)(containers.get(i));
									if (word.equals(dC.word)) columnNumber[i] = nColumns;
								  }
								isANumber = false;
							  }
							if (isANumber) 
							  {
								startData = lineNo; // so next if () will start extracting numbers							
								break;
							  }
						  }
					  }
					if (lineNo == startData)
					  {
						// Just before we break out, test all columns needed
						// column titles
						for (i=0;i<containers.size();i++)
						  {
							if (columnNumber[i]==0)
							  {
								// no column for this word
								dC = (dataContainer)(containers.get(i));
								if (dC.mandatory)
								  {
									EH.die_err("ASCII file " + filename+" has no "+ dC.word + " data");
								  }
								else columnNumber[i]=-1; // flags non-necessity
							  }
						  }
					  }
				  }
				else
				  {
					// Plain ascii 3 column useXYZformat
					if (lineNo==1) 
					  {
						nColumns=3;
						startData = lineNo;
						for(i=0;i<3;i++) columnNumber[i]=(i+1);
					  }
				  }
				// data reading common to both formats
				if (lineNo >= startData)
				  {
					// actual data
					st = new StringTokenizer(s);
					int col=1;
					while (st.hasMoreTokens()) 
					  {
						try
						  {
							val = Double.parseDouble(st.nextToken());
						  }
						catch (NumberFormatException e)
						  {
							EH.die_err("Line "+lineNo+" of file " + filename+" has non-numbers");
						  }
						// add data to lists
						for (i=0;i<containers.size();i++)
						  {
							if (col == columnNumber[i])  // will skip optional quantities with columnNumber[i]=-1
							  {
								dC = (dataContainer)(containers.get(i));
								dC.data.add(val);
							  }
						  }
						col++;
					  }
				  }
			  }
		  }
		in.close();
	  }
	catch (Exception e) 
	  {
		EH.die_err("Can't open or decipher file "+filename);
	  }

	// trim containers
	for (i=0;i<containers.size();i++)
	  {
		((dataContainer)(containers.get(i))).data.trimToSize();
	  }

	return success;
  }
	
  public static class dataContainer
  {
	/** key word for data
	 */
	String word;
	/**
	   list of data under this keyword
	*/
	DoubleArrayList data;
	/**
	   flag to indicate whether mandatory or not
	*/
	boolean mandatory = true;
	public dataContainer(String word,DoubleArrayList data)
	{
	  this.word = word;
	  this.data = data;
	}

	public dataContainer(String word,DoubleArrayList data,boolean mandatory)
	{
	  this.word = word;
	  this.data = data;
	  this.mandatory = mandatory;
	}
  } 

  /**
	 Class of useful info to return about parsing information
  */
  public static class ParsingInfo
  {
	boolean isLAS2 = false;
	HashMap <String,String>unitsMap;
  }

}
