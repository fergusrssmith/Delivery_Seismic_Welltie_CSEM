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
package script;

import org.w3c.dom.*;
import java.net.URL;
import java.util.Hashtable;
import java.util.HashSet;
import java.util.Vector;
import java.util.Stack;
import java.util.regex.*;
import java.util.Arrays;
import java.util.Enumeration;
import java.io.*;
import model.*;
import javax.swing.tree.DefaultTreeModel;

/**
  * Class which generates and writes out script files, for linear and 
  * bsub jobs.
  */
public final class ScriptGenerator {
    // types of jobs in the array
    private static final int SOLO_TYPE = 0;
    private static final int START_TYPE = 1;
    private static final int MIDDLE_TYPE = 2;
    private static final int END_TYPE = 3;
    
    /* key -> String for a command name, 
     * value -> HashSet containing Strings of commonly handled command line  
     *          options. */
    private final Hashtable commandOpts = new Hashtable();
    // state variable which tracks the number of scripts currently being
    // generated
    private int numScripts = 1;
    // state variable which holds the last type of command encountered
    private int lastType = SOLO_TYPE;
    // character to use for piping (switched briefly to '>' for P_file)
    private char pipeChar = '|';
    // table of all parameters defined for the script and their values
    // key -> parameter name, value -> parameter value
    private Hashtable parameterTable;
    // a stack of all the current Scatters being employed
    private Stack scatters = new Stack();
    // all warnings gathered over the course of the script generation
    private Vector warnings = new Vector();
    // base of all files before .sh, .bsub.sh, or [index info].[index].sh
    private String fileBase = null;
    // end of file base, the start of all job names for bsub
    private String scriptBase = null;
    // String which holds the text of the script which will be the script to
    // run to run the job linearly and not on lsf / bsub
    private String masterScriptLinear = null;
    // String which holds the text of the script which will run with lsf / bsub
    private String masterScriptBsub = null;
    private int depth = 0;
    private String baseParamText = "";
    // one-level commands (unix commands or script text)
    private HashSet shallowStarts;
    private HashSet shallowMiddles;
    private HashSet shallowEnds;
    private HashSet shallowSolos;
    // full commands of arbitrary heirarchy
    private HashSet fullStarts;
    private HashSet fullMiddles;
    private HashSet fullEnds;
    private HashSet fullSolos;

    private DEditor editor;  // editor (if any) this generator is running from
    // needed to see if running as an applet.
    private String xmlDir = null; // location of xml file that the script is
    // generating.

    /**
      * Constructor for ScriptGenerators.  Doesn't do much, except sets up
      * the commandOpts Hashtable appropriately.
      */
    public ScriptGenerator(DEditor owner) {
	this.editor = owner;
	String[] temp = new String[]{"unix_command_P","script_text_P"};
	shallowStarts = new HashSet(Arrays.asList(temp));
	temp = new String[]{"P_unix_command_P","P_script_text_P"};
	shallowMiddles = new HashSet(Arrays.asList(temp));
	temp = new String[]{"P_unix_command","P_script_text"};
	shallowEnds = new HashSet(Arrays.asList(temp));
	temp = new String[]{"cp","rm","unix_command","script_text"};
	shallowSolos = new HashSet(Arrays.asList(temp));
	temp = new String[]{"segyread","bhpread"};
	fullStarts = new HashSet(Arrays.asList(temp));
	temp = new String[]{"sushw","bhploadhdr","bhpintegrate","bhpstorehdr",
			    "bhphorizon","suhilb","suclean","suswapbytes",
			    "suchw","sugain","suwind","bhptd","suresamp"};
	fullMiddles = new HashSet(Arrays.asList(temp));
	temp = new String[]{"surange","suximage","bhpwritecube"};
	fullEnds = new HashSet(Arrays.asList(temp));
	temp = new String[]{"bhpio","bhptranspose"};
	fullSolos = new HashSet(Arrays.asList(temp));
	temp = new String[]{"filename","pathlist", "verbose","delete"};
	HashSet hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhpio",hs);
	temp = new String[]{"filename","pathlist","verbose","pkey"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhptranspose",hs);
	temp = new String[]{"tape","buff","verbose","vblock","hfile","bfile", 
			    "over","format","conv","ns","trmin","trmax",
			    "endian","errmax"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("segyread",hs);
	temp = new String[]{"filename","pathlist","rule","maxtraces",
			    "request","endian","foldfile","display",
			    "interpolate","verbose"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhpread",hs);
	temp = new String[]{"infile","bias","scalar","null","gridtype",
			    "extrap","verbose"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhploadhdr",hs);
	temp = new String[]{"mode","Z","A","L","H","P1","P2","P3","type",
			    "limit","enforce","nval"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhpintegrate",hs);
	temp = new String[]{"ofile","bias","scalar","compress"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhpstorehdr",hs);
	temp = new String[]{"phdr","shdr","file","gridtype","extrap","bias",
			    "scalar","null"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhphorizon",hs);
	temp = new String[]{"verbose","nsegy","imach","omach"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("suclean",hs);
	temp = new String[]{"panel","tpow","epow","gpow","agc","gagc",
			    "wagc","trap","clip","qclip","qbal","pbal",
			    "mbal","scale","bias","jon","verbose","tmpdir"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("sugain",hs);
	temp = new String[]{"verbose","key","min","max","abs","j","s",
			    "count","reject","accept","dt","tmin","tmax",
			    "itmin","itmax","nt"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("suwind",hs);
	temp = new String[]{"verbose","tmpdir","n2","d1","d2","f1","f2"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("suximage",hs);
	temp = new String[]{"filename","init","pathlist","bin","action","rule",
			    "binhdr","maxens","prealloc","stripe","size",
			    "endian","stdin_endian","transpose","vkey","units",
			    "verbose"};
	hs = new HashSet(Arrays.asList(temp));
	commandOpts.put("bhpwritecube",hs);
	temp = new String[]{"format","ns"};
	hs = new HashSet(Arrays.asList(temp));
    }
    
    /**
      * Retrieve all warnings for this script generation.  (Reset on 
      * every call to either public generateScript(...).
      * @return a Vector of Strings containing warning messages
      */
    public Vector getWarnings() {
	return warnings;
    }

    /**
      * Resets important state variables every time either public 
      * generateScript(...) is called.
      */
    private void reset() {
	baseParamText = "";
	masterScriptLinear = "";
	masterScriptBsub = "";
	numScripts = 1;
	lastType = SOLO_TYPE;
	pipeChar = '|';
	parameterTable = new Hashtable();
	scatters = new Stack();
	warnings = new Vector();
	scriptBase = null;
	fileBase = null;
	xmlDir = null;
	depth = 0;
    }

    /**
      * Generates a script (and possible sub-scripts) given a base name 
      * from which to derive more script names, and a root Node to start 
      * at (the node should be called "script").  This is the default
      * behavior option and generates a non-linear script if possible.
      * @param dirName the directory path up to the directory the script 
      * exists in
      * @param scriptName the base name for the script to create, like 
      * "script"
      * @param root the root Node (under the Document node) for the XML 
      * tree to be turned into a script
      */
    public void generateScript(String dirName, String scriptName,
			       String xmlDir, Node root) 
	throws MalformedScriptException {
	generateScript(dirName,scriptName,xmlDir, root,false);
    }

    /**
      * Generates a script (and possible sub-scripts) given a base name 
      * from which to derive more script names, and a root Node to start 
      * at (the node should be called "script").  The script will be linear
      * if linear is true, and non-linear if possible otherwise
      * @param dirName the directory path up to the directory the script 
      * exists in
      * @param scriptName the base name for the script to create, like 
      * "script"
      * @param root the root Node (under the Document node) for the XML 
      * tree to be turned into a script
      * @param linear whether or not the script generation should force a 
      * linear script
      */
    private void generateScript(String dirName, String scriptName, 
				String xmlDir, Node root, boolean linear) 
	throws MalformedScriptException {
	reset();
	this.xmlDir = xmlDir;
	scriptBase = scriptName;
	fileBase = dirName + scriptBase;
	String baseName = "1";
	NodeList commands = root.getChildNodes();
	int start = fillParameters(commands,xmlDir);
	String scriptTextStart = "#!/bin/sh\n\n";
	scriptTextStart += baseParamText + "\n";
	masterScriptLinear = scriptTextStart;
	masterScriptLinear += "sh " + fileBase + baseName 
	    + "a.1.sh\n";
	masterScriptBsub = scriptTextStart;
	masterScriptBsub += "bsub -J \"" + scriptBase + baseName + "a"; 
	masterScriptBsub +="\" $* " + fileBase + baseName 
	    + "a.1.sh\n";
	generateScript(baseName,new String[]{baseParamText},commands,start,
	'a',linear);
	writeScriptFile(fileBase + ".sh",masterScriptLinear);
	writeScriptFile(fileBase + ".bsub.sh",masterScriptBsub);
	//chmod(fileBase + "*","775");
	DEditDebug.println(this,6,"chmod call");
	chmod(fileBase+"*.sh","775");
    }

    /**
      * The major script generation function called by the other 
      * generateScript() methods.  It generates all scripts except for master
      * scripts through the magic of for loops and recursion.  Master scripts 
      * are scripts that call other scripts, with "bsub jobarray" commands 
      * and the like.   Where to generate those???
      * @param baseName the base name for the script to create, like 
      * "script1" or "script1b4d3".  Multiplying the numbers should give 
      * you a number == numScripts
      * @param paramText the parameter text to put at the top of every script
      * @param commands the list of Nodes representing commands that will 
      * form this and all scripts  This particular script will contain 
      * commands in the list indexed from startIndex to when it a) hits a 
      * scatter, but it will then command the scatter to be built and resume
      * with a different name when it completes,  b) hits a gather (down a
      * level after a previous scatter), or c) reaches the end of commands.
      * @param startIndex the index into commands from which to start
      * generating the script
      * @param lastChar the character to append to the end of baseName for 
      * scriptNames.  It gets ++ treatment on scatters.
      * @param linear if true, ignore scatter/gathers (and therefore generate
      * just one script).  This means the value supplied in the top-level 
      * parameters will be used for any scatter parameters
      */
    private char generateScript(String baseName, String[] paramText, 
				NodeList commands, int startIndex,
				char lastChar, boolean linear) 
	throws MalformedScriptException {
	// go item by item and create strings for commands, then
	// append them to the return string
	boolean fromGather = false;
	lastType = SOLO_TYPE;	
	String[] scriptText = new String[numScripts];
	for (int i = 0; i < scriptText.length; i++) {
	    scriptText[i] = "#!/bin/sh\n\n" + paramText[i] + "\n";
	    DEditDebug.println(this,3,i + ":" + paramText[i]);
	}
	for (int i = startIndex; i < commands.getLength(); i++) {
	    Node command = commands.item(i);
	    if (command.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    String commandName = command.getNodeName();
	    NodeList children = command.getChildNodes();
	    String[] append = null;
	    if (shallowStarts.contains(commandName)) {
		append = buildShallowCommand(command,commandName,START_TYPE);
	    } else if (shallowMiddles.contains(commandName)) {
		append = buildShallowCommand(command,commandName,MIDDLE_TYPE);
	    } else if (shallowEnds.contains(commandName)) {
		append = buildShallowCommand(command,commandName,END_TYPE);
	    } else if (shallowSolos.contains(commandName)) {
		append = buildShallowCommand(command,commandName,SOLO_TYPE);
	    } else if (fullStarts.contains(commandName)) {
		append = buildCommand(children,commandName,START_TYPE);
	    } else if (fullMiddles.contains(commandName)) {
		append = buildCommand(children,commandName,MIDDLE_TYPE);
	    } else if (fullEnds.contains(commandName)) {
		append = buildCommand(children,commandName,END_TYPE);
	    } else if (fullSolos.contains(commandName)) {
		append = buildCommand(children,commandName,SOLO_TYPE);
	    } else if (commandName.equals("delivery")) {
		append = buildDeliveryCommand(children);
	    } else if (commandName.equals("file_P")) {
		append = buildPipeFromFile(command);
	    } else if (commandName.equals("scatter")) {
		if (linear == false) {
		    for (int j = 0; j < scriptText.length; j++) {
			String writeFileName = fileBase + baseName 
			    + lastChar + "." + (j+1);
			DEditDebug.println(this,4,"writing " 
					   + writeFileName + ".sh");
			writeScriptFile(writeFileName + ".sh",
					scriptText[j]);
		    }
		    char prev = lastChar;
		    if (lastChar == 'y') {
			lastChar = 'A';
		    } else if (lastChar == 'Y') {
			// over 50 sections...make new indexing scheme if 
			// more are required !!!
			String err = "Too many gather/scatters at one level ("
			    + " change code in ScripGenerator to allow more).";
			throw new MalformedScriptException(err);
		    } else if (!Character.isLetter(lastChar)) {
			throw new MalformedScriptException("Indexing error");
		    } else {
		        lastChar = (char)((int)lastChar + 1);
		    }
		    Scatter scatter = generateScatter(baseName,children,
						      paramText,commands,i+1,
						      prev,lastChar);
		    warnings.addAll(scatter.getWarnings());
		    prev = lastChar;
		    lastChar = (char)((int)lastChar + 1);
		    
		    for (int j = 0; j < numScripts; j++) {
			masterScriptLinear += "sh " + fileBase 
			    + baseName + lastChar + "." + (j+1)
			    + ".sh\n";
		    }
		    String bsubString = "bsub -w \'done(" + scriptBase 
			+ scatter.dependency + ")\'";
		    bsubString += " -J \"" + scriptBase + baseName + lastChar;
		    if (numScripts > 1) {
			bsubString += "[1-" + numScripts + "]";
		    }
		    bsubString += "\" $* " + fileBase + baseName + lastChar;
		    if (numScripts > 1) {
			bsubString += ".\\$LSB_JOBINDEX.sh\n";
		    } else {
			bsubString += ".1.sh\n";
		    }
		    masterScriptBsub += bsubString;
		    int numCommands = scatter.commandsConsumed;
		    i += numCommands;
		    
		    // continue
		    for (int j = 0; j < numScripts; j++) {
			scriptText[j] = "#!/bin/sh\n\n" + paramText[j] + "\n";
		    }
		}
		append = null;
	    } else if (commandName.equals("deliveryAnalyser")) {
		append = buildDeliveryAnalyserCommand(children);
	    } else if (commandName.equals("P_file")) {
		append = buildPipeToFile(command);
	    } else if (commandName.equals("gather")) {
		if (linear == false) {
		    if (lastType != END_TYPE && lastType != SOLO_TYPE) {
			String err 
			    = "gather cannot appear in the middle of a pipe.";
			throw new MalformedScriptException(err);
		    }
		    if (scatters.empty()) {
			throw new MalformedScriptException("gather but no "
							   + "scatter.");
		    }
		    Scatter scatter = (Scatter)scatters.peek();
		    scatter.commandsConsumed = i + 1 - startIndex;
		    fromGather = true;
		    break; // ???
		}
	    } else {
		if (command.getNodeType() != Node.TEXT_NODE) {
		    throw new MalformedScriptException("Unrecognized command "
						       + commandName);
		}
	    }
	    
	    if (append != null) {
		NamedNodeMap nnm = command.getAttributes();
		if (nnm != null && nnm.getLength() > 0) {
		    Node background = nnm.getNamedItem("run_in_background");
		    if (background != null) {
			String val = background.getNodeValue();
			if (val.equals("true")) {
			    for (int j = 0; j < append.length; j++) {
				append[j] += " &";
			    }
			}
		    }
		}
		
		for (int j = 0; j < append.length; j++) {
		    switch(lastType) {
		    case END_TYPE:append[j] = pipeChar + " " + append[j]
				      + "\n"; break;
		    case MIDDLE_TYPE:append[j] = pipeChar + " " + append[j] 
					 + " \\\n"; break;
		    case SOLO_TYPE:append[j] = "\n" + append[j] + "\n"; break;
		    case START_TYPE:append[j] = "\n" + append[j] + " \\\n"; 
			break;
		    }
		    scriptText[j] += append[j];
		}
		pipeChar = '|';
	    }
	}

	if (!fromGather && !scatters.empty()) {
	    throw new MalformedScriptException(scatters.size() 
					       + " closing gather(s) needed.");
	}

        // write finished script and chmod it
	FileWriter fw;
	for (int i = 0; i < numScripts; i++) {
	    String filename = fileBase + baseName + lastChar + "." 
		+ (i+1) + ".sh";
	    writeScriptFile(filename,scriptText[i]);
	}
	return lastChar;
    }
    
    /**
      * Fills in the parameter hashtable and returns an int which represents 
      * how many nodes into NodeList commands actual commands begin.  This is 
      * called recursively for includes.
      * @param commands the list of include, parameters, and command nodes 
      * of which the first two are of interest here
      * @param dir the String for the path of the directory in which the 
      * xml file being looked at resigns (this can change on recursive calls 
      * caused by included files with includes)
      * @return an int which represents how many nodes into NodeList commands 
      * actual commands begin
      */
    private int fillParameters(NodeList commands, String dir) 
	throws MalformedScriptException {
	int start;
	Node paramNode = null;
	for (start = 0; start < commands.getLength(); start++) {
	    Node temp = commands.item(start);
	    if (temp != null && temp.getNodeType() == Node.ELEMENT_NODE) { 
		String name = temp.getNodeName();
		if (name.equals("parameters")) {
		    paramNode = temp;
		    start++;
		    break;
		} else if (name.equals("include")) {
		    String[] val = getOptionValue(temp,"include","filename");
		    String filename = val[0];
		    DEditDebug.println(this,6,"include " + filename);
		    if (!filename.startsWith("/")) {
			filename = dir + filename;
		    }
		    Document d = editor.createDOMforXML(filename);
		    if (d != null) {
			Node n = d.getDocumentElement();
			NodeList nl = n.getChildNodes();
			int lastDir = Math.max(0,filename.lastIndexOf("/"));
			dir = filename.substring(0,lastDir+1);
			fillParameters(nl,dir);
			DEditDebug.println(this,6,"included " 
					   + filename);
		    }
		}
	    }
	}
	
	if (paramNode == null 
	    || !paramNode.getNodeName().equals("parameters")) {
	    throw new MalformedScriptException("missing parameters group");
	}
	NodeList parameters = paramNode.getChildNodes();
	for (int i = 0; i < parameters.getLength(); i++) {
	    Node parameter = parameters.item(i);
	    if (parameter.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    NodeList paramChilds = parameter.getChildNodes();
	    Node nameChild = null;
	    int j = 0;
	    for (; j < paramChilds.getLength(); j++) {
		Node temp = paramChilds.item(j);
		if (temp.getNodeType() == Node.ELEMENT_NODE) {
		    nameChild = temp;
		    DEditDebug.println(this,3,"found on " + j + " of " 
				       + paramChilds.getLength());
		    j++;
		    break;
		}
	    }
	    if (nameChild == null || !nameChild.getNodeName().equals("name")) {
		throw new MalformedScriptException("parameter missing name");
	    }
	    Node nameTextChild = nameChild.getFirstChild();
	    if (nameTextChild == null 
		|| nameTextChild.getNodeType() != Node.TEXT_NODE) {
		throw new MalformedScriptException("name has no value");
	    }
	    String name = nameTextChild.getNodeValue();
	    Node valueChild = null;
	    for (; j < paramChilds.getLength(); j++) {
		Node temp = paramChilds.item(j);
		if (temp.getNodeType() == Node.ELEMENT_NODE) {
		    valueChild = temp;
		    DEditDebug.println(this,3,"found on " + j + " of " 
				       + paramChilds.getLength());
		     break;
		}
	    }
	    if (valueChild == null 
		|| !valueChild.getNodeName().equals("value")) {
		throw new MalformedScriptException("parameter "+ name 
						   + " missing value");
	    }
	    
	    String value = getOptionValue(valueChild,"parameter","value")[0];
	    baseParamText += name + "=" + value + "\n";
	    parameterTable.put(name,value);
	}
	return start;
    }

    /**
      * Checks to see that it is legal to go from the type of the previous 
      * command to the next command.  Throws a MalformedScriptException if it 
      * isn't.
      * @param oldType the type of the previous command 
      * @param newType the type of the next command 
      */
    private static void typeCheck(int oldType, int newType) 
	throws MalformedScriptException {
	switch (oldType) {
	case SOLO_TYPE:
	case END_TYPE:if (newType == START_TYPE 
			  || newType == SOLO_TYPE) return;break;
	case MIDDLE_TYPE:
	case START_TYPE: if (newType == END_TYPE 
			     || newType == MIDDLE_TYPE) return;break;
	default: throw new MalformedScriptException("Previous type is invalid "
						    + oldType);
	}
	throw new MalformedScriptException("Pipe flow error: " + oldType 
					   + " " + newType);
    }

    /**
      * Builds one of the "shallow" commands.  That is, commands which are 
      * consist solely of a leaf node and have no children.
      * @param command the Node representing the command
      * @param name the name of the command, .equals(command.getNodeName())
      * @param type the type of the command to make (START_TYPE, MIDDLE_TYPE,
      * END_TYPE, and SOLO_TYPE)
      * @return an array of "shallow" commands as Strings
      */
    private String[] buildShallowCommand(Node command, String name, int type) 
	throws MalformedScriptException {
	typeCheck(lastType,type);
	String[] ret = getOptionValue(command,name,"text"); 
	// no further checking ???
	if (name.equals("cp") || name.equals("rm")) {
	    prepend(name + " ",ret);
	}
	lastType = type;
	return ret;
    }

    /**
      * Builds a command given the Nodes it will have as command line options,
      * the name of the command, and the type of command (in terms of placing 
      * in the flow of a pipe).  It defers to handleSpecialCase if the name of 
      * a Node in children is not in the HashSet corresponding to 
      * commandOpts.get(name).
      * @param children the Nodes which will build the command line options 
      * for the command named "name".
      * @param name the name of the command to build
      * @param type the type of the command in terms of flow
      * @return an array of commands named "name" as Strings
      */
    private String[] buildCommand(NodeList children, String name, int type) 
	throws MalformedScriptException {
	typeCheck(lastType,type);
	String[] ret = expand(name,numScripts);
	HashSet basic = (HashSet)commandOpts.get(name);
	for (int i = 0; i < children.getLength(); i++) {
	    Node cmdlnOption = children.item(i);
	    if (cmdlnOption.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    String optName = cmdlnOption.getNodeName();
	    if (basic != null && basic.contains(optName)) {
		append(ret," "+ optName + "=");
		append(ret,getOptionValue(cmdlnOption,name,optName));
	    } else {
		append(ret,handleSpecialCase(cmdlnOption,name,optName));
	    }
	}
	lastType = type;
	return ret;
    }

    /**
      * Handles special cases not covered by "basic" behavior ("'top-level 
      * node name'='that node's value'" type options).
      * @param cmdlnOption the Node corresponding to the complex command line 
      * option
      * @param name the name of the command that has this option
      * @param optName the name of the "complex" option
      * @return an array of the text for a command line option for command 
      * "name" as Strings
      */
    private String[] handleSpecialCase(Node cmdlnOption, String name, 
				       String optName) 
	throws MalformedScriptException {
	if (name.equals("bhpread")) {
	    return buildBHPReadCommandComponent(cmdlnOption,optName);
	} else if (name.equals("sushw")) {
	    return buildSUShwCommandComponent(cmdlnOption,optName);
	} else if (name.equals("suchw")) {
	    return buildSUChwCommandComponent(cmdlnOption,optName);
	} else if (name.equals("bhploadhdr")) {
	    return buildBHPLoadHdrCommandComponent(cmdlnOption,optName);
	} else if (name.equals("bhpstorehdr")) {
	    return buildBHPStoreHdrCommandComponent(cmdlnOption,optName);
	} else if (name.equals("bhpintegrate")) {
	    return buildBHPIntegrateCommandComponent(cmdlnOption,optName);
	} else if (name.equals("bhpwritecube")) {
	    return buildBHPWriteCubeCommandComponent(cmdlnOption,optName);
	} else {
	    String err = "Illegal command line option " + optName 
		+ " in " + name + " command.";
	    throw new MalformedScriptException(err);
	}
    }

    /**
      * Constructs a mid-pipe InversionAnalyser.sh given the provided children.
      * @param children the Nodes from which to construct the 
      * DeliveryAnalyser.sh command.
      * @return an array of Strings which are InverstionAnalyser commands.
      */
    private String[] buildDeliveryAnalyserCommand(NodeList children)
	throws MalformedScriptException {
	if (lastType != START_TYPE && lastType != MIDDLE_TYPE) {
	    throw new MalformedScriptException("Must run deliveryAnalyser" 
					       + "in the middle of a pipe.");
	}
	String[] ret = expand("deliveryAnalyser",numScripts);
	String[] listArr = new String[]{"BHPcommand","Raftery","filter",
					"histogram","layer","mean","quantile",
					"stddev","p","p_ascii","p_ascii2","s"};
	Vector listArgs = new Vector(Arrays.asList(listArr));
	for (int i = 0; i < children.getLength(); i++) {
	    Node cmdlnOption = children.item(i);
	    if (cmdlnOption.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    String optName = cmdlnOption.getNodeName();
	    if (optName.equals("f") || optName.equals("p_ascii_full")) {
		// options with one '-' and no arguments
		optName = optName.replace('_','-');
		append(ret," -" + optName);
	    } else if (optName.equals("full") || optName.equals("full_means")
		       || optName.equals("full_stddevs")
		       || optName.equals("wait")) {
		// options with two '-' and no arguments
		optName = optName.replace('_','-');
		append(ret," --" + optName);
	    } else if (optName.equals("i") || optName.equals("v")) {
		// options with one '-' and one argument
		append(ret," -" + optName + " ");
		append(ret,getOptionValue(cmdlnOption,"deliveryAnalyser",
					  optName));
	    } else if (optName.equals("realisation_header_word")
		       || optName.equals("full_quantile")) {
		// options with two '-' and one argument
		optName = optName.replace('_','-');
		append(ret," --" + optName + " ");
		append(ret,getOptionValue(cmdlnOption,"deliveryAnalyser",
					  optName));
	    } else if (listArgs.contains(optName)) {
		// complicated arguments
		NodeList kids = cmdlnOption.getChildNodes();
		String optNames[] = {" --BHPcommand "," --Raftery ",
				     " --filter "," --histogram "," --layer ",
				     " --mean "," --quantile "," --stddev ",
				     " -p "," -p-ascii "," -p-ascii2 ",
				     " -s "};
		String[] propArr = {"property","N"};
		String[][] names = {
		    {"filename"}, // BHPCommand
		    propArr, // Raftery
		    {"property","comparator","comparison_value",
		     "layer_no"}, //filter
		    propArr, // histogram
		    {"t_start_in_seconds","t_end_in_seconds"}, // layer
		    propArr, // mean
		    {"property","N","Q"}, // quantile
		    propArr, // stddev
		    propArr, // p
		    propArr, // p-ascii
		    {"property1","property2","N"}, // p-ascii2
		    {"t_start_in_seconds","t_end_in_seconds",
		     "reflection_coefficient","wavelet_filename"} // s
		};
		int index = listArgs.indexOf(optName);
		append(ret,optNames[index]);
		if (kids == null || kids.getLength() < names[index].length) {
		    String err = "Not enough children for " + optNames[index];
		    throw new MalformedScriptException(err + ".");
		}
		int j = 0;
		for (int k = 0; k < names[index].length; k++) {
		    String name = names[index][k];
		    for (; j < kids.getLength(); j++) {
			Node temp = kids.item(j);
			if (temp.getNodeName().equals(name)) {
			    boolean escape = name.equals("comparator");
			    if (escape) {
				append(ret,"'");
			    }
			    append(ret,getOptionValue(temp,"deliveryAnalyser",
						      optNames[index] + "->" 
						      + name));
			    if (escape) {
				append(ret,"'");
			    }
			    break;
			}
		    }
		    j++; 
		    if (j > kids.getLength()) {
			throw new MalformedScriptException(optNames[index] 
							   + " missing " 
							   + name +".");
		    }
		    append(ret," ");
		}
	    } else {
		throw new MalformedScriptException("Illegal command line "
						   + "option " + optName 
						   + " in deliveryAnalyser.");
	    }
	}
	
	lastType = MIDDLE_TYPE;
	return ret;
    }

    /**
      * Constructs a beginning pipe seismicinverter.sh given the provided 
      * children.
      * @param children the Nodes from which to construct the 
      * seismicinverter.sh command.
      * @return an array of Strings which are SeimsicInverter commands.
      */
    private String[] buildDeliveryCommand(NodeList children) 
	throws MalformedScriptException {
	if (lastType != END_TYPE && lastType != SOLO_TYPE) {
	    throw new MalformedScriptException("Can't run delivery" 
					       + "in the middle of a pipe.");
	}
	String[] ret = expand("delivery",numScripts);
	String[] listArr = new String[]{"BHPcommand","make_synthetic","S"};
	Vector listArgs = new Vector(Arrays.asList(listArr));
	for (int i = 0; i < children.getLength(); i++) {
	    Node cmdlnOption = children.item(i);
	    if (cmdlnOption.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    String optName = cmdlnOption.getNodeName();
	    if (optName.equals("ID") || optName.equals("II") 
		|| optName.equals("IS") || optName.equals("PD") 
		|| optName.equals("RWS")|| optName.equals("p")) {
		// all options with one '-' and no parameters
		append(ret," -" + optName);
	    } else if (optName.equals("Gassman_diagnostics")
		       || optName.equals("Gelman_diagnostics")
		       || optName.equals("Newton")
		       || optName.equals("Npre_decimation")
		       || optName.equals("Raftery_auto_decimation")
		       || optName.equals("Raftery_diagnostics")
		       || optName.equals("Xsections")
		       || optName.equals("additive_times")
		       || optName.equals("check_prior")
		       || optName.equals("dump_optimal_seismic")
		       || optName.equals("lock_depths")
		       || optName.equals("mode_profiles")
		       || optName.equals("no_decimation")
		       || optName.equals("no_smart_mode_finding_starting_points")
		       || optName.equals("p_fluid")
		       || optName.equals("parametric_depths")
		       || optName.equals("suppress_model_dump")
		       || optName.equals("times_only")
		       || optName.equals("truncate_synthetics")) {
		// all options with two '-' and no parameters
		optName = optName.replace('_','-');
		append(ret," --" + optName);
	    } else if (optName.equals("N") || optName.equals("m")
		       || optName.equals("o") || optName.equals("v")) {
		// all options with one '-' and one parameter
		append(ret," -" + optName + " ");
		append(ret,getOptionValue(cmdlnOption,"delivery",
					  optName));
	    } else if (optName.equals("IF") 
		       || optName.equals("far_stack_filename")
		       || optName.equals("far_stack_wavelet_filename")
		       || optName.equals("near_stack_filename")
		       || optName.equals("near_stack_wavelet_filename")
		       || optName.equals("trace_filter")) {
		// all options with two '-' and one parameter
		optName = optName.replace('_','-');
		append(ret," --" + optName + " ");
		append(ret,getOptionValue(cmdlnOption,"delivery",
					  optName));
	    } else if (optName.equals("filename")) {
		// the filename parameter
		append(ret," ");
		append(ret,getOptionValue(cmdlnOption,"delivery",
					  optName));
	    } else if (listArgs.contains(optName)) {
		// special complicated options
		NodeList kids = cmdlnOption.getChildNodes();
		String optNames[] = {" --BHPcommand "," --make-synthetic ",
				     " -S "};
		String[][] names = {
		    {"filename"}, // BHPcommand 
		    {"freq_in_Hz","sample_interval_in_ms",
		     "t_start_in_seconds","t_end_in_seconds"},// make-synthetic
		    {"t_start_in_seconds","t_end_in_seconds"}
		};
		int index = listArgs.indexOf(optName);
		append(ret,optNames[index]);
		if (kids == null || kids.getLength() < names[index].length) {
		    String err = "Not enough children for " + optNames[index];
		    throw new MalformedScriptException(err + ".");
		}
		int j = 0;
		for (int k = 0; k < names[index].length; k++) {
		    String name = names[index][k];
		    for (; j < kids.getLength(); j++) {
			Node temp = kids.item(j);
			if (temp.getNodeName().equals(name)) {
			    append(ret,getOptionValue(temp,"delivery",
						      optNames[index] + "->" 
						      + name));
			    break;
			}
		    }
		    j++; 
		    if (j > kids.getLength()) {
			throw new MalformedScriptException(optNames[index] 
							   + " missing " 
							   + name +".");
		    }
		    append(ret," ");
		}
	    } else {
		throw new MalformedScriptException("Illegal command line "
						   + "option " + optName 
						   + " in delivery");
	    }
	}
	
	lastType = SOLO_TYPE;
	return ret;
	
    }

    /**
      * Builds a command via "cat" which turns a file into a data pipe to 
      * be fed as standard input to another process.
      * @param command the Node from which to generate the command
      * @return an array of "cat file" commands as Strings
      */
    private String[] buildPipeFromFile(Node command)  
	throws MalformedScriptException {
	if (lastType != END_TYPE && lastType != SOLO_TYPE) {
	    throw new MalformedScriptException("file_P must be at the "
					       + "start of a pipe.");
	}
	String[] ret = null;
	ret = getOptionValue(command,"file_P","filename");
	lastType = START_TYPE;
	prepend("cat ",ret);
	return ret;
    }

    /**
      * Builds one of the "complex" command line options for bhpread (actually,
      * if optName.equals("keys") it builds two) and returns it as an array of 
      * Strings (which may very by script index via a scatter).
      * @param cmdlnOption a Node corresponding to a complex option
      * @param optName the name of the option, 
      * .equals(cmdlnOption.getNodeName())
      * @return an array of one of the "complex" options as Strings
      */
    private String[] buildBHPReadCommandComponent(Node cmdlnOption,
						  String optName) 
	throws MalformedScriptException {
	if (optName.equals("keys")) {
	    NodeList keys = cmdlnOption.getChildNodes();
	    if (keys == null || keys.getLength() == 0) {
		throw new MalformedScriptException("Keys for bhpread specified"
						   + " but none provided");
	    }
	    Node[] elemKeys = pareList(keys);
	    int length = elemKeys.length;
	    String[][] keyNames = new String[length][numScripts];
	    String[][] keyListRanges = new String[length][numScripts];
	    for (int j = 0; j < length; j++) {
		Node keyPair = elemKeys[j];
		NodeList keyPairKids = keyPair.getChildNodes();
		Node keyName = null;
		Node keyRange = null;
		int k = 0;
		for (; k < keyPairKids.getLength(); k++) {
		    Node temp = keyPairKids.item(k);
		    if (temp != null 
			&& temp.getNodeType() == Node.ELEMENT_NODE
			&& temp.getNodeName().equals("name")) {
			keyName = temp;
			k++;
			break;
		    }
		}
		for (; k < keyPairKids.getLength(); k++) {
		    Node temp = keyPairKids.item(k);
		    if (temp != null
			&& temp.getNodeType() == Node.ELEMENT_NODE
			&& temp.getNodeName().equals("range")) {
			keyRange = temp;
		    }
		}
		if (keyName == null || keyRange == null) {
		    throw new MalformedScriptException("key pair missing "
						       + " name and "
						       + "range.");
		}
		// more checking???
		String[] names = getOptionValue(keyName,"bhpread",
						"key->name");
		String[] ranges = getOptionValue(keyRange,"bhpread",
						 "key->range");
		// !!! may be able to check range validity here if need be 
		// - but could only be invalid with list-from-file scatters, so
		// in that case it requires the user to give correct data
		keyNames[j] = names;
		keyListRanges[j] = ranges;
	    }
	    String[] ret = buildSeparatedList(",",keyNames," keys=");
	    append(ret,buildSeparatedList(":",keyListRanges," keylist="));
	    return ret;
	} else if (optName.equals("properties")
		   || optName.equals("horizons")) {
	    NodeList nodes = cmdlnOption.getChildNodes();
	    return buildSeparatedList(",",nodes," " + optName + "=",
					   "bhpread",optName);
	} else {
	    String err = "Illegal command line option " + optName 
		+ " in bhpread command";
	    throw new MalformedScriptException(err);
	}
    }

    /**
      * Builds and returns a "complex" command line option for sushw, in this 
      * case it is always "operators".
      * @param cmdlnOption the Node representing the "complex" command line 
      * option
      * @param optName the name of the "complex" command line option, 
      * .equals(cmdlnOption.getNodeName())
      * @return an array of "operators" options as Strings
      */
    private String[] buildSUShwCommandComponent(Node cmdlnOption,
						String optName) 
	throws MalformedScriptException {
	if (optName.equals("operators")) {
	    NodeList operators = cmdlnOption.getChildNodes();
	    if (operators == null || operators.getLength() == 0) {
		throw new MalformedScriptException("Operators for sushw "
						   +"specified "
						   + "but none provided.");
	    }
	    String[] ret = expand("",numScripts);
	    Node[] ops = pareList(operators);
	    int length = ops.length;
	    String[][] key = new String[length][numScripts];
	    String[][] a = new String[length][numScripts];
	    String[][] b = new String[length][numScripts];
	    String[][] c = new String[length][numScripts];
	    String[][] d = new String[length][numScripts];
	    String[][] j = new String[length][numScripts];
	    String[][][] types = new String[][][]{key,a,b,c,d,j};
	    String[] typeNames = new String[]{"key","a","b","c","d","j"};
	    for (int n = 0; n < length; n++) {
		Node operator = ops[n];
		NodeList opElems = operator.getChildNodes();
		if (opElems == null || opElems.getLength() == 0) {
		    throw new MalformedScriptException("Operator for sushw "
						       + "specified but none "
						       + "provided.");
		}
		int m = 0;
		for (int k = 0; k < opElems.getLength(); k++) {
		    Node temp = opElems.item(k);
		    if (temp.getNodeType() == Node.ELEMENT_NODE) {
			types[m][n] = getOptionValue(temp,"sushw",
						     "operator->"
						     +typeNames[m]);
			m++;
		    }
		}
	    }
	    for (int n = 0; n < types.length; n++) {
		String[][] type = types[n];
		append(ret,buildSeparatedList(",",type," "+typeNames[n]+"="));
	    }
	    return ret;
	} else {
	    throw new MalformedScriptException("sushw has a top-level node "
					       + "other than operators.");
	}
    }
   
    /**
      * Builds and returns a "complex" command line option for suchw, in this 
      * case it is always "operators".
      * @param cmdlnOption the Node representing the "complex" command line 
      * option
      * @param optName the name of the "complex" command line option, 
      * .equals(cmdlnOption.getNodeName())
      * @return an array of "operators" options as Strings
      */
    private String[] buildSUChwCommandComponent(Node cmdlnOption, 
						String optName) 
	throws MalformedScriptException {
	if (optName.equals("operators")) {
	    NodeList operators = cmdlnOption.getChildNodes();
	    if (operators == null || operators.getLength() == 0) {
		throw new MalformedScriptException("Operators for suchw "
						   +"specified "
						   + "but none provided.");
	    }
	    String[] ret = expand("",numScripts);
	    int length = operators.getLength();
	    Node[] ops = pareList(operators);
	    length = ops.length;
	    String[][] key1 = new String[length][numScripts];
	    String[][] key2 = new String[length][numScripts];
	    String[][] key3 = new String[length][numScripts];
	    String[][] a = new String[length][numScripts];
	    String[][] b = new String[length][numScripts];
	    String[][] c = new String[length][numScripts];
	    String[][] d = new String[length][numScripts];
	    String[][][] types = new String[][][]{key1,key2,key3,a,b,c,d};
	    String[] typeNames = new String[]{"key1","key2","key3","a","b","c",
					      "d"};
	    for (int n = 0; n < length; n++) {
		Node operator = ops[n];
		NodeList opElems = operator.getChildNodes();
		if (opElems == null || opElems.getLength() == 0) {
		    throw new MalformedScriptException("Operator for suchw "
						       + "specified but none "
						       + "provided.");
		}
		int m = 0;
		for (int k = 0; k < opElems.getLength(); k++) {
		    Node temp = opElems.item(k);
		    if (temp.getNodeType() == Node.ELEMENT_NODE) {
			types[m][n] = getOptionValue(temp,"suscw",
						     "operator->"
						     +typeNames[m]);
			m++;
		    }
		}
	    }
	    for (int n = 0; n < types.length; n++) {
		String[][] type = types[n];
		append(ret,buildSeparatedList(",",type," "+typeNames[n]+"="));
	    }
	    return ret;
	} else {
	    throw new MalformedScriptException("suchw has a top-level node "
					       +"other than operators.");
	}
    }

    /**
      * Builds and returns an array of the "keys" command line option for
      * bhploadhdr.
      * @param cmdlnOption the Node for "keys"
      * @param optName .equals("keys")
      * @return an array of "keys" command line options as Strings
      */
    private String[] buildBHPLoadHdrCommandComponent(Node cmdlnOption,
						     String optName) 
	throws MalformedScriptException {
	if (optName.equals("keys")) {
	    NodeList keysList = cmdlnOption.getChildNodes();
	    Node[] keys = pareList(keysList);
	    if (keys == null || keys.length < 3 
		|| keys.length > 4) {
		throw new MalformedScriptException("Command bhploadhdr "
						   + "needs 3 or 4 keys.");
	    }
	    return buildSeparatedList(",",keys," key=","bhploadhdr","key");
	} else {
	    String err = "Illegal command " + optName + " specified"
		+ " for command bhploadhdr.";
	    throw new MalformedScriptException(err);
	}
    }
    
    /**
      * Builds and returns the "search" option for a bhpintegrate command.
      * @param cmdlnOption the Node for the "search" option
      * @param optName .equals("search")
      * @return an array of "search" options as Strings
      */
    private String[] buildBHPIntegrateCommandComponent(Node cmdlnOption,
						       String optName) 
	throws MalformedScriptException {
	if (optName.equals("search")) {
	    NodeList searchValsNodes = cmdlnOption.getChildNodes();
	    Node[] searchVals = pareList(searchValsNodes);
	    if (searchVals == null || searchVals.length <= 0 
		|| searchVals.length > 4) {
		String err = "Search values improperly specified for "
		    + "command bhpintegrate";
		throw new MalformedScriptException(err);
	    }
	    return buildSeparatedList(",",searchVals," search=",
					  "bhpintegrate",optName);
	} else {
	    throw new MalformedScriptException("Illegal option " + optName
					       + " specified for command "
					       + "bhpintegrate.");
	}
    }

    /**
      * Builds and returns an array of the "keys" command line option for
      * bhpstorehdr.
      * @param cmdlnOption the Node for "keys"
      * @param optName .equals("keys")
      * @return an array of "keys" command line options as Strings
      */
    private String[] buildBHPStoreHdrCommandComponent(Node cmdlnOption,
						      String optName) 
	throws MalformedScriptException {
	if (optName.equals("keys")) {
	    NodeList keysNodes = cmdlnOption.getChildNodes();
	    Node[] keys = pareList(keysNodes);
	    if (keys == null || keys.length <= 0 
		|| keys.length > 3) {
		throw new MalformedScriptException("Command bhpstorehdr "
						   +"needs up to 3 keys.");
	    }
	    return buildSeparatedList(",",keys," key=","bhpstorehdr","keys");
	} else {
	    throw new MalformedScriptException("Illegal command " + optName
					       + " specified for command "
					       + "bhpstorehdr.");
	}
    }
    
    /**
      * Builds and returns a "complex" command line option for bhpwritecube.
      * @param cmdlnOption a Node corresponding to a "complex" command line 
      * option
      * @param optName the name of the "complex" command line option,
      * .equals(cmdlnOption.getNodeName())
      * @return an array of "complex" command line options as Strings
      */
    private String[] buildBHPWriteCubeCommandComponent(Node cmdlnOption, 
						       String optName) 
	throws MalformedScriptException {
	if (optName.equals("keys")) {
	    NodeList keysNodes = cmdlnOption.getChildNodes();
	    if (keysNodes == null || keysNodes.getLength() == 0) {
		throw new MalformedScriptException("Keys for bhpwritecube "
						   + "specified but none "
						   + "provided.");
	    }
	    Node[] keys = pareList(keysNodes);
	    int length = keys.length;
	    if (length > 5) {
		throw new MalformedScriptException("bhpwritecube accepts "
						   + "only up to 5 keys.");
	    }
	    String[] ret = expand("",numScripts);
	    for (int n = 0; n < length; n++) {
		Node key = keys[n];
		NodeList keyElemsNodes = key.getChildNodes();
		String name = "key" + (n+1);
		append(ret,buildSeparatedList(",",keyElemsNodes,
						   " " + name +  "=",
						   "bhpwritecube",name));
	    }
	    return ret;
	} else if (optName.equals("properties")
		   || optName.equals("horizons")) {
	    NodeList nodes = cmdlnOption.getChildNodes();
	    if (nodes == null || nodes.getLength() == 0) {
		throw new MalformedScriptException(optName + " specified "
						   + "in bhpread but none "
						   + "provided");
	    }
	    return buildSeparatedList(",",nodes," " + optName + "=",
					   "bhwritecube",optName);
	} else {
	    throw new MalformedScriptException("Illegal option " + optName
					       + " specified for command "
					       + "bhpwritecube.");
	}
    }

    /**
      * Builds a pipe to a file which expects input from standard in and 
      * terminates a pipe flow.
      * @param command the Node to derive the file name from
      * @return an array of the name of the file (might be different from each 
      * other if in a scatter)
      */
    private String[] buildPipeToFile(Node command) 
	throws MalformedScriptException {
	if (lastType != MIDDLE_TYPE && lastType != START_TYPE) {
	    throw new MalformedScriptException("P_file must be at the "
					       + "end of a pipe.");
	}
	String[] ret = null;
	ret = getOptionValue(command,"P_file","filename");
	pipeChar = '>';
	lastType = END_TYPE;
	return ret;
    }

    /**
      * Generates a scatter of a script based on a Scatter object derived from 
      * the options specified by the Nodes in children, then continues on 
      * with commands from startIndex+1, terminating at the gather nodes.
      * Returns the generated scatter which tells the script that called 
      * generateScatter where to "pick up."
      * @param baseName the baseName to use for the scatter script files
      * @param children the Nodes to generate a new Scatter object from
      * @param paramText the parameter text to be placed at the top of every 
      * script file
      * @param commands the list of Nodes which specify commands
      * @param startIndex the current index into commands where the "scatter" 
      * node was hit.
      * @param prev the previous character that was the last character of 
      * the name of a script, == lastChar-- unless lastChar == 'A', then 
      * it's 'y'.  Used for job array dependency
      * purposes.
      * @param lastChar the current last character to append to the name of 
      * the job and script files
      * @return the Scatter object generated and completed
      */
    private Scatter generateScatter(String baseName, NodeList children,
				    String[] paramText,NodeList commands,
				    int startIndex, char prev, char lastChar) 
	throws MalformedScriptException {
	if (lastType != END_TYPE && lastType != SOLO_TYPE) {
	    String err 
		= "scatter cannot appear in the middle of a pipe.";
	    throw new MalformedScriptException(err);
	}
	depth++;
	int maxThreads = -1;
	int currThreads = 1;
	Vector params = new Vector();
	for (int i = 0; i < children.getLength(); i++) {
	    Node scatterOption = children.item(i);
	    if (scatterOption.getNodeType() != Node.ELEMENT_NODE) {
		continue;
	    }
	    String optName = scatterOption.getNodeName();
	    // hashtable???
	    if (optName.equals("max_number_threads")) {
		String[] threadNum = getOptionValue(scatterOption,"scatter",
						  optName);
		maxThreads = Integer.parseInt(threadNum[0]);
	    } else if (optName.equals("scatter_parameter")) {
		NodeList paramOpts = scatterOption.getChildNodes();
		String name = "name";
		int min = 0, max = 0, increment = 0;
		int min_number_per_thread = -1;
		String listfile = null;
		for (int j = 0; j < paramOpts.getLength(); j++) {
		    Node paramOption = paramOpts.item(j);
		    if (paramOption.getNodeType() != Node.ELEMENT_NODE) {
			continue;
		    }
		    String paramName = paramOption.getNodeName();
		    if (paramName.equals("name")) {
			String[] names = getOptionValue(paramOption,
							"scatter->parameter",
							paramName);
			name = names[0];
		    } else if (paramName.equals("range")) {
			NodeList rangeOpts = paramOption.getChildNodes();
			for (int k = 0; k < rangeOpts.getLength(); k++) {
			    Node rangeOption = rangeOpts.item(k);
			    if (rangeOption.getNodeType() 
				!= Node.ELEMENT_NODE) {
				continue;
			    }
			    paramName = rangeOption.getNodeName();
			    if (paramName.equals("min")) {
				String[] minStr 
				    = getOptionValue(rangeOption,
						     "scatter->parameter",
						     paramName);
				min = Integer.parseInt(minStr[0]);
			    } else if (paramName.equals("max")) {
				String[] maxStr 
				    = getOptionValue(rangeOption,
						     "scatter->parameter",
						     paramName);
				max = Integer.parseInt(maxStr[0]);
			    } else if (paramName.equals("increment")) {
				String[] incStr
				    = getOptionValue(rangeOption,
						     "scatter->parameter",
						     paramName);
				increment = Integer.parseInt(incStr[0]);
			    } else if (paramName.equals("min_number"
							+"_per_thread")) {
				String[] mnptStr 
				    = getOptionValue(rangeOption,
						     "scatter->parameter",
						     paramName);
				min_number_per_thread 
				    = Integer.parseInt(mnptStr[0]);
			    }
			}
		    } else if (paramName.equals("listfile")) {
			String[] files = getOptionValue(paramOption,
							"scatter->parameter",
							paramName);
			listfile = files[0];
		    }
		}
		AScatterParam sp;
		if (listfile == null) {
		    sp  = new RangeScatterParam(name,min,max,increment,
						min_number_per_thread);
		} else {
		    if (!listfile.startsWith("/")) {
			listfile = xmlDir + listfile;
		    }
		    sp = new ListScatterParam(name,listfile,editor);
		}
		params.add(sp);
	    }
	}
	AScatterParam[] scatParams = new AScatterParam[params.size()];
	params.toArray(scatParams);
	if (!scatters.empty()) {
	    Scatter lastScatter = (Scatter)scatters.peek();
	    if (lastScatter.maxThreads == -1 || maxThreads == -1) {
		maxThreads = -1;
	    } else {
		maxThreads *= lastScatter.maxThreads;
	    }
	    currThreads = lastScatter.currThreads;
	}
	parameterTable.put("index" + depth,"meaningless_value"+depth);
	Scatter scatter = new Scatter(maxThreads,currThreads,scatParams);
	scatters.push(scatter);
	int oldNumScripts = numScripts;
	numScripts = scatter.currThreads;
	String scriptName = baseName + lastChar + scatter.mulFact;
	String[] newParamText = new String[numScripts];	
	for (int i = 0; i < oldNumScripts; i++) {
	    String newBase = scriptName + "a.";
	    for (int j = 0; j < scatter.mulFact; j++) {
		int index = i*scatter.mulFact+j;
		masterScriptLinear += "sh " + fileBase + newBase  
		    + (index+1) + ".sh\n";
		newParamText[index] = paramText[i] + "index" + depth + "=" 
		    + (j+1) + "\n";
	    }
	}
	String range = "[1-" + numScripts + "]";
	String bsubString = "bsub -w \'done(" + scriptBase + baseName 
	    + prev + ")\'";
	bsubString += " -J \"" + scriptBase + scriptName + 'a' + range;
	bsubString += "\" $* " + fileBase + scriptName + 'a'
	    + ".\\$LSB_JOBINDEX.sh\n";
	masterScriptBsub += bsubString;
	char last = generateScript(scriptName,newParamText,commands,
				   startIndex,'a',false);
	scatter.dependency = scriptName + last;
	scatters.pop();
	numScripts = oldNumScripts;
	// resume after scatter
	parameterTable.remove("index" + depth);
	depth--;
	lastType = SOLO_TYPE;
	return scatter;	
    }
    
    // The following are helper static functions which assist when dealing 
    // with combining Strings and arrays of Strings
    private static void append(String[] partA, String[] partB) 
	throws MalformedScriptException {
	if (partA.length != partB.length) {
	    if (partB.length == 1) {
		append(partA,partB[0]);
		//DEditDebug.println("",6,"append " + partB[0]);
		return;
 	    } else {
		throw new MalformedScriptException("String append arrays of "
						   + "unequal size.");
	    }
	}
	for (int i = 0; i < partA.length; i++) {
	    partA[i] += partB[i];
	}
    }

    private static void append(String[] partA, String partB) {
	for (int i = 0; i < partA.length; i++) {
	    partA[i] += partB;
	}
    }

    private static void prepend(String partA, String[] partB) {
	for (int i = 0; i < partB.length; i++) {
	    partB[i] = partA + partB[i];
	}
    }

    private static String[] expand(String str, int size) {
	String[] ret = new String[size];
	for (int i = 0; i < size; i++) {
	    ret[i] = str;
	}
	return ret;
    }

    /**
      * Constructs and returns an array of separated lists following the 
      * pattern "start"+value+"separator"+value+...+value+"separator"+value.
      * @param separator the String used to separate values in the list
      * @param nl the NodeList containing the Node's with the values
      * @param start the String to use at the start of the list (be sure to 
      * include '=' character in this string if it is required, like key=...)
      * @param commandName the name of the command for which this 
      * list is being built
      * @param optionName the name of the command line option within 
      * "commandName" for which this list is being built
      * @return an array of separated lists
      */
    private String[] buildSeparatedList(String separator, NodeList nl,
					String start, String commandName,
					String optionName) 
	throws MalformedScriptException {
	if (nl == null || nl.getLength() == 0) {
	    throw new MalformedScriptException(optionName + " specified for "
					       + commandName + " but none " 
					       + "provided.");
	}
	Node[] nodes = pareList(nl);
	return buildSeparatedList(separator,nodes,start,commandName,
				  optionName);
    }
    
    /**
      * Constructs and returns an array of separated lists following the 
      * pattern "start"+value+"separator"+value+...+value+"separator"+value.
      * @param separator the String used to separate values in the list
      * @param nodes Node's with the values for the list
      * @param start the String to use at the start of the list (be sure to 
      * include '=' character in this string if it is required, like key=...)
      * @param commandName the name of the command for which this 
      * list is being built
      * @param optionName the name of the command line option within 
      * "commandName" for which this list is being built
      * @return an array of separated lists
      */
    private String[] buildSeparatedList(String separator, Node[] nodes,
					String start, String commandName,
					String optionName) 
	throws MalformedScriptException {
	if (nodes == null || nodes.length == 0) {
	    throw new MalformedScriptException(optionName + " specified for "
					       + commandName + " but none " 
					       + "provided.");
	}
	if (start == null) {
	    start = "";
	}
	String[][] vals = new String[nodes.length][numScripts];
	for (int i = 0; i < nodes.length; i++) {
	    vals[i] = getOptionValue(nodes[i],commandName,optionName);
	}
	return buildSeparatedList(separator,vals,start);
    }
    
    /**
      * Constructs and returns an array of separated lists following the 
      * pattern "start"+value+"separator"+value+...+value+"separator"+value.
      * @param separator the String used to separate values in the list
      * @param values the values to put in the lists
      * @param start the String to use at the start of the list (be sure to 
      * @return an array of separated lists
      */
    private String[] buildSeparatedList(String separator, String[][] values, 
					String start) 
	throws MalformedScriptException {
	String[] ret = expand(start,numScripts);
	append(ret,values[0]);
	for (int j = 1; j < values.length; j++) {
	    append(ret,separator);
	    append(ret,values[j]);
	}
	return ret;
    }

    /**
      * Returns the value from Node node with the name optionName.
      * The command name is provided for supplying error information.
      * @param node the Node whose name is optionName and whose value will 
      * be returned
      * @param commandName the Name of the command being processed
      * @param optionName the name of the option of the command being processed
      * @return the values for the option with name optionName as Strings 
      * (they may differ if Node node is in a scatter)
      */
    private String[] getOptionValue(Node node, String commandName,
				    String optionName) 
	throws MalformedScriptException {
	NamedNodeMap nnm = node.getAttributes();
	if (nnm != null && nnm.getLength() > 0) {
	    Node subNode = nnm.getNamedItem("substitution_string");
	    if (subNode != null) {
		String val = subNode.getNodeValue();
		if (val != null && !val.equals("")) {
		    String[] vals = expand("",numScripts);
		    // check against parameters
		    int subIndex = val.indexOf("${");
		    String substr = val;
		    while (subIndex != -1) {
			append(vals,substr.substring(0,subIndex));
			substr = substr.substring(subIndex+2);
			int endIndex = substr.indexOf("}");
			if (endIndex == -1) {
			    String err = "Referenced parameter " + substr 
				+ " missing '}'";
			    throw new MalformedScriptException(err);
							       
			}
			String paramName = substr.substring(0,endIndex);
			append(vals,resolve(paramName,commandName+"="
					    +optionName));
			substr = substr.substring(endIndex+1);
			subIndex = substr.indexOf("${");
		    }
		    append(vals,substr);
		    return vals;
		} else {
		    DEditDebug.println(this,3,">"+val+"<");
		}
	    } else {
		DEditDebug.println(this,3,"no substitution_string attr"
				   + " for " + optionName); 
	    }
	} else {
	    DEditDebug.println(this,3,"no attributes for " + optionName);
	}
	Node commandText = node.getFirstChild();
	if (commandText.getNodeType() != Node.TEXT_NODE) {
	    throw new MalformedScriptException("Command " + commandName
					       + ", option " + optionName +
					       " has no value.");
	}
	return expand(commandText.getNodeValue(),numScripts);
    }

    /**
      * Resolves a value found inside ${} in the script.  If it is inside 
      * a scatter as a parameter or is a normal parameter used in a scatter 
      * node (which won't appear in the script text), it replaces it with 
      * the real value.  Otherwise it just checks to make sure it is in 
      * the parameters block.
      * @param paramName the name of the parameter to resolve
      * @param commandOpt a String which specifies the source of the 
      * call to resolve()
      * @return an array of Strings which give the appropriate value for 
      * paramName in each of the scripts.
      */
    private String[] resolve(String paramName,String commandOpt) 
	throws MalformedScriptException {
	// bhpread=key->range
	Object value = parameterTable.get(paramName);
	if (value == null) {
	    String err = "No match for referenced parameter "+paramName+".";
	    throw new MalformedScriptException(err);
	}
	
	if (!scatters.empty()) {
	    ResolvedScatterParam chosen = null;
	    Stack holder = new Stack();
	    int divFact = 1;
	    int modFact = 1;
	    while (!scatters.empty()) {
		Scatter curr = (Scatter)scatters.pop();
		holder.push(curr);
		ResolvedScatterParam[] params = curr.params;
		for (int i = 0; i < params.length; i++) {
		    if (params[i].name.equals(paramName)) {
			chosen = params[i];
			break;
		    }
		    divFact *= params[i].getMulFact();
		}
		if (chosen != null) {
		    modFact = chosen.getMulFact();
		    break;
		}
	    }
	    while (!holder.empty()) {
		scatters.push(holder.pop());
	    }
	    if (chosen != null) {
		String[] ret = new String[numScripts];
		DEditDebug.println(this,5,"length " + ret.length + " divFact "
				   + divFact + " modFact " + modFact);
		for (int i = 0; i < ret.length; i++) {
		    int index = i / divFact;
		    index = index % modFact;
		    DEditDebug.println(this,5,i + " index " + index + ":"
				       + divFact + " " + modFact);
		    ret[i] = chosen.getValue(index, commandOpt);
		}
		return ret;
	    }
	}
	if (commandOpt.startsWith("scatter")) {
	    return new String[]{(String)parameterTable.get(paramName)};
	}
	return new String[]{"${"+paramName+"}"}; // or {(String)value)}; ???
    }

    /** 
      * This helper function takes a Node list and pares it down to just 
      * an array of Node's of type Node.ELEMENT_NODE.  This is extremely 
      * useful for getting rid of unwanted nodes between elements.
      * @param nl the NodeList to pare down to Element Node's
      * @return an array of Node's from NodeList which are Element Node's.
      */
    private static Node[] pareList(NodeList nl) {
	Vector elemNodes = new Vector();
	for (int i = 0; i < nl.getLength(); i++) {
	    Node temp = nl.item(i);
	    if (temp.getNodeType() == Node.ELEMENT_NODE) {
		elemNodes.add(temp);
	    }
	}
	Node[] ret = new Node[elemNodes.size()];
	elemNodes.toArray(ret);
	return ret;
    }

    /**
      * A useful helper function which performs chmod on the filename provided 
      * setting it to the permissions provided.  User beware: there is no 
      * checking on either provided String.
      * @param filename the name of the file to chmod
      * @param permissions new permissions for the file (you could put any 
      * other options in this String as well)
      * @return the return value of the Process.waitFor() where the Process is
      * generated from Runtime.exec("chmod " + permissions + " " + filename).  
      * In other words, should be the return value from system call chmod().
      */
    public int chmod(String filename, String permissions) 
	throws MalformedScriptException {
	try {
	    String command = "chmod " + permissions + " " + filename;
	    if (editor != null && editor.isApplet()) {
		return editor.executeRemote(command);
	    } 
	    String tempScriptName 
		= "reallylongandarbitraryandyoudneverthinktouseit.sh";
	    FileWriter fw = new FileWriter(tempScriptName);
	    fw.write(command);
	    fw.close();
	    Runtime rt = Runtime.getRuntime();
	    // chmod
	    
	    DEditDebug.println("",6,command);
	    Process p = rt.exec("sh " + tempScriptName);
	    int res = p.waitFor();
	    if (res != 0) {
		System.err.println("chmod error for " + filename 
				   + ": " + res);
	    }
	    p = rt.exec("rm " + tempScriptName);
	    res = p.waitFor();
	    if (res != 0) {
		System.err.println("error removing temporary script " 
				   + filename + ": " + res);
	    }
	    return res;
	} catch (IOException ioe) {
	    throw new MalformedScriptException("Couldn't write script: " 
					       + ioe.toString());
	} catch (InterruptedException ie) {
	    throw new MalformedScriptException("Command interrupted: " 
					       + ie.toString());
	}
    }
    
    /**
      * Writes a script file, remotely or locally depending on which 
      * mode the program is running (applet or app), with the given 
      * filename and file contents.
      * @param filename name of the file to write
      * @param contents contents of the file to write as a String
      */
    private void writeScriptFile(String filename, String contents) 
	throws MalformedScriptException {
	try {
	    if (editor != null && editor.isApplet()) {
		editor.writeRemoteFile(filename,contents);
		//chmod(filename,"775");
	    } else {
		FileWriter fw = new FileWriter(filename);
		fw.write(contents);
		fw.close();
		//chmod(filename,"775");
	    }
	} catch (Exception e) {
	    throw new MalformedScriptException("Error writing " + filename 
					       + ":" + e.toString());
	}
    }

    /**
      * Runs the ScriptGenerator with the provided arguments.
      * arg 1: name of XML file to generate script(s) from.
      * arg 2: directory to put scripts in
      * arg 3: base name for scripts
      * arg 4&5 (optional): -pl DEBUG_LEVEL
      * @param args the command line arguments
      */
    public static void main(String[] args) {
	String filename = null;
	String dirName = null;
	String scriptName = null;
	if (args.length > 2) {
	    filename = args[0];
	    dirName = args[1];
	    scriptName = args[2];
	}
	if (args.length > 4) {
	    if (args[3].equals("-pl")) {
		int debugLevel = Integer.parseInt(args[4]);
		DEditDebug.setPriority(debugLevel);
	    }
	}
	DefaultTreeModel m = new DefaultTreeModel(DEditNode.EMPTY_NODE); 
	if (filename != null) {
	    Document d = DEditor.createDOMforLocalXML(filename);
	    if (d != null) {
		ScriptGenerator sg = new ScriptGenerator(null);
		try {
		    Node n = d.getDocumentElement();
		    DEditNode newNode = new DEditNode(n,m);
		    int index = Math.max(0,filename.lastIndexOf("/"));
		    String xmlDir = filename.substring(0,index+1);
		    sg.generateScript(dirName,scriptName,xmlDir,n);
		    Vector v = sg.getWarnings();
		    if (v == null) {
			return;
		    }
		    for (int i = 0; i < v.size(); i++) {
			System.err.println(v.get(i));
		    }
		} catch (MalformedScriptException mse) {
		    System.err.println(mse.toString());
		}
	    }
	}
    }

}
