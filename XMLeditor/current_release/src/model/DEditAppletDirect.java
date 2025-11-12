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
package model;

import javax.swing.*;
import java.awt.event.*;
import java.lang.reflect.Method;

/**
  * Class which can be used to run the DEditor as an applet.  It automatically
  * launches the editor when the applet is init()ed.  To set parameters and/or 
  * not launch automatically, create a subclass and override the constructor
  * and the function init().  Be sure that getParameter() will return the 
  * values you want.
  */
public class DEditAppletDirect extends JApplet {

    // default configuration parameters 
    private String curDirString = "";
    private String configFile 
	= "/home/ahglim/nprg/xmledit/beta/configFiles/scriptconfig.xml";
    private String filename = null;
    private String schemaFilename 
	= "/home/ahglim/nprg/xmledit/beta/xsd/script.xsd";
    private boolean allowNewFile = true;
    private boolean allowOpenFile = true;
    private boolean allowLoadSchema = false;
    private boolean relaxed = false;
    private DEditPlugin plugin;
    private String pluginName = null;
    
    /** 
      * Creates the applet and makes it visible.
      */
    public DEditAppletDirect() {
	super();
	this.setVisible(true);
	
    }
    
    /**
      * Initializes the applet by starting the editor.
      */
    public void init() {
	startEditor();
    }

    /**
      * Now that all the parameters are gathered, creates the editor and 
      * starts it.
      */
    protected void startEditor() {
	// get a file to open if there is one specified
	filename = getParameter("filename");
	if (filename != null && filename.equals("")) {
	    filename = null;
	}
	// build a plugin if there is one provided
	pluginName = getParameter("plugin");
	if (pluginName != null) {
	    try {
		String flags = getParameter("plugin_flags"); 
		String[] leftOvers;
		if (flags == null) {
		    leftOvers = new String[]{};
		} else {
		    leftOvers = flags.split(" ");
		}
		Class c = Class.forName(pluginName);
		if (c != null) {
		    String[] tempArr = new String[0];
		    Method m=c.getDeclaredMethod("buildPlugin",
						 new Class[]{String[].class});
		    if (m != null) {
			plugin =(DEditPlugin)m.invoke(null, 
						      new Object[]{leftOvers});
		    }
		}
	    } catch (Exception e) {
		System.err.println("load exception error: " + e.toString());
		e.printStackTrace();
	    }
	}
	// get all the other parmeters
	String temp = getParameter("configFile");
	if (temp != null) configFile = temp;
	if (configFile != null && configFile.equals("")) {
	    configFile = null;
	}
	temp = getParameter("schema_filename");
	if (temp != null) schemaFilename = temp;
	temp = getParameter("allow_new_file");
	if (temp != null) allowNewFile = Boolean.valueOf(temp).booleanValue();
	temp = getParameter("allow_open_file");
	if (temp != null) allowOpenFile = Boolean.valueOf(temp).booleanValue();
	temp = getParameter("relaxed");
	if (temp != null) relaxed = Boolean.valueOf(temp).booleanValue();
	System.out.println("Schema file name: " + schemaFilename);
	DEditor editor = new DEditor(curDirString,configFile,filename,
				     schemaFilename,allowNewFile,allowOpenFile,
				     allowLoadSchema,relaxed,plugin,this);
	editor.start();
    }
}
