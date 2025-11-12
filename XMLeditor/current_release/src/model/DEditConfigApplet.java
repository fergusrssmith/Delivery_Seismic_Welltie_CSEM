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
import javax.swing.event.*;
import java.awt.event.*;
import java.lang.reflect.Method;
import java.util.Hashtable;
import java.awt.GridLayout;

/**
  * Class which can be used to run the DEditor as an applet.  It is a work
  * in progress which will allow total configurability (although a PARAM in 
  * the html may restrict what can be configured).
  */
public class DEditConfigApplet extends DEditAppletDirect {

    private JButton startButton;
    private Hashtable params;
    private final JComboBox presets;
    private static final String bhpBase = "/home/ahglim/nprg/xmledit/beta";
    private static final String wtrBase = "/home/williamr/shigoto/beta";
    private String base;

    /** 
      * Creates the applet with it's GUI components.
      */
    public DEditConfigApplet() {
	super();
	this.base = bhpBase;
	this.params = new Hashtable();
	this.startButton = new JButton("Launch XML Editor");
	startButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    String selected = (String)presets.getSelectedItem();
		    setParametersBy(selected);
		    startEditor();
		}
	    });
	String[] presetNames = new String[]{"wavelet_id",
					    "amplitude_extraction",
					    "delivery_editor",
					    "script_editor"};
	this.presets = new JComboBox(presetNames);
	this.presets.addItemListener(new ItemListener() {
		public void itemStateChanged(ItemEvent ie) {
		    String selected = (String)presets.getSelectedItem();
		    setParametersBy(selected);
		}
	    });
	JPanel appPanel = new JPanel(new GridLayout(2,1));
	appPanel.add(presets);
	appPanel.add(startButton);
	this.getContentPane().add(appPanel);
	this.setVisible(true);
    }
  
    /** 
      * Initializes the applet by grabbing the base directory for files.
      */ 
    public void init() {
	String tempBase = getParameter("remoteBaseDirectory");
	if (tempBase != null) {
	    this.base = tempBase;
	}
    }

    /**
      * Overrides JApplet method so that it can set parameters through
      * the applet window.
      * @param name the name of the parameter to get
      * @return the value for that paramter name
      */
    public String getParameter(String name) {
	//String defVal = super.getParameter(name);
	if (params.containsKey(name)) {
	    return (String)params.get(name);
	} else {
	    return super.getParameter(name);// defVal;
	} 
    }

    /**
      * Sets the value of a parameter so that future calls to getParameter()
      * will use this value.
      * @param name the name of the parameter to set
      * @param value the value to set param "name" to
      */
    public void setParameter(String name, String value) {
	params.put(name,value);
    } 

    /**
      * Sets all major patterns according to some presets.  There should
      * be a cleaner way of doing this than the current code.  Note that
      * the preset names should match those in the JComboBox.
      * @param preset the name of a preset group of parameters.  Current 
      * valid choices include "script_editor","delivery_editor",
      * "amplitude_extraction", and "wavelet_id".
      */
    public void setParametersBy(String preset) {
	if (preset.equals("script_editor")) {
	    setParameter("configFile",base+"/configFiles/scriptconfig.xml");
	    setParameter("schema_filename","xsd/script.xsd");
	    setParameter("filename","");
	    setParameter("allow_open_file","true");
	    setParameter("allow_new_file","true");
	    setParameter("plugin_flags","");
	    setParameter("plugin","script.ScriptPlugin");
	} else if (preset.equals("delivery_editor")) {
	    setParameter("configFile",base+"/configFiles/deliveryConfig.xml");
	    setParameter("schema_filename","xsd/Delivery.xsd");
	    setParameter("filename","");
	    setParameter("allow_open_file","true");
	    setParameter("allow_new_file","true");
	    setParameter("plugin_flags","");
	    setParameter("plugin","delivery.DeliveryPlugin");
	} else if (preset.equals("amplitude_extraction")) {
	    setParameter("configFile",base+"/configFiles/scriptconfig.xml");
	    setParameter("schema_filename","xsd/script.xsd");
	    setParameter("filename",base+"/sampleXML/east_bhpintegrate.xml");
	    setParameter("allow_open_file","false");
	    setParameter("allow_new_file","false");
	    setParameter("plugin_flags","-r script-hs="+base+"/plugins/help/script/amplitude_extraction/amplitude_extraction.hs");
	    setParameter("plugin","script.ScriptPlugin");
	} else if (preset.equals("wavelet_id")) {
setParameter("configFile",base+"/configFiles/scriptconfig.xml");
	    setParameter("schema_filename","xsd/script.xsd");
	    setParameter("filename",base+"/sampleXML/wavelet2.xml");
	    setParameter("allow_open_file","false");
	    setParameter("allow_new_file","false");
	    setParameter("plugin_flags","-r script-hs="+base
			 +"/plugins/help/script/wavelet_id/wavelet_id.hs");
	    setParameter("plugin","script.ScriptPlugin");
	} else {
	    // no match, no setting
	}
    }
}
