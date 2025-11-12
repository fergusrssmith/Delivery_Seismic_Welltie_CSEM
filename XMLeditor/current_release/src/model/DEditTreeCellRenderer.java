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

import org.w3c.dom.Node;
import org.w3c.dom.Attr;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.NodeList;
import java.awt.Component;
import java.awt.Font;
import java.awt.Color;
import javax.swing.tree.*;
import javax.swing.JTree;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Pattern;

/**
  * Class which renders the tree specially if need be.  This should
  * eventually base it's methods in part on a configuration file or something
  * for a given xsd. 
  */
public class DEditTreeCellRenderer extends DefaultTreeCellRenderer {

    // Singleton for a basic "no-frills" RenderInfo to use if no other 
    // ones are specified
    private static final RenderInfo BASIC 
	= RenderInfo.makeRenderInfo("dv.0.0.0.0.255");
    /**
      * Stores the RenderInfos keyed by "attribute_name=attribute_value".
      */
    private Hashtable attHt; 
    /**
      * Stores the Icons keyed by node name
      */
    private Hashtable iconHt;
    // the plugin used by the DEditor that this DEditTreeCellRenderer renders
    // trees for
    private DEditPlugin plugin = null;

    /**
      * Creates a DEditTreeCellRenderer for use with a JTree with DEditNodes
      * in it.  Constructs RenderInfos for customized rendering.
      * @param configPairs an array of arrays of 2 strings, of the format
      * "attribute_name=attribute_value" and "disp-type.int.int.int.int.int". 
      * See README for more information.
      * @param iconPairs like configPairs but store tag name and icon file 
      * location
      */
    public DEditTreeCellRenderer(String[][] configPairs,String[][] iconPairs) {
	attHt = new Hashtable();
	iconHt = new Hashtable();
	if (configPairs != null ) {
	    for (int i = 0; i < configPairs.length; i++) {
		if (configPairs[i] != null && configPairs[i].length == 2) {
		    attHt.put(configPairs[i][0],
			      RenderInfo.makeRenderInfo(configPairs[i][1]));
		}
	    }
	} 
	if (iconPairs != null) {
	    ClassLoader cl = getClass().getClassLoader();
	    Pattern p = Pattern.compile("[,]");
	    for (int i = 0; i < iconPairs.length; i++) {
	       if (iconPairs[i] != null) {
		   String iconName = iconPairs[i][0];
		   DEditDebug.println(this,3,iconName + " " + iconPairs[i][1]);
		   ImageIcon icon 
		       = new ImageIcon(cl.getResource(iconName));
		   if (icon != null) {
		       String nodeNameList = iconPairs[i][1];
		       String[] nodeNameStrs = p.split(nodeNameList);
		       for (int j = 0; j < nodeNameStrs.length; j++) {
			   iconHt.put(nodeNameStrs[j],icon);
		       }
		   } 
	       }
	   }
	}
    }
    
    /**
      * Major function in class, uses inputs to render component of tree
      * associated with Object value.  See TreeCellRenderer javadoc for 
      * meaning of parameters.
      */
    public Component getTreeCellRendererComponent(JTree tree, 
						  Object value, 
						  boolean selected, 
						  boolean expanded,
						  boolean leaf,
						  int row,
						  boolean hasFocus) {
	
	super.getTreeCellRendererComponent(tree,value,selected,expanded,
					   leaf,row,hasFocus);
	this.setFont(this.getFont().deriveFont(Font.PLAIN));
	this.setForeground(Color.BLACK);
	if (!(value instanceof DEditNode)) {
	    return this;
	}
	DEditNode denode = (DEditNode)value;
	if (denode.getNode() == null) {
	    return this;
	}
	Icon special = (Icon)iconHt.get(denode.getNode().getNodeName());
	if (special != null) {
	    DEditDebug.println(this,2,"->" + special.getIconHeight() + " " 
			       + special.getIconWidth());
	    setIcon(special);
	}
	String dispName = denode.getDisplayName();
	Node node = denode.getNode();
	String nodeName = node.getNodeName();
	String nodeValue = "";
	String attributeValues = "";
	if (leaf) {
	    // look for text based children, which indicates that a special
	    // rendering method should be used
	    NodeList nlist = node.getChildNodes();
	    for (int i = 0; i < nlist.getLength(); i++) {
		Node child = nlist.item(i);
		if (child.getNodeType() == Node.TEXT_NODE) {
		    String childVal = child.getNodeValue();
		    if (!childVal.trim().equals("")) {
			nodeValue += " " + childVal;
		    }
		}
	    }
	}
	Vector renderInfos = new Vector();
	// check attributes and treat some specially
	NamedNodeMap nnm = node.getAttributes();
	if(nnm != null) {
	    int len = nnm.getLength();
	    Attr attr;
	    for (int j = 0; j < len; j++) {
		attr = (Attr)nnm.item(j);
		int newStyle = 0;
		String aname = attr.getNodeName();
		String avalue = attr.getNodeValue();
		if (!avalue.equals("")) {
		    RenderInfo ri = (RenderInfo)attHt.get(aname+"="+avalue);
		    if (ri != null) {
			DEditDebug.println(this,1,aname+"="+avalue);
			renderInfos.add(ri);
		    }
		    ri = (RenderInfo)attHt.get(aname+"=");
		    if (ri != null) {
			renderInfos.add(ri);
			attributeValues += " " + avalue;
			DEditDebug.println(this,3,avalue);
		    }
		}
	    }
	}
	if (renderInfos.size() == 0) {
	    renderInfos.add(BASIC);
	}
	RenderInfo[] synth = new RenderInfo[renderInfos.size()];
	renderInfos.toArray(synth);
	RenderInfo toRender = RenderInfo.synthesize(synth);
	this.setForeground(toRender.textColor);
	this.setFont(this.getFont().deriveFont(toRender.styleBits));
	String forLabel = toRender.getDisplayString(dispName,nodeName,
						    attributeValues,nodeValue);
	this.setText(forLabel);
	Component ret = this;
	//DEditDebug.println(this,6,this.getText());
	if (plugin != null) {
	    ret = plugin.getTreeCellRendererComponent(this,tree,denode,
						      selected,expanded,
						      leaf,row,hasFocus);
	}
	//DEditDebug.println(this,6,this.getText());
	return ret;
    }

    /**
      * Sets the plugin that this DEditTreeCellRenderer uses to augment its
      * display of certain DEditNodes.
      * @param plugin the plugin to use to augment display
      */
    public void setPlugin(DEditPlugin plugin) {
	this.plugin = plugin;
    }
}

