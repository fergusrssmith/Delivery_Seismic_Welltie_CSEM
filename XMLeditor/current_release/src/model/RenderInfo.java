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

import java.awt.Color;
import java.util.regex.Pattern;

/**
  * This class contains information about how to render a DEditNode 
  * in the JTree.
  */
public class RenderInfo {

    private String display; // specifies what to show

    /**
      * The bits to be used for font style, an | of some combination of 
      * Font.PLAIN (0), Font.BOLD (1), and Font.ITALIC (2).  Hence, 0-3.
      */
    public int styleBits;
    /**
      * The color to use to display the text.
      */
    public Color textColor;
    
    /**
      * Private constructor.  Use RenderInfo.makeRenderInfo(String).
      * @param display how to display text describing the node
      * @param style the style bits for this RenderInfo
      * @param tColor the Color for this RenderInfo
      */
    private RenderInfo(String display, int style, Color tColor) {
	this.styleBits = style;
	this.textColor = tColor;
	this.display = display;
    }
    
    /**
      * Combines an array of RenderInfo's into one RenderInfo
      * which can be used to display a DEditNode.
      * @param ris the RenderInfos to join into one
      * @return RenderInfo a combination of the RenderInfo's in ris
      */
    public static RenderInfo synthesize(RenderInfo[] ris) {
	int newRed = 0, newGreen = 0, newBlue = 0, newAlpha = 0;
	int styleB = 0;
	String display = "";
	boolean hasN = false, hasD = false, hasA = false, hasV = false;
	for (int i = 0; i < ris.length; i++) {
	    newRed += ris[i].textColor.getRed();
	    newGreen += ris[i].textColor.getGreen();
	    newBlue += ris[i].textColor.getBlue();
	    newAlpha += ris[i].textColor.getAlpha();
	    styleB = (styleB | ris[i].styleBits);
	    hasN = hasN || (ris[i].display.indexOf("n") != -1);
	    hasD = hasD || (ris[i].display.indexOf("d") != -1);
	    hasA = hasA || (ris[i].display.indexOf("a") != -1);
	    hasV = hasV || (ris[i].display.indexOf("v") != -1);
	}
	newRed = newRed / ris.length;
	newGreen = newGreen / ris.length;
	newBlue = newBlue / ris.length;
	newAlpha = newAlpha / ris.length;
	if (hasN) {
	    display += "n";
	}
	if (hasD) {
	    display += "d";
	}
	if (hasA) {
	    display += "a";
	}
	if (hasV) {
	    display += "v";
	}
	return new RenderInfo(display,styleB,
			      new Color(newRed,newGreen,newBlue,newAlpha));
    }
    
    /**
      * Returns a String to use to display a node given a bunch of 
      * info about the node.
      * @param displayName the result of DEditNode.getDisplayName()
      * @param nodeName the result of DEditNode.getNode().getNodeValue();
      * @param attributeVals a String representing DEditNode.getNode()'s 
      * attributes
      * @param nodeValue a String representing DEditNode.getNodes()'s children 
      * of type TEXT_NODE
      * @return a String to use to display a node given a bunch of 
      * info about the node
      */
    public String getDisplayString(String displayName, String nodeName,
				   String attributeVals, String nodeValue) {
	String ret = "";
	if (display.indexOf("d") != -1) {
	    ret += displayName;
	} else if (display.indexOf("n") != -1) {
	    ret += nodeName;
	} 
	if (display.indexOf("v") != -1 && !nodeValue.equals("")) {
	    ret += ":" + nodeValue;
	}
	if (display.indexOf("a") != -1) {
	    ret += " > " + attributeVals; 
	}
	return ret;
    }

    /**
      * Creates a RenderInfo from a string.  It defaults to plain black
      * text with a shown value if there is some error.  See README for 
      * how to construct the string.
      * @param str the String to make the RenderInfo from
      * @return a RenderInfo constructed from a string.
      */
    public static RenderInfo makeRenderInfo(String str) {
	Pattern p = Pattern.compile("[.]");
	String display = "";
	int bits = 0, r = 0, g = 0, b = 0, a = 255;
	String[] paramStrs = p.split(str);
	if (paramStrs.length >= 6) {
	    display = paramStrs[0];
	    bits = Integer.parseInt(paramStrs[1]);
	    r = Integer.parseInt(paramStrs[2]);
	    g = Integer.parseInt(paramStrs[3]);
	    b = Integer.parseInt(paramStrs[4]);
	    a = Integer.parseInt(paramStrs[5]);
	} 
	Color color = new Color(r,g,b,a);
	return new RenderInfo(display,bits,color);
    }
}
