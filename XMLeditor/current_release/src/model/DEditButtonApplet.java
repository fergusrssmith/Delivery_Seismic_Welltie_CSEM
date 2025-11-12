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
  * Class which can be used to run the DEditor as an applet.  It creates a 
  * button which launches the editor.  It subclasses DEditAppletDirect so
  * that all the configuration properties can be set up there with
  * startEditor().
  */
public class DEditButtonApplet extends DEditAppletDirect {

    private JButton startButton;
    /** 
      * Creates the applet with it's button.
      */
    public DEditButtonApplet() {
	super();
	startButton = new JButton("");
	startButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    startEditor();
		}
	    });
	this.getContentPane().add(startButton);
	this.setVisible(true);	
    }
    
    /**
      * Initializes the applet by filling in the text for the button.
      */
    public void init() {
	String buttonText = getParameter("button_text");
	if (buttonText == null) buttonText = "Launch XML Editor...";
	startButton.setText(buttonText);
    }
}

